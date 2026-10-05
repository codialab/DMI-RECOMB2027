import importlib.util
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

MODULE_PATH = Path(__file__).parents[1] / "scripts" / "dmi_bridge_pl2d_hypothesis_freeze_core_v1.py"
spec = importlib.util.spec_from_file_location("pl2dh", MODULE_PATH)
m = importlib.util.module_from_spec(spec)
assert spec.loader
spec.loader.exec_module(m)


def test_primary_contract_is_single_feature_reaction_level_and_fixed():
    assert m.PRIMARY["unit"] == "reaction"
    assert m.PRIMARY["predictor"] == "mean_sign_magnitude_eta2_across_evaluations"
    assert m.PRIMARY["response"] == "mean_directionally_useful_fraction_across_defined_evaluations"
    assert m.PRIMARY["minimum_finite_reactions"] == 700
    assert m.PRIMARY["minimum_spearman_rho"] == 0.50
    assert m.PRIMARY["minimum_positive_contexts"] == 12
    assert m.BOOTSTRAP["replicates"] == 5000
    assert m.BOOTSTRAP["seed"] == 20260929


def test_development_snapshot_is_frozen():
    s = m.DEVELOPMENT_SNAPSHOT
    assert s["finite_primary_reactions"] == 2455
    assert s["primary_reaction_level_rho"] == pytest.approx(0.8570601823118907, abs=1e-15)
    assert s["primary_partial_rho_controlling_non_tie_pair_fraction"] == pytest.approx(0.8356721975298953, abs=1e-15)
    assert s["eta2_vs_mean_information_advantage_rho"] == pytest.approx(0.3730769973633655, abs=1e-15)
    assert s["positive_algorithm_rna_contexts"] == 16


def _split_rows():
    rows = []
    for i in range(4180):
        rows.append({
            "reaction_id": f"R{i:04d}",
            "split_hash": f"{i:064x}",
            "split_rank": str(i),
            "analysis_split": "CONFIRMATION_HOLDOUT" if i < 1045 else "DEVELOPMENT",
            "sign_magnitude_eta2": "POISON",
            "directionally_useful_fraction": "POISON",
        })
    return rows


def test_confirmation_registry_is_identity_only_and_exact():
    out = m.build_confirmation_registry(_split_rows())
    assert len(out) == 1045
    assert [int(r["split_rank"]) for r in out] == list(range(1045))
    assert all(set(r) == {"reaction_id", "split_hash", "split_rank", "analysis_split"} for r in out)
    assert all("POISON" not in str(r) for r in out)


def test_confirmation_registry_rejects_split_count_change():
    rows = _split_rows()
    rows[-1]["analysis_split"] = "CONFIRMATION_HOLDOUT"
    with pytest.raises(ValueError):
        m.build_confirmation_registry(rows)


def test_authoritative_numerical_reference_is_frozen():
    assert "pandas.DataFrame.groupby" in m.NUMERICAL_REFERENCE["reaction_mean"]
    assert "scipy.stats.spearmanr" in m.NUMERICAL_REFERENCE["spearman"]
    assert "three scipy.stats.spearmanr" in m.NUMERICAL_REFERENCE["partial_spearman"]


def test_reaction_aggregation_uses_pandas_groupby_mean_not_fsum():
    rows = [
        {"reaction_id": "R1", "sign_magnitude_eta2": 1e16, "directionally_useful_fraction": 0.0, "non_tie_pair_fraction": 0.0, "mean_information_advantage": 0.0, "directional_entropy3": 0.0, "dominant_sign_mass": 1.0, "n_supported_sign_states": 1.0},
        {"reaction_id": "R1", "sign_magnitude_eta2": 1.0, "directionally_useful_fraction": 0.5, "non_tie_pair_fraction": 0.5, "mean_information_advantage": 1.0, "directional_entropy3": 0.5, "dominant_sign_mass": 0.5, "n_supported_sign_states": 2.0},
        {"reaction_id": "R1", "sign_magnitude_eta2": -1e16, "directionally_useful_fraction": 1.0, "non_tie_pair_fraction": 1.0, "mean_information_advantage": 2.0, "directional_entropy3": 1.0, "dominant_sign_mass": 0.0, "n_supported_sign_states": 3.0},
    ]
    out = m.aggregate_reaction_rows(rows)
    assert len(out) == 1
    # pandas groupby.mean gives 0.0 here; math.fsum(values)/3 would give 1/3.
    assert out[0]["mean_sign_magnitude_eta2"] == 0.0
    expected = pd.DataFrame.from_records(rows).groupby("reaction_id", sort=True)["sign_magnitude_eta2"].mean().iloc[0]
    assert out[0]["mean_sign_magnitude_eta2"] == float(expected)


def test_aggregate_predictor_does_not_condition_on_response_availability():
    rows = [
        {"reaction_id": "R1", "sign_magnitude_eta2": "1", "directionally_useful_fraction": "", "non_tie_pair_fraction": "0", "mean_information_advantage": "", "directional_entropy3": "", "dominant_sign_mass": "", "n_supported_sign_states": ""},
        {"reaction_id": "R1", "sign_magnitude_eta2": "3", "directionally_useful_fraction": "0.5", "non_tie_pair_fraction": "0.5", "mean_information_advantage": "2", "directional_entropy3": "0.5", "dominant_sign_mass": "0.5", "n_supported_sign_states": "2"},
        {"reaction_id": "R1", "sign_magnitude_eta2": "5", "directionally_useful_fraction": "1", "non_tie_pair_fraction": "1", "mean_information_advantage": "4", "directional_entropy3": "1", "dominant_sign_mass": "0", "n_supported_sign_states": "3"},
    ]
    out = m.aggregate_reaction_rows(rows)
    assert len(out) == 1
    r = out[0]
    assert r["mean_sign_magnitude_eta2"] == pytest.approx(3.0)
    assert r["mean_directionally_useful_fraction"] == pytest.approx(0.75)
    assert r["mean_non_tie_pair_fraction"] == pytest.approx(0.5)
    assert r["mean_information_advantage"] == pytest.approx(3.0)
    assert r["n_finite_eta2_evaluations"] == 3
    assert r["n_defined_response_evaluations"] == 2


def test_spearman_handles_ties_deterministically():
    assert m.spearman([1, 1, 2, 3], [1, 2, 3, 4]) == pytest.approx(0.9486832980505138)
    assert m.spearman([1, 2, 3], [3, 2, 1]) == pytest.approx(-1.0)


def test_spearman_matches_scipy_reference_exactly():
    x = [0.0, 0.0, 1.0, 2.0, 2.0, 4.0, 7.0]
    y = [3.0, 1.0, 1.0, 2.0, 2.0, 8.0, 5.0]
    assert m.spearman(x, y) == float(spearmanr(x, y).statistic)


def test_partial_spearman_one_control_is_finite_and_directional():
    x = [1, 2, 3, 4, 5, 6]
    y = [1, 3, 2, 5, 4, 6]
    z = [1, 1, 2, 2, 3, 3]
    r = m.partial_spearman_one_control(x, y, z)
    assert math.isfinite(r)
    assert r > 0


def test_primary_rho_uses_only_finite_reaction_pairs():
    rows = [
        {"mean_sign_magnitude_eta2": 1.0, "mean_directionally_useful_fraction": 1.0},
        {"mean_sign_magnitude_eta2": 2.0, "mean_directionally_useful_fraction": 2.0},
        {"mean_sign_magnitude_eta2": None, "mean_directionally_useful_fraction": 0.0},
        {"mean_sign_magnitude_eta2": 3.0, "mean_directionally_useful_fraction": math.nan},
    ]
    n, rho = m.primary_rho(rows)
    assert n == 2
    assert rho == pytest.approx(1.0)


def test_confirmation_status_support_effect_and_context_gates():
    assert m.classify_confirmation(finite_reactions=699, overall_rho=0.9, evaluable_contexts=16, positive_contexts=16) == "PL2D_INSUFFICIENT_SUPPORT"
    assert m.classify_confirmation(finite_reactions=800, overall_rho=0.49, evaluable_contexts=16, positive_contexts=16) == "PL2D_NOT_CONFIRMED"
    assert m.classify_confirmation(finite_reactions=800, overall_rho=0.8, evaluable_contexts=15, positive_contexts=15) == "PL2D_INSUFFICIENT_SUPPORT"
    assert m.classify_confirmation(finite_reactions=800, overall_rho=0.8, evaluable_contexts=16, positive_contexts=11) == "PL2D_NOT_CONFIRMED"
    assert m.classify_confirmation(finite_reactions=800, overall_rho=0.8, evaluable_contexts=16, positive_contexts=12) == "PL2D_CONFIRMED"


def test_supportive_features_are_not_primary():
    assert set(m.SUPPORTIVE_FEATURES) == {"directional_entropy3", "dominant_sign_mass", "n_supported_sign_states"}
    assert m.SUPPORTIVE_FEATURES["dominant_sign_mass"] == "NEGATIVE"
    assert "directional_entropy3" not in m.PRIMARY["predictor"]


def test_pinned_pl2c_identity_and_split():
    assert m.EXPECTED_PL2C_MANIFEST_SHA256 == "6abcf149cd6adb294c0931585cec28f992b453451a5b5f03e10b3ffb2f9b9c39"
    assert m.EXPECTED_SPLIT_SHA256 == "7dd8676508416a13f88af184110ef40e0f8d92782cd5d62b42e7006af5316e2f"
    assert m.EXPECTED_CONFIRMATION_REGISTRY_SHA256 == "02a38afe3cb25b8c85a7774c1f20118ad9601ba74f0df0fbff021ba2494bd2e5"
    assert m.EXPECTED_SPLIT_COUNTS == {"DEVELOPMENT": 3135, "CONFIRMATION_HOLDOUT": 1045}
    assert len(m.EXPECTED_PL2C_ARTIFACT_SHA256) == 9
