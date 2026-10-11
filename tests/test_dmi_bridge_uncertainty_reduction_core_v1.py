import math
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import dmi_bridge_uncertainty_reduction_core_v1 as core


def direct_summary(delta, weights, truth, direction, lam):
    signed = np.asarray(delta, dtype=float)
    w = np.asarray(weights, dtype=float)
    w = w / w.sum()
    sign = np.where(signed > core.TIE_TOL, 1, np.where(signed < -core.TIE_TOL, -1, 0))
    if direction is not None:
        w = w * np.exp(lam * direction * sign)
        w = w / w.sum()
    mag = np.abs(signed)
    mean = float(np.dot(w, mag))
    return mean, float(np.dot(w, (mag - mean) ** 2)), w


def test_variance_mixture_identity_matches_direct_reweighting_and_support():
    delta = np.array([-5.0, -1.0, 0.0, 0.5, 3.0, 9.0])
    weights = np.array([0.10, 0.15, 0.2, 0.1, 0.25, 0.2])
    truth = 2.0
    result = core.analyze_distribution(delta, weights, truth)
    for lam in core.LAMBDA_GRID:
        for cue, direction in (("correct", 1), ("incorrect", -1)):
            expected_mean, expected_var, expected_weights = direct_summary(delta, weights, truth, direction, lam)
            got = result["updates"][format(lam, ".17g")][cue]
            assert got["mean"] == pytest.approx(expected_mean, abs=1e-14)
            assert got["variance"] == pytest.approx(expected_var, abs=1e-13)
            assert got["support_count_unchanged"] is True
            assert np.array_equal(expected_weights > 0, np.asarray(weights) > 0)


def test_zero_lambda_exactly_reuses_baseline_uncertainty_and_crps():
    out = core.analyze_distribution([-4, -2, 0, 1, 3], [0.1, 0.2, 0.1, 0.3, 0.3], 2.0)
    baseline = out["baseline"]
    for cue in ("correct", "incorrect"):
        update = out["updates"]["0"][cue]
        for metric in ("mean", "variance", "sd", "iqr", "central90_width", "crps"):
            assert update[metric] == baseline[metric]
        assert update["absolute_variance_reduction"] == 0.0
        assert update["relative_variance_reduction"] == 0.0


def test_truth_ties_have_no_correct_or_incorrect_cue():
    out = core.analyze_distribution([-1, 0, 1], [0.3, 0.4, 0.3], 0.5e-12)
    assert out["truth_direction"] == 0
    assert all(not rows for rows in out["updates"].values())


def test_weighted_inverse_cdf_quantiles_match_pl1_convention():
    values = np.array([0.0, 1.0, 2.0, 3.0])
    weights = np.array([0.2, 0.3, 0.1, 0.4])
    assert core.weighted_quantile_sorted(values, weights, 0.0) == 0.0
    assert core.weighted_quantile_sorted(values, weights, 0.25) == 1.0
    assert core.weighted_quantile_sorted(values, weights, 0.5) == 1.0
    assert core.weighted_quantile_sorted(values, weights, 0.6) == 2.0


def test_relative_variance_threshold_excludes_degenerate_case_but_keeps_absolute_change():
    out = core.analyze_distribution([2.0, 2.0, 2.0], [0.2, 0.3, 0.5], 2.0)
    assert out["relative_variance_eligible"] is False
    assert math.isnan(out["updates"]["0.25"]["correct"]["relative_variance_reduction"])
    assert out["updates"]["0.25"]["correct"]["absolute_variance_reduction"] == 0.0


@pytest.mark.parametrize("seed", [2, 19, 404])
def test_weighted_discrete_crps_matches_quadratic_reference(seed):
    rng = np.random.default_rng(seed)
    values = np.sort(rng.uniform(0, 5, size=11))
    weights = rng.uniform(0.01, 1.0, size=11)
    weights /= weights.sum()
    truth = 1.37
    quadratic = float(np.dot(weights, np.abs(values - truth)))
    quadratic -= 0.5 * float(np.sum(weights[:, None] * weights[None, :] * np.abs(values[:, None] - values[None, :])))
    assert core.crps_sorted(values, weights, truth) == pytest.approx(quadratic, abs=1e-14)


def test_crps_rejects_unsorted_values():
    with pytest.raises(ValueError, match="sorted"):
        core.crps_sorted([1, 0], [0.5, 0.5], 0.5)


def test_cue_crps_reused_mixture_matches_direct():
    rng = np.random.default_rng(17)
    delta = rng.normal(size=41)
    weights = rng.random(41)
    truth = 0.8
    result = core.analyze_distribution(delta, weights, truth)
    signs = np.where(delta > 1e-12, 1, np.where(delta < -1e-12, -1, 0))
    for lam in (0.0, 0.25, 0.5, 1.0):
        for cue in ("correct", "incorrect"):
            direction = result["truth_direction"] * (1 if cue == "correct" else -1)
            posterior = weights * np.exp(lam * direction * signs)
            posterior /= posterior.sum()
            order = np.argsort(np.abs(delta), kind="mergesort")
            direct = core.crps_sorted(np.abs(delta)[order], posterior[order], abs(truth))
            assert result["updates"][format(lam, ".17g")][cue]["crps"] == pytest.approx(direct, abs=1e-12)


def test_a1_a2_key_match_is_exact_and_rejects_duplicates_or_mismatch():
    assert core.exact_key_match([("p2", "R2"), ("p1", "R1")], [("p1", "R1"), ("p2", "R2")]) == [("p1", "R1"), ("p2", "R2")]
    with pytest.raises(ValueError, match="unique"):
        core.exact_key_match([("p1",), ("p1",)], [("p1",), ("p2",)])
    with pytest.raises(ValueError, match="populations"):
        core.exact_key_match([("p1",)], [("p2",)])


def test_crossed_mouse_bootstrap_preserves_product_multiplicity():
    # Seed produces a reproducible draw with duplicated identities in both cohorts.
    rng = np.random.default_rng(1)
    mult = core.crossed_mouse_multiplicities(["c1", "c2", "c3"], ["g1", "g2", "g3"], rng)
    assert sum(mult.values()) == 9
    assert all(value > 0 for value in mult.values())
    # Multiplicity is the product of independently sampled cohort counts.
    assert all(value in {1, 2, 3, 4, 6, 9} for value in mult.values())


def _model_rows(n=40):
    rng = np.random.default_rng(10)
    rows = []
    for i in range(n):
        eta = i / n
        entropy = rng.normal()
        correct = 0.4 * eta + 0.1 * entropy
        incorrect = -0.2 * eta + 0.05 * entropy
        rows.append({
            "reaction_id": f"R{i:03d}",
            "analysis_split": "DEVELOPMENT" if i < 30 else "CONFIRMATION_HOLDOUT",
            "sign_magnitude_eta2": eta,
            "directional_entropy3": entropy,
            "dominant_sign_mass": 0.5 + eta / 2,
            "non_tie_coverage": 0.8,
            "n_supported_sign_states": 2 + i % 2,
            "correct_relative_variance_reduction": correct,
            "incorrect_relative_variance_reduction": incorrect,
            "observed_usefulness": float(i % 2),
        })
    return rows


def test_correct_and_incorrect_prospective_targets_are_independently_reported():
    rows = _model_rows()
    correct = core.fit_prospective_models(rows, target="correct_relative_variance_reduction")
    wrong = core.fit_prospective_models(rows, target="incorrect_relative_variance_reduction")
    assert {row["target"] for row in correct} == {"correct_relative_variance_reduction"}
    assert {row["target"] for row in wrong} == {"incorrect_relative_variance_reduction"}
    assert [row["model"] for row in correct] == [row["model"] for row in wrong]
    assert all(row["n_confirmation"] == 10 for row in correct + wrong)
    assert all("usefulness" not in row["predictors"] for row in correct + wrong)


def test_prospective_models_ignore_truth_derived_usefulness_columns():
    rows = _model_rows()
    baseline = core.fit_prospective_models(rows, target="correct_relative_variance_reduction")
    changed = [dict(row, observed_usefulness=10000.0 - i) for i, row in enumerate(rows)]
    actual = core.fit_prospective_models(changed, target="correct_relative_variance_reduction")
    assert [row["mae"] for row in actual] == pytest.approx([row["mae"] for row in baseline])


def test_prospective_rejects_pooled_or_unknown_cue_target():
    with pytest.raises(ValueError, match="distinguish"):
        core.fit_prospective_models(_model_rows(), target="pooled_relative_variance_reduction")


def test_explanatory_partial_association_is_separate_and_finite():
    rows = _model_rows()
    result = core.explanatory_partial_rank_association(
        [row["correct_relative_variance_reduction"] for row in rows],
        [row["observed_usefulness"] for row in rows],
        [row["sign_magnitude_eta2"] for row in rows],
    )
    assert result["n"] == len(rows)
    assert math.isfinite(result["rho_eta2_usefulness_adjusted"])


def test_declared_prediction_combinations_are_baseline_descriptors_only():
    specs = core.prospective_model_specs()
    assert specs["eta2"] == ("sign_magnitude_eta2",)
    assert specs["eta2_plus_entropy"] == ("sign_magnitude_eta2", "directional_entropy3")
    assert set(specs["all_geometry"]) == {
        "sign_magnitude_eta2", "directional_entropy3", "dominant_sign_mass",
        "non_tie_coverage", "n_supported_sign_states",
    }
    assert all("usefulness" not in col for cols in specs.values() for col in cols)
