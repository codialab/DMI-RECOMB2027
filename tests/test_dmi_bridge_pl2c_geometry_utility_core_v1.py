import math

import numpy as np
import pytest

from scripts import dmi_bridge_pl2c_geometry_utility_core_v1 as core


def test_reaction_split_is_order_invariant_and_exact_quarter():
    ids = [f"R{i:04d}" for i in range(20)]
    first = core.freeze_reaction_split(ids)
    second = core.freeze_reaction_split(list(reversed(ids)))
    assert first == second
    assert sum(row["analysis_split"] == "CONFIRMATION_HOLDOUT" for row in first) == 5
    assert sum(row["analysis_split"] == "DEVELOPMENT" for row in first) == 15
    assert core.split_registry_sha256(first) == core.split_registry_sha256(second)


def test_frozen_4180_split_counts_are_exact():
    ids = [f"R{i:05d}" for i in range(core.EXPECTED_REACTIONS)]
    rows = core.freeze_reaction_split(ids)
    assert sum(row["analysis_split"] == "CONFIRMATION_HOLDOUT" for row in rows) == 1045
    assert sum(row["analysis_split"] == "DEVELOPMENT" for row in rows) == 3135


def test_reaction_split_rejects_duplicate_ids():
    with pytest.raises(ValueError, match="unique"):
        core.freeze_reaction_split(["R1", "R1", "R2", "R3"])


def test_information_advantage_identity_and_directionally_useful_status():
    result = core.case_responses(
        correct_gain=4e-6,
        wrong_gain=-2e-6,
        random_gain=1e-6,
        truth_status="NON_TIE",
    )
    assert result["information_advantage"] == pytest.approx(3e-6)
    assert result["information_advantage"] == pytest.approx(0.5 * (4e-6 - (-2e-6)))
    assert result["correct_gain_status"] == "POSITIVE"
    assert result["information_advantage_status"] == "POSITIVE"
    assert result["directionally_useful"] is True


def test_correct_gain_alone_is_not_enough_for_directionally_useful():
    result = core.case_responses(
        correct_gain=2e-6,
        wrong_gain=4e-6,
        random_gain=3e-6,
        truth_status="NON_TIE",
    )
    assert result["correct_gain_status"] == "POSITIVE"
    assert result["information_advantage_status"] == "NEGATIVE"
    assert result["directionally_useful"] is False


def test_random_sign_identity_fails_closed():
    with pytest.raises(ValueError, match="random-sign"):
        core.case_responses(
            correct_gain=2e-6,
            wrong_gain=-2e-6,
            random_gain=1e-3,
            truth_status="NON_TIE",
        )


def test_truth_tie_has_no_directional_response():
    result = core.case_responses(
        correct_gain=math.nan,
        wrong_gain=math.nan,
        random_gain=math.nan,
        truth_status="TRUTH_TIE",
    )
    assert result["response_status"] == "TRUTH_TIE"
    assert math.isnan(result["information_advantage"])
    assert result["directionally_useful"] is False


def test_distinct_pair_aggregation_is_equal_weight_and_retains_ties():
    rows = [
        {"truth_pair_id": "P1", "truth_status": "NON_TIE", "correct_gain": 4e-6, "wrong_gain": -2e-6, "random_gain": 1e-6},
        {"truth_pair_id": "P2", "truth_status": "NON_TIE", "correct_gain": -1e-6, "wrong_gain": -3e-6, "random_gain": -2e-6},
        {"truth_pair_id": "P3", "truth_status": "TRUTH_TIE", "correct_gain": math.nan, "wrong_gain": math.nan, "random_gain": math.nan},
    ]
    out = core.aggregate_distinct_pair_responses(rows)
    assert out["n_distinct_pairs"] == 3
    assert out["n_non_tie_pairs"] == 2
    assert out["n_truth_tie_pairs"] == 1
    assert out["non_tie_pair_fraction"] == pytest.approx(2.0 / 3.0)
    assert out["mean_correct_gain"] == pytest.approx(1.5e-6)
    assert out["mean_information_advantage"] == pytest.approx(2e-6)
    assert out["positive_correct_gain_fraction"] == pytest.approx(0.5)
    assert out["positive_information_advantage_fraction"] == pytest.approx(1.0)
    assert out["directionally_useful_fraction"] == pytest.approx(0.5)


def test_aggregation_rejects_duplicate_truth_pair_ids():
    rows = [
        {"truth_pair_id": "P1", "truth_status": "TRUTH_TIE"},
        {"truth_pair_id": "P1", "truth_status": "TRUTH_TIE"},
    ]
    with pytest.raises(ValueError, match="unique"):
        core.aggregate_distinct_pair_responses(rows)


def test_all_truth_ties_are_retained_but_directional_response_is_undefined():
    out = core.aggregate_distinct_pair_responses([
        {"truth_pair_id": "P1", "truth_status": "TRUTH_TIE"},
        {"truth_pair_id": "P2", "truth_status": "TRUTH_TIE"},
    ])
    assert out["aggregate_status"] == "NO_DIRECTIONAL_TRUTH"
    assert out["n_truth_tie_pairs"] == 2
    assert out["non_tie_pair_fraction"] == 0.0
    assert math.isnan(out["directionally_useful_fraction"])


def test_spearman_uses_average_ranks_and_finite_pairs_only():
    result = core.spearman_finite([1, 2, 2, 4, math.nan], [10, 20, 20, 40, 50])
    assert result["finite_count"] == 4
    assert result["nonfinite_count"] == 1
    assert result["rho"] == pytest.approx(1.0)


def test_linear_quartiles_and_labels_are_frozen():
    edges = core.linear_quartile_edges([0, 1, 2, 3, 4])
    assert edges == pytest.approx((1.0, 2.0, 3.0))
    assert core.quartile_label(0.5, edges) == "Q1"
    assert core.quartile_label(1.5, edges) == "Q2"
    assert core.quartile_label(2.5, edges) == "Q3"
    assert core.quartile_label(3.5, edges) == "Q4"
    assert core.quartile_label(math.nan, edges) == "NONFINITE"


def test_contract_pins_sources_features_and_confirmation_boundary():
    contract = core.analysis_contract()
    assert contract["source_identity"]["pl2b_manifest_sha256"] == core.PL2B_MANIFEST_SHA256
    assert contract["source_identity"]["pl2b_case_logical_sha256"] == core.PL2B_CASE_LOGICAL_SHA256
    assert contract["reaction_split"]["expected_development"] == 3135
    assert contract["reaction_split"]["expected_confirmation"] == 1045
    assert contract["reaction_split"]["confirmation_outcomes_must_not_be_analyzed_in_pl2c"] is True
    assert contract["primary_features"] == list(core.PRIMARY_FEATURES)
    assert contract["pl2d_confirmation_status_on_success"] == "UNTOUCHED"
    assert "confirmation_outcome_analysis" in contract["forbidden"]
