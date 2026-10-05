import math
import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import dmi_bridge_pl2_sign_only_core_v1 as core


def test_joint_delta_distribution_cartesian_not_index_paired():
    d, w = core.joint_delta_distribution([3.0, 5.0], [1.0, 4.0], [0.25, 0.75], [0.6, 0.4])
    np.testing.assert_allclose(d, [2.0, -1.0, 4.0, 1.0])
    np.testing.assert_allclose(w, [0.15, 0.10, 0.45, 0.30])
    assert math.isclose(float(w.sum()), 1.0)


def test_sign_group_statistics_partition_and_baseline():
    d = np.array([-3.0, -1.0, 0.0, 2.0])
    w = np.array([0.1, 0.2, 0.3, 0.4])
    s = core.sign_group_statistics(d, w)
    assert s.negative_mass == pytest.approx(0.3)
    assert s.tie_mass == pytest.approx(0.3)
    assert s.positive_mass == pytest.approx(0.4)
    assert s.baseline_abs_mean == pytest.approx(1.3)
    assert s.baseline_signed_mean == pytest.approx(0.3)


def test_lambda_zero_is_exact_baseline_for_correct_and_wrong():
    d = np.array([-4.0, -2.0, 1.0, 3.0])
    w = np.array([0.1, 0.2, 0.3, 0.4])
    s = core.sign_group_statistics(d, w)
    p1 = core.posterior_from_sign_groups(s, observed_direction=1, lambda_total=0.0)
    pm = core.posterior_from_sign_groups(s, observed_direction=-1, lambda_total=0.0)
    assert p1["magnitude_mean"] == pytest.approx(s.baseline_abs_mean)
    assert pm["magnitude_mean"] == pytest.approx(s.baseline_abs_mean)
    assert p1["signed_mean"] == pytest.approx(s.baseline_signed_mean)
    assert pm["signed_mean"] == pytest.approx(s.baseline_signed_mean)
    assert p1["joint_ess"] == pytest.approx(s.baseline_ess)
    assert pm["joint_ess"] == pytest.approx(s.baseline_ess)


def test_correct_direction_can_improve_magnitude_without_using_truth_magnitude():
    # Positive branch magnitude is near the truth; negative branch is much larger.
    d = np.array([-10.0, -8.0, 2.0, 3.0])
    w = np.full(4, 0.25)
    out = core.evaluate_truth_case(
        delta_b=d,
        baseline_weights=w,
        truth_delta_b=2.5,
        lambda_grid=[1.0],
        reliability_grid=[0.0, 0.5, 1.0],
    )
    row = out["lambda_results"][0]
    assert row["correct_abs_error_gain"] > 0.0
    assert row["wrong_abs_error_gain"] < 0.0
    assert row["correct"]["favorable_mass"] > out["baseline_positive_mass"]


def test_random_sign_is_exact_q_half_expected_gain():
    d = [-7.0, -2.0, 1.0, 4.0]
    out = core.evaluate_truth_case(
        delta_b=d,
        baseline_weights=[0.2, 0.3, 0.1, 0.4],
        truth_delta_b=3.0,
        lambda_grid=[0.5],
    )
    row = out["lambda_results"][0]
    expected = 0.5 * (row["correct_abs_error_gain"] + row["wrong_abs_error_gain"])
    assert row["random_sign_expected_abs_error_gain"] == pytest.approx(expected)


def test_reliability_endpoints_match_wrong_and_correct():
    assert core.reliability_expected_gain(correct_gain=2.0, wrong_gain=-1.0, reliability=0.0) == -1.0
    assert core.reliability_expected_gain(correct_gain=2.0, wrong_gain=-1.0, reliability=1.0) == 2.0
    assert core.reliability_expected_gain(correct_gain=2.0, wrong_gain=-1.0, reliability=0.5) == 0.5
    assert core.break_even_reliability(correct_gain=2.0, wrong_gain=-1.0) == pytest.approx(1.0 / 3.0)


def test_break_even_nan_when_no_unique_crossing_or_outside_unit_interval():
    assert math.isnan(core.break_even_reliability(correct_gain=1.0, wrong_gain=1.0))
    assert math.isnan(core.break_even_reliability(correct_gain=2.0, wrong_gain=1.0))


def test_truth_tie_is_retained_without_invented_direction():
    out = core.evaluate_truth_case(
        delta_b=[-1.0, 0.0, 1.0],
        baseline_weights=[0.25, 0.5, 0.25],
        truth_delta_b=0.0,
    )
    assert out["status"] == "TRUTH_TIE"
    assert out["truth_direction"] == 0
    assert out["lambda_results"] == []


def test_posterior_matches_direct_reweighting():
    d = np.array([-4.0, -1.0, 0.0, 2.0, 8.0])
    w = core.normalize_weights([0.05, 0.15, 0.2, 0.25, 0.35])
    lam = 0.7
    s = core.sign_group_statistics(d, w)
    got = core.posterior_from_sign_groups(s, observed_direction=1, lambda_total=lam)
    score = np.where(d > core.TIE_TOL, 1.0, np.where(d < -core.TIE_TOL, -1.0, 0.0))
    direct = w * np.exp(lam * score)
    direct /= direct.sum()
    assert got["magnitude_mean"] == pytest.approx(float(np.dot(direct, np.abs(d))))
    assert got["signed_mean"] == pytest.approx(float(np.dot(direct, d)))
    assert got["joint_ess"] == pytest.approx(float(1.0 / np.sum(direct * direct)))


def test_contract_separates_pl2_from_historical_tanh_operator():
    c = core.analysis_contract()
    assert c["weak_information"] == "binary_direction_only"
    assert c["primary_target"] == "absolute_between_condition_difference_magnitude_E_abs_delta_B"
    assert c["primary_lambda"] == 0.25
    assert c["truth_selection"]["reaction_B_values_must_not_select_truth_candidates"] is True
    assert "historical_tanh_operator_relabeling" in c["forbidden"]


def test_gain_classification_is_frozen_and_symmetric():
    assert core.classify_gain(2e-12) == "IMPROVED"
    assert core.classify_gain(-2e-12) == "HARMED"
    assert core.classify_gain(0.5e-12) == "TIED"
    assert core.classify_gain(-0.5e-12) == "TIED"
    assert core.classify_gain(float("nan")) == "UNDEFINED"


def test_contract_uses_distinct_truth_pairs_and_compact_reliability_storage():
    c = core.analysis_contract()
    assert c["gain_classification"]["tolerance"] == 1e-12
    assert c["primary_benchmark_unit"] == "distinct_held_out_truth_pair_once"
    assert "do_not_duplicate_weight_primary_summaries" in c["quantile_alias_policy"]
    assert c["reliability_storage"] == "derive_from_correct_and_wrong_endpoints_without_materializing_q_rows"
