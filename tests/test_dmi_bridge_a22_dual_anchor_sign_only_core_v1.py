import math

import numpy as np
import pytest

from scripts import dmi_bridge_pl2_sign_only_core_v1 as pl2
from scripts.dmi_bridge_a22_dual_anchor_sign_only_core_v1 import (
    EXPECTED_CASE_ROWS,
    analysis_contract,
    common_target_reactions,
    evaluate_truth_case,
    frozen_pl2_semantics,
    remove_holdout_and_normalize,
    validate_frozen_pl2_semantics,
)


def test_frozen_pl2_semantics_are_exact():
    validate_frozen_pl2_semantics()
    observed = frozen_pl2_semantics()
    assert observed["tie_tol"] == 1e-12
    assert observed["gain_tol"] == 1e-12
    assert observed["lambda_grid"] == (0.0, 0.25, 0.5, 1.0)
    assert observed["primary_lambda"] == 0.25
    assert observed["reliability_grid"] == (0.0, 0.25, 0.5, 0.75, 1.0)


def test_common_target_reactions_excludes_both_strong_anchors():
    reactions = [f"R{i}" for i in range(4179)] + ["HEX1", "LDH_L"]
    targets = common_target_reactions(reactions)
    assert len(targets) == 4179
    assert "HEX1" not in targets
    assert "LDH_L" not in targets


def test_holdout_removal_renormalizes_without_new_threshold():
    out = remove_holdout_and_normalize([0.05, 0.15, 0.30, 0.50], 2)
    np.testing.assert_allclose(out.weights, [0.05 / 0.70, 0.15 / 0.70, 0.50 / 0.70])
    assert math.isclose(out.post_holdout_mass, 0.70)
    assert out.positive_weight_candidates_pre == 4
    assert out.positive_weight_candidates_post == 3
    assert out.post_holdout_ess > 1.0


def test_holdout_zero_remaining_mass_fails_closed():
    with pytest.raises(ValueError, match="zero post-holdout mass"):
        remove_holdout_and_normalize([0.0, 1.0, 0.0], 1)


def test_a22_truth_case_delegates_exactly_to_frozen_pl2_core():
    kwargs = {
        "delta_b": [-3.0, -1.0, 0.0, 2.0, 4.0],
        "baseline_weights": [0.1, 0.2, 0.1, 0.3, 0.3],
        "truth_delta_b": 2.5,
    }
    direct = pl2.evaluate_truth_case(**kwargs)
    wrapped = evaluate_truth_case(**kwargs)

    def canon(value):
        if isinstance(value, float) and math.isnan(value):
            return "NaN"
        if isinstance(value, dict):
            return {key: canon(item) for key, item in value.items()}
        if isinstance(value, list):
            return [canon(item) for item in value]
        return value

    assert canon(wrapped) == canon(direct)


def test_random_sign_control_is_exact_endpoint_average():
    result = evaluate_truth_case(
        delta_b=[-2.0, -1.0, 1.0, 4.0],
        baseline_weights=[0.1, 0.2, 0.3, 0.4],
        truth_delta_b=3.0,
    )
    primary = next(row for row in result["lambda_results"] if row["lambda_total"] == 0.25)
    expected = 0.5 * (
        primary["correct_abs_error_gain"] + primary["wrong_abs_error_gain"]
    )
    assert primary["random_sign_expected_abs_error_gain"] == expected


def test_contract_freezes_matched_population_and_expected_case_count():
    contract = analysis_contract()
    assert contract["qualified_baseline"] == "A2-L"
    assert contract["matched_truth_pair_count"] == 836
    assert contract["target_reaction_count"] == 4179
    assert contract["expected_case_rows_if_support_gate_passes"] == 3_493_644
    assert EXPECTED_CASE_ROWS == 3_493_644
    assert contract["new_ess_threshold"] is False
    assert contract["a2_g_used"] is False
    assert contract["geometry_utility_analysis"] is False
