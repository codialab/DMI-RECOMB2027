import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('candidate_sensitivity', ROOT / 'analyses/reviewer_checks/candidate_sensitivity.py')
cs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cs)


def test_vectorized_descriptors_match_frozen_pl1_core():
    sys.path.insert(0, '/home/pty/work/project_MRI_brain_tumor/scripts')
    import dmi_bridge_pl1_predictability_core_v1 as pl1
    left = np.array([[0.0, 1.0, 3.0], [1.0, 2.0, 2.5], [0.0, -1.0, 2.0]])
    right = np.array([[0.0, 2.0, 1.0], [1.0, 1.0, 0.5]])
    wl = np.array([0.2, 0.3, 0.5]); wr = np.array([0.7, 0.3])
    got = cs.descriptor_block(left, right, wl, wr)
    for j in range(left.shape[1]):
        expected = pl1.target_landscape_metrics(
            anchor_ct2a=left[:, j], anchor_gl261=right[:, j],
            target_ct2a=left[:, j], target_gl261=right[:, j],
            strong_weights_ct2a=wl, strong_weights_gl261=wr)
        assert got['eta2'][j] == pytest.approx(expected['sign_magnitude_eta2'], abs=1e-12)
        assert got['H_dir'][j] == pytest.approx(expected['directional_entropy3'], abs=1e-12)
        assert got['D_dir'][j] == pytest.approx(expected['dominant_sign_mass'], abs=1e-12)


def test_vectorized_weak_cue_gains_match_frozen_scalar_core():
    sys.path.insert(0, '/home/pty/work/project_MRI_brain_tumor/scripts')
    import dmi_bridge_pl2_sign_only_core_v1 as pl2
    left = np.array([[1.0, 3.0], [2.0, 1.0], [4.0, 2.0]])
    right = np.array([[0.0, 1.0], [2.0, 0.5]])
    wl = np.array([0.2, 0.3, 0.5]); wr = np.array([0.6, 0.4])
    truth = np.array([1.5, -2.0])
    got = cs.vector_outcomes(left, right, wl, wr, truth)
    for j in range(2):
        delta = (left[:, j, None] - right[None, :, j]).reshape(-1)
        joint = (wl[:, None] * wr[None, :]).reshape(-1)
        expected = pl2.evaluate_truth_case(delta_b=delta, baseline_weights=joint,
                                           truth_delta_b=truth[j], lambda_grid=(0.25,))
        row = expected['lambda_results'][0]
        assert got['correct_gain'][j] == pytest.approx(row['correct_abs_error_gain'], abs=1e-12)
        assert got['wrong_gain'][j] == pytest.approx(row['wrong_abs_error_gain'], abs=1e-12)


def test_truth_exclusion_operates_on_whole_vectors_and_renormalizes():
    weights = np.array([0.1, 0.2, 0.3, 0.4])
    sample_ids = np.array([10, 11, 12, 13])
    keep = sample_ids != 12
    remaining = weights[keep] / weights[keep].sum()
    assert sample_ids[keep].tolist() == [10, 11, 13]
    assert remaining.sum() == pytest.approx(1.0)
    assert remaining.tolist() == pytest.approx([1/7, 2/7, 4/7])


def test_global_pool_temperature_operator_reproduces_target_ess():
    sys.path.insert(0, '/home/pty/work/project_MRI_brain_tumor/scripts')
    import dmi_bridge_a20_dual_anchor_core_v1 as a20
    distance = np.linspace(0.0, 0.2, 320) ** 2
    for target in (10.0, 20.0, 40.0):
        result = a20.solve_temperature_for_ess(distance, target_ess=target)
        assert result.achieved_ess == pytest.approx(target, abs=1e-9)
        assert result.weights.sum() == pytest.approx(1.0, abs=1e-14)
