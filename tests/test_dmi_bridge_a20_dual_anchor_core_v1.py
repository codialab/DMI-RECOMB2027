import numpy as np

from scripts.dmi_bridge_a20_dual_anchor_core_v1 import (
    average_ranks,
    candidate_rank_coordinate,
    dual_squared_distance,
    effective_sample_size,
    kernel_weights,
    mouse_rank_coordinate,
    solve_temperature_for_ess,
)


def test_average_ranks_ties():
    x = [10.0, 20.0, 20.0, 40.0]
    got = average_ranks(x)
    np.testing.assert_allclose(got, [1.0, 2.5, 2.5, 4.0])


def test_rank_coordinate_conventions():
    np.testing.assert_allclose(
        candidate_rank_coordinate([1.0, 2.0, 3.0, 4.0]),
        [0.125, 0.375, 0.625, 0.875],
    )
    np.testing.assert_allclose(
        mouse_rank_coordinate([1.0, 2.0, 3.0, 4.0, 5.0]),
        [0.0, 0.25, 0.5, 0.75, 1.0],
    )


def test_dual_distance_is_sum_of_squared_coordinate_residuals():
    q1 = np.array([0.0, 0.5, 1.0])
    q2 = np.array([1.0, 0.5, 0.0])
    got = dual_squared_distance(q1, q2, 0.5, 0.5)
    np.testing.assert_allclose(got, [0.5, 0.0, 0.5])


def test_kernel_is_normalized_and_ess_is_bounded():
    w = kernel_weights([0.0, 0.1, 0.2, 0.3], 0.1)
    assert np.isclose(w.sum(), 1.0)
    assert 1.0 <= effective_sample_size(w) <= 4.0


def test_temperature_solver_hits_target_ess():
    d = np.linspace(0.0, 1.0, 100) ** 2
    out = solve_temperature_for_ess(d, target_ess=20.0)
    assert np.isclose(out.weights.sum(), 1.0)
    assert abs(out.achieved_ess - 20.0) < 1e-8
    assert 1e-12 < out.temperature <= 1.0


def test_temperature_solver_expands_archived_upper_bracket():
    d = np.linspace(0.0, 100.0, 100) ** 2
    out = solve_temperature_for_ess(d, target_ess=20.0)
    assert np.isclose(out.weights.sum(), 1.0)
    assert abs(out.achieved_ess - 20.0) < 1e-8
    assert 1.0 < out.temperature <= 1e6


def test_temperature_solver_returns_final_hi_not_geometric_mean():
    d = np.linspace(0.0, 1.0, 100) ** 2
    out = solve_temperature_for_ess(d, target_ess=20.0, iterations=7)

    lo = 1e-12
    hi = 1.0
    for _ in range(7):
        mid = np.sqrt(lo * hi)
        ess_mid = effective_sample_size(kernel_weights(d, mid))
        if ess_mid < 20.0:
            lo = mid
        else:
            hi = mid

    assert out.temperature == hi
    np.testing.assert_allclose(out.weights, kernel_weights(d, hi), rtol=0, atol=0)
