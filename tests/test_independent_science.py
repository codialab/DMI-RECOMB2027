import numpy as np
import pytest

from analyses.reviewer_checks import independent_science as science


def test_sign_only_posterior_tilts_directions_and_keeps_ties_neutral():
    delta = np.array([-2.0, 0.0, 1.0])
    weights = np.array([0.2, 0.3, 0.5])
    correct = science.posterior_magnitude(delta, weights, 1)
    incorrect = science.posterior_magnitude(delta, weights, -1)
    factors = np.exp(np.array([-0.25, 0.0, 0.25]))
    expected = np.sum(weights * factors * np.abs(delta)) / np.sum(weights * factors)
    assert correct == pytest.approx(expected, abs=1e-15)
    assert incorrect == pytest.approx(
        np.sum(weights * factors[::-1] * np.abs(delta)) / np.sum(weights * factors[::-1]), abs=1e-15
    )
    assert correct != incorrect


def test_product_distribution_is_full_cartesian_with_condition_weights():
    delta, weights = science.product_delta_distribution(
        [1.0, 3.0], [0.0, 2.0], [0.25, 0.75], [0.6, 0.4]
    )
    assert delta.tolist() == [1.0, -1.0, 3.0, 1.0]
    assert weights.tolist() == pytest.approx([0.15, 0.10, 0.45, 0.30])
    assert weights.sum() == pytest.approx(1.0)


def test_level2a_rederives_error_gain_and_usefulness_from_predictions():
    row = science.endpoint_from_predictions(2.0, 1.0, 2.5, 0.5)
    assert row["baseline_error"] == pytest.approx(1.0)
    assert row["correct_gain"] == pytest.approx(0.5)
    assert row["incorrect_gain"] == pytest.approx(-0.5)
    assert row["information_advantage"] == pytest.approx(0.5)
    assert row["directionally_useful"] is True
    tie = science.endpoint_from_predictions(1e-12, 0.5, 0.4, 0.6)
    assert tie["truth_status"] == "TRUTH_TIE"
    assert np.isnan(tie["correct_gain"])


def test_heldout_exclusion_removes_an_entire_state_and_renormalizes():
    vectors = np.array([[1.0, 11.0], [2.0, 22.0], [3.0, 33.0], [4.0, 44.0]])
    kept, weights = science.exclude_whole_vector(vectors, [10, 11, 12, 13], 12,
                                                  [0.1, 0.2, 0.3, 0.4])
    assert kept.tolist() == [[1.0, 11.0], [2.0, 22.0], [4.0, 44.0]]
    assert weights.tolist() == pytest.approx([1 / 7, 2 / 7, 4 / 7])


def test_heldout_exclusion_rejects_ambiguous_identity():
    with pytest.raises(ValueError, match="must be unique"):
        science.exclude_whole_vector(np.ones((3, 2)), [1, 1, 2], 1, [0.2, 0.3, 0.5])


def test_level2b_reconstructs_predictions_after_per_side_whole_vector_exclusion():
    left = np.array([[0.0, 3.0], [1.0, 2.0], [4.0, 1.0]])
    right = np.array([[0.0, 1.0], [2.0, 0.0], [3.0, 2.0]])
    got = science.level2b_predictions(
        left, right, [0, 1, 2], [0, 1, 2], 1, 0,
        [0.2, 0.3, 0.5], [0.6, 0.1, 0.3], [2.0, -1.0]
    )
    kept_left, wl = science.exclude_whole_vector(left, [0, 1, 2], 1, [0.2, 0.3, 0.5])
    kept_right, wr = science.exclude_whole_vector(right, [0, 1, 2], 0, [0.6, 0.1, 0.3])
    for j in range(2):
        delta, joint = science.product_delta_distribution(kept_left[:, j], kept_right[:, j], wl, wr)
        assert got["baseline"][j] == pytest.approx(np.dot(joint, np.abs(delta)), abs=1e-15)
        cue = 1 if [2.0, -1.0][j] > 0 else -1
        assert got["correct"][j] == pytest.approx(science.posterior_magnitude(delta, joint, cue), abs=1e-15)
        assert got["incorrect"][j] == pytest.approx(science.posterior_magnitude(delta, joint, -cue), abs=1e-15)


def test_eta2_direction_conditioned_variance_including_ties_and_degenerate_case():
    delta = np.array([-3.0, -1.0, 0.0, 0.0, 2.0])
    weights = np.array([0.1, 0.2, 0.15, 0.15, 0.4])
    eta = science.eta2_directional_magnitude(delta, weights)
    assert 0 <= eta <= 1
    assert np.isnan(science.eta2_directional_magnitude(np.zeros(5), weights))


def test_eta2_geometry_block_uses_full_candidate_support():
    left = np.array([[0.0, 0.0], [2.0, 2.0]])
    right = np.array([[0.0, 1.0], [1.0, 0.0]])
    wl, wr = [0.3, 0.7], [0.4, 0.6]
    got = science.eta2_geometry_block(left, right, wl, wr)
    expected = []
    for j in range(2):
        d, w = science.product_delta_distribution(left[:, j], right[:, j], wl, wr)
        expected.append(science.eta2_directional_magnitude(d, w))
    assert got == pytest.approx(expected, abs=1e-15)


def test_geometry_parallel_matches_serial_with_multiple_batches():
    rng = np.random.default_rng(17)
    mats = rng.normal(size=(4, 20, 37))
    evals = np.array([[0, 1], [1, 2], [2, 3], [3, 0], [0, 2], [1, 3]])
    weights = rng.random((len(evals), 4, 20))
    serial = science.reconstruct_geometry_parallel(mats, evals, weights, np.arange(31), workers=1, batch_size=2)
    parallel = science.reconstruct_geometry_parallel(mats, evals, weights, np.arange(31), workers=2, batch_size=2)
    assert [x[0] for x in serial] == [x[0] for x in parallel] == list(range(len(evals)))
    for left, right in zip(serial, parallel):
        np.testing.assert_array_equal(left[1], right[1])
        np.testing.assert_array_equal(left[2], right[2])


def test_geometry_worker_exception_is_not_silently_accepted():
    mats = np.ones((2, 20, 3))
    evals = np.array([[0, 1]])
    weights = np.ones((1, 4, 20))
    with pytest.raises(RuntimeError, match="geometry worker pool failed"):
        science.reconstruct_geometry_parallel(mats, evals, weights, np.array([4]), workers=2, batch_size=1)


def test_case_aggregation_uses_distinct_non_tie_truths_and_pair_weighting():
    rows = [
        {"anchor_setting":"A1", "evaluation_id":"e1", "reaction_id":"r1", "truth_pair_id":"p1",
         "truth_delta":2., "baseline_prediction":1., "correct_prediction":2., "incorrect_prediction":.5},
        {"anchor_setting":"A1", "evaluation_id":"e1", "reaction_id":"r1", "truth_pair_id":"p2",
         "truth_delta":-1., "baseline_prediction":1., "correct_prediction":.5, "incorrect_prediction":2.},
        {"anchor_setting":"A1", "evaluation_id":"e1", "reaction_id":"r1", "truth_pair_id":"pt",
         "truth_delta":0., "baseline_prediction":1., "correct_prediction":1., "incorrect_prediction":1.},
        {"anchor_setting":"A1", "evaluation_id":"e1", "reaction_id":"r2", "truth_pair_id":"p3",
         "truth_delta":1., "baseline_prediction":.5, "correct_prediction":1., "incorrect_prediction":0.},
    ]
    cases, pairs = science.aggregate_case_predictions(rows)
    assert len(cases) == 4
    assert len(pairs) == 2
    p1 = pairs.set_index("reaction_id").loc["r1"]
    assert p1.n_distinct_non_tie_truths == 2
    assert p1.correct_gain_mean == pytest.approx(.25)
    assert p1.incorrect_gain_mean == pytest.approx(-.75)
    assert p1.directional_usefulness_fraction == pytest.approx(.5)


def test_case_aggregation_rejects_duplicate_and_missing_identities():
    row = {"anchor_setting":"A1", "evaluation_id":"e", "reaction_id":"r", "truth_pair_id":"p",
           "truth_delta":1., "baseline_prediction":1., "correct_prediction":1., "incorrect_prediction":1.}
    with pytest.raises(ValueError, match="unique"):
        science.aggregate_case_predictions([row, row])
    row["truth_pair_id"] = None
    with pytest.raises(ValueError, match="complete and unique"):
        science.aggregate_case_predictions([row])


def test_subnormal_weights_and_zero_lambda_are_stable():
    d = np.array([-2., 1., 3.])
    w = np.array([np.nextafter(0., 1.), .5, .5])
    assert science.posterior_magnitude(d, w, 1, lam=0.) == pytest.approx(2.)
    assert np.isfinite(science.posterior_magnitude(d, w, 1, lam=.25))
    assert np.isfinite(science.posterior_magnitude(d, [0., .5, .5], -1, lam=.25))


def test_geometry_serial_and_parallel_identity_with_nan_reactions():
    mats = np.zeros((2, 20, 4))
    evals = np.array([[0, 1], [1, 0]])
    weights = np.full((2, 4, 20), 1/20)
    rxns = np.arange(4)
    serial = science.reconstruct_geometry_parallel(mats, evals, weights, rxns, workers=1)
    parallel = science.reconstruct_geometry_parallel(mats, evals, weights, rxns, workers=2)
    for a, b in zip(serial, parallel):
        assert a[0] == b[0]
        np.testing.assert_array_equal(np.isnan(a[1]), np.isnan(b[1]))
        np.testing.assert_array_equal(np.isnan(a[2]), np.isnan(b[2]))
