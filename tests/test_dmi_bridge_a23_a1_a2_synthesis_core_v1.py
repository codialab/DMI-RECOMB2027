import numpy as np
import pytest

from scripts.dmi_bridge_a23_a1_a2_synthesis_core_v1 import (
    average_ranks,
    directional_specificity,
    paired_delta,
    require_exact_key_match,
    spearman_rho,
)


def test_exact_key_match_sorts_and_accepts_reordered_input():
    left = [("p2", "R2"), ("p1", "R1")]
    right = [("p1", "R1"), ("p2", "R2")]
    assert require_exact_key_match(left, right) == [("p1", "R1"), ("p2", "R2")]


def test_exact_key_match_rejects_duplicates():
    with pytest.raises(ValueError):
        require_exact_key_match([("p1",), ("p1",)], [("p1",), ("p2",)])


def test_exact_key_match_rejects_population_difference():
    with pytest.raises(ValueError):
        require_exact_key_match([("p1",)], [("p2",)])


def test_paired_delta_is_a2_minus_a1():
    np.testing.assert_allclose(
        paired_delta([1.0, 4.0, -2.0], [1.5, 3.0, 1.0]),
        [0.5, -1.0, 3.0],
    )


def test_average_ranks_uses_average_ties():
    np.testing.assert_allclose(average_ranks([3.0, 1.0, 1.0, 7.0]), [3.0, 1.5, 1.5, 4.0])


def test_spearman_matches_perfect_monotone_relationship():
    assert np.isclose(spearman_rho([1, 2, 3, 4], [10, 20, 30, 40]), 1.0)
    assert np.isclose(spearman_rho([1, 2, 3, 4], [40, 30, 20, 10]), -1.0)


def test_directional_specificity_is_correct_minus_wrong():
    np.testing.assert_allclose(
        directional_specificity([0.3, -0.1, 0.0], [0.1, -0.4, 0.2]),
        [0.2, 0.3, -0.2],
    )
