"""Numerical core for DMI-BRIDGE-PL2 sign-only truth-referenced utility.

This module contains no repository I/O.  It evaluates a weak sign-only prior on a
fixed finite candidate-pair distribution using sufficient statistics for the
negative/tie/positive branches.  The primary target is the magnitude
``|Delta B|`` rather than the signed contrast.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Sequence

import numpy as np

SCHEMA = "bridge.pl2.sign_only_utility.core.v1"
TIE_TOL = 1e-12
LAMBDA_GRID = (0.0, 0.25, 0.5, 1.0)
RELIABILITY_GRID = (0.0, 0.25, 0.5, 0.75, 1.0)
PRIMARY_LAMBDA = 0.25
GAIN_TOL = 1e-12


def _finite_1d(values: Sequence[float], *, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    if arr.size == 0:
        raise ValueError(f"{name} must be non-empty")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    return arr


def normalize_weights(weights: Sequence[float]) -> np.ndarray:
    w = _finite_1d(weights, name="weights")
    if np.any(w < 0.0):
        raise ValueError("weights must be non-negative")
    total = float(np.sum(w))
    if not math.isfinite(total) or total <= 0.0:
        raise ValueError("weights must have positive finite mass")
    return w / total


def effective_sample_size(weights: Sequence[float]) -> float:
    w = normalize_weights(weights)
    return float(1.0 / np.sum(w * w))


def joint_delta_distribution(
    ct2a_values: Sequence[float],
    gl261_values: Sequence[float],
    ct2a_weights: Sequence[float],
    gl261_weights: Sequence[float],
) -> tuple[np.ndarray, np.ndarray]:
    """Return Cartesian ``CT2A - GL261`` deltas and product weights."""
    left = _finite_1d(ct2a_values, name="ct2a_values")
    right = _finite_1d(gl261_values, name="gl261_values")
    wl = normalize_weights(ct2a_weights)
    wr = normalize_weights(gl261_weights)
    if left.size != wl.size or right.size != wr.size:
        raise ValueError("value and weight lengths must match within each condition")
    delta = (left[:, None] - right[None, :]).reshape(-1)
    joint = (wl[:, None] * wr[None, :]).reshape(-1)
    return delta, joint


@dataclass(frozen=True)
class SignGroupStatistics:
    negative_mass: float
    tie_mass: float
    positive_mass: float
    negative_abs_sum: float
    tie_abs_sum: float
    positive_abs_sum: float
    negative_signed_sum: float
    tie_signed_sum: float
    positive_signed_sum: float
    negative_weight_sq_sum: float
    tie_weight_sq_sum: float
    positive_weight_sq_sum: float

    @property
    def total_mass(self) -> float:
        return self.negative_mass + self.tie_mass + self.positive_mass

    @property
    def baseline_abs_mean(self) -> float:
        return self.negative_abs_sum + self.tie_abs_sum + self.positive_abs_sum

    @property
    def baseline_signed_mean(self) -> float:
        return self.negative_signed_sum + self.tie_signed_sum + self.positive_signed_sum

    @property
    def baseline_ess(self) -> float:
        denom = self.negative_weight_sq_sum + self.tie_weight_sq_sum + self.positive_weight_sq_sum
        if denom <= 0.0:
            return float("nan")
        return float(1.0 / denom)


def sign_group_statistics(
    delta: Sequence[float], weights: Sequence[float], *, tie_tol: float = TIE_TOL
) -> SignGroupStatistics:
    d = _finite_1d(delta, name="delta")
    w = normalize_weights(weights)
    if d.size != w.size:
        raise ValueError("delta and weights must have equal length")
    neg = d < -tie_tol
    pos = d > tie_tol
    tie = ~(neg | pos)
    mag = np.abs(d)

    def _sum(mask: np.ndarray, values: np.ndarray) -> float:
        return float(np.dot(w[mask], values[mask])) if np.any(mask) else 0.0

    def _mass(mask: np.ndarray) -> float:
        return float(np.sum(w[mask]))

    def _sq(mask: np.ndarray) -> float:
        return float(np.sum(w[mask] ** 2))

    out = SignGroupStatistics(
        negative_mass=_mass(neg),
        tie_mass=_mass(tie),
        positive_mass=_mass(pos),
        negative_abs_sum=_sum(neg, mag),
        tie_abs_sum=_sum(tie, mag),
        positive_abs_sum=_sum(pos, mag),
        negative_signed_sum=_sum(neg, d),
        tie_signed_sum=_sum(tie, d),
        positive_signed_sum=_sum(pos, d),
        negative_weight_sq_sum=_sq(neg),
        tie_weight_sq_sum=_sq(tie),
        positive_weight_sq_sum=_sq(pos),
    )
    if not math.isclose(out.total_mass, 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("sign-group masses do not sum to one")
    return out


def _direction(value: float, *, tie_tol: float = TIE_TOL) -> int:
    if not math.isfinite(value):
        raise ValueError("truth delta must be finite")
    if value > tie_tol:
        return 1
    if value < -tie_tol:
        return -1
    return 0


def posterior_from_sign_groups(
    stats: SignGroupStatistics,
    *,
    observed_direction: int,
    lambda_total: float,
) -> dict[str, float]:
    """Apply ``exp(lambda * s)`` where ``s`` is sign agreement in {-1,0,+1}."""
    if observed_direction not in (-1, 1):
        raise ValueError("observed_direction must be -1 or +1")
    if not math.isfinite(lambda_total) or lambda_total < 0.0:
        raise ValueError("lambda_total must be finite and non-negative")

    hi = math.exp(lambda_total)
    lo = math.exp(-lambda_total)
    if observed_direction == 1:
        f_neg, f_tie, f_pos = lo, 1.0, hi
        fav_mass0, unfav_mass0 = stats.positive_mass, stats.negative_mass
    else:
        f_neg, f_tie, f_pos = hi, 1.0, lo
        fav_mass0, unfav_mass0 = stats.negative_mass, stats.positive_mass

    z = f_neg * stats.negative_mass + f_tie * stats.tie_mass + f_pos * stats.positive_mass
    if not math.isfinite(z) or z <= 0.0:
        raise ValueError("posterior normalization is not positive finite")

    abs_mean = (
        f_neg * stats.negative_abs_sum
        + f_tie * stats.tie_abs_sum
        + f_pos * stats.positive_abs_sum
    ) / z
    signed_mean = (
        f_neg * stats.negative_signed_sum
        + f_tie * stats.tie_signed_sum
        + f_pos * stats.positive_signed_sum
    ) / z
    sq = (
        (f_neg * f_neg) * stats.negative_weight_sq_sum
        + (f_tie * f_tie) * stats.tie_weight_sq_sum
        + (f_pos * f_pos) * stats.positive_weight_sq_sum
    ) / (z * z)
    ess = float(1.0 / sq) if sq > 0.0 else float("nan")

    favorable_mass = (hi * fav_mass0) / z
    unfavorable_mass = (lo * unfav_mass0) / z
    tie_mass = stats.tie_mass / z
    return {
        "normalization": float(z),
        "magnitude_mean": float(abs_mean),
        "signed_mean": float(signed_mean),
        "joint_ess": ess,
        "favorable_mass": float(favorable_mass),
        "tie_mass": float(tie_mass),
        "unfavorable_mass": float(unfavorable_mass),
    }


def classify_gain(gain: float, *, tol: float = GAIN_TOL) -> str:
    """Classify an absolute-error gain without outcome-dependent threshold fitting."""
    if not math.isfinite(gain):
        return "UNDEFINED"
    if not math.isfinite(tol) or tol < 0.0:
        raise ValueError("tol must be finite and non-negative")
    if gain > tol:
        return "IMPROVED"
    if gain < -tol:
        return "HARMED"
    return "TIED"


def reliability_expected_gain(
    *, correct_gain: float, wrong_gain: float, reliability: float
) -> float:
    """Expected gain when the supplied binary direction is correct with probability q."""
    if not all(math.isfinite(x) for x in (correct_gain, wrong_gain, reliability)):
        raise ValueError("reliability inputs must be finite")
    if reliability < 0.0 or reliability > 1.0:
        raise ValueError("reliability must lie in [0, 1]")
    return float(reliability * correct_gain + (1.0 - reliability) * wrong_gain)


def break_even_reliability(*, correct_gain: float, wrong_gain: float, tol: float = 1e-15) -> float:
    """Return q where expected gain crosses zero, or NaN when no unique crossing exists."""
    if not all(math.isfinite(x) for x in (correct_gain, wrong_gain)):
        return float("nan")
    slope = correct_gain - wrong_gain
    if abs(slope) <= tol:
        return float("nan")
    q = -wrong_gain / slope
    if q < 0.0 or q > 1.0:
        return float("nan")
    return float(q)


def evaluate_truth_case(
    *,
    delta_b: Sequence[float],
    baseline_weights: Sequence[float],
    truth_delta_b: float,
    lambda_grid: Iterable[float] = LAMBDA_GRID,
    reliability_grid: Iterable[float] = RELIABILITY_GRID,
    tie_tol: float = TIE_TOL,
) -> dict[str, object]:
    """Evaluate correct/wrong sign-only utility for one held-out truth case.

    The primary estimate is ``E[|Delta B|]`` and the primary error is absolute
    error to ``|truth_delta_b|``.  A truth tie is retained but is not assigned a
    correct or wrong binary direction.
    """
    stats = sign_group_statistics(delta_b, baseline_weights, tie_tol=tie_tol)
    truth_direction = _direction(float(truth_delta_b), tie_tol=tie_tol)
    truth_magnitude = abs(float(truth_delta_b))
    baseline_estimate = stats.baseline_abs_mean
    baseline_error = abs(baseline_estimate - truth_magnitude)

    result: dict[str, object] = {
        "truth_direction": truth_direction,
        "truth_magnitude": truth_magnitude,
        "baseline_magnitude_estimate": baseline_estimate,
        "baseline_signed_estimate": stats.baseline_signed_mean,
        "baseline_abs_error": baseline_error,
        "baseline_joint_ess": stats.baseline_ess,
        "baseline_negative_mass": stats.negative_mass,
        "baseline_tie_mass": stats.tie_mass,
        "baseline_positive_mass": stats.positive_mass,
        "status": "TRUTH_TIE" if truth_direction == 0 else "OK",
        "lambda_results": [],
    }
    if truth_direction == 0:
        return result

    rows: list[dict[str, object]] = []
    for lam in lambda_grid:
        lam = float(lam)
        correct = posterior_from_sign_groups(
            stats, observed_direction=truth_direction, lambda_total=lam
        )
        wrong = posterior_from_sign_groups(
            stats, observed_direction=-truth_direction, lambda_total=lam
        )
        correct_error = abs(float(correct["magnitude_mean"]) - truth_magnitude)
        wrong_error = abs(float(wrong["magnitude_mean"]) - truth_magnitude)
        correct_gain = baseline_error - correct_error
        wrong_gain = baseline_error - wrong_error
        reliability_rows = [
            {
                "reliability": float(q),
                "expected_abs_error_gain": reliability_expected_gain(
                    correct_gain=correct_gain, wrong_gain=wrong_gain, reliability=float(q)
                ),
            }
            for q in reliability_grid
        ]
        rows.append(
            {
                "lambda_total": lam,
                "correct": correct,
                "wrong": wrong,
                "correct_abs_error": correct_error,
                "wrong_abs_error": wrong_error,
                "correct_abs_error_gain": correct_gain,
                "wrong_abs_error_gain": wrong_gain,
                "random_sign_expected_abs_error_gain": reliability_expected_gain(
                    correct_gain=correct_gain, wrong_gain=wrong_gain, reliability=0.5
                ),
                "break_even_reliability": break_even_reliability(
                    correct_gain=correct_gain, wrong_gain=wrong_gain
                ),
                "reliability_curve": reliability_rows,
            }
        )
    result["lambda_results"] = rows
    return result


def analysis_contract() -> dict[str, object]:
    return {
        "schema": "bridge.pl2.sign_only_utility.v1",
        "stage_role": "TRUTH_REFERENCED_NETWORK_WIDE_SIGN_ONLY_UTILITY",
        "strong_anchor_reaction": "HEX1",
        "primary_target": "absolute_between_condition_difference_magnitude_E_abs_delta_B",
        "weak_information": "binary_direction_only",
        "primary_lambda": PRIMARY_LAMBDA,
        "lambda_grid": list(LAMBDA_GRID),
        "reliability_grid": list(RELIABILITY_GRID),
        "random_sign_null": "q_equals_0.5_exact_expectation_of_correct_and_wrong_direction_errors",
        "gain_classification": {
            "metric": "baseline_abs_error_minus_posterior_abs_error",
            "tolerance": GAIN_TOL,
            "positive": "IMPROVED",
            "within_tolerance": "TIED",
            "negative": "HARMED",
        },
        "primary_benchmark_unit": "distinct_held_out_truth_pair_once",
        "quantile_alias_policy": "preserve_all_aliases_for_audit_but_do_not_duplicate_weight_primary_summaries",
        "reliability_storage": "derive_from_correct_and_wrong_endpoints_without_materializing_q_rows",
        "pl1_role": "frozen_outcome_blind_predictors_only",
        "truth_selection": {
            "basis": "strong_anchor_HEX1_coordinate_and_strong_anchor_weights_only",
            "quantiles": [0.10, 0.50, 0.90],
            "collapse_duplicate_candidate_selections": True,
            "pair_only_shared_quantile_aliases": True,
            "remove_truth_candidate_from_each_condition_before_inference": True,
            "reaction_B_values_must_not_select_truth_candidates": True,
        },
        "forbidden": [
            "PL2_outcome_driven_reaction_selection",
            "PL2_outcome_driven_lambda_selection",
            "PL2_outcome_driven_predictor_threshold_fitting",
            "gene_score_use_in_primary_PL2",
            "historical_tanh_operator_relabeling",
            "solver_or_optimization",
            "FVA",
            "new_sampling",
            "model_reconstruction",
            "new_strong_anchor_weight_fitting",
        ],
    }
