#!/usr/bin/env python3
"""Vectorized accelerator for DMI-BRIDGE-PL2B sign-only utility.

The scientific authority remains ``evaluate_truth_case`` in
``dmi_bridge_pl2_sign_only_core_v1``.  This module applies the same sufficient-
statistics equations to a block of reaction columns.
"""
from __future__ import annotations

from collections.abc import Sequence
import math

import numpy as np

try:
    from scripts import dmi_bridge_pl2_sign_only_core_v1 as core
except ModuleNotFoundError:  # direct script execution
    import dmi_bridge_pl2_sign_only_core_v1 as core


def _matrix(values, *, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 2 or min(arr.shape) <= 0:
        raise ValueError(f"{name} must be a non-empty candidate-by-reaction matrix")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr


def _vector(values, *, name: str, size: int | None = None) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError(f"{name} must be a non-empty vector")
    if size is not None and arr.size != size:
        raise ValueError(f"{name} has the wrong length")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr


def _posterior(
    *,
    masses: np.ndarray,
    abs_sums: np.ndarray,
    signed_sums: np.ndarray,
    weight_sq_sums: np.ndarray,
    observed_direction: int,
    lambda_total: float,
) -> dict[str, np.ndarray]:
    if observed_direction not in (-1, 1):
        raise ValueError("observed_direction must be -1 or +1")
    if not math.isfinite(lambda_total) or lambda_total < 0.0:
        raise ValueError("lambda_total must be finite and non-negative")
    hi = math.exp(lambda_total)
    lo = math.exp(-lambda_total)
    factors = np.asarray((lo, 1.0, hi) if observed_direction == 1 else (hi, 1.0, lo))
    z = np.sum(factors[:, None] * masses, axis=0)
    if np.any(~np.isfinite(z)) or np.any(z <= 0.0):
        raise ValueError("posterior normalization is not positive finite")
    magnitude = np.sum(factors[:, None] * abs_sums, axis=0) / z
    signed = np.sum(factors[:, None] * signed_sums, axis=0) / z
    denom = np.sum((factors[:, None] ** 2) * weight_sq_sums, axis=0) / (z * z)
    ess = np.full(z.shape, np.nan)
    ok = denom > 0.0
    ess[ok] = 1.0 / denom[ok]
    return {
        "normalization": z,
        "magnitude_mean": magnitude,
        "signed_mean": signed,
        "joint_ess": ess,
    }


def evaluate_truth_block(
    *,
    ct2a_values,
    gl261_values,
    ct2a_weights: Sequence[float],
    gl261_weights: Sequence[float],
    truth_ct2a: Sequence[float],
    truth_gl261: Sequence[float],
    lambda_grid: Sequence[float] = core.LAMBDA_GRID,
    tie_tol: float = core.TIE_TOL,
) -> dict[str, object]:
    """Evaluate a block of reaction-B columns using Cartesian held-out support."""
    left = _matrix(ct2a_values, name="ct2a_values")
    right = _matrix(gl261_values, name="gl261_values")
    if left.shape[1] != right.shape[1]:
        raise ValueError("condition matrices must have the same reaction count")
    wl = core.normalize_weights(ct2a_weights)
    wr = core.normalize_weights(gl261_weights)
    if left.shape[0] != wl.size or right.shape[0] != wr.size:
        raise ValueError("candidate and weight dimensions do not match")
    truth_left = _vector(truth_ct2a, name="truth_ct2a", size=left.shape[1])
    truth_right = _vector(truth_gl261, name="truth_gl261", size=left.shape[1])

    joint = (wl[:, None] * wr[None, :]).reshape(-1)
    delta = (left[:, None, :] - right[None, :, :]).reshape(-1, left.shape[1])
    negative = delta < -tie_tol
    positive = delta > tie_tol
    tie = ~(negative | positive)
    masks = (negative, tie, positive)
    magnitude = np.abs(delta)
    masses = np.stack([np.sum(joint[:, None] * mask, axis=0) for mask in masks])
    if not np.allclose(np.sum(masses, axis=0), 1.0, rtol=0.0, atol=1e-12):
        raise ValueError("sign-group masses do not sum to one")
    abs_sums = np.stack([np.sum(joint[:, None] * magnitude * mask, axis=0) for mask in masks])
    signed_sums = np.stack([np.sum(joint[:, None] * delta * mask, axis=0) for mask in masks])
    weight_sq_sums = np.stack(
        [np.sum((joint[:, None] ** 2) * mask, axis=0) for mask in masks]
    )

    truth_delta = truth_left - truth_right
    truth_direction = np.where(
        truth_delta > tie_tol, 1, np.where(truth_delta < -tie_tol, -1, 0)
    ).astype(np.int8)
    truth_magnitude = np.abs(truth_delta)
    baseline_magnitude = np.sum(abs_sums, axis=0)
    baseline_signed = np.sum(signed_sums, axis=0)
    baseline_ess = np.full(left.shape[1], 1.0 / np.sum(joint * joint))
    baseline_error = np.abs(baseline_magnitude - truth_magnitude)

    result: dict[str, object] = {
        "truth_delta_b": truth_delta,
        "truth_direction": truth_direction,
        "truth_magnitude": truth_magnitude,
        "baseline_magnitude_estimate": baseline_magnitude,
        "baseline_signed_estimate": baseline_signed,
        "baseline_abs_error": baseline_error,
        "baseline_joint_ess": baseline_ess,
        "baseline_negative_mass": masses[0],
        "baseline_tie_mass": masses[1],
        "baseline_positive_mass": masses[2],
        "lambda_results": {},
        "cartesian_support": int(joint.size),
    }
    lambda_results: dict[float, dict[str, np.ndarray]] = {}
    non_tie = truth_direction != 0
    for value in lambda_grid:
        lam = float(value)
        positive_post = _posterior(
            masses=masses,
            abs_sums=abs_sums,
            signed_sums=signed_sums,
            weight_sq_sums=weight_sq_sums,
            observed_direction=1,
            lambda_total=lam,
        )
        negative_post = _posterior(
            masses=masses,
            abs_sums=abs_sums,
            signed_sums=signed_sums,
            weight_sq_sums=weight_sq_sums,
            observed_direction=-1,
            lambda_total=lam,
        )
        correct_is_positive = truth_direction == 1
        correct_mag = np.where(
            correct_is_positive, positive_post["magnitude_mean"], negative_post["magnitude_mean"]
        )
        wrong_mag = np.where(
            correct_is_positive, negative_post["magnitude_mean"], positive_post["magnitude_mean"]
        )
        correct_ess = np.where(
            correct_is_positive, positive_post["joint_ess"], negative_post["joint_ess"]
        )
        wrong_ess = np.where(
            correct_is_positive, negative_post["joint_ess"], positive_post["joint_ess"]
        )
        correct_signed = np.where(
            correct_is_positive, positive_post["signed_mean"], negative_post["signed_mean"]
        )
        wrong_signed = np.where(
            correct_is_positive, negative_post["signed_mean"], positive_post["signed_mean"]
        )
        correct_error = np.abs(correct_mag - truth_magnitude)
        wrong_error = np.abs(wrong_mag - truth_magnitude)
        correct_gain = baseline_error - correct_error
        wrong_gain = baseline_error - wrong_error
        for array in (
            correct_mag, wrong_mag, correct_ess, wrong_ess, correct_signed, wrong_signed,
            correct_error, wrong_error, correct_gain, wrong_gain,
        ):
            array[~non_tie] = np.nan
        slope = correct_gain - wrong_gain
        break_even = np.full(left.shape[1], np.nan)
        crossing = non_tie & np.isfinite(slope) & (np.abs(slope) > 1e-15)
        candidates = np.full(left.shape[1], np.nan)
        candidates[crossing] = -wrong_gain[crossing] / slope[crossing]
        inside = crossing & (candidates >= 0.0) & (candidates <= 1.0)
        break_even[inside] = candidates[inside]
        lambda_results[lam] = {
            "correct_magnitude_estimate": correct_mag,
            "wrong_magnitude_estimate": wrong_mag,
            "correct_signed_estimate": correct_signed,
            "wrong_signed_estimate": wrong_signed,
            "correct_absolute_error": correct_error,
            "wrong_absolute_error": wrong_error,
            "correct_absolute_error_gain": correct_gain,
            "wrong_absolute_error_gain": wrong_gain,
            "correct_posterior_ess": correct_ess,
            "wrong_posterior_ess": wrong_ess,
            "random_sign_expected_gain": 0.5 * (correct_gain + wrong_gain),
            "break_even_reliability": break_even,
        }
    result["lambda_results"] = lambda_results
    return result
