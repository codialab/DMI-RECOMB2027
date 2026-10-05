#!/usr/bin/env python3
"""Vectorized execution accelerator for DMI-BRIDGE-PL1.

The scalar scientific authority remains target_landscape_metrics() in
``dmi_bridge_pl1_predictability_core_v1``.  This module evaluates the same finite
20x20 Cartesian geometry for a block of reaction-B columns.  Repository code must
validate this accelerator against the scalar authority on deterministic sentinel
reactions before publication.
"""
from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from scripts.dmi_bridge_pl1_predictability_core_v1 import TIE_TOL, normalize_weights


def _as_matrix(values: np.ndarray | Sequence[Sequence[float]], *, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 2 or arr.shape[0] == 0 or arr.shape[1] == 0:
        raise ValueError(f"{name} must be a non-empty 2-D candidate-by-reaction matrix")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def _weighted_mean_variance_matrix(x: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = np.sum(w[:, None] * x, axis=0)
    var = np.sum(w[:, None] * (x - mean[None, :]) ** 2, axis=0)
    return mean, np.maximum(var, 0.0)


def _weighted_quantiles_matrix(x: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    order = np.argsort(x, axis=0, kind="stable")
    xs = np.take_along_axis(x, order, axis=0)
    ww = np.take_along_axis(np.broadcast_to(w[:, None], x.shape), order, axis=0)
    cs = np.cumsum(ww, axis=0)
    cols = np.arange(x.shape[1])
    out = []
    for q in (0.10, 0.50, 0.90):
        idx = np.argmax(cs >= q, axis=0)
        out.append(xs[idx, cols])
    return out[0], out[1], out[2]


def _conditional_summary_matrix(
    abs_delta: np.ndarray,
    joint: np.ndarray,
    mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    k = abs_delta.shape[1]
    weighted = joint[:, None] * mask
    mass = np.sum(weighted, axis=0)
    ok = mass > 0.0

    denom = np.sum((joint[:, None] ** 2) * mask, axis=0)
    ess = np.full(k, np.nan)
    ess[ok] = mass[ok] ** 2 / denom[ok]

    mean = np.full(k, np.nan)
    numer = np.sum(weighted * abs_delta, axis=0)
    mean[ok] = numer[ok] / mass[ok]

    order = np.argsort(abs_delta, axis=0, kind="stable")
    xs = np.take_along_axis(abs_delta, order, axis=0)
    ww = np.take_along_axis(weighted, order, axis=0)
    cs = np.cumsum(ww, axis=0)
    cols = np.arange(k)
    quantiles = []
    for q in (0.10, 0.50, 0.90):
        idx = np.argmax(cs >= q * mass[None, :], axis=0)
        values = xs[idx, cols].copy()
        values[~ok] = np.nan
        quantiles.append(values)
    q10, q50, q90 = quantiles
    return mass, ess, mean, q10, q50, q90, q90 - q10


def _safe_contraction_array(strong: np.ndarray, reference: np.ndarray, *, tol: float) -> np.ndarray:
    out = np.full(strong.shape, np.nan)
    ok = np.isfinite(reference) & (reference > tol)
    out[ok] = 1.0 - strong[ok] / reference[ok]
    return out


def target_landscape_metrics_block(
    *,
    anchor_ct2a: Sequence[float],
    anchor_gl261: Sequence[float],
    target_ct2a: np.ndarray | Sequence[Sequence[float]],
    target_gl261: np.ndarray | Sequence[Sequence[float]],
    strong_weights_ct2a: Sequence[float],
    strong_weights_gl261: Sequence[float],
    tie_tol: float = TIE_TOL,
) -> dict[str, np.ndarray]:
    """Evaluate scalar PL1 definitions for a reaction block.

    ``target_ct2a`` and ``target_gl261`` have shape candidates x reactions.  The
    returned dictionary contains one vector per scalar authority field.  Floating
    reductions are vectorized, so last-bit roundoff may differ from separate BLAS
    scalar calls; callers must enforce deterministic scalar-parity tolerances.
    """
    ac = np.asarray(anchor_ct2a, dtype=float)
    ag = np.asarray(anchor_gl261, dtype=float)
    if ac.ndim != 1 or ag.ndim != 1 or ac.size == 0 or ag.size == 0:
        raise ValueError("anchor vectors must be non-empty 1-D vectors")
    if not np.all(np.isfinite(ac)) or not np.all(np.isfinite(ag)):
        raise ValueError("anchor vectors contain non-finite values")
    bc = _as_matrix(target_ct2a, name="target_ct2a")
    bg = _as_matrix(target_gl261, name="target_gl261")
    if bc.shape[1] != bg.shape[1] or bc.shape[0] != ac.size or bg.shape[0] != ag.size:
        raise ValueError("anchor/target candidate or reaction dimensions do not match")

    wc = normalize_weights(strong_weights_ct2a)
    wg = normalize_weights(strong_weights_gl261)
    if wc.size != ac.size or wg.size != ag.size:
        raise ValueError("strong weights do not match candidate counts")

    joint = (wc[:, None] * wg[None, :]).reshape(-1)
    joint = joint / np.sum(joint)
    uniform = np.full(joint.size, 1.0 / joint.size)
    delta_a = (ac[:, None] - ag[None, :]).reshape(-1)
    delta_b = (bc[:, None, :] - bg[None, :, :]).reshape(-1, bc.shape[1])
    abs_b = np.abs(delta_b)
    k = delta_b.shape[1]

    mean_b, var_b = _weighted_mean_variance_matrix(delta_b, joint)
    mean_u, var_u = _weighted_mean_variance_matrix(delta_b, uniform)
    sd_b = np.sqrt(var_b)
    sd_u = np.sqrt(var_u)
    q10, q50, q90 = _weighted_quantiles_matrix(delta_b, joint)
    uq10, uq50, uq90 = _weighted_quantiles_matrix(delta_b, uniform)
    width80 = q90 - q10
    uniform_width80 = uq90 - uq10

    pos = delta_b > tie_tol
    neg = delta_b < -tie_tol
    tie = ~(pos | neg)
    neg_s = _conditional_summary_matrix(abs_b, joint, neg)
    tie_s = _conditional_summary_matrix(abs_b, joint, tie)
    pos_s = _conditional_summary_matrix(abs_b, joint, pos)

    masses = np.stack((neg_s[0], tie_s[0], pos_s[0]), axis=0)
    supported = masses > 0.0
    n_supported = np.count_nonzero(supported, axis=0)
    log_mass = np.zeros_like(masses)
    log_mass[supported] = np.log(masses[supported])
    raw_entropy = -np.sum(np.where(supported, masses * log_mass, 0.0), axis=0) / math.log(3.0)
    entropy = np.zeros(k)
    entropy[n_supported > 1] = raw_entropy[n_supported > 1]

    abs_mean, abs_var = _weighted_mean_variance_matrix(abs_b, joint)
    sign_numeric = np.where(pos, 1.0, np.where(neg, -1.0, 0.0))
    sign_mean, sign_var = _weighted_mean_variance_matrix(sign_numeric, joint)
    sign_abs_cov = np.sum(
        joint[:, None] * (sign_numeric - sign_mean[None, :]) * (abs_b - abs_mean[None, :]), axis=0
    )
    sign_mag_corr = np.full(k, np.nan)
    corr_ok = (sign_var > tie_tol) & (abs_var > tie_tol)
    sign_mag_corr[corr_ok] = sign_abs_cov[corr_ok] / np.sqrt(sign_var[corr_ok] * abs_var[corr_ok])

    between = np.zeros(k)
    for summary in (neg_s, tie_s, pos_s):
        mass, group_mean = summary[0], summary[2]
        ok = mass > 0.0
        between[ok] += mass[ok] * (group_mean[ok] - abs_mean[ok]) ** 2
    eta2 = np.full(k, np.nan)
    eta_ok = abs_var > tie_tol
    eta2[eta_ok] = np.clip(between[eta_ok] / abs_var[eta_ok], 0.0, 1.0)

    mean_a = float(np.sum(joint * delta_a))
    var_a = float(np.sum(joint * (delta_a - mean_a) ** 2))
    var_a = max(var_a, 0.0)
    cov_ab = np.sum(joint[:, None] * (delta_a - mean_a)[:, None] * (delta_b - mean_b[None, :]), axis=0)
    corr_ab = np.full(k, np.nan)
    ab_ok = (var_a > tie_tol) & (var_b > tie_tol)
    corr_ab[ab_ok] = cov_ab[ab_ok] / np.sqrt(var_a * var_b[ab_ok])

    abs_a = np.abs(delta_a)
    abs_a_mean = float(np.sum(joint * abs_a))
    abs_a_var = float(np.sum(joint * (abs_a - abs_a_mean) ** 2))
    abs_a_var = max(abs_a_var, 0.0)
    cov_abs = np.sum(joint[:, None] * (abs_a - abs_a_mean)[:, None] * (abs_b - abs_mean[None, :]), axis=0)
    corr_abs = np.full(k, np.nan)
    abs_ok = (abs_a_var > tie_tol) & (abs_var > tie_tol)
    corr_abs[abs_ok] = cov_abs[abs_ok] / np.sqrt(abs_a_var * abs_var[abs_ok])

    return {
        "cartesian_support": np.full(k, delta_b.shape[0], dtype=float),
        "joint_ess": np.full(k, 1.0 / np.dot(joint, joint), dtype=float),
        "max_joint_mass": np.full(k, float(np.max(joint)), dtype=float),
        "delta_a_mean": np.full(k, mean_a, dtype=float),
        "delta_a_sd": np.full(k, math.sqrt(var_a), dtype=float),
        "delta_b_mean": mean_b,
        "delta_b_sd": sd_b,
        "delta_b_q10": q10,
        "delta_b_q50": q50,
        "delta_b_q90": q90,
        "delta_b_width80": width80,
        "uniform_delta_b_mean": mean_u,
        "uniform_delta_b_sd": sd_u,
        "uniform_delta_b_q10": uq10,
        "uniform_delta_b_q50": uq50,
        "uniform_delta_b_q90": uq90,
        "uniform_delta_b_width80": uniform_width80,
        "anchor_sd_contraction": _safe_contraction_array(sd_b, sd_u, tol=tie_tol),
        "anchor_width80_contraction": _safe_contraction_array(width80, uniform_width80, tol=tie_tol),
        "p_negative": neg_s[0],
        "p_tie": tie_s[0],
        "p_positive": pos_s[0],
        "dominant_sign_mass": np.max(masses, axis=0),
        "directional_entropy3": entropy,
        "n_supported_sign_states": n_supported.astype(float),
        "negative_ess": neg_s[1],
        "tie_ess": tie_s[1],
        "positive_ess": pos_s[1],
        "negative_abs_mean": neg_s[2],
        "negative_abs_q10": neg_s[3],
        "negative_abs_q50": neg_s[4],
        "negative_abs_q90": neg_s[5],
        "negative_abs_width80": neg_s[6],
        "positive_abs_mean": pos_s[2],
        "positive_abs_q10": pos_s[3],
        "positive_abs_q50": pos_s[4],
        "positive_abs_q90": pos_s[5],
        "positive_abs_width80": pos_s[6],
        "negative_width80_contraction_vs_all": _safe_contraction_array(neg_s[6], width80, tol=tie_tol),
        "positive_width80_contraction_vs_all": _safe_contraction_array(pos_s[6], width80, tol=tie_tol),
        "sign_magnitude_eta2": eta2,
        "sign_magnitude_correlation": sign_mag_corr,
        "delta_anchor_target_correlation": corr_ab,
        "abs_delta_anchor_target_correlation": corr_abs,
        "target_delta_degenerate": (var_b <= tie_tol) & (width80 <= tie_tol),
        "target_abs_magnitude_degenerate": abs_var <= tie_tol,
        "anchor_delta_degenerate": np.full(k, var_a <= tie_tol, dtype=bool),
        "sign_magnitude_status": np.where(np.isfinite(eta2), "OK", "DEGENERATE"),
        "anchor_target_coupling_status": np.where(np.isfinite(corr_ab), "OK", "DEGENERATE"),
    }
