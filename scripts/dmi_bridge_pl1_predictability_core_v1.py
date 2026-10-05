#!/usr/bin/env python3
"""DMI-BRIDGE-PL1: outcome-blind network-wide predictability landscape primitives.

PL1 characterizes the finite candidate geometry available after the qualified strong
HEX1 anchor. It does not inject a weak prior, use hidden/synthetic truth, use gene
scores, choose a reaction B, or run a solver/FVA/sampler/reconstruction.

Repository-bound code must scan all admissible finite reactions except HEX1 and keep
same-pathway/proximal reactions in the primary analysis. Pathway relation is annotation
only; pathway exclusion is sensitivity-only.
"""
from __future__ import annotations

import hashlib
import json
import math
from typing import Any, Sequence

import numpy as np

SCHEMA = "bridge.pl1.predictability_landscape.v1"
COMPLETE_STATUS = "PL1_PREDICTABILITY_LANDSCAPE_COMPLETE"
STRONG_ANCHOR_REACTION = "HEX1"
TIE_TOL = 1e-12
QUANTILES = (0.10, 0.50, 0.90)


def canonical_json_bytes(obj: Any) -> bytes:
    return (json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=True) + "\n").encode("ascii")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _finite_1d(values: Sequence[float], *, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError(f"{name} must be a non-empty 1-D vector")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def normalize_weights(weights: Sequence[float]) -> np.ndarray:
    w = _finite_1d(weights, name="weights")
    if np.any(w < 0):
        raise ValueError("weights must be non-negative")
    total = float(np.sum(w))
    if not math.isfinite(total) or total <= 0:
        raise ValueError("weights must have positive finite mass")
    return w / total


def effective_sample_size(weights: Sequence[float]) -> float:
    w = normalize_weights(weights)
    return float(1.0 / np.dot(w, w))


def weighted_mean_variance(values: Sequence[float], weights: Sequence[float]) -> tuple[float, float]:
    x = _finite_1d(values, name="values")
    w = normalize_weights(weights)
    if x.size != w.size:
        raise ValueError("values and weights must have equal length")
    mean = float(np.dot(w, x))
    var = float(np.dot(w, (x - mean) ** 2))
    return mean, max(var, 0.0)


def weighted_quantile(values: Sequence[float], weights: Sequence[float], q: float) -> float:
    """Stable weighted inverse-CDF quantile: first value whose cumulative mass reaches q."""
    if not (0.0 <= q <= 1.0):
        raise ValueError("q must lie in [0,1]")
    x = _finite_1d(values, name="values")
    w = normalize_weights(weights)
    if x.size != w.size:
        raise ValueError("values and weights must have equal length")
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
    cs = np.cumsum(w[order])
    if q <= 0.0:
        return float(xs[0])
    idx = int(np.searchsorted(cs, q, side="left"))
    return float(xs[min(idx, xs.size - 1)])


def weighted_correlation(
    x: Sequence[float],
    y: Sequence[float],
    weights: Sequence[float],
    *,
    tol: float = TIE_TOL,
) -> float:
    a = _finite_1d(x, name="x")
    b = _finite_1d(y, name="y")
    w = normalize_weights(weights)
    if a.size != b.size or a.size != w.size:
        raise ValueError("x, y, and weights must have equal length")
    ma = float(np.dot(w, a))
    mb = float(np.dot(w, b))
    ac = a - ma
    bc = b - mb
    va = float(np.dot(w, ac * ac))
    vb = float(np.dot(w, bc * bc))
    if va <= tol or vb <= tol:
        return float("nan")
    cov = float(np.dot(w, ac * bc))
    return float(cov / math.sqrt(va * vb))


def joint_delta_distribution(
    ct2a_values: Sequence[float],
    gl261_values: Sequence[float],
    ct2a_weights: Sequence[float],
    gl261_weights: Sequence[float],
) -> tuple[np.ndarray, np.ndarray]:
    c = _finite_1d(ct2a_values, name="ct2a_values")
    g = _finite_1d(gl261_values, name="gl261_values")
    wc = normalize_weights(ct2a_weights)
    wg = normalize_weights(gl261_weights)
    if c.size != wc.size or g.size != wg.size:
        raise ValueError("value/weight dimensions do not match")
    delta = (c[:, None] - g[None, :]).reshape(-1)
    joint = (wc[:, None] * wg[None, :]).reshape(-1)
    return delta, normalize_weights(joint)


def _conditional_summary(values: np.ndarray, weights: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    mass = float(np.sum(weights[mask]))
    if mass <= 0:
        return {
            "mass": 0.0,
            "ess": float("nan"),
            "mean": float("nan"),
            "q10": float("nan"),
            "q50": float("nan"),
            "q90": float("nan"),
            "width80": float("nan"),
        }
    x = values[mask]
    w = normalize_weights(weights[mask])
    mean, _ = weighted_mean_variance(x, w)
    q10 = weighted_quantile(x, w, 0.10)
    q50 = weighted_quantile(x, w, 0.50)
    q90 = weighted_quantile(x, w, 0.90)
    return {
        "mass": mass,
        "ess": effective_sample_size(w),
        "mean": mean,
        "q10": q10,
        "q50": q50,
        "q90": q90,
        "width80": q90 - q10,
    }


def sign_magnitude_eta2(delta: Sequence[float], weights: Sequence[float], *, tol: float = TIE_TOL) -> float:
    """Weighted eta^2 for sign category {-1,0,+1} explaining |delta| variation."""
    d = _finite_1d(delta, name="delta")
    w = normalize_weights(weights)
    if d.size != w.size:
        raise ValueError("delta and weights must have equal length")
    sign = np.where(d > tol, 1, np.where(d < -tol, -1, 0))
    mag = np.abs(d)
    overall = float(np.dot(w, mag))
    total = float(np.dot(w, (mag - overall) ** 2))
    if total <= tol:
        return float("nan")
    between = 0.0
    for group in (-1, 0, 1):
        mask = sign == group
        mass = float(np.sum(w[mask]))
        if mass <= 0:
            continue
        mean_g = float(np.dot(w[mask], mag[mask]) / mass)
        between += mass * (mean_g - overall) ** 2
    eta2 = between / total
    return float(min(1.0, max(0.0, eta2)))


def _safe_contraction(strong_value: float, reference_value: float, *, tol: float = TIE_TOL) -> float:
    if not math.isfinite(reference_value) or reference_value <= tol:
        return float("nan")
    return float(1.0 - strong_value / reference_value)


def target_landscape_metrics(
    *,
    anchor_ct2a: Sequence[float],
    anchor_gl261: Sequence[float],
    target_ct2a: Sequence[float],
    target_gl261: Sequence[float],
    strong_weights_ct2a: Sequence[float],
    strong_weights_gl261: Sequence[float],
    tie_tol: float = TIE_TOL,
) -> dict[str, float | str | bool | int]:
    """Compute PL1 geometry for one reaction B in one exact context/mouse contrast.

    The strong weights are normalized within the exact reconstruction/RNA stratum.
    The reference distribution is uniform on the exact same finite candidate support;
    therefore contraction fields describe weight-induced concentration, not FVA or
    feasible-set contraction.
    """
    ac = _finite_1d(anchor_ct2a, name="anchor_ct2a")
    ag = _finite_1d(anchor_gl261, name="anchor_gl261")
    bc = _finite_1d(target_ct2a, name="target_ct2a")
    bg = _finite_1d(target_gl261, name="target_gl261")
    if ac.size != bc.size or ag.size != bg.size:
        raise ValueError("anchor/target candidate counts must match within each condition")

    wc = normalize_weights(strong_weights_ct2a)
    wg = normalize_weights(strong_weights_gl261)
    if wc.size != ac.size or wg.size != ag.size:
        raise ValueError("strong weights do not match candidate counts")

    delta_a, joint = joint_delta_distribution(ac, ag, wc, wg)
    delta_b, joint_b = joint_delta_distribution(bc, bg, wc, wg)
    if not np.allclose(joint, joint_b, rtol=0.0, atol=1e-15):
        raise ValueError("anchor/target joint weights differ")

    uc = np.full(ac.size, 1.0 / ac.size)
    ug = np.full(ag.size, 1.0 / ag.size)
    delta_b_uniform, uniform_joint = joint_delta_distribution(bc, bg, uc, ug)

    mean_b, var_b = weighted_mean_variance(delta_b, joint)
    mean_u, var_u = weighted_mean_variance(delta_b_uniform, uniform_joint)
    sd_b = math.sqrt(var_b)
    sd_u = math.sqrt(var_u)

    q10 = weighted_quantile(delta_b, joint, 0.10)
    q50 = weighted_quantile(delta_b, joint, 0.50)
    q90 = weighted_quantile(delta_b, joint, 0.90)
    uq10 = weighted_quantile(delta_b_uniform, uniform_joint, 0.10)
    uq50 = weighted_quantile(delta_b_uniform, uniform_joint, 0.50)
    uq90 = weighted_quantile(delta_b_uniform, uniform_joint, 0.90)
    width80 = q90 - q10
    uniform_width80 = uq90 - uq10

    pos = delta_b > tie_tol
    neg = delta_b < -tie_tol
    tie = ~(pos | neg)
    pos_s = _conditional_summary(np.abs(delta_b), joint, pos)
    neg_s = _conditional_summary(np.abs(delta_b), joint, neg)
    tie_s = _conditional_summary(np.abs(delta_b), joint, tie)

    masses = np.asarray([neg_s["mass"], tie_s["mass"], pos_s["mass"]], dtype=float)
    nz = masses[masses > 0]
    entropy = 0.0 if nz.size <= 1 else float(-np.sum(nz * np.log(nz)) / math.log(3.0))

    sign_numeric = np.where(pos, 1.0, np.where(neg, -1.0, 0.0))
    eta2 = sign_magnitude_eta2(delta_b, joint, tol=tie_tol)
    sign_mag_corr = weighted_correlation(sign_numeric, np.abs(delta_b), joint, tol=tie_tol)
    corr_ab = weighted_correlation(delta_a, delta_b, joint, tol=tie_tol)
    corr_abs_ab = weighted_correlation(np.abs(delta_a), np.abs(delta_b), joint, tol=tie_tol)

    mean_a, var_a = weighted_mean_variance(delta_a, joint)
    _, var_abs_b = weighted_mean_variance(np.abs(delta_b), joint)

    return {
        "cartesian_support": int(delta_b.size),
        "joint_ess": effective_sample_size(joint),
        "max_joint_mass": float(np.max(joint)),
        "delta_a_mean": mean_a,
        "delta_a_sd": math.sqrt(var_a),
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
        "anchor_sd_contraction": _safe_contraction(sd_b, sd_u, tol=tie_tol),
        "anchor_width80_contraction": _safe_contraction(width80, uniform_width80, tol=tie_tol),
        "p_negative": neg_s["mass"],
        "p_tie": tie_s["mass"],
        "p_positive": pos_s["mass"],
        "dominant_sign_mass": float(np.max(masses)),
        "directional_entropy3": entropy,
        "n_supported_sign_states": int(np.count_nonzero(masses > 0.0)),
        "negative_ess": neg_s["ess"],
        "tie_ess": tie_s["ess"],
        "positive_ess": pos_s["ess"],
        "negative_abs_mean": neg_s["mean"],
        "negative_abs_q10": neg_s["q10"],
        "negative_abs_q50": neg_s["q50"],
        "negative_abs_q90": neg_s["q90"],
        "negative_abs_width80": neg_s["width80"],
        "positive_abs_mean": pos_s["mean"],
        "positive_abs_q10": pos_s["q10"],
        "positive_abs_q50": pos_s["q50"],
        "positive_abs_q90": pos_s["q90"],
        "positive_abs_width80": pos_s["width80"],
        "negative_width80_contraction_vs_all": _safe_contraction(neg_s["width80"], width80, tol=tie_tol),
        "positive_width80_contraction_vs_all": _safe_contraction(pos_s["width80"], width80, tol=tie_tol),
        "sign_magnitude_eta2": eta2,
        "sign_magnitude_correlation": sign_mag_corr,
        "delta_anchor_target_correlation": corr_ab,
        "abs_delta_anchor_target_correlation": corr_abs_ab,
        "target_delta_degenerate": bool(var_b <= tie_tol and width80 <= tie_tol),
        "target_abs_magnitude_degenerate": bool(var_abs_b <= tie_tol),
        "anchor_delta_degenerate": bool(var_a <= tie_tol),
        "sign_magnitude_status": "DEGENERATE" if not math.isfinite(eta2) else "OK",
        "anchor_target_coupling_status": "DEGENERATE" if not math.isfinite(corr_ab) else "OK",
    }


def analysis_contract() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "status_on_success": COMPLETE_STATUS,
        "strong_anchor_reaction": STRONG_ANCHOR_REACTION,
        "scientific_question": "under_what_preexisting_strong_anchor_geometry_can_sign_information_about_B_inform_B_magnitude",
        "stage_role": "OUTCOME_BLIND_NETWORK_WIDE_PREDICTABILITY_LANDSCAPE",
        "primary_scan": {
            "scope": "all_admissible_finite_flux_cache_reactions_except_HEX1",
            "preselect_LDH_or_TCA": False,
            "exclude_same_pathway_or_proximal_reactions": False,
            "pathway_relation_is_annotation": True,
            "pathway_exclusion_is_sensitivity_only": True,
        },
        "evaluation_unit": "exact_reconstruction_method_x_RNA_context_x_CT2A_mouse_x_GL261_mouse",
        "candidate_geometry": {
            "use_existing_strong_anchor_ensembles": True,
            "condition_on_exact_stratum": True,
            "reference": "uniform_weights_on_the_exact_same_finite_candidate_support",
            "finite_support_not_FVA": True,
            "anchor_contraction_means_weight_induced_concentration_not_feasible_set_contraction": True,
        },
        "metrics": [
            "anchor_weight_induced_sd_and_q10_q90_width_contraction",
            "positive_tie_negative_mass_and_entropy",
            "joint_and_sign_specific_ESS",
            "within_sign_absolute_delta_concentration",
            "sign_magnitude_eta2_and_correlation",
            "delta_HEX1_delta_B_coupling",
            "pathway_and_proximity_annotation",
        ],
        "forbidden_in_PL1": [
            "hidden_or_synthetic_truth",
            "gene_score_concordance",
            "weak_prior_injection",
            "correct_wrong_or_random_direction_performance",
            "reaction_B_selection_or_ranking",
            "bridgeability_threshold_fitting",
            "solver_or_optimization",
            "FVA",
            "flux_sampling",
            "model_reconstruction",
            "new_weight_fitting",
        ],
        "next_stage": "PL2_truth_referenced_sign_only_utility_after_PL1_is_frozen",
    }
