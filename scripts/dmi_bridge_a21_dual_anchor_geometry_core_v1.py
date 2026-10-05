#!/usr/bin/env python3
"""Pure numerical helpers for DMI-BRIDGE-A2.1.

A2.1 reuses the frozen PL1 geometry definitions under the qualified A2-L
(Vmax + Vlac) weights. HEX1 remains the common anchor coordinate used for
A1-versus-A2 comparable anchor-target coupling. The lactate production
coordinate is added only as a separately named supportive coupling feature.
"""
from __future__ import annotations

import math
from typing import Any, Sequence

import numpy as np

from scripts.dmi_bridge_pl1_predictability_core_v1 import (
    TIE_TOL,
    joint_delta_distribution,
    target_landscape_metrics,
    weighted_correlation,
    weighted_mean_variance,
)

SCHEMA = "bridge.a21.dual_anchor_geometry.v1"
COMPLETE_STATUS = "BRIDGE_A21_DUAL_ANCHOR_GEOMETRY_FROZEN"
QUALIFIED_ARM = "A2-L"
COMMON_ANCHOR_REACTION = "HEX1"
SECOND_ANCHOR_REACTION = "LDH_L"

A1_COMPARABLE_FEATURES = (
    "sign_magnitude_eta2",
    "sign_magnitude_correlation",
    "directional_entropy3",
    "dominant_sign_mass",
    "n_supported_sign_states",
    "anchor_width80_contraction",
    "joint_ess",
    "abs_delta_anchor_target_correlation",
)

A2_ONLY_SUPPORTIVE_FEATURES = (
    "delta_lactate_target_correlation",
    "abs_delta_lactate_target_correlation",
)


def _finite_1d(values: Sequence[float], *, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1 or arr.size == 0 or not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be a non-empty finite 1-D vector")
    return arr


def lactate_production(ldh_l_flux: Sequence[float]) -> np.ndarray:
    """A2-L candidate observable: max(-LDH_L, 0)."""
    raw = _finite_1d(ldh_l_flux, name="ldh_l_flux")
    return np.maximum(-raw, 0.0)


def target_landscape_metrics_a21(
    *,
    hex1_ct2a: Sequence[float],
    hex1_gl261: Sequence[float],
    lactate_ct2a: Sequence[float],
    lactate_gl261: Sequence[float],
    target_ct2a: Sequence[float],
    target_gl261: Sequence[float],
    strong_weights_ct2a: Sequence[float],
    strong_weights_gl261: Sequence[float],
    tie_tol: float = TIE_TOL,
) -> dict[str, float | str | bool | int]:
    """Compute PL1-compatible geometry plus separately named lactate coupling.

    All original PL1 fields are returned unchanged in meaning. In particular,
    ``abs_delta_anchor_target_correlation`` continues to mean the absolute
    HEX1-contrast versus absolute target-contrast coupling, which preserves the
    A1 comparator semantics. A2-L's second anchor is represented by additional
    lactate-specific fields; no composite dual-anchor coupling score is created.
    """
    hc = _finite_1d(hex1_ct2a, name="hex1_ct2a")
    hg = _finite_1d(hex1_gl261, name="hex1_gl261")
    lc = _finite_1d(lactate_ct2a, name="lactate_ct2a")
    lg = _finite_1d(lactate_gl261, name="lactate_gl261")
    bc = _finite_1d(target_ct2a, name="target_ct2a")
    bg = _finite_1d(target_gl261, name="target_gl261")
    if hc.size != lc.size or hc.size != bc.size:
        raise ValueError("CT2A anchor/target candidate counts must match")
    if hg.size != lg.size or hg.size != bg.size:
        raise ValueError("GL261 anchor/target candidate counts must match")

    out = target_landscape_metrics(
        anchor_ct2a=hc,
        anchor_gl261=hg,
        target_ct2a=bc,
        target_gl261=bg,
        strong_weights_ct2a=strong_weights_ct2a,
        strong_weights_gl261=strong_weights_gl261,
        tie_tol=tie_tol,
    )

    delta_l, joint_l = joint_delta_distribution(
        lc, lg, strong_weights_ct2a, strong_weights_gl261
    )
    delta_b, joint_b = joint_delta_distribution(
        bc, bg, strong_weights_ct2a, strong_weights_gl261
    )
    if not np.allclose(joint_l, joint_b, rtol=0.0, atol=1e-15):
        raise ValueError("lactate/target joint weights differ")

    mean_l, var_l = weighted_mean_variance(delta_l, joint_l)
    corr_l = weighted_correlation(delta_l, delta_b, joint_l, tol=tie_tol)
    corr_abs_l = weighted_correlation(
        np.abs(delta_l), np.abs(delta_b), joint_l, tol=tie_tol
    )

    # Explicit aliases make the preserved PL1 common-anchor semantics auditable.
    out.update(
        {
            "delta_hex1_mean": out["delta_a_mean"],
            "delta_hex1_sd": out["delta_a_sd"],
            "delta_hex1_target_correlation": out["delta_anchor_target_correlation"],
            "abs_delta_hex1_target_correlation": out["abs_delta_anchor_target_correlation"],
            "hex1_delta_degenerate": out["anchor_delta_degenerate"],
            "hex1_target_coupling_status": out["anchor_target_coupling_status"],
            "delta_lactate_mean": mean_l,
            "delta_lactate_sd": math.sqrt(var_l),
            "delta_lactate_target_correlation": corr_l,
            "abs_delta_lactate_target_correlation": corr_abs_l,
            "lactate_delta_degenerate": bool(var_l <= tie_tol),
            "lactate_target_coupling_status": "DEGENERATE"
            if not math.isfinite(corr_l)
            else "OK",
        }
    )
    return out


def analysis_contract() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "status_on_success": COMPLETE_STATUS,
        "qualified_input_arm": QUALIFIED_ARM,
        "strong_anchor_set": [COMMON_ANCHOR_REACTION, SECOND_ANCHOR_REACTION],
        "second_anchor_observable": "lac_production=max(-LDH_L,0)",
        "a2_g_used": False,
        "stage_role": "OUTCOME_BLIND_DUAL_ANCHOR_PREDICTABILITY_LANDSCAPE",
        "primary_scan": {
            "scope": "all_admissible_finite_flux_cache_reactions_except_HEX1_and_LDH_L",
            "preselect_reactions": False,
            "exclude_same_pathway_or_proximal_reactions": False,
            "pathway_relation_is_annotation": True,
            "pathway_exclusion_is_sensitivity_only": True,
        },
        "reference": "uniform_weights_on_the_exact_same_finite_candidate_support",
        "weight_semantics": "A2-L_weights_conditioned_and_normalized_within_exact_method_x_RNA_context",
        "a1_comparable_features": list(A1_COMPARABLE_FEATURES),
        "common_anchor_coupling_semantics": "HEX1_contrast_preserved_exactly_for_A1_A2_comparability",
        "a2_only_supportive_features": list(A2_ONLY_SUPPORTIVE_FEATURES),
        "forbidden": [
            "weak_prior",
            "truth",
            "PL2_outcomes",
            "gene_scores",
            "new_weight_fitting",
            "solver",
            "FVA",
            "sampling",
            "reconstruction",
            "reaction_ranking",
            "composite_dual_anchor_bridgeability_score",
        ],
    }
