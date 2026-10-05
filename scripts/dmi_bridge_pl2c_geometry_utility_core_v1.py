#!/usr/bin/env python3
"""Frozen numerical primitives for DMI-BRIDGE-PL2C development analysis.

This module contains only outcome definitions, the outcome-blind reaction split,
equal-weight distinct-pair aggregation, and deterministic descriptive utilities.
Repository admission, PL1/PL2B joins, and publication live in the PL2C adapter.
"""
from __future__ import annotations

from collections.abc import Iterable, Sequence
import hashlib
import math
from typing import Any

import numpy as np

SCHEMA = "bridge.pl2c.geometry_utility_development.v1"
COMPLETE_STATUS = "PL2C_DEVELOPMENT_GEOMETRY_UTILITY_COMPLETE"
PL2D_STATUS = "UNTOUCHED"

PL2A_MANIFEST_SHA256 = "8ae3eb60eb23614e379e38e3aa3cb72ec16c6bdb69aa9c74cc94499e31f1c312"
PL2B_MANIFEST_SHA256 = "7d9ca97a3f4ea5aa8b0b983e09ee83fc3c82f326d1b1cbf9c6137cfadf25539d"
PL2B_CASE_LOGICAL_SHA256 = "c3670992a0ad33698af6c7620429bed48e0ecaa98adb3c87e648eefb1ba70363"
PL1_FEATURE_LOGICAL_SHA256 = "46aa9f0410f89d03b020e01fefa517348faae1d55d9bfe8387ac04b39bb9ac93"

PRIMARY_LAMBDA = 0.25
GAIN_TOL = 1e-12
REACTION_SPLIT_SALT = "DMI-BRIDGE-PL2C-REACTION-HOLDOUT-V1"
EXPECTED_REACTIONS = 4180
EXPECTED_DEVELOPMENT_REACTIONS = 3135
EXPECTED_CONFIRMATION_REACTIONS = 1045

PRIMARY_FEATURES = (
    "sign_magnitude_eta2",
    "sign_magnitude_correlation",
    "directional_entropy3",
    "dominant_sign_mass",
    "n_supported_sign_states",
    "anchor_width80_contraction",
    "joint_ess",
    "abs_delta_anchor_target_correlation",
)


def _finite(value: float, *, name: str) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite")
    return out


def reaction_split_hash(reaction_id: str) -> str:
    rid = str(reaction_id)
    if not rid:
        raise ValueError("reaction_id must be non-empty")
    return hashlib.sha256((REACTION_SPLIT_SALT + "\0" + rid).encode("utf-8")).hexdigest()


def freeze_reaction_split(reaction_ids: Sequence[str]) -> list[dict[str, Any]]:
    """Assign exactly one quarter of frozen reactions to confirmation, outcome-blind.

    Assignment depends only on reaction identity and is invariant to caller order.
    For the frozen 4,180-reaction PL1 registry this gives 1,045 confirmation and
    3,135 development reactions.
    """
    ids = [str(value) for value in reaction_ids]
    if not ids or any(not value for value in ids):
        raise ValueError("reaction_ids must be non-empty strings")
    if len(set(ids)) != len(ids):
        raise ValueError("reaction_ids must be unique")
    ranked = sorted((reaction_split_hash(rid), rid) for rid in ids)
    n_confirmation = len(ranked) // 4
    confirmation = {rid for _, rid in ranked[:n_confirmation]}
    records = [
        {
            "reaction_id": rid,
            "split_hash": digest,
            "split_rank": rank,
            "analysis_split": "CONFIRMATION_HOLDOUT" if rid in confirmation else "DEVELOPMENT",
        }
        for rank, (digest, rid) in enumerate(ranked)
    ]
    return records


def split_registry_sha256(records: Sequence[dict[str, Any]]) -> str:
    canonical = "reaction_id\tsplit_hash\tsplit_rank\tanalysis_split\n" + "".join(
        f"{row['reaction_id']}\t{row['split_hash']}\t{int(row['split_rank'])}\t{row['analysis_split']}\n"
        for row in records
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def classify_signed(value: float, *, tol: float = GAIN_TOL) -> str:
    if not math.isfinite(float(value)):
        return "UNDEFINED"
    if not math.isfinite(float(tol)) or tol < 0.0:
        raise ValueError("tol must be finite and non-negative")
    if value > tol:
        return "POSITIVE"
    if value < -tol:
        return "NEGATIVE"
    return "TIED"


def case_responses(
    *,
    correct_gain: float,
    wrong_gain: float,
    random_gain: float,
    truth_status: str,
    tol: float = GAIN_TOL,
    identity_tol: float = 5e-12,
) -> dict[str, Any]:
    """Return frozen PL2C responses for one distinct PL2B truth-pair case."""
    if truth_status == "TRUTH_TIE":
        return {
            "response_status": "TRUTH_TIE",
            "correct_gain": math.nan,
            "random_gain": math.nan,
            "information_advantage": math.nan,
            "correct_gain_status": "UNDEFINED",
            "information_advantage_status": "UNDEFINED",
            "directionally_useful": False,
        }
    if truth_status != "NON_TIE":
        raise ValueError(f"unexpected truth_status: {truth_status}")
    gc = _finite(correct_gain, name="correct_gain")
    gw = _finite(wrong_gain, name="wrong_gain")
    gr = _finite(random_gain, name="random_gain")
    expected_random = 0.5 * (gc + gw)
    if abs(gr - expected_random) > identity_tol:
        raise ValueError("random-sign gain identity failed")
    advantage = gc - gr
    half_endpoint_gap = 0.5 * (gc - gw)
    if abs(advantage - half_endpoint_gap) > identity_tol:
        raise ValueError("information-advantage identity failed")
    correct_status = classify_signed(gc, tol=tol)
    advantage_status = classify_signed(advantage, tol=tol)
    return {
        "response_status": "DEFINED",
        "correct_gain": gc,
        "random_gain": gr,
        "information_advantage": advantage,
        "correct_gain_status": correct_status,
        "information_advantage_status": advantage_status,
        "directionally_useful": correct_status == "POSITIVE" and advantage_status == "POSITIVE",
    }


def _median(values: Sequence[float]) -> float:
    return float(np.median(np.asarray(values, dtype=float)))


def aggregate_distinct_pair_responses(records: Iterable[dict[str, Any]]) -> dict[str, Any]:
    """Aggregate one evaluation/reaction across distinct truth pairs with equal weight."""
    rows = list(records)
    if not rows:
        raise ValueError("at least one distinct-pair row is required")
    pair_ids = [str(row["truth_pair_id"]) for row in rows]
    if len(set(pair_ids)) != len(pair_ids):
        raise ValueError("truth_pair_id must be unique within an evaluation/reaction aggregate")

    defined: list[dict[str, Any]] = []
    tie_count = 0
    for row in rows:
        response = case_responses(
            correct_gain=float(row.get("correct_gain", math.nan)),
            wrong_gain=float(row.get("wrong_gain", math.nan)),
            random_gain=float(row.get("random_gain", math.nan)),
            truth_status=str(row["truth_status"]),
        )
        if response["response_status"] == "TRUTH_TIE":
            tie_count += 1
        else:
            defined.append(response)

    out: dict[str, Any] = {
        "n_distinct_pairs": len(rows),
        "n_non_tie_pairs": len(defined),
        "n_truth_tie_pairs": tie_count,
        "non_tie_pair_fraction": len(defined) / float(len(rows)),
    }
    if not defined:
        out.update({
            "aggregate_status": "NO_DIRECTIONAL_TRUTH",
            "mean_correct_gain": math.nan,
            "median_correct_gain": math.nan,
            "mean_information_advantage": math.nan,
            "median_information_advantage": math.nan,
            "positive_correct_gain_fraction": math.nan,
            "positive_information_advantage_fraction": math.nan,
            "directionally_useful_fraction": math.nan,
        })
        return out

    correct = [float(row["correct_gain"]) for row in defined]
    advantage = [float(row["information_advantage"]) for row in defined]
    n = float(len(defined))
    out.update({
        "aggregate_status": "DEFINED",
        "mean_correct_gain": float(np.mean(correct)),
        "median_correct_gain": _median(correct),
        "mean_information_advantage": float(np.mean(advantage)),
        "median_information_advantage": _median(advantage),
        "positive_correct_gain_fraction": sum(value > GAIN_TOL for value in correct) / n,
        "positive_information_advantage_fraction": sum(value > GAIN_TOL for value in advantage) / n,
        "directionally_useful_fraction": sum(bool(row["directionally_useful"]) for row in defined) / n,
    })
    return out


def _average_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.size, dtype=float)
    i = 0
    while i < values.size:
        j = i + 1
        while j < values.size and values[order[j]] == values[order[i]]:
            j += 1
        average = 0.5 * ((i + 1) + j)
        ranks[order[i:j]] = average
        i = j
    return ranks


def spearman_finite(x: Sequence[float], y: Sequence[float]) -> dict[str, float | int]:
    """Deterministic descriptive Spearman correlation using finite paired values only."""
    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    if xa.shape != ya.shape or xa.ndim != 1:
        raise ValueError("x and y must be same-length 1-D vectors")
    mask = np.isfinite(xa) & np.isfinite(ya)
    xv = xa[mask]
    yv = ya[mask]
    n = int(xv.size)
    if n < 2:
        return {"rho": math.nan, "finite_count": n, "nonfinite_count": int(xa.size - n)}
    xr = _average_ranks(xv)
    yr = _average_ranks(yv)
    xsd = float(np.std(xr))
    ysd = float(np.std(yr))
    rho = math.nan if xsd == 0.0 or ysd == 0.0 else float(np.corrcoef(xr, yr)[0, 1])
    return {"rho": rho, "finite_count": n, "nonfinite_count": int(xa.size - n)}


def linear_quartile_edges(values: Sequence[float]) -> tuple[float, float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        raise ValueError("quartiles require at least one finite value")
    q = np.quantile(arr, [0.25, 0.5, 0.75], method="linear")
    return float(q[0]), float(q[1]), float(q[2])


def quartile_label(value: float, edges: Sequence[float]) -> str:
    x = float(value)
    if not math.isfinite(x):
        return "NONFINITE"
    q25, q50, q75 = (float(v) for v in edges)
    if x <= q25:
        return "Q1"
    if x <= q50:
        return "Q2"
    if x <= q75:
        return "Q3"
    return "Q4"


def analysis_contract() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "status_on_success": COMPLETE_STATUS,
        "stage_role": "DEVELOPMENT_ONLY_FROZEN_PL1_GEOMETRY_TO_FROZEN_PL2B_UTILITY",
        "source_identity": {
            "pl2a_manifest_sha256": PL2A_MANIFEST_SHA256,
            "pl2b_manifest_sha256": PL2B_MANIFEST_SHA256,
            "pl2b_case_logical_sha256": PL2B_CASE_LOGICAL_SHA256,
            "pl1_feature_logical_sha256": PL1_FEATURE_LOGICAL_SHA256,
        },
        "reaction_split": {
            "salt": REACTION_SPLIT_SALT,
            "method": "sha256_rank_then_first_quarter_confirmation",
            "expected_total": EXPECTED_REACTIONS,
            "expected_development": EXPECTED_DEVELOPMENT_REACTIONS,
            "expected_confirmation": EXPECTED_CONFIRMATION_REACTIONS,
            "outcomes_must_not_define_split": True,
            "confirmation_outcomes_must_not_be_analyzed_in_pl2c": True,
        },
        "primary_lambda": PRIMARY_LAMBDA,
        "gain_tolerance": GAIN_TOL,
        "primary_response": "directionally_useful_fraction_after_distinct_pair_equal_weight_aggregation",
        "required_coverage_diagnostic": "non_tie_pair_fraction",
        "primary_continuous_secondary_response": "mean_information_advantage",
        "information_advantage": "g_correct_minus_g_random_equals_half_correct_minus_wrong_gap",
        "directionally_useful": "g_correct_gt_tol_AND_information_advantage_gt_tol",
        "predictor_unit": "evaluation_id_x_reaction_id",
        "truth_pair_aggregation": "equal_weight_distinct_truth_pairs_never_alias_weighted",
        "primary_features": list(PRIMARY_FEATURES),
        "missingness": "preserve_nonfinite_PL1_values_and_report_finite_denominators_no_zero_imputation",
        "permitted": [
            "development_only_univariate_spearman",
            "algorithm_x_rna_context_stability_summaries",
            "outcome_blind_linear_quartile_bins",
            "fixed_eta2_x_entropy_4x4_map",
            "labeled_pathway_proximity_sensitivity",
        ],
        "forbidden": [
            "confirmation_outcome_analysis",
            "reaction_ranking",
            "outcome_driven_feature_selection",
            "outcome_driven_threshold_optimization",
            "lambda_tuning",
            "unrestricted_machine_learning",
            "biological_p_values_from_reaction_or_mouse_cross_pseudoreplication",
            "solver_or_optimization",
            "FVA",
            "new_flux_sampling",
            "reconstruction",
            "new_strong_anchor_fitting",
            "new_flux_weights",
            "PL1_PL2A_or_PL2B_mutation",
        ],
        "pl2d_confirmation_status_on_success": PL2D_STATUS,
    }
