#!/usr/bin/env python3
"""Frozen primitives for DMI-BRIDGE PL2D one-shot confirmation."""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

try:
    from scripts import dmi_bridge_pl2d_hypothesis_freeze_core_v1 as hcore
except ModuleNotFoundError:
    import dmi_bridge_pl2d_hypothesis_freeze_core_v1 as hcore

SCHEMA = "bridge.pl2d.confirmation.v1"
EXPECTED_CONFIRMATION_REACTIONS = 1045
EXPECTED_DEVELOPMENT_REACTIONS = 3135
EXPECTED_EVALUATIONS = 370
EXPECTED_CONTEXTS = 16
EXPECTED_EVALUATION_REACTION_ROWS = EXPECTED_CONFIRMATION_REACTIONS * EXPECTED_EVALUATIONS
EXPECTED_EVALUABLE_TRUTH_PAIRS = 836
EXPECTED_CONFIRMATION_CASE_ROWS = EXPECTED_CONFIRMATION_REACTIONS * EXPECTED_EVALUABLE_TRUTH_PAIRS
EXPECTED_DEVELOPMENT_CASE_ROWS_TO_SKIP = EXPECTED_DEVELOPMENT_REACTIONS * EXPECTED_EVALUABLE_TRUTH_PAIRS

PL2A_MANIFEST_SHA256 = "8ae3eb60eb23614e379e38e3aa3cb72ec16c6bdb69aa9c74cc94499e31f1c312"
PL2B_MANIFEST_SHA256 = "7d9ca97a3f4ea5aa8b0b983e09ee83fc3c82f326d1b1cbf9c6137cfadf25539d"
PL2B_CASE_LOGICAL_SHA256 = "c3670992a0ad33698af6c7620429bed48e0ecaa98adb3c87e648eefb1ba70363"
PL1_FEATURE_LOGICAL_SHA256 = "46aa9f0410f89d03b020e01fefa517348faae1d55d9bfe8387ac04b39bb9ac93"
PL2C_MANIFEST_SHA256 = "6abcf149cd6adb294c0931585cec28f992b453451a5b5f03e10b3ffb2f9b9c39"
PL2DH_MANIFEST_SHA256 = "a1dd76d6b40e741e07efe64c34bc39eaad289a46bf6c870917a4f59cab688427"
PL2DH_CONFIRMATION_REGISTRY_SHA256 = "02a38afe3cb25b8c85a7774c1f20118ad9601ba74f0df0fbff021ba2494bd2e5"
PL2DHA_MANIFEST_SHA256 = "fb7df9e208da2e6c9841d43a3b572038e59943a4221535705eef9f5dc0ed6a3a"
PL2DHA_CONTRACT_SHA256 = "08baccb983fa7690d1e10dd965b95d42f0223a798e797dbf0b0c1a830491cfb4"

GAIN_TOL = 1e-12
PRIMARY_MINIMUM_FINITE_REACTIONS = 700
PRIMARY_MINIMUM_FINITE_PER_CONTEXT = 30
PRIMARY_MINIMUM_RHO = 0.50
PRIMARY_MINIMUM_POSITIVE_CONTEXTS = 12
BOOTSTRAP_REPLICATES = 5000
BOOTSTRAP_SEED = 20260929
BOOTSTRAP_Q = (0.025, 0.975)

PRIMARY_FEATURE = "sign_magnitude_eta2"
SUPPORTIVE_FEATURE_DIRECTIONS = {
    "directional_entropy3": "POSITIVE",
    "dominant_sign_mass": "NEGATIVE",
    "n_supported_sign_states": "POSITIVE",
}


@dataclass(frozen=True)
class PrimaryResult:
    status: str
    finite_reactions: int
    overall_rho: float
    context_count: int
    context_min_finite_reactions: int
    positive_contexts: int


def classify_confirmation(
    *,
    finite_reactions: int,
    overall_rho: float,
    context_finite_counts: Sequence[int],
    context_rhos: Sequence[float],
) -> PrimaryResult:
    """Apply only the frozen amended PL2D primary gate."""
    counts = [int(v) for v in context_finite_counts]
    rhos = [float(v) for v in context_rhos]
    positive = sum(math.isfinite(v) and v > 0.0 for v in rhos)
    min_count = min(counts) if counts else 0
    insufficient = (
        int(finite_reactions) < PRIMARY_MINIMUM_FINITE_REACTIONS
        or len(counts) != EXPECTED_CONTEXTS
        or len(rhos) != EXPECTED_CONTEXTS
        or any(v < PRIMARY_MINIMUM_FINITE_PER_CONTEXT for v in counts)
        or any(not math.isfinite(v) for v in rhos)
        or not math.isfinite(float(overall_rho))
    )
    if insufficient:
        status = "PL2D_INSUFFICIENT_SUPPORT"
    elif float(overall_rho) < PRIMARY_MINIMUM_RHO or positive < PRIMARY_MINIMUM_POSITIVE_CONTEXTS:
        status = "PL2D_NOT_CONFIRMED"
    else:
        status = "PL2D_CONFIRMED"
    return PrimaryResult(
        status=status,
        finite_reactions=int(finite_reactions),
        overall_rho=float(overall_rho),
        context_count=len(counts),
        context_min_finite_reactions=min_count,
        positive_contexts=positive,
    )


def primary_association(reaction_rows: Sequence[Mapping[str, object]]) -> tuple[int, float]:
    pairs = []
    for row in reaction_rows:
        x = row.get("mean_sign_magnitude_eta2")
        y = row.get("mean_directionally_useful_fraction")
        if x is not None and y is not None and math.isfinite(float(x)) and math.isfinite(float(y)):
            pairs.append((float(x), float(y)))
    if not pairs:
        return 0, math.nan
    return len(pairs), hcore.spearman([a for a, _ in pairs], [b for _, b in pairs])


def partial_primary_controlling_coverage(reaction_rows: Sequence[Mapping[str, object]]) -> tuple[int, float]:
    triples = []
    for row in reaction_rows:
        values = (
            row.get("mean_sign_magnitude_eta2"),
            row.get("mean_directionally_useful_fraction"),
            row.get("mean_non_tie_pair_fraction"),
        )
        if all(v is not None and math.isfinite(float(v)) for v in values):
            triples.append(tuple(float(v) for v in values))
    if not triples:
        return 0, math.nan
    return len(triples), hcore.partial_spearman_one_control(
        [v[0] for v in triples], [v[1] for v in triples], [v[2] for v in triples]
    )


def bootstrap_primary(
    reaction_rows: Sequence[Mapping[str, object]],
    *,
    replicates: int = BOOTSTRAP_REPLICATES,
    seed: int = BOOTSTRAP_SEED,
) -> dict[str, object]:
    """Frozen supportive reaction bootstrap over all confirmation reaction rows."""
    if len(reaction_rows) != EXPECTED_CONFIRMATION_REACTIONS:
        raise ValueError(f"expected {EXPECTED_CONFIRMATION_REACTIONS} confirmation reactions")
    if int(replicates) <= 0:
        raise ValueError("replicates must be positive")
    x = np.asarray([
        math.nan if r.get("mean_sign_magnitude_eta2") is None else float(r["mean_sign_magnitude_eta2"])
        for r in reaction_rows
    ], dtype=float)
    y = np.asarray([
        math.nan if r.get("mean_directionally_useful_fraction") is None else float(r["mean_directionally_useful_fraction"])
        for r in reaction_rows
    ], dtype=float)
    rng = np.random.Generator(np.random.PCG64(int(seed)))
    rho = np.full(int(replicates), math.nan, dtype=float)
    n_finite = np.zeros(int(replicates), dtype=int)
    for i in range(int(replicates)):
        idx = rng.integers(0, len(reaction_rows), size=len(reaction_rows), endpoint=False)
        xv = x[idx]
        yv = y[idx]
        mask = np.isfinite(xv) & np.isfinite(yv)
        n_finite[i] = int(mask.sum())
        if n_finite[i] >= 2:
            rho[i] = hcore.spearman(xv[mask], yv[mask])
    finite = rho[np.isfinite(rho)]
    if finite.size:
        lo, hi = np.quantile(finite, BOOTSTRAP_Q, method="linear")
        median = np.quantile(finite, 0.5, method="linear")
    else:
        lo = hi = median = math.nan
    return {
        "replicates": int(replicates),
        "seed": int(seed),
        "rng": "numpy.random.Generator(PCG64)",
        "rho": rho,
        "finite_pair_counts": n_finite,
        "finite_replicates": int(finite.size),
        "nonfinite_replicates": int(replicates - finite.size),
        "rho_q025": float(lo),
        "rho_median": float(median),
        "rho_q975": float(hi),
    }


def contract_constants() -> dict[str, object]:
    return {
        "schema": SCHEMA,
        "confirmation_reactions": EXPECTED_CONFIRMATION_REACTIONS,
        "confirmation_evaluation_reaction_rows": EXPECTED_EVALUATION_REACTION_ROWS,
        "confirmation_case_rows": EXPECTED_CONFIRMATION_CASE_ROWS,
        "development_case_rows_skipped_before_outcome_parse": EXPECTED_DEVELOPMENT_CASE_ROWS_TO_SKIP,
        "primary_feature": PRIMARY_FEATURE,
        "primary_minimum_finite_reactions": PRIMARY_MINIMUM_FINITE_REACTIONS,
        "primary_minimum_finite_reactions_per_context": PRIMARY_MINIMUM_FINITE_PER_CONTEXT,
        "primary_minimum_rho": PRIMARY_MINIMUM_RHO,
        "primary_minimum_positive_contexts": PRIMARY_MINIMUM_POSITIVE_CONTEXTS,
        "context_count": EXPECTED_CONTEXTS,
        "gain_tolerance": GAIN_TOL,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "supportive_feature_directions": dict(SUPPORTIVE_FEATURE_DIRECTIONS),
    }
