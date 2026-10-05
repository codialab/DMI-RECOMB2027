#!/usr/bin/env python3
"""DMI-BRIDGE PL2D-H support-gate amendment primitives.

This module amends only the per-context support floor before confirmation opening.
It contains no code for reading confirmation predictors or outcomes.
"""
from __future__ import annotations

import math
from fractions import Fraction
from typing import Mapping, Sequence

SCHEMA = "bridge.pl2d_h.support_gate_amendment.v1"
STATUS = "PL2D_H_SUPPORT_GATE_AMENDED_READY_FOR_CONFIRMATION"
ORIGINAL_STATUS = "PL2D_HYPOTHESIS_FROZEN_READY_FOR_CONFIRMATION"
ORIGINAL_MANIFEST_SHA256 = "a1dd76d6b40e741e07efe64c34bc39eaad289a46bf6c870917a4f59cab688427"
ORIGINAL_ARTIFACT_SHA256 = {
    "BRIDGEPL2DH_SOURCE_AUDIT.json": "21f8ea38f331922034a869a75864b33da72455a66d7e403725cb16016d7a05c9",
    "BRIDGEPL2DH_DEVELOPMENT_SNAPSHOT.json": "ebd60bbdd480ec83fe233e068a398776b4bd4ff964cc4bc6bf8c615c66cc28b8",
    "BRIDGEPL2DH_CONFIRMATION_REGISTRY.tsv": "02a38afe3cb25b8c85a7774c1f20118ad9601ba74f0df0fbff021ba2494bd2e5",
    "BRIDGEPL2DH_CONFIRMATION_CONTRACT.json": "ac7fa325491ab96402fc774bc3c3f3843a884fa711c8bf1104eb056390bb2054",
}

DEVELOPMENT_REACTIONS = 3135
CONFIRMATION_REACTIONS = 1045
CONTEXT_COUNT = 16
ORIGINAL_MINIMUM_FINITE_REACTIONS_PER_CONTEXT = 100
AMENDED_MINIMUM_FINITE_REACTIONS_PER_CONTEXT = 30
RETENTION_NUMERATOR = 7
RETENTION_DENOMINATOR = 10

FROZEN_DEVELOPMENT_CONTEXT_SUPPORT = (
    ("CORDA", "training_samples=setx1,setx2", 1847),
    ("CORDA", "training_samples=setx1,setx2,setx3", 1815),
    ("CORDA", "training_samples=setx1,setx3", 1854),
    ("CORDA", "training_samples=setx2,setx3", 1841),
    ("GIMME", "training_samples=setx1,setx2", 1235),
    ("GIMME", "training_samples=setx1,setx2,setx3", 1264),
    ("GIMME", "training_samples=setx1,setx3", 1237),
    ("GIMME", "training_samples=setx2,setx3", 1283),
    ("RIPTiDe", "training_samples=setx1,setx2", 129),
    ("RIPTiDe", "training_samples=setx1,setx2,setx3", 159),
    ("RIPTiDe", "training_samples=setx1,setx3", 134),
    ("RIPTiDe", "training_samples=setx2,setx3", 162),
    ("iMAT", "training_samples=setx1,setx2", 1729),
    ("iMAT", "training_samples=setx1,setx2,setx3", 1737),
    ("iMAT", "training_samples=setx1,setx3", 1701),
    ("iMAT", "training_samples=setx2,setx3", 1726),
)

UNCHANGED_PRIMARY = {
    "unit": "reaction",
    "predictor": "mean_sign_magnitude_eta2_across_evaluations",
    "response": "mean_directionally_useful_fraction_across_defined_evaluations",
    "expected_direction": "POSITIVE",
    "minimum_finite_reactions": 700,
    "minimum_spearman_rho": 0.50,
    "context_count": 16,
    "minimum_positive_contexts": 12,
}


def expected_proportional_support(min_development_support: int) -> Fraction:
    """Exact proportional support under the frozen 3:1 development/holdout split."""
    if min_development_support < 0:
        raise ValueError("development support must be non-negative")
    return Fraction(min_development_support * CONFIRMATION_REACTIONS, DEVELOPMENT_REACTIONS)


def derive_amended_context_floor(context_support: Sequence[int]) -> dict[str, object]:
    """Derive the frozen 30-pair floor from development support only."""
    if len(context_support) != CONTEXT_COUNT:
        raise ValueError(f"expected {CONTEXT_COUNT} development context counts")
    counts = [int(v) for v in context_support]
    if any(v < 0 or v > DEVELOPMENT_REACTIONS for v in counts):
        raise ValueError("invalid development context support")
    minimum = min(counts)
    expected = expected_proportional_support(minimum)
    retained = expected * Fraction(RETENTION_NUMERATOR, RETENTION_DENOMINATOR)
    floor = retained.numerator // retained.denominator
    return {
        "minimum_development_context_support": minimum,
        "proportional_confirmation_support_numerator": expected.numerator,
        "proportional_confirmation_support_denominator": expected.denominator,
        "proportional_confirmation_support": float(expected),
        "retention_fraction_numerator": RETENTION_NUMERATOR,
        "retention_fraction_denominator": RETENTION_DENOMINATOR,
        "derived_context_floor": floor,
    }


def verify_frozen_development_contexts(rows: Sequence[Mapping[str, object]]) -> list[int]:
    """Verify the original development snapshot's context support inventory."""
    observed = []
    for row in rows:
        observed.append((str(row["algorithm"]), str(row["rna_context_key"]), int(row["finite_reactions"])))
    observed.sort()
    expected = sorted(FROZEN_DEVELOPMENT_CONTEXT_SUPPORT)
    if observed != expected:
        raise ValueError("development context support does not match frozen PL2D-H snapshot")
    return [count for _, _, count in FROZEN_DEVELOPMENT_CONTEXT_SUPPORT]


def amended_primary(original_primary: Mapping[str, object]) -> dict[str, object]:
    """Return the effective primary contract with exactly one changed numeric field."""
    result = dict(original_primary)
    required = dict(UNCHANGED_PRIMARY)
    for key, expected in required.items():
        if result.get(key) != expected:
            raise ValueError(f"original primary field changed unexpectedly: {key}")
    if result.get("minimum_finite_reactions_per_context") != ORIGINAL_MINIMUM_FINITE_REACTIONS_PER_CONTEXT:
        raise ValueError("original per-context support floor is not the qualified value 100")
    result["minimum_finite_reactions_per_context"] = AMENDED_MINIMUM_FINITE_REACTIONS_PER_CONTEXT
    return result


def classify_confirmation(
    *,
    finite_reactions: int,
    overall_rho: float,
    context_finite_counts: Sequence[int],
    positive_contexts: int,
) -> str:
    """Apply the amended support gate and otherwise unchanged primary effect gates."""
    if int(finite_reactions) < UNCHANGED_PRIMARY["minimum_finite_reactions"]:
        return "PL2D_INSUFFICIENT_SUPPORT"
    if len(context_finite_counts) != UNCHANGED_PRIMARY["context_count"]:
        return "PL2D_INSUFFICIENT_SUPPORT"
    if any(int(v) < AMENDED_MINIMUM_FINITE_REACTIONS_PER_CONTEXT for v in context_finite_counts):
        return "PL2D_INSUFFICIENT_SUPPORT"
    if not math.isfinite(float(overall_rho)):
        return "PL2D_INSUFFICIENT_SUPPORT"
    if float(overall_rho) < UNCHANGED_PRIMARY["minimum_spearman_rho"]:
        return "PL2D_NOT_CONFIRMED"
    if int(positive_contexts) < UNCHANGED_PRIMARY["minimum_positive_contexts"]:
        return "PL2D_NOT_CONFIRMED"
    return "PL2D_CONFIRMED"
