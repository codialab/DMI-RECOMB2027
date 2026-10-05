"""Pure helpers for DMI-BRIDGE-A2.2 matched sign-only utility.

A2.2 deliberately delegates the weak-prior numerical operator to the frozen PL2
core.  This module adds only the dual-anchor holdout/support and common-target
contracts needed to run the same operator on frozen A2-L weights.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np

try:
    from scripts import dmi_bridge_pl2_sign_only_core_v1 as pl2
except ModuleNotFoundError:  # direct script execution
    import dmi_bridge_pl2_sign_only_core_v1 as pl2

SCHEMA = "bridge.a22.dual_anchor_sign_only_utility.core.v1"
STRONG_ANCHORS = ("HEX1", "LDH_L")
EXPECTED_CACHE_REACTIONS = 4181
EXPECTED_TARGET_REACTIONS = 4179
EXPECTED_MATCHED_TRUTH_PAIRS = 836
EXPECTED_ALL_DISTINCT_TRUTH_PAIRS = 885
EXPECTED_TRUTH_PAIR_ALIASES = 1200
EXPECTED_CASE_ROWS = EXPECTED_MATCHED_TRUTH_PAIRS * EXPECTED_TARGET_REACTIONS


@dataclass(frozen=True)
class HoldoutResult:
    weights: np.ndarray
    pre_holdout_mass: float
    post_holdout_mass: float
    positive_weight_candidates_pre: int
    positive_weight_candidates_post: int
    pre_holdout_ess: float
    post_holdout_ess: float


def frozen_pl2_semantics() -> dict[str, object]:
    """Return the exact PL2 numerical constants A2.2 is allowed to use."""
    return {
        "tie_tol": pl2.TIE_TOL,
        "gain_tol": pl2.GAIN_TOL,
        "lambda_grid": tuple(pl2.LAMBDA_GRID),
        "primary_lambda": pl2.PRIMARY_LAMBDA,
        "reliability_grid": tuple(pl2.RELIABILITY_GRID),
        "random_sign_null": "q=0.5 exact expectation of correct/wrong endpoint gains",
    }


def validate_frozen_pl2_semantics() -> None:
    expected = {
        "tie_tol": 1e-12,
        "gain_tol": 1e-12,
        "lambda_grid": (0.0, 0.25, 0.5, 1.0),
        "primary_lambda": 0.25,
        "reliability_grid": (0.0, 0.25, 0.5, 0.75, 1.0),
    }
    observed = frozen_pl2_semantics()
    for key, value in expected.items():
        if observed[key] != value:
            raise RuntimeError(f"frozen PL2 semantic mismatch: {key}")


def common_target_reactions(reaction_ids: Sequence[str]) -> list[str]:
    """Return the exact A1/A2 common B universe, excluding both strong anchors."""
    ids = [str(value) for value in reaction_ids]
    if len(ids) != EXPECTED_CACHE_REACTIONS or len(set(ids)) != len(ids):
        raise ValueError("reaction inventory must contain 4,181 unique reactions")
    for anchor in STRONG_ANCHORS:
        if ids.count(anchor) != 1:
            raise ValueError(f"strong anchor identity is not unique: {anchor}")
    targets = [reaction for reaction in ids if reaction not in STRONG_ANCHORS]
    if len(targets) != EXPECTED_TARGET_REACTIONS:
        raise RuntimeError("A2.2 target reaction count mismatch")
    return targets


def remove_holdout_and_normalize(
    weights: Sequence[float], holdout_index: int
) -> HoldoutResult:
    """Remove one frozen truth candidate and renormalize the remaining support."""
    w = np.asarray(weights, dtype=float)
    if w.ndim != 1 or w.size < 2 or not np.all(np.isfinite(w)):
        raise ValueError("weights must be a finite one-dimensional array with >=2 entries")
    if np.any(w < 0.0):
        raise ValueError("weights must be non-negative")
    if holdout_index < 0 or holdout_index >= w.size:
        raise IndexError("holdout_index is outside the candidate vector")
    total = float(np.sum(w))
    if not math.isfinite(total) or total <= 0.0:
        raise ValueError("weights must have positive finite mass")
    p = w / total
    positive_pre = int(np.count_nonzero(p > 0.0))
    pre_ess = float(1.0 / np.sum(p * p))
    kept = np.delete(p, holdout_index)
    # Summing retained mass avoids cancellation when the truth has nearly all
    # weight; it is the mass actually used in the 19-candidate renormalization.
    post_mass = float(np.sum(kept))
    positive_post = int(np.count_nonzero(kept > 0.0))
    if not math.isfinite(post_mass) or post_mass <= 0.0:
        raise ValueError("zero post-holdout mass under frozen A2-L weights")
    post = kept / post_mass
    if not np.all(np.isfinite(post)) or np.any(post < 0.0):
        raise ValueError("invalid post-holdout normalized weights")
    if not math.isclose(float(np.sum(post)), 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("post-holdout weights do not normalize to one")
    post_ess = float(1.0 / np.sum(post * post))
    return HoldoutResult(
        weights=post,
        pre_holdout_mass=1.0,
        post_holdout_mass=post_mass,
        positive_weight_candidates_pre=positive_pre,
        positive_weight_candidates_post=positive_post,
        pre_holdout_ess=pre_ess,
        post_holdout_ess=post_ess,
    )


def evaluate_truth_case(**kwargs):
    """Delegate exactly to the frozen PL2 sign-only truth-case operator."""
    validate_frozen_pl2_semantics()
    return pl2.evaluate_truth_case(**kwargs)


def analysis_contract() -> dict[str, object]:
    validate_frozen_pl2_semantics()
    return {
        "schema": "bridge.a22.dual_anchor_sign_only_utility.v1",
        "stage_role": "POST_FREEZE_MATCHED_SIGN_ONLY_UTILITY_ON_A2L",
        "qualified_baseline": "A2-L",
        "strong_anchors": list(STRONG_ANCHORS),
        "target_reaction_count": EXPECTED_TARGET_REACTIONS,
        "matched_truth_pair_count": EXPECTED_MATCHED_TRUTH_PAIRS,
        "all_frozen_truth_pair_count": EXPECTED_ALL_DISTINCT_TRUTH_PAIRS,
        "truth_pair_alias_count": EXPECTED_TRUTH_PAIR_ALIASES,
        "expected_case_rows_if_support_gate_passes": EXPECTED_CASE_ROWS,
        "truth_identity_policy": "reuse original PL2A identities; no A2-L reselection",
        "primary_population": "original A1-evaluable 836 distinct truth pairs only",
        "a1_non_evaluable_policy": "audit support only; no A2.2 v1 reaction-B outcomes",
        "post_holdout_policy": "remove exact frozen truth candidate then renormalize A2-L stratum weights",
        "new_ess_threshold": False,
        "tie_tol": pl2.TIE_TOL,
        "gain_tol": pl2.GAIN_TOL,
        "lambda_grid": list(pl2.LAMBDA_GRID),
        "primary_lambda": pl2.PRIMARY_LAMBDA,
        "reliability_grid": list(pl2.RELIABILITY_GRID),
        "random_sign_control": "q=0.5 exact analytic expectation; no Monte Carlo",
        "primary_target": "absolute_between_condition_difference_magnitude_E_abs_delta_B",
        "a2_g_used": False,
        "geometry_utility_analysis": False,
        "a1_a2_synthesis": False,
        "untouched_confirmation_claim": False,
    }
