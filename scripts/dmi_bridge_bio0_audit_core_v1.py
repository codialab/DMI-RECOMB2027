#!/usr/bin/env python3
"""DMI-BRIDGE-BIO-0: outcome-free audit primitives for BIO-1 concordance decomposition.

This module defines the scientific contract for the updated DMI master roadmap.
BIO-0 may reconcile candidate identities, existing weight operators, support, scoring
provenance, and gene-score provenance.  It MUST NOT compute new gene--flux
concordance outcomes or use them to choose operators, lambdas, panels, or controls.

The primary BIO-1 decomposition is, for reconstruction/RNA context k,

    p0(k, v) = pi0(k) p0(v | k)

with four matched arms:

    baseline       = pi0(k) p0(v | k)
    within_context = pi0(k) p1(v | k)
    context_only   = pi1(k) p0(v | k)
    full_update    = pi1(k) p1(v | k)

where pi1 and p1(v|k) are derived from the SAME full weak-update distribution.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from typing import Any, Iterable, Mapping, Sequence

SCHEMA = "bridge.bio0.audit.v1"
READY_STATUS = "BIO0_READY_FOR_BIO1_FOUR_ARM_DECOMPOSITION"
BLOCKED_STATUS = "BIO0_BLOCKED_INPUT_OR_SUPPORT"
ROADMAP_REFERENCE = "DMI_SCOPE_SHIFT_BRIDGE_MASTER_ROADMAP_20260928"

PRIMARY_WEAK_OPERATOR = "historical_lactate_direction"
PRIMARY_LAMBDA = 0.5
SENSITIVITY_LAMBDAS = (0.25, 1.0)
STRATUM_FIELDS = ("algorithm", "rna_context_key")

FOUR_ARMS = {
    "strong_anchor_baseline": {
        "context_weights": "baseline",
        "within_context_weights": "baseline",
    },
    "within_context_weak_update": {
        "context_weights": "baseline",
        "within_context_weights": "updated",
    },
    "context_only_update": {
        "context_weights": "updated",
        "within_context_weights": "baseline",
    },
    "full_weak_update": {
        "context_weights": "updated",
        "within_context_weights": "updated",
    },
}

BIO1_FREEZE = {
    "primary_weak_operator": PRIMARY_WEAK_OPERATOR,
    "primary_lambda": PRIMARY_LAMBDA,
    "sensitivity_lambdas": list(SENSITIVITY_LAMBDAS),
    "stratum_fields": list(STRATUM_FIELDS),
    "primary_estimand": "delta_rho_updated_minus_strong_anchor_baseline",
    "concordance_role": "EXPLORATORY_BIOLOGICAL_CONSISTENCY",
    "untouched_confirmation_set": False,
    "controls": ["reversed_direction", "historical_deterministic_randomized_direction"],
    "mouse_uncertainty": {
        "scheme": "paired_within_tumor_mouse_bootstrap",
        "replicates": 10000,
        "seed": 20260918,
        "ct2a_draws": 5,
        "gl261_draws": 5,
        "same_draws_across_arms": True,
        "tumor_aggregation": "median",
        "bootstrap_pathways": False,
    },
    "forbidden_claims": [
        "independent_biological_validation",
        "absolute_flux_accuracy",
        "retrospective_confirmation_of_C4",
        "rescue_of_C5",
    ],
}


@dataclass(frozen=True)
class WeightRow:
    candidate_id: str
    stratum_id: str
    weight: float


@dataclass(frozen=True)
class DecompositionDiagnostics:
    baseline_zero_mass_strata: tuple[str, ...]
    updated_zero_mass_strata: tuple[str, ...]
    nonderivable_within_context_strata: tuple[str, ...]
    support_identical: bool


def canonical_json_bytes(obj: Any) -> bytes:
    return (json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=True) + "\n").encode("ascii")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _finite_nonnegative(value: float, *, name: str) -> float:
    x = float(value)
    if not math.isfinite(x) or x < 0:
        raise ValueError(f"{name} must be finite and non-negative, got {value!r}")
    return x


def normalize_rows(rows: Iterable[WeightRow], *, tol: float = 1e-15) -> list[WeightRow]:
    items = list(rows)
    if not items:
        raise ValueError("empty weight distribution")
    seen: set[str] = set()
    total = 0.0
    for row in items:
        if not row.candidate_id or not row.stratum_id:
            raise ValueError("candidate_id and stratum_id are required")
        if row.candidate_id in seen:
            raise ValueError(f"duplicate candidate_id: {row.candidate_id}")
        seen.add(row.candidate_id)
        total += _finite_nonnegative(row.weight, name=f"weight[{row.candidate_id}]")
    if total <= tol:
        raise ValueError("distribution has no positive mass")
    out = [WeightRow(r.candidate_id, r.stratum_id, r.weight / total) for r in items]
    out.sort(key=lambda r: (r.stratum_id, r.candidate_id))
    return out


def distribution_fingerprint(rows: Sequence[WeightRow]) -> str:
    payload = [
        {"candidate_id": r.candidate_id, "stratum_id": r.stratum_id, "weight": format(r.weight, ".17g")}
        for r in sorted(rows, key=lambda x: (x.stratum_id, x.candidate_id))
    ]
    return sha256_bytes(canonical_json_bytes(payload))


def _as_map(rows: Sequence[WeightRow]) -> dict[str, WeightRow]:
    return {row.candidate_id: row for row in rows}


def _masses(rows: Sequence[WeightRow]) -> dict[str, float]:
    out: dict[str, float] = {}
    for row in rows:
        out[row.stratum_id] = out.get(row.stratum_id, 0.0) + row.weight
    return out


def _conditionals(rows: Sequence[WeightRow], masses: Mapping[str, float], *, tol: float) -> dict[str, float | None]:
    out: dict[str, float | None] = {}
    for row in rows:
        mass = masses[row.stratum_id]
        # A tiny positive mass is still valid support.  The decomposition must
        # preserve the exact finite distribution; using the normalization
        # tolerance here silently discards real (if numerically tiny) strata.
        out[row.candidate_id] = None if mass == 0.0 else row.weight / mass
    return out


def four_arm_decomposition(
    baseline_rows: Iterable[WeightRow],
    full_update_rows: Iterable[WeightRow],
    *,
    tol: float = 1e-15,
) -> tuple[dict[str, list[WeightRow]], DecompositionDiagnostics]:
    """Create the matched four-arm weight decomposition without any outcome access.

    Candidate support and candidate->stratum membership must be identical between the
    baseline and full-update distributions.  If a baseline-positive stratum has zero
    full-update mass, p1(v|k) is undefined and the within-context arm is not
    constructible; fail closed instead of borrowing support or assigning a conditional.
    """
    q = normalize_rows(baseline_rows, tol=tol)
    p = normalize_rows(full_update_rows, tol=tol)
    qmap = _as_map(q)
    pmap = _as_map(p)
    support_identical = set(qmap) == set(pmap)
    if not support_identical:
        missing_q = sorted(set(pmap) - set(qmap))
        missing_p = sorted(set(qmap) - set(pmap))
        raise ValueError(f"candidate support mismatch: only_updated={missing_q[:3]}, only_baseline={missing_p[:3]}")
    for cid in qmap:
        if qmap[cid].stratum_id != pmap[cid].stratum_id:
            raise ValueError(f"stratum mismatch for candidate {cid}")

    qmass = _masses(q)
    pmass = _masses(p)
    strata = sorted(set(qmass) | set(pmass))
    baseline_zero = tuple(k for k in strata if qmass.get(k, 0.0) == 0.0)
    updated_zero = tuple(k for k in strata if pmass.get(k, 0.0) == 0.0)
    nonderivable = tuple(k for k in strata if qmass.get(k, 0.0) > 0.0 and pmass.get(k, 0.0) == 0.0)
    if nonderivable:
        raise ValueError(
            "updated distribution has zero mass in baseline-positive strata; "
            f"within-context update is undefined for {list(nonderivable)}"
        )
    # A full update cannot invent mass in a baseline-zero stratum under a matched-support
    # reweighting experiment.  Such a case is a different support intervention.
    invented = [k for k in strata if qmass.get(k, 0.0) == 0.0 and pmass.get(k, 0.0) > 0.0]
    if invented:
        raise ValueError(f"full update invents support in baseline-zero strata: {invented}")

    qcond = _conditionals(q, qmass, tol=tol)
    pcond = _conditionals(p, pmass, tol=tol)

    baseline: list[WeightRow] = []
    within: list[WeightRow] = []
    context_only: list[WeightRow] = []
    full: list[WeightRow] = []
    for cid in sorted(qmap, key=lambda c: (qmap[c].stratum_id, c)):
        k = qmap[cid].stratum_id
        qc = qcond[cid]
        pc = pcond[cid]
        # q-positive strata are guaranteed p-positive above.
        if qc is None or pc is None:
            # Both q and p have zero mass for this stratum. Candidate rows may exist but
            # all four arms carry zero weight; no conditional is invented.
            b = w = c = f = 0.0
        else:
            b = qmass[k] * qc
            w = qmass[k] * pc
            c = pmass[k] * qc
            f = pmass[k] * pc
        baseline.append(WeightRow(cid, k, b))
        within.append(WeightRow(cid, k, w))
        context_only.append(WeightRow(cid, k, c))
        full.append(WeightRow(cid, k, f))

    arms = {
        "strong_anchor_baseline": normalize_rows(baseline, tol=tol),
        "within_context_weak_update": normalize_rows(within, tol=tol),
        "context_only_update": normalize_rows(context_only, tol=tol),
        "full_weak_update": normalize_rows(full, tol=tol),
    }
    diagnostics = DecompositionDiagnostics(
        baseline_zero_mass_strata=baseline_zero,
        updated_zero_mass_strata=updated_zero,
        nonderivable_within_context_strata=nonderivable,
        support_identical=True,
    )
    return arms, diagnostics


def stratum_masses(rows: Sequence[WeightRow]) -> dict[str, float]:
    return _masses(normalize_rows(rows))


def gene_score_provenance_classification(
    *,
    used_in_model_reconstruction: bool,
    used_in_weak_prior_construction: bool,
    used_in_outcome_guided_selection: bool,
    external_to_flux_model_generation: bool,
) -> str:
    """Return a conservative provenance label; never labels BIO as confirmation."""
    if used_in_model_reconstruction or used_in_weak_prior_construction or used_in_outcome_guided_selection:
        return "OVERLAPPING_NOT_INDEPENDENT"
    if external_to_flux_model_generation:
        return "EXTERNAL_EXPLORATORY_NOT_UNTOUCHED"
    return "PROVENANCE_UNRESOLVED"


def build_analysis_contract(
    *,
    candidate_panel_id: str,
    baseline_operator_id: str,
    weak_operator_id: str,
    reaction_panel_id: str,
    flux_score_id: str,
    gene_score_id: str,
    gene_score_provenance: str,
    controls: Sequence[str],
) -> dict[str, Any]:
    required = {
        "candidate_panel_id": candidate_panel_id,
        "baseline_operator_id": baseline_operator_id,
        "weak_operator_id": weak_operator_id,
        "reaction_panel_id": reaction_panel_id,
        "flux_score_id": flux_score_id,
        "gene_score_id": gene_score_id,
        "gene_score_provenance": gene_score_provenance,
    }
    missing = [key for key, value in required.items() if not str(value).strip()]
    if missing:
        raise ValueError(f"cannot freeze BIO-1 contract with unresolved fields: {missing}")
    controls_set = tuple(sorted(set(controls)))
    for needed in ("reversed_direction", "historical_deterministic_randomized_direction"):
        if needed not in controls_set:
            raise ValueError(f"required control missing: {needed}")
    return {
        "schema_version": SCHEMA,
        "status": READY_STATUS,
        "roadmap_reference": ROADMAP_REFERENCE,
        "candidate_panel_id": candidate_panel_id,
        "baseline_operator_id": baseline_operator_id,
        "weak_operator_id": weak_operator_id,
        "reaction_panel_id": reaction_panel_id,
        "flux_score_id": flux_score_id,
        "gene_score_id": gene_score_id,
        "gene_score_provenance": gene_score_provenance,
        "four_arm_decomposition": FOUR_ARMS,
        "bio1_freeze": BIO1_FREEZE,
        "controls": list(controls_set),
        "outcome_accessed_during_bio0": False,
        "new_solver_or_sampling_invoked": False,
    }
