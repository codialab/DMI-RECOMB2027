#!/usr/bin/env python3
"""Solver-free qualification helpers for DMI-BRIDGE-BIO-0R.

BIO-0R is an additive provenance qualification layer.  It must not calculate
new gene--flux concordance outcomes, alter frozen BIO-0 four-arm weights, or
select BIO-1 operators/parameters from biological outcomes.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from typing import Iterable, Mapping, Sequence

import numpy as np

SCHEMA = "bridge.bio0r.qualification.v1"
READY_STATUS = "BIO0R_QUALIFIED_READY_FOR_BIO1"
BLOCKED_STATUS = "BIO0R_BLOCKED_PROVENANCE_OR_REPRODUCTION"
PRIMARY_LAMBDA = 0.5
SENSITIVITY_LAMBDAS = (0.25, 1.0)
RANDOM_SEEDS = tuple(range(20260926, 20260942))
TOL = 1e-12


def canonical_json_bytes(obj: object) -> bytes:
    return (json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=True) + "\n").encode("ascii")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _finite(values: Sequence[float], *, name: str) -> np.ndarray:
    a = np.asarray(values, dtype=float)
    if a.ndim != 1 or len(a) == 0 or not np.isfinite(a).all():
        raise ValueError(f"{name} must be a non-empty finite 1-D vector")
    return a


def average_ranks(values: Sequence[float]) -> np.ndarray:
    """One-based average ranks with deterministic tie handling."""
    x = _finite(values, name="values")
    order = np.argsort(x, kind="stable")
    out = np.empty(len(x), dtype=float)
    i = 0
    while i < len(x):
        j = i + 1
        while j < len(x) and x[order[j]] == x[order[i]]:
            j += 1
        rank = ((i + 1) + j) / 2.0
        out[order[i:j]] = rank
        i = j
    return out


def candidate_midrank_coordinate(values: Sequence[float]) -> np.ndarray:
    """Historical candidate coordinate: (average_rank - 0.5) / N."""
    x = _finite(values, name="candidate values")
    return (average_ranks(x) - 0.5) / float(len(x))


def mouse_rank_coordinate(values: Sequence[float]) -> np.ndarray:
    """Historical mouse coordinate: (average_rank - 1) / (N - 1)."""
    x = _finite(values, name="mouse values")
    if len(x) < 2:
        raise ValueError("mouse rank coordinate requires at least two values")
    return (average_ranks(x) - 1.0) / float(len(x) - 1)


def normalize(weights: Sequence[float]) -> np.ndarray:
    w = _finite(weights, name="weights")
    if (w < 0).any() or float(w.sum()) <= 0.0:
        raise ValueError("weights must be non-negative with positive mass")
    return w / float(w.sum())


def effective_sample_size(weights: Sequence[float]) -> float:
    q = normalize(weights)
    return float(1.0 / np.square(q).sum())


def rank_kernel_weights(
    candidate_coordinate: Sequence[float],
    mouse_coordinate: float,
    temperature: float,
) -> np.ndarray:
    c = _finite(candidate_coordinate, name="candidate_coordinate")
    m = float(mouse_coordinate)
    t = float(temperature)
    if not math.isfinite(m) or not math.isfinite(t) or t <= 0.0:
        raise ValueError("mouse_coordinate must be finite and temperature positive")
    dist = np.square(c - m)
    z = -(dist - float(dist.min())) / t
    z -= float(z.max())
    return normalize(np.exp(z))


def assert_close_vector(actual: Sequence[float], expected: Sequence[float], *, atol: float = TOL) -> None:
    a = _finite(actual, name="actual")
    e = _finite(expected, name="expected")
    if a.shape != e.shape or not np.allclose(a, e, rtol=0.0, atol=atol):
        diff = float(np.max(np.abs(a - e))) if a.shape == e.shape else math.inf
        raise ValueError(f"vector mismatch; max_abs_diff={diff}")


def identity_overlap(old_ids: Iterable[str], current_ids: Iterable[str]) -> dict[str, int]:
    old = set(map(str, old_ids))
    cur = set(map(str, current_ids))
    return {
        "historical": len(old),
        "current": len(cur),
        "overlap": len(old & cur),
        "historical_only": len(old - cur),
        "current_only": len(cur - old),
    }


def weak_lactate_score(q_vlac: Sequence[float], tumor: str) -> np.ndarray:
    q = _finite(q_vlac, name="q_Vlac")
    if tumor == "CT2A":
        d = 1.0
    elif tumor == "GL261":
        d = -1.0
    else:
        raise ValueError(f"unknown tumor: {tumor}")
    return d * (2.0 * q - 1.0)


def tilted_weights(base_weights: Sequence[float], score: Sequence[float], lam: float) -> np.ndarray:
    q = normalize(base_weights)
    s = _finite(score, name="score")
    if len(q) != len(s):
        raise ValueError("base_weights and score length mismatch")
    lam = float(lam)
    if not math.isfinite(lam) or lam < 0:
        raise ValueError("lambda must be finite and non-negative")
    if lam == 0.0:
        return q.copy()
    z = lam * s
    z -= float(z.max())
    return normalize(q * np.exp(z))


def random_sign_vector(seed: int, n: int) -> np.ndarray:
    if n <= 0:
        raise ValueError("n must be positive")
    rng = np.random.default_rng(int(seed))
    return rng.choice(np.asarray([-1, 1], dtype=np.int8), size=int(n), replace=True)


def random_sign_fingerprint(signs: Sequence[int]) -> str:
    s = np.asarray(signs, dtype=np.int8)
    if s.ndim != 1 or len(s) == 0 or not np.isin(s, [-1, 1]).all():
        raise ValueError("random signs must be a non-empty vector in {-1,+1}")
    payload = {"schema": "bridge.bio0r.random_signs.v1", "signs": [int(x) for x in s]}
    return sha256_bytes(canonical_json_bytes(payload))


def support_signature(weights: Sequence[float], strata: Sequence[str]) -> dict[str, object]:
    q = _finite(weights, name="weights")
    if len(q) != len(strata):
        raise ValueError("weights/strata length mismatch")
    labels = np.asarray([str(x) for x in strata], dtype=object)
    masses: dict[str, float] = {}
    for k in sorted(set(labels.tolist())):
        masses[k] = float(q[labels == k].sum())
    return {
        "positive_candidate_count": int(np.count_nonzero(q > 0.0)),
        "zero_candidate_count": int(np.count_nonzero(q == 0.0)),
        "positive_strata": int(sum(v > 0.0 for v in masses.values())),
        "zero_strata": int(sum(v == 0.0 for v in masses.values())),
        "stratum_masses": masses,
    }


def assert_matched_support(base_weights: Sequence[float], updated_weights: Sequence[float], strata: Sequence[str]) -> None:
    b = _finite(base_weights, name="base_weights")
    u = _finite(updated_weights, name="updated_weights")
    if b.shape != u.shape:
        raise ValueError("baseline/update shape mismatch")
    bs = support_signature(b, strata)
    us = support_signature(u, strata)
    if not np.array_equal(b > 0.0, u > 0.0):
        raise ValueError("candidate positive-weight support differs")
    if bs["positive_strata"] != us["positive_strata"] or bs["zero_strata"] != us["zero_strata"]:
        raise ValueError("stratum support differs")
    if any((bs["stratum_masses"][k] == 0.0) != (us["stratum_masses"][k] == 0.0) for k in bs["stratum_masses"]):
        raise ValueError("literal-zero stratum support differs")


@dataclass(frozen=True)
class QualificationSummary:
    historical_overlap_verified: bool
    current_vmax_semantics_verified: bool
    lactate_semantics_verified: bool
    random_family_verified: bool
    all_operator_support_verified: bool
    baseline_full_reproduction_verified: bool
    ensemble_projection_cache_identity_verified: bool
    historical_reproduction_contract_verified: bool

    @property
    def ready(self) -> bool:
        return all(self.__dict__.values())

    @property
    def status(self) -> str:
        return READY_STATUS if self.ready else BLOCKED_STATUS
