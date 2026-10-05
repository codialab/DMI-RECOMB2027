"""Pure numerical helpers for DMI-BRIDGE-A2.0.

This module contains no repository discovery and no scientific-outcome access.
The repository adapter is responsible for verifying source identities and the
archived Stage-12 operator before calling these helpers.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class KernelResult:
    weights: np.ndarray
    temperature: float
    achieved_ess: float


def average_ranks(values: Sequence[float]) -> np.ndarray:
    """Return one-based average ranks with stable tie handling."""
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or x.size == 0 or not np.all(np.isfinite(x)):
        raise ValueError("values must be a non-empty finite 1-D array")

    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(x.size, dtype=float)

    i = 0
    while i < x.size:
        j = i + 1
        while j < x.size and x[order[j]] == x[order[i]]:
            j += 1
        avg = 0.5 * ((i + 1) + j)
        ranks[order[i:j]] = avg
        i = j
    return ranks


def candidate_rank_coordinate(values: Sequence[float]) -> np.ndarray:
    """Archived candidate convention: (average_rank - 0.5) / n."""
    r = average_ranks(values)
    n = r.size
    return (r - 0.5) / float(n)


def mouse_rank_coordinate(values: Sequence[float]) -> np.ndarray:
    """Archived mouse convention: (average_rank - 1) / (n - 1)."""
    r = average_ranks(values)
    n = r.size
    if n < 2:
        raise ValueError("at least two mouse values are required")
    return (r - 1.0) / float(n - 1)


def dual_squared_distance(
    candidate_coord_1: Sequence[float],
    candidate_coord_2: Sequence[float],
    mouse_coord_1: float,
    mouse_coord_2: float,
) -> np.ndarray:
    """Equal-weight 2-D squared rank distance."""
    q1 = np.asarray(candidate_coord_1, dtype=float)
    q2 = np.asarray(candidate_coord_2, dtype=float)
    if q1.shape != q2.shape or q1.ndim != 1 or q1.size == 0:
        raise ValueError("candidate coordinates must be non-empty 1-D arrays of equal length")
    if not (np.all(np.isfinite(q1)) and np.all(np.isfinite(q2))):
        raise ValueError("candidate coordinates must be finite")
    if not (math.isfinite(mouse_coord_1) and math.isfinite(mouse_coord_2)):
        raise ValueError("mouse coordinates must be finite")
    return (q1 - float(mouse_coord_1)) ** 2 + (q2 - float(mouse_coord_2)) ** 2


def effective_sample_size(weights: Sequence[float]) -> float:
    w = np.asarray(weights, dtype=float)
    if w.ndim != 1 or w.size == 0 or np.any(w < 0) or not np.all(np.isfinite(w)):
        raise ValueError("weights must be a non-empty finite non-negative 1-D array")
    s = float(w.sum())
    if s <= 0:
        raise ValueError("weights must have positive sum")
    p = w / s
    return 1.0 / float(np.dot(p, p))


def kernel_weights(distance: Sequence[float], temperature: float) -> np.ndarray:
    d = np.asarray(distance, dtype=float)
    if d.ndim != 1 or d.size == 0 or not np.all(np.isfinite(d)):
        raise ValueError("distance must be a non-empty finite 1-D array")
    if not math.isfinite(temperature) or temperature <= 0:
        raise ValueError("temperature must be finite and positive")

    z = -(d - float(np.min(d))) / float(temperature)
    z -= float(np.max(z))
    w = np.exp(z)
    s = float(w.sum())
    if not math.isfinite(s) or s <= 0:
        raise ValueError("kernel normalization failed")
    return w / s


def solve_temperature_for_ess(
    distance: Sequence[float],
    target_ess: float = 20.0,
    *,
    low: float = 1e-12,
    high: float = 1.0,
    max_high: float = 1e6,
    iterations: int = 80,
) -> KernelResult:
    """Solve using the verified archived Stage-12 temperature procedure.

    Start at [1e-12, 1], expand the upper bracket by factors of ten until the
    target ESS is bracketed (maximum 1e6), run 80 geometric-mean bisection
    steps, and return the weights evaluated at the final upper bracket.
    """
    d = np.asarray(distance, dtype=float)
    if d.ndim != 1 or d.size == 0 or not np.all(np.isfinite(d)):
        raise ValueError("distance must be a non-empty finite 1-D array")
    if not (1.0 <= target_ess <= float(d.size)):
        raise ValueError("target_ess must lie in [1, n]")
    if not (0 < low < high <= max_high) or iterations <= 0:
        raise ValueError("invalid temperature search contract")

    lo = float(low)
    hi = float(high)
    cap = float(max_high)

    ess_lo = effective_sample_size(kernel_weights(d, lo))
    ess_hi = effective_sample_size(kernel_weights(d, hi))
    if target_ess < ess_lo - 1e-10:
        raise ValueError("target ESS is not bracketed by temperature bounds")

    while ess_hi < target_ess - 1e-10 and hi < cap:
        hi = min(10.0 * hi, cap)
        ess_hi = effective_sample_size(kernel_weights(d, hi))

    if target_ess > ess_hi + 1e-10:
        raise ValueError("target ESS is not bracketed by temperature bounds")

    for _ in range(iterations):
        mid = math.sqrt(lo * hi)
        ess_mid = effective_sample_size(kernel_weights(d, mid))
        if ess_mid < target_ess:
            lo = mid
        else:
            hi = mid

    t = hi
    w = kernel_weights(d, t)
    return KernelResult(
        weights=w,
        temperature=t,
        achieved_ess=effective_sample_size(w),
    )
