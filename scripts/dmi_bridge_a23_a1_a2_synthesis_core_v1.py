"""Pure helpers for DMI-BRIDGE-A2.3 matched A1/A2 synthesis.

No repository discovery, outcome selection, or scientific production lives here.
The repository adapter must first recover and bind the frozen A1 PL2C/PL2D
estimator contract.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np


def require_unique_keys(keys: Iterable[tuple], *, label: str) -> list[tuple]:
    out = list(keys)
    if len(out) != len(set(out)):
        raise ValueError(f"{label} contains duplicate keys")
    return out


def require_exact_key_match(
    left_keys: Iterable[tuple],
    right_keys: Iterable[tuple],
    *,
    left_label: str = "left",
    right_label: str = "right",
) -> list[tuple]:
    left = require_unique_keys(left_keys, label=left_label)
    right = require_unique_keys(right_keys, label=right_label)
    if set(left) != set(right):
        missing_right = len(set(left) - set(right))
        missing_left = len(set(right) - set(left))
        raise ValueError(
            f"key mismatch: missing from {right_label}={missing_right}, "
            f"missing from {left_label}={missing_left}"
        )
    return sorted(left)


def paired_delta(a1: Sequence[float], a2: Sequence[float]) -> np.ndarray:
    x = np.asarray(a1, dtype=float)
    y = np.asarray(a2, dtype=float)
    if x.shape != y.shape or x.ndim != 1:
        raise ValueError("a1 and a2 must be equal-length 1-D arrays")
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        raise ValueError("paired values must be finite")
    return y - x


def average_ranks(values: Sequence[float]) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or x.size == 0 or not np.all(np.isfinite(x)):
        raise ValueError("values must be a non-empty finite 1-D array")
    order = np.argsort(x, kind="mergesort")
    out = np.empty(x.size, dtype=float)
    i = 0
    while i < x.size:
        j = i + 1
        while j < x.size and x[order[j]] == x[order[i]]:
            j += 1
        out[order[i:j]] = 0.5 * ((i + 1) + j)
        i = j
    return out


def spearman_rho(x: Sequence[float], y: Sequence[float]) -> float:
    a = np.asarray(x, dtype=float)
    b = np.asarray(y, dtype=float)
    if a.shape != b.shape or a.ndim != 1 or a.size < 2:
        raise ValueError("x and y must be equal-length 1-D arrays with n>=2")
    mask = np.isfinite(a) & np.isfinite(b)
    a = a[mask]
    b = b[mask]
    if a.size < 2:
        return float("nan")
    ra = average_ranks(a)
    rb = average_ranks(b)
    ra = ra - float(np.mean(ra))
    rb = rb - float(np.mean(rb))
    den = float(np.sqrt(np.dot(ra, ra) * np.dot(rb, rb)))
    if den == 0:
        return float("nan")
    return float(np.dot(ra, rb) / den)


def directional_specificity(
    correct_gain: Sequence[float],
    wrong_gain: Sequence[float],
) -> np.ndarray:
    """Return paired correct-minus-wrong gain on exact matched cases."""
    c = np.asarray(correct_gain, dtype=float)
    w = np.asarray(wrong_gain, dtype=float)
    if c.shape != w.shape or c.ndim != 1:
        raise ValueError("gain arrays must be equal-length 1-D arrays")
    if not (np.all(np.isfinite(c)) and np.all(np.isfinite(w))):
        raise ValueError("gain arrays must be finite")
    return c - w
