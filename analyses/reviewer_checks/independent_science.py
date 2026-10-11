"""Independent calculations for the frozen DMI directional-gain endpoints.

This module intentionally does not import the production PL2/A22 numerical
cores.  Level 2A starts from predictions; Level 2B starts from candidate
vectors and already-frozen conditional weights.
"""
from __future__ import annotations

import math
import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import shared_memory
from typing import Sequence

import numpy as np

TIE_TOL = 1e-12
GAIN_TOL = 1e-12
PRIMARY_LAMBDA = 0.25


def _normalized(values: Sequence[float]) -> np.ndarray:
    w = np.asarray(values, dtype=np.float64)
    if w.ndim != 1 or not np.isfinite(w).all() or (w < 0).any():
        raise ValueError("weights must be a finite nonnegative vector")
    total = math.fsum(float(x) for x in w)
    if total <= 0 or not math.isfinite(total):
        raise ValueError("weights must have positive finite mass")
    return w / total


def posterior_magnitude(delta: Sequence[float], weights: Sequence[float],
                        cue_direction: int, lam: float = PRIMARY_LAMBDA) -> float:
    """E[|delta|] under the frozen sign-only exponential tilt."""
    d = np.asarray(delta, dtype=np.float64)
    w = _normalized(weights)
    if d.ndim != 1 or d.size != w.size or not np.isfinite(d).all():
        raise ValueError("delta must be finite and match the weights")
    if cue_direction not in (-1, 1) or not math.isfinite(lam) or lam < 0:
        raise ValueError("cue direction must be +/-1 and lambda finite/nonnegative")
    factors = np.ones(d.size, dtype=np.float64)
    positive = d > TIE_TOL
    negative = d < -TIE_TOL
    factors[positive] = math.exp(lam * cue_direction)
    factors[negative] = math.exp(-lam * cue_direction)
    tilted = w * factors
    z = math.fsum(float(x) for x in tilted)
    if z <= 0 or not math.isfinite(z):
        raise ValueError("posterior normalizer is not positive finite")
    return math.fsum(float(a) * abs(float(b)) for a, b in zip(tilted, d)) / z


def product_delta_distribution(left: Sequence[float], right: Sequence[float],
                               left_weights: Sequence[float],
                               right_weights: Sequence[float]) -> tuple[np.ndarray, np.ndarray]:
    """Construct all left-minus-right states using normalized product weights."""
    a = np.asarray(left, dtype=np.float64)
    b = np.asarray(right, dtype=np.float64)
    wa, wb = _normalized(left_weights), _normalized(right_weights)
    if a.ndim != 1 or b.ndim != 1 or a.size != wa.size or b.size != wb.size:
        raise ValueError("candidate values and weights have inconsistent dimensions")
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("candidate values must be finite")
    return (a[:, None] - b[None, :]).reshape(-1), (wa[:, None] * wb[None, :]).reshape(-1)


def endpoint_from_predictions(truth: float, baseline: float, correct: float,
                              incorrect: float) -> dict[str, float | bool | str]:
    """Level 2A: derive gains from predictions and truth, never gain columns."""
    if not math.isfinite(float(truth)) or not math.isfinite(float(baseline)):
        raise ValueError("truth and baseline prediction must be finite")
    if baseline < 0:
        raise ValueError("baseline magnitude prediction must be nonnegative")
    base_error = abs(abs(truth) - baseline)
    if abs(truth) <= TIE_TOL:
        return {"truth_status": "TRUTH_TIE", "baseline_error": base_error,
                "correct_gain": math.nan, "incorrect_gain": math.nan,
                "information_advantage": math.nan, "directionally_useful": False}
    if not math.isfinite(float(correct)) or not math.isfinite(float(incorrect)):
        raise ValueError("non-tie cue predictions must be finite")
    if min(correct, incorrect) < 0:
        raise ValueError("cue magnitude predictions must be nonnegative")
    correct_error = abs(abs(truth) - correct)
    incorrect_error = abs(abs(truth) - incorrect)
    gc, gi = base_error - correct_error, base_error - incorrect_error
    advantage = 0.5 * (gc - gi)
    return {"truth_status": "NON_TIE", "baseline_error": base_error,
            "correct_gain": gc, "incorrect_gain": gi,
            "information_advantage": advantage,
            "directionally_useful": gc > GAIN_TOL and advantage > GAIN_TOL}


def aggregate_case_predictions(rows):
    """Aggregate case rows after independently deriving gains from predictions.

    Each row must identify one (anchor, evaluation, reaction, truth-pair) case
    and provide truth_delta plus baseline/correct/incorrect magnitude predictions.
    Pair aggregation is the scientific unit used by Figure 3: distinct non-tie
    truth pairs are averaged within each reaction-evaluation pair first.
    """
    import pandas as pd
    frame = pd.DataFrame(rows).copy()
    keys = ["anchor_setting", "evaluation_id", "reaction_id", "truth_pair_id"]
    required = keys + ["truth_delta", "baseline_prediction", "correct_prediction", "incorrect_prediction"]
    missing = [name for name in required if name not in frame]
    if missing:
        raise ValueError(f"case input missing columns: {missing}")
    if frame[keys].isna().any().any() or frame.duplicated(keys).any():
        raise ValueError("case identities must be complete and unique")
    computed = [endpoint_from_predictions(r.truth_delta, r.baseline_prediction,
                                          r.correct_prediction, r.incorrect_prediction)
                for r in frame.itertuples(index=False)]
    frame["truth_status"] = [x["truth_status"] for x in computed]
    for name in ("baseline_error", "correct_gain", "incorrect_gain", "information_advantage"):
        frame[name] = [x.get(name, math.nan) for x in computed]
    frame["directionally_useful"] = [x["directionally_useful"] for x in computed]
    eligible = frame.loc[frame.truth_status.eq("NON_TIE")].copy()
    pair_keys = ["anchor_setting", "evaluation_id", "reaction_id"]
    pair = eligible.groupby(pair_keys, sort=True, as_index=False).agg(
        n_distinct_non_tie_truths=("truth_pair_id", "nunique"),
        correct_gain_mean=("correct_gain", "mean"),
        incorrect_gain_mean=("incorrect_gain", "mean"),
        information_advantage_mean=("information_advantage", "mean"),
        directional_usefulness_fraction=("directionally_useful", "mean"),
    )
    return frame, pair


def exclude_whole_vector(candidate_matrix: np.ndarray, sample_indices: Sequence[int],
                         heldout_index: int, weights: Sequence[float]) -> tuple[np.ndarray, np.ndarray]:
    """Remove one held-out state across every reaction and renormalize its stratum."""
    matrix = np.asarray(candidate_matrix, dtype=np.float64)
    indices = np.asarray(sample_indices, dtype=np.int64)
    w = _normalized(weights)
    if matrix.ndim != 2 or matrix.shape[0] != indices.size or indices.size != w.size:
        raise ValueError("candidate matrix, indices, and weights must align by row")
    if len(set(indices.tolist())) != indices.size:
        raise ValueError("candidate sample indices must be unique")
    keep = indices != int(heldout_index)
    if keep.sum() != indices.size - 1:
        raise ValueError("held-out candidate identity must occur exactly once")
    post = w[keep]
    mass = math.fsum(float(x) for x in post)
    if mass <= 0 or not math.isfinite(mass):
        raise ValueError("held-out exclusion leaves no positive mass")
    return matrix[keep], post / mass


def level2b_predictions(left_vectors: np.ndarray, right_vectors: np.ndarray,
                        left_indices: Sequence[int], right_indices: Sequence[int],
                        heldout_left: int, heldout_right: int,
                        left_weights: Sequence[float], right_weights: Sequence[float],
                        truth_delta: Sequence[float],
                        lam: float = PRIMARY_LAMBDA,
                        reaction_block_size: int = 256) -> dict[str, np.ndarray]:
    """Reconstruct baseline/cue magnitude means conditional on frozen weights.

    Each matrix row is one candidate's complete reaction vector.  Held-out
    states are removed before the Cartesian product is formed.
    """
    left, wl = exclude_whole_vector(left_vectors, left_indices, heldout_left, left_weights)
    right, wr = exclude_whole_vector(right_vectors, right_indices, heldout_right, right_weights)
    truth = np.asarray(truth_delta, dtype=np.float64)
    if left.shape[1] != right.shape[1] or truth.ndim != 1 or truth.size != left.shape[1]:
        raise ValueError("left, right, and truth reaction axes must match")
    if not np.isfinite(truth).all() or not math.isfinite(lam) or lam < 0:
        raise ValueError("truth must be finite and lambda nonnegative")
    if reaction_block_size < 1:
        raise ValueError("reaction_block_size must be positive")
    joint = wl[:, None] * wr[None, :]
    n = truth.size
    baseline = np.empty(n, dtype=np.float64)
    correct = np.empty(n, dtype=np.float64)
    incorrect = np.empty(n, dtype=np.float64)
    truth_sign = np.where(truth > TIE_TOL, 1, np.where(truth < -TIE_TOL, -1, 0))
    for start in range(0, n, reaction_block_size):
        stop = min(start + reaction_block_size, n)
        delta = left[:, None, start:stop] - right[None, :, start:stop]
        magnitude = np.abs(delta)
        baseline[start:stop] = np.sum(joint[:, :, None] * magnitude, axis=(0, 1))
        candidate_sign = np.where(delta > TIE_TOL, 1, np.where(delta < -TIE_TOL, -1, 0))
        aligned = candidate_sign * truth_sign[None, None, start:stop]
        correct_factor = np.where(aligned == 0, 1.0,
                                  np.where(aligned > 0, math.exp(lam), math.exp(-lam)))
        incorrect_factor = np.where(aligned == 0, 1.0,
                                    np.where(aligned > 0, math.exp(-lam), math.exp(lam)))

        def tilted(factor: np.ndarray) -> np.ndarray:
            weighted = joint[:, :, None] * factor
            z = np.sum(weighted, axis=(0, 1))
            if not np.isfinite(z).all() or (z <= 0).any():
                raise ValueError("posterior normalizer is not positive finite")
            return np.sum(weighted * magnitude, axis=(0, 1)) / z

        correct[start:stop] = tilted(correct_factor)
        incorrect[start:stop] = tilted(incorrect_factor)
    return {"baseline": baseline, "correct": correct, "incorrect": incorrect}


def eta2_directional_magnitude(delta: Sequence[float], weights: Sequence[float],
                                tol: float = TIE_TOL) -> float:
    """Weighted direction-conditioned variance fraction on the supplied support."""
    d = np.asarray(delta, dtype=np.float64)
    w = _normalized(weights)
    if d.ndim != 1 or d.size != w.size or not np.isfinite(d).all():
        raise ValueError("delta must be finite and match the weights")
    if not math.isfinite(tol) or tol < 0:
        raise ValueError("tie tolerance must be finite and nonnegative")
    mag = np.abs(d)
    mean = math.fsum(float(a) * float(b) for a, b in zip(w, mag))
    total = math.fsum(float(a) * (float(b) - mean) ** 2 for a, b in zip(w, mag))
    if total <= tol:
        return math.nan
    sign = np.where(d > tol, 1, np.where(d < -tol, -1, 0))
    between = 0.0
    for group in (-1, 0, 1):
        mask = sign == group
        mass = math.fsum(float(x) for x in w[mask])
        if mass == 0.0:
            continue
        conditional_mean = math.fsum(float(a) * float(b) for a, b in zip(w[mask], mag[mask])) / mass
        between += mass * (conditional_mean - mean) ** 2
    return min(1.0, max(0.0, between / total))


def eta2_geometry_block(left: np.ndarray, right: np.ndarray,
                        left_weights: Sequence[float], right_weights: Sequence[float],
                        tie_tol: float = TIE_TOL) -> np.ndarray:
    """Independent 20x20 geometry η² for a reaction block (no truth holdout)."""
    a, b = np.asarray(left, dtype=np.float64), np.asarray(right, dtype=np.float64)
    wa, wb = _normalized(left_weights), _normalized(right_weights)
    if a.ndim != 2 or b.ndim != 2 or a.shape[0] != wa.size or b.shape[0] != wb.size or a.shape[1] != b.shape[1]:
        raise ValueError("candidate matrices and weights must align on their reaction axes")
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("candidate matrices must be finite")
    delta = a[:, None, :] - b[None, :, :]
    weight = (wa[:, None] * wb[None, :])[:, :, None]
    magnitude = np.abs(delta)
    mean = np.sum(weight * magnitude, axis=(0, 1))
    total = np.sum(weight * (magnitude - mean[None, None, :]) ** 2, axis=(0, 1))
    sign = np.where(delta > tie_tol, 1, np.where(delta < -tie_tol, -1, 0))
    between = np.zeros(a.shape[1], dtype=np.float64)
    for group in (-1, 0, 1):
        mask = sign == group
        mass = np.sum(weight * mask, axis=(0, 1))
        numerator = np.sum(weight * magnitude * mask, axis=(0, 1))
        conditional = np.divide(numerator, mass, out=np.zeros_like(numerator), where=mass > 0)
        between += mass * (conditional - mean) ** 2
    output = np.divide(between, total, out=np.full_like(total, np.nan), where=total > tie_tol)
    return np.clip(output, 0.0, 1.0, out=output)


def resolve_workers(workers: int | str, *, memory_budget_mb: int = 2048,
                    estimated_worker_mb: int = 512) -> int:
    """Resolve an explicit count or conservative CPU/memory-aware auto count."""
    cpus = os.cpu_count() or 1
    if workers == "auto":
        return max(1, min(cpus, memory_budget_mb // max(1, estimated_worker_mb)))
    try:
        count = int(workers)
    except (TypeError, ValueError) as exc:
        raise ValueError("workers must be a positive integer or 'auto'") from exc
    if count < 1:
        raise ValueError("workers must be a positive integer")
    return min(count, cpus)


_GEOM_MATS = None
_GEOM_EVALS = None
_GEOM_WEIGHTS = None
_GEOM_RXNS = None
_GEOM_SHM = None
_THREAD_LIMITER = None


def _geometry_worker_init(shm_name: str, mats_shape: tuple[int, ...], mats_dtype: str,
                          evals: np.ndarray, weights: np.ndarray, reactions: np.ndarray) -> None:
    global _GEOM_MATS, _GEOM_EVALS, _GEOM_WEIGHTS, _GEOM_RXNS, _GEOM_SHM, _THREAD_LIMITER
    _GEOM_SHM = shared_memory.SharedMemory(name=shm_name)
    _GEOM_MATS = np.ndarray(mats_shape, dtype=np.dtype(mats_dtype), buffer=_GEOM_SHM.buf)
    _GEOM_MATS.flags.writeable = False
    _GEOM_EVALS, _GEOM_WEIGHTS, _GEOM_RXNS = evals, weights, reactions
    try:
        from threadpoolctl import threadpool_limits
        _THREAD_LIMITER = threadpool_limits(limits=1)
    except ImportError:
        # Geometry uses NumPy reductions rather than threaded BLAS operations.
        _THREAD_LIMITER = None


def _geometry_worker_batch(task: tuple[int, ...]) -> list[tuple[int, np.ndarray, np.ndarray]]:
    if _GEOM_MATS is None:
        raise RuntimeError("geometry worker was not initialized")
    results = []
    for eval_index in task:
        left_idx, right_idx = _GEOM_EVALS[eval_index]
        left = _GEOM_MATS[left_idx][:, _GEOM_RXNS]
        right = _GEOM_MATS[right_idx][:, _GEOM_RXNS]
        w = _GEOM_WEIGHTS[eval_index]
        a1 = eta2_geometry_block(left, right, w[0], w[1])
        a2 = eta2_geometry_block(left, right, w[2], w[3])
        results.append((eval_index, a1, a2))
    return results


def reconstruct_geometry_parallel(mats: np.ndarray, evaluation_ensemble_indices: np.ndarray,
                                  evaluation_weights: np.ndarray, reaction_indices: np.ndarray,
                                  workers: int | str = 1, batch_size: int = 4,
                                  memory_budget_mb: int = 2048) -> list[tuple[int, np.ndarray, np.ndarray]]:
    """Recompute η² for independent evaluations, sharing the candidate cache read-only.

    ``evaluation_weights`` has shape (n_evaluations, 4, 20), ordered A1-CT2A,
    A1-GL261, A2-CT2A, A2-GL261.  Geometry uses the full 20x20 support.
    """
    global _GEOM_MATS, _GEOM_EVALS, _GEOM_WEIGHTS, _GEOM_RXNS
    data = np.asarray(mats)
    evals = np.asarray(evaluation_ensemble_indices, dtype=np.int64)
    weights = np.asarray(evaluation_weights, dtype=np.float64)
    reactions = np.asarray(reaction_indices, dtype=np.int64)
    if data.ndim != 3 or evals.ndim != 2 or evals.shape[1] != 2:
        raise ValueError("candidate cache/evaluation ensemble indices have invalid shapes")
    if weights.shape != (evals.shape[0], 4, data.shape[1]) or reactions.ndim != 1:
        raise ValueError("evaluation weights or reaction indices have invalid shapes")
    if batch_size < 1 or reactions.size == 0 or not np.isfinite(data).all():
        raise ValueError("batch size, reaction support, or candidate values are invalid")
    count = resolve_workers(workers, memory_budget_mb=memory_budget_mb)
    if count == 1:
        old = (_GEOM_MATS, _GEOM_EVALS, _GEOM_WEIGHTS, _GEOM_RXNS)
        try:
            _GEOM_MATS, _GEOM_EVALS, _GEOM_WEIGHTS, _GEOM_RXNS = data, evals, weights, reactions
            return _geometry_worker_batch(tuple(range(evals.shape[0])))
        finally:
            _GEOM_MATS, _GEOM_EVALS, _GEOM_WEIGHTS, _GEOM_RXNS = old
    shm = shared_memory.SharedMemory(create=True, size=data.nbytes)
    try:
        shared = np.ndarray(data.shape, dtype=data.dtype, buffer=shm.buf)
        shared[:] = data
        shared.flags.writeable = False
        tasks = [tuple(range(i, min(i + batch_size, evals.shape[0])))
                 for i in range(0, evals.shape[0], batch_size)]
        context = multiprocessing.get_context("spawn")
        try:
            pool = ProcessPoolExecutor(
                max_workers=count, mp_context=context,
                initializer=_geometry_worker_init,
                initargs=(shm.name, data.shape, data.dtype.str, evals, weights, reactions),
            )
        except OSError:
            # Pool creation can be unavailable in constrained hosts.  The
            # deterministic serial path is a safe fallback for that condition.
            old = (_GEOM_MATS, _GEOM_EVALS, _GEOM_WEIGHTS, _GEOM_RXNS)
            try:
                _GEOM_MATS, _GEOM_EVALS, _GEOM_WEIGHTS, _GEOM_RXNS = data, evals, weights, reactions
                return _geometry_worker_batch(tuple(range(evals.shape[0])))
            finally:
                _GEOM_MATS, _GEOM_EVALS, _GEOM_WEIGHTS, _GEOM_RXNS = old
        with pool:
            output = []
            try:
                for block in pool.map(_geometry_worker_batch, tasks, chunksize=1):
                    output.extend(block)
            except Exception as exc:
                raise RuntimeError(f"geometry worker pool failed: {exc}") from exc
            return output
    finally:
        shm.close()
        shm.unlink()
