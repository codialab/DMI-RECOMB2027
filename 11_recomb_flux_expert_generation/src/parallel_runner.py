"""Process-safe outer parallelism for canonical Stage-11 manifest jobs.

Each manifest row is an independent cache identity, so it is safe to run
multiple rows concurrently.  Workers are spawned (rather than forked) to avoid
inheriting solver/model state from the parent process.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import pandas as pd

from .core import STAGE, load_configuration
from .ensembles import generate_flux_ensemble, valid_ensemble
from .jobs import (
    assert_production_manifests_authorized,
    fit_and_persist_expert,
    parent_reaction_universe,
    stage4_reaction_universe,
)
from .projections import create_projection, valid_projection


_MANIFEST_BY_PHASE = {
    "fit": "expert_jobs.tsv",
    "ensemble": "ensemble_jobs.tsv",
    "projection": "projection_jobs.tsv",
}

# Per-process immutable context populated by _worker_init().
_WORKER_PHASE: str | None = None
_WORKER_PATHS: Mapping[str, Path] | None = None
_WORKER_CONFIG: Mapping[str, Any] | None = None
_WORKER_REACTION_UNIVERSE: Sequence[str] | None = None


def _record_id(phase: str, record: Mapping[str, Any]) -> str:
    return str({
        "fit": record.get("expert_id", record.get("expert_hash", "unknown")),
        "ensemble": record.get("ensemble_hash", "unknown"),
        "projection": record.get("projection_hash", "unknown"),
    }[phase])


def _worker_init(
    phase: str,
    paths: Mapping[str, Path],
    config: Mapping[str, Any],
    reaction_universe: Sequence[str] | None,
) -> None:
    """Initialize one spawned worker without sharing live solver/model objects."""
    global _WORKER_PHASE, _WORKER_PATHS, _WORKER_CONFIG, _WORKER_REACTION_UNIVERSE
    _WORKER_PHASE = phase
    _WORKER_PATHS = paths
    _WORKER_CONFIG = config
    _WORKER_REACTION_UNIVERSE = reaction_universe


def _worker_run(record: Mapping[str, Any]) -> tuple[str, str]:
    """Run one manifest row and return (identity, status)."""
    if _WORKER_PHASE is None or _WORKER_PATHS is None or _WORKER_CONFIG is None:
        raise RuntimeError("Stage-11 worker was not initialized")

    phase = _WORKER_PHASE
    identity = _record_id(phase, record)

    if phase == "fit":
        return identity, str(fit_and_persist_expert(record, _WORKER_PATHS, _WORKER_CONFIG))

    if _WORKER_REACTION_UNIVERSE is None:
        raise RuntimeError(f"Reaction universe was not initialized for phase={phase}")

    if phase == "ensemble":
        was_valid = valid_ensemble(record)
        generate_flux_ensemble(record, _WORKER_REACTION_UNIVERSE, _WORKER_CONFIG)
        return identity, "skipped_valid" if was_valid else "completed"

    if phase == "projection":
        was_valid = valid_projection(record)
        create_projection(record, _WORKER_REACTION_UNIVERSE)
        return identity, "skipped_valid" if was_valid else "completed"

    raise ValueError(f"Unknown Stage-11 phase: {phase}")


def _planned_records(phase: str, limit: int) -> list[dict[str, Any]]:
    manifest = STAGE / "manifests" / _MANIFEST_BY_PHASE[phase]
    jobs = pd.read_csv(manifest, sep="\t")
    if "execution_status" in jobs.columns:
        jobs = jobs.loc[jobs.execution_status.eq("planned")]
    if limit:
        jobs = jobs.head(limit)
    return jobs.to_dict("records")


def _reaction_universe_for_phase(
    phase: str, paths: Mapping[str, Path]
) -> Sequence[str] | None:
    if phase == "ensemble":
        return tuple(parent_reaction_universe(paths))
    if phase == "projection":
        return tuple(stage4_reaction_universe(paths))
    return None


def run_phase(phase: str, *, workers: int = 1, limit: int = 0) -> int:
    """Run one canonical phase; return a shell-style exit code.

    workers=1 preserves serial execution.  workers>1 runs independent manifest
    rows in spawned processes.  Failures are collected while other independent
    jobs finish; the overall return code is non-zero if any row fails.
    """
    if phase not in _MANIFEST_BY_PHASE:
        raise ValueError(f"Unknown Stage-11 phase: {phase}")
    if workers < 1:
        raise ValueError("--workers must be >= 1")
    if limit < 0:
        raise ValueError("--limit must be >= 0")

    # Qualification workers use their dedicated qualification-only entrypoints.
    # The canonical production runner must never trust execution_status=planned
    # from an old manifest without revalidating the current eight-case barrier
    # and the exact final-manifest snapshot authorized after qualification.
    assert_production_manifests_authorized()

    paths, config = load_configuration()
    records = _planned_records(phase, limit)
    if not records:
        print(f"phase={phase} planned_jobs=0 workers=0", flush=True)
        return 0

    worker_count = min(int(workers), len(records))
    inner_sampling_processes = int(config.get("sampling", {}).get("processes", 1))
    if phase == "ensemble" and worker_count > 1 and inner_sampling_processes > 1:
        raise RuntimeError(
            "Refusing nested multiprocessing: outer --workers > 1 while "
            "sampling.processes > 1. Keep sampling.processes: 1 when using "
            "parallel Stage-11 ensemble jobs."
        )

    reaction_universe = _reaction_universe_for_phase(phase, paths)
    gurobi_threads = int(config.get("gurobi", {}).get("Threads", 1))
    print(
        f"phase={phase} planned_jobs={len(records)} workers={worker_count} "
        f"gurobi_threads_per_fit={gurobi_threads} "
        f"sampling_processes_per_ensemble={inner_sampling_processes}",
        flush=True,
    )

    # Serial mode uses the exact same worker function for behavioral parity.
    if worker_count == 1:
        _worker_init(phase, paths, config, reaction_universe)
        failures = 0
        for index, record in enumerate(records, start=1):
            identity = _record_id(phase, record)
            try:
                _, status = _worker_run(record)
                print(f"[{index}/{len(records)}] {identity} {status}", flush=True)
            except Exception as exc:
                failures += 1
                print(
                    f"[{index}/{len(records)}] {identity} FAILED "
                    f"{type(exc).__name__}: {exc}",
                    file=sys.stderr,
                    flush=True,
                )
        if failures:
            print(f"phase={phase} failures={failures}/{len(records)}", file=sys.stderr, flush=True)
            return 1
        return 0

    failures: list[tuple[str, str]] = []
    completed = 0
    # 'spawn' is deliberate: do not fork a process that may already have
    # imported Gurobi/COBRA or initialized native numerical libraries.
    context = get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=worker_count,
        mp_context=context,
        initializer=_worker_init,
        initargs=(phase, paths, config, reaction_universe),
    ) as executor:
        futures = {
            executor.submit(_worker_run, record): _record_id(phase, record)
            for record in records
        }
        for future in as_completed(futures):
            identity = futures[future]
            completed += 1
            try:
                _, status = future.result()
                print(f"[{completed}/{len(records)}] {identity} {status}", flush=True)
            except Exception as exc:
                message = f"{type(exc).__name__}: {exc}"
                failures.append((identity, message))
                print(
                    f"[{completed}/{len(records)}] {identity} FAILED {message}",
                    file=sys.stderr,
                    flush=True,
                )

    if failures:
        print(f"phase={phase} failures={len(failures)}/{len(records)}", file=sys.stderr, flush=True)
        for identity, message in failures:
            print(f"  {identity}: {message}", file=sys.stderr, flush=True)
        return 1

    return 0
