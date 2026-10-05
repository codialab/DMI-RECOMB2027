"""Deterministic derived projections of a persisted full-parent q_k(v)."""
from __future__ import annotations
import json
import os
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence
import numpy as np
import pandas as pd
from .cache_identity import PROJECTION_CACHE_SCHEMA
from .core import STAGE, atomic_json, atomic_table, file_sha256, reaction_universe_hash
from .ensembles import valid_ensemble
from .expert_artifacts import marker_valid


def projection_directory(record: Mapping[str, Any]) -> Path:
    return STAGE / "flux_ensembles" / str(record["algorithm"]) / str(record["ensemble_hash"]) / "projections" / str(record["projection_name"]) / str(record["projection_hash"])


def _source_path(record: Mapping[str, Any]) -> Path:
    return STAGE / "flux_ensembles" / str(record["algorithm"]) / str(record["ensemble_hash"]) / "full_parent_fluxes/reaction_samples.parquet"


def valid_projection(record: Mapping[str, Any]) -> bool:
    directory = projection_directory(record)
    if not marker_valid(directory, {"ensemble_hash": str(record["ensemble_hash"]), "projection_hash": str(record["projection_hash"])}, {"reaction_samples.parquet", "metadata.json"}, PROJECTION_CACHE_SCHEMA): return False
    try:
        metadata = json.loads((directory / "metadata.json").read_text()); projected = pd.read_parquet(directory / "reaction_samples.parquet"); source = _source_path(record)
        ids = list(map(str, metadata["ordered_projection_reactions"])); full = pd.read_parquet(source, columns=["sample_index", *ids])
        return (metadata.get("cache_schema") == PROJECTION_CACHE_SCHEMA and metadata.get("source_full_ensemble_checksum") == file_sha256(source)
            and list(projected.columns) == ["sample_index", *ids] and projected.equals(full)
            and projected.iloc[:, 1:].dtypes.map(lambda dtype: np.dtype(dtype) == np.dtype("float32")).all()
            and metadata.get("reaction_universe_sha256") == reaction_universe_hash(ids))
    except Exception: return False


def create_projection(record: Mapping[str, Any], reaction_ids: Sequence[str]) -> Path:
    if not valid_ensemble(record): raise RuntimeError("Valid full-parent ensemble is required before projection")
    if valid_projection(record): return projection_directory(record)
    ids = list(map(str, reaction_ids))
    if len(ids) != len(set(ids)): raise ValueError("Projection reaction universe has duplicates")
    final = projection_directory(record); final.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{record['projection_hash'][:16]}.", dir=final.parent))
    try:
        source = _source_path(record); stored = pd.read_parquet(source)
        if list(stored.columns[:1]) != ["sample_index"]: raise RuntimeError("Full ensemble lacks sample_index")
        missing = set(ids) - set(stored.columns)
        if missing: raise RuntimeError(f"Projection reactions absent from full parent ensemble: {sorted(missing)[:5]}")
        projected = stored.loc[:, ["sample_index", *ids]].copy()
        if not projected.iloc[:, 1:].dtypes.map(lambda dtype: np.dtype(dtype) == np.dtype("float32")).all(): raise RuntimeError("Projection source must remain float32")
        atomic_table(temporary / "reaction_samples.parquet", projected); source_checksum = file_sha256(source)
        atomic_json(temporary / "metadata.json", dict(record) | {"cache_schema": PROJECTION_CACHE_SCHEMA, "source_full_ensemble_checksum": source_checksum,
            "reaction_universe_sha256": reaction_universe_hash(ids), "ordered_projection_reactions": ids, "reaction_count": len(ids),
            "storage_dtype": "float32", "calculation_dtype": "float64", "created_utc": datetime.now(timezone.utc).isoformat()})
        checksums = {str(p.relative_to(temporary)): file_sha256(p) for p in temporary.rglob("*") if p.is_file()}
        atomic_json(temporary / "DONE.json", {"cache_schema": PROJECTION_CACHE_SCHEMA, "ensemble_hash": record["ensemble_hash"], "projection_hash": record["projection_hash"], "qc_passed": True, "checksums": checksums})
        if final.exists():
            target = STAGE / "qc/stale_projections" / f"{final.name}.{time.time_ns()}"; target.parent.mkdir(parents=True, exist_ok=True); os.replace(final, target)
        os.replace(temporary, final); return final
    except Exception:
        target = STAGE / "qc/failed_projections" / f"{record['projection_hash']}.{time.time_ns()}"; target.parent.mkdir(parents=True, exist_ok=True); os.replace(temporary, target); raise
