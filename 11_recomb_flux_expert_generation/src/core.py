"""Stable IO, hashing, GPR, and flux QC reused from Stage-10 ``src/core.py``.

The functions below preserve the validated semantics while Stage 11 changes
the artifact hierarchy from contexts to RNA-derived flux experts.
"""
from __future__ import annotations

import ast
import hashlib
import json
import lzma
import math
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
STAGE = Path(__file__).resolve().parents[1]


def _is_parquet_name(path: Path) -> bool:
    text = str(path)
    return text.endswith(".parquet") or text.endswith(".parquet.xz")


def _parquet_base_and_xz(path: Path) -> tuple[Path, Path]:
    text = str(path)
    if text.endswith(".parquet.xz"):
        return Path(text[:-3]), path
    if text.endswith(".parquet"):
        return path, Path(text + ".xz")
    raise ValueError(f"Not a Parquet path: {path}")


def resolve_parquet_path(path: str | os.PathLike[str] | Path, *, prefer_xz: bool = True) -> Path:
    """Resolve either ``foo.parquet`` or ``foo.parquet.xz``.

    When both exist, prefer the XZ-wrapped representation by default.  Callers
    may continue passing the historical ``.parquet`` path during migration.
    """
    requested = Path(path)
    if not _is_parquet_name(requested):
        return requested
    base, xz = _parquet_base_and_xz(requested)
    order = (xz, base) if prefer_xz else (base, xz)
    for candidate in order:
        if candidate.is_file():
            return candidate
    return requested


def parquet_exists(path: str | os.PathLike[str] | Path) -> bool:
    requested = Path(path)
    if not _is_parquet_name(requested):
        return requested.is_file()
    base, xz = _parquet_base_and_xz(requested)
    return base.is_file() or xz.is_file()


_PANDAS_READ_PARQUET = pd.read_parquet


def read_parquet_compat(path, *args, **kwargs):
    """Read ordinary Parquet or an outer XZ-wrapped Parquet file.

    ``.parquet.xz`` is an XZ container around the complete Parquet byte stream;
    it is not a native Parquet codec.  XZ input is therefore expanded into a
    temporary ``.parquet`` file before delegating to pandas/pyarrow.
    """
    if not isinstance(path, (str, os.PathLike, Path)):
        return _PANDAS_READ_PARQUET(path, *args, **kwargs)
    requested = Path(path)
    if not _is_parquet_name(requested):
        return _PANDAS_READ_PARQUET(path, *args, **kwargs)
    physical = resolve_parquet_path(requested)
    if not str(physical).endswith(".parquet.xz"):
        return _PANDAS_READ_PARQUET(physical, *args, **kwargs)
    with tempfile.NamedTemporaryFile("wb", suffix=".parquet", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        with lzma.open(physical, "rb") as source, temporary.open("wb") as destination:
            shutil.copyfileobj(source, destination, length=8 * 1024 * 1024)
        return _PANDAS_READ_PARQUET(temporary, *args, **kwargs)
    finally:
        temporary.unlink(missing_ok=True)


# Compatibility for existing code paths that call pandas directly.  Importing
# src.core makes historical ``pd.read_parquet('foo.parquet')`` calls accept a
# migrated ``foo.parquet.xz`` without changing scientific code fingerprints.
pd.read_parquet = read_parquet_compat


def parquet_payload_sha256(path: str | os.PathLike[str] | Path, block_size: int = 1024 * 1024) -> str:
    """SHA256 of the underlying Parquet bytes, independent of outer XZ wrapping."""
    requested = Path(path)
    physical = resolve_parquet_path(requested)
    digest = hashlib.sha256()
    if str(physical).endswith(".parquet.xz"):
        with lzma.open(physical, "rb") as handle:
            for block in iter(lambda: handle.read(block_size), b""):
                digest.update(block)
    else:
        with physical.open("rb") as handle:
            for block in iter(lambda: handle.read(block_size), b""):
                digest.update(block)
    return digest.hexdigest()


def required_artifacts_present(required: set[str], checksums: Mapping[str, str]) -> bool:
    """Treat ``x.parquet`` and ``x.parquet.xz`` as the same logical artifact."""
    names = set(checksums)
    for name in required:
        if name in names:
            continue
        if name.endswith(".parquet") and f"{name}.xz" in names:
            continue
        return False
    return True


def stored_checksum_matches(directory: Path, name: str, checksum: str) -> bool:
    """Validate marker checksums across a storage-only Parquet→XZ migration.

    Historical DONE markers name/hash the raw Parquet payload.  If only the
    XZ wrapper remains, compare the marker hash against the decompressed
    payload.  New markers that explicitly name ``.parquet.xz`` continue to
    validate the physical compressed file byte-for-byte.
    """
    physical = directory / name
    if physical.is_file():
        return file_sha256(physical) == checksum
    if name.endswith(".parquet"):
        xz = directory / f"{name}.xz"
        if xz.is_file():
            return parquet_payload_sha256(xz) == checksum
    return False


def stable_hash(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def file_sha256(path: Path, block_size: int = 1024 * 1024) -> str:
    path = Path(path)
    if not path.is_file() and str(path).endswith(".parquet") and Path(str(path) + ".xz").is_file():
        # Preserve historical checksums/metadata across storage-only XZ wrapping.
        return parquet_payload_sha256(path, block_size=block_size)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def seed_for(master_seed: int, *parts: object) -> int:
    digest = hashlib.sha256("|".join(map(str, (master_seed, *parts))).encode()).digest()
    return int.from_bytes(digest[:4], "big") & 0x7FFFFFFF


def load_yaml(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected YAML mapping: {path}")
    return value


def load_configuration() -> tuple[dict[str, Path], dict[str, Any]]:
    raw = load_yaml(STAGE / "config/input_paths.yaml")
    configured = Path(raw.pop("project_root"))
    project = (ROOT / configured).resolve() if not configured.is_absolute() else configured.resolve()
    if project != ROOT.resolve():
        raise RuntimeError(f"Configured project root is not this repository: {project}")
    paths: dict[str, Path] = {}
    for key, value in raw.items():
        if key == "schema_version":
            continue
        path = Path(value)
        resolved = (project / path).resolve(strict=True) if not path.is_absolute() else path.resolve(strict=True)
        if not resolved.is_relative_to(project):
            raise RuntimeError(f"Input escapes repository boundary: {key}={resolved}")
        paths[key] = resolved
    return paths, load_yaml(STAGE / "provenance/master_config.yaml")


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False) as handle:
        json.dump(value, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n"); handle.flush(); os.fsync(handle.fileno())
        temporary = Path(handle.name)
    os.replace(temporary, path)


def atomic_table(path: Path, frame: pd.DataFrame) -> None:
    """Atomically write tables; Parquet outputs are outer-XZ wrapped at preset 9."""
    requested = Path(path)
    requested.parent.mkdir(parents=True, exist_ok=True)
    if _is_parquet_name(requested):
        base, final = _parquet_base_and_xz(requested)
        with tempfile.NamedTemporaryFile("wb", dir=requested.parent, delete=False, suffix=".parquet") as handle:
            parquet_tmp = Path(handle.name)
        with tempfile.NamedTemporaryFile("wb", dir=requested.parent, delete=False, suffix=".parquet.xz") as handle:
            xz_tmp = Path(handle.name)
        try:
            frame.to_parquet(parquet_tmp, index=False, compression="zstd")
            with parquet_tmp.open("rb") as source, lzma.open(xz_tmp, "wb", format=lzma.FORMAT_XZ, preset=9) as destination:
                shutil.copyfileobj(source, destination, length=8 * 1024 * 1024)
            os.replace(xz_tmp, final)
            if base != final:
                base.unlink(missing_ok=True)
        finally:
            parquet_tmp.unlink(missing_ok=True)
            xz_tmp.unlink(missing_ok=True)
        return
    with tempfile.NamedTemporaryFile("wb", dir=requested.parent, delete=False, suffix=requested.suffix) as handle:
        temporary = Path(handle.name)
    try:
        frame.to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, requested)
    finally:
        temporary.unlink(missing_ok=True)


def _combine(values: np.ndarray, op: str) -> float:
    finite = values[np.isfinite(values)]
    if not finite.size:
        return np.nan
    if op == "MIN":
        return float(finite.min()) if finite.size == values.size else np.nan
    if op == "MAX": return float(finite.max())
    if op == "SUM": return float(finite.sum())
    if op == "MEAN": return float(finite.mean())
    raise ValueError(op)


def evaluate_gpr(rule: str | ast.AST | None, gene_values: Mapping[str, float], and_operator: str = "MIN", or_operator: str = "MAX") -> float:
    if rule is None: return np.nan
    node = ast.parse(rule, mode="eval").body if isinstance(rule, str) and rule.strip() else rule
    if node is None: return np.nan
    def visit(item: ast.AST) -> float:
        if isinstance(item, ast.Name): return float(gene_values.get(item.id, np.nan))
        if isinstance(item, ast.BoolOp):
            return _combine(np.asarray([visit(x) for x in item.values]), and_operator if isinstance(item.op, ast.And) else or_operator)
        raise TypeError(f"Unsupported GPR syntax: {type(item).__name__}")
    return visit(node)


def reaction_universe_hash(reaction_ids: Sequence[str]) -> str:
    ids = [str(x) for x in reaction_ids]
    if len(ids) != len(set(ids)): raise ValueError("Reaction universe contains duplicate IDs")
    return stable_hash(ids)


def validate_samples(samples: pd.DataFrame, model, config: Mapping[str, float]) -> dict[str, Any]:
    ids = [reaction.id for reaction in model.reactions]
    if list(samples.columns) != ids:
        raise ValueError("QC samples must be in the expert model reaction order")
    x = samples.to_numpy(dtype=np.float64)
    if not np.isfinite(x).all(): raise RuntimeError("Samples contain NaN or Inf")
    import cobra
    S = np.asarray(cobra.util.array.create_stoichiometric_matrix(model, array_type="dense"), dtype=np.float64)
    residual = float(np.max(np.abs(S @ x.T))) if S.size else 0.0
    lb = np.asarray([r.lower_bound for r in model.reactions], dtype=np.float64)
    ub = np.asarray([r.upper_bound for r in model.reactions], dtype=np.float64)
    violation = max(float(np.max(lb - x)), float(np.max(x - ub)), 0.0)
    if residual > float(config["steady_state_tolerance"]): raise RuntimeError(f"Mass-balance residual {residual}")
    if violation > float(config["bound_tolerance"]): raise RuntimeError(f"Bound violation {violation}")
    f32 = x.astype(np.float32).astype(np.float64)
    return {"successful_vectors": len(samples), "largest_mass_balance_residual": residual,
            "largest_bound_violation": violation, "float32_max_abs_error": float(np.max(np.abs(x - f32))),
            "float32_median_abs_error": float(np.median(np.abs(x - f32)))}


def summarize_samples(samples: pd.DataFrame, tolerance: float = 1e-9) -> pd.DataFrame:
    values = samples.astype(np.float64); absolute = values.abs()
    return pd.DataFrame({"reaction_id": values.columns, "R1_median_signed_flux": values.median().to_numpy(),
        "R2_mean_signed_flux": values.mean().to_numpy(), "R3_median_absolute_flux": absolute.median().to_numpy(),
        "R4_mean_absolute_flux": absolute.mean().to_numpy(), "R5_absolute_median_flux": values.median().abs().to_numpy(),
        "R6_probability_active": (absolute > tolerance).mean().to_numpy(),
        "R7_directionality_probability": ((values > tolerance).mean() - (values < -tolerance).mean()).to_numpy(),
        "sampling_sd": values.std(ddof=1).to_numpy(), "retained_vectors": len(values)})
