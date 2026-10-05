"""Deterministic storage helpers for the DMI-BRIDGE PL1 artifacts."""
from __future__ import annotations

import csv
import hashlib
import json
import lzma
import os
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np

XZ_NPY_PRESET = 9
XZ_TSV_PRESET = 0
FEATURE_PARTS = 4
MAX_PART_BYTES = 90 * 1024 * 1024
FEATURE_NAME = "BRIDGEPL1_REACTION_EVALUATION_FEATURES"
FEATURE_MANIFEST = FEATURE_NAME + ".parts.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_replace(path: Path, writer) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.unlink(missing_ok=True)
    try:
        writer(tmp)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def compress_file_xz(source: Path, destination: Path, *, preset: int = XZ_NPY_PRESET) -> str:
    def write(tmp: Path) -> None:
        with source.open("rb") as src, lzma.open(tmp, "wb", preset=preset, format=lzma.FORMAT_XZ) as dst:
            shutil.copyfileobj(src, dst, length=1 << 20)

    _atomic_replace(destination, write)
    return sha256_file(destination)


def decompress_file_xz(source: Path, destination: Path) -> None:
    def write(tmp: Path) -> None:
        with lzma.open(source, "rb", format=lzma.FORMAT_XZ) as src, tmp.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=1 << 20)

    _atomic_replace(destination, write)


def compress_file_xz_verified(source: Path, destination: Path, *, preset: int = XZ_NPY_PRESET) -> str:
    """Atomically compress a raw artifact and verify its exact byte stream."""
    source_sha = sha256_file(source)
    compressed_sha = compress_file_xz(source, destination, preset=preset)
    with tempfile.TemporaryDirectory(prefix="xz-verify-") as td:
        raw = Path(td) / source.name
        decompress_file_xz(destination, raw)
        if sha256_file(raw) != source_sha:
            destination.unlink(missing_ok=True)
            raise RuntimeError(f"XZ round-trip mismatch: {source}")
    return compressed_sha


def save_npy_xz(path: Path, array: np.ndarray) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="pl1-npy-") as td:
        raw = Path(td) / "array.npy"
        np.save(raw, array, allow_pickle=False)
        return compress_file_xz(raw, path, preset=XZ_NPY_PRESET)


def load_npy_xz(path: Path, mmap_mode: str | None = None) -> np.ndarray:
    if mmap_mode is not None:
        raise ValueError("memory mapping compressed NPY is unsupported; use open_npy_xz_work")
    with tempfile.TemporaryDirectory(prefix="pl1-npy-load-") as td:
        raw = Path(td) / "array.npy"
        decompress_file_xz(path, raw)
        return np.load(raw, allow_pickle=False).copy()


class NpyXZWork:
    """External temporary mmap for a canonical compressed NPY checkpoint."""

    def __init__(self, compressed: Path, legacy: Path | None = None, *, shape=None, dtype=None):
        self.compressed = compressed
        self.legacy = legacy
        self._tmpdir = Path(tempfile.mkdtemp(prefix="pl1-npy-work-"))
        self.raw = self._tmpdir / "array.npy"
        if compressed.is_file() and legacy is not None and legacy.is_file():
            with tempfile.TemporaryDirectory(prefix="pl1-npy-compare-") as td:
                candidate = Path(td) / "compressed.npy"
                decompress_file_xz(compressed, candidate)
                if sha256_file(candidate) != sha256_file(legacy):
                    self.close(remove_temp=True)
                    raise RuntimeError(f"conflicting compressed and legacy NPY files: {compressed}")
        if compressed.is_file():
            decompress_file_xz(compressed, self.raw)
        elif legacy is not None and legacy.is_file():
            shutil.copyfile(legacy, self.raw)
        elif shape is not None and dtype is not None:
            np.lib.format.open_memmap(self.raw, mode="w+", dtype=dtype, shape=shape).flush()
        else:
            self.close(remove_temp=True)
            raise FileNotFoundError(compressed)
        self.array = np.load(self.raw, mmap_mode="r+", allow_pickle=False)

    def checkpoint(self) -> str:
        self.array.flush()
        return compress_file_xz(self.raw, self.compressed, preset=XZ_NPY_PRESET)

    def close(self, *, remove_temp: bool = True) -> None:
        if getattr(self, "array", None) is not None:
            self.array.flush()
            self.array = None
        if remove_temp:
            shutil.rmtree(self._tmpdir, ignore_errors=True)


def _part_name(index: int) -> str:
    return f"{FEATURE_NAME}.part-{index:03d}.tsv.xz"


def _manifest_path(directory: Path) -> Path:
    return directory / FEATURE_MANIFEST


def _write_json(path: Path, value: dict) -> None:
    _atomic_replace(path, lambda tmp: tmp.write_bytes(json.dumps(value, sort_keys=True, indent=2).encode("utf-8") + b"\n"))


def _row_text(fields: list[str], row: dict) -> bytes:
    from io import StringIO
    buf = StringIO()
    writer = csv.DictWriter(buf, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
    writer.writerow({key: str(row[key]) for key in fields})
    return buf.getvalue().encode("utf-8")


def write_partitioned_tsv(
    directory: Path,
    fields: list[str],
    rows: Iterable[dict],
    *,
    total_rows: int,
    source_monolithic_sha256: str | None = None,
    schema: str = "bridge.pl1.partitioned_tsv.v1",
) -> dict:
    if total_rows < FEATURE_PARTS:
        raise ValueError("total_rows must be at least the number of parts")
    directory.mkdir(parents=True, exist_ok=True)
    boundaries = [(total_rows * i) // FEATURE_PARTS for i in range(FEATURE_PARTS + 1)]
    paths = [directory / _part_name(i) for i in range(FEATURE_PARTS)]
    for path in paths:
        path.unlink(missing_ok=True)
    handles = [lzma.open(path, "wt", encoding="utf-8", newline="", preset=XZ_TSV_PRESET, format=lzma.FORMAT_XZ) for path in paths]
    writers = [csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise") for handle in handles]
    try:
        for writer in writers:
            writer.writeheader()
        counts = [0] * FEATURE_PARTS
        seen = 0
        for row in rows:
            if seen >= total_rows:
                raise ValueError("row iterator produced more rows than declared")
            part = min(FEATURE_PARTS - 1, (seen * FEATURE_PARTS) // total_rows)
            writers[part].writerow({key: str(row[key]) for key in fields})
            counts[part] += 1
            seen += 1
        if seen != total_rows:
            raise ValueError(f"row iterator produced {seen} rows; expected {total_rows}")
    finally:
        for handle in handles:
            handle.close()
    manifest = validate_partitioned_tsv(directory, fields=fields)
    manifest.update({"schema": schema, "logical_artifact": FEATURE_NAME, "format": "tsv.xz.parts", "source_monolithic_sha256": source_monolithic_sha256})
    _write_json(_manifest_path(directory), manifest)
    return manifest


def _canonical_parts(manifest: dict) -> list[str]:
    parts = manifest.get("parts")
    if not isinstance(parts, list) or len(parts) != FEATURE_PARTS:
        raise RuntimeError("invalid PL1 feature parts manifest")
    expected = [_part_name(i) for i in range(FEATURE_PARTS)]
    names = [p.get("filename") for p in parts]
    if names != expected:
        raise RuntimeError("PL1 feature parts are not in canonical numeric order")
    return names


def validate_partitioned_tsv(directory: Path, *, fields: list[str] | None = None) -> dict:
    manifest_path = _manifest_path(directory)
    old_manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else None
    expected_names = _canonical_parts(old_manifest) if old_manifest else [_part_name(i) for i in range(FEATURE_PARTS)]
    actual_part_names = sorted(path.name for path in directory.glob(f"{FEATURE_NAME}.part-*.tsv.xz"))
    if actual_part_names != sorted(expected_names):
        raise RuntimeError("PL1 feature parts include missing or unmanifested files")
    header: list[str] | None = None
    part_records = []
    logical = hashlib.sha256()
    total = 0
    for index, name in enumerate(expected_names):
        path = directory / name
        if not path.is_file():
            raise RuntimeError(f"missing PL1 feature part: {path}")
        if path.stat().st_size >= MAX_PART_BYTES:
            raise RuntimeError(f"PL1 feature part exceeds safety ceiling: {path}")
        count = 0
        with lzma.open(path, "rt", encoding="utf-8", newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            try:
                current = next(reader)
            except StopIteration as exc:
                raise RuntimeError(f"empty PL1 feature part: {path}") from exc
            if header is None:
                header = current
                logical.update(("\t".join(header) + "\n").encode("utf-8"))
            elif current != header:
                raise RuntimeError(f"PL1 feature schema mismatch in {path}")
            if fields is not None and current != fields:
                raise RuntimeError("PL1 feature schema does not match expected fields")
            for row in reader:
                if len(row) != len(header):
                    raise RuntimeError(f"PL1 feature row width mismatch in {path}")
                logical.update(("\t".join(row) + "\n").encode("utf-8"))
                count += 1
        part_records.append({"filename": name, "rows": count, "bytes": path.stat().st_size, "sha256": sha256_file(path)})
        total += count
    result = {"parts": part_records, "number_of_parts": FEATURE_PARTS, "total_data_rows": total, "column_names": header, "logical_content_sha256": logical.hexdigest()}
    if old_manifest:
        for key in ("number_of_parts", "total_data_rows", "column_names", "logical_content_sha256"):
            if old_manifest.get(key) != result[key]:
                raise RuntimeError(f"PL1 feature manifest mismatch: {key}")
        if old_manifest.get("parts") != result["parts"]:
            raise RuntimeError("PL1 feature part hash/size/row manifest mismatch")
    return result


def discover_feature_artifact(directory: Path) -> tuple[str, object]:
    manifest = _manifest_path(directory)
    legacy = directory / (FEATURE_NAME + ".tsv.xz")
    if manifest.is_file():
        checked = validate_partitioned_tsv(directory)
        if legacy.is_file():
            raise RuntimeError("legacy and partitioned PL1 feature artifacts coexist")
        return "parts", checked
    if legacy.is_file():
        return "legacy", legacy
    raise FileNotFoundError(f"no PL1 feature artifact in {directory}")


def iter_feature_rows(directory: Path) -> Iterator[dict[str, str]]:
    kind, value = discover_feature_artifact(directory)
    if kind == "legacy":
        with lzma.open(value, "rt", encoding="utf-8", newline="") as fh:
            yield from csv.DictReader(fh, delimiter="\t")
        return
    manifest = json.loads(_manifest_path(directory).read_text())
    for name in _canonical_parts(manifest):
        with lzma.open(directory / name, "rt", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            yield from reader


def logical_legacy_tsv_digest(path: Path) -> tuple[list[str], int, str]:
    logical = hashlib.sha256()
    rows = 0
    with lzma.open(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        logical.update(("\t".join(header) + "\n").encode("utf-8"))
        for row in reader:
            logical.update(("\t".join(row) + "\n").encode("utf-8"))
            rows += 1
    return header, rows, logical.hexdigest()
