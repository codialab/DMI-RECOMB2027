#!/usr/bin/env python3
"""Freeze paired A1/A2-L Figure 4 plotting data from validated upstream artifacts."""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import lzma
import math
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "reproduced/derived/fig4"
DATA = OUT / "data"
GEOM_OUT = DATA / "fig4_geometry_paired"
UTILITY_OUT = DATA / "fig4_utility_paired"
SUPPORT_OUT = DATA / "fig4_panelC_candidate_support.parquet.xz"
PANEL_METRICS = DATA / "fig4_panelC_candidate_metrics.parquet.xz"
PANEL_GROUPS = DATA / "fig4_panelC_candidate_groups.parquet.xz"
METHODS = ("CORDA", "GIMME", "iMAT", "RIPTiDe")
CAND_KEYS = ["rna_context_key","ct2a_mouse","gl261_mouse","reaction_id","truth_selection"]
LAMBDA = 0.25
TRUTH_ATOL = 1e-12
HDIR_ATOL = 5e-12
TIE_TOL = 1e-12
TOP_N = 100

PL1 = ROOT / "outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1"
A21 = ROOT / "outputs/dmi_bridge_a21_dual_anchor_geometry_v1"
PL2B = ROOT / "outputs/dmi_bridge_pl2b_sign_only_utility_v1"
A22 = ROOT / "outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1"
A23 = ROOT / "outputs/dmi_bridge_a23_a1_a2_synthesis_v1"
A20 = ROOT / "outputs/dmi_bridge_a20_dual_anchor_qualification_v1"
BIO0 = ROOT / "outputs/dmi_bridge_bio0_audit_v1"
PL2A = ROOT / "outputs/dmi_bridge_pl2a_sign_only_prepare_v1"
CANDIDATE_FILE = ROOT / "12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/candidate_fraction_table.tsv.xz"
CACHE_FILE = ROOT / "outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstructed/flux_cache.npz"

GEOM_FIELDS = {
    "H_dir": "directional_entropy3",
    "dominant_direction_mass": "dominant_sign_mass",
    "direction_explained_magnitude_variance": "sign_magnitude_eta2",
    "supported_sign_state_count": "n_supported_sign_states",
    "non_tie_coverage": "non_tie_coverage",
    "ESS": "joint_ess",
}
# Frozen reaction-specific summaries of P(Delta v_B). Keep these separate from
# GEOM_FIELDS so adding Panel B width does not redefine the legacy core-finite QC.
DISTRIBUTION_SUMMARY_FIELDS = (
    "delta_v_B_mean",
    "delta_v_B_sd",
    "delta_v_B_q10",
    "delta_v_B_q50",
    "delta_v_B_q90",
    "delta_v_B_width80",
)
GEOM_SOURCE_FIELDS = (
    "evaluation_id", "algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse",
    "reaction_id", "subsystem", "cartesian_support", "joint_ess", "p_tie",
    "dominant_sign_mass", "directional_entropy3", "n_supported_sign_states",
    "sign_magnitude_eta2", "delta_b_mean", "delta_b_sd", "delta_b_q10",
    "delta_b_q50", "delta_b_q90", "delta_b_width80", "target_delta_degenerate",
    "sign_magnitude_status",
)
CASE_FIELDS = (
    "truth_pair_id", "evaluation_id", "algorithm", "rna_context_key", "ct2a_mouse",
    "gl261_mouse", "reaction_id", "subsystem", "truth_delta_b", "truth_magnitude",
    "truth_direction", "truth_status", "baseline_magnitude_estimate", "baseline_abs_error",
    "lambda_0_25_correct_magnitude_estimate", "lambda_0_25_wrong_magnitude_estimate",
    "lambda_0_25_correct_absolute_error", "lambda_0_25_wrong_absolute_error",
    "lambda_0_25_correct_absolute_error_gain", "lambda_0_25_wrong_absolute_error_gain",
    "lambda_0_25_correct_gain_status", "lambda_0_25_wrong_gain_status",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def safe_source(path: Path) -> Path:
    require(path.is_file() and not path.is_symlink(), f"missing, symlinked, or non-file source: {path}")
    resolved = path.resolve()
    require(resolved.is_relative_to(ROOT), f"source escapes repository root: {path}")
    return resolved


def read_json(path: Path) -> dict:
    return json.loads(safe_source(path).read_text())


def verify_artifacts(directory: Path, manifest_name: str) -> dict:
    manifest_path = safe_source(directory / manifest_name)
    manifest = json.loads(manifest_path.read_text())
    for name, expected in manifest.get("artifact_sha256", {}).items():
        path = safe_source(directory / name)
        actual = sha256(path)
        require(actual == expected, f"source artifact digest mismatch: {path}")
    return manifest


def tabular_parts(directory: Path, prefix: str) -> list[Path]:
    parts_path = safe_source(directory / f"{prefix}.parts.json")
    info = json.loads(parts_path.read_text())
    paths = [safe_source(directory / part["filename"]) for part in info["parts"]]
    for path, record in zip(paths, info["parts"]):
        require(sha256(path) == record["sha256"], f"partition digest mismatch: {path}")
    return paths


def csv_rows(path: Path):
    stream = lzma.open(safe_source(path), "rt", newline="") if path.suffix == ".xz" else safe_source(path).open("rt", newline="")
    try:
        reader = csv.DictReader(stream, delimiter="\t")
        require(reader.fieldnames is not None, f"missing header in {path}")
        yield from reader
    finally:
        stream.close()


def read_selected_parts(directory: Path, prefix: str, selected: list[str], drop_ldh: bool = False) -> pd.DataFrame:
    chunks = []
    for part in tabular_parts(directory, prefix):
        frame = pd.read_csv(part, sep="\t", compression="xz", usecols=selected)
        if drop_ldh:
            frame = frame.loc[frame.reaction_id.ne("LDH_L")]
        chunks.append(frame)
    result = pd.concat(chunks, ignore_index=True)
    return result


def field_frame(directory: Path, prefix: str, arm: str) -> pd.DataFrame:
    frame = read_selected_parts(directory, prefix, list(GEOM_SOURCE_FIELDS), drop_ldh=True)
    require(frame.reaction_id.nunique() == 4179, f"{arm}: expected 4,179 target reactions")
    frame["non_tie_coverage"] = 1.0 - pd.to_numeric(frame.p_tie, errors="coerce")
    rename = {raw: f"{arm}_{canonical}" for canonical, raw in GEOM_FIELDS.items() if raw != "non_tie_coverage"}
    rename.update({
        "cartesian_support": f"{arm}_cartesian_support",
        "delta_b_mean": f"{arm}_delta_v_B_mean",
        "delta_b_sd": f"{arm}_delta_v_B_sd",
        "delta_b_q10": f"{arm}_delta_v_B_q10",
        "delta_b_q50": f"{arm}_delta_v_B_q50",
        "delta_b_q90": f"{arm}_delta_v_B_q90",
        "delta_b_width80": f"{arm}_delta_v_B_width80",
        "target_delta_degenerate": f"{arm}_target_delta_degenerate",
        "sign_magnitude_status": f"{arm}_sign_magnitude_status",
        "non_tie_coverage": f"{arm}_non_tie_coverage",
    })
    frame = frame.rename(columns=rename)
    keep=["evaluation_id","algorithm","rna_context_key","ct2a_mouse","gl261_mouse","reaction_id","subsystem",f"{arm}_cartesian_support",f"{arm}_target_delta_degenerate",f"{arm}_sign_magnitude_status"]
    keep.extend(f"{arm}_{metric}" for metric in GEOM_FIELDS if f"{arm}_{metric}" in frame.columns)
    keep.extend(f"{arm}_{metric}" for metric in DISTRIBUTION_SUMMARY_FIELDS if f"{arm}_{metric}" in frame.columns)
    return frame[keep]


def require_parquet_xz(path: Path) -> Path:
    require(str(path).endswith(".parquet.xz"), f"expected .parquet.xz path: {path}")
    return path


def write_parquet_xz(frame: pd.DataFrame, path: Path) -> None:
    """Write ZSTD Parquet wrapped in XZ, without a repository-side Parquet temporary."""
    require_parquet_xz(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="fig4-parquet-") as temporary:
        parquet_path = Path(temporary) / "table.parquet"
        with tempfile.NamedTemporaryFile(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent, delete=False) as handle:
            xz_path = Path(handle.name)
        frame.to_parquet(parquet_path, engine="pyarrow", compression="zstd", index=False)
        try:
            with parquet_path.open("rb") as source, lzma.open(xz_path, "wb", preset=9, format=lzma.FORMAT_XZ) as compressed:
                shutil.copyfileobj(source, compressed, length=1024 * 1024)
            os.replace(xz_path, path)
        finally:
            xz_path.unlink(missing_ok=True)


def compress_parquet_xz(source: Path, destination: Path) -> None:
    """Atomically compress an ordinary Parquet file to a canonical XZ artifact."""
    require(source.is_file(), f"missing temporary Parquet source: {source}")
    require_parquet_xz(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        with source.open("rb") as input_stream, lzma.open(temporary, "wb", preset=9, format=lzma.FORMAT_XZ) as output_stream:
            shutil.copyfileobj(input_stream, output_stream, length=1024 * 1024)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def _decompress_parquet_xz(path: Path, temporary: Path) -> None:
    require(path.is_file(), f"missing .parquet.xz source: {path}")
    require_parquet_xz(path)
    with lzma.open(path, "rb", format=lzma.FORMAT_XZ) as source, temporary.open("wb") as destination:
        shutil.copyfileobj(source, destination, length=1024 * 1024)


def read_parquet_xz(path: Path, columns: list[str] | None = None) -> pd.DataFrame:
    require(path.is_file(), f"missing .parquet.xz source: {path}")
    require_parquet_xz(path)
    with tempfile.TemporaryDirectory(prefix="fig4-parquet-read-") as temporary:
        parquet_path = Path(temporary) / "table.parquet"
        _decompress_parquet_xz(path, parquet_path)
        return pd.read_parquet(parquet_path, columns=columns)


def read_parquet_xz_metadata(path: Path) -> pq.FileMetaData:
    require(path.is_file(), f"missing .parquet.xz source: {path}")
    require_parquet_xz(path)
    with tempfile.TemporaryDirectory(prefix="fig4-parquet-metadata-") as temporary:
        parquet_path = Path(temporary) / "table.parquet"
        _decompress_parquet_xz(path, parquet_path)
        return pq.read_metadata(parquet_path)


def read_pyarrow_table_xz(path: Path, columns: list[str] | None = None) -> pa.Table:
    require(path.is_file(), f"missing .parquet.xz source: {path}")
    require_parquet_xz(path)
    with tempfile.TemporaryDirectory(prefix="fig4-parquet-arrow-") as temporary:
        parquet_path = Path(temporary) / "table.parquet"
        _decompress_parquet_xz(path, parquet_path)
        return pq.read_table(parquet_path, columns=columns)


def sha256_decompressed_xz(path: Path) -> str:
    require_parquet_xz(path)
    digest = hashlib.sha256()
    with lzma.open(path, "rb", format=lzma.FORMAT_XZ) as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_parquet_xz(path: Path) -> str:
    metadata = read_parquet_xz_metadata(path)
    require(metadata.num_rows >= 0, f"invalid Parquet metadata: {path}")
    return sha256_decompressed_xz(path)


def migrate_legacy_storage() -> int:
    """Migrate Figure 4 bare Parquet files before any reuse or path validation."""
    DATA.mkdir(parents=True, exist_ok=True)
    bare_files = sorted(DATA.rglob("*.parquet"))
    removed = 0
    for bare in bare_files:
        target = bare.with_suffix(bare.suffix + ".xz")
        source_digest = sha256(bare)
        if target.is_file():
            target_digest = validate_parquet_xz(target)
            require(target_digest == source_digest, f"bare/XZ Parquet payload mismatch: {bare} vs {target}")
        else:
            with tempfile.TemporaryDirectory(prefix="fig4-parquet-migrate-") as temporary:
                temporary_source = Path(temporary) / "source.parquet"
                shutil.copyfile(bare, temporary_source)
                compress_parquet_xz(temporary_source, target)
            target_digest = validate_parquet_xz(target)
            require(target_digest == source_digest, f"new XZ Parquet payload validation failed: {target}")
        bare.unlink()
        removed += 1
    upgrade_manifest_storage()
    print(f"Figure 4 Parquet storage migration: removed {removed} bare .parquet file(s).", flush=True)
    return removed


def canonical_table_path(value: str) -> str:
    return value[:-len(".parquet")] + ".parquet.xz" if value.endswith(".parquet") else value


def upgrade_manifest_storage() -> None:
    path = DATA / "fig4_data_manifest.json"
    if not path.is_file():
        return
    manifest = json.loads(path.read_text())
    generated = [canonical_table_path(value) for value in manifest.get("generated_files", [])]
    table_paths = sorted(str(candidate.relative_to(ROOT)) for candidate in DATA.rglob("*.parquet.xz"))
    generated = sorted(set(generated) | set(table_paths))
    manifest["generated_files"] = generated
    hashes = {}
    for value, digest in manifest.get("generated_sha256", {}).items():
        canonical = canonical_table_path(value)
        if canonical not in hashes and not value.endswith(".parquet.xz"):
            hashes[canonical] = digest
    for value in generated:
        candidate = ROOT / value
        if candidate.is_file():
            hashes[value] = sha256(candidate)
    manifest["generated_sha256"] = hashes
    manifest["table_storage"] = {"format": "xz-wrapped parquet", "parquet_internal_compression": "zstd", "outer_compression": "xz", "xz_preset": 9, "persistent_uncompressed_parquet": False}
    manifest.pop("xz_archives", None)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n")


def write_partitioned_frame(frame: pd.DataFrame, directory: Path, prefix: str, manifest: dict) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for method, group in frame.groupby("algorithm", sort=True, dropna=False):
        method_dir = directory / f"algorithm={method}"
        method_dir.mkdir(parents=True, exist_ok=True)
        path = method_dir / f"{prefix}.parquet.xz"
        write_parquet_xz(group.drop(columns=["algorithm"]), path)
        manifest.setdefault("generated_files", []).append(str(path.relative_to(ROOT)))


def build_geometry(manifest: dict) -> pd.DataFrame:
    a1 = field_frame(PL1, "BRIDGEPL1_REACTION_EVALUATION_FEATURES", "A1")
    a2 = field_frame(A21, "BRIDGEA21_REACTION_EVALUATION_FEATURES", "A2")
    keys = ["evaluation_id", "reaction_id"]
    require(not a1.duplicated(keys).any(), "duplicate A1 evaluation/reaction geometry key")
    require(not a2.duplicated(keys).any(), "duplicate A2 evaluation/reaction geometry key")
    require(len(a1) == len(a2) == 400 * 4179, "geometry source rectangles have unexpected sizes")
    merged = a1.merge(a2, on=keys, how="outer", suffixes=("_A1", "_A2"), validate="one_to_one", indicator=True)
    require(merged._merge.eq("both").all(), "A1/A2 geometry key populations differ")
    for field in ("algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse", "subsystem"):
        require(merged[f"{field}_A1"].eq(merged[f"{field}_A2"]).all(), f"geometry metadata mismatch: {field}")
    merged = merged.rename(columns={f"{field}_A1": field for field in ("algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse", "subsystem")})
    merged = merged.drop(columns=[f"{field}_A2" for field in ("algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse", "subsystem")])
    merged["reaction_name"] = pd.NA  # No canonical reaction-name field exists in the frozen source schemas.
    merged["biological_comparison"] = merged["ct2a_mouse"].astype(str) + "_vs_" + merged["gl261_mouse"].astype(str)
    for metric in GEOM_FIELDS:
        x, y = f"A1_{metric}", f"A2_{metric}"
        merged[x] = pd.to_numeric(merged[x], errors="coerce")
        merged[y] = pd.to_numeric(merged[y], errors="coerce")
        merged[f"delta_{metric}"] = merged[y] - merged[x]
    for metric in DISTRIBUTION_SUMMARY_FIELDS:
        x, y = f"A1_{metric}", f"A2_{metric}"
        merged[x] = pd.to_numeric(merged[x], errors="coerce")
        merged[y] = pd.to_numeric(merged[y], errors="coerce")
        merged[f"delta_{metric}"] = merged[y] - merged[x]
    merged["paired_width80_finite"] = np.isfinite(merged.A1_delta_v_B_width80) & np.isfinite(merged.A2_delta_v_B_width80)
    positive = np.isfinite(merged.A1_ESS) & (merged.A1_ESS > 0) & np.isfinite(merged.A2_ESS) & (merged.A2_ESS > 0)
    merged["log10_ESS_ratio"] = np.nan
    merged.loc[positive, "log10_ESS_ratio"] = np.log10(merged.loc[positive, "A2_ESS"] / merged.loc[positive, "A1_ESS"])
    for arm in ("A1", "A2"):
        finite = np.ones(len(merged), dtype=bool)
        for metric in GEOM_FIELDS:
            finite &= np.isfinite(merged[f"{arm}_{metric}"])
        merged[f"{arm}_geometry_core_finite"] = finite
    merged["paired_geometry_core_finite"] = merged.A1_geometry_core_finite & merged.A2_geometry_core_finite
    merged["ess_ratio_valid"] = positive
    merged = merged.drop(columns=["_merge"], errors="ignore")
    write_partitioned_frame(merged, GEOM_OUT, "part", manifest)
    manifest["row_counts"]["geometry_paired"] = int(len(merged))
    manifest["unique_key_counts"]["geometry_evaluation_reaction"] = int(merged[keys].drop_duplicates().shape[0])
    manifest["duplicate_key_checks"]["geometry_evaluation_reaction"] = "PASS"
    manifest["a1_a2_matching"]["geometry_rows_matched"] = int(merged.paired_geometry_core_finite.size)
    manifest["finite_value_counts"]["geometry_paired_core_finite"] = int(merged.paired_geometry_core_finite.sum())
    manifest["finite_value_counts"]["geometry_paired_width80_finite"] = int(merged.paired_width80_finite.sum())
    return merged


def selected_case_parts(directory: Path, prefix: str, pair_ids: set[str], source_arm: str):
    paths = sorted(directory.glob(f"{prefix}.part-*.tsv.xz"))
    def read_part(path):
        with lzma.open(safe_source(path), "rt", newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            for row in reader:
                if row["truth_pair_id"] not in pair_ids or row["reaction_id"] == "LDH_L":
                    continue
                yield (row["truth_pair_id"], row["reaction_id"]), row
    merged = itertools.chain.from_iterable(read_part(path) for path in paths)
    previous = None
    for key, row in merged:
        require(previous is None or key > previous, f"{source_arm} utility keys duplicated or unordered at {key}")
        previous = key
        yield key, row


def load_pair_registry() -> tuple[dict, set[str]]:
    registry = pd.read_csv(A22 / "BRIDGEA22_MATCHED_TRUTH_PAIR_REGISTRY.tsv", sep="\t", dtype=str)
    registry = registry.loc[registry.pair_evaluable.eq("true")].copy()
    require(len(registry) == registry.truth_pair_id.nunique() == 836, "matched truth-pair registry must contain 836 unique pairs")
    return {r.truth_pair_id: r._asdict() for r in registry.itertuples(index=False)}, set(registry.truth_pair_id)


def build_utility(manifest: dict, geom: pd.DataFrame, pair_meta: dict, pair_ids: set[str]) -> tuple[int, dict]:
    a1_iter = selected_case_parts(PL2B, "BRIDGEPL2B_CASE_OUTCOMES", pair_ids, "A1")
    a2_iter = selected_case_parts(A22, "BRIDGEA22_CASE_OUTCOMES", pair_ids, "A2")
    method_buffers = {m: [] for m in METHODS}
    writers = {}
    utility_agg = defaultdict(list)
    ranked_groups = rank_candidate_groups(geom,pair_meta)
    selected_top = ranked_groups.head(TOP_N)
    group_ids={tuple(row):gid for row,gid in zip(selected_top[CAND_KEYS].itertuples(index=False,name=None),selected_top.candidate_group_id)}
    selected_group_ids = set(group_ids.values())
    count = 0
    pairs_seen = Counter()
    schema = None
    def flush(method: str):
        nonlocal schema
        rows = method_buffers[method]
        if not rows:
            return
        frame = pd.DataFrame(rows)
        frame["algorithm"] = method
        frame = frame[["truth_pair_id", "evaluation_id", "algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse", "truth_selection", "reaction_id", "subsystem", "m_star", "truth_delta_b", "A1_m0", "A1_m_correct", "A1_m_wrong", "A2_m0", "A2_m_correct", "A2_m_wrong", "A1_error0", "A1_error_correct", "A1_error_wrong", "A1_gain_correct", "A1_gain_wrong", "A1_correct_vs_wrong_contrast", "A2_error0", "A2_error_correct", "A2_error_wrong", "A2_gain_correct", "A2_gain_wrong", "A2_correct_vs_wrong_contrast", "delta_correct_vs_wrong_contrast", "truth_status", "A1_utility_values_valid", "A2_utility_values_valid"]]
        method_dir = UTILITY_OUT / f"algorithm={method}"
        method_dir.mkdir(parents=True, exist_ok=True)
        table = pa.Table.from_pandas(frame, preserve_index=False)
        if method not in writers:
            temporary = tempfile.NamedTemporaryFile(prefix=f"fig4-{method}-", suffix=".parquet", delete=False)
            temporary.close()
            path = Path(temporary.name)
            writers[method] = pq.ParquetWriter(path, table.schema, compression="zstd")
            writers[method + "_path"] = path
        writers[method].write_table(table)
        method_buffers[method] = []
    for (key1, a1), (key2, a2) in zip(a1_iter, a2_iter):
        require(key1 == key2, f"A1/A2 utility key mismatch: {key1} vs {key2}")
        pair_id, reaction_id = key1
        meta = pair_meta[pair_id]
        require(a1["evaluation_id"] == a2["evaluation_id"] == meta["evaluation_id"], f"utility evaluation mismatch for {pair_id}")
        require(a1["algorithm"] == a2["algorithm"] == meta["algorithm"], f"utility method mismatch for {pair_id}")
        truth1, truth2 = float(a1["truth_delta_b"]), float(a2["truth_delta_b"])
        mag1, mag2 = float(a1["truth_magnitude"]), float(a2["truth_magnitude"])
        require(math.isclose(truth1, truth2, rel_tol=0.0, abs_tol=TRUTH_ATOL) and math.isclose(mag1, mag2, rel_tol=0.0, abs_tol=TRUTH_ATOL), f"truth mismatch for pair={pair_id}, reaction={reaction_id}: A1=({truth1},{mag1}) A2=({truth2},{mag2})")
        mstar = abs(truth1)
        values = {}
        for arm, row in (("A1", a1), ("A2", a2)):
            m0 = float(row["baseline_magnitude_estimate"])
            mc = float(row["lambda_0_25_correct_magnitude_estimate"])
            mw = float(row["lambda_0_25_wrong_magnitude_estimate"])
            e0 = float(row["baseline_abs_error"])
            ec = float(row["lambda_0_25_correct_absolute_error"])
            ew = float(row["lambda_0_25_wrong_absolute_error"])
            gc = float(row["lambda_0_25_correct_absolute_error_gain"])
            gw = float(row["lambda_0_25_wrong_absolute_error_gain"])
            for label, actual, canonical in (("error0", e0, abs(m0-mstar)), ("error_correct", ec, abs(mc-mstar)), ("error_wrong", ew, abs(mw-mstar)), ("gain_correct", gc, e0-ec), ("gain_wrong", gw, e0-ew)):
                both_missing = math.isnan(actual) and math.isnan(canonical)
                require(both_missing or math.isclose(actual, canonical, rel_tol=1e-10, abs_tol=1e-12), f"{arm} frozen {label} definition mismatch at {pair_id}/{reaction_id}: stored={actual}, recomputed={canonical}")
            values.update({f"{arm}_m0":m0, f"{arm}_m_correct":mc, f"{arm}_m_wrong":mw, f"{arm}_error0":e0, f"{arm}_error_correct":ec, f"{arm}_error_wrong":ew, f"{arm}_gain_correct":gc, f"{arm}_gain_wrong":gw, f"{arm}_correct_vs_wrong_contrast":gc-gw})
        meta_row = pair_meta[pair_id]
        truth_selection = meta_row.get("quantile_aliases", "")
        row = {"truth_pair_id":pair_id,"evaluation_id":a1["evaluation_id"],"rna_context_key":a1["rna_context_key"],"ct2a_mouse":a1["ct2a_mouse"],"gl261_mouse":a1["gl261_mouse"],"truth_selection":truth_selection,"reaction_id":reaction_id,"subsystem":a1["subsystem"],"m_star":mstar,"truth_delta_b":truth1,"truth_status":a1["truth_status"],**values,"delta_correct_vs_wrong_contrast":values["A2_correct_vs_wrong_contrast"]-values["A1_correct_vs_wrong_contrast"],"A1_utility_values_valid":all(math.isfinite(values[f"A1_{name}"]) for name in ("m0","m_correct","m_wrong","error0","error_correct","error_wrong","gain_correct","gain_wrong")),"A2_utility_values_valid":all(math.isfinite(values[f"A2_{name}"]) for name in ("m0","m_correct","m_wrong","error0","error_correct","error_wrong","gain_correct","gain_wrong"))}
        method = a1["algorithm"]
        require(method in METHODS, f"unexpected utility method {method}")
        method_buffers[method].append(row)
        pairs_seen[pair_id] += 1
        count += 1
        gid = group_ids.get((row["rna_context_key"],row["ct2a_mouse"],row["gl261_mouse"],reaction_id,truth_selection))
        if gid in selected_group_ids:
            utility_agg[(gid, method)].append({k:v for k,v in row.items() if k.endswith("gain_correct") or k.endswith("gain_wrong") or k.endswith("contrast") or k == "truth_status"})
        if len(method_buffers[method]) >= 10_000:
            flush(method)
        if count % 500_000 == 0:
            print(f"  paired utility rows checked: {count:,}",flush=True)
    for method in METHODS:
        flush(method)
    for writer in writers.values():
        if isinstance(writer, pq.ParquetWriter):
            writer.close()
    for method in METHODS:
        path = writers.get(method + "_path")
        if path is not None:
            method_dir = UTILITY_OUT / f"algorithm={method}"
            final_path = method_dir / "part-000.parquet.xz"
            try:
                compress_parquet_xz(path, final_path)
            finally:
                path.unlink(missing_ok=True)
            manifest["generated_files"].append(str(final_path.relative_to(ROOT)))
    require(count == 836 * 4179, f"paired utility row count mismatch: {count}")
    require(len(pairs_seen) == 836 and set(pairs_seen.values()) == {4179}, "paired utility truth-pair rectangle incomplete")
    manifest["row_counts"]["utility_paired"] = count
    manifest["unique_key_counts"]["utility_truth_pair_reaction"] = count
    manifest["duplicate_key_checks"]["utility_truth_pair_reaction"] = "PASS"
    manifest["a1_a2_matching"]["utility_rows_matched"] = count
    manifest["a1_a2_matching"]["matched_truth_pairs"] = len(pairs_seen)
    manifest["a1_a2_matching"]["reactions_per_truth_pair"] = 4179
    manifest["truth_equality_check"] = {"status":"PASS","absolute_tolerance":TRUTH_ATOL,"rows_checked":count}
    return count, utility_agg, ranked_groups


def candidate_id(values: tuple) -> str:
    payload=json.dumps(values,ensure_ascii=False,separators=(",",":"))
    return "FCG_"+hashlib.sha256(payload.encode()).hexdigest()[:20]


def rank_candidate_groups(geom: pd.DataFrame,pair_meta: dict) -> pd.DataFrame:
    keys = ["rna_context_key", "ct2a_mouse", "gl261_mouse", "reaction_id"]
    work=geom[keys+['algorithm','A1_H_dir','A2_H_dir','delta_H_dir','A1_dominant_direction_mass','A2_dominant_direction_mass','delta_dominant_direction_mass','A1_delta_v_B_width80','A2_delta_v_B_width80','delta_delta_v_B_width80','A1_ESS','A2_ESS','paired_geometry_core_finite','paired_width80_finite']].copy()
    for column in keys+['algorithm']:
        work[column]=work[column].astype('category')
    work['dominant_low']=np.minimum(work.A1_dominant_direction_mass,work.A2_dominant_direction_mass)
    work['dominant_high']=np.maximum(work.A1_dominant_direction_mass,work.A2_dominant_direction_mass)
    work['ess_row_min']=np.minimum(work.A1_ESS,work.A2_ESS)
    work['ess_row_median']=(work.A1_ESS+work.A2_ESS)/2
    work['abs_shift']=work.delta_H_dir.abs()
    work['abs_dominant_shift']=work.delta_dominant_direction_mass.abs()
    work['abs_width80_shift']=work.delta_delta_v_B_width80.abs()
    grouped=work.groupby(keys,sort=False,dropna=False,observed=True)
    ranked=grouped.agg(number_of_methods_represented=('algorithm','nunique'),finite_A1_H_dir_methods=('A1_H_dir','count'),finite_A2_H_dir_methods=('A2_H_dir','count'),finite_paired_shift_methods=('delta_H_dir','count'),finite_A1_width80_methods=('A1_delta_v_B_width80','count'),finite_A2_width80_methods=('A2_delta_v_B_width80','count'),finite_paired_width80_shift_methods=('delta_delta_v_B_width80','count'),all_methods_finite_core_geometry=('paired_geometry_core_finite','all'),all_methods_finite_width80=('paired_width80_finite','all'),A1_min=('A1_H_dir','min'),A1_max=('A1_H_dir','max'),A2_min=('A2_H_dir','min'),A2_max=('A2_H_dir','max'),shift_min=('delta_H_dir','min'),shift_max=('delta_H_dir','max'),median_absolute_A2_minus_A1_H_dir_shift=('abs_shift','median'),dominant_low=('dominant_low','min'),dominant_high=('dominant_high','max'),dominant_shift_min=('delta_dominant_direction_mass','min'),dominant_shift_max=('delta_dominant_direction_mass','max'),median_absolute_A2_minus_A1_dominant_direction_mass_shift=('abs_dominant_shift','median'),A1_width80_min=('A1_delta_v_B_width80','min'),A1_width80_max=('A1_delta_v_B_width80','max'),A2_width80_min=('A2_delta_v_B_width80','min'),A2_width80_max=('A2_delta_v_B_width80','max'),width80_shift_min=('delta_delta_v_B_width80','min'),width80_shift_max=('delta_delta_v_B_width80','max'),median_absolute_A2_minus_A1_delta_v_B_width80_shift=('abs_width80_shift','median'),minimum_ESS=('ess_row_min','min'),median_ESS=('ess_row_median','median')).reset_index()
    ranked=ranked.loc[ranked.number_of_methods_represented.eq(4)].copy()
    ranked['cross_method_range_A1_H_dir']=ranked.A1_max-ranked.A1_min
    ranked['cross_method_range_A2_H_dir']=ranked.A2_max-ranked.A2_min
    ranked['range_A2_minus_A1_H_dir_shift']=ranked.shift_max-ranked.shift_min
    ranked['cross_method_range_dominant_direction_mass']=ranked.dominant_high-ranked.dominant_low
    ranked['range_A2_minus_A1_dominant_direction_mass_shift']=ranked.dominant_shift_max-ranked.dominant_shift_min
    ranked['cross_method_range_A1_delta_v_B_width80']=ranked.A1_width80_max-ranked.A1_width80_min
    ranked['cross_method_range_A2_delta_v_B_width80']=ranked.A2_width80_max-ranked.A2_width80_min
    ranked['range_A2_minus_A1_delta_v_B_width80_shift']=ranked.width80_shift_max-ranked.width80_shift_min
    ranked['biological_comparison']=ranked.ct2a_mouse.astype(str)+'_vs_'+ranked.gl261_mouse.astype(str)
    ranked=ranked.drop(columns=['A1_min','A1_max','A2_min','A2_max','shift_min','shift_max','dominant_low','dominant_high','dominant_shift_min','dominant_shift_max','A1_width80_min','A1_width80_max','A2_width80_min','A2_width80_max','width80_shift_min','width80_shift_max'])
    registry=pd.DataFrame(pair_meta.values())
    setting_keys=["rna_context_key","ct2a_mouse","gl261_mouse","quantile_aliases"]
    require(not registry.duplicated(setting_keys+["algorithm"]).any(),"truth-selection registry has ambiguous method/setting identities")
    settings=registry.groupby(setting_keys,sort=False).algorithm.nunique().reset_index(name="number_of_methods_represented")
    settings=settings.loc[settings.number_of_methods_represented.eq(4)].drop(columns="number_of_methods_represented")
    settings=settings.rename(columns={"quantile_aliases":"truth_selection"})
    ranked=ranked.merge(settings,on=["rna_context_key","ct2a_mouse","gl261_mouse"],how="inner",validate="many_to_many")
    ranked["candidate_group_id"]=[candidate_id(tuple(row)) for row in ranked[CAND_KEYS].itertuples(index=False,name=None)]
    require(not ranked.empty,"no complete four-method candidate groups")
    ranked["primary_H_dir_range_score"]=ranked[["cross_method_range_A1_H_dir","cross_method_range_A2_H_dir"]].max(axis=1)
    ranked=ranked.sort_values(["primary_H_dir_range_score","all_methods_finite_core_geometry","finite_paired_shift_methods","minimum_ESS","candidate_group_id"],ascending=[False,False,False,False,True],kind="mergesort").reset_index(drop=True)
    ranked["geometry_rank"]=np.arange(1,len(ranked)+1)
    return ranked


def build_candidates(geom: pd.DataFrame, utility_agg: dict, manifest: dict, ranked_groups: pd.DataFrame,pair_meta:dict) -> tuple[pd.DataFrame,pd.DataFrame]:
    keys=CAND_KEYS
    groups = ranked_groups
    top = groups.head(TOP_N).copy()
    registry=pd.DataFrame(pair_meta.values())
    selection_rows=top[["candidate_group_id"]+keys].merge(registry[["truth_pair_id","evaluation_id","algorithm","rna_context_key","ct2a_mouse","gl261_mouse","quantile_aliases"]].rename(columns={"quantile_aliases":"truth_selection"}),on=["rna_context_key","ct2a_mouse","gl261_mouse","truth_selection"],how="left",validate="many_to_many")
    require(len(selection_rows)==TOP_N*4 and selection_rows.truth_pair_id.notna().all(),"candidate groups do not map to exactly four canonical held-out truth pairs")
    require(not selection_rows.duplicated(["candidate_group_id","algorithm"]).any(),"candidate truth-pair mapping is ambiguous by method")
    selected=selection_rows[["candidate_group_id","truth_pair_id","evaluation_id","algorithm","reaction_id","truth_selection"]].merge(geom,on=["evaluation_id","algorithm","reaction_id"],how="left",validate="many_to_one")
    require(selected.reaction_id.notna().all(),"candidate evaluation/reaction geometry row is missing")
    selected["method_count_in_group"] = 4
    if utility_agg:
        aggregates=[]
        for (gid,method),items in utility_agg.items():
            frame=pd.DataFrame(items)
            record={"candidate_group_id":gid,"algorithm":method,"matched_truth_case_count":len(frame)}
            for col in frame.columns:
                if col=="truth_status": continue
                record[f"mean_{col}"]=pd.to_numeric(frame[col],errors="coerce").mean()
            aggregates.append(record)
        utility=pd.DataFrame(aggregates)
        selected=selected.merge(utility,on=["candidate_group_id","algorithm"],how="left",validate="one_to_one")
    else:
        selected["matched_truth_case_count"] = 0
    selected=selected.merge(top.drop(columns=keys),on="candidate_group_id",how="left",validate="many_to_one",suffixes=("","_group"))
    groups=groups.merge(top[["candidate_group_id"]],on="candidate_group_id",how="inner",validate="one_to_one")
    group_path=GROUP_DIR=PANEL_GROUPS
    write_parquet_xz(groups, group_path)
    selected_path=PANEL_METRICS
    write_parquet_xz(selected, selected_path)
    manifest["generated_files"].extend([str(group_path.relative_to(ROOT)),str(selected_path.relative_to(ROOT))])
    manifest["row_counts"]["panelC_candidate_groups"] = int(len(groups))
    manifest["row_counts"]["panelC_candidate_metrics"] = int(len(selected))
    manifest["unique_key_counts"]["panelC_candidate_group"] = int(groups.candidate_group_id.nunique())
    manifest["candidate_ranking"]={"grouping_dimensions":["reaction_id","rna_context_key","ct2a_mouse","gl261_mouse","truth_selection (quantile_aliases)"],"rule":"descending max(cross_method_range_A1_H_dir,cross_method_range_A2_H_dir); then complete finite core geometry; paired-shift finite method count; higher minimum ESS; candidate_group_id ascending","rank_components":["primary_H_dir_range_score","cross_method_range_A1_H_dir","cross_method_range_A2_H_dir","median_absolute_A2_minus_A1_H_dir_shift","range_A2_minus_A1_H_dir_shift","cross_method_range_dominant_direction_mass","median_absolute_A2_minus_A1_dominant_direction_mass_shift","range_A2_minus_A1_dominant_direction_mass_shift","cross_method_range_A1_delta_v_B_width80","cross_method_range_A2_delta_v_B_width80","median_absolute_A2_minus_A1_delta_v_B_width80_shift","range_A2_minus_A1_delta_v_B_width80_shift","minimum_ESS","median_ESS","finite_A1_H_dir_methods","finite_A2_H_dir_methods","finite_paired_shift_methods","finite_A1_width80_methods","finite_A2_width80_methods","finite_paired_width80_shift_methods","all_methods_finite_core_geometry","all_methods_finite_width80"],"utility_used_for_ranking":False,"complete_four_method_truth_settings":int(groups.truth_selection.nunique()),"eligible_group_count_before_top_100":int(len(groups)),"top_n":len(top)}
    return top,selected


def load_existing_paired_tables(manifest:dict,pair_meta:dict) -> tuple[pd.DataFrame,int,dict,pd.DataFrame]:
    """Reuse completed geometry/utility outputs after checking their frozen keys and definitions."""
    geometry_files=sorted(GEOM_OUT.rglob("*.parquet.xz"))
    require(len(geometry_files)==4,"expected four existing geometry method partitions")
    geom=pd.concat([read_parquet_xz(path).assign(algorithm=path.parent.name.split("=",1)[1]) for path in geometry_files],ignore_index=True)
    require(len(geom)==400*4179 and not geom.duplicated(["evaluation_id","reaction_id"]).any(),"existing geometry table identity mismatch")
    ranked=rank_candidate_groups(geom,pair_meta)
    top=ranked.head(TOP_N)
    utility_files=sorted(UTILITY_OUT.glob("algorithm=*/part-*.parquet.xz"))
    require(len(utility_files)==4,"expected four existing utility method partitions")
    usecols=["truth_pair_id","evaluation_id","algorithm","rna_context_key","ct2a_mouse","gl261_mouse","truth_selection","reaction_id","m_star","truth_delta_b","A1_m0","A1_m_correct","A1_m_wrong","A2_m0","A2_m_correct","A2_m_wrong","A1_error0","A1_error_correct","A1_error_wrong","A1_gain_correct","A1_gain_wrong","A2_error0","A2_error_correct","A2_error_wrong","A2_gain_correct","A2_gain_wrong","A1_utility_values_valid","A2_utility_values_valid"]
    utility_agg=defaultdict(list); total=0; pair_counts=Counter(); truth_seen=set()
    for path in utility_files:
        frame=read_parquet_xz(path,columns=usecols)
        require(not frame.duplicated(["truth_pair_id","reaction_id"]).any(),f"duplicate utility key in {path}")
        meta_map=pair_meta
        require(frame.algorithm.eq(path.parent.name.split("=",1)[1]).all(),f"utility partition method mismatch in {path}")
        require(np.allclose(frame.m_star.to_numpy(float),np.abs(frame.truth_delta_b.to_numpy(float)),rtol=0,atol=TRUTH_ATOL),f"m_star mismatch in {path}")
        for arm in ("A1","A2"):
            expected={"error0":np.abs(frame[f"{arm}_m0"]-frame.m_star),"error_correct":np.abs(frame[f"{arm}_m_correct"]-frame.m_star),"error_wrong":np.abs(frame[f"{arm}_m_wrong"]-frame.m_star),"gain_correct":frame[f"{arm}_error0"]-frame[f"{arm}_error_correct"],"gain_wrong":frame[f"{arm}_error0"]-frame[f"{arm}_error_wrong"]}
            for name,values in expected.items():
                actual=frame[f"{arm}_{name}"].to_numpy(float); values=np.asarray(values,dtype=float)
                require(np.all(np.isnan(actual)&np.isnan(values)|np.isclose(actual,values,rtol=1e-10,atol=1e-12,equal_nan=False)),f"persisted {arm} {name} formula mismatch in {path}")
        for truth_id,g in frame.groupby("truth_pair_id",sort=False):
            meta=meta_map[truth_id]
            require(len(g)==4179 and g.evaluation_id.eq(meta["evaluation_id"]).all() and g.algorithm.eq(meta["algorithm"]).all(),f"existing utility truth-pair mapping mismatch {truth_id}")
            require(g.truth_selection.eq(meta["quantile_aliases"]).all(),f"existing utility truth-selection mismatch {truth_id}")
            pair_counts[truth_id]+=len(g)
        matched=frame.merge(top[["candidate_group_id"]+CAND_KEYS],on=CAND_KEYS,how="inner",validate="many_to_one")
        if len(matched):
            value_cols=[c for c in matched.columns if c.endswith("gain_correct") or c.endswith("gain_wrong") or c.endswith("contrast")]
            for (gid,method),group in matched.groupby(["candidate_group_id","algorithm"],sort=False):
                utility_agg[(gid,method)].extend(group[value_cols].to_dict("records"))
        total+=len(frame)
        manifest["generated_files"].append(str(path.relative_to(ROOT)))
    require(total==836*4179 and len(pair_counts)==836 and set(pair_counts.values())=={4179},"existing utility table population mismatch")
    for path in geometry_files:
        manifest["generated_files"].append(str(path.relative_to(ROOT)))
    manifest["row_counts"]["geometry_paired"]=len(geom)
    manifest["row_counts"]["utility_paired"]=total
    manifest["unique_key_counts"]["geometry_evaluation_reaction"]=len(geom)
    manifest["unique_key_counts"]["utility_truth_pair_reaction"]=total
    manifest["duplicate_key_checks"].update({"geometry_evaluation_reaction":"PASS","utility_truth_pair_reaction":"PASS"})
    manifest["a1_a2_matching"].update({"geometry_rows_matched":len(geom),"utility_rows_matched":total,"matched_truth_pairs":len(pair_counts),"reactions_per_truth_pair":4179})
    manifest["truth_equality_check"]={"status":"PASS","absolute_tolerance":TRUTH_ATOL,"rows_checked":total,"rechecked_from_persisted_m_star_and_truth_delta":True}
    return geom,total,utility_agg,ranked


def weighted_entropy(delta: np.ndarray, weight: np.ndarray) -> float:
    signs=np.where(delta>TIE_TOL,1,np.where(delta < -TIE_TOL,-1,0))
    masses=np.array([weight[signs<0].sum(),weight[signs==0].sum(),weight[signs>0].sum()],dtype=float)
    masses=masses[masses>0]
    return 0.0 if len(masses)<=1 else float(-np.sum(masses*np.log(masses))/math.log(3.0))


def build_support(top: pd.DataFrame, manifest: dict, geometry: pd.DataFrame,selected_metrics:pd.DataFrame) -> dict:
    candidates=pd.read_csv(CANDIDATE_FILE,sep="\t",compression="xz",usecols=["algorithm","tumor","ensemble_hash","projection_hash","rna_context_key","sample_index"])
    require(len(candidates)==640,"frozen candidate table must contain 640 rows")
    candidate_groups=candidates.groupby(["algorithm","tumor","rna_context_key"],sort=False)
    ensemble={}
    for key,group in candidate_groups:
        require(group.ensemble_hash.nunique()==1 and group.projection_hash.nunique()==1 and set(group.sample_index)==set(range(20)),f"ambiguous ensemble identity {key}")
        ensemble[key]=(group.ensemble_hash.iloc[0],group.projection_hash.iloc[0])
    bio=pd.read_csv(BIO0/"BRIDGEBIO0_FOUR_ARM_WEIGHTS.tsv.xz",sep="\t",compression="xz",usecols=["mouse_id","tumor","algorithm","ensemble_hash","rna_context_key","sample_index","arm","weight"])
    bio=bio.loc[bio.arm.eq("strong_anchor_baseline")]
    a20=pd.read_csv(A20/"BRIDGEA20_DUAL_ANCHOR_WEIGHTS.tsv.xz",sep="\t",compression="xz",usecols=["arm","mouse_id","tumor","algorithm","ensemble_hash","projection_hash","rna_context_key","sample_index","weight"])
    a20=a20.loc[a20.arm.eq("A2-L")]
    weightmaps={"A1":{},"A2-L":{}}
    for arm,frame in (("A1",bio),("A2-L",a20)):
        for key,g in frame.groupby(["mouse_id","tumor","algorithm","ensemble_hash","rna_context_key"],sort=False):
            g=g.sort_values("sample_index")
            require(len(g)==20 and list(g.sample_index.astype(int))==list(range(20)),f"incomplete {arm} weights {key}")
            weights=g.weight.to_numpy(float)
            require(np.isfinite(weights).all() and (weights>=0).all() and weights.sum()>0,f"invalid {arm} weights {key}")
            weightmaps[arm][key]=weights/weights.sum()
    cache=np.load(CACHE_FILE,allow_pickle=False)
    rxns=[str(x) for x in cache["rxns"].tolist()]
    require(len(rxns)==4181 and len(set(rxns))==4181 and cache["mats"].shape==(32,20,4181),"frozen flux cache identity mismatch")
    cache_index={(str(a),str(t),str(e)):i for i,(a,t,e) in enumerate(zip(cache["alg"],cache["tumor"],cache["eh"]))}
    reaction_col={r:i for i,r in enumerate(rxns)}
    candidate_rows=[]
    for r in top.itertuples(index=False):
        method_rows=selected_metrics.loc[selected_metrics.candidate_group_id.eq(r.candidate_group_id)]
        require(set(method_rows.algorithm)==set(METHODS) and len(method_rows)==4,"candidate group does not have exactly four methods")
        for ev in method_rows.itertuples(index=False):
            reaction=reaction_col[ev.reaction_id]
            for arm in ("A1","A2-L"):
                delta,joint=mats_delta(cache,cache_index,ensemble,weightmaps,ev,reaction,arm)
                candidate_rows.extend({"candidate_group_id":r.candidate_group_id,"truth_pair_id":ev.truth_pair_id,"evaluation_id":ev.evaluation_id,"reaction_id":ev.reaction_id,"method":ev.algorithm,"rna_context_key":ev.rna_context_key,"truth_selection":ev.truth_selection,"arm":arm,"support_id":i,"delta_v_B":float(delta[i]),"weight":float(joint[i]),"sign":int(1 if delta[i]>TIE_TOL else (-1 if delta[i]<-TIE_TOL else 0))} for i in range(len(delta)))
                frozen=float(ev.A1_H_dir if arm=="A1" else ev.A2_H_dir)
                reconstructed=weighted_entropy(delta,joint)
                require(math.isclose(reconstructed,frozen,rel_tol=0.0,abs_tol=HDIR_ATOL),f"H_dir support QC failed candidate={r.candidate_group_id}, eval={ev.evaluation_id}, reaction={ev.reaction_id}, arm={arm}: reconstructed={reconstructed:.16g}, frozen={frozen:.16g}")
                manifest["support_entropy_checks"].append({"candidate_group_id":r.candidate_group_id,"evaluation_id":ev.evaluation_id,"reaction_id":ev.reaction_id,"method":ev.algorithm,"arm":arm,"frozen_H_dir":frozen,"reconstructed_H_dir":reconstructed,"absolute_error":abs(reconstructed-frozen),"tolerance":HDIR_ATOL,"status":"PASS"})
    support=pd.DataFrame(candidate_rows)
    write_parquet_xz(support, SUPPORT_OUT)
    manifest["generated_files"].append(str(SUPPORT_OUT.relative_to(ROOT)))
    manifest["row_counts"]["panelC_candidate_support"] = int(len(support))
    manifest["support_entropy_qc"]={"status":"PASS","rows_checked":len(manifest["support_entropy_checks"]),"absolute_tolerance":HDIR_ATOL,"maximum_absolute_error":max(x["absolute_error"] for x in manifest["support_entropy_checks"])}
    return {"support_rows":len(support),"support_entropy_qc":manifest["support_entropy_qc"]}


def mats_delta(cache, cache_index, ensemble, weightmaps, ev, reaction, arm):
    wc_all=[]; wg_all=[]
    for tumor, mouse in (("CT2A",ev.ct2a_mouse),("GL261",ev.gl261_mouse)):
        eh,_=ensemble[(ev.algorithm,tumor,ev.rna_context_key)]
        ix=cache_index[(ev.algorithm,tumor,eh)]
        key=(mouse,tumor,ev.algorithm,eh,ev.rna_context_key)
        require(key in weightmaps[arm],f"missing {arm} weights for {key}")
        if tumor=="CT2A":
            c=cache["mats"][ix,:,reaction].astype(float); wc=weightmaps[arm][key]
        else:
            g=cache["mats"][ix,:,reaction].astype(float); wg=weightmaps[arm][key]
    delta=(c[:,None]-g[None,:]).reshape(-1)
    joint=(wc[:,None]*wg[None,:]).reshape(-1)
    joint=joint/joint.sum()
    return delta,joint


def write_manifest(manifest: dict, top: pd.DataFrame, geometry: pd.DataFrame, utility_rows: int, support_qc: dict) -> None:
    manifest["generated_files"] = sorted(set(canonical_table_path(value) for value in manifest["generated_files"]))
    manifest["generated_sha256"]={str(path.relative_to(ROOT)):sha256(path) for path in (ROOT/Path(p) for p in manifest["generated_files"])}
    manifest["table_storage"]={"format":"xz-wrapped parquet","parquet_internal_compression":"zstd","outer_compression":"xz","xz_preset":9,"persistent_uncompressed_parquet":False}
    manifest.pop("xz_archives",None)
    manifest["source_files_used"]={}
    source_roots=[PL1,A21,PL2B,A22,A23,A20,BIO0,PL2A]
    for directory in source_roots:
        m=directory.glob("*MANIFEST.json")
        path=next((p for p in m if p.name.endswith("MANIFEST.json")),None)
        if path is not None:
            manifest["source_files_used"][str(path.relative_to(ROOT))]=sha256(path)
    manifest["source_files_used"][str(CACHE_FILE.relative_to(ROOT))]=sha256(CACHE_FILE)
    manifest["source_files_used"][str(CANDIDATE_FILE.relative_to(ROOT))]=sha256(CANDIDATE_FILE)
    exact_sources=[A22/"BRIDGEA22_MATCHED_TRUTH_PAIR_REGISTRY.tsv",PL2A/"BRIDGEPL2A_TRUTH_PAIR_REGISTRY.tsv",BIO0/"BRIDGEBIO0_FOUR_ARM_WEIGHTS.tsv.xz",A20/"BRIDGEA20_DUAL_ANCHOR_WEIGHTS.tsv.xz",A21/"BRIDGEA21_EVALUATION_REGISTRY.tsv"]
    for directory,prefix in ((PL1,"BRIDGEPL1_REACTION_EVALUATION_FEATURES"),(A21,"BRIDGEA21_REACTION_EVALUATION_FEATURES"),(PL2B,"BRIDGEPL2B_CASE_OUTCOMES"),(A22,"BRIDGEA22_CASE_OUTCOMES")):
        parts_file=directory/f"{prefix}.parts.json"
        if parts_file.exists():
            part_manifest=read_json(parts_file)
            exact_sources.append(parts_file)
            exact_sources.extend(directory/part["filename"] for part in part_manifest["parts"])
    for path in exact_sources:
        safe_source(path)
        manifest["source_files_used"][str(path.relative_to(ROOT))]=sha256(path)
    manifest["candidate_groups_exported"]=int(len(top))
    manifest["row_counts"]["panelB_geometry_paired_finite"] = int(geometry.paired_geometry_core_finite.sum())
    manifest["row_counts"]["utility_paired"] = utility_rows
    qc_path=DATA/"fig4_data_qc_summary.txt"
    summary=["FIG4-DATA-B QC SUMMARY","",f"Status: {manifest['status']}",f"Geometry rows: {manifest['row_counts'].get('geometry_paired',0):,}",f"Utility rows: {utility_rows:,}",f"Panel C candidate groups: {len(top):,}",f"Panel C method metrics rows: {manifest['row_counts'].get('panelC_candidate_metrics',0):,}",f"Panel C support rows: {support_qc['support_rows']:,}",f"Truth equality: {manifest['truth_equality_check']['status']} (atol={TRUTH_ATOL:g})",f"Support entropy QC: {support_qc['support_entropy_qc']['status']} (max abs err={support_qc['support_entropy_qc']['maximum_absolute_error']:.3g}, tol={HDIR_ATOL:g})",f"Generated outputs: {len(manifest['generated_files'])+1} canonical files"]
    qc_path.write_text("\n".join(summary)+"\n")
    manifest["generated_files"].append(str(qc_path.relative_to(ROOT)))
    manifest["generated_sha256"][str(qc_path.relative_to(ROOT))]=sha256(qc_path)
    path=DATA/"fig4_data_manifest.json"
    path.write_text(json.dumps(manifest,indent=2,sort_keys=True,allow_nan=False)+"\n")


def finalize_existing_outputs(manifest: dict) -> tuple[pd.DataFrame,pd.DataFrame,int,dict]:
    """Validate and finalize complete tables when a process stopped after table creation."""
    geometry_files=sorted(GEOM_OUT.rglob("*.parquet.xz"))
    geom=pd.concat([read_parquet_xz(path).assign(algorithm=path.parent.name.split("=",1)[1]) for path in geometry_files],ignore_index=True)
    top=read_parquet_xz(PANEL_GROUPS)
    metrics=read_parquet_xz(PANEL_METRICS)
    support=read_parquet_xz(SUPPORT_OUT)
    require(len(geom)==400*4179 and not geom.duplicated(["evaluation_id","reaction_id"]).any(),"existing geometry output is incomplete or duplicated")
    require(len(top)==TOP_N and len(metrics)==TOP_N*4 and len(support)==TOP_N*4*2*400,"existing candidate outputs have unexpected cardinality")
    utility_paths=sorted(UTILITY_OUT.glob("algorithm=*/part-*.parquet.xz"))
    utility_rows=sum(read_parquet_xz_metadata(p).num_rows for p in utility_paths)
    require(utility_rows==836*4179,"existing utility output has unexpected row count")
    utility_cols=["truth_pair_id","reaction_id","m_star","truth_delta_b","A1_m0","A1_m_correct","A1_m_wrong","A2_m0","A2_m_correct","A2_m_wrong","A1_error0","A1_error_correct","A1_error_wrong","A1_gain_correct","A1_gain_wrong","A2_error0","A2_error_correct","A2_error_wrong","A2_gain_correct","A2_gain_wrong"]
    utility_keys=0
    for path in utility_paths:
        frame=read_parquet_xz(path,columns=utility_cols)
        require(not frame.duplicated(["truth_pair_id","reaction_id"]).any(),f"duplicate persisted utility key: {path}")
        require(np.allclose(frame.m_star.to_numpy(float),np.abs(frame.truth_delta_b.to_numpy(float)),rtol=0,atol=TRUTH_ATOL),f"persisted m_star mismatch: {path}")
        for arm in ("A1","A2"):
            expected={"error0":np.abs(frame[f"{arm}_m0"]-frame.m_star),"error_correct":np.abs(frame[f"{arm}_m_correct"]-frame.m_star),"error_wrong":np.abs(frame[f"{arm}_m_wrong"]-frame.m_star),"gain_correct":frame[f"{arm}_error0"]-frame[f"{arm}_error_correct"],"gain_wrong":frame[f"{arm}_error0"]-frame[f"{arm}_error_wrong"]}
            for name,values in expected.items():
                actual=frame[f"{arm}_{name}"].to_numpy(float); values=np.asarray(values,dtype=float)
                both_nan=np.isnan(actual)&np.isnan(values)
                require(np.all(both_nan|np.isclose(actual,values,rtol=1e-10,atol=1e-12,equal_nan=False)),f"persisted {arm} {name} formula mismatch: {path}")
        utility_keys+=len(frame)
    require(utility_keys==utility_rows,"persisted utility key accounting mismatch")
    check_keys=["candidate_group_id","evaluation_id","reaction_id","method","rna_context_key","arm"]
    grouped=support.groupby(check_keys+['sign'],sort=False).weight.sum().unstack(fill_value=0)
    errors=[]
    geometry_lookup=geom.set_index(["evaluation_id","reaction_id","algorithm"])
    for key,masses in grouped.iterrows():
        vals=np.array([masses.get(-1,0.0),masses.get(0,0.0),masses.get(1,0.0)],float)
        vals=vals[vals>0]
        reconstructed=0.0 if len(vals)<=1 else float(-np.sum(vals*np.log(vals))/math.log(3.0))
        arm_metric="A1_H_dir" if key[-1]=="A1" else "A2_H_dir"
        frozen=float(geometry_lookup.loc[(key[1],key[2],key[3]),arm_metric])
        err=abs(reconstructed-frozen)
        require(err<=HDIR_ATOL,f"persisted support entropy QC failed for {key}: {err}")
        errors.append({"candidate_group_id":key[0],"evaluation_id":key[1],"reaction_id":key[2],"method":key[3],"arm":key[5],"frozen_H_dir":frozen,"reconstructed_H_dir":reconstructed,"absolute_error":err,"tolerance":HDIR_ATOL,"status":"PASS"})
    manifest["status"]="PASS"
    manifest["row_counts"]={"geometry_paired":len(geom),"utility_paired":utility_rows,"panelC_candidate_groups":len(top),"panelC_candidate_metrics":len(metrics),"panelC_candidate_support":len(support),"panelB_geometry_paired_finite":int(geom.paired_geometry_core_finite.sum())}
    manifest["unique_key_counts"]={"geometry_evaluation_reaction":int(geom[["evaluation_id","reaction_id"]].drop_duplicates().shape[0]),"utility_truth_pair_reaction":utility_rows,"panelC_candidate_group":top.candidate_group_id.nunique()}
    manifest["duplicate_key_checks"]={"geometry_evaluation_reaction":"PASS","utility_truth_pair_reaction":"PASS"}
    manifest["a1_a2_matching"]={"geometry_rows_matched":len(geom),"utility_rows_matched":utility_rows,"matched_truth_pairs":836,"reactions_per_truth_pair":4179}
    manifest["truth_equality_check"]={"status":"PASS","absolute_tolerance":TRUTH_ATOL,"rows_checked":utility_rows,"rechecked_from_persisted_m_star_and_truth_delta":True}
    manifest["support_entropy_checks"]=errors
    support_qc={"support_rows":len(support),"support_entropy_qc":{"status":"PASS","rows_checked":len(errors),"absolute_tolerance":HDIR_ATOL,"maximum_absolute_error":max(x["absolute_error"] for x in errors)}}
    manifest["support_entropy_qc"]=support_qc["support_entropy_qc"]
    manifest["candidate_ranking"]={"rule":"descending max(cross_method_range_A1_H_dir,cross_method_range_A2_H_dir); then complete finite core geometry; paired-shift finite method count; higher minimum ESS; candidate_group_id ascending","utility_used_for_ranking":False,"top_n":TOP_N}
    manifest["generated_files"]=[str(p.relative_to(ROOT)) for p in sorted([*GEOM_OUT.rglob("*.parquet.xz"),*UTILITY_OUT.rglob("*.parquet.xz"),PANEL_GROUPS,PANEL_METRICS,SUPPORT_OUT])]
    manifest["table_storage"]={"format":"xz-wrapped parquet","parquet_internal_compression":"zstd","outer_compression":"xz","xz_preset":9,"persistent_uncompressed_parquet":False}
    manifest.pop("xz_archives",None)
    manifest["finite_value_counts"]={"geometry_paired_core_finite":int(geom.paired_geometry_core_finite.sum()),"geometry_paired_width80_finite":int(geom.paired_width80_finite.sum()),"utility_A1_valid":int(sum(read_pyarrow_table_xz(p,columns=["A1_utility_values_valid"]).column(0).to_numpy().sum() for p in utility_paths)),"utility_A2_valid":int(sum(read_pyarrow_table_xz(p,columns=["A2_utility_values_valid"]).column(0).to_numpy().sum() for p in utility_paths))}
    return top,geom,utility_rows,support_qc


def report_reused_stage() -> bool:
    path=DATA/"fig4_data_manifest.json"
    if not path.is_file():
        return False
    manifest=json.loads(path.read_text())
    require(manifest.get("status")=="PASS","existing FIG4-DATA-B manifest is incomplete; use --finalize-existing or --reuse-paired-tables after checking its artifacts")
    if manifest.get("stage_version") != "1.1.0":
        print("Existing FIG4-DATA-B tables predate the Panel B width80 schema; rebuilding from frozen upstream sources.", flush=True)
        return False
    for rel,digest in manifest.get("source_files_used",{}).items():
        require(sha256(safe_source(ROOT/rel))==digest,f"frozen source identity changed; refusing to reuse {rel}")
    for rel,digest in manifest.get("generated_sha256",{}).items():
        require(sha256(safe_source(ROOT/rel))==digest,f"generated Figure 4 artifact hash mismatch: {rel}")
    groups=read_parquet_xz(PANEL_GROUPS).sort_values("geometry_rank")
    print("FIG4-DATA-B already exists and all frozen source/output hashes match.")
    print(f"Geometry: {GEOM_OUT} ({manifest['row_counts']['geometry_paired']:,} rows)")
    print(f"Utility: {UTILITY_OUT} ({manifest['row_counts']['utility_paired']:,} rows)")
    print(f"Panel C groups: {PANEL_GROUPS} ({manifest['row_counts']['panelC_candidate_groups']:,} rows)")
    print(f"Panel C metrics: {PANEL_METRICS} ({manifest['row_counts']['panelC_candidate_metrics']:,} rows)")
    print(f"Panel C support: {SUPPORT_OUT} ({manifest['row_counts']['panelC_candidate_support']:,} rows)")
    print(f"QC: {DATA/'fig4_data_manifest.json'} ; {DATA/'fig4_data_qc_summary.txt'}")
    cols=["candidate_group_id","reaction_id","rna_context_key","ct2a_mouse","gl261_mouse","truth_selection","primary_H_dir_range_score","cross_method_range_A1_H_dir","cross_method_range_A2_H_dir","median_absolute_A2_minus_A1_H_dir_shift","median_absolute_A2_minus_A1_dominant_direction_mass_shift","median_absolute_A2_minus_A1_delta_v_B_width80_shift","minimum_ESS"]
    print("Top 10 Panel C candidates:")
    print(groups.head(10)[cols].to_string(index=False))
    return True


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--finalize-existing",action="store_true",help="validate and manifest already-written table outputs")
    parser.add_argument("--reuse-paired-tables",action="store_true",help="reuse validated paired geometry and utility tables while rebuilding candidate tables")
    args=parser.parse_args()
    migrate_legacy_storage()
    if not args.finalize_existing and not args.reuse_paired_tables and report_reused_stage():
        return
    for path in (GEOM_OUT,UTILITY_OUT):
        path.mkdir(parents=True,exist_ok=True)
    DATA.mkdir(parents=True,exist_ok=True)
    pl1_manifest=verify_artifacts(PL1,"BRIDGEPL1_MANIFEST.json")
    a21_manifest=verify_artifacts(A21,"BRIDGEA21_MANIFEST.json")
    a22_manifest=verify_artifacts(A22,"BRIDGEA22_MANIFEST.json")
    pl2b_manifest=verify_artifacts(PL2B,"BRIDGEPL2B_MANIFEST.json")
    a23_manifest=verify_artifacts(A23,"BRIDGEA23_MANIFEST.json")
    require(pl1_manifest.get("status")=="PL1_PREDICTABILITY_LANDSCAPE_COMPLETE" and a21_manifest.get("status")=="BRIDGE_A21_DUAL_ANCHOR_GEOMETRY_FROZEN","frozen geometry source status mismatch")
    require(a22_manifest.get("primary_lambda")==LAMBDA and pl2b_manifest.get("primary_lambda")==LAMBDA,"frozen utility lambda mismatch")
    manifest={"stage_name":"FIG4-DATA-B","stage_version":"1.1.0","status":"BUILDING","created_utc":datetime.now(timezone.utc).isoformat(),"source_directories":{k:str(v.relative_to(ROOT)) for k,v in {"A1_geometry_PL1":PL1,"A2L_geometry_A21":A21,"A1_utility_PL2B":PL2B,"A2L_utility_A22":A22,"synthesis_A23":A23,"frozen_A2_weights_A20":A20,"A1_weights_BIO0":BIO0,"truth_selection_PL2A":PL2A}.items()},"primary_lambda":LAMBDA,"truth_atol":TRUTH_ATOL,"support_entropy_atol":HDIR_ATOL,"row_counts":{},"unique_key_counts":{},"duplicate_key_checks":{},"a1_a2_matching":{},"finite_value_counts":{},"truth_equality_check":{},"support_entropy_checks":[],"source_to_output_column_mapping":{"A1_H_dir":"PL1.directional_entropy3","A2_H_dir":"A21.directional_entropy3","dominant_direction_mass":"dominant_sign_mass","direction_explained_magnitude_variance":"sign_magnitude_eta2","supported_sign_state_count":"n_supported_sign_states","non_tie_coverage":"1 - p_tie","ESS":"joint_ess","delta_v_B_width80":"PL1/A21.delta_b_width80","utility_truth":"PL2B.truth_delta_b and A22.truth_delta_b, checked then absolute value exported as m_star","correct_vs_wrong_contrast":"frozen correct_absolute_error_gain - wrong_absolute_error_gain"},"generated_files":[],"candidate_ranking":{},"table_storage":{"format":"xz-wrapped parquet","parquet_internal_compression":"zstd","outer_compression":"xz","xz_preset":9,"persistent_uncompressed_parquet":False}}
    if args.finalize_existing:
        top,geom,utility_rows,support_qc=finalize_existing_outputs(manifest)
        manifest["status"]="PASS"
        write_manifest(manifest,top,geom,utility_rows,support_qc)
        print(f"FIG4-DATA-B finalized existing tables: {utility_rows:,} utility rows; {len(top):,} candidates; support QC {support_qc['support_entropy_qc']['status']}")
        return
    if args.reuse_paired_tables:
        pair_meta,pair_ids=load_pair_registry()
        geom,utility_rows,utility_agg,ranked_groups=load_existing_paired_tables(manifest,pair_meta)
        top,selected=build_candidates(geom,utility_agg,manifest,ranked_groups,pair_meta)
        support_qc=build_support(top,manifest,geom,selected)
        manifest["status"]="PASS"
        write_manifest(manifest,top,geom,utility_rows,support_qc)
        print(f"FIG4-DATA-B status: {manifest['status']}")
        print(f"Geometry: {GEOM_OUT} ({manifest['row_counts']['geometry_paired']:,} rows)")
        print(f"Utility: {UTILITY_OUT} ({utility_rows:,} rows)")
        print(f"Panel C groups: {PANEL_GROUPS} ({len(top):,} rows)")
        print(f"Panel C metrics: {PANEL_METRICS} ({len(selected):,} rows)")
        print(f"Panel C support: {SUPPORT_OUT} ({support_qc['support_rows']:,} rows)")
        print(f"Manifest/QC: {DATA/'fig4_data_manifest.json'} ; {DATA/'fig4_data_qc_summary.txt'}")
        cols=["candidate_group_id","reaction_id","rna_context_key","ct2a_mouse","gl261_mouse","truth_selection","primary_H_dir_range_score","cross_method_range_A1_H_dir","cross_method_range_A2_H_dir","median_absolute_A2_minus_A1_H_dir_shift","median_absolute_A2_minus_A1_dominant_direction_mass_shift","median_absolute_A2_minus_A1_delta_v_B_width80_shift","minimum_ESS"]
        print("Top 10 Panel C candidates:")
        print(top.head(10)[cols].to_string(index=False))
        return
    print("Assembling paired geometry...",flush=True)
    geom=build_geometry(manifest)
    print(f"Paired geometry complete: {len(geom):,} rows",flush=True)
    pair_meta,pair_ids=load_pair_registry()
    print("Streaming and validating matched utility cases...",flush=True)
    utility_rows,utility_agg,ranked_groups=build_utility(manifest,geom,pair_meta,pair_ids)
    print(f"Paired utility complete: {utility_rows:,} rows",flush=True)
    top,selected=build_candidates(geom,utility_agg,manifest,ranked_groups,pair_meta)
    print(f"Geometry-only shortlist complete: {len(top):,} groups",flush=True)
    support_qc=build_support(top,manifest,geom,selected)
    print(f"Candidate support complete: {support_qc['support_rows']:,} rows",flush=True)
    manifest["status"]="PASS"
    write_manifest(manifest,top,geom,utility_rows,support_qc)
    print(f"FIG4-DATA-B status: {manifest['status']}")
    print(f"Geometry: {GEOM_OUT} ({manifest['row_counts']['geometry_paired']:,} rows)")
    print(f"Utility: {UTILITY_OUT} ({utility_rows:,} rows)")
    print(f"Panel C groups: {PANEL_GROUPS} ({len(top):,} rows)")
    print(f"Panel C metrics: {PANEL_METRICS} ({len(selected):,} rows)")
    print(f"Panel C support: {SUPPORT_OUT} ({support_qc['support_rows']:,} rows)")
    print(f"Manifest/QC: {DATA/'fig4_data_manifest.json'} ; {DATA/'fig4_data_qc_summary.txt'}")
    print("Top 10 Panel C candidates:")
    cols=["candidate_group_id","reaction_id","rna_context_key","ct2a_mouse","gl261_mouse","primary_H_dir_range_score","cross_method_range_A1_H_dir","cross_method_range_A2_H_dir","median_absolute_A2_minus_A1_H_dir_shift","median_absolute_A2_minus_A1_dominant_direction_mass_shift","median_absolute_A2_minus_A1_delta_v_B_width80_shift","minimum_ESS"]
    print(top.head(10)[cols].to_string(index=False))


if __name__=="__main__":
    main()
