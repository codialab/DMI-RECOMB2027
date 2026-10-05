#!/usr/bin/env python3
"""Aggregate frozen matched A1/A2-L synthetic-truth gains for Figure 3."""
from __future__ import annotations

import hashlib
import io
import json
import lzma
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
TABLE_DIR = ROOT / "reproduced/derived/fig3"
GAIN_TOL = 1e-12
IDENTITY_TOL = 5e-12
SOURCES = {
    "A1": ("outputs/dmi_bridge_pl2b_sign_only_utility_v1", "BRIDGEPL2B", "BRIDGEPL2B_MANIFEST.json", "7d9ca97a3f4ea5aa8b0b983e09ee83fc3c82f326d1b1cbf9c6137cfadf25539d", "PL2B_SIGN_ONLY_UTILITY_COMPLETE"),
    "A2-L": ("outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1", "BRIDGEA22", "BRIDGEA22_MANIFEST.json", "15b237151005538e5077a352868f873ef00466f4f5ea49307b781a8a8cd0e7cb", "BRIDGE_A22_SIGN_ONLY_UTILITY_FROZEN"),
}
A23_MANIFEST = ("outputs/dmi_bridge_a23_a1_a2_synthesis_v1/BRIDGEA23_MANIFEST.json", "b8e57074f3df23fe7708145bc01f51a19f946545b1c16e61628c2678dbcabd18")
OUT_NAMES = (
    "fig3_gain_reaction_evaluation.parquet.xz",
    "fig3_gain_cue_direction.parquet.xz",
    "fig3_gain_summary.parquet.xz",
    "fig3_gain_definition_audit.parquet.xz",
    "fig3_directional_advantage_summary.parquet.xz",
)
XZ_PRESET = 9
QUANTILES = (0.05, 0.25, 0.50, 0.75, 0.95)
CASE_COLS = [
    "truth_pair_id", "evaluation_id", "reaction_id", "truth_direction", "truth_status",
    "lambda_0_25_correct_absolute_error_gain", "lambda_0_25_wrong_absolute_error_gain",
    "lambda_0_25_random_sign_expected_gain",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def checked_repo_file(relative: str) -> Path:
    path = ROOT / relative
    resolved = path.resolve(strict=True)
    if path.is_symlink() or not resolved.is_relative_to(ROOT) or not resolved.is_file():
        raise RuntimeError(f"unsafe or missing repository input: {relative}")
    return resolved


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def load_admitted_inputs() -> tuple[set[str], dict[str, str], set[str], dict[str, str]]:
    matched_path = checked_repo_file("outputs/dmi_bridge_a23_a1_a2_synthesis_v1/BRIDGEA23_MATCHED_POPULATION_AUDIT.tsv")
    matched = read_tsv(matched_path)
    matched = matched.loc[matched["matched"].str.lower().eq("true")]
    if len(matched) != 836 or matched["truth_pair_id"].nunique() != 836:
        raise RuntimeError("A2.3 matched population is not exactly 836 distinct A1-evaluable truth pairs")
    pair_to_eval = dict(zip(matched["truth_pair_id"], matched["evaluation_id"]))

    registry_path = checked_repo_file("outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_REACTION_REGISTRY.tsv")
    reaction_registry = read_tsv(registry_path)
    reaction_ids = set(reaction_registry["reaction_id"])
    if len(reaction_ids) != 4179 or "LDH_L" in reaction_ids:
        raise RuntimeError("frozen A1/A2-L common reaction universe changed")

    source_hashes: dict[str, str] = {}
    for anchor, (directory, prefix, manifest_name, expected_hash, expected_status) in SOURCES.items():
        relative = f"{directory}/{manifest_name}"
        manifest_path = checked_repo_file(relative)
        digest = sha256(manifest_path)
        if digest != expected_hash:
            raise RuntimeError(f"frozen {anchor} manifest hash mismatch")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("status") != expected_status:
            raise RuntimeError(f"unexpected frozen {anchor} status")
        source_hashes[relative] = digest
        for name, expected in manifest.get("artifact_sha256", {}).items():
            artifact = checked_repo_file(f"{directory}/{name}")
            if sha256(artifact) != expected:
                raise RuntimeError(f"frozen {anchor} artifact hash mismatch: {name}")
        parts_path = checked_repo_file(f"{directory}/{prefix}_CASE_OUTCOMES.parts.json")
        parts = json.loads(parts_path.read_text(encoding="utf-8"))
        if parts.get("number_of_parts") != len(parts.get("parts", [])):
            raise RuntimeError(f"{anchor} case parts manifest is inconsistent")
        if any(Path(record["filename"]).name != record["filename"] for record in parts["parts"]):
            raise RuntimeError("unsafe case-part path")

    a23_relative, a23_expected = A23_MANIFEST
    a23_path = checked_repo_file(a23_relative)
    if sha256(a23_path) != a23_expected:
        raise RuntimeError("frozen A2.3 matched-population manifest hash mismatch")
    source_hashes[a23_relative] = a23_expected
    return set(pair_to_eval), pair_to_eval, reaction_ids, source_hashes


def quantile_summary(values: np.ndarray) -> dict:
    quantiles = np.quantile(values, QUANTILES)
    return {f"q{int(q * 100):02d}": float(value) for q, value in zip(QUANTILES, quantiles)}


def aggregate_anchor(anchor: str, pair_ids: set[str], pair_to_eval: dict[str, str], reaction_ids: set[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    directory, prefix, manifest_name, _, _ = SOURCES[anchor]
    parts = sorted((ROOT / directory).glob(f"{prefix}_CASE_OUTCOMES.part-*.tsv.xz"))
    if not parts:
        raise RuntimeError(f"no frozen case parts found for {anchor}")

    chunks: list[pd.DataFrame] = []
    seen_case_count = 0
    dropped_wrong_population = 0
    excluded_ties = 0
    for part in parts:
        for chunk in pd.read_csv(part, sep="\t", compression="xz", usecols=CASE_COLS, dtype=str, chunksize=100_000):
            original_rows = len(chunk)
            chunk = chunk.loc[chunk["truth_pair_id"].isin(pair_ids) & chunk["reaction_id"].isin(reaction_ids)].copy()
            dropped_wrong_population += original_rows - len(chunk)
            if chunk.empty:
                continue
            if not (chunk["evaluation_id"].to_numpy() == chunk["truth_pair_id"].map(pair_to_eval).to_numpy()).all():
                raise RuntimeError(f"{anchor} evaluation_id disagrees with matched truth registry")
            seen_case_count += len(chunk)
            for col in CASE_COLS[-3:]:
                chunk[col] = pd.to_numeric(chunk[col], errors="coerce")
            status_ok = chunk["truth_status"].isin(["TRUTH_TIE", "NON_TIE"])
            if not status_ok.all():
                raise RuntimeError(f"{anchor} contains an unknown truth_status")
            ties = chunk["truth_status"].eq("TRUTH_TIE")
            excluded_ties += int(ties.sum())
            valid = chunk.loc[~ties].copy()
            if valid.empty:
                continue
            values = valid[CASE_COLS[-3:]].to_numpy(dtype=float)
            if not np.isfinite(values).all():
                raise RuntimeError(f"{anchor} contains a non-finite gain for a non-tie matched truth")
            correct = valid[CASE_COLS[-3]].to_numpy(dtype=float)
            wrong = valid[CASE_COLS[-2]].to_numpy(dtype=float)
            random = valid[CASE_COLS[-1]].to_numpy(dtype=float)
            if (np.abs(random - 0.5 * (correct + wrong)) > IDENTITY_TOL).any():
                raise RuntimeError(f"{anchor} frozen random-sign gain identity failed")
            valid["g_info"] = correct - random
            if (np.abs(valid["g_info"].to_numpy() - 0.5 * (correct - wrong)) > IDENTITY_TOL).any():
                raise RuntimeError(f"{anchor} frozen information-advantage identity failed")
            if not valid["truth_direction"].isin(["POSITIVE", "NEGATIVE"]).all():
                raise RuntimeError(f"{anchor} truth direction has an unexpected encoding")
            chunks.append(valid[["truth_pair_id", "evaluation_id", "reaction_id", "truth_direction", "g_info", CASE_COLS[-3], CASE_COLS[-2]]].rename(columns={CASE_COLS[-3]: "correct_gain", CASE_COLS[-2]: "wrong_gain"}))

    cases = pd.concat(chunks, ignore_index=True)
    if seen_case_count != 836 * 4179:
        raise RuntimeError(f"{anchor} matched case population is {seen_case_count}, expected {836 * 4179}")
    if cases.duplicated(["truth_pair_id", "reaction_id"]).any():
        raise RuntimeError(f"{anchor} has duplicate matched synthetic-truth/reaction keys")
    per_pair = cases.groupby("truth_pair_id", sort=True)["reaction_id"].nunique()
    if len(per_pair) != 836 or not per_pair.gt(0).all():
        raise RuntimeError(f"{anchor} valid non-tie cases do not retain all matched truth pairs")

    groups = ["evaluation_id", "reaction_id"]
    wide = cases.groupby(groups, sort=True).agg(
        g_info_mean=("g_info", "mean"),
        correct_gain_mean=("correct_gain", "mean"),
        wrong_gain_mean=("wrong_gain", "mean"),
        n_valid_synthetic_truths=("truth_pair_id", "size"),
        n_truth_direction_positive=("truth_direction", lambda x: int(x.eq("POSITIVE").sum())),
        n_truth_direction_negative=("truth_direction", lambda x: int(x.eq("NEGATIVE").sum())),
    ).reset_index()
    wide.insert(0, "anchor_setting", anchor)
    wide["directional_advantage"] = wide["correct_gain_mean"] - wide["wrong_gain_mean"]
    wide = wide.sort_values(["anchor_setting", "evaluation_id", "reaction_id"], kind="stable").reset_index(drop=True)
    if wide.duplicated(["anchor_setting", "evaluation_id", "reaction_id"]).any():
        raise RuntimeError(f"{anchor} reaction-evaluation aggregate keys are not unique")
    if wide[["g_info_mean", "correct_gain_mean", "wrong_gain_mean"]].isna().any().any():
        raise RuntimeError(f"{anchor} aggregate gain table has missing plotted values")
    if not np.isfinite(wide[["g_info_mean", "correct_gain_mean", "wrong_gain_mean", "directional_advantage"]].to_numpy(dtype=float)).all():
        raise RuntimeError(f"{anchor} aggregate gain table contains non-finite values")
    if (np.abs(wide["directional_advantage"] - (wide["correct_gain_mean"] - wide["wrong_gain_mean"])) > GAIN_TOL).any():
        raise RuntimeError(f"{anchor} directional-advantage identity failed")

    long_parts = []
    for cue_direction, column in (("correct", "correct_gain"), ("wrong", "wrong_gain")):
        arm = cases[groups + [column]].groupby(groups, sort=True)[column].mean().rename("gain_mean").reset_index()
        arm.insert(0, "anchor_setting", anchor)
        arm.insert(3, "cue_direction", cue_direction)
        long_parts.append(arm)
    long = pd.concat(long_parts, ignore_index=True).sort_values(["anchor_setting", "evaluation_id", "reaction_id", "cue_direction"], kind="stable").reset_index(drop=True)
    if long.duplicated(["anchor_setting", "evaluation_id", "reaction_id", "cue_direction"]).any():
        raise RuntimeError(f"{anchor} cue endpoint aggregate keys are not unique")

    display_label = "A1" if anchor == "A1" else "A2"
    audit_rows = []
    audit_measures = (
        ("g_info", "g_info", "g_info_mean"),
        ("correct_cue_gain", "correct_gain", "correct_gain_mean"),
        ("wrong_cue_gain", "wrong_gain", "wrong_gain_mean"),
    )
    for gain_measure, case_column, reaction_evaluation_column in audit_measures:
        for definition, source in (("truth_level", cases[case_column]), ("reaction_evaluation_level", wide[reaction_evaluation_column])):
            values = source.to_numpy(dtype=float)
            audit_rows.append({
                "anchor_setting": anchor,
                "display_label": display_label,
                "gain_measure": gain_measure,
                "aggregation_definition": definition,
                "n_observations": int(len(values)),
                "mean_gain": float(values.mean()),
                "median_gain": float(np.median(values)),
                "standard_deviation": float(values.std(ddof=1)) if len(values) > 1 else 0.0,
                **quantile_summary(values),
            })
    audit = pd.DataFrame(audit_rows)
    # Compare the two existing marginal aggregates. Do not pair or pool individual
    # truth-level correct/wrong gains into a truth-level directional advantage.
    truth_contrast = float(cases["correct_gain"].mean() - cases["wrong_gain"].mean())
    reaction_eval_advantage = float(wide["directional_advantage"].mean())
    audit_rows = [{
        "anchor_setting": anchor,
        "display_label": display_label,
        "truth_level_correct_mean_minus_wrong_mean": truth_contrast,
        "reaction_evaluation_directional_advantage_mean": reaction_eval_advantage,
        "difference_between_aggregation_means": reaction_eval_advantage - truth_contrast,
        "truth_level_directional_advantage_calculated": False,
    }]
    qc = {
        "source_manifest": f"{directory}/{manifest_name}",
        "matched_truth_pairs": 836,
        "matched_case_rows_including_ties": seen_case_count,
        "excluded_truth_ties": excluded_ties,
        "valid_synthetic_truth_case_rows": len(cases),
        "rows_outside_common_population_filtered": dropped_wrong_population,
        "reaction_evaluation_rows": len(wide),
        "cue_endpoint_rows": len(long),
        "missing_aggregate_values": 0,
    }
    qc["directional_advantage_identity_max_abs_error"] = float(np.max(np.abs(wide["directional_advantage"] - (wide["correct_gain_mean"] - wide["wrong_gain_mean"]))))
    return wide, long, audit, pd.DataFrame(audit_rows), qc


def main() -> None:
    print("Verifying frozen A1/A2-L manifests and case artifacts...", flush=True)
    pair_ids, pair_to_eval, reaction_ids, source_hashes = load_admitted_inputs()
    print("Frozen inputs verified; aggregating matched synthetic-truth gains...", flush=True)
    wide_frames = []
    long_frames = []
    audit_frames = []
    directional_audit_frames = []
    directional_summary_rows = []
    anchor_qc = {}
    for anchor in SOURCES:
        print(f"Aggregating {anchor}...", flush=True)
        wide, long, audit, directional_audit, qc = aggregate_anchor(anchor, pair_ids, pair_to_eval, reaction_ids)
        wide_frames.append(wide)
        long_frames.append(long)
        audit_frames.append(audit)
        directional_audit_frames.append(directional_audit)
        values = wide["directional_advantage"].to_numpy(dtype=float)
        directional_summary_rows.append({
            "anchor_setting": anchor, "display_label": "A1" if anchor == "A1" else "A2",
            "n_observations": int(len(values)), "n_reactions": int(wide["reaction_id"].nunique()),
            "n_evaluations": int(wide["evaluation_id"].nunique()),
            "mean": float(values.mean()), "median": float(np.median(values)),
            "standard_deviation": float(values.std(ddof=1)) if len(values) > 1 else 0.0,
            "q05": float(np.quantile(values, .05)), "q95": float(np.quantile(values, .95)),
            "positive_fraction": float(np.mean(values > 0)), "zero_fraction": float(np.mean(values == 0)),
            "negative_fraction": float(np.mean(values < 0)), "missing_count": int((~np.isfinite(values)).sum()),
        })
        anchor_qc[anchor] = qc
    wide = pd.concat(wide_frames, ignore_index=True).sort_values(["anchor_setting", "evaluation_id", "reaction_id"], kind="stable").reset_index(drop=True)
    cue = pd.concat(long_frames, ignore_index=True).sort_values(["anchor_setting", "evaluation_id", "reaction_id", "cue_direction"], kind="stable").reset_index(drop=True)
    definition_audit = pd.concat(audit_frames, ignore_index=True).sort_values(["anchor_setting", "aggregation_definition"], kind="stable").reset_index(drop=True)
    directional_audit = pd.concat(directional_audit_frames, ignore_index=True).sort_values("anchor_setting", kind="stable").reset_index(drop=True)
    directional_summary = pd.DataFrame(directional_summary_rows)

    summary_rows = []
    for anchor in SOURCES:
        subset = wide.loc[wide["anchor_setting"].eq(anchor)]
        for metric, column in (("g_info", "g_info_mean"), ("correct", "correct_gain_mean"), ("wrong", "wrong_gain_mean")):
            values = subset[column].to_numpy(dtype=float)
            summary_rows.append({"anchor_setting": anchor, "display_label": "A1" if anchor == "A1" else "A2", "gain_measure": metric, "cue_direction": "information" if metric == "g_info" else metric, "n_reaction_evaluations": len(values), "mean_gain": float(values.mean()), "median_gain": float(np.median(values)), "sd_gain": float(values.std(ddof=1)) if len(values) > 1 else 0.0, **quantile_summary(values), "missing_gain_count": int((~np.isfinite(values)).sum())})
    summary = pd.DataFrame(summary_rows).sort_values(["anchor_setting", "gain_measure"], kind="stable").reset_index(drop=True)

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    frames = {OUT_NAMES[0]: wide, OUT_NAMES[1]: cue, OUT_NAMES[2]: summary, OUT_NAMES[3]: definition_audit, OUT_NAMES[4]: directional_summary}
    output_hashes = {}
    logical_hashes = {}
    sizes = {}
    for name, frame in frames.items():
        buffer = io.BytesIO()
        frame.to_parquet(buffer, index=False, engine="pyarrow", compression="zstd")
        logical_bytes = buffer.getvalue()
        compressed_bytes = lzma.compress(logical_bytes, format=lzma.FORMAT_XZ, check=lzma.CHECK_CRC64, preset=XZ_PRESET)
        path = TABLE_DIR / name
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_bytes(compressed_bytes)
        temporary.replace(path)
        output_hashes[name] = hashlib.sha256(compressed_bytes).hexdigest()
        logical_hashes[name] = hashlib.sha256(logical_bytes).hexdigest()
        sizes[name] = {"before_xz_bytes": len(logical_bytes), "after_xz_bytes": len(compressed_bytes)}

    # Read back every compressed logical table and validate its schema/row count before
    # removing only the superseded Figure 3 uncompressed derived copies.
    for name, frame in frames.items():
        compressed_path = TABLE_DIR / name
        with lzma.open(compressed_path, "rb", format=lzma.FORMAT_XZ) as source:
            restored = pd.read_parquet(io.BytesIO(source.read()), engine="pyarrow")
        if list(restored.columns) != list(frame.columns) or len(restored) != len(frame):
            raise RuntimeError(f"compressed Parquet round-trip mismatch: {name}")
    metadata = {
        "schema": "recomb.fig3.gain_distribution_metadata.v2",
        "aggregation_unit": "anchor_setting x reaction_id x evaluation_id",
        "g_info_definition": "correct_absolute_error_gain - random_sign_expected_gain; verified equal to 0.5 * (correct_gain - wrong_gain)",
        "directional_advantage_definition": "correct_gain_mean - wrong_gain_mean, calculated only after each gain is averaged within reaction_id x evaluation_id; identity checked within 1e-12",
        "directional_advantage_aggregation_audit": directional_audit.to_dict(orient="records"),
        "aggregation": "arithmetic mean across valid distinct matched synthetic truths within each reaction x evaluation",
        "compression": {"method": "XZ wrapper around Parquet bytes", "format": "XZ", "preset": XZ_PRESET, "check": "CRC64", "parquet_engine": "pyarrow", "parquet_internal_codec": "zstd"},
        "gain_definition_audit": {"truth_level": "valid non-tie synthetic truth gain values", "reaction_evaluation_level": "one aggregated gain per observed reaction x evaluation", "gain_measures": {"g_info": "correct gain minus random-sign expected gain", "correct_cue_gain": "correct-direction absolute-error gain", "wrong_cue_gain": "wrong-direction absolute-error gain"}, "quantile_probabilities": list(QUANTILES), "quantiles_are_descriptive_not_confidence_intervals": True},
        "summary_statistics": {"standard_deviation_ddof": 1, "quantile_interpolation": "linear"},
        "cue_direction_definition": "correct or wrong endpoint direction; aggregated into a separate long table",
        "truth_direction_preservation": "positive and negative direction counts among valid truths are retained in the wide table; truth ties are separately counted in QC",
        "inclusion": "836 distinct truth pairs admitted by frozen A1 evaluation QC, joined by exact truth_pair_id; common 4,179 reaction universe",
        "exclusion": "TRUTH_TIE cases excluded from gains; A1-non-evaluable truth pairs excluded by the frozen matched registry; the extra A1 LDH_L target is excluded to match the A2-L reaction universe",
        "source_manifest_sha256": source_hashes,
        "anchor_qc": anchor_qc,
        "row_counts": {name: int(len(frame)) for name, frame in frames.items()},
        "missingness": {name: {column: int(frame[column].isna().sum()) for column in frame.columns} for name, frame in frames.items()},
        "output_sha256": output_hashes,
        "uncompressed_logical_parquet_sha256": logical_hashes,
        "file_sizes": sizes,
        "outputs": list(OUT_NAMES),
    }
    metadata_path = TABLE_DIR / "fig3_gain_distribution_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    for old_name in ("fig3_gain_reaction_evaluation.parquet", "fig3_gain_cue_direction.parquet", "fig3_gain_summary.parquet"):
        old_path = TABLE_DIR / old_name
        if old_path.exists() and old_path.is_file() and not old_path.is_symlink():
            old_path.unlink()
    print(json.dumps({"status": "PASS", "output_dir": str(TABLE_DIR.relative_to(ROOT)), "row_counts": metadata["row_counts"], "output_sha256": output_hashes, "metadata": metadata_path.name}, sort_keys=True))


if __name__ == "__main__":
    main()
