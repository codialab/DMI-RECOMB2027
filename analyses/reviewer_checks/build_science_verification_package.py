#!/usr/bin/env python3
"""Build a deterministic, compact reviewer bundle from frozen source inputs.

The input repository is read-only. Output is written only to the caller's
explicit destination, which should be outside the public Git checkout.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import lzma
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
import pandas as pd

BASE_COMMIT = "af18245db1ac32918c7596d3d6660f253ae02639"
A1_REL = Path(".cache/recomb2027-dmi/frozen/outputs/dmi_bridge_pl2b_sign_only_utility_v1")
PL1_REL = Path(".cache/recomb2027-dmi/frozen/outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_MANIFEST.json")
CACHE_REL = Path(".cache/recomb2027-dmi/frozen/outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstructed/flux_cache.npz")
CATALOG_REL = Path(".cache/recomb2027-dmi/frozen/12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/candidate_fraction_table.tsv.xz")
CASE_FIELDS = [
    "truth_pair_id", "evaluation_id", "algorithm", "rna_context_key",
    "ct2a_mouse", "gl261_mouse", "reaction_id", "truth_delta_b",
    "baseline_magnitude_estimate", "lambda_0_25_correct_magnitude_estimate",
    "lambda_0_25_wrong_magnitude_estimate", "baseline_abs_error",
    "lambda_0_25_correct_absolute_error",
    "lambda_0_25_correct_absolute_error_gain",
    "lambda_0_25_wrong_absolute_error_gain",
    "lambda_0_25_random_sign_expected_gain",
]
STRING_FIELDS = CASE_FIELDS[:7]
FLOAT_FIELDS = CASE_FIELDS[7:]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def safe_file(root: Path, relative: str | Path) -> Path:
    path = root / relative
    if path.is_symlink() or not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"source path escapes source root or is a symlink: {relative}")
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def export_case_parts(source_root: Path, output: Path, source_manifest: dict) -> dict:
    base = safe_file(source_root, A1_REL / "BRIDGEPL2B_CASE_OUTCOMES.parts.json").parent
    parts_manifest = json.loads((base / "BRIDGEPL2B_CASE_OUTCOMES.parts.json").read_text())
    parts = parts_manifest.get("parts", [])
    if len(parts) != 128:
        raise ValueError("A1 source must contain the frozen 128-part case collection")
    expected_hashes = source_manifest["artifact_sha256"]
    output.mkdir()
    total = 0
    compact_parts = []
    for record in parts:
        name = record["filename"]
        source = pd.read_csv(base / record["filename"], sep="\t", compression="xz",
                             usecols=CASE_FIELDS,
                             dtype={name: str for name in STRING_FIELDS})
        src_path = safe_file(source_root, A1_REL / name)
        if digest(src_path) != record["sha256"] or expected_hashes.get(name) != record["sha256"]:
            raise ValueError(f"frozen source hash mismatch: {name}")
        if len(source) != int(record["rows"]):
            raise ValueError(f"source row count mismatch: {name}")
        compact_name = name.replace("BRIDGEPL2B_CASE_OUTCOMES", "A1_CASE_PREDICTIONS")
        compact_path = output / compact_name
        # Project the source's numeric lexemes directly. Reformatting parsed
        # floats can change the result of pandas' default parser by one ULP.
        with lzma.open(src_path, "rt", encoding="utf-8", newline="") as src_stream, \
             lzma.open(compact_path, "wt", encoding="utf-8", newline="",
                       format=lzma.FORMAT_XZ, check=lzma.CHECK_CRC64, preset=6) as sink:
            reader = csv.DictReader(src_stream, delimiter="\t")
            if reader.fieldnames is None or not set(CASE_FIELDS).issubset(reader.fieldnames):
                raise ValueError(f"source case schema lacks required columns: {name}")
            writer = csv.DictWriter(sink, fieldnames=CASE_FIELDS, delimiter="\t",
                                    lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            written = 0
            for row in reader:
                writer.writerow(row)
                written += 1
            if written != len(source):
                raise ValueError(f"streamed projection row count mismatch: {name}")
        exported = pd.read_csv(compact_path, sep="\t", compression="xz", usecols=CASE_FIELDS,
                               dtype={name: str for name in STRING_FIELDS})
        if len(exported) != len(source):
            raise ValueError(f"round-trip row count mismatch: {compact_name}")
        for field in STRING_FIELDS:
            if not source[field].equals(exported[field]):
                raise ValueError(f"round-trip identity mismatch in {field}: {compact_name}")
        for field in FLOAT_FIELDS:
            a = source[field].to_numpy(dtype=np.float64)
            b = exported[field].to_numpy(dtype=np.float64)
            if not np.array_equal(a.view(np.uint64), b.view(np.uint64), equal_nan=True):
                raise ValueError(f"round-trip float64 mismatch in {field}: {compact_name}")
        compact_parts.append({"filename": compact_name, "sha256": digest(compact_path),
                              "bytes": compact_path.stat().st_size, "rows": len(source),
                              "source_filename": name, "source_sha256": record["sha256"]})
        total += len(source)
    if total != 3_494_480:
        raise ValueError(f"unexpected A1 raw case row count: {total}")
    return {"schema": "recomb2027.compact_case_predictions.partitioned_tsv_xz.v1",
            "number_of_parts": len(compact_parts), "total_data_rows": total,
            "columns": CASE_FIELDS, "float_encoding": "original decimal lexemes preserved; xz lossless; default pandas parse parity",
            "parts": compact_parts}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True,
                        help="Read-only root containing frozen master artifacts")
    parser.add_argument("--public-root", type=Path, required=True,
                        help="Reviewer checkout at the Stage 5 base commit")
    parser.add_argument("--output", type=Path, required=True,
                        help="New local package directory outside the public checkout")
    args = parser.parse_args()
    source_root = args.source_root.resolve()
    public_root = args.public_root.resolve()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite package output: {output}")
    if output.is_relative_to(public_root):
        raise ValueError("package output must be outside the public Git checkout")
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=public_root, text=True).strip()
    ancestry = subprocess.run(["git", "merge-base", "--is-ancestor", BASE_COMMIT, "HEAD"],
                              cwd=public_root, check=False)
    if head != BASE_COMMIT and ancestry.returncode != 0:
        raise ValueError(f"public checkout must descend from Stage 5 base {BASE_COMMIT}")
    a1 = source_root / A1_REL
    source_manifest_path = safe_file(source_root, A1_REL / "BRIDGEPL2B_MANIFEST.json")
    source_manifest = json.loads(source_manifest_path.read_text())
    if source_manifest.get("status") != "PL2B_SIGN_ONLY_UTILITY_COMPLETE":
        raise ValueError("A1 PL2B manifest is not in the expected completed state")
    output.mkdir(parents=True)
    try:
        files = {}
        case_parts_dir = output / "a1_case_predictions"
        case_parts_manifest = export_case_parts(source_root, case_parts_dir, source_manifest)
        parts_manifest_path = output / "A1_CASE_PREDICTIONS.parts.json"
        parts_manifest_path.write_text(json.dumps(case_parts_manifest, indent=2, sort_keys=True) + "\n")
        files[parts_manifest_path.name] = {"sha256": digest(parts_manifest_path),
                                           "bytes": parts_manifest_path.stat().st_size,
                                           "rows": case_parts_manifest["total_data_rows"],
                                           "parts": case_parts_manifest["number_of_parts"]}
        for record in case_parts_manifest["parts"]:
            path = case_parts_dir / record["filename"]
            files[str(Path("a1_case_predictions") / record["filename"])] = {
                "sha256": record["sha256"], "bytes": record["bytes"], "rows": record["rows"]}

        registry_rel = A1_REL / "BRIDGEPL2B_DISTINCT_TRUTH_PAIR_REGISTRY.tsv"
        registry = safe_file(source_root, registry_rel)
        shutil.copyfile(registry, output / registry.name)
        files[registry.name] = {"sha256": digest(registry), "bytes": registry.stat().st_size,
                                "rows": 885, "evaluable_rows": 836}

        pl1 = safe_file(source_root, PL1_REL)
        cache = safe_file(source_root, CACHE_REL)
        catalog = safe_file(source_root, CATALOG_REL)
        for src in (pl1, cache, catalog):
            shutil.copyfile(src, output / src.name)
            files[src.name] = {"sha256": digest(src), "bytes": src.stat().st_size}
        rights_source = safe_file(public_root, "manifests/provenance/VERIFICATION_PACKAGE_RIGHTS.tsv")
        shutil.copyfile(rights_source, output / rights_source.name)
        files[rights_source.name] = {"sha256": digest(rights_source), "bytes": rights_source.stat().st_size}

        public_inputs = {}
        manifest_json = json.loads((public_root / "manifests/FIGURE_INPUTS.json").read_text())
        for item in manifest_json["inputs"]:
            public_inputs[item["destination_path"]] = item["sha256"]
        required_public = [
            "outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_MANIFEST.json",
            "outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_CASE_OUTCOMES.parts.json",
            "outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_MATCHED_TRUTH_PAIR_REGISTRY.tsv",
            "data/figure_inputs/s1/candidate_weights.tsv",
            "outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_MANIFEST.json",
            "outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_EVALUATION_REGISTRY.tsv",
            "outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_MANIFEST.json",
            "outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_REACTION_SPLIT.tsv",
            "data/figure_inputs/fig2/figure2_DG_points.csv.xz",
            "data/figure_inputs/fig3/fig3_gain_reaction_evaluation.parquet.xz",
            *[f"data/figure_inputs/fig4/fig4_geometry_paired/algorithm={method}/part.parquet.xz"
              for method in ("CORDA", "GIMME", "RIPTiDe", "iMAT")],
        ]
        checkout_hashes = {}
        for rel in required_public:
            path = public_root / rel
            if not path.is_file():
                raise FileNotFoundError(f"required public checkout input missing: {rel}")
            # Preserve each checkout dependency's committed, currently validated identity.
            indexed = public_inputs.get(rel)
            if indexed and digest(path) != indexed:
                raise ValueError(f"public figure input hash does not match FIGURE_INPUTS: {rel}")
            checkout_hashes[rel] = digest(path)

        implementation_commit = None if head == BASE_COMMIT else head
        manifest = {
            "schema": "recomb2027.science_verification_package.v1",
            "package_version": "1.0.0-local-review",
            "release_ready": False,
            "implementation_commit": implementation_commit,
            "implementation_commit_status": ("PENDING_FINAL_REVIEW_COMMIT" if implementation_commit is None else "PINNED"),
            "implementation_base_commit": BASE_COMMIT,
            "public_checkout_base_commit": BASE_COMMIT,
            "implementation_source_sha256": {
                "analyses/reviewer_checks/science_verification.py": digest(public_root / "analyses/reviewer_checks/science_verification.py"),
                "reproduce.py": digest(public_root / "reproduce.py"),
                "analyses/reviewer_checks/build_science_verification_package.py": digest(Path(__file__).resolve()),
            },
            "scope": "Level 2A, Level 2B conditional posterior, Figure 2-4 geometry verification",
            "redistribution_status": "PENDING_AUTHOR_APPROVAL; see reviewer checkout RIGHTS_REVIEW.tsv",
            "rights_review_file": "VERIFICATION_PACKAGE_RIGHTS.tsv",
            "source_root_is_not_required_at_runtime": True,
            "source_manifests": {
                "a1_manifest_sha256": digest(source_manifest_path),
                "a1_parts_manifest_sha256": digest(a1 / "BRIDGEPL2B_CASE_OUTCOMES.parts.json"),
                "pl1_manifest_sha256": digest(pl1),
            },
            "case_schema": {"columns": CASE_FIELDS, "string_columns": STRING_FIELDS,
                            "float64_columns": FLOAT_FIELDS, "encoding": "original decimal lexemes preserved; xz lossless; default pandas parse parity",
                            "partition_count": 128, "rows": 3_494_480},
            "files": files,
            "required_public_checkout_inputs": checkout_hashes,
            "scientific_expectations": {
                "level2a_numerical_parity": "PASS",
                "level2b_posterior_parity": "PASS conditional on frozen S1 weights",
                "figure3_endpoint_parity": "PASS",
                "geometry_numerical_parity": "PASS",
                "exact_usefulness_classification_parity": "one A1 cutoff difference retained",
                "independent_spearman_parity": "not exact; report floating-point rank instability",
            },
        }
        (output / "science_verification_package.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps({"status": "PASS", "output": str(output),
                          "payload_files": len(files),
                          "payload_bytes": sum(path.stat().st_size for path in output.rglob("*") if path.is_file()),
                          "case_rows": case_parts_manifest["total_data_rows"],
                          "case_parts": case_parts_manifest["number_of_parts"],
                          "manifest_sha256": digest(output / "science_verification_package.json")}, indent=2))
    except Exception:
        shutil.rmtree(output)
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
