#!/usr/bin/env python3
"""One-shot DMI-BRIDGE PL2D confirmation under the frozen amended contract."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import lzma
import math
import os
import platform
import sqlite3
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

try:
    from scripts import dmi_bridge_pl2d_confirmation_core_v1 as core
    from scripts import dmi_bridge_pl2d_hypothesis_freeze_core_v1 as hcore
    from scripts import dmi_bridge_pl2d_h_support_gate_amendment_v1 as amendment
    from scripts import dmi_bridge_pl2c_geometry_utility_core_v1 as pl2c
    from scripts import dmi_bridge_pl2c_geometry_utility_development_v1 as pl2c_adapter
    from scripts import dmi_bridge_pl2b_production_v1 as pl2b
    from scripts import dmi_bridge_pl2a_prepare_v1 as pl2a
except ModuleNotFoundError:
    import dmi_bridge_pl2d_confirmation_core_v1 as core
    import dmi_bridge_pl2d_hypothesis_freeze_core_v1 as hcore
    import dmi_bridge_pl2d_h_support_gate_amendment_v1 as amendment
    import dmi_bridge_pl2c_geometry_utility_core_v1 as pl2c
    import dmi_bridge_pl2c_geometry_utility_development_v1 as pl2c_adapter
    import dmi_bridge_pl2b_production_v1 as pl2b
    import dmi_bridge_pl2a_prepare_v1 as pl2a

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs/dmi_bridge_pl2d_confirmation_v1"
H = ROOT / "outputs/dmi_bridge_pl2d_hypothesis_freeze_v1"
HA = ROOT / "outputs/dmi_bridge_pl2d_h_support_gate_amendment_v1"
C = ROOT / "outputs/dmi_bridge_pl2c_geometry_utility_development_v1"
B = ROOT / "outputs/dmi_bridge_pl2b_sign_only_utility_v1"
A = ROOT / "outputs/dmi_bridge_pl2a_sign_only_prepare_v1"
P = ROOT / pl2a.PL1_REL
PATCH = ROOT / "dmi_bridge_pl2d_confirmation_foundation_v1.patch"
DOC = ROOT / "docs/DMI_BRIDGE_PL2D_CONFIRMATION_CONTRACT.md"
CORE = ROOT / "scripts/dmi_bridge_pl2d_confirmation_core_v1.py"
ADAPTER = Path(__file__).resolve()
PATCH_SHA = "40630a291e8e8b3c483aaebb8ad1f7b238b527a4f9a1e9ef6d2ee7dbb23264be"
FILES = (
    "BRIDGEPL2D_SOURCE_AUDIT.json",
    "BRIDGEPL2D_CONFIRMATION_EVALUATION_REACTION.tsv.xz",
    "BRIDGEPL2D_CONFIRMATION_REACTION.tsv",
    "BRIDGEPL2D_CONTEXT_RESULTS.tsv",
    "BRIDGEPL2D_PRIMARY_RESULT.json",
    "BRIDGEPL2D_SUPPORTIVE_RESULTS.json",
    "BRIDGEPL2D_BOOTSTRAP.tsv.xz",
    "BRIDGEPL2D_QC.json",
)
MEANS = tuple(hcore.REACTION_MEAN_SOURCES)
AGG = (
    "n_distinct_pairs", "n_non_tie_pairs", "n_truth_tie_pairs", "non_tie_pair_fraction",
    "aggregate_status", "mean_correct_gain", "median_correct_gain",
    "mean_information_advantage", "median_information_advantage",
    "positive_correct_gain_fraction", "positive_information_advantage_fraction",
    "directionally_useful_fraction",
)
PREDICTORS = ("sign_magnitude_eta2", "directional_entropy3", "dominant_sign_mass", "n_supported_sign_states")
ANNOTATIONS = ("proximal_set_member", "same_subsystem_as_HEX1")
EVAL_FIELDS = ("evaluation_id", "reaction_id", "algorithm", "rna_context_key", *ANNOTATIONS, *PREDICTORS, *AGG)
REACTION_FIELDS = ("reaction_id", *hcore.REACTION_MEAN_OUTPUTS.values(), "n_finite_eta2_evaluations", "n_defined_response_evaluations", *ANNOTATIONS)
CONTEXT_FIELDS = ("algorithm", "rna_context_key", "total_confirmation_reactions", "finite_primary_reaction_pairs", "rho", "rho_finite", "support_gate_ge_30", "evaluable", "direction")


def safe_file(path: Path) -> Path:
    resolved = path.resolve(strict=True)
    if path.is_symlink() or not resolved.is_relative_to(ROOT) or not resolved.is_file():
        raise RuntimeError(f"unsafe or missing repository file: {path}")
    return resolved


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with safe_file(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(safe_file(path).read_text(encoding="utf-8"))


def write_json(path: Path, value: dict) -> None:
    def finite_json(item):
        if isinstance(item, dict):
            return {key: finite_json(v) for key, v in item.items()}
        if isinstance(item, (list, tuple)):
            return [finite_json(v) for v in item]
        if isinstance(item, (float, np.floating)) and not math.isfinite(float(item)):
            return None
        return item
    path.write_bytes((json.dumps(finite_json(value), sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode("ascii"))


def write_tsv(path: Path, fields: tuple[str, ...], rows, *, compressed: bool = False) -> None:
    opener = (lambda: lzma.open(path, "wt", encoding="utf-8", newline="", preset=0, format=lzma.FORMAT_XZ)) if compressed else (lambda: path.open("w", encoding="utf-8", newline=""))
    with opener() as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: pl2b.format_value(row.get(field)) for field in fields})


def verify_bundle(base: Path, manifest_name: str, expected_sha: str) -> tuple[dict, dict[str, str]]:
    if sha(base / manifest_name) != expected_sha:
        raise RuntimeError(f"manifest hash mismatch: {manifest_name}")
    manifest = read_json(base / manifest_name)
    hashes = manifest.get("artifact_sha256")
    if not isinstance(hashes, dict):
        raise RuntimeError(f"missing artifact hashes: {manifest_name}")
    for name, digest in hashes.items():
        if Path(name).name != name or sha(base / name) != digest:
            raise RuntimeError(f"artifact hash mismatch: {name}")
    return manifest, hashes


def verify_source_admission() -> dict:
    """Hash-only admission: no confirmation row is decompressed or parsed here."""
    if sha(PATCH) != PATCH_SHA:
        raise RuntimeError("PL2D foundation patch mismatch")
    original, original_hashes = verify_bundle(H, "BRIDGEPL2DH_MANIFEST.json", core.PL2DH_MANIFEST_SHA256)
    amended, amended_hashes = verify_bundle(HA, "BRIDGEPL2DHA_MANIFEST.json", core.PL2DHA_MANIFEST_SHA256)
    c, c_hashes = verify_bundle(C, "BRIDGEPL2C_MANIFEST.json", core.PL2C_MANIFEST_SHA256)
    b, b_hashes = verify_bundle(B, "BRIDGEPL2B_MANIFEST.json", core.PL2B_MANIFEST_SHA256)
    a, a_hashes = verify_bundle(A, "BRIDGEPL2A_MANIFEST.json", core.PL2A_MANIFEST_SHA256)
    if sha(H / "BRIDGEPL2DH_CONFIRMATION_REGISTRY.tsv") != core.PL2DH_CONFIRMATION_REGISTRY_SHA256:
        raise RuntimeError("confirmation registry hash mismatch")
    if sha(HA / "BRIDGEPL2DHA_AMENDED_CONFIRMATION_CONTRACT.json") != core.PL2DHA_CONTRACT_SHA256:
        raise RuntimeError("amended contract hash mismatch")
    amendment.verify_source()
    old_contract = read_json(H / "BRIDGEPL2DH_CONFIRMATION_CONTRACT.json")
    new_contract = read_json(HA / "BRIDGEPL2DHA_AMENDED_CONFIRMATION_CONTRACT.json")
    differences = amendment.contract_diff(old_contract, new_contract)
    if amended.get("contract_semantic_diff") != differences or amended.get("original_pl2d_h_manifest_sha256") != core.PL2DH_MANIFEST_SHA256:
        raise RuntimeError("amendment provenance mismatch")
    if original.get("status") != hcore.STATUS or amended.get("pl2d_confirmation_status") != "NOT_RUN":
        raise RuntimeError("PL2D-H/amendment status mismatch")
    if c.get("status") != hcore.EXPECTED_PL2C_STATUS or c.get("split_registry_sha256") != hcore.EXPECTED_SPLIT_SHA256:
        raise RuntimeError("PL2C identity mismatch")
    if b.get("case_logical_content_sha256") != core.PL2B_CASE_LOGICAL_SHA256 or b.get("counts", {}).get("case_rows") != 3_494_480:
        raise RuntimeError("PL2B logical identity or case count mismatch")
    if len(b_hashes) != 135:
        raise RuntimeError("PL2B artifact inventory mismatch")
    b_parts = read_json(B / "BRIDGEPL2B_CASE_OUTCOMES.parts.json")
    if b_parts.get("logical_content_sha256") != core.PL2B_CASE_LOGICAL_SHA256 or len(b_parts.get("parts", [])) != 128:
        raise RuntimeError("PL2B part manifest mismatch")
    for part in b_parts["parts"]:
        name = part["filename"]
        if b_hashes.get(name) != part["sha256"] or (B / name).stat().st_size != part["bytes"]:
            raise RuntimeError(f"PL2B part identity mismatch: {name}")
    p_manifest = read_json(P / "BRIDGEPL1_MANIFEST.json")
    p_hashes = p_manifest.get("artifact_sha256", {})
    for name, digest in p_hashes.items():
        if Path(name).name != name or sha(P / name) != digest:
            raise RuntimeError(f"PL1 artifact hash mismatch: {name}")
    p_parts = read_json(P / "BRIDGEPL1_REACTION_EVALUATION_FEATURES.parts.json")
    if p_parts.get("logical_content_sha256") != core.PL1_FEATURE_LOGICAL_SHA256 or len(p_parts.get("parts", [])) != 4:
        raise RuntimeError("PL1 feature logical identity mismatch")
    for part in p_parts["parts"]:
        name = part["filename"]
        if p_hashes.get(name) != part["sha256"] or (P / name).stat().st_size != part["bytes"]:
            raise RuntimeError(f"PL1 part identity mismatch: {name}")
    if a.get("status") != "PL2A_READY_FOR_PL2B_SIGN_ONLY_UTILITY":
        raise RuntimeError("PL2A status mismatch")
    versions = {"python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__, "scipy": scipy.__version__}
    return {
        "schema": "bridge.pl2d.source_audit.v1", "status": "PASS",
        "pl2d_h_manifest_sha256": core.PL2DH_MANIFEST_SHA256,
        "pl2d_h_artifact_sha256": original_hashes,
        "pl2d_h_amendment_manifest_sha256": core.PL2DHA_MANIFEST_SHA256,
        "pl2d_h_amendment_artifact_sha256": amended_hashes,
        "pl2c_manifest_sha256": core.PL2C_MANIFEST_SHA256, "pl2c_artifact_sha256": c_hashes,
        "pl2b_manifest_sha256": core.PL2B_MANIFEST_SHA256, "pl2b_artifact_sha256": b_hashes,
        "pl2b_case_logical_content_sha256": core.PL2B_CASE_LOGICAL_SHA256,
        "pl2a_manifest_sha256": core.PL2A_MANIFEST_SHA256, "pl2a_artifact_sha256": a_hashes,
        "pl1_manifest_sha256": sha(P / "BRIDGEPL1_MANIFEST.json"), "pl1_artifact_sha256": p_hashes,
        "pl1_feature_logical_content_sha256": core.PL1_FEATURE_LOGICAL_SHA256,
        "foundation_patch_sha256": sha(PATCH), "contract_sha256": sha(DOC),
        "core_sha256": sha(CORE), "adapter_sha256": sha(ADAPTER),
        "contract_semantic_diff": differences, "versions": versions,
    }


def confirmation_population() -> set[str]:
    with safe_file(H / "BRIDGEPL2DH_CONFIRMATION_REGISTRY.tsv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    ids = [row["reaction_id"] for row in rows]
    ranks = [int(row["split_rank"]) for row in rows]
    if len(rows) != 1045 or len(set(ids)) != 1045 or ranks != list(range(1045)) or any(row["analysis_split"] != "CONFIRMATION_HOLDOUT" for row in rows):
        raise RuntimeError("frozen confirmation population mismatch")
    return set(ids)


def pair_registry() -> tuple[dict, dict[str, set[str]]]:
    with safe_file(B / "BRIDGEPL2B_DISTINCT_TRUTH_PAIR_REGISTRY.tsv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    pairs = {row["truth_pair_id"]: row for row in rows}
    if len(rows) != 885 or len(pairs) != 885:
        raise RuntimeError("distinct truth-pair registry mismatch")
    expected: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        if row["pair_evaluable"] == "true":
            expected[row["evaluation_id"]].add(row["truth_pair_id"])
    if len(expected) != 370 or sum(map(len, expected.values())) != 836:
        raise RuntimeError("evaluable truth-pair registry mismatch")
    return pairs, expected


def stream_cases(stage: Path, confirmation: set[str], pairs: dict, expected: dict[str, set[str]], qc: dict) -> sqlite3.Connection:
    parts = read_json(B / "BRIDGEPL2B_CASE_OUTCOMES.parts.json")["parts"]
    db = sqlite3.connect(stage / "confirmation_pairs.work.sqlite")
    db.execute("CREATE TABLE cases (evaluation_id TEXT, reaction_id TEXT, truth_pair_id TEXT, truth_status TEXT, correct_gain REAL, wrong_gain REAL, random_gain REAL, PRIMARY KEY(evaluation_id,reaction_id,truth_pair_id)) WITHOUT ROWID")
    counters = Counter()
    logical = hashlib.sha256()
    batch = []
    for part_index, part in enumerate(parts):
        with lzma.open(safe_file(B / part["filename"]), "rb", format=lzma.FORMAT_XZ) as handle:
            header_line = handle.readline()
            header = header_line.rstrip(b"\n").decode().split("\t")
            if header != pl2b.CASE_FIELDS:
                raise RuntimeError("PL2B case header mismatch")
            if part_index == 0:
                logical.update(header_line)
            index = {name: header.index(name) for name in ("reaction_id", "evaluation_id", "truth_pair_id", "algorithm", "rna_context_key", "truth_status", "lambda_0_25_correct_absolute_error_gain", "lambda_0_25_wrong_absolute_error_gain", "lambda_0_25_random_sign_expected_gain")}
            for line in handle:
                logical.update(line)
                fields = line.rstrip(b"\n").decode().split("\t")
                if len(fields) != len(header):
                    raise RuntimeError("PL2B case row width mismatch")
                rid = fields[index["reaction_id"]]
                if rid not in confirmation:
                    counters["development_skipped_before_outcome_parse"] += 1
                    continue
                counters["confirmation_cases"] += 1
                eid, pid = fields[index["evaluation_id"]], fields[index["truth_pair_id"]]
                pair = pairs.get(pid)
                if pair is None or pair["pair_evaluable"] != "true" or pair["evaluation_id"] != eid or pair["algorithm"] != fields[index["algorithm"]] or pair["rna_context_key"] != fields[index["rna_context_key"]] or pid not in expected[eid]:
                    raise RuntimeError("confirmation truth-pair identity mismatch")
                status = fields[index["truth_status"]]
                if status == "TRUTH_TIE":
                    counters["truth_tie_cases"] += 1
                    gains = (None, None, None)
                elif status == "NON_TIE":
                    counters["non_tie_cases"] += 1
                    gains = tuple(float(fields[index[f"lambda_0_25_{suffix}"]]) for suffix in ("correct_absolute_error_gain", "wrong_absolute_error_gain", "random_sign_expected_gain"))
                    response = pl2c.case_responses(correct_gain=gains[0], wrong_gain=gains[1], random_gain=gains[2], truth_status=status)
                    qc["maximum_random_sign_identity_deviation"] = max(qc["maximum_random_sign_identity_deviation"], abs(gains[2] - 0.5 * (gains[0] + gains[1])))
                    qc["maximum_information_advantage_identity_deviation"] = max(qc["maximum_information_advantage_identity_deviation"], abs(response["information_advantage"] - 0.5 * (gains[0] - gains[1])))
                else:
                    raise RuntimeError("invalid truth status")
                batch.append((eid, rid, pid, status, *gains))
                if len(batch) >= 5000:
                    db.executemany("INSERT INTO cases VALUES (?,?,?,?,?,?,?)", batch)
                    batch.clear()
    if batch:
        db.executemany("INSERT INTO cases VALUES (?,?,?,?,?,?,?)", batch)
    db.commit()
    if logical.hexdigest() != core.PL2B_CASE_LOGICAL_SHA256:
        raise RuntimeError("PL2B logical content mismatch")
    if counters["confirmation_cases"] != 873620 or counters["development_skipped_before_outcome_parse"] != 2620860 or db.execute("SELECT COUNT(*) FROM cases").fetchone()[0] != 873620:
        raise RuntimeError("PL2B confirmation/development accounting mismatch")
    qc["counts"].update(counters)
    qc["unique_confirmation_cases"] = True
    qc["no_quantile_alias_weighting"] = True
    return db


def eligible_pl1_key(key: tuple[str, str], confirmation: set[str], evaluable: dict[str, set[str]]) -> bool:
    """Apply PL2C's frozen evaluable-evaluation restriction before predictor conversion."""
    return key[1] in confirmation and key[0] in evaluable


def make_evaluation_rows(stage: Path, db: sqlite3.Connection, confirmation: set[str], expected: dict[str, set[str]], qc: dict) -> pd.DataFrame:
    parts = read_json(P / "BRIDGEPL1_REACTION_EVALUATION_FEATURES.parts.json")["parts"]
    cursor = db.execute("SELECT * FROM cases ORDER BY evaluation_id,reaction_id,truth_pair_id")
    current = cursor.fetchone()
    records = []
    logical = hashlib.sha256()
    previous = None
    all_rows = 0
    with lzma.open(stage / "BRIDGEPL2D_CONFIRMATION_EVALUATION_REACTION.tsv.xz", "wt", encoding="utf-8", newline="", preset=0, format=lzma.FORMAT_XZ) as out:
        writer = csv.DictWriter(out, fieldnames=EVAL_FIELDS, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        for part_index, part in enumerate(parts):
            with lzma.open(safe_file(P / part["filename"]), "rb", format=lzma.FORMAT_XZ) as handle:
                header_line = handle.readline()
                header = header_line.rstrip(b"\n").decode().split("\t")
                if part_index == 0:
                    logical.update(header_line)
                elif header != first_header:
                    raise RuntimeError("PL1 feature header mismatch")
                else:
                    pass
                if part_index == 0:
                    first_header = header
                index = {name: header.index(name) for name in EVAL_FIELDS if name not in AGG}
                for line in handle:
                    logical.update(line)
                    all_rows += 1
                    fields = line.rstrip(b"\n").decode().split("\t")
                    if len(fields) != len(header):
                        raise RuntimeError("PL1 feature row width mismatch")
                    key = (fields[index["evaluation_id"]], fields[index["reaction_id"]])
                    if previous is not None and key <= previous:
                        raise RuntimeError("PL1 key order or uniqueness mismatch")
                    previous = key
                    if not eligible_pl1_key(key, confirmation, expected):
                        continue
                    if current is None or (current[0], current[1]) != key:
                        raise RuntimeError("PL1 confirmation join missing case aggregate")
                    case_rows = []
                    while current is not None and (current[0], current[1]) == key:
                        case_rows.append({"truth_pair_id": current[2], "truth_status": current[3], "correct_gain": math.nan if current[4] is None else current[4], "wrong_gain": math.nan if current[5] is None else current[5], "random_gain": math.nan if current[6] is None else current[6]})
                        current = cursor.fetchone()
                    if {row["truth_pair_id"] for row in case_rows} != expected[key[0]]:
                        raise RuntimeError("truth-pair coverage mismatch")
                    aggregate = pl2c.aggregate_distinct_pair_responses(case_rows)
                    row = {name: fields[index[name]] for name in index}
                    row.update(aggregate)
                    writer.writerow({name: pl2b.format_value(row.get(name)) for name in EVAL_FIELDS})
                    records.append(row)
    if logical.hexdigest() != core.PL1_FEATURE_LOGICAL_SHA256 or all_rows != 1672000:
        raise RuntimeError("PL1 feature logical content or row count mismatch")
    if current is not None or len(records) != core.EXPECTED_EVALUATION_REACTION_ROWS:
        raise RuntimeError("PL1 confirmation join cardinality mismatch")
    qc["counts"]["pl1_all_rows"] = all_rows
    qc["counts"]["confirmation_evaluation_reaction_rows"] = len(records)
    qc["predictor_aggregation_independent_of_response_availability"] = True
    return pd.DataFrame.from_records(records)


def association(rows: list[dict], x: str, y: str) -> dict:
    pairs = [(r.get(x), r.get(y)) for r in rows]
    pairs = [(float(a), float(b)) for a, b in pairs if a is not None and b is not None and math.isfinite(float(a)) and math.isfinite(float(b))]
    return {"finite_reactions": len(pairs), "rho": hcore.spearman([a for a, _ in pairs], [b for _, b in pairs]) if len(pairs) >= 2 else math.nan}


def reaction_rows(frame: pd.DataFrame) -> list[dict]:
    rows = hcore.aggregate_reaction_frame(frame)
    if len(rows) != 1045:
        raise RuntimeError("reaction aggregation count mismatch")
    annotations = {}
    for rid, subset in frame.groupby("reaction_id", sort=True):
        annotations[rid] = {}
        for name in ANNOTATIONS:
            vals = set(subset[name])
            if len(vals) != 1 or next(iter(vals)) not in ("true", "false", "True", "False"):
                raise RuntimeError(f"reaction annotation varies: {rid} {name}")
            annotations[rid][name] = next(iter(vals)).lower()
    for row in rows:
        row.update(annotations[row["reaction_id"]])
    return rows


def context_rows(frame: pd.DataFrame, expected_contexts: set[tuple[str, str]]) -> list[dict]:
    rows = []
    for key, subset in frame.groupby(["algorithm", "rna_context_key"], sort=True):
        reaction = hcore.aggregate_reaction_frame(subset)
        result = association(reaction, "mean_sign_magnitude_eta2", "mean_directionally_useful_fraction")
        rho, finite = result["rho"], result["finite_reactions"]
        rows.append({"algorithm": key[0], "rna_context_key": key[1], "total_confirmation_reactions": len(reaction), "finite_primary_reaction_pairs": finite, "rho": rho, "rho_finite": math.isfinite(rho), "support_gate_ge_30": finite >= 30, "evaluable": finite >= 30 and math.isfinite(rho), "direction": "UNDEFINED" if not math.isfinite(rho) else "POSITIVE" if rho > 0 else "NEGATIVE" if rho < 0 else "ZERO"})
    if len(rows) != 16 or {(r["algorithm"], r["rna_context_key"]) for r in rows} != expected_contexts or any(r["total_confirmation_reactions"] != 1045 for r in rows):
        raise RuntimeError("context population mismatch")
    return rows


def supportive_results(rows: list[dict], bootstrap: dict) -> dict:
    x, y = "mean_sign_magnitude_eta2", "mean_directionally_useful_fraction"
    triples = core.partial_primary_controlling_coverage(rows)
    sensitivity = {}
    for name, excluded in (("exclude_proximal", ("proximal_set_member",)), ("exclude_same_subsystem_as_HEX1", ("same_subsystem_as_HEX1",)), ("exclude_either", ANNOTATIONS)):
        sensitivity[name] = association([r for r in rows if not any(r[f] == "true" for f in excluded)], x, y)
    return {
        "schema": "bridge.pl2d.supportive.v1",
        "eta2_vs_mean_information_advantage": association(rows, x, "mean_information_advantage"),
        "partial_eta2_vs_usefulness_controlling_non_tie_coverage": {"finite_reactions": triples[0], "rho": triples[1]},
        "directional_entropy3_vs_usefulness": association(rows, "mean_directional_entropy3", y),
        "dominant_sign_mass_vs_usefulness": association(rows, "mean_dominant_sign_mass", y),
        "n_supported_sign_states_vs_usefulness": association(rows, "mean_n_supported_sign_states", y),
        "pathway_sensitivity": sensitivity,
        "bootstrap": {key: bootstrap[key] for key in ("replicates", "seed", "rng", "finite_replicates", "nonfinite_replicates", "rho_q025", "rho_median", "rho_q975")},
        "bootstrap_interpretation": "descriptive computational stability across reaction benchmark units",
    }


def validate_existing(output: Path, audit: dict) -> dict:
    if output.is_symlink() or not output.is_dir() or {p.name for p in output.iterdir()} != set(FILES) | {"BRIDGEPL2D_MANIFEST.json"}:
        raise RuntimeError("existing PL2D output inventory mismatch")
    manifest = read_json(output / "BRIDGEPL2D_MANIFEST.json")
    if manifest.get("source_audit_sha256") != hashlib.sha256(pl2b.json_bytes(audit)).hexdigest() or manifest.get("status") not in ("PL2D_CONFIRMED", "PL2D_NOT_CONFIRMED", "PL2D_INSUFFICIENT_SUPPORT"):
        raise RuntimeError("existing PL2D output source/status mismatch")
    if set(manifest.get("artifact_sha256", {})) != set(FILES):
        raise RuntimeError("existing PL2D output hash inventory mismatch")
    for name, digest in manifest["artifact_sha256"].items():
        if sha(output / name) != digest:
            raise RuntimeError(f"existing PL2D output mutation: {name}")
    primary = read_json(output / "BRIDGEPL2D_PRIMARY_RESULT.json")
    if primary.get("status") != manifest["status"]:
        raise RuntimeError("existing PL2D primary/manifest status mismatch")
    return manifest


def produce(output: Path = OUTPUT) -> dict:
    outputs = ROOT / "outputs"
    if output.is_symlink() or not output.resolve().is_relative_to(outputs.resolve()) or output.resolve() == outputs.resolve():
        raise RuntimeError("PL2D output must be an in-repository outputs directory")
    audit = verify_source_admission()
    if output.exists():
        return {"status": "NO_OP_EXISTING_VALIDATED_PL2D_CONFIRMATION", "manifest": validate_existing(output, audit)}
    # The existing validator recomputes PL1/PL2B logical digests. Run it only
    # after the hash-only admission above has opened the one-shot source gate.
    pl2c_adapter.verify_sources()
    confirmation = confirmation_population()
    pairs, expected = pair_registry()
    snapshot = read_json(H / "BRIDGEPL2DH_DEVELOPMENT_SNAPSHOT.json")
    expected_contexts = {(r["algorithm"], r["rna_context_key"]) for r in snapshot["algorithm_rna_contexts"]}
    qc = {"schema": "bridge.pl2d.qc.v1", "status": "PASS", "counts": {}, "maximum_random_sign_identity_deviation": 0.0, "maximum_information_advantage_identity_deviation": 0.0}
    outputs.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="dmi_bridge_pl2d_stage_", dir=outputs) as temp:
        stage = Path(temp) / output.name
        stage.mkdir()
        write_json(stage / FILES[0], audit)
        db = stream_cases(stage, confirmation, pairs, expected, qc)
        frame = make_evaluation_rows(stage, db, confirmation, expected, qc)
        db.close()
        (stage / "confirmation_pairs.work.sqlite").unlink()
        reactions = reaction_rows(frame)
        contexts = context_rows(frame, expected_contexts)
        write_tsv(stage / "BRIDGEPL2D_CONFIRMATION_REACTION.tsv", REACTION_FIELDS, reactions)
        write_tsv(stage / "BRIDGEPL2D_CONTEXT_RESULTS.tsv", CONTEXT_FIELDS, contexts)
        finite, rho = core.primary_association(reactions)
        primary = core.classify_confirmation(finite_reactions=finite, overall_rho=rho, context_finite_counts=[r["finite_primary_reaction_pairs"] for r in contexts], context_rhos=[r["rho"] for r in contexts])
        primary_json = {"schema": "bridge.pl2d.primary.v1", **primary.__dict__, "context_results": contexts, "rule": core.contract_constants()}
        write_json(stage / "BRIDGEPL2D_PRIMARY_RESULT.json", primary_json)
        primary_hash = sha(stage / "BRIDGEPL2D_PRIMARY_RESULT.json")
        bootstrap = core.bootstrap_primary(reactions)
        write_tsv(stage / "BRIDGEPL2D_BOOTSTRAP.tsv.xz", ("replicate_index", "finite_reaction_pair_count", "rho"), ({"replicate_index": i, "finite_reaction_pair_count": int(bootstrap["finite_pair_counts"][i]), "rho": float(bootstrap["rho"][i])} for i in range(5000)), compressed=True)
        write_json(stage / "BRIDGEPL2D_SUPPORTIVE_RESULTS.json", supportive_results(reactions, bootstrap))
        if sha(stage / "BRIDGEPL2D_PRIMARY_RESULT.json") != primary_hash:
            raise RuntimeError("primary result changed after supportive analysis")
        qc["counts"].update({"confirmation_reactions": len(confirmation), "reaction_rows": len(reactions), "context_rows": len(contexts), "bootstrap_rows": 5000})
        qc.update({"primary_result_frozen_before_supportive": True, "source_admission_passed_before_confirmation_opening": True, "numerical_reference_versions": audit["versions"], "no_development_reaction_in_outputs": set(frame["reaction_id"]) == confirmation, "complete_pl1_join": True, "pathway_annotations_invariant": True, "bootstrap_replicate_count": 5000})
        if qc["maximum_random_sign_identity_deviation"] > 5e-12 or qc["maximum_information_advantage_identity_deviation"] > 5e-12 or not qc["no_development_reaction_in_outputs"]:
            raise RuntimeError("confirmation QC failure")
        write_json(stage / "BRIDGEPL2D_QC.json", qc)
        hashes = {name: sha(stage / name) for name in FILES}
        manifest = {"schema": core.SCHEMA, "status": primary.status, "source_audit_sha256": hashes["BRIDGEPL2D_SOURCE_AUDIT.json"], "artifact_sha256": hashes, "confirmation_set_opened": True, "confirmation_reaction_count": 1045, "post_opening_threshold_changes": False, "post_opening_feature_changes": False, "post_opening_lambda_changes": False, "post_opening_population_changes": False, "post_opening_status_rule_changes": False, **{name: False for name in ("solver_invoked", "fva_invoked", "sampling_invoked", "reconstruction_invoked", "new_weights_computed", "model_fitting_invoked", "reaction_selection_invoked")}}
        write_json(stage / "BRIDGEPL2D_MANIFEST.json", manifest)
        verify_source_admission()
        os.replace(stage, output)
    return {"status": primary.status, "manifest": manifest}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    result = produce(args.output)
    print(json.dumps({"status": result["status"], "primary_status": result["manifest"]["status"]}, sort_keys=True))


if __name__ == "__main__":
    main()
