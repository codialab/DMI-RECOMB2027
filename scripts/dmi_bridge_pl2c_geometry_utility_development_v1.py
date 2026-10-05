#!/usr/bin/env python3
"""Development-only PL1 geometry versus PL2B directional utility."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import lzma
import math
import os
import sqlite3
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

try:
    from scripts import dmi_bridge_pl2c_geometry_utility_core_v1 as core
    from scripts import dmi_bridge_pl2b_production_v1 as pl2b
    from scripts import dmi_bridge_pl2a_prepare_v1 as pl2a
    from scripts import dmi_bridge_pl1_storage_v1 as pl1_storage
except ModuleNotFoundError:  # direct script execution
    import dmi_bridge_pl2c_geometry_utility_core_v1 as core
    import dmi_bridge_pl2b_production_v1 as pl2b
    import dmi_bridge_pl2a_prepare_v1 as pl2a
    import dmi_bridge_pl1_storage_v1 as pl1_storage

ROOT = Path(__file__).resolve().parents[1]
OUTPUT_REL = Path("outputs/dmi_bridge_pl2c_geometry_utility_development_v1")
PATCH = Path("dmi_bridge_pl2c_geometry_utility_foundation_v1.patch")
PATCH_SHA256 = "e84ef61e4c0d77cc8282be1f5c8e8718b8d517be78ee51be4b4805533673b73c"
CONTRACT = Path("docs/DMI_BRIDGE_PL2C_GEOMETRY_UTILITY_DEVELOPMENT_CONTRACT.md")
SPLIT_NAME = "BRIDGEPL2C_REACTION_SPLIT.tsv"
DATA_NAME = "BRIDGEPL2C_DEVELOPMENT_EVALUATION_REACTION.tsv.xz"
RESPONSES = ("directionally_useful_fraction", "mean_information_advantage", "mean_correct_gain", "non_tie_pair_fraction")
AGG_FIELDS = ("n_distinct_pairs", "n_non_tie_pairs", "n_truth_tie_pairs", "non_tie_pair_fraction", "aggregate_status", "mean_correct_gain", "median_correct_gain", "mean_information_advantage", "median_information_advantage", "positive_correct_gain_fraction", "positive_information_advantage_fraction", "directionally_useful_fraction")
AUDIT_FIELDS = ("subsystem", "same_subsystem_as_HEX1", "proximal_set_member", "target_delta_degenerate", "target_abs_magnitude_degenerate", "anchor_delta_degenerate", "sign_magnitude_status", "anchor_target_coupling_status")
DATA_FIELDS = ("evaluation_id", "reaction_id", "algorithm", "rna_context_key", *AUDIT_FIELDS, *core.PRIMARY_FEATURES, *AGG_FIELDS)
EXPECTED_SPLIT_SHA256 = "7dd8676508416a13f88af184110ef40e0f8d92782cd5d62b42e7006af5316e2f"


def _file(path: Path) -> Path:
    resolved = path.resolve(strict=True)
    if not resolved.is_relative_to(ROOT) or path.is_symlink() or not resolved.is_file():
        raise RuntimeError(f"source is not a regular in-repository file: {path}")
    return resolved


def _hash(path: Path) -> str:
    digest = hashlib.sha256()
    with _file(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json(path: Path) -> dict:
    return json.loads(_file(path).read_text(encoding="utf-8"))


def _write_json(path: Path, data: dict) -> None:
    path.write_bytes(pl2b.json_bytes(data))


def _write_tsv(path: Path, fields: tuple[str, ...] | list[str], rows) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: pl2b.format_value(row.get(field)) for field in fields})


def _float(value: str) -> float:
    return float(value) if value else math.nan


def _bool(value: str) -> bool:
    if value in ("true", "True"):
        return True
    if value in ("false", "False"):
        return False
    raise RuntimeError(f"invalid frozen Boolean annotation: {value}")


def verify_sources() -> dict:
    if _hash(ROOT / PATCH) != PATCH_SHA256:
        raise RuntimeError("PL2C foundation patch hash mismatch")
    admission = pl2b.verify_admission(ROOT)
    pl2a_base = ROOT / pl2b.PL2A_REL
    pl2b_base = ROOT / pl2b.OUTPUT_REL
    if _hash(pl2a_base / "BRIDGEPL2A_MANIFEST.json") != core.PL2A_MANIFEST_SHA256:
        raise RuntimeError("PL2A manifest mismatch")
    if _hash(pl2b_base / "BRIDGEPL2B_MANIFEST.json") != core.PL2B_MANIFEST_SHA256:
        raise RuntimeError("PL2B manifest mismatch")
    manifest = _json(pl2b_base / "BRIDGEPL2B_MANIFEST.json")
    if manifest.get("status") != pl2b.COMPLETE_STATUS or manifest.get("counts", {}).get("case_rows") != pl2b.EXPECTED_CASE_ROWS:
        raise RuntimeError("PL2B status or case count mismatch")
    if manifest.get("case_logical_content_sha256") != core.PL2B_CASE_LOGICAL_SHA256:
        raise RuntimeError("PL2B logical digest mismatch")
    artifacts = manifest.get("artifact_sha256", {})
    if len(artifacts) != 135:
        raise RuntimeError("PL2B artifact inventory mismatch")
    for name, expected in artifacts.items():
        if Path(name).name != name or _hash(pl2b_base / name) != expected:
            raise RuntimeError(f"PL2B artifact hash mismatch: {name}")
    parts_name = "BRIDGEPL2B_CASE_OUTCOMES.parts.json"
    if _hash(pl2b_base / parts_name) != manifest.get("case_parts_manifest_sha256"):
        raise RuntimeError("PL2B parts manifest hash mismatch")
    parts = _json(pl2b_base / parts_name)
    if len(parts.get("parts", [])) != pl2b.PART_COUNT:
        raise RuntimeError("PL2B case part inventory mismatch")
    for part in parts["parts"]:
        name = part["filename"]
        if part["sha256"] != artifacts.get(name) or (pl2b_base / name).stat().st_size != part["bytes"]:
            raise RuntimeError(f"PL2B individual part mismatch: {name}")
    checked, _ = pl2b.inspect_case_parts(pl2b_base, expected_manifest=parts, collect_summary=False)
    if checked["logical_content_sha256"] != core.PL2B_CASE_LOGICAL_SHA256:
        raise RuntimeError("PL2B recomputed logical digest mismatch")
    if admission["pl1_consumed_artifacts"]["feature_logical_content_sha256"] != core.PL1_FEATURE_LOGICAL_SHA256:
        raise RuntimeError("PL1 feature logical digest mismatch")
    return {"schema": "bridge.pl2c.source_audit.v1", "status": "PASS", "pl2a_manifest_sha256": core.PL2A_MANIFEST_SHA256, "pl2b_manifest_sha256": core.PL2B_MANIFEST_SHA256, "pl2b_case_logical_content_sha256": checked["logical_content_sha256"], "pl2b_case_parts_manifest_sha256": manifest["case_parts_manifest_sha256"], "pl1_feature_logical_content_sha256": core.PL1_FEATURE_LOGICAL_SHA256, "pl2c_foundation_patch_sha256": PATCH_SHA256, "pl2c_contract_sha256": _hash(ROOT / CONTRACT), "upstream_admission": admission, "pl2b_artifact_sha256": artifacts}


def freeze_split(stage: Path) -> tuple[set[str], dict[str, dict], str]:
    registry = pl2b.read_tsv(ROOT / pl2b.PL2A_REL / "BRIDGEPL2A_REACTION_REGISTRY.tsv")
    ids = [row["reaction_id"] for row in registry]
    if len(ids) != core.EXPECTED_REACTIONS or len(set(ids)) != len(ids) or "HEX1" in ids:
        raise RuntimeError("frozen non-HEX1 reaction inventory mismatch")
    split = core.freeze_reaction_split(ids)
    digest = core.split_registry_sha256(split)
    if digest != EXPECTED_SPLIT_SHA256:
        raise RuntimeError("frozen reaction split hash mismatch")
    _write_tsv(stage / SPLIT_NAME, ["reaction_id", "split_hash", "split_rank", "analysis_split"], split)
    if _hash(stage / SPLIT_NAME) != digest:
        raise RuntimeError("published split registry hash mismatch")
    development = {row["reaction_id"] for row in split if row["analysis_split"] == "DEVELOPMENT"}
    if len(development) != core.EXPECTED_DEVELOPMENT_REACTIONS:
        raise RuntimeError("development split count mismatch")
    return development, {row["reaction_id"]: row for row in registry}, digest


def iter_development_cases(parts_dir: Path, parts: dict, development: set[str], pair_registry: dict, reactions: dict, counters: Counter, qc: dict):
    """Check reaction identity before touching any case outcome field."""
    for part in parts["parts"]:
        with lzma.open(_file(parts_dir / part["filename"]), "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != pl2b.CASE_FIELDS:
                raise RuntimeError("PL2B case header mismatch")
            for row in reader:
                rid = row["reaction_id"]
                if rid not in reactions:
                    raise RuntimeError("PL2B case reaction outside frozen inventory")
                if rid not in development:
                    counters["confirmation_skipped"] += 1
                    continue
                counters["development"] += 1
                pair_id = row["truth_pair_id"]
                pair = pair_registry.get(pair_id)
                if pair is None or pair["pair_evaluable"] != "true" or pair["evaluation_id"] != row["evaluation_id"]:
                    raise RuntimeError("case truth-pair identity mismatch")
                if row["algorithm"] != pair["algorithm"] or row["rna_context_key"] != pair["rna_context_key"]:
                    raise RuntimeError("case context identity mismatch")
                if any(row[field] != reactions[rid][field] for field in ("subsystem", "same_subsystem_as_HEX1", "proximal_set_member")):
                    raise RuntimeError("case reaction annotation mismatch")
                status = row["truth_status"]
                if status == "TRUTH_TIE":
                    counters["truth_tie"] += 1
                    gains = (None, None, None)
                elif status == "NON_TIE":
                    counters["non_tie"] += 1
                    gains = tuple(_float(row[f"lambda_0_25_{suffix}"]) for suffix in ("correct_absolute_error_gain", "wrong_absolute_error_gain", "random_sign_expected_gain"))
                    response = core.case_responses(correct_gain=gains[0], wrong_gain=gains[1], random_gain=gains[2], truth_status=status)
                    qc["maximum_random_sign_identity_deviation"] = max(qc["maximum_random_sign_identity_deviation"], abs(gains[2] - 0.5 * (gains[0] + gains[1])))
                    qc["maximum_information_advantage_identity_deviation"] = max(qc["maximum_information_advantage_identity_deviation"], abs(response["information_advantage"] - 0.5 * (gains[0] - gains[1])))
                    counters["directionally_useful_cases"] += int(response["directionally_useful"])
                else:
                    raise RuntimeError("unexpected truth status")
                yield (row["evaluation_id"], rid, pair_id, status, *gains)


def _materialize_development_pairs(stage: Path, development: set[str], reactions: dict, counters: Counter, qc: dict) -> tuple[sqlite3.Connection, dict]:
    pl2b_base = ROOT / pl2b.OUTPUT_REL
    parts = _json(pl2b_base / "BRIDGEPL2B_CASE_OUTCOMES.parts.json")
    pair_rows = pl2b.read_tsv(pl2b_base / "BRIDGEPL2B_DISTINCT_TRUTH_PAIR_REGISTRY.tsv")
    pair_registry = {row["truth_pair_id"]: row for row in pair_rows}
    if len(pair_rows) != 885 or len(pair_registry) != 885 or sum(row["pair_evaluable"] == "true" for row in pair_rows) != 836:
        raise RuntimeError("distinct truth-pair registry mismatch")
    expected_by_eval: dict[str, set[str]] = defaultdict(set)
    for pair in pair_rows:
        if pair["pair_evaluable"] == "true":
            expected_by_eval[pair["evaluation_id"]].add(pair["truth_pair_id"])
    if len(expected_by_eval) != 370:
        raise RuntimeError("evaluable evaluation count mismatch")
    db = sqlite3.connect(stage / "development_pairs.work.sqlite")
    db.execute("PRAGMA journal_mode=OFF")
    db.execute("PRAGMA synchronous=OFF")
    db.execute("CREATE TABLE cases (evaluation_id TEXT, reaction_id TEXT, truth_pair_id TEXT, truth_status TEXT, correct_gain REAL, wrong_gain REAL, random_gain REAL, PRIMARY KEY (evaluation_id, reaction_id, truth_pair_id)) WITHOUT ROWID")
    stream = iter_development_cases(pl2b_base, parts, development, pair_registry, reactions, counters, qc)
    buffer = []
    for record in stream:
        buffer.append(record)
        if len(buffer) == 5000:
            db.executemany("INSERT INTO cases VALUES (?,?,?,?,?,?,?)", buffer)
            buffer.clear()
    if buffer:
        db.executemany("INSERT INTO cases VALUES (?,?,?,?,?,?,?)", buffer)
    db.commit()
    if counters["development"] != 2_620_860 or counters["confirmation_skipped"] != 873_620 or sum((counters["development"], counters["confirmation_skipped"])) != 3_494_480:
        raise RuntimeError("development/confirmation case row accounting mismatch")
    if db.execute("SELECT COUNT(*) FROM cases").fetchone()[0] != counters["development"]:
        raise RuntimeError("duplicate or missing development cases")
    return db, expected_by_eval


def _summary(values: np.ndarray, mask: np.ndarray) -> dict:
    good = values[mask & np.isfinite(values)]
    return {"mean": float(np.mean(good)) if good.size else math.nan, "median": float(np.median(good)) if good.size else math.nan, "finite_count": int(good.size)}


def _analysis_rows(stage: Path, db: sqlite3.Connection, expected_by_eval: dict, development: set[str], counters: Counter) -> tuple[dict, dict]:
    n = len(expected_by_eval) * len(development)
    arrays = {field: np.full(n, math.nan, dtype=float) for field in (*core.PRIMARY_FEATURES, *RESPONSES)}
    context = np.empty(n, dtype=object)
    proximal = np.zeros(n, dtype=bool)
    same = np.zeros(n, dtype=bool)
    eta_degenerate = np.zeros(n, dtype=bool)
    all_predictor_values = {field: [] for field in core.PRIMARY_FEATURES}
    cursor = db.execute("SELECT evaluation_id,reaction_id,truth_pair_id,truth_status,correct_gain,wrong_gain,random_gain FROM cases ORDER BY evaluation_id,reaction_id,truth_pair_id")
    current = cursor.fetchone()
    previous_pl1 = None
    output_count = 0
    group_pair_counts = Counter()
    with lzma.open(stage / DATA_NAME, "wt", encoding="utf-8", newline="", preset=0, format=lzma.FORMAT_XZ) as handle:
        writer = csv.DictWriter(handle, fieldnames=DATA_FIELDS, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        for source in pl1_storage.iter_feature_rows(ROOT / pl2a.PL1_REL):
            key = (source["evaluation_id"], source["reaction_id"])
            if previous_pl1 is not None and key <= previous_pl1:
                raise RuntimeError("PL1 predictor key order or uniqueness mismatch")
            previous_pl1 = key
            counters["pl1_all_rows"] += 1
            if key[1] not in development:
                continue
            counters["pl1_development_rows"] += 1
            feature_values = {field: _float(source[field]) for field in core.PRIMARY_FEATURES}
            for field in core.PRIMARY_FEATURES:
                all_predictor_values[field].append(feature_values[field])
            if key[0] not in expected_by_eval:
                continue
            if current is None or (current[0], current[1]) != key:
                raise RuntimeError("PL1 predictor coverage missing a development aggregate")
            pair_records = []
            while current is not None and (current[0], current[1]) == key:
                pair_records.append({"truth_pair_id": current[2], "truth_status": current[3], "correct_gain": math.nan if current[4] is None else current[4], "wrong_gain": math.nan if current[5] is None else current[5], "random_gain": math.nan if current[6] is None else current[6]})
                current = cursor.fetchone()
            if {record["truth_pair_id"] for record in pair_records} != expected_by_eval[key[0]]:
                raise RuntimeError("evaluation/reaction truth-pair coverage mismatch")
            agg = core.aggregate_distinct_pair_responses(pair_records)
            group_pair_counts[agg["n_distinct_pairs"]] += 1
            counters["aggregates_defined"] += agg["aggregate_status"] == "DEFINED"
            counters["aggregates_no_directional_truth"] += agg["aggregate_status"] == "NO_DIRECTIONAL_TRUTH"
            row = {field: source[field] for field in ("evaluation_id", "reaction_id", "algorithm", "rna_context_key", *AUDIT_FIELDS)}
            row.update(feature_values)
            row.update(agg)
            writer.writerow({field: pl2b.format_value(row.get(field)) for field in DATA_FIELDS})
            for field in core.PRIMARY_FEATURES:
                arrays[field][output_count] = feature_values[field]
            for field in RESPONSES:
                arrays[field][output_count] = agg[field]
            context[output_count] = (source["algorithm"], source["rna_context_key"])
            proximal[output_count] = _bool(source["proximal_set_member"])
            same[output_count] = _bool(source["same_subsystem_as_HEX1"])
            eta_degenerate[output_count] = source["sign_magnitude_status"] != "OK"
            output_count += 1
    if current is not None or output_count != n or counters["pl1_all_rows"] != 1_672_000 or counters["pl1_development_rows"] != 1_254_000:
        raise RuntimeError("PL1 join or aggregate row accounting mismatch")
    if sum(group_pair_counts[k] * k for k in group_pair_counts) != 2_620_860:
        raise RuntimeError("aggregate distinct pair accounting mismatch")
    arrays = {key: value[:output_count] for key, value in arrays.items()}
    meta = {"context": context[:output_count], "proximal": proximal[:output_count], "same": same[:output_count], "eta_degenerate": eta_degenerate[:output_count], "all_predictor_values": all_predictor_values, "group_pair_counts": dict(sorted(group_pair_counts.items()))}
    return arrays, meta


def _summaries(stage: Path, arrays: dict, meta: dict) -> dict:
    n = len(arrays[RESPONSES[0]])
    all_mask = np.ones(n, dtype=bool)
    scopes = {"PRIMARY_ALL_DEVELOPMENT": all_mask, "SENSITIVITY_EXCLUDE_PROXIMAL": ~meta["proximal"], "SENSITIVITY_EXCLUDE_SAME_SUBSYSTEM": ~meta["same"], "SENSITIVITY_EXCLUDE_EITHER": ~(meta["proximal"] | meta["same"])}
    assoc = []
    primary_rho = {}
    for scope, mask in scopes.items():
        for feature in core.PRIMARY_FEATURES:
            for response in RESPONSES:
                result = core.spearman_finite(arrays[feature][mask], arrays[response][mask])
                assoc.append({"scope": scope, "feature": feature, "response": response, **result})
                if scope == "PRIMARY_ALL_DEVELOPMENT":
                    primary_rho[(feature, response)] = result["rho"]
    _write_tsv(stage / "BRIDGEPL2C_FEATURE_ASSOCIATIONS.tsv", ["scope", "feature", "response", "rho", "finite_count", "nonfinite_count"], assoc)
    context_rows = []
    keys = sorted(set(meta["context"]))
    for algorithm, rna_key in keys:
        mask = np.fromiter((item == (algorithm, rna_key) for item in meta["context"]), dtype=bool, count=n)
        for feature in core.PRIMARY_FEATURES:
            for response in RESPONSES:
                result = core.spearman_finite(arrays[feature][mask], arrays[response][mask])
                rho, overall = result["rho"], primary_rho[(feature, response)]
                direction = "UNDEFINED" if not math.isfinite(rho) else ("POSITIVE" if rho > 0 else "NEGATIVE" if rho < 0 else "ZERO")
                stable = math.isfinite(rho) and math.isfinite(overall) and ((rho > 0) == (overall > 0)) and ((rho < 0) == (overall < 0))
                context_rows.append({"algorithm": algorithm, "rna_context_key": rna_key, "feature": feature, "response": response, **result, "direction": direction, "direction_matches_overall": stable})
    _write_tsv(stage / "BRIDGEPL2C_CONTEXT_ASSOCIATIONS.tsv", ["algorithm", "rna_context_key", "feature", "response", "rho", "finite_count", "nonfinite_count", "direction", "direction_matches_overall"], context_rows)
    edges = {feature: core.linear_quartile_edges(meta["all_predictor_values"][feature]) for feature in core.PRIMARY_FEATURES if feature != "n_supported_sign_states"}
    quartile_rows = []
    for feature in core.PRIMARY_FEATURES:
        if feature == "n_supported_sign_states":
            levels = sorted({int(x) for x in meta["all_predictor_values"][feature] if math.isfinite(x)})
            labels = [str(int(v)) if math.isfinite(v) else "NONFINITE" for v in arrays[feature]]
            bins = [str(v) for v in levels] + ["NONFINITE"]
        else:
            labels = [core.quartile_label(v, edges[feature]) for v in arrays[feature]]
            if feature == "sign_magnitude_eta2":
                labels = ["DEGENERATE" if meta["eta_degenerate"][i] else label for i, label in enumerate(labels)]
            bins = ["Q1", "Q2", "Q3", "Q4", "NONFINITE"] + (["DEGENERATE"] if feature == "sign_magnitude_eta2" else [])
        label_arr = np.asarray(labels, dtype=object)
        for bin_name in bins:
            mask = label_arr == bin_name
            row = {"feature": feature, "bin": bin_name, "cut_q25": edges.get(feature, (math.nan,)*3)[0], "cut_q50": edges.get(feature, (math.nan,)*3)[1], "cut_q75": edges.get(feature, (math.nan,)*3)[2], "evaluation_reaction_rows": int(mask.sum()), "defined_directional_rows": int(np.sum(mask & np.isfinite(arrays["directionally_useful_fraction"])))}
            for response in RESPONSES:
                summary = _summary(arrays[response], mask)
                row[f"{response}_mean"] = summary["mean"]
                row[f"{response}_median"] = summary["median"]
            quartile_rows.append(row)
    quartile_fields = ["feature", "bin", "cut_q25", "cut_q50", "cut_q75", "evaluation_reaction_rows", "defined_directional_rows", *[f"{response}_{stat}" for response in RESPONSES for stat in ("mean", "median")]]
    _write_tsv(stage / "BRIDGEPL2C_FEATURE_QUARTILES.tsv", quartile_fields, quartile_rows)
    eta_labels = np.asarray(["DEGENERATE" if meta["eta_degenerate"][i] else core.quartile_label(v, edges["sign_magnitude_eta2"]) for i, v in enumerate(arrays["sign_magnitude_eta2"])], dtype=object)
    entropy_labels = np.asarray([core.quartile_label(v, edges["directional_entropy3"]) for v in arrays["directional_entropy3"]], dtype=object)
    map_rows = []
    for eta in ("Q1", "Q2", "Q3", "Q4", "NONFINITE", "DEGENERATE"):
        for entropy in ("Q1", "Q2", "Q3", "Q4", "NONFINITE"):
            mask = (eta_labels == eta) & (entropy_labels == entropy)
            if not mask.any() and (eta not in ("Q1", "Q2", "Q3", "Q4") or entropy == "NONFINITE"):
                continue
            row = {"eta2_bin": eta, "entropy3_bin": entropy, "evaluation_reaction_rows": int(mask.sum()), "defined_directional_rows": int(np.sum(mask & np.isfinite(arrays["directionally_useful_fraction"])))}
            for response in RESPONSES:
                summary = _summary(arrays[response], mask)
                row[f"{response}_mean"] = summary["mean"]
                row[f"{response}_median"] = summary["median"]
            map_rows.append(row)
    map_fields = ["eta2_bin", "entropy3_bin", "evaluation_reaction_rows", "defined_directional_rows", *[f"{response}_{stat}" for response in RESPONSES for stat in ("mean", "median")]]
    _write_tsv(stage / "BRIDGEPL2C_ETA2_ENTROPY_MAP.tsv", map_fields, map_rows)
    stability = {}
    for feature in core.PRIMARY_FEATURES:
        for response in RESPONSES:
            rows = [row for row in context_rows if row["feature"] == feature and row["response"] == response]
            stability[f"{feature}:{response}"] = {"contexts": len(rows), "positive": sum(row["direction"] == "POSITIVE" for row in rows), "negative": sum(row["direction"] == "NEGATIVE" for row in rows), "zero": sum(row["direction"] == "ZERO" for row in rows), "undefined": sum(row["direction"] == "UNDEFINED" for row in rows), "matches_overall": sum(row["direction_matches_overall"] for row in rows)}
    return {"quartile_edges_from_all_development_pl1": {k: list(v) for k, v in edges.items()}, "context_stability": stability, "primary_coverage": _summary(arrays["non_tie_pair_fraction"], all_mask), "primary_utility": _summary(arrays["directionally_useful_fraction"], all_mask), "map_rows": len(map_rows)}


def validate_existing_output(output: Path, fingerprint: dict) -> dict:
    manifest_path = output / "BRIDGEPL2C_MANIFEST.json"
    if not manifest_path.is_file():
        raise RuntimeError("existing PL2C output lacks manifest")
    manifest = _json(manifest_path)
    if manifest.get("status") != core.COMPLETE_STATUS or manifest.get("fingerprint") != fingerprint or manifest.get("pl2d_confirmation_status") != core.PL2D_STATUS:
        raise RuntimeError("existing PL2C output identity/status mismatch")
    expected = manifest.get("artifact_sha256", {})
    if len(expected) != 9 or set(path.name for path in output.iterdir()) != set(expected) | {"BRIDGEPL2C_MANIFEST.json"}:
        raise RuntimeError("existing PL2C artifact inventory mismatch")
    for name, digest in expected.items():
        if _hash(output / name) != digest:
            raise RuntimeError(f"existing PL2C artifact hash mismatch: {name}")
    if expected[SPLIT_NAME] != EXPECTED_SPLIT_SHA256:
        raise RuntimeError("existing PL2C split hash mismatch")
    return manifest


def produce(output: Path) -> dict:
    output = output.resolve()
    if not output.is_relative_to(ROOT / "outputs") or output == ROOT / "outputs":
        raise RuntimeError("PL2C output must remain under repository outputs")
    source_audit = verify_sources()
    fingerprint = {"schema": core.SCHEMA, "source_sha256": {key: source_audit[key] for key in ("pl2a_manifest_sha256", "pl2b_manifest_sha256", "pl2b_case_logical_content_sha256", "pl1_feature_logical_content_sha256", "pl2c_foundation_patch_sha256", "pl2c_contract_sha256")}, "implementation_sha256": {name: _hash(ROOT / name) for name in ("scripts/dmi_bridge_pl2c_geometry_utility_core_v1.py", "scripts/dmi_bridge_pl2c_geometry_utility_development_v1.py")}}
    if output.exists():
        return {"status": "NO_OP_EXISTING_IDENTICAL_PL2C", "manifest": validate_existing_output(output, fingerprint)}
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="dmi_bridge_pl2c_publish_", dir=output.parent) as temporary:
        stage = Path(temporary)
        development, reactions, split_digest = freeze_split(stage)
        _write_json(stage / "BRIDGEPL2C_SOURCE_AUDIT.json", source_audit)
        contract = core.analysis_contract()
        contract.update({"split_registry_sha256": split_digest, "quartile_population": "all_1254000_development_PL1_predictor_rows_before_response_conditioning", "pathway_sensitivity": ["exclude_proximal", "exclude_same_subsystem", "exclude_either"]})
        _write_json(stage / "BRIDGEPL2C_ANALYSIS_CONTRACT.json", contract)
        counters = Counter()
        qc = {"maximum_random_sign_identity_deviation": 0.0, "maximum_information_advantage_identity_deviation": 0.0}
        db, expected_by_eval = _materialize_development_pairs(stage, development, reactions, counters, qc)
        try:
            arrays, meta = _analysis_rows(stage, db, expected_by_eval, development, counters)
        finally:
            db.close()
            (stage / "development_pairs.work.sqlite").unlink()
        qc.update(_summaries(stage, arrays, meta))
        qc.update({"schema": "bridge.pl2c.qc.v1", "status": "PASS", "counts": dict(sorted(counters.items())), "aggregate_pair_multiplicity": meta["group_pair_counts"], "split_registry_sha256": split_digest, "confirmation_outcomes_analyzed": False, "pl2d_confirmation_status": core.PL2D_STATUS, "forbidden_work": {name: False for name in ("solver_invoked", "optimization_invoked", "fva_invoked", "sampling_invoked", "reconstruction_invoked", "new_weights_generated", "lambda_tuning", "feature_selection", "reaction_ranking", "confirmation_outcome_analysis")}})
        _write_json(stage / "BRIDGEPL2C_QC.json", qc)
        names = ("BRIDGEPL2C_SOURCE_AUDIT.json", "BRIDGEPL2C_ANALYSIS_CONTRACT.json", SPLIT_NAME, DATA_NAME, "BRIDGEPL2C_FEATURE_ASSOCIATIONS.tsv", "BRIDGEPL2C_CONTEXT_ASSOCIATIONS.tsv", "BRIDGEPL2C_FEATURE_QUARTILES.tsv", "BRIDGEPL2C_ETA2_ENTROPY_MAP.tsv", "BRIDGEPL2C_QC.json")
        hashes = {name: _hash(stage / name) for name in names}
        manifest = {"schema": core.SCHEMA, "status": core.COMPLETE_STATUS, "fingerprint": fingerprint, "artifact_sha256": hashes, "counts": dict(sorted(counters.items())), "reaction_split_counts": {"DEVELOPMENT": 3135, "CONFIRMATION_HOLDOUT": 1045}, "split_registry_sha256": split_digest, "confirmation_outcomes_analyzed": False, "pl2d_confirmation_status": core.PL2D_STATUS}
        _write_json(stage / "BRIDGEPL2C_MANIFEST.json", manifest)
        os.replace(stage, output)
    validate_existing_output(output, fingerprint)
    return {"status": core.COMPLETE_STATUS, "manifest": manifest}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / OUTPUT_REL)
    args = parser.parse_args()
    result = produce(args.output)
    print(json.dumps({"status": result["status"], "counts": result["manifest"]["counts"], "split_registry_sha256": result["manifest"]["split_registry_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    main()
