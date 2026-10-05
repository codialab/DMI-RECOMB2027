#!/usr/bin/env python3
"""Prepare the outcome-blind DMI-BRIDGE-PL2A registry.

This adapter deliberately stops before any reaction-B truth delta, posterior
error, gain, or PL2 scientific outcome is computed.  Truth identities are
selected from HEX1 and the frozen strong-anchor weights only.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import lzma
import math
import os
import shutil
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
try:
    from scripts import dmi_bridge_pl1_storage_v1 as pl1_storage
except ModuleNotFoundError:  # direct script execution
    import dmi_bridge_pl1_storage_v1 as pl1_storage

ROOT = Path(__file__).resolve().parents[1]
OUTPUT_REL = Path("outputs/dmi_bridge_pl2a_sign_only_prepare_v1")
PL1_REL = Path("outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1")
PL1_MANIFEST_REL = PL1_REL / "BRIDGEPL1_MANIFEST.json"
EXPECTED_HISTORICAL_PL1_MANIFEST_SHA256 = "7465c8bc8602b4b1b7a623a2f4c9322c77f5d5ad1c070e7bab5aa6668c411e9c"
EXPECTED_LIVE_PL1_MANIFEST_SHA256 = "b0f81620ff24f5791dd964a8063416784db86a358b0464b8c68a504e84e5f529"
PATCH_REL = Path("dmi_bridge_pl2a_sign_only_foundation_v1.patch")
EXPECTED_PATCH_SHA256 = "8b402487fa1bdc90b683d8761c0789bf00c15f3e80ede5d37d424a4e8095225f"
QUANTILES = (0.10, 0.50, 0.90)
METHODS = ("CORDA", "GIMME", "RIPTiDe", "iMAT")
TUMORS = ("CT2A", "GL261")
STRONG_ANCHOR = "HEX1"
REACTION_BLOCK_SIZE = 64
NON_EVALUABLE_STATUS = "NON_EVALUABLE_ZERO_POST_HOLDOUT_MASS under the frozen normalized-weight representation"
EVALUABLE_STATUS = "EVALUABLE"
PL1_FEATURE_FIELDS = (
    "sign_magnitude_eta2", "directional_entropy3", "dominant_sign_mass",
    "negative_width80_contraction_vs_all", "positive_width80_contraction_vs_all",
    "anchor_sd_contraction", "anchor_width80_contraction",
    "delta_anchor_target_correlation", "abs_delta_anchor_target_correlation",
    "joint_ess", "negative_ess", "tie_ess", "positive_ess",
    "target_delta_degenerate", "target_abs_magnitude_degenerate",
    "anchor_delta_degenerate", "sign_magnitude_status", "anchor_target_coupling_status",
)
PL1_JOIN_FIELDS = ("evaluation_id", "reaction_id")
PROXIMAL_SUBSYSTEMS = {
    "Glycolysis/gluconeogenesis", "Pyruvate metabolism",
    "Citric acid cycle", "Glutamate metabolism",
}

PINNED_SOURCES = {
    "candidate_table": (Path("12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/candidate_fraction_table.tsv.xz"), "421c2c1d196fbd088a4e92e2ecc998ad578aab814a025fff8e1c843597fed1da"),
    "bio0_manifest": (Path("outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_MANIFEST.json"), "40df8378876820bb616fdbdc8ef8a0f77d43e36f464104fca6cae9a16b39d7c1"),
    "bio0_contract": (Path("outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_ANALYSIS_CONTRACT.json"), "5e991aa58275e1545628bd9b3e9014b20c1fe799bfcc9810a0d6d712a28c8e2a"),
    "bio0_source_audit": (Path("outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_SOURCE_AUDIT.json"), "acf214cf3bbdb9a0737673e7808696da83f7b52cef599d43f660c21ba5b59787"),
    "bio0_support_audit": (Path("outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_SUPPORT_AUDIT.tsv"), "340c6ef821a093d4d7949b6488cdba2fa6cf9c3300b54b9a8b057016339de69e"),
    "reaction_panel": (Path("outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_REACTION_PANEL.tsv"), "f6e881ef949b2a519f0ac731d0dd8b08dba0d44f40e24ca94e18a47017907f4d"),
    "strong_weights": (Path("outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_FOUR_ARM_WEIGHTS.tsv.xz"), "843e17b7012a6fa212d1e9ac9e1dc6e309e854d94ab5c4b3f9e316b589768f4f"),
    "bio0r_manifest": (Path("outputs/dmi_bridge_bio0r_qualification_v1/BRIDGEBIO0R_MANIFEST.json"), "d155b63f59db3479b1dfd2d87c0b9ca470358c0166496abf288a649f04890cea"),
    "bio0r_identity_audit": (Path("outputs/dmi_bridge_bio0r_qualification_v1/BRIDGEBIO0R_IDENTITY_AUDIT.json"), "5b816176e2532ff076defd5fab103afa4210d9162990c85812209a7f6fa6f4a1"),
    "bio0r_source_audit": (Path("outputs/dmi_bridge_bio0r_qualification_v1/BRIDGEBIO0R_SOURCE_AUDIT.json"), "33ab67e8d8d12b84cdee1fc43b36bdba6b3680e9451bc01ed071e5fecd1c6fc4"),
    "flux_cache": (Path("outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstructed/flux_cache.npz"), "f9644da4ebd875cf5f0b34cc58b20faf0a0457d7c7a44635e8a010d8896c3dbb"),
    "pl1_context_summary": (PL1_REL / "BRIDGEPL1_CONTEXT_REACTION_SUMMARY.tsv.xz", "6f181d5354f07b5ea173f47453587e72734cdca455c72e36af106be5181784c3"),
    "pl1_reaction_summary": (PL1_REL / "BRIDGEPL1_REACTION_SUMMARY.tsv", "e4c9daed2923cff1b357a1e0f7e974485453f555317e7b5e9f68d2ff4303bf10"),
    "pl1_sensitivity_summary": (PL1_REL / "BRIDGEPL1_SENSITIVITY_SUMMARY.tsv", "98715eec90bbc150c94b37ea8ef3cff1fde553ee571b33b81667c433f097a78f"),
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = lzma.open if path.suffix == ".xz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field)) for field in fields})


def format_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, (float, np.floating)):
        return "nan" if not math.isfinite(float(value)) else format(float(value), ".17g")
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def stable_json_hash(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")).hexdigest()


def verify_pl1_artifacts(live_manifest: dict) -> dict[str, object]:
    artifact_hashes = live_manifest.get("artifact_sha256", {})
    required_names = [
        "BRIDGEPL1_EVALUATION_REGISTRY.tsv",
        "BRIDGEPL1_REACTION_EVALUATION_FEATURES.parts.json",
        *[f"BRIDGEPL1_REACTION_EVALUATION_FEATURES.part-{i:03d}.tsv.xz" for i in range(4)],
    ]
    for name in required_names:
        if name not in artifact_hashes:
            raise RuntimeError(f"PL1 live manifest lacks consumed artifact hash: {name}")
        path = ROOT / PL1_REL / name
        if not path.is_file() or sha256_file(path) != artifact_hashes[name]:
            raise RuntimeError(f"PL1 consumed artifact hash mismatch: {name}")
    feature_manifest = pl1_storage.validate_partitioned_tsv(ROOT / PL1_REL)
    feature_manifest_name = "BRIDGEPL1_REACTION_EVALUATION_FEATURES.parts.json"
    if sha256_file(ROOT / PL1_REL / feature_manifest_name) != artifact_hashes[feature_manifest_name]:
        raise RuntimeError("PL1 feature manifest hash mismatch")
    if feature_manifest["number_of_parts"] != 4 or feature_manifest["total_data_rows"] != 1672000:
        raise RuntimeError("PL1 feature partition accounting mismatch")
    manifest_parts = {part["filename"]: part for part in feature_manifest["parts"]}
    for name in required_names[2:]:
        record = manifest_parts.get(name)
        if record is None or record["sha256"] != artifact_hashes[name]:
            raise RuntimeError(f"PL1 feature-part manifest/hash mismatch: {name}")
    return {
        "evaluation_registry_sha256": artifact_hashes["BRIDGEPL1_EVALUATION_REGISTRY.tsv"],
        "feature_manifest_sha256": artifact_hashes[feature_manifest_name],
        "feature_part_sha256": {name: artifact_hashes[name] for name in required_names[2:]},
        "feature_part_records": feature_manifest["parts"],
        "feature_logical_content_sha256": feature_manifest["logical_content_sha256"],
        "feature_schema": feature_manifest["column_names"],
        "feature_rows": feature_manifest["total_data_rows"],
    }


def verify_sources() -> dict[str, str]:
    observed: dict[str, str] = {}
    patch_hash = sha256_file(ROOT / PATCH_REL)
    if patch_hash != EXPECTED_PATCH_SHA256:
        raise RuntimeError("foundation patch hash mismatch")
    for name, (relative, expected) in PINNED_SOURCES.items():
        path = ROOT / relative
        if not path.is_file():
            raise RuntimeError(f"missing pinned source: {relative}")
        actual = sha256_file(path)
        if actual != expected:
            raise RuntimeError(f"pinned source hash mismatch: {relative}")
        observed[name] = actual
    live_manifest_hash = sha256_file(ROOT / PL1_MANIFEST_REL)
    if live_manifest_hash != EXPECTED_LIVE_PL1_MANIFEST_SHA256:
        raise RuntimeError("live PL1 manifest changed after reconciliation audit")
    observed["pl1_manifest_live"] = live_manifest_hash
    observed["foundation_patch"] = patch_hash
    return observed


def load_and_validate_inputs() -> tuple[dict, list[dict], list[dict], dict, np.ndarray, dict]:
    manifest = json.loads((ROOT / PL1_MANIFEST_REL).read_text(encoding="utf-8"))
    source_hashes = verify_sources()
    pl1_artifact_audit = verify_pl1_artifacts(manifest)
    if manifest.get("status") != "PL1_PREDICTABILITY_LANDSCAPE_COMPLETE":
        raise RuntimeError("PL1 status is not complete")
    if manifest.get("pl2_status") != "NOT_RUN":
        raise RuntimeError("PL1 pl2_status is not NOT_RUN")
    if manifest.get("schema") != "bridge.pl1.predictability_landscape.v1":
        raise RuntimeError("PL1 schema mismatch")
    counts = manifest.get("row_counts", {})
    expected_counts = {"evaluation_registry": 400, "reaction_evaluation_features": 1672000, "context_reaction_summary": 66880, "reaction_summary": 4180}
    if any(counts.get(key) != value for key, value in expected_counts.items()):
        raise RuntimeError("PL1 row-count contract mismatch")
    candidates = read_tsv(ROOT / PINNED_SOURCES["candidate_table"][0])
    weights = read_tsv(ROOT / PINNED_SOURCES["strong_weights"][0])
    panel_rows = read_tsv(ROOT / PINNED_SOURCES["reaction_panel"][0])
    if len(candidates) != 640 or len({tuple(r[k] for k in ("algorithm", "tumor", "ensemble_hash", "sample_index")) for r in candidates}) != 640:
        raise RuntimeError("candidate panel is not exactly 640 unique identities")
    strata = defaultdict(list)
    for row in candidates:
        strata[(row["algorithm"], row["tumor"], row["ensemble_hash"])].append(row)
    if len(strata) != 32 or any(len(rows) != 20 for rows in strata.values()):
        raise RuntimeError("candidate panel is not exactly 32 strata x 20 candidates")
    if len(weights) != 12800:
        raise RuntimeError("BIO-0 four-arm weight row count is not 12800")
    strong = [row for row in weights if row["arm"] == "strong_anchor_baseline"]
    if len(strong) != 3200:
        raise RuntimeError("strong-anchor row count is not 3200")
    by_weight_key = defaultdict(list)
    for row in strong:
        key = (row["mouse_id"], row["tumor"], row["algorithm"], row["ensemble_hash"], row["rna_context_key"])
        weight = float(row["weight"])
        if not math.isfinite(weight) or weight < 0:
            raise RuntimeError("strong-anchor weights are not finite non-negative")
        by_weight_key[key].append((int(row["sample_index"]), weight))
    if len(by_weight_key) != 160 or any(len(v) != 20 for v in by_weight_key.values()):
        raise RuntimeError("strong-anchor support is not exactly 160 x 20")
    cache = np.load(ROOT / PINNED_SOURCES["flux_cache"][0], allow_pickle=False)
    rxns = [str(x) for x in cache["rxns"].tolist()]
    if len(rxns) != 4181 or len(set(rxns)) != 4181 or rxns.count(STRONG_ANCHOR) != 1:
        raise RuntimeError("cache reaction identity is invalid")
    if cache["mats"].shape != (32, 20, 4181) or not np.isfinite(cache["mats"]).all():
        raise RuntimeError("cache shape or finiteness mismatch")
    cache_keys = {(str(a), str(t), str(e)) for a, t, e in zip(cache["alg"], cache["tumor"], cache["eh"])}
    if cache_keys != set(strata):
        raise RuntimeError("candidate/cache strata do not match")
    return manifest, candidates, strong, {"panel": panel_rows, "weights": by_weight_key, "pl1_artifact_audit": pl1_artifact_audit}, cache["mats"], {"rxns": rxns, "hashes": source_hashes, "pl1_artifact_audit": pl1_artifact_audit, "alg": [str(x) for x in cache["alg"].tolist()], "tumor": [str(x) for x in cache["tumor"].tolist()], "eh": [str(x) for x in cache["eh"].tolist()]}


def candidate_identity(row: dict[str, str]) -> str:
    fields = ("algorithm", "tumor", "ensemble_hash", "projection_hash", "rna_context_key", "rna_training_samples", "sample_index")
    return stable_json_hash({"schema": "bridge.pl2a.candidate_identity.v1", "identity": {field: str(row[field]) for field in fields}})


def select_truth(candidates: list[dict[str, str]], strong: list[dict[str, str]], by_data: dict, mats: np.ndarray, rxns: list[str], cache_info: dict) -> tuple[list[dict[str, object]], dict[tuple[str, str, str, str], dict[str, dict[str, object]]]]:
    candidate_map = {(r["algorithm"], r["tumor"], r["ensemble_hash"], int(r["sample_index"])): r for r in candidates}
    cache_row = {(algorithm, tumor, ensemble_hash): index for index, (algorithm, tumor, ensemble_hash) in enumerate(zip(cache_info["alg"], cache_info["tumor"], cache_info["eh"]))}
    if len(cache_row) != 32:
        raise RuntimeError("cache row map construction failed")
    selected: list[dict[str, object]] = []
    alias_map: dict[tuple[str, str, str, str], dict[str, dict[str, object]]] = defaultdict(dict)
    for key in sorted(by_data["weights"]):
        mouse, tumor, algorithm, ensemble_hash, rna = key
        rows = [candidate_map[(algorithm, tumor, ensemble_hash, sample)] for sample, _ in by_data["weights"][key]]
        raw_by_sample = dict(by_data["weights"][key])
        cache_values = mats[cache_row[(algorithm, tumor, ensemble_hash)], :, rxns.index(STRONG_ANCHOR)]
        rows.sort(key=lambda row: (float(cache_values[int(row["sample_index"])]), candidate_identity(row)))
        weights = np.asarray([raw_by_sample[int(row["sample_index"])] for row in rows], dtype=float)
        normalized = weights / weights.sum()
        cumulative = np.cumsum(normalized)
        chosen_by_candidate: dict[str, list[str]] = defaultdict(list)
        chosen_rows: dict[str, dict[str, object]] = {}
        for quantile in QUANTILES:
            index = min(int(np.searchsorted(cumulative, quantile, side="left")), len(rows) - 1)
            row = rows[index]
            cid = candidate_identity(row)
            remaining = float(1.0 - normalized[index])
            positive_before = int(np.count_nonzero(normalized > 0.0))
            without = np.delete(normalized, index)
            positive_after = int(np.count_nonzero(without > 0.0))
            status = EVALUABLE_STATUS if remaining > 0.0 else NON_EVALUABLE_STATUS
            post_ess = float("nan") if remaining <= 0.0 else float(1.0 / np.sum((without / remaining) ** 2))
            record = {
                "mouse_id": mouse, "tumor": tumor, "algorithm": algorithm,
                "ensemble_hash": ensemble_hash, "projection_hash": row["projection_hash"],
                "rna_context_key": rna, "rna_training_samples": row["rna_training_samples"],
                "requested_quantile": f"{quantile:.2f}", "candidate_id": cid,
                "sample_index": int(row["sample_index"]),
                "HEX1": float(cache_values[int(row["sample_index"])]),
                "pre_holdout_normalized_weight": float(normalized[index]),
                "pre_holdout_mass": 1.0, "post_holdout_mass": remaining,
                "positive_weight_candidates_pre": positive_before,
                "positive_weight_candidates_post": positive_after,
                "pre_holdout_ess": float(1.0 / np.sum(normalized ** 2)),
                "post_holdout_ess": post_ess,
                "sole_positive_weight_candidate": positive_before == 1,
                "holdout_support_status": status,
            }
            selected.append(record)
            chosen_by_candidate[cid].append(f"{quantile:.2f}")
            chosen_rows[cid] = record
        for cid, aliases in chosen_by_candidate.items():
            for alias in aliases:
                alias_map[(mouse, tumor, algorithm, rna)][alias] = {**chosen_rows[cid], "collapsed_quantile_aliases": ",".join(aliases)}
    if len(selected) != 480:
        raise RuntimeError("truth selection count is not 480")
    if sum(row["holdout_support_status"] == NON_EVALUABLE_STATUS for row in selected) != 18:
        raise RuntimeError("zero post-holdout truth selection count is not 18")
    return selected, alias_map


def build_evaluation_registry() -> list[dict[str, str]]:
    return read_tsv(ROOT / PL1_REL / "BRIDGEPL1_EVALUATION_REGISTRY.tsv")


def build_pair_registry(evaluations: list[dict[str, str]], aliases: dict) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for evaluation in evaluations:
        for quantile in ("0.10", "0.50", "0.90"):
            left = aliases[(evaluation["ct2a_mouse"], "CT2A", evaluation["algorithm"], evaluation["rna_context_key"])][quantile]
            right = aliases[(evaluation["gl261_mouse"], "GL261", evaluation["algorithm"], evaluation["rna_context_key"])][quantile]
            left_valid = left["post_holdout_mass"] > 0.0
            right_valid = right["post_holdout_mass"] > 0.0
            status = "VALID" if left_valid and right_valid else NON_EVALUABLE_STATUS
            rows.append({
                "evaluation_id": evaluation["evaluation_id"], "algorithm": evaluation["algorithm"],
                "rna_context_key": evaluation["rna_context_key"], "ct2a_mouse": evaluation["ct2a_mouse"],
                "gl261_mouse": evaluation["gl261_mouse"], "quantile_alias": quantile,
                "ct2a_candidate_id": left["candidate_id"], "gl261_candidate_id": right["candidate_id"],
                "ct2a_sample_index": left["sample_index"], "gl261_sample_index": right["sample_index"],
                "ct2a_post_holdout_mass": left["post_holdout_mass"], "gl261_post_holdout_mass": right["post_holdout_mass"],
                "ct2a_post_holdout_ess": left["post_holdout_ess"], "gl261_post_holdout_ess": right["post_holdout_ess"],
                "ct2a_holdout_status": left["holdout_support_status"], "gl261_holdout_status": right["holdout_support_status"],
                "pair_evaluable": left_valid and right_valid, "pair_status": status,
                "non_evaluable_reason": "" if left_valid and right_valid else "zero post-holdout mass on one or both condition distributions under frozen normalized weights",
            })
    if len(rows) != 1200:
        raise RuntimeError("truth-pair alias count is not 1200")
    if sum(bool(row["pair_evaluable"]) for row in rows) != 1110:
        raise RuntimeError("evaluable truth-pair alias count is not 1110")
    return rows


def audit_feature_join(evaluations: list[dict[str, str]], rxns: list[str]) -> dict[str, object]:
    expected = {(row["evaluation_id"], reaction) for row in evaluations for reaction in rxns if reaction != STRONG_ANCHOR}
    seen: set[tuple[str, str]] = set()
    feature_fields: list[str] | None = None
    part_dir = ROOT / PL1_REL
    parts = sorted(part_dir.glob("BRIDGEPL1_REACTION_EVALUATION_FEATURES.part-*.tsv.xz"))
    for part in parts:
        with lzma.open(part, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            feature_fields = reader.fieldnames
            for row in reader:
                key = (row["evaluation_id"], row["reaction_id"])
                if key in seen:
                    raise RuntimeError("duplicate PL1 predictor join key")
                seen.add(key)
    missing = expected - seen
    extra = seen - expected
    if missing or extra or len(seen) != 1672000 or not feature_fields:
        raise RuntimeError("PL1 predictor join coverage mismatch")
    missing_fields = [field for field in PL1_FEATURE_FIELDS if field not in feature_fields]
    if missing_fields:
        raise RuntimeError(f"missing PL1 predictor fields: {missing_fields}")
    return {"join_fields": list(PL1_JOIN_FIELDS), "coverage_rows": len(seen), "expected_rows": len(expected), "coverage_fraction": 1.0, "feature_fields": list(PL1_FEATURE_FIELDS), "missing_fields": [], "part_count": len(parts)}


def build_reaction_registry(panel_rows: list[dict[str, str]], rxns: list[str]) -> list[dict[str, object]]:
    panel = {row["reaction_id"]: row["subsystem"] for row in panel_rows}
    reactions = sorted(set(rxns) - {STRONG_ANCHOR})
    if len(reactions) != 4180:
        raise RuntimeError("reaction registry is not exactly 4180 reactions")
    anchor_subsystem = panel.get(STRONG_ANCHOR, "UNMAPPED_PATHWAY")
    return [{"reaction_id": reaction, "subsystem": panel.get(reaction, "UNMAPPED_PATHWAY"), "same_subsystem_as_HEX1": panel.get(reaction, "UNMAPPED_PATHWAY") == anchor_subsystem, "proximal_set_member": panel.get(reaction, "UNMAPPED_PATHWAY") in PROXIMAL_SUBSYSTEMS} for reaction in reactions]


def json_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True) + "\n").encode("utf-8")


def prepare(output: Path) -> dict[str, object]:
    if output.exists():
        manifest_path = output / "BRIDGEPL2A_MANIFEST.json"
        if not manifest_path.is_file():
            raise RuntimeError("refusing to overwrite existing PL2A directory without manifest")
        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
        if existing.get("status") != "PL2A_READY_FOR_PL2B_SIGN_ONLY_UTILITY":
            raise RuntimeError("refusing to overwrite differing PL2A output")
        for name, expected in existing.get("artifact_sha256", {}).items():
            path = output / name
            if not path.is_file() or sha256_file(path) != expected:
                raise RuntimeError("refusing to accept differing PL2A output artifact")
        if existing.get("counts", {}).get("truth_pair_aliases_total") != 1200 or existing.get("counts", {}).get("truth_pair_aliases_evaluable") != 1110 or existing.get("counts", {}).get("truth_pair_aliases_non_evaluable") != 90:
            raise RuntimeError("refusing to accept differing PL2A output counts")
        live_manifest = json.loads((ROOT / PL1_MANIFEST_REL).read_text(encoding="utf-8"))
        verify_sources()
        verify_pl1_artifacts(live_manifest)
        return {"status": "NO_OP_EXISTING_IDENTICAL_CANDIDATE", "manifest": existing}
    manifest, candidates, strong, data, mats, cache_info = load_and_validate_inputs()
    evaluations = build_evaluation_registry()
    selected, aliases = select_truth(candidates, strong, data, mats, cache_info["rxns"], cache_info)
    pairs = build_pair_registry(evaluations, aliases)
    reaction_registry = build_reaction_registry(data["panel"], cache_info["rxns"])
    join_audit = audit_feature_join(evaluations, cache_info["rxns"])
    source_hashes = cache_info["hashes"]
    source_hashes["pl1_manifest_historical_contract"] = EXPECTED_HISTORICAL_PL1_MANIFEST_SHA256
    summary = {
        "truth_selection_rows": len(selected), "truth_selection_zero_post_holdout_mass_rows": 18,
        "truth_pair_aliases_total": len(pairs), "truth_pair_aliases_evaluable": sum(bool(row["pair_evaluable"]) for row in pairs),
        "truth_pair_aliases_non_evaluable": sum(not bool(row["pair_evaluable"]) for row in pairs),
    }
    artifact_names = ["BRIDGEPL2A_SOURCE_AUDIT.json", "BRIDGEPL2A_ANALYSIS_CONTRACT.json", "BRIDGEPL2A_REACTION_REGISTRY.tsv", "BRIDGEPL2A_TRUTH_SELECTION.tsv", "BRIDGEPL2A_TRUTH_PAIR_REGISTRY.tsv", "BRIDGEPL2A_PL1_PREDICTOR_CONTRACT.json", "BRIDGEPL2A_EXECUTION_BUDGET.json"]
    with tempfile.TemporaryDirectory(prefix="dmi_bridge_pl2a_", dir=str(output.parent)) as temp_name:
        temp = Path(temp_name)
        source_audit = {
            "schema": "bridge.pl2a.source_audit.v1", "status": "PL2A_SOURCE_AUDIT_COMPLETE",
            "foundation_patch_sha256": EXPECTED_PATCH_SHA256,
            "pl1_manifest_historical_sha256": EXPECTED_HISTORICAL_PL1_MANIFEST_SHA256,
            "pl1_manifest_live_sha256": source_hashes["pl1_manifest_live"],
            "manifest_provenance_treatment": "live manifest is recorded as a provenance/storage successor; historical hash is not silently repinned",
            "source_sha256": source_hashes, "pl1_counts": manifest["row_counts"], "summary": summary,
            "pl1_consumed_artifacts": cache_info["pl1_artifact_audit"],
        }
        analysis_contract = {
            "schema": "bridge.pl2a.sign_only_prepare.v1", "stage_role": "OUTCOME_BLIND_PL2_PREPARATION",
            "strong_anchor_reaction": STRONG_ANCHOR, "candidate_scope": "all 4180 finite cache reactions except HEX1",
            "truth_selection": {"basis": "HEX1 coordinate and strong_anchor_baseline weights only", "quantiles": list(QUANTILES), "candidate_identity_fields": ["algorithm", "tumor", "ensemble_hash", "projection_hash", "rna_context_key", "rna_training_samples", "sample_index"], "collapse_duplicate_selected_candidates": True, "matched_aliases_only": True, "holdout": "remove selected candidate from its own condition inference pool", "zero_mass_status": NON_EVALUABLE_STATUS},
            "lambda_total": {"primary": 0.25, "sensitivity": [0.5, 1.0], "lambda_zero_qc": 0.0},
            "reliability_grid": [0.0, 0.25, 0.5, 0.75, 1.0],
            "primary_endpoint_deferred_to_pl2b": "absolute_error_of_E_abs_delta_B", "historical_bridge_outcomes_used": False,
            "forbidden_in_pl2a": ["reaction_B_truth_deltas", "correct_wrong_gains", "random_sign_gains", "reliability_curves", "solver", "optimization", "FVA", "sampling", "reconstruction", "new_strong_anchor_weights"],
        }
        predictor_contract = {"schema": "bridge.pl2a.pl1_predictor_contract.v1", "join_fields": list(PL1_JOIN_FIELDS), "coverage": join_audit, "primary_predictor": "sign_magnitude_eta2", "secondary_predictors": list(PL1_FEATURE_FIELDS[1:]), "values_frozen_before_truth_holdout": True, "degenerate_and_nonfinite_values": "preserved explicitly; no PL2-derived bins or composite score", "pl1_consumed_artifacts": cache_info["pl1_artifact_audit"]}
        blocks = int(math.ceil(4180 / REACTION_BLOCK_SIZE))
        evaluable = summary["truth_pair_aliases_evaluable"]
        execution_budget = {"schema": "bridge.pl2a.execution_budget.v1", "reaction_count": 4180, "truth_pair_aliases_total": 1200, "truth_pair_aliases_evaluable": evaluable, "truth_pair_aliases_non_evaluable": 90, "lambda_arms": [0.0, 0.25, 0.5, 1.0], "reliability_values": 5, "potential_cases_per_lambda_all_aliases": 4180 * 1200, "evaluable_cases_per_lambda": 4180 * evaluable, "potential_cases_all_lambda_arms_all_aliases": 4180 * 1200 * 4, "evaluable_cases_all_lambda_arms": 4180 * evaluable * 4, "evaluable_reliability_rows": 4180 * evaluable * 4 * 5, "reaction_block_size": REACTION_BLOCK_SIZE, "reaction_blocks": blocks, "resumable_shard_key": ["algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse", "quantile_alias", "reaction_block"], "resumable_shard_count_upper_bound": 16 * 25 * 3 * blocks}
        write_json = lambda name, value: (temp / name).write_bytes(json_bytes(value))
        write_json("BRIDGEPL2A_SOURCE_AUDIT.json", source_audit)
        write_json("BRIDGEPL2A_ANALYSIS_CONTRACT.json", analysis_contract)
        write_tsv(temp / "BRIDGEPL2A_REACTION_REGISTRY.tsv", ["reaction_id", "subsystem", "same_subsystem_as_HEX1", "proximal_set_member"], reaction_registry)
        write_tsv(temp / "BRIDGEPL2A_TRUTH_SELECTION.tsv", ["mouse_id", "tumor", "algorithm", "ensemble_hash", "projection_hash", "rna_context_key", "rna_training_samples", "requested_quantile", "collapsed_quantile_aliases", "candidate_id", "sample_index", "HEX1", "pre_holdout_normalized_weight", "pre_holdout_mass", "post_holdout_mass", "positive_weight_candidates_pre", "positive_weight_candidates_post", "pre_holdout_ess", "post_holdout_ess", "sole_positive_weight_candidate", "holdout_support_status"], [{**row, "collapsed_quantile_aliases": next((r["collapsed_quantile_aliases"] for r in aliases[(row["mouse_id"], row["tumor"], row["algorithm"], row["rna_context_key"])].values() if r["candidate_id"] == row["candidate_id"]), "")} for row in selected])
        write_tsv(temp / "BRIDGEPL2A_TRUTH_PAIR_REGISTRY.tsv", ["evaluation_id", "algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse", "quantile_alias", "ct2a_candidate_id", "gl261_candidate_id", "ct2a_sample_index", "gl261_sample_index", "ct2a_post_holdout_mass", "gl261_post_holdout_mass", "ct2a_post_holdout_ess", "gl261_post_holdout_ess", "ct2a_holdout_status", "gl261_holdout_status", "pair_evaluable", "pair_status", "non_evaluable_reason"], pairs)
        write_json("BRIDGEPL2A_PL1_PREDICTOR_CONTRACT.json", predictor_contract)
        write_json("BRIDGEPL2A_EXECUTION_BUDGET.json", execution_budget)
        artifact_hashes = {name: sha256_file(temp / name) for name in artifact_names}
        final_manifest = {"schema": "bridge.pl2a.prepare.manifest.v1", "status": "PL2A_READY_FOR_PL2B_SIGN_ONLY_UTILITY", "pl2_status": "NOT_RUN", "pl1_manifest_historical_sha256": EXPECTED_HISTORICAL_PL1_MANIFEST_SHA256, "pl1_manifest_live_sha256": source_hashes["pl1_manifest_live"], "pl1_manifest_identity_treatment": "historical contract identity preserved; live manifest recorded as provenance/storage successor", "counts": {**summary, "reaction_registry": len(reaction_registry), "pl1_predictor_join_rows": join_audit["coverage_rows"]}, "artifact_sha256": artifact_hashes, "pl2_scientific_outcomes_computed": False, "reaction_B_truth_deltas_computed": False, "correct_wrong_gains_computed": False, "historical_bridge_outcomes_used": False, "solver_invoked": False, "optimization_invoked": False, "fva_invoked": False, "sampling_invoked": False, "reconstruction_invoked": False, "new_strong_anchor_weights_generated": False}
        (temp / "BRIDGEPL2A_MANIFEST.json").write_bytes(json_bytes(final_manifest))
        output.parent.mkdir(parents=True, exist_ok=True)
        os.replace(temp, output)
    return {"status": final_manifest["status"], "manifest": final_manifest}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / OUTPUT_REL)
    args = parser.parse_args()
    result = prepare(args.output)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
