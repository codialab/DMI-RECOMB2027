#!/usr/bin/env python3
"""Terminal, read-only-source matched A1/A2 BRIDGE synthesis."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import lzma
import math
import platform
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

from scripts import dmi_bridge_pl2c_geometry_utility_core_v1 as pl2c
from scripts import dmi_bridge_pl2d_confirmation_core_v1 as dcore
from scripts import dmi_bridge_pl2d_hypothesis_freeze_core_v1 as hcore

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/dmi_bridge_a23_a1_a2_synthesis_v1"
DIRS = {
    "A20": ROOT / "outputs/dmi_bridge_a20_dual_anchor_qualification_v1",
    "A21": ROOT / "outputs/dmi_bridge_a21_dual_anchor_geometry_v1",
    "A22": ROOT / "outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1",
    "PL1": ROOT / "outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1",
    "PL2A": ROOT / "outputs/dmi_bridge_pl2a_sign_only_prepare_v1",
    "PL2B": ROOT / "outputs/dmi_bridge_pl2b_sign_only_utility_v1",
    "PL2C": ROOT / "outputs/dmi_bridge_pl2c_geometry_utility_development_v1",
    "PL2DH": ROOT / "outputs/dmi_bridge_pl2d_hypothesis_freeze_v1",
    "PL2DHA": ROOT / "outputs/dmi_bridge_pl2d_h_support_gate_amendment_v1",
    "PL2D": ROOT / "outputs/dmi_bridge_pl2d_confirmation_v1",
    "M1": ROOT / "outputs/dmi_bridge_m1_manuscript_evidence_v1",
}
MANIFESTS = {
    "A20": ("BRIDGEA20_MANIFEST.json", "01f8004a0bf9d5358d7ce7281e827cc1a1fbe00a6750a5c063a1edad2eedf380", "BRIDGE_A20_PARTIAL_QUALIFICATION"),
    "A21": ("BRIDGEA21_MANIFEST.json", "db9948d7b99447a695cf8aa21429dde6d642e1fd0632bcce52643deac59e6feb", "BRIDGE_A21_DUAL_ANCHOR_GEOMETRY_FROZEN"),
    "A22": ("BRIDGEA22_MANIFEST.json", "15b237151005538e5077a352868f873ef00466f4f5ea49307b781a8a8cd0e7cb", "BRIDGE_A22_SIGN_ONLY_UTILITY_FROZEN"),
    "PL1": ("BRIDGEPL1_MANIFEST.json", "b0f81620ff24f5791dd964a8063416784db86a358b0464b8c68a504e84e5f529", "PL1_PREDICTABILITY_LANDSCAPE_COMPLETE"),
    "PL2A": ("BRIDGEPL2A_MANIFEST.json", "8ae3eb60eb23614e379e38e3aa3cb72ec16c6bdb69aa9c74cc94499e31f1c312", "PL2A_READY_FOR_PL2B_SIGN_ONLY_UTILITY"),
    "PL2B": ("BRIDGEPL2B_MANIFEST.json", "7d9ca97a3f4ea5aa8b0b983e09ee83fc3c82f326d1b1cbf9c6137cfadf25539d", "PL2B_SIGN_ONLY_UTILITY_COMPLETE"),
    "PL2C": ("BRIDGEPL2C_MANIFEST.json", "6abcf149cd6adb294c0931585cec28f992b453451a5b5f03e10b3ffb2f9b9c39", "PL2C_DEVELOPMENT_GEOMETRY_UTILITY_COMPLETE"),
    "PL2DH": ("BRIDGEPL2DH_MANIFEST.json", "a1dd76d6b40e741e07efe64c34bc39eaad289a46bf6c870917a4f59cab688427", "PL2D_HYPOTHESIS_FROZEN_READY_FOR_CONFIRMATION"),
    "PL2DHA": ("BRIDGEPL2DHA_MANIFEST.json", "fb7df9e208da2e6c9841d43a3b572038e59943a4221535705eef9f5dc0ed6a3a", "PL2D_H_SUPPORT_GATE_AMENDED_READY_FOR_CONFIRMATION"),
    "PL2D": ("BRIDGEPL2D_MANIFEST.json", "9f9eb730473e7ad209c1d5ec3340808ab614087b8ae1121d929eff61be12b699", "PL2D_CONFIRMED"),
    "M1": ("BRIDGEM1_MANIFEST.json", "201ba9f4fc2f74290d991cd294aec558c5e5dd23282877e871a5c0669a44cdde", "BRIDGE_M1_MANUSCRIPT_EVIDENCE_FROZEN"),
}
LOGICAL = {"A21": "50bc3e47239e4b22d0931246f4fa0d5d126a5b69226d47cb342cb8a5b3ee99d2", "A22": "70321eb7d322f06a3c1e29c048d08f9b4c55d291763b5194f6d626325ed786f4", "PL1": "46aa9f0410f89d03b020e01fefa517348faae1d55d9bfe8387ac04b39bb9ac93", "PL2B": "c3670992a0ad33698af6c7620429bed48e0ecaa98adb3c87e648eefb1ba70363"}
FILES = (
    "BRIDGEA23_SOURCE_AUDIT.json", "BRIDGEA23_ANALYSIS_CONTRACT.json",
    "BRIDGEA23_MATCHED_POPULATION_AUDIT.tsv", "BRIDGEA23_GEOMETRY_COMPARISON.tsv",
    "BRIDGEA23_UTILITY_COMPARISON.tsv", "BRIDGEA23_GEOMETRY_UTILITY_PERSISTENCE.tsv",
    "BRIDGEA23_CONTROL_SUMMARY.tsv", "BRIDGEA23_SENSITIVITY_SUMMARY.tsv", "BRIDGEA23_QC.json",
)
FEATURES = pl2c.PRIMARY_FEATURES
EXPECTED_ROWS = {
    "BRIDGEA23_MATCHED_POPULATION_AUDIT.tsv": 836,
    "BRIDGEA23_GEOMETRY_COMPARISON.tsv": 24,
    "BRIDGEA23_UTILITY_COMPARISON.tsv": 3,
    "BRIDGEA23_GEOMETRY_UTILITY_PERSISTENCE.tsv": 42,
    "BRIDGEA23_CONTROL_SUMMARY.tsv": 9,
    "BRIDGEA23_SENSITIVITY_SUMMARY.tsv": 16,
}


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def safe_file(path: Path) -> Path:
    resolved = path.resolve(strict=True)
    require(not path.is_symlink() and resolved.is_relative_to(ROOT) and resolved.is_file(), f"unsafe source: {path}")
    return resolved


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with safe_file(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def checked(path: Path, expected: str) -> None:
    require(sha(path) == expected, f"SHA256 mismatch: {path.relative_to(ROOT)}")


def read_json(path: Path) -> dict:
    return json.loads(safe_file(path).read_text(encoding="utf-8"))


def rows(path: Path):
    with (lzma.open(safe_file(path), "rt", encoding="utf-8", newline="") if path.suffix == ".xz" else safe_file(path).open("r", encoding="utf-8", newline="")) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_json(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def fmt(value):
    if value is None:
        return ""
    if isinstance(value, (float, np.floating)):
        return format(float(value), ".17g") if math.isfinite(float(value)) else ""
    return str(value)


def write_tsv(path: Path, records: list[dict]) -> None:
    require(bool(records), f"empty output: {path.name}")
    fields = tuple(records[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        for row in records:
            writer.writerow({field: fmt(row.get(field)) for field in fields})


def parts(stage: str, prefix: str) -> list[Path]:
    info = read_json(DIRS[stage] / f"{prefix}.parts.json")
    require(info["number_of_parts"] == len(info["parts"]) and info["total_data_rows"] in (3493644, 3494480, 1671600, 1672000), f"{stage} part count")
    out = []
    for record in info["parts"]:
        name = record["filename"]
        require(Path(name).name == name, "unsafe part name")
        path = DIRS[stage] / name
        checked(path, record["sha256"])
        out.append(path)
    return out


def admit() -> tuple[dict, dict]:
    manifests = {}
    hashes = {}
    for stage, (name, digest, status) in MANIFESTS.items():
        path = DIRS[stage] / name
        checked(path, digest)
        manifest = read_json(path)
        require(manifest["status"] == status, f"{stage} status changed")
        manifests[stage], hashes[str(path.relative_to(ROOT))] = manifest, digest
        for artifact, expected in manifest.get("artifact_sha256", {}).items():
            require(Path(artifact).name == artifact, "unsafe artifact name")
            checked(DIRS[stage] / artifact, expected)
    require(manifests["A20"]["qualification"] == {"A2-L": "QUALIFIED", "A2-G": "NOT_QUALIFIED"}, "A2 qualification changed")
    require(manifests["A21"]["feature_logical_content_sha256"] == LOGICAL["A21"] and manifests["A21"]["a2_g_used"] is False, "A2.1 logical identity changed")
    require(manifests["A22"]["case_logical_content_sha256"] == LOGICAL["A22"] and manifests["A22"]["a2_g_used"] is False and manifests["A22"]["a1_vs_a2_synthesis_computed"] is False, "A2.2 logical identity changed")
    require(manifests["PL2B"]["case_logical_content_sha256"] == LOGICAL["PL2B"], "A1 case identity changed")
    require(manifests["PL2B"]["fingerprint"]["admission"]["pl1_consumed_artifacts"]["feature_logical_content_sha256"] == LOGICAL["PL1"], "A1 PL1 feature identity changed")
    require(manifests["PL2C"]["split_registry_sha256"] == "7dd8676508416a13f88af184110ef40e0f8d92782cd5d62b42e7006af5316e2f", "split changed")
    m1_sources = manifests["M1"]["source_sha256"]
    for stage in ("PL2B", "PL2C", "PL2DH", "PL2DHA", "PL2D"):
        p = DIRS[stage] / MANIFESTS[stage][0]
        require(m1_sources[str(p.relative_to(ROOT))] == MANIFESTS[stage][1], f"M1 did not bind {stage}")
    for stage, prefix in (("A22", "BRIDGEA22_CASE_OUTCOMES"), ("PL2B", "BRIDGEPL2B_CASE_OUTCOMES"), ("A21", "BRIDGEA21_REACTION_EVALUATION_FEATURES"), ("PL1", "BRIDGEPL1_REACTION_EVALUATION_FEATURES")):
        # The manifest hashes above bind the parts manifests and every part.
        if stage in ("A22", "PL2B"):
            parts(stage, prefix)
    contracts = {stage: read_json(DIRS[stage] / name) for stage, name in (("PL2B", "BRIDGEPL2B_ANALYSIS_CONTRACT.json"), ("PL2C", "BRIDGEPL2C_ANALYSIS_CONTRACT.json"), ("PL2DH", "BRIDGEPL2DH_CONFIRMATION_CONTRACT.json"), ("PL2DHA", "BRIDGEPL2DHA_AMENDED_CONFIRMATION_CONTRACT.json"), ("A21", "BRIDGEA21_ANALYSIS_CONTRACT.json"), ("A22", "BRIDGEA22_ANALYSIS_CONTRACT.json"))}
    require(contracts["PL2B"]["primary_lambda"] == contracts["PL2C"]["primary_lambda"] == contracts["A22"]["primary_lambda"] == .25, "lambda changed")
    require(contracts["PL2B"]["gain_classification"]["tolerance"] == contracts["PL2C"]["gain_tolerance"] == contracts["A22"]["gain_tol"] == 1e-12, "gain tolerance changed")
    require(contracts["PL2C"]["primary_features"] == list(FEATURES) == contracts["A21"]["a1_comparable_features"], "features changed")
    require(contracts["PL2DHA"]["primary"]["predictor"] == "mean_sign_magnitude_eta2_across_evaluations" and contracts["PL2DHA"]["primary"]["minimum_finite_reactions_per_context"] == 30, "effective PL2D contract changed")
    require(contracts["A21"]["common_anchor_coupling_semantics"] == "HEX1_contrast_preserved_exactly_for_A1_A2_comparability", "HEX1 coupling changed")
    qc = read_json(DIRS["A22"] / "BRIDGEA22_QC.json")
    require(qc["status"] == "PASS" and qc["case_logical_content_sha256"] == LOGICAL["A22"], "A2.2 QC failed")
    return manifests, {"schema": "bridge.a23.source_audit.v1", "manifest_sha256": hashes, "logical_sha256": LOGICAL, "m1_source_audit_sha256": sha(DIRS["M1"] / "BRIDGEM1_SOURCE_AUDIT.json"), "foundation_patch_sha256": sha(ROOT / "dmi_bridge_a23_a1_a2_synthesis_foundation_v1.patch"), "status": "PASS"}


def identity(manifests: dict):
    a2_pairs = list(rows(DIRS["A22"] / "BRIDGEA22_MATCHED_TRUTH_PAIR_REGISTRY.tsv"))
    a1_pairs = [r for r in rows(DIRS["PL2B"] / "BRIDGEPL2B_DISTINCT_TRUTH_PAIR_REGISTRY.tsv") if r["pair_evaluable"] == "true"]
    require(len(a2_pairs) == len(a1_pairs) == 836 and a2_pairs == a1_pairs, "truth-pair registry changed")
    pair_by_id = {r["truth_pair_id"]: r for r in a2_pairs}
    require(len(pair_by_id) == 836, "duplicate truth pair")
    reactions = list(rows(DIRS["A22"] / "BRIDGEA22_REACTION_REGISTRY.tsv"))
    rid = [r["reaction_id"] for r in reactions]
    split = list(rows(DIRS["PL2C"] / "BRIDGEPL2C_REACTION_SPLIT.tsv"))
    split_by_id = {r["reaction_id"]: r["analysis_split"] for r in split}
    require(len(rid) == len(set(rid)) == 4179 and len(split_by_id) == 4180 and set(rid) == set(split_by_id) - {"LDH_L"} and split_by_id["LDH_L"] == "DEVELOPMENT", "reaction population changed")
    require(Counter(split_by_id[r] for r in rid) == {"DEVELOPMENT": 3134, "CONFIRMATION_HOLDOUT": 1045}, "matched split changed")
    a1_eval = list(rows(DIRS["PL1"] / "BRIDGEPL1_EVALUATION_REGISTRY.tsv"))
    a2_eval = list(rows(DIRS["A21"] / "BRIDGEA21_EVALUATION_REGISTRY.tsv"))
    common = ("evaluation_id", "algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse", "cartesian_support", "joint_weight_sum", "status")
    require(len(a1_eval) == len(a2_eval) == 400 and [{k:r[k] for k in common} for r in a1_eval] == [{k:r[k] for k in common} for r in a2_eval], "evaluation identity changed")
    evals = [r["evaluation_id"] for r in a1_eval]
    require(len(set(evals)) == 400, "duplicate evaluation")
    eval_by_id = {r["evaluation_id"]: r for r in a1_eval}
    require({r["evaluation_id"] for r in a2_pairs} <= set(evals) and len({r["evaluation_id"] for r in a2_pairs}) == 370, "truth evaluation population changed")
    return pair_by_id, reactions, split_by_id, evals, eval_by_id


def iter_parts(stage: str, prefix: str):
    for path in parts(stage, prefix):
        yield from rows(path)


def case_stream(stage: str, prefix: str, pair_by_id: dict, rid_set: set, expected: int, audit: dict):
    previous = None
    count = 0
    skipped = 0
    by_pair = Counter()
    for row in iter_parts(stage, prefix):
        key = (row["truth_pair_id"], row["reaction_id"])
        require(previous is None or key > previous, f"{stage} case key order or duplicate")
        previous = key
        if stage == "PL2B" and row["reaction_id"] == "LDH_L":
            skipped += 1
            continue
        require(row["truth_pair_id"] in pair_by_id and row["reaction_id"] in rid_set, f"{stage} unexpected case key")
        require(row["evaluation_id"] == pair_by_id[row["truth_pair_id"]]["evaluation_id"], f"{stage} pair evaluation mismatch")
        by_pair[row["truth_pair_id"]] += 1
        count += 1
        yield key, row
    require(count == expected and (skipped == 836 if stage == "PL2B" else skipped == 0), f"{stage} case row count")
    require(len(by_pair) == 836 and set(by_pair.values()) == {4179}, f"{stage} case rectangle incomplete")
    audit[stage] = by_pair


def check_case_population(pair_by_id: dict, rid: list[str]) -> list[dict]:
    audit = {}
    streams = [case_stream("PL2B", "BRIDGEPL2B_CASE_OUTCOMES", pair_by_id, set(rid), 3493644, audit), case_stream("A22", "BRIDGEA22_CASE_OUTCOMES", pair_by_id, set(rid), 3493644, audit)]
    from itertools import zip_longest
    for a, b in zip_longest(*streams):
        require(a is not None and b is not None and a[0] == b[0], "A1/A2 case key mismatch")
        require(all(a[1][field] == b[1][field] for field in ("truth_status", "truth_direction", "truth_delta_b", "truth_magnitude")), "truth semantics mismatch")
    return [{"truth_pair_id": pid, "evaluation_id": pair_by_id[pid]["evaluation_id"], "a1_case_rows": audit["PL2B"][pid], "a2_case_rows": audit["A22"][pid], "matched": "true"} for pid in sorted(pair_by_id)]


def geometry(stage: str, prefix: str, evals: list[str], rid: list[str]) -> np.ndarray:
    n_r = len(rid)
    e_index = {v: i for i, v in enumerate(evals)}
    r_index = {v: i for i, v in enumerate(rid)}
    result = np.full((len(evals) * n_r, len(FEATURES)), np.nan, dtype=float)
    seen = np.zeros(result.shape[0], dtype=bool)
    count = 0
    for row in iter_parts(stage, prefix):
        if stage == "PL1" and row["reaction_id"] == "LDH_L":
            continue
        require(row["evaluation_id"] in e_index and row["reaction_id"] in r_index, f"{stage} geometry key outside rectangle")
        ix = e_index[row["evaluation_id"]] * n_r + r_index[row["reaction_id"]]
        require(not seen[ix], f"{stage} duplicate geometry key")
        seen[ix] = True
        if stage == "A21":
            require(row["abs_delta_anchor_target_correlation"] == row["abs_delta_hex1_target_correlation"], "A2 coupling changed from HEX1")
        result[ix] = [float(row[f]) if row[f] else math.nan for f in FEATURES]
        count += 1
    require(count == 1671600 and bool(seen.all()), f"{stage} geometry rectangle incomplete")
    return result.reshape((len(evals), n_r, len(FEATURES)))


def finite_summary(a1: np.ndarray, a2: np.ndarray) -> dict:
    x, y = np.asarray(a1, dtype=float), np.asarray(a2, dtype=float)
    require(x.shape == y.shape, "paired summary shape mismatch")
    mask = np.isfinite(x) & np.isfinite(y)
    delta = y[mask] - x[mask]
    return {"candidate_rows": int(x.size), "finite_paired_rows": int(mask.sum()), "a1_mean": float(np.mean(x[mask])) if mask.any() else None, "a2_mean": float(np.mean(y[mask])) if mask.any() else None, "paired_delta_mean": float(np.mean(delta)) if mask.any() else None, "paired_delta_median": float(np.median(delta)) if mask.any() else None}


def geometry_summary(a1: np.ndarray, a2: np.ndarray, split: dict, rid: list[str], evals: list[str], truth_evals: set[str]) -> list[dict]:
    out = []
    truth_idx = [i for i, e in enumerate(evals) if e in truth_evals]
    for j, feature in enumerate(FEATURES):
        out.append({"scope": "ALL_EVALUATION_REACTION", "feature": feature, "unit": "evaluation_x_reaction", **finite_summary(a1[:, :, j], a2[:, :, j])})
        for name in ("DEVELOPMENT", "CONFIRMATION_HOLDOUT"):
            ids = [i for i, r in enumerate(rid) if split[r] == name]
            x = pd.DataFrame(a1[np.ix_(truth_idx, ids, [j])][:, :, 0]).mean(axis=0).to_numpy()
            y = pd.DataFrame(a2[np.ix_(truth_idx, ids, [j])][:, :, 0]).mean(axis=0).to_numpy()
            out.append({"scope": name, "feature": feature, "unit": "reaction_mean_across_truth_evaluations", **finite_summary(x, y)})
    return out


def utility_and_aggregates(pair_by_id: dict, rid: list[str], split: dict, evals: list[str]):
    from itertools import zip_longest
    n_r = len(rid)
    e_index = {v: i for i, v in enumerate(evals)}
    r_index = {v: i for i, v in enumerate(rid)}
    size = len(evals) * n_r
    counts = np.zeros(size, dtype=np.int32)
    useful = [np.zeros(size, dtype=np.int32), np.zeros(size, dtype=np.int32)]
    advantage = [np.zeros(size, dtype=float), np.zeros(size, dtype=float)]
    sums = {scope: {arm: np.zeros(4, dtype=float) for arm in ("correct", "wrong", "random", "specificity")} for scope in ("ALL", "DEVELOPMENT", "CONFIRMATION_HOLDOUT")}
    nontie = Counter()
    ties = Counter()
    audit = {}
    streams = [case_stream("PL2B", "BRIDGEPL2B_CASE_OUTCOMES", pair_by_id, set(rid), 3493644, audit), case_stream("A22", "BRIDGEA22_CASE_OUTCOMES", pair_by_id, set(rid), 3493644, audit)]
    for a, b in zip_longest(*streams):
        require(a is not None and b is not None and a[0] == b[0], "case key mismatch during utility")
        ra, rb = a[1], b[1]
        require(ra["truth_status"] == rb["truth_status"], "truth tie mismatch")
        scope = split[ra["reaction_id"]]
        scopes = ("ALL", scope)
        idx = e_index[ra["evaluation_id"]] * n_r + r_index[ra["reaction_id"]]
        if ra["truth_status"] == "TRUTH_TIE":
            for s in scopes: ties[s] += 1
            continue
        require(ra["truth_status"] == "NON_TIE", "unknown truth status")
        counts[idx] += 1
        gains = []
        for row, arm in ((ra, 0), (rb, 1)):
            c = float(row["lambda_0_25_correct_absolute_error_gain"])
            w = float(row["lambda_0_25_wrong_absolute_error_gain"])
            z = float(row["lambda_0_25_random_sign_expected_gain"])
            require(all(math.isfinite(v) for v in (c, w, z)) and abs(z - .5 * (c + w)) <= 5e-12, "random-sign identity failed")
            gains.append((c, w, z, c - w))
            useful[arm][idx] += int(c > pl2c.GAIN_TOL and c - z > pl2c.GAIN_TOL)
            advantage[arm][idx] += c - z
        for s in scopes:
            nontie[s] += 1
            for j, name in enumerate(("correct", "wrong", "random", "specificity")):
                v1, v2 = gains[0][j], gains[1][j]
                sums[s][name] += (v1, v2, v2 - v1, float(v1 > pl2c.GAIN_TOL) - float(v2 > pl2c.GAIN_TOL))
    require(nontie["ALL"] + ties["ALL"] == 3493644 and nontie["ALL"] == 1440313, "A2.2 tie population changed")
    out = []
    controls = []
    for scope in ("ALL", "DEVELOPMENT", "CONFIRMATION_HOLDOUT"):
        n = nontie[scope]
        c, w, z, d = (sums[scope][name] for name in ("correct", "wrong", "random", "specificity"))
        out.append({"scope": scope, "all_case_keys": n + ties[scope], "non_tie_pairs": n, "truth_ties": ties[scope], "unit": "distinct_truth_pair_x_reaction", "a1_correct_mean_gain": c[0]/n, "a2_correct_mean_gain": c[1]/n, "paired_correct_delta": c[2]/n, "a1_specificity_mean": d[0]/n, "a2_specificity_mean": d[1]/n, "paired_specificity_delta": d[2]/n})
        for name, v in (("correct_q1", c), ("wrong_q0", w), ("analytic_random_q0_5", z)):
            controls.append({"scope": scope, "arm": name, "unit": "non_tie_distinct_truth_pair_x_reaction", "denominator": n, "a1_mean_gain": v[0]/n, "a2_mean_gain": v[1]/n, "paired_delta_mean": v[2]/n})
    return out, controls, counts.reshape((len(evals), n_r)), [x.reshape((len(evals), n_r)) for x in useful], [x.reshape((len(evals), n_r)) for x in advantage]


def reaction_frame(geom: np.ndarray, counts: np.ndarray, useful: np.ndarray, advantage: np.ndarray, nfrac: np.ndarray, rid: list[str], evals: list[str], truth_evals: set[str], indices: list[int]) -> pd.DataFrame:
    eidx = [i for i, e in enumerate(evals) if e in truth_evals]
    g = geom[np.ix_(eidx, indices, range(len(FEATURES)))]
    c = counts[np.ix_(eidx, indices)]
    u = useful[np.ix_(eidx, indices)]
    a = advantage[np.ix_(eidx, indices)]
    response = np.divide(u, c, out=np.full(c.shape, np.nan), where=c > 0)
    mean_adv = np.divide(a, c, out=np.full(c.shape, np.nan), where=c > 0)
    frame = pd.DataFrame({"reaction_id": np.tile([rid[i] for i in indices], len(eidx)), "sign_magnitude_eta2": g[:, :, 0].ravel(), "directionally_useful_fraction": response.ravel(), "non_tie_pair_fraction": nfrac[np.ix_(eidx, indices)].ravel(), "mean_information_advantage": mean_adv.ravel(), "directional_entropy3": g[:, :, 2].ravel(), "dominant_sign_mass": g[:, :, 3].ravel(), "n_supported_sign_states": g[:, :, 4].ravel()})
    return frame


def association(rows_: list[dict], x: str, y: str) -> tuple[int, float]:
    pairs = [(float(r[x]), float(r[y])) for r in rows_ if r.get(x) is not None and r.get(y) is not None and math.isfinite(float(r[x])) and math.isfinite(float(r[y]))]
    return len(pairs), hcore.spearman([v[0] for v in pairs], [v[1] for v in pairs])


def persistent(geometry: list[np.ndarray], counts: np.ndarray, useful: list[np.ndarray], advantage: list[np.ndarray], pair_by_id: dict, rid: list[str], split: dict, evals: list[str], eval_by_id: dict, reactions: list[dict]):
    # Counts per evaluation are frozen by the matched registry, independent of outcomes.
    pair_count = Counter(r["evaluation_id"] for r in pair_by_id.values())
    truth_evals = set(pair_count)
    eidx = [i for i, e in enumerate(evals) if e in truth_evals]
    nfrac = np.full(counts.shape, np.nan, dtype=float)
    for i in eidx:
        nfrac[i] = counts[i] / pair_count[evals[i]]
    ids = {name: [i for i, r in enumerate(rid) if split[r] == name] for name in ("DEVELOPMENT", "CONFIRMATION_HOLDOUT")}
    react_rows = {}
    frames = {}
    for arm in (0, 1):
        for name, rx_idx in ids.items():
            frame = reaction_frame(geometry[arm], counts, useful[arm], advantage[arm], nfrac, rid, evals, truth_evals, rx_idx)
            frames[arm, name] = frame
            react_rows[arm, name] = hcore.aggregate_reaction_frame(frame)
    out = []
    # PL2C evaluation x reaction development association, on the frozen feature set.
    for j, feature in enumerate(FEATURES):
        x1 = geometry[0][np.ix_(eidx, ids["DEVELOPMENT"], [j])][:, :, 0].ravel()
        x2 = geometry[1][np.ix_(eidx, ids["DEVELOPMENT"], [j])][:, :, 0].ravel()
        y1 = frames[0, "DEVELOPMENT"]["directionally_useful_fraction"].to_numpy()
        y2 = frames[1, "DEVELOPMENT"]["directionally_useful_fraction"].to_numpy()
        r1, r2 = pl2c.spearman_finite(x1, y1), pl2c.spearman_finite(x2, y2)
        out.append({"level": "PL2C_EVALUATION_REACTION", "split": "DEVELOPMENT", "context": "ALL", "feature": feature, "a1_rho": r1["rho"], "a2_rho": r2["rho"], "a1_finite": r1["finite_count"], "a2_finite": r2["finite_count"], "rho_a2_minus_a1": r2["rho"] - r1["rho"]})
    for name in ("DEVELOPMENT", "CONFIRMATION_HOLDOUT"):
        a1n, a1rho = association(react_rows[0, name], "mean_sign_magnitude_eta2", "mean_directionally_useful_fraction")
        a2n, a2rho = association(react_rows[1, name], "mean_sign_magnitude_eta2", "mean_directionally_useful_fraction")
        out.append({"level": "PL2D_REACTION", "split": name, "context": "ALL", "feature": "sign_magnitude_eta2", "a1_rho": a1rho, "a2_rho": a2rho, "a1_finite": a1n, "a2_finite": a2n, "rho_a2_minus_a1": a2rho-a1rho})
        for context in sorted({(r["algorithm"], r["rna_context_key"]) for r in eval_by_id.values()}):
            ctx_e = {e for e in truth_evals if (eval_by_id[e]["algorithm"], eval_by_id[e]["rna_context_key"]) == context}
            context_rows = []
            for arm in (0, 1):
                sub = frames[arm, name].iloc[np.flatnonzero(np.repeat(np.array([evals[i] in ctx_e for i in eidx]), len(ids[name])))]
                rr = hcore.aggregate_reaction_frame(sub)
                context_rows.append(association(rr, "mean_sign_magnitude_eta2", "mean_directionally_useful_fraction"))
            (n1, rho1), (n2, rho2) = context_rows
            out.append({"level": "PL2D_CONTEXT_REACTION", "split": name, "context": context[0] + "|" + context[1], "feature": "sign_magnitude_eta2", "a1_rho": rho1, "a2_rho": rho2, "a1_finite": n1, "a2_finite": n2, "rho_a2_minus_a1": rho2-rho1})
    sensitivity = []
    annotation = {r["reaction_id"]: r for r in reactions}
    for name in ("DEVELOPMENT", "CONFIRMATION_HOLDOUT"):
        for label, flags in (("exclude_proximal", ("proximal_set_member",)), ("exclude_same_subsystem", ("same_subsystem_as_HEX1",)), ("exclude_either", ("proximal_set_member", "same_subsystem_as_HEX1"))):
            selected = [r for r in react_rows[0, name] if not any(annotation[r["reaction_id"]][flag].lower() == "true" for flag in flags)]
            selected2 = [r for r in react_rows[1, name] if not any(annotation[r["reaction_id"]][flag].lower() == "true" for flag in flags)]
            n1, v1 = association(selected, "mean_sign_magnitude_eta2", "mean_directionally_useful_fraction")
            n2, v2 = association(selected2, "mean_sign_magnitude_eta2", "mean_directionally_useful_fraction")
            sensitivity.append({"split": name, "analysis": label, "a1_rho": v1, "a2_rho": v2, "a1_finite": n1, "a2_finite": n2})
        for label, x, y in (("eta2_information_advantage", "mean_sign_magnitude_eta2", "mean_information_advantage"), ("entropy_usefulness", "mean_directional_entropy3", "mean_directionally_useful_fraction"), ("dominant_mass_usefulness", "mean_dominant_sign_mass", "mean_directionally_useful_fraction"), ("sign_states_usefulness", "mean_n_supported_sign_states", "mean_directionally_useful_fraction")):
            n1, v1 = association(react_rows[0, name], x, y)
            n2, v2 = association(react_rows[1, name], x, y)
            sensitivity.append({"split": name, "analysis": label, "a1_rho": v1, "a2_rho": v2, "a1_finite": n1, "a2_finite": n2})
        n1, v1 = dcore.partial_primary_controlling_coverage(react_rows[0, name])
        n2, v2 = dcore.partial_primary_controlling_coverage(react_rows[1, name])
        sensitivity.append({"split": name, "analysis": "partial_eta2_coverage", "a1_rho": v1, "a2_rho": v2, "a1_finite": n1, "a2_finite": n2})
    return out, sensitivity, react_rows


def verify_a1_parity(rows_: list[dict]) -> None:
    frozen = list(rows(DIRS["PL2D"] / "BRIDGEPL2D_CONFIRMATION_REACTION.tsv"))
    actual = {r["reaction_id"]: r for r in rows_}
    require(len(frozen) == len(actual) == 1045, "confirmation reaction count mismatch")
    for row in frozen:
        got = actual[row["reaction_id"]]
        for field in hcore.REACTION_MEAN_OUTPUTS.values():
            expected = float(row[field]) if row[field] else math.nan
            value = got[field] if got[field] is not None else math.nan
            require((math.isnan(expected) and math.isnan(value)) or math.isclose(expected, value, rel_tol=0, abs_tol=1e-12), f"A1 parity failed: {row['reaction_id']} {field}")


def existing_noop(audit: dict) -> bool:
    if not OUT.exists():
        return False
    require(OUT.is_dir() and not OUT.is_symlink(), "conflicting A2.3 output path")
    manifest_path = OUT / "BRIDGEA23_MANIFEST.json"
    require(manifest_path.is_file() and not manifest_path.is_symlink(), "partial A2.3 output exists")
    manifest = read_json(manifest_path)
    require(manifest["status"] == "BRIDGE_A23_A1_A2_SYNTHESIS_FROZEN" and manifest["source_sha256"] == audit["manifest_sha256"] and manifest["implementation_sha256"] == sha(Path(__file__)), "conflicting A2.3 identity")
    require(set(p.name for p in OUT.iterdir()) == set(FILES) | {manifest_path.name}, "conflicting A2.3 files")
    for name, digest in manifest["artifact_sha256"].items():
        checked(OUT / name, digest)
    return True


def run() -> dict:
    manifests, source = admit()
    if existing_noop(source):
        return {"status": "BRIDGE_A23_A1_A2_SYNTHESIS_FROZEN", "rerun": "NO_OP"}
    pairs, reactions, split, evals, eval_by_id = identity(manifests)
    rid = [r["reaction_id"] for r in reactions]
    matched = check_case_population(pairs, rid)
    # Only after exact population admission may numerical synthesis be accessed.
    geometry_data = [geometry("PL1", "BRIDGEPL1_REACTION_EVALUATION_FEATURES", evals, rid), geometry("A21", "BRIDGEA21_REACTION_EVALUATION_FEATURES", evals, rid)]
    utility, controls, counts, useful, advantage = utility_and_aggregates(pairs, rid, split, evals)
    truth_evals = {p["evaluation_id"] for p in pairs.values()}
    geo_summary = geometry_summary(*geometry_data, split, rid, evals, truth_evals)
    persistence, sensitivity, react_rows = persistent(geometry_data, counts, useful, advantage, pairs, rid, split, evals, eval_by_id, reactions)
    verify_a1_parity(react_rows[0, "CONFIRMATION_HOLDOUT"])
    bootstrap = dcore.bootstrap_primary(react_rows[1, "CONFIRMATION_HOLDOUT"])
    primary_rows = [r for r in persistence if r["level"] == "PL2D_REACTION"]
    qc = {"schema": "bridge.a23.qc.v1", "status": "PASS", "matched_truth_pairs": 836, "audit_only_truth_pairs_excluded": 49, "reactions": 4179, "matched_case_keys": 3493644, "geometry_keys": 1671600, "development_reactions": 3134, "confirmation_reactions": 1045, "primary_finite_reactions": {r["split"]: {"a1": r["a1_finite"], "a2": r["a2_finite"]} for r in primary_rows}, "a2_joint_ess_range_from_frozen_qc": [read_json(DIRS["A21"] / "BRIDGEA21_QC.json")["conditioned_joint_ess_min"], read_json(DIRS["A21"] / "BRIDGEA21_QC.json")["conditioned_joint_ess_max"]], "a1_confirmation_parity": "PASS", "a2_g_used": False, "bootstrap": {k: bootstrap[k] for k in ("replicates", "seed", "rng", "finite_replicates", "nonfinite_replicates", "rho_q025", "rho_median", "rho_q975")}}
    contract = {"schema": "bridge.a23.analysis_contract.v1", "role": "post_freeze_descriptive_robustness_not_independent_confirmation", "primary_lambda": .25, "primary_reliability_q": 1, "gain_tolerance": 1e-12, "truth_tie_tolerance": 1e-12, "response": "fraction_non_tie_distinct_pairs_with_correct_gain_gt_tol_and_information_advantage_gt_tol", "predictor": "reaction_mean_sign_magnitude_eta2", "development_reactions": 3134, "confirmation_reactions": 1045, "features": FEATURES, "primary_numerical_reference": "frozen_PL2D_H_pandas_groupby_mean_and_scipy_spearman", "paired_delta": "A2_minus_A1", "bootstrap_role": "supportive_computational_stability_only", "versions": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__}}
    require(not OUT.exists(), "A2.3 output appeared during synthesis")
    OUT.mkdir(parents=False)
    write_json(OUT / FILES[0], source)
    write_json(OUT / FILES[1], contract)
    for name, data in zip(FILES[2:8], (matched, geo_summary, utility, persistence, controls, sensitivity)):
        require(len(data) == EXPECTED_ROWS[name], f"{name} row count changed")
        write_tsv(OUT / name, data)
    write_json(OUT / FILES[8], qc)
    artifact_sha256 = {name: sha(OUT / name) for name in FILES}
    manifest = {"schema": "bridge.a23.manifest.v1", "status": "BRIDGE_A23_A1_A2_SYNTHESIS_FROZEN", "source_sha256": source["manifest_sha256"], "implementation_sha256": sha(Path(__file__)), "artifact_sha256": artifact_sha256, "row_counts": EXPECTED_ROWS, "case_population": {"truth_pairs": 836, "reactions": 4179, "case_keys": 3493644}, "a2_g_used": False, "new_production_performed": False}
    write_json(OUT / "BRIDGEA23_MANIFEST.json", manifest)
    return {"status": manifest["status"], "output": str(OUT.relative_to(ROOT)), "row_counts": EXPECTED_ROWS}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    print(json.dumps(run(), sort_keys=True))


if __name__ == "__main__":
    main()
