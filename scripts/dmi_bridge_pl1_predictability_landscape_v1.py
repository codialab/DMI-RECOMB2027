#!/usr/bin/env python3
"""Build the outcome-blind, network-wide PL1 feature landscape.

This adapter only joins the frozen candidate/cache/weight records and applies the
scalar functions in dmi_bridge_pl1_predictability_core_v1.  It has no model or
solver dependency.  The publication directory is assembled beside the final
directory and is never partially published.
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

from scripts.dmi_bridge_pl1_predictability_batch_v1 import target_landscape_metrics_block
from scripts.dmi_bridge_pl1_predictability_core_v1 import (
    COMPLETE_STATUS,
    STRONG_ANCHOR_REACTION,
    analysis_contract,
    target_landscape_metrics,
)
from scripts import dmi_bridge_pl1_storage_v1 as storage

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/dmi_bridge_pl1_predictability_landscape_v1"
CORE_SCHEMA = "bridge.pl1.predictability_landscape.v1"
METHODS = ("CORDA", "GIMME", "RIPTiDe", "iMAT")
TUMORS = ("CT2A", "GL261")
REACTION_BLOCK_SIZE = 64
XZ_PRESET = 0
SCALAR_PARITY_RTOL = 5e-12
SCALAR_PARITY_ATOL = 5e-12
REACTION_STATS = (
    "delta_b_mean", "delta_b_sd", "delta_b_q10", "delta_b_q50", "delta_b_q90",
    "delta_b_width80", "uniform_delta_b_mean", "uniform_delta_b_sd",
    "uniform_delta_b_q10", "uniform_delta_b_q50", "uniform_delta_b_q90",
    "uniform_delta_b_width80", "anchor_sd_contraction", "anchor_width80_contraction",
    "p_negative", "p_tie", "p_positive", "dominant_sign_mass", "directional_entropy3",
    "joint_ess", "negative_ess", "tie_ess", "positive_ess", "negative_abs_mean",
    "negative_abs_q10", "negative_abs_q50", "negative_abs_q90", "negative_abs_width80",
    "positive_abs_mean", "positive_abs_q10", "positive_abs_q50", "positive_abs_q90",
    "positive_abs_width80", "negative_width80_contraction_vs_all",
    "positive_width80_contraction_vs_all", "sign_magnitude_eta2",
    "sign_magnitude_correlation", "delta_anchor_target_correlation",
    "abs_delta_anchor_target_correlation",
)
DEGENERACY_COUNT_FIELDS = (
    "target_delta_degenerate_count",
    "target_abs_magnitude_degenerate_count",
    "anchor_delta_degenerate_count",
    "sign_magnitude_degenerate_count",
    "anchor_target_coupling_degenerate_count",
)

PINNED = {
    "outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_MANIFEST.json": "40df8378876820bb616fdbdc8ef8a0f77d43e36f464104fca6cae9a16b39d7c1",
    "outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_ANALYSIS_CONTRACT.json": "5e991aa58275e1545628bd9b3e9014b20c1fe799bfcc9810a0d6d712a28c8e2a",
    "outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_SOURCE_AUDIT.json": "acf214cf3bbdb9a0737673e7808696da83f7b52cef599d43f660c21ba5b59787",
    "outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_SUPPORT_AUDIT.tsv": "340c6ef821a093d4d7949b6488cdba2fa6cf9c3300b54b9a8b057016339de69e",
    "outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_REACTION_PANEL.tsv": "f6e881ef949b2a519f0ac731d0dd8b08dba0d44f40e24ca94e18a47017907f4d",
    "outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_FOUR_ARM_WEIGHTS.tsv.xz": "843e17b7012a6fa212d1e9ac9e1dc6e309e854d94ab5c4b3f9e316b589768f4f",
    "outputs/dmi_bridge_bio0r_qualification_v1/BRIDGEBIO0R_MANIFEST.json": "d155b63f59db3479b1dfd2d87c0b9ca470358c0166496abf288a649f04890cea",
    "outputs/dmi_bridge_bio0r_qualification_v1/BRIDGEBIO0R_IDENTITY_AUDIT.json": "5b816176e2532ff076defd5fab103afa4210d9162990c85812209a7f6fa6f4a1",
    "outputs/dmi_bridge_bio0r_qualification_v1/BRIDGEBIO0R_SOURCE_AUDIT.json": "33ab67e8d8d12b84cdee1fc43b36bdba6b3680e9451bc01ed071e5fecd1c6fc4",
    "12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/candidate_fraction_table.tsv.xz": "421c2c1d196fbd088a4e92e2ecc998ad578aab814a025fff8e1c843597fed1da",
    "outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstructed/flux_cache.npz": "f9644da4ebd875cf5f0b34cc58b20faf0a0457d7c7a44635e8a010d8896c3dbb",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def source_hashes() -> dict[str, str]:
    got = {}
    for rel, expected in PINNED.items():
        p = ROOT / rel
        if not p.is_file():
            raise RuntimeError(f"missing pinned source: {rel}")
        actual = sha256(p)
        if actual != expected:
            raise RuntimeError(f"pinned source hash mismatch: {rel}")
        got[rel] = actual
    return got


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = lzma.open if path.suffix == ".xz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as fh:
        out = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        out.writeheader()
        for row in rows:
            out.writerow({k: f(row[k]) for k in fields})


def write_xz_tsv(path: Path, fields: list[str], rows: Iterable[dict]) -> None:
    # Deterministic streaming XZ avoids a large uncompressed sibling and the
    # expensive default preset. Compression level is an execution detail only.
    with lzma.open(path, "wt", encoding="utf-8", newline="", preset=XZ_PRESET) as fh:
        out = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        out.writeheader()
        for row in rows:
            out.writerow({k: f(row[k]) for k in fields})


def write_partitioned_feature_tsv(directory: Path, fields: list[str], rows: Iterable[dict], total_rows: int) -> dict:
    return storage.write_partitioned_tsv(directory, fields, rows, total_rows=total_rows)


def open_xz_tsv(path: Path, fields: list[str]):
    fh = lzma.open(path, "wt", encoding="utf-8", newline="", preset=XZ_PRESET)
    out = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
    out.writeheader()
    return fh, out


def f(value: object) -> str:
    if isinstance(value, (float, np.floating)):
        return "nan" if not math.isfinite(float(value)) else format(float(value), ".17g")
    return str(value)


def qstats(values: list[float]) -> dict[str, object]:
    a = np.asarray(values, dtype=float)
    finite = np.isfinite(a)
    x = a[finite]
    return {
        "median": float(np.median(x)) if x.size else float("nan"),
        "q10": float(np.quantile(x, .1, method="linear")) if x.size else float("nan"),
        "q90": float(np.quantile(x, .9, method="linear")) if x.size else float("nan"),
        "finite_count": int(x.size),
        "nonfinite_count": int(np.count_nonzero(~finite)),
    }


def _check_statuses() -> None:
    bio0 = json.loads((ROOT / "outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_MANIFEST.json").read_text())
    bio0r = json.loads((ROOT / "outputs/dmi_bridge_bio0r_qualification_v1/BRIDGEBIO0R_MANIFEST.json").read_text())
    if bio0.get("status") != "BIO0_READY_FOR_BIO1_FOUR_ARM_DECOMPOSITION":
        raise RuntimeError("BIO-0 is not currently admitted")
    if bio0r.get("status") != "BIO0R_QUALIFIED_READY_FOR_BIO1":
        raise RuntimeError("BIO-0R is not currently qualified")
    if not bio0r.get("qualification_gates", {}).get("current_vmax_semantics_verified"):
        raise RuntimeError("current Vmax semantics are not verified")
    if not bio0r.get("qualification_gates", {}).get("ensemble_projection_cache_identity_verified"):
        raise RuntimeError("ensemble/cache identity is not verified")


def load_inputs() -> tuple[list[dict], dict, dict, np.ndarray, np.ndarray, dict[str, str]]:
    hashes = source_hashes()
    _check_statuses()
    candidates = read_tsv(ROOT / "12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/candidate_fraction_table.tsv.xz")
    weights = read_tsv(ROOT / "outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_FOUR_ARM_WEIGHTS.tsv.xz")
    panel_rows = read_tsv(ROOT / "outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_REACTION_PANEL.tsv")
    panel = {r["reaction_id"]: r["subsystem"] for r in panel_rows}
    if len(panel) != 3893:
        raise RuntimeError("reaction panel cardinality mismatch")
    z = np.load(ROOT / "outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstructed/flux_cache.npz", allow_pickle=False)
    rxns, mats = z["rxns"].tolist(), z["mats"]
    if len(rxns) != 4181 or len(set(rxns)) != len(rxns) or rxns.count(STRONG_ANCHOR_REACTION) != 1:
        raise RuntimeError("cache reaction inventory is invalid")
    if mats.shape != (32, 20, 4181) or not np.all(np.isfinite(mats)):
        raise RuntimeError("cache matrix shape or finiteness mismatch")
    cache = {(str(a), str(t), str(e)): i for i, (a, t, e) in enumerate(zip(z["alg"], z["tumor"], z["eh"]))}
    if len(cache) != 32:
        raise RuntimeError("cache stratum identities are not unique")
    ckey = lambda r: (r["algorithm"], r["tumor"], r["ensemble_hash"])
    cgroups = defaultdict(list)
    for r in candidates:
        cgroups[ckey(r)].append(r)
    if len(candidates) != 640 or len(cgroups) != 32 or any(len(v) != 20 for v in cgroups.values()):
        raise RuntimeError("candidate inventory is not exactly 32x20")
    by_key = {}
    for key, rows in cgroups.items():
        rows.sort(key=lambda r: int(r["sample_index"]))
        if [int(r["sample_index"]) for r in rows] != list(range(20)):
            raise RuntimeError("candidate sample indices are incomplete")
        by_key[key] = rows
        if key not in cache:
            raise RuntimeError("candidate/cache identity join is incomplete")
    if set(cache) != set(by_key):
        raise RuntimeError("cache contains an unjoined identity")
    strong = [r for r in weights if r["arm"] == "strong_anchor_baseline"]
    if len(strong) != 3200:
        raise RuntimeError("strong-anchor weight row count mismatch")
    wm = defaultdict(list)
    for r in strong:
        key = (r["mouse_id"], r["tumor"], r["algorithm"], r["ensemble_hash"], r["rna_context_key"])
        wm[key].append((int(r["sample_index"]), float(r["weight"])))
    if len(wm) != 160 or any(len(v) != 20 for v in wm.values()):
        raise RuntimeError("strong-anchor support is incomplete")
    for key, rows in wm.items():
        if sorted(i for i, _ in rows) != list(range(20)) or any(not math.isfinite(w) or w < 0 for _, w in rows):
            raise RuntimeError("invalid strong-anchor support")
    return candidates, by_key, wm, np.asarray(rxns), mats, cache, hashes


def _summary_rows(groups: dict[tuple, list[dict]], keys: list[str], prefix: str = "") -> list[dict]:
    rows = []
    for key in sorted(groups):
        vals = groups[key]
        row = dict(zip(keys, key))
        for metric in REACTION_STATS:
            s = qstats([float(v[metric]) for v in vals])
            row[f"{prefix}{metric}_median"] = s["median"]
            row[f"{prefix}{metric}_q10"] = s["q10"]
            row[f"{prefix}{metric}_q90"] = s["q90"]
            row[f"{prefix}{metric}_finite_count"] = s["finite_count"]
            row[f"{prefix}{metric}_degenerate_count"] = s["degenerate_count"]
        rows.append(row)
    return rows


def _stats_row(row: dict, vectors: np.ndarray) -> dict:
    for j, metric in enumerate(REACTION_STATS):
        s = qstats(vectors[:, j].tolist())
        for suffix, value in s.items():
            row[f"{metric}_{suffix}"] = value
    return row


def _summary_arrays(values: np.ndarray) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """Vectorized finite-value summaries over axis 0 for reaction x metric arrays."""
    out = {}
    for j, metric in enumerate(REACTION_STATS):
        x = np.asarray(values[:, :, j], dtype=float)
        finite = np.isfinite(x)
        count = np.sum(finite, axis=0)
        with np.errstate(invalid="ignore"), np.testing.suppress_warnings() as sup:
            sup.filter(RuntimeWarning)
            q10 = np.nanquantile(x, 0.10, axis=0, method="linear")
            med = np.nanmedian(x, axis=0)
            q90 = np.nanquantile(x, 0.90, axis=0, method="linear")
        out[metric] = (med, q10, q90, count, x.shape[0] - count)
    return out


def _summary_row_from_arrays(
    base: dict,
    summaries: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    rpos: int,
) -> dict:
    row = dict(base)
    for metric in REACTION_STATS:
        med, q10, q90, count, deg = summaries[metric]
        row[f"{metric}_median"] = med[rpos]
        row[f"{metric}_q10"] = q10[rpos]
        row[f"{metric}_q90"] = q90[rpos]
        row[f"{metric}_finite_count"] = int(count[rpos])
        row[f"{metric}_nonfinite_count"] = int(deg[rpos])
    return row


def _degeneracy_count_arrays(bool_values: np.ndarray, status_values: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "target_delta_degenerate_count": np.sum(bool_values[:, :, 0] != 0, axis=0),
        "target_abs_magnitude_degenerate_count": np.sum(bool_values[:, :, 1] != 0, axis=0),
        "anchor_delta_degenerate_count": np.sum(bool_values[:, :, 2] != 0, axis=0),
        "sign_magnitude_degenerate_count": np.sum(status_values[:, :, 0] == 0, axis=0),
        "anchor_target_coupling_degenerate_count": np.sum(status_values[:, :, 1] == 0, axis=0),
    }


def _assert_scalar_parity(
    *,
    batch: dict[str, np.ndarray],
    local_pos: int,
    anchor_ct2a: np.ndarray,
    anchor_gl261: np.ndarray,
    target_ct2a: np.ndarray,
    target_gl261: np.ndarray,
    strong_weights_ct2a: np.ndarray,
    strong_weights_gl261: np.ndarray,
) -> None:
    scalar = target_landscape_metrics(
        anchor_ct2a=anchor_ct2a,
        anchor_gl261=anchor_gl261,
        target_ct2a=target_ct2a,
        target_gl261=target_gl261,
        strong_weights_ct2a=strong_weights_ct2a,
        strong_weights_gl261=strong_weights_gl261,
    )
    for name, expected in scalar.items():
        got = batch[name][local_pos]
        if isinstance(expected, str):
            if str(got) != expected:
                raise RuntimeError(f"PL1 batch/scalar status mismatch for {name}")
        elif isinstance(expected, bool):
            if bool(got) is not expected:
                raise RuntimeError(f"PL1 batch/scalar boolean mismatch for {name}")
        elif isinstance(expected, int):
            if int(got) != expected:
                raise RuntimeError(f"PL1 batch/scalar integer mismatch for {name}")
        elif math.isnan(float(expected)):
            if not math.isnan(float(got)):
                raise RuntimeError(f"PL1 batch/scalar NaN mismatch for {name}")
        elif not np.isclose(float(got), float(expected), rtol=SCALAR_PARITY_RTOL, atol=SCALAR_PARITY_ATOL):
            raise RuntimeError(f"PL1 batch/scalar numeric mismatch for {name}: {got} vs {expected}")


def build(out: Path = OUT) -> dict:
    candidates, by_key, wm, rxns, mats, cache_index, hashes = load_inputs()
    panel_rows = read_tsv(ROOT / "outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_REACTION_PANEL.tsv")
    panel = {r["reaction_id"]: r["subsystem"] for r in panel_rows}
    hex_idx = int(np.flatnonzero(rxns == STRONG_ANCHOR_REACTION)[0])
    primary_rxns = sorted(r for r in rxns.tolist() if r != STRONG_ANCHOR_REACTION)
    proximal = {"Glycolysis/gluconeogenesis", "Pyruvate metabolism", "Citric acid cycle", "Glutamate metabolism"}
    registry: list[dict] = []
    rxn_index = {r: i for i, r in enumerate(rxns.tolist())}
    probe = target_landscape_metrics(anchor_ct2a=[0], anchor_gl261=[0], target_ct2a=[0], target_gl261=[0], strong_weights_ct2a=[1], strong_weights_gl261=[1])
    metric_fields = [k for k, v in probe.items() if isinstance(v, (int, float, np.number)) and not isinstance(v, bool)]
    status_fields = ["target_delta_degenerate", "target_abs_magnitude_degenerate", "anchor_delta_degenerate", "sign_magnitude_status", "anchor_target_coupling_status"]
    metric_pos = {m: i for i, m in enumerate(metric_fields)}
    stat_pos = {m: i for i, m in enumerate(REACTION_STATS)}
    work_files = []
    def mapped(shape, dtype):
        fd, name = tempfile.mkstemp(prefix="pl1-metrics-", dir="/tmp")
        os.close(fd)
        work_files.append(Path(name))
        return np.memmap(name, mode="w+", dtype=dtype, shape=shape)
    all_metrics = mapped((400, len(primary_rxns), len(metric_fields)), "f8")
    bool_metrics = mapped((400, len(primary_rxns), 3), "u1")
    status_metrics = mapped((400, len(primary_rxns), 2), "u1")
    context_labels: list[tuple[str, str]] = []
    feature_fields = ["evaluation_id", "algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse", "reaction_id", "subsystem", "same_subsystem_as_HEX1", "proximal_set_member"] + metric_fields + status_fields
    stratum_keys = sorted(by_key, key=lambda k: (k[0], k[1], by_key[k][0]["rna_context_key"]))
    eval_id = 0
    for key in stratum_keys:
        method, tumor, ensemble = key
        rows = by_key[key]
        rna = rows[0]["rna_context_key"]
        if any(r["rna_context_key"] != rna for r in rows):
            raise RuntimeError("candidate stratum has mixed RNA identities")
        other = "GL261" if tumor == "CT2A" else "CT2A"
        other_keys = [k for k in stratum_keys if k[0] == method and k[1] == other and by_key[k][0]["rna_context_key"] == rna]
        if len(other_keys) != 1:
            raise RuntimeError("paired tumor stratum is missing")
        other_key = other_keys[0]
        if tumor != "CT2A":
            continue
        ckey, gkey = key, other_key
        ci, gi = cache_index[ckey], cache_index[gkey]
        crows, grows = by_key[ckey], by_key[gkey]
        cmice = sorted({k[0] for k in wm if k[1] == "CT2A" and k[2] == method and k[3] == ensemble and k[4] == rna})
        gmice = sorted({k[0] for k in wm if k[1] == "GL261" and k[2] == method and k[3] == grows[0]["ensemble_hash"] and k[4] == rna})
        if len(cmice) != 5 or len(gmice) != 5:
            raise RuntimeError("expected five mice per tumor")
        context_labels.append((method, rna))
        for cmouse in cmice:
            wc = np.asarray([w for _, w in sorted(wm[(cmouse, "CT2A", method, ensemble, rna)])], dtype=float)
            for gmouse in gmice:
                ge = grows[0]["ensemble_hash"]
                wg = np.asarray([w for _, w in sorted(wm[(gmouse, "GL261", method, ge, rna)])], dtype=float)
                eid = f"E{eval_id:04d}"
                eval_id += 1
                registry.append({"evaluation_id": eid, "algorithm": method, "rna_context_key": rna, "ct2a_mouse": cmouse, "gl261_mouse": gmouse, "cartesian_support": 400, "joint_weight_sum": 1.0, "status": "ADMITTED"})
                anchor_c = mats[ci, :, hex_idx]
                anchor_g = mats[gi, :, hex_idx]
                for start in range(0, len(primary_rxns), REACTION_BLOCK_SIZE):
                    stop = min(start + REACTION_BLOCK_SIZE, len(primary_rxns))
                    block_rxns = primary_rxns[start:stop]
                    cols = [rxn_index[r] for r in block_rxns]
                    target_c = mats[ci][:, cols]
                    target_g = mats[gi][:, cols]
                    block = target_landscape_metrics_block(
                        anchor_ct2a=anchor_c,
                        anchor_gl261=anchor_g,
                        target_ct2a=target_c,
                        target_gl261=target_g,
                        strong_weights_ct2a=wc,
                        strong_weights_gl261=wg,
                    )
                    for name in metric_fields:
                        all_metrics[eval_id - 1, start:stop, metric_pos[name]] = np.asarray(block[name], dtype=float)
                    for j, name in enumerate(status_fields[:3]):
                        bool_metrics[eval_id - 1, start:stop, j] = np.asarray(block[name], dtype=bool)
                    for j, name in enumerate(status_fields[3:]):
                        status_metrics[eval_id - 1, start:stop, j] = np.asarray(block[name] == np.asarray("OK"), dtype=np.uint8)
                    # One deterministic sentinel per block on the first mouse pair
                    # of each context preserves the scalar core as scientific authority.
                    if cmouse == cmice[0] and gmouse == gmice[0]:
                        sentinel = 0
                        _assert_scalar_parity(
                            batch=block,
                            local_pos=sentinel,
                            anchor_ct2a=anchor_c,
                            anchor_gl261=anchor_g,
                            target_ct2a=target_c[:, sentinel],
                            target_gl261=target_g[:, sentinel],
                            strong_weights_ct2a=wc,
                            strong_weights_gl261=wg,
                        )
    if eval_id != 400:
        raise RuntimeError(f"unexpected PL1 evaluation count: {eval_id}")
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and not out.is_dir():
        raise RuntimeError("PL1 output path exists and is not a directory")
    tmp = Path(tempfile.mkdtemp(prefix=out.name + ".tmp.", dir=out.parent))
    try:
        reg_fields = list(registry[0])
        write_tsv(tmp / "BRIDGEPL1_EVALUATION_REGISTRY.tsv", reg_fields, registry)
        def feature_rows():
            for e, reg in enumerate(registry):
                for rpos, rid in enumerate(primary_rxns):
                    subsystem = panel.get(rid, "UNMAPPED_PATHWAY")
                    row = {"evaluation_id": reg["evaluation_id"], "algorithm": reg["algorithm"], "rna_context_key": reg["rna_context_key"], "ct2a_mouse": reg["ct2a_mouse"], "gl261_mouse": reg["gl261_mouse"], "reaction_id": rid, "subsystem": subsystem, "same_subsystem_as_HEX1": subsystem == panel.get(STRONG_ANCHOR_REACTION), "proximal_set_member": subsystem in proximal}
                    row.update({m: all_metrics[e, rpos, metric_pos[m]] for m in metric_fields})
                    row.update({m: bool_metrics[e, rpos, j] for j, m in enumerate(status_fields[:3])})
                    row.update({m: ("OK" if status_metrics[e, rpos, j] else "DEGENERATE") for j, m in enumerate(status_fields[3:])})
                    yield {k: f(row[k]) for k in feature_fields}
        feature_manifest = write_partitioned_feature_tsv(tmp, feature_fields, feature_rows(), 1_672_000)
        cfields = ["algorithm", "rna_context_key", "reaction_id"] + [f"{m}_{s}" for m in REACTION_STATS for s in ("median", "q10", "q90", "finite_count", "nonfinite_count")] + list(DEGENERACY_COUNT_FIELDS)
        cfh, cwriter = open_xz_tsv(tmp / "BRIDGEPL1_CONTEXT_REACTION_SUMMARY.tsv.xz", cfields)
        try:
            stat_indices = [metric_pos[m] for m in REACTION_STATS]
            for cpos, (method, rna) in enumerate(context_labels):
                start = cpos * 25
                values = np.asarray(all_metrics[start:start + 25, :, :][:, :, stat_indices])
                summaries = _summary_arrays(values)
                degeneracy = _degeneracy_count_arrays(
                    np.asarray(bool_metrics[start:start + 25, :, :]),
                    np.asarray(status_metrics[start:start + 25, :, :]),
                )
                for rpos, rid in enumerate(primary_rxns):
                    row = _summary_row_from_arrays(
                        {"algorithm": method, "rna_context_key": rna, "reaction_id": rid}, summaries, rpos
                    )
                    for name in DEGENERACY_COUNT_FIELDS:
                        row[name] = int(degeneracy[name][rpos])
                    cwriter.writerow({k: f(row[k]) for k in cfields})
        finally:
            cfh.close()

        stat_indices = [metric_pos[m] for m in REACTION_STATS]
        reaction_summaries = _summary_arrays(np.asarray(all_metrics[:, :, :][:, :, stat_indices]))
        reaction_degeneracy = _degeneracy_count_arrays(np.asarray(bool_metrics), np.asarray(status_metrics))
        rrows = []
        for rpos, rid in enumerate(primary_rxns):
            row = _summary_row_from_arrays({"reaction_id": rid}, reaction_summaries, rpos)
            for name in DEGENERACY_COUNT_FIELDS:
                row[name] = int(reaction_degeneracy[name][rpos])
            rrows.append(row)
        rfields = ["reaction_id"] + [f"{m}_{s}" for m in REACTION_STATS for s in ("median", "q10", "q90", "finite_count", "nonfinite_count")] + list(DEGENERACY_COUNT_FIELDS)
        write_tsv(tmp / "BRIDGEPL1_REACTION_SUMMARY.tsv", rfields, rrows)
        pathway_positions = defaultdict(list)
        for rpos, rid in enumerate(primary_rxns):
            pathway_positions[panel.get(rid, "UNMAPPED_PATHWAY")].append(rpos)
        pathway_vectors = defaultdict(list)
        for rpos, rid in enumerate(primary_rxns):
            subsystem = panel.get(rid, "UNMAPPED_PATHWAY")
            pathway_vectors[subsystem].append([r[f"{m}_median"] for m in REACTION_STATS for r in [rrows[rpos]]])
        prows = [_stats_row({"subsystem": subsystem}, np.asarray(vectors, dtype=float)) for subsystem, vectors in sorted(pathway_vectors.items())]
        pfields = ["subsystem"] + [f"{m}_{s}" for m in REACTION_STATS for s in ("median", "q10", "q90", "finite_count", "nonfinite_count")]
        write_tsv(tmp / "BRIDGEPL1_PATHWAY_SUMMARY.tsv", pfields, prows)
        sensitivity = []
        for scope, predicate in (("primary", lambda r: True), ("exclude_same_subsystem", lambda r: not r["same_subsystem_as_HEX1"]), ("exclude_proximal_set", lambda r: not r["proximal_set_member"])):
            positions = [i for i, rid in enumerate(primary_rxns) if predicate({"same_subsystem_as_HEX1": panel.get(rid, "UNMAPPED_PATHWAY") == panel.get(STRONG_ANCHOR_REACTION), "proximal_set_member": panel.get(rid, "UNMAPPED_PATHWAY") in proximal})]
            selected = np.asarray([[rrows[i][f"{metric}_median"] for metric in REACTION_STATS] for i in positions], dtype=float)
            row = {"scope": scope, "reaction_count": len(positions), "feature_count": 400 * len(positions)}
            for metric in REACTION_STATS:
                s = qstats(selected[:, stat_pos[metric]].tolist())
                for k, v in s.items(): row[f"{metric}_{k}"] = v
            sensitivity.append(row)
        sfields = list(sensitivity[0])
        write_tsv(tmp / "BRIDGEPL1_SENSITIVITY_SUMMARY.tsv", sfields, sensitivity)
        manifest = {"schema": CORE_SCHEMA, "status": COMPLETE_STATUS, "pl2_status": "NOT_RUN", "row_counts": {"evaluation_registry": 400, "reaction_evaluation_features": 1672000, "context_reaction_summary": 66880, "reaction_summary": 4180, "pathway_summary": len(prows)}, "source_sha256": hashes, "analysis_contract": analysis_contract(), "reaction_inventory": {"cache_reactions": 4181, "primary_reactions": 4180, "hex1_index": hex_idx, "panel_reactions": 3893, "unmapped_reactions": 288}, "execution": {"reaction_block_size": REACTION_BLOCK_SIZE, "xz_preset": XZ_PRESET, "feature_storage": "partitioned_tsv_v1", "scalar_authority": "target_landscape_metrics", "scalar_parity_rtol": SCALAR_PARITY_RTOL, "scalar_parity_atol": SCALAR_PARITY_ATOL}, "publication": "PL1_PREDICTABILITY_LANDSCAPE_COMPLETE"}
        artifact_files = ["BRIDGEPL1_SOURCE_AUDIT.json", "BRIDGEPL1_EVALUATION_REGISTRY.tsv", storage.FEATURE_MANIFEST] + [p["filename"] for p in feature_manifest["parts"]] + ["BRIDGEPL1_CONTEXT_REACTION_SUMMARY.tsv.xz", "BRIDGEPL1_REACTION_SUMMARY.tsv", "BRIDGEPL1_PATHWAY_SUMMARY.tsv", "BRIDGEPL1_SENSITIVITY_SUMMARY.tsv"]
        (tmp / "BRIDGEPL1_SOURCE_AUDIT.json").write_bytes(json.dumps({"schema": "bridge.pl1.source_audit.v1", "status": COMPLETE_STATUS, "source_sha256": hashes, "counts": manifest["row_counts"]}, sort_keys=True, indent=2, ensure_ascii=True).encode("ascii") + b"\n")
        manifest["artifact_sha256"] = {name: sha256(tmp / name) for name in artifact_files}
        (tmp / "BRIDGEPL1_MANIFEST.json").write_bytes(json.dumps(manifest, sort_keys=True, indent=2, ensure_ascii=True).encode("ascii") + b"\n")
        if source_hashes() != hashes:
            raise RuntimeError("source changed during PL1 build")
        existing = list(out.iterdir()) if out.exists() else []
        if existing:
            names = {p.name for p in existing}
            new_names = {p.name for p in tmp.iterdir()}
            if names != new_names or any((out / n).read_bytes() != (tmp / n).read_bytes() for n in names):
                raise RuntimeError("refusing to overwrite differing PL1 outputs")
            shutil.rmtree(tmp)
            for p in work_files:
                p.unlink(missing_ok=True)
            return manifest
        if out.exists():
            out.rmdir()
        tmp.rename(out)
        for p in work_files:
            p.unlink(missing_ok=True)
        return manifest
    except Exception:
        shutil.rmtree(tmp, ignore_errors=True)
        for p in work_files:
            p.unlink(missing_ok=True)
        raise


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=OUT)
    args = ap.parse_args()
    print(json.dumps(build(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
