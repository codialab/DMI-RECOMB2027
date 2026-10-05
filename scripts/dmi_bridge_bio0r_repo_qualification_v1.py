#!/usr/bin/env python3
"""Repository-bound, outcome-free DMI-BRIDGE-BIO-0R qualification.

This adapter independently reproduces provenance-critical weight semantics on
the frozen current panel.  It deliberately does not import GEM, solver,
sampling, or concordance/scoring implementations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import stat
import sys
import tarfile
import tempfile
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import dmi_bridge_bio0r_qualification_core_v1 as core
from dmi_bridge_bio0_audit_core_v1 import WeightRow, distribution_fingerprint


ROOT = Path(__file__).resolve().parents[1]
OUT_REL = "outputs/dmi_bridge_bio0r_qualification_v1"
PATCH_REL = "dmi_bridge_bio0r_qualification_v1.patch"
PATCH_SHA256 = "9aa50d571f7fbea4a74f157c721ce794e7cf6aaa5e7c85b5f9ea1d9d97392ac2"
PANEL_ID = "9c27c8e75aa128a4a4cb15547f23557845f37d35a7323e875af5a1711a527f42"
NUMERIC_ATOL = 1e-14
NUMERIC_RTOL = 1e-12

BIO0_REL = Path("outputs/dmi_bridge_bio0_audit_v1")
BIO0_FILES = {
    "BRIDGEBIO0_SOURCE_AUDIT.json": "acf214cf3bbdb9a0737673e7808696da83f7b52cef599d43f660c21ba5b59787",
    "BRIDGEBIO0_OPERATOR_REGISTRY.tsv": "c843e7289910665b2a2e54976fbc7c30294c570fe089d3946d8e2019cdb2fec0",
    "BRIDGEBIO0_GENE_SCORE_PROVENANCE.tsv": "b7b743ee087ca3166b96771b6ff883c1bcfa30680a1b6fe5817927a81ecb0100",
    "BRIDGEBIO0_REACTION_PANEL.tsv": "f6e881ef949b2a519f0ac731d0dd8b08dba0d44f40e24ca94e18a47017907f4d",
    "BRIDGEBIO0_SUPPORT_AUDIT.tsv": "340c6ef821a093d4d7949b6488cdba2fa6cf9c3300b54b9a8b057016339de69e",
    "BRIDGEBIO0_ANALYSIS_CONTRACT.json": "5e991aa58275e1545628bd9b3e9014b20c1fe799bfcc9810a0d6d712a28c8e2a",
    "BRIDGEBIO0_FOUR_ARM_WEIGHTS.tsv.xz": "843e17b7012a6fa212d1e9ac9e1dc6e309e854d94ab5c4b3f9e316b589768f4f",
    "BRIDGEBIO0_MANIFEST.json": "40df8378876820bb616fdbdc8ef8a0f77d43e36f464104fca6cae9a16b39d7c1",
}

EXTRA_SOURCES = {
    "bio0r_patch": (PATCH_REL, PATCH_SHA256),
    "sa_preflight": ("13_recomb_identifiability/outputs/sa_u1_u2_followup_20260926/PREFLIGHT.json", "dc4450e3994b2d08ccfc15c6ba218552db2b5bd9be9c26298414c8d9ce36c2cc"),
    "sa_runner": ("13_recomb_identifiability/scripts/run_sa_u1_u2_followup.py", "70823e2c2490d40732e47c5d68078afd10dae33aca6ddb5b016332a59bef8ade"),
    "sa_core": ("13_recomb_identifiability/scripts/sa_u1_u2_core.py", "971f7530eee3eb1bc268f5f1d8a67e2329469a29a2376389c66496a42d612349"),
    "candidate_builder": ("outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/build_current_candidates.py", "c7d566dd612e2957cbe0d758e6cc257738b062fcf8cd82bde42ba21941b25893"),
    "cache_reconstructor": ("outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstruct_flux_cache.py", "70c77a76a863a7ca800b7c3610175417b2720c5cbc3d4c1b6db331ec3cef1d32"),
    "reproduction_runner": ("outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reproduction/run_six_arm_reproduction.py", "c013b4aba776290558dc2ed1d4140a7d43cf4a526ab6c8da8d6feb32b70df373"),
    "historical_analysis_manifest": ("outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/results_current/analysis_manifest.json", "12dcc9bf71885b402f510d989baf2d5204fb9cafd2de8c83ec41c29b4b8c1459"),
    "coordinate_verification": ("outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/results_current/coordinate_reconstruction_verification.json", "64a915ad0899c8e03d1b59cc2e717bc2eb971a21efbffe179b902fcbefdcb5f2"),
    "corrected_archive": ("13_recomb_identifiability/reference_inputs/sa_prime49_corrected_20260926/DMI_SA_six_arm_corrected_PRIME49_20260926.tar.xz", "bcf96da38f6fe38ba652b22b1976dc564468dbf2006c9025449a158295a6f806"),
    "rerun_archive": ("13_recomb_identifiability/reference_inputs/sa_prime49_corrected_20260926/DMI_SA_stage12_exact_rerun_ready_20260926.tar.xz", "3bb7bcd90f59500fb29fea70e6783e5b3226a6de92ac39039a3e3cd40a9bef62"),
}

ARCHIVE_MEMBERS = {
    "corrected_archive": {
        "./corrected_lambda_prime49.py": "80d76b52d0568358e4cac2654d5e8c920920797686c2a25d29b97ff6d6f6393d",
        "./analysis_manifest.json": "49b20601efd7f806a8b3bc09d69ab9d626b5972d2df81333386b50f1d27c2d7b",
        "./stage12_calibration.json": "5af467aaaf1403e76c79fd4aab247b85d5426d6158019974a3cb7139432746be",
        "./rerun_six_arm_stage12_exact.py": "0f22fdc623a95b598b6eb22343a2f06fe79b9eea610af3804ac31670f0a5718d",
    },
    "rerun_archive": {
        "corrected_rerun_audit.json": "702865890811e4d8e0295f38d22da0ececfa54a8e01a0db642bc0d3f3ffa878a",
        "rerun_six_arm_stage12_exact.py": "0f22fdc623a95b598b6eb22343a2f06fe79b9eea610af3804ac31670f0a5718d",
        "expected_gate_output.txt": "a63b7eed04f7ae66817c066043325f686869e4d2cd37afcf47045ff03f448a9e",
    },
}

EXPECTED_SIGN_HASHES = {
    20260926: "3c46cfd6c63da4d229e76ae88530bd961f34466542085a1f77d18ccf10ace675",
    20260927: "75a15a25c9a8f16415bfd138b2abdbd8556849e919486922a67c7013c47e480e",
    20260928: "ce5813fd0d9d797c3a6e130466c73d3a2b6fd6291cc6d8f82b01d4ee2952ffdc",
    20260929: "2faad87595edbb61eff0f49179c5aa84b3f0640ac41f2d24918c04d414bc5be1",
    20260930: "5c7000b4f333f01e6f018547f948273084c3adaf0dc0b9c9017a6b545629fd90",
    20260931: "655a17e35286c6cb1ee71d2ef30ba6d235c638e84abe393b912da98eefdc1341",
    20260932: "d6393b44c06a2687a5872696590b03977b493bce6aad342c79e3d9306ba68b88",
    20260933: "26001bb542362c1b9b501f3669b79a0388b002d5d506ec6d1accb0185d739252",
    20260934: "49b5e09e5198bbdfb31ab72483351a3400b470de751b375c7594feb240f00cbf",
    20260935: "dec47d12434cebfae0c065ab586b714078e05d1b0d5d7cd40f72c965e3f679ce",
    20260936: "783770d6195ea82c18d45670765dd353a5a2d25a58ddd3f0750f51c994444810",
    20260937: "a1017e8917edf001bdcc69a24461ab8190f033765bc2a32a1dad59a17b58a0d2",
    20260938: "24131184c4d48a4734488110cb27f2e81b8948e6bf7d616be9078a2d2246abf9",
    20260939: "5f2e6b8a615d883e7037d7353ebc3afdc4cb05d6ff82cb6165e5f757b7ce8498",
    20260940: "3583588a402212726b1d0a018c78810d556e6291429a1a8e9d27e513c000f4d7",
    20260941: "8ffc3d2eb65caf22bd1361844b4462831e58765bdc5c28c047af66584ce68f95",
}
ORDER_HASH = "13f7cfb15a446e0b1e61d75bf1152ab51737c6567c3a5b55f92c8728e953c4de"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def secure_file(repo: Path, path: Path) -> Path:
    candidate = path if path.is_absolute() else repo / path
    resolved = candidate.resolve(strict=True)
    resolved.relative_to(repo)
    mode = os.lstat(candidate).st_mode
    require(stat.S_ISREG(mode) and not stat.S_ISLNK(mode), f"not a regular non-symlink file: {candidate}")
    return resolved


def verified_source(repo: Path, relative_path: str, expected: str) -> dict[str, Any]:
    path = secure_file(repo, Path(relative_path))
    observed = sha256_file(path)
    require(observed == expected, f"source hash mismatch: {relative_path}: {observed} != {expected}")
    return {"relative_path": relative_path, "sha256": observed, "bytes": path.stat().st_size}


def canonical_hash(obj: object) -> str:
    return core.sha256_bytes(core.canonical_json_bytes(obj))


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"expected JSON object: {path}")
    return value


def source_audit(repo: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    frozen: dict[str, Any] = {}
    for name, expected in BIO0_FILES.items():
        frozen[name] = verified_source(repo, str(BIO0_REL / name), expected)

    manifest = read_json(repo / BIO0_REL / "BRIDGEBIO0_MANIFEST.json")
    require(manifest.get("status") == "BIO0_READY_FOR_BIO1_FOUR_ARM_DECOMPOSITION", "BIO-0 status mismatch")
    require(manifest.get("row_counts") == {
        "operator_registry": 13,
        "gene_score_provenance": 5,
        "reaction_panel": 3893,
        "support_audit": 160,
        "four_arm_weights": 12800,
    }, "BIO-0 row accounting mismatch")

    frozen_source_doc = read_json(repo / BIO0_REL / "BRIDGEBIO0_SOURCE_AUDIT.json")
    require(frozen_source_doc.get("candidate_panel", {}).get("id") == PANEL_ID, "BIO-0 panel ID mismatch")
    upstream = {}
    declared = frozen_source_doc.get("sources")
    require(isinstance(declared, dict) and len(declared) == 25, "expected 25 BIO-0 source records")
    for role, item in sorted(declared.items()):
        require(isinstance(item, dict) and set(item) == {"relative_path", "sha256"}, f"bad BIO-0 source record: {role}")
        upstream[role] = verified_source(repo, str(item["relative_path"]), str(item["sha256"]))

    extra = {role: verified_source(repo, rel, expected) for role, (rel, expected) in EXTRA_SOURCES.items()}
    members: dict[str, Any] = {}
    for archive_role, expected_members in ARCHIVE_MEMBERS.items():
        archive_path = repo / EXTRA_SOURCES[archive_role][0]
        member_records = {}
        with tarfile.open(archive_path) as archive:
            for member, expected in expected_members.items():
                handle = archive.extractfile(member)
                require(handle is not None, f"archive member missing: {archive_role}:{member}")
                payload = handle.read()
                observed = hashlib.sha256(payload).hexdigest()
                require(observed == expected, f"archive member hash mismatch: {archive_role}:{member}")
                member_records[member] = {"sha256": observed, "bytes": len(payload)}
        members[archive_role] = member_records
    return {
        "schema": "bridge.bio0r.source_audit.v1",
        "frozen_bio0_outputs": frozen,
        "frozen_bio0_sources_reverified": upstream,
        "additional_sources": extra,
        "archive_members": members,
        "bio0_status": manifest["status"],
        "bio0_row_counts": manifest["row_counts"],
        "operator_reproduction_only": True,
        "new_scientific_weight_generation": False,
        "concordance_outcome_accessed": False,
    }, frozen_source_doc


def require_columns(frame: pd.DataFrame, required: Iterable[str], name: str) -> None:
    missing = set(required) - set(frame.columns)
    require(not missing, f"{name} missing columns: {sorted(missing)}")


def panel_fingerprint(candidate: pd.DataFrame, source_doc: dict[str, Any]) -> str:
    fields = ["algorithm", "tumor", "ensemble_hash", "projection_hash", "rna_context_key", "rna_training_samples", "sample_index"]
    payload = {
        "schema": "bridge3a.candidate_panel.v1",
        "identity_fields": fields,
        "source_sha256": {
            "stage12_candidates": source_doc["sources"]["candidate_panel"]["sha256"],
            "stage12_provenance": source_doc["sources"]["candidate_provenance"]["sha256"],
        },
        "candidates": candidate[fields].astype(str).sort_values(fields, kind="stable").to_dict("records"),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()).hexdigest()


def identity_overlap_audit(candidate: pd.DataFrame, historical: pd.DataFrame, current: pd.DataFrame) -> dict[str, Any]:
    fields = ["algorithm", "tumor", "ensemble_hash", "projection_hash", "rna_context_key", "sample_index"]
    require_columns(historical, [*fields, "observed_targets"], "historical Vmax table")
    require_columns(current, fields, "current Vmax table")
    hist_vmax = historical.loc[historical["observed_targets"].astype(str).eq("Vmax")]
    old = set(map(tuple, hist_vmax[fields].astype(str).drop_duplicates().to_numpy()))
    new = set(map(tuple, current[fields].astype(str).drop_duplicates().to_numpy()))
    overlap = old & new
    old_only = old - new
    new_only = new - old
    by_algorithm = lambda values: {a: sum(row[0] == a for row in values) for a in sorted({row[0] for row in values})}
    result = {
        "identity_fields": fields,
        "historical_total_rows": int(len(historical)),
        "historical_observed_target_counts": {str(k): int(v) for k, v in historical["observed_targets"].value_counts().sort_index().items()},
        "historical_vmax_rows": int(len(hist_vmax)),
        "historical_unique": len(old),
        "current_rows": int(len(current)),
        "current_unique": len(new),
        "overlap": len(overlap),
        "historical_only": len(old_only),
        "current_only": len(new_only),
        "overlap_by_algorithm": by_algorithm(overlap),
        "historical_only_by_algorithm": by_algorithm(old_only),
        "current_only_by_algorithm": by_algorithm(new_only),
        "identity_remapping_used": False,
    }
    require(result["historical_observed_target_counts"] == {"Vglx": 3200, "Vlac": 3200, "Vmax": 3200}, "historical observed-target accounting mismatch")
    require((result["historical_unique"], result["current_unique"], result["overlap"], result["historical_only"], result["current_only"]) == (640, 640, 480, 160, 160), "historical/current overlap mismatch")
    require(result["overlap_by_algorithm"] == {"CORDA": 160, "GIMME": 160, "iMAT": 160}, "non-RIPTiDe overlap mismatch")
    require(result["historical_only_by_algorithm"] == {"RIPTiDe": 160} and result["current_only_by_algorithm"] == {"RIPTiDe": 160}, "identity turnover is not exclusively RIPTiDe")
    panel_ids = set(map(tuple, candidate[fields].astype(str).drop_duplicates().to_numpy()))
    require(panel_ids == new, "current Vmax identities differ from current candidate panel")
    return result


def solve_kernel(distance: np.ndarray, target_ess: float = 20.0) -> tuple[np.ndarray, float]:
    distance = np.asarray(distance, float)
    distance -= float(distance.min())

    def make(temperature: float) -> np.ndarray:
        z = -distance / temperature
        z -= float(z.max())
        return core.normalize(np.exp(z))

    lo, hi = 1e-12, 1.0
    while core.effective_sample_size(make(hi)) < target_ess and hi < 1e6:
        hi *= 10.0
    require(core.effective_sample_size(make(hi)) >= target_ess, "failed to bracket ESS-20 temperature")
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        if core.effective_sample_size(make(mid)) < target_ess:
            lo = mid
        else:
            hi = mid
    return make(hi), hi


def reproduce_coordinates_and_vmax(candidate: pd.DataFrame, observables: pd.DataFrame, current: pd.DataFrame, dmi: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    key = ["algorithm", "tumor", "ensemble_hash", "sample_index"]
    require_columns(candidate, [*key, "projection_hash", "rna_context_key", "EX_glc__D_e", "LDH_L", "lac_production"], "candidate table")
    require_columns(observables, [*key, "Vmax", "Vlac", "q_Vmax", "q_Vlac"], "candidate observables")
    require_columns(current, ["mouse_id", *key, "projection_hash", "rna_context_key", "Vmax_coordinate", "q_Vmax", "weight", "target_ess", "achieved_ess"], "current Vmax weights")
    require_columns(dmi, ["mouse_id", "tumor_type", "Vmax", "Vlac"], "DMI table")
    require(len(candidate) == 640 and not candidate.duplicated(key).any(), "candidate table must contain 640 unique keys")
    require(dmi.groupby("tumor_type").mouse_id.nunique().to_dict() == {"CT2A": 5, "GL261": 5}, "DMI mouse accounting mismatch")

    coords = candidate[key + ["projection_hash", "rna_context_key", "EX_glc__D_e", "LDH_L", "lac_production"]].merge(observables, on=key, validate="one_to_one")
    coords["Vmax_recomputed"] = np.maximum(-coords["EX_glc__D_e"].to_numpy(float), 0.0)
    coords["lac_production_recomputed"] = np.maximum(-coords["LDH_L"].to_numpy(float), 0.0)
    coords["q_Vmax_recomputed"] = np.nan
    coords["q_Vlac_recomputed"] = np.nan
    for _, indices in coords.groupby("tumor", sort=True).groups.items():
        require(len(indices) == 320, "each tumor must have 320 candidates")
        coords.loc[indices, "q_Vmax_recomputed"] = core.candidate_midrank_coordinate(coords.loc[indices, "Vmax_recomputed"].to_numpy(float))
        coords.loc[indices, "q_Vlac_recomputed"] = core.candidate_midrank_coordinate(coords.loc[indices, "lac_production_recomputed"].to_numpy(float))

    raw_vmax_diff = float(np.max(np.abs(coords["Vmax_recomputed"] - coords["Vmax"])))
    raw_lac_diff = float(np.max(np.abs(coords["lac_production_recomputed"] - coords["lac_production"])))
    observable_lac_diff = float(np.max(np.abs(coords["lac_production_recomputed"] - coords["Vlac"])))
    qv_diff = float(np.max(np.abs(coords["q_Vmax_recomputed"] - coords["q_Vmax"])))
    ql_diff = float(np.max(np.abs(coords["q_Vlac_recomputed"] - coords["q_Vlac"])))
    require(raw_vmax_diff == raw_lac_diff == observable_lac_diff == qv_diff == ql_diff == 0.0, "raw/rank coordinate reproduction mismatch")

    mouse_q: dict[str, float] = {}
    for _, group in dmi.groupby("tumor_type", sort=True):
        ranks = core.mouse_rank_coordinate(group["Vmax"].to_numpy(float))
        mouse_q.update({str(mouse): float(value) for mouse, value in zip(group["mouse_id"], ranks)})

    joined = current.merge(coords[key + ["q_Vmax_recomputed", "q_Vlac_recomputed"]], on=key, validate="many_to_one")
    rows = []
    temperatures = {}
    max_weight_diff = max_ess_diff = max_coordinate_diff = max_q_weight_diff = 0.0
    for mouse, group in joined.groupby("mouse_id", sort=True):
        require(len(group) == 320 and group["tumor"].nunique() == 1, f"bad current Vmax group: {mouse}")
        weight, temperature = solve_kernel(np.square(group["q_Vmax_recomputed"].to_numpy(float) - mouse_q[str(mouse)]))
        achieved = core.effective_sample_size(weight)
        temperatures[str(mouse)] = float(temperature)
        for i, row in enumerate(group.itertuples(index=False)):
            saved_weight = float(row.weight)
            rec = row._asdict()
            rec.update({
                "mouse_Vmax_coordinate": mouse_q[str(mouse)],
                "kernel_temperature": float(temperature),
                "recomputed_weight": float(weight[i]),
                "recomputed_ess": achieved,
                "coordinate_abs_diff": abs(float(row.Vmax_coordinate) - float(row.q_Vmax_recomputed)),
                "weight_abs_diff": abs(saved_weight - float(weight[i])),
                "ess_abs_diff": abs(float(row.achieved_ess) - achieved),
                "q_Vmax_is_weight": float(row.q_Vmax) == saved_weight,
            })
            rows.append(rec)
        max_coordinate_diff = max(max_coordinate_diff, float(np.max(np.abs(group["Vmax_coordinate"] - group["q_Vmax_recomputed"]))))
        max_q_weight_diff = max(max_q_weight_diff, float(np.max(np.abs(group["q_Vmax"] - group["weight"]))))
        max_weight_diff = max(max_weight_diff, float(np.max(np.abs(group["weight"].to_numpy(float) - weight))))
        max_ess_diff = max(max_ess_diff, abs(achieved - 20.0), float(np.max(np.abs(group["achieved_ess"].to_numpy(float) - achieved))))
    require(max_coordinate_diff == 0.0 and max_q_weight_diff == 0.0, "saved Vmax field semantics mismatch")
    require(max_weight_diff <= NUMERIC_ATOL and max_ess_diff <= NUMERIC_ATOL, "Vmax kernel reproduction exceeds frozen tolerance")

    lactate = coords[[*key, "projection_hash", "rna_context_key", "LDH_L", "lac_production", "lac_production_recomputed", "Vlac", "q_Vlac", "q_Vlac_recomputed"]].copy()
    lactate["direction"] = lactate["tumor"].map({"CT2A": 1.0, "GL261": -1.0})
    require(lactate["direction"].notna().all(), "unknown tumor direction")
    lactate["historical_score"] = lactate["direction"] * (2.0 * lactate["q_Vlac_recomputed"] - 1.0)
    summary = {
        "raw_vmax_max_abs_diff": raw_vmax_diff,
        "raw_lactate_max_abs_diff": raw_lac_diff,
        "observable_lactate_max_abs_diff": observable_lac_diff,
        "candidate_vmax_coordinate_max_abs_diff": qv_diff,
        "candidate_lactate_coordinate_max_abs_diff": ql_diff,
        "saved_vmax_coordinate_max_abs_diff": max_coordinate_diff,
        "saved_q_vmax_vs_weight_max_abs_diff": max_q_weight_diff,
        "weight_max_abs_diff": max_weight_diff,
        "ess_max_abs_diff": max_ess_diff,
        "temperatures": temperatures,
        "numeric_tolerance": {"atol": NUMERIC_ATOL, "rtol": NUMERIC_RTOL, "source": "REPRODUCTION_COMPARISON.json"},
    }
    return pd.DataFrame(rows), lactate, summary


def preserve_stratum_mass_tilt(base: np.ndarray, score: np.ndarray, strata: np.ndarray, lam: float) -> np.ndarray:
    base = core.normalize(base)
    score = np.asarray(score, float)
    strata = np.asarray(strata, object)
    require(len(base) == len(score) == len(strata), "stratum tilt length mismatch")
    out = np.zeros_like(base)
    for label in pd.unique(strata):
        mask = strata == label
        mass = float(base[mask].sum())
        if mass == 0.0:
            continue
        conditional = base[mask] / mass
        out[mask] = mass * core.tilted_weights(conditional, score[mask], lam)
    return core.normalize(out)


def historical_global_tilt(base: np.ndarray, score: np.ndarray, lam: float) -> np.ndarray:
    """Reproduce the archived/BIO-0 global tilt without pre-normalizing q.

    The saved Vmax vector is already normalized.  The producing code multiplied
    that serialized vector directly and normalized only the result; retaining
    that operation order matters for deterministic float fingerprints.
    """
    base = np.asarray(base, float)
    score = np.asarray(score, float)
    require(base.ndim == score.ndim == 1 and len(base) == len(score), "global tilt length mismatch")
    require(np.isfinite(base).all() and np.isfinite(score).all() and (base >= 0.0).all() and float(base.sum()) > 0.0, "invalid global tilt inputs")
    z = float(lam) * score
    z -= float(z.max())
    updated = base * np.exp(z)
    return updated / float(updated.sum())


def random_registry(candidate: pd.DataFrame) -> tuple[pd.DataFrame, dict[int, np.ndarray]]:
    fields = ["algorithm", "tumor", "ensemble_hash", "projection_hash", "rna_context_key", "sample_index"]
    payload = {"schema": "bridge.bio0r.random_candidate_order.v1", "identity_fields": fields, "candidates": candidate[fields].astype(str).to_dict("records")}
    order_hash = canonical_hash(payload)
    require(order_hash == ORDER_HASH, "historical random candidate ordering mismatch")
    rows, signs_by_seed = [], {}
    for index, seed in enumerate(core.RANDOM_SEEDS):
        signs = core.random_sign_vector(seed, len(candidate))
        fingerprint = core.random_sign_fingerprint(signs)
        require(fingerprint == EXPECTED_SIGN_HASHES[seed], f"random sign fingerprint mismatch for seed {seed}")
        signs_by_seed[seed] = signs
        rows.append({
            "control_index": index,
            "seed": seed,
            "candidate_order_sha256": order_hash,
            "sign_vector_sha256": fingerprint,
            "minus_count": int(np.count_nonzero(signs == -1)),
            "plus_count": int(np.count_nonzero(signs == 1)),
            "lambda": 0.5,
            "score_formula": "sign[candidate_file_position]*tumor_direction",
            "normalization": "within_algorithm_rna_context_mass_preserving",
        })
    return pd.DataFrame(rows), signs_by_seed


def candidate_id(row: Any) -> str:
    return "|".join(map(str, (row.algorithm, row.tumor, row.ensemble_hash, int(row.sample_index))))


def weight_rows(group: pd.DataFrame, column: str) -> list[WeightRow]:
    return [WeightRow(candidate_id(row), str(row.stratum_id), float(getattr(row, column))) for row in group.itertuples(index=False)]


def reconstruct_bio0(candidate: pd.DataFrame, current: pd.DataFrame, lactate: pd.DataFrame, frozen: pd.DataFrame, signs_by_seed: dict[int, np.ndarray]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    key = ["algorithm", "tumor", "ensemble_hash", "sample_index"]
    current_lac = current.merge(lactate[key + ["q_Vlac_recomputed"]], on=key, validate="many_to_one")
    position = {tuple(row): i for i, row in enumerate(candidate[key].astype(str).to_numpy())}
    require(len(position) == 640, "candidate-position map is not one-to-one")
    reproduction_rows = []
    support_rows = []
    rec_fingerprints, frozen_fingerprints = {}, {}
    max_baseline_diff = max_full_diff = 0.0
    max_random_mass_diff = 0.0

    frozen_base = frozen.loc[frozen["arm"].eq("strong_anchor_baseline")]
    frozen_full = frozen.loc[frozen["arm"].eq("full_weak_update")]
    require(len(frozen_base) == len(frozen_full) == 3200, "frozen BIO-0 arm row accounting mismatch")

    rec_groups = []
    for mouse, group0 in current_lac.groupby("mouse_id", sort=True):
        group = group0.copy()
        tumor = str(group["tumor"].iloc[0])
        direction = 1.0 if tumor == "CT2A" else -1.0
        group["stratum_id"] = group["algorithm"].astype(str) + "::" + group["rna_context_key"].astype(str)
        base = group["weight"].to_numpy(float)
        score = direction * (2.0 * group["q_Vlac_recomputed"].to_numpy(float) - 1.0)
        full = historical_global_tilt(base, score, 0.5)
        group["baseline_recomputed"] = base
        group["full_recomputed"] = full
        rec_groups.append(group)

        variants = []
        for lam in (0.25, 0.5, 1.0):
            variants.append(("correct_lactate_direction", lam, historical_global_tilt(base, score, lam)))
            variants.append(("reversed_lactate_direction", lam, historical_global_tilt(base, -score, lam)))
        candidate_positions = np.asarray([position[tuple(values)] for values in group[key].astype(str).to_numpy()], dtype=int)
        strata = group["stratum_id"].to_numpy(object)
        for index, seed in enumerate(core.RANDOM_SEEDS):
            random_score = signs_by_seed[seed][candidate_positions].astype(float) * direction
            variants.append((f"random_direction_{index:02d}", 0.5, preserve_stratum_mass_tilt(base, random_score, strata, 0.5)))

        for operator, lam, updated in variants:
            core.assert_matched_support(base, updated, strata)
            for label in sorted(set(strata)):
                mask = strata == label
                baseline_mass = float(base[mask].sum())
                updated_mass = float(updated[mask].sum())
                baseline_support = int(np.count_nonzero(base[mask] > 0.0))
                updated_support = int(np.count_nonzero(updated[mask] > 0.0))
                if operator.startswith("random_direction"):
                    max_random_mass_diff = max(max_random_mass_diff, abs(baseline_mass - updated_mass))
                support_rows.append({
                    "mouse_id": mouse,
                    "tumor": tumor,
                    "stratum_id": label,
                    "operator": operator,
                    "lambda": lam,
                    "candidate_rows": int(mask.sum()),
                    "baseline_mass": baseline_mass,
                    "updated_mass": updated_mass,
                    "baseline_literal_zero": baseline_mass == 0.0,
                    "updated_literal_zero": updated_mass == 0.0,
                    "baseline_positive_candidate_support": baseline_support,
                    "updated_positive_candidate_support": updated_support,
                    "support_differs": baseline_support != updated_support or not np.array_equal(base[mask] > 0.0, updated[mask] > 0.0),
                    "support_invented": bool(np.any((base[mask] == 0.0) & (updated[mask] > 0.0))),
                    "support_lost": bool(np.any((base[mask] > 0.0) & (updated[mask] == 0.0))),
                })

    reconstructed = pd.concat(rec_groups, ignore_index=True)
    join_fields = ["mouse_id", "tumor", "algorithm", "ensemble_hash", "sample_index", "rna_context_key", "stratum_id"]
    comparison = reconstructed.merge(frozen_base[join_fields + ["weight"]].rename(columns={"weight": "baseline_frozen"}), on=join_fields, validate="one_to_one")
    comparison = comparison.merge(frozen_full[join_fields + ["weight"]].rename(columns={"weight": "full_frozen"}), on=join_fields, validate="one_to_one")
    comparison["baseline_abs_diff"] = np.abs(comparison["baseline_recomputed"] - comparison["baseline_frozen"])
    comparison["full_abs_diff"] = np.abs(comparison["full_recomputed"] - comparison["full_frozen"])
    max_baseline_diff = float(comparison["baseline_abs_diff"].max())
    max_full_diff = float(comparison["full_abs_diff"].max())
    require(np.allclose(comparison["baseline_recomputed"], comparison["baseline_frozen"], atol=NUMERIC_ATOL, rtol=NUMERIC_RTOL), "BIO-0 baseline reproduction mismatch")
    require(np.allclose(comparison["full_recomputed"], comparison["full_frozen"], atol=NUMERIC_ATOL, rtol=NUMERIC_RTOL), "BIO-0 full reproduction mismatch")

    for mouse, group in comparison.groupby("mouse_id", sort=True):
        rec_fingerprints[str(mouse)] = {
            "strong_anchor_baseline": distribution_fingerprint(weight_rows(group, "baseline_recomputed")),
            "full_weak_update": distribution_fingerprint(weight_rows(group, "full_recomputed")),
        }
        frozen_fingerprints[str(mouse)] = {
            "strong_anchor_baseline": distribution_fingerprint(weight_rows(group, "baseline_frozen")),
            "full_weak_update": distribution_fingerprint(weight_rows(group, "full_frozen")),
        }

    global_fp = {}
    for label, fps in (("reconstructed", rec_fingerprints), ("frozen", frozen_fingerprints)):
        global_fp[label] = {}
        for arm in ("strong_anchor_baseline", "full_weak_update"):
            schema = "bridge.bio0r.distribution_fingerprints.v1" if label == "reconstructed" else "bridge.bio0r.frozen_distribution_fingerprints.v1"
            global_fp[label][arm] = canonical_hash({"schema": schema, "arm": arm, "per_mouse": {m: fps[m][arm] for m in sorted(fps)}})

    require(global_fp["reconstructed"] == {
        "strong_anchor_baseline": "31a9acd56a1e065e127af7b986968cd9f26c30d48d5a8dd084cbdd97deba0c67",
        "full_weak_update": "144b6faecbcdd8b5daf73137675942a18705d76c7e878a986d176d0bb21d2c89",
    }, "reconstructed distribution fingerprints changed")
    require(global_fp["frozen"] == {
        "strong_anchor_baseline": "60ad32b694c9f16a2915ca26b14832178600f2fbc3cd81c5ac4cd7755b555b23",
        "full_weak_update": "5b51e46bea3fed97238e23dcdc3b06d51dc92cd4107323030c8504f793082094",
    }, "frozen distribution fingerprints changed")

    support = pd.DataFrame(support_rows)
    require(len(support) == 3520, "operator support row accounting mismatch")
    require(not support[["support_differs", "support_invented", "support_lost", "baseline_literal_zero", "updated_literal_zero"]].any().any(), "operator support qualification failed")
    require(max_random_mass_diff <= 1e-12, "random control failed stratum mass preservation")
    summary = {
        "baseline_max_abs_diff": max_baseline_diff,
        "full_max_abs_diff": max_full_diff,
        "max_random_stratum_mass_diff": max_random_mass_diff,
        "reconstructed_per_mouse_fingerprints": rec_fingerprints,
        "frozen_per_mouse_fingerprints": frozen_fingerprints,
        "global_fingerprints": global_fp,
    }
    keep = join_fields + ["baseline_recomputed", "baseline_frozen", "baseline_abs_diff", "full_recomputed", "full_frozen", "full_abs_diff"]
    return comparison[keep], support, summary


def identity_files_audit(repo: Path, candidate: pd.DataFrame, provenance: dict[str, Any], cache_path: Path) -> dict[str, Any]:
    rows = candidate[["algorithm", "tumor", "ensemble_hash", "projection_hash", "full_parent_path", "projection_path"]].drop_duplicates()
    require(len(rows) == 32, "expected 32 ensemble/projection identity pairs")
    require(set(rows["ensemble_hash"].astype(str)) == set(map(str, provenance.get("stage11_ensemble_hashes", []))), "ensemble identities differ from Stage-12 provenance")
    require(set(rows["projection_hash"].astype(str)) == set(map(str, provenance.get("stage11_projection_hashes", []))), "projection identities differ from Stage-12 provenance")
    records = []
    for row in rows.sort_values(["algorithm", "tumor", "ensemble_hash"], kind="stable").itertuples(index=False):
        full = secure_file(repo, Path(row.full_parent_path))
        projection = secure_file(repo, Path(row.projection_path))
        ensemble_dir = full.parents[1]
        projection_dir = projection.parent
        ensemble_done_path = secure_file(repo, ensemble_dir / "DONE.json")
        ensemble_meta_path = secure_file(repo, ensemble_dir / "metadata.json")
        projection_done_path = secure_file(repo, projection_dir / "DONE.json")
        projection_meta_path = secure_file(repo, projection_dir / "metadata.json")
        ensemble_done = read_json(ensemble_done_path)
        ensemble_meta = read_json(ensemble_meta_path)
        projection_done = read_json(projection_done_path)
        projection_meta = read_json(projection_meta_path)
        require(str(row.ensemble_hash) == ensemble_dir.name == str(ensemble_done.get("ensemble_hash")) == str(ensemble_meta.get("ensemble_hash")), "ensemble identity mismatch")
        require(str(row.projection_hash) == projection_dir.name == str(projection_done.get("projection_hash")) == str(projection_meta.get("projection_hash")), "projection identity mismatch")
        require(str(row.ensemble_hash) == str(projection_done.get("ensemble_hash")) == str(projection_meta.get("ensemble_hash")), "projection-to-ensemble linkage mismatch")
        require(ensemble_done.get("qc_passed") is True and projection_done.get("qc_passed") is True, "ensemble/projection QC marker is not PASS")
        full_sha = sha256_file(full)
        projection_sha = sha256_file(projection)
        require(full_sha == ensemble_done.get("checksums", {}).get("full_parent_fluxes/reaction_samples.parquet.xz"), "full-parent packaged checksum mismatch")
        require(sha256_file(ensemble_meta_path) == ensemble_done.get("checksums", {}).get("metadata.json"), "ensemble metadata checksum mismatch")
        require(projection_sha == projection_done.get("checksums", {}).get("reaction_samples.parquet.xz"), "projection packaged checksum mismatch")
        require(sha256_file(projection_meta_path) == projection_done.get("checksums", {}).get("metadata.json"), "projection metadata checksum mismatch")
        records.append({
            "algorithm": str(row.algorithm), "tumor": str(row.tumor),
            "ensemble_hash": str(row.ensemble_hash), "projection_hash": str(row.projection_hash),
            "full_parent_path": str(full.relative_to(repo)), "full_parent_file_sha256": full_sha,
            "projection_path": str(projection.relative_to(repo)), "projection_file_sha256": projection_sha,
            "ensemble_done_sha256": sha256_file(ensemble_done_path),
            "ensemble_metadata_sha256": sha256_file(ensemble_meta_path),
            "projection_done_sha256": sha256_file(projection_done_path),
            "projection_metadata_sha256": sha256_file(projection_meta_path),
            "projection_source_full_ensemble_checksum": projection_meta.get("source_full_ensemble_checksum"),
            "qc_passed": True,
        })

    cache = np.load(cache_path, allow_pickle=False)
    shapes = {name: list(cache[name].shape) for name in cache.files}
    require(shapes == {"rxns": [4181], "alg": [32], "tumor": [32], "eh": [32], "mats": [32, 20, 4181], "qc": [32, 3]}, "flux-cache shape mismatch")
    cache_keys = {(str(a), str(t), str(e), str(i)) for a, t, e in zip(cache["alg"], cache["tumor"], cache["eh"]) for i in range(20)}
    fields = ["algorithm", "tumor", "ensemble_hash", "sample_index"]
    candidate_keys = set(map(tuple, candidate[fields].astype(str).to_numpy()))
    require(len(cache_keys) == len(candidate_keys) == 640 and cache_keys == candidate_keys, "flux-cache candidate linkage mismatch")
    return {
        "ensemble_projection_records": records,
        "n_ensembles": 32,
        "n_projections": 32,
        "cache_shapes": shapes,
        "cache_key_fields": fields,
        "cache_unique_keys": 640,
        "cache_overlap": 640,
        "cache_only": 0,
        "candidate_only": 0,
        "ordinal_or_hash_remapping_used": False,
    }


def historical_pass_audit(repo: Path, source_doc: dict[str, Any]) -> dict[str, Any]:
    source = source_doc["sources"]
    comparison = read_json(repo / source["reproduction"]["relative_path"])
    require(comparison.get("overall_reproduction_status") == "PASS" and comparison.get("scientific_reproduction_pass") is True, "historical reproduction is not PASS")
    require(comparison.get("reconstructed_cache_sha256") == source["flux_cache"]["sha256"], "historical reproduction cache hash mismatch")
    require(comparison.get("numeric_tolerance") == {"atol": 1e-14, "description": "numpy allclose with equal_nan=True", "rtol": 1e-12}, "historical reproduction tolerance mismatch")
    comparisons = comparison.get("comparisons")
    require(isinstance(comparisons, dict) and len(comparisons) == 11, "expected 11 historical reproduction comparisons")
    for name, record in comparisons.items():
        require(record.get("comparison_status") == "PASS", f"historical comparison not PASS: {name}")
        require(record.get("key_string_equality") is True, f"historical comparison key mismatch: {name}")
        if record.get("numeric_equality_within_tolerance") is not None:
            require(record.get("numeric_equality_within_tolerance") is True, f"historical comparison numerical mismatch: {name}")
        if record.get("historical_row_count") is not None or record.get("reproduced_row_count") is not None:
            require(record.get("historical_row_count") == record.get("reproduced_row_count"), f"historical comparison row mismatch: {name}")

    preflight = read_json(repo / EXTRA_SOURCES["sa_preflight"][0])
    require(preflight.get("status") == "PASS_CANONICAL_CURRENT_PANEL_INPUTS_ONLY_V2", "current preflight is not PASS")
    require(preflight.get("panel") == {"n_candidates": 640, "methods": ["CORDA", "GIMME", "iMAT", "RIPTiDe"], "tumors": ["CT2A", "GL261"], "ensembles_per_method_tumor": 4, "vectors_per_ensemble": 20, "n_ensemble_hashes": 32, "n_projection_hashes": 32}, "current preflight panel mismatch")
    require(preflight.get("historical_reproduction") == {"overall_reproduction_status": "PASS", "scientific_reproduction_pass": True, "comparisons": 11}, "preflight historical reproduction mismatch")

    vmax_manifest = read_json(repo / source["current_vmax_manifest"]["relative_path"])
    require(vmax_manifest == {"operator": "archived corrected_lambda_prime49.py", "target_ess": 20.0, "cache_sha256": source["flux_cache"]["sha256"], "n_rows": 3200}, "current Vmax manifest mismatch")

    contract = read_json(repo / source["sa_contract"]["relative_path"])
    require(contract.get("operator") == "corrected_lambda_prime49.py", "SA operator mismatch")
    require(contract.get("secondary_score") == "direction*(2*q_Vlac-1)" and contract.get("lambda_primary") == 0.5, "SA weak-operator contract mismatch")
    require(contract.get("lambda_grid") == [0.0, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0], "SA lambda grid mismatch")
    require(contract.get("stale_stage12b_weights") == "DO_NOT_USE" and contract.get("solver_jobs") is False, "SA stale/solver contract mismatch")

    baseline = pd.read_csv(repo / source["sa_baseline"]["relative_path"], sep="\t", usecols=["arm", "lambda", "n_overlap", "abs_error", "pass"])
    expected_arms = {"unweighted", "Vmax_only", "LDH_only", "Vmax_plus_weak_LDH", "Vmax_plus_reversed_LDH", "two_quantitative_anchors", "within_stratum_weak", "within_stratum_reversed"}
    require(len(baseline) == 8 and set(baseline["arm"]) == expected_arms and baseline["pass"].eq(True).all(), "baseline reproduction PASS rows mismatch")
    require(baseline["lambda"].eq(0.5).all() and baseline["n_overlap"].eq(49).all(), "baseline reproduction schema mismatch")
    referenced = baseline[baseline["abs_error"].notna()]
    require((referenced["abs_error"] <= 1e-10).all(), "baseline reproduction error exceeds producer gate")

    random_provenance = pd.read_csv(repo / source["sa_random"]["relative_path"], sep="\t", usecols=["control_index", "seed"])
    require(random_provenance["control_index"].tolist() == list(range(16)), "random-control indices mismatch")
    require(random_provenance["seed"].tolist() == list(core.RANDOM_SEEDS), "random-control seeds mismatch")

    corrected_archive = repo / EXTRA_SOURCES["corrected_archive"][0]
    rerun_archive = repo / EXTRA_SOURCES["rerun_archive"][0]
    with tarfile.open(corrected_archive) as archive:
        archived_manifest = json.loads(archive.extractfile("./analysis_manifest.json").read())
        archived_source = archive.extractfile("./corrected_lambda_prime49.py").read().decode("utf-8")
    require(archived_manifest.get("gene_score_id") == "PRIME_SIGNED__COLLECTRI__TF_EXPRESSION_ROBUST_Z__MEDIAN" and archived_manifest.get("n_gene_subsystems") == 49, "archived PRIME contract mismatch")
    require(archived_manifest.get("seed") == 20260918 and archived_manifest.get("bootstrap_replicates") == 10000, "archived bootstrap contract mismatch")
    for token in ("TARGET_ESS=20.0", "LAMBDAS=[0.25,0.5,1.0]", "s=1.0 if tumor=='CT2A' else -1.0", "h=2*ql-1"):
        require(token in archived_source, f"archived operator source lacks exact token: {token}")
    with tarfile.open(rerun_archive) as archive:
        rerun_audit = json.loads(archive.extractfile("corrected_rerun_audit.json").read())
    require(rerun_audit.get("status") == "BLOCKED_ONLY_ON_EXACT_PRIME_VECTOR" and rerun_audit.get("withdrawn_previous_validation") is True, "superseded rerun audit semantics changed")
    require(rerun_audit.get("required_input", {}).get("required_gene_score_id") == "PRIME_SIGNED__COLLECTRI__TF_EXPRESSION_ROBUST_Z__MEDIAN", "superseded rerun audit PRIME requirement mismatch")
    return {
        "reproduction_overall_status": "PASS",
        "scientific_reproduction_pass": True,
        "comparison_count": 11,
        "preflight_status": preflight["status"],
        "baseline_pass_rows": 8,
        "random_provenance_rows": 16,
        "current_vmax_manifest_verified": True,
        "analysis_contract_verified": True,
        "archived_operator_contract_verified": True,
        "superseded_blocked_audit_treated_as_pass": False,
        "outcome_columns_used_for_parameter_selection": False,
    }


def write_tsv(frame: pd.DataFrame, path: Path, sort_fields: list[str]) -> None:
    frame.sort_values(sort_fields, kind="stable").to_csv(path, sep="\t", index=False, lineterminator="\n")


def build(repo: Path, output_dir: Path) -> dict[str, Any]:
    repo = repo.resolve(strict=True)
    require(repo == ROOT.resolve(strict=True), "--repo-root must be the native repository root")
    output_dir = output_dir.resolve()
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    source_report, frozen_source_doc = source_audit(repo)
    sources = frozen_source_doc["sources"]
    candidate = pd.read_csv(repo / sources["candidate_panel"]["relative_path"], sep="\t", compression="xz")
    provenance = read_json(repo / sources["candidate_provenance"]["relative_path"])
    historical = pd.read_csv(repo / sources["historical_vmax"]["relative_path"], sep="\t", compression="xz")
    current = pd.read_csv(repo / sources["current_vmax"]["relative_path"], sep="\t", compression="xz")
    observables = pd.read_csv(repo / sources["observables"]["relative_path"])
    dmi = pd.read_csv(repo / sources["dmi"]["relative_path"])
    frozen_weights = pd.read_csv(repo / BIO0_REL / "BRIDGEBIO0_FOUR_ARM_WEIGHTS.tsv.xz", sep="\t", compression="xz")

    require(panel_fingerprint(candidate, frozen_source_doc) == PANEL_ID, "current candidate panel fingerprint mismatch")
    require(len(candidate) == 640 and set(candidate["algorithm"]) == {"CORDA", "GIMME", "iMAT", "RIPTiDe"}, "current candidate panel membership mismatch")
    identity_overlap = identity_overlap_audit(candidate, historical, current)
    vmax, lactate, vmax_summary = reproduce_coordinates_and_vmax(candidate, observables, current, dmi)
    random_controls, signs = random_registry(candidate)
    bio0_reproduction, support, reproduction_summary = reconstruct_bio0(candidate, current, lactate, frozen_weights, signs)
    identity_files = identity_files_audit(repo, candidate, provenance, repo / sources["flux_cache"]["relative_path"])
    pass_audit = historical_pass_audit(repo, frozen_source_doc)

    identity_report = {
        "schema": "bridge.bio0r.identity_audit.v1",
        "candidate_panel_id": PANEL_ID,
        "candidate_rows": 640,
        "historical_current_vmax": identity_overlap,
        **identity_files,
    }
    source_report["historical_contract_pass_audit"] = pass_audit
    source_report["vmax_reproduction_summary"] = vmax_summary
    source_report["bio0_reproduction_summary"] = reproduction_summary
    source_report["status"] = core.READY_STATUS

    manifest = {
        "schema": core.SCHEMA,
        "status": core.READY_STATUS,
        "operator_reproduction_only": True,
        "new_scientific_weight_generation": False,
        "concordance_outcome_computed": False,
        "candidate_panel_id": PANEL_ID,
        "qualification_gates": {
            "historical_overlap_verified": True,
            "current_vmax_semantics_verified": True,
            "lactate_semantics_verified": True,
            "random_family_verified": True,
            "all_operator_support_verified": True,
            "baseline_full_reproduction_verified": True,
            "ensemble_projection_cache_identity_verified": True,
            "historical_reproduction_contract_verified": True,
        },
        "row_counts": {
            "vmax_reproduction": 3200,
            "lactate_operator_audit": 640,
            "random_control_registry": 16,
            "operator_support_audit": 3520,
            "bio0_reproduction": 3200,
        },
        "outputs": {},
    }

    with tempfile.TemporaryDirectory(prefix="dmi_bridge_bio0r_", dir=str(output_dir.parent)) as temp_name:
        stage = Path(temp_name)
        (stage / "BRIDGEBIO0R_SOURCE_AUDIT.json").write_bytes(core.canonical_json_bytes(source_report))
        (stage / "BRIDGEBIO0R_IDENTITY_AUDIT.json").write_bytes(core.canonical_json_bytes(identity_report))
        write_tsv(vmax, stage / "BRIDGEBIO0R_VMAX_REPRODUCTION.tsv", ["mouse_id", "algorithm", "rna_context_key", "ensemble_hash", "sample_index"])
        write_tsv(lactate, stage / "BRIDGEBIO0R_LACTATE_OPERATOR_AUDIT.tsv", ["tumor", "algorithm", "rna_context_key", "ensemble_hash", "sample_index"])
        write_tsv(random_controls, stage / "BRIDGEBIO0R_RANDOM_CONTROL_REGISTRY.tsv", ["control_index"])
        write_tsv(support, stage / "BRIDGEBIO0R_OPERATOR_SUPPORT_AUDIT.tsv", ["mouse_id", "operator", "lambda", "stratum_id"])
        write_tsv(bio0_reproduction, stage / "BRIDGEBIO0R_BIO0_REPRODUCTION.tsv", ["mouse_id", "algorithm", "rna_context_key", "ensemble_hash", "sample_index"])
        for path in sorted(stage.iterdir()):
            manifest["outputs"][path.name] = {"sha256": sha256_file(path), "bytes": path.stat().st_size}
        (stage / "BRIDGEBIO0R_MANIFEST.json").write_bytes(core.canonical_json_bytes(manifest))

        if output_dir.exists():
            existing = {path.name: path for path in output_dir.iterdir() if path.is_file()}
            staged = {path.name: path for path in stage.iterdir() if path.is_file()}
            require(set(existing) == set(staged), f"refusing differing existing BIO-0R inventory: {output_dir}")
            require(all(existing[name].read_bytes() == staged[name].read_bytes() for name in staged), f"refusing to overwrite differing existing BIO-0R output: {output_dir}")
        else:
            output_dir.mkdir(parents=False)
            for path in stage.iterdir():
                shutil.copy2(path, output_dir / path.name)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path, default=ROOT / OUT_REL)
    args = parser.parse_args()
    manifest = build(args.repo_root, args.output_dir)
    print(json.dumps({"status": manifest["status"], "output_dir": str(args.output_dir.resolve()), "row_counts": manifest["row_counts"]}, sort_keys=True))


if __name__ == "__main__":
    main()
