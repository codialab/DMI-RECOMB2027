#!/usr/bin/env python3
"""Qualify the frozen A2-L and A2-G strong-anchor proposals without outcomes."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import lzma
import os
from pathlib import Path
import stat
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import dmi_bridge_a20_dual_anchor_core_v1 as core
import dmi_bridge_bio0r_repo_qualification_v1 as bio0r


ROOT = Path(__file__).resolve().parents[1]
OUT = Path("outputs/dmi_bridge_a20_dual_anchor_qualification_v1")
PANEL_ID = "9c27c8e75aa128a4a4cb15547f23557845f37d35a7323e875af5a1711a527f42"
TARGET_ESS = 20.0
ATOL = 1e-14
KEY = ["algorithm", "tumor", "ensemble_hash", "projection_hash", "rna_context_key", "sample_index"]
SOURCES = {
    "m1_manifest": ("outputs/dmi_bridge_m1_manuscript_evidence_v1/BRIDGEM1_MANIFEST.json", "201ba9f4fc2f74290d991cd294aec558c5e5dd23282877e871a5c0669a44cdde"),
    "m1_source_audit": ("outputs/dmi_bridge_m1_manuscript_evidence_v1/BRIDGEM1_SOURCE_AUDIT.json", "07c16a60e79969c2a751cf1e7e471fdce88510723cec716b5bc3fe4b100ea073"),
    "candidate_panel": ("12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/candidate_fraction_table.tsv.xz", "421c2c1d196fbd088a4e92e2ecc998ad578aab814a025fff8e1c843597fed1da"),
    "candidate_provenance": ("12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/DMI_FRACTION_ENSEMBLE_PROVENANCE.json", "fd05e5febba47a92005619db38c53b3ba4896d81eb498226d87c61f1e9536168"),
    "ensemble_selection": ("12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/ensemble_selection.tsv", "ffaa8bedba3bc1e4a11f2cdb930de7d0b128936afb88f325125380a5586e8c70"),
    "dmi": ("pre-simulation_constraint_analysis/results/run_20260824T164247Z_corrected_analysis_A/mouse_level_dmi_input.csv", "ae3cbad6d422405780d0dc44e0e31f95770ebc204d1f80f2ccdada12bb119082"),
    "historical_weights": ("12_recomb_method_comparison/outputs_stage12b/dmi_single_anchor_weights.tsv.xz", "b3357587b70329fca7d77b99e4b6c6487feccb54ac7f251d8a0a8437353a8030"),
    "vmax_manifest": ("13_recomb_identifiability/outputs/sa_u1_u2_followup_20260926/current_vmax_base_weights_manifest.json", "2abf9e7c4d3c922bda8c81426001c1f27c6e43f1761e9c9073e0d8ca078bd5bc"),
    "bio0r_manifest": ("outputs/dmi_bridge_bio0r_qualification_v1/BRIDGEBIO0R_MANIFEST.json", "d155b63f59db3479b1dfd2d87c0b9ca470358c0166496abf288a649f04890cea"),
    "bio0r_source_audit": ("outputs/dmi_bridge_bio0r_qualification_v1/BRIDGEBIO0R_SOURCE_AUDIT.json", "33ab67e8d8d12b84cdee1fc43b36bdba6b3680e9451bc01ed071e5fecd1c6fc4"),
    "bio0r_vmax_reproduction": ("outputs/dmi_bridge_bio0r_qualification_v1/BRIDGEBIO0R_VMAX_REPRODUCTION.tsv", "4843397b2f2b777fe8d3a1070812068bd6d434c2fd4888b8024f79d14cccd295"),
    "bio0r_vmax_producer": ("scripts/dmi_bridge_bio0r_repo_qualification_v1.py", "de50fc12ab6382a69d0242beedd799e427fa00e026fcb3e9618f1f2ce66e2b1d"),
    "flux_cache": ("outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstructed/flux_cache.npz", "f9644da4ebd875cf5f0b34cc58b20faf0a0457d7c7a44635e8a010d8896c3dbb"),
    "current_observables": ("outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/current_candidate_observables.csv", "946cd613d673c9aa9367bce705ad5d05562c6a174acc830fbbd36cb395aa5504"),
    "archived_kernel": ("outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstruct_flux_cache.py", "70c77a76a863a7ca800b7c3610175417b2720c5cbc3d4c1b6db331ec3cef1d32"),
    "operator_audit": ("outputs/dmi_bridge3a_registry_v1/BRIDGE3A_OPERATOR_AUDIT.json", "37deb5887ebe3f74976e96fb2a7e2614f952d065d537aeeb74a66cb2fa52f8f8"),
    "bio0_registry": ("outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_OPERATOR_REGISTRY.tsv", "c843e7289910665b2a2e54976fbc7c30294c570fe089d3946d8e2019cdb2fec0"),
    "dmi_mapping": ("analysis/brain_glioma_gem/config/dmi_to_gem_mapping.tsv", "3002fa055a5efaccc27044c3a158351a3f845e5a31c8e8688e34728e6575b489"),
}
VMAX_SHA = "91121bf1251354505cde5f5ae783dc8d14bc9a22ecdf8daa06be1e5d82e3385d"


def require(condition: bool, reason: str) -> None:
    if not condition:
        raise RuntimeError(reason)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(obj: object) -> bytes:
    return (json.dumps(obj, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def safe_source(repo: Path, relative: str, expected: str) -> tuple[Path, int]:
    path = repo / relative
    mode = os.lstat(path).st_mode
    require(stat.S_ISREG(mode) and not stat.S_ISLNK(mode), f"not a regular source: {relative}")
    resolved = path.resolve(strict=True)
    resolved.relative_to(repo)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    require(digest.hexdigest() == expected, f"source SHA256 mismatch: {relative}")
    return resolved, path.stat().st_size


def source_gate(repo: Path) -> tuple[dict[str, Path], pd.DataFrame]:
    actual_root = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=repo, text=True).strip()).resolve()
    require(actual_root == repo, "repository root mismatch")
    paths: dict[str, Path] = {}
    rows: list[dict[str, object]] = []
    for role, (relative, expected) in SOURCES.items():
        paths[role], size = safe_source(repo, relative, expected)
        rows.append(dict(role=role, relative_path=relative, sha256=expected, bytes=size, status="PASS"))

    m1 = json.loads(paths["m1_manifest"].read_text())
    bio = json.loads(paths["bio0r_manifest"].read_text())
    vmax_manifest = json.loads(paths["vmax_manifest"].read_text())
    require(m1.get("status") == "BRIDGE_M1_MANUSCRIPT_EVIDENCE_FROZEN", "M1 is not frozen")
    require(bio.get("status") == "BIO0R_QUALIFIED_READY_FOR_BIO1", "BIO0R is not qualified")
    require(bio.get("candidate_panel_id") == PANEL_ID and bio.get("qualification_gates", {}).get("current_vmax_semantics_verified") is True, "BIO0R panel/operator gate failed")
    require(vmax_manifest == {"operator": "archived corrected_lambda_prime49.py", "target_ess": 20.0, "cache_sha256": SOURCES["flux_cache"][1], "n_rows": 3200}, "current Vmax baseline manifest mismatch")

    receipt = json.loads(paths["bio0r_source_audit"].read_text())
    binding = receipt["frozen_bio0_sources_reverified"]["current_vmax"]
    require(binding["sha256"] == VMAX_SHA, "BIO0R current Vmax binding mismatch")
    relative = str(binding["relative_path"])
    paths["current_vmax"], size = safe_source(repo, relative, VMAX_SHA)
    require(size == binding["bytes"], "current Vmax byte count mismatch")
    rows.append(dict(role="current_vmax", relative_path=relative, sha256=VMAX_SHA, bytes=size, status="PASS"))
    return paths, pd.DataFrame(rows)


def panel_fingerprint(panel: pd.DataFrame) -> str:
    fields = ["algorithm", "tumor", "ensemble_hash", "projection_hash", "rna_context_key", "rna_training_samples", "sample_index"]
    payload = {
        "schema": "bridge3a.candidate_panel.v1", "identity_fields": fields,
        "source_sha256": {"stage12_candidates": SOURCES["candidate_panel"][1], "stage12_provenance": SOURCES["candidate_provenance"][1]},
        "candidates": panel[fields].astype(str).sort_values(fields, kind="stable").to_dict("records"),
    }
    return sha(json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode())


def load_panel(paths: dict[str, Path]) -> tuple[pd.DataFrame, dict[str, object]]:
    panel = pd.read_csv(paths["candidate_panel"], sep="\t", compression="xz")
    require(len(panel) == 640 and not panel.duplicated(KEY).any(), "candidate count or key mismatch")
    require(panel_fingerprint(panel) == PANEL_ID, "candidate panel fingerprint mismatch")
    require(panel.groupby(["algorithm", "tumor"]).size().to_dict() == {(a, t): 80 for a in ("CORDA", "GIMME", "iMAT", "RIPTiDe") for t in ("CT2A", "GL261")}, "candidate method/tumor accounting mismatch")
    prov = json.loads(paths["candidate_provenance"].read_text())
    require(set(panel["ensemble_hash"].astype(str)) == set(map(str, prov["stage11_ensemble_hashes"])), "candidate/provenance ensemble mismatch")
    require(set(panel["projection_hash"].astype(str)) == set(map(str, prov["stage11_projection_hashes"])), "candidate/provenance projection mismatch")
    panel = panel.sort_values(KEY, kind="stable").reset_index(drop=True)

    with np.load(paths["flux_cache"], allow_pickle=False) as cache:
        require(cache["mats"].shape == (32, 20, 4181), "flux-cache shape mismatch")
        rxns = list(map(str, cache["rxns"]))
        require(len(rxns) == len(set(rxns)) == 4181, "flux-cache reaction inventory mismatch")
        needed = ("LDH_L", "HEX1", "AKGDm")
        require(set(needed).issubset(rxns), "required reaction absent from flux cache")
        indices = [rxns.index(r) for r in needed]
        records = []
        mats = cache["mats"]
        for j, (alg, tumor, ensemble) in enumerate(zip(cache["alg"], cache["tumor"], cache["eh"])):
            for sample in range(20):
                records.append((str(alg), str(tumor), str(ensemble), sample, *map(float, mats[j, sample, indices])))
    flux = pd.DataFrame(records, columns=["algorithm", "tumor", "ensemble_hash", "sample_index", *needed])
    short_key = ["algorithm", "tumor", "ensemble_hash", "sample_index"]
    require(len(flux) == 640 and not flux.duplicated(short_key).any(), "flux-cache candidate duplication")
    p = panel.merge(flux, on=short_key, how="left", validate="one_to_one", suffixes=("_panel", "_cache"), indicator=True)
    require(p["_merge"].eq("both").all(), "flux-cache candidate linkage mismatch")
    require(np.isfinite(p[["LDH_L_cache", "HEX1", "AKGDm"]].to_numpy(float)).all(), "nonfinite flux-cache anchor reaction")
    # The flux cache is carbon-normalized for scoring; the frozen quantitative
    # Vlac operator uses raw LDH_L from the candidate panel.
    lac = np.maximum(-p["LDH_L_panel"].to_numpy(float), 0.0)
    require(np.allclose(lac, p["lac_production"].to_numpy(float), atol=ATOL, rtol=0), "candidate lactate semantics mismatch")
    p["candidate_Vmax"] = np.maximum(-p["EX_glc__D_e"].to_numpy(float), 0.0)
    p["candidate_Vlac"] = lac
    for tumor, idx in p.groupby("tumor", sort=True).groups.items():
        require(len(idx) == 320, f"candidate count mismatch: {tumor}")
        for name in ("Vmax", "Vlac"):
            p.loc[idx, "q_" + name] = core.candidate_rank_coordinate(p.loc[idx, "candidate_" + name].to_numpy(float))
    return p.sort_values(KEY, kind="stable").reset_index(drop=True), {"cache_candidate_keys": 640, "cache_reactions": 4181, "panel_id": PANEL_ID}


def mouse_coordinates(dmi: pd.DataFrame) -> pd.DataFrame:
    require(len(dmi) == 10 and not dmi.mouse_id.duplicated().any(), "DMI mouse inventory mismatch")
    require(dmi.groupby("tumor_type").size().to_dict() == {"CT2A": 5, "GL261": 5}, "DMI tumor accounting mismatch")
    rows = []
    for tumor, group in dmi.groupby("tumor_type", sort=True):
        group = group.sort_values("mouse_id", kind="stable")
        for observable in ("Vmax", "Vlac", "Vglx"):
            require(np.isfinite(group[observable].to_numpy(float)).all(), f"nonfinite DMI {observable}")
            q = core.mouse_rank_coordinate(group[observable].to_numpy(float))
            for mouse, coordinate in zip(group.mouse_id, q):
                rows.append(dict(mouse_id=str(mouse), tumor=str(tumor), observable=observable, coordinate=float(coordinate)))
    return pd.DataFrame(rows).sort_values(["tumor", "mouse_id", "observable"], kind="stable").reset_index(drop=True)


def check_historical_lactate(panel: pd.DataFrame, historical: pd.DataFrame) -> dict[str, object]:
    old = historical.loc[historical.observed_targets.eq("Vlac")].copy()
    require(len(old) == 3200, "historical Vlac row count mismatch")
    unique = old.drop_duplicates(KEY)
    require(len(unique) == 640 and old.groupby(KEY).size().eq(5).all(), "historical Vlac candidate accounting mismatch")
    require(old.groupby("tumor").mouse_id.nunique().to_dict() == {"CT2A": 5, "GL261": 5}, "historical Vlac mouse accounting mismatch")
    require(np.isfinite(unique[["Vlac", "q_Vlac"]].to_numpy(float)).all(), "historical Vlac nonfinite")
    rank_diff = 0.0
    for _, idx in unique.groupby("tumor", sort=True).groups.items():
        expected = core.candidate_rank_coordinate(unique.loc[idx, "Vlac"].to_numpy(float))
        rank_diff = max(rank_diff, float(np.max(np.abs(expected - unique.loc[idx, "q_Vlac"].to_numpy(float)))))
    require(rank_diff <= ATOL, "historical Vlac rank semantics mismatch")
    joined = panel[KEY + ["candidate_Vlac"]].merge(unique[KEY + ["Vlac"]], on=KEY, how="inner", validate="one_to_one")
    require(len(joined) == 480 and joined.algorithm.value_counts().to_dict() == {"CORDA": 160, "GIMME": 160, "iMAT": 160}, "historical/current Vlac overlap mismatch")
    raw_diff = float(np.max(np.abs(joined.candidate_Vlac.to_numpy(float) - joined.Vlac.to_numpy(float))))
    require(raw_diff <= ATOL, "historical/current Vlac raw observable mismatch")
    return {"historical_candidates": 640, "current_overlap_candidates": 480, "overlap_by_algorithm": {"CORDA": 160, "GIMME": 160, "iMAT": 160}, "historical_rank_max_abs_diff": rank_diff, "overlap_lactate_max_abs_diff": raw_diff}


def archived_kernel(path: Path):
    # The pinned archived module executes its workflow only under __main__.
    spec = importlib.util.spec_from_file_location("a20_pinned_archived_kernel", path)
    require(spec is not None and spec.loader is not None, "cannot load archived kernel")
    module = importlib.util.module_from_spec(spec)
    prior = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = prior
    return module.kernel_from_distance, module.ess


def qualify(paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    panel, panel_qc = load_panel(paths)
    dmi = pd.read_csv(paths["dmi"])
    mice = mouse_coordinates(dmi)
    historical = pd.read_csv(paths["historical_weights"], sep="\t", compression="xz")
    lactate_qc = check_historical_lactate(panel, historical)

    current = pd.read_csv(paths["current_vmax"], sep="\t", compression="xz")
    require(len(current) == 3200 and not current.duplicated(["mouse_id", *KEY]).any(), "current Vmax row/key mismatch")
    require(current.groupby("mouse_id").size().eq(320).all(), "current Vmax per-mouse count mismatch")
    require(set(map(tuple, current[KEY].drop_duplicates().to_numpy())) == set(map(tuple, panel[KEY].to_numpy())), "current Vmax panel identity mismatch")
    observed = pd.read_csv(paths["current_observables"])
    short_key = ["algorithm", "tumor", "ensemble_hash", "sample_index"]
    check = panel[short_key + ["candidate_Vmax", "candidate_Vlac", "q_Vmax", "q_Vlac"]].merge(observed, on=short_key, validate="one_to_one", suffixes=("_new", "_old"))
    require(len(check) == 640, "current observable row count mismatch")
    obs_diffs = {}
    for new, old in (("candidate_Vmax", "Vmax"), ("candidate_Vlac", "Vlac"), ("q_Vmax_new", "q_Vmax_old"), ("q_Vlac_new", "q_Vlac_old")):
        diff = float(np.max(np.abs(check[new].to_numpy(float) - check[old].to_numpy(float))))
        require(diff <= ATOL, f"current observable mismatch: {new}/{old}")
        obs_diffs[new] = diff

    archived_solve, archived_ess = archived_kernel(paths["archived_kernel"])
    qmouse = mice.pivot(index="mouse_id", columns="observable", values="coordinate")
    base_diff = 0.0
    base_ess_diff = 0.0
    helper_ess_diff = 0.0
    reg_weight_diff = 0.0
    reg_temp_diff = 0.0
    a2_rows = []
    summary_rows = []
    for mouse, existing in current.groupby("mouse_id", sort=True):
        tumor = str(existing.tumor.iloc[0])
        require(existing.tumor.nunique() == 1, "mixed tumor in current Vmax group")
        block = existing.merge(panel[KEY + ["candidate_Vmax", "candidate_Vlac", "q_Vmax", "q_Vlac"]], on=KEY, how="left", validate="one_to_one", suffixes=("_saved", ""), indicator=True)
        require(len(block) == 320 and block["_merge"].eq("both").all(), "current Vmax candidate join mismatch")
        # Retain the pinned Vmax table order through the solver: floating sums
        # change in the last bit if the same 320 candidates are reordered.
        block = block.reset_index(drop=True)
        qv = block.q_Vmax.to_numpy(float)
        ql = block.q_Vlac.to_numpy(float)
        tv = float(qmouse.loc[mouse, "Vmax"])
        tl = float(qmouse.loc[mouse, "Vlac"])
        vdist = np.square(qv - tv)
        archived_vmax, _ = archived_solve(vdist)
        baseline_weights, _ = bio0r.solve_kernel(vdist.copy())
        base = core.solve_temperature_for_ess(vdist)
        base_diff = max(base_diff, float(np.max(np.abs(baseline_weights - block.weight.to_numpy(float)))))
        producing_ess = bio0r.core.effective_sample_size(baseline_weights)
        base_ess_diff = max(base_ess_diff, abs(producing_ess - TARGET_ESS), float(np.max(np.abs(block.achieved_ess.to_numpy(float) - producing_ess))))
        helper_ess_diff = max(helper_ess_diff, abs(base.achieved_ess - producing_ess), float(np.max(np.abs(base.weights - archived_vmax))))
        require(np.allclose(block.Vmax_coordinate.to_numpy(float), qv, atol=ATOL, rtol=0), "saved Vmax coordinate mismatch")
        require(np.allclose(block.q_Vmax_saved.to_numpy(float), block.weight.to_numpy(float), atol=ATOL, rtol=0), "saved q_Vmax/weight semantics mismatch")
        require(np.allclose(block.target_ess.to_numpy(float), TARGET_ESS, atol=ATOL, rtol=0), "saved target ESS mismatch")

        dist = core.dual_squared_distance(qv, ql, tv, tl)
        result = core.solve_temperature_for_ess(dist)
        archived_weights, archived_temperature = archived_solve(dist)
        reg_weight_diff = max(reg_weight_diff, float(np.max(np.abs(result.weights - archived_weights))))
        reg_temp_diff = max(reg_temp_diff, abs(result.temperature - archived_temperature))
        require(np.isfinite(result.weights).all() and (result.weights > 0).all(), "A2-L support or finiteness failure")
        require(abs(float(result.weights.sum()) - 1.0) <= ATOL and abs(result.achieved_ess - TARGET_ESS) <= 1e-10, "A2-L normalization/ESS failure")
        require(abs(result.achieved_ess - archived_ess(archived_weights)) <= 1e-10, "archived ESS regression mismatch")
        for i, row in enumerate(block.itertuples(index=False)):
            a2_rows.append({"arm": "A2-L", "mouse_id": str(mouse), "tumor": tumor, **{k: getattr(row, k) for k in KEY}, "q_Vmax": qv[i], "q_Vlac": ql[i], "mouse_q_Vmax": tv, "mouse_q_Vlac": tl, "distance": float(dist[i]), "weight": float(result.weights[i])})
        summary_rows.append({"arm": "A2-L", "mouse_id": str(mouse), "tumor": tumor, "candidate_rows": len(block), "positive_weight_rows": int(np.count_nonzero(result.weights > 0)), "weight_sum": float(result.weights.sum()), "target_ess": TARGET_ESS, "achieved_ess": result.achieved_ess, "temperature": result.temperature, "max_weight": float(result.weights.max()), "archived_weight_max_abs_diff": float(np.max(np.abs(result.weights - archived_weights))), "archived_temperature_abs_diff": abs(result.temperature - archived_temperature)})

    require(base_diff <= ATOL and base_ess_diff <= ATOL, f"current Vmax baseline reproduction mismatch: weight={base_diff:.17g}, ess={base_ess_diff:.17g}")
    require(reg_weight_diff <= ATOL and reg_temp_diff <= ATOL, "archived temperature regression failed")
    weights = pd.DataFrame(a2_rows).sort_values(["arm", "tumor", "mouse_id", *KEY], kind="stable").reset_index(drop=True)
    summary = pd.DataFrame(summary_rows).sort_values(["arm", "tumor", "mouse_id"], kind="stable").reset_index(drop=True)
    require(len(weights) == 3200 and len(summary) == 10, "A2-L output accounting mismatch")

    audit = json.loads(paths["operator_audit"].read_text())
    g = audit["operators"]["G"]["S"]
    require(audit["candidate_panel_id"] == PANEL_ID and g["status"] == "UNAVAILABLE" and "No independent current-panel Vglx quantitative operator exists" in g["reason"], "A2-G operator audit unexpectedly changed")
    registry = pd.read_csv(paths["bio0_registry"], sep="\t")
    require(((registry.quantity == "historical_vglx_pyr") & (registry.coordinate == "PYRt2m")).any(), "historical Vglx registry mismatch")
    require(((registry.quantity == "current_vglx_akgdm") & (registry.coordinate == "AKGDm") & (registry.bio1_status == "NOT_BIO1")).any(), "AKGDm sensitivity registry mismatch")
    mapping = pd.read_csv(paths["dmi_mapping"], sep="\t")
    require(((mapping.dmi_parameter == "Vglx") & (mapping.candidate_reaction_ids == "AKGDm") & (mapping.active_by_default.astype(str).str.lower() == "false")).any(), "Vglx/AKGDm sensitivity mapping missing")
    reason = "Archived Stage-12 Vglx uses PYRt2m; current AKGDm is a sensitivity association without an independently qualified current-panel quantitative Vglx operator (frozen BRIDGE3A G/S status UNAVAILABLE). Exact archived Vglx semantics cannot be reproduced through AKGDm."
    arms = pd.DataFrame([
        {"arm": "A2-L", "strong_measurements": "Vmax;Vlac", "reaction_associations": "HEX1;LDH_L", "candidate_observables": "max(-EX_glc__D_e,0);max(-LDH_L,0)", "qualification_status": "QUALIFIED", "reason": "Current Vmax and historical Vlac semantics reproduced; archived ESS-20 operator regression passed", "weight_rows": 3200},
        {"arm": "A2-G", "strong_measurements": "Vmax;Vglx", "reaction_associations": "HEX1;AKGDm", "candidate_observables": "max(-EX_glc__D_e,0);AKGDm sensitivity association", "qualification_status": "NOT_QUALIFIED", "reason": reason, "weight_rows": 0},
    ])
    qc = {"panel": panel_qc, "current_observable_max_abs_diffs": obs_diffs, "lactate": lactate_qc, "vmax_weight_max_abs_diff": base_diff, "vmax_ess_max_abs_diff": base_ess_diff, "helper_vs_archived_ess_max_abs_diff": helper_ess_diff, "archived_temperature_regression": {"representative_real_panel_distances": 20, "max_weight_abs_diff": reg_weight_diff, "max_temperature_abs_diff": reg_temp_diff, "archived_source_sha256": SOURCES["archived_kernel"][1]}, "a2_g_reason": reason}
    return arms, mice, weights, summary, qc


def tsv(frame: pd.DataFrame, compressed: bool = False) -> bytes:
    raw = frame.to_csv(sep="\t", index=False, lineterminator="\n", float_format="%.17g").encode("utf-8")
    return lzma.compress(raw, preset=6, format=lzma.FORMAT_XZ) if compressed else raw


def build(repo: Path) -> tuple[dict[str, bytes], dict[str, object]]:
    paths, sources = source_gate(repo)
    arms, mice, weights, summary, qc = qualify(paths)
    status = "BRIDGE_A20_PARTIAL_QUALIFICATION"
    status_doc = {"schema": "bridge.a20.status.v1", "status": status, "arm_status": {"A2-L": "QUALIFIED", "A2-G": "NOT_QUALIFIED"}, "a2_g_reason": qc["a2_g_reason"], "a2_g_outcome_accessed": False, "replacement_anchor_selected": False, "analysis_label": "post-freeze robustness/generalization", "a2_1_a2_2_a2_3_computed": False}
    blobs = {
        "BRIDGEA20_SOURCE_AUDIT.tsv": tsv(sources),
        "BRIDGEA20_ANCHOR_REGISTRY.tsv": tsv(arms),
        "BRIDGEA20_MOUSE_COORDINATES.tsv": tsv(mice),
        "BRIDGEA20_DUAL_ANCHOR_WEIGHTS.tsv.xz": tsv(weights, compressed=True),
        "BRIDGEA20_WEIGHT_SUMMARY.tsv": tsv(summary),
        "BRIDGEA20_STATUS.json": canonical(status_doc),
    }
    manifest = {"schema": "bridge.a20.manifest.v1", "status": status, "analysis_label": status_doc["analysis_label"], "panel_id": PANEL_ID, "source_sha256": {row.role: {"relative_path": row.relative_path, "sha256": row.sha256} for row in sources.itertuples()}, "artifact_sha256": {name: sha(payload) for name, payload in blobs.items()}, "row_counts": {"source_audit": len(sources), "anchor_registry": len(arms), "mouse_coordinates": len(mice), "dual_anchor_weights": len(weights), "weight_summary": len(summary), "A2-L_weights": 3200, "A2-G_weights": 0}, "qualification": {"A2-L": "QUALIFIED", "A2-G": "NOT_QUALIFIED"}, "diagnostics": qc, "scientific_outcomes_computed": False, "a2_g_outcome_accessed": False, "replacement_anchor_selected": False}
    blobs["BRIDGEA20_MANIFEST.json"] = canonical(manifest)
    return blobs, manifest


def publish(repo: Path, blobs: dict[str, bytes]) -> str:
    out = repo / OUT
    if out.exists():
        require(out.is_dir() and not out.is_symlink(), "A2 output path is not a normal directory")
        existing = {p.name for p in out.iterdir()}
        require(existing == set(blobs), "existing A2 output inventory differs")
        for name, expected in blobs.items():
            path = out / name
            require(path.is_file() and not path.is_symlink() and path.read_bytes() == expected, f"existing A2 output differs: {name}")
        return "NO_OP_IDENTICAL_FROZEN_RERUN"
    out.mkdir(parents=False, exist_ok=False)
    for name, payload in blobs.items():
        with (out / name).open("xb") as stream:
            stream.write(payload)
    return "WRITTEN"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    args = parser.parse_args()
    repo = args.repo_root.resolve(strict=True)
    blobs, manifest = build(repo)
    write_status = publish(repo, blobs)
    print(json.dumps({"status": manifest["status"], "write_status": write_status, "qualification": manifest["qualification"], "row_counts": manifest["row_counts"], "max_ess_error": manifest["diagnostics"]["vmax_ess_max_abs_diff"], "max_archived_weight_diff": manifest["diagnostics"]["archived_temperature_regression"]["max_weight_abs_diff"]}, sort_keys=True))


if __name__ == "__main__":
    main()
