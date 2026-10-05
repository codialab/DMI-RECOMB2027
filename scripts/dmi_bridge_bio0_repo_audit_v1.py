#!/usr/bin/env python3
"""Repository-bound, outcome-free DMI-BRIDGE-BIO-0 audit."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from dmi_bridge_bio0_audit_core_v1 import (
    BIO1_FREEZE,
    FOUR_ARMS,
    READY_STATUS,
    WeightRow,
    build_analysis_contract,
    canonical_json_bytes,
    distribution_fingerprint,
    four_arm_decomposition,
)

ROOT = Path(__file__).resolve().parents[1]
PANEL_ID = "9c27c8e75aa128a4a4cb15547f23557845f37d35a7323e875af5a1711a527f42"
FLUX_PANEL_ID = "313b44a46d271dbb471a16ce598ec6bae67074f4e95ec53c3ebca22b2e5cadca"
GENE_ID = "PRIME_SIGNED__COLLECTRI__TF_EXPRESSION_ROBUST_Z__MEDIAN"
FLUX_ID = "FLUX__SUBSYSTEM__ABSOLUTE_STATE_FIRST__R5__NA__N0__C1__GA0__PB0__ANDNA__ORNA__W0__MEDIAN"

SOURCES = {
    "patch": ("dmi_bridge_bio0_audit_v1.patch", "ae7fb3caa8371d16648bb85b8bf1f477830a928ccde8e32c2755299d344e2963"),
    "candidate_panel": ("12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/candidate_fraction_table.tsv.xz", "421c2c1d196fbd088a4e92e2ecc998ad578aab814a025fff8e1c843597fed1da"),
    "candidate_provenance": ("12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/DMI_FRACTION_ENSEMBLE_PROVENANCE.json", "fd05e5febba47a92005619db38c53b3ba4896d81eb498226d87c61f1e9536168"),
    "ensemble_selection": ("12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/ensemble_selection.tsv", "ffaa8bedba3bc1e4a11f2cdb930de7d0b128936afb88f325125380a5586e8c70"),
    "ensemble_weights": ("12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/ensemble_candidate_weights.tsv.xz", "994856d63d1b80898ef549f91952c68c344b28e54bc08197d2291854968cf2ba"),
    "dmi": ("pre-simulation_constraint_analysis/results/run_20260824T164247Z_corrected_analysis_A/mouse_level_dmi_input.csv", "ae3cbad6d422405780d0dc44e0e31f95770ebc204d1f80f2ccdada12bb119082"),
    "historical_vmax": ("12_recomb_method_comparison/outputs_stage12b/dmi_single_anchor_weights.tsv.xz", "b3357587b70329fca7d77b99e4b6c6487feccb54ac7f251d8a0a8437353a8030"),
    "current_vmax": ("13_recomb_identifiability/outputs/sa_u1_u2_followup_20260926/current_vmax_base_weights.tsv.xz", "91121bf1251354505cde5f5ae783dc8d14bc9a22ecdf8daa06be1e5d82e3385d"),
    "current_vmax_manifest": ("13_recomb_identifiability/outputs/sa_u1_u2_followup_20260926/current_vmax_base_weights_manifest.json", "2abf9e7c4d3c922bda8c81426001c1f27c6e43f1761e9c9073e0d8ca078bd5bc"),
    "observables": ("outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/current_candidate_observables.csv", "946cd613d673c9aa9367bce705ad5d05562c6a174acc830fbbd36cb395aa5504"),
    "flux_cache": ("outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstructed/flux_cache.npz", "f9644da4ebd875cf5f0b34cc58b20faf0a0457d7c7a44635e8a010d8896c3dbb"),
    "reproduction": ("outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstructed/REPRODUCTION_COMPARISON.json", "f2525f3ad285c52694f056c6f2a4ab04b44b27cd48456bc29bb7251c28920d79"),
    "sa_contract": ("13_recomb_identifiability/outputs/sa_u1_u2_followup_20260926/ANALYSIS_CONTRACT.json", "ff67c20120eb86e0a1f2730f8af73802ef4c60c4d15a28fbf955020d0be0f8d7"),
    "sa_baseline": ("13_recomb_identifiability/outputs/sa_u1_u2_followup_20260926/BASELINE_REPRODUCTION.tsv", "662a56a67fed85cd8905440070e477ce12fdf95587e3af79d18a27d5b076e773"),
    "sa_random": ("13_recomb_identifiability/outputs/sa_u1_u2_followup_20260926/random_direction_controls.tsv", "610e0c5d99f13c3f570832d82041a4be05a01347666bbbf760e37d5e698f2a1b"),
    "mapping": ("analysis/brain_glioma_gem/config/dmi_to_gem_mapping.tsv", "3002fa055a5efaccc27044c3a158351a3f845e5a31c8e8688e34728e6575b489"),
    "mapping_full": ("analysis/brain_glioma_gem/05_dmi_gem_mapping/dmi_to_full_model_mapping.tsv", "4dd29f9a7b3a0d436f98ecf7e83951b53808138821e15438364f11c81afd1d64"),
    "flux_spec": ("1_pathway_flux_scoring/analysis/pathway_flux_postprocessing/frozen_flux_scoring_spec.json", "a306da97301457b5bd95a761b6a901144828214bea16b6aee985ff5ef15787ab"),
    "score_manifest": ("4_scoring_definition_benchmark/score_manifest.csv.xz", "a802ee300fb228eefa4205e02cdcbe3f0693aa97872519e60a2056f60005e85e"),
    "subsystem_memberships": ("4_scoring_definition_benchmark/results/subsystem_memberships.csv.xz", "54f681f4331da34dd4aaf0df611d8286e7e36e0f92c8e95ead38f023b5fc8b32"),
    "reaction_metadata": ("13_recomb_identifiability/exports/dmi_sa_offline_feasibility/reaction_metadata/reaction_metadata.tsv", "6b7dad70307458fa6634958e4fe817a13c9ceec69e330c9ef5f9627b8902de8f"),
    "stage4_metadata": ("dmi_sa_r1_cup_reuse_20260926/selected_model/ensembles/stage4_projection_4747/metadata.json", "f5542cd608a625f098358c7f152ca6039402c3b473a1eee47f0423aedc9f9f53"),
    "gene_score": ("9_recomb_gene_score_validation_combined/02_khalsa_to_mikolajewicz_external/results/mikolajewicz_gene_effects_canonical.tsv.xz", "c9dabb2b2c3fa174b7bc2d797e7ff1ad6debda10d644e707c7e6be8764cfa9f1"),
    "gene_readme": ("9_recomb_gene_score_validation_combined/02_khalsa_to_mikolajewicz_external/README.md", "88b0e2079cacd7b3e64ba22e979a4514137b9610765fdb1c655a8cd25bb15cb0"),
    "gene_sample_qc": ("9_recomb_gene_score_validation_combined/02_khalsa_to_mikolajewicz_external/inputs_frozen/mikolajewicz_sample_qc.tsv", "b8c96d91c1cb3f1139e983b1e550a9ea440aacbb0c33530eeadb620f6cd1f1b5"),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def source_audit(repo: Path) -> dict:
    records = {}
    for role, (rel, expected) in SOURCES.items():
        p = (repo / rel).resolve(strict=True)
        p.relative_to(repo.resolve())
        got = sha256(p)
        if got != expected:
            raise RuntimeError(f"source hash mismatch for {role}: {got} != {expected}")
        records[role] = {"relative_path": rel, "sha256": got}
    return records


def panel_payload(c: pd.DataFrame) -> dict:
    fields = ["algorithm", "tumor", "ensemble_hash", "projection_hash", "rna_context_key", "rna_training_samples", "sample_index"]
    return {"schema": "bridge3a.candidate_panel.v1", "identity_fields": fields,
            "source_sha256": {"stage12_candidates": SOURCES["candidate_panel"][1], "stage12_provenance": SOURCES["candidate_provenance"][1]},
            "candidates": c[fields].astype(str).sort_values(fields, kind="stable").to_dict("records")}


def make_weights(c: pd.DataFrame, obs: pd.DataFrame, qtab: pd.DataFrame, dmi: pd.DataFrame):
    key = ["mouse_id", "algorithm", "tumor", "ensemble_hash", "sample_index"]
    candidate_key = ["algorithm", "tumor", "ensemble_hash", "sample_index"]
    joined = qtab.merge(obs[candidate_key + ["q_Vlac"]], on=candidate_key, how="inner", validate="many_to_one")
    if len(joined) != 3200:
        raise RuntimeError("current Vmax table does not join exactly to 3200 current mouse-candidate rows")
    records = []
    support = []
    for mouse, g in joined.groupby("mouse_id", sort=True):
        tumor = str(g.tumor.iloc[0])
        q = g.weight.to_numpy(float)
        if not np.isfinite(q).all() or (q < 0).any() or q.sum() <= 0:
            raise RuntimeError(f"invalid baseline weights for {mouse}")
        direction = 1.0 if tumor == "CT2A" else -1.0
        score = direction * (2.0 * g.q_Vlac.to_numpy(float) - 1.0)
        z = 0.5 * score
        z -= z.max()
        p = q * np.exp(z)
        p /= p.sum()
        strata = (g.algorithm.astype(str) + "::" + g.rna_context_key.astype(str)).to_numpy()
        rows0 = [WeightRow("|".join(map(str, row)), str(st), float(w)) for row, st, w in zip(g[candidate_key].itertuples(index=False, name=None), strata, q)]
        rows1 = [WeightRow(r.candidate_id, r.stratum_id, float(w)) for r, w in zip(rows0, p)]
        arms, diag = four_arm_decomposition(rows0, rows1)
        if diag.baseline_zero_mass_strata or diag.updated_zero_mass_strata or diag.nonderivable_within_context_strata:
            raise RuntimeError(f"unexpected literal-zero support in {mouse}: {diag}")
        amap = {name: {r.candidate_id: r.weight for r in rows} for name, rows in arms.items()}
        for row, st in zip(g[candidate_key].itertuples(index=False, name=None), strata):
            cid = "|".join(map(str, row))
            for arm in FOUR_ARMS:
                records.append({"mouse_id": str(mouse), "tumor": tumor, "algorithm": row[0], "ensemble_hash": row[2], "sample_index": int(row[3]), "rna_context_key": str(st.split("::", 1)[1]), "stratum_id": st, "arm": arm, "weight": float(amap[arm][cid])})
        for st in sorted(set(strata)):
            mask = strata == st
            support.append({"mouse_id": str(mouse), "tumor": tumor, "stratum_id": st, "candidate_rows": int(mask.sum()), "baseline_mass": float(q[mask].sum()), "updated_mass": float(p[mask].sum()), "baseline_positive": bool(q[mask].sum() > 0.0), "updated_positive": bool(p[mask].sum() > 0.0), "nonderivable": False, "support_invented": False, **{f"{arm}_mass": float(sum(amap[arm][r] for r in ["|".join(map(str, x)) for x in g.loc[mask, candidate_key].itertuples(index=False, name=None)])) for arm in FOUR_ARMS}})
    return pd.DataFrame(records), pd.DataFrame(support)


def build(repo: Path, out: Path) -> None:
    source_records = source_audit(repo)
    c = pd.read_csv(repo / SOURCES["candidate_panel"][0], sep="\t", compression="xz")
    prov = json.loads((repo / SOURCES["candidate_provenance"][0]).read_text())
    payload = panel_payload(c)
    panel_id = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()).hexdigest()
    if panel_id != PANEL_ID or len(c) != 640 or c.duplicated(["algorithm", "tumor", "ensemble_hash", "sample_index"]).any():
        raise RuntimeError("candidate panel identity/count gate failed")
    if set(c.algorithm) != {"CORDA", "GIMME", "iMAT", "RIPTiDe"} or set(c.tumor) != {"CT2A", "GL261"}:
        raise RuntimeError("candidate method/tumor gate failed")
    if not all(v == 4 for v in c.groupby(["algorithm", "tumor"]).rna_context_key.nunique()) or not all(v == 20 for v in c.groupby(["algorithm", "tumor", "ensemble_hash"]).size()):
        raise RuntimeError("candidate context/vector count gate failed")
    dmi = pd.read_csv(repo / SOURCES["dmi"][0])
    if dmi.groupby("tumor_type").mouse_id.nunique().to_dict() != {"CT2A": 5, "GL261": 5}:
        raise RuntimeError("DMI mouse gate failed")
    obs = pd.read_csv(repo / SOURCES["observables"][0], dtype={"sample_index": str})
    qtab = pd.read_csv(repo / SOURCES["current_vmax"][0], sep="\t", compression="xz", dtype={"sample_index": str})
    identity_cols = ["algorithm", "tumor", "ensemble_hash", "projection_hash", "rna_context_key", "sample_index"]
    panel_identity = {tuple(row) for row in c[identity_cols].astype(str).drop_duplicates().itertuples(index=False, name=None)}
    vmax_identity = {tuple(row) for row in qtab[identity_cols].astype(str).drop_duplicates().itertuples(index=False, name=None)}
    if panel_identity != vmax_identity or len(vmax_identity) != 640:
        raise RuntimeError("current Vmax cache is not identity-compatible with the 640-candidate panel")
    weights, support = make_weights(c, obs, qtab, dmi)
    prime = pd.read_csv(repo / SOURCES["gene_score"][0], sep="\t")
    prime = prime[prime.gene_score_id.eq(GENE_ID)].dropna(subset=["subsystem", "rna_effect"])
    if len(prime) != 49 or prime.subsystem.duplicated().any():
        raise RuntimeError("PRIME-49 gate failed")
    cache = np.load(repo / SOURCES["flux_cache"][0])
    stage4 = set(map(str, json.loads((repo / SOURCES["stage4_metadata"][0]).read_text())["ordered_projection_reactions"]))
    meta = pd.read_csv(repo / SOURCES["reaction_metadata"][0], sep="\t")["reaction_id subsystem".split()].dropna().drop_duplicates()
    panel = meta[meta.reaction_id.astype(str).isin(set(map(str, cache["rxns"]))) & meta.reaction_id.astype(str).isin(stage4) & meta.subsystem.astype(str).isin(set(prime.subsystem))].drop_duplicates().sort_values(["subsystem", "reaction_id"])
    panel_bytes = ("subsystem\treaction_id\n" + "".join(f"{r.subsystem}\t{r.reaction_id}\n" for r in panel.itertuples())).encode()
    if len(panel) != 3893 or panel.subsystem.nunique() != 49 or hashlib.sha256(panel_bytes).hexdigest() != FLUX_PANEL_ID:
        raise RuntimeError("reaction panel gate failed")
    registry = pd.DataFrame([
        ["Vmax_strong", "EX_glc__D_e", "within_tumor_midrank", "ESS20 exponential squared-rank kernel", "BIO1_ALLOWED"],
        ["LDH_net", "LDH_L", "raw_signed", "negative is pyruvate-to-lactate", "PROVENANCE_ONLY"],
        ["Vlac_lactate", "max(-LDH_L,0)", "q_Vlac=(rank-.5)/320", "primary historical weak coordinate", "BIO1_ALLOWED"],
        ["weak_reversed", "max(-LDH_L,0)", "same as primary", "score sign reversed", "CONTROL_ONLY"],
        ["random_direction_family", "candidate identity", "16 seeded sign vectors", "seeds 20260926..20260941", "CONTROL_ONLY"],
        ["lactate_transport", "L_LACt2r", "reversible", "negative exports", "NOT_BIO1"],
        ["lactate_exchange", "EX_lac__L_e", "positive secretion", "export interface", "NOT_BIO1"],
        ["historical_f_lac", "max(-LDH_L,0)/(max(-LDH_L,0)+max(PDHm,0))", "allocation fraction", "historical proxy only", "NOT_BIO1"],
        ["historical_f_lac_plus_pcm", "max(-LDH_L,0)/(max(-LDH_L,0)+max(PDHm,0)+max(PCm,0))", "allocation fraction", "sensitivity only", "NOT_BIO1"],
        ["historical_vglx_pyr", "PYRt2m", "forward", "historical proxy", "NOT_BIO1"],
        ["historical_vglx_pdh_pc", "PDHm+PCm", "forward sum", "historical proxy", "NOT_BIO1"],
        ["current_vglx_akgdm", "AKGDm", "forward association", "V_glx->AKGDm; no numerical conversion", "NOT_BIO1"],
        ["pdhm_reference", "PDHm", "forward", "unanchored/held-out reference", "NOT_BIO1"],
    ], columns=["quantity", "coordinate", "sign_or_transform", "role", "bio1_status"])
    gene_prov = pd.DataFrame([
        ["Khalsa_RNA", "seta1,seta2,seta3,setb1,setb2,setb3", "Stage-11 contextualization", "overlap with flux-side reconstruction", "not independent of flux generation"],
        ["DMI_mice", "5 CT2A + 5 GL261", "Vmax baseline and weak direction", "same biological DMI support", "operator provenance, not gene validation"],
        ["Mikolajewicz_pseudobulks", "4 CT2A + 3 GL261", "PRIME gene score", "previously used in Stage-7/SA concordance", "OVERLAPPING_NOT_INDEPENDENT"],
        ["CollecTRI_PRIME_definition", "regulatory network and robust-z expression", "gene-score construction", "defines the selected endpoint", "exploratory endpoint provenance"],
        ["Historical_SA_outcome_use", "frozen PRIME-49 analysis", "operator/endpoint motivation", "prior inspection before BIO", "not untouched confirmation"],
    ], columns=["dataset", "samples_or_cohort", "role", "overlap", "interpretation"])
    contract = build_analysis_contract(candidate_panel_id=PANEL_ID, baseline_operator_id="vmax_ess20_current_panel", weak_operator_id="historical_lactate_exp_tilt_lambda_0.5", reaction_panel_id=FLUX_PANEL_ID, flux_score_id=FLUX_ID, gene_score_id=GENE_ID, gene_score_provenance="OVERLAPPING_NOT_INDEPENDENT", controls=["reversed_direction", "historical_deterministic_randomized_direction"])
    contract["support_audit"] = {"literal_zero_strata": True, "baseline_zero": 0, "updated_zero": 0, "nonderivable": 0, "invented": 0, "rows": 160}
    contract["candidate_panel_fingerprint"] = PANEL_ID
    contract["reaction_panel_fingerprint"] = FLUX_PANEL_ID
    contract["outcome_accessed_during_bio0"] = False
    contract["new_solver_or_sampling_invoked"] = False
    source_audit_doc = {"schema": "bridge.bio0.source_audit.v1", "status": READY_STATUS, "sources": source_records, "candidate_panel": {"id": PANEL_ID, "rows": 640, "methods": sorted(c.algorithm.unique()), "tumors": sorted(c.tumor.unique()), "contexts_per_method_tumor": 4, "vectors_per_ensemble": 20}, "historical_vmax": {"rows_vmax": 3200, "identity_overlap_current": 480, "current_only": 160, "historical_only": 160, "status": "PRESENT_DO_NOT_USE_FOR_CURRENT_PANEL"}, "flux_score_id": FLUX_ID, "gene_score_id": GENE_ID, "gene_score_provenance": "OVERLAPPING_NOT_INDEPENDENT", "new_outcome_generated": False}
    manifest = {"schema": "bridge.bio0.manifest.v1", "status": READY_STATUS, "outputs": {}, "row_counts": {"operator_registry": len(registry), "gene_score_provenance": len(gene_prov), "reaction_panel": len(panel), "support_audit": len(support), "four_arm_weights": len(weights)}}
    with tempfile.TemporaryDirectory(prefix="dmi_bridge_bio0_", dir=str(out.parent)) as td:
        stage = Path(td)
        (stage / "BRIDGEBIO0_SOURCE_AUDIT.json").write_bytes(canonical_json_bytes(source_audit_doc))
        registry.to_csv(stage / "BRIDGEBIO0_OPERATOR_REGISTRY.tsv", sep="\t", index=False)
        gene_prov.to_csv(stage / "BRIDGEBIO0_GENE_SCORE_PROVENANCE.tsv", sep="\t", index=False)
        panel.to_csv(stage / "BRIDGEBIO0_REACTION_PANEL.tsv", sep="\t", index=False)
        support.sort_values(["mouse_id", "stratum_id"]).to_csv(stage / "BRIDGEBIO0_SUPPORT_AUDIT.tsv", sep="\t", index=False)
        weights.sort_values(["mouse_id", "stratum_id", "algorithm", "ensemble_hash", "sample_index", "arm"]).to_csv(stage / "BRIDGEBIO0_FOUR_ARM_WEIGHTS.tsv.xz", sep="\t", index=False, compression={"method": "xz", "preset": 9})
        (stage / "BRIDGEBIO0_ANALYSIS_CONTRACT.json").write_bytes(canonical_json_bytes(contract))
        for p in sorted(stage.iterdir()):
            manifest["outputs"][p.name] = {"sha256": sha256(p), "bytes": p.stat().st_size}
        (stage / "BRIDGEBIO0_MANIFEST.json").write_bytes(canonical_json_bytes(manifest))
        if out.exists():
            existing = {p.name: p for p in out.iterdir() if p.is_file()}
            staged = {p.name: p for p in stage.iterdir()}
            if set(existing) != set(staged) or any(existing[n].read_bytes() != staged[n].read_bytes() for n in staged):
                raise RuntimeError(f"refusing to overwrite differing existing output: {out}")
        else:
            out.mkdir(parents=True)
            for p in stage.iterdir():
                shutil.copy2(p, out / p.name)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=ROOT)
    ap.add_argument("--output-dir", type=Path, default=ROOT / "outputs/dmi_bridge_bio0_audit_v1")
    args = ap.parse_args()
    repo = args.repo_root.resolve(strict=True)
    if repo != ROOT.resolve():
        raise RuntimeError("--repo-root must be the native repository root")
    build(repo, args.output_dir.resolve())
    print(json.dumps({"status": READY_STATUS, "output_dir": str(args.output_dir.resolve())}, sort_keys=True))


if __name__ == "__main__":
    main()
