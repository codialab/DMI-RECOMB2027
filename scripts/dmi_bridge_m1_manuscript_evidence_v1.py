#!/usr/bin/env python3
"""Publish the frozen DMI-BRIDGE-M1 manuscript evidence without scientific analysis."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import tempfile
from pathlib import Path

try:
    from scripts import dmi_bridge_m1_manuscript_core_v1 as core
except ModuleNotFoundError:
    import dmi_bridge_m1_manuscript_core_v1 as core

ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"
DEFAULT_OUTPUT = OUTPUTS / "dmi_bridge_m1_manuscript_evidence_v1"
PATCH = ROOT / "dmi_bridge_m1_manuscript_evidence_foundation_v1.patch"
CONTRACT = ROOT / "docs/DMI_BRIDGE_M1_MANUSCRIPT_EVIDENCE_CONSOLIDATION.md"
CORE = ROOT / "scripts/dmi_bridge_m1_manuscript_core_v1.py"
PATCH_SHA256 = "20dc51e5df3a0d0fd419431faab7356557f082c0a1f8733d35e39fc6ab1d1ac1"
SOURCE_BUNDLES = (
    ("PL2B", "dmi_bridge_pl2b_sign_only_utility_v1", "BRIDGEPL2B_MANIFEST.json", core.PL2B_MANIFEST_SHA256),
    ("PL2C", "dmi_bridge_pl2c_geometry_utility_development_v1", "BRIDGEPL2C_MANIFEST.json", core.PL2C_MANIFEST_SHA256),
    ("PL2D_H", "dmi_bridge_pl2d_hypothesis_freeze_v1", "BRIDGEPL2DH_MANIFEST.json", core.PL2DH_MANIFEST_SHA256),
    ("PL2D_H_AMENDMENT", "dmi_bridge_pl2d_h_support_gate_amendment_v1", "BRIDGEPL2DHA_MANIFEST.json", core.PL2DHA_MANIFEST_SHA256),
    ("PL2D", "dmi_bridge_pl2d_confirmation_v1", "BRIDGEPL2D_MANIFEST.json", core.PL2D_MANIFEST_SHA256),
)
PINNED = {
    "PL2B": {"BRIDGEPL2B_CONTEXT_SUMMARY.tsv": core.PL2B_CONTEXT_SUMMARY_SHA256},
    "PL2C": {},
    "PL2D_H": {"BRIDGEPL2DH_DEVELOPMENT_SNAPSHOT.json": core.PL2DH_DEVELOPMENT_SNAPSHOT_SHA256},
    "PL2D_H_AMENDMENT": {},
    "PL2D": {
        "BRIDGEPL2D_PRIMARY_RESULT.json": core.PL2D_PRIMARY_SHA256,
        "BRIDGEPL2D_SUPPORTIVE_RESULTS.json": core.PL2D_SUPPORTIVE_SHA256,
        "BRIDGEPL2D_CONTEXT_RESULTS.tsv": core.PL2D_CONTEXT_SHA256,
        "BRIDGEPL2D_QC.json": core.PL2D_QC_SHA256,
    },
}
UTILITY_FIELDS = (
    "reliability_q", "total_distinct_reaction_truth_pair_cases", "non_tie_count",
    "truth_tie_count", "improved_count", "improved_fraction_of_non_ties",
    "tied_count", "tied_fraction_of_non_ties", "harmed_count",
    "harmed_fraction_of_non_ties", "gain_median", "gain_q10", "gain_q90",
)
CONTEXT_FIELDS = (
    "algorithm", "rna_context_key", "finite_primary_reaction_pairs", "rho",
    "direction", "evaluable", "support_gate_ge_30",
)
ARTIFACTS = (
    "BRIDGEM1_SOURCE_AUDIT.json", "BRIDGEM1_CLAIM_LEDGER.tsv",
    "BRIDGEM1_RESULTS_FACTS.tsv", "BRIDGEM1_FIGURE_PANEL_PLAN.tsv",
    "BRIDGEM1_PL2B_UTILITY.tsv", "BRIDGEM1_DEVELOPMENT_CONFIRMATION.tsv",
    "BRIDGEM1_CONFIRMATION_CONTEXTS.tsv", "BRIDGEM1_LIMITATIONS.tsv",
)
FALSE_FLAGS = (
    "new_scientific_analysis_performed", "new_inferential_statistics_computed",
    "upstream_outputs_modified", "development_confirmation_pooled",
    "reaction_ranking_performed", "model_fitting_invoked", "solver_invoked",
    "fva_invoked", "sampling_invoked", "reconstruction_invoked",
)


def safe_file(path: Path) -> Path:
    resolved = path.resolve(strict=True)
    if path.is_symlink() or not resolved.is_file() or not resolved.is_relative_to(ROOT):
        raise RuntimeError(f"unsafe or missing repository file: {path}")
    return resolved


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with safe_file(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def json_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def tsv_bytes(fields: tuple[str, ...], rows: list[dict[str, object]]) -> bytes:
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "") for field in fields})
    return out.getvalue().encode("utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with safe_file(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def source_path(stage: str, name: str) -> str:
    return f"outputs/{next(sub for key, sub, _, _ in SOURCE_BUNDLES if key == stage)}/{name}"


def verify_sources() -> dict[str, object]:
    if sha(PATCH) != PATCH_SHA256:
        raise RuntimeError("M1 patch SHA256 mismatch")
    source_hashes: dict[str, str] = {}
    manifests: dict[str, dict] = {}
    for stage, sub, manifest_name, expected in SOURCE_BUNDLES:
        base = OUTPUTS / sub
        manifest_path = base / manifest_name
        if sha(manifest_path) != expected:
            raise RuntimeError(f"{stage} manifest SHA256 mismatch")
        source_hashes[source_path(stage, manifest_name)] = expected
        manifest = json.loads(safe_file(manifest_path).read_text(encoding="utf-8"))
        declared = manifest.get("artifact_sha256")
        if not isinstance(declared, dict) or not declared:
            raise RuntimeError(f"{stage} artifact inventory missing")
        for name, expected_artifact in sorted(declared.items()):
            if Path(name).name != name or sha(base / name) != expected_artifact:
                raise RuntimeError(f"{stage} artifact SHA256 mismatch: {name}")
            source_hashes[source_path(stage, name)] = expected_artifact
        for name, expected_artifact in PINNED[stage].items():
            if declared.get(name) != expected_artifact:
                raise RuntimeError(f"{stage} pinned artifact mismatch: {name}")
        manifests[stage] = manifest
    if manifests["PL2D"].get("status") != "PL2D_CONFIRMED":
        raise RuntimeError("PL2D is not confirmed")
    return {
        "schema": "bridge.m1.source_audit.v1", "status": "PASS",
        "source_sha256": dict(sorted(source_hashes.items())),
        "manifest_artifact_counts": {
            stage: len(manifests[stage]["artifact_sha256"])
            for stage, _, _, _ in SOURCE_BUNDLES
        },
        "pl2d_status": "PL2D_CONFIRMED",
        "patch_sha256": PATCH_SHA256,
        "m1_core_sha256": sha(CORE),
        "adapter_sha256": sha(Path(__file__)),
        "contract_document_sha256": sha(CONTRACT),
    }


def frozen_rows() -> tuple[list[dict], list[dict], list[dict], dict, dict, dict]:
    b = read_tsv(ROOT / source_path("PL2B", "BRIDGEPL2B_CONTEXT_SUMMARY.tsv"))
    utility = [r for r in b if r["summary_scope"] == "OVERALL" and r["lambda_total"] == "0.25" and r["lambda_role"] == "PRIMARY" and r["reliability_q"] in ("0", "0.5", "1")]
    if len(utility) != 3 or {r["reliability_q"] for r in utility} != {"0", "0.5", "1"}:
        raise RuntimeError("PL2B primary utility row inventory mismatch")
    utility.sort(key=lambda r: float(r["reliability_q"]))
    core.validate_pl2b_display_rows(utility)
    for row in utility:
        if row["fraction_denominator"] != "NON_TIE_CASES" or row["quantile_alias_weighting"] != "DISTINCT_TRUTH_PAIR_ONCE":
            raise RuntimeError("PL2B denominator or alias weighting mismatch")
        if row["algorithm"] != "ALL" or row["rna_context_key"] != "ALL":
            raise RuntimeError("PL2B overall population mismatch")
        n = int(row["non_tie_count"])
        if int(row["total_distinct_reaction_truth_pair_cases"]) != n + int(row["truth_tie_count"]):
            raise RuntimeError("PL2B case accounting mismatch")
        if n != sum(int(row[f"{kind}_count"]) for kind in ("improved", "tied", "harmed")):
            raise RuntimeError("PL2B outcome accounting mismatch")
        for kind in ("improved", "tied", "harmed"):
            if not math.isclose(float(row[f"{kind}_fraction_of_non_ties"]), int(row[f"{kind}_count"]) / n, abs_tol=5e-16):
                raise RuntimeError("PL2B fraction/count mismatch")
    h = json.loads(safe_file(ROOT / source_path("PL2D_H", "BRIDGEPL2DH_DEVELOPMENT_SNAPSHOT.json")).read_text())
    dev = {"finite_primary_reactions": h["finite_primary_reactions"], "eta2_usefulness_rho": h["primary_reaction_level_rho"], "positive_contexts": h["positive_algorithm_rna_contexts"]}
    core.validate_development_snapshot(dev)
    if len(h["algorithm_rna_contexts"]) != 16:
        raise RuntimeError("development context count mismatch")
    p = json.loads(safe_file(ROOT / source_path("PL2D", "BRIDGEPL2D_PRIMARY_RESULT.json")).read_text())
    core.validate_confirmation_primary(p)
    s = json.loads(safe_file(ROOT / source_path("PL2D", "BRIDGEPL2D_SUPPORTIVE_RESULTS.json")).read_text())
    core.validate_supportive(s)
    qc = json.loads(safe_file(ROOT / source_path("PL2D", "BRIDGEPL2D_QC.json")).read_text())
    core.validate_confirmation_qc(qc)
    contexts = read_tsv(ROOT / source_path("PL2D", "BRIDGEPL2D_CONTEXT_RESULTS.tsv"))
    expected_keys = {(r["algorithm"], r["rna_context_key"]) for r in h["algorithm_rna_contexts"]}
    if len(contexts) != 16 or {(r["algorithm"], r["rna_context_key"]) for r in contexts} != expected_keys:
        raise RuntimeError("confirmation context inventory mismatch")
    if len({(r["algorithm"], r["rna_context_key"]) for r in contexts}) != 16 or {r["algorithm"] for r in contexts} != {"GIMME", "iMAT", "CORDA", "RIPTiDe"}:
        raise RuntimeError("confirmation algorithm inventory mismatch")
    if sum(r["evaluable"] == "true" for r in contexts) != 16 or sum(r["direction"] == "POSITIVE" and float(r["rho"]) > 0 for r in contexts) != 16 or min(int(r["finite_primary_reaction_pairs"]) for r in contexts) != 34:
        raise RuntimeError("confirmation context gates mismatch")
    if any(r["support_gate_ge_30"] != "true" for r in contexts):
        raise RuntimeError("confirmation support gate mismatch")
    weak = [r for r in contexts if r["algorithm"] == "RIPTiDe" and r["rna_context_key"] == "training_samples=setx2,setx3"]
    if len(weak) != 1 or weak[0]["finite_primary_reaction_pairs"] != "39" or not math.isclose(float(weak[0]["rho"]), 0.27561031478816606, abs_tol=1e-15):
        raise RuntimeError("weak RIPTiDe context missing or changed")
    if len(p["context_results"]) != 16 or any(
        (a["algorithm"], a["rna_context_key"], a["finite_primary_reaction_pairs"]) !=
        (b["algorithm"], b["rna_context_key"], int(b["finite_primary_reaction_pairs"]))
        or not math.isclose(a["rho"], float(b["rho"]), abs_tol=1e-15)
        for a, b in zip(p["context_results"], contexts)
    ):
        raise RuntimeError("PL2D primary/context disagreement")
    return utility, contexts, h, p, s, qc


def build_artifacts() -> dict[str, bytes]:
    audit = verify_sources()
    utility, contexts, development, primary, supportive, qc = frozen_rows()
    blobs: dict[str, bytes] = {"BRIDGEM1_SOURCE_AUDIT.json": json_bytes(audit)}
    utility_rows = []
    for row in utility:
        q = float(row["reliability_q"])
        utility_rows.append({
            **{field: row[field] for field in UTILITY_FIELDS},
            "arm_label": "RANDOM_SIGN_CONTROL" if q == 0.5 else f"RELIABILITY_Q_{q:.2f}",
            "fraction_denominator": row["fraction_denominator"],
            "quantile_alias_weighting": row["quantile_alias_weighting"],
            **{f"{kind}_manuscript_percent": f"{100 * round(float(row[f'{kind}_fraction_of_non_ties']), 4):.2f}%" for kind in ("improved", "tied", "harmed")},
        })
    blobs["BRIDGEM1_PL2B_UTILITY.tsv"] = tsv_bytes((*UTILITY_FIELDS, "arm_label", "fraction_denominator", "quantile_alias_weighting", "improved_manuscript_percent", "tied_manuscript_percent", "harmed_manuscript_percent"), utility_rows)
    dc_rows = [
        {"cohort": "DEVELOPMENT", "status": development["status"], "finite_reactions": development["finite_primary_reactions"], "eta2_usefulness_rho": development["primary_reaction_level_rho"], "positive_contexts": development["positive_algorithm_rna_contexts"], "context_count": len(development["algorithm_rna_contexts"]), "minimum_finite_context_count": "", "prespecified_reaction_holdout": "false", "bootstrap_median": "", "bootstrap_interval_lower": "", "bootstrap_interval_upper": "", "bootstrap_interpretation": "", "source_artifact": source_path("PL2D_H", "BRIDGEPL2DH_DEVELOPMENT_SNAPSHOT.json")},
        {"cohort": "CONFIRMATION", "status": primary["status"], "finite_reactions": primary["finite_reactions"], "eta2_usefulness_rho": primary["overall_rho"], "positive_contexts": primary["positive_contexts"], "context_count": primary["context_count"], "minimum_finite_context_count": primary["context_min_finite_reactions"], "prespecified_reaction_holdout": "true", "bootstrap_median": supportive["bootstrap"]["rho_median"], "bootstrap_interval_lower": supportive["bootstrap"]["rho_q025"], "bootstrap_interval_upper": supportive["bootstrap"]["rho_q975"], "bootstrap_interpretation": "computational reaction-level stability, not biological replication", "source_artifact": source_path("PL2D", "BRIDGEPL2D_PRIMARY_RESULT.json")},
    ]
    dc_fields = tuple(dc_rows[0])
    blobs["BRIDGEM1_DEVELOPMENT_CONFIRMATION.tsv"] = tsv_bytes(dc_fields, dc_rows)
    blobs["BRIDGEM1_CONFIRMATION_CONTEXTS.tsv"] = tsv_bytes(CONTEXT_FIELDS, [{f: r[f] for f in CONTEXT_FIELDS} for r in contexts])
    limitation_values = (
        ("confirmation_cases", qc["counts"]["confirmation_cases"], "cases", "counts.confirmation_cases"),
        ("non_tie_cases", qc["counts"]["non_tie_cases"], "cases", "counts.non_tie_cases"),
        ("truth_tie_cases", qc["counts"]["truth_tie_cases"], "cases", "counts.truth_tie_cases"),
        ("confirmation_evaluation_reaction_rows", qc["counts"]["confirmation_evaluation_reaction_rows"], "rows", "counts.confirmation_evaluation_reaction_rows"),
        ("defined_directional_response_rows", core.EXPECTED_DEFINED_ROWS, "rows", "aggregate_status=DEFINED"),
        ("no_directional_truth_rows", core.EXPECTED_NO_DIRECTIONAL_ROWS, "rows", "aggregate_status=NO_DIRECTIONAL_TRUTH"),
        ("mean_non_tie_pair_fraction", "0.42593732488469338", "evaluation_reaction_rows", "non_tie_pair_fraction"),
    )
    if core.EXPECTED_DEFINED_ROWS + core.EXPECTED_NO_DIRECTIONAL_ROWS != qc["counts"]["confirmation_evaluation_reaction_rows"]:
        raise RuntimeError("coverage accounting mismatch")
    limitations = [{"metric": k, "exact_value": v, "manuscript_display": "0.425937" if k == "mean_non_tie_pair_fraction" else v, "unit_or_denominator": unit, "source_artifact": source_path("PL2D", "BRIDGEPL2D_QC.json" if field.startswith("counts.") else "BRIDGEPL2D_CONFIRMATION_EVALUATION_REACTION.tsv.xz"), "source_field": field, "interpretation": "Directional usefulness is conditional on a defined binary directional truth; a substantial fraction of benchmark cases are truth ties. Computational reactions/evaluations are not independent biological replicates."} for k, v, unit, field in limitation_values]
    blobs["BRIDGEM1_LIMITATIONS.tsv"] = tsv_bytes(tuple(limitations[0]), limitations)
    claim_sources = {
        "C1_SIGN_ONLY_UTILITY_NOT_UNIVERSAL": ("PL2B", source_path("PL2B", "BRIDGEPL2B_CONTEXT_SUMMARY.tsv"), "q=0/0.5/1 improved/tied/harmed percentages; lambda=0.25"),
        "C2_GEOMETRY_PREDICTS_DIRECTIONAL_USEFULNESS": ("PL2D_H;PL2D", source_path("PL2D_H", "BRIDGEPL2DH_DEVELOPMENT_SNAPSHOT.json") + ";" + source_path("PL2D", "BRIDGEPL2D_PRIMARY_RESULT.json"), "development rho=0.8570601823118907;n=2455;confirmation rho=0.8351435823183292;n=838;16/16 positive"),
        "C3_OCCURRENCE_STRONGER_THAN_GAIN_MAGNITUDE": ("PL2D", source_path("PL2D", "BRIDGEPL2D_PRIMARY_RESULT.json") + ";" + source_path("PL2D", "BRIDGEPL2D_SUPPORTIVE_RESULTS.json"), "usefulness rho=0.8351435823183292;information-advantage rho=0.34542078458252307"),
        "C4_CONTEXT_AND_PATHWAY_STABILITY": ("PL2D", source_path("PL2D", "BRIDGEPL2D_CONTEXT_RESULTS.tsv") + ";" + source_path("PL2D", "BRIDGEPL2D_SUPPORTIVE_RESULTS.json"), "16/16 positive;minimum n=34;pathway rho=0.8386984413243371/0.8378450478802439/0.8386984413243371"),
        "L1_DIRECTIONAL_TRUTH_COVERAGE": ("PL2D", source_path("PL2D", "BRIDGEPL2D_QC.json"), "873620 cases;359099 non-ties;514521 truth ties;164690 defined rows;221960 undefined rows"),
    }
    claims = []
    for spec in core.claim_specs():
        stage, artifact, exact = claim_sources[spec["claim_id"]]
        claims.append({"claim_id": spec["claim_id"], "claim_level": spec["claim_level"], "manuscript_safe_statement": spec["statement"], "evidence_stage": stage, "source_artifact": artifact, "exact_supporting_values": exact, "manuscript_display_values": exact, "interpretation_boundary": spec["boundary"], "prohibited_overclaim": spec["boundary"]})
    blobs["BRIDGEM1_CLAIM_LEDGER.tsv"] = tsv_bytes(tuple(claims[0]), claims)
    panel_details = {
        "A": ("conceptual schematic", "HEX1-conditioned strong anchor; conditional weak direction utility", "No quantitative recovery or causal interpretation"),
        "B": ("stacked improved/tied/harmed bars", "lambda=0.25; q=0,0.5,1; q=0.5 RANDOM_SIGN_CONTROL; percentages among non-ties", "Do not label random-sign control as baseline or hide ties/harm"),
        "C": ("separate development and confirmation facets", "DEVELOPMENT n=2455; CONFIRMATION prespecified reaction holdout n=838; separate rhos", "Do not pool cohorts or describe bootstrap as biological replication"),
        "D": ("16-point context dot/forest-style display with finite n", "all four algorithms and all four RNA contexts; retain RIPTiDe setx2,setx3 n=39 rho=0.27561031478816606", "Do not omit or reorder by rho"),
        "E": ("small specificity and pathway robustness comparison", "primary usefulness rho=0.8351435823183292; information advantage rho=0.34542078458252307; pathway sensitivity values", "No new significance comparison or confirmation classifier"),
    }
    panels = []
    panel_artifacts = {
        "A": "BRIDGEM1_CLAIM_LEDGER.tsv", "B": "BRIDGEM1_PL2B_UTILITY.tsv",
        "C": "BRIDGEM1_DEVELOPMENT_CONFIRMATION.tsv", "D": "BRIDGEM1_CONFIRMATION_CONTEXTS.tsv",
        "E": "BRIDGEM1_RESULTS_FACTS.tsv",
    }
    for spec in core.figure_panels():
        visual, annotations, prohibited = panel_details[spec["panel"]]
        panels.append({"panel": spec["panel"], "title": spec["title"], "scientific_message": spec["message"], "source_stage": spec["source"], "source_artifact": panel_artifacts[spec["panel"]], "required_metrics": annotations, "visual_form": visual, "mandatory_annotations": annotations, "prohibited_visual_interpretation": prohibited})
    blobs["BRIDGEM1_FIGURE_PANEL_PLAN.tsv"] = tsv_bytes(tuple(panels[0]), panels)
    facts = []
    def add(fid: str, stage: str, role: str, metric: str, exact: object, display: object, unit: str, source: str, field: str, caveat: str = "") -> None:
        facts.append({"fact_id": fid, "stage": stage, "analysis_role": role, "metric": metric, "exact_value": exact, "manuscript_display": display, "unit_or_denominator": unit, "source_file": source, "source_field": field, "caveat": caveat})
    for row in utility:
        q = float(row["reliability_q"])
        prefix = f"PL2B_Q{q:.2f}"
        for field in UTILITY_FIELDS:
            display = f"{100 * round(float(row[field]), 4):.2f}%" if field.endswith("fraction_of_non_ties") else row[field]
            add(f"{prefix}_{field.upper()}", "PL2B", "PRIMARY_CONTEXT", field, row[field], display, "non-tie cases" if field.endswith("fraction_of_non_ties") else "source units", source_path("PL2B", "BRIDGEPL2B_CONTEXT_SUMMARY.tsv"), field, "q=0.5 is RANDOM_SIGN_CONTROL" if q == 0.5 else "")
    for cohort, row in zip(("DEVELOPMENT", "CONFIRMATION"), dc_rows):
        for field in ("finite_reactions", "eta2_usefulness_rho", "positive_contexts", "context_count", "minimum_finite_context_count"):
            if row[field] != "":
                add(f"{cohort}_{field.upper()}", "PL2D_H" if cohort == "DEVELOPMENT" else "PL2D", "PRIMARY", field, row[field], row[field], "reactions" if "reactions" in field or "count" in field else "rho or contexts", row["source_artifact"], field, "Separate populations; no pooling")
    for field in ("rho_median", "rho_q025", "rho_q975"):
        add(f"CONFIRMATION_BOOTSTRAP_{field.upper()}", "PL2D", "SUPPORTIVE", field, supportive["bootstrap"][field], supportive["bootstrap"][field], "rho", source_path("PL2D", "BRIDGEPL2D_SUPPORTIVE_RESULTS.json"), f"bootstrap.{field}", "Computational reaction-level stability, not biological replication")
    supportive_metrics = (
        ("information_advantage_rho", "eta2_vs_mean_information_advantage"),
        ("partial_coverage_rho", "partial_eta2_vs_usefulness_controlling_non_tie_coverage"),
        ("directional_entropy_rho", "directional_entropy3_vs_usefulness"),
        ("dominant_sign_mass_rho", "dominant_sign_mass_vs_usefulness"),
        ("supported_sign_state_count_rho", "n_supported_sign_states_vs_usefulness"),
    )
    for metric, key in supportive_metrics:
        add(f"CONFIRMATION_{metric.upper()}", "PL2D", "SUPPORTIVE", metric, supportive[key]["rho"], supportive[key]["rho"], "rho", source_path("PL2D", "BRIDGEPL2D_SUPPORTIVE_RESULTS.json"), f"{key}.rho", "Supportive, not a confirmation gate")
    for key, result in supportive["pathway_sensitivity"].items():
        add(f"PATHWAY_{key.upper()}_RHO", "PL2D", "SUPPORTIVE", key + "_rho", result["rho"], result["rho"], "rho", source_path("PL2D", "BRIDGEPL2D_SUPPORTIVE_RESULTS.json"), f"pathway_sensitivity.{key}.rho", "Supportive pathway sensitivity")
        add(f"PATHWAY_{key.upper()}_N", "PL2D", "SUPPORTIVE", key + "_finite_reactions", result["finite_reactions"], result["finite_reactions"], "reactions", source_path("PL2D", "BRIDGEPL2D_SUPPORTIVE_RESULTS.json"), f"pathway_sensitivity.{key}.finite_reactions")
    for index, row in enumerate(contexts, 1):
        for field in ("rho", "finite_primary_reaction_pairs"):
            add(f"CONTEXT_{index:02d}_{field.upper()}", "PL2D", "PRIMARY_CONTEXT", field, row[field], row[field], "rho" if field == "rho" else "reactions", source_path("PL2D", "BRIDGEPL2D_CONTEXT_RESULTS.tsv"), field, row["algorithm"] + ";" + row["rna_context_key"])
    for row in limitations:
        add(f"LIMITATION_{row['metric'].upper()}", "PL2D", "LIMITATION", row["metric"], row["exact_value"], row["manuscript_display"], row["unit_or_denominator"], row["source_artifact"], row["source_field"], "Directional truth is conditional; computational units are not biological replicates")
    blobs["BRIDGEM1_RESULTS_FACTS.tsv"] = tsv_bytes(tuple(facts[0]), facts)
    if set(blobs) != set(ARTIFACTS):
        raise RuntimeError("M1 artifact inventory mismatch")
    manifest = {
        "schema": core.SCHEMA, "status": core.TERMINAL_STATUS,
        "source_sha256": audit["source_sha256"],
        "patch_sha256": PATCH_SHA256,
        "m1_core_sha256": audit["m1_core_sha256"],
        "adapter_sha256": audit["adapter_sha256"],
        "contract_document_sha256": audit["contract_document_sha256"],
        "artifact_sha256": {name: hashlib.sha256(blobs[name]).hexdigest() for name in ARTIFACTS},
        "row_counts": {"claim_ledger": len(claims), "results_facts": len(facts), "figure_panel_plan": len(panels), "pl2b_utility": len(utility_rows), "development_confirmation": len(dc_rows), "confirmation_contexts": len(contexts), "limitations": len(limitations)},
        **{flag: False for flag in FALSE_FLAGS},
    }
    blobs["BRIDGEM1_MANIFEST.json"] = json_bytes(manifest)
    return blobs


def produce(output: Path = DEFAULT_OUTPUT) -> dict[str, object]:
    if output.is_symlink() or output.resolve() == OUTPUTS.resolve() or not output.resolve().is_relative_to(OUTPUTS.resolve()) or output.parent.resolve() != OUTPUTS.resolve():
        raise RuntimeError("M1 output must be a direct in-repository outputs directory")
    blobs = build_artifacts()
    if output.exists():
        if not output.is_dir() or {p.name for p in output.iterdir()} != set(blobs):
            raise RuntimeError("existing M1 output inventory mismatch")
        for name, expected in blobs.items():
            if safe_file(output / name).read_bytes() != expected:
                raise RuntimeError(f"existing M1 output mutation: {name}")
        return {"status": "NO_OP_EXISTING_IDENTICAL_BRIDGE_M1", "manifest": json.loads(blobs["BRIDGEM1_MANIFEST.json"])}
    with tempfile.TemporaryDirectory(prefix=".dmi_bridge_m1_", dir=OUTPUTS) as temp:
        stage = Path(temp) / "bundle"
        stage.mkdir()
        for name, data in blobs.items():
            (stage / name).write_bytes(data)
        os.replace(stage, output)
    return {"status": core.TERMINAL_STATUS, "manifest": json.loads(blobs["BRIDGEM1_MANIFEST.json"])}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = produce(args.output)
    print(json.dumps({"status": result["status"], "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
