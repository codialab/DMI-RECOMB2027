#!/usr/bin/env python3
"""Publish frozen A2 robustness facts as deterministic manuscript evidence."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import tempfile
from pathlib import Path

try:
    from scripts import dmi_bridge_m2_a2_robustness_core_v1 as core
except ModuleNotFoundError:
    import dmi_bridge_m2_a2_robustness_core_v1 as core

ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"
DEFAULT_OUTPUT = OUTPUTS / "dmi_bridge_m2_a2_robustness_evidence_v1"
PATCH = ROOT / "dmi_bridge_m2_a2_robustness_evidence_foundation_v1.patch"
CONTRACT = ROOT / "docs/DMI_BRIDGE_M2_A2_ROBUSTNESS_EVIDENCE_CONSOLIDATION.md"
CORE = ROOT / "scripts/dmi_bridge_m2_a2_robustness_core_v1.py"
PATCH_SHA256 = "1467481a6a0ccc77c1e7e84f748e440f7c426a24cca69a4bfc7207e45c8e6d63"
SOURCES = {
    "M1": ("dmi_bridge_m1_manuscript_evidence_v1", "BRIDGEM1_MANIFEST.json", core.M1_MANIFEST_SHA256, "BRIDGE_M1_MANUSCRIPT_EVIDENCE_FROZEN"),
    "A20": ("dmi_bridge_a20_dual_anchor_qualification_v1", "BRIDGEA20_MANIFEST.json", core.A20_MANIFEST_SHA256, "BRIDGE_A20_PARTIAL_QUALIFICATION"),
    "A21": ("dmi_bridge_a21_dual_anchor_geometry_v1", "BRIDGEA21_MANIFEST.json", core.A21_MANIFEST_SHA256, "BRIDGE_A21_DUAL_ANCHOR_GEOMETRY_FROZEN"),
    "A22": ("dmi_bridge_a22_dual_anchor_sign_only_utility_v1", "BRIDGEA22_MANIFEST.json", core.A22_MANIFEST_SHA256, "BRIDGE_A22_SIGN_ONLY_UTILITY_FROZEN"),
    "A23": ("dmi_bridge_a23_a1_a2_synthesis_v1", "BRIDGEA23_MANIFEST.json", core.A23_MANIFEST_SHA256, "BRIDGE_A23_A1_A2_SYNTHESIS_FROZEN"),
}
ARTIFACTS = (
    "BRIDGEM2_SOURCE_AUDIT.json", "BRIDGEM2_CLAIM_LEDGER.tsv",
    "BRIDGEM2_RESULTS_FACTS.tsv", "BRIDGEM2_A1_A2_UTILITY.tsv",
    "BRIDGEM2_A1_A2_GEOMETRY.tsv", "BRIDGEM2_A1_A2_PERSISTENCE.tsv",
    "BRIDGEM2_CONTROLS.tsv", "BRIDGEM2_LIMITATIONS.tsv",
    "BRIDGEM2_FIGURE_PANEL_PLAN.tsv",
)
FALSE_FLAGS = (
    "new_scientific_analysis_performed", "new_inferential_statistics_computed",
    "upstream_outputs_modified", "solver_invoked", "fva_invoked",
    "sampling_invoked", "reconstruction_invoked", "reaction_ranking_performed",
    "figure_rendered", "a2_g_used",
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


def checked(path: Path, expected: str) -> None:
    if sha(path) != expected:
        raise RuntimeError(f"SHA256 mismatch: {path}")


def json_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def tsv_bytes(fields: tuple[str, ...], rows: list[dict]) -> bytes:
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
    writer.writeheader()
    writer.writerows(rows)
    return out.getvalue().encode("utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with safe_file(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def source_path(stage: str, name: str) -> str:
    return f"outputs/{SOURCES[stage][0]}/{name}"


def verify_sources() -> tuple[dict, dict]:
    checked(PATCH, PATCH_SHA256)
    manifests = {}
    hashes = {}
    artifact_counts = {}
    source_counts = {}
    for stage, (directory, name, digest, status) in SOURCES.items():
        path = OUTPUTS / directory / name
        checked(path, digest)
        hashes[source_path(stage, name)] = digest
        manifest = json.loads(safe_file(path).read_text(encoding="utf-8"))
        if manifest.get("status") != status:
            raise RuntimeError(f"{stage} terminal status mismatch")
        declared = manifest.get("artifact_sha256")
        if not isinstance(declared, dict) or not declared:
            raise RuntimeError(f"{stage} artifact inventory missing")
        for artifact, expected in declared.items():
            if Path(artifact).name != artifact:
                raise RuntimeError(f"{stage} unsafe artifact name")
            checked(OUTPUTS / directory / artifact, expected)
            hashes[source_path(stage, artifact)] = expected
        artifact_counts[stage] = len(declared)
        source_map = manifest.get("source_sha256", {})
        if not isinstance(source_map, dict):
            raise RuntimeError(f"{stage} source inventory invalid")
        for role, entry in source_map.items():
            if stage == "A20":
                relative, expected = entry["relative_path"], entry["sha256"]
            else:
                relative, expected = role, entry
            if Path(relative).is_absolute() or ".." in Path(relative).parts:
                raise RuntimeError(f"{stage} unsafe source path")
            checked(ROOT / relative, expected)
            hashes[relative] = expected
        source_counts[stage] = len(source_map)
        manifests[stage] = manifest
    if artifact_counts != {"M1": 8, "A20": 6, "A21": 25, "A22": 137, "A23": 9}:
        raise RuntimeError("manifest artifact inventory changed")
    if source_counts != {"M1": 164, "A20": 19, "A21": 39, "A22": 0, "A23": 11}:
        raise RuntimeError("manifest source inventory changed")
    a20 = manifests["A20"]
    status = json.loads(safe_file(ROOT / source_path("A20", "BRIDGEA20_STATUS.json")).read_text())
    if (a20.get("qualification") != {"A2-G": "NOT_QUALIFIED", "A2-L": "QUALIFIED"}
            or a20.get("a2_g_outcome_accessed") is not False
            or status.get("arm_status") != a20["qualification"]
            or status.get("a2_g_outcome_accessed") is not False):
        raise RuntimeError("A2 qualification mismatch")
    if (manifests["A21"].get("feature_logical_content_sha256") != "50bc3e47239e4b22d0931246f4fa0d5d126a5b69226d47cb342cb8a5b3ee99d2"
            or manifests["A21"].get("a2_g_used") is not False
            or manifests["A22"].get("case_logical_content_sha256") != "70321eb7d322f06a3c1e29c048d08f9b4c55d291763b5194f6d626325ed786f4"
            or manifests["A22"].get("a2_g_used") is not False):
        raise RuntimeError("A2 logical identity mismatch")
    for name, digest in core.A23_ARTIFACT_SHA256.items():
        checked(ROOT / source_path("A23", name), digest)
    core.validate_a23_manifest(manifests["A23"])
    return manifests, {
        "schema": "bridge.m2.source_audit.v1", "status": "PASS",
        "source_sha256": dict(sorted(hashes.items())),
        "manifest_sha256": {stage: spec[2] for stage, spec in SOURCES.items()},
        "manifest_artifact_counts": artifact_counts,
        "manifest_source_counts": source_counts,
        "patch_sha256": sha(PATCH), "core_sha256": sha(CORE),
        "adapter_sha256": sha(Path(__file__)), "contract_document_sha256": sha(CONTRACT),
    }


def frozen_rows(manifests: dict) -> tuple[dict, list[dict], list[dict], list[dict], list[dict], dict]:
    def rows(name: str) -> list[dict]:
        return read_tsv(ROOT / source_path("A23", name))

    qc = json.loads(safe_file(ROOT / source_path("A23", "BRIDGEA23_QC.json")).read_text())
    utility = rows("BRIDGEA23_UTILITY_COMPARISON.tsv")
    controls = rows("BRIDGEA23_CONTROL_SUMMARY.tsv")
    geometry = rows("BRIDGEA23_GEOMETRY_COMPARISON.tsv")
    persistence = rows("BRIDGEA23_GEOMETRY_UTILITY_PERSISTENCE.tsv")
    population = rows("BRIDGEA23_MATCHED_POPULATION_AUDIT.tsv")
    if (len(utility), len(controls), len(geometry), len(persistence), len(population)) != (3, 9, 24, 42, 836):
        raise RuntimeError("A2.3 source row inventory mismatch")
    if {r["scope"] for r in utility} != {"ALL", "DEVELOPMENT", "CONFIRMATION_HOLDOUT"}:
        raise RuntimeError("utility scope inventory mismatch")
    all_utility = next(r for r in utility if r["scope"] == "ALL")
    core.validate_utility_all(all_utility)
    core.validate_control_rows(controls)
    core.validate_persistence_rows(persistence)
    core.validate_qc(qc)
    if qc.get("audit_only_truth_pairs_excluded") != 49 or qc.get("a2_g_used") is not False:
        raise RuntimeError("A2.3 QC exclusion/arm mismatch")
    features = (
        "sign_magnitude_eta2", "sign_magnitude_correlation", "directional_entropy3",
        "dominant_sign_mass", "n_supported_sign_states", "anchor_width80_contraction",
        "joint_ess", "abs_delta_anchor_target_correlation",
    )
    scopes = ("ALL_EVALUATION_REACTION", "DEVELOPMENT", "CONFIRMATION_HOLDOUT")
    if {(r["feature"], r["scope"]) for r in geometry} != {(f, s) for f in features for s in scopes}:
        raise RuntimeError("geometry feature/scope inventory mismatch")
    for row in geometry:
        expected = 1671600 if row["scope"] == "ALL_EVALUATION_REACTION" else 3134 if row["scope"] == "DEVELOPMENT" else 1045
        if int(row["candidate_rows"]) != expected or not 0 <= int(row["finite_paired_rows"]) <= expected:
            raise RuntimeError("geometry denominator mismatch")
    eta2 = next(r for r in geometry if r["feature"] == "sign_magnitude_eta2" and r["scope"] == "ALL_EVALUATION_REACTION")
    if int(eta2["finite_paired_rows"]) != core.EXPECTED_ETA2_FINITE_PAIRED:
        raise RuntimeError("eta2 finite denominator mismatch")
    return qc, utility, controls, geometry, persistence, manifests["A23"]


def build_artifacts() -> dict[str, bytes]:
    manifests, audit = verify_sources()
    qc, utility, controls, geometry, persistence, a23 = frozen_rows(manifests)
    blobs = {"BRIDGEM2_SOURCE_AUDIT.json": json_bytes(audit)}
    def table(name: str, records: list[dict]) -> None:
        blobs[name] = tsv_bytes(tuple(records[0]), records)

    table("BRIDGEM2_A1_A2_UTILITY.tsv", utility)
    table("BRIDGEM2_CONTROLS.tsv", controls)
    table("BRIDGEM2_A1_A2_GEOMETRY.tsv", geometry)
    primary_persistence = [r for r in persistence if r["level"] == "PL2D_REACTION" and r["context"] == "ALL" and r["feature"] == "sign_magnitude_eta2"]
    primary_persistence.sort(key=lambda r: ("DEVELOPMENT", "CONFIRMATION_HOLDOUT").index(r["split"]))
    table("BRIDGEM2_A1_A2_PERSISTENCE.tsv", primary_persistence)
    claim_sources = {
        "R1_MATCHED_CORRECT_DIRECTION_UTILITY": ("BRIDGEA23_UTILITY_COMPARISON.tsv", "ALL: a1_correct_mean_gain;a2_correct_mean_gain;paired_correct_delta"),
        "R2_DIRECTIONAL_SPECIFICITY": ("BRIDGEA23_UTILITY_COMPARISON.tsv;BRIDGEA23_CONTROL_SUMMARY.tsv", "ALL: a1_specificity_mean;a2_specificity_mean;paired_specificity_delta;wrong_q0;analytic_random_q0_5"),
        "R3_GEOMETRY_UTILITY_PERSISTENCE": ("BRIDGEA23_GEOMETRY_UTILITY_PERSISTENCE.tsv;BRIDGEA23_QC.json", "PL2D_REACTION/ALL/sign_magnitude_eta2: DEVELOPMENT;CONFIRMATION_HOLDOUT;bootstrap"),
        "R4_RESIDUAL_GEOMETRY_SHIFT": ("BRIDGEA23_GEOMETRY_COMPARISON.tsv", "eight features x three scopes;candidate_rows;finite_paired_rows"),
        "L1_POST_FREEZE_STATUS": ("BRIDGEA23_ANALYSIS_CONTRACT.json", "post-freeze robustness/generalization"),
        "L2_A2L_ONLY": ("BRIDGEA20_STATUS.json", "A2-L QUALIFIED;A2-G NOT_QUALIFIED;a2_g_outcome_accessed=false"),
    }
    claims = []
    for spec in core.claim_specs():
        artifact, fields = claim_sources[spec["claim_id"]]
        boundary = spec["boundary"]
        if spec["claim_id"] == "R1_MATCHED_CORRECT_DIRECTION_UTILITY":
            boundary += " No new p-value or significance claim."
        if spec["claim_id"] == "R2_DIRECTIONAL_SPECIFICITY":
            boundary += " Do not imply information gain from random data."
        if spec["claim_id"] == "R3_GEOMETRY_UTILITY_PERSISTENCE":
            boundary += " Do not claim A2 improves the predictor; bootstrap is computational reaction-level stability only."
        claims.append({"claim_id": spec["claim_id"], "claim_level": spec["claim_level"], "manuscript_safe_statement": spec["statement"], "source_artifacts": artifact, "source_fields": fields, "interpretation_boundary": boundary})
    table("BRIDGEM2_CLAIM_LEDGER.tsv", claims)
    limitations = [
        ("L1_POST_FREEZE", "A2 is post-freeze robustness/generalization, not untouched confirmation.", "BRIDGEA23_ANALYSIS_CONTRACT.json"),
        ("L2_A2L_ONLY", "Only A2-L qualified; A2-G failed its operator/provenance gate and supplies no outcome.", "BRIDGEA20_STATUS.json"),
        ("L3_EXCLUDED_PAIRS", "49 A1-non-evaluable truth pairs remain outside the primary matched population.", "BRIDGEA23_QC.json"),
        ("L4_TRUTH_TIES", "Directional gains are conditional on 1440313 non-tie case keys; 2053331 truth ties remain.", "BRIDGEA23_UTILITY_COMPARISON.tsv"),
        ("L5_FINITE_GEOMETRY", "Finite support varies by descriptor; eta2 has 573882 finite paired rows of 1671600 candidates; undefined is not zero.", "BRIDGEA23_GEOMETRY_COMPARISON.tsv"),
        ("L6_COMPUTATIONAL_UNITS", "Computational reactions and evaluations are not biological replicates; bootstrap is reaction-level stability only.", "BRIDGEA23_QC.json"),
    ]
    table("BRIDGEM2_LIMITATIONS.tsv", [{"limitation_id": k, "statement": v, "source_artifact": s} for k, v, s in limitations])
    panels = [
        ("A", "Matched design", "A1 HEX1 versus A2-L HEX1+LDH_L; 836 pairs and 4179 weak targets", "BRIDGEM2_RESULTS_FACTS.tsv", "Do not infer A2-G outcomes"),
        ("B", "Utility and specificity", "Correct, wrong, q=0.5 random-sign gains and correct-minus-wrong contrast; 1440313 non-ties", "BRIDGEM2_A1_A2_UTILITY.tsv;BRIDGEM2_CONTROLS.tsv", "Do not portray random sign as directional information"),
        ("C", "Geometry-to-utility persistence", "Separate development and confirmation rhos and finite reaction counts", "BRIDGEM2_A1_A2_PERSISTENCE.tsv", "No significance stars or independent-confirmation claim"),
        ("D", "Residual geometry", "All eight original descriptors across three scopes with finite denominators", "BRIDGEM2_A1_A2_GEOMETRY.tsv", "No universal more-constrained or more-ambiguous label"),
        ("E", "Provenance and finite support", "Post-freeze A2-L only; A2-G gate; 49 exclusions; ties; finite support; computational units", "BRIDGEM2_LIMITATIONS.tsv", "No main-versus-supplement placement based on favorability"),
    ]
    table("BRIDGEM2_FIGURE_PANEL_PLAN.tsv", [{"panel": p, "title": t, "required_content": c, "source_artifacts": s, "interpretation_boundary": b} for p, t, c, s, b in panels])
    facts = []
    def add(fid: str, value: object, source: str, field: str, unit: str = "", caveat: str = "") -> None:
        facts.append({"fact_id": fid, "exact_value": value, "manuscript_display": value, "unit_or_denominator": unit, "source_artifact": source, "source_field": field, "caveat": caveat})
    population = a23["case_population"]
    for fid, value, field, unit in (
        ("MATCHED_TRUTH_PAIRS", population["truth_pairs"], "case_population.truth_pairs", "truth pairs"),
        ("REACTIONS", population["reactions"], "case_population.reactions", "reactions"),
        ("MATCHED_CASE_KEYS", population["case_keys"], "case_population.case_keys", "case keys"),
        ("NON_TIE_CASE_KEYS", next(r for r in utility if r["scope"] == "ALL")["non_tie_pairs"], "ALL.non_tie_pairs", "case keys"),
        ("TRUTH_TIE_CASE_KEYS", next(r for r in utility if r["scope"] == "ALL")["truth_ties"], "ALL.truth_ties", "case keys"),
        ("GEOMETRY_KEYS", qc["geometry_keys"], "geometry_keys", "evaluation x reaction keys"),
        ("ETA2_FINITE_PAIRED", eta2_finite(geometry), "ALL_EVALUATION_REACTION.sign_magnitude_eta2.finite_paired_rows", "paired keys"),
        ("EXCLUDED_TRUTH_PAIRS", qc["audit_only_truth_pairs_excluded"], "audit_only_truth_pairs_excluded", "truth pairs"),
    ):
        add(fid, value, "BRIDGEA23_MANIFEST.json" if field.startswith("case_population") else "BRIDGEA23_UTILITY_COMPARISON.tsv" if field.startswith("ALL.") else "BRIDGEA23_GEOMETRY_COMPARISON.tsv" if fid == "ETA2_FINITE_PAIRED" else "BRIDGEA23_QC.json", field, unit)
    all_row = next(r for r in utility if r["scope"] == "ALL")
    for field in core.EXPECTED_ALL_UTILITY:
        add(field.upper(), all_row[field], "BRIDGEA23_UTILITY_COMPARISON.tsv", f"ALL.{field}", "mean gain on non-ties", "Descriptive matched robustness")
    for arm in ("wrong_q0", "analytic_random_q0_5"):
        row = next(r for r in controls if r["scope"] == "ALL" and r["arm"] == arm)
        for field in ("a1_mean_gain", "a2_mean_gain"):
            add(f"{arm.upper()}_{field.upper()}", row[field], "BRIDGEA23_CONTROL_SUMMARY.tsv", f"ALL.{arm}.{field}", "mean gain on non-ties")
    for row in primary_persistence:
        for field in ("a1_rho", "a2_rho", "a1_finite", "a2_finite"):
            add(f"{row['split']}_{field.upper()}", row[field], "BRIDGEA23_GEOMETRY_UTILITY_PERSISTENCE.tsv", f"PL2D_REACTION.{row['split']}.ALL.sign_magnitude_eta2.{field}", "rho" if field.endswith("rho") else "reactions")
    for field in ("replicates", "finite_replicates", "rho_median", "rho_q025", "rho_q975"):
        add(f"A2_BOOTSTRAP_{field.upper()}", qc["bootstrap"][field], "BRIDGEA23_QC.json", f"bootstrap.{field}", "replicates" if "replicates" in field else "rho", "Computational reaction-level stability only")
    if len(facts) != 31:
        raise RuntimeError("M2 results fact inventory mismatch")
    table("BRIDGEM2_RESULTS_FACTS.tsv", facts)
    if set(blobs) != set(ARTIFACTS):
        raise RuntimeError("M2 artifact inventory mismatch")
    manifest = {
        "schema": core.SCHEMA, "status": core.TERMINAL_STATUS,
        "source_sha256": audit["source_sha256"],
        "patch_sha256": audit["patch_sha256"], "core_sha256": audit["core_sha256"],
        "adapter_sha256": audit["adapter_sha256"],
        "contract_document_sha256": audit["contract_document_sha256"],
        "artifact_sha256": {name: hashlib.sha256(blobs[name]).hexdigest() for name in ARTIFACTS},
        "row_counts": {"claim_ledger": len(claims), "results_facts": len(facts), "a1_a2_utility": len(utility), "a1_a2_geometry": len(geometry), "a1_a2_persistence": len(primary_persistence), "controls": len(controls), "limitations": len(limitations), "figure_panel_plan": len(panels)},
        **{flag: False for flag in FALSE_FLAGS},
    }
    blobs["BRIDGEM2_MANIFEST.json"] = json_bytes(manifest)
    return blobs


def eta2_finite(geometry: list[dict]) -> str:
    return next(r["finite_paired_rows"] for r in geometry if r["scope"] == "ALL_EVALUATION_REACTION" and r["feature"] == "sign_magnitude_eta2")


def produce(output: Path = DEFAULT_OUTPUT) -> dict:
    if output.is_symlink() or output.resolve() != DEFAULT_OUTPUT.resolve() or output.parent.resolve() != OUTPUTS.resolve():
        raise RuntimeError("M2 output must be the canonical direct in-repository outputs directory")
    blobs = build_artifacts()
    if output.exists():
        if not output.is_dir() or {p.name for p in output.iterdir()} != set(blobs):
            raise RuntimeError("existing M2 output inventory mismatch")
        for name, expected in blobs.items():
            if safe_file(output / name).read_bytes() != expected:
                raise RuntimeError(f"existing M2 output mutation: {name}")
        return {"status": "NO_OP_EXISTING_IDENTICAL_BRIDGE_M2", "manifest": json.loads(blobs["BRIDGEM2_MANIFEST.json"])}
    with tempfile.TemporaryDirectory(prefix=".dmi_bridge_m2_", dir=OUTPUTS) as temp:
        stage = Path(temp) / "bundle"
        stage.mkdir()
        for name, data in blobs.items():
            (stage / name).write_bytes(data)
        os.replace(stage, output)
    return {"status": core.TERMINAL_STATUS, "manifest": json.loads(blobs["BRIDGEM2_MANIFEST.json"])}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = produce(args.output)
    print(json.dumps({"status": result["status"], "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
