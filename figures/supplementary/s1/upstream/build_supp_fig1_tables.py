#!/usr/bin/env python3
"""Build S1 diagnostics from frozen candidate-weighting artifacts only."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "reproduced/derived/s1"
DATA = OUT / "data"
METHODS = ("iMAT", "GIMME", "CORDA", "RIPTiDe")
CONTEXTS = (
    "training_samples=setx1,setx2",
    "training_samples=setx1,setx3",
    "training_samples=setx2,setx3",
    "training_samples=setx1,setx2,setx3",
)
TUMORS = ("CT2A", "GL261")
ANCHORS = ("A1", "A2")
KEY = ["anchor", "method", "rna_context_key", "tumor", "mouse_id", "candidate_index"]
TOL = 1e-10

SOURCES = {
    "bio0_manifest": "outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_MANIFEST.json",
    "a1_weights": "outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_FOUR_ARM_WEIGHTS.tsv.xz",
    "fig4_manifest": "data/figure_inputs/fig4/fig4_data_manifest.json",
    "a20_manifest": "outputs/dmi_bridge_a20_dual_anchor_qualification_v1/BRIDGEA20_MANIFEST.json",
    "a2_weights": "outputs/dmi_bridge_a20_dual_anchor_qualification_v1/BRIDGEA20_DUAL_ANCHOR_WEIGHTS.tsv.xz",
    "pl1_manifest": "outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_MANIFEST.json",
    "evaluation_registry": "outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_EVALUATION_REGISTRY.tsv",
    "pl2a_manifest": "outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_MANIFEST.json",
    "truth_selection": "outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_TRUTH_SELECTION.tsv",
    "truth_pair_registry": "outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_TRUTH_PAIR_REGISTRY.tsv",
    "pl2b_manifest": "outputs/dmi_bridge_pl2b_sign_only_utility_v1/BRIDGEPL2B_MANIFEST.json",
    "a1_matched_pair_registry": "outputs/dmi_bridge_pl2b_sign_only_utility_v1/BRIDGEPL2B_DISTINCT_TRUTH_PAIR_REGISTRY.tsv",
    "a22_manifest": "outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_MANIFEST.json",
    "matched_pair_registry": "outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_MATCHED_TRUTH_PAIR_REGISTRY.tsv",
    "holdout_support": "outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_TRUTH_SUPPORT_AUDIT.tsv",
}


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def source_path(relative: str) -> Path:
    path = ROOT / relative
    require(path.is_file() and not path.is_symlink(), f"missing or non-regular frozen input: {relative}")
    resolved = path.resolve(strict=True)
    require(resolved.is_relative_to(ROOT), f"input escapes repository: {relative}")
    return resolved


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_manifest(manifest_key: str, selected_keys: list[str]) -> dict:
    manifest_path = source_path(SOURCES[manifest_key])
    manifest = json.loads(manifest_path.read_text())
    inventory = manifest.get("artifact_sha256", {})
    for key in selected_keys:
        relative = SOURCES[key]
        name = Path(relative).name
        expected = inventory.get(name)
        require(expected is not None, f"{name} is not checksum-listed in {manifest_path}")
        require(sha256(source_path(relative)) == expected, f"frozen artifact checksum mismatch: {relative}")
    return manifest


def effective_n(weights: np.ndarray) -> float:
    total = float(weights.sum())
    require(total > 0 and np.isfinite(total), "weight vector has no finite positive mass")
    normalized = weights / total
    return float(1.0 / np.dot(normalized, normalized))


def load_weights() -> pd.DataFrame:
    a1_path = source_path(SOURCES["a1_weights"])
    a2_path = source_path(SOURCES["a2_weights"])
    a1 = pd.read_csv(a1_path, sep="\t", compression="xz")
    a1 = a1.loc[a1.arm.eq("strong_anchor_baseline")].copy()
    require(len(a1) == 3200, "A1 strong_anchor_baseline should have 3,200 candidate-weight rows")
    a1 = a1.rename(columns={"algorithm": "method", "sample_index": "candidate_index"})
    a1["anchor"] = "A1"
    a2 = pd.read_csv(a2_path, sep="\t", compression="xz")
    a2 = a2.loc[a2.arm.eq("A2-L")].copy()
    require(len(a2) == 3200, "A2-L frozen weights should have 3,200 candidate-weight rows")
    a2 = a2.rename(columns={"algorithm": "method", "sample_index": "candidate_index"})
    a2["anchor"] = "A2"
    frame = pd.concat([
        a1[["anchor", "method", "rna_context_key", "tumor", "mouse_id", "candidate_index", "weight"]],
        a2[["anchor", "method", "rna_context_key", "tumor", "mouse_id", "candidate_index", "weight"]],
    ], ignore_index=True)
    frame["candidate_index"] = pd.to_numeric(frame.candidate_index, errors="raise").astype(int)
    frame["global_weight"] = pd.to_numeric(frame.weight, errors="coerce")
    require(frame["global_weight"].notna().all() and np.isfinite(frame.global_weight.to_numpy()).all(), "nonfinite/missing candidate weights")
    require((frame.global_weight >= 0).all(), "negative candidate weight")
    require(not frame.duplicated(KEY).any(), "duplicate candidate-weight keys")
    require(set(frame.anchor) == set(ANCHORS), "unexpected anchor values")
    require(set(frame.method) == set(METHODS), "unexpected methods in candidate weights")
    require(set(frame.rna_context_key) == set(CONTEXTS), "unexpected RNA contexts in candidate weights")
    require(set(frame.tumor) == set(TUMORS), "unexpected tumors in candidate weights")
    expected_mice = {"CT2A": {"C1", "C2", "C3", "C4", "C5"}, "GL261": {"G1", "G2", "G3", "G4", "G5"}}
    require({t: set(frame.loc[frame.tumor.eq(t), "mouse_id"]) for t in TUMORS} == expected_mice,
            "unexpected DMI mouse inventory")
    require(frame.groupby(KEY[:-1], observed=True).size().eq(20).all(), "a conditional stratum does not contain exactly 20 candidates")
    require(len(frame) == 6400 and frame.groupby(KEY[:-1], observed=True).ngroups == 320, "expected 6,400 rows across 320 conditional strata")
    global_sums = frame.groupby(["anchor", "tumor", "mouse_id"], observed=True).global_weight.sum()
    require(np.allclose(global_sums.to_numpy(), 1.0, atol=TOL, rtol=0), "global 320-candidate weights are not normalized")
    frame["conditional_weight"] = frame.global_weight / frame.groupby(KEY[:-1], observed=True).global_weight.transform("sum")
    frame = frame.drop(columns="weight")
    require(np.allclose(frame.groupby(KEY[:-1], observed=True).conditional_weight.sum().to_numpy(), 1.0, atol=TOL, rtol=0),
            "conditional 20-candidate weights are not normalized")
    return frame.sort_values(KEY, kind="stable").reset_index(drop=True)


def build() -> dict:
    DATA.mkdir(parents=True, exist_ok=True)
    manifests = {
        "bio0": json.loads(source_path(SOURCES["bio0_manifest"]).read_text()),
        "a20": verify_manifest("a20_manifest", ["a2_weights"]),
        "pl1": verify_manifest("pl1_manifest", ["evaluation_registry"]),
        "pl2a": verify_manifest("pl2a_manifest", ["truth_selection", "truth_pair_registry"]),
        "pl2b": verify_manifest("pl2b_manifest", ["a1_matched_pair_registry"]),
        "a22": verify_manifest("a22_manifest", ["matched_pair_registry", "holdout_support"]),
    }
    fig4_manifest = json.loads(source_path(SOURCES["fig4_manifest"]).read_text())
    fig4_a1_hash = fig4_manifest.get("source_files_used", {}).get(SOURCES["a1_weights"])
    require(fig4_a1_hash is not None and sha256(source_path(SOURCES["a1_weights"])) == fig4_a1_hash,
            "A1 weight checksum does not match frozen Figure 4 source provenance")
    manifests["fig4"] = fig4_manifest
    weights = load_weights()
    weights.to_csv(DATA / "candidate_weights.tsv", sep="\t", index=False, float_format="%.17g")

    # One row per conditional context. ESS uses the 20 normalized candidate weights.
    stratum_rows = []
    for keys, group in weights.groupby(KEY[:-1], sort=True, observed=True):
        anchor, method, context, tumor, mouse = keys
        w = group.conditional_weight.to_numpy(float)
        stratum_rows.append({"anchor": anchor, "method": method, "rna_context_key": context,
                             "tumor": tumor, "mouse_id": mouse, "candidate_count": len(w),
                             "weight_sum": float(w.sum()), "conditional_ess": effective_n(w),
                             "max_candidate_weight": float(w.max()), "positive_candidate_count": int((w > 0).sum()),
                             "zero_candidate_count": int((w == 0).sum())})
    strata = pd.DataFrame(stratum_rows).sort_values(KEY[:-1], kind="stable").reset_index(drop=True)
    require(strata.groupby("anchor", observed=True).size().to_dict() == {"A1": 160, "A2": 160},
            "expected 160 conditional strata per anchor setting")
    require(strata.groupby(["anchor", "tumor", "mouse_id"], observed=True).size().eq(16).all(),
            "each tumor/mouse/anchor pool should contain 16 method × RNA context strata")
    strata.to_csv(DATA / "conditional_stratum_summary.tsv", sep="\t", index=False, float_format="%.17g")

    # Global ESS is over the 320 method × RNA-context × candidate states for each tumor/mouse/anchor.
    global_rows = []
    for keys, group in weights.groupby(["anchor", "tumor", "mouse_id"], sort=True, observed=True):
        anchor, tumor, mouse = keys
        global_rows.append({"anchor": anchor, "tumor": tumor, "mouse_id": mouse,
                            "candidate_count": len(group), "weight_sum": float(group.global_weight.sum()),
                            "global_ess": effective_n(group.global_weight.to_numpy(float)), "ess_target": 20.0})
    global_ess = pd.DataFrame(global_rows)
    require(len(global_ess) == 20 and global_ess.candidate_count.eq(320).all(), "global candidate pool should contain 320 states per anchor/tumor/mouse")
    require(np.allclose(global_ess.global_ess.to_numpy(), 20.0, atol=1e-8, rtol=0),
            "global ESS does not meet the frozen approximately-20 target")
    global_ess.to_csv(DATA / "global_ess_summary.tsv", sep="\t", index=False, float_format="%.17g")

    # Product support uses the product weights over two 20-state condition distributions.
    # Its ESS is the product of the two condition-specific ESS values (400 Cartesian pairs).
    registry = pd.read_csv(source_path(SOURCES["evaluation_registry"]), sep="\t")
    require(len(registry) == 400 and not registry.evaluation_id.duplicated().any(), "PL1 evaluation registry is not the frozen 400-evaluation population")
    require(registry.status.eq("ADMITTED").all(), "PL1 registry includes non-admitted evaluations")
    ess_lookup = strata.set_index(["anchor", "method", "rna_context_key", "tumor", "mouse_id"]).conditional_ess
    contrast_rows = []
    for row in registry.itertuples(index=False):
        for anchor in ANCHORS:
            left = float(ess_lookup.loc[(anchor, row.algorithm, row.rna_context_key, "CT2A", row.ct2a_mouse)])
            right = float(ess_lookup.loc[(anchor, row.algorithm, row.rna_context_key, "GL261", row.gl261_mouse)])
            contrast_rows.append({"evaluation_id": row.evaluation_id, "anchor": anchor, "method": row.algorithm,
                                  "rna_context_key": row.rna_context_key, "ct2a_mouse": row.ct2a_mouse,
                                  "gl261_mouse": row.gl261_mouse, "cartesian_candidate_pairs": 400,
                                  "contrast_product_ess": left * right,
                                  "ct2a_conditional_ess": left, "gl261_conditional_ess": right})
    contrast = pd.DataFrame(contrast_rows)
    require(len(contrast) == 800 and not contrast.duplicated(["evaluation_id", "anchor"]).any(), "contrast support table should have two anchor rows per evaluation")
    contrast.to_csv(DATA / "contrast_product_support.tsv", sep="\t", index=False, float_format="%.17g")

    selection = pd.read_csv(source_path(SOURCES["truth_selection"]), sep="\t")
    require(len(selection) > 0 and selection.holdout_support_status.notna().all()
            and selection.holdout_support_status.str.startswith(("EVALUABLE", "NON_EVALUABLE")).all(),
            "unexpected holdout selection status")
    selection["post_holdout_ess"] = pd.to_numeric(selection.post_holdout_ess, errors="coerce")
    selection["post_holdout_mass"] = pd.to_numeric(selection.post_holdout_mass, errors="coerce")
    selection.to_csv(DATA / "post_holdout_candidate_support.tsv", sep="\t", index=False, float_format="%.17g")
    pairs = pd.read_csv(source_path(SOURCES["matched_pair_registry"]), sep="\t")
    require(len(pairs) == 836 and pairs.truth_pair_id.nunique() == 836, "A22 primary matched truth population must be exactly 836 pairs")
    require(not pairs.truth_pair_id.duplicated().any(), "duplicate matched truth-pair IDs")
    a1_pairs = pd.read_csv(source_path(SOURCES["a1_matched_pair_registry"]), sep="\t")
    a1_pairs = a1_pairs.loc[a1_pairs.truth_pair_id.isin(pairs.truth_pair_id)].copy()
    require(len(a1_pairs) == 836 and set(a1_pairs.truth_pair_id) == set(pairs.truth_pair_id),
            "A1 registry does not contain the exact fixed 836-pair population")
    support_rows = []
    for arm, pair_frame in (("A1", a1_pairs), ("A2", pairs)):
        for row in pair_frame.itertuples(index=False):
            for condition, mouse_col, ess_col, mass_col, status_col in (
                ("CT2A", "ct2a_mouse", "ct2a_post_holdout_ess", "ct2a_post_holdout_mass", "ct2a_holdout_status"),
                ("GL261", "gl261_mouse", "gl261_post_holdout_ess", "gl261_post_holdout_mass", "gl261_holdout_status"),
            ):
                support_rows.append({"anchor": arm, "truth_pair_id": row.truth_pair_id,
                                     "evaluation_id": row.evaluation_id, "method": row.algorithm,
                                     "rna_context_key": row.rna_context_key, "condition": condition,
                                     "mouse_id": getattr(row, mouse_col), "post_holdout_ess": getattr(row, ess_col),
                                     "post_holdout_mass": getattr(row, mass_col),
                                     "holdout_status": getattr(row, status_col)})
    pair_support = pd.DataFrame(support_rows)
    require(len(pair_support) == 3344 and set(pair_support.truth_pair_id) == set(pairs.truth_pair_id),
            "post-holdout support does not contain both conditions and anchors for the fixed 836 pairs")
    require(pair_support.groupby(["anchor", "truth_pair_id"], observed=True).condition.nunique().eq(2).all(),
            "each matched truth pair must have CT2A and GL261 support records for both anchors")
    require(not pair_support.duplicated(["anchor", "truth_pair_id", "condition"]).any(), "duplicate condition support per matched truth pair")
    pair_support["post_holdout_ess"] = pd.to_numeric(pair_support.post_holdout_ess, errors="coerce")
    pair_support["post_holdout_mass"] = pd.to_numeric(pair_support.post_holdout_mass, errors="coerce")
    require(np.isfinite(pair_support[["post_holdout_ess", "post_holdout_mass"]].to_numpy(float)).all(),
            "nonfinite post-holdout support values in the fixed truth population")
    require((pair_support.post_holdout_ess >= 1.0 - 1e-3).all() and (pair_support.post_holdout_ess <= 20.0 + 1e-8).all(),
            "post-holdout ESS outside the valid 20-candidate support range")
    require(pair_support.holdout_status.eq("EVALUABLE").all(), "fixed matched population contains non-evaluable holdout support")
    pair_support.to_csv(DATA / "post_holdout_pair_support.tsv", sep="\t", index=False, float_format="%.17g")

    # Summarize documented exclusions without deleting source rows from the saved diagnostics.
    exclusions = selection.groupby("holdout_support_status", dropna=False).size().to_dict()
    require(len(pairs) == 836, "fixed A1/A2 comparison population mismatch")
    inputs = {}
    for key, relative in SOURCES.items():
        path = source_path(relative)
        inputs[relative] = {"sha256": sha256(path), "bytes": path.stat().st_size,
                            "frozen_manifest_verified": key not in {"bio0_manifest", "a20_manifest", "pl1_manifest", "pl2a_manifest", "a22_manifest"}}
    summary = {
        "status": "PASS", "methods": list(METHODS), "rna_contexts": list(CONTEXTS), "tumors": list(TUMORS),
        "anchors": {"A1": "1-Strong-Anchor (HEX1)", "A2": "2-Strong-Anchors (HEX1 + LDH_L)"},
        "source_row_counts": {"A1_candidate_weights_before_arm_filter": 12800, "A1_strong_anchor_baseline": 3200,
                              "A2_weights": 3200, "PL1_evaluation_registry": 400,
                              "PL2A_truth_selection": int(len(selection)), "A22_matched_truth_pairs": 836,
                              "A1_A2_post_holdout_support_rows": int(len(pair_support))},
        "output_row_counts": {"candidate_weights": len(weights), "conditional_strata": len(strata),
                              "global_ess": len(global_ess), "contrast_product_support": len(contrast),
                              "post_holdout_candidate_support": len(selection), "post_holdout_pair_support": len(pair_support)},
        "keys": {"candidate_weights": len(weights), "conditional_strata": len(strata), "truth_pairs": len(pairs)},
        "weight_checks": {"duplicate_keys": 0, "nonfinite": 0, "negative": 0,
                          "conditional_normalization_atol": TOL,
                          "max_abs_conditional_normalization_error": float(np.max(np.abs(weights.groupby(KEY[:-1]).conditional_weight.sum().to_numpy() - 1))),
                          "max_abs_global_normalization_error": float(np.max(np.abs(weights.groupby(["anchor", "tumor", "mouse_id"]).global_weight.sum().to_numpy() - 1)))},
        "ess_ranges": {"conditional_min": float(strata.conditional_ess.min()), "conditional_max": float(strata.conditional_ess.max()),
                       "global_min": float(global_ess.global_ess.min()), "global_max": float(global_ess.global_ess.max()),
                       "contrast_product_min": float(contrast.contrast_product_ess.min()), "contrast_product_max": float(contrast.contrast_product_ess.max()),
                       "post_holdout_pair_min": float(pair_support.post_holdout_ess.min()), "post_holdout_pair_max": float(pair_support.post_holdout_ess.max())},
        "concentrated_strata": {"conditional_ess_le_2": int((strata.conditional_ess <= 2).sum()),
                                "singleton_positive_support": int((strata.positive_candidate_count <= 1).sum()),
                                "max_weight_ge_0_5": int((strata.max_candidate_weight >= 0.5).sum())},
        "holdout_status_counts": {str(k): int(v) for k, v in exclusions.items()},
        "input_files": inputs,
        "a1_weight_provenance": {"sha256": fig4_a1_hash, "manifest": SOURCES["fig4_manifest"]},
        "production_vs_exploratory": "Inputs are frozen BIO0/A20/PL1/PL2A/A22 production outputs. No secondary-anchor exploratory outputs are used.",
    }
    (DATA / "validation_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    summary = build()
    print(json.dumps({"status": summary["status"], "output_row_counts": summary["output_row_counts"],
                      "ess_ranges": summary["ess_ranges"], "concentrated_strata": summary["concentrated_strata"]}, indent=2))


if __name__ == "__main__":
    main()
