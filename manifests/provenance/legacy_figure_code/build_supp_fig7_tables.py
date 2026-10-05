#!/usr/bin/env python3
"""Build deterministic S7 tables from frozen Figure 4C artifacts."""
from __future__ import annotations

import hashlib
import json
import lzma
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DATA = OUT / "data"
FIG4_DATA = ROOT / "figures/fig4/data"
FIG4_MANIFEST_PATH = FIG4_DATA / "fig4_data_manifest.json"
GROUPS_PATH = FIG4_DATA / "fig4_panelC_candidate_groups.parquet.xz"
METRICS_PATH = FIG4_DATA / "fig4_panelC_candidate_metrics.parquet.xz"
SUPPORT_PATH = FIG4_DATA / "fig4_panelC_candidate_support.parquet.xz"
NOTEBOOK_PATH = ROOT / "figures/fig4/fig4_panelB_distributions_updated.ipynb"
FIG4_BUILDER_PATH = ROOT / "figures/fig4/build_fig4_tables.py"
BIO0_WEIGHTS_PATH = ROOT / "outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_FOUR_ARM_WEIGHTS.tsv.xz"
BIO0_MANIFEST_PATH = BIO0_WEIGHTS_PATH.parent / "BRIDGEBIO0_MANIFEST.json"
A20_WEIGHTS_PATH = ROOT / "outputs/dmi_bridge_a20_dual_anchor_qualification_v1/BRIDGEA20_DUAL_ANCHOR_WEIGHTS.tsv.xz"
A20_MANIFEST_PATH = A20_WEIGHTS_PATH.parent / "BRIDGEA20_MANIFEST.json"
FLUX_CACHE_PATH = ROOT / "outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstructed/flux_cache.npz"
CANDIDATE_TABLE_PATH = ROOT / "12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/candidate_fraction_table.tsv.xz"

METHODS = ("CORDA", "GIMME", "iMAT", "RIPTiDe")
ARMS = {"A1": "1-Strong-Anchor (HEX1)", "A2-L": "2-Strong-Anchors (HEX1 + LDH_L)"}
SHIFT_METRIC = "median_absolute_A2_minus_A1_H_dir_shift"
TIE_TOL = 1e-12
WEIGHT_TOL = 1e-10


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_path(path: Path) -> Path:
    require(path.is_file() and not path.is_symlink(), f"missing or non-regular frozen source: {path.relative_to(ROOT)}")
    resolved = path.resolve(strict=True)
    require(resolved.is_relative_to(ROOT), f"source escapes repository root: {path.relative_to(ROOT)}")
    return resolved


def read_parquet_xz(path: Path) -> pd.DataFrame:
    with tempfile.TemporaryDirectory(prefix="supp-s7-parquet-") as tempdir:
        parquet_path = Path(tempdir) / "source.parquet"
        with lzma.open(path, "rb", format=lzma.FORMAT_XZ) as source, parquet_path.open("wb") as destination:
            shutil.copyfileobj(source, destination, length=1 << 20)
        return pd.read_parquet(parquet_path)


def verify_input_provenance() -> tuple[dict, dict[str, str]]:
    fig4_manifest = json.loads(source_path(FIG4_MANIFEST_PATH).read_text())
    require(fig4_manifest.get("status") == "PASS", "Figure 4 data manifest is not PASS")
    generated = fig4_manifest.get("generated_sha256", {})
    listed_generated = {
        GROUPS_PATH.relative_to(ROOT).as_posix(): GROUPS_PATH,
        METRICS_PATH.relative_to(ROOT).as_posix(): METRICS_PATH,
        SUPPORT_PATH.relative_to(ROOT).as_posix(): SUPPORT_PATH,
    }
    input_hashes: dict[str, str] = {}
    for relative, path in listed_generated.items():
        expected = generated.get(relative)
        require(expected is not None, f"Figure 4 manifest has no generated checksum for {relative}")
        actual = sha256(source_path(path))
        require(actual == expected, f"Figure 4 generated-artifact checksum mismatch: {relative}")
        input_hashes[relative] = actual

    source_used = fig4_manifest.get("source_files_used", {})
    for path in (BIO0_WEIGHTS_PATH, A20_WEIGHTS_PATH, FLUX_CACHE_PATH, CANDIDATE_TABLE_PATH):
        relative = path.relative_to(ROOT).as_posix()
        expected = source_used.get(relative)
        require(expected is not None, f"Figure 4 manifest has no source checksum for {relative}")
        actual = sha256(source_path(path))
        require(actual == expected, f"Figure 4 source checksum mismatch: {relative}")
        input_hashes[relative] = actual

    a20_manifest = json.loads(source_path(A20_MANIFEST_PATH).read_text())
    a20_expected = a20_manifest.get("artifact_sha256", {}).get(A20_WEIGHTS_PATH.name)
    require(a20_expected == input_hashes[A20_WEIGHTS_PATH.relative_to(ROOT).as_posix()],
            "A20 manifest checksum mismatch for frozen A2 weights")
    input_hashes[FIG4_MANIFEST_PATH.relative_to(ROOT).as_posix()] = sha256(source_path(FIG4_MANIFEST_PATH))
    input_hashes[BIO0_MANIFEST_PATH.relative_to(ROOT).as_posix()] = sha256(source_path(BIO0_MANIFEST_PATH))
    input_hashes[A20_MANIFEST_PATH.relative_to(ROOT).as_posix()] = sha256(source_path(A20_MANIFEST_PATH))
    input_hashes[FIG4_BUILDER_PATH.relative_to(ROOT).as_posix()] = sha256(source_path(FIG4_BUILDER_PATH))
    input_hashes[NOTEBOOK_PATH.relative_to(ROOT).as_posix()] = sha256(source_path(NOTEBOOK_PATH))
    return fig4_manifest, input_hashes


def select_examples(groups: pd.DataFrame) -> pd.DataFrame:
    required = {
        "candidate_group_id", "reaction_id", "rna_context_key", "ct2a_mouse", "gl261_mouse",
        "truth_selection", "geometry_rank", SHIFT_METRIC,
    }
    require(required.issubset(groups.columns), f"candidate group table lacks columns: {sorted(required - set(groups.columns))}")
    require(len(groups) == 100 and not groups.candidate_group_id.duplicated().any(),
            "expected 100 unique frozen Figure 4C candidate groups")
    main = groups.loc[groups.candidate_group_id.eq("FCG_ed8b1873c8dbc03fe5eb")].copy()
    require(len(main) == 1 and main.iloc[0].reaction_id == "FACOAL204", "frozen Main Figure 4C FACOAL204 group is absent or ambiguous")

    # Fix selection from precomputed geometry summaries before reading support rows or plotting.
    pool = groups.loc[~groups.candidate_group_id.eq(main.iloc[0].candidate_group_id)].copy()
    pool[SHIFT_METRIC] = pd.to_numeric(pool[SHIFT_METRIC], errors="coerce")
    require(np.isfinite(pool[SHIFT_METRIC].to_numpy(float)).all(), "nonfinite frozen geometry shift in candidate selection pool")
    pool = pool.sort_values([SHIFT_METRIC, "geometry_rank", "candidate_group_id"], ascending=[False, True, True], kind="mergesort")
    strong = pool.iloc[[0]].copy()
    remaining = pool.loc[~pool.candidate_group_id.isin(strong.candidate_group_id)]
    weak = remaining.sort_values([SHIFT_METRIC, "geometry_rank", "candidate_group_id"], ascending=[True, True, True], kind="mergesort").iloc[[0]].copy()

    main["example_order"] = 1
    main["example_label"] = "Main Figure 4C example"
    main["selection_basis"] = "Existing Main Figure 4C example; the repository labels it outcome-blind, but does not record a frozen rule specific to its selection."
    strong["example_order"] = 2
    strong["example_label"] = "Largest frozen A1-to-A2 geometry change"
    strong["selection_basis"] = f"Descriptive post hoc selection: largest {SHIFT_METRIC} among the other 99 frozen candidates; ties resolved by geometry_rank then candidate_group_id."
    weak["example_order"] = 3
    weak["example_label"] = "No median directional-entropy shift"
    weak["selection_basis"] = f"Descriptive post hoc selection: smallest {SHIFT_METRIC} among the remaining frozen candidates; ties resolved by geometry_rank then candidate_group_id. A zero entropy shift does not imply the full flux-magnitude distributions are unchanged."
    selected = pd.concat([main, strong, weak], ignore_index=True).sort_values("example_order", kind="stable").reset_index(drop=True)
    require(selected.candidate_group_id.nunique() == 3, "S7 examples must be distinct candidate groups")
    return selected


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    fig4_manifest, input_hashes = verify_input_provenance()
    groups = read_parquet_xz(source_path(GROUPS_PATH))
    selected = select_examples(groups)

    metrics_all = read_parquet_xz(source_path(METRICS_PATH))
    metrics = metrics_all.loc[metrics_all.candidate_group_id.isin(selected.candidate_group_id)].copy()
    require(len(metrics) == 12 and set(metrics.algorithm) == set(METHODS),
            "each selected example must have exactly four frozen method metrics")
    require(not metrics.duplicated(["candidate_group_id", "algorithm"]).any(), "duplicate selected example/method metric rows")

    selection_columns = [
        "example_order", "example_label", "selection_basis", "candidate_group_id", "geometry_rank",
        "reaction_id", "rna_context_key", "ct2a_mouse", "gl261_mouse", "truth_selection", SHIFT_METRIC,
    ]
    selection = selected[selection_columns].merge(
        metrics[["candidate_group_id", "algorithm", "evaluation_id", "truth_pair_id", "subsystem"]],
        on="candidate_group_id", how="left", validate="one_to_many",
    ).rename(columns={"algorithm": "method"})
    require(len(selection) == 12 and selection.evaluation_id.notna().all(), "selected example method/evaluation provenance is incomplete")
    selection = selection.sort_values(["example_order", "method"], kind="stable").reset_index(drop=True)

    # Distribution support is loaded only after the deterministic geometry-metric selection is fixed.
    support_all = read_parquet_xz(source_path(SUPPORT_PATH))
    support = support_all.loc[support_all.candidate_group_id.isin(selected.candidate_group_id)].copy()
    expected_rows = 3 * len(METHODS) * 2 * 400
    require(len(support) == expected_rows, f"expected {expected_rows} frozen selected support rows; found {len(support)}")
    require(set(support.arm) == set(ARMS) and set(support.method) == set(METHODS), "unexpected method or anchor arm in selected support")
    require(not support.duplicated(["candidate_group_id", "method", "arm", "support_id"]).any(), "duplicate weighted support state")
    require(support.groupby(["candidate_group_id", "method", "arm"], observed=True).size().eq(400).all(),
            "each selected example/method/arm must contain 400 support states")
    require(support.groupby(["candidate_group_id", "method", "arm"], observed=True).support_id.nunique().eq(400).all(),
            "support IDs are not unique within each selected example/method/arm")
    support_ids = support.groupby(["candidate_group_id", "method", "arm"], observed=True).support_id.agg(
        lambda values: frozenset(values.astype(int)))
    for gid, method in support[["candidate_group_id", "method"]].drop_duplicates().itertuples(index=False, name=None):
        require(support_ids.loc[(gid, method, "A1")] == support_ids.loc[(gid, method, "A2-L")],
                f"A1/A2 Cartesian support IDs are not paired for {gid}/{method}")
    require(np.isfinite(support[["delta_v_B", "weight"]].to_numpy(float)).all(), "nonfinite weighted distribution values")
    require((support.weight >= 0).all(), "negative candidate-state weight")
    weight_sums = support.groupby(["candidate_group_id", "method", "arm"], observed=True).weight.sum()
    require(np.allclose(weight_sums.to_numpy(float), 1.0, rtol=0, atol=WEIGHT_TOL), "frozen product weights are not normalized")
    expected_sign = np.where(support.delta_v_B.to_numpy(float) > TIE_TOL, 1,
                             np.where(support.delta_v_B.to_numpy(float) < -TIE_TOL, -1, 0))
    require(np.array_equal(support.sign.to_numpy(int), expected_sign), "frozen sign labels disagree with delta_v_B and Figure 4C tie tolerance")

    # Check that support identity matches the frozen metric/evaluation records for both arms.
    identity = support[["candidate_group_id", "method", "evaluation_id", "reaction_id", "rna_context_key"]].drop_duplicates()
    require(len(identity) == 12 and not identity.duplicated(["candidate_group_id", "method"]).any(),
            "support identity grain is not one row per example and method")
    compare = identity.merge(
        metrics[["candidate_group_id", "algorithm", "evaluation_id", "reaction_id", "rna_context_key"]],
        left_on=["candidate_group_id", "method", "evaluation_id", "reaction_id", "rna_context_key"],
        right_on=["candidate_group_id", "algorithm", "evaluation_id", "reaction_id", "rna_context_key"],
        how="left", validate="one_to_one",
    )
    require(compare.algorithm.notna().all(), "support evaluation or reaction identity does not match frozen candidate metrics")

    group_lookup = selected.set_index("candidate_group_id")
    method_meta = selection[["candidate_group_id", "method", "evaluation_id", "example_order", "example_label",
                             "reaction_id", "rna_context_key", "ct2a_mouse", "gl261_mouse"]]
    support = support.merge(method_meta, on=["candidate_group_id", "method", "evaluation_id", "reaction_id", "rna_context_key"],
                            how="left", validate="many_to_one")
    require(support.example_order.notna().all(), "support rows did not map to selected example metadata")
    support["example_order"] = support.candidate_group_id.map(group_lookup.example_order)
    support["example_label"] = support.candidate_group_id.map(group_lookup.example_label)
    support["anchor"] = support.arm.map({"A1": "A1", "A2-L": "A2"})
    support["anchor_definition"] = support.arm.map(ARMS)
    support = support.sort_values(["example_order", "method", "arm", "support_id"], kind="stable").reset_index(drop=True)
    plot_path = DATA / "supp_fig7_weighted_distributions.tsv.gz"
    support[["example_order", "example_label", "candidate_group_id", "reaction_id", "method", "evaluation_id",
             "rna_context_key", "ct2a_mouse", "gl261_mouse", "truth_selection", "arm", "anchor", "anchor_definition",
             "support_id", "delta_v_B", "weight", "sign"]].to_csv(
        plot_path, sep="\t", index=False, compression="gzip", float_format="%.17g")
    selection_path = DATA / "supp_fig7_example_selection.tsv"
    selection.to_csv(selection_path, sep="\t", index=False, float_format="%.17g")

    mass_rows = []
    for keys, frame in support.groupby(["candidate_group_id", "method", "arm"], sort=True, observed=True):
        gid, method, arm = keys
        mass_rows.append({"candidate_group_id": gid, "method": method, "arm": arm,
                          "support_states": int(len(frame)), "weight_sum": float(frame.weight.sum()),
                          "weighted_delta_min": float(frame.delta_v_B.min()), "weighted_delta_max": float(frame.delta_v_B.max()),
                          "effective_support": float(1.0 / np.square(frame.weight.to_numpy(float)).sum())})
    support_summary = pd.DataFrame(mass_rows)
    summary_path = DATA / "supp_fig7_support_summary.tsv"
    support_summary.to_csv(summary_path, sep="\t", index=False, float_format="%.17g")

    output = {
        "status": "PASS",
        "scientific_definition": {
            "A1": "1-Strong-Anchor (HEX1)",
            "A2": "2-Strong-Anchors (HEX1 + LDH_L)",
            "delta_v_B": "CT2A candidate flux minus GL261 candidate flux, as constructed by the authoritative Figure 4C mats_delta implementation",
            "flux_units": "mmol gDW^-1 h^-1",
            "production_distribution": "Figure 4C frozen 400-state Cartesian support with product weights; no independent weighting redesign",
        },
        "selection": {
            "pool": "100 unique groups from the frozen Figure 4C candidate_groups table (the frozen manifest reports 100 eligible groups before top-100 export)",
            "main_example": "FACOAL204 is retained because it is the existing Main Figure 4C example; no formal frozen criterion specific to its selection was found",
            "additional_examples": f"After excluding FACOAL204, select the largest {SHIFT_METRIC}; then from the remaining groups select the smallest {SHIFT_METRIC}. Ties use ascending geometry_rank then candidate_group_id. Criteria use frozen geometry summaries and are applied before reading support distributions or plotting.",
            "representativeness": "Descriptive post hoc examples only; no statistical representativeness claim",
            "examples": selection[selection.example_order.isin([1, 2, 3])][["example_order", "example_label", "candidate_group_id", "geometry_rank", "reaction_id", "rna_context_key", "ct2a_mouse", "gl261_mouse", SHIFT_METRIC, "selection_basis"]].drop_duplicates("candidate_group_id").to_dict(orient="records"),
        },
        "dimensions": {
            "candidate_groups_in_pool": int(len(groups)), "candidate_groups_displayed": 3,
            "selected_method_evaluations": int(len(selection)), "methods": list(METHODS),
            "support_rows": int(len(support)), "rows_per_example_method_arm": 400,
            "example_method_arm_distributions": int(support.groupby(["candidate_group_id", "method", "arm"], observed=True).ngroups),
            "unique_support_states": int(support[["candidate_group_id", "method", "arm", "support_id"]].drop_duplicates().shape[0]),
        },
        "integrity": {
            "duplicate_support_keys": 0, "nonfinite_delta_or_weight": 0, "negative_weights": 0,
            "weight_sum_atol": WEIGHT_TOL, "maximum_absolute_weight_sum_error": float(np.max(np.abs(weight_sums.to_numpy(float) - 1.0))),
            "sign_labels_match_delta_v_B": True,
            "A1_A2_cartesian_support_ids_matched": True,
            "all_method_evaluations_match_frozen_metrics": True,
            "mass_summary_rows": int(len(support_summary)),
            "effective_support_range": [float(support_summary.effective_support.min()), float(support_summary.effective_support.max())],
        },
        "exclusions": {"candidate groups": {"count": int(len(groups) - 3), "reason": "Only the retained Main Figure 4C example and two groups selected by the declared deterministic frozen-summary rule are displayed."}},
        "source_hashes": input_hashes,
        "figure4_source_manifest_sha256": sha256(FIG4_MANIFEST_PATH),
        "figure4_candidate_ranking_rule": fig4_manifest.get("candidate_ranking", {}).get("rule"),
        "output_hashes": {},
    }
    output["output_hashes"] = {
        path.relative_to(ROOT).as_posix(): sha256(path)
        for path in (plot_path, selection_path, summary_path)
    }
    manifest_path = DATA / "supp_fig7_manifest.json"
    manifest_path.write_text(json.dumps(output, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"status": output["status"], "selected_examples": output["selection"]["examples"],
                      "dimensions": output["dimensions"], "integrity": output["integrity"]}, indent=2))


if __name__ == "__main__":
    main()
