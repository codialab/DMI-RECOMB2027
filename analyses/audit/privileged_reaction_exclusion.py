#!/usr/bin/env python3
"""Read-only frozen-row audit; never executes production code or changes inputs.

Requires explicit access to --source-root. Writes only aggregate audit evidence
inside this Git repository. Rights-pending source tables are not redistributed.
"""
import argparse
import hashlib
import io
import json
import lzma
import os
from pathlib import Path
import platform
import subprocess

import numpy as np
import pandas as pd
import pyarrow
import pyarrow.parquet as pq
import scipy
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
EX = "EX_glc__D_e"
PRIVILEGED = {"HEX1", EX, "LDH_L"}
EXPECTED_SOURCE_HEAD = "16379a169f45ec556274953eb66ae12293e1e985"
CACHE = "outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstructed/flux_cache.npz"
A23 = "outputs/dmi_bridge_a23_a1_a2_synthesis_v1/"


def git(root, *args):
    return subprocess.check_output(
        ["git", "-C", str(root), *args],
        env={**os.environ, "GIT_OPTIONAL_LOCKS": "0"}, text=True,
    ).strip()


def source_state(root):
    return {
        "head": git(root, "rev-parse", "HEAD"),
        "tracked_status": git(root, "status", "--short", "--untracked-files=no"),
        "tracked_diff": git(root, "diff", "--raw", "HEAD"),
        "index_listing_sha256": hashlib.sha256(
            git(root, "ls-files", "--stage").encode()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--replace-audit-outputs", action="store_true",
                        help="Replace only this script's two destination audit summaries")
    args = parser.parse_args()
    source = args.source_root.resolve(strict=True)
    assert Path(git(source, "rev-parse", "--show-toplevel")) == source
    assert Path(git(ROOT, "rev-parse", "--show-toplevel")) == ROOT
    assert source != ROOT
    before = source_state(source)
    assert before["head"] == EXPECTED_SOURCE_HEAD
    assert before["tracked_status"] == before["tracked_diff"] == ""
    evidence = {"audit_date": "2026-10-06", "source_state_before": before,
                "starting_destination_commit": git(ROOT, "rev-parse", "HEAD"),
                "scope": "Frozen-row filtering and aggregates only; no production execution",
                "versions": {"python": platform.python_version(), "numpy": np.__version__,
                             "pandas": pd.__version__, "scipy": scipy.__version__,
                             "pyarrow": pyarrow.__version__},
                "sources": [], "membership": []}
    expected = {x["source_path"]: x["sha256"] for x in json.loads(
        (ROOT / "manifests/FIGURE_INPUTS.json").read_text())["inputs"]}
    expected[CACHE] = "f9644da4ebd875cf5f0b34cc58b20faf0a0457d7c7a44635e8a010d8896c3dbb"
    expected[A23 + "BRIDGEA23_MANIFEST.json"] = "b8e57074f3df23fe7708145bc01f51a19f946545b1c16e61628c2678dbcabd18"
    seen = set()

    def path(rel):
        p = (source / rel).resolve(strict=True)
        p.relative_to(source)  # Reject any outside-repository symlink target.
        assert p.is_file()
        if rel not in seen:
            h = hashlib.sha256()
            with p.open("rb") as f:
                for block in iter(lambda: f.read(1024 * 1024), b""):
                    h.update(block)
            digest = h.hexdigest()
            if rel in expected:
                assert digest == expected[rel], f"Frozen identity mismatch: {rel}"
            evidence["sources"].append({"source_path": rel, "sha256": digest,
                                        "bytes": p.stat().st_size, "copied": False,
                                        "expected_hash_verified": rel in expected})
            seen.add(rel)
        return p

    def tsv(rel, **kwargs):
        return pd.read_csv(path(rel), sep="\t", float_precision="round_trip", **kwargs)

    def parquet(rel, columns=None):
        # In-memory decompression: no temporary source/data copies on disk.
        with lzma.open(path(rel), "rb") as f:
            buf = io.BytesIO(f.read())
        return pq.read_table(buf, columns=columns).to_pandas()

    def membership(rel, ids, exact=None):
        ids = pd.Series(ids, dtype="object")
        universe = set(ids)
        if exact is not None:
            assert universe == exact, rel
        assert EX not in universe, rel
        evidence["membership"].append({"source_path": rel, "rows": len(ids),
            "unique_reactions": len(universe),
            "privileged_row_counts": {r: int(ids.eq(r).sum()) for r in sorted(PRIVILEGED)}})
        return universe

    with np.load(path(CACHE), allow_pickle=False) as cache:
        rxns = cache["rxns"].astype(str)
    universe = membership(CACHE, rxns)
    assert len(rxns) == len(universe) == 4181
    a = universe - {"HEX1", "LDH_L"}
    b = universe - {EX, "LDH_L"}
    c = universe - PRIVILEGED
    assert a == c and len(a) == 4179 and len(b) == 4180
    evidence["policies"] = {"A_count": len(a), "B_count": len(b), "C_count": len(c),
        "A_equals_C": a == c, "B_restores": sorted(b - a),
        "EX_present_on_candidate_axis": EX in universe}

    stage4_rel = "dmi_sa_r1_cup_reuse_20260926/selected_model/ensembles/stage4_projection_4747/metadata.json"
    stage4 = json.loads(path(stage4_rel).read_text())["ordered_projection_reactions"]
    offline = "13_recomb_identifiability/exports/dmi_sa_offline_feasibility/"
    subsystems = set(pd.read_csv(path(offline + "validation/rna_effects_used.csv"),
                                 usecols=["subsystem"])["subsystem"].dropna())
    annotation = tsv(offline + "reaction_metadata/reaction_metadata.tsv",
                     usecols=["reaction_id", "subsystem"]).dropna().drop_duplicates()
    eligible = set(annotation.loc[annotation.subsystem.isin(subsystems), "reaction_id"])
    recovered = [r for r in stage4 if r in eligible]
    assert recovered == list(rxns)
    evidence["upstream_axis_selection"] = {
        "stage4_count": len(stage4), "stage4_contains_EX": EX in stage4,
        "rna_subsystem_count": len(subsystems),
        "EX_subsystems": sorted(set(annotation.loc[annotation.reaction_id.eq(EX), "subsystem"])),
        "EX_eligible": EX in eligible, "recovered_axis_exact_match": True}
    path("outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstruct_flux_cache.py")

    for rel, exact in [
        ("outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_REACTION_SUMMARY.tsv", universe - {"HEX1"}),
        ("outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_REACTION_REGISTRY.tsv", universe - {"HEX1"}),
        ("outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_REACTION_SPLIT.tsv", universe - {"HEX1"}),
        ("outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_REACTION_SUMMARY.tsv.xz", a),
        ("outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_REACTION_REGISTRY.tsv", a),
        ("outputs/dmi_bridge_pl2d_confirmation_v1/BRIDGEPL2D_CONFIRMATION_REACTION.tsv", None),
    ]:
        membership(rel, tsv(rel, usecols=["reaction_id"]).reaction_id, exact)

    # The admitted A23 manifests bind its summary records to the common registry.
    manifest_rel = A23 + "BRIDGEA23_MANIFEST.json"
    manifest = json.loads(path(manifest_rel).read_text())
    for name, digest in manifest["artifact_sha256"].items():
        expected[A23 + name] = digest
    assert manifest["new_production_performed"] is False
    evidence["A23_case_population"] = manifest["case_population"]
    results = []

    def row(section, anchor, scope, metric, unit, x, y, nx, ny, source_rel,
            method="recomputed_from_filtered_frozen_rows"):
        assert nx == ny and (x == y or (np.isnan(x) and np.isnan(y)))
        results.append(dict(section=section, anchor_setting=anchor, scope=scope,
            metric=metric, aggregation_unit=unit, current_value=float(x),
            policy_c_value=float(y), signed_change=float(y-x), absolute_change=float(abs(y-x)),
            current_n=int(nx), policy_c_n=int(ny), removed_rows=int(nx-ny),
            calculation_method=method, source_path=source_rel))

    def summaries(section, scope, frame, columns, unit, source_rel):
        filtered = frame.loc[frame.reaction_id.isin(c)]
        assert len(filtered) == len(frame)
        for arm, col, metric in columns:
            x = frame[col].to_numpy(float)
            y = filtered[col].to_numpy(float)
            x, y = x[np.isfinite(x)], y[np.isfinite(y)]
            assert np.array_equal(x, y)
            for stat, fn in [("mean", np.mean), ("median", np.median)]:
                row(section, arm, scope, metric+"_"+stat, unit,
                    fn(x), fn(y), len(x), len(y), source_rel)

    fig2_rel = "figures/fig2/figure2_DG_points.csv.xz"
    fig2 = pd.read_csv(path(fig2_rel), float_precision="round_trip")
    membership(fig2_rel, fig2.reaction_id)
    expected_rho = {"D": (2454, .8571160663367106), "E": (2485, .8339440056039364),
                    "F": (838, .8351435823183291), "G": (843, .8493657186460156)}
    for panel, frame in fig2.groupby("panel", sort=True):
        arm = "A1" if panel in ("D", "F") else "A2"
        filtered = frame.loc[frame.reaction_id.isin(c)]
        xcol = "x_direction_explained_magnitude_variance"
        ycol = "y_directional_usefulness_fraction"
        rho = float(spearmanr(frame[xcol], frame[ycol]).statistic)
        rho_c = float(spearmanr(filtered[xcol], filtered[ycol]).statistic)
        assert len(frame) == expected_rho[panel][0]
        assert abs(rho - expected_rho[panel][1]) <= 1e-6
        row("Figure2", arm, panel, "development_rho" if panel in ("D", "E") else "confirmation_rho",
            "finite_plotted_reaction", rho, rho_c, len(frame), len(filtered), fig2_rel)
        summaries("Figure2", panel, frame, [(arm, ycol, "directional_usefulness_fraction"),
            (arm, "y_mean_g_info", "information_advantage")],
            "equal_mean_of_finite_plotted_reaction_means", fig2_rel)
    evidence["fig2_baseline_parity_tolerance"] = 1e-6

    fig3_rel = "figures/fig3/tables/fig3_gain_reaction_evaluation.parquet.xz"
    fig3 = parquet(fig3_rel, ["anchor_setting", "reaction_id", "g_info_mean",
        "correct_gain_mean", "wrong_gain_mean", "directional_advantage"])
    membership(fig3_rel, fig3.reaction_id)
    for arm, frame in fig3.groupby("anchor_setting", sort=True):
        assert len(frame) == 660237
        public_arm = "A1" if arm == "A1" else "A2"
        summaries("Figure3", "matched_non_tie_evaluation_reaction", frame,
            [(public_arm, "correct_gain_mean", "correct_cue_gain"),
             (public_arm, "wrong_gain_mean", "wrong_cue_gain"),
             (public_arm, "directional_advantage", "correct_minus_wrong_separation"),
             (public_arm, "g_info_mean", "information_advantage")],
            "equal_weight_evaluation_reaction_after_within_evaluation_distinct_truth_mean", fig3_rel)
        filtered = frame.loc[frame.reaction_id.isin(c)]
        for col in ["correct_gain_mean", "wrong_gain_mean", "directional_advantage"]:
            x, y = frame[col].to_numpy(float), filtered[col].to_numpy(float)
            for q in [.05, .25, .75, .95]:
                row("Figure3", public_arm, "matched_non_tie_evaluation_reaction", col+f"_q{int(q*100)}",
                    "evaluation_reaction", np.quantile(x,q), np.quantile(y,q),len(x),len(y),fig3_rel)
            # This is an ECDF diagnostic, NOT the truth-level usefulness predicate.
            row("Figure3", public_arm, "matched_non_tie_evaluation_reaction", col+"_strict_positive_fraction",
                "evaluation_reaction", np.mean(x>0), np.mean(y>0),len(x),len(y),fig3_rel)
    del fig3
    cue_rel = "figures/fig3/tables/fig3_gain_cue_direction.parquet.xz"
    membership(cue_rel, parquet(cue_rel, ["reaction_id"]).reaction_id)

    metrics = ["H_dir", "dominant_direction_mass", "delta_v_B_width80",
               "direction_explained_magnitude_variance", "supported_sign_state_count", "non_tie_coverage"]
    frames = []
    for method in ["CORDA", "GIMME", "RIPTiDe", "iMAT"]:
        rel = f"figures/fig4/data/fig4_geometry_paired/algorithm={method}/part.parquet.xz"
        cols = ["evaluation_id", "reaction_id"] + [f"{arm}_{m}" for arm in ("A1","A2") for m in metrics+["ESS"]]
        frame = parquet(rel, cols)
        membership(rel, frame.reaction_id, a)
        assert len(frame) == 417900
        frames.append(frame)
    fig4 = pd.concat(frames, ignore_index=True)
    assert len(fig4) == 1671600
    geometry_source = "figures/fig4/data/fig4_geometry_paired/algorithm=*/part.parquet.xz"
    for metric in metrics:
        paired = fig4.loc[np.isfinite(fig4["A1_"+metric]) & np.isfinite(fig4["A2_"+metric])].copy()
        paired["paired_delta"] = paired["A2_"+metric] - paired["A1_"+metric]
        summaries("Figure4", "paired_finite_all_evaluations", paired,
            [("A1", "A1_"+metric, metric), ("A2", "A2_"+metric, metric),
             ("A2_minus_A1", "paired_delta", metric)], "evaluation_reaction", geometry_source)
    # ESS is a joint-ensemble property: do not count reactions as independent ESS observations.
    for arm in ["A1", "A2"]:
        ess = fig4.groupby("evaluation_id")[arm+"_ESS"]
        assert ess.nunique().max() == 1
        values = ess.first().to_numpy(float)
        assert len(values) == 400
        row("Figure4", arm, "all_evaluations", "joint_ESS_mean", "evaluation",
            np.mean(values), np.mean(values), 400, 400, geometry_source,
            "recomputed_from_identical_filtered_evaluation_sets")
    del fig4, frames
    for filename in ["fig4_panelC_candidate_groups.parquet.xz", "fig4_panelC_candidate_metrics.parquet.xz"]:
        rel = "figures/fig4/data/"+filename
        frame = parquet(rel, ["reaction_id"])
        membership(rel, frame.reaction_id)
        assert frame.reaction_id.isin(c).all()
    evidence["figure_distribution_arrays_unchanged"] = True
    evidence["figure4_example_candidates_unchanged"] = True

    # Every A23 summary input has exactly the same admitted rows under A and C.
    # Record its existing statistics without pretending to recompute absent
    # reaction-level fields from already aggregated numbers.
    for filename, section, fields in [
        ("BRIDGEA23_GEOMETRY_COMPARISON.tsv", "A23_geometry", ["a1_mean","a2_mean","paired_delta_mean","paired_delta_median"]),
        ("BRIDGEA23_UTILITY_COMPARISON.tsv", "A23_case_pooled_utility", ["a1_correct_mean_gain","a2_correct_mean_gain","a1_specificity_mean","a2_specificity_mean","paired_correct_delta","paired_specificity_delta"]),
        ("BRIDGEA23_CONTROL_SUMMARY.tsv", "A23_case_pooled_control", ["a1_mean_gain","a2_mean_gain","paired_delta_mean"]),
        ("BRIDGEA23_GEOMETRY_UTILITY_PERSISTENCE.tsv", "A23_association", ["a1_rho","a2_rho","rho_a2_minus_a1"]),
    ]:
        rel = A23+filename
        frame = tsv(rel)
        for _, record in frame.iterrows():
            scope = ":".join(str(record[k]) for k in ["scope","level","split","context","feature","arm"] if k in record)
            for col in fields:
                arm = "A1" if col.startswith("a1_") else "A2" if col.startswith("a2_") else "A2_minus_A1"
                count = next((int(record[k]) for k in ["finite_paired_rows","non_tie_pairs","denominator", "a1_finite" if arm=="A1" else "a2_finite"] if k in record), 0)
                row(section, arm, scope, col, str(record.get("unit", "source_defined_finite_population")),
                    record[col], record[col], count,count,rel,"frozen_statistic_preserved_by_proven_identical_input_set")

    m1_rel = "outputs/dmi_bridge_m1_manuscript_evidence_v1/BRIDGEM1_DEVELOPMENT_CONFIRMATION.tsv"
    m1 = tsv(m1_rel)
    evidence["original_A1_M1_development_confirmation_records"] = m1.astype(object).where(m1.notna(), None).to_dict(orient="records")
    evidence["original_A1_population_note"] = "PL2 population 4180; no EX row; original development result retains LDH_L, unlike matched final Figure2"
    for _, record in m1.iterrows():
        if record["cohort"] != "CONFIRMATION":
            # The older 4180-reaction development result includes LDH_L.
            # It is not the frozen 4179 baseline requested for this comparison.
            continue
        for col in ["eta2_usefulness_rho", "positive_contexts", "context_count",
                    "minimum_finite_context_count", "bootstrap_median", "bootstrap_interval_lower", "bootstrap_interval_upper"]:
            if pd.notna(record[col]):
                row("M1_original_A1", "A1", record["cohort"], col,
                    "original_source_defined_reaction_population", record[col], record[col],
                    record["finite_reactions"],record["finite_reactions"],m1_rel,
                    "frozen_statistic_preserved_by_proven_identical_input_set_no_bootstrap_rerun")

    after = source_state(source)
    assert after == before
    evidence["source_state_after"] = after
    output = ROOT / "tables/audit/privileged_reaction_exclusion_sensitivity.tsv"
    evidence_path = ROOT / "manifests/provenance/PRIVILEGED_REACTION_EXCLUSION_EVIDENCE.json"
    table_bytes = pd.DataFrame(results).to_csv(sep="\t", index=False, float_format="%.17g", na_rep="NA").encode()
    evidence["sensitivity_table_sha256"] = hashlib.sha256(table_bytes).hexdigest()
    evidence["sensitivity_rows"] = len(results)
    evidence["all_absolute_changes_zero"] = all(r["absolute_change"] == 0 for r in results)
    assert evidence["all_absolute_changes_zero"]
    evidence_bytes = (json.dumps(evidence, indent=2, allow_nan=False)+"\n").encode()
    # An interrupted audit may have written the identical table. Never replace
    # existing outputs with different bytes, and never touch frozen inputs.
    for p, data in [(output, table_bytes), (evidence_path, evidence_bytes)]:
        if p.exists() and not args.replace_audit_outputs:
            assert p.read_bytes() == data, f"Existing audit output differs: {p}; use --replace-audit-outputs"
        else:
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(data)
    print(json.dumps({"policies": evidence["policies"], "sensitivity_rows": len(results),
                      "all_absolute_changes_zero": True, "source_unchanged": before==after}))


if __name__ == "__main__":
    main()
