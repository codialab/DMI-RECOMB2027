#!/usr/bin/env python3
"""Build deterministic S5 plotting tables from frozen PL2/PL2B/A22 outputs."""
from __future__ import annotations

import hashlib
import io
import json
import lzma
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures/supp_fig5"
DATA = OUT / "data"
FIG3_TABLES = ROOT / "figures/fig3/tables"
SOURCES = {
    "A1": {
        "directory": "outputs/dmi_bridge_pl2b_sign_only_utility_v1",
        "prefix": "BRIDGEPL2B",
        "manifest": "BRIDGEPL2B_MANIFEST.json",
        "expected_manifest_sha256": "7d9ca97a3f4ea5aa8b0b983e09ee83fc3c82f326d1b1cbf9c6137cfadf25539d",
        "expected_status": "PL2B_SIGN_ONLY_UTILITY_COMPLETE",
    },
    "A2": {
        "directory": "outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1",
        "prefix": "BRIDGEA22",
        "manifest": "BRIDGEA22_MANIFEST.json",
        "expected_manifest_sha256": "15b237151005538e5077a352868f873ef00466f4f5ea49307b781a8a8cd0e7cb",
        "expected_status": "BRIDGE_A22_SIGN_ONLY_UTILITY_FROZEN",
    },
}
FIG3_HASHES = {
    "fig3_gain_reaction_evaluation.parquet.xz": "7a83a31d29317f7b0d906fff73691428d98427d93d2fdebd4ef91294f8a46676",
    "fig3_gain_cue_direction.parquet.xz": "3c571ad4b97bb451d165c00fd464ccb51b0b5de8c8e69a49d3b423d78eeee30d",
    "fig3_gain_summary.parquet.xz": "0529bab84c6ad4cfc4a3c9e0d07ad418afcab3be52d2254183400dcd679af482",
    "fig3_gain_definition_audit.parquet.xz": "32370fc711e923d1476464e7e16e62b8f1f44cee7e1ae0de967cf695de0f4533",
    "fig3_directional_advantage_summary.parquet.xz": "d2ec6a86e344877f878b63e9ad03c5de0823f821b7527b8e4191f1ab5d39afab",
}
FIG3_METADATA_SHA256 = "d6ebfaa08aec0255104463bd327c55a74f24a1d074589d563782cb5953ab14a7"
A23_MANIFEST = "outputs/dmi_bridge_a23_a1_a2_synthesis_v1/BRIDGEA23_MANIFEST.json"
A23_MANIFEST_SHA256 = "b8e57074f3df23fe7708145bc01f51a19f946545b1c16e61628c2678dbcabd18"
EVALUATION_REGISTRY = "outputs/dmi_bridge_pl1_predictability_landscape_v1/BRIDGEPL1_EVALUATION_REGISTRY.tsv"
EVALUATION_REGISTRY_SHA256 = "aa468dcec11dea764682bd3d7eb14c03c003063e74d92080d0bfaabb366290a3"
LAMBDA = 0.25
Q_ENDPOINTS = {0.0: "wrong cue (q=0)", 1.0: "correct cue (q=1)"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def repo_file(relative: str) -> Path:
    path = ROOT / relative
    resolved = path.resolve(strict=True)
    if path.is_symlink() or not resolved.is_relative_to(ROOT) or not resolved.is_file():
        raise RuntimeError(f"Unsafe or missing repository input: {relative}")
    return resolved


def read_parquet_xz(path: Path) -> pd.DataFrame:
    with lzma.open(path, "rb", format=lzma.FORMAT_XZ) as source:
        return pd.read_parquet(io.BytesIO(source.read()), engine="pyarrow")


def verify_sources() -> tuple[dict[str, dict], dict[str, str]]:
    manifests: dict[str, dict] = {}
    hashes: dict[str, str] = {}
    for anchor, spec in SOURCES.items():
        manifest_rel = f"{spec['directory']}/{spec['manifest']}"
        manifest_path = repo_file(manifest_rel)
        digest = sha256(manifest_path)
        if digest != spec["expected_manifest_sha256"]:
            raise RuntimeError(f"Frozen {anchor} manifest checksum mismatch")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("status") != spec["expected_status"]:
            raise RuntimeError(f"Unexpected frozen {anchor} status")
        if manifest.get("primary_lambda") != LAMBDA:
            raise RuntimeError(f"Unexpected primary lambda in {anchor} manifest")
        context_name = f"{spec['prefix']}_CONTEXT_SUMMARY.tsv"
        context_rel = f"{spec['directory']}/{context_name}"
        context_path = repo_file(context_rel)
        expected_context = manifest.get("artifact_sha256", {}).get(context_name)
        if expected_context != sha256(context_path):
            raise RuntimeError(f"Frozen {anchor} context summary checksum mismatch")
        manifests[anchor] = manifest
        hashes[manifest_rel] = digest
        hashes[context_rel] = expected_context

    metadata_path = repo_file("figures/fig3/tables/fig3_gain_distribution_metadata.json")
    if sha256(metadata_path) != FIG3_METADATA_SHA256:
        raise RuntimeError("Figure 3 metadata checksum mismatch")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    hashes["figures/fig3/tables/fig3_gain_distribution_metadata.json"] = sha256(metadata_path)
    expected_sources = {
        "outputs/dmi_bridge_pl2b_sign_only_utility_v1/BRIDGEPL2B_MANIFEST.json": SOURCES["A1"]["expected_manifest_sha256"],
        "outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_MANIFEST.json": SOURCES["A2"]["expected_manifest_sha256"],
        A23_MANIFEST: A23_MANIFEST_SHA256,
    }
    if metadata.get("source_manifest_sha256") != expected_sources:
        raise RuntimeError("Figure 3 table source-manifest provenance mismatch")
    for name, expected in FIG3_HASHES.items():
        path = repo_file(f"figures/fig3/tables/{name}")
        actual = sha256(path)
        if actual != expected or metadata.get("output_sha256", {}).get(name) != expected:
            raise RuntimeError(f"Figure 3 table checksum mismatch: {name}")
        hashes[f"figures/fig3/tables/{name}"] = actual

    a23 = repo_file(A23_MANIFEST)
    if sha256(a23) != A23_MANIFEST_SHA256:
        raise RuntimeError("Frozen A1/A2 matched-population manifest checksum mismatch")
    hashes[A23_MANIFEST] = sha256(a23)
    registry = repo_file(EVALUATION_REGISTRY)
    if sha256(registry) != EVALUATION_REGISTRY_SHA256:
        raise RuntimeError("Frozen evaluation registry checksum mismatch")
    hashes[EVALUATION_REGISTRY] = sha256(registry)
    pl2_fingerprint = manifests["A1"]["fingerprint"]["admission"]["pl1_consumed_artifacts"]
    if pl2_fingerprint.get("evaluation_registry_sha256") != EVALUATION_REGISTRY_SHA256:
        raise RuntimeError("Evaluation registry does not match the PL2B production fingerprint")
    return manifests, hashes


def build_context_gain_table() -> tuple[pd.DataFrame, dict]:
    cue_path = repo_file("figures/fig3/tables/fig3_gain_cue_direction.parquet.xz")
    cue = read_parquet_xz(cue_path)
    wide_path = repo_file("figures/fig3/tables/fig3_gain_reaction_evaluation.parquet.xz")
    wide = read_parquet_xz(wide_path)
    wide_columns = {
        "anchor_setting", "evaluation_id", "reaction_id", "g_info_mean",
        "correct_gain_mean", "wrong_gain_mean", "directional_advantage",
    }
    if not wide_columns.issubset(wide.columns):
        raise RuntimeError(f"Figure 3 reaction-evaluation table lacks columns: {sorted(wide_columns - set(wide.columns))}")
    if wide.duplicated(["anchor_setting", "evaluation_id", "reaction_id"]).any():
        raise RuntimeError("Duplicate Figure 3 reaction-evaluation keys")
    wide_values = wide[["g_info_mean", "correct_gain_mean", "wrong_gain_mean", "directional_advantage"]].to_numpy(dtype=float)
    if not np.isfinite(wide_values).all():
        raise RuntimeError("Nonfinite values in frozen Figure 3 reaction-evaluation table")
    g_info_error = float(np.max(np.abs(wide["g_info_mean"] - 0.5 * (wide["correct_gain_mean"] - wide["wrong_gain_mean"]))))
    direction_error = float(np.max(np.abs(wide["directional_advantage"] - (wide["correct_gain_mean"] - wide["wrong_gain_mean"]))))
    if g_info_error > 5e-12 or direction_error > 1e-12:
        raise RuntimeError("Frozen Figure 3 gain identities failed the S5 audit")
    needed = {"anchor_setting", "evaluation_id", "reaction_id", "cue_direction", "gain_mean"}
    if not needed.issubset(cue.columns):
        raise RuntimeError(f"Figure 3 cue table lacks columns: {sorted(needed - set(cue.columns))}")
    cue["gain_mean"] = pd.to_numeric(cue["gain_mean"], errors="coerce")
    if not np.isfinite(cue["gain_mean"].to_numpy(dtype=float)).all():
        raise RuntimeError("Nonfinite cue gain in the frozen Figure 3 table")
    if cue.duplicated(["anchor_setting", "evaluation_id", "reaction_id", "cue_direction"]).any():
        raise RuntimeError("Duplicate Figure 3 cue keys")
    if set(cue["anchor_setting"].unique()) != {"A1", "A2-L"}:
        raise RuntimeError("Unexpected anchor arms in Figure 3 cue table")
    if set(cue["cue_direction"].unique()) != {"correct", "wrong"}:
        raise RuntimeError("Unexpected cue directions in Figure 3 cue table")

    registry = pd.read_csv(repo_file(EVALUATION_REGISTRY), sep="\t", dtype=str)
    if registry["evaluation_id"].duplicated().any():
        raise RuntimeError("Duplicate evaluation_id in the frozen evaluation registry")
    cue = cue.merge(
        registry[["evaluation_id", "algorithm", "rna_context_key"]],
        on="evaluation_id", how="left", validate="many_to_one", indicator=True,
    )
    if not cue["_merge"].eq("both").all():
        raise RuntimeError("Unmapped evaluation_id in Figure 3 gain table")
    cue = cue.drop(columns="_merge")
    if set(cue["algorithm"].unique()) != {"GIMME", "iMAT", "CORDA", "RIPTiDe"}:
        raise RuntimeError("Figure 3 gain table contains a method outside the current primary panel")
    cue["anchor_setting"] = cue["anchor_setting"].replace({"A2-L": "A2"})
    grouped = cue.groupby(
        ["anchor_setting", "algorithm", "rna_context_key", "cue_direction"], sort=True
    )["gain_mean"]
    context = grouped.agg(
        n_reaction_evaluations="size", mean_gain="mean", median_gain="median",
        q05=lambda values: values.quantile(0.05),
        q95=lambda values: values.quantile(0.95),
    ).reset_index()
    context["gain_unit"] = "mmol gDW^-1 h^-1"
    context = context.sort_values(
        ["algorithm", "rna_context_key", "anchor_setting", "cue_direction"], kind="stable"
    ).reset_index(drop=True)
    if context.duplicated(["anchor_setting", "algorithm", "rna_context_key", "cue_direction"]).any():
        raise RuntimeError("Duplicate context gain summary keys")
    validation = {
        "input_rows": int(len(cue)),
        "unique_evaluations": int(cue["evaluation_id"].nunique()),
        "unique_reactions": int(cue["reaction_id"].nunique()),
        "anchor_rows": cue.groupby("anchor_setting").size().astype(int).to_dict(),
        "cue_endpoint_rows": {
            f"{anchor}|{direction}": int(count)
            for (anchor, direction), count in cue.groupby(["anchor_setting", "cue_direction"]).size().items()
        },
        "output_context_summaries": int(len(context)),
        "reaction_evaluation_rows": int(len(wide)),
        "maximum_abs_g_info_identity_error": g_info_error,
        "maximum_abs_directional_advantage_identity_error": direction_error,
        "context_summary_grain": "anchor x algorithm x RNA context x cue; reaction-evaluation gains summarized descriptively",
        "quantile_interpretation": "descriptive q05/q95 across reaction-evaluation rows; not confidence intervals",
    }
    return context, validation


def build_outcome_tables(manifests: dict[str, dict]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    outcome_rows: list[dict] = []
    coverage_rows: list[dict] = []
    input_rows: dict[str, int] = {}
    for anchor, spec in SOURCES.items():
        path = repo_file(f"{spec['directory']}/{spec['prefix']}_CONTEXT_SUMMARY.tsv")
        frame = pd.read_csv(path, sep="\t")
        input_rows[anchor] = int(len(frame))
        overall = frame.loc[
            frame["summary_scope"].eq("OVERALL")
            & frame["algorithm"].eq("ALL")
            & frame["lambda_total"].eq(LAMBDA)
            & frame["lambda_role"].eq("PRIMARY")
            & frame["reliability_q"].isin(Q_ENDPOINTS)
        ].copy()
        if set(overall["reliability_q"].astype(float)) != set(Q_ENDPOINTS):
            raise RuntimeError(f"Missing primary q endpoints for {anchor}")
        if overall["reliability_q"].duplicated().any():
            raise RuntimeError(f"Duplicate overall primary q rows for {anchor}")
        for row in overall.to_dict(orient="records"):
            non_ties = int(row["non_tie_count"])
            total = int(row["total_distinct_reaction_truth_pair_cases"])
            improved, tied, harmed = (int(row[key]) for key in ("improved_count", "tied_count", "harmed_count"))
            if improved + tied + harmed != non_ties:
                raise RuntimeError(f"Improved/tied/harmed cases do not sum to the non-tie denominator for {anchor}")
            if int(row["truth_tie_count"]) + non_ties != total:
                raise RuntimeError(f"Truth tie and non-tie populations do not sum for {anchor}")
            q = float(row["reliability_q"])
            outcome_rows.append({
                "anchor_setting": anchor,
                "display_label": "A1 — 1-Strong-Anchor (HEX1)" if anchor == "A1" else "A2 — 2-Strong-Anchors (HEX1 + LDH_L)",
                "lambda_total": LAMBDA,
                "reliability_q": q,
                "cue_label": Q_ENDPOINTS[q],
                "non_tie_denominator": non_ties,
                "improved_count": improved,
                "improved_fraction": float(row["improved_fraction_of_non_ties"]),
                "gain_tie_count": tied,
                "gain_tie_fraction": float(row["tied_fraction_of_non_ties"]),
                "harmed_count": harmed,
                "harmed_fraction": float(row["harmed_fraction_of_non_ties"]),
                "fraction_unit": "proportion",
            })
            coverage_rows.append({
                "anchor_setting": anchor,
                "display_label": "A1 — 1-Strong-Anchor (HEX1)" if anchor == "A1" else "A2 — 2-Strong-Anchors (HEX1 + LDH_L)",
                "lambda_total": LAMBDA,
                "total_cases": total,
                "truth_tie_count": int(row["truth_tie_count"]),
                "truth_tie_fraction": int(row["truth_tie_count"]) / total,
                "non_tie_count": non_ties,
                "non_tie_fraction": non_ties / total,
            })
    outcomes = pd.DataFrame(outcome_rows).sort_values(["anchor_setting", "reliability_q"]).reset_index(drop=True)
    coverage = pd.DataFrame(coverage_rows).drop_duplicates("anchor_setting").sort_values("anchor_setting").reset_index(drop=True)
    for anchor, group in outcomes.groupby("anchor_setting"):
        total = group[["improved_fraction", "gain_tie_fraction", "harmed_fraction"]].sum(axis=1)
        if not np.allclose(total, 1.0, rtol=0, atol=1e-12):
            raise RuntimeError(f"Outcome proportions do not sum to one for {anchor}")
    return outcomes, coverage, {"context_summary_input_rows": input_rows}


def main() -> None:
    manifests, source_hashes = verify_sources()
    context_gains, gain_validation = build_context_gain_table()
    outcomes, coverage, outcome_validation = build_outcome_tables(manifests)
    DATA.mkdir(parents=True, exist_ok=True)
    outputs = {
        "supp_fig5_context_gains.tsv": context_gains,
        "supp_fig5_outcome_composition.tsv": outcomes,
        "supp_fig5_truth_tie_coverage.tsv": coverage,
    }
    output_hashes: dict[str, str] = {}
    for name, frame in outputs.items():
        path = DATA / name
        frame.to_csv(path, sep="\t", index=False, float_format="%.12g")
        output_hashes[name] = sha256(path)
    metadata = {
        "schema": "recomb.supp_fig5.build_manifest.v1",
        "status": "PASS",
        "inputs_sha256": source_hashes,
        "output_sha256": output_hashes,
        "output_rows": {name: int(len(frame)) for name, frame in outputs.items()},
        "parameters": {"primary_lambda": LAMBDA, "reliability_endpoints": {str(k): v for k, v in Q_ENDPOINTS.items()}},
        "matched_gain_population": {
            "matched_truth_pairs": 836,
            "reaction_universe": 4179,
            "input_gain_table_rows_per_anchor": 660237,
            "source_population": "Figure 3 compact matched A1/A2 gain tables",
        },
        "full_source_population_note": "PL2B A1 has 4,180 reactions including the additional LDH_L reaction; A22 A2 has 4,179. Overall source-composition tables preserve their source denominators and are not labeled as a matched A1/A2 comparison.",
        "gain_validation": gain_validation,
        "outcome_validation": outcome_validation,
    }
    (DATA / "supp_fig5_build_manifest.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": "PASS", "output_rows": metadata["output_rows"], "output_sha256": output_hashes}, sort_keys=True))


if __name__ == "__main__":
    main()
