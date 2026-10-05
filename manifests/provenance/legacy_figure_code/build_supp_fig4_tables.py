#!/usr/bin/env python3
"""Build compact S4 plotting tables from frozen PL2C/PL2D production artifacts."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures/supp_fig4"
DATA = OUT / "data"
PL2C = ROOT / "outputs/dmi_bridge_pl2c_geometry_utility_development_v1"
PL2DH = ROOT / "outputs/dmi_bridge_pl2d_hypothesis_freeze_v1"
PL2D = ROOT / "outputs/dmi_bridge_pl2d_confirmation_v1"
M1 = ROOT / "outputs/dmi_bridge_m1_manuscript_evidence_v1"
METHOD_ORDER = ("GIMME", "iMAT", "CORDA", "RIPTiDe")
CONTEXT_ORDER = (
    "training_samples=setx1,setx2",
    "training_samples=setx1,setx2,setx3",
    "training_samples=setx1,setx3",
    "training_samples=setx2,setx3",
)
EXPECTED = {
    "development": {"split_n": 3135, "finite_n": 2455, "rho": 0.8570601823118907},
    "confirmation": {"split_n": 1045, "finite_n": 838, "rho": 0.8351435823183292},
    "bootstrap": {"replicates": 5000, "finite_replicates": 5000, "seed": 20260929,
                  "rho_q025": 0.8015629118742916, "rho_median": 0.8354193419051692,
                  "rho_q975": 0.8631275423376382},
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def finite_numeric(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series.replace("", np.nan), errors="coerce")
    return values.where(np.isfinite(values), np.nan)


def close(actual: float, expected: float, label: str, tol: float = 1e-12) -> None:
    if not math.isfinite(float(actual)) or abs(float(actual) - expected) > tol:
        raise RuntimeError(f"{label} mismatch: {actual!r} != {expected!r}")


def main() -> None:
    sys.path.insert(0, str(ROOT))
    from scripts import dmi_bridge_m1_manuscript_core_v1 as m1core
    from scripts import dmi_bridge_pl2d_hypothesis_freeze_core_v1 as hcore

    # Verify the canonical manifests and every artifact listed by the PL2C,
    # PL2D-H, and PL2D production manifests before using production values.
    manifest_specs = (
        (PL2C, "BRIDGEPL2C_MANIFEST.json", m1core.PL2C_MANIFEST_SHA256),
        (PL2DH, "BRIDGEPL2DH_MANIFEST.json", m1core.PL2DH_MANIFEST_SHA256),
        (PL2D, "BRIDGEPL2D_MANIFEST.json", m1core.PL2D_MANIFEST_SHA256),
    )
    source_audit = {"status": "PASS", "bundles": {}}
    for base, name, expected_manifest_hash in manifest_specs:
        manifest_path = base / name
        actual_manifest_hash = sha256(manifest_path)
        if actual_manifest_hash != expected_manifest_hash:
            raise RuntimeError(f"Production manifest hash mismatch: {manifest_path.relative_to(ROOT)}")
        manifest = read_json(manifest_path)
        artifacts = manifest.get("artifact_sha256")
        if not isinstance(artifacts, dict) or not artifacts:
            raise RuntimeError(f"Manifest has no artifact inventory: {manifest_path.relative_to(ROOT)}")
        artifact_hashes = {}
        for filename, expected_hash in sorted(artifacts.items()):
            if Path(filename).name != filename:
                raise RuntimeError(f"Unsafe manifest artifact path: {filename}")
            path = base / filename
            actual_hash = sha256(path)
            if actual_hash != expected_hash:
                raise RuntimeError(f"Production artifact hash mismatch: {path.relative_to(ROOT)}")
            artifact_hashes[filename] = actual_hash
        source_audit["bundles"][name] = {
            "manifest_sha256": actual_manifest_hash, "status": manifest.get("status"),
            "artifact_sha256": artifact_hashes,
        }
    if source_audit["bundles"]["BRIDGEPL2D_MANIFEST.json"]["status"] != "PL2D_CONFIRMED":
        raise RuntimeError("PL2D production manifest is not confirmed")
    # M1 is a frozen derived summary, used only as a parity cross-check.
    m1_manifest_path = M1 / "BRIDGEM1_MANIFEST.json"
    m1_manifest = read_json(m1_manifest_path)
    if m1_manifest.get("status") != "BRIDGE_M1_MANUSCRIPT_EVIDENCE_FROZEN":
        raise RuntimeError("M1 derived evidence bundle is not frozen")
    for filename, expected_hash in m1_manifest.get("artifact_sha256", {}).items():
        if Path(filename).name != filename or sha256(M1 / filename) != expected_hash:
            raise RuntimeError(f"M1 derived artifact hash mismatch: {filename}")

    split_path = PL2C / "BRIDGEPL2C_REACTION_SPLIT.tsv"
    dev_path = PL2C / "BRIDGEPL2C_DEVELOPMENT_EVALUATION_REACTION.tsv.xz"
    snapshot_path = PL2DH / "BRIDGEPL2DH_DEVELOPMENT_SNAPSHOT.json"
    primary_path = PL2D / "BRIDGEPL2D_PRIMARY_RESULT.json"
    context_path = PL2D / "BRIDGEPL2D_CONTEXT_RESULTS.tsv"
    supportive_path = PL2D / "BRIDGEPL2D_SUPPORTIVE_RESULTS.json"
    bootstrap_path = PL2D / "BRIDGEPL2D_BOOTSTRAP.tsv.xz"
    qc_path = PL2D / "BRIDGEPL2D_QC.json"

    split = read_tsv(split_path)
    if set(split.columns) != {"reaction_id", "split_hash", "split_rank", "analysis_split"}:
        raise RuntimeError(f"Unexpected split registry columns: {list(split.columns)}")
    if split["reaction_id"].duplicated().any() or split["reaction_id"].eq("").any():
        raise RuntimeError("Reaction split has missing or duplicate identifiers")
    split_sets = {name: set(split.loc[split["analysis_split"].eq(name), "reaction_id"])
                  for name in ("DEVELOPMENT", "CONFIRMATION_HOLDOUT")}
    if len(split_sets["DEVELOPMENT"]) != 3135 or len(split_sets["CONFIRMATION_HOLDOUT"]) != 1045:
        raise RuntimeError("Frozen reaction split counts do not match PL2C contract")
    if split_sets["DEVELOPMENT"] & split_sets["CONFIRMATION_HOLDOUT"]:
        raise RuntimeError("Development and confirmation reaction populations overlap")

    # Development PL2C rows are evaluation x reaction; reproduce the frozen
    # pandas reaction-mean aggregation and SciPy Spearman definition.
    dev = pd.read_csv(dev_path, sep="\t", compression="xz", usecols=[
        "evaluation_id", "reaction_id", "sign_magnitude_eta2", "directionally_useful_fraction"
    ])
    dev_input_rows = int(len(dev))
    if dev_input_rows != 3135 * 370 or dev["evaluation_id"].nunique() != 370:
        raise RuntimeError("Unexpected PL2C development evaluation-reaction population")
    if dev["evaluation_id"].isna().any() or dev["reaction_id"].isna().any():
        raise RuntimeError("PL2C development table has missing observation keys")
    if set(dev["reaction_id"].unique()) != split_sets["DEVELOPMENT"]:
        raise RuntimeError("PL2C evaluation-reaction table does not match the frozen development split")
    if dev.duplicated(["evaluation_id", "reaction_id"]).any():
        raise RuntimeError("Duplicate PL2C evaluation-reaction keys")
    for col in ("sign_magnitude_eta2", "directionally_useful_fraction"):
        dev[col] = finite_numeric(dev[col])
        if ((dev[col].dropna() < -1e-12) | (dev[col].dropna() > 1 + 1e-12)).any():
            raise RuntimeError(f"PL2C {col} has values outside its valid [0, 1] range")
    dev_reaction = dev.groupby("reaction_id", sort=True)[
        ["sign_magnitude_eta2", "directionally_useful_fraction"]
    ].mean()
    dev_pair = dev_reaction.dropna(subset=["sign_magnitude_eta2", "directionally_useful_fraction"])
    dev_rho = hcore.spearman(dev_pair["sign_magnitude_eta2"], dev_pair["directionally_useful_fraction"])
    if len(dev_pair) != EXPECTED["development"]["finite_n"]:
        raise RuntimeError(f"Development finite pair count mismatch: {len(dev_pair)}")
    close(dev_rho, EXPECTED["development"]["rho"], "development rho")

    snapshot = read_json(snapshot_path)
    m1core.validate_development_snapshot({
        "finite_primary_reactions": snapshot["finite_primary_reactions"],
        "eta2_usefulness_rho": snapshot["primary_reaction_level_rho"],
        "positive_contexts": snapshot["positive_algorithm_rna_contexts"],
    })
    if snapshot["development_reactions"] != len(split_sets["DEVELOPMENT"]):
        raise RuntimeError("Development snapshot split size mismatch")
    if snapshot["finite_primary_reactions"] != len(dev_pair):
        raise RuntimeError("Development snapshot finite count disagrees with recomputation")
    close(snapshot["primary_reaction_level_rho"], dev_rho, "snapshot development rho")

    primary = read_json(primary_path)
    supportive = read_json(supportive_path)
    qc = read_json(qc_path)
    m1core.validate_confirmation_primary(primary)
    m1core.validate_supportive(supportive)
    m1core.validate_confirmation_qc(qc)
    if primary["finite_reactions"] != EXPECTED["confirmation"]["finite_n"]:
        raise RuntimeError("Unexpected confirmation finite pair count")
    close(primary["overall_rho"], EXPECTED["confirmation"]["rho"], "confirmation rho")
    if primary["status"] != "PL2D_CONFIRMED":
        raise RuntimeError("PL2D primary status is not confirmed")

    # Round-trip float parsing preserves exact serialized values and their tie structure.
    confirmation = pd.read_csv(PL2D / "BRIDGEPL2D_CONFIRMATION_REACTION.tsv", sep="\t", float_precision="round_trip")
    if len(confirmation) != 1045 or confirmation["reaction_id"].duplicated().any():
        raise RuntimeError("Unexpected PL2D reaction table row count or duplicate keys")
    if set(confirmation["reaction_id"]) != split_sets["CONFIRMATION_HOLDOUT"]:
        raise RuntimeError("PL2D confirmation reaction table does not match frozen holdout")
    x = finite_numeric(confirmation["mean_sign_magnitude_eta2"])
    y = finite_numeric(confirmation["mean_directionally_useful_fraction"])
    for label, values in (("eta2", x), ("directionally useful fraction", y)):
        if ((values.dropna() < -1e-12) | (values.dropna() > 1 + 1e-12)).any():
            raise RuntimeError(f"PL2D {label} has values outside its valid [0, 1] range")
    confirmation_pair = x.notna() & y.notna()
    recomputed_confirmation_rho = hcore.spearman(x[confirmation_pair], y[confirmation_pair])
    if int(confirmation_pair.sum()) != EXPECTED["confirmation"]["finite_n"]:
        raise RuntimeError("Confirmation finite pair count differs from frozen expected count")
    close(recomputed_confirmation_rho, primary["overall_rho"], "recomputed confirmation rho")

    context = read_tsv(context_path)
    expected_context_keys = {(m, c) for m in METHOD_ORDER for c in CONTEXT_ORDER}
    observed_context_keys = set(zip(context["algorithm"], context["rna_context_key"]))
    if len(context) != 16 or observed_context_keys != expected_context_keys:
        raise RuntimeError("Context result inventory does not contain exactly the frozen 16 contexts")
    context["rho"] = pd.to_numeric(context["rho"], errors="raise")
    context_n = pd.to_numeric(context["finite_primary_reaction_pairs"], errors="raise").astype(int)
    if not np.isfinite(context["rho"]).all() or (context_n < 30).any():
        raise RuntimeError("Context correlation is nonfinite or below the frozen support gate")
    if (context["rho"] > 0).sum() != 16:
        raise RuntimeError("Expected all 16 context correlations to be positive")
    primary_context = {(r["algorithm"], r["rna_context_key"]): r for r in primary["context_results"]}
    for row in context.to_dict("records"):
        key = (row["algorithm"], row["rna_context_key"])
        frozen = primary_context.get(key)
        if frozen is None or int(row["finite_primary_reaction_pairs"]) != int(frozen["finite_primary_reaction_pairs"]):
            raise RuntimeError(f"Context table and primary result disagree for {key}")
        close(float(row["rho"]), float(frozen["rho"]), f"context rho {key}")
    context["method_order"] = context["algorithm"].map({m: i for i, m in enumerate(METHOD_ORDER)})
    context["context_order"] = context["rna_context_key"].map({c: i for i, c in enumerate(CONTEXT_ORDER)})
    context = context.sort_values(["method_order", "context_order"]).copy()
    context["context_label"] = context["rna_context_key"].str.removeprefix("training_samples=").str.replace(",", " + ", regex=False)
    context.rename(columns={"finite_primary_reaction_pairs": "finite_n", "total_confirmation_reactions": "split_n"}, inplace=True)
    context["finite_n"] = pd.to_numeric(context["finite_n"], errors="raise").astype(int)
    context["rho"] = pd.to_numeric(context["rho"], errors="raise")
    context_out = context[["algorithm", "context_label", "rna_context_key", "split_n", "finite_n", "rho", "direction", "evaluable"]]

    # The raw replicate table itself is the figure's bootstrap distribution;
    # validate against PL2D's frozen metadata and recomputed percentile summary.
    bootstrap = pd.read_csv(bootstrap_path, sep="\t", compression="xz")
    bmeta = supportive["bootstrap"]
    if len(bootstrap) != EXPECTED["bootstrap"]["replicates"] or bootstrap["replicate_index"].duplicated().any():
        raise RuntimeError("Bootstrap replicate table has unexpected size or duplicate indices")
    if bootstrap["replicate_index"].min() != 0 or bootstrap["replicate_index"].max() != 4999:
        raise RuntimeError("Bootstrap replicate indices do not span 0 through 4999")
    if not np.isfinite(bootstrap["rho"].to_numpy(float)).all():
        raise RuntimeError("Unexpected nonfinite bootstrap replicate")
    if ((bootstrap["rho"] < -1 - 1e-12) | (bootstrap["rho"] > 1 + 1e-12)).any():
        raise RuntimeError("Bootstrap rho outside the valid [-1, 1] range")
    if int(bmeta["replicates"]) != 5000 or int(bmeta["finite_replicates"]) != 5000 or int(bmeta["seed"]) != 20260929:
        raise RuntimeError("Bootstrap provenance metadata mismatch")
    if bmeta["rng"] != "numpy.random.Generator(PCG64)":
        raise RuntimeError("Bootstrap RNG metadata mismatch")
    quantiles = np.quantile(bootstrap["rho"].to_numpy(float), [0.025, 0.5, 0.975], method="linear")
    for actual, key in zip(quantiles, ("rho_q025", "rho_median", "rho_q975")):
        close(actual, float(bmeta[key]), f"bootstrap {key}")
        close(actual, EXPECTED["bootstrap"][key], f"expected bootstrap {key}")

    summary = pd.DataFrame([
        {"cohort": "DEVELOPMENT", "split_n": 3135, "finite_n": int(len(dev_pair)), "rho": dev_rho,
         "status": "FROZEN_VALUES_REPRODUCED", "predictor": "mean_sign_magnitude_eta2",
         "response": "mean_directionally_useful_fraction"},
        {"cohort": "CONFIRMATION", "split_n": 1045, "finite_n": int(confirmation_pair.sum()), "rho": recomputed_confirmation_rho,
         "status": primary["status"], "predictor": "mean_sign_magnitude_eta2",
         "response": "mean_directionally_useful_fraction"},
    ])
    m1_dc = read_tsv(M1 / "BRIDGEM1_DEVELOPMENT_CONFIRMATION.tsv").set_index("cohort")
    for cohort, n, rho in (("DEVELOPMENT", 2455, dev_rho), ("CONFIRMATION", 838, recomputed_confirmation_rho)):
        if int(m1_dc.loc[cohort, "finite_reactions"]) != n:
            raise RuntimeError(f"M1 summary finite count mismatch for {cohort}")
        close(float(m1_dc.loc[cohort, "eta2_usefulness_rho"]), rho, f"M1 summary rho {cohort}")
    DATA.mkdir(parents=True, exist_ok=True)
    summary.to_csv(DATA / "supp_fig4_correlations.tsv", sep="\t", index=False, float_format="%.17g")
    context_out.to_csv(DATA / "supp_fig4_contexts.tsv", sep="\t", index=False, float_format="%.17g")
    with gzip.open(DATA / "supp_fig4_confirmation_bootstrap.tsv.gz", "wt", encoding="utf-8", newline="") as f:
        bootstrap[["replicate_index", "finite_reaction_pair_count", "rho"]].to_csv(f, sep="\t", index=False, float_format="%.17g")

    inputs = [split_path, dev_path, snapshot_path, primary_path, context_path, supportive_path, bootstrap_path, qc_path,
              PL2C / "BRIDGEPL2C_MANIFEST.json", PL2D / "BRIDGEPL2D_MANIFEST.json",
              PL2DH / "BRIDGEPL2DH_MANIFEST.json", M1 / "BRIDGEM1_MANIFEST.json",
              M1 / "BRIDGEM1_DEVELOPMENT_CONFIRMATION.tsv"]
    hashes = {str(p.relative_to(ROOT)): sha256(p) for p in inputs}
    validation = {
        "status": "PASS", "development": {"split_reactions": 3135, "input_evaluation_reaction_rows": dev_input_rows,
            "unique_reactions": int(dev["reaction_id"].nunique()), "finite_reaction_pairs": int(len(dev_pair)), "rho": float(dev_rho),
            "excluded_nonfinite_pairs": int(len(dev_reaction) - len(dev_pair))},
        "confirmation": {"split_reactions": 1045, "reaction_rows": int(len(confirmation)),
            "finite_reaction_pairs": int(confirmation_pair.sum()), "rho": float(recomputed_confirmation_rho),
            "excluded_nonfinite_pairs": int((~confirmation_pair).sum()), "status": primary["status"]},
        "contexts": {"rows": int(len(context_out)), "finite_pair_min": int(context_out["finite_n"].min()),
            "finite_pair_max": int(context_out["finite_n"].max()), "positive_correlations": int((context_out["rho"] > 0).sum())},
        "bootstrap": {"replicates": int(len(bootstrap)), "finite_replicates": int(np.isfinite(bootstrap["rho"]).sum()),
            "seed": 20260929, "rng": bmeta["rng"], "interval": "2.5th and 97.5th percentiles; NumPy quantile method=linear",
            "q025": float(quantiles[0]), "median": float(quantiles[1]), "q975": float(quantiles[2]),
            "interpretation": "computational reaction-level stability, not biological replication"},
        "exclusions": [
            {"rule": "pairwise finite predictor and response required for Spearman correlation", "development_count": int(len(dev_reaction) - len(dev_pair)),
             "confirmation_count": int((~confirmation_pair).sum()), "reason": "correlation undefined when either value is nonfinite"},
            {"rule": "no exclusions from the confirmation bootstrap distribution", "count": 0, "reason": "all 5,000 frozen replicate rho values are finite"}],
        "development_confirmation_overlap": 0,
        "source_audit": source_audit,
        "input_sha256": hashes,
    }
    (DATA / "supp_fig4_validation_summary.json").write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n")
    out_hashes = {p.name: sha256(p) for p in (DATA / "supp_fig4_correlations.tsv", DATA / "supp_fig4_contexts.tsv",
                                               DATA / "supp_fig4_confirmation_bootstrap.tsv.gz", DATA / "supp_fig4_validation_summary.json")}
    manifest = {"status": "PASS", "schema": "supp_fig4.tables.v1", "inputs": hashes,
                "output_sha256": out_hashes, "development_confirmation_pooled": False,
                "production_pipeline_rerun": False}
    (DATA / "supp_fig4_build_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "PASS", "development": {"n": len(dev_pair), "rho": dev_rho},
                      "confirmation": {"n": int(confirmation_pair.sum()), "rho": recomputed_confirmation_rho},
                      "contexts": len(context_out), "bootstrap_replicates": len(bootstrap)}, indent=2))


if __name__ == "__main__":
    main()
