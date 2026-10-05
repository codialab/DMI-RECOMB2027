#!/usr/bin/env python3
"""Publish the development-only DMI-BRIDGE PL2D-H hypothesis freeze."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

try:
    from scripts import dmi_bridge_pl2d_hypothesis_freeze_core_v1 as core
except ModuleNotFoundError:  # direct script execution
    import dmi_bridge_pl2d_hypothesis_freeze_core_v1 as core


ROOT = Path(__file__).resolve().parents[1]
PL2C_REL = Path("outputs/dmi_bridge_pl2c_geometry_utility_development_v1")
OUTPUT_REL = Path("outputs/dmi_bridge_pl2d_hypothesis_freeze_v1")
PATCH_REL = Path("dmi_bridge_pl2d_hypothesis_freeze_foundation_v2.patch")
PATCH_SHA256 = "7d74dc0f0590e5d06e6c77135f5e97fb4bbd04eb27457629d02cec8c728c57f8"
CONTRACT_REL = Path("docs/DMI_BRIDGE_PL2D_HYPOTHESIS_FREEZE_CONTRACT.md")
MANIFEST_NAME = "BRIDGEPL2C_MANIFEST.json"
SPLIT_NAME = "BRIDGEPL2C_REACTION_SPLIT.tsv"
DEVELOPMENT_NAME = "BRIDGEPL2C_DEVELOPMENT_EVALUATION_REACTION.tsv.xz"
ARTIFACT_NAMES = (
    "BRIDGEPL2DH_SOURCE_AUDIT.json",
    "BRIDGEPL2DH_DEVELOPMENT_SNAPSHOT.json",
    "BRIDGEPL2DH_CONFIRMATION_REGISTRY.tsv",
    "BRIDGEPL2DH_CONFIRMATION_CONTRACT.json",
)
NUMERIC_COLUMNS = core.REACTION_MEAN_SOURCES
DEVELOPMENT_COLUMNS = (
    "evaluation_id", "reaction_id", "algorithm", "rna_context_key", *NUMERIC_COLUMNS
)


def _regular_repo_file(path: Path) -> Path:
    resolved = path.resolve(strict=True)
    if not resolved.is_relative_to(ROOT) or path.is_symlink() or not resolved.is_file():
        raise RuntimeError(f"source is not a regular in-repository file: {path}")
    return resolved


def _source_file(name: str) -> Path:
    if name not in {MANIFEST_NAME, *core.EXPECTED_PL2C_ARTIFACT_SHA256}:
        raise RuntimeError(f"PL2C source not allowlisted: {name}")
    base = ROOT / PL2C_REL
    if base.is_symlink() or base.resolve(strict=True) != base:
        raise RuntimeError("canonical PL2C bundle is redirected")
    return _regular_repo_file(base / name)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with _regular_repo_file(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_bytes(value: dict) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def _versions() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }


def verify_pl2c() -> tuple[dict, dict]:
    patch_hash = _sha256(ROOT / PATCH_REL)
    if patch_hash != PATCH_SHA256:
        raise RuntimeError("PL2D-H foundation patch SHA256 mismatch")
    manifest_path = _source_file(MANIFEST_NAME)
    manifest_hash = _sha256(manifest_path)
    if manifest_hash != core.EXPECTED_PL2C_MANIFEST_SHA256:
        raise RuntimeError("canonical PL2C manifest SHA256 mismatch")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != core.EXPECTED_PL2C_STATUS:
        raise RuntimeError("canonical PL2C status mismatch")
    if manifest.get("pl2d_confirmation_status") != "UNTOUCHED":
        raise RuntimeError("PL2C confirmation is not untouched")
    if manifest.get("confirmation_outcomes_analyzed") is not False:
        raise RuntimeError("PL2C confirmation outcomes flag is not false")
    if manifest.get("reaction_split_counts") != core.EXPECTED_SPLIT_COUNTS:
        raise RuntimeError("canonical PL2C reaction split counts mismatch")
    if manifest.get("split_registry_sha256") != core.EXPECTED_SPLIT_SHA256:
        raise RuntimeError("canonical PL2C split registry hash mismatch")
    if manifest.get("artifact_sha256") != core.EXPECTED_PL2C_ARTIFACT_SHA256:
        raise RuntimeError("canonical PL2C artifact inventory mismatch")
    for name, expected in core.EXPECTED_PL2C_ARTIFACT_SHA256.items():
        if _sha256(_source_file(name)) != expected:
            raise RuntimeError(f"canonical PL2C artifact SHA256 mismatch: {name}")
    audit = {
        "schema": "bridge.pl2d_h.source_audit.v1",
        "status": "PASS",
        "pl2d_h_foundation_patch_sha256": patch_hash,
        "pl2c_manifest_sha256": manifest_hash,
        "pl2c_artifact_sha256": core.EXPECTED_PL2C_ARTIFACT_SHA256,
        "pl2c_split_registry_sha256": core.EXPECTED_SPLIT_SHA256,
        "pl2c_status": manifest["status"],
        "pl2c_confirmation_status": manifest["pl2d_confirmation_status"],
        "pl2c_confirmation_outcomes_analyzed": False,
        "reaction_split_counts": core.EXPECTED_SPLIT_COUNTS,
        "parsed_pl2c_files": [MANIFEST_NAME, DEVELOPMENT_NAME, SPLIT_NAME],
        "other_pl2c_artifacts": "SHA256_CHECK_ONLY",
        "confirmation_predictor_values_accessed": False,
        "confirmation_outcomes_accessed": False,
        "confirmation_geometry_utility_analysis_invoked": False,
        "numerical_versions": _versions(),
    }
    return manifest, audit


def _split_and_registry() -> tuple[set[str], bytes]:
    with _source_file(SPLIT_NAME).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or set(rows[0]) != {"reaction_id", "split_hash", "split_rank", "analysis_split"}:
        raise RuntimeError("PL2C split registry fields mismatch")
    registry = core.build_confirmation_registry(rows)
    blob = core.confirmation_registry_tsv_bytes(registry)
    if core.sha256_bytes(blob) != core.EXPECTED_CONFIRMATION_REGISTRY_SHA256:
        raise RuntimeError("confirmation identity registry SHA256 mismatch")
    development = {r["reaction_id"] for r in rows if r["analysis_split"] == "DEVELOPMENT"}
    if len(development) != core.EXPECTED_SPLIT_COUNTS["DEVELOPMENT"]:
        raise RuntimeError("development reaction identity count mismatch")
    return development, blob


def _association(rows: list[dict[str, object]], x: str, y: str) -> tuple[int, float]:
    pairs = [(r[x], r[y]) for r in rows if r[x] is not None and r[y] is not None]
    return len(pairs), core.spearman([a for a, _ in pairs], [b for _, b in pairs])


def _close(actual: float, expected: float, label: str) -> None:
    if not math.isfinite(actual) or not math.isclose(actual, expected, rel_tol=0.0, abs_tol=1e-12):
        raise RuntimeError(f"frozen development snapshot mismatch: {label}: {actual!r}")


def reproduce_development(development: set[str]) -> dict:
    frame = pd.read_csv(
        _source_file(DEVELOPMENT_NAME), sep="\t", usecols=list(DEVELOPMENT_COLUMNS),
        compression="xz", dtype=str, keep_default_na=False,
    )
    if len(frame) != core.EXPECTED_DEVELOPMENT_ROWS:
        raise RuntimeError("PL2C development row count mismatch")
    if set(frame["reaction_id"]) != development:
        raise RuntimeError("PL2C development table has wrong reaction identities")
    if frame["evaluation_id"].nunique() != core.EXPECTED_EVALUATIONS:
        raise RuntimeError("PL2C development evaluation count mismatch")
    if not (frame.groupby("reaction_id", sort=True).size() == core.EXPECTED_EVALUATIONS).all():
        raise RuntimeError("PL2C development evaluation coverage mismatch")
    if frame[["algorithm", "rna_context_key"]].drop_duplicates().shape[0] != core.EXPECTED_CONTEXTS:
        raise RuntimeError("PL2C development context count mismatch")

    reaction_rows = core.aggregate_reaction_frame(frame)
    observed: dict[str, object] = {
        "schema": "bridge.pl2d_h.development_snapshot.v1",
        "source": DEVELOPMENT_NAME,
        "development_reactions": len(development),
        "development_evaluation_reaction_rows": len(frame),
        "numerical_reference": core.NUMERICAL_REFERENCE,
        "numerical_versions": _versions(),
    }
    n, rho = core.primary_rho(reaction_rows)
    observed["finite_primary_reactions"] = n
    observed["primary_reaction_level_rho"] = rho
    if n != core.DEVELOPMENT_SNAPSHOT["finite_primary_reactions"]:
        raise RuntimeError("frozen primary reaction support mismatch")
    _close(rho, core.DEVELOPMENT_SNAPSHOT["primary_reaction_level_rho"], "primary rho")

    triple = [r for r in reaction_rows if all(r[k] is not None for k in (
        "mean_sign_magnitude_eta2", "mean_directionally_useful_fraction", "mean_non_tie_pair_fraction"
    ))]
    partial = core.partial_spearman_one_control(
        [r["mean_sign_magnitude_eta2"] for r in triple],
        [r["mean_directionally_useful_fraction"] for r in triple],
        [r["mean_non_tie_pair_fraction"] for r in triple],
    )
    observed["finite_partial_reactions"] = len(triple)
    observed["primary_partial_rho_controlling_non_tie_pair_fraction"] = partial
    _close(partial, core.DEVELOPMENT_SNAPSHOT["primary_partial_rho_controlling_non_tie_pair_fraction"], "partial rho")

    comparisons = (
        ("eta2_vs_mean_information_advantage_rho", "mean_sign_magnitude_eta2", "mean_information_advantage", "finite_information_advantage_reactions"),
        ("supportive_directional_entropy3_rho", "mean_directional_entropy3", "mean_directionally_useful_fraction", "finite_directional_entropy3_reactions"),
        ("supportive_dominant_sign_mass_rho", "mean_dominant_sign_mass", "mean_directionally_useful_fraction", "finite_dominant_sign_mass_reactions"),
        ("supportive_n_supported_sign_states_rho", "mean_n_supported_sign_states", "mean_directionally_useful_fraction", "finite_n_supported_sign_states_reactions"),
    )
    for name, x, y, count_name in comparisons:
        count, value = _association(reaction_rows, x, y)
        observed[name] = value
        observed[count_name] = count
        _close(value, core.DEVELOPMENT_SNAPSHOT[name], name)
        expected_count = n if name == "eta2_vs_mean_information_advantage_rho" else core.DEVELOPMENT_SNAPSHOT["supportive_finite_reactions"]
        if count != expected_count:
            raise RuntimeError(f"frozen reaction support mismatch: {name}")

    contexts = []
    for (algorithm, rna_context_key), subset in frame.groupby(["algorithm", "rna_context_key"], sort=True):
        within = core.aggregate_reaction_frame(subset)
        count, value = _association(within, "mean_sign_magnitude_eta2", "mean_directionally_useful_fraction")
        contexts.append({"algorithm": algorithm, "rna_context_key": rna_context_key, "finite_reactions": count, "rho": value})
    observed["algorithm_rna_contexts"] = contexts
    observed["positive_algorithm_rna_contexts"] = sum(math.isfinite(row["rho"]) and row["rho"] > 0 for row in contexts)
    if len(contexts) != core.EXPECTED_CONTEXTS or observed["positive_algorithm_rna_contexts"] != core.DEVELOPMENT_SNAPSHOT["positive_algorithm_rna_contexts"]:
        raise RuntimeError("frozen development context directions mismatch")
    observed["status"] = "FROZEN_VALUES_REPRODUCED"
    return observed


def _contract() -> dict:
    return {
        "schema": "bridge.pl2d_h.confirmation_contract.v1",
        "status": "FROZEN_BEFORE_CONFIRMATION",
        "pl2d_confirmation_status": "NOT_RUN",
        "primary_feature": "sign_magnitude_eta2",
        "primary_unit": "reaction",
        "reaction_split": "CONFIRMATION_HOLDOUT",
        "confirmation_reaction_count": 1045,
        "predictor_aggregation": "mean of all finite sign_magnitude_eta2 evaluation rows, independent of response availability",
        "primary_response_aggregation": "mean of finite directionally_useful_fraction values from rows with defined directional response",
        "coverage_aggregation": "mean finite non_tie_pair_fraction across evaluation rows",
        "secondary_magnitude_response_aggregation": "mean finite mean_information_advantage from rows with defined directional response",
        "primary_statistic": "reaction-level Spearman rho: mean_sign_magnitude_eta2 vs mean_directionally_useful_fraction",
        "primary": core.PRIMARY,
        "primary_minimum_rho": 0.50,
        "primary_minimum_positive_contexts": 12,
        "primary_minimum_finite_reactions": 700,
        "status_rule": {
            "PL2D_INSUFFICIENT_SUPPORT": "fewer than 700 finite reaction pairs or fewer than 16 evaluable contexts with at least 100 finite pairs each",
            "PL2D_NOT_CONFIRMED": "support adequate but rho < 0.50 or fewer than 12 positive contexts",
            "PL2D_CONFIRMED": "all primary support and effect requirements pass",
        },
        "supportive_analyses": [
            "eta2 vs mean information advantage: expected positive",
            "partial Spearman eta2 vs directional usefulness controlling mean non-tie coverage: expected positive",
            "directional entropy vs directional usefulness: expected positive",
            "dominant sign mass vs directional usefulness: expected negative",
            "supported sign-state count vs directional usefulness: expected positive",
            "algorithm x RNA context-specific associations",
            "proximal-reaction exclusion",
            "same-subsystem-as-HEX1 exclusion",
            "union pathway exclusion",
        ],
        "bootstrap": core.BOOTSTRAP | {"unit": "confirmation reaction ID", "resampling": "reaction IDs with replacement", "gate": False},
        "numerical_reference": core.NUMERICAL_REFERENCE,
        "interpretation_limits": [
            "PL2C development statistics are discovery evidence, not confirmation",
            "bootstrap is computational stability across reaction benchmark units, not biological replication or an independent biological confidence interval",
            "supportive analyses cannot rescue PL2D_NOT_CONFIRMED",
            "PL2C eta2 x entropy quartile map and entropy-Q4 group are not confirmation rules",
            "PL2D confirmation is a computational holdout, not independent biological validation",
        ],
        "confirmation_predictor_values_accessed": False,
        "confirmation_outcomes_accessed": False,
        "confirmation_geometry_utility_analysis_invoked": False,
    }


def _fingerprint() -> dict:
    return {
        "foundation_patch_sha256": _sha256(ROOT / PATCH_REL),
        "pl2c_manifest_sha256": core.EXPECTED_PL2C_MANIFEST_SHA256,
        "pl2c_artifact_sha256": core.EXPECTED_PL2C_ARTIFACT_SHA256,
        "implementation_sha256": {
            "scripts/dmi_bridge_pl2d_hypothesis_freeze_core_v1.py": _sha256(ROOT / "scripts/dmi_bridge_pl2d_hypothesis_freeze_core_v1.py"),
            "scripts/dmi_bridge_pl2d_hypothesis_freeze_v1.py": _sha256(ROOT / "scripts/dmi_bridge_pl2d_hypothesis_freeze_v1.py"),
        },
        "contract_document_sha256": _sha256(ROOT / CONTRACT_REL),
        "numerical_versions": _versions(),
    }


def validate_existing_output(output: Path, fingerprint: dict) -> dict:
    if output.is_symlink() or not output.is_dir():
        raise RuntimeError("existing PL2D-H output is not a regular directory")
    if {p.name for p in output.iterdir()} != set(ARTIFACT_NAMES) | {"BRIDGEPL2DH_MANIFEST.json"}:
        raise RuntimeError("existing PL2D-H artifact inventory mismatch")
    manifest = json.loads(_regular_repo_file(output / "BRIDGEPL2DH_MANIFEST.json").read_text(encoding="utf-8"))
    if manifest.get("status") != core.STATUS or manifest.get("fingerprint") != fingerprint:
        raise RuntimeError("existing PL2D-H status or fingerprint mismatch")
    if manifest.get("pl2d_confirmation_status") != "NOT_RUN" or manifest.get("confirmation_predictor_values_accessed") is not False or manifest.get("confirmation_outcomes_accessed") is not False or manifest.get("confirmation_geometry_utility_analysis_invoked") is not False:
        raise RuntimeError("existing PL2D-H confirmation firewall mismatch")
    hashes = manifest.get("artifact_sha256")
    if not isinstance(hashes, dict) or set(hashes) != set(ARTIFACT_NAMES):
        raise RuntimeError("existing PL2D-H artifact hashes missing")
    for name, digest in hashes.items():
        if _sha256(output / name) != digest:
            raise RuntimeError(f"existing PL2D-H artifact mutation: {name}")
    if hashes["BRIDGEPL2DH_CONFIRMATION_REGISTRY.tsv"] != core.EXPECTED_CONFIRMATION_REGISTRY_SHA256:
        raise RuntimeError("existing PL2D-H registry identity mismatch")
    return manifest


def produce(output: Path) -> dict:
    outputs = (ROOT / "outputs").resolve(strict=True)
    if output.is_symlink():
        raise RuntimeError("PL2D-H output may not be a symlink")
    output = output.resolve()
    if not output.is_relative_to(outputs) or output == outputs:
        raise RuntimeError("PL2D-H output must remain under repository outputs")
    _, audit = verify_pl2c()
    fingerprint = _fingerprint()
    if output.exists():
        return {"status": "NO_OP_EXISTING_IDENTICAL_PL2D_H", "manifest": validate_existing_output(output, fingerprint)}
    development, registry_blob = _split_and_registry()
    snapshot = reproduce_development(development)
    contract = _contract()
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="dmi_bridge_pl2d_h_publish_", dir=output.parent) as temporary:
        stage = Path(temporary) / output.name
        stage.mkdir()
        blobs = {
            "BRIDGEPL2DH_SOURCE_AUDIT.json": _json_bytes(audit),
            "BRIDGEPL2DH_DEVELOPMENT_SNAPSHOT.json": _json_bytes(snapshot),
            "BRIDGEPL2DH_CONFIRMATION_REGISTRY.tsv": registry_blob,
            "BRIDGEPL2DH_CONFIRMATION_CONTRACT.json": _json_bytes(contract),
        }
        for name, blob in blobs.items():
            (stage / name).write_bytes(blob)
        hashes = {name: core.sha256_bytes(blob) for name, blob in blobs.items()}
        manifest = {
            "schema": core.SCHEMA,
            "status": core.STATUS,
            "fingerprint": fingerprint,
            "artifact_sha256": hashes,
            "pl2c_manifest_sha256": core.EXPECTED_PL2C_MANIFEST_SHA256,
            "pl2c_artifact_sha256": core.EXPECTED_PL2C_ARTIFACT_SHA256,
            "pl2c_split_registry_sha256": core.EXPECTED_SPLIT_SHA256,
            "confirmation_registry_sha256": hashes["BRIDGEPL2DH_CONFIRMATION_REGISTRY.tsv"],
            "confirmation_reaction_count": 1045,
            "confirmation_predictor_values_accessed": False,
            "confirmation_outcomes_accessed": False,
            "confirmation_geometry_utility_analysis_invoked": False,
            "pl2d_confirmation_status": "NOT_RUN",
            "primary_feature": "sign_magnitude_eta2",
            "primary_unit": "reaction",
            "primary_minimum_rho": 0.50,
            "primary_minimum_positive_contexts": 12,
            "primary_minimum_finite_reactions": 700,
            "numerical_versions": _versions(),
        }
        (stage / "BRIDGEPL2DH_MANIFEST.json").write_bytes(_json_bytes(manifest))
        os.replace(stage, output)
    validate_existing_output(output, fingerprint)
    return {"status": core.STATUS, "manifest": manifest}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / OUTPUT_REL)
    args = parser.parse_args()
    result = produce(args.output)
    print(json.dumps({
        "status": result["status"],
        "confirmation_reaction_count": result["manifest"]["confirmation_reaction_count"],
        "confirmation_registry_sha256": result["manifest"]["confirmation_registry_sha256"],
        "pl2d_confirmation_status": result["manifest"]["pl2d_confirmation_status"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
