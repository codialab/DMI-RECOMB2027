#!/usr/bin/env python3
"""Publish the development-only PL2D-H support-gate amendment."""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import tempfile
from pathlib import Path

try:
    from scripts import dmi_bridge_pl2d_h_support_gate_amendment_core_v1 as core
except ModuleNotFoundError:  # direct script execution
    import dmi_bridge_pl2d_h_support_gate_amendment_core_v1 as core


ROOT = Path(__file__).resolve().parents[1]
SOURCE_REL = Path("outputs/dmi_bridge_pl2d_hypothesis_freeze_v1")
OUTPUT_REL = Path("outputs/dmi_bridge_pl2d_h_support_gate_amendment_v1")
PATCH_REL = Path("dmi_bridge_pl2d_h_support_gate_amendment_foundation_v1.patch")
DOC_REL = Path("docs/DMI_BRIDGE_PL2D_H_SUPPORT_GATE_AMENDMENT_V1.md")
CORE_REL = Path("scripts/dmi_bridge_pl2d_h_support_gate_amendment_core_v1.py")
ADAPTER_REL = Path("scripts/dmi_bridge_pl2d_h_support_gate_amendment_v1.py")
PATCH_SHA256 = "b7c53dd5e7848fb8a57d4a56d7804652cae445e943287be1cd232bef412741a3"
SOURCE_MANIFEST = "BRIDGEPL2DH_MANIFEST.json"
SOURCE_FILES = {SOURCE_MANIFEST, *core.ORIGINAL_ARTIFACT_SHA256}
ARTIFACTS = (
    "BRIDGEPL2DHA_SOURCE_AUDIT.json",
    "BRIDGEPL2DHA_SUPPORT_CALIBRATION.json",
    "BRIDGEPL2DHA_AMENDED_CONFIRMATION_CONTRACT.json",
)
MANIFEST = "BRIDGEPL2DHA_MANIFEST.json"
INSUFFICIENT_RULE = (
    "fewer than 700 finite overall reaction pairs, fewer than 16 context results, "
    "any of the 16 contexts with fewer than 30 finite reaction pairs, or non-finite overall primary Spearman rho"
)
ORIGINAL_INSUFFICIENT_RULE = (
    "fewer than 700 finite reaction pairs or fewer than 16 evaluable contexts with at least 100 finite pairs each"
)


def _file(path: Path) -> Path:
    resolved = path.resolve(strict=True)
    if path.is_symlink() or not resolved.is_relative_to(ROOT) or not resolved.is_file():
        raise RuntimeError(f"not a regular in-repository file: {path}")
    return resolved


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with _file(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _source_file(name: str) -> Path:
    if name not in SOURCE_FILES:
        raise RuntimeError(f"PL2D-H source not allowlisted: {name}")
    base = ROOT / SOURCE_REL
    if base.is_symlink() or base.resolve(strict=True) != base:
        raise RuntimeError("original PL2D-H bundle is redirected")
    return _file(base / name)


def _json_bytes(value: dict) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def _source_hashes() -> dict[str, str]:
    expected = {SOURCE_MANIFEST: core.ORIGINAL_MANIFEST_SHA256, **core.ORIGINAL_ARTIFACT_SHA256}
    observed = {name: _sha256(_source_file(name)) for name in sorted(expected)}
    if observed != expected:
        raise RuntimeError("original PL2D-H source hash mismatch")
    if {p.name for p in (ROOT / SOURCE_REL).iterdir()} != SOURCE_FILES:
        raise RuntimeError("original PL2D-H source inventory mismatch")
    return observed


def verify_source() -> tuple[dict, dict, dict, dict[str, str]]:
    hashes = _source_hashes()
    manifest = json.loads(_source_file(SOURCE_MANIFEST).read_text(encoding="utf-8"))
    if manifest.get("status") != core.ORIGINAL_STATUS or manifest.get("pl2d_confirmation_status") != "NOT_RUN":
        raise RuntimeError("original PL2D-H status mismatch")
    firewall = ("confirmation_predictor_values_accessed", "confirmation_outcomes_accessed", "confirmation_geometry_utility_analysis_invoked")
    if any(manifest.get(key) is not False for key in firewall):
        raise RuntimeError("original PL2D-H confirmation firewall mismatch")
    if manifest.get("artifact_sha256") != core.ORIGINAL_ARTIFACT_SHA256:
        raise RuntimeError("original PL2D-H artifact inventory mismatch")
    audit = json.loads(_source_file("BRIDGEPL2DH_SOURCE_AUDIT.json").read_text(encoding="utf-8"))
    if audit.get("status") != "PASS" or any(audit.get(key) is not False for key in firewall):
        raise RuntimeError("original PL2D-H source audit firewall mismatch")
    snapshot = json.loads(_source_file("BRIDGEPL2DH_DEVELOPMENT_SNAPSHOT.json").read_text(encoding="utf-8"))
    if snapshot.get("development_reactions") != core.DEVELOPMENT_REACTIONS:
        raise RuntimeError("development split size mismatch")
    core.verify_frozen_development_contexts(snapshot["algorithm_rna_contexts"])
    if manifest.get("confirmation_reaction_count") != core.CONFIRMATION_REACTIONS:
        raise RuntimeError("confirmation split size mismatch")
    contract = json.loads(_source_file("BRIDGEPL2DH_CONFIRMATION_CONTRACT.json").read_text(encoding="utf-8"))
    if contract.get("pl2d_confirmation_status") != "NOT_RUN" or any(contract.get(key) is not False for key in firewall):
        raise RuntimeError("original contract firewall mismatch")
    return manifest, snapshot, contract, hashes


def contract_diff(original: dict, amended: dict) -> list[str]:
    """Validate the whole effective contract, including every unchanged nested field."""
    expected = copy.deepcopy(original)
    if expected.get("status_rule", {}).get("PL2D_INSUFFICIENT_SUPPORT") != ORIGINAL_INSUFFICIENT_RULE:
        raise RuntimeError("original status rule mismatch")
    expected["primary"] = core.amended_primary(expected["primary"])
    expected["status_rule"]["PL2D_INSUFFICIENT_SUPPORT"] = INSUFFICIENT_RULE
    expected["supersedes_support_gate_from_manifest_sha256"] = core.ORIGINAL_MANIFEST_SHA256
    expected["amendment_scope"] = "minimum_finite_reactions_per_context_only"
    if amended != expected:
        raise RuntimeError("amended contract changes fields beyond the support gate")
    return ["primary.minimum_finite_reactions_per_context", "status_rule.PL2D_INSUFFICIENT_SUPPORT"]


def amended_contract(original: dict) -> tuple[dict, list[str]]:
    amended = copy.deepcopy(original)
    amended["primary"] = core.amended_primary(original["primary"])
    amended["status_rule"]["PL2D_INSUFFICIENT_SUPPORT"] = INSUFFICIENT_RULE
    amended["supersedes_support_gate_from_manifest_sha256"] = core.ORIGINAL_MANIFEST_SHA256
    amended["amendment_scope"] = "minimum_finite_reactions_per_context_only"
    return amended, contract_diff(original, amended)


def _implementation_hashes() -> dict[str, str]:
    paths = {"amendment_patch_sha256": PATCH_REL, "amendment_core_sha256": CORE_REL,
             "amendment_adapter_sha256": ADAPTER_REL, "amendment_contract_document_sha256": DOC_REL}
    hashes = {key: _sha256(ROOT / path) for key, path in paths.items()}
    if hashes["amendment_patch_sha256"] != PATCH_SHA256:
        raise RuntimeError("amendment patch hash mismatch")
    return hashes


def _payloads() -> tuple[dict[str, bytes], dict]:
    original_manifest, snapshot, original_contract, source_hashes = verify_source()
    implementation_hashes = _implementation_hashes()
    counts = core.verify_frozen_development_contexts(snapshot["algorithm_rna_contexts"])
    derivation = core.derive_amended_context_floor(counts)
    if (derivation["minimum_development_context_support"] != 129 or
            derivation["proportional_confirmation_support_numerator"] != 43 or
            derivation["proportional_confirmation_support_denominator"] != 1 or
            derivation["derived_context_floor"] != core.AMENDED_MINIMUM_FINITE_REACTIONS_PER_CONTEXT or
            core.ORIGINAL_MINIMUM_FINITE_REACTIONS_PER_CONTEXT <= 43):
        raise RuntimeError("development support calibration mismatch")
    contract, differences = amended_contract(original_contract)
    calibration = {
        "schema": core.SCHEMA,
        "development_context_support": [
            {"algorithm": algorithm, "rna_context_key": context, "finite_reactions": count}
            for algorithm, context, count in core.FROZEN_DEVELOPMENT_CONTEXT_SUPPORT
        ],
        "development_reaction_count": core.DEVELOPMENT_REACTIONS,
        "confirmation_reaction_count": core.CONFIRMATION_REACTIONS,
        **derivation,
        "retention_fraction": 0.70,
        "original_context_floor": core.ORIGINAL_MINIMUM_FINITE_REACTIONS_PER_CONTEXT,
        "amended_context_floor": core.AMENDED_MINIMUM_FINITE_REACTIONS_PER_CONTEXT,
        "confirmation_support_values_inspected": False,
        "confirmation_support_statement": "No confirmation support value was inspected.",
    }
    audit = {
        "schema": core.SCHEMA,
        "status": "PASS",
        "original_pl2d_h_manifest_sha256": source_hashes[SOURCE_MANIFEST],
        "original_pl2d_h_artifact_sha256": core.ORIGINAL_ARTIFACT_SHA256,
        **implementation_hashes,
        "parsed_original_files": [SOURCE_MANIFEST, "BRIDGEPL2DH_SOURCE_AUDIT.json", "BRIDGEPL2DH_DEVELOPMENT_SNAPSHOT.json", "BRIDGEPL2DH_CONFIRMATION_CONTRACT.json"],
        "confirmation_registry_hashed_not_parsed": True,
        "confirmation_predictor_values_accessed": False,
        "confirmation_outcomes_accessed": False,
        "confirmation_support_counts_accessed": False,
        "confirmation_geometry_utility_analysis_invoked": False,
    }
    blobs = {
        "BRIDGEPL2DHA_SOURCE_AUDIT.json": _json_bytes(audit),
        "BRIDGEPL2DHA_SUPPORT_CALIBRATION.json": _json_bytes(calibration),
        "BRIDGEPL2DHA_AMENDED_CONFIRMATION_CONTRACT.json": _json_bytes(contract),
    }
    manifest = {
        "schema": core.SCHEMA,
        "status": core.STATUS,
        "artifact_sha256": {name: hashlib.sha256(blob).hexdigest() for name, blob in blobs.items()},
        "original_pl2d_h_manifest_sha256": source_hashes[SOURCE_MANIFEST],
        "original_pl2d_h_artifact_sha256": core.ORIGINAL_ARTIFACT_SHA256,
        "implementation_sha256": implementation_hashes,
        "contract_semantic_diff": differences,
        "original_pl2d_h_preserved": True,
        "pl2d_confirmation_status": "NOT_RUN",
        "confirmation_predictor_values_accessed": False,
        "confirmation_outcomes_accessed": False,
        "confirmation_support_counts_accessed": False,
        "confirmation_geometry_utility_analysis_invoked": False,
        "original_minimum_finite_reactions_per_context": 100,
        "amended_minimum_finite_reactions_per_context": 30,
        "primary_minimum_finite_reactions": original_contract["primary"]["minimum_finite_reactions"],
        "primary_minimum_rho": original_contract["primary"]["minimum_spearman_rho"],
        "primary_minimum_positive_contexts": original_contract["primary"]["minimum_positive_contexts"],
        "primary_context_count": original_contract["primary"]["context_count"],
    }
    blobs[MANIFEST] = _json_bytes(manifest)
    if _source_hashes() != source_hashes:
        raise RuntimeError("original PL2D-H source mutated during amendment")
    return blobs, manifest


def produce(output: Path) -> dict:
    outputs = (ROOT / "outputs").resolve(strict=True)
    if output.is_symlink():
        raise RuntimeError("amendment output may not be a symlink")
    output = output.resolve()
    if not output.is_relative_to(outputs) or output == outputs or output.is_relative_to(ROOT / SOURCE_REL):
        raise RuntimeError("amendment output must be a separate in-repository outputs path")
    blobs, manifest = _payloads()
    if output.exists():
        if not output.is_dir() or {p.name for p in output.iterdir()} != set(blobs):
            raise RuntimeError("existing amendment output inventory mismatch")
        if any(_file(output / name).read_bytes() != blob for name, blob in blobs.items()):
            raise RuntimeError("existing amendment output mutated or incompatible")
        return {"status": "NO_OP_EXISTING_IDENTICAL_PL2D_H_AMENDMENT", "manifest": manifest}
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="dmi_bridge_pl2d_ha_publish_", dir=output.parent) as temporary:
        stage = Path(temporary) / output.name
        stage.mkdir()
        for name, blob in blobs.items():
            (stage / name).write_bytes(blob)
        _source_hashes()
        os.replace(stage, output)
    return {"status": core.STATUS, "manifest": manifest}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / OUTPUT_REL)
    args = parser.parse_args()
    result = produce(args.output)
    print(json.dumps({"status": result["status"], "pl2d_confirmation_status": result["manifest"]["pl2d_confirmation_status"]}, sort_keys=True))


if __name__ == "__main__":
    main()
