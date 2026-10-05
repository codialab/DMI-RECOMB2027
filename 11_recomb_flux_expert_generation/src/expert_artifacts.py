"""Portable, checksum-protected, atomically promoted expert artifacts."""
from __future__ import annotations
import json
import lzma
import os
import shutil
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence
import cobra
import numpy as np
import pandas as pd
from .cache_identity import EXPERT_CACHE_SCHEMA
from .core import (STAGE, atomic_json, atomic_table, file_sha256, parquet_exists,
                   reaction_universe_hash, required_artifacts_present, stored_checksum_matches)
from .models import enforce_gurobi, model_hash


def expert_directory(record: Mapping[str, Any]) -> Path:
    return STAGE / "expert_artifacts" / str(record["algorithm"]) / str(record["expert_hash"])


def expert_relative_locator(record: Mapping[str, Any]) -> str:
    return str(Path("expert_artifacts") / str(record["algorithm"]) / str(record["expert_hash"]))


def _atomic_model(path: Path, model) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("wb", dir=path.parent, delete=False) as handle: temporary = Path(handle.name)
    try:
        with lzma.open(temporary, "wt") as handle: handle.write(cobra.io.to_json(model))
        os.replace(temporary, path)
    finally: temporary.unlink(missing_ok=True)


def load_expert_model(path: Path, *, seed: int | None = None,
                      activate_solver: bool = True):
    """Load a persisted model, activating Gurobi only for solve-capable callers."""
    configuration = cobra.Configuration()
    original_solver = configuration.solver
    try:
        # COBRApy constructs a solver while deserializing.  Cache-integrity
        # validation is structural and must remain usable when a transient
        # Gurobi license is unavailable.
        configuration.solver = "glpk"
        with lzma.open(path, "rt") as handle:
            model = cobra.io.from_json(handle.read())
    finally:
        configuration.solver = original_solver
    if activate_solver:
        enforce_gurobi(model, seed=seed)
    return model


def marker_valid(directory: Path, hashes: Mapping[str, str], required: set[str], schema: int) -> bool:
    marker = directory / "DONE.json"
    if not marker.is_file(): return False
    try: done = json.loads(marker.read_text())
    except Exception: return False
    if done.get("cache_schema") != schema or not done.get("qc_passed") or any(done.get(k) != v for k, v in hashes.items()): return False
    checksums = done.get("checksums", {})
    return required_artifacts_present(required, checksums) and all(
        stored_checksum_matches(directory, name, checksum) for name, checksum in checksums.items()
    )


def valid_expert(record: Mapping[str, Any]) -> bool:
    directory = expert_directory(record); required = {"expert_model.json.xz", "metadata.json", "algorithm_state.json", "qc/dmi_flux_capability.tsv", "qc/model_qualification.json"}
    if not marker_valid(directory, {"expert_hash": str(record["expert_hash"])}, required, EXPERT_CACHE_SCHEMA): return False
    try:
        metadata = json.loads((directory / "metadata.json").read_text()); model = load_expert_model(directory / "expert_model.json.xz", activate_solver=False)
        qualification = json.loads((directory / "qc/model_qualification.json").read_text())
        model_ids = [reaction.id for reaction in model.reactions]; parent_ids = list(map(str, metadata["ordered_parent_reactions"])); removed = list(map(str, metadata["contextualization_removed_reactions"]))
        base_valid = (metadata.get("cache_schema") == EXPERT_CACHE_SCHEMA and metadata.get("expert_hash") == record["expert_hash"]
            and metadata.get("artifact_relative_locator") == expert_relative_locator(record)
            and metadata.get("model_hash") == model_hash(model) and metadata.get("ordered_expert_reactions") == model_ids
            and len(parent_ids) == len(set(parent_ids)) and set(model_ids) <= set(parent_ids)
            and set(removed) == set(parent_ids) - set(model_ids)
            and metadata.get("parent_reaction_universe_sha256") == reaction_universe_hash(parent_ids)
            and qualification.get("passed") is True
            and qualification.get("artifact_model_hash") == metadata.get("model_hash")
            and qualification.get("expert_hash") == record["expert_hash"]
            and qualification.get("qualification_identity") == metadata.get("qualification_identity")
            and qualification.get("parent_model_unchanged") is True
            and qualification.get("reload_equivalence", {}).get("passed") is True
            and (record.get("generation_spec_identity") is None or qualification.get("generation_spec_identity") == record.get("generation_spec_identity")))
        if not base_valid:
            return False
        if str(record.get("algorithm")) == "GIMME" and metadata.get("mapping_audit_sha256"):
            relative = metadata.get("gimme_mapping_audit_relative_path")
            audit = directory / str(relative) if relative else None
            return bool(
                audit is not None and audit.is_file()
                and metadata.get("gimme_mapping_audit_file_sha256") == file_sha256(audit)
                and qualification.get("gimme_mapping_audit_sha256") == metadata.get("mapping_audit_sha256")
            )
        return True
    except Exception: return False


def valid_generated_expert(record: Mapping[str, Any]) -> bool:
    """Verify an immutable generated artifact independently of qualification rules."""
    if not valid_expert(record):
        return False
    try:
        metadata = json.loads((expert_directory(record) / "metadata.json").read_text())
        generation = record.get("generation_spec_identity")
        return bool(
            metadata.get("artifact_model_hash") == metadata.get("model_hash")
            and (generation is None or metadata.get("generation_spec_identity") == generation)
        )
    except Exception:
        return False


_QUALIFICATION_METADATA_FILES = {
    "DONE.json",
    "metadata.json",
    "qc/model_qualification.json",
}


def immutable_generated_checksums(directory: Path) -> dict[str, str]:
    """Hash every generated artifact that qualification refresh must preserve."""
    return {
        str(path.relative_to(directory)): file_sha256(path)
        for path in sorted(directory.rglob("*"))
        if path.is_file()
        and str(path.relative_to(directory)) not in _QUALIFICATION_METADATA_FILES
    }


def refresh_expert_qualification(
    record: Mapping[str, Any], model_qualification: Mapping[str, Any]
) -> dict[str, Any]:
    """Atomically refresh qualification metadata without regenerating the expert.

    The existing checksum-valid artifact is copied byte-for-byte.  Only the
    qualification report, its pointer in metadata, and the completion-marker
    checksums are rewritten.  Every other file must retain its SHA-256 digest.
    """
    final = expert_directory(record)
    if not valid_expert(record):
        raise RuntimeError(
            "Qualification refresh requires a checksum-valid persisted expert"
        )
    metadata_path = final / "metadata.json"
    qualification_path = final / "qc/model_qualification.json"
    metadata = json.loads(metadata_path.read_text())
    previous_qualification = json.loads(qualification_path.read_text())
    generation_identity = record.get("generation_spec_identity")
    if generation_identity is None or str(metadata.get("generation_spec_identity")) != str(
        generation_identity
    ):
        raise RuntimeError(
            "Qualification refresh requires the exact current generation identity"
        )
    artifact_hash = str(metadata.get("model_hash", ""))
    if not artifact_hash or str(metadata.get("artifact_model_hash")) != artifact_hash:
        raise RuntimeError("Persisted expert semantic model hashes are inconsistent")
    required = {
        "passed": True,
        "expert_hash": str(record["expert_hash"]),
        "generation_spec_identity": str(generation_identity),
        "artifact_model_hash": artifact_hash,
    }
    mismatched = [
        key
        for key, expected in required.items()
        if model_qualification.get(key) != expected
    ]
    if mismatched:
        raise RuntimeError(
            "Replacement qualification is not bound to the current generated expert: "
            f"{mismatched}"
        )
    qualification_identity = str(model_qualification.get("qualification_identity", ""))
    if not qualification_identity:
        raise RuntimeError("Replacement qualification identity is absent")

    immutable_before = immutable_generated_checksums(final)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{str(record['expert_hash'])[:16]}.requalify.", dir=final.parent)
    )
    try:
        shutil.copytree(final, temporary, dirs_exist_ok=True, copy_function=shutil.copy2)
        refreshed_metadata = dict(metadata)
        refreshed_metadata["qualification_identity"] = qualification_identity
        atomic_json(
            temporary / "qc/model_qualification.json", _json_safe(model_qualification)
        )
        atomic_json(temporary / "metadata.json", refreshed_metadata)

        marker = json.loads((temporary / "DONE.json").read_text())
        marker["checksums"] = {
            str(path.relative_to(temporary)): file_sha256(path)
            for path in sorted(temporary.rglob("*"))
            if path.is_file() and path.name != "DONE.json"
        }
        atomic_json(temporary / "DONE.json", marker)
        immutable_after = immutable_generated_checksums(temporary)
        if immutable_after != immutable_before:
            raise RuntimeError(
                "Qualification refresh changed generated scientific artifact bytes"
            )
        if previous_qualification.get("artifact_model_hash") != artifact_hash:
            raise RuntimeError(
                "Existing qualification was not bound to the persisted semantic model"
            )

        _archive_existing(final, "stale_experts")
        os.replace(temporary, final)
        if not valid_expert(record):
            raise RuntimeError("Refreshed expert artifact failed checksum validation")
        immutable_promoted = immutable_generated_checksums(final)
        if immutable_promoted != immutable_before:
            raise RuntimeError(
                "Promoted qualification refresh changed generated scientific artifacts"
            )
        return {
            "status": "requalified_no_refit",
            "previous_qualification_identity": previous_qualification.get(
                "qualification_identity"
            ),
            "qualification_identity": qualification_identity,
            "immutable_generated_checksums_before": immutable_before,
            "immutable_generated_checksums_after": immutable_promoted,
        }
    except Exception:
        if temporary.exists():
            _archive_existing(temporary, "failed_experts")
        raise


def _json_safe(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)): return value
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, Mapping): return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)): return [_json_safe(v) for v in value]
    return repr(value)


def _archive_existing(path: Path, category: str) -> None:
    if path.exists():
        destination = STAGE / "qc" / category / f"{path.name}.{time.time_ns()}"; destination.parent.mkdir(parents=True, exist_ok=True); os.replace(path, destination)


def persist_expert(record: Mapping[str, Any], result, provenance: Mapping[str, Any], parent_reactions: Sequence[str] | None = None,
                   dmi_capability: pd.DataFrame | None = None,
                   model_qualification: Mapping[str, Any] | None = None) -> Path:
    final = expert_directory(record); final.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{record['expert_hash'][:16]}.", dir=final.parent))
    try:
        parent = list(map(str, parent_reactions)) if parent_reactions is not None else [reaction.id for reaction in result.model.reactions]; expert = [reaction.id for reaction in result.model.reactions]
        if dmi_capability is None: dmi_capability = pd.DataFrame(columns=["reaction_id"])
        if model_qualification is None or model_qualification.get("passed") is not True:
            raise RuntimeError("A passing identity-bound model qualification is required")
        if len(parent) != len(set(parent)) or not set(expert) <= set(parent): raise ValueError("Expert reaction IDs do not form a unique parent subset")
        removed = [rid for rid in parent if rid not in set(expert)]
        _atomic_model(temporary / "expert_model.json.xz", result.model)
        persisted_hash = model_hash(result.model)
        reloaded = load_expert_model(temporary / "expert_model.json.xz")
        reloaded_hash = model_hash(reloaded)
        reload_equivalence = {
            "passed": bool(
                reloaded_hash == persisted_hash
                and [reaction.id for reaction in reloaded.reactions] == expert
                and all(
                    reloaded.reactions.get_by_id(rid).bounds == result.model.reactions.get_by_id(rid).bounds
                    for rid in expert
                )
            ),
            "persisted_model_hash": persisted_hash,
            "reloaded_model_hash": reloaded_hash,
            "reaction_order_equal": [reaction.id for reaction in reloaded.reactions] == expert,
        }
        if not reload_equivalence["passed"]:
            raise RuntimeError("Serialized/reloaded expert model is not semantically equivalent")
        algorithm_state = dict(result.algorithm_state)
        mapping_rows = algorithm_state.pop("_dedicated_gimme_mapping_audit", None)
        if mapping_rows is not None:
            audit_relative = "qc/gimme_irreversible_mapping.tsv"
            audit_path = temporary / audit_relative
            atomic_table(audit_path, pd.DataFrame(_json_safe(mapping_rows)))
            audit_file_sha256 = file_sha256(audit_path)
            audit_content_sha256 = str(algorithm_state.get("mapping_audit_sha256", ""))
            if not audit_content_sha256:
                raise RuntimeError("GIMME mapping audit lacks a stable content identity")
            algorithm_state.update({
                "gimme_mapping_audit_relative_path": audit_relative,
                "gimme_mapping_audit_file_sha256": audit_file_sha256,
            })
            result.metadata.update({
                "gimme_mapping_audit_relative_path": audit_relative,
                "gimme_mapping_audit_file_sha256": audit_file_sha256,
            })
            model_qualification = dict(model_qualification) | {
                "gimme_mapping_audit_sha256": audit_content_sha256,
                "gimme_mapping_audit_relative_path": audit_relative,
                "gimme_mapping_audit_file_sha256": audit_file_sha256,
            }
        qualification = _json_safe(dict(model_qualification) | {"reload_equivalence": reload_equivalence})
        if qualification.get("artifact_model_hash") != persisted_hash:
            raise RuntimeError("Model qualification is not bound to the persisted model hash")
        atomic_json(temporary / "algorithm_state.json", _json_safe(algorithm_state))
        if result.native_samples is not None:
            frame = result.native_samples.astype(np.float32).copy(); frame.insert(0, "sample_index", np.arange(len(frame), dtype=np.int64)); atomic_table(temporary / "native_samples.parquet", frame)
        atomic_table(temporary / "qc/dmi_flux_capability.tsv", dmi_capability)
        atomic_json(temporary / "qc/model_qualification.json", qualification)
        metadata = dict(record) | dict(result.metadata) | dict(provenance) | {"cache_schema": EXPERT_CACHE_SCHEMA,
            "artifact_kind": result.artifact_kind, "artifact_relative_locator": expert_relative_locator(record), "model_hash": persisted_hash,
            "artifact_model_hash": persisted_hash, "qualification_identity": qualification["qualification_identity"],
            "native_samples_persisted": result.native_samples is not None, "ordered_parent_reactions": parent,
            "ordered_expert_reactions": expert, "contextualization_removed_reactions": removed,
            "parent_reaction_universe_sha256": reaction_universe_hash(parent), "created_utc": datetime.now(timezone.utc).isoformat()}
        metadata.pop("expert_artifact", None); atomic_json(temporary / "metadata.json", metadata)
        checksums = {str(p.relative_to(temporary)): file_sha256(p) for p in temporary.rglob("*") if p.is_file()}
        atomic_json(temporary / "DONE.json", {"cache_schema": EXPERT_CACHE_SCHEMA, "expert_hash": record["expert_hash"], "expert_id": record["expert_id"], "qc_passed": True, "checksums": checksums})
        _archive_existing(final, "stale_experts"); os.replace(temporary, final); return final
    except Exception:
        _archive_existing(temporary, "failed_experts"); raise


def load_expert_state(record: Mapping[str, Any], *, seed: int | None = None) -> tuple[Any, dict[str, Any], pd.DataFrame | None, dict[str, Any]]:
    directory = expert_directory(record)
    if not valid_expert(record): raise RuntimeError(f"Expert cache is absent, corrupt, legacy, or stale: {directory}")
    model = load_expert_model(directory / "expert_model.json.xz", seed=seed); state = json.loads((directory / "algorithm_state.json").read_text()); metadata = json.loads((directory / "metadata.json").read_text())
    native_path = directory / "native_samples.parquet"
    native = pd.read_parquet(native_path).drop(columns="sample_index").astype(np.float64) if parquet_exists(native_path) else None
    return model, state, native, metadata
