"""Fail-closed helpers for bounded Stage-11 full-GEM qualification."""
from __future__ import annotations

import json
import ctypes
import sys
import os
import signal
import subprocess
import time
from pathlib import Path
from typing import Any, Mapping, Sequence

import pandas as pd
import numpy as np

from .cache_identity import ensemble_identity, projection_identity
from .core import STAGE, file_sha256, stable_hash
from .expert_artifacts import expert_directory
from .jobs import stage4_reaction_universe

QUALIFICATION_PROTOCOL = "khalsa_all_mikolajewicz_validation"
QUALIFICATION_FOLD = "all"
QUALIFICATION_TUMOR = "CT2A"
QUALIFICATION_SAMPLES = ("seta1", "seta2", "seta3")
QUALIFICATION_VARIANT = "rna_only_primary"
QUALIFICATION_VECTORS = 2
FULL_GEM_QUALIFICATION_SCHEMA = 5
RIPTIDE_NATIVE_SET_BOUNDS_QUALIFICATION_POLICY = "riptide_native_set_bounds_qualification"
RIPTIDE_NATIVE_SET_BOUNDS_QUALIFICATION_POLICY_VERSION = 1


def _reaction_signature(reaction) -> dict[str, Any]:
    """Stable structural fields that RIPTiDe qualification must not alter."""
    return {
        "stoichiometry": sorted((metabolite.id, float(value)) for metabolite, value in reaction.metabolites.items()),
        "gpr": str(reaction.gene_reaction_rule),
        "boundary": bool(reaction.boundary),
    }


def _single_positive_objective(model, *, label: str) -> tuple[str, float]:
    """Resolve one biological reaction objective without trusting a display string."""
    from cobra.util.solver import linear_reaction_coefficients

    coefficients = {
        reaction.id: float(coefficient)
        for reaction, coefficient in linear_reaction_coefficients(model).items()
        if float(coefficient) != 0.0
    }
    if len(coefficients) != 1:
        raise RuntimeError(f"{label} must have exactly one nonzero reaction objective: {coefficients}")
    reaction_id, coefficient = next(iter(coefficients.items()))
    if not np.isfinite(coefficient) or coefficient <= 0.0:
        raise RuntimeError(f"{label} objective coefficient must be finite and positive: {coefficients}")
    if not np.isclose(coefficient, 1.0, rtol=0.0, atol=0.0):
        raise RuntimeError(f"{label} objective coefficient must be exactly one: {coefficients}")
    if str(model.objective.direction) != "max":
        raise RuntimeError(f"{label} objective direction must be max")
    return reaction_id, coefficient


def validate_riptide_native_set_bounds_qualification(
    model,
    conditioned_parent,
    native_samples: pd.DataFrame | None,
    state: Mapping[str, Any],
    metadata: Mapping[str, Any],
    record: Mapping[str, Any],
    config: Mapping[str, Any],
    objective: Mapping[str, Any],
    required_after_fit: Sequence[str],
) -> dict[str, Any]:
    """Validate canonical RIPTiDe ``set_bounds`` provenance without refitting.

    RIPTiDe applies its fraction to a transient post-pruning/pre-``set_bounds``
    model.  That denominator is not persisted, so the final-model ratio below
    is a necessary consistency check, never a reconstruction of that transient
    quantity.
    """
    from .cache_identity import FIT_SEMANTICS_VERSION, algorithm_fit_fingerprint

    tolerance = float(config["validation"]["bound_tolerance"])
    objective_tolerance = float(config["validation"].get("objective_tolerance", 1.0e-8))
    if str(record.get("algorithm")) != "RIPTiDe":
        raise RuntimeError("RIPTiDe set-bounds qualification invoked for another algorithm")
    try:
        spec = json.loads(str(record["expert_spec"]))
    except (KeyError, TypeError, json.JSONDecodeError) as exc:
        raise RuntimeError("RIPTiDe qualification lacks a readable authoritative expert specification") from exc
    if str(spec.get("algorithm")) != "RIPTiDe":
        raise RuntimeError("RIPTiDe expert specification has a mismatched algorithm")
    if int(spec.get("algorithm_semantics_version", -1)) != FIT_SEMANTICS_VERSION["RIPTiDe"]:
        raise RuntimeError("RIPTiDe expert specification has a stale fit semantics version")
    current_fit_fingerprint = algorithm_fit_fingerprint("RIPTiDe")
    if str(spec.get("algorithm_fit_fingerprint")) != current_fit_fingerprint:
        raise RuntimeError("RIPTiDe expert specification has a stale fit fingerprint")

    parameters = dict(spec.get("fit_parameters") or {})
    riptide_config = dict(config["contextualization"]["riptide"])
    configured_fraction = float(riptide_config["fraction"])
    persisted_fraction = float(state.get("fraction_of_optimum", np.nan))
    if not (parameters.get("set_bounds") is True and metadata.get("set_bounds") is True and state.get("set_bounds") is True):
        raise RuntimeError("RIPTiDe set-bounds provenance is absent or false")
    if list(parameters.get("tasks") or []) or list(parameters.get("exclude") or []):
        raise RuntimeError("RIPTiDe generation specification is not canonical tasks=[], exclude=[]")
    if not np.isclose(float(parameters.get("fraction", np.nan)), configured_fraction, rtol=0.0, atol=tolerance):
        raise RuntimeError("RIPTiDe generation fraction differs from configuration")
    if not np.isfinite(persisted_fraction) or not np.isclose(persisted_fraction, configured_fraction, rtol=0.0, atol=tolerance):
        raise RuntimeError("RIPTiDe persisted native fraction differs from configuration")

    role_entries = list((spec.get("requirement_policy") or {}).get("objective_reactions") or [])
    if len(role_entries) != 1:
        raise RuntimeError("RIPTiDe generation specification must declare exactly one objective role")
    configured_objective = str(config["objective"]["biomass_reaction"])
    policy_objective = str(role_entries[0].get("reaction_id", ""))
    state_objective = str(state.get("objective_reaction", ""))
    parent_objective, _ = _single_positive_objective(conditioned_parent, label="Conditioned parent")
    model_objective, _ = _single_positive_objective(model, label="RIPTiDe expert")
    identities = {configured_objective, policy_objective, state_objective, parent_objective, model_objective}
    if len(identities) != 1 or not configured_objective:
        raise RuntimeError(f"RIPTiDe objective identity is ambiguous: {sorted(identities)}")
    objective_id = configured_objective
    if objective_id not in required_after_fit:
        raise RuntimeError("RIPTiDe objective is not required after fitting")
    if objective_id not in model.reactions or objective_id not in conditioned_parent.reactions:
        raise RuntimeError("RIPTiDe objective is absent from expert or conditioned parent")
    expert_objective = model.reactions.get_by_id(objective_id)
    parent_objective_reaction = conditioned_parent.reactions.get_by_id(objective_id)
    if expert_objective.boundary or parent_objective_reaction.boundary:
        raise RuntimeError("RIPTiDe biological objective cannot be a boundary reaction")
    if _reaction_signature(expert_objective) != _reaction_signature(parent_objective_reaction):
        raise RuntimeError("RIPTiDe objective structure changed during fitting")

    required = set(map(str, required_after_fit))
    model_ids = {reaction.id for reaction in model.reactions}
    missing = sorted(required - model_ids)
    if missing:
        raise RuntimeError(f"Required reactions absent after RIPTiDe fitting: {missing}")
    if native_samples is None or native_samples.empty:
        raise RuntimeError("RIPTiDe native samples are required for set-bounds qualification")
    expected_samples = int(parameters.get("fit_native_samples", riptide_config["fit_native_samples"]))
    if len(native_samples) != expected_samples:
        raise RuntimeError(f"RIPTiDe native sample count differs: {len(native_samples)} != {expected_samples}")

    parent_ids = {reaction.id for reaction in conditioned_parent.reactions}
    retained_audit: list[dict[str, Any]] = []
    for reaction in model.reactions:
        rid = reaction.id
        if rid not in parent_ids:
            raise RuntimeError(f"Retained RIPTiDe reaction is absent from conditioned parent: {rid}")
        if rid not in native_samples.columns:
            raise RuntimeError(f"Native RIPTiDe samples lack retained reaction coordinate: {rid}")
        parent_reaction = conditioned_parent.reactions.get_by_id(rid)
        if _reaction_signature(reaction) != _reaction_signature(parent_reaction):
            raise RuntimeError(f"Retained RIPTiDe reaction structure changed: {rid}")
        lower, upper = float(reaction.lower_bound), float(reaction.upper_bound)
        parent_lower, parent_upper = float(parent_reaction.lower_bound), float(parent_reaction.upper_bound)
        values = native_samples[rid].to_numpy(dtype=float, copy=False)
        if not (np.isfinite([lower, upper, parent_lower, parent_upper]).all() and np.isfinite(values).all()):
            raise RuntimeError(f"RIPTiDe bounds/native samples are non-finite: {rid}")
        if lower > upper:
            raise RuntimeError(f"RIPTiDe bounds are reversed: {rid}")
        if lower < parent_lower - tolerance or upper > parent_upper + tolerance:
            raise RuntimeError(f"RIPTiDe bounds widen conditioned-parent interval: {rid}")
        sample_lower, sample_upper = float(np.min(values)), float(np.max(values))
        parent_unchanged = np.isclose(lower, parent_lower, rtol=0.0, atol=tolerance) and np.isclose(upper, parent_upper, rtol=0.0, atol=tolerance)
        extrema_match = np.isclose(lower, sample_lower, rtol=0.0, atol=tolerance) and np.isclose(upper, sample_upper, rtol=0.0, atol=tolerance)
        if not (parent_unchanged or extrema_match):
            raise RuntimeError(f"RIPTiDe fitted bounds lack native-extrema provenance: {rid}")
        retained_audit.append({"reaction_id": rid, "parent_bounds": [parent_lower, parent_upper], "expert_bounds": [lower, upper], "native_sample_bounds": [sample_lower, sample_upper], "parent_bounds_retained": parent_unchanged, "native_extrema_match": extrema_match, "boundary": bool(reaction.boundary)})

    if str(objective.get("status")) != "optimal":
        raise RuntimeError("RIPTiDe contextual objective solve is not optimal")
    contextual_bmax = float(objective.get("objective_value", np.nan))
    if not np.isfinite(contextual_bmax) or contextual_bmax <= objective_tolerance:
        raise RuntimeError("RIPTiDe contextual objective capacity is non-finite or nonpositive")
    lower, upper = float(expert_objective.lower_bound), float(expert_objective.upper_bound)
    samples = native_samples[objective_id].to_numpy(dtype=float, copy=False)
    sample_lower, sample_upper = float(np.min(samples)), float(np.max(samples))
    if lower < 0.0:
        raise RuntimeError("RIPTiDe objective lower bound is negative")
    if not np.isclose(upper, contextual_bmax, rtol=0.0, atol=tolerance):
        raise RuntimeError("RIPTiDe objective upper bound does not match contextual Bmax")
    required_minimum = configured_fraction * contextual_bmax
    shortfall = required_minimum - lower
    if shortfall > tolerance:
        raise RuntimeError("RIPTiDe objective lower bound violates persisted native-fraction consistency")
    return {
        "passed": True,
        "policy_name": RIPTIDE_NATIVE_SET_BOUNDS_QUALIFICATION_POLICY,
        "policy_version": RIPTIDE_NATIVE_SET_BOUNDS_QUALIFICATION_POLICY_VERSION,
        "objective_id": objective_id,
        "objective_enabled_evidence": "current RIPTiDe fit fingerprint binds adapter objective=True invocation",
        "canonical_set_bounds": True,
        "configured_fraction": configured_fraction,
        "persisted_fraction_of_optimum": persisted_fraction,
        "native_sample_count": int(len(native_samples)),
        "objective_parent_bounds": [float(parent_objective_reaction.lower_bound), float(parent_objective_reaction.upper_bound)],
        "objective_expert_bounds": [lower, upper],
        "objective_native_sample_bounds": [sample_lower, sample_upper],
        "contextual_bmax": contextual_bmax,
        "persisted_consistency_required_minimum": required_minimum,
        "objective_lower_to_contextual_bmax_ratio": lower / contextual_bmax,
        "persisted_consistency_shortfall": shortfall,
        "bound_tolerance": tolerance,
        "historical_pre_set_bounds_bmax": None,
        "historical_pre_set_bounds_bmax_interpretation": "not persisted; final-Bmax fraction is necessary consistency only",
        "required_reactions_after_fit": sorted(required),
        "retained_reaction_bound_audit": retained_audit,
    }


def full_gem_qualification_fingerprint() -> str:
    """Fingerprint the controller/worker/helper implementation that issued a qualification."""
    return stable_hash({
        "schema": FULL_GEM_QUALIFICATION_SCHEMA,
        "qualification_helpers_sha256": file_sha256(STAGE / "src/qualification.py"),
        "qualification_controller_sha256": file_sha256(STAGE / "scripts/qualify_full_gem.py"),
        "qualification_worker_sha256": file_sha256(STAGE / "scripts/run_full_gem_qualification_worker.py"),
        "jobs_orchestration_sha256": file_sha256(STAGE / "src/jobs.py"),
    })


def _optional_text(value: Any) -> str | None:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    return None if text.lower() in {"", "nan", "none"} else text


def validate_qualification_context(record: Mapping[str, Any]) -> None:
    """Fail closed on ambiguous RNA-only/RNA+DMI biological contexts."""
    required = ("algorithm", "condition_type", "tumor", "context_id", "rna_training_slice_hash",
                "protocol", "fold", "expert_variant", "expert_hash",
                "generation_spec_identity", "artifact_model_hash")
    missing = [key for key in required if not record.get(key)]
    if missing:
        raise ValueError(f"Qualification context lacks required keys: {missing}")
    tumor, condition, context, mouse = (str(record["tumor"]), str(record["condition_type"]),
                                        str(record["context_id"]), record.get("dmi_mouse_id"))
    if condition == "rna_only":
        if mouse not in (None, ""):
            raise ValueError("RNA-only qualification cannot carry a DMI mouse")
        if not context.startswith(f"{tumor}__rna_only__"):
            raise ValueError("RNA-only context does not agree with tumor")
    elif condition == "rna_dmi":
        if mouse in (None, "") or str(mouse) not in {"C1", "C2", "C3", "C4", "C5", "G1", "G2", "G3", "G4", "G5"}:
            raise ValueError("RNA+DMI qualification requires a recognized mouse")
        if not str(mouse).startswith("C" if tumor == "CT2A" else "G") or context != str(mouse):
            raise ValueError("RNA+DMI mouse context is inconsistent with tumor")
    else:
        raise ValueError(f"Unknown qualification condition type: {condition}")


def qualification_matches_reference(reference: Mapping[str, Any], qualification: Mapping[str, Any]) -> bool:
    """Bind a passing qualification to the exact current pre-fit expert reference."""
    try:
        validate_qualification_context(qualification)
    except ValueError:
        return False
    exact_keys = (
        "algorithm", "condition_type", "tumor", "context_id", "rna_training_slice_hash",
        "protocol", "fold", "expert_variant", "expert_hash", "generation_spec_identity",
    )
    if any(str(qualification.get(key)) != str(reference.get(key)) for key in exact_keys):
        return False
    return bool(
        qualification.get("qualification_status") == "passed"
        and _optional_text(qualification.get("dmi_mouse_id")) == _optional_text(reference.get("dmi_mouse_id"))
    )


def currently_qualified_expert(record: Mapping[str, Any], qualification: Mapping[str, Any]) -> bool:
    """Require a current passing qualification for the exact immutable artifact/context."""
    return bool(
        qualification_matches_reference(record, qualification)
        and qualification.get("artifact_model_hash") == record.get("artifact_model_hash")
        and qualification.get("qualification_identity") == record.get("qualification_identity")
    )


def select_qualification_record(references: pd.DataFrame, algorithm: str, *, tumor: str = QUALIFICATION_TUMOR) -> dict[str, Any]:
    """Resolve exactly one intended administrative reference; never choose by row order."""
    if tumor not in {"CT2A", "GL261"}:
        raise ValueError(f"Unsupported qualification tumor: {tumor}")
    selected = references.loc[
        references["algorithm"].eq(algorithm)
        & references["protocol"].eq(QUALIFICATION_PROTOCOL)
        & references["fold"].astype(str).eq(QUALIFICATION_FOLD)
        & references["tumor"].eq(tumor)
        & references["expert_variant"].eq(QUALIFICATION_VARIANT)
    ]
    if len(selected) != 1:
        raise RuntimeError(
            "qualification_selection_error: expected exactly one reference for "
            f"{algorithm}, found {len(selected)}"
        )
    record = selected.iloc[0].to_dict()
    if not str(record.get("expert_hash", "")):
        raise RuntimeError("qualification_selection_error: selected reference has no expert identity")
    return record


def qualification_record(algorithm: str, config: Mapping[str, Any], *, tumor: str = QUALIFICATION_TUMOR) -> tuple[dict[str, Any], dict[str, Any]]:
    references = pd.read_csv(STAGE / "manifests/expert_references.tsv", sep="\t")
    expert = select_qualification_record(references, algorithm, tumor=tumor)
    if pd.isna(expert.get("dmi_mouse_id")):
        expert["dmi_mouse_id"] = None
    # A qualification may exercise a blocked adapter, but never changes production manifests.
    expert["execution_status"] = "planned"
    ensemble_hash, seed, ensemble_spec = ensemble_identity(
        str(expert["expert_hash"]), algorithm, config, QUALIFICATION_VECTORS
    )
    ensemble = expert | {
        "ensemble_hash": ensemble_hash,
        "ensemble_seed": seed,
        "ensemble_spec": json.dumps(ensemble_spec, sort_keys=True, separators=(",", ":")),
        "requested_vectors": QUALIFICATION_VECTORS,
        "ensemble_replicate": 0,
        "sampling_thinning": int(config["sampling"]["thinning"]),
        "sampling_processes": int(config["sampling"]["processes"]),
    }
    return expert, ensemble


def projection_record(ensemble: Mapping[str, Any], paths: Mapping[str, Path]) -> tuple[dict[str, Any], list[str]]:
    reactions = stage4_reaction_universe(paths)
    digest, spec = projection_identity(str(ensemble["ensemble_hash"]), "stage4_4747", reactions)
    return dict(ensemble) | {
        "projection_name": "stage4_4747",
        "projection_hash": digest,
        "projection_spec": json.dumps(spec, sort_keys=True, separators=(",", ":")),
    }, reactions


def tree_checksums(root: Path) -> dict[str, str]:
    if not root.exists():
        return {}
    return {
        str(path.relative_to(root)): file_sha256(path)
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def qualification_smoke_artifact_identity(
    record: Mapping[str, Any],
    ensemble: Mapping[str, Any],
    paths: Mapping[str, Path],
    config: Mapping[str, Any],
) -> tuple[str, dict[str, Any]]:
    """Hash the exact two-vector smoke artifacts without any production job."""
    from .cache_identity import sampling_semantics_fingerprint
    from .ensembles import ensemble_directory, valid_ensemble
    from .projections import projection_directory, valid_projection

    if int(ensemble.get("requested_vectors", -1)) != QUALIFICATION_VECTORS:
        raise RuntimeError("Qualification smoke identity requires exactly two vectors")
    projection, _ = projection_record(ensemble, paths)
    if not valid_ensemble(ensemble) or not valid_projection(projection):
        raise RuntimeError("Qualification smoke artifacts are absent, stale, or invalid")
    ensemble_root = ensemble_directory(ensemble)
    projection_root = projection_directory(projection)
    spec = {
        "scope": "qualification_only_two_vector_smoke",
        "production": False,
        "algorithm": str(record["algorithm"]),
        "tumor": str(record["tumor"]),
        "expert_hash": str(record["expert_hash"]),
        "generation_spec_identity": str(record["generation_spec_identity"]),
        "sampling_semantics_fingerprint": sampling_semantics_fingerprint(
            str(record["algorithm"]), config
        ),
        "ensemble_hash": str(ensemble["ensemble_hash"]),
        "ensemble_files": tree_checksums(ensemble_root),
        "projection_hash": str(projection["projection_hash"]),
        "projection_files": tree_checksums(projection_root),
        "requested_vectors": QUALIFICATION_VECTORS,
    }
    return stable_hash(spec), spec


def terminate_process_group(process: subprocess.Popen, grace_seconds: float = 20.0) -> bool:
    """Terminate one qualification group and return whether any member survived."""
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        process.communicate(timeout=grace_seconds)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.communicate()
    deadline = time.monotonic() + grace_seconds
    while time.monotonic() < deadline:
        try:
            while True:
                waited, _ = os.waitpid(-process.pid, os.WNOHANG)
                if waited == 0:
                    break
        except ChildProcessError:
            pass
        try:
            os.killpg(process.pid, 0)
        except ProcessLookupError:
            return False
        time.sleep(0.05)
    return True


def run_bounded_phase(command: Sequence[str], cwd: Path, timeout_seconds: int) -> dict[str, Any]:
    if sys.platform.startswith("linux"):
        ctypes.CDLL(None).prctl(36, 1, 0, 0, 0)  # PR_SET_CHILD_SUBREAPER
    started = time.monotonic()
    process = subprocess.Popen(
        list(command), cwd=cwd, start_new_session=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    try:
        output, _ = process.communicate(timeout=timeout_seconds)
        survivors = False
        try:
            os.killpg(process.pid, 0)
            survivors = True
        except ProcessLookupError:
            pass
        result = {
            "status": "passed" if process.returncode == 0 else "qualification_failed_validation",
            "returncode": process.returncode,
            "elapsed_seconds": time.monotonic() - started,
            "timed_out": False,
            "process_group": process.pid,
            "workers_survived": survivors,
            "output": output[-12000:],
        }
        try:
            payload = json.loads(output.rstrip().splitlines()[-1])
            if isinstance(payload, dict) and isinstance(payload.get("details"), dict):
                result["details"] = payload["details"]
        except (IndexError, json.JSONDecodeError):
            pass
        return result
    except subprocess.TimeoutExpired:
        survivors = terminate_process_group(process)
        return {
            "status": "qualification_failed_timeout",
            "returncode": process.returncode,
            "elapsed_seconds": time.monotonic() - started,
            "timed_out": True,
            "process_group": process.pid,
            "workers_survived": survivors,
            "output": "",
        }


def expert_model_checksum(record: Mapping[str, Any]) -> str:
    return file_sha256(expert_directory(record) / "expert_model.json.xz")
