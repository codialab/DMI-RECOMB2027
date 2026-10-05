"""Full-parent Stage-11 ensemble generation, strict mapping, and QC."""
from __future__ import annotations
import json
import os
import tempfile
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence
import numpy as np
import pandas as pd
import scipy.linalg
from cobra.sampling import OptGPSampler, hr_sampler
from .cache_identity import ENSEMBLE_CACHE_SCHEMA, ensemble_implementation_fingerprint
from .core import STAGE, atomic_json, atomic_table, file_sha256, reaction_universe_hash, validate_samples
from .expert_artifacts import expert_directory, load_expert_state, marker_valid
from .models import apply_biomass_floor, model_hash
from .qc import flux_capability_table
from .riptide_native import native_resample

_COBRA_NULLSPACE = hr_sampler.nullspace
_SVD_FALLBACK_USED = False
_NULLSPACE_ATOL = 1.0e-13
_NULLSPACE_RTOL = 0.0


def _robust_nullspace(matrix: np.ndarray, atol: float = 1e-13, rtol: float = 0.0) -> np.ndarray:
    global _SVD_FALLBACK_USED
    use_atol, use_rtol = max(float(atol), _NULLSPACE_ATOL), max(float(rtol), _NULLSPACE_RTOL)
    try: return _COBRA_NULLSPACE(matrix, atol=use_atol, rtol=use_rtol)
    except np.linalg.LinAlgError:
        _, singular_values, vh = scipy.linalg.svd(np.atleast_2d(matrix), full_matrices=True, lapack_driver="gesvd")
        threshold = max(use_atol, use_rtol * singular_values[0]); _SVD_FALLBACK_USED = True
        return vh[(singular_values >= threshold).sum():].conj().T
hr_sampler.nullspace = _robust_nullspace


def ensemble_directory(record: Mapping[str, Any]) -> Path:
    return STAGE / "flux_ensembles" / str(record["algorithm"]) / str(record["ensemble_hash"])


def valid_ensemble(record: Mapping[str, Any]) -> bool:
    directory = ensemble_directory(record); required = {"full_parent_fluxes/reaction_samples.parquet", "metadata.json", "qc/qc.json", "qc/dmi_flux_capability_generation.tsv"}
    try:
        spec = json.loads(str(record["ensemble_spec"]))
        current = ensemble_implementation_fingerprint(str(record["algorithm"]))
        if spec.get("sampling_implementation_fingerprint") != current:
            return False
    except Exception:
        return False
    if not marker_valid(directory, {"expert_hash": str(record["expert_hash"]), "ensemble_hash": str(record["ensemble_hash"])}, required, ENSEMBLE_CACHE_SCHEMA): return False
    try:
        metadata = json.loads((directory / "metadata.json").read_text()); stored = pd.read_parquet(directory / "full_parent_fluxes/reaction_samples.parquet")
        parent = list(map(str, metadata["ordered_parent_reactions"])); requested = int(metadata["requested_vectors"])
        return (metadata.get("cache_schema") == ENSEMBLE_CACHE_SCHEMA and metadata.get("expert_hash") == record["expert_hash"]
            and metadata.get("ensemble_hash") == record["ensemble_hash"] and len(stored) == requested
            and list(stored.columns) == ["sample_index", *parent] and stored["sample_index"].tolist() == list(range(requested))
            and stored.iloc[:, 1:].dtypes.map(lambda dtype: np.dtype(dtype) == np.dtype("float32")).all()
            and metadata.get("full_parent_reaction_universe_sha256") == reaction_universe_hash(parent))
    except Exception: return False


def map_full_parent(samples: pd.DataFrame, parent_reactions: Sequence[str], *, contextualization_removed_reactions: Sequence[str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    parent = list(map(str, parent_reactions)); present = list(map(str, samples.columns)); removed = list(map(str, contextualization_removed_reactions))
    if len(parent) != len(set(parent)): raise ValueError("Parent reaction universe contains duplicates")
    if len(present) != len(set(present)): raise ValueError("Expert samples contain duplicate reaction IDs")
    if len(removed) != len(set(removed)): raise ValueError("Contextualization-removed reaction list contains duplicates")
    extra = set(present) - set(parent)
    if extra: raise ValueError(f"Expert output contains non-parent reaction IDs: {sorted(extra)[:5]}")
    expected_removed = set(parent) - set(present)
    if set(removed) != expected_removed: raise ValueError("Documented contextualization-removed reactions differ from parent minus sampled expert reactions")
    if set(present) & set(removed): raise ValueError("A sampled reaction is also documented as contextualization-removed")
    full = pd.DataFrame(0.0, index=samples.index, columns=parent, dtype=np.float64)
    full.loc[:, present] = samples.loc[:, present].to_numpy(dtype=np.float64)
    return full, {"parent_reaction_count": len(parent), "expert_reaction_count": len(present),
        "contextualization_removed_reactions": removed, "contextualization_removed_count": len(removed),
        "parent_reaction_universe_sha256": reaction_universe_hash(parent)}


def _common_samples(model, record: Mapping[str, Any], config: Mapping[str, Any]) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame]:
    global _SVD_FALLBACK_USED, _NULLSPACE_ATOL, _NULLSPACE_RTOL
    geometry = config["sampling"]["geometry"]; _NULLSPACE_ATOL = float(geometry["nullspace_atol"]); _NULLSPACE_RTOL = float(geometry["nullspace_rtol"])
    before = model_hash(model); floor = apply_biomass_floor(
        model, config["objective"]["biomass_reaction"], float(config["objective"]["ensemble_objective_fraction"]),
        tolerance=float(config["validation"].get("objective_tolerance", 1.0e-8)),
        direction=str(config["objective"].get("qualification_direction", "max")),
    )
    removed = record["contextualization_removed_reactions"]; parent = record["ordered_parent_reactions"]
    capability = flux_capability_table(
        model, parent, removed, float(config["validation"]["zero_tolerance"]),
        biological_objective_id=config["objective"]["biomass_reaction"],
        physiological_fraction=float(config["objective"]["ensemble_objective_fraction"]),
        objective_tolerance=float(config["validation"].get("objective_tolerance", 1.0e-8)),
        objective_direction=str(config["objective"].get("qualification_direction", "max")),
    )
    model.objective = model.problem.Objective(0, direction="min"); _SVD_FALLBACK_USED = False
    sampler = OptGPSampler(model, thinning=int(record["sampling_thinning"]), processes=int(record["sampling_processes"]), seed=int(record["ensemble_seed"]))
    samples = sampler.sample(int(record["requested_vectors"])).astype(np.float64).iloc[:int(record["requested_vectors"])].copy()
    if len(samples) != int(record["requested_vectors"]): raise RuntimeError("OptGP did not retain exactly requested vectors")
    qc = validate_samples(samples, model, config["validation"])
    return samples, {"sampler": "COBRApy OptGPSampler", "pre_sampling_model_hash": before, "state_model_hash": model_hash(model),
        "objective_floor": floor, "nullspace_svd_fallback": _SVD_FALLBACK_USED, "qc": qc}, capability


def _riptide_samples(model, state: Mapping[str, Any], record: Mapping[str, Any], config: Mapping[str, Any]) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame]:
    """Generate native RIPTiDe samples under the same paper-wide biomass floor.

    RIPTiDe's persisted ``fraction_of_optimum`` is part of its native sampling
    state, but the RECOMB comparison additionally declares a common biomass
    floor for every algorithm.  Apply that floor to the fresh in-memory expert
    copy before entering the private/native GapSplit adapter so RIPTiDe is not
    sampled from a looser feasible space than the common-OptGP experts.
    """
    before = model_hash(model)
    floor = apply_biomass_floor(
        model,
        config["objective"]["biomass_reaction"],
        float(config["objective"]["ensemble_objective_fraction"]),
        tolerance=float(config["validation"].get("objective_tolerance", 1.0e-8)),
        direction=str(config["objective"].get("qualification_direction", "max")),
    )
    removed = record["contextualization_removed_reactions"]
    parent = record["ordered_parent_reactions"]
    capability = flux_capability_table(
        model, parent, removed, float(config["validation"]["zero_tolerance"]),
        biological_objective_id=config["objective"]["biomass_reaction"],
        physiological_fraction=float(config["objective"]["ensemble_objective_fraction"]),
        objective_tolerance=float(config["validation"].get("objective_tolerance", 1.0e-8)),
        objective_direction=str(config["objective"].get("qualification_direction", "max")),
    )
    expected_native_fraction = float(config["contextualization"]["riptide"]["fraction"])
    native_fraction = float(state["fraction_of_optimum"])
    fraction_tolerance = float(config["validation"]["bound_tolerance"])
    if not np.isclose(native_fraction, expected_native_fraction, rtol=0.0, atol=fraction_tolerance):
        raise RuntimeError(
            "Persisted RIPTiDe fraction_of_optimum differs from the configured value: "
            f"{native_fraction} != {expected_native_fraction}"
        )
    samples = native_resample(
        model,
        state,
        int(record["requested_vectors"]),
        int(record["ensemble_seed"]),
    )
    qc = validate_samples(samples, model, config["validation"])
    return samples, {
        "sampler": "RIPTiDe native GapSplit",
        "pre_sampling_model_hash": before,
        "state_model_hash": model_hash(model),
        "objective_floor": floor,
        "native_fraction_of_optimum": native_fraction,
        "configured_native_fraction_of_optimum": expected_native_fraction,
        "qc": qc,
    }, capability


def validate_persisted_biomass_floor(biomass_values: Sequence[float], biomass_floor: float, configured_bound_tolerance: float) -> dict[str, Any]:
    """Validate persisted float32 biomass coordinates against the model floor."""
    values = np.asarray(biomass_values, dtype=np.float32)
    if not values.size or not np.isfinite(values).all():
        raise RuntimeError("Persisted biomass values are empty or non-finite")
    floor = float(biomass_floor)
    spacing = abs(float(np.spacing(np.float32(floor))))
    tolerance = float(configured_bound_tolerance) + spacing
    minimum = float(np.min(values))
    result = {"biomass_floor": floor, "persisted_minimum_biomass": minimum, "configured_bound_tolerance": float(configured_bound_tolerance), "float32_floor_spacing": spacing, "numerical_tolerance": tolerance, "passed": bool(minimum >= floor - tolerance)}
    if not result["passed"]:
        raise RuntimeError(f"Persisted float32 biomass minimum {minimum} is below floor {floor} minus tolerance {tolerance}")
    return result


def _archive_existing(path: Path, category: str) -> None:
    if path.exists():
        target = STAGE / "qc" / category / f"{path.name}.{time.time_ns()}"; target.parent.mkdir(parents=True, exist_ok=True); os.replace(path, target)


def generate_flux_ensemble(record: Mapping[str, Any], parent_reactions: Sequence[str], config: Mapping[str, Any]) -> Path:
    if valid_ensemble(record): return ensemble_directory(record)
    final = ensemble_directory(record); final.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{record['ensemble_hash'][:16]}.", dir=final.parent)); started = time.time(); phase = "load_expert"
    try:
        expert_path = expert_directory(record); expert_file = expert_path / "expert_model.json.xz"; checksum_before = file_sha256(expert_file)
        model, state, native, metadata = load_expert_state(record, seed=int(record["ensemble_seed"])); model_ids = [reaction.id for reaction in model.reactions]
        if list(metadata["ordered_parent_reactions"]) != list(map(str, parent_reactions)): raise RuntimeError("Requested parent universe differs from persisted expert provenance")
        record = dict(record) | {"ordered_parent_reactions": metadata["ordered_parent_reactions"], "contextualization_removed_reactions": metadata["contextualization_removed_reactions"]}
        phase = "sampling_geometry_and_vectors"
        if metadata["artifact_kind"] in {"context_specific_model", "rna_bounded_model"}: samples, sampling, capability = _common_samples(model, record, config)
        elif metadata["artifact_kind"] == "riptide_native_expert":
            samples, sampling, capability = _riptide_samples(model, state, record, config)
        else: raise RuntimeError(f"No valid ensemble generator for {metadata['artifact_kind']}")
        if list(samples.columns) != model_ids: raise RuntimeError("Sampled columns do not exactly match the persisted expert model reaction order")
        if file_sha256(expert_file) != checksum_before: raise RuntimeError("Flux generation mutated persisted expert artifact")
        phase = "full_parent_mapping"; full, mapping = map_full_parent(samples, parent_reactions, contextualization_removed_reactions=metadata["contextualization_removed_reactions"])
        storage = full.astype(np.float32)
        objective_id = str(config["objective"]["biomass_reaction"])
        if objective_id not in storage.columns:
            raise RuntimeError(f"Biomass reaction absent from full-parent ensemble: {objective_id}")
        floor = sampling.get("objective_floor", {})
        floor["persisted_float32_validation"] = validate_persisted_biomass_floor(
            storage[objective_id].to_numpy(copy=False),
            float(floor["biomass_floor"]),
            float(config["validation"]["bound_tolerance"]),
        )
        sampling["objective_floor"] = floor
        storage.insert(0, "sample_index", np.arange(len(storage), dtype=np.int64)); atomic_table(temporary / "full_parent_fluxes/reaction_samples.parquet", storage)
        atomic_table(temporary / "qc/dmi_flux_capability_generation.tsv", capability)
        atomic_json(temporary / "qc/qc.json", {"qc_passed": True, "sampling": sampling, "full_parent_mapping": mapping, "calculation_dtype": "float64", "storage_dtype": "float32"})
        outmeta = dict(record) | {"cache_schema": ENSEMBLE_CACHE_SCHEMA, "source_expert_model_hash": metadata["model_hash"], "created_utc": datetime.now(timezone.utc).isoformat(),
            "runtime_seconds": time.time()-started, "full_parent_reaction_universe_sha256": mapping["parent_reaction_universe_sha256"],
            "full_parent_reaction_count": len(parent_reactions), "ordered_parent_reactions": list(parent_reactions), "artifact_kind": metadata["artifact_kind"]}
        atomic_json(temporary / "metadata.json", outmeta); checksums = {str(p.relative_to(temporary)): file_sha256(p) for p in temporary.rglob("*") if p.is_file()}
        atomic_json(temporary / "DONE.json", {"cache_schema": ENSEMBLE_CACHE_SCHEMA, "expert_hash": record["expert_hash"], "ensemble_hash": record["ensemble_hash"], "qc_passed": True, "checksums": checksums})
        _archive_existing(final, "stale_ensembles"); os.replace(temporary, final); return final
    except Exception as exc:
        atomic_json(temporary / "failure.json", {"phase": phase, "elapsed_seconds": time.time()-started, "error": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc()})
        _archive_existing(temporary, "failed_ensembles"); raise
