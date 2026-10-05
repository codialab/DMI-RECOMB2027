"""Stage-11 manifest construction and schema-2 execution orchestration."""
from __future__ import annotations
import json
import math
import numbers
from pathlib import Path
from typing import Any, Mapping
import pandas as pd
import yaml
from .adapters import REGISTRY, effective_status, four_method_barrier_validated
from .cache_identity import (FIT_SEMANTICS_VERSION, algorithm_fit_fingerprint, expert_identity,
    ensemble_identity, projection_identity, qualification_identity)
from .core import STAGE, atomic_json, atomic_table, file_sha256, load_configuration, stable_hash
from .expert_artifacts import expert_directory, persist_expert, valid_expert
from .expert_fitting import fit_expert
from .models import (ANCHOR_ABLATION_VARIANT, PRIMARY_EXPERT_VARIANT, apply_frozen_medium,
    algorithm_required_reactions, closed_source_negative_control, load_expression, load_model,
    model_hash, postfit_required_reactions, protected_reactions_for_variant,
    qualify_context_model, qualify_riptide_context_model, requirement_policy_for_variant, resolved_frozen_medium_routes,
    validate_environmental_protection, validate_model)
from .qc import flux_capability_table


PRODUCTION_AUTHORIZATION_SCHEMA = 2
PRODUCTION_AUTHORIZATION_FILE = "PRODUCTION_AUTHORIZATION.json"
PRODUCTION_MANIFEST_FILES = (
    "expert_references.tsv",
    "expert_jobs.tsv",
    "ensemble_references.tsv",
    "ensemble_jobs.tsv",
    "projection_jobs.tsv",
)
_CANONICAL_FOLD_COLUMNS = ("protocol", "fold", "role", "sample", "tumor", "cohort")


def _production_authorization_path() -> Path:
    return STAGE / "manifests" / PRODUCTION_AUTHORIZATION_FILE


def _current_qualification_implementation_fingerprint() -> str:
    from .qualification import full_gem_qualification_fingerprint
    return full_gem_qualification_fingerprint()


def _authorization_code_hashes() -> dict[str, str]:
    source_dir = Path(__file__).resolve().parent
    return {
        "jobs_sha256": file_sha256(Path(__file__).resolve()),
        "adapters_sha256": file_sha256(source_dir / "adapters.py"),
        "parallel_runner_sha256": file_sha256(source_dir / "parallel_runner.py"),
    }


def _production_manifest_hashes() -> dict[str, str]:
    manifest_dir = STAGE / "manifests"
    return {
        name: file_sha256(manifest_dir / name)
        for name in PRODUCTION_MANIFEST_FILES
    }


def _hash_required_regular_file(path: Path, *, label: str) -> str:
    """Hash one configured regular file, refusing ambiguous authorization inputs."""
    candidate = Path(path)
    if not candidate.is_file():
        raise RuntimeError(f"Stage-11 authorization input is missing or not a regular file: {label}={candidate}")
    try:
        return file_sha256(candidate)
    except OSError as exc:
        raise RuntimeError(f"Stage-11 authorization input is unreadable: {label}={candidate}") from exc


def _snapshot_from_canonical_paths(paths: Mapping[str, Path]) -> dict[str, Any]:
    """Hash the complete mapping returned by the canonical configuration loader."""
    return {
        "input_paths_yaml_sha256": _hash_required_regular_file(
            STAGE / "config/input_paths.yaml", label="config/input_paths.yaml"
        ),
        "master_config_yaml_sha256": _hash_required_regular_file(
            STAGE / "provenance/master_config.yaml", label="provenance/master_config.yaml"
        ),
        "resolved_inputs_sha256": {
            str(key): _hash_required_regular_file(Path(path), label=str(key))
            for key, path in sorted(paths.items(), key=lambda item: str(item[0]))
        },
    }


def _canonical_authorization_state() -> tuple[dict[str, Path], dict[str, Any], pd.DataFrame, dict[str, Any]]:
    """Load the only paths/config/folds that may produce an authorization receipt."""
    from .folds import canonical_folds

    paths, config = load_configuration()
    folds = canonical_folds(
        pd.read_csv(paths["frozen_loo_manifest"], sep="\t"),
        pd.read_csv(paths["khalsa_manifest"], sep="\t"),
        pd.read_csv(paths["mikolajewicz_gene_scores"], sep="\t")[["sample"]],
    )
    return paths, config, folds, _snapshot_from_canonical_paths(paths)


def _production_input_hashes() -> dict[str, Any]:
    """Reconstruct the complete current authoritative input/config snapshot."""
    return _canonical_authorization_state()[3]


def _fold_scalar(value: Any) -> list[str]:
    """Stable, dtype-independent scalar representation for canonical fold records."""
    try:
        if pd.isna(value):
            return ["na", ""]
    except (TypeError, ValueError):
        pass
    if isinstance(value, bool):
        return ["value", "true" if value else "false"]
    if isinstance(value, numbers.Integral):
        return ["value", str(int(value))]
    if isinstance(value, numbers.Real):
        numeric = float(value)
        if math.isfinite(numeric) and numeric.is_integer():
            return ["value", str(int(numeric))]
    return ["value", str(value)]


def canonical_fold_fingerprint(folds: pd.DataFrame) -> str:
    """Fingerprint fold semantics while ignoring DataFrame index/dtype presentation."""
    missing = [column for column in _CANONICAL_FOLD_COLUMNS if column not in folds.columns]
    if missing:
        raise ValueError(f"Fold state lacks canonical columns: {missing}")
    records = [
        {column: _fold_scalar(value) for column, value in row.items()}
        for row in folds.loc[:, _CANONICAL_FOLD_COLUMNS].to_dict("records")
    ]
    records.sort(key=lambda row: tuple(tuple(row[column]) for column in _CANONICAL_FOLD_COLUMNS))
    return stable_hash({"columns": list(_CANONICAL_FOLD_COLUMNS), "records": records})


def _canonical_paths_match(supplied: Mapping[str, Path], canonical: Mapping[str, Path]) -> bool:
    if {str(key) for key in supplied} != {str(key) for key in canonical}:
        return False
    try:
        return all(
            Path(supplied[key]).resolve(strict=True) == Path(canonical[key]).resolve(strict=True)
            for key in canonical
        )
    except (OSError, RuntimeError):
        return False


def _assert_authorized_final_state(
    supplied_folds: pd.DataFrame,
    supplied_paths: Mapping[str, Path],
    supplied_config: Mapping[str, Any],
) -> tuple[dict[str, Path], dict[str, Any], pd.DataFrame, dict[str, Any], str]:
    """Reject caller state that is not the fresh canonical authorization state."""
    paths, config, folds, snapshot = _canonical_authorization_state()
    failures = []
    if not _canonical_paths_match(supplied_paths, paths):
        failures.append("paths")
    if stable_hash(supplied_config) != stable_hash(config):
        failures.append("config")
    canonical_fingerprint = canonical_fold_fingerprint(folds)
    try:
        supplied_fingerprint = canonical_fold_fingerprint(supplied_folds)
    except Exception as exc:
        raise RuntimeError("Stage-11 authorized-final manifest state is noncanonical: invalid caller folds") from exc
    if supplied_fingerprint != canonical_fingerprint:
        failures.append("folds")
    if failures:
        raise RuntimeError(
            "Stage-11 authorized-final manifest state is stale or noncanonical "
            f"(mismatched: {failures}); rebuild from fresh canonical inputs"
        )
    return paths, config, folds, snapshot, canonical_fingerprint


def _invalidate_production_authorization() -> None:
    """Fail closed before any manifest rebuild or qualification-only rewrite."""
    _production_authorization_path().unlink(missing_ok=True)


def _write_production_authorization(
    expected_snapshot: Mapping[str, Any], expected_fold_fingerprint: str,
) -> None:
    """Bind final production authorization to code, inputs, config, and manifest bytes."""
    _, _, folds, current_snapshot = _canonical_authorization_state()
    current_fold_fingerprint = canonical_fold_fingerprint(folds)
    if current_snapshot != dict(expected_snapshot) or current_fold_fingerprint != expected_fold_fingerprint:
        raise RuntimeError(
            "Stage-11 production authorization is stale: authoritative input/config or fold state changed during manifest build"
        )
    payload = {
        "schema": PRODUCTION_AUTHORIZATION_SCHEMA,
        "authorized": True,
        "qualification_implementation_fingerprint": (
            _current_qualification_implementation_fingerprint()
        ),
        "authorization_code": _authorization_code_hashes(),
        "authoritative_inputs": current_snapshot,
        "canonical_folds_fingerprint": current_fold_fingerprint,
        "manifest_sha256": _production_manifest_hashes(),
    }
    atomic_json(_production_authorization_path(), payload)


def assert_production_manifests_authorized() -> None:
    """Refuse execution of stale, prequalification, or modified production manifests."""
    path = _production_authorization_path()
    try:
        receipt = json.loads(path.read_text())
    except Exception as exc:
        raise RuntimeError(
            "Stage-11 production is not authorized: missing or unreadable "
            f"{PRODUCTION_AUTHORIZATION_FILE}. Rebuild final manifests only after "
            "all eight model/smoke qualifications pass."
        ) from exc

    if receipt.get("schema") != PRODUCTION_AUTHORIZATION_SCHEMA or receipt.get("authorized") is not True:
        raise RuntimeError("Stage-11 production authorization receipt is invalid or obsolete")
    if not four_method_barrier_validated():
        raise RuntimeError(
            "Stage-11 production is blocked: the current four-method/eight-case "
            "qualification barrier is not satisfied"
        )
    if receipt.get("qualification_implementation_fingerprint") != _current_qualification_implementation_fingerprint():
        raise RuntimeError(
            "Stage-11 production authorization is stale: qualification implementation changed"
        )
    if receipt.get("authorization_code") != _authorization_code_hashes():
        raise RuntimeError(
            "Stage-11 production authorization is stale: authorization/runner code changed"
        )

    try:
        _, _, current_folds, current_inputs = _canonical_authorization_state()
        current_fold_fingerprint = canonical_fold_fingerprint(current_folds)
    except Exception as exc:
        raise RuntimeError(
            "Stage-11 production authorization is stale: authoritative input/config "
            "snapshot cannot be reconstructed"
        ) from exc
    if receipt.get("authoritative_inputs") != current_inputs:
        raise RuntimeError(
            "Stage-11 production authorization is stale: authoritative input/config "
            "bytes no longer match the authorized final-manifest snapshot"
        )
    if receipt.get("canonical_folds_fingerprint") != current_fold_fingerprint:
        raise RuntimeError(
            "Stage-11 production authorization is stale: canonical fold state no longer matches the authorized snapshot"
        )

    try:
        current_hashes = _production_manifest_hashes()
    except Exception as exc:
        raise RuntimeError(
            "Stage-11 production authorization is stale: one or more authorized manifests are missing"
        ) from exc
    if receipt.get("manifest_sha256") != current_hashes:
        raise RuntimeError(
            "Stage-11 production authorization is stale: manifest bytes no longer match "
            "the authorized final-manifest snapshot"
        )


def parent_reaction_universe(paths: Mapping[str, Path]) -> list[str]:
    ids = [reaction.id for reaction in load_model(paths["model"]).reactions]
    if len(ids) != len(set(ids)): raise ValueError("Parent model has duplicate reaction IDs")
    return ids


def stage4_reaction_universe(paths: Mapping[str, Path]) -> list[str]:
    values = pd.read_csv(paths["stage4_canonical_reactions"])["reaction_id"].astype(str).tolist()
    if len(values) != len(set(values)): raise ValueError("Stage-4 universe has duplicates")
    return values


def _current_full_parent_qualification_reference(folds: pd.DataFrame, algorithm: str, tumor: str,
                                                 paths: Mapping[str, Path], config: Mapping[str, Any]) -> dict[str, Any] | None:
    """Build the canonical current full-parent reference without consulting stale manifests."""
    protocol = "khalsa_all_mikolajewicz_validation"
    selected = folds.loc[
        folds.protocol.eq(protocol) & folds.fold.astype(str).eq("all")
        & folds.role.eq("train") & folds.tumor.eq(tumor)
    ]
    samples = sorted(selected["sample"].astype(str))
    if not samples:
        return None
    expert_id, expert_hash, spec = expert_identity(
        algorithm, tumor, samples, paths, config, PRIMARY_EXPERT_VARIANT
    )
    return {
        "protocol": protocol, "fold": "all", "algorithm": algorithm, "tumor": tumor,
        "expert_variant": PRIMARY_EXPERT_VARIANT, "condition_type": "rna_only",
        "context_id": f"{tumor}__rna_only__all", "dmi_mouse_id": None,
        "rna_training_slice_hash": spec["expression_slice_sha256"],
        "generation_spec_identity": spec["generation_spec_identity"],
        "rna_training_samples": ";".join(samples), "expert_id": expert_id, "expert_hash": expert_hash,
    }


def build_expert_references(folds: pd.DataFrame, paths: Mapping[str, Path], config: Mapping[str, Any], *, include_anchor_ablation: bool = False) -> pd.DataFrame:
    variants = [PRIMARY_EXPERT_VARIANT] + ([ANCHOR_ABLATION_VARIANT] if include_anchor_ablation else [])
    qualification_references = {
        (algorithm, tumor): _current_full_parent_qualification_reference(folds, algorithm, tumor, paths, config)
        for algorithm in config["algorithms"] for tumor in ("CT2A", "GL261")
    }
    global_authorized = four_method_barrier_validated(qualification_references)
    rows = []
    for protocol in config["protocols"]:
        for fold in sorted(folds.loc[folds.protocol.eq(protocol), "fold"].unique()):
            block = folds.loc[folds.protocol.eq(protocol) & folds.fold.eq(fold) & folds.role.eq("train")]
            for algorithm in config["algorithms"]:
                for tumor in ("CT2A", "GL261"):
                    samples = sorted(block.loc[block.tumor.eq(tumor), "sample"].astype(str))
                    for variant in variants:
                        expert_id, expert_hash, spec = expert_identity(algorithm, tumor, samples, paths, config, variant)
                        context_id = f"{tumor}__rna_only__all" if protocol == "khalsa_all_mikolajewicz_validation" and str(fold) == "all" else f"{tumor}__rna_only__{expert_hash[:16]}"
                        status = effective_status(
                            algorithm, tumor=tumor,
                            expected_record=qualification_references[(algorithm, tumor)],
                            global_authorized=global_authorized,
                        )
                        is_primary = variant == PRIMARY_EXPERT_VARIANT
                        planned = status == "primary_runnable" and (is_primary or include_anchor_ablation)
                        rows.append({
                            "protocol": protocol, "fold": fold, "algorithm": algorithm, "tumor": tumor,
                            "expert_variant": variant, "condition_type": "rna_only", "context_id": context_id,
                            "dmi_mouse_id": None, "rna_training_slice_hash": spec["expression_slice_sha256"],
                            "generation_spec_identity": spec["generation_spec_identity"],
                            "manifest_role": "primary" if is_primary else "optional_ablation",
                            "rna_training_samples": ";".join(samples), "expert_id": expert_id,
                            "expert_hash": expert_hash, "expert_spec": json.dumps(spec, sort_keys=True, separators=(",", ":")),
                            "design_status": REGISTRY[algorithm].design_status, "effective_status": status,
                            "execution_status": "planned" if planned else "blocked",
                            "blocker": "none" if planned else (REGISTRY[algorithm].blocker or status),
                        })
    return pd.DataFrame(rows).sort_values(["protocol", "fold", "algorithm", "tumor", "expert_variant"]).reset_index(drop=True)


def canonical_expert_jobs(references: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, block in references.groupby("expert_hash", sort=True):
        first = block.iloc[0].to_dict(); rows.append({key: first[key] for key in first if key not in {"protocol", "fold"}} |
            {"reference_count": len(block), "administrative_references": ";".join(sorted(f"{r.protocol}:{r.fold}" for r in block.itertuples()))})
    return pd.DataFrame(rows).sort_values("expert_id").reset_index(drop=True)


def build_ensemble_references(expert_references: pd.DataFrame, config: Mapping[str, Any], requested_vectors: int | None = None, ensemble_replicate: int = 0) -> pd.DataFrame:
    requested = int(config["sampling"]["retained_vectors"] if requested_vectors is None else requested_vectors); rows = []
    for row in expert_references.to_dict("records"):
        if row["execution_status"] != "planned":
            rows.append(row | {"ensemble_hash": "", "ensemble_seed": 0, "ensemble_spec": "", "requested_vectors": 0,
                "ensemble_replicate": ensemble_replicate, "sampling_thinning": 0, "sampling_processes": 0}); continue
        digest, seed, spec = ensemble_identity(row["expert_hash"], row["algorithm"], config, requested, ensemble_replicate)
        rows.append(row | {"ensemble_hash": digest, "ensemble_seed": seed, "ensemble_spec": json.dumps(spec, sort_keys=True, separators=(",", ":")),
            "requested_vectors": requested, "ensemble_replicate": ensemble_replicate, "sampling_thinning": int(config["sampling"]["thinning"]),
            "sampling_processes": int(config["sampling"]["processes"])})
    return pd.DataFrame(rows)


def canonical_ensemble_jobs(references: pd.DataFrame) -> pd.DataFrame:
    planned = references.loc[references.execution_status.eq("planned")]; rows = []
    for _, block in planned.groupby("ensemble_hash", sort=True):
        first = block.iloc[0].to_dict(); rows.append({key: first[key] for key in first if key not in {"protocol", "fold"}} |
            {"reference_count": len(block), "administrative_references": ";".join(sorted(f"{r.protocol}:{r.fold}" for r in block.itertuples()))})
    return pd.DataFrame(rows)


def build_projection_jobs(ensemble_jobs: pd.DataFrame, paths: Mapping[str, Path]) -> pd.DataFrame:
    universe = stage4_reaction_universe(paths); rows = []
    for row in ensemble_jobs.to_dict("records"):
        digest, spec = projection_identity(row["ensemble_hash"], "stage4_4747", universe)
        rows.append(row | {"projection_name": "stage4_4747", "projection_hash": digest, "projection_spec": json.dumps(spec, sort_keys=True, separators=(",", ":"))})
    return pd.DataFrame(rows)


def _record_training_samples(record: Mapping[str, Any]) -> list[str]:
    raw = record.get("rna_training_samples", "")
    if isinstance(raw, str):
        return sorted(sample for sample in raw.split(";") if sample)
    return sorted(map(str, raw))


def _assert_current_expert_manifest_identity(
    record: Mapping[str, Any],
    paths: Mapping[str, Path],
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """Re-derive a real manifest row from current authoritative inputs before reuse/fit."""
    spec = json.loads(record["expert_spec"])
    algorithm = str(record["algorithm"])
    tumor = str(record["tumor"])
    variant = str(record["expert_variant"])
    samples = _record_training_samples(record)
    current_id, current_hash, current_spec = expert_identity(
        algorithm, tumor, samples, paths, config, variant
    )
    expected_record = {
        "expert_id": current_id,
        "expert_hash": current_hash,
        "generation_spec_identity": current_hash,
        "rna_training_slice_hash": current_spec["expression_slice_sha256"],
    }
    mismatches = [
        key for key, expected in expected_record.items()
        if str(record.get(key)) != str(expected)
    ]
    if spec != current_spec:
        mismatches.append("expert_spec")
    if mismatches:
        raise RuntimeError(
            "Stale expert manifest identity after authoritative input/config change; "
            f"rebuild manifests before fitting/reuse (mismatched: {sorted(set(mismatches))})"
        )
    return spec


def fit_and_persist_expert(record: Mapping[str, Any], paths: Mapping[str, Path], config: Mapping[str, Any]) -> str:
    if record.get("execution_status") == "blocked": return "blocked"
    if str(record.get("condition_type", "rna_only")) == "rna_dmi" and not record.get("dmi_mapping_fingerprint"):
        raise RuntimeError("RNA+DMI fitting requires an approved mapping fingerprint")
    # A checksum-valid cached expert can be reused by legacy/minimal callers that
    # do not carry a full expert_spec.  Real Stage-11 manifest rows do carry the
    # spec and must still pass the current fitting-fingerprint guard below, so a
    # stale manifest cannot silently reuse an expert built under old semantics.
    cached_valid = valid_expert(record)
    if cached_valid and "expert_spec" not in record:
        return "skipped_valid"
    spec = json.loads(record["expert_spec"])
    algorithm = str(record["algorithm"])
    current_fingerprint = algorithm_fit_fingerprint(algorithm)
    if spec.get("algorithm_fit_fingerprint") != current_fingerprint or spec.get("algorithm_semantics_version") != FIT_SEMANTICS_VERSION[algorithm]:
        raise RuntimeError("Stale expert manifest identity after fitting-code change; rebuild manifests before fitting")
    spec = _assert_current_expert_manifest_identity(record, paths, config)
    variant = str(spec["expert_variant"]); protected = list(spec["protected_reactions"])
    if cached_valid: return "skipped_valid"
    parent = load_model(paths["model"], seed=int(spec["fit_seed"])); final_a = yaml.safe_load(paths["final_a_config"].read_text())
    expected = protected_reactions_for_variant(final_a, variant)
    if protected != expected: raise RuntimeError("Manifest protection policy differs from authoritative medium-derived policy")
    expected_policy = requirement_policy_for_variant(final_a, variant, objective_id=str(config["objective"]["biomass_reaction"]))
    if spec.get("requirement_policy") != expected_policy:
        raise RuntimeError("Manifest requirement policy differs from authoritative role policy")
    algorithm_required = algorithm_required_reactions(expected_policy)
    if spec.get("algorithm_required_reactions") != algorithm_required:
        raise RuntimeError("Manifest algorithm requirements differ from authoritative role policy")
    validate_environmental_protection(parent, final_a, protected)
    expected_source_routes = resolved_frozen_medium_routes(parent, final_a)
    apply_frozen_medium(parent, final_a)
    from .boundary_qc import (assert_no_unexpected_boundary_sources,
        assert_source_routes_within_resolved_medium, open_boundary_source_ids)
    assert_no_unexpected_boundary_sources(
        parent, expected_source_routes, context="Frozen-medium parent"
    )
    parent_source_ids = open_boundary_source_ids(parent)
    objective_id = str(config["objective"]["biomass_reaction"])
    objective_tolerance = float(config["validation"].get("objective_tolerance", 1.0e-8))
    parent_admission = validate_model(
        parent, objective_id, algorithm_required,
        objective_tolerance=objective_tolerance,
        objective_direction=str(config["objective"].get("qualification_direction", "max")),
    )
    parent_negative_control = closed_source_negative_control(
        parent, objective_id, tolerance=objective_tolerance,
        direction=str(config["objective"].get("qualification_direction", "max")),
    )
    conditioned_parent_hash = model_hash(parent)
    expression = load_expression(paths["model_gene_expression"], str(record["rna_training_samples"]).split(";"))
    if algorithm == "iMAT":
        from .imat_canonical import fit_imat_expert
        result = fit_imat_expert(
            parent, expression, config, fit_seed=int(spec["fit_seed"]),
            protected_reactions=algorithm_required,
        )
    elif algorithm == "GIMME":
        from .gimme_repair import fit_gimme_expert
        result = fit_gimme_expert(parent, expression, config, fit_seed=int(spec["fit_seed"]),
                                  protected_reactions=algorithm_required)
    elif algorithm == "CORDA":
        from .corda_canonical import fit_corda_expert
        result = fit_corda_expert(
            parent, expression, config, fit_seed=int(spec["fit_seed"]),
            protected_reactions=algorithm_required,
        )
    elif algorithm == "RIPTiDe":
        from .riptide_repair import fit_riptide_expert
        configured_bounds = config["contextualization"]["riptide"].get("set_bounds")
        if configured_bounds is not True:
            raise RuntimeError("Canonical RIPTiDe requires set_bounds=True")
        result = fit_riptide_expert(parent, expression, config, fit_seed=int(spec["fit_seed"]),
                                    protected_reactions=algorithm_required, set_bounds=True)
    else:
        result = fit_expert(algorithm, parent, expression, config, fit_seed=int(spec["fit_seed"]),
                            protected_reactions=algorithm_required, riptide_set_bounds=None)
    if model_hash(parent) != conditioned_parent_hash:
        raise RuntimeError(f"{algorithm} mutated the conditioned parent model")
    required_after_fit = postfit_required_reactions(parent, algorithm_required, algorithm, objective_id)
    if algorithm == "RIPTiDe":
        model_qualification = qualify_riptide_context_model(
            result.model, parent, objective_id, required_after_fit, expected_source_routes,
            config["validation"],
            objective_direction=str(config["objective"].get("qualification_direction", "max")),
        )
        from .qualification import validate_riptide_native_set_bounds_qualification
        model_qualification["riptide_native_set_bounds_qualification"] = (
            validate_riptide_native_set_bounds_qualification(
                result.model, parent, result.native_samples, result.algorithm_state,
                result.metadata, record, config, model_qualification["objective"],
                required_after_fit,
            )
        )
    else:
        model_qualification = qualify_context_model(
            result.model, objective_id, required_after_fit, expected_source_routes,
            config["validation"],
            objective_direction=str(config["objective"].get("qualification_direction", "max")),
        )
    if algorithm == "GIMME":
        parent_capacity = float(parent_admission["objective_value"])
        contextual_capacity = float(model_qualification["objective"]["objective_value"])
        required_fraction = float(config["contextualization"]["objective_fraction"])
        minimum_capacity = required_fraction * parent_capacity
        if contextual_capacity + objective_tolerance < minimum_capacity:
            raise RuntimeError(
                "GIMME contextualized model violates the configured parent-objective preservation contract: "
                f"contextual_capacity={contextual_capacity:.12g}, "
                f"required_minimum={minimum_capacity:.12g}, "
                f"parent_capacity={parent_capacity:.12g}, fraction={required_fraction:.12g}. "
                "Do not sample this expert; diagnose the GIMME extraction/objective mapping first."
            )
        native_biomass = float(result.metadata["native_biomass_flux"])
        native_passed = native_biomass + objective_tolerance >= minimum_capacity
        if not native_passed:
            raise RuntimeError(
                "GIMME native biological biomass flux violates the configured parent-objective contract"
            )
        model_qualification["gimme_parent_objective_contract"] = {
            "parent_capacity": parent_capacity,
            "required_fraction": required_fraction,
            "required_minimum": minimum_capacity,
            "native_solver_status": result.metadata["native_solver_status"],
            "native_biomass_flux": native_biomass,
            "gimme_optimization_objective_value": result.metadata["gimme_optimization_objective_value"],
            "native_objective_contract_passed": native_passed,
            "retained_reaction_count": result.metadata["retained_reaction_count"],
            "retained_reaction_indices_sha256": result.metadata["retained_reaction_indices_sha256"],
            "corrected_asymmetric_reversible_pair_count": result.metadata["corrected_asymmetric_reversible_pair_count"],
            "contextual_capacity": contextual_capacity,
            "reconstructed_objective_contract_passed": True,
            "mapping_audit_sha256": result.metadata["mapping_audit_sha256"],
            "mapping_audit_relative_path": result.metadata["mapping_audit_relative_path"],
            "mapping_audit_summary": result.metadata["mapping_audit_summary"],
            "biomass_objective_mapping": result.metadata["biomass_objective_mapping"],
            "passed": True,
        }
    assert_no_unexpected_boundary_sources(
        result.model,
        parent_source_ids,
        context=f"{algorithm} contextualized expert",
    )
    assert_source_routes_within_resolved_medium(
        result.model, expected_source_routes, context=f"{algorithm} contextualized expert"
    )
    parent_ids = [reaction.id for reaction in parent.reactions]; removed = [rid for rid in parent_ids if rid not in {r.id for r in result.model.reactions}]
    capability = flux_capability_table(
        result.model, parent_ids, removed, float(config["validation"]["zero_tolerance"]),
        biological_objective_id=objective_id,
        physiological_fraction=float(config["objective"]["ensemble_objective_fraction"]),
        objective_tolerance=float(config["validation"].get("objective_tolerance", 1.0e-8)),
    )
    artifact_hash = model_hash(result.model)
    if algorithm == "RIPTiDe":
        model_qualification["authoritative_contextual_expert"] = {
            "artifact_model_hash": artifact_hash,
            "reaction_count": len(result.model.reactions),
            "reaction_order_sha256": model_qualification["reaction_order_sha256"],
            "bounds_sha256": model_qualification["bounds_sha256"],
            "generation_spec_identity": str(spec["generation_spec_identity"]),
            "algorithm_fit_fingerprint": str(spec["algorithm_fit_fingerprint"]),
            "fit_seed": int(spec["fit_seed"]),
            "frozen_medium_sha256": str(spec["frozen_medium_sha256"]),
            "expression_slice_sha256": str(spec["expression_slice_sha256"]),
        }
    qualification_hash, qualification_spec = qualification_identity(
        str(spec["generation_spec_identity"]), artifact_hash, config, expected_source_routes,
        algorithm=algorithm,
    )
    model_qualification.update({
        "expert_hash": str(record["expert_hash"]),
        "generation_spec_identity": str(spec["generation_spec_identity"]),
        "artifact_model_hash": artifact_hash,
        "qualification_identity": qualification_hash,
        "qualification_spec": qualification_spec,
        "parent_model_unchanged": True,
        "parent_admission": parent_admission,
        "parent_closed_source_negative_control": parent_negative_control,
    })
    persist_expert(record, result, {"parent_model_sha256": file_sha256(paths["model"]), "dmi_measurements_used": False,
        "expert_variant": variant, "resolved_protected_reactions": protected,
        "requirement_policy": expected_policy, "algorithm_required_reactions": algorithm_required,
        "required_reactions_after_fit": required_after_fit,
        "environment_semantics": "frozen parent bounds; extracting algorithms may prune unused allowed boundary exchanges"}, parent_ids, capability, model_qualification)
    return "completed"


def acquire_expert_for_qualification(
    record: Mapping[str, Any], paths: Mapping[str, Path], config: Mapping[str, Any]
) -> tuple[str, str, bool]:
    """Resolve a current-generation expert without refitting for stale QC alone."""
    _assert_current_expert_manifest_identity(record, paths, config)
    if valid_expert(record):
        return "skipped_valid", "reused_current", True
    fit_status = fit_and_persist_expert(record, paths, config)
    if fit_status != "completed" or not valid_expert(record):
        raise RuntimeError(f"Expert persistence failed: {fit_status}")
    return fit_status, "refit_current_generation", False


def write_manifests(
    folds: pd.DataFrame,
    paths: Mapping[str, Path],
    config: Mapping[str, Any],
    *,
    include_anchor_ablation: bool = False,
    require_primary_authorization: bool = False,
) -> dict[str, pd.DataFrame]:
    # Never allow a receipt from an older snapshot to survive any rebuild.
    # A fresh authorization receipt is written only after all final manifests
    # have been atomically replaced below.
    _invalidate_production_authorization()
    authorization_snapshot = None
    authorization_fold_fingerprint = None
    if require_primary_authorization:
        # From this point onward, final manifests must be derived only from the
        # freshly reloaded canonical state.  Caller-owned objects were checked
        # for equality above, but must not remain part of the authorized build.
        (
            paths,
            config,
            folds,
            authorization_snapshot,
            authorization_fold_fingerprint,
        ) = _assert_authorized_final_state(folds, paths, config)
    refs = build_expert_references(folds, paths, config, include_anchor_ablation=include_anchor_ablation)
    if require_primary_authorization:
        from .adapters import PRIMARY_ALGORITHMS
        primary = refs.loc[refs["algorithm"].isin(PRIMARY_ALGORITHMS)]
        failures = primary.loc[
            ~primary["effective_status"].eq("primary_runnable")
            | ~primary["execution_status"].eq("planned")
        ]
        if len(failures):
            cases = sorted(set(zip(failures["algorithm"], failures["tumor"])))
            raise RuntimeError(
                "Final production manifest is blocked by unauthorized or stale primary cases: "
                f"{cases}"
            )
    experts = canonical_expert_jobs(refs)
    ensemble_refs = build_ensemble_references(refs, config); ensembles = canonical_ensemble_jobs(ensemble_refs); projections = build_projection_jobs(ensembles, paths)
    outputs = {"expert_references": refs, "expert_jobs": experts, "ensemble_references": ensemble_refs, "ensemble_jobs": ensembles, "projection_jobs": projections}
    for name, frame in outputs.items(): atomic_table(STAGE / "manifests" / f"{name}.tsv", frame)
    if require_primary_authorization:
        assert authorization_snapshot is not None
        assert authorization_fold_fingerprint is not None
        _write_production_authorization(
            authorization_snapshot, authorization_fold_fingerprint
        )
    return outputs
