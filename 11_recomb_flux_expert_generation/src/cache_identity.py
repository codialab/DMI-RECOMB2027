"""Scientific cache identities for the portable Stage-11 schema-2 artifacts."""
from __future__ import annotations
import hashlib
import inspect
from functools import lru_cache
from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Mapping, Sequence
import numpy as np
import pandas as pd
import yaml
from .adapters import REGISTRY
from .core import STAGE, file_sha256, seed_for, stable_hash, validate_samples
from .models import (apply_biomass_floor, apply_frozen_medium, expression_thresholds,
    load_expression, load_model, model_hash, protected_reactions_for_variant,
    requirement_policy_for_variant, algorithm_required_reactions, reaction_expression,
    validate_environmental_protection)

EXPERT_CACHE_SCHEMA = 4
ENSEMBLE_CACHE_SCHEMA = 3
PROJECTION_CACHE_SCHEMA = 3
SCORING_CACHE_SCHEMA = 2
FIT_SEMANTICS_VERSION = {"GIMME": 4, "iMAT": 9, "INIT_tINIT": 2, "CORDA": 5,
    "FASTCORE_FASTCORMICS": 2, "E-Flux": 2, "RIPTiDe": 6, "RegrEx": 2}


def _version(name: str) -> str:
    try: return version(name)
    except PackageNotFoundError: return "not_installed"


def _source_hash(*objects: object) -> str:
    sources = []
    for item in objects:
        try: sources.append(inspect.getsource(item))
        except (OSError, TypeError): sources.append(repr(item))
    return stable_hash(sources)


def _medium_application_fingerprint() -> str:
    """Fingerprint the full orientation-aware medium/source semantics."""
    from . import boundary_qc
    return _source_hash(
        apply_frozen_medium,
        requirement_policy_for_variant,
        algorithm_required_reactions,
        validate_environmental_protection,
        boundary_qc._boundary_coefficient,
        boundary_qc.boundary_source_direction,
        boundary_qc.close_boundary_source_direction,
        boundary_qc.assert_resolved_source_routes,
        boundary_qc.assert_source_routes_within_resolved_medium,
    )


def _external_algorithm_source(algorithm: str) -> dict[str, str]:
    targets = {
        "GIMME": ("troppo.methods.reconstruction.gimme", "GIMME"),
        "iMAT": ("troppo.methods.reconstruction.imat", "IMAT"),
        "INIT_tINIT": ("troppo.methods.reconstruction.tINIT", "tINIT"),
        "CORDA": ("corda", "CORDA"),
        "FASTCORE_FASTCORMICS": ("troppo.methods.reconstruction.fastcore", "FASTcore"),
        "RIPTiDe": ("riptide", "contextualize"),
    }
    if algorithm not in targets: return {"entry_point_source_sha256": "implemented_or_blocked_locally"}
    module_name, object_name = targets[algorithm]
    try:
        item = getattr(import_module(module_name), object_name)
        return {"entry_point": f"{module_name}.{object_name}", "entry_point_signature": str(inspect.signature(item)),
                "entry_point_source_sha256": _source_hash(item)}
    except Exception as exc:
        return {"entry_point": f"{module_name}.{object_name}", "entry_point_source_sha256": f"unavailable:{type(exc).__name__}"}


def algorithm_fit_fingerprint(algorithm: str) -> str:
    # Keep the pre-existing shared fingerprint calculation for unaffected
    # algorithms. iMAT, CORDA, and RIPTiDe use isolated modules so their fitting
    # changes do not invalidate GIMME or other unrelated methods.
    if algorithm == "iMAT":
        from . import imat_canonical
        return stable_hash({
            "repair_module": "src/imat_canonical.py",
            "complete_implementation_file_sha256": file_sha256(STAGE / "src/imat_canonical.py"),
            "external": {
                "troppo_version": _version("troppo"),
                "cobamp_version": _version("cobamp"),
                "troppo_imat": _external_algorithm_source("iMAT"),
            },
        })
    if algorithm == "CORDA":
        return stable_hash({
            "repair_module": "src/corda_canonical.py",
            "complete_implementation_file_sha256": file_sha256(
                STAGE / "src/corda_canonical.py"
            ),
            "external": {
                "corda_version": _version("corda"),
                "corda": _external_algorithm_source("CORDA"),
            },
        })
    if algorithm == "GIMME":
        return stable_hash({
            "repair_module": "src/gimme_repair.py",
            "complete_implementation_file_sha256": file_sha256(STAGE / "src/gimme_repair.py"),
            "external": _external_algorithm_source(algorithm),
        })
    if algorithm == "RIPTiDe":
        from . import riptide_repair
        return stable_hash({
            "repair_module": "src/riptide_repair.py",
            "complete_implementation_file_sha256": file_sha256(
                STAGE / "src/riptide_repair.py"
            ),
            "external": _external_algorithm_source(algorithm),
        })

    from . import expert_fitting
    helpers = {
        "GIMME": (expert_fitting.fit_expert, expert_fitting.gimme_expression_vector, expert_fitting._matrix, expert_fitting._reduce),
        "INIT_tINIT": (expert_fitting.fit_expert, expert_fitting._matrix, expert_fitting._reduce),
        "FASTCORE_FASTCORMICS": (expert_fitting.fit_expert, expert_fitting._numpy_compatible_fastcore, expert_fitting._matrix, expert_fitting._reduce),
        "E-Flux": (expert_fitting.fit_expert, expert_fitting._fit_eflux, expert_fitting._eflux_bounds),
        "RegrEx": (expert_fitting.fit_expert,),
    }[algorithm]
    return stable_hash({"expert_fitting_file_sha256": file_sha256(STAGE / "src/expert_fitting.py"),
        "algorithm_helpers_sha256": _source_hash(*helpers), "external": _external_algorithm_source(algorithm)})


@lru_cache(maxsize=None)
def _expression_slice_hash(path_text: str, samples: tuple[str, ...], content_sha256: str) -> str:
    frame = pd.read_csv(Path(path_text), sep="\t", dtype={"model_gene_id": str})
    missing = set(samples) - set(frame.columns)
    if missing: raise KeyError(f"Expression samples absent: {sorted(missing)}")
    selected = frame[["model_gene_id", *samples]].sort_values("model_gene_id")
    if not selected.model_gene_id.is_unique: raise ValueError("Expression model_gene_id must be unique")
    digest = hashlib.sha256(); digest.update("\0".join(selected.model_gene_id.astype(str)).encode()); digest.update("\0".join(samples).encode())
    digest.update(np.ascontiguousarray(selected[list(samples)].to_numpy(dtype="<f8")).tobytes()); return digest.hexdigest()


def expression_slice_hash(path_text: str, samples: tuple[str, ...]) -> str:
    """Hash one RNA training slice while invalidating the cache on file changes."""
    path = Path(path_text)
    return _expression_slice_hash(path_text, samples, file_sha256(path))


@lru_cache(maxsize=None)
def _parent_semantic_hash(path_text: str, content_sha256: str) -> str:
    return model_hash(load_model(Path(path_text)))


@lru_cache(maxsize=None)
def _medium_parent_hash(path_text: str, content_sha256: str, medium_sha256: str, medium_yaml: str) -> str:
    model = load_model(Path(path_text)); apply_frozen_medium(model, yaml.safe_load(medium_yaml)); return model_hash(model)


def fit_parameters(algorithm: str, config: Mapping[str, Any]) -> dict[str, Any]:
    c = config["contextualization"]
    common = {"pooled_statistic": c["pooled_statistic"]}
    values = {
      "GIMME": {"objective_fraction": c["objective_fraction"], "expression_threshold": c["gimme"]["expression_threshold"], "numerical_epsilon": c["epsilon"]},
      "iMAT": {"lower_quantile": c["lower_quantile"], "upper_quantile": c["upper_quantile"],
               "epsilon": c["imat"]["epsilon"], "activity_tolerance": c["imat"]["tolerance"],
               "mip_gap": c["imat"]["mip_gap"],
               "biomass_floor_fraction": c["imat"]["biomass_floor_fraction"],
               "imat_properties_overrides": ["exp_vector", "exp_thresholds", "tolerance", "epsilon", "core"],
               "core_reactions": [str(config["objective"]["biomass_reaction"])],
               "biomass_preservation_policy": "biomass_only_core_union_with_hard_parent_bmax_fraction_constraint"},
      "INIT_tINIT": {"lower_quantile": c["lower_quantile"], **c["tinit"]},
      "CORDA": {"lower_quantile": c["lower_quantile"], "upper_quantile": c["upper_quantile"],
                "algorithm_parameters": "corda_package_defaults",
                "manually_high_confidence_reactions": [str(config["objective"]["biomass_reaction"])],
                "biomass_preservation_policy": "corda_high_confidence"},
      "FASTCORE_FASTCORMICS": {"upper_quantile": c["upper_quantile"], "epsilon": c["epsilon"]},
      "E-Flux": dict(c["eflux"]), "RIPTiDe": dict(c["riptide"]), "RegrEx": dict(c["regrex"]),
    }
    return common | values[algorithm]


def expert_identity(algorithm: str, tumor: str, training_samples: Sequence[str], paths: Mapping[str, Path],
                    config: Mapping[str, Any], expert_variant: str | None = None) -> tuple[str, str, dict[str, Any]]:
    samples = tuple(sorted(map(str, training_samples))); variant = expert_variant or str(config["primary_expert_variant"])
    final_a = yaml.safe_load(paths["final_a_config"].read_text())
    role_policy = requirement_policy_for_variant(final_a, variant, objective_id=str(config["objective"]["biomass_reaction"]))
    protected = protected_reactions_for_variant(final_a, variant)
    algorithm_required = algorithm_required_reactions(role_policy)
    parent_sha = file_sha256(paths["model"]); medium = final_a["medium"]; medium_hash = stable_hash(medium)
    medium_yaml = yaml.safe_dump({"medium": medium}, sort_keys=True)
    spec = {"schema": EXPERT_CACHE_SCHEMA, "algorithm": algorithm, "artifact_kind": REGISTRY[algorithm].artifact_kind,
      "algorithm_package": REGISTRY[algorithm].package, "algorithm_package_version": _version(REGISTRY[algorithm].package),
      "algorithm_semantics_version": FIT_SEMANTICS_VERSION[algorithm], "algorithm_fit_fingerprint": algorithm_fit_fingerprint(algorithm),
      "tumor": tumor, "rna_training_samples": list(samples), "expression_slice_sha256": expression_slice_hash(str(paths["model_gene_expression"]), samples),
      "expression_preprocessing_sha256": _source_hash(load_expression, reaction_expression, expression_thresholds),
      "parent_model_sha256": parent_sha, "parent_model_semantic_hash": _parent_semantic_hash(str(paths["model"]), parent_sha),
      "frozen_medium": medium, "frozen_medium_sha256": medium_hash,
      "medium_conditioned_parent_model_hash": _medium_parent_hash(str(paths["model"]), parent_sha, medium_hash, medium_yaml),
      "medium_application_sha256": _medium_application_fingerprint(),
      "expert_variant": variant, "protected_reactions": protected,
      "requirement_policy": role_policy, "algorithm_required_reactions": algorithm_required,
      "protection_policy": config["contextualization"]["protection_policy"],
      "fit_parameters": fit_parameters(algorithm, config), "gurobi": config["gurobi"],
      "fit_seed": seed_for(int(config["master_seed"]), "expert_fit", algorithm, tumor, stable_hash(samples), variant)}
    # This identity is computable before fitting. The resulting model hash is
    # deliberately persisted separately and never feeds this specification.
    digest = stable_hash(spec)
    spec = spec | {"generation_spec_identity": digest}
    return f"{algorithm}__{tumor}__{variant}__{digest[:16]}", digest, spec


def ensemble_implementation_fingerprint(algorithm: str) -> str:
    payload = {"ensembles_file_sha256": file_sha256(STAGE / "src/ensembles.py"),
        "biomass_floor_sha256": _source_hash(apply_biomass_floor), "validation_sha256": _source_hash(validate_samples),
        "ensemble_method": REGISTRY[algorithm].ensemble_method}
    if algorithm == "RIPTiDe":
        from . import riptide_native
        payload["riptide_native_sampling_sha256"] = _source_hash(
            riptide_native.api_description, riptide_native.native_resample
        )
    return stable_hash(payload)


def qualification_identity(generation_spec_identity: str, artifact_model_hash: str,
                           config: Mapping[str, Any], expected_source_routes: Mapping[str, Any],
                           *, algorithm: str | None = None) -> tuple[str, dict[str, Any]]:
    """Identity for validating an immutable model, not for generating it."""
    from . import boundary_qc, qc
    from .models import (closed_source_negative_control, evaluate_biological_objective,
        qualify_context_model)
    spec = {
        "generation_spec_identity": str(generation_spec_identity),
        "artifact_model_hash": str(artifact_model_hash),
        "biological_objective": {"reaction": config["objective"]["biomass_reaction"],
                                 "direction": config["objective"].get("qualification_direction", "max")},
        "expected_source_routes": expected_source_routes,
        "objective_tolerance": float(config["validation"].get("objective_tolerance", 1.0e-8)),
        "physiological_fraction": float(config["objective"]["ensemble_objective_fraction"]),
        "numerical_tolerances": dict(config["validation"]),
        "qualification_code": _source_hash(boundary_qc.assert_resolved_source_routes,
                                             boundary_qc.close_boundary_source_direction_for_negative_control,
                                             evaluate_biological_objective,
                                             qc.flux_capability_table,
                                             closed_source_negative_control,
                                             qualify_context_model),
    }
    if algorithm == "RIPTiDe":
        from .qualification import (
            RIPTIDE_NATIVE_SET_BOUNDS_QUALIFICATION_POLICY,
            RIPTIDE_NATIVE_SET_BOUNDS_QUALIFICATION_POLICY_VERSION,
            validate_riptide_native_set_bounds_qualification,
        )
        from .models import (
            RIPTIDE_CLOSED_SOURCE_QC_POLICY_NAME,
            RIPTIDE_CLOSED_SOURCE_QC_POLICY_VERSION,
            qualify_riptide_context_model,
            riptide_closed_source_negative_control,
        )
        spec["riptide_closed_source_policy"] = {
            "policy_name": RIPTIDE_CLOSED_SOURCE_QC_POLICY_NAME,
            "policy_version": RIPTIDE_CLOSED_SOURCE_QC_POLICY_VERSION,
            "canonical_set_bounds": True,
            "strict_control_first": True,
            "infeasible_fallback": "relax_only_riptide_compulsory_bounds_that_exclude_parent_allowed_zero_toward_zero",
            "source_routes_must_remain_closed": True,
            "fallback_requires_optimal_zero_biomass": True,
            "code_sha256": _source_hash(
                riptide_closed_source_negative_control, qualify_riptide_context_model
            ),
        }
        spec["riptide_native_set_bounds_qualification"] = {
            "policy_name": RIPTIDE_NATIVE_SET_BOUNDS_QUALIFICATION_POLICY,
            "policy_version": RIPTIDE_NATIVE_SET_BOUNDS_QUALIFICATION_POLICY_VERSION,
            "objective_enabled": True,
            "canonical_set_bounds": True,
            "native_extrema_provenance_required": True,
            "final_bmax_fraction_is_consistency_only": True,
            "code_sha256": _source_hash(validate_riptide_native_set_bounds_qualification),
        }
    return stable_hash(spec), spec


def sampling_semantics_fingerprint(algorithm: str, config: Mapping[str, Any]) -> str:
    """Identity of sampling semantics, independent of any production artifact."""
    sampling = config["sampling"]
    payload = {
        "algorithm": str(algorithm),
        "ensemble_implementation_fingerprint": ensemble_implementation_fingerprint(algorithm),
        "ensemble_method": REGISTRY[algorithm].ensemble_method,
        "objective": {
            "reaction": config["objective"]["biomass_reaction"],
            "contextual_model_fraction": float(config["objective"]["ensemble_objective_fraction"]),
            "preprocessing": "fresh_copy_capacity_then_lower_bound_floor",
        },
        "sampling": {
            key: sampling[key]
            for key in ("sampler", "thinning", "processes", "calculation_dtype", "raw_storage_dtype", "geometry")
        },
        "validation_tolerances": config["validation"],
        "riptide_native_fraction": (
            float(config["contextualization"]["riptide"]["fraction"])
            if algorithm == "RIPTiDe" else None
        ),
    }
    return stable_hash(payload)


def smoke_run_identity(artifact_model_hash: str, algorithm: str, config: Mapping[str, Any], sample_count: int) -> tuple[str, dict[str, Any]]:
    """Sampler smoke provenance deliberately independent of biological qualification."""
    sampler = "RIPTiDe native GapSplit" if REGISTRY[str(algorithm)].ensemble_method == "native_gapsplit" else config["sampling"]["sampler"]
    spec = {
        "artifact_model_hash": str(artifact_model_hash), "algorithm": str(algorithm),
        "sample_count": int(sample_count), "sampler": sampler,
        "sampling_semantics_fingerprint": sampling_semantics_fingerprint(algorithm, config),
        "sampling": {key: config["sampling"][key] for key in ("thinning", "processes", "geometry")},
        "smoke_test_only": True, "sampling_convergence_assessed": False,
        "flux_distribution_stability_assessed": False,
    }
    return stable_hash(spec), spec


def ensemble_identity(expert_hash: str, algorithm: str, config: Mapping[str, Any], requested_vectors: int,
                      ensemble_replicate: int = 0) -> tuple[str, int, dict[str, Any]]:
    sampling = config["sampling"]; seed = seed_for(int(config["master_seed"]), "ensemble", expert_hash, ensemble_replicate)
    spec = {"schema": ENSEMBLE_CACHE_SCHEMA, "expert_hash": expert_hash, "algorithm": algorithm,
      "ensemble_method": REGISTRY[algorithm].ensemble_method, "requested_vectors": int(requested_vectors),
      "ensemble_replicate": int(ensemble_replicate), "seed": seed,
      "objective": {"reaction": config["objective"]["biomass_reaction"], "fraction": float(config["objective"]["ensemble_objective_fraction"]),
                    "preprocessing": "fresh_copy_capacity_then_lower_bound_floor"},
      "sampling": {k: sampling[k] for k in ("sampler", "thinning", "processes", "calculation_dtype", "raw_storage_dtype")},
      "geometry": sampling["geometry"], "validation_tolerances": config["validation"], "gurobi": config["gurobi"],
      "cobra_version": _version("cobra"), "gurobi_version": _version("gurobipy"),
      "sampling_implementation_fingerprint": ensemble_implementation_fingerprint(algorithm)}
    return stable_hash(spec), seed, spec


def projection_identity(ensemble_hash: str, projection_name: str, reaction_ids: Sequence[str]) -> tuple[str, dict[str, Any]]:
    spec = {"schema": PROJECTION_CACHE_SCHEMA, "ensemble_hash": ensemble_hash, "projection_name": projection_name,
            "reaction_universe_sha256": stable_hash(list(map(str, reaction_ids))), "storage_dtype": "float32",
            "projection_implementation_sha256": file_sha256(STAGE / "src/projections.py")}
    return stable_hash(spec), spec


def scoring_identity(projection_hash: str, paths: Mapping[str, Path], config: Mapping[str, Any]) -> tuple[str, dict[str, Any]]:
    spec = {"schema": SCORING_CACHE_SCHEMA, "projection_hash": projection_hash, "scoring": config["scoring"],
      "stage4_memberships_sha256": file_sha256(paths["stage4_memberships"]), "scoring_code_sha256": file_sha256(STAGE / "src/scoring.py"),
      "scope": "stage4_reaction_summary_only_not_sparseflux_weighting"}
    return stable_hash(spec), spec
