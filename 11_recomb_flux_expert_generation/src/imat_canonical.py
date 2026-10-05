"""Canonical single-run Troppo iMAT adapter for Stage 11.

Expression consistency uses the configured iMAT epsilon. Tumor growth is an
independent phenotype constraint: the one iMAT MILP must carry at least the
configured fraction of the exact conditioned parent's biomass maximum.
`BIOMASS_reaction` is the sole explicit iMAT core; no other project requirement
is promoted to core status. Extraction unions that core with Troppo's native
active-reaction set.
"""
from __future__ import annotations

import inspect
import math
import time
from importlib.metadata import version
from typing import Any, Sequence

import numpy as np
import pandas as pd

from .expert_fitting import ExpertResult, _matrix
from .models import (
    enforce_gurobi,
    evaluate_biological_objective,
    expression_thresholds,
    reaction_expression,
)


def imat_expression_vector(scores: pd.Series) -> np.ndarray:
    """Encode unmapped reactions with Troppo iMAT's conventional -1 sentinel."""
    values = scores.fillna(-1.0).to_numpy(dtype=np.float64)
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError("iMAT expression must be a finite one-dimensional vector")
    if (values < -1.0).any():
        raise ValueError("iMAT expression contains values below the unmapped sentinel -1")
    return values


def validate_retained_indices(values: Any, reaction_count: int) -> np.ndarray:
    """Semantically validate and canonicalize Troppo's retained indices.

    NumPy integer scalar/index types are accepted.  Floating values are accepted
    only when every value is exactly integral.  Boolean, duplicate, negative,
    non-finite, multidimensional, and out-of-range values fail closed.
    """
    array = np.asarray(values)
    if array.ndim != 1:
        raise RuntimeError("Troppo iMAT returned a non-one-dimensional index collection")
    if array.dtype.kind == "b":
        raise RuntimeError("Troppo iMAT returned boolean values instead of reaction indices")
    if array.dtype.kind in "iu":
        canonical = array.astype(np.int64, copy=False)
    elif array.dtype.kind == "f":
        if not np.isfinite(array).all() or not np.equal(array, np.rint(array)).all():
            raise RuntimeError("Troppo iMAT returned non-integer-valued reaction indices")
        canonical = array.astype(np.int64)
    else:
        try:
            numeric = array.astype(np.float64)
        except (TypeError, ValueError) as exc:
            raise RuntimeError("Troppo iMAT returned malformed reaction indices") from exc
        if not np.isfinite(numeric).all() or not np.equal(numeric, np.rint(numeric)).all():
            raise RuntimeError("Troppo iMAT returned non-integer-valued reaction indices")
        canonical = numeric.astype(np.int64)
    if len(np.unique(canonical)) != len(canonical):
        raise RuntimeError("Troppo iMAT returned duplicate reaction indices")
    if (canonical < 0).any() or (canonical >= int(reaction_count)).any():
        raise RuntimeError("Troppo iMAT returned negative or out-of-range reaction indices")
    return canonical


def actual_cobamp_backend() -> dict[str, str]:
    """Return the backend Troppo's solver-less GenericLinearSystem will use."""
    from cobamp.core import linear_systems

    backend = str(linear_systems.get_default_solver()).upper()
    interface = linear_systems.get_solver_interfaces().get(backend)
    module = getattr(interface, "__name__", "")
    if backend != "GUROBI" or "gurobi" not in module.lower():
        raise RuntimeError(
            "Canonical Troppo iMAT requires Cobamp's actual backend to be GUROBI; "
            f"resolved backend={backend!r}, interface={module!r}"
        )
    return {"backend": backend, "interface": module}


def _native_gurobi_model(optimizer: Any, linear_system: Any) -> Any:
    """Resolve the one native model shared by Cobamp's generated objects."""
    candidates = []
    for label, owner in (("optimizer", optimizer), ("linear_system", linear_system)):
        optlang_model = getattr(owner, "model", None)
        native = getattr(optlang_model, "problem", None)
        if native is not None:
            candidates.append((label, native))
    if not candidates:
        raise RuntimeError("Cobamp did not expose the generated native solver model")
    native = candidates[0][1]
    if any(candidate is not native for _, candidate in candidates[1:]):
        labels = [label for label, _ in candidates]
        raise RuntimeError(
            "Cobamp returned ambiguous native solver models through "
            f"{labels}; refusing to configure an uncertain optimization model"
        )
    if not hasattr(native, "Params"):
        raise RuntimeError("Cobamp's generated native solver model has no Gurobi Params")
    return native


def _parameter_matches(name: str, expected: Any, actual: Any) -> bool:
    """Compare native Gurobi readbacks according to their parameter type."""
    if name == "OutputFlag":
        try:
            expected_float = float(expected)
            actual_float = float(actual)
        except (TypeError, ValueError, OverflowError):
            return False
        return (
            math.isfinite(expected_float)
            and math.isfinite(actual_float)
            and expected_float.is_integer()
            and actual_float.is_integer()
            and int(expected_float) in (0, 1)
            and int(actual_float) in (0, 1)
            and bool(int(actual_float)) == bool(int(expected_float))
        )
    if name in {"Threads", "Seed"}:
        if isinstance(expected, bool) or isinstance(actual, bool):
            return False
        try:
            expected_float = float(expected)
            actual_float = float(actual)
        except (TypeError, ValueError, OverflowError):
            return False
        return (
            math.isfinite(expected_float)
            and math.isfinite(actual_float)
            and expected_float.is_integer()
            and actual_float.is_integer()
            and int(actual_float) == int(expected_float)
        )
    try:
        expected_float = float(expected)
        actual_float = float(actual)
    except (TypeError, ValueError, OverflowError):
        return False
    return bool(
        math.isfinite(expected_float)
        and math.isfinite(actual_float)
        and math.isclose(
            actual_float, expected_float, rel_tol=1.0e-12, abs_tol=1.0e-15
        )
    )


def _configure_cobamp_gurobi(
    optimizer: Any,
    linear_system: Any,
    config: dict[str, Any],
    fit_seed: int,
) -> dict[str, Any]:
    """Configure and verify the actual Cobamp/Optlang Gurobi model.

    Cobamp initializes each GenericLinearSystem with its own solver
    configuration, including 1e-9 feasibility/optimality/integrality
    tolerances.  Process-global ``gurobipy.setParam`` calls therefore do not
    prove what the native model will use.  Apply Stage-11's declared policy to
    the native model after Troppo constructs it and read the values back before
    the one canonical iMAT optimization.
    """
    configured = dict(config.get("gurobi", {}))
    requested: dict[str, Any] = {
        "OutputFlag": 0,
        "Seed": int(fit_seed),
    }
    for name in ("Threads", "MIPGap", "FeasibilityTol", "OptimalityTol", "TimeLimit"):
        if name in configured and configured[name] is not None:
            requested[name] = configured[name]

    # The validated iMAT reconstruction gap is deliberately looser than the
    # historical project-wide setting. Keep it local to this adapter so an
    # iMAT runtime choice cannot silently change other algorithms.
    imat_mip_gap = config["contextualization"]["imat"].get("mip_gap")
    if imat_mip_gap is not None:
        requested["MIPGap"] = float(imat_mip_gap)

    native = _native_gurobi_model(optimizer, linear_system)
    for name, value in requested.items():
        setattr(native.Params, name, value)

    observed_names = (
        "OutputFlag", "Threads", "Seed", "MIPGap", "FeasibilityTol",
        "OptimalityTol", "IntFeasTol", "TimeLimit",
    )
    observed = {name: getattr(native.Params, name) for name in observed_names}
    for name, expected in requested.items():
        actual = observed[name]
        if not _parameter_matches(name, expected, actual):
            raise RuntimeError(
                f"Cobamp/Gurobi parameter {name} did not stick: "
                f"requested={expected!r}, observed={actual!r}"
            )
    return {"requested": requested, "observed": observed}


def _reduce_to_active(model, retained: np.ndarray):
    """Extract exactly the canonical raw-retained/core-union topology."""
    keep_ids = {model.reactions[int(index)].id for index in retained}
    result = model.copy()
    result.remove_reactions(
        [reaction for reaction in result.reactions if reaction.id not in keep_ids],
        remove_orphans=True,
    )
    return result


def _add_biomass_floor_constraint(
    linear_system: Any,
    biomass_index: int,
    biomass_floor: float,
) -> dict[str, Any]:
    """Add and verify the temporary phenotype row in the one iMAT MILP."""
    variable_count = len(linear_system.model.variables)
    if not 0 <= int(biomass_index) < variable_count:
        raise RuntimeError("Biomass index is outside the generated iMAT variables")
    if not math.isfinite(float(biomass_floor)) or float(biomass_floor) <= 0.0:
        raise RuntimeError(f"iMAT biomass floor is not finite and positive: {biomass_floor}")
    row = np.zeros((1, variable_count), dtype=np.float64)
    row[0, int(biomass_index)] = 1.0
    constraints = linear_system.add_rows_to_model(
        row,
        [float(biomass_floor)],
        [None],
        only_nonzero=True,
        names=["stage11_parent_biomass_floor"],
    )
    constraint = constraints[0]
    observed_lower = float(constraint.lb)
    if not math.isclose(
        observed_lower, float(biomass_floor), rel_tol=1.0e-12, abs_tol=1.0e-15
    ):
        raise RuntimeError(
            "Generated iMAT biomass floor did not stick: "
            f"requested={biomass_floor}, observed={observed_lower}"
        )
    return {
        "name": str(constraint.name),
        "biomass_index": int(biomass_index),
        "lower_bound": observed_lower,
        "upper_bound": None if constraint.ub is None else float(constraint.ub),
    }


def _validate_extracted_biomass_floor(
    extracted_bmax: float,
    biomass_floor: float,
    objective_tolerance: float,
) -> bool:
    """Fail closed when extraction loses the native iMAT phenotype floor."""
    passed = bool(float(extracted_bmax) + float(objective_tolerance) >= float(biomass_floor))
    if not passed:
        raise RuntimeError(
            "iMAT extraction inconsistency: the native constrained solution met "
            "the parent-Bmax phenotype floor but the extracted topology failed to "
            f"preserve it: extracted_bmax={extracted_bmax}, "
            f"floor={biomass_floor}, tolerance={objective_tolerance}"
        )
    return passed


def _imat_properties_defaults(properties_class: type) -> dict[str, Any]:
    """Fail closed if the pinned Troppo iMAT defaults drift unexpectedly."""
    signature = inspect.signature(properties_class)
    expected = {"core": None, "tolerance": 1.0e-8, "epsilon": 1}
    observed: dict[str, Any] = {}
    for name, value in expected.items():
        parameter = signature.parameters.get(name)
        if parameter is None or parameter.default is inspect.Parameter.empty:
            raise RuntimeError(f"Troppo IMATProperties no longer defines default {name!r}")
        observed[name] = parameter.default
        if name == "core":
            if parameter.default is not None:
                raise RuntimeError(f"Unexpected Troppo iMAT core default: {parameter.default!r}")
        elif not math.isclose(float(parameter.default), float(value), rel_tol=0.0, abs_tol=0.0):
            raise RuntimeError(
                f"Unexpected Troppo iMAT {name} default: {parameter.default!r} != {value!r}"
            )
    return observed


def fit_imat_expert(
    model,
    expression: pd.DataFrame,
    config: dict[str, Any],
    *,
    fit_seed: int,
    protected_reactions: Sequence[str],
) -> ExpertResult:
    """Run canonical Troppo iMAT exactly once and construct its context GEM."""
    from troppo.methods.reconstruction.imat import IMAT, IMATProperties

    started = time.perf_counter()
    enforce_gurobi(model, seed=fit_seed)
    backend = actual_cobamp_backend()
    solver_parameters: dict[str, Any] = {}
    declared_requirements = list(dict.fromkeys(map(str, protected_reactions)))
    reaction_ids = [reaction.id for reaction in model.reactions]
    missing_declared = sorted(set(declared_requirements) - set(reaction_ids))
    if missing_declared:
        raise KeyError(f"Declared project requirements absent from the parent: {missing_declared}")
    biomass_id = str(config["objective"]["biomass_reaction"])
    if biomass_id not in reaction_ids:
        raise KeyError(f"Biomass reaction absent from parent: {biomass_id}")
    boundary_ids = {reaction.id for reaction in model.boundary}
    if biomass_id in boundary_ids:
        raise RuntimeError(f"Biomass reaction may not be a boundary reaction: {biomass_id}")
    biomass_index = reaction_ids.index(biomass_id)

    contextualization = config["contextualization"]
    imat_config = contextualization["imat"]
    epsilon = float(imat_config["epsilon"])
    activity_tolerance = float(imat_config["tolerance"])
    biomass_fraction = float(imat_config["biomass_floor_fraction"])
    if epsilon <= 0.0 or not math.isfinite(epsilon):
        raise ValueError(f"iMAT epsilon must be finite and positive: {epsilon}")
    if activity_tolerance != 1.0e-8:
        raise ValueError(
            "Canonical Stage-11 iMAT extraction tolerance must remain 1e-8: "
            f"{activity_tolerance}"
        )
    if not 0.0 < biomass_fraction <= 1.0:
        raise ValueError("iMAT biomass floor fraction must be in (0, 1]")
    objective_tolerance = float(config["validation"].get("objective_tolerance", 1.0e-8))
    parent_objective = evaluate_biological_objective(
        model,
        biomass_id,
        tolerance=objective_tolerance,
        direction=str(config["objective"].get("qualification_direction", "max")),
    )
    parent_bmax = float(parent_objective["capacity"])
    biomass_floor = biomass_fraction * parent_bmax
    scores = reaction_expression(model, expression, contextualization["pooled_statistic"])
    low, high = expression_thresholds(
        scores,
        contextualization["lower_quantile"],
        contextualization["upper_quantile"],
    )
    vector = imat_expression_vector(scores)
    if len(vector) != len(reaction_ids):
        raise RuntimeError("iMAT reaction-expression vector length does not match the parent")
    property_defaults = _imat_properties_defaults(IMATProperties)

    S, lb, ub = _matrix(model)
    # Biomass is the sole explicit iMAT core.  It is independent of expression
    # evidence and is unioned with the raw Troppo topology at extraction.
    core = np.asarray([biomass_index], dtype=np.int64)
    properties = IMATProperties(
        exp_vector=vector,
        exp_thresholds=(float(low), float(high)),
        tolerance=activity_tolerance,
        epsilon=epsilon,
        core=core.tolist(),
    )
    floor_constraint: dict[str, Any] = {}
    class _ConfiguredIMAT(IMAT):
        def generate_imat_problem(self, *args, **kwargs):
            lso, linear_system = super().generate_imat_problem(*args, **kwargs)
            floor_constraint.update(
                _add_biomass_floor_constraint(
                    linear_system, biomass_index, biomass_floor
                )
            )
            solver_parameters.update(
                _configure_cobamp_gurobi(lso, linear_system, config, fit_seed)
            )
            return lso, linear_system

    imat = _ConfiguredIMAT(S, lb, ub, properties)
    raw_keep = imat.run()  # The one and only iMAT optimization for this expert fit.
    solution = getattr(imat, "sol", None)
    if solution is None:
        raise RuntimeError("Troppo iMAT returned retained indices without a stored solution")
    try:
        status = str(solution.status()).lower()
    except Exception as exc:
        raise RuntimeError("Troppo iMAT solution status is unavailable") from exc
    if status != "optimal":
        raise RuntimeError(f"Troppo iMAT solution is not optimal: {status}")
    try:
        native_values = np.asarray(solution.x(), dtype=np.float64)
        imat_objective = float(solution.objective_value())
    except Exception as exc:
        raise RuntimeError("Troppo iMAT native solution values are unavailable") from exc
    if len(native_values) < len(reaction_ids) or not np.isfinite(native_values).all():
        raise RuntimeError("Troppo iMAT returned malformed or non-finite native values")
    native_biomass = float(native_values[biomass_index])
    native_floor_passed = bool(
        native_biomass + objective_tolerance >= biomass_floor
    )
    if not native_floor_passed:
        raise RuntimeError(
            "Troppo iMAT native biomass violates the hard parent-phenotype floor: "
            f"native={native_biomass}, floor={biomass_floor}, "
            f"tolerance={objective_tolerance}"
        )
    raw_retained = validate_retained_indices(raw_keep, len(reaction_ids))
    raw_retained_biomass = bool(biomass_index in raw_retained)
    final_retained = np.union1d(raw_retained, core).astype(np.int64, copy=False)
    if biomass_index not in final_retained:
        raise RuntimeError("Canonical iMAT final retained topology omitted biomass core")
    context_model = _reduce_to_active(model, final_retained)
    enforce_gurobi(context_model, seed=fit_seed)
    try:
        extracted_objective = evaluate_biological_objective(
            context_model,
            biomass_id,
            tolerance=objective_tolerance,
            direction=str(config["objective"].get("qualification_direction", "max")),
        )
    except Exception as exc:
        raise RuntimeError(
            "iMAT extraction inconsistency: the native constrained solution meets "
            "the biomass floor but the extracted model has no positive biomass capacity"
        ) from exc
    extracted_bmax = float(extracted_objective["capacity"])
    extracted_biomass_floor_passed = _validate_extracted_biomass_floor(
        extracted_bmax, biomass_floor, objective_tolerance
    )

    elapsed = time.perf_counter() - started
    metadata = {
        "algorithm": "iMAT",
        "implementation": "canonical Troppo 0.0.7 iMAT single-run active-reaction extraction with a hard parent-biomass phenotype constraint",
        "troppo_version": version("troppo"),
        "cobamp_version": version("cobamp"),
        "cobamp_solver_backend": backend["backend"],
        "cobamp_solver_interface": backend["interface"],
        "cobamp_gurobi_parameters_requested": solver_parameters.get("requested", {}),
        "cobamp_gurobi_parameters_observed": solver_parameters.get("observed", {}),
        # Backward-compatible alias for readers predating the explicit
        # requested/observed distinction.
        "cobamp_gurobi_parameters": solver_parameters.get("observed", {}),
        "threshold_low": float(low),
        "threshold_high": float(high),
        "imat_properties_overrides": ["exp_vector", "exp_thresholds", "tolerance", "epsilon", "core"],
        "imat_properties_defaults": property_defaults,
        "epsilon": epsilon,
        "activity_tolerance": activity_tolerance,
        "core_reactions": [biomass_id],
        "boundary_core_reactions": [],
        "declared_project_requirements_not_forced": declared_requirements,
        "forced_core_reactions": [biomass_id],
        "contextualization_biomass_forcing": True,
        "biomass_preservation_policy": "biomass_only_core_union_with_hard_parent_bmax_fraction_constraint",
        "parent_bmax": parent_bmax,
        "biomass_floor_fraction": biomass_fraction,
        "biomass_floor": biomass_floor,
        "native_biomass_flux": native_biomass,
        "native_biomass_floor_passed": native_floor_passed,
        "biomass_floor_constraint": floor_constraint,
        "imat_objective_value": imat_objective,
        "raw_retained_biomass": raw_retained_biomass,
        "raw_retained_reaction_count": len(raw_retained),
        "final_retained_reaction_count": len(final_retained),
        "extracted_bmax": extracted_bmax,
        "extracted_bmax_parent_fraction": extracted_bmax / parent_bmax,
        "extracted_biomass_floor_passed": extracted_biomass_floor_passed,
        "canonical_imat_invocations": 1,
        "troppo_solution_status": status,
        "dmi_measurements_used": False,
        "training_samples": list(expression.columns),
        "fit_seed": int(fit_seed),
        "runtime_seconds": elapsed,
        "reactions_retained": len(context_model.reactions),
        "reactions_removed": len(model.reactions) - len(context_model.reactions),
    }
    return ExpertResult(
        context_model,
        "context_specific_model",
        None,
        {
            "retained_indices": final_retained.tolist(),
            "troppo_active_indices": raw_retained.tolist(),
            "core_indices": core.tolist(),
            "final_retained_indices": final_retained.tolist(),
            "raw_retained_biomass": raw_retained_biomass,
            "parent_bmax": parent_bmax,
            "biomass_floor": biomass_floor,
            "native_biomass_flux": native_biomass,
            "imat_objective_value": imat_objective,
            "canonical_imat_invocations": 1,
            "troppo_solution_status": status,
        },
        metadata,
    )
