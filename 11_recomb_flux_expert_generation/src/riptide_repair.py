"""Algorithm-scoped Stage-11 RIPTiDe repair.

The frozen cancer medium is imposed on the parent model before this adapter is
called.  Allowed exchange reactions therefore describe environmental bounds;
they are deliberately *not* converted into RIPTiDe metabolic tasks or exclude
entries.
"""
from __future__ import annotations

import time
from typing import Any, Sequence

import numpy as np
import pandas as pd

from .expert_fitting import ExpertResult
from .models import enforce_gurobi


def _apply_measurement_observability_guard(
    model,
    rcfg: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    """
    Prevent experimentally observed reactions from being structurally pruned.

    This uses only reaction identity/direction, not numerical DMI values.

    Positive:
        v >= +epsilon

    Negative:
        v <= -epsilon
    """

    guard_cfg = rcfg.get("measurement_observability", {})

    if not bool(guard_cfg.get("enabled", False)):
        return {}

    epsilon = float(guard_cfg.get("epsilon", 1.0e-6))

    if not np.isfinite(epsilon) or epsilon <= 0.0:
        raise ValueError(
            "RIPTiDe measurement-observability epsilon "
            f"must be finite and > 0; got {epsilon}"
        )

    reactions = guard_cfg.get(
        "reactions",
        {
            "LDH_L": "negative",
            "PDHm": "positive",
        },
    )

    if not isinstance(reactions, dict) or not reactions:
        raise ValueError(
            "RIPTiDe measurement_observability.reactions "
            "must be a non-empty mapping."
        )

    applied: dict[str, dict[str, Any]] = {}

    for reaction_id, direction in reactions.items():
        reaction_id = str(reaction_id)
        direction = str(direction).lower()

        try:
            reaction = model.reactions.get_by_id(reaction_id)
        except KeyError as exc:
            raise KeyError(
                f"Measurement-observability reaction absent "
                f"from conditioned parent: {reaction_id}"
            ) from exc

        old_lb = float(reaction.lower_bound)
        old_ub = float(reaction.upper_bound)

        if direction == "positive":
            if old_ub < epsilon:
                raise RuntimeError(
                    f"{reaction_id} cannot carry required positive "
                    f"observability flux: bounds={reaction.bounds}, "
                    f"epsilon={epsilon}"
                )

            reaction.lower_bound = max(old_lb, epsilon)

        elif direction == "negative":
            if old_lb > -epsilon:
                raise RuntimeError(
                    f"{reaction_id} cannot carry required negative "
                    f"observability flux: bounds={reaction.bounds}, "
                    f"epsilon={epsilon}"
                )

            reaction.upper_bound = min(old_ub, -epsilon)

        else:
            raise ValueError(
                f"Unknown direction for {reaction_id}: {direction!r}; "
                "expected 'positive' or 'negative'"
            )

        if reaction.lower_bound > reaction.upper_bound:
            raise RuntimeError(
                f"Observability guard made {reaction_id} infeasible: "
                f"{reaction.bounds}"
            )

        applied[reaction_id] = {
            "direction": direction,
            "epsilon": epsilon,
            "original_lower_bound": old_lb,
            "original_upper_bound": old_ub,
            "guarded_lower_bound": float(reaction.lower_bound),
            "guarded_upper_bound": float(reaction.upper_bound),
        }

    # Fail immediately if the very small structural constraints themselves
    # make the parent infeasible.
    solution = model.optimize()

    if solution.status != "optimal":
        raise RuntimeError(
            "Conditioned parent became infeasible after applying "
            "measurement-observability guard."
        )

    return applied


def fit_riptide_expert(
    model,
    expression: pd.DataFrame,
    config: dict[str, Any],
    *,
    fit_seed: int,
    protected_reactions: Sequence[str],
    set_bounds: bool,
) -> ExpertResult:
    """Fit RIPTiDe using frozen parent bounds as the environmental definition."""
    import riptide

    started = time.perf_counter()
    enforce_gurobi(model, seed=fit_seed)
    required = list(dict.fromkeys(map(str, protected_reactions)))
    parent_ids = {reaction.id for reaction in model.reactions}
    missing = set(required) - parent_ids
    if missing:
        raise KeyError(f"Required reactions absent from parent: {sorted(missing)}")

    rcfg = config["contextualization"]["riptide"]
    if rcfg.get("set_bounds") is not True or list(rcfg.get("tasks", [])) or list(rcfg.get("exclude", [])):
        raise RuntimeError("Canonical RIPTiDe requires set_bounds=True, tasks=[], and exclude=[]")
    if set_bounds is not True:
        raise RuntimeError("RIPTiDe adapter received set_bounds other than true")
    objective_id = str(config["objective"]["biomass_reaction"])
    boundary_ids = {reaction.id for reaction in model.boundary}
    nonboundary_required = set(required) - boundary_ids
    unsupported_structural = sorted(nonboundary_required - {"ATPM", objective_id})
    if unsupported_structural:
        raise RuntimeError(
            "RIPTiDe non-task repair does not implement optional structural anchor protection; "
            f"unsupported protected reactions: {unsupported_structural}"
        )
    transcriptome = expression.median(axis=1).to_dict()

    native_input = model.copy()
    native_input.objective = objective_id

    measurement_guard = _apply_measurement_observability_guard(
        native_input,
        rcfg,
    )

    native = riptide.contextualize(
        native_input,
        transcriptome=transcriptome,
        samples=int(rcfg["fit_native_samples"]),
        silent=True,
        prune=bool(rcfg["prune"]),
        fraction=float(rcfg["fraction"]),
        objective=True,
        set_bounds=bool(set_bounds),
        gpr=bool(rcfg["gpr"]),

        # Deliberately remain empty:
        # RIPTiDe tasks use a positive minimum-flux convention,
        # which is unsuitable for the negative LDH_L production direction.
        tasks=[],
        exclude=[],

        task_frac=float(rcfg["task_fraction"]),
    )

    flux = getattr(native, "flux_samples", None)
    frame = None if isinstance(flux, str) or flux is None else pd.DataFrame(flux).astype(np.float64)
    state = {
        key: getattr(native, key, None)
        for key in (
            "transcriptome",
            "minimization_coefficients",
            "maximization_coefficients",
            "pruned",
            "fraction_of_optimum",
            "metabolic_tasks",
            "concordance",
            "gpr_integration",
            "percent_of_mapping",
            "additional_parameters",
            "fraction_bounds",
            "maxfit",
        )
    }
    if not isinstance(state["metabolic_tasks"], (list, tuple, set)):
        state["metabolic_tasks"] = []
    state.update(
        {
            "objective_reaction": objective_id,
            "task_fraction": float(rcfg["task_fraction"]),
            "minimum_flux": 1.0e-8,
            "native_api": "RIPTiDe GapSplit",
            "set_bounds": bool(set_bounds),
            "measurement_observability_guard": measurement_guard,
        }
    )

    result_model = native.model
    for reaction_id, info in measurement_guard.items():

        if reaction_id not in result_model.reactions:
            raise RuntimeError(
                f"RIPTiDe pruned measurement-observable reaction "
                f"{reaction_id} despite observability guard."
            )

        reaction = result_model.reactions.get_by_id(reaction_id)

        epsilon = float(info["epsilon"])
        direction = str(info["direction"])

        if direction == "positive":
            if reaction.upper_bound < epsilon:
                raise RuntimeError(
                    f"{reaction_id} lost positive flux capability "
                    f"after RIPTiDe: bounds={reaction.bounds}"
                )

        elif direction == "negative":
            if reaction.lower_bound > -epsilon:
                raise RuntimeError(
                    f"{reaction_id} lost negative flux capability "
                    f"after RIPTiDe: bounds={reaction.bounds}"
                )
    enforce_gurobi(result_model, seed=fit_seed)
    metadata = {
        "algorithm": "RIPTiDe",
        "threshold_low": np.nan,
        "threshold_high": np.nan,
        "protected_reactions": required,
        "environment_reactions": sorted(set(required) & boundary_ids),
        "structural_protected_reactions": [rid for rid in required if rid not in boundary_ids],
        "training_samples": list(expression.columns),
        "fit_seed": int(fit_seed),
        "dmi_measurements_used": False,
        "measurement_observability_guard_used": bool(measurement_guard),
        "measurement_observability_guard": measurement_guard,
        "measurement_observability_semantics": (
            "reaction_identity_and_direction_only_no_numeric_dmi_values"
        ),
        "protected_as_native_tasks": [],
        "native_exclude": [],
        "environment_handling": "frozen_medium_bounds_only",
        "native_sampling": "GapSplit",
        "set_bounds": bool(set_bounds),
        "contextualization_biomass_forcing": False,
        "gpr": bool(rcfg["gpr"]),
        "fit_native_samples": int(rcfg["fit_native_samples"]),
        "runtime_seconds": time.perf_counter() - started,
        "reactions_retained": len(result_model.reactions),
        "reactions_removed": len(model.reactions) - len(result_model.reactions),
    }
    return ExpertResult(result_model, "riptide_native_expert", frame, state, metadata)
