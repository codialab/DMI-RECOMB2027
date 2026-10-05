"""Isolated Stage-11 repair for Troppo 0.0.7 GIMME.

Troppo/Cobamp's irreversible conversion can reverse the forward/reverse
capacities of asymmetric reversible reactions.  This adapter audits the
mapping against the frozen COBRA parent, normalizes only that representational
defect, validates the biomass-objective mapping independently, and then runs
the published GIMME optimization exactly once.
"""
from __future__ import annotations

import time
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from .core import stable_hash
from .expert_fitting import ExpertResult, _matrix, gimme_expression_vector
from .models import enforce_gurobi, evaluate_biological_objective, reaction_expression


_MAPPING_TOLERANCE = 1.0e-10


def _column_orientation(candidate: np.ndarray, original: np.ndarray) -> int:
    if np.allclose(candidate, original, rtol=0.0, atol=_MAPPING_TOLERANCE):
        return 1
    if np.allclose(candidate, -original, rtol=0.0, atol=_MAPPING_TOLERANCE):
        return -1
    return 0


def _mapped_indices(mapping_value: Any) -> tuple[int, ...]:
    values = mapping_value if isinstance(mapping_value, (tuple, list)) else (mapping_value,)
    array = np.asarray(values)
    if array.ndim != 1 or array.dtype.kind == "b":
        raise RuntimeError("GIMME returned a malformed irreversible reaction mapping")
    try:
        numeric = array.astype(np.float64)
    except (TypeError, ValueError) as exc:
        raise RuntimeError("GIMME returned a non-numeric irreversible reaction mapping") from exc
    if not np.isfinite(numeric).all() or not np.equal(numeric, np.rint(numeric)).all():
        raise RuntimeError("GIMME returned a non-integer-valued irreversible reaction mapping")
    canonical = tuple(int(value) for value in numeric)
    if len(set(canonical)) != len(canonical) or any(value < 0 for value in canonical):
        raise RuntimeError("GIMME returned duplicate or negative irreversible indices")
    return canonical


def _validate_active_indices(values: Any, reaction_count: int) -> np.ndarray:
    array = np.asarray(values)
    if array.ndim != 1 or array.dtype.kind == "b":
        raise RuntimeError("GIMME returned malformed active reaction indices")
    try:
        numeric = array.astype(np.float64)
    except (TypeError, ValueError) as exc:
        raise RuntimeError("GIMME returned non-numeric active reaction indices") from exc
    if not np.isfinite(numeric).all() or not np.equal(numeric, np.rint(numeric)).all():
        raise RuntimeError("GIMME returned non-integer-valued active reaction indices")
    canonical = numeric.astype(np.int64)
    if len(np.unique(canonical)) != len(canonical):
        raise RuntimeError("GIMME returned duplicate active reaction indices")
    if (canonical < 0).any() or (canonical >= int(reaction_count)).any():
        raise RuntimeError("GIMME returned negative or out-of-range active reaction indices")
    return canonical


def audit_and_repair_irreversible_mapping(worker, model) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Audit every original-to-irreversible reaction mapping and repair known swaps."""
    original_S, original_lb, original_ub = _matrix(model)
    irreversible_S = np.asarray(worker.gm.get_stoichiometric_matrix(), dtype=np.float64)
    irreversible_count = irreversible_S.shape[1]
    rows: list[dict[str, Any]] = []
    corrected = 0
    reversible_pairs = 0
    for original_index, reaction in enumerate(model.reactions):
        if original_index not in worker.gm.mapping:
            raise RuntimeError(f"Missing irreversible mapping for {reaction.id}")
        mapped = _mapped_indices(worker.gm.mapping[original_index])
        if any(index >= irreversible_count for index in mapped):
            raise RuntimeError(f"Out-of-range irreversible mapping for {reaction.id}: {mapped}")
        orientations = tuple(
            _column_orientation(irreversible_S[:, index], original_S[:, original_index])
            for index in mapped
        )
        if any(value == 0 for value in orientations):
            raise RuntimeError(f"Unrecognized irreversible column orientation for {reaction.id}")
        pre_bounds = tuple(tuple(map(float, worker.gm.get_reaction_bounds(index))) for index in mapped)
        action = "unchanged"
        forward_index = mapped[orientations.index(1)] if 1 in orientations else None
        reverse_index = mapped[orientations.index(-1)] if -1 in orientations else None
        if len(mapped) == 2:
            reversible_pairs += 1
            if sorted(orientations) != [-1, 1]:
                raise RuntimeError(f"Unrecognized reversible-pair orientation for {reaction.id}")
            expected_forward = max(float(original_ub[original_index]), 0.0)
            expected_reverse = max(-float(original_lb[original_index]), 0.0)
            forward_bounds = tuple(map(float, worker.gm.get_reaction_bounds(forward_index)))
            reverse_bounds = tuple(map(float, worker.gm.get_reaction_bounds(reverse_index)))
            current = (forward_bounds[1], reverse_bounds[1])
            expected = (expected_forward, expected_reverse)
            swapped = (expected_reverse, expected_forward)
            lower_bounds_zero = np.allclose(
                [forward_bounds[0], reverse_bounds[0]], [0.0, 0.0], rtol=0.0, atol=_MAPPING_TOLERANCE
            )
            if lower_bounds_zero and np.allclose(current, expected, rtol=0.0, atol=_MAPPING_TOLERANCE):
                pass
            elif lower_bounds_zero and np.allclose(current, swapped, rtol=0.0, atol=_MAPPING_TOLERANCE):
                worker.gm.set_reaction_bounds(forward_index, lb=0.0, ub=expected_forward, temporary=False)
                worker.gm.set_reaction_bounds(reverse_index, lb=0.0, ub=expected_reverse, temporary=False)
                corrected += 1
                action = "corrected_swapped_capacities"
            else:
                raise RuntimeError(
                    "Unrecognized irreversible bound mapping for "
                    f"{reaction.id}: observed={current}, expected={expected}, swapped={swapped}"
                )
        elif len(mapped) == 1:
            orientation = orientations[0]
            expected_lb = max(float(original_lb[original_index]), 0.0) if orientation == 1 else max(-float(original_ub[original_index]), 0.0)
            expected_ub = max(float(original_ub[original_index]), 0.0) if orientation == 1 else max(-float(original_lb[original_index]), 0.0)
            observed = pre_bounds[0]
            if not np.allclose(observed, (expected_lb, expected_ub), rtol=0.0, atol=_MAPPING_TOLERANCE):
                raise RuntimeError(
                    f"Unrecognized one-way irreversible mapping for {reaction.id}: "
                    f"observed={observed}, expected={(expected_lb, expected_ub)}"
                )
        else:
            raise RuntimeError(f"Unexpected irreversible mapping arity for {reaction.id}: {mapped}")
        post_bounds = tuple(tuple(map(float, worker.gm.get_reaction_bounds(index))) for index in mapped)
        rows.append({
            "original_index": int(original_index),
            "reaction_id": reaction.id,
            "original_lower_bound": float(original_lb[original_index]),
            "original_upper_bound": float(original_ub[original_index]),
            "irreversible_indices": list(mapped),
            "orientations": list(orientations),
            "forward_index": forward_index,
            "reverse_index": reverse_index,
            "pre_bounds": [list(item) for item in pre_bounds],
            "post_bounds": [list(item) for item in post_bounds],
            "action": action,
        })
    return rows, {
        "original_reaction_count": len(model.reactions),
        "irreversible_reaction_count": irreversible_count,
        "reversible_pair_count": reversible_pairs,
        "corrected_asymmetric_reversible_pair_count": corrected,
        "unchanged_mapping_count": len(model.reactions) - corrected,
    }


def validate_biomass_objective_mapping(worker, model, objective_id: str) -> dict[str, Any]:
    """Prove how the original biological objective maps into irreversible space."""
    original_index = model.reactions.index(objective_id)
    mapped = _mapped_indices(worker.gm.mapping[original_index])
    original_S, _, _ = _matrix(model)
    irreversible_S = np.asarray(worker.gm.get_stoichiometric_matrix(), dtype=np.float64)
    orientations = tuple(
        _column_orientation(irreversible_S[:, index], original_S[:, original_index])
        for index in mapped
    )
    if any(value == 0 for value in orientations):
        raise RuntimeError("Biomass objective has an unrecognized irreversible orientation")
    expected_coefficients = {index: float(orientation) for index, orientation in zip(mapped, orientations)}
    adjusted = worker.gm._GIMMEModel__adjust_objective_to_irreversible({original_index: 1.0})
    adjusted = {int(key): float(value) for key, value in adjusted.items()}
    if len(mapped) == 1:
        if orientations != (1,) or adjusted != expected_coefficients or mapped[0] != original_index:
            raise RuntimeError(
                "Biomass objective index is interpreted inconsistently between reversible and irreversible models"
            )
    else:
        # A multi-index mapping is valid only when it is the expected +/- pair
        # for this same reaction and Troppo's coefficient signs preserve net flux.
        if sorted(orientations) != [-1, 1] or adjusted != expected_coefficients:
            raise RuntimeError(
                "Unexplained or semantically inconsistent multi-index biomass objective mapping"
            )
    return {
        "original_reaction_id": objective_id,
        "original_index": int(original_index),
        "irreversible_indices": list(mapped),
        "orientations": list(orientations),
        "expected_irreversible_coefficients": {str(key): value for key, value in expected_coefficients.items()},
        "troppo_adjusted_coefficients": {str(key): value for key, value in adjusted.items()},
        "multi_index_mapping": len(mapped) > 1,
        "validated": True,
    }


def native_biomass_flux(worker, objective_mapping: Mapping[str, Any]) -> float:
    """Reconstruct biological biomass flux from the validated native mapping.

    Troppo's GIMMESolution persists original-reaction values after the raw
    irreversible optimizer has been reverted.  The current biomass objective
    has one validated forward irreversible variable, so that value is exactly
    its biological flux.  A future multi-index objective must expose signed raw
    primals before this adapter can accept it.
    """
    indices = objective_mapping["irreversible_indices"]
    orientations = objective_mapping["orientations"]
    if len(indices) != 1 or list(orientations) != [1]:
        raise RuntimeError(
            "A semantically valid multi-index biomass mapping cannot be reconstructed "
            "from Troppo's collapsed GIMMESolution"
        )
    value = worker.sol.var_values().get(objective_mapping["original_reaction_id"])
    if value is None or not np.isfinite(float(value)):
        raise RuntimeError("Native GIMME biological biomass flux is unavailable")
    return float(value)


def _reduce(model, retained: np.ndarray, required: Sequence[str]):
    keep_ids = {model.reactions[int(index)].id for index in retained}
    keep_ids.update(map(str, required))
    result = model.copy()
    result.remove_reactions(
        [reaction for reaction in result.reactions if reaction.id not in keep_ids],
        remove_orphans=True,
    )
    return result


def fit_gimme_expert(
    model,
    expression: pd.DataFrame,
    config: dict[str, Any],
    *,
    fit_seed: int,
    protected_reactions: Sequence[str],
) -> ExpertResult:
    """Run one repaired GIMME fit and enforce native and reconstructed contracts."""
    from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties

    started = time.perf_counter()
    enforce_gurobi(model, seed=fit_seed)
    required = list(dict.fromkeys(map(str, protected_reactions)))
    reaction_ids = [reaction.id for reaction in model.reactions]
    missing = sorted(set(required) - set(reaction_ids))
    if missing:
        raise KeyError(f"Required structural reactions are absent from the parent: {missing}")
    contextualization = config["contextualization"]
    objective_id = str(config["objective"]["biomass_reaction"])
    tolerance = float(config["validation"].get("objective_tolerance", 1.0e-8))
    direction = str(config["objective"].get("qualification_direction", "max"))
    parent_capacity = float(
        evaluate_biological_objective(model, objective_id, direction=direction, tolerance=tolerance)["capacity"]
    )
    fraction = float(contextualization["objective_fraction"])
    required_minimum = fraction * parent_capacity
    scores = reaction_expression(model, expression, contextualization["pooled_statistic"])
    vector = gimme_expression_vector(scores)
    threshold = float(contextualization["gimme"]["expression_threshold"])
    S, lb, ub = _matrix(model)
    objective_index = reaction_ids.index(objective_id)
    properties = GIMMEProperties(
        vector,
        [{objective_index: 1.0}],
        obj_frac=fraction,
        flux_threshold=threshold,
        solver="GUROBI",
        reaction_ids=reaction_ids,
        metabolite_ids=[metabolite.id for metabolite in model.metabolites],
    )
    worker = GIMME(S, lb, ub, properties)

    mapping_rows, mapping_summary = audit_and_repair_irreversible_mapping(worker, model)
    objective_mapping = validate_biomass_objective_mapping(worker, model, objective_id)
    raw_keep = worker.run()  # Exactly one repaired production GIMME optimization.
    solution = getattr(worker, "sol", None)
    if solution is None:
        raise RuntimeError("GIMME returned active reactions without a native solution")
    try:
        native_status = str(solution.status()).lower()
        penalty_objective = float(solution.objective_value())
    except Exception as exc:
        raise RuntimeError("GIMME native solution status/objective is unavailable") from exc
    if native_status != "optimal" or not np.isfinite(penalty_objective):
        raise RuntimeError(f"GIMME native optimization is not valid: status={native_status}")
    biological_flux = float(native_biomass_flux(worker, objective_mapping))
    native_passed = biological_flux + tolerance >= required_minimum
    if not native_passed:
        raise RuntimeError(
            "Native GIMME biological biomass flux violates the configured parent-objective contract: "
            f"native_biomass_flux={biological_flux:.12g}, required_minimum={required_minimum:.12g}, "
            f"parent_capacity={parent_capacity:.12g}, fraction={fraction:.12g}"
        )

    retained = _validate_active_indices(raw_keep, len(reaction_ids))
    structural = np.asarray([reaction_ids.index(rid) for rid in required], dtype=np.int64)
    union = np.union1d(retained, structural).astype(np.int64)
    context_model = _reduce(model, union, required)
    enforce_gurobi(context_model, seed=fit_seed)
    contextual_capacity = float(
        evaluate_biological_objective(
            context_model, objective_id, direction=direction, tolerance=tolerance
        )["capacity"]
    )
    reconstructed_passed = contextual_capacity + tolerance >= required_minimum
    if not reconstructed_passed:
        native_values = solution.var_values()
        native_support = {
            rid for rid, value in native_values.items() if abs(float(value)) > tolerance
        }
        missing_support = sorted(native_support - {reaction.id for reaction in context_model.reactions})
        raise RuntimeError(
            "Reconstructed GIMME model violates the configured parent-objective contract: "
            f"contextual_capacity={contextual_capacity:.12g}, required_minimum={required_minimum:.12g}, "
            f"native_flux_support_lost={missing_support[:50]}"
        )

    audit_hash = stable_hash(mapping_rows)
    diagnostics = {
        "parent_capacity": parent_capacity,
        "objective_fraction": fraction,
        "required_minimum": required_minimum,
        "native_solver_status": native_status,
        "native_biomass_flux": biological_flux,
        "gimme_optimization_objective_value": penalty_objective,
        "native_objective_contract_passed": native_passed,
        "retained_reaction_count": int(len(union)),
        "retained_reaction_indices_sha256": stable_hash(union.tolist()),
        "corrected_asymmetric_reversible_pair_count": mapping_summary[
            "corrected_asymmetric_reversible_pair_count"
        ],
        "reconstructed_contextual_capacity": contextual_capacity,
        "reconstructed_objective_contract_passed": reconstructed_passed,
        "mapping_audit_sha256": audit_hash,
        "mapping_audit_relative_path": "qc/gimme_irreversible_mapping.tsv",
        "mapping_audit_summary": mapping_summary,
        "biomass_objective_mapping": objective_mapping,
        "gimme_optimization_runs": 1,
    }
    return ExpertResult(
        context_model,
        "context_specific_model",
        None,
        {
            "retained_indices": union.tolist(),
            "_dedicated_gimme_mapping_audit": mapping_rows,
            **diagnostics,
        },
        {
            "algorithm": "GIMME",
            "implementation": "Troppo 0.0.7 GIMME with audited irreversible-capacity compatibility repair",
            "gimme_input_semantics": "native_reaction_expression",
            "gimme_expression_threshold": threshold,
            "numerical_epsilon": float(contextualization["epsilon"]),
            "fit_objective_fraction": fraction,
            "dmi_measurements_used": False,
            "training_samples": list(expression.columns),
            "fit_seed": int(fit_seed),
            "runtime_seconds": time.perf_counter() - started,
            "reactions_retained": len(context_model.reactions),
            "reactions_removed": len(model.reactions) - len(context_model.reactions),
            **diagnostics,
        },
    )
