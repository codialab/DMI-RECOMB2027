"""Parent-model handling reused from Stage-10 ``src/models.py``."""
from __future__ import annotations
import ast
import io
import lzma
from pathlib import Path
from typing import Any, Sequence
import cobra
import numpy as np
import pandas as pd
from .boundary_qc import (active_boundary_source_fluxes, assert_resolved_source_routes,
    assert_source_routes_within_resolved_medium, boundary_source_direction,
    boundary_source_routes, close_boundary_source_direction,
    close_boundary_source_direction_for_negative_control, unexpected_boundary_sources)
from .core import evaluate_gpr, stable_hash

PRIMARY_EXPERT_VARIANT = "rna_only_primary"
RIPTIDE_CLOSED_SOURCE_QC_POLICY_NAME = "riptide_set_bounds_closed_source_two_stage"
RIPTIDE_CLOSED_SOURCE_QC_POLICY_VERSION = 1


class RiptideClosedSourceQCError(RuntimeError):
    """Fail-closed RIPTiDe qualification error with machine-readable evidence."""

    def __init__(self, message: str, diagnostics: dict[str, Any]):
        super().__init__(message)
        self.diagnostics = diagnostics


ANCHOR_ABLATION_VARIANT = "rna_with_anchor_protection_ablation"
TECHNICAL_REACTIONS = ("ATPM", "BIOMASS_reaction")
HISTORICAL_DMI_ANCHORS = ("HEX1", "r0354", "r0355", "SBTR", "CBPPer", "UGALGTg", "LDH_L")
DMI_OBSERVABLE_GROUPS = {
    "glucose": ("EX_glc__D_e", "HEX1", "r0354", "r0355", "SBTR", "CBPPer", "UGALGTg"),
    "lactate": ("L_LACt2r", "EX_lac__L_e", "LDH_L"),
    "glx": ("PYRt2m", "PDHm", "PCm", "AKGDm"),
}


def enforce_gurobi(model, *, seed: int | None = None) -> dict[str, Any]:
    model.solver = "gurobi"
    if "gurobi" not in model.solver.interface.__name__.lower(): raise RuntimeError("COBRApy did not activate Gurobi")
    model.solver.configuration.verbosity = 0; model.solver.configuration.timeout = None
    params = model.solver.problem.Params; params.Threads = 1; params.OutputFlag = 0
    if seed is not None: params.Seed = int(seed)
    return {"interface": model.solver.interface.__name__, "Threads": int(params.Threads), "Seed": int(params.Seed),
            "MIPGap": float(params.MIPGap), "FeasibilityTol": float(params.FeasibilityTol),
            "OptimalityTol": float(params.OptimalityTol), "TimeLimit": float(params.TimeLimit)}


def load_model(path: Path, *, seed: int | None = None):
    with lzma.open(path, "rt") as handle: model = cobra.io.read_sbml_model(io.StringIO(handle.read()))
    enforce_gurobi(model, seed=seed); return model


def model_hash(model) -> str:
    objective = {r.id: float(v) for r, v in cobra.util.solver.linear_reaction_coefficients(model).items()}
    metabolite_constraints = {m.id for m in model.metabolites}
    extras = sorted((c.name, float(c.lb) if c.lb is not None else None, float(c.ub) if c.ub is not None else None, str(c.expression))
                    for c in model.solver.constraints if c.name not in metabolite_constraints)
    return stable_hash({"reactions": [{"id": r.id, "lb": float(r.lower_bound), "ub": float(r.upper_bound), "gpr": r.gene_reaction_rule,
        "stoichiometry": sorted((m.id, float(v)) for m, v in r.metabolites.items())} for r in model.reactions], "objective": objective,
        "extra_constraints": extras})


def frozen_medium_reactions(final_a: dict[str, Any]) -> list[str]:
    medium = final_a["medium"]
    return sorted(set(map(str, medium["uptake_bounds"])) | set(map(str, medium.get("secretion_only", []))))


def resolved_frozen_medium_routes(model, final_a: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Resolve configured medium magnitudes into stoichiometry-aware source routes."""
    routes: dict[str, dict[str, Any]] = {}
    for rid, item in final_a["medium"]["uptake_bounds"].items():
        if rid not in model.reactions:
            raise KeyError(f"Frozen medium reaction absent: {rid}")
        reaction = model.reactions.get_by_id(rid)
        if not reaction.boundary:
            raise RuntimeError(f"Frozen medium reaction is not a boundary reaction: {rid}")
        capacity = float(item["magnitude"])
        if not np.isfinite(capacity) or capacity <= 0.0:
            raise ValueError(f"Frozen medium source capacity must be finite and positive: {rid}")
        routes[str(rid)] = {
            "source_direction": boundary_source_direction(reaction),
            "source_capacity": capacity,
        }
    return dict(sorted(routes.items()))


def requirement_policy_for_variant(
    final_a: dict[str, Any], variant: str, *, objective_id: str = "BIOMASS_reaction"
) -> dict[str, list[dict[str, Any]]]:
    """Build the role-tagged RNA requirement policy with explicit provenance."""
    if variant not in {PRIMARY_EXPERT_VARIANT, ANCHOR_ABLATION_VARIANT}:
        raise ValueError(f"Unknown expert variant: {variant}")
    environment = [
        {"reaction_id": str(rid), "role": "environment_reactions",
         "provenance_source": "resolved_frozen_medium", "condition_type": "rna_only",
         "required_bound_or_constraint": dict(item)}
        for rid, item in final_a["medium"]["uptake_bounds"].items()
    ]
    environment.extend(
        {"reaction_id": str(rid), "role": "environment_reactions",
         "provenance_source": "resolved_frozen_medium", "condition_type": "rna_only",
         "required_bound_or_constraint": {"secretion_only": True}}
        for rid in final_a["medium"].get("secretion_only", [])
    )
    policy: dict[str, list[dict[str, Any]]] = {
        "environment_reactions": environment,
        "structural_required_reactions": [
            {"reaction_id": "ATPM", "role": "structural_required_reactions",
             "provenance_source": "stage11_technical_infrastructure", "condition_type": "rna_only",
             "required_bound_or_constraint": None}
        ],
        "measurement_required_reactions": [],
        "objective_reactions": [
            {"reaction_id": str(objective_id), "role": "objective_reactions",
             "provenance_source": "master_config.objective", "condition_type": "rna_only",
             "required_bound_or_constraint": {"direction": "max"}}
        ],
        "algorithm_specific_requirements": [],
    }
    if variant == ANCHOR_ABLATION_VARIANT:
        policy["measurement_required_reactions"] = [
            {"reaction_id": rid, "role": "measurement_required_reactions",
             "provenance_source": "historical_anchor_ablation", "condition_type": "rna_only_anchor_ablation",
             "required_bound_or_constraint": None}
            for rid in HISTORICAL_DMI_ANCHORS
        ]
    return policy


def requirement_ids(policy: dict[str, Sequence[dict[str, Any]]], roles: Sequence[str]) -> list[str]:
    """Resolve role entries, rejecting unclassified reactions rather than guessing."""
    values: list[str] = []
    for role in roles:
        if role not in policy:
            raise ValueError(f"Requirement policy lacks role: {role}")
        for entry in policy[role]:
            if entry.get("role") != role or not entry.get("reaction_id") or not entry.get("provenance_source"):
                raise ValueError(f"Unclassified or malformed requirement entry for role {role}: {entry}")
            values.append(str(entry["reaction_id"]))
    return list(dict.fromkeys(values))


def algorithm_required_reactions(policy: dict[str, Sequence[dict[str, Any]]]) -> list[str]:
    """Roles eligible for structural retention; environment is bounds-only."""
    return requirement_ids(policy, (
        "structural_required_reactions", "measurement_required_reactions",
        "objective_reactions", "algorithm_specific_requirements",
    ))


def protected_reactions_for_variant(final_a: dict[str, Any], variant: str) -> list[str]:
    # Compatibility view only. Adapter behavior must be derived from roles.
    policy = requirement_policy_for_variant(final_a, variant)
    return sorted(requirement_ids(policy, tuple(policy)))


def postfit_required_reactions(parent, protected: Sequence[str], algorithm: str, objective_id: str) -> list[str]:
    """Resolve which protected IDs must remain structurally present after fitting.

    The frozen-medium boundary set constrains the parent before contextualization,
    but extracting methods may legitimately prune an allowed exchange that carries
    zero flux.  Canonical iMAT and CORDA do not structurally protect project
    requirement roles; biomass viability is a downstream qualification result.
    RIPTiDe uses the frozen bounds as environment and requires its biomass objective
    to survive; medium exchanges are not converted into tasks.  Other adapters keep
    the legacy full protected-set retention contract.
    """
    values = list(dict.fromkeys(map(str, protected)))
    ids = {reaction.id for reaction in parent.reactions}
    missing = set(values) - ids
    if missing:
        raise KeyError(f"Protected reactions absent from parent: {sorted(missing)}")
    boundary = {reaction.id for reaction in parent.boundary}
    if algorithm in {"iMAT", "CORDA"}:
        # Canonical iMAT and CORDA decide topology from expression confidence.
        # Do not turn project requirement roles into structural or post-fit rescue.
        # Biomass viability is evaluated independently after contextualization.
        return []
    if algorithm == "RIPTiDe":
        if objective_id not in ids:
            raise KeyError(f"RIPTiDe objective absent from parent: {objective_id}")
        return [str(objective_id)]
    return values


def validate_environmental_protection(model, final_a: dict[str, Any], protected: Sequence[str]) -> None:
    ids = {reaction.id for reaction in model.reactions}; missing = set(protected) - ids
    if missing: raise KeyError(f"Protected reactions absent from parent: {sorted(missing)}")
    boundary = {reaction.id for reaction in model.boundary}
    nonboundary = set(frozen_medium_reactions(final_a)) - boundary
    if nonboundary: raise RuntimeError(f"Frozen-medium reactions are not parent boundary reactions: {sorted(nonboundary)}")


def apply_frozen_medium(model, final_a: dict[str, Any]) -> None:
    medium = final_a["medium"]
    expected = resolved_frozen_medium_routes(model, final_a)
    configured = set(expected)
    secretion_only = set(map(str, medium.get("secretion_only", [])))
    overlap = configured & secretion_only
    if overlap:
        raise ValueError(f"Frozen medium cannot declare source and secretion-only: {sorted(overlap)}")
    for rid in secretion_only:
        if rid not in model.reactions:
            raise KeyError(f"Frozen secretion reaction absent: {rid}")
        if not model.reactions.get_by_id(rid).boundary:
            raise RuntimeError(f"Frozen secretion reaction is not a boundary reaction: {rid}")
    # Reuse the Analysis-A invariant over every boundary reaction, not exchanges
    # only. Closing the source sign deliberately leaves legitimate drain flux.
    for reaction in model.boundary:
        close_boundary_source_direction(reaction)
    for rid, route in expected.items():
        reaction = model.reactions.get_by_id(rid)
        if route["source_direction"] == "negative":
            reaction.lower_bound = -float(route["source_capacity"])
        else:
            reaction.upper_bound = float(route["source_capacity"])
    assert_resolved_source_routes(model, expected, context="Frozen medium")


def load_expression(path: Path, training_samples: Sequence[str]) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype={"model_gene_id": str}); missing = set(training_samples) - set(frame.columns)
    if missing: raise KeyError(f"Expression samples absent: {sorted(missing)}")
    result = frame.set_index("model_gene_id")[list(training_samples)].astype(np.float64)
    if not np.isfinite(result.to_numpy()).all() or (result < 0).any().any(): raise ValueError("Expression must be finite and nonnegative")
    return result


def reaction_expression(model, expression: pd.DataFrame, pooled_statistic: str = "median", *, and_operator: str = "MIN", or_operator: str = "MAX") -> pd.Series:
    pooled = expression.median(axis=1) if pooled_statistic == "median" else expression.mean(axis=1)
    return pd.Series({r.id: evaluate_gpr(r.gpr.body, pooled.to_dict(), and_operator, or_operator) for r in model.reactions}, dtype=np.float64)


def reaction_expression_with_evidence(model, expression: pd.DataFrame, pooled_statistic: str = "median", *,
                                      and_operator: str = "MIN", or_operator: str = "SUM") -> tuple[pd.Series, dict[str, str]]:
    """Evaluate GPRs without silently dropping unmapped genes from OR expressions."""
    pooled = expression.median(axis=1) if pooled_statistic == "median" else expression.mean(axis=1)
    values = pooled.to_dict()

    def combine(numbers: list[float], operator: str) -> float:
        array = np.asarray(numbers, dtype=float)
        if not np.isfinite(array).all():
            return np.nan
        if operator == "MIN": return float(array.min())
        if operator == "MAX": return float(array.max())
        if operator == "SUM": return float(array.sum())
        raise ValueError(operator)

    def one(rule) -> tuple[float, str]:
        node = rule.body if hasattr(rule, "body") else rule
        if node is None:
            return np.nan, "unmapped"

        def visit(item: ast.AST) -> tuple[float, int, int]:
            if isinstance(item, ast.Name):
                value = values.get(item.id, np.nan)
                return float(value), int(np.isfinite(value)), 1
            if isinstance(item, ast.BoolOp):
                parts = [visit(child) for child in item.values]
                operator = and_operator if isinstance(item.op, ast.And) else or_operator
                return combine([part[0] for part in parts], operator), sum(part[1] for part in parts), sum(part[2] for part in parts)
            raise TypeError(f"Unsupported GPR syntax: {type(item).__name__}")

        value, mapped, total = visit(node)
        if total == 0 or mapped == 0: return np.nan, "unmapped"
        if mapped < total: return np.nan, "partial_mapping"
        return value, "fully_mapped_zero" if value == 0.0 else "fully_mapped_positive"

    evaluated = {reaction.id: one(reaction.gpr.body) for reaction in model.reactions}
    return (pd.Series({rid: value for rid, (value, _) in evaluated.items()}, dtype=np.float64),
            {rid: status for rid, (_, status) in evaluated.items()})


def expression_thresholds(scores: pd.Series, lower: float, upper: float) -> tuple[float, float]:
    finite = scores[np.isfinite(scores) & (scores > 0)]
    if finite.empty: raise ValueError("No positive reaction expression in training data")
    lo, hi = finite.quantile([lower, upper]).tolist(); return float(lo), float(hi)


def evaluate_biological_objective(
    model,
    objective_id: str,
    *,
    direction: str = "max",
    tolerance: float = 1.0e-8,
) -> dict[str, Any]:
    """Independently qualify the experimental phenotype before a floor is made."""
    if direction != "max":
        raise ValueError(f"Unsupported biological objective direction: {direction}")
    if objective_id not in model.reactions:
        raise KeyError(f"Biological objective reaction absent: {objective_id}")
    with model:
        model.objective = objective_id
        model.objective_direction = direction
        solution = model.optimize()
        capacity = float(solution.objective_value) if solution.objective_value is not None else np.nan
    if solution.status != "optimal":
        raise RuntimeError(f"Biological objective is not feasible: {solution.status}")
    if not np.isfinite(capacity):
        raise RuntimeError(f"Biological objective capacity is non-finite: {capacity}")
    if capacity <= float(tolerance):
        raise RuntimeError(f"Biological objective capacity is not positive above tolerance: {capacity}")
    return {"status": solution.status, "objective": objective_id, "direction": direction,
            "capacity": capacity, "tolerance": float(tolerance)}


def closed_source_negative_control(
    model,
    objective_id: str,
    *,
    tolerance: float = 1.0e-8,
    direction: str = "max",
) -> dict[str, Any]:
    """Require objective capacity to collapse after every source sign is closed."""
    if direction != "max":
        raise ValueError(f"Unsupported biological objective direction: {direction}")
    control = model.copy()
    for reaction in control.boundary:
        close_boundary_source_direction_for_negative_control(reaction)
    control.objective = objective_id
    control.objective_direction = direction
    solution = control.optimize()
    residual = float(solution.objective_value) if solution.objective_value is not None else np.nan
    result = {
        "status": solution.status,
        "objective": objective_id,
        "residual_biomass": residual,
        "tolerance": float(tolerance),
        "passed": bool(solution.status == "optimal" and np.isfinite(residual) and abs(residual) <= float(tolerance)),
    }
    if not result["passed"]:
        raise RuntimeError(
            "Closed-source negative control retained objective capacity: "
            f"status={solution.status}, residual={residual}, tolerance={tolerance}"
        )
    return result


def riptide_closed_source_negative_control(
    model,
    conditioned_parent,
    objective_id: str,
    *,
    tolerance: float = 1.0e-8,
    direction: str = "max",
) -> dict[str, Any]:
    """RIPTiDe-specific source-leak control for canonical ``set_bounds=True`` models.

    The ordinary closed-source control remains authoritative when it is feasible.
    If closing all source directions makes the RIPTiDe expert infeasible, a
    disposable second control relaxes *only* compulsory-flux bounds introduced
    relative to the conditioned parent: bounds that exclude zero even though the
    corresponding parent reaction admitted zero.  The relaxation extends such a
    bound only to zero, never beyond zero, so it cannot open a nutrient source or
    add a new reaction direction.  Biomass must then optimize to zero.

    This is QC-only.  Neither the persisted expert nor the conditioned parent is
    mutated.
    """
    if direction != "max":
        raise ValueError(f"Unsupported biological objective direction: {direction}")
    if objective_id not in model.reactions:
        raise KeyError(f"Biological objective reaction absent: {objective_id}")

    diagnostics: dict[str, Any] = {
        "policy_name": RIPTIDE_CLOSED_SOURCE_QC_POLICY_NAME,
        "policy_version": RIPTIDE_CLOSED_SOURCE_QC_POLICY_VERSION,
        "objective": objective_id,
        "tolerance": float(tolerance),
        "strict_status": None,
        "strict_residual_biomass": None,
        "strict_source_routes_after_closure": {},
        "fallback_invoked": False,
        "fallback_solve_performed": False,
        "fallback_status": None,
        "fallback_residual_biomass": None,
        "qualifying_bound_count": 0,
        "qualifying_reaction_count": 0,
        "nonqualifying_compulsory_bound_count": 0,
        "relaxed_riptide_compulsory_bound_count": 0,
        "relaxed_riptide_compulsory_reaction_count": 0,
        "relaxed_riptide_compulsory_reaction_ids": [],
        "relaxed_riptide_compulsory_bounds": [],
        "nonqualifying_compulsory_bounds": [],
        "remaining_forced_nonzero_bounds": [],
        "remaining_forced_nonzero_bound_count": 0,
        "remaining_forced_nonzero_reaction_count": 0,
        "missing_conditioned_parent_reaction_ids": [],
        "source_routes_after_fallback": None,
        "open_source_route_count_after_fallback": None,
        "no_source_routes_open_after_fallback": None,
        "passed": False,
        "failure_reason": None,
        # Compatibility keys retained for current qualification readers.
        "status": None,
        "residual_biomass": None,
        "fallback_applied": False,
        "policy": RIPTIDE_CLOSED_SOURCE_QC_POLICY_NAME,
    }

    def fail(reason: str, message: str) -> None:
        diagnostics["failure_reason"] = reason
        raise RiptideClosedSourceQCError(message, dict(diagnostics))

    def finite_numeric_or_none(value) -> float | None:
        if value is None:
            return None
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            return None
        return numeric if np.isfinite(numeric) else None

    def optimize_status_and_objective(candidate) -> tuple[str, float | None]:
        """Solve status-first so infeasible controls never request primal fluxes.

        COBRApy's high-level ``Model.optimize()`` constructs a full ``Solution``
        and therefore requests primal reaction values.  Gurobi has no primal
        ``X`` values after an infeasible solve, which can raise before this QC
        can inspect the expected ``infeasible`` status and enter Stage B.  The
        optlang solver API exposes status without constructing fluxes, so use it
        when available and read the objective only after an optimal solve.
        """
        solver = getattr(candidate, "solver", None)
        solver_optimize = getattr(solver, "optimize", None)
        if callable(solver_optimize):
            raw_status = solver_optimize()
            status = str(
                raw_status if raw_status is not None else getattr(solver, "status", "unknown")
            ).strip().lower()
            if status != "optimal":
                return status, None
            objective = getattr(solver, "objective", None)
            return status, finite_numeric_or_none(getattr(objective, "value", None))

        # Test doubles and non-COBRA callers may expose only the high-level API.
        solution = candidate.optimize()
        status = str(getattr(solution, "status", "unknown")).strip().lower()
        return status, finite_numeric_or_none(getattr(solution, "objective_value", None))

    control = model.copy()
    for reaction in control.boundary:
        close_boundary_source_direction_for_negative_control(reaction)
    control.objective = objective_id
    control.objective_direction = direction

    strict_status, strict_residual = optimize_status_and_objective(control)
    diagnostics.update({
        "strict_status": strict_status,
        "strict_residual_biomass": strict_residual,
        "strict_source_routes_after_closure": boundary_source_routes(control, tolerance=0.0),
    })
    if strict_status == "optimal":
        diagnostics.update({"status": strict_status, "residual_biomass": strict_residual})
        if strict_residual is None:
            fail(
                "strict_optimal_nonfinite_objective",
                "RIPTiDe strict closed-source solve reported optimal with a non-finite biomass objective",
            )
        if abs(strict_residual) > float(tolerance):
            fail(
                "strict_source_independent_biomass_capacity",
                "RIPTiDe strict closed-source negative control retained biomass capacity: "
                f"status={strict_status}, residual={strict_residual}, tolerance={tolerance}",
            )
        diagnostics["passed"] = True
        return diagnostics

    if strict_status != "infeasible":
        diagnostics.update({"status": strict_status, "residual_biomass": strict_residual})
        fail(
            "strict_status_not_optimal_or_infeasible",
            "RIPTiDe strict closed-source negative control failed before fallback: "
            f"status={strict_status}, residual={strict_residual}, tolerance={tolerance}",
        )

    diagnostics["fallback_invoked"] = True
    relaxed: list[dict[str, Any]] = []
    nonqualifying: list[dict[str, Any]] = []
    parent_ids = {reaction.id for reaction in conditioned_parent.reactions}
    missing_from_parent = sorted(
        reaction.id for reaction in control.reactions if reaction.id not in parent_ids
    )
    diagnostics["missing_conditioned_parent_reaction_ids"] = missing_from_parent
    if missing_from_parent:
        fail(
            "missing_conditioned_parent_reactions",
            "RIPTiDe fallback contains reactions absent from the conditioned parent: "
            f"{missing_from_parent}",
        )

    for reaction in control.reactions:
        parent_reaction = conditioned_parent.reactions.get_by_id(reaction.id)
        lower, upper = map(float, reaction.bounds)
        parent_lower, parent_upper = map(float, parent_reaction.bounds)
        new_lower, new_upper = lower, upper
        relaxed_sides: list[str] = []
        if lower > 0.0 and parent_lower <= 0.0:
            new_lower = 0.0
            relaxed_sides.append("lower_bound")
        elif lower > 0.0:
            nonqualifying.append({
                "reaction_id": reaction.id,
                "bound_side": "lower_bound",
                "strict_closed_source_expert_bound": lower,
                "conditioned_parent_bound": parent_lower,
                "reason": "conditioned_parent_also_excluded_zero",
            })
        if upper < 0.0 and parent_upper >= 0.0:
            new_upper = 0.0
            relaxed_sides.append("upper_bound")
        elif upper < 0.0:
            nonqualifying.append({
                "reaction_id": reaction.id,
                "bound_side": "upper_bound",
                "strict_closed_source_expert_bound": upper,
                "conditioned_parent_bound": parent_upper,
                "reason": "conditioned_parent_also_excluded_zero",
            })
        if relaxed_sides:
            reaction.bounds = (new_lower, new_upper)
            relaxed.append({
                "reaction": reaction.id,
                "reaction_id": reaction.id,
                "strict_closed_source_bounds": [lower, upper],
                "conditioned_parent_bounds": [parent_lower, parent_upper],
                "fallback_bounds": [new_lower, new_upper],
                "relaxed_bound_sides": relaxed_sides,
            })

    relaxed_ids = [item["reaction_id"] for item in relaxed]
    relaxed_bound_count = sum(len(item["relaxed_bound_sides"]) for item in relaxed)
    diagnostics.update({
        "qualifying_bound_count": relaxed_bound_count,
        "qualifying_reaction_count": len(relaxed),
        "nonqualifying_compulsory_bound_count": len(nonqualifying),
        "relaxed_riptide_compulsory_bound_count": relaxed_bound_count,
        "relaxed_riptide_compulsory_reaction_count": len(relaxed),
        "relaxed_riptide_compulsory_reaction_ids": relaxed_ids,
        "relaxed_riptide_compulsory_bounds": relaxed,
        "nonqualifying_compulsory_bounds": nonqualifying,
        "fallback_applied": bool(relaxed),
    })
    remaining_forced = []
    for reaction in control.reactions:
        forced_sides = []
        if float(reaction.lower_bound) > 0.0:
            forced_sides.append("lower_bound")
        if float(reaction.upper_bound) < 0.0:
            forced_sides.append("upper_bound")
        if forced_sides:
            remaining_forced.append({
                "reaction_id": reaction.id,
                "bounds": list(map(float, reaction.bounds)),
                "forced_bound_sides": forced_sides,
            })
    diagnostics.update({
        "remaining_forced_nonzero_bounds": remaining_forced,
        "remaining_forced_nonzero_bound_count": sum(
            len(item["forced_bound_sides"]) for item in remaining_forced
        ),
        "remaining_forced_nonzero_reaction_count": len(remaining_forced),
    })
    if not relaxed:
        fail(
            "no_riptide_introduced_compulsory_bounds",
            "RIPTiDe closed-source model is infeasible but no RIPTiDe-introduced compulsory-flux "
            "bound excluding parent-allowed zero was found; refusing to reinterpret infeasibility as a pass",
        )

    reopened_sources = boundary_source_routes(control, tolerance=0.0)
    diagnostics.update({
        "source_routes_after_fallback": reopened_sources,
        "open_source_route_count_after_fallback": len(reopened_sources),
        "no_source_routes_open_after_fallback": not bool(reopened_sources),
    })
    if reopened_sources:
        fail(
            "fallback_reopened_source_routes",
            "RIPTiDe negative-control fallback reopened source routes: "
            f"{sorted(reopened_sources)}",
        )

    fallback_status, fallback_residual = optimize_status_and_objective(control)
    diagnostics.update({
        "fallback_solve_performed": True,
        "fallback_status": fallback_status,
        "fallback_residual_biomass": fallback_residual,
        "status": fallback_status,
        "residual_biomass": fallback_residual,
    })
    passed = bool(
        fallback_status == "optimal"
        and fallback_residual is not None
        and abs(fallback_residual) <= float(tolerance)
    )
    if not passed:
        fail(
            "fallback_did_not_prove_zero_biomass_capacity",
            "RIPTiDe closed-source negative-control fallback did not prove zero biomass capacity: "
            f"strict_status={strict_status}, fallback_status={fallback_status}, "
            f"fallback_residual={fallback_residual}, tolerance={tolerance}, "
            f"relaxed_bounds={relaxed_bound_count}",
        )
    diagnostics["passed"] = True
    return diagnostics


def qualify_context_model(
    model,
    objective_id: str,
    required: Sequence[str],
    expected_source_routes: dict[str, dict[str, Any]],
    validation: dict[str, Any],
    *,
    objective_direction: str = "max",
) -> dict[str, Any]:
    """Independently admit a contextual model for persistence and sampling."""
    objective_tolerance = float(validation.get("objective_tolerance", 1.0e-8))
    assert_source_routes_within_resolved_medium(
        model, expected_source_routes, context="Contextual-model qualification",
        tolerance=float(validation["bound_tolerance"]),
    )
    qualified = validate_model(
        model, objective_id, required, objective_tolerance=objective_tolerance,
        objective_direction=objective_direction,
    )
    with model:
        model.objective = objective_id
        model.objective_direction = objective_direction
        solution = model.optimize()
    if solution.status != "optimal":
        raise RuntimeError(f"Qualification representative solve failed: {solution.status}")
    fluxes = solution.fluxes.reindex([reaction.id for reaction in model.reactions])
    from .core import validate_samples
    numerical = validate_samples(pd.DataFrame([fluxes.to_numpy()], columns=fluxes.index), model, validation)
    unexpected = unexpected_boundary_sources(
        model, expected_source_routes, tolerance=float(validation["zero_tolerance"])
    )
    if unexpected:
        raise RuntimeError(f"Contextual model contains undeclared sources: {unexpected}")
    negative_control = closed_source_negative_control(
        model, objective_id, tolerance=objective_tolerance, direction=objective_direction
    )
    return {
        "passed": True,
        "model_hash": model_hash(model),
        "reaction_count": len(model.reactions),
        "objective": qualified,
        "actual_source_routes": boundary_source_routes(model),
        "undeclared_source_count": len(unexpected),
        "active_source_fluxes": active_boundary_source_fluxes(
            model, solution.fluxes, tolerance=float(validation["zero_tolerance"])
        ),
        "required_reactions": list(map(str, required)),
        "structural_reactions_retained": True,
        "mass_balance": numerical,
        "bounds_passed": True,
        "closed_source_negative_control": negative_control,
    }


def qualify_riptide_context_model(
    model,
    conditioned_parent,
    objective_id: str,
    required: Sequence[str],
    expected_source_routes: dict[str, dict[str, Any]],
    validation: dict[str, Any],
    *,
    objective_direction: str = "max",
) -> dict[str, Any]:
    """Qualify canonical RIPTiDe while preserving generic qualification semantics.

    All ordinary structural, environmental, numerical, and objective checks are
    identical to :func:`qualify_context_model`; only the closed-source negative
    control is RIPTiDe-aware.
    """
    objective_tolerance = float(validation.get("objective_tolerance", 1.0e-8))
    assert_source_routes_within_resolved_medium(
        model, expected_source_routes, context="RIPTiDe contextual-model qualification",
        tolerance=float(validation["bound_tolerance"]),
    )
    qualified = validate_model(
        model, objective_id, required, objective_tolerance=objective_tolerance,
        objective_direction=objective_direction,
    )
    with model:
        model.objective = objective_id
        model.objective_direction = objective_direction
        solution = model.optimize()
    if solution.status != "optimal":
        raise RuntimeError(f"RIPTiDe qualification representative solve failed: {solution.status}")
    fluxes = solution.fluxes.reindex([reaction.id for reaction in model.reactions])
    from .core import validate_samples
    numerical = validate_samples(pd.DataFrame([fluxes.to_numpy()], columns=fluxes.index), model, validation)
    unexpected = unexpected_boundary_sources(
        model, expected_source_routes, tolerance=float(validation["zero_tolerance"])
    )
    if unexpected:
        raise RuntimeError(f"RIPTiDe contextual model contains undeclared sources: {unexpected}")
    negative_control = riptide_closed_source_negative_control(
        model, conditioned_parent, objective_id, tolerance=objective_tolerance,
        direction=objective_direction,
    )
    return {
        "passed": True,
        "model_hash": model_hash(model),
        "reaction_count": len(model.reactions),
        "reaction_order_sha256": stable_hash([reaction.id for reaction in model.reactions]),
        "bounds_sha256": stable_hash([
            [reaction.id, float(reaction.lower_bound), float(reaction.upper_bound)]
            for reaction in model.reactions
        ]),
        "objective": qualified,
        "actual_source_routes": boundary_source_routes(model),
        "undeclared_source_count": len(unexpected),
        "active_source_fluxes": active_boundary_source_fluxes(
            model, solution.fluxes, tolerance=float(validation["zero_tolerance"])
        ),
        "required_reactions": list(map(str, required)),
        "structural_reactions_retained": True,
        "mass_balance": numerical,
        "bounds_passed": True,
        "closed_source_negative_control": negative_control,
    }



def apply_biomass_floor(model, biomass_id: str, fraction: float, *, tolerance: float = 1.0e-8,
                        direction: str = "max") -> dict[str, float]:
    if not 0.0 < float(fraction) <= 1.0:
        raise ValueError("Biomass floor fraction must be in (0, 1]")
    qualified = evaluate_biological_objective(model, biomass_id, tolerance=tolerance, direction=direction)
    reaction = model.reactions.get_by_id(biomass_id); capacity = float(qualified["capacity"])
    previous_lower_bound = float(reaction.lower_bound)
    requested_floor = float(fraction) * capacity
    reaction.lower_bound = max(previous_lower_bound, requested_floor)
    effective_floor = float(reaction.lower_bound)
    return {
        # Compatibility keys retained for current readers.
        "biomass_capacity": capacity,
        "biomass_floor": effective_floor,
        "fraction": float(fraction),
        # Explicit requested-versus-effective sampling-floor provenance.
        "contextualized_biomass_maximum": capacity,
        "requested_common_sampling_fraction": float(fraction),
        "requested_common_sampling_floor": requested_floor,
        "preexisting_contextualized_biomass_lower_bound": previous_lower_bound,
        "effective_sampling_biomass_lower_bound": effective_floor,
    }


def validate_model(model, objective_id: str, required: Sequence[str], *, objective_tolerance: float = 1.0e-8,
                   objective_direction: str = "max") -> dict[str, Any]:
    missing = set(required) - {r.id for r in model.reactions}
    if missing: raise KeyError(f"Required reactions absent: {sorted(missing)}")
    bad = [r.id for r in model.reactions if r.lower_bound > r.upper_bound]
    if bad: raise ValueError(f"LB > UB: {bad[:10]}")
    qualified = evaluate_biological_objective(model, objective_id, tolerance=objective_tolerance, direction=objective_direction)
    return {"status": qualified["status"], "objective": objective_id,
            "objective_value": qualified["capacity"], "objective_direction": qualified["direction"],
            "objective_tolerance": objective_tolerance, "model_hash": model_hash(model),
            "reactions": len(model.reactions)}
