"""Scientific QC for Stage-11 expert and ensemble feasibility."""
from __future__ import annotations
from typing import Mapping, Sequence
import numpy as np
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis
from .models import (DMI_OBSERVABLE_GROUPS, apply_biomass_floor,
    enforce_gurobi, evaluate_biological_objective)


def dmi_observable_records() -> list[tuple[str, str]]:
    return [(group, reaction) for group, reactions in DMI_OBSERVABLE_GROUPS.items() for reaction in reactions]


def flux_capability_table(model, parent_reactions: Sequence[str], contextualization_removed_reactions: Sequence[str],
                          tolerance: float, *, biological_objective_id: str | None = None,
                          physiological_fraction: float | None = None,
                          objective_tolerance: float = 1.0e-8,
                          objective_direction: str = "max") -> pd.DataFrame:
    """Classify observables under the production physiological state when supplied.

    Callers that omit an objective receive only the explicitly structural,
    objective-free diagnostic retained for backwards-compatible toy tests.
    """
    parent = set(map(str, parent_reactions)); expert = {reaction.id for reaction in model.reactions}; removed = set(map(str, contextualization_removed_reactions))
    if removed != parent - expert: raise ValueError("Contextualization-removed reaction accounting differs from parent minus expert")
    query_present = [rid for _, rid in dmi_observable_records() if rid in expert]
    ranges = pd.DataFrame(columns=["minimum", "maximum"])
    objective_capacity = np.nan
    physiological_state = biological_objective_id is not None
    if query_present:
        fresh = model.copy(); enforce_gurobi(fresh)
        if biological_objective_id is not None:
            if physiological_fraction is None:
                raise ValueError("physiological_fraction is required with a biological objective")
            qualified = evaluate_biological_objective(
                fresh, biological_objective_id, tolerance=objective_tolerance, direction=objective_direction
            )
            objective_capacity = float(qualified["capacity"])
            apply_biomass_floor(
                fresh, biological_objective_id, float(physiological_fraction),
                tolerance=objective_tolerance, direction=objective_direction,
            )
        # Bounds, including the admitted objective floor, define the FVA state.
        fresh.objective = fresh.problem.Objective(0, direction="min")
        ranges = flux_variability_analysis(fresh, reaction_list=query_present, fraction_of_optimum=0.0, processes=1)
    rows = []
    for group, rid in dmi_observable_records():
        parent_present = rid in parent; expert_present = rid in expert; contextualization_removed = rid in removed
        minimum = float(ranges.loc[rid, "minimum"]) if expert_present else np.nan
        maximum = float(ranges.loc[rid, "maximum"]) if expert_present else np.nan
        nonzero = bool(expert_present and max(abs(minimum), abs(maximum)) > float(tolerance))
        rows.append({"observable_group": group, "reaction_id": rid, "parent_present": parent_present,
            "expert_present": expert_present, "contextualization_removed": contextualization_removed,
            "structurally_absent_from_parent": not parent_present, "fva_minimum": minimum, "fva_maximum": maximum,
            "structural_capability": expert_present, "nonzero_flux_possible": nonzero,
            # Measured-scale and reproduction require a measurement constraint;
            # an RNA-only FVA must not manufacture either conclusion.
            "measured_scale_flux_possible": None, "phenotype_reproduced": None,
            "physiological_state": physiological_state, "objective_capacity": objective_capacity,
            "flux_capable": nonzero, "retained_but_flux_blocked": bool(expert_present and not nonzero)})
    return pd.DataFrame(rows)
