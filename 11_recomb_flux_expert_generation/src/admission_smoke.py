"""Non-generating, identity-bound production-admission solve smoke."""
from __future__ import annotations

from typing import Any, Mapping
import math
import pandas as pd

from .core import stable_hash, validate_samples
from .models import apply_biomass_floor, enforce_gurobi
from .boundary_qc import (active_boundary_source_fluxes,
    assert_source_routes_within_resolved_medium, unexpected_boundary_sources)

ADMISSION_SMOKE_SCHEMA = 1
COMMON_FLOOR_FRACTION = 0.90


def admission_smoke_spec(*, record: Mapping[str, Any], semantic_model_hash: str,
                         qualification_identity: str, controller_fingerprint: str,
                         fit_fingerprint: str, sampling_fingerprint: str) -> dict[str, Any]:
    """Stable identity inputs for one disposable solve-only admission check."""
    return {
        "schema": ADMISSION_SMOKE_SCHEMA, "scope": "qualification_only_solve_only",
        "production": False, "algorithm": str(record["algorithm"]),
        "tumor": str(record["tumor"]),
        "generation_spec_identity": str(record["generation_spec_identity"]),
        "semantic_model_hash": str(semantic_model_hash),
        "qualification_identity": str(qualification_identity),
        "full_controller_fingerprint": str(controller_fingerprint),
        "fit_fingerprint": str(fit_fingerprint),
        "sampling_fingerprint": str(sampling_fingerprint),
        "common_floor_fraction": COMMON_FLOOR_FRACTION,
        "implementation": "src/admission_smoke.py",
    }


def admission_smoke_identity(spec: Mapping[str, Any]) -> str:
    return stable_hash(dict(spec))


def run_admission_smoke(model, *, objective_id: str, expected_routes: Mapping[str, Any],
                        validation: Mapping[str, Any], seed: int | None = None) -> dict[str, Any]:
    """Solve one disposable copy; neither sample nor persist a scientific model."""
    smoke = model.copy()
    solver = enforce_gurobi(smoke, seed=seed)
    if int(solver["Threads"]) != 1:
        raise RuntimeError("Admission smoke requires Gurobi Threads=1 after readback")
    assert_source_routes_within_resolved_medium(smoke, expected_routes,
                                                 context="admission smoke model")
    unexpected = unexpected_boundary_sources(smoke, expected_routes)
    if unexpected:
        raise RuntimeError(f"Admission smoke found undeclared open source routes: {sorted(unexpected)}")
    floor = apply_biomass_floor(smoke, objective_id, COMMON_FLOOR_FRACTION,
                                tolerance=1.0e-8, direction="max")
    smoke.objective = objective_id
    smoke.objective_direction = "max"
    solution = smoke.optimize()
    if str(solution.status) != "optimal":
        raise RuntimeError(f"Admission smoke solve status is not exactly optimal: {solution.status}")
    objective = float(solution.objective_value)
    biomass = float(solution.fluxes[objective_id])
    if not math.isfinite(objective) or objective <= 1.0e-8:
        raise RuntimeError("Admission smoke objective is non-finite or nonpositive")
    effective_floor = float(floor["effective_biomass_lower_bound"])
    if biomass < effective_floor - 1.0e-7:
        raise RuntimeError("Admission smoke biomass flux is below the effective common floor")
    samples = pd.DataFrame([solution.fluxes.reindex([r.id for r in smoke.reactions]).to_numpy()],
                           columns=[r.id for r in smoke.reactions])
    qc = validate_samples(samples, smoke, {
        "steady_state_tolerance": 1.0e-7, "bound_tolerance": 1.0e-7,
    })
    active = active_boundary_source_fluxes(smoke, solution.fluxes, tolerance=1.0e-9)
    undeclared_active = sorted(set(active) - set(expected_routes))
    if undeclared_active:
        raise RuntimeError(f"Admission smoke has undeclared active source routes: {undeclared_active}")
    return {
        "solver": solver, "solver_status": str(solution.status), "contextual_bmax": floor["biomass_capacity"],
        "existing_biomass_lower_bound": floor["preexisting_biomass_lower_bound"],
        "requested_common_floor": floor["requested_biomass_floor"],
        "effective_biomass_lower_bound": effective_floor, "achieved_biomass_flux": biomass,
        "objective_value": objective, "sampling_qc": qc, "active_declared_source_routes": active,
        "undeclared_open_source_routes": unexpected, "undeclared_active_source_routes": undeclared_active,
        "disposable_model": True, "scientific_model_serialized": False,
    }
