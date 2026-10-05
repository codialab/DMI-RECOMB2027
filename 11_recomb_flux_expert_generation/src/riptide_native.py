"""Narrow, fail-closed adapter for persisted RIPTiDe native GapSplit sampling.

RIPTiDe 3.4.81 exposes contextualization publicly but does not expose a public
resampling method.  Private access is isolated here and source/signature gated.
"""
from __future__ import annotations
import inspect
import random
from importlib.metadata import version
from typing import Any, Mapping
import numpy as np
import pandas as pd
from .core import stable_hash
from .models import enforce_gurobi, model_hash


def api_description() -> dict[str, Any]:
    import riptide
    public = getattr(riptide, "contextualize")
    private = getattr(riptide, "_constrain_for_sampling", None)
    public_resamplers = sorted(name for name in dir(riptide) if not name.startswith("_") and "sampl" in name.lower())
    return {"riptide_version": version("riptide"), "public_contextualize_signature": str(inspect.signature(public)),
        "public_native_resamplers": public_resamplers, "private_sampling_available": private is not None,
        "private_sampling_signature": str(inspect.signature(private)) if private else None,
        "public_contextualize_source_sha256": stable_hash(inspect.getsource(public)),
        "private_sampling_source_sha256": stable_hash(inspect.getsource(private)) if private else None}


def qualification_fingerprint() -> str:
    return stable_hash({
        "api": api_description(),
        "required_contextualization_set_bounds": True,
        "adapter_source": inspect.getsource(native_resample),
    })


def native_resample(model, state: Mapping[str, Any], requested: int, seed: int) -> pd.DataFrame:
    import riptide
    api = api_description()
    if api["public_native_resamplers"]:
        raise RuntimeError("A public RIPTiDe sampling API now exists; adapter requires revalidation")
    helper = getattr(riptide, "_constrain_for_sampling", None)
    if helper is None: raise RuntimeError("Qualified RIPTiDe private native sampling helper is unavailable")
    random.seed(int(seed)); np.random.seed(int(seed))
    fresh = model.copy(); enforce_gurobi(fresh, seed=seed); before = model_hash(fresh)
    samples, _ = helper(model=fresh, max_coefficients=state["maximization_coefficients"], sampling_depth=int(requested),
        objective=state["objective_reaction"], obj_frac=float(state["fraction_of_optimum"]), tasks=list(state.get("metabolic_tasks") or []),
        task_frac=float(state["task_fraction"]), min_flux=float(state["minimum_flux"]), silent=True)
    if isinstance(samples, str): raise RuntimeError(f"RIPTiDe native sampling failed: {samples}")
    result = pd.DataFrame(samples).astype(np.float64)
    if len(result) != int(requested): raise RuntimeError(f"RIPTiDe returned {len(result)} native vectors; requested {requested}")
    if list(result.columns) != [reaction.id for reaction in fresh.reactions]: raise RuntimeError("RIPTiDe native reaction coordinates changed")
    if model_hash(fresh) != before: raise RuntimeError("RIPTiDe native resampling mutated the persisted expert copy")
    return result
