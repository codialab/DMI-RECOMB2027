"""RNA fitting adapters refactored from the validated Stage-10 implementations.

Primary adapters receive an already resolved technical/environmental protection
set.  They have no DMI measurement input.
"""
from __future__ import annotations
import time
from dataclasses import dataclass
from typing import Any, Sequence
import cobra
import numpy as np
import pandas as pd
from .models import (enforce_gurobi, expression_thresholds, reaction_expression,
    reaction_expression_with_evidence)


@dataclass
class ExpertResult:
    model: Any | None
    artifact_kind: str
    native_samples: pd.DataFrame | None
    algorithm_state: dict[str, Any]
    metadata: dict[str, Any]


def _matrix(model):
    return (np.asarray(cobra.util.array.create_stoichiometric_matrix(model, array_type="dense"), dtype=np.float64),
            np.asarray([r.lower_bound for r in model.reactions], dtype=np.float64),
            np.asarray([r.upper_bound for r in model.reactions], dtype=np.float64))


def _reduce(model, keep_indices, required):
    keep = {model.reactions[int(index)].id for index in np.asarray(keep_indices).ravel()}; keep.update(required)
    result = model.copy(); result.remove_reactions([r for r in result.reactions if r.id not in keep], remove_orphans=True); return result


def _numpy_compatible_fastcore(base_class):
    class CompatibleFASTcore(base_class):
        def fastcore(self):
            flipped = False; singleton = False; J, A, P, irreversible = self.preprocessing()
            while J.size > 0:
                P = np.setdiff1d(P, A); support = self.findSparseMode(J, P, singleton); A = np.union1d(A, support).astype(int)
                if np.intersect1d(J, A).size > 0: J = np.setdiff1d(J, A); flipped = False
                else:
                    JiRev = np.setdiff1d(J[0] if singleton else J, irreversible)
                    if flipped or JiRev.size == 0:
                        if singleton: return sorted(np.union1d(A, J))
                        flipped = False; singleton = True
                    else: self.reverse_irreversible_reactions_in_reverse_direction(JiRev); flipped = True
            return sorted(A)
    return CompatibleFASTcore


def gimme_expression_vector(scores: pd.Series) -> np.ndarray:
    values = scores.fillna(0.0).to_numpy(dtype=np.float64)
    if not np.isfinite(values).all() or (values < 0).any(): raise ValueError("GIMME expression must be finite and nonnegative")
    return values


def _eflux_bounds(model, expression: pd.DataFrame, config: dict[str, Any], protected_reactions: Sequence[str] = ()) -> tuple[dict[str, float], dict[str, Any]]:
    c = config["contextualization"]["eflux"]
    if c["formulation"] != "colijn_2009_unit_normalized": raise RuntimeError("Only the documented E-Flux formulation is implemented")
    scores, evidence = reaction_expression_with_evidence(
        model, expression, config["contextualization"]["pooled_statistic"],
        and_operator="MIN", or_operator="SUM",
    )
    finite = scores[np.isfinite(scores) & (scores > 0)]
    if finite.empty: raise RuntimeError("E-Flux has no mapped positive reaction expression")
    maximum = float(finite.max()); protected = set(protected_reactions); bounds: dict[str, float] = {}; mapped = []
    for reaction in model.reactions:
        if reaction.boundary or reaction.id in protected: continue
        value = scores.get(reaction.id, np.nan)
        # Partial GPR evidence uses the documented unmapped fallback. In
        # particular, a mapped zero in an OR rule cannot erase an unmapped arm.
        cap = 1.0 if not np.isfinite(value) else max(0.0, float(value) / maximum)
        bounds[reaction.id] = cap
        if np.isfinite(value): mapped.append(reaction.id)
    counts = pd.Series([evidence[rid] for rid in bounds]).value_counts().to_dict()
    return bounds, {"formulation": c["formulation"], "gpr_and": "MIN", "gpr_or": "SUM", "normalization_denominator": maximum,
        "unmapped_internal_cap": 1.0, "boundary_policy": "preserve_frozen_medium_bounds", "protected_policy": "preserve_technical_environmental_bounds",
        "mapped_reactions": len(mapped), "evidence_status_counts": {str(key): int(value) for key, value in counts.items()}}


def _fit_eflux(model, expression: pd.DataFrame, config: dict[str, Any], protected_reactions: Sequence[str]) -> ExpertResult:
    result = model.copy(); caps, metadata = _eflux_bounds(result, expression, config, protected_reactions)
    for rid, cap in caps.items():
        reaction = result.reactions.get_by_id(rid); old_lb, old_ub = float(reaction.lower_bound), float(reaction.upper_bound)
        reaction.lower_bound = max(old_lb, -cap); reaction.upper_bound = min(old_ub, cap)
        if reaction.lower_bound > reaction.upper_bound: raise RuntimeError(f"E-Flux made {rid} invalid")
    return ExpertResult(result, "rna_bounded_model", None, {"reaction_caps": caps}, metadata)


def _fit_riptide(model, expression: pd.DataFrame, config: dict[str, Any], protected_reactions: Sequence[str], *, set_bounds: bool) -> ExpertResult:
    raise RuntimeError("Legacy RIPTiDe adapter is disabled; use src.riptide_repair.fit_riptide_expert")


def fit_expert(algorithm: str, model, expression: pd.DataFrame, config: dict[str, Any], *, fit_seed: int,
               protected_reactions: Sequence[str], riptide_set_bounds: bool | None = None) -> ExpertResult:
    started = time.perf_counter(); enforce_gurobi(model, seed=fit_seed); c = config["contextualization"]
    required = list(dict.fromkeys(map(str, protected_reactions))); missing = set(required) - {r.id for r in model.reactions}
    if missing: raise KeyError(f"Required reactions absent from parent: {sorted(missing)}")
    scores = reaction_expression(model, expression, c["pooled_statistic"])
    needs_thresholds = algorithm in {"GIMME", "iMAT", "INIT_tINIT", "CORDA", "FASTCORE_FASTCORMICS"}
    low, high = expression_thresholds(scores, c["lower_quantile"], c["upper_quantile"]) if needs_thresholds else (np.nan, np.nan)
    filled = scores.fillna(0.0); common = {"algorithm": algorithm, "threshold_low": low, "threshold_high": high,
        "protected_reactions": required, "training_samples": list(expression.columns), "fit_seed": fit_seed, "dmi_measurements_used": False}
    if algorithm == "E-Flux": result = _fit_eflux(model, expression, config, required)
    elif algorithm == "RIPTiDe":
        raise RuntimeError("Legacy RIPTiDe route is disabled; use src.riptide_repair")
    elif algorithm == "RegrEx": raise RuntimeError("blocked_pending_validation: validated RegrEx-LAD/AOS implementation is unavailable")
    elif algorithm == "CORDA":
        from corda import CORDA, reaction_confidence
        pool = expression.median(axis=1); positive = pool[pool > 0]; glo, ghi = positive.quantile([c["lower_quantile"], c["upper_quantile"]])
        gene_conf = {g: (-1 if v <= 0 else 0 if v < glo else 2 if v < ghi else 3) for g, v in pool.items()}
        confidence = {r.id: int(reaction_confidence(r, gene_conf)) for r in model.reactions}
        # Preserve the validated confidence mapping, but use corda 0.5.1's
        # constructor defaults for algorithm parameters (n=3, penalty_factor=100,
        # support=5) rather than the former project-specific n=1 override. Project
        # requirement roles are not promoted to high confidence.
        worker = CORDA(model.copy(), confidence); worker.build()
        result = ExpertResult(worker.cobra_model(name="stage11 CORDA"), "context_specific_model", None, {"reaction_confidence": confidence},
            {"stochastic": False, "corda_algorithm_parameters": "package_defaults",
             "declared_project_requirements_not_forced": required,
             "forced_confidence_reactions": [], "contextualization_biomass_forcing": False})
    else:
        import gurobipy as gp
        gp.setParam("Threads", 1); gp.setParam("Seed", int(fit_seed)); gp.setParam("OutputFlag", 0)
        S, lb, ub = _matrix(model); rids = [r.id for r in model.reactions]; mids = [m.id for m in model.metabolites]; force = [rids.index(rid) for rid in required]
        if algorithm == "GIMME":
            from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties
            threshold = float(c["gimme"]["expression_threshold"]); objective = {rids.index(config["objective"]["biomass_reaction"]): 1.0}
            props = GIMMEProperties(gimme_expression_vector(scores), [objective], obj_frac=float(c["objective_fraction"]), flux_threshold=threshold, solver="GUROBI", reaction_ids=rids, metabolite_ids=mids)
            keep = np.union1d(GIMME(S, lb, ub, props).run(), force); meta = {"gimme_input_semantics": "native_reaction_expression", "gimme_expression_threshold": threshold,
                "numerical_epsilon": float(c["epsilon"]), "fit_objective_fraction": float(c["objective_fraction"])}
        elif algorithm == "iMAT":
            raise RuntimeError("Legacy iMAT route is disabled; use src.imat_repair")
        elif algorithm == "INIT_tINIT":
            from troppo.methods.reconstruction.tINIT import tINIT, tINITProperties
            p = c["tinit"]; centered = np.log1p(filled.to_numpy()) - np.log1p(low)
            raw = tINIT(S, lb, ub, tINITProperties(centered.tolist(), essential_reactions=force, production_weight=float(p["production_weight"]), allow_excretion=bool(p["allow_excretion"]), no_reverse_loops=bool(p["no_reverse_loops"]), solver="GUROBI")).run()
            if raw is None or len(raw) == 0: raise RuntimeError("tINIT did not return a qualified completed reconstruction")
            keep = np.union1d(raw, force); meta = {"task_weight_semantics": p}
        elif algorithm == "FASTCORE_FASTCORMICS":
            from cobra.flux_analysis import fastcc
            from troppo.methods.reconstruction.fastcore import FASTcore, FastcoreProperties
            consistent = fastcc(model.copy(), flux_threshold=float(c["epsilon"])); ids = [r.id for r in consistent.reactions]
            absent = set(required) - set(ids)
            if absent: raise RuntimeError(f"FASTCC removed required infrastructure: {sorted(absent)}")
            S, lb, ub = _matrix(consistent); force = [ids.index(rid) for rid in required]
            core = np.union1d(np.where(filled.reindex(ids).fillna(0).to_numpy() >= high)[0], force).astype(int).tolist()
            keep = _numpy_compatible_fastcore(FASTcore)(S, lb, ub, FastcoreProperties(core, float(c["epsilon"]), "GUROBI")).run()
            result = ExpertResult(_reduce(consistent, keep, required), "context_specific_model", None, {"fastcc_preprocessed": True}, {"numpy2_compatibility_repair": True})
            result.metadata = common | result.metadata | {"runtime_seconds": time.perf_counter() - started}; return result
        else: raise ValueError(f"Unknown algorithm: {algorithm}")
        result = ExpertResult(_reduce(model, keep, required), "context_specific_model", None, {"retained_indices": np.asarray(keep).tolist()}, meta)
    enforce_gurobi(result.model, seed=fit_seed)
    result.metadata = common | result.metadata | {"runtime_seconds": time.perf_counter() - started, "reactions_retained": len(result.model.reactions), "reactions_removed": len(model.reactions) - len(result.model.reactions)}
    return result
