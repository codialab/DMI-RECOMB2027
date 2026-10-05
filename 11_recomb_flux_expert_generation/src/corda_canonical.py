"""Canonical CORDA adapter for Stage 11.

The RNA-to-gene-to-reaction confidence mapping and installed CORDA defaults
match the previously canonical shared adapter. Biomass alone is promoted
through CORDA's native high-confidence class immediately before reconstruction.
No reaction is added to the reconstructed model afterward.
"""
from __future__ import annotations

import time
from typing import Any, Sequence

import pandas as pd

from .expert_fitting import ExpertResult
from .models import enforce_gurobi, expression_thresholds, reaction_expression


def fit_corda_expert(
    model,
    expression: pd.DataFrame,
    config: dict[str, Any],
    *,
    fit_seed: int,
    protected_reactions: Sequence[str],
) -> ExpertResult:
    """Run CORDA with the native RNA confidence map plus biomass-only class 3."""
    from corda import CORDA, reaction_confidence

    started = time.perf_counter()
    enforce_gurobi(model, seed=fit_seed)
    contextualization = config["contextualization"]
    required = list(dict.fromkeys(map(str, protected_reactions)))
    reaction_ids = [reaction.id for reaction in model.reactions]
    missing = sorted(set(required) - set(reaction_ids))
    if missing:
        raise KeyError(f"Required reactions absent from parent: {missing}")

    biomass_id = str(config["objective"]["biomass_reaction"])
    if biomass_id not in reaction_ids:
        raise KeyError(f"Biomass reaction absent from parent: {biomass_id}")
    boundary_ids = {reaction.id for reaction in model.boundary}
    if biomass_id in boundary_ids:
        raise RuntimeError(f"Biomass reaction may not be a boundary reaction: {biomass_id}")

    scores = reaction_expression(
        model, expression, contextualization["pooled_statistic"]
    )
    low, high = expression_thresholds(
        scores,
        contextualization["lower_quantile"],
        contextualization["upper_quantile"],
    )
    pool = expression.median(axis=1)
    positive = pool[pool > 0]
    gene_low, gene_high = positive.quantile([
        contextualization["lower_quantile"],
        contextualization["upper_quantile"],
    ])
    gene_confidence = {
        gene: (
            -1 if value <= 0
            else 0 if value < gene_low
            else 2 if value < gene_high
            else 3
        )
        for gene, value in pool.items()
    }
    rna_confidence = {
        reaction.id: int(reaction_confidence(reaction, gene_confidence))
        for reaction in model.reactions
    }
    confidence = dict(rna_confidence)
    confidence[biomass_id] = 3

    worker = CORDA(model.copy(), confidence)
    worker.build()
    context_model = worker.cobra_model(name="stage11 CORDA")
    enforce_gurobi(context_model, seed=fit_seed)

    metadata = {
        "algorithm": "CORDA",
        "threshold_low": float(low),
        "threshold_high": float(high),
        "protected_reactions": required,
        "training_samples": list(expression.columns),
        "fit_seed": int(fit_seed),
        "dmi_measurements_used": False,
        "stochastic": False,
        "corda_algorithm_parameters": "package_defaults",
        "declared_project_requirements_not_forced": [
            reaction_id for reaction_id in required if reaction_id != biomass_id
        ],
        "forced_confidence_reactions": [biomass_id],
        "manually_high_confidence_reactions": [biomass_id],
        "contextualization_biomass_forcing": True,
        "biomass_preservation_policy": "corda_high_confidence",
        "runtime_seconds": time.perf_counter() - started,
        "reactions_retained": len(context_model.reactions),
        "reactions_removed": len(model.reactions) - len(context_model.reactions),
    }
    return ExpertResult(
        context_model,
        "context_specific_model",
        None,
        {
            "reaction_confidence": confidence,
            "rna_reaction_confidence": rna_confidence,
        },
        metadata,
    )
