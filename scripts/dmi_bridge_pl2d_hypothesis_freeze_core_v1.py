#!/usr/bin/env python3
"""DMI-BRIDGE PL2D-H: frozen one-shot confirmation definitions.

This module contains only deterministic hypothesis-freeze primitives.  It must not
open PL2D confirmation predictor values or outcomes.
"""
from __future__ import annotations

import hashlib
import math
from collections import defaultdict
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr as scipy_spearmanr

SCHEMA = "bridge.pl2d.hypothesis_freeze.v1"
STATUS = "PL2D_HYPOTHESIS_FROZEN_READY_FOR_CONFIRMATION"
EXPECTED_PL2C_STATUS = "PL2C_DEVELOPMENT_GEOMETRY_UTILITY_COMPLETE"
EXPECTED_PL2C_MANIFEST_SHA256 = "6abcf149cd6adb294c0931585cec28f992b453451a5b5f03e10b3ffb2f9b9c39"
EXPECTED_SPLIT_SHA256 = "7dd8676508416a13f88af184110ef40e0f8d92782cd5d62b42e7006af5316e2f"
EXPECTED_CONFIRMATION_REGISTRY_SHA256 = "02a38afe3cb25b8c85a7774c1f20118ad9601ba74f0df0fbff021ba2494bd2e5"
EXPECTED_SPLIT_COUNTS = {"DEVELOPMENT": 3135, "CONFIRMATION_HOLDOUT": 1045}
EXPECTED_DEVELOPMENT_ROWS = 1_159_950
EXPECTED_EVALUATIONS = 370
EXPECTED_CONTEXTS = 16

EXPECTED_PL2C_ARTIFACT_SHA256 = {
    "BRIDGEPL2C_ANALYSIS_CONTRACT.json": "7b3525d06b73e79bf24d621462603edd38ecb88fd07f1b2a6891a372b79d14ee",
    "BRIDGEPL2C_CONTEXT_ASSOCIATIONS.tsv": "7d5566520de95706cadf3c9cf624bf6dd2d8c05c996fa9e07b7b8a73483bcdc6",
    "BRIDGEPL2C_DEVELOPMENT_EVALUATION_REACTION.tsv.xz": "31c140edf99187a818e1db1a0816a4fb43368674ec96863d37ff5f0a0b3af752",
    "BRIDGEPL2C_ETA2_ENTROPY_MAP.tsv": "d0e0c0c92feeda4a625e71c0bf26fd04db55d9d2dc653f89e9114efc9f27a565",
    "BRIDGEPL2C_FEATURE_ASSOCIATIONS.tsv": "4eb7f41e140aae52095619ff67d895558d911e30dd6985a1dd2b7446492b0dcb",
    "BRIDGEPL2C_FEATURE_QUARTILES.tsv": "51def33bd51138cfd6aa05ae9a014f906ec830791462c3d59a4f1f462a165655",
    "BRIDGEPL2C_QC.json": "ae7f69daa906b62eb148712412c898c5f0549a873083432af3a6dc2d5c3fcde4",
    "BRIDGEPL2C_REACTION_SPLIT.tsv": EXPECTED_SPLIT_SHA256,
    "BRIDGEPL2C_SOURCE_AUDIT.json": "2f5b5a19ac05703290b3006a737df774ed9f9f1feeed035db86435fec4e96819",
}

PRIMARY = {
    "unit": "reaction",
    "predictor": "mean_sign_magnitude_eta2_across_evaluations",
    "response": "mean_directionally_useful_fraction_across_defined_evaluations",
    "expected_direction": "POSITIVE",
    "minimum_finite_reactions": 700,
    "minimum_spearman_rho": 0.50,
    "context_count": 16,
    "minimum_finite_reactions_per_context": 100,
    "minimum_positive_contexts": 12,
}

BOOTSTRAP = {
    "role": "SUPPORTIVE_COMPUTATIONAL_STABILITY_ONLY",
    "replicates": 5000,
    "seed": 20260929,
    "interval": "percentile_95",
    "lower_quantile": 0.025,
    "upper_quantile": 0.975,
}

DEVELOPMENT_SNAPSHOT = {
    "finite_primary_reactions": 2455,
    "primary_reaction_level_rho": 0.8570601823118907,
    "primary_partial_rho_controlling_non_tie_pair_fraction": 0.8356721975298953,
    "eta2_vs_mean_information_advantage_rho": 0.3730769973633655,
    "positive_algorithm_rna_contexts": 16,
    "supportive_directional_entropy3_rho": 0.8102865680476654,
    "supportive_dominant_sign_mass_rho": -0.7966029638193884,
    "supportive_n_supported_sign_states_rho": 0.7964386893537818,
    "supportive_finite_reactions": 2716,
}

SUPPORTIVE_FEATURES = {
    "directional_entropy3": "POSITIVE",
    "dominant_sign_mass": "NEGATIVE",
    "n_supported_sign_states": "POSITIVE",
}

NUMERICAL_REFERENCE = {
    "reaction_mean": "pandas.DataFrame.groupby(reaction_id, sort=True).mean() on finite numeric values",
    "spearman": "scipy.stats.spearmanr on finite paired vectors (average ranks for exact ties)",
    "partial_spearman": "closed-form one-control partial correlation from three scipy.stats.spearmanr coefficients",
    "nonfinite_policy": "pandas.to_numeric(errors=coerce); +/-inf -> NaN; group means skip NaN",
}

REACTION_MEAN_SOURCES = (
    "sign_magnitude_eta2",
    "directionally_useful_fraction",
    "non_tie_pair_fraction",
    "mean_information_advantage",
    "directional_entropy3",
    "dominant_sign_mass",
    "n_supported_sign_states",
)

REACTION_MEAN_OUTPUTS = {
    "sign_magnitude_eta2": "mean_sign_magnitude_eta2",
    "directionally_useful_fraction": "mean_directionally_useful_fraction",
    "non_tie_pair_fraction": "mean_non_tie_pair_fraction",
    "mean_information_advantage": "mean_information_advantage",
    "directional_entropy3": "mean_directional_entropy3",
    "dominant_sign_mass": "mean_dominant_sign_mass",
    "n_supported_sign_states": "mean_n_supported_sign_states",
}


def split_identity_row(row: Mapping[str, str]) -> dict[str, str]:
    """Return only metadata permitted in the PL2D-H confirmation registry."""
    required = ("reaction_id", "split_hash", "split_rank", "analysis_split")
    missing = [k for k in required if k not in row]
    if missing:
        raise ValueError(f"split row missing required fields: {missing}")
    return {k: str(row[k]) for k in required}


def build_confirmation_registry(rows: Iterable[Mapping[str, str]]) -> list[dict[str, str]]:
    """Build the identity-only holdout registry without any predictor/outcome fields."""
    out = []
    counts = defaultdict(int)
    seen_reactions: set[str] = set()
    seen_ranks: set[int] = set()
    for raw in rows:
        row = split_identity_row(raw)
        split = row["analysis_split"]
        counts[split] += 1
        reaction_id = row["reaction_id"]
        if reaction_id in seen_reactions:
            raise ValueError(f"duplicate reaction_id in split registry: {reaction_id}")
        seen_reactions.add(reaction_id)
        try:
            rank = int(row["split_rank"])
        except ValueError as exc:
            raise ValueError(f"invalid split_rank for {reaction_id}: {row['split_rank']!r}") from exc
        if rank in seen_ranks:
            raise ValueError(f"duplicate split_rank: {rank}")
        seen_ranks.add(rank)
        if split == "CONFIRMATION_HOLDOUT":
            out.append(row)
    if dict(counts) != EXPECTED_SPLIT_COUNTS:
        raise ValueError(f"unexpected split counts: {dict(counts)}")
    if len(out) != EXPECTED_SPLIT_COUNTS["CONFIRMATION_HOLDOUT"]:
        raise ValueError(f"unexpected confirmation registry size: {len(out)}")
    out.sort(key=lambda r: (int(r["split_rank"]), r["reaction_id"]))
    expected_ranks = list(range(EXPECTED_SPLIT_COUNTS["CONFIRMATION_HOLDOUT"]))
    observed_ranks = [int(r["split_rank"]) for r in out]
    if observed_ranks != expected_ranks:
        raise ValueError("confirmation holdout is not exactly the first 1,045 frozen split ranks")
    return out


def _finite_float(value: object) -> float | None:
    if value is None or value == "":
        return None
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def spearman(x: Sequence[float], y: Sequence[float]) -> float:
    """Spearman rho using the frozen SciPy numerical reference."""
    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    if xa.shape != ya.shape or xa.ndim != 1 or xa.size < 2:
        return math.nan
    mask = np.isfinite(xa) & np.isfinite(ya)
    xv = xa[mask]
    yv = ya[mask]
    if xv.size < 2:
        return math.nan
    return float(scipy_spearmanr(xv, yv).statistic)


def partial_spearman_one_control(x: Sequence[float], y: Sequence[float], z: Sequence[float]) -> float:
    """One-control partial Spearman from three frozen SciPy Spearman coefficients."""
    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    za = np.asarray(z, dtype=float)
    if xa.shape != ya.shape or xa.shape != za.shape or xa.ndim != 1 or xa.size < 3:
        return math.nan
    mask = np.isfinite(xa) & np.isfinite(ya) & np.isfinite(za)
    xv = xa[mask]
    yv = ya[mask]
    zv = za[mask]
    if xv.size < 3:
        return math.nan
    rxy = float(scipy_spearmanr(xv, yv).statistic)
    rxz = float(scipy_spearmanr(xv, zv).statistic)
    ryz = float(scipy_spearmanr(yv, zv).statistic)
    denom2 = (1.0 - rxz * rxz) * (1.0 - ryz * ryz)
    if not math.isfinite(denom2) or denom2 <= 0.0:
        return math.nan
    return float((rxy - rxz * ryz) / math.sqrt(denom2))


def aggregate_reaction_frame(frame: pd.DataFrame) -> list[dict[str, object]]:
    """Aggregate canonical PL2C rows to reaction level with pandas groupby means.

    This is the authoritative numerical reference for PL2D-H/PL2D reaction means.
    The use of pandas grouped means is intentional: replacing it with ``math.fsum``
    or another mathematically equivalent summation changes some 1e-16-scale means,
    which changes exact tie/rank identities and no longer reproduces the frozen PL2C
    development snapshot.
    """
    required = {"reaction_id", *REACTION_MEAN_SOURCES}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"reaction aggregation missing required columns: {missing}")
    if frame.empty:
        return []

    work = frame.loc[:, ["reaction_id", *REACTION_MEAN_SOURCES]].copy()
    if work["reaction_id"].isna().any():
        raise ValueError("reaction_id must be non-null")
    work["reaction_id"] = work["reaction_id"].astype(str)
    if (work["reaction_id"] == "").any():
        raise ValueError("reaction_id must be non-empty")

    for column in REACTION_MEAN_SOURCES:
        work[column] = pd.to_numeric(work[column], errors="coerce")
        work.loc[~np.isfinite(work[column].to_numpy(dtype=float)), column] = np.nan

    grouped = work.groupby("reaction_id", sort=True, observed=True)
    means = grouped[list(REACTION_MEAN_SOURCES)].mean()
    counts = grouped[list(REACTION_MEAN_SOURCES)].count()

    out: list[dict[str, object]] = []
    for reaction_id in means.index:
        row: dict[str, object] = {"reaction_id": str(reaction_id)}
        for source in REACTION_MEAN_SOURCES:
            value = float(means.at[reaction_id, source])
            row[REACTION_MEAN_OUTPUTS[source]] = value if math.isfinite(value) else None
        row["n_finite_eta2_evaluations"] = int(counts.at[reaction_id, "sign_magnitude_eta2"])
        row["n_defined_response_evaluations"] = int(counts.at[reaction_id, "directionally_useful_fraction"])
        out.append(row)
    return out


def aggregate_reaction_rows(rows: Iterable[Mapping[str, object]]) -> list[dict[str, object]]:
    """Row-iterator wrapper around the authoritative pandas groupby implementation."""
    records = list(rows)
    if not records:
        return []
    return aggregate_reaction_frame(pd.DataFrame.from_records(records))


def primary_rho(reaction_rows: Iterable[Mapping[str, object]]) -> tuple[int, float]:
    x: list[float] = []
    y: list[float] = []
    for row in reaction_rows:
        a = _finite_float(row.get("mean_sign_magnitude_eta2"))
        b = _finite_float(row.get("mean_directionally_useful_fraction"))
        if a is not None and b is not None:
            x.append(a)
            y.append(b)
    return len(x), spearman(x, y)


def classify_confirmation(*, finite_reactions: int, overall_rho: float, evaluable_contexts: int, positive_contexts: int) -> str:
    """Apply the frozen PL2D primary confirmation gates."""
    if finite_reactions < PRIMARY["minimum_finite_reactions"]:
        return "PL2D_INSUFFICIENT_SUPPORT"
    if evaluable_contexts < PRIMARY["context_count"]:
        return "PL2D_INSUFFICIENT_SUPPORT"
    if not math.isfinite(overall_rho):
        return "PL2D_INSUFFICIENT_SUPPORT"
    if overall_rho < PRIMARY["minimum_spearman_rho"]:
        return "PL2D_NOT_CONFIRMED"
    if positive_contexts < PRIMARY["minimum_positive_contexts"]:
        return "PL2D_NOT_CONFIRMED"
    return "PL2D_CONFIRMED"


def confirmation_registry_tsv_bytes(rows: Sequence[Mapping[str, str]]) -> bytes:
    header = "reaction_id\tsplit_hash\tsplit_rank\tanalysis_split\n"
    body = "".join(
        f"{r['reaction_id']}\t{r['split_hash']}\t{r['split_rank']}\t{r['analysis_split']}\n"
        for r in rows
    )
    return (header + body).encode("utf-8")


def sha256_bytes(blob: bytes) -> str:
    return hashlib.sha256(blob).hexdigest()
