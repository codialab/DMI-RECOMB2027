"""Fold construction reused unchanged in behavior from Stage-10 ``src/folds.py``."""
from __future__ import annotations
import pandas as pd


def canonical_folds(frozen: pd.DataFrame, khalsa_manifest: pd.DataFrame | None = None, external_samples: pd.DataFrame | None = None) -> pd.DataFrame:
    required = {"design", "fold", "contextualization_CT2A", "contextualization_GL261", "gene_CT2A", "gene_GL261"}
    if not required <= set(frozen): raise ValueError(f"Frozen LOO manifest lacks {sorted(required-set(frozen))}")
    label_map = {} if khalsa_manifest is None else khalsa_manifest.set_index("sample_title")["tumor_model"].to_dict(); rows = []
    for row in frozen.drop_duplicates(["design", "fold"]).itertuples(index=False):
        protocol = "khalsa_LOO_A" if str(row.design).startswith("LOO_A") else "khalsa_LOO_B"
        for role, samples in (("train", sorted(str(row.contextualization_CT2A).split(";") + str(row.contextualization_GL261).split(";"))), ("validation", sorted(str(row.gene_CT2A).split(";") + str(row.gene_GL261).split(";")))):
            for sample in samples:
                tumor = label_map.get(sample, "CT2A" if sample.startswith("seta") else "GL261")
                rows.append({"protocol": protocol, "fold": row.fold, "role": role, "sample": sample, "tumor": tumor, "cohort": "Khalsa"})
    for sample, tumor in sorted(label_map.items()): rows.append({"protocol": "khalsa_all_mikolajewicz_validation", "fold": "all", "role": "train", "sample": sample, "tumor": tumor, "cohort": "Khalsa"})
    if external_samples is not None:
        for sample in sorted(external_samples["sample"].astype(str).unique()): rows.append({"protocol": "khalsa_all_mikolajewicz_validation", "fold": "all", "role": "validation", "sample": sample, "tumor": "CT2A" if sample.startswith("CT2A") else "GL261", "cohort": "Mikolajewicz"})
    result = pd.DataFrame(rows).drop_duplicates().sort_values(["protocol", "fold", "role", "sample"]).reset_index(drop=True); assert_no_leakage(result); return result


def assert_no_leakage(folds: pd.DataFrame) -> None:
    for (protocol, fold), block in folds.groupby(["protocol", "fold"]):
        overlap = set(block.loc[block.role.eq("train"), "sample"]) & set(block.loc[block.role.eq("validation"), "sample"])
        if overlap: raise ValueError(f"Leakage in {protocol}/{fold}: {sorted(overlap)}")
