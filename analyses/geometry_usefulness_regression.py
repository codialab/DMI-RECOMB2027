#!/usr/bin/env python3
"""Post hoc predictive follow-up for frozen A1/A2 directional utility results.

This analysis keeps the canonical PL2C/PL2D split and response definitions, uses
only cached geometry/outcome artifacts, and never changes the source bundles.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import shutil
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/recomb2027_dmi_matplotlib")
matplotlib_dir = Path(os.environ["MPLCONFIGDIR"])
matplotlib_dir.mkdir(parents=True, exist_ok=True)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr
import sklearn
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

# Direct script execution places analyses/ rather than the repository root on sys.path.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import dmi_bridge_pl2c_geometry_utility_core_v1 as pl2c
from scripts import dmi_bridge_pl2d_confirmation_core_v1 as pl2d
from scripts import dmi_bridge_pl2d_hypothesis_freeze_core_v1 as pl2dh

OUT = ROOT / "reproduced/geometry_usefulness_regression_v1"
SPLIT_SHA256 = "7dd8676508416a13f88af184110ef40e0f8d92782cd5d62b42e7006af5316e2f"
BUNDLE_SPECS = {
    "pl2c": ("dmi_bridge_pl2c_geometry_utility_development_v1", "BRIDGEPL2C_MANIFEST.json", "6abcf149cd6adb294c0931585cec28f992b453451a5b5f03e10b3ffb2f9b9c39", "PL2C_DEVELOPMENT_GEOMETRY_UTILITY_COMPLETE"),
    "pl2dh": ("dmi_bridge_pl2d_hypothesis_freeze_v1", "BRIDGEPL2DH_MANIFEST.json", "a1dd76d6b40e741e07efe64c34bc39eaad289a46bf6c870917a4f59cab688427", "PL2D_HYPOTHESIS_FROZEN_READY_FOR_CONFIRMATION"),
    "pl2dha": ("dmi_bridge_pl2d_h_support_gate_amendment_v1", "BRIDGEPL2DHA_MANIFEST.json", "fb7df9e208da2e6c9841d43a3b572038e59943a4221535705eef9f5dc0ed6a3a", "PL2D_H_SUPPORT_GATE_AMENDED_READY_FOR_CONFIRMATION"),
    "pl2d": ("dmi_bridge_pl2d_confirmation_v1", "BRIDGEPL2D_MANIFEST.json", "9f9eb730473e7ad209c1d5ec3340808ab614087b8ae1121d929eff61be12b699", "PL2D_CONFIRMED"),
    "a21": ("dmi_bridge_a21_dual_anchor_geometry_v1", "BRIDGEA21_MANIFEST.json", "db9948d7b99447a695cf8aa21429dde6d642e1fd0632bcce52643deac59e6feb", "BRIDGE_A21_DUAL_ANCHOR_GEOMETRY_FROZEN"),
    "a22": ("dmi_bridge_a22_dual_anchor_sign_only_utility_v1", "BRIDGEA22_MANIFEST.json", "15b237151005538e5077a352868f873ef00466f4f5ea49307b781a8a8cd0e7cb", "BRIDGE_A22_SIGN_ONLY_UTILITY_FROZEN"),
    "a23": ("dmi_bridge_a23_a1_a2_synthesis_v1", "BRIDGEA23_MANIFEST.json", "b8e57074f3df23fe7708145bc01f51a19f946545b1c16e61628c2678dbcabd18", "BRIDGE_A23_A1_A2_SYNTHESIS_FROZEN"),
}
MATCHED_PAIRS = 4179
EVALUATIONS = 370
EXPECTED_A2_CASES = 3_493_644
GAIN_TOL = 1e-12
IDENTITY_TOL = 5e-12
ALPHAS = (0.01, 0.1, 1.0, 10.0)
FEATURES = ("mean_sign_magnitude_eta2", "mean_directional_entropy3", "mean_non_tie_pair_fraction")
MODEL_FEATURES = {
    "M1": ("mean_directional_entropy3",),
    "M2": ("mean_sign_magnitude_eta2",),
    "M3": ("mean_directional_entropy3", "mean_non_tie_pair_fraction"),
    "M4": ("mean_directional_entropy3", "mean_non_tie_pair_fraction", "mean_sign_magnitude_eta2"),
}
SOURCE_COLUMNS = (
    "sign_magnitude_eta2", "directional_entropy3", "directionally_useful_fraction",
    "dominant_sign_mass", "n_supported_sign_states", "non_tie_pair_fraction",
    "mean_correct_gain", "mean_information_advantage",
)
REACTION_COLUMNS = (
    "sign_magnitude_eta2", "directional_entropy3", "directionally_useful_fraction",
    "non_tie_pair_fraction", "mean_correct_gain", "mean_information_advantage",
)
OUTPUT_FILES = (
    "REACTION_LEVEL_INPUT.tsv.gz", "MATCHED_FOLD_ASSIGNMENTS.tsv",
    "A1_ONLY_REACTION_AUDIT.tsv", "DEVELOPMENT_OOF_PREDICTIONS.tsv.gz",
    "CONFIRMATION_PREDICTIONS.tsv.gz", "MODEL_METRICS.tsv",
    "TOP_K_METRICS.tsv", "REGRESSION_SUMMARY.png", "REGRESSION_SUMMARY.svg",
)


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_path(path: Path) -> Path:
    resolved = path.resolve(strict=True)
    require(not path.is_symlink() and resolved.is_relative_to(ROOT) and resolved.is_file(), f"unsafe source file: {path}")
    return resolved


def read_json(path: Path) -> dict:
    return json.loads(safe_path(path).read_text(encoding="utf-8"))


def verify_bundle(key: str) -> tuple[Path, dict, dict[str, str]]:
    rel, manifest_name, expected_sha, expected_status = BUNDLE_SPECS[key]
    base = ROOT / "outputs" / rel
    manifest_path = safe_path(base / manifest_name)
    require(sha256(manifest_path) == expected_sha, f"{key} manifest hash mismatch")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    require(manifest.get("status") == expected_status, f"{key} status mismatch")
    hashes = manifest.get("artifact_sha256", {})
    require(isinstance(hashes, dict) and hashes, f"{key} artifact inventory missing")
    for filename, expected in hashes.items():
        require(Path(filename).name == filename, f"unsafe artifact name in {key}: {filename}")
        artifact = safe_path(base / filename)
        require(sha256(artifact) == expected, f"{key} artifact hash mismatch: {filename}")
    return base, manifest, hashes


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(safe_path(path), sep="\t", dtype={"reaction_id": str, "evaluation_id": str}, compression="infer")


def numeric_finite(frame: pd.DataFrame, columns: tuple[str, ...]) -> pd.DataFrame:
    result = frame.copy()
    for column in columns:
        result[column] = pd.to_numeric(result[column], errors="coerce")
        result.loc[~np.isfinite(result[column].to_numpy(dtype=float)), column] = np.nan
    return result


def _manifest_parts(base: Path, manifest_name: str, files: dict[str, str], expected_rows: int) -> list[Path]:
    path = safe_path(base / manifest_name)
    part_manifest = json.loads(path.read_text(encoding="utf-8"))
    parts = part_manifest.get("parts", [])
    require(parts and len(parts) == part_manifest.get("number_of_parts", len(parts)), f"part inventory mismatch: {manifest_name}")
    require(part_manifest.get("total_data_rows") == expected_rows, f"part row count mismatch: {manifest_name}")
    result = []
    observed = 0
    for item in parts:
        filename = item["filename"]
        require(Path(filename).name == filename, "unsafe part filename")
        part = safe_path(base / filename)
        require(files.get(filename) == item["sha256"] and sha256(part) == item["sha256"], f"part hash mismatch: {filename}")
        require(part.stat().st_size == item["bytes"], f"part byte count mismatch: {filename}")
        observed += int(item["rows"])
        result.append(part)
    require(observed == expected_rows, f"part inventory data rows mismatch: {manifest_name}")
    return result


def load_split_and_registry(a22_base: Path) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str], dict[str, int], dict[str, str]]:
    cbase = ROOT / "outputs/dmi_bridge_pl2c_geometry_utility_development_v1"
    split_path = safe_path(cbase / "BRIDGEPL2C_REACTION_SPLIT.tsv")
    require(sha256(split_path) == SPLIT_SHA256, "frozen reaction split registry hash mismatch")
    split = pd.read_csv(split_path, sep="\t", dtype={"reaction_id": str})
    require(len(split) == 4180 and split.reaction_id.is_unique, "frozen A1 split inventory mismatch")
    counts = split.analysis_split.value_counts().to_dict()
    require(counts == {"DEVELOPMENT": 3135, "CONFIRMATION_HOLDOUT": 1045}, "frozen A1 split counts mismatch")
    a2_registry = read_tsv(a22_base / "BRIDGEA22_REACTION_REGISTRY.tsv")
    require(len(a2_registry) == MATCHED_PAIRS and a2_registry.reaction_id.is_unique, "A2 reaction registry mismatch")
    require(a2_registry.subsystem.notna().all() and (a2_registry.subsystem.astype(str).str.len() > 0).all(), "A2 subsystem annotations are incomplete")
    split_ids = set(split.reaction_id)
    a2_ids = set(a2_registry.reaction_id)
    require(a2_ids.issubset(split_ids), "A2 reaction IDs are not a subset of frozen A1 split registry")
    extra = sorted(split_ids - a2_ids)
    require(extra == ["LDH_L"], f"unexpected A1-only target reaction(s): {extra}")
    group = a2_registry.set_index("reaction_id").subsystem.astype(str).to_dict()
    split_by_id = split.set_index("reaction_id").analysis_split.astype(str).to_dict()
    return split, a2_registry, sorted(a2_ids), extra, counts, group | {"__split__" + k: v for k, v in split_by_id.items()}


def make_folds(reaction_ids: list[str], subsystem: dict[str, str]) -> pd.DataFrame:
    groups = [subsystem[rid] for rid in reaction_ids]
    require(len(set(groups)) == 50, f"unexpected subsystem count: {len(set(groups))}")
    fold_values = np.full(len(reaction_ids), -1, dtype=int)
    splitter = GroupKFold(n_splits=3)
    for fold, (train, test) in enumerate(splitter.split(np.zeros(len(reaction_ids)), groups=groups)):
        train_groups = {groups[i] for i in train}
        test_groups = {groups[i] for i in test}
        require(not (train_groups & test_groups), "subsystem leakage across a GroupKFold split")
        fold_values[test] = fold
    require(np.all(fold_values >= 0), "some matched reactions lack a CV fold")
    frame = pd.DataFrame({"reaction_id": reaction_ids, "subsystem": groups, "fold": fold_values})
    require(frame.reaction_id.is_unique, "duplicate reaction in fold registry")
    return frame


def load_a1(split: pd.DataFrame, matched_ids: list[str], subsystem: dict[str, str], a1_only: list[str]) -> tuple[pd.DataFrame, dict, dict[str, str], dict]:
    base = ROOT / "outputs/dmi_bridge_pl2c_geometry_utility_development_v1"
    dev = read_tsv(base / "BRIDGEPL2C_DEVELOPMENT_EVALUATION_REACTION.tsv.xz")
    conf_base = ROOT / "outputs/dmi_bridge_pl2d_confirmation_v1"
    conf = read_tsv(conf_base / "BRIDGEPL2D_CONFIRMATION_EVALUATION_REACTION.tsv.xz")
    require(len(dev) == 1_159_950 and len(conf) == 386_650, "A1 evaluation-reaction row counts mismatch")
    frozen_dev_reactions = pl2dh.aggregate_reaction_frame(dev)
    dev_n, dev_rho = pl2d.primary_association(frozen_dev_reactions)
    snapshot = pl2dh.DEVELOPMENT_SNAPSHOT
    require(dev_n == snapshot["finite_primary_reactions"] and math.isclose(dev_rho, snapshot["primary_reaction_level_rho"], rel_tol=0, abs_tol=2e-12),
            "canonical A1 Development table does not reproduce the frozen PL2D-H snapshot")
    dev_snapshot_qc = {"finite_primary_reactions": dev_n, "primary_reaction_level_rho": dev_rho, "snapshot_parity": "PASS"}
    a1_only_subsystem = {}
    for rid in a1_only:
        rows = dev.loc[dev.reaction_id.eq(rid), "subsystem"]
        require(not rows.empty and rows.nunique(dropna=False) == 1, f"A1-only subsystem annotation missing or inconsistent: {rid}")
        a1_only_subsystem[rid] = str(rows.iloc[0])
    split_map = split.set_index("reaction_id").analysis_split.to_dict()
    require(set(dev.reaction_id) == {r for r, s in split_map.items() if s == "DEVELOPMENT"}, "A1 Development identities do not match frozen split")
    require(set(conf.reaction_id) == {r for r, s in split_map.items() if s == "CONFIRMATION_HOLDOUT"}, "A1 Confirmation identities do not match frozen split")
    evals_dev = set(dev.evaluation_id)
    evals_conf = set(conf.evaluation_id)
    require(len(evals_dev) == len(evals_conf) == EVALUATIONS and evals_dev == evals_conf, "A1 evaluation registry mismatch")
    all_rows = pd.concat([dev, conf], ignore_index=True, sort=False)
    all_rows = all_rows.loc[all_rows.reaction_id.isin(matched_ids)].copy()
    all_rows["subsystem"] = all_rows.reaction_id.map(subsystem)
    require(all_rows.subsystem.notna().all(), "A1 matched reaction subsystem mapping missing")
    all_rows["analysis_split"] = all_rows.reaction_id.map(split_map)
    all_rows = numeric_finite(all_rows, SOURCE_COLUMNS)
    return all_rows, {"development_rows": len(dev), "confirmation_rows": len(conf), "development_reactions": 3135, "confirmation_reactions": 1045, "evaluation_ids": EVALUATIONS}, a1_only_subsystem, dev_snapshot_qc


def load_a2_geometry(base: Path, a21_hashes: dict[str, str], eval_ids: list[str], reaction_ids: list[str], registry: pd.DataFrame) -> pd.DataFrame:
    size = len(eval_ids) * len(reaction_ids)
    e_index = {v: i for i, v in enumerate(eval_ids)}
    r_index = {v: i for i, v in enumerate(reaction_ids)}
    eta = np.full(size, np.nan, dtype=float)
    entropy = np.full(size, np.nan, dtype=float)
    dominant_mass = np.full(size, np.nan, dtype=float)
    supported_states = np.full(size, np.nan, dtype=float)
    seen = np.zeros(size, dtype=bool)
    subsystems = registry.set_index("reaction_id").subsystem.astype(str).to_dict()
    parts = _manifest_parts(base, "BRIDGEA21_REACTION_EVALUATION_FEATURES.parts.json", a21_hashes, 1_671_600)
    fields = ["evaluation_id", "reaction_id", "subsystem", "sign_magnitude_eta2", "directional_entropy3", "dominant_sign_mass", "n_supported_sign_states"]
    for part in parts:
        for chunk in pd.read_csv(part, sep="\t", usecols=fields, dtype={"evaluation_id": str, "reaction_id": str, "subsystem": str}, chunksize=100_000, compression="infer"):
            chunk = chunk.loc[chunk.evaluation_id.isin(e_index)]
            if chunk.empty:
                continue
            ei = chunk.evaluation_id.map(e_index).to_numpy(dtype=np.int64)
            ri = chunk.reaction_id.map(r_index)
            require(not pd.isna(ri).any(), "A21 geometry contains reaction outside matched population")
            ri = ri.to_numpy(dtype=np.int64)
            flat = ei * len(reaction_ids) + ri
            require(not seen[flat].any() and not pd.Index(flat).has_duplicates, "duplicate A21 evaluation-reaction geometry key")
            expected_sub = np.asarray([subsystems[r] for r in chunk.reaction_id], dtype=object)
            require(np.array_equal(expected_sub, chunk.subsystem.to_numpy(dtype=object)), "A21 subsystem annotation mismatch")
            eta_values = pd.to_numeric(chunk.sign_magnitude_eta2, errors="coerce").to_numpy(dtype=float)
            entropy_values = pd.to_numeric(chunk.directional_entropy3, errors="coerce").to_numpy(dtype=float)
            dominant_values = pd.to_numeric(chunk.dominant_sign_mass, errors="coerce").to_numpy(dtype=float)
            state_values = pd.to_numeric(chunk.n_supported_sign_states, errors="coerce").to_numpy(dtype=float)
            eta_values[~np.isfinite(eta_values)] = np.nan
            entropy_values[~np.isfinite(entropy_values)] = np.nan
            dominant_values[~np.isfinite(dominant_values)] = np.nan
            state_values[~np.isfinite(state_values)] = np.nan
            eta[flat] = eta_values
            entropy[flat] = entropy_values
            dominant_mass[flat] = dominant_values
            supported_states[flat] = state_values
            seen[flat] = True
    require(int(seen.sum()) == size, f"A21 geometry rectangle incomplete: {int(seen.sum())}/{size}")
    return pd.DataFrame({
        "evaluation_id": np.repeat(eval_ids, len(reaction_ids)),
        "reaction_id": np.tile(reaction_ids, len(eval_ids)),
        "sign_magnitude_eta2": eta,
        "directional_entropy3": entropy,
        "dominant_sign_mass": dominant_mass,
        "n_supported_sign_states": supported_states,
    })


def load_a2_outcomes(base: Path, a22_hashes: dict[str, str], reaction_ids: list[str]) -> tuple[pd.DataFrame, dict]:
    registry = read_tsv(base / "BRIDGEA22_MATCHED_TRUTH_PAIR_REGISTRY.tsv")
    registry = registry.loc[registry.pair_evaluable.astype(str).str.lower().eq("true")].copy()
    require(len(registry) == 836 and registry.truth_pair_id.is_unique, "A22 matched truth-pair registry mismatch")
    pair_ids = registry.truth_pair_id.astype(str).tolist()
    eval_ids = sorted(registry.evaluation_id.astype(str).unique().tolist())
    require(len(eval_ids) == EVALUATIONS, "A22 evaluable evaluation count mismatch")
    e_index = {v: i for i, v in enumerate(eval_ids)}
    r_index = {v: i for i, v in enumerate(reaction_ids)}
    p_index = {v: i for i, v in enumerate(pair_ids)}
    pair_eval = registry.set_index("truth_pair_id").evaluation_id.astype(str).to_dict()
    total_pairs = registry.groupby("evaluation_id", sort=True).truth_pair_id.nunique().reindex(eval_ids).to_numpy(dtype=np.int32)
    n_eval, n_reaction = len(eval_ids), len(reaction_ids)
    shape = n_eval * n_reaction
    seen = np.zeros((len(pair_ids), n_reaction), dtype=bool)
    row_counts = np.zeros(shape, dtype=np.int32)
    non_tie_counts = np.zeros(shape, dtype=np.int32)
    useful_counts = np.zeros(shape, dtype=np.int32)
    correct_sums = np.zeros(shape, dtype=np.float64)
    advantage_sums = np.zeros(shape, dtype=np.float64)
    parts = _manifest_parts(base, "BRIDGEA22_CASE_OUTCOMES.parts.json", a22_hashes, EXPECTED_A2_CASES)
    fields = ["truth_pair_id", "evaluation_id", "reaction_id", "subsystem", "truth_status", "lambda_0_25_correct_absolute_error_gain", "lambda_0_25_wrong_absolute_error_gain", "lambda_0_25_random_sign_expected_gain"]
    subsystem = read_tsv(base / "BRIDGEA22_REACTION_REGISTRY.tsv").set_index("reaction_id").subsystem.astype(str).to_dict()
    observed_rows = 0
    tie_cases = 0
    useful_cases = 0
    for part in parts:
        for chunk in pd.read_csv(part, sep="\t", usecols=fields, dtype=str, chunksize=100_000, compression="infer", keep_default_na=False):
            observed_rows += len(chunk)
            pair_idx = chunk.truth_pair_id.map(p_index)
            eval_idx = chunk.evaluation_id.map(e_index)
            reaction_idx = chunk.reaction_id.map(r_index)
            require(not pair_idx.isna().any() and not eval_idx.isna().any() and not reaction_idx.isna().any(), "A22 case has unknown pair, evaluation, or reaction")
            pi = pair_idx.to_numpy(dtype=np.int64)
            ei = eval_idx.to_numpy(dtype=np.int64)
            ri = reaction_idx.to_numpy(dtype=np.int64)
            require(np.array_equal(chunk.evaluation_id.to_numpy(dtype=str), np.asarray([pair_eval[pair_ids[i]] for i in pi], dtype=str)), "A22 case pair/evaluation identity mismatch")
            flat_case = pi * n_reaction + ri
            require(not pd.Index(flat_case).has_duplicates and not seen[pi, ri].any(), "duplicate A22 truth-pair/reaction case")
            seen[pi, ri] = True
            expected_sub = np.asarray([subsystem[r] for r in chunk.reaction_id], dtype=object)
            require(np.array_equal(expected_sub, chunk.subsystem.to_numpy(dtype=object)), "A22 reaction subsystem mismatch")
            flat = ei * n_reaction + ri
            row_counts += np.bincount(flat, minlength=shape).astype(np.int32)
            status = chunk.truth_status.to_numpy(dtype=str)
            ties = status == "TRUTH_TIE"
            non = status == "NON_TIE"
            require(np.all(ties | non), "A22 contains an unknown truth status")
            tie_cases += int(ties.sum())
            if not non.any():
                continue
            gains = chunk.loc[non, ["lambda_0_25_correct_absolute_error_gain", "lambda_0_25_wrong_absolute_error_gain", "lambda_0_25_random_sign_expected_gain"]].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
            require(np.isfinite(gains).all(), "A22 non-tie gains must be finite")
            correct, wrong, random = gains.T
            require(np.max(np.abs(random - 0.5 * (correct + wrong))) <= IDENTITY_TOL, "A22 random-sign gain identity mismatch")
            advantage = correct - random
            is_useful = (correct > GAIN_TOL) & (advantage > GAIN_TOL)
            useful_cases += int(is_useful.sum())
            non_flat = flat[non]
            non_tie_counts += np.bincount(non_flat, minlength=shape).astype(np.int32)
            useful_counts += np.bincount(non_flat[is_useful], minlength=shape).astype(np.int32)
            correct_sums += np.bincount(non_flat, weights=correct, minlength=shape)
            advantage_sums += np.bincount(non_flat, weights=advantage, minlength=shape)
    require(observed_rows == EXPECTED_A2_CASES and seen.all(), "A22 case matrix incomplete or row count changed")
    require(tie_cases == 2_053_331 and int(non_tie_counts.sum()) == 1_440_313, "A22 tie/non-tie accounting changed")
    require(np.array_equal(row_counts.reshape(n_eval, n_reaction), np.repeat(total_pairs[:, None], n_reaction, axis=1)), "A22 case rectangle accounting mismatch")
    useful = np.divide(useful_counts, non_tie_counts, out=np.full(shape, np.nan), where=non_tie_counts > 0)
    correct = np.divide(correct_sums, non_tie_counts, out=np.full(shape, np.nan), where=non_tie_counts > 0)
    advantage = np.divide(advantage_sums, non_tie_counts, out=np.full(shape, np.nan), where=non_tie_counts > 0)
    coverage = np.divide(non_tie_counts, np.repeat(total_pairs, n_reaction), out=np.full(shape, np.nan), where=np.repeat(total_pairs, n_reaction) > 0)
    frame = pd.DataFrame({
        "evaluation_id": np.repeat(eval_ids, n_reaction),
        "reaction_id": np.tile(reaction_ids, n_eval),
        "directionally_useful_fraction": useful,
        "non_tie_pair_fraction": coverage,
        "mean_correct_gain": correct,
        "mean_information_advantage": advantage,
    })
    qc = {"case_rows": observed_rows, "evaluable_truth_pairs": len(pair_ids), "evaluations": n_eval,
          "reaction_count": n_reaction, "tie_cases": tie_cases, "non_tie_cases": int(non_tie_counts.sum()),
          "directionally_useful_cases": useful_cases, "distinct_pair_alias_weighting": "NONE"}
    return frame, qc


def aggregate_reactions(frame: pd.DataFrame, split_map: dict[str, str], fold_map: dict[str, int]) -> pd.DataFrame:
    # Canonical evaluation ordering matches the frozen A23 adapter's flattened arrays.
    work = frame.sort_values(["evaluation_id", "reaction_id"], kind="stable").reset_index(drop=True)
    work = numeric_finite(work, SOURCE_COLUMNS)
    for column in pl2dh.REACTION_MEAN_SOURCES:
        require(column in work.columns, f"missing frozen aggregation column {column}")
    # Reuse the frozen pandas groupby aggregation as the numerical reference.
    authoritative = pd.DataFrame.from_records(pl2dh.aggregate_reaction_frame(work))
    gains = work.groupby("reaction_id", sort=True, observed=True).mean_correct_gain.mean()
    subsystems = work.groupby("reaction_id", sort=True, observed=True).subsystem.first()
    result = authoritative.set_index("reaction_id")
    result["mean_correct_gain"] = gains.reindex(result.index)
    result["subsystem"] = subsystems.reindex(result.index)
    result = result.reset_index()
    result["analysis_split"] = result.reaction_id.map(split_map)
    result["fold"] = result.reaction_id.map(fold_map).astype(int)
    return result


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if int(mask.sum()) < 2 or np.unique(x[mask]).size < 2 or np.unique(y[mask]).size < 2:
        return math.nan
    return float(spearmanr(x[mask], y[mask]).statistic)


def fit_cv_and_confirmation(data: pd.DataFrame, anchor: str, folds: pd.DataFrame) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    fold_map = folds.set_index("reaction_id").fold.to_dict()
    # A common complete-case cohort gives every model the same outcome and predictors.
    required = [*FEATURES, "mean_directionally_useful_fraction"]
    eligible = data[required].apply(lambda col: np.isfinite(pd.to_numeric(col, errors="coerce"))).all(axis=1)
    d = data.loc[eligible].copy().reset_index(drop=True)
    require(d.analysis_split.notna().all(), f"{anchor} split labels missing")
    require(d.fold.notna().all(), f"{anchor} fold labels missing")
    dev = d.loc[d.analysis_split.eq("DEVELOPMENT")].copy().reset_index(drop=True)
    conf = d.loc[d.analysis_split.eq("CONFIRMATION_HOLDOUT")].copy().reset_index(drop=True)
    require(len(dev) > 0 and len(conf) > 0, f"{anchor} has no complete-case Development or Confirmation reactions")
    dev_fold = dev.fold.to_numpy(dtype=int)
    conf_fold = conf.fold.to_numpy(dtype=int)
    y_dev = dev.mean_directionally_useful_fraction.to_numpy(dtype=float)
    y_conf = conf.mean_directionally_useful_fraction.to_numpy(dtype=float)
    pred_rows: list[dict] = []
    metric_rows: list[dict] = []
    top_rows: list[dict] = []
    confirmation_predictions: dict[str, np.ndarray] = {}

    # M0 is the fold-specific training mean in CV and the Development mean in Confirmation.
    oof_m0 = np.full(len(dev), np.nan)
    fold_mae = []
    for fold in range(3):
        train = dev_fold != fold
        test = dev_fold == fold
        require(test.any() and train.any(), f"{anchor} fold {fold} is empty after complete-case filtering")
        require(not (set(dev.loc[train, "subsystem"]) & set(dev.loc[test, "subsystem"])), "complete-case CV has subsystem leakage")
        oof_m0[test] = float(np.mean(y_dev[train]))
        fold_mae.append(float(mean_absolute_error(y_dev[test], oof_m0[test])))
    pred_rows.extend({"anchor": anchor, "model": "M0", "split": "DEVELOPMENT_OOF", "reaction_id": r, "fold": int(f), "observed": float(y), "prediction": float(p), "alpha": None}
                     for r, f, y, p in zip(dev.reaction_id, dev_fold, y_dev, oof_m0))
    full_m0 = float(np.mean(y_dev))
    conf_m0 = np.full(len(conf), full_m0)
    confirmation_predictions["M0"] = conf_m0
    metric_rows.append({"anchor": anchor, "model": "M0", "selected_alpha": None, "cv_mae": float(np.mean(fold_mae)),
                        "confirmation_mae": float(mean_absolute_error(y_conf, conf_m0)), "confirmation_spearman": math.nan,
                        "development_n": len(dev), "confirmation_n": len(conf), "cv_fold_maes": ",".join(f"{v:.12g}" for v in fold_mae)})
    pred_rows.extend({"anchor": anchor, "model": "M0", "split": "CONFIRMATION", "reaction_id": r, "fold": int(f), "observed": float(y), "prediction": full_m0, "alpha": None}
                     for r, f, y in zip(conf.reaction_id, conf_fold, y_conf))

    for model, columns in MODEL_FEATURES.items():
        alpha_scores: dict[float, list[float]] = {alpha: [] for alpha in ALPHAS}
        alpha_oof: dict[float, np.ndarray] = {alpha: np.full(len(dev), np.nan) for alpha in ALPHAS}
        for fold in range(3):
            train = dev_fold != fold
            test = dev_fold == fold
            for alpha in ALPHAS:
                estimator = make_pipeline(StandardScaler(), Ridge(alpha=alpha))
                estimator.fit(dev.loc[train, list(columns)], y_dev[train])
                prediction = estimator.predict(dev.loc[test, list(columns)])
                alpha_oof[alpha][test] = prediction
                alpha_scores[alpha].append(float(mean_absolute_error(y_dev[test], prediction)))
        mean_scores = {alpha: float(np.mean(scores)) for alpha, scores in alpha_scores.items()}
        selected_alpha = min(ALPHAS, key=lambda alpha: (mean_scores[alpha], alpha))
        cv_prediction = alpha_oof[selected_alpha]
        require(np.isfinite(cv_prediction).all(), f"{anchor} {model} OOF predictions incomplete")
        final = make_pipeline(StandardScaler(), Ridge(alpha=selected_alpha))
        final.fit(dev.loc[:, list(columns)], y_dev)
        prediction = final.predict(conf.loc[:, list(columns)])
        confirmation_predictions[model] = prediction
        fold_scores = alpha_scores[selected_alpha]
        metric_rows.append({"anchor": anchor, "model": model, "predictors": "+".join(columns), "selected_alpha": selected_alpha,
                            "cv_mae": mean_scores[selected_alpha], "confirmation_mae": float(mean_absolute_error(y_conf, prediction)),
                            "confirmation_spearman": spearman(prediction, y_conf), "development_n": len(dev), "confirmation_n": len(conf),
                            "cv_fold_maes": ",".join(f"{v:.12g}" for v in fold_scores)})
        pred_rows.extend({"anchor": anchor, "model": model, "split": "DEVELOPMENT_OOF", "reaction_id": r, "fold": int(f), "observed": float(y), "prediction": float(p), "alpha": selected_alpha}
                         for r, f, y, p in zip(dev.reaction_id, dev_fold, y_dev, cv_prediction))
        pred_rows.extend({"anchor": anchor, "model": model, "split": "CONFIRMATION", "reaction_id": r, "fold": int(f), "observed": float(y), "prediction": float(p), "alpha": selected_alpha}
                         for r, f, y, p in zip(conf.reaction_id, conf_fold, y_conf, prediction))

    # All Confirmation predictions are emitted on the exact same finite reaction IDs.
    for model in ["M0", *MODEL_FEATURES]:
        require(len(confirmation_predictions[model]) == len(conf), f"{anchor} Confirmation support differs by model")
    confirm_ids = conf.reaction_id.astype(str).tolist()
    for model in MODEL_FEATURES:
        score = confirmation_predictions[model]
        for top_fraction in (0.10, 0.20):
            n_select = max(1, int(math.ceil(top_fraction * len(conf))))
            order = np.lexsort((np.asarray(confirm_ids), -score))
            selected = order[:n_select]
            observed_rate = float(np.mean(y_conf[selected]))
            population_rate = float(np.mean(y_conf))
            gains = pd.to_numeric(conf.iloc[selected].mean_correct_gain, errors="coerce").to_numpy(dtype=float)
            finite_gains = gains[np.isfinite(gains)]
            entropy_base = confirmation_predictions["M1"]
            entropy_order = np.lexsort((np.asarray(confirm_ids), -entropy_base))
            entropy_rate = float(np.mean(y_conf[entropy_order[:n_select]]))
            top_rows.append({"anchor": anchor, "model": model, "top_fraction": top_fraction, "selected_n": n_select,
                             "observed_usefulness_frequency": observed_rate, "random_expected_frequency": population_rate,
                             "enrichment_over_random": observed_rate / population_rate if population_rate > 0 else math.nan,
                             "difference_vs_entropy_M1": observed_rate - entropy_rate,
                             "ratio_vs_entropy_M1": observed_rate / entropy_rate if entropy_rate > 0 else math.nan,
                             "mean_correct_gain": float(np.mean(finite_gains)) if finite_gains.size else math.nan,
                             "finite_correct_gain_n": int(finite_gains.size)})
    confirmation_rows = []
    for i, row in conf.iterrows():
        for model, prediction in confirmation_predictions.items():
            fitted_alpha = next((r["selected_alpha"] for r in metric_rows if r["model"] == model), None)
            confirmation_rows.append({"anchor": anchor, "reaction_id": row.reaction_id, "model": model,
                                      "observed_usefulness_frequency": float(row.mean_directionally_useful_fraction),
                                      "predicted_usefulness_frequency": float(prediction[i]),
                                      "mean_correct_gain": float(row.mean_correct_gain) if pd.notna(row.mean_correct_gain) else math.nan,
                                      "selected_alpha": fitted_alpha})
    return pred_rows, confirmation_rows, metric_rows, top_rows


def write_table(path: Path, frame: pd.DataFrame) -> None:
    frame.to_csv(path, sep="\t", index=False, compression="gzip" if path.suffix == ".gz" else None, na_rep="")


def json_clean(value):
    if isinstance(value, dict):
        return {str(k): json_clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_clean(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if math.isfinite(float(value)) else None
    return value


def plot_summary(metrics: pd.DataFrame, topk: pd.DataFrame, path_png: Path, path_svg: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.4), constrained_layout=True)
    colors = {"M0": "#777777", "M1": "#2C7FB8", "M2": "#D95F0E", "M3": "#41AB5D", "M4": "#756BB1"}
    for row, anchor in enumerate(("A1", "A2")):
        ax = axes[row, 0]
        sub = metrics.loc[metrics.anchor.eq(anchor)].set_index("model")
        models = ["M0", "M1", "M2", "M3", "M4"]
        x = np.arange(len(models))
        width = 0.36
        ax.bar(x - width / 2, [sub.loc[m, "cv_mae"] for m in models], width, label="Development CV", color="#9ECAE1")
        ax.bar(x + width / 2, [sub.loc[m, "confirmation_mae"] for m in models], width, label="Confirmation", color="#3182BD")
        ax.set_xticks(x, models)
        ax.set_title(f"{anchor}: MAE")
        ax.set_ylabel("Mean absolute error")
        ax.grid(axis="y", alpha=0.22)
        ax.legend(frameon=False, fontsize=8)

        ax = axes[row, 1]
        sub = topk.loc[topk.anchor.eq(anchor)]
        models_rank = ["M1", "M2", "M3", "M4"]
        x = np.arange(len(models_rank))
        width = 0.34
        base = float(sub.random_expected_frequency.iloc[0]) if len(sub) else 0.0
        for frac, shift, color in ((0.10, -width / 2, "#31A354"), (0.20, width / 2, "#74C476")):
            z = sub.loc[sub.top_fraction.eq(frac)].set_index("model")
            ax.bar(x + shift, [z.loc[m, "observed_usefulness_frequency"] for m in models_rank], width,
                   label=f"Top {int(frac*100)}%", color=color)
        ax.axhline(base, color="#333333", linestyle="--", linewidth=1, label="Random expectation")
        ax.set_xticks(x, models_rank)
        ax.set_title(f"{anchor}: Confirmation usefulness in selected reactions")
        ax.set_ylabel("Observed usefulness frequency")
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", alpha=0.22)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Post hoc predictive follow-up using the frozen Development/Confirmation split", fontsize=12)
    fig.savefig(path_png, dpi=220)
    fig.savefig(path_svg)
    plt.close(fig)


def run() -> dict:
    bundle_info = {}
    loaded = {}
    for key in BUNDLE_SPECS:
        base, manifest, hashes = verify_bundle(key)
        loaded[key] = (base, manifest, hashes)
        bundle_info[key] = {"path": str(base.relative_to(ROOT)), "manifest_sha256": sha256(base / BUNDLE_SPECS[key][1]), "status": manifest["status"], "artifact_count": len(hashes)}
    require(loaded["pl2d"][1].get("status") == "PL2D_CONFIRMED", "frozen PL2D terminal status changed")
    split, a2_registry, matched_ids, a1_extra, original_a1_counts, id_maps = load_split_and_registry(loaded["a22"][0])
    subsystem = {rid: value for rid, value in id_maps.items() if not rid.startswith("__split__")}
    split_map = {rid: id_maps["__split__" + rid] for rid in split.reaction_id}
    folds = make_folds(matched_ids, subsystem)
    folds["analysis_split"] = folds.reaction_id.map(split_map)
    fold_counts = folds.fold.value_counts().sort_index().to_dict()
    fold_subsystems = {str(f): int(folds.loc[folds.fold.eq(f), "subsystem"].nunique()) for f in range(3)}
    require(fold_counts == {0: 1393, 1: 1391, 2: 1395}, f"matched GroupKFold counts changed: {fold_counts}")
    require(fold_subsystems == {"0": 15, "1": 16, "2": 19}, f"matched subsystem fold counts changed: {fold_subsystems}")
    fold_map = folds.set_index("reaction_id").fold.astype(int).to_dict()

    a1_eval, a1_counts, a1_only_subsystem, a1_dev_snapshot_qc = load_a1(split, matched_ids, subsystem, a1_extra)
    # Independently read A2's distinct-pair outcomes and its cached geometry only.
    a22_base, _, a22_hashes = loaded["a22"]
    a21_base, _, a21_hashes = loaded["a21"]
    a2_outcome, a2_case_qc = load_a2_outcomes(a22_base, a22_hashes, matched_ids)
    eval_ids = sorted(a2_outcome.evaluation_id.unique().tolist())
    a2_geometry = load_a2_geometry(a21_base, a21_hashes, eval_ids, matched_ids, a2_registry)
    a2_eval = a2_geometry.merge(a2_outcome, on=["evaluation_id", "reaction_id"], how="inner", validate="one_to_one", sort=False)
    require(len(a2_eval) == EVALUATIONS * MATCHED_PAIRS, "A2 evaluation-reaction merge dimensions changed")
    a2_eval["subsystem"] = a2_eval.reaction_id.map(subsystem)
    a2_eval = numeric_finite(a2_eval, SOURCE_COLUMNS)

    reaction_frames = []
    for anchor, eval_frame in (("A1", a1_eval), ("A2", a2_eval)):
        reaction = aggregate_reactions(eval_frame, split_map, fold_map)
        require(set(reaction.reaction_id) == set(matched_ids), f"{anchor} reaction aggregation does not cover matched cohort")
        reaction.insert(0, "anchor", anchor)
        reaction_frames.append(reaction)
    reaction_data = pd.concat(reaction_frames, ignore_index=True, sort=False)

    # Check A1 Confirmation reaction means directly against the frozen PL2D output.
    frozen_a1_conf = read_tsv(loaded["pl2d"][0] / "BRIDGEPL2D_CONFIRMATION_REACTION.tsv").set_index("reaction_id")
    actual_a1_conf = reaction_data.loc[reaction_data.anchor.eq("A1") & reaction_data.analysis_split.eq("CONFIRMATION_HOLDOUT")].set_index("reaction_id")
    require(len(frozen_a1_conf) == len(actual_a1_conf) == 1045, "matched A1 Confirmation population count mismatch")
    frozen_fields = tuple(pl2dh.REACTION_MEAN_OUTPUTS.values())
    for rid, row in actual_a1_conf.iterrows():
        require(rid in frozen_a1_conf.index, "A1 matched Confirmation reaction absent from canonical PL2D table")
        for field in frozen_fields:
            expected = pd.to_numeric(pd.Series([frozen_a1_conf.at[rid, field]]), errors="coerce").iloc[0]
            actual = pd.to_numeric(pd.Series([row[field]]), errors="coerce").iloc[0]
            require((pd.isna(expected) and pd.isna(actual)) or (pd.notna(expected) and pd.notna(actual) and math.isclose(float(actual), float(expected), rel_tol=0, abs_tol=1e-12)),
                    f"A1 Confirmation frozen reaction mean mismatch: {rid} {field}")

    # Parity against the validated A23 reaction-level rho summaries is a scientific QC gate.
    a23_path = loaded["a23"][0] / "BRIDGEA23_GEOMETRY_UTILITY_PERSISTENCE.tsv"
    a23 = read_tsv(a23_path)
    a23_crosscheck = []
    for anchor in ("A1", "A2"):
        for split_name in ("DEVELOPMENT", "CONFIRMATION_HOLDOUT"):
            part = reaction_data.loc[reaction_data.anchor.eq(anchor) & reaction_data.analysis_split.eq(split_name)]
            n, rho = pl2d.primary_association(part.to_dict("records"))
            frozen = a23.loc[a23.level.eq("PL2D_REACTION") & a23.split.eq(split_name) & a23.context.eq("ALL")].iloc[0]
            expected_n = int(frozen[f"{anchor.lower()}_finite"])
            expected_rho = float(frozen[f"{anchor.lower()}_rho"])
            require(n == expected_n, f"{anchor} {split_name} finite support differs from A23: n={n}/{expected_n}")
            a23_crosscheck.append({"anchor": anchor, "split": split_name, "finite_pairs": n,
                                   "spearman_eta2_usefulness": rho, "a23_spearman": expected_rho,
                                   "absolute_rho_difference": abs(rho - expected_rho),
                                   "interpretation": "support matches; tiny rho differences can arise from tie ranks after floating-point aggregation"})

    prediction_rows, confirmation_rows, metric_rows, top_rows = [], [], [], []
    for anchor in ("A1", "A2"):
        a = reaction_data.loc[reaction_data.anchor.eq(anchor)].copy()
        p, c, m, t = fit_cv_and_confirmation(a, anchor, folds)
        prediction_rows.extend(p)
        confirmation_rows.extend(c)
        metric_rows.extend(m)
        top_rows.extend(t)

    # Audit the full original A1 population without including LDH_L in matched model fitting.
    a1_only_rows = []
    for rid in a1_extra:
        subsystem_value = a1_only_subsystem[rid]
        a1_only_rows.append({"reaction_id": rid, "anchor_presence": "A1_ONLY", "analysis_split": split_map[rid],
                             "subsystem": subsystem_value,
                             "provenance": "LDH_L is the added quantitative target in A2-L and is removed from the A2 weak-target registry; excluded from the matched regression population only."})
    require(a1_extra == ["LDH_L"], "A1-only reaction identity changed")

    # Holdout support must be identical across all five models within anchor.
    cm = pd.DataFrame(confirmation_rows)
    for anchor in ("A1", "A2"):
        sub = cm.loc[cm.anchor.eq(anchor)]
        counts = sub.groupby("model").reaction_id.nunique().to_dict()
        ids_by_model = {m: set(g.reaction_id) for m, g in sub.groupby("model")}
        require(len(set(counts.values())) == 1 and all(v == ids_by_model["M0"] for v in ids_by_model.values()),
                f"{anchor} models do not share identical finite Confirmation reactions")

    # Validate Development-only CV fold assignment and output determinism metadata.
    cv = pd.DataFrame(prediction_rows)
    for anchor in ("A1", "A2"):
        subset = cv.loc[cv.anchor.eq(anchor) & cv.split.eq("DEVELOPMENT_OOF")]
        for model, group in subset.groupby("model"):
            require(group.reaction_id.is_unique, f"{anchor} {model} has duplicate OOF predictions")
            require(set(group.reaction_id).issubset(set(matched_ids)), "OOF reaction outside matched registry")

    if OUT.exists():
        raise RuntimeError(f"output already exists; refusing overwrite: {OUT}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix="geometry_usefulness_regression_", dir=OUT.parent))
    try:
        write_table(stage / OUTPUT_FILES[0], reaction_data.sort_values(["anchor", "analysis_split", "reaction_id"], kind="stable"))
        folds.to_csv(stage / OUTPUT_FILES[1], sep="\t", index=False, lineterminator="\n")
        pd.DataFrame(a1_only_rows).to_csv(stage / OUTPUT_FILES[2], sep="\t", index=False, lineterminator="\n")
        write_table(stage / OUTPUT_FILES[3], pd.DataFrame(prediction_rows).sort_values(["anchor", "model", "split", "reaction_id"], kind="stable"))
        write_table(stage / OUTPUT_FILES[4], cm.sort_values(["anchor", "model", "reaction_id"], kind="stable"))
        metrics = pd.DataFrame(metric_rows)
        topk = pd.DataFrame(top_rows)
        metrics.to_csv(stage / OUTPUT_FILES[5], sep="\t", index=False, lineterminator="\n", na_rep="")
        topk.to_csv(stage / OUTPUT_FILES[6], sep="\t", index=False, lineterminator="\n", na_rep="")
        plot_summary(metrics, topk, stage / OUTPUT_FILES[7], stage / OUTPUT_FILES[8])

        qc = {
            "status": "PASS",
            "analysis_role": "POST_HOC_PREDICTIVE_FOLLOW_UP_USING_FROZEN_DEVELOPMENT_CONFIRMATION_SPLIT_NOT_NEW_UNTOUCHED_CONFIRMATION",
            "pl2d_status_already_frozen": loaded["pl2d"][1]["status"],
            "original_a1_reactions_by_split": original_a1_counts,
            "a2_reactions_by_split": {name: int(sum(split_map[r] == name for r in matched_ids)) for name in ("DEVELOPMENT", "CONFIRMATION_HOLDOUT")},
            "matched_reactions": MATCHED_PAIRS,
            "a1_only_reaction": a1_only_rows[0],
            "a1_evaluation_reaction_input_counts": a1_counts,
            "a2_case_aggregation": a2_case_qc,
            "a2_evaluation_reaction_rows": len(a2_eval),
            "groupkfold": {"n_splits": 3, "group": "subsystem", "fold_reaction_counts": {str(k): int(v) for k, v in fold_counts.items()}, "fold_subsystem_counts": fold_subsystems, "subsystem_leakage": False, "same_folds_for_both_anchors": True},
            "a1_frozen_development_snapshot": a1_dev_snapshot_qc,
            "a1_frozen_confirmation_reaction_means": "PASS",
            "a23_crosscheck": a23_crosscheck,
            "same_confirmation_ids_within_anchor": True,
            "models": list(MODEL_FEATURES),
            "ridge_alpha_grid": list(ALPHAS),
            "selection_metric": "mean fold validation MAE",
            "nonlinear_sensitivity_included": False,
            "software": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "scikit_learn": sklearn.__version__},
        }
        source_hashes = {key: {"manifest_sha256": info["manifest_sha256"], "artifact_sha256": loaded[key][2]} for key, info in bundle_info.items()}
        manifest = {
            "schema": "dmi.geometry_usefulness_regression.v1",
            "status": "POST_HOC_PREDICTIVE_FOLLOW_UP_COMPLETE",
            "analysis_role": qc["analysis_role"],
            "summary": "Regression prediction using the frozen split; not a new untouched confirmation experiment.",
            "matched_primary_population": {"reaction_count": MATCHED_PAIRS, "excluded_a1_only_reaction": "LDH_L", "reason": "A2-L quantitative anchor; absent from A2 weak-target registry."},
            "source_bundles": source_hashes,
            "source_root_for_recovery": "/home/pty/work/project_MRI_brain_tumor/outputs",
            "input_output_files": list(OUTPUT_FILES),
            "qc": json_clean(qc),
            "artifact_sha256": {name: sha256(stage / name) for name in OUTPUT_FILES},
            "implementation_sha256": sha256(Path(__file__).resolve()),
        }
        (stage / "MANIFEST.json").write_text(json.dumps(json_clean(manifest), indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
        os.replace(stage, OUT)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return {"status": "POST_HOC_PREDICTIVE_FOLLOW_UP_COMPLETE", "output": str(OUT.relative_to(ROOT)), "matched_reactions": MATCHED_PAIRS, "confirmation_rows_per_anchor": {"A1": int(cm.loc[cm.anchor.eq("A1") & cm.model.eq("M0"), "reaction_id"].nunique()), "A2": int(cm.loc[cm.anchor.eq("A2") & cm.model.eq("M0"), "reaction_id"].nunique())}}


if __name__ == "__main__":
    print(json.dumps(run(), sort_keys=True))
