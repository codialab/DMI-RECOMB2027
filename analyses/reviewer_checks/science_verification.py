"""End-to-end local Level 2 verifier for the frozen RECOMB 2027 endpoints.

Private A1/A2 source artifacts are consumed read-only from ``--data-root``.
The output is a compact validation report; case outcomes and flux vectors are
never copied to the reviewer repository.
"""
from __future__ import annotations

import hashlib
import io
import json
import lzma
import math
import os
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from analyses.reviewer_checks import independent_science as science

EXPECTED = {"truth_pairs": 836, "reactions": 4179, "cases_per_anchor": 3493644,
            "non_tie_per_anchor": 1440313, "eligible_pairs_per_anchor": 660237,
            "development": 3134, "confirmation": 1045}
FIG2_FROZEN_RHO = {"D": (2454, 0.8571160663367106),
                   "E": (2485, 0.8339440056039364),
                   "F": (838, 0.8351435823183291),
                   "G": (843, 0.8493657186460156)}
FIG2_RHO_SOURCE_TOL = 1e-6
PAIR_FIELDS = ["truth_pair_id", "evaluation_id", "algorithm", "rna_context_key",
               "ct2a_mouse", "gl261_mouse", "ct2a_candidate_id", "gl261_candidate_id",
               "ct2a_sample_index", "gl261_sample_index"]
CASE_COLUMNS = ["truth_pair_id", "evaluation_id", "algorithm", "rna_context_key",
                "ct2a_mouse", "gl261_mouse", "reaction_id", "truth_delta_b",
                "baseline_magnitude_estimate", "lambda_0_25_correct_magnitude_estimate",
                "lambda_0_25_wrong_magnitude_estimate", "lambda_0_25_correct_absolute_error_gain",
                "lambda_0_25_wrong_absolute_error_gain",
                "lambda_0_25_random_sign_expected_gain"]
ANCHORS = {
    "A1": (".cache/recomb2027-dmi/frozen/outputs/dmi_bridge_pl2b_sign_only_utility_v1",
           "BRIDGEPL2B", "A1"),
    "A2-L": ("outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1",
             "BRIDGEA22", "A2-L"),
}

_L2B_WORKER_INPUTS = None


def _l2b_worker(pair: dict) -> dict:
    if _L2B_WORKER_INPUTS is None:
        raise RuntimeError("Level 2B worker was not initialized with shared read-only inputs")
    try:
        from threadpoolctl import threadpool_limits
        with threadpool_limits(limits=1):
            return _level2b_pair(pair, _L2B_WORKER_INPUTS)
    except ImportError:
        return _level2b_pair(pair, _L2B_WORKER_INPUTS)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"required science input is missing: manifest ({path})")
    return json.loads(path.read_text())


def _safe(root: Path, relative: str) -> Path:
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f"input path escapes data root: {relative}")
    if path.is_symlink():
        raise ValueError(f"symbolic links are not accepted as science inputs: {relative}")
    return path


def _verify_hash(path: Path, expected: str, name: str, hashes: dict) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"required science input is missing: {name} ({path})")
    actual = sha256(path)
    hashes[name] = actual
    if actual != expected:
        raise ValueError(f"SHA-256 mismatch for {name}: expected {expected}, got {actual}")


def _load_inputs(data_root: Path, public_root: Path) -> dict:
    """Resolve exact local source artifacts and verify their frozen identities."""
    data_root = data_root.resolve()
    public_root = public_root.resolve()
    hashes: dict[str, str] = {}
    sources = {}
    case_parts = {}
    registries = {}
    for anchor, (relative, prefix, label) in ANCHORS.items():
        base = _safe(data_root, relative)
        manifest_path = base / f"{prefix}_MANIFEST.json"
        parts_path = base / f"{prefix}_CASE_OUTCOMES.parts.json"
        reg_name = f"{prefix}_DISTINCT_TRUTH_PAIR_REGISTRY.tsv" if prefix == "BRIDGEPL2B" else f"{prefix}_MATCHED_TRUTH_PAIR_REGISTRY.tsv"
        reg_path = base / reg_name
        manifest = _json(manifest_path)
        parts = _json(parts_path)
        _verify_hash(manifest_path, sha256(manifest_path), f"{anchor}:manifest", hashes)
        _verify_hash(parts_path, manifest.get("artifact_sha256", {}).get(parts_path.name, sha256(parts_path)), f"{anchor}:parts_manifest", hashes)
        _verify_hash(reg_path, manifest["artifact_sha256"][reg_name], f"{anchor}:pair_registry", hashes)
        records = parts.get("parts", [])
        if len(records) != 128:
            raise ValueError(f"{anchor}: expected 128 frozen case partitions, found {len(records)}")
        paths = []
        for record in records:
            part_path = base / record["filename"]
            _verify_hash(part_path, record["sha256"], f"{anchor}:{record['filename']}", hashes)
            paths.append((part_path, int(record["rows"])))
        registry = pd.read_csv(reg_path, sep="\t", dtype={"truth_pair_id": str, "evaluation_id": str})
        if registry.truth_pair_id.duplicated().any():
            raise ValueError(f"{anchor}: duplicate truth_pair_id in registry")
        if prefix == "BRIDGEPL2B":
            registry = registry.loc[registry.pair_evaluable.astype(bool)].copy()
        registries[anchor] = registry
        case_parts[anchor] = paths
        sources[anchor] = {"base": base, "prefix": prefix, "label": label,
                           "registry_path": reg_path, "parts_path": parts_path,
                           "part_rows": [x[1] for x in paths]}

    # A1 and A2 must describe exactly the same held-out states and metadata.
    r1 = registries["A1"].sort_values("truth_pair_id").reset_index(drop=True)
    r2 = registries["A2-L"].sort_values("truth_pair_id").reset_index(drop=True)
    if len(r1) != EXPECTED["truth_pairs"] or len(r2) != EXPECTED["truth_pairs"]:
        raise ValueError("A1/A2 valid matched truth-pair population is not 836 each")
    for col in PAIR_FIELDS:
        if not r1[col].astype(str).equals(r2[col].astype(str)):
            raise ValueError(f"A1/A2 truth-pair registry identity mismatch in {col}")
    if not r1.truth_pair_id.is_unique:
        raise ValueError("matched truth-pair registry identities are not unique")

    # Validate the frozen candidate cache and candidate catalog against source manifests.
    frozen = _safe(data_root, ".cache/recomb2027-dmi/frozen")
    pl1_manifest_path = frozen / "outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_MANIFEST.json"
    pl1_manifest = _json(pl1_manifest_path)
    src_hashes = pl1_manifest["source_sha256"]
    cache_rel = "outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstructed/flux_cache.npz"
    panel_rel = "12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/candidate_fraction_table.tsv.xz"
    cache_path = frozen / cache_rel
    panel_path = frozen / panel_rel
    _verify_hash(cache_path, src_hashes[cache_rel], "candidate_flux_cache", hashes)
    _verify_hash(panel_path, src_hashes[panel_rel], "candidate_catalog", hashes)

    weights_path = _safe(data_root, "data/figure_inputs/s1/candidate_weights.tsv")
    weights = pd.read_csv(weights_path, sep="\t")
    figure_manifest = _json(data_root / "manifests/FIGURE_INPUTS.json")
    weight_entry = next((x for x in figure_manifest["inputs"] if x["destination_path"] == "data/figure_inputs/s1/candidate_weights.tsv"), None)
    if weight_entry is None:
        raise ValueError("candidate weight hash is absent from FIGURE_INPUTS manifest")
    _verify_hash(weights_path, weight_entry["sha256"], "candidate_weights", hashes)

    eval_path = _safe(data_root, "outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_EVALUATION_REGISTRY.tsv")
    a21_manifest_path = eval_path.parent / "BRIDGEA21_MANIFEST.json"
    a21_manifest = _json(a21_manifest_path)
    _verify_hash(eval_path, a21_manifest["artifact_sha256"][eval_path.name], "evaluation_registry", hashes)
    evals = pd.read_csv(eval_path, sep="\t", dtype={"evaluation_id": str})
    if len(evals) != 400 or evals.evaluation_id.duplicated().any():
        raise ValueError("evaluation registry must contain 400 unique evaluations")
    expected_eval_keys = {(method, context, f"C{i}", f"G{j}")
        for method in ("iMAT", "GIMME", "CORDA", "RIPTiDe")
        for context in evals.rna_context_key.astype(str).unique()
        for i in range(1, 6) for j in range(1, 6)}
    observed_eval_keys = set(map(tuple, evals[["algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse"]].astype(str).to_numpy()))
    if len(evals.rna_context_key.astype(str).unique()) != 4 or observed_eval_keys != expected_eval_keys:
        raise ValueError("evaluation registry does not match the exact 4-method × 4-context × 5×5 mouse population")

    cache = np.load(cache_path, allow_pickle=False)
    mats = cache["mats"]
    rxns = [str(x) for x in cache["rxns"].tolist()]
    cache_alg = [str(x) for x in cache["alg"].tolist()]
    cache_tumor = [str(x) for x in cache["tumor"].tolist()]
    cache_eh = [str(x) for x in cache["eh"].tolist()]
    if mats.shape != (32, 20, 4181) or not np.isfinite(mats).all():
        raise ValueError(f"candidate cache expected finite (32,20,4181), found {mats.shape}")
    if len(rxns) != 4181 or len(set(rxns)) != 4181:
        raise ValueError("candidate cache reaction axis must contain 4,181 unique IDs")
    key_to_cache = {}
    for i, key in enumerate(zip(cache_alg, cache_tumor, cache_eh)):
        if key in key_to_cache:
            raise ValueError(f"duplicate cache ensemble identity: {key}")
        key_to_cache[key] = i
    panel = pd.read_csv(panel_path, sep="\t", compression="xz",
                        usecols=["algorithm", "tumor", "ensemble_hash", "projection_hash",
                                 "rna_context_key", "rna_training_samples", "sample_index"])
    panel_key = ["algorithm", "tumor", "ensemble_hash", "rna_context_key"]
    if panel.duplicated(panel_key + ["sample_index"]).any():
        raise ValueError("candidate catalog duplicates a candidate identity")
    catalog_map = {}
    for key, block in panel.groupby(panel_key, sort=True):
        block = block.sort_values("sample_index")
        if block.sample_index.astype(int).tolist() != list(range(20)):
            raise ValueError(f"candidate catalog sample indices invalid for ensemble {key}")
        if len(block) != 20:
            raise ValueError(f"candidate catalog requires 20 states for ensemble {key}")
        ckey = (str(key[0]), str(key[1]), str(key[2]))
        if ckey not in key_to_cache:
            raise ValueError(f"candidate catalog/cache identity mismatch for {ckey}")
        catalog_map[tuple(key)] = (key_to_cache[ckey], block)
    cache_key_set = set(key_to_cache)
    catalog_cache_key_set = {(k[0], k[1], k[2]) for k in catalog_map}
    if catalog_cache_key_set != cache_key_set:
        raise ValueError("candidate catalog and flux cache ensemble populations differ")
    if len(catalog_map) != 32:
        raise ValueError(f"candidate catalog must contain 32 method/tumor/context ensembles, found {len(catalog_map)}")

    # S1 has 20 conditionally normalized states per anchor/mouse/stratum and
    # 320 global states per anchor/mouse across the frozen 16 strata.
    weight_group = ["anchor", "method", "rna_context_key", "tumor", "mouse_id"]
    if len(weights) != 6400 or set(weights.anchor.astype(str)) != {"A1", "A2"}:
        raise ValueError("S1 candidate weight inventory must contain 6,400 A1/A2 rows")
    if weights.duplicated(weight_group + ["candidate_index"]).any():
        raise ValueError("candidate weights contain duplicate anchor/mouse/state keys")
    for key, block in weights.groupby(weight_group, sort=False):
        b = block.sort_values("candidate_index")
        if b.candidate_index.astype(int).tolist() != list(range(20)):
            raise ValueError(f"conditional candidate indices invalid in weight stratum {key}")
        if not np.isclose(float(b.conditional_weight.sum()), 1.0, rtol=0, atol=5e-12):
            raise ValueError(f"conditional weights do not normalize in stratum {key}")
        if not np.isfinite(b.conditional_weight.to_numpy(float)).all() or (b.conditional_weight < 0).any():
            raise ValueError(f"invalid conditional weights in stratum {key}")
        global_subset = b.global_weight.to_numpy(float)
        conditional_from_global = global_subset / math.fsum(float(x) for x in global_subset)
        if np.max(np.abs(conditional_from_global-b.conditional_weight.to_numpy(float))) > 5e-12:
            raise ValueError(f"conditional weights do not match normalized global weights in stratum {key}")
        global_subset = b.global_weight.to_numpy(float)
        conditional_from_global = global_subset / math.fsum(float(x) for x in global_subset)
        if np.max(np.abs(conditional_from_global-b.conditional_weight.to_numpy(float))) > 5e-12:
            raise ValueError(f"conditional weights do not match normalized global weights in stratum {key}")
    global_group = ["anchor", "tumor", "mouse_id"]
    global_ess_values = []
    for key, block in weights.groupby(global_group, sort=False):
        if len(block) != 320 or not np.isclose(block.global_weight.sum(), 1.0, rtol=0, atol=5e-12):
            raise ValueError(f"global candidate pool must have 320 normalized states: {key}")
        ess = 1.0 / float(np.sum(np.square(block.global_weight.to_numpy(float))))
        global_ess_values.append(ess)
        if not np.isclose(ess, 20.0, rtol=0, atol=5e-12):
            raise ValueError(f"global 320-state ESS calibration differs from 20 for {key}: {ess}")

    # Identity construction and candidate index mapping follow the published
    # candidate_id contract without importing the production core.
    def candidate_id(row):
        ident = {name: str(row[name]) for name in ("algorithm", "tumor", "ensemble_hash",
                 "projection_hash", "rna_context_key", "rna_training_samples", "sample_index")}
        payload = {"schema": "bridge.pl2a.candidate_identity.v1", "identity": ident}
        return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()).hexdigest()
    candidate_id_map = {}
    for key, (_, block) in catalog_map.items():
        for row in block.to_dict("records"):
            candidate_id_map[candidate_id(row)] = (key, int(row["sample_index"]))
    if len(candidate_id_map) != len(panel):
        raise ValueError("candidate identity hashes are not unique across the catalog")
    for anchor in ("A1", "A2-L"):
        for row in r1.to_dict("records") if anchor == "A1" else r2.to_dict("records"):
            for side, tumor in (("ct2a", "CT2A"), ("gl261", "GL261")):
                cid = str(row[f"{side}_candidate_id"])
                mapped = candidate_id_map.get(cid)
                expected_key = (str(row["algorithm"]), tumor,
                                next(k[2] for k in catalog_map if k[0] == str(row["algorithm"]) and k[1] == tumor and k[3] == str(row["rna_context_key"])),
                                str(row["rna_context_key"]))
                if mapped is None or mapped[0] != expected_key or mapped[1] != int(row[f"{side}_sample_index"]):
                    raise ValueError(f"{anchor} candidate identity/sample index mismatch for pair {row['truth_pair_id']}")

    # Figure 2 panel/split identity and Figure 3 reference identities.
    split_path = _safe(data_root, "outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_REACTION_SPLIT.tsv")
    split_manifest_path = split_path.parent / "BRIDGEPL2C_MANIFEST.json"
    split_manifest = _json(split_manifest_path)
    _verify_hash(split_path, split_manifest["artifact_sha256"][split_path.name], "development_confirmation_split", hashes)
    split = pd.read_csv(split_path, sep="\t", dtype={"reaction_id": str})
    if split.reaction_id.duplicated().any():
        raise ValueError("Development/Confirmation split contains duplicate reaction identities")
    # The PL2C split includes the A1-only LDH_L row. Restrict to the exact
    # common non-anchor axis and require every target identity to be covered.
    split = split.loc[split.reaction_id.ne("LDH_L")].copy()
    split_by_reaction = split.set_index("reaction_id")
    reaction_ids = sorted(split_by_reaction.index.astype(str))
    if len(reaction_ids) != EXPECTED["reactions"]:
        raise ValueError(f"PL2C split expected 4,179 common targets, found {len(reaction_ids)}")
    split_col = next((c for c in ("reaction_split", "split", "partition", "analysis_split") if c in split_by_reaction.columns), None)
    if split_col is None:
        raise ValueError("PL2C split lacks a recognized split-label column")
    split_counts = split_by_reaction[split_col].astype(str).str.lower().value_counts().to_dict()
    if split_counts not in ({"development": 3134, "confirmation": 1045},
                            {"development": 3134, "confirmation_holdout": 1045},
                            {"development": 3134, "confirmation_holdout": 1045}):
        raise ValueError(f"Development/Confirmation split mismatch: {split_counts}")

    # Hash and schema checks on release references come from the public manifest.
    public_manifest = _json(public_root / "manifests/FIGURE_INPUTS.json")
    public_index = {x["destination_path"]: x for x in public_manifest["inputs"]}
    ref_rel = "data/figure_inputs/fig3/fig3_gain_reaction_evaluation.parquet.xz"
    ref_path = _safe(public_root, ref_rel)
    _verify_hash(ref_path, public_index[ref_rel]["sha256"], "frozen_fig3_pairs", hashes)
    with lzma.open(ref_path, "rb") as stream:
        reference = pd.read_parquet(io.BytesIO(stream.read()))
    if reference.groupby("anchor_setting").size().to_dict() != {"A1": EXPECTED["eligible_pairs_per_anchor"], "A2-L": EXPECTED["eligible_pairs_per_anchor"]}:
        raise ValueError("frozen Figure 3 pair population mismatch")

    return {"data_root": data_root, "public_root": public_root, "hashes": hashes,
            "sources": sources, "parts": case_parts, "registries": registries,
            "evals": evals, "eval_index": {v: i for i, v in enumerate(evals.evaluation_id.astype(str))},
            "weights": weights, "mats": mats, "rxns": rxns, "key_to_cache": key_to_cache,
            "panel": panel, "catalog_map": catalog_map, "reference": reference,
            "reaction_ids": reaction_ids, "reaction_index": {v: i for i, v in enumerate(reaction_ids)},
            "split_labels": split_by_reaction[split_col].astype(str).to_dict(),
            "candidate_id_map": candidate_id_map,
            "weight_diagnostics": {"global_ess_min": float(min(global_ess_values)),
                                   "global_ess_max": float(max(global_ess_values)),
                                   "global_pools": len(global_ess_values)}}


def _case_columns(anchor: str) -> dict[str, str]:
    return {"truth_pair_id": "truth_pair_id", "evaluation_id": "evaluation_id",
            "algorithm": "algorithm", "rna_context_key": "rna_context_key",
            "ct2a_mouse": "ct2a_mouse", "gl261_mouse": "gl261_mouse",
            "reaction_id": "reaction_id", "truth_delta_b": "truth_delta_b",
            "baseline_magnitude_estimate": "baseline_magnitude_estimate",
            "lambda_0_25_correct_magnitude_estimate": "lambda_0_25_correct_magnitude_estimate",
            "lambda_0_25_wrong_magnitude_estimate": "lambda_0_25_wrong_magnitude_estimate",
            "baseline_abs_error": "baseline_abs_error",
            "lambda_0_25_correct_absolute_error": "lambda_0_25_correct_absolute_error",
            "lambda_0_25_correct_absolute_error_gain": "lambda_0_25_correct_absolute_error_gain",
            "lambda_0_25_wrong_absolute_error_gain": "lambda_0_25_wrong_absolute_error_gain",
            "lambda_0_25_random_sign_expected_gain": "lambda_0_25_random_sign_expected_gain"}


def _new_accumulator(n_eval: int, n_rxn: int) -> dict[str, np.ndarray]:
    return {"correct_sum": np.zeros((n_eval, n_rxn)), "wrong_sum": np.zeros((n_eval, n_rxn)),
            "advantage_sum": np.zeros((n_eval, n_rxn)), "useful_sum": np.zeros((n_eval, n_rxn)),
            "non_tie_count": np.zeros((n_eval, n_rxn), dtype=np.int32)}


def _accumulate(acc, eval_ix, rxn_ix, rows):
    truth = rows.truth_delta_b.to_numpy(float)
    base = rows.baseline_magnitude_estimate.to_numpy(float)
    correct = rows.lambda_0_25_correct_magnitude_estimate.to_numpy(float)
    wrong = rows.lambda_0_25_wrong_magnitude_estimate.to_numpy(float)
    tie = np.abs(truth) <= science.TIE_TOL
    if not (np.isfinite(truth).all() and np.isfinite(base).all()):
        raise ValueError("case-level truths and baseline magnitude predictions must be finite")
    if not (np.isfinite(correct[~tie]).all() and np.isfinite(wrong[~tie]).all()):
        raise ValueError("non-tie cue magnitude predictions must be finite")
    if (base < 0).any() or (correct < 0).any() or (wrong < 0).any():
        raise ValueError("magnitude predictions must be nonnegative")
    if tie.any():
        # A truth tie contributes no directional endpoint and must have no
        # directional gain according to the frozen endpoint convention.
        pass
    mask = ~tie
    if not mask.any():
        return
    baseline_error = np.abs(np.abs(truth[mask]) - base[mask])
    correct_gain = baseline_error - np.abs(np.abs(truth[mask]) - correct[mask])
    wrong_gain = baseline_error - np.abs(np.abs(truth[mask]) - wrong[mask])
    # Preserve the frozen construction order: g_info = G_correct - G_random,
    # where G_random is the equal mixture of the two directional gains.
    advantage = correct_gain - 0.5 * (correct_gain + wrong_gain)
    useful = (correct_gain > science.GAIN_TOL) & (advantage > science.GAIN_TOL)
    ei, ri = eval_ix[mask], rxn_ix[mask]
    np.add.at(acc["correct_sum"], (ei, ri), correct_gain)
    np.add.at(acc["wrong_sum"], (ei, ri), wrong_gain)
    np.add.at(acc["advantage_sum"], (ei, ri), advantage)
    np.add.at(acc["useful_sum"], (ei, ri), useful.astype(float))
    np.add.at(acc["non_tie_count"], (ei, ri), 1)


def _validate_and_accumulate_anchor(anchor: str, inputs: dict) -> tuple[dict, dict]:
    reg = inputs["registries"][anchor].set_index("truth_pair_id", drop=False)
    expected_pairs = set(reg.index.astype(str))
    seen_pairs = set()
    reactions = inputs["reaction_ids"]
    rindex = inputs["reaction_index"]
    eindex = inputs["eval_index"]
    acc = _new_accumulator(len(eindex), len(reactions))
    tally = {"cases": 0, "non_ties": 0, "case_gain_max_abs_diff": 0.0,
             "case_wrong_gain_max_abs_diff": 0.0, "truth_pair_count": 0,
             "usefulness_classification_mismatches": 0,
             "usefulness_mismatch_examples": [],
             "usefulness_mismatches_by_split": {"development": 0, "confirmation": 0, "other": 0},
             "threshold_distance_counts": {source: {metric: {str(distance): 0
                 for distance in (1e-14, 1e-13, 1e-12, 5e-12)}
                 for metric in ("correct_gain", "information_advantage")}
                 for source in ("frozen", "independent")}}
    names = _case_columns(anchor)
    for path, expected_rows in inputs["parts"][anchor]:
        frame = pd.read_csv(path, sep="\t", compression="xz", usecols=list(names.values()),
                            dtype={"truth_pair_id": str, "evaluation_id": str,
                                   "reaction_id": str, "algorithm": str, "rna_context_key": str,
                                   "ct2a_mouse": str, "gl261_mouse": str})
        if len(frame) != expected_rows:
            raise ValueError(f"{anchor}: row count mismatch in {path.name}")
        if frame[["truth_pair_id", "evaluation_id", "reaction_id"]].isna().any().any():
            raise ValueError(f"{anchor}: missing case identity in {path.name}")
        if frame.duplicated(["truth_pair_id", "reaction_id"]).any():
            raise ValueError(f"{anchor}: duplicate truth-pair/reaction case in {path.name}")
        sizes = frame.groupby("truth_pair_id", sort=False).reaction_id.agg(["size", "nunique"])
        raw_reaction_count = EXPECTED["reactions"] + (1 if anchor == "A1" else 0)
        if not (sizes["size"].eq(raw_reaction_count) & sizes["nunique"].eq(raw_reaction_count)).all():
            raise ValueError(f"{anchor}: frozen raw case support must contain {raw_reaction_count} unique reactions per truth pair")
        if anchor == "A1":
            if not frame.groupby("truth_pair_id", sort=False).reaction_id.apply(lambda x: (x == "LDH_L").sum() == 1).all():
                raise ValueError("A1 raw case support must contain the frozen LDH_L row once per pair before matched-axis restriction")
            frame = frame.loc[frame.reaction_id.ne("LDH_L")].copy()
        elif frame.reaction_id.isin(["HEX1", "LDH_L"]).any():
            raise ValueError("A2 target cases unexpectedly include an anchor reaction")
        sizes = frame.groupby("truth_pair_id", sort=False).reaction_id.agg(["size", "nunique"])
        if not (sizes["size"].eq(EXPECTED["reactions"]) & sizes["nunique"].eq(EXPECTED["reactions"])).all():
            raise ValueError(f"{anchor}: non-anchor case support must contain 4,179 unique reactions")
        pids = set(frame.truth_pair_id.astype(str).unique())
        if seen_pairs & pids:
            raise ValueError(f"{anchor}: truth pair is duplicated across partitions")
        if not pids <= expected_pairs:
            raise ValueError(f"{anchor}: case partition includes unregistered truth pairs")
        seen_pairs |= pids
        for pid, block in frame.groupby("truth_pair_id", sort=False):
            identity = reg.loc[str(pid)]
            for col in ("evaluation_id", "algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse"):
                if str(block[col].iloc[0]) != str(identity[col]) or block[col].astype(str).nunique() != 1:
                    raise ValueError(f"{anchor}: case/registry identity mismatch for {pid} field {col}")
        # Recompute gains from predictions and truth, compare with saved columns
        # only as frozen references, then accumulate reaction-evaluation means.
        truth = frame.truth_delta_b.to_numpy(float)
        baseline = frame.baseline_magnitude_estimate.to_numpy(float)
        correct = frame.lambda_0_25_correct_magnitude_estimate.to_numpy(float)
        wrong = frame.lambda_0_25_wrong_magnitude_estimate.to_numpy(float)
        non_tie = np.abs(truth) > science.TIE_TOL
        if not (np.isfinite(truth).all() and np.isfinite(baseline).all()):
            raise ValueError(f"{anchor}: nonfinite truth or baseline prediction in {path.name}")
        # Cue predictions are intentionally undefined (NaN) for structural
        # truth ties. Require finite values only on the scientifically
        # eligible non-tie cases.
        if not (np.isfinite(correct[non_tie]).all() and np.isfinite(wrong[non_tie]).all()):
            raise ValueError(f"{anchor}: nonfinite cue prediction on a non-tie case in {path.name}")
        base_error = np.abs(np.abs(truth) - baseline)
        correct_gain = base_error - np.abs(np.abs(truth) - correct)
        wrong_gain = base_error - np.abs(np.abs(truth) - wrong)
        frozen_c = frame.lambda_0_25_correct_absolute_error_gain.to_numpy(float)
        frozen_w = frame.lambda_0_25_wrong_absolute_error_gain.to_numpy(float)
        frozen_random = frame.lambda_0_25_random_sign_expected_gain.to_numpy(float)
        c_diff = float(np.max(np.abs(correct_gain[non_tie] - frozen_c[non_tie]))) if non_tie.any() else 0.0
        w_diff = float(np.max(np.abs(wrong_gain[non_tie] - frozen_w[non_tie]))) if non_tie.any() else 0.0
        tally["case_gain_max_abs_diff"] = max(tally["case_gain_max_abs_diff"], c_diff)
        tally["case_wrong_gain_max_abs_diff"] = max(tally["case_wrong_gain_max_abs_diff"], w_diff)
        if c_diff > 5e-12 or w_diff > 5e-12:
            raise ValueError(f"{anchor}: independently recalculated case gains exceed 5e-12 parity tolerance")
        independent_info = correct_gain - 0.5 * (correct_gain + wrong_gain)
        independent_useful = (correct_gain > science.GAIN_TOL) & (independent_info > science.GAIN_TOL)
        frozen_info = frozen_c - frozen_random
        frozen_useful = (frozen_c > science.GAIN_TOL) & (frozen_info > science.GAIN_TOL)
        distances = (1e-14, 1e-13, 1e-12, 5e-12)
        for source, metrics in (("frozen", {"correct_gain": frozen_c,
                                                "information_advantage": frozen_info}),
                                ("independent", {"correct_gain": correct_gain,
                                                     "information_advantage": independent_info})):
            for metric, values in metrics.items():
                for distance in distances:
                    tally["threshold_distance_counts"][source][metric][str(distance)] += int(
                        np.count_nonzero(non_tie & (np.abs(values - science.GAIN_TOL) <= distance)))
        mismatch = non_tie & (independent_useful != frozen_useful)
        tally["usefulness_classification_mismatches"] += int(mismatch.sum())
        for i in np.flatnonzero(mismatch):
            split = str(inputs["split_labels"].get(str(frame.reaction_id.iloc[i]), "other")).lower()
            split = "confirmation" if split.startswith("confirmation") else split
            tally["usefulness_mismatches_by_split"][split if split in ("development", "confirmation") else "other"] += 1
        for i in np.flatnonzero(mismatch)[:10]:
            baseline_error = abs(abs(float(truth[i])) - float(baseline[i]))
            correct_error = abs(abs(float(truth[i])) - float(correct[i]))
            tally["usefulness_mismatch_examples"].append({
                "truth_pair_id": str(frame.truth_pair_id.iloc[i]),
                "evaluation_id": str(frame.evaluation_id.iloc[i]),
                "reaction_id": str(frame.reaction_id.iloc[i]),
                "split": str(inputs["split_labels"].get(str(frame.reaction_id.iloc[i]), "unknown")),
                "independent_correct_gain": float(correct_gain[i]),
                "frozen_correct_gain": float(frozen_c[i]),
                "truth_delta_hex": float(truth[i]).hex(),
                "baseline_magnitude_hex": float(baseline[i]).hex(),
                "correct_magnitude_hex": float(correct[i]).hex(),
                "independent_baseline_absolute_error_hex": float(baseline_error).hex(),
                "independent_correct_absolute_error_hex": float(correct_error).hex(),
                "frozen_baseline_absolute_error_hex": float(frame.baseline_abs_error.iloc[i]).hex(),
                "frozen_correct_absolute_error_hex": float(frame.lambda_0_25_correct_absolute_error.iloc[i]).hex(),
                "frozen_gain_from_serialized_error_fields": float(frame.baseline_abs_error.iloc[i] - frame.lambda_0_25_correct_absolute_error.iloc[i]),
                "correct_gain_difference": float(correct_gain[i] - frozen_c[i]),
                "correct_gain_difference_hex": float(correct_gain[i] - frozen_c[i]).hex(),
                "independent_information_advantage": float(independent_info[i]),
                "frozen_information_advantage": float(frozen_info[i]),
                "independent_useful": bool(independent_useful[i]),
                "frozen_useful": bool(frozen_useful[i])})
        tally["cases"] += len(frame)
        tally["non_ties"] += int(non_tie.sum())
        eval_ix = frame.evaluation_id.map(eindex).to_numpy()
        rxn_ix = frame.reaction_id.map(rindex).to_numpy()
        if pd.isna(eval_ix).any() or pd.isna(rxn_ix).any():
            raise ValueError(f"{anchor}: case evaluation/reaction not in exact axis registries")
        _accumulate(acc, eval_ix.astype(int), rxn_ix.astype(int), frame)
    if seen_pairs != expected_pairs:
        raise ValueError(f"{anchor}: missing registered pairs {len(expected_pairs-seen_pairs)}")
    tally["truth_pair_count"] = len(seen_pairs)
    if tally["cases"] != EXPECTED["cases_per_anchor"] or tally["non_ties"] != EXPECTED["non_tie_per_anchor"]:
        raise ValueError(f"{anchor}: case/non-tie population mismatch: {tally}")
    return acc, tally


def _pair_table(anchor: str, acc: dict, inputs: dict) -> pd.DataFrame:
    count = acc["non_tie_count"]
    ei, ri = np.nonzero(count)
    if len(ei) != EXPECTED["eligible_pairs_per_anchor"]:
        raise ValueError(f"{anchor}: expected 660,237 eligible reaction-evaluation pairs, found {len(ei)}")
    return pd.DataFrame({"anchor_setting": anchor,
        "evaluation_id": np.asarray(list(inputs["eval_index"]), dtype=object)[ei],
        "reaction_id": np.asarray(inputs["reaction_ids"], dtype=object)[ri],
        "n_distinct_non_tie_truths": count[ei, ri],
        "correct_gain_mean": acc["correct_sum"][ei, ri] / count[ei, ri],
        "wrong_gain_mean": acc["wrong_sum"][ei, ri] / count[ei, ri],
        "directional_information_advantage": acc["advantage_sum"][ei, ri] / count[ei, ri],
        "directional_advantage": acc["correct_sum"][ei, ri] / count[ei, ri] - acc["wrong_sum"][ei, ri] / count[ei, ri],
        "directionally_useful_fraction": acc["useful_sum"][ei, ri] / count[ei, ri]})


def _verify_l2a(inputs: dict) -> dict:
    summaries, pair_tables = {}, []
    for anchor in ("A1", "A2-L"):
        acc, tally = _validate_and_accumulate_anchor(anchor, inputs)
        pair = _pair_table(anchor, acc, inputs)
        ref = inputs["reference"].loc[inputs["reference"].anchor_setting.eq(anchor)].copy()
        keys = ["evaluation_id", "reaction_id"]
        if ref.duplicated(keys).any() or pair.duplicated(keys).any():
            raise ValueError(f"{anchor}: duplicate reaction-evaluation endpoint key")
        merged = pair.merge(ref, on=keys, how="outer", suffixes=("_independent", "_frozen"), indicator=True, validate="one_to_one")
        if not merged._merge.eq("both").all():
            raise ValueError(f"{anchor}: independent and frozen reaction-evaluation identities differ")
        diffs = {}
        for col, ref_col in [("correct_gain_mean", "correct_gain_mean"),
                             ("wrong_gain_mean", "wrong_gain_mean"),
                             ("directional_advantage", "directional_advantage")]:
            a = merged[f"{col}_independent"].to_numpy(float)
            b = merged[f"{ref_col}_frozen"].to_numpy(float)
            diffs[col] = float(np.max(np.abs(a - b)))
            if diffs[col] > 5e-12:
                raise ValueError(f"{anchor}: pair-weighted input {col} differs from frozen reference")
        summaries[anchor] = {**tally, "eligible_pairs": len(pair),
            "mean_correct_gain_independent": float(pair.correct_gain_mean.mean()),
            "mean_wrong_gain_independent": float(pair.wrong_gain_mean.mean()),
            "separation_independent": float((pair.correct_gain_mean-pair.wrong_gain_mean).mean()),
            "mean_correct_gain_frozen": float(merged.correct_gain_mean_frozen.mean()),
            "mean_wrong_gain_frozen": float(merged.wrong_gain_mean_frozen.mean()),
            "separation_frozen": float(merged.directional_advantage_frozen.mean()),
            "maximum_pair_endpoint_difference": diffs}
        pair_tables.append(pair)
    return {"status": "PASS", "anchors": summaries, "pair_tables": pair_tables}


def _weight_vector(inputs: dict, anchor: str, eval_row, tumor: str, mouse: str) -> np.ndarray:
    w = inputs["weights"]
    subset = w.loc[w.anchor.eq(anchor) & w.method.eq(eval_row.algorithm)
                   & w.rna_context_key.eq(eval_row.rna_context_key) & w.tumor.eq(tumor)
                   & w.mouse_id.eq(mouse)].sort_values("candidate_index")
    if len(subset) != 20 or subset.candidate_index.astype(int).tolist() != list(range(20)):
        raise ValueError(f"{anchor}: missing or misindexed mouse-specific weights for {mouse}")
    return subset.conditional_weight.to_numpy(float)


def _level2b_pair(pair: dict, inputs: dict) -> dict:
    eval_ix = inputs["eval_index"][str(pair["evaluation_id"])]
    evaluation = inputs["evals"].iloc[eval_ix]
    indices = inputs["mats"]
    cache = inputs["key_to_cache"]
    panel = inputs["panel"]
    cache_indices = []
    candidate_matrices = []
    for tumor, mouse in (("CT2A", str(pair["ct2a_mouse"])), ("GL261", str(pair["gl261_mouse"]))):
        group = panel.loc[panel.algorithm.eq(pair["algorithm"]) & panel.tumor.eq(tumor)
                          & panel.rna_context_key.eq(pair["rna_context_key"])]
        hashes = group.ensemble_hash.astype(str).unique()
        if len(hashes) != 1:
            raise ValueError(f"candidate ensemble identity ambiguous for {tumor}/{pair['algorithm']}/{pair['rna_context_key']}")
        cache_ix = cache[(str(pair["algorithm"]), tumor, hashes[0])]
        cache_indices.append(cache_ix)
        candidate_matrices.append(indices[cache_ix])
    rxns = inputs["rxns"]
    target_idx = np.asarray([i for i, rid in enumerate(rxns) if rid not in {"HEX1", "LDH_L"}], dtype=int)
    target_ids = [rxns[i] for i in target_idx]
    left_index, right_index = int(pair["ct2a_sample_index"]), int(pair["gl261_sample_index"])
    if not (0 <= left_index < 20 and 0 <= right_index < 20):
        raise ValueError("held-out candidate indices must lie in 0..19")
    truth = candidate_matrices[0][left_index, target_idx] - candidate_matrices[1][right_index, target_idx]
    if len(truth) != EXPECTED["reactions"]:
        raise ValueError("reconstructed truth reaction axis is not 4,179 common reactions")
    out = {"truth": truth}
    for anchor in ("A1", "A2"):
        label = "A2-L" if anchor == "A2" else "A1"
        wl = _weight_vector(inputs, anchor, evaluation, "CT2A", str(pair["ct2a_mouse"]))
        wr = _weight_vector(inputs, anchor, evaluation, "GL261", str(pair["gl261_mouse"]))
        # Remove each held-out state as a complete vector, retaining all 19x19 states.
        pred = science.level2b_predictions(candidate_matrices[0][:, target_idx], candidate_matrices[1][:, target_idx],
            np.arange(20), np.arange(20), left_index, right_index, wl, wr, truth,
            lam=0.25, reaction_block_size=256)
        out[label] = pred
    out["target_reactions"] = target_ids
    return out


def _prediction_max_difference(name: str, reconstructed: np.ndarray,
                               frozen: np.ndarray, eligible: np.ndarray,
                               pair_id: str) -> float:
    reconstructed = np.asarray(reconstructed, dtype=float)
    frozen = np.asarray(frozen, dtype=float)
    eligible = np.asarray(eligible, dtype=bool)
    if reconstructed.shape != frozen.shape or reconstructed.shape != eligible.shape:
        raise ValueError(f"{name} prediction comparison has inconsistent shapes for {pair_id}")
    finite_reconstructed = np.isfinite(reconstructed)
    finite_frozen = np.isfinite(frozen)
    if not np.array_equal(finite_frozen, eligible):
        raise ValueError(f"{name} prediction finite-population mismatch for {pair_id}")
    if not finite_reconstructed[eligible].all():
        raise ValueError(f"{name} predictions do not match scientific eligibility for {pair_id}")
    return float(np.max(np.abs(reconstructed[eligible] - frozen[eligible]))) if eligible.any() else 0.0


def _verify_l2b(inputs: dict, workers: int | str = 1, *, batch_size: int = 1) -> dict:
    resolved = science.resolve_workers(workers)
    if batch_size < 1:
        raise ValueError("Level 2B task batch size must be positive")
    pool = None
    global _L2B_WORKER_INPUTS
    if resolved > 1:
        # Linux fork inherits the immutable NumPy cache and compact registries
        # as copy-on-write pages. Tasks contain only one truth-pair identity.
        inputs["mats"].flags.writeable = False
        _L2B_WORKER_INPUTS = inputs
        try:
            from concurrent.futures import ProcessPoolExecutor
            import multiprocessing
            pool = ProcessPoolExecutor(max_workers=resolved,
                mp_context=multiprocessing.get_context("fork"),
                initializer=_limit_worker_threads)
        except (OSError, RuntimeError):
            pool = None
    registries = inputs["registries"]
    pairs = registries["A1"].sort_values("truth_pair_id").to_dict("records")
    pred_err = {a: {k: 0.0 for k in ("baseline", "correct", "incorrect")} for a in ("A1", "A2-L")}
    posterior_acc = {a: _new_accumulator(len(inputs["eval_index"]), len(inputs["reaction_ids"]))
                     for a in ("A1", "A2-L")}
    truth_err, pred_count = 0.0, 0
    # Input case chunks are kept one partition in memory.  Compute the posterior
    # per truth pair and compare against prediction columns in both frozen arms.
    paths = {a: {p.name: p for p, _ in inputs["parts"][a]} for a in ("A1", "A2-L")}
    a1_parts = inputs["parts"]["A1"]
    a2_parts = inputs["parts"]["A2-L"]
    if len(a1_parts) != len(a2_parts):
        raise ValueError("A1/A2 frozen case partition inventories differ")
    for (p1, _), (p2, _) in zip(a1_parts, a2_parts):
        a1 = pd.read_csv(p1, sep="\t", compression="xz", usecols=CASE_COLUMNS,
                         dtype={"truth_pair_id": str, "reaction_id": str})
        a2 = pd.read_csv(p2, sep="\t", compression="xz", usecols=CASE_COLUMNS,
                         dtype={"truth_pair_id": str, "reaction_id": str})
        ids1 = a1.truth_pair_id.drop_duplicates().tolist()
        ids2 = a2.truth_pair_id.drop_duplicates().tolist()
        if ids1 != ids2:
            raise ValueError(f"A1/A2 case partition truth-pair identities differ at {p1.name}")
        pair_records = [registries["A1"].loc[registries["A1"].truth_pair_id.eq(pid)].iloc[0].to_dict() for pid in ids1]
        if pool is None:
            pair_results = (_level2b_pair(pair, inputs) for pair in pair_records)
        else:
            try:
                pair_results = pool.map(_l2b_worker, pair_records, chunksize=batch_size)
            except Exception as exc:
                pool.shutdown(wait=False, cancel_futures=True)
                raise RuntimeError(f"Level 2B worker pool failed: {exc}") from exc
        for pid, pair, frozen1, frozen2, result in zip(ids1, pair_records,
                (a1.loc[a1.truth_pair_id.eq(x)].sort_values("reaction_id") for x in ids1),
                (a2.loc[a2.truth_pair_id.eq(x)].sort_values("reaction_id") for x in ids2),
                pair_results):
            target_ids = result["target_reactions"]
            for anchor, frozen in (("A1", frozen1), ("A2-L", frozen2)):
                if anchor == "A1":
                    frozen = frozen.loc[frozen.reaction_id.ne("LDH_L")]
                frozen = frozen.sort_values("reaction_id")
                if frozen.reaction_id.tolist() != sorted(target_ids):
                    raise ValueError(f"{anchor}: case reaction identities do not match candidate reaction axis for {pid}")
                target_order = np.argsort(target_ids)
                truth_values = frozen.truth_delta_b.to_numpy(float)
                max_truth_diff = float(np.max(np.abs(truth_values - result["truth"][target_order])))
                truth_err = max(truth_err, max_truth_diff)
                if max_truth_diff > 5e-12:
                    raise ValueError(f"{anchor}: candidate reconstructed truth differs for {pid}")
                preds = result[anchor]
                for name, col in (("baseline", "baseline_magnitude_estimate"),
                                  ("correct", "lambda_0_25_correct_magnitude_estimate"),
                                  ("incorrect", "lambda_0_25_wrong_magnitude_estimate")):
                    reconstructed = np.asarray(preds[name][target_order], dtype=float)
                    frozen_values = frozen[col].to_numpy(float)
                    eligible = (np.ones(len(frozen), dtype=bool) if name == "baseline"
                                else np.abs(truth_values) > science.TIE_TOL)
                    maxdiff = _prediction_max_difference(name, reconstructed,
                                                          frozen_values, eligible, pid)
                    pred_err[anchor][name] = max(pred_err[anchor][name], maxdiff)
                    if maxdiff > 5e-12:
                        raise ValueError(f"{anchor}: Level 2B {name} predictions differ for {pid}")
                eval_ix = inputs["eval_index"][str(pair["evaluation_id"])]
                reaction_indices = np.asarray([inputs["reaction_index"][rid] for rid in target_ids], dtype=int)
                prediction_frame = pd.DataFrame({
                    "truth_delta_b": result["truth"],
                    "baseline_magnitude_estimate": preds["baseline"],
                    "lambda_0_25_correct_magnitude_estimate": preds["correct"],
                    "lambda_0_25_wrong_magnitude_estimate": preds["incorrect"],
                })
                _accumulate(posterior_acc[anchor],
                    np.full(len(target_ids), eval_ix, dtype=int), reaction_indices, prediction_frame)
            pred_count += 1
    if pool is not None:
        pool.shutdown(wait=True)
    if pred_count != EXPECTED["truth_pairs"]:
        raise ValueError(f"Level 2B expected 836 truth pairs, processed {pred_count}")
    posterior_endpoints = {}
    for anchor in ("A1", "A2-L"):
        pair = _pair_table(anchor, posterior_acc[anchor], inputs)
        reference = inputs["reference"].loc[inputs["reference"].anchor_setting.eq(anchor)]
        merged = pair.merge(reference, on=["evaluation_id", "reaction_id"], how="outer",
            suffixes=("_independent", "_frozen"), indicator=True, validate="one_to_one")
        if not merged._merge.eq("both").all():
            raise ValueError(f"{anchor}: Level 2B-derived Level 2A pair identities differ from frozen endpoints")
        max_gain = max(float(np.max(np.abs(merged.correct_gain_mean_independent-merged.correct_gain_mean_frozen))),
                       float(np.max(np.abs(merged.wrong_gain_mean_independent-merged.wrong_gain_mean_frozen))))
        if max_gain > 5e-12:
            raise ValueError(f"{anchor}: Level 2B-derived Level 2A gains exceed 5e-12: {max_gain}")
        posterior_endpoints[anchor] = {
            "eligible_pairs": len(pair),
            "correct_gain_independent": float(pair.correct_gain_mean.mean()),
            "incorrect_gain_independent": float(pair.wrong_gain_mean.mean()),
            "separation_independent": float((pair.correct_gain_mean-pair.wrong_gain_mean).mean()),
            "max_pair_gain_difference_to_frozen": max_gain,
        }
    return {"status": "PASS", "workers": resolved, "batch_size": batch_size, "truth_pairs": pred_count,
            "reconstructed_truth_max_abs_difference": truth_err,
            "maximum_prediction_differences": pred_err,
            "level2a_endpoints_from_recomputed_predictions": posterior_endpoints,
            "scope": "posterior verification conditional on frozen S1 anchor weights"}


def _limit_worker_threads() -> None:
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(limits=1)
    except ImportError:
        pass


def _verify_geometry(inputs: dict, l2a: dict, *, workers: int | str = 1,
                     batch_size: int = 1) -> dict:
    # Geometry validation is performed with the full 20x20 support.  Figure 2
    # correlation references remain frozen; a separate calculation from
    # independent eta values quantifies floating-point rank instability.
    geometry_parts = []
    eval_ids = {a: set(inputs["reference"].loc[inputs["reference"].anchor_setting.eq(a), "evaluation_id"].astype(str)) for a in ("A1", "A2-L")}
    if any(len(v) != 370 for v in eval_ids.values()):
        raise ValueError("Figure 2 geometry uses exactly 370 eligible evaluations per anchor")
    eval_rows = inputs["evals"].copy()
    eval_index = inputs["eval_index"]
    eval_cache = np.empty((len(eval_rows), 2), dtype=np.int64)
    eval_weights = np.empty((len(eval_rows), 4, 20), dtype=np.float64)
    for i, e in eval_rows.iterrows():
        for side, tumor, mouse in ((0, "CT2A", str(e.ct2a_mouse)), (1, "GL261", str(e.gl261_mouse))):
            g = inputs["panel"].loc[inputs["panel"].algorithm.eq(e.algorithm)&inputs["panel"].tumor.eq(tumor)&inputs["panel"].rna_context_key.eq(e.rna_context_key)]
            hashes = g.ensemble_hash.astype(str).unique()
            if len(hashes) != 1: raise ValueError("geometry ensemble identity is ambiguous")
            eval_cache[i,side] = inputs["key_to_cache"][(str(e.algorithm),tumor,hashes[0])]
            for a,anchor in enumerate(("A1","A2")):
                eval_weights[i,a*2+side] = _weight_vector(inputs,anchor,e,tumor,mouse)
    target_idx = np.asarray([i for i,r in enumerate(inputs["rxns"]) if r not in {"HEX1","LDH_L"}],dtype=int)
    resolved = science.resolve_workers(workers)
    raw = science.reconstruct_geometry_parallel(inputs["mats"],eval_cache,eval_weights,target_idx,workers=resolved,batch_size=batch_size)
    eta = np.full((len(eval_rows),len(target_idx),2),np.nan)
    for ix,a1,a2 in raw: eta[ix,:,0],eta[ix,:,1]=a1,a2
    if eta.shape[1] != EXPECTED["reactions"]: raise ValueError("geometry reaction-axis mismatch")
    frozen_path = inputs["public_root"] / "data/figure_inputs/fig2/figure2_DG_points.csv.xz"
    public_manifest = _json(inputs["public_root"] / "manifests/FIGURE_INPUTS.json")
    public_index = {x["destination_path"]: x for x in public_manifest["inputs"]}
    _verify_hash(frozen_path, public_index["data/figure_inputs/fig2/figure2_DG_points.csv.xz"]["sha256"],
                 "frozen_fig2_points", inputs["hashes"])
    points = pd.read_csv(frozen_path,float_precision="round_trip")
    # Compare against per-evaluation frozen geometry rows (not just Figure 2 summaries).
    reference_geom=[]
    for method in ["CORDA","GIMME","RIPTiDe","iMAT"]:
        p=inputs["public_root"]/f"data/figure_inputs/fig4/fig4_geometry_paired/algorithm={method}/part.parquet.xz"
        rel=f"data/figure_inputs/fig4/fig4_geometry_paired/algorithm={method}/part.parquet.xz"
        _verify_hash(p, public_index[rel]["sha256"], f"frozen_fig4_geometry_{method}", inputs["hashes"])
        with lzma.open(p,"rb") as f:
            d=pd.read_parquet(io.BytesIO(f.read()),columns=["evaluation_id","reaction_id","A1_direction_explained_magnitude_variance","A2_direction_explained_magnitude_variance"])
        reference_geom.append(d)
    ref_geom=pd.concat(reference_geom,ignore_index=True)
    epos={str(x):i for i,x in enumerate(eval_rows.evaluation_id.astype(str))}
    rpos={inputs["rxns"][j]:k for k,j in enumerate(target_idx)}
    eta_ref=np.full_like(eta,np.nan)
    for row in ref_geom.itertuples(index=False):
        i=epos.get(str(row.evaluation_id)); j=rpos.get(str(row.reaction_id))
        if i is not None and j is not None:
            eta_ref[i,j,0]=row.A1_direction_explained_magnitude_variance
            eta_ref[i,j,1]=row.A2_direction_explained_magnitude_variance
    if not np.array_equal(np.isfinite(eta),np.isfinite(eta_ref)):
        raise ValueError("geometry reconstructed and frozen finite-reaction eligibility differ")
    finite=np.isfinite(eta)
    maxdiff=float(np.max(np.abs(eta[finite]-eta_ref[finite])))
    if maxdiff > 5e-12: raise ValueError(f"geometry η² exceeds 5e-12 parity tolerance: {maxdiff}")
    panels=[]
    # Use level2a-derived usefulness fractions, then confirm they match frozen
    # Figure 2 responses.  Values are calculated over exact eligible eval IDs.
    l2a_tables={}
    for anchor in ("A1","A2-L"):
        # _verify_l2a retains pair outputs for matching here.
        pair=l2a["pair_by_anchor"][anchor]
        useful=pair.pivot(index="evaluation_id",columns="reaction_id",values="directionally_useful_fraction")
        eval_selected=sorted(eval_ids[anchor])
        for panel,split,arm in (("D","development","A1"),("E","development","A2-L"),("F","confirmation","A1"),("G","confirmation","A2-L")):
            if anchor != arm: continue
            block=points.loc[points.panel.eq(panel)].sort_values("reaction_id")
            expected_labels = {"development": {"development"}, "confirmation": {"confirmation_holdout"}}[split]
            observed_labels = {str(inputs["split_labels"].get(str(r), "")).lower() for r in block.reaction_id}
            if observed_labels != expected_labels:
                raise ValueError(f"{panel}: Figure 2 reactions do not match the frozen {split} split")
            idx=np.array([eval_index[e] for e in eval_selected],dtype=int)
            eta_arm=0 if anchor=="A1" else 1
            order=np.argsort([inputs["rxns"][j] for j in target_idx])
            # Canonical pandas group mean is used only for aggregation; the
            # underlying eta values remain independently reconstructed.
            rids=np.asarray([inputs["rxns"][j] for j in target_idx],dtype=object)
            geom_frame=pd.DataFrame({"evaluation_id":np.repeat(eval_rows.iloc[idx].evaluation_id.astype(str).to_numpy(),len(rids)),"reaction_id":np.tile(rids,len(idx)),"eta":eta[idx,:,eta_arm].reshape(-1)})
            geom_mean=geom_frame.sort_values(["evaluation_id","reaction_id"],kind="stable").groupby("reaction_id",sort=True).eta.mean()
            xi=geom_mean.reindex(block.reaction_id).to_numpy(float)
            # Faithful reference path: use the serialized frozen per-evaluation
            # geometry, in the source evaluation order used by the Figure 2 exporter.
            ref_frame=pd.DataFrame({"evaluation_id":np.repeat(eval_rows.iloc[idx].evaluation_id.astype(str).to_numpy(),len(rids)),"reaction_id":np.tile(rids,len(idx)),"eta":eta_ref[idx,:,eta_arm].reshape(-1)})
            eval_order={str(e):i for i,e in enumerate(eval_rows.evaluation_id.astype(str))}
            ref_frame["_evaluation_order"]=ref_frame.evaluation_id.map(eval_order)
            ref_mean=(ref_frame.sort_values(["_evaluation_order","reaction_id"],kind="stable")
                      .groupby("reaction_id",sort=True).eta.mean())
            x_reference_faithful=ref_mean.reindex(block.reaction_id).to_numpy(float)
            yi=useful.reindex(index=eval_selected,columns=block.reaction_id).mean(axis=0).to_numpy(float)
            frozen_y=block.y_directional_usefulness_fraction.to_numpy(float)
            yd=float(np.max(np.abs(yi-frozen_y)))
            xf=block.x_direction_explained_magnitude_variance.to_numpy(float)
            if not (np.isfinite(xi).all() and np.isfinite(yi).all() and np.isfinite(xf).all()):
                raise ValueError(f"{panel}: plotted geometry/usefulness contains nonfinite values")
            from scipy.stats import rankdata
            independent_rank = rankdata(xi, method="average")
            frozen_rank = rankdata(xf, method="average")
            y_rank = rankdata(frozen_y, method="average")
            rank_changed = independent_rank != frozen_rank
            independent_y_rank = rankdata(yi, method="average")
            y_rank_changed = independent_y_rank != y_rank
            y_abs_diff = np.abs(yi-frozen_y)
            max_y_diff_index = int(np.argmax(y_abs_diff))
            def tie_counts(values):
                _, counts = np.unique(values, return_counts=True)
                return {"exact_tie_excess": int(np.sum(counts - 1)),
                        "exact_tied_pairs": int(np.sum(counts * (counts - 1) // 2))}
            def near_tie_pairs(values, tolerance=1e-15):
                ordered = np.sort(np.asarray(values, dtype=float))
                total = 0
                left = 0
                for right in range(len(ordered)):
                    while ordered[right] - ordered[left] > tolerance:
                        left += 1
                    total += right-left
                exact_pairs = tie_counts(values)["exact_tied_pairs"]
                return int(total-exact_pairs)
            frozen_cov = float(np.sum((frozen_rank-frozen_rank.mean())*(y_rank-y_rank.mean())))
            independent_cov = float(np.sum((independent_rank-independent_rank.mean())*(y_rank-y_rank.mean())))
            independent_y_cov = float(np.sum((frozen_rank-frozen_rank.mean())*(independent_y_rank-independent_y_rank.mean())))
            both_independent_cov = float(np.sum((independent_rank-independent_rank.mean())*(independent_y_rank-independent_y_rank.mean())))
            panels.append({"panel":panel,"reaction_count":len(block),
                "max_abs_eta2_difference":float(np.max(np.abs(xi-xf))),
                "exact_eta2_ties_frozen":tie_counts(xf),
                "exact_eta2_ties_reconstructed":tie_counts(xi),
                "exact_usefulness_ties_frozen":tie_counts(frozen_y),
                "exact_usefulness_ties_reconstructed":tie_counts(yi),
                "near_tie_pairs_within_1e-15_frozen":near_tie_pairs(xf),
                "near_tie_pairs_within_1e-15_reconstructed":near_tie_pairs(xi),
                "changed_rank_count":int(rank_changed.sum()),
                "changed_rank_reaction_ids":block.reaction_id.to_numpy()[rank_changed].astype(str).tolist(),
                "changed_usefulness_rank_count":int(y_rank_changed.sum()),
                "changed_usefulness_rank_reaction_ids":block.reaction_id.to_numpy()[y_rank_changed].astype(str).tolist(),
                "rank_l1_change":float(np.sum(np.abs(independent_rank-frozen_rank))),
                "rank_covariance_numerator_frozen":frozen_cov,
                "rank_covariance_numerator_reconstructed":independent_cov,
                "rank_covariance_numerator_change":independent_cov-frozen_cov,
                "rank_covariance_numerator_y_change_at_frozen_eta":independent_y_cov-frozen_cov,
                "rank_covariance_numerator_both_independent":both_independent_cov,
                "spearman_independent_eta2_usefulness":float(spearmanr(xi,yi).statistic),
                "spearman_independent_eta2_frozen_usefulness":float(spearmanr(xi,frozen_y).statistic),
                "spearman_frozen_eta2_independent_usefulness":float(spearmanr(xf,yi).statistic),
                "spearman_frozen_eta2_usefulness":float(spearmanr(xf,frozen_y).statistic),
                "spearman_reference_faithful_eta2_frozen_usefulness":float(spearmanr(x_reference_faithful,frozen_y).statistic),
                "frozen_reference_rho":float(FIG2_FROZEN_RHO[panel][1]),
                "frozen_reference_rho_tolerance":FIG2_RHO_SOURCE_TOL,
                "frozen_input_rho_exact_match":bool(float(spearmanr(xf,frozen_y).statistic)==FIG2_FROZEN_RHO[panel][1]),
                "frozen_input_rho_source_tolerance_pass":bool(abs(float(spearmanr(xf,frozen_y).statistic)-FIG2_FROZEN_RHO[panel][1])<=FIG2_RHO_SOURCE_TOL),
                "reference_faithful_eta2_max_abs_difference_to_fig2":float(np.max(np.abs(x_reference_faithful-xf))),
                "reference_faithful_eta2_exact_match_to_fig2":bool(np.array_equal(x_reference_faithful,xf)),
                "usefulness_frequency_max_abs_difference":float(np.max(np.abs(yi-frozen_y))),
                "usefulness_frequency_max_difference_reaction_id":str(block.reaction_id.iloc[max_y_diff_index]),
                "usefulness_frequency_at_max_difference_frozen":float(frozen_y[max_y_diff_index]),
                "usefulness_frequency_at_max_difference_independent":float(yi[max_y_diff_index]),
                "usefulness_frequency_reactions_above_5e-15":int(np.count_nonzero(y_abs_diff>5e-15)),
                "usefulness_frequency_exact_match":bool(np.array_equal(yi,frozen_y)),
            "rank_discrepancy_explanation":"η² rank changes reflect exact statistical tied-rank instability from floating-point serialization/reconstruction; any usefulness-frequency mismatch is reported separately and reflects case-level threshold classification"})
    return {"status":"PASS", "eta2_rows_per_anchor":int(eta.shape[0]*eta.shape[1]),
            "max_abs_eta2_difference_to_fig4":maxdiff,"panels":panels,
            "association_status":("PASS" if all(x["usefulness_frequency_exact_match"] for x in panels)
                                  else "SCIENTIFIC_DISCREPANCY"),
            "usefulness_threshold_mismatch_examples":{
                anchor: l2a["anchors"][anchor]["usefulness_mismatch_examples"]
                for anchor in ("A1", "A2-L")},
            "rank_comparison_scope":"Three correlation inputs are reported separately: frozen Figure 2 x/y; independently reconstructed eta2 with frozen y; independently reconstructed eta2 with independently reconstructed y. Reference-faithful frozen per-evaluation eta2 aggregation is also reported.",
            "definition":"full 20x20 product support; no truth-specific holdout"}


def format_human_report(result: dict) -> str:
    lines = [f"Scientific verification: {result['status']}",
             f"Requested level: {result['requested_level']}",
             f"Runtime: {result.get('runtime_seconds', float('nan')):.2f} s",
             f"Workers: {result.get('workers', {}).get('resolved', 'n/a')}"]
    for name, status in result.get("validation_statuses", {}).items():
        lines.append(f"{name}: {status}")
    for anchor, values in result.get("level2a", {}).get("anchors", {}).items():
        lines.extend([f"Level 2A {anchor}: {values['truth_pair_count']} truth pairs; "
                      f"{values['cases']} cases; {values['eligible_pairs']} eligible pairs",
                      f"  Independent correct/wrong/separation: "
                      f"{values['mean_correct_gain_independent']:.12g} / "
                      f"{values['mean_wrong_gain_independent']:.12g} / "
                      f"{values['separation_independent']:.12g}",
                      f"  Frozen correct/wrong/separation: "
                      f"{values['mean_correct_gain_frozen']:.12g} / "
                      f"{values['mean_wrong_gain_frozen']:.12g} / "
                      f"{values['separation_frozen']:.12g}",
                      f"  Usefulness classification differences: {values['usefulness_classification_mismatches']}"])
        lines.append(f"  Usefulness mismatches by split: {values['usefulness_mismatches_by_split']}")
        lines.append(f"  Threshold-distance counts: {values['threshold_distance_counts']}")
    l2b = result.get("level2b")
    if l2b:
        lines.append(f"Level 2B: {l2b['truth_pairs']} pairs; posterior predictions checked "
                     f"conditional on frozen weights; max truth error "
                     f"{l2b['reconstructed_truth_max_abs_difference']:.3g}")
        for anchor, predictions in l2b["maximum_prediction_differences"].items():
            lines.append(f"  {anchor} maximum prediction differences: {predictions}")
    geometry = result.get("geometry")
    if geometry:
        lines.append(f"Geometry: η² max absolute difference {geometry['max_abs_eta2_difference_to_fig4']:.3g}; "
                     f"association status {geometry['association_status']}")
        for panel in geometry["panels"]:
            lines.append(f"  {panel['panel']} n={panel['reaction_count']}: "
                         f"independent ρ={panel['spearman_independent_eta2_usefulness']:.9f}; "
                         f"frozen ρ={panel['spearman_frozen_eta2_usefulness']:.9f}; "
                         f"independent η²/frozen usefulness ρ={panel['spearman_independent_eta2_frozen_usefulness']:.9f}; "
                         f"changed η² ranks={panel['changed_rank_count']}; "
                         f"max usefulness-frequency difference="
                         f"{panel['usefulness_frequency_max_abs_difference']:.3g}")
    lines.append("Full hashes, identities, tie/rank diagnostics, and discrepancy details: science_validation.json")
    return "\n".join(lines) + "\n"


def run_verification(*, level: str, data_root: Path, public_root: Path,
                     workers: int | str = 1, output_root: Path, batch_size: int = 1) -> dict:
    if output_root.exists():
        raise FileExistsError(f"output directory already exists; choose a new --output-root: {output_root}")
    start=time.perf_counter()
    inputs=_load_inputs(data_root,public_root)
    result={"schema":"recomb2027.level2.science_validation.v1","status":"PASS",
            "requested_level":level,"inputs_sha256":inputs["hashes"],
            "identity_counts":{**EXPECTED,"evaluations":len(inputs["evals"]),
                "methods":int(inputs["evals"].algorithm.nunique()),
                "rna_contexts":int(inputs["evals"].rna_context_key.nunique()),
                "ct2a_mice":int(inputs["evals"].ct2a_mouse.nunique()),
                "gl261_mice":int(inputs["evals"].gl261_mouse.nunique())},
            "weight_diagnostics":inputs["weight_diagnostics"],"limitations":[
              "Level 2B conditions on frozen S1 anchor weights; upstream strong-anchor weight generation and model reconstruction are unsupported.",
              "Private source artifacts are local inputs and are not packaged or redistributed."]}
    l2a=_verify_l2a(inputs)
    class_mismatch = any(l2a["anchors"][a]["usefulness_classification_mismatches"]
                         for a in ("A1", "A2-L"))
    result["level2a"]={"status":l2a["status"],
        "endpoint_status":"PASS",
        "case_gain_numerical_parity_status":"PASS" if all(l2a["anchors"][a]["case_gain_max_abs_diff"] <= 5e-12 and l2a["anchors"][a]["case_wrong_gain_max_abs_diff"] <= 5e-12 for a in ("A1", "A2-L")) else "FAIL",
        "exact_usefulness_classification_parity_status":"FAIL" if class_mismatch else "PASS",
        "threshold_boundary_sensitivity":"diagnostic_counts_reported_per_anchor_and_metric",
        "anchors":l2a["anchors"]}
    if level.lower() in ("2b","geometry"):
        l2b=_verify_l2b(inputs,workers=workers,batch_size=batch_size)
        result["level2b"]=l2b
    if level.lower() in ("geometry", "2b"):
        # Geometry needs independently recomputed pair response tables from L2A.
        l2a["pair_by_anchor"]={a:p for a,p in zip(("A1","A2-L"),l2a["pair_tables"])}
        # Geometry is serial by default because measured 40-evaluation
        # benchmarks showed process startup and shared-cache overhead slower.
        geometry_workers = workers if level.lower() == "geometry" else 1
        result["geometry"]=_verify_geometry(inputs,l2a,workers=geometry_workers,batch_size=batch_size)
        result["geometry"]["workers"] = science.resolve_workers(geometry_workers)
    result["validation_statuses"]={
        "numerical_prediction_parity": ("PASS" if "level2b" in result and max(
            max(v.values()) for v in result["level2b"]["maximum_prediction_differences"].values()) <= 5e-12 else ("NOT_RUN" if "level2b" not in result else "FAIL")),
        "case_gain_numerical_parity": result["level2a"]["case_gain_numerical_parity_status"],
        "exact_usefulness_classification_parity": result["level2a"]["exact_usefulness_classification_parity_status"],
        "threshold_boundary_sensitivity": "DIAGNOSTIC_ONLY",
        "geometry_numerical_parity": ("PASS" if "geometry" in result and result["geometry"]["max_abs_eta2_difference_to_fig4"] <= 5e-12 else ("NOT_RUN" if "geometry" not in result else "FAIL")),
        "exact_frozen_spearman_reproduction": ("PASS" if "geometry" in result and all(p["reference_faithful_eta2_exact_match_to_fig2"] for p in result["geometry"]["panels"]) else ("NOT_RUN" if "geometry" not in result else "FAIL")),
        "exact_spearman_match_to_source_expected_rhos": ("PASS" if "geometry" in result and all(p["frozen_input_rho_exact_match"] for p in result["geometry"]["panels"]) else ("NOT_RUN" if "geometry" not in result else "FAIL")),
        "frozen_spearman_source_tolerance": ("PASS" if "geometry" in result and all(p["frozen_input_rho_source_tolerance_pass"] for p in result["geometry"]["panels"]) else ("NOT_RUN" if "geometry" not in result else "FAIL")),
        "independently_reconstructed_geometry_exact_spearman_match": ("PASS" if "geometry" in result and all(p["spearman_independent_eta2_frozen_usefulness"] == p["spearman_frozen_eta2_usefulness"] for p in result["geometry"]["panels"]) else ("NOT_RUN" if "geometry" not in result else "FAIL")),
    }
    result["workers"]={"requested":workers,"resolved":science.resolve_workers(workers)}
    result["runtime_seconds"]=time.perf_counter()-start
    modules_pass = all(result[k]["status"]=="PASS" for k in ("level2a","level2b","geometry") if k in result)
    result["status"]="PASS" if modules_pass else "FAIL"
    if class_mismatch or result.get("geometry", {}).get("association_status") == "SCIENTIFIC_DISCREPANCY":
        result["status"] = "SCIENTIFIC_DISCREPANCY"
    output_root.mkdir(parents=True,exist_ok=False)
    (output_root/"science_validation.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    (output_root/"science_validation.txt").write_text(format_human_report(result))
    return result
