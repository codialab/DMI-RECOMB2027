#!/usr/bin/env python3
"""Additive uncertainty analysis for the frozen RECOMB directional-cue study.

Frozen artifacts may be staged under .cache/recomb2027-dmi/frozen with their
repository-relative paths. The empirical command refuses to run unless every
required artifact is present and matches its pinned SHA-256.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import lzma
import math
import os
from pathlib import Path
import sys
import time
from collections import Counter, defaultdict

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
RESTORE_ROOT = ROOT / ".cache/recomb2027-dmi/frozen"
OUT = ROOT / "outputs/dmi_bridge_uncertainty_reduction_v1"
LINEAGE = ROOT / "manifests/provenance/FROZEN_LINEAGE.json"
EXPECTED = {
    "candidate_table": ("12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/candidate_fraction_table.tsv.xz", "421c2c1d196fbd088a4e92e2ecc998ad578aab814a025fff8e1c843597fed1da"),
    "flux_cache": ("outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstructed/flux_cache.npz", "f9644da4ebd875cf5f0b34cc58b20faf0a0457d7c7a44635e8a010d8896c3dbb"),
    "a1_weights": ("outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_FOUR_ARM_WEIGHTS.tsv.xz", "843e17b7012a6fa212d1e9ac9e1dc6e309e854d94ab5c4b3f9e316b589768f4f"),
    "a20_manifest": ("outputs/dmi_bridge_a20_dual_anchor_qualification_v1/BRIDGEA20_MANIFEST.json", "01f8004a0bf9d5358d7ce7281e827cc1a1fbe00a6750a5c063a1edad2eedf380"),
    "pl1_manifest": ("outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_MANIFEST.json", "b0f81620ff24f5791dd964a8063416784db86a358b0464b8c68a504e84e5f529"),
    "pl2a_manifest": ("outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_MANIFEST.json", "8ae3eb60eb23614e379e38e3aa3cb72ec16c6bdb69aa9c74cc94499e31f1c312"),
    "pl2b_manifest": ("outputs/dmi_bridge_pl2b_sign_only_utility_v1/BRIDGEPL2B_MANIFEST.json", "7d9ca97a3f4ea5aa8b0b983e09ee83fc3c82f326d1b1cbf9c6137cfadf25539d"),
    "a20_weights": ("outputs/dmi_bridge_a20_dual_anchor_qualification_v1/BRIDGEA20_DUAL_ANCHOR_WEIGHTS.tsv.xz", "b153e2fde45dca9573c71c548783e704313f056d0980f64af60180ad90ebc1f4"),
}
OUTPUT_SCHEMA = "bridge.uncertainty_reduction.v1"
CASE_FIELDS = (
    "anchor", "truth_pair_id", "evaluation_id", "algorithm", "rna_context_key",
    "ct2a_mouse", "gl261_mouse", "reaction_id", "truth_status", "truth_delta",
    "truth_magnitude", "baseline_mean_magnitude", "baseline_variance", "baseline_sd",
    "baseline_iqr", "baseline_central90_width", "baseline_crps",
    "relative_variance_eligible", "relative_variance_threshold", "support_count",
    "prediction_gain_correct", "prediction_gain_incorrect", "prediction_status_correct",
    "prediction_status_incorrect", "directionally_useful",
) + tuple(
    f"{metric}_{cue}_lambda_{lam:g}"
    for lam in (0.0, 0.25, 0.5, 1.0)
    for cue in ("correct", "incorrect")
    for metric in ("variance", "sd", "iqr", "central90_width", "absolute_variance_reduction", "relative_variance_reduction", "crps", "crps_improvement", "variance_change_status")
) + tuple(
    f"prediction_gain_{cue}_lambda_{lam:g}"
    for lam in (0.0, 0.25, 0.5, 1.0)
    for cue in ("correct", "incorrect")
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _candidate_paths(relative: str):
    yield ROOT / relative
    yield RESTORE_ROOT / relative


def resolve_file(relative: str, expected_hash: str | None = None) -> tuple[Path | None, str | None]:
    """Resolve a regular repository-contained input and optionally verify it."""
    errors = []
    for path in _candidate_paths(relative):
        if not path.exists():
            continue
        if path.is_symlink():
            errors.append(f"symlink rejected: {path}")
            continue
        resolved = path.resolve(strict=True)
        if not resolved.is_relative_to(ROOT) or not resolved.is_file():
            errors.append(f"path outside repository or not regular: {path}")
            continue
        digest = sha256_file(resolved)
        if expected_hash is not None and digest != expected_hash:
            errors.append(f"SHA-256 mismatch: {relative} expected {expected_hash} got {digest}")
            continue
        return resolved, None
    return None, "; ".join(errors) if errors else None


def required_artifacts() -> dict[str, str]:
    """Return every frozen artifact hash needed by the empirical pipeline."""
    requirements = {path: digest for path, digest in EXPECTED.values()}
    lineage = json.loads(LINEAGE.read_text(encoding="utf-8"))
    for stage in (
        "dmi_bridge_a20_dual_anchor_qualification_v1",
        "dmi_bridge_pl1_predictability_landscape_v1_resumable_v1",
        "dmi_bridge_pl2a_sign_only_prepare_v1",
        "dmi_bridge_pl2b_sign_only_utility_v1",
    ):
        record = lineage[stage]
        requirements[record["source_path"]] = record["source_sha256"]
        for name, digest in record.get("artifact_sha256", {}).items():
            requirements[str(Path(record["source_path"]).parent / name)] = digest
    # These frozen output bundles are already in this publication checkout and
    # bind the A2 cases, A2 geometry, and reaction split used by the analysis.
    for relative in (
        "outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_MANIFEST.json",
        "outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_MANIFEST.json",
        "outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_REACTION_SPLIT.tsv",
    ):
        if relative.endswith("BRIDGEA21_MANIFEST.json"):
            requirements[relative] = "db9948d7b99447a695cf8aa21429dde6d642e1fd0632bcce52643deac59e6feb"
        elif relative.endswith("BRIDGEA22_MANIFEST.json"):
            requirements[relative] = "15b237151005538e5077a352868f873ef00466f4f5ea49307b781a8a8cd0e7cb"
        else:
            requirements[relative] = "7dd8676508416a13f88af184110ef40e0f8d92782cd5d62b42e7006af5316e2f"
    # Pin every file that the empirical run later consumes from these two
    # already-published bundles, rather than treating a valid manifest alone
    # as sufficient preflight evidence.
    for relative, digest in (
        ("outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_MANIFEST.json", requirements["outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_MANIFEST.json"]),
        ("outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_MANIFEST.json", requirements["outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_MANIFEST.json"]),
    ):
        path = ROOT / relative
        if path.is_file() and not path.is_symlink() and sha256_file(path) == digest:
            bundle = json.loads(path.read_text(encoding="utf-8"))
            base = str(Path(relative).parent)
            for name, artifact_digest in bundle.get("artifact_sha256", {}).items():
                requirements[str(Path(base) / name)] = artifact_digest
    return requirements


def preflight() -> dict[str, object]:
    requirements = required_artifacts()
    available: dict[str, str] = {}
    missing: list[dict[str, str]] = []
    mismatched: list[dict[str, str]] = []
    for relative, expected_hash in sorted(requirements.items()):
        path, error = resolve_file(relative, expected_hash)
        if path is not None:
            available[relative] = expected_hash
        elif error and "SHA-256 mismatch" in error:
            mismatched.append({"path": relative, "expected_sha256": expected_hash, "detail": error})
        else:
            missing.append({"path": relative, "expected_sha256": expected_hash})
    return {
        "schema": OUTPUT_SCHEMA,
        "status": "READY" if not missing and not mismatched else "BLOCKED_MISSING_FROZEN_INPUTS",
        "repository_root": str(ROOT),
        "restore_root": str(RESTORE_ROOT),
        "required_count": len(requirements),
        "available_count": len(available),
        "missing": missing,
        "mismatched": mismatched,
        "available_sha256": available,
    }


def _open_text(path: Path):
    if path.suffix == ".xz":
        return lzma.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with _open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _manifest_bundle(stage: str, manifest_relative: str, expected_manifest_sha: str):
    path, error = resolve_file(manifest_relative, expected_manifest_sha)
    if path is None:
        raise RuntimeError(f"frozen input unavailable: {manifest_relative}: {error or 'missing'}")
    manifest = json.loads(path.read_text(encoding="utf-8"))
    base = str(Path(manifest_relative).parent)
    artifact_paths: dict[str, Path] = {}
    for name, digest in manifest.get("artifact_sha256", {}).items():
        relative = str(Path(base) / name)
        resolved, problem = resolve_file(relative, digest)
        if resolved is None:
            raise RuntimeError(f"{stage} artifact unavailable or hash mismatch: {relative}: {problem or 'missing'}")
        artifact_paths[name] = resolved
    return manifest, artifact_paths


def _support_specs_a2(pairs, candidate_by_id, weights_rows):
    from scripts import dmi_bridge_a22_dual_anchor_sign_only_core_v1 as a22
    pools = defaultdict(dict)
    for row in weights_rows:
        if row["arm"] != "A2-L":
            raise RuntimeError("A20 bundle contains unexpected arm")
        key = (row["mouse_id"], row["tumor"], row["algorithm"], row["ensemble_hash"], row["projection_hash"], row["rna_context_key"])
        sample = int(row["sample_index"])
        if sample in pools[key]:
            raise RuntimeError("duplicate A2-L sample weight")
        pools[key][sample] = float(row["weight"])
    specs = []
    for pair in pairs:
        if not pair["pair_evaluable"]:
            continue
        condition = {}
        for prefix, tumor in (("ct2a", "CT2A"), ("gl261", "GL261")):
            candidate = candidate_by_id[str(pair[f"{prefix}_candidate_id"])]
            key = (str(pair[f"{prefix}_mouse"]), tumor, pair["algorithm"], candidate["ensemble_hash"], candidate["projection_hash"], pair["rna_context_key"])
            pool = pools[key]
            samples = sorted(pool)
            holdout = int(pair[f"{prefix}_sample_index"])
            loc = samples.index(holdout)
            held = a22.remove_holdout_and_normalize([pool[s] for s in samples], loc)
            condition[prefix] = {
                "held_samples": np.asarray(samples[:loc] + samples[loc + 1:], dtype=int),
                "held_weights": held.weights,
                "truth_sample": holdout,
                "cache_key": (candidate["algorithm"], tumor, candidate["ensemble_hash"]),
            }
        specs.append({"pair": pair, **condition})
    specs.sort(key=lambda value: str(value["pair"]["truth_pair_id"]))
    if len(specs) != 836:
        raise RuntimeError(f"A2 truth support count mismatch: {len(specs)}")
    return specs


def _case_part_reader(parts, index):
    part_path = parts[index]
    handle = _open_text(part_path)
    return handle, csv.DictReader(handle, delimiter="\t")


def _load_existing_case_part(path, expected_rows, accumulator):
    if not path.is_file() or path.is_symlink():
        return False
    with _open_text(path) as handle:
        count = 0
        for row in csv.DictReader(handle, delimiter="\t"):
            pair = {key: row[key] for key in (
                "truth_pair_id", "evaluation_id", "algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse"
            )}
            row["relative_variance_eligible"] = row["relative_variance_eligible"].lower() == "true"
            row["directionally_useful"] = row["directionally_useful"].lower() == "true"
            accumulator.add(pair, row["reaction_id"], row)
            count += 1
    if count != expected_rows:
        raise RuntimeError(f"partial output part has {count} rows, expected {expected_rows}: {path}")
    return True


def _expected_cache_row(rows, anchor, reaction_id):
    row = next(rows, None)
    while anchor == "A1" and row is not None and row["reaction_id"] == "LDH_L":
        row = next(rows, None)
    if row is None:
        raise RuntimeError("cached frozen case outcomes ended early")
    return row


class CaseAccumulator:
    """Compact evaluation-by-reaction sums for equal-weight pair aggregation."""

    UNIT_FIELDS = (
        "relative_correct", "relative_incorrect", "absolute_correct", "absolute_incorrect",
        "gain_correct", "gain_incorrect", "crps_correct", "crps_incorrect", "usefulness",
        "baseline_sd", "baseline_iqr", "baseline_central90_width",
        "sd_correct", "sd_incorrect", "iqr_correct", "iqr_incorrect",
        "central90_correct", "central90_incorrect",
    )

    def __init__(self, eval_ids, reaction_ids):
        self.eval_index = {key: i for i, key in enumerate(eval_ids)}
        self.reaction_index = {key: i for i, key in enumerate(reaction_ids)}
        shape = (len(eval_ids), len(reaction_ids))
        self.sums = {name: np.zeros(shape, dtype=np.float64) for name in self.UNIT_FIELDS}
        self.counts = {name: np.zeros(shape, dtype=np.uint16) for name in self.UNIT_FIELDS}
        self.truth_counts = np.zeros(shape, dtype=np.uint16)
        self.case_summaries = Counter()
        self.stratum_summaries = Counter()

    def add(self, pair, reaction_id, row):
        self.case_summaries[("case_rows",)] += 1
        ei = self.eval_index[str(pair["evaluation_id"])]
        ri = self.reaction_index[reaction_id]
        if row["truth_status"] == "TRUTH_TIE":
            self.case_summaries[("truth_ties",)] += 1
            return
        self.truth_counts[ei, ri] += 1
        self.case_summaries[("non_ties",)] += 1
        eligible = row["relative_variance_eligible"]
        if isinstance(eligible, str):
            eligible = eligible.strip().lower() in {"true", "1"}
        if not bool(eligible):
            self.case_summaries[("relative_variance_excluded_near_zero",)] += 1
        for cue in ("correct", "incorrect"):
            gain = float(row[f"prediction_gain_{cue}"])
            gain_status = "IMPROVED" if gain > 1e-12 else ("WORSENED" if gain < -1e-12 else "TIED")
            uncertainty = row[f"variance_change_status_{cue}_lambda_0.25"]
            self.case_summaries[(cue, uncertainty)] += 1
            self.case_summaries[(cue, uncertainty, gain_status)] += 1
            for lam in (0.0, 0.25, 0.5, 1.0):
                lam_key = f"{lam:g}"
                status = row[f"variance_change_status_{cue}_lambda_{lam_key}"]
                if status and status != "nan":
                    self.case_summaries[(cue, f"lambda_{lam_key}", status)] += 1
                gain_value = float(row[f"prediction_gain_{cue}_lambda_{lam_key}"])
                sensitivity_gain_status = "IMPROVED" if gain_value > 1e-12 else ("WORSENED" if gain_value < -1e-12 else "TIED")
                self.case_summaries[(cue, f"lambda_{lam_key}", "prediction_gain", sensitivity_gain_status)] += 1
                if status and status != "nan":
                    self.stratum_summaries[(str(pair["algorithm"]), str(pair["rna_context_key"]), cue,
                                            f"lambda_{lam_key}", status, sensitivity_gain_status)] += 1
        for field, source in (
            ("relative_correct", "relative_variance_reduction_correct_lambda_0.25"),
            ("relative_incorrect", "relative_variance_reduction_incorrect_lambda_0.25"),
            ("absolute_correct", "absolute_variance_reduction_correct_lambda_0.25"),
            ("absolute_incorrect", "absolute_variance_reduction_incorrect_lambda_0.25"),
            ("gain_correct", "prediction_gain_correct"),
            ("gain_incorrect", "prediction_gain_incorrect"),
            ("crps_correct", "crps_improvement_correct_lambda_0.25"),
            ("crps_incorrect", "crps_improvement_incorrect_lambda_0.25"),
            ("usefulness", "directionally_useful"),
            ("baseline_sd", "baseline_sd"),
            ("baseline_iqr", "baseline_iqr"),
            ("baseline_central90_width", "baseline_central90_width"),
            ("sd_correct", "sd_correct_lambda_0.25"),
            ("sd_incorrect", "sd_incorrect_lambda_0.25"),
            ("iqr_correct", "iqr_correct_lambda_0.25"),
            ("iqr_incorrect", "iqr_incorrect_lambda_0.25"),
            ("central90_correct", "central90_width_correct_lambda_0.25"),
            ("central90_incorrect", "central90_width_incorrect_lambda_0.25"),
        ):
            value = row[source]
            try:
                value = float(value)
            except (TypeError, ValueError):
                continue
            if not math.isfinite(value):
                continue
            self.sums[field][ei, ri] += value
            self.counts[field][ei, ri] += 1

    def reaction_rows(self, anchor, features, split_map):
        output = []
        for reaction, ri in self.reaction_index.items():
            row = {"anchor": anchor, "reaction_id": reaction, "analysis_split": split_map[reaction]}
            for field in self.UNIT_FIELDS:
                vals = np.divide(self.sums[field][:, ri], self.counts[field][:, ri], out=np.full(len(self.eval_index), np.nan), where=self.counts[field][:, ri] > 0)
                finite = np.isfinite(vals)
                row[field] = float(np.mean(vals[finite])) if np.any(finite) else math.nan
            row.update(features.get(reaction, {}))
            output.append(row)
        return output


def _case_cache_values(cache_row, anchor):
    prefixes = {0.25: "lambda_0_25", 0.5: "lambda_0_5", 1.0: "lambda_1_0"}
    gain_by_lambda = {}
    for lam, prefix in prefixes.items():
        gain_by_lambda[(lam, "correct")] = float(cache_row[f"{prefix}_correct_absolute_error_gain"])
        gain_by_lambda[(lam, "incorrect")] = float(cache_row[f"{prefix}_wrong_absolute_error_gain"])
    gain_by_lambda[(0.0, "correct")] = gain_by_lambda[(0.0, "incorrect")] = 0.0
    return {
        "baseline": float(cache_row["baseline_magnitude_estimate"]),
        "correct_gain": gain_by_lambda[(0.25, "correct")],
        "incorrect_gain": gain_by_lambda[(0.25, "incorrect")],
        "random_gain": float(cache_row["lambda_0_25_random_sign_expected_gain"]),
        "gain_by_lambda": gain_by_lambda,
        "truth_status": cache_row["truth_status"],
    }


def _support_cache_position(support, cache_index):
    """Accept frozen A1 positions and A2 keys without mutating their schemas."""
    if "cache_index" in support:
        return int(support["cache_index"])
    return int(cache_index[support["cache_key"]])


def _analyze_anchor(anchor, specs, reaction_ids, rxn_index, mats, cache_index,
                    cached_parts, output_dir, evaluation_ids):
    from scripts import dmi_bridge_uncertainty_reduction_core_v1 as core
    from scripts import dmi_bridge_pl2_sign_only_core_v1 as pl2
    from scripts import dmi_bridge_pl2c_geometry_utility_core_v1 as pl2c

    part_manifest = json.loads(cached_parts["manifest_path"].read_text(encoding="utf-8"))
    cached_paths = [cached_parts["artifacts"][rec["filename"]] for rec in part_manifest["parts"]]
    accumulator = CaseAccumulator(evaluation_ids, reaction_ids)
    pair_ranges = [((len(specs) * i) // len(cached_paths), (len(specs) * (i + 1)) // len(cached_paths)) for i in range(len(cached_paths))]
    case_records = []
    illustration_written = (output_dir / "ILLUSTRATIVE_DISTRIBUTION.tsv").exists()
    max_baseline_dev, max_gain_dev, total_cases = 0.0, 0.0, 0
    for part_idx, (first, stop) in enumerate(pair_ranges):
        final_path = output_dir / f"CASE_METRICS_{anchor}.part-{part_idx:03d}.tsv.xz"
        checkpoint_path = output_dir / "RUN_STATE.json"
        checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        checkpoint_record = checkpoint.get("completed_parts", {}).get(final_path.name)
        if final_path.exists():
            if not checkpoint_record or checkpoint_record.get("sha256") != sha256_file(final_path):
                raise RuntimeError(f"existing output part lacks a matching checkpoint: {final_path}")
            expected_rows = (stop - first) * len(reaction_ids)
            _load_existing_case_part(final_path, expected_rows, accumulator)
            cached_handle, _ = _case_part_reader(cached_paths, part_idx)
            cached_handle.close()
            case_records.append({"path": final_path, "sha256": checkpoint_record["sha256"], "rows": expected_rows})
            max_baseline_dev = max(max_baseline_dev, float(checkpoint_record.get("maximum_baseline_mean_deviation", 0.0)))
            max_gain_dev = max(max_gain_dev, float(checkpoint_record.get("maximum_prediction_gain_deviation", 0.0)))
            total_cases += expected_rows
            print(f"{anchor}: resumed partition {part_idx + 1}/{len(pair_ranges)}", flush=True)
            continue
        cached_handle, cached_rows = _case_part_reader(cached_paths, part_idx)
        temp_path = output_dir / f"CASE_METRICS_{anchor}.part-{part_idx:03d}.tsv.xz.tmp"
        with lzma.open(temp_path, "wt", encoding="utf-8", newline="", preset=3) as handle:
            writer = csv.DictWriter(handle, fieldnames=CASE_FIELDS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            for spec in specs[first:stop]:
                pair = spec["pair"]
                left, right = spec["ct2a"], spec["gl261"]
                li = _support_cache_position(left, cache_index)
                ri = _support_cache_position(right, cache_index)
                left_values = mats[li, left["held_samples"]]
                right_values = mats[ri, right["held_samples"]]
                truth_left = int(left["truth_sample"])
                truth_right = int(right["truth_sample"])
                for reaction in reaction_ids:
                    cached = _expected_cache_row(cached_rows, anchor, reaction)
                    if cached["truth_pair_id"] != str(pair["truth_pair_id"]) or cached["reaction_id"] != reaction:
                        raise RuntimeError(f"cached case order/key mismatch for {anchor}: {pair['truth_pair_id']} {reaction}")
                    col = rxn_index[reaction]
                    delta = (left_values[:, col, None] - right_values[None, :, col]).reshape(-1)
                    joint = (left["held_weights"][:, None] * right["held_weights"][None, :]).reshape(-1)
                    truth_delta = float(mats[li, truth_left, col] - mats[ri, truth_right, col])
                    dist = core.analyze_distribution(delta, joint, truth_delta)
                    if anchor == "A1" and not illustration_written and dist["truth_direction"] != 0:
                        order = np.argsort(np.abs(delta), kind="mergesort")
                        sign = np.where(delta > 1e-12, 1, np.where(delta < -1e-12, -1, 0))
                        cue_weights = {}
                        for cue, direction in (("correct", dist["truth_direction"]), ("incorrect", -dist["truth_direction"])):
                            posterior = joint * np.exp(0.25 * direction * sign)
                            cue_weights[cue] = posterior / posterior.sum()
                        _write_tsv(output_dir / "ILLUSTRATIVE_DISTRIBUTION.tsv", [
                            {"magnitude": float(abs(delta[idx])), "baseline_weight": float(joint[idx]),
                             "correct_lambda_0_25_weight": float(cue_weights["correct"][idx]),
                             "incorrect_lambda_0_25_weight": float(cue_weights["incorrect"][idx]),
                             "truth_magnitude": abs(truth_delta), "truth_pair_id": pair["truth_pair_id"],
                             "reaction_id": reaction}
                            for idx in order
                        ])
                        illustration_written = True
                    prediction = pl2.evaluate_truth_case(delta_b=delta, baseline_weights=joint, truth_delta_b=truth_delta)
                    old = _case_cache_values(cached, anchor)
                    baseline_dev = abs(float(dist["baseline"]["mean"]) - old["baseline"])
                    max_baseline_dev = max(max_baseline_dev, baseline_dev)
                    if baseline_dev > 5e-12:
                        raise RuntimeError(f"frozen baseline mean parity failure: {baseline_dev}")
                    truth_status = "TRUTH_TIE" if dist["truth_direction"] == 0 else "NON_TIE"
                    if old["truth_status"] != truth_status:
                        raise RuntimeError("frozen truth tie status mismatch")
                    pred_by_lam = {format(float(row["lambda_total"]), ".17g"): row for row in prediction["lambda_results"]}
                    row = {
                        "anchor": anchor,
                        "truth_pair_id": pair["truth_pair_id"],
                        "evaluation_id": pair["evaluation_id"],
                        "algorithm": pair["algorithm"],
                        "rna_context_key": pair["rna_context_key"],
                        "ct2a_mouse": pair["ct2a_mouse"],
                        "gl261_mouse": pair["gl261_mouse"],
                        "reaction_id": reaction,
                        "truth_status": truth_status,
                        "truth_delta": truth_delta,
                        "truth_magnitude": abs(truth_delta),
                        "baseline_mean_magnitude": dist["baseline"]["mean"],
                        "baseline_variance": dist["baseline"]["variance"],
                        "baseline_sd": dist["baseline"]["sd"],
                        "baseline_iqr": dist["baseline"]["iqr"],
                        "baseline_central90_width": dist["baseline"]["central90_width"],
                        "baseline_crps": dist["baseline"]["crps"],
                        "relative_variance_eligible": dist["relative_variance_eligible"],
                        "relative_variance_threshold": dist["relative_variance_threshold"],
                        "support_count": dist["support_count"],
                        "prediction_gain_correct": old["correct_gain"],
                        "prediction_gain_incorrect": old["incorrect_gain"],
                        "prediction_status_correct": pl2.classify_gain(old["correct_gain"]),
                        "prediction_status_incorrect": pl2.classify_gain(old["incorrect_gain"]),
                        "directionally_useful": pl2c.case_responses(
                            correct_gain=old["correct_gain"], wrong_gain=old["incorrect_gain"],
                            random_gain=old["random_gain"], truth_status=truth_status,
                        )["directionally_useful"],
                    }
                    for lam in (0.0, 0.25, 0.5, 1.0):
                        lam_key = format(lam, ".17g")
                        for cue in ("correct", "incorrect"):
                            suffix = "correct" if cue == "correct" else "wrong"
                            distribution = dist["updates"].get(lam_key, {}).get(cue)
                            if distribution is None:
                                for metric in ("variance", "sd", "iqr", "central90_width", "absolute_variance_reduction", "relative_variance_reduction", "crps", "crps_improvement", "variance_change_status"):
                                    row[f"{metric}_{cue}_lambda_{lam:g}"] = math.nan
                                continue
                            for metric in ("variance", "sd", "iqr", "central90_width", "absolute_variance_reduction", "relative_variance_reduction", "crps", "crps_improvement", "variance_change_status"):
                                row[f"{metric}_{cue}_lambda_{lam:g}"] = distribution.get(metric, math.nan)
                            if truth_status == "NON_TIE" and lam == 0.25:
                                predicted_gain = float(pred_by_lam[lam_key][f"{suffix}_abs_error_gain"])
                                max_gain_dev = max(max_gain_dev, abs(predicted_gain - old[f"{cue}_gain"]))
                                if abs(predicted_gain - old[f"{cue}_gain"]) > 5e-12:
                                    raise RuntimeError("frozen prediction gain parity failure")
                            if truth_status == "NON_TIE":
                                predicted_gain = float(pred_by_lam[lam_key][f"{suffix}_abs_error_gain"])
                                cached_gain = float(old["gain_by_lambda"][(lam, cue)])
                                if abs(predicted_gain - cached_gain) > 5e-12:
                                    raise RuntimeError(f"frozen prediction gain parity failure at lambda={lam:g}")
                                row[f"prediction_gain_{cue}_lambda_{lam:g}"] = predicted_gain
                            else:
                                row[f"prediction_gain_{cue}_lambda_{lam:g}"] = math.nan
                    accumulator.add(pair, reaction, row)
                    writer.writerow(row)
                    total_cases += 1
        cached_handle.close()
        os.replace(temp_path, final_path)
        case_records.append({"path": final_path, "sha256": sha256_file(final_path), "rows": (stop - first) * len(reaction_ids)})
        checkpoint["completed_parts"][final_path.name] = {
            "sha256": case_records[-1]["sha256"], "rows": case_records[-1]["rows"],
            "maximum_baseline_mean_deviation": max_baseline_dev,
            "maximum_prediction_gain_deviation": max_gain_dev,
        }
        _write_json_atomic(checkpoint_path, checkpoint)
        print(f"{anchor}: partition {part_idx + 1}/{len(pair_ranges)} complete", flush=True)
    return accumulator, case_records, {"case_rows": total_cases, "maximum_baseline_mean_deviation": max_baseline_dev, "maximum_prediction_gain_deviation": max_gain_dev}


def _feature_means(parts, reaction_ids, fields):
    sums = {reaction: {field: 0.0 for field in fields} for reaction in reaction_ids}
    counts = {reaction: {field: 0 for field in fields} for reaction in reaction_ids}
    for path in parts:
        with _open_text(path) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                reaction = row["reaction_id"]
                if reaction not in sums:
                    continue
                for field in fields:
                    value = float(row.get(field, "nan"))
                    if math.isfinite(value):
                        sums[reaction][field] += value
                        counts[reaction][field] += 1
                negative = float(row.get("p_negative", "nan"))
                positive = float(row.get("p_positive", "nan"))
                if math.isfinite(negative) and math.isfinite(positive):
                    sums[reaction].setdefault("non_tie_coverage", 0.0)
                    counts[reaction].setdefault("non_tie_coverage", 0)
                    sums[reaction]["non_tie_coverage"] += negative + positive
                    counts[reaction]["non_tie_coverage"] += 1
    result = {}
    for reaction in reaction_ids:
        result[reaction] = {
            "sign_magnitude_eta2": _finite_mean(sums[reaction]["sign_magnitude_eta2"], counts[reaction]["sign_magnitude_eta2"]),
            "directional_entropy3": _finite_mean(sums[reaction]["directional_entropy3"], counts[reaction]["directional_entropy3"]),
            "dominant_sign_mass": _finite_mean(sums[reaction]["dominant_sign_mass"], counts[reaction]["dominant_sign_mass"]),
            "n_supported_sign_states": _finite_mean(sums[reaction]["n_supported_sign_states"], counts[reaction]["n_supported_sign_states"]),
            "non_tie_coverage": _finite_mean(sums[reaction].get("non_tie_coverage", 0.0), counts[reaction].get("non_tie_coverage", 0)),
        }
    return result


def _finite_mean(total, n):
    return float(total / n) if n else math.nan


def _feature_parts(manifest, artifact_paths, prefix):
    info_name = f"{prefix}_REACTION_EVALUATION_FEATURES.parts.json"
    info_path = artifact_paths[info_name]
    info = json.loads(info_path.read_text(encoding="utf-8"))
    return [artifact_paths[part["filename"]] for part in info["parts"]]


def _save_unit_table(anchor, accumulator, eval_rows, reaction_ids, path):
    with lzma.open(path, "wt", encoding="utf-8", newline="", preset=3) as handle:
        fields = ("anchor", "evaluation_id", "reaction_id", "algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse", "n_non_tie_pairs", *accumulator.UNIT_FIELDS)
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for eval_index, evaluation in enumerate(eval_rows):
            for reaction_index, reaction in enumerate(reaction_ids):
                n = int(accumulator.truth_counts[eval_index, reaction_index])
                if n == 0:
                    continue
                row = {
                    "anchor": anchor,
                    "evaluation_id": evaluation["evaluation_id"],
                    "reaction_id": reaction,
                    "algorithm": evaluation["algorithm"],
                    "rna_context_key": evaluation["rna_context_key"],
                    "ct2a_mouse": evaluation["ct2a_mouse"],
                    "gl261_mouse": evaluation["gl261_mouse"],
                    "n_non_tie_pairs": n,
                }
                for field in accumulator.UNIT_FIELDS:
                    count = int(accumulator.counts[field][eval_index, reaction_index])
                    row[field] = accumulator.sums[field][eval_index, reaction_index] / count if count else math.nan
                writer.writerow(row)


def _write_tsv(path, records):
    if not records:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in records:
            writer.writerow({key: ("" if isinstance(value, float) and not math.isfinite(value) else value) for key, value in row.items()})


def _write_json_atomic(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _plot_results(out, case_paths, reaction_rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    r = pd.DataFrame(reaction_rows)
    case_frames = []
    rng = np.random.default_rng(20261011)
    columns = [
        "anchor", "relative_variance_reduction_correct_lambda_0.25",
        "relative_variance_reduction_incorrect_lambda_0.25", "variance_change_status_correct_lambda_0.25",
        "variance_change_status_incorrect_lambda_0.25", "prediction_gain_correct", "prediction_gain_incorrect",
    ]
    for anchor, parts in case_paths.items():
        for path in parts:
            sample = pd.read_csv(path, sep="\t", compression="xz", usecols=columns)
            case_frames.append(sample.sample(n=min(100, len(sample)), random_state=int(rng.integers(0, 2**31 - 1))))
    cases = pd.concat(case_frames, ignore_index=True)
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    ax = axes[0, 0]
    illustration = pd.read_csv(out / "ILLUSTRATIVE_DISTRIBUTION.tsv", sep="\t")
    values = illustration.magnitude.to_numpy(float)
    bins = np.linspace(float(values.min()), float(values.max()), 21)
    colors = {"baseline_weight": "#555555", "correct_lambda_0_25_weight": "#2171b5", "incorrect_lambda_0_25_weight": "#cb181d"}
    labels = {"baseline_weight": "Anchor baseline", "correct_lambda_0_25_weight": "Correct cue", "incorrect_lambda_0_25_weight": "Incorrect cue"}
    for column in colors:
        weights = illustration[column].to_numpy(float)
        ax.hist(values, bins=bins, weights=weights, histtype="step", linewidth=1.8, color=colors[column], label=labels[column])
        mean = float(np.dot(values, weights) / weights.sum())
        sd = float(np.sqrt(np.dot(weights, (values - mean) ** 2) / weights.sum()))
        ax.axvline(mean, color=colors[column], linewidth=.8, alpha=.8)
        ax.plot([], [], color=colors[column], label=f"μ={mean:.3g}, SD={sd:.3g}")
    ax.axvline(float(illustration.truth_magnitude.iloc[0]), color="black", linestyle="--", linewidth=1, label="Held-out truth")
    ax.set_xlabel("Magnitude |Δv|")
    ax.set_ylabel("Weighted probability per bin")
    ax.legend(frameon=False, fontsize=7)
    ax.set_title("A. Same support, reweighted mixture")
    ax = axes[0, 1]
    groups, labels = [], []
    for anchor in ("A1", "A2"):
        for cue, column in (("correct", "relative_variance_reduction_correct_lambda_0.25"), ("incorrect", "relative_variance_reduction_incorrect_lambda_0.25")):
            values = cases.loc[(cases.anchor == anchor) & cases[column].notna(), column].to_numpy()
            groups.append(values)
            labels.append(f"{anchor}\n{cue}")
    ax.boxplot(groups, labels=labels, showfliers=False)
    ax.axhline(0, color="black", linewidth=.7)
    ax.set_ylabel("Relative variance reduction")
    ax.set_title("B. Correct and incorrect cue uncertainty changes")
    ax = axes[1, 0]
    colors = {"BOTH_IMPROVE": "#238b45", "OVERCONFIDENCE_RISK": "#cb181d",
              "ACCURACY_IMPROVES_DISPERSION_INCREASES": "#2171b5", "BOTH_WORSEN": "#636363",
              "NEGLIGIBLE_UNCERTAINTY": "#bdbdbd"}
    categories = []
    plot_points = []
    for _, row in cases.iterrows():
        for cue in ("correct", "incorrect"):
            x = row[f"relative_variance_reduction_{cue}_lambda_0.25"]
            y = row[f"prediction_gain_{cue}"]
            status = row[f"variance_change_status_{cue}_lambda_0.25"]
            if not (math.isfinite(float(x)) and math.isfinite(float(y))):
                continue
            gain = "IMPROVES" if y > 1e-12 else ("WORSENS" if y < -1e-12 else "TIED")
            if status == "NEGLIGIBLE" or gain == "TIED":
                category = "NEGLIGIBLE_UNCERTAINTY"
            elif status == "DECREASED" and gain == "IMPROVES":
                category = "BOTH_IMPROVE"
            elif status == "DECREASED":
                category = "OVERCONFIDENCE_RISK"
            elif gain == "IMPROVES":
                category = "ACCURACY_IMPROVES_DISPERSION_INCREASES"
            else:
                category = "BOTH_WORSEN"
            plot_points.append((x, y, cue))
            categories.append(category)
    plot_points = np.asarray(plot_points, dtype=object)
    for category, color in colors.items():
        mask = np.asarray(categories) == category
        for cue, marker in (("correct", "o"), ("incorrect", "x")):
            sub = mask & (plot_points[:, 2] == cue)
            if np.any(sub):
                ax.scatter(plot_points[sub, 0].astype(float), plot_points[sub, 1].astype(float),
                           s=8, alpha=.35, color=color, marker=marker, label=category.replace("_", " ") if cue == "correct" else None)
    ax.axvline(0, color="black", linewidth=.7)
    ax.axhline(1e-12, color="gray", linewidth=.6)
    ax.set_xlabel("Relative variance reduction (○ correct, × incorrect cue)")
    ax.set_ylabel("Existing prediction-error gain")
    ax.legend(frameon=False, fontsize=6, loc="best")
    ax.set_title("C. Four uncertainty–accuracy outcomes")
    ax = axes[1, 1]
    for anchor, color in (("A1", "#2171b5"), ("A2", "#cb181d")):
        subset = r.loc[r.anchor.eq(anchor)]
        ax.scatter(subset["sign_magnitude_eta2"], subset["relative_correct"], s=8, alpha=.4,
                   marker="o", color=color, label=f"{anchor} correct")
        ax.scatter(subset["sign_magnitude_eta2"], subset["relative_incorrect"], s=8, alpha=.35,
                   marker="x", color=color, label=f"{anchor} incorrect")
    ax.set_xlabel(r"Baseline $\eta^2_{dir}$")
    ax.set_ylabel("Relative variance reduction")
    ax.legend(frameon=False, fontsize=7)
    ax.set_title("D. η²_dir versus uncertainty reduction")
    fig.savefig(out / "FIGURE_UNCERTAINTY_REDUCTION.png", dpi=220)
    fig.savefig(out / "FIGURE_UNCERTAINTY_REDUCTION.svg")
    plt.close(fig)


def _descriptive_summary(reaction_rows):
    import pandas as pd
    frame = pd.DataFrame(reaction_rows)
    rows = []
    for keys, subset in frame.groupby(["anchor", "analysis_split"], dropna=False):
        anchor, split = keys
        for outcome in ("relative_correct", "relative_incorrect", "absolute_correct", "absolute_incorrect",
                        "baseline_sd", "baseline_iqr", "baseline_central90_width", "sd_correct", "sd_incorrect",
                        "iqr_correct", "iqr_incorrect", "central90_correct", "central90_incorrect",
                        "gain_correct", "gain_incorrect", "crps_correct", "crps_incorrect"):
            values = pd.to_numeric(subset[outcome], errors="coerce").dropna()
            rows.append({"anchor": anchor, "analysis_split": split, "outcome": outcome,
                         "n_reactions": int(len(values)), "mean": float(values.mean()) if len(values) else math.nan,
                         "median": float(values.median()) if len(values) else math.nan,
                         "sd_across_reaction_means": float(values.std(ddof=1)) if len(values) > 1 else math.nan})
    return rows


def _explanatory_associations(reaction_rows):
    from scripts.dmi_bridge_uncertainty_reduction_core_v1 import explanatory_partial_rank_association
    import pandas as pd
    frame = pd.DataFrame(reaction_rows)
    records = []
    for anchor, subset in frame.groupby("anchor"):
        for cue, outcome, gain in (
            ("correct", "relative_correct", "gain_correct"),
            ("incorrect", "relative_incorrect", "gain_incorrect"),
        ):
            values = explanatory_partial_rank_association(
                pd.to_numeric(subset[outcome], errors="coerce"),
                pd.to_numeric(subset["usefulness"], errors="coerce"),
                pd.to_numeric(subset["sign_magnitude_eta2"], errors="coerce"),
            )
            records.append({"anchor": anchor, "cue": cue, "adjustment": "observed_usefulness",
                            "outcome": "relative_variance_reduction", **values})
            # Keep the established prediction-usefulness association separate.
            values = explanatory_partial_rank_association(
                pd.to_numeric(subset[gain], errors="coerce"),
                pd.to_numeric(subset["usefulness"], errors="coerce"),
                pd.to_numeric(subset["sign_magnitude_eta2"], errors="coerce"),
            )
            records.append({"anchor": anchor, "cue": cue, "adjustment": "observed_usefulness",
                            "outcome": "prediction_error_gain", **values})
    return records


def _mouse_cluster_bootstrap(summary_paths, *, replicates=5000, seed=20260929):
    """Cross-resample mouse identities, sharing each draw across A1 and A2."""
    import pandas as pd
    rng = np.random.default_rng(seed)
    columns = ("relative_correct", "relative_incorrect", "absolute_correct", "absolute_incorrect",
               "gain_correct", "gain_incorrect", "crps_correct", "crps_incorrect")
    pair_values = {}
    mouse_ids = {"CT2A": set(), "GL261": set()}
    for anchor, path in summary_paths.items():
        pair_values[anchor] = {}
        accum = defaultdict(lambda: {"sum": np.zeros(len(columns)), "n": np.zeros(len(columns), dtype=np.int64)})
        for chunk in pd.read_csv(path, sep="\t", compression="xz", chunksize=50000):
            for tumor, column in (("CT2A", "ct2a_mouse"), ("GL261", "gl261_mouse")):
                mouse_ids[tumor].update(chunk[column].dropna().astype(str).unique())
            chunk["_pair"] = list(zip(chunk.ct2a_mouse.astype(str), chunk.gl261_mouse.astype(str)))
            for pair, group in chunk.groupby("_pair", sort=False):
                vals = group[list(columns)].apply(pd.to_numeric, errors="coerce").to_numpy(float)
                finite = np.isfinite(vals)
                entry = accum[pair]
                entry["sum"] += np.where(finite, vals, 0.0).sum(axis=0)
                entry["n"] += finite.sum(axis=0)
        pair_values[anchor] = accum
    mice_left, mice_right = sorted(mouse_ids["CT2A"]), sorted(mouse_ids["GL261"])
    if len(mice_left) != 5 or len(mice_right) != 5:
        raise RuntimeError(f"expected five mice per cohort for crossed bootstrap; got {len(mice_left)}, {len(mice_right)}")
    boot = {anchor: {column: [] for column in columns} for anchor in summary_paths}
    paired_boot = {column: [] for column in columns} if {"A1", "A2"}.issubset(summary_paths) else None
    for _ in range(replicates):
        left_draw = rng.choice(mice_left, size=5, replace=True)
        right_draw = rng.choice(mice_right, size=5, replace=True)
        lc = Counter(left_draw.tolist()); rc = Counter(right_draw.tolist())
        weights = {pair: lc[pair[0]] * rc[pair[1]] for pair in pair_values["A1"]}
        draw_values = {}
        for anchor in summary_paths:
            draw_values[anchor] = {}
            for index, column in enumerate(columns):
                numerator = denominator = 0.0
                for pair, count in weights.items():
                    entry = pair_values[anchor].get(pair)
                    if count and entry and entry["n"][index]:
                        weight = count * int(entry["n"][index])
                        numerator += count * float(entry["sum"][index])
                        denominator += weight
                value = numerator / denominator if denominator else math.nan
                boot[anchor][column].append(value)
                draw_values[anchor][column] = value
        if paired_boot is not None:
            for column in columns:
                left, right = draw_values["A1"][column], draw_values["A2"][column]
                paired_boot[column].append(right - left if math.isfinite(left) and math.isfinite(right) else math.nan)
    records = []
    for anchor in summary_paths:
        for column in columns:
            values = np.asarray(boot[anchor][column], dtype=float)
            values = values[np.isfinite(values)]
            records.append({"anchor": anchor, "outcome": column, "replicates": replicates,
                            "seed": seed, "mean": float(values.mean()) if len(values) else math.nan,
                            "ci95_low": float(np.quantile(values, .025)) if len(values) else math.nan,
                            "ci95_high": float(np.quantile(values, .975)) if len(values) else math.nan,
                            "finite_replicates": int(len(values)), "mouse_ids_ct2a": len(mice_left),
                            "mouse_ids_gl261": len(mice_right), "method": "crossed_mouse_bootstrap",
                            "mouse_cohort": "", "omitted_mouse": ""})
            for cohort, mice, index in (("CT2A", mice_left, 0), ("GL261", mice_right, 1)):
                for omitted in mice:
                    numerator = denominator = 0.0
                    for pair, entry in pair_values[anchor].items():
                        if pair[index] == omitted or not entry["n"][columns.index(column)]:
                            continue
                        n = int(entry["n"][columns.index(column)])
                        numerator += float(entry["sum"][columns.index(column)])
                        denominator += n
                    records.append({"anchor": anchor, "outcome": column, "replicates": 0,
                                    "seed": seed, "mean": numerator / denominator if denominator else math.nan,
                                    "ci95_low": math.nan, "ci95_high": math.nan, "finite_replicates": 0,
                                    "mouse_ids_ct2a": len(mice_left), "mouse_ids_gl261": len(mice_right),
                                    "method": "leave_one_mouse_out", "mouse_cohort": cohort,
                                    "omitted_mouse": omitted})
    if paired_boot is not None:
        for column, draws in paired_boot.items():
            values = np.asarray(draws, dtype=float)
            values = values[np.isfinite(values)]
            records.append({"anchor": "A2_minus_A1", "outcome": column, "replicates": replicates,
                            "seed": seed, "mean": float(values.mean()) if len(values) else math.nan,
                            "ci95_low": float(np.quantile(values, .025)) if len(values) else math.nan,
                            "ci95_high": float(np.quantile(values, .975)) if len(values) else math.nan,
                            "finite_replicates": int(len(values)), "mouse_ids_ct2a": len(mice_left),
                            "mouse_ids_gl261": len(mice_right), "method": "crossed_mouse_bootstrap_paired_difference",
                            "mouse_cohort": "", "omitted_mouse": ""})
    return records


def benchmark(cases: int = 64, repeats: int = 2) -> dict[str, object]:
    """Synthetic performance benchmark; values are never empirical findings."""
    import resource
    import dmi_bridge_uncertainty_reduction_core_v1 as core

    rng = np.random.default_rng(20261011)
    elapsed = []
    for _ in range(repeats):
        start = time.perf_counter()
        for _case in range(cases):
            delta = rng.normal(size=361)
            weights = rng.random(361)
            weights /= weights.sum()
            core.analyze_distribution(delta, weights, float(rng.normal()))
        elapsed.append(time.perf_counter() - start)
    best = min(elapsed)
    return {
        "status": "SYNTHETIC_BENCHMARK_ONLY",
        "seed": 20261011,
        "synthetic_cases": cases,
        "repeats": repeats,
        "seconds": elapsed,
        "best_cases_per_second": cases / best if best > 0 else math.nan,
        "estimated_case_count_two_anchors": 2 * 836 * 4179,
        "extrapolated_seconds_at_best_rate": (2 * 836 * 4179) * best / cases,
        "peak_rss_kib": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
        "interpretation": "synthetic throughput projection only; no frozen biological inputs were accessed",
    }


def run_empirical(output_dir: Path = OUT):
    """Run the full empirical analysis after every frozen input passes preflight."""
    import pandas as pd
    from scripts import dmi_bridge_pl1_predictability_core_v1 as pl1
    from scripts import dmi_bridge_pl2_sign_only_core_v1 as pl2
    from scripts import dmi_bridge_pl2a_prepare_v1 as pl2a
    from scripts import dmi_bridge_pl2b_production_v1 as pl2b
    from scripts import dmi_bridge_pl2c_geometry_utility_core_v1 as pl2c
    import dmi_bridge_uncertainty_reduction_core_v1 as core

    gate = preflight()
    if gate["status"] != "READY":
        raise RuntimeError("empirical analysis blocked: frozen inputs are missing or mismatched")
    output_dir.mkdir(parents=True, exist_ok=True)
    state_path = output_dir / "RUN_STATE.json"
    if state_path.exists():
        state = json.loads(state_path.read_text(encoding="utf-8"))
        if state.get("input_hashes") != gate["available_sha256"] or state.get("schema") != OUTPUT_SCHEMA:
            raise RuntimeError("existing output directory belongs to different inputs or schema")
    elif any(output_dir.iterdir()):
        raise RuntimeError(f"nonempty output directory has no resumable state: {output_dir}")
    else:
        state = {"schema": OUTPUT_SCHEMA, "input_hashes": gate["available_sha256"], "completed_parts": {}}
        _write_json_atomic(state_path, state)
    def source(name):
        relative, digest = EXPECTED[name]
        resolved, error = resolve_file(relative, digest)
        if resolved is None:
            raise RuntimeError(error or f"missing {relative}")
        return resolved

    pl1_manifest, pl1_artifacts = _manifest_bundle("PL1", EXPECTED["pl1_manifest"][0], EXPECTED["pl1_manifest"][1])
    pl2a_manifest, pl2a_artifacts = _manifest_bundle("PL2A", EXPECTED["pl2a_manifest"][0], EXPECTED["pl2a_manifest"][1])
    pl2b_manifest, pl2b_artifacts = _manifest_bundle("PL2B", EXPECTED["pl2b_manifest"][0], EXPECTED["pl2b_manifest"][1])
    a20_manifest, a20_artifacts = _manifest_bundle("A20", EXPECTED["a20_manifest"][0], EXPECTED["a20_manifest"][1])
    a21_manifest, a21_artifacts = _manifest_bundle("A21", "outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_MANIFEST.json", "db9948d7b99447a695cf8aa21429dde6d642e1fd0632bcce52643deac59e6feb")
    a22_manifest, a22_artifacts = _manifest_bundle("A22", "outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_MANIFEST.json", "15b237151005538e5077a352868f873ef00466f4f5ea49307b781a8a8cd0e7cb")
    a20_weights = read_tsv(source("a20_weights"))
    candidates = read_tsv(source("candidate_table"))
    bio0_weights = read_tsv(source("a1_weights"))
    pairs_alias = read_tsv(pl2a_artifacts["BRIDGEPL2A_TRUTH_PAIR_REGISTRY.tsv"])
    pairs, _aliases = pl2b.collapse_truth_pairs(pairs_alias)
    candidate_by_id = {pl2a.candidate_identity(row): row for row in candidates}
    with np.load(source("flux_cache"), allow_pickle=False) as cache:
        mats = cache["mats"]
        rxns = [str(value) for value in cache["rxns"].tolist()]
        algs, tumors, ensembles = [str(x) for x in cache["alg"].tolist()], [str(x) for x in cache["tumor"].tolist()], [str(x) for x in cache["eh"].tolist()]
    cache_info = {"rxns": rxns, "alg": algs, "tumor": tumors, "eh": ensembles}
    cache_index = {(a, t, e): i for i, (a, t, e) in enumerate(zip(algs, tumors, ensembles))}
    reaction_rows = read_tsv(pl2a_artifacts["BRIDGEPL2A_REACTION_REGISTRY.tsv"])
    reaction_ids = sorted(row["reaction_id"] for row in reaction_rows if row["reaction_id"] not in {"HEX1", "LDH_L"})
    if len(reaction_ids) != 4179 or len(set(reaction_ids)) != 4179:
        raise RuntimeError("common reaction population mismatch")
    rxn_index = {reaction: i for i, reaction in enumerate(rxns)}
    # PL2A deliberately consumes and preserves the pre-truth PL1 evaluation registry.
    eval_rows = read_tsv(pl1_artifacts["BRIDGEPL1_EVALUATION_REGISTRY.tsv"])
    eval_ids = [row["evaluation_id"] for row in eval_rows]
    split_path, split_error = resolve_file("outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_REACTION_SPLIT.tsv", "7dd8676508416a13f88af184110ef40e0f8d92782cd5d62b42e7006af5316e2f")
    if split_path is None:
        raise RuntimeError(split_error or "missing frozen reaction split")
    split_map = {row["reaction_id"]: row["analysis_split"] for row in read_tsv(split_path)}
    if len(split_map) != 4180:
        raise RuntimeError("frozen split reaction count mismatch")
    # PL2A's A1 holdout uses the frozen BIO-0 strong-anchor baseline arm.
    a1_pools = defaultdict(list)
    for row in bio0_weights:
        if row["arm"] == "strong_anchor_baseline":
            key = (row["mouse_id"], row["tumor"], row["algorithm"], row["ensemble_hash"], row["rna_context_key"])
            a1_pools[key].append((int(row["sample_index"]), float(row["weight"])))
    a1_specs, a1_qc = pl2b.build_holdout_specs(pairs, candidates, {"weights": a1_pools}, mats, cache_info)
    a2_specs = _support_specs_a2(pairs, candidate_by_id, a20_weights)
    common_a1 = {str(spec["pair"]["truth_pair_id"]) for spec in a1_specs}
    common_a2 = {str(spec["pair"]["truth_pair_id"]) for spec in a2_specs}
    core.exact_key_match([(key,) for key in common_a1], [(key,) for key in common_a2])
    a1_case_parts = {"manifest_path": pl2b_artifacts["BRIDGEPL2B_CASE_OUTCOMES.parts.json"], "artifacts": pl2b_artifacts}
    a2_case_parts = {"manifest_path": a22_artifacts["BRIDGEA22_CASE_OUTCOMES.parts.json"], "artifacts": a22_artifacts}
    case_paths = {}
    accumulators = {}
    qc = {"status": "PASS", "a1_holdout_qc": a1_qc, "matched_truth_pairs": len(common_a1), "common_reactions": len(reaction_ids), "per_anchor": {}}
    for anchor, specs, cached_parts in (("A1", a1_specs, a1_case_parts), ("A2", a2_specs, a2_case_parts)):
        accumulator, records, anchor_qc = _analyze_anchor(
            anchor, specs, reaction_ids, rxn_index, mats, cache_index,
            cached_parts, output_dir, eval_ids,
        )
        accumulators[anchor] = accumulator
        case_paths[anchor] = [record["path"] for record in records]
        qc["per_anchor"][anchor] = {**anchor_qc, "case_parts": [{"name": x["path"].name, "rows": x["rows"], "sha256": x["sha256"]} for x in records]}
        unit_path = output_dir / f"REACTION_EVALUATION_SUMMARY_{anchor}.tsv.xz"
        _save_unit_table(anchor, accumulator, eval_rows, reaction_ids, unit_path)
        qc["per_anchor"][anchor]["reaction_evaluation_rows"] = int(np.count_nonzero(accumulator.truth_counts))
    feature_fields = ("sign_magnitude_eta2", "directional_entropy3", "dominant_sign_mass", "n_supported_sign_states", "p_negative", "p_positive")
    a1_feature_paths = _feature_parts(pl1_manifest, pl1_artifacts, "BRIDGEPL1")
    a2_feature_paths = _feature_parts(a21_manifest, a21_artifacts, "BRIDGEA21")
    a1_features = _feature_means(a1_feature_paths, reaction_ids, feature_fields)
    a2_features_raw = _feature_means(a2_feature_paths, reaction_ids, feature_fields)
    for table in (a1_features, a2_features_raw):
        for item in table.values():
            item["non_tie_coverage"] = item.pop("p_negative") + item.pop("p_positive") if math.isfinite(item.get("p_negative", math.nan)) and math.isfinite(item.get("p_positive", math.nan)) else math.nan
    reaction_rows = []
    for anchor, features in (("A1", a1_features), ("A2", a2_features_raw)):
        reaction_rows.extend(accumulators[anchor].reaction_rows(anchor, features, split_map))
    # Prospective predictions require a shared finite population across A1/A2
    # and both cue targets, and use only outcome-blind baseline descriptors.
    frames = {anchor: {row["reaction_id"]: row for row in reaction_rows if row["anchor"] == anchor} for anchor in ("A1", "A2")}
    common_reactions = []
    for reaction in reaction_ids:
        candidates_rows = [frames[anchor][reaction] for anchor in ("A1", "A2")]
        required = ["relative_correct", "relative_incorrect", "sign_magnitude_eta2", "directional_entropy3", "dominant_sign_mass", "non_tie_coverage", "n_supported_sign_states"]
        if all(math.isfinite(float(row.get(field, math.nan))) for row in candidates_rows for field in required):
            if candidates_rows[0]["analysis_split"] == candidates_rows[1]["analysis_split"]:
                common_reactions.append(reaction)
    models = []
    for anchor in ("A1", "A2"):
        data_rows = []
        for reaction in common_reactions:
            row = frames[anchor][reaction]
            data_rows.append({
                "reaction_id": reaction,
                "analysis_split": row["analysis_split"],
                "sign_magnitude_eta2": row["sign_magnitude_eta2"],
                "directional_entropy3": row["directional_entropy3"],
                "dominant_sign_mass": row["dominant_sign_mass"],
                "non_tie_coverage": row["non_tie_coverage"],
                "n_supported_sign_states": row["n_supported_sign_states"],
                "correct_relative_variance_reduction": row["relative_correct"],
                "incorrect_relative_variance_reduction": row["relative_incorrect"],
            })
        for target in ("correct_relative_variance_reduction", "incorrect_relative_variance_reduction"):
            for result in core.fit_prospective_models(data_rows, target=target):
                result["anchor"] = anchor
                models.append(result)
    _write_tsv(output_dir / "PROSPECTIVE_GEOMETRY_PREDICTION.tsv", models)
    _write_tsv(output_dir / "REACTION_LEVEL_PLOTTING_TABLE.tsv", reaction_rows)
    _write_tsv(output_dir / "CASE_LEVEL_OUTCOME_COUNTS.tsv", [
        {"anchor": anchor, "aggregation_level": "matched_truth_case_reaction",
         "outcome": "|".join(map(str, key)), "count": int(value)}
        for anchor, acc in accumulators.items() for key, value in acc.case_summaries.items()
    ])
    _write_tsv(output_dir / "SAMPLE_COUNTS.tsv", [
        {"anchor": anchor, "case_rows": int(acc.case_summaries[("case_rows",)]),
         "truth_tie_rows": int(acc.case_summaries[("truth_ties",)]),
         "non_tie_case_rows": int(acc.case_summaries[("non_ties",)]),
         "near_zero_baseline_variance_rows": int(acc.case_summaries[("relative_variance_excluded_near_zero",)]),
         "relative_variance_eligible_rows": int(acc.case_summaries[("non_ties",)] - acc.case_summaries[("relative_variance_excluded_near_zero",)]),
         "nonfinite_inputs_excluded": 0, "nonfinite_policy": "fail-fast; no case is silently dropped"}
        for anchor, acc in accumulators.items()
    ])
    _write_tsv(output_dir / "CASE_STRATUM_OUTCOME_COUNTS.tsv", [
        {"anchor": anchor, "algorithm": key[0], "rna_context_key": key[1], "cue": key[2],
         "lambda": key[3], "uncertainty_status": key[4], "prediction_status": key[5], "count": int(value)}
        for anchor, acc in accumulators.items() for key, value in acc.stratum_summaries.items()
    ])
    _write_tsv(output_dir / "SENSITIVITY_SUMMARY.tsv", [
        {"anchor": anchor, "cue": key[0], "lambda": key[1], "endpoint": "variance_change",
         "status": key[2], "count": int(value)}
        for anchor, acc in accumulators.items() for key, value in acc.case_summaries.items()
        if len(key) == 3 and str(key[1]).startswith("lambda_")
    ] + [
        {"anchor": anchor, "cue": key[0], "lambda": key[1], "endpoint": key[2],
         "status": key[3], "count": int(value)}
        for anchor, acc in accumulators.items() for key, value in acc.case_summaries.items()
        if len(key) == 4 and str(key[1]).startswith("lambda_")
    ])
    _write_tsv(output_dir / "DESCRIPTIVE_SUMMARY.tsv", _descriptive_summary(reaction_rows))
    _write_tsv(output_dir / "EXPLANATORY_USEFULNESS_ASSOCIATION.tsv", _explanatory_associations(reaction_rows))
    bootstrap_rows = _mouse_cluster_bootstrap(
        {anchor: output_dir / f"REACTION_EVALUATION_SUMMARY_{anchor}.tsv.xz" for anchor in ("A1", "A2")}
    )
    _write_tsv(output_dir / "MOUSE_CLUSTER_BOOTSTRAP.tsv", bootstrap_rows)
    _plot_results(output_dir, case_paths, reaction_rows)
    empirical_summary = _empirical_report_summary(accumulators, reaction_rows, models, common_reactions)
    report = {
        "schema": OUTPUT_SCHEMA,
        "status": "COMPLETE",
        "case_rows_per_anchor": 836 * 4179,
        "matched_truth_pairs": 836,
        "reaction_population": 4179,
        "relative_variance_threshold": "V0 > 100*eps*max(E[M^2], E[M]^2)",
        "reaction_aggregation": "equal-weight truth pairs within evaluation/reaction, then equal-weight eligible evaluation/reaction units",
        "prospective_targets": ["correct_relative_variance_reduction", "incorrect_relative_variance_reduction"],
        "prospective_predictors": list(core.prospective_model_specs()),
        "prospective_finite_common_reactions": len(common_reactions),
        "prospective_split_counts": dict(Counter(split_map[r] for r in common_reactions)),
        "prospective_population": "identical complete-case reaction set and frozen Development/Confirmation labels across A1/A2 and both separate cue targets",
        "accuracy_definition": "existing PL2 correct/wrong expected-magnitude absolute-error gain",
        "crps_definition": "weighted empirical CRPS; positive improvement is baseline minus update",
        "analysis_limit": "synthetic-truth ensemble analysis; not experimental calibration",
        "bootstrap": "5000 independently crossed CT2A and GL261 mouse resamples (5 mice each), shared draws across A1/A2; reaction/evaluation units remain clustered within sampled mouse pairs",
        **empirical_summary,
    }
    (output_dir / "SCIENTIFIC_REPORT.md").write_text(_render_report(report, models), encoding="utf-8")
    manifest = {
        "schema": OUTPUT_SCHEMA,
        "status": "UNCERTAINTY_REDUCTION_ANALYSIS_COMPLETE",
        "input_hashes": gate["available_sha256"],
        "matched_truth_pairs": 836,
        "common_reactions": 4179,
        "qc": qc,
        "report": report,
        "outputs": {p.name: sha256_file(p) for p in output_dir.iterdir() if p.is_file()},
    }
    (output_dir / "MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def _empirical_report_summary(accumulators, reaction_rows, models, common_reactions):
    import pandas as pd
    frequency = {}
    groups = {}
    for anchor, accumulator in accumulators.items():
        total = int(accumulator.case_summaries[("non_ties",)])
        excluded = int(accumulator.case_summaries[("relative_variance_excluded_near_zero",)])
        frequency[anchor] = {}
        groups[anchor] = {}
        for cue in ("correct", "incorrect"):
            counts = {status.lower(): int(accumulator.case_summaries[(cue, status)]) for status in ("DECREASED", "INCREASED", "NEGLIGIBLE")}
            frequency[anchor][cue] = {
                "decreased": counts["decreased"], "increased": counts["increased"],
                "negligible": counts["negligible"], "non_tie": total,
                "decreased_fraction": counts["decreased"] / total if total else math.nan,
                "relative_excluded": excluded,
            }
            both_improve = overconfidence = accuracy_only = both_worsen = negligible = 0
            for u_status in ("DECREASED", "INCREASED", "NEGLIGIBLE"):
                for p_status in ("IMPROVED", "WORSENED", "TIED"):
                    n = int(accumulator.case_summaries[(cue, u_status, p_status)])
                    if u_status == "NEGLIGIBLE" or p_status == "TIED":
                        negligible += n
                    elif u_status == "DECREASED" and p_status == "IMPROVED":
                        both_improve += n
                    elif u_status == "DECREASED":
                        overconfidence += n
                    elif p_status == "IMPROVED":
                        accuracy_only += n
                    else:
                        both_worsen += n
            groups[anchor][cue] = {
                "both_improve": both_improve, "overconfidence_risk": overconfidence,
                "accuracy_improves_dispersion_increases": accuracy_only,
                "both_worsen": both_worsen, "negligible_or_tied": negligible,
            }
    reaction_frame = pd.DataFrame(reaction_rows)
    reaction_means = {}
    for anchor, subset in reaction_frame.groupby("anchor"):
        reaction_means[anchor] = {}
        for cue, column in (("correct", "relative_correct"), ("incorrect", "relative_incorrect")):
            values = pd.to_numeric(subset[column], errors="coerce").dropna()
            reaction_means[anchor][cue] = {
                "finite_reactions": int(len(values)),
                "mean_equal_reaction": float(values.mean()) if len(values) else math.nan,
                "median_equal_reaction": float(values.median()) if len(values) else math.nan,
            }
    model_index = {(row.get("anchor"), row.get("target"), row.get("model")): row for row in models}
    return {
        "uncertainty_frequencies": frequency,
        "accuracy_uncertainty_groups": groups,
        "reaction_level_relative_reduction": reaction_means,
        "model_index": {
            f"{anchor}|{target}|{name}": model_index.get((anchor, target, name), {})
            for anchor in ("A1", "A2")
            for target in ("correct_relative_variance_reduction", "incorrect_relative_variance_reduction")
            for name in ("development_mean_baseline", "eta2", "all_geometry")
        },
    }


def _title_verdict(summary):
    frequency = summary["uncertainty_frequencies"]
    model_index = summary["model_index"]
    broad_reduction = all(
        frequency[anchor]["correct"]["decreased_fraction"] > 0.5
        for anchor in ("A1", "A2")
    )
    eta2_predictive = all(
        model_index[f"{anchor}|correct_relative_variance_reduction|eta2"].get("mae", math.inf)
        < model_index[f"{anchor}|correct_relative_variance_reduction|development_mean_baseline"].get("mae", -math.inf)
        for anchor in ("A1", "A2")
    )
    any_geometry_predictive = all(
        model_index[f"{anchor}|correct_relative_variance_reduction|all_geometry"].get("mae", math.inf)
        < model_index[f"{anchor}|correct_relative_variance_reduction|development_mean_baseline"].get("mae", -math.inf)
        for anchor in ("A1", "A2")
    )
    if broad_reduction and eta2_predictive:
        return "A. Strongly supported", "Correct cues reduce variance in a majority of eligible cases under both anchors, and η²_dir alone predicts the correct-cue reduction better than the Development-mean baseline on Confirmation for both anchors. Accuracy gains and uncertainty changes remain separate outcomes."
    observed_reduction = any(frequency[anchor]["correct"]["decreased"] > 0 for anchor in ("A1", "A2"))
    if observed_reduction and any_geometry_predictive:
        return "B. Partially supported", "Variance reduction occurs in identifiable cases and the prespecified baseline geometry set predicts correct-cue reduction better than the Development-mean baseline on Confirmation under both anchors; reduction is not sufficiently consistent for the strongest wording."
    return "C. Not supported", "The empirical criteria for common correct-cue variance reduction and replicated Confirmation prediction were not both met; prediction-error improvement remains a separate possible conclusion."


def _render_report(summary, models):
    title_verdict, verdict_reason = _title_verdict(summary)
    frequency_lines = []
    for anchor, cues in summary["uncertainty_frequencies"].items():
        for cue, count in cues.items():
            frequency_lines.append(
                f"- {anchor} {cue}: decreased {count['decreased']}/{count['non_tie']} "
                f"({count['decreased_fraction']:.1%}), increased {count['increased']}, "
                f"negligible {count['negligible']}; near-zero relative exclusions {count['relative_excluded']}."
            )
    outcome_lines = [
        f"- {anchor} {cue}: both improve={counts['both_improve']}; overconfidence risk={counts['overconfidence_risk']}; "
        f"accuracy improves as dispersion increases={counts['accuracy_improves_dispersion_increases']}; "
        f"both worsen={counts['both_worsen']}; negligible/tied={counts['negligible_or_tied']}."
        for anchor, cues in summary["accuracy_uncertainty_groups"].items()
        for cue, counts in cues.items()
    ]
    model_lines = [
        f"- {row['anchor']} {row['target']} {row['model']}: Confirmation MAE={row.get('mae', math.nan):.5g}, "
        f"RMSE={row.get('rmse', math.nan):.5g}, R²={row.get('r2', math.nan):.4g}, "
        f"MAE improvement vs η²-only={row.get('mae_improvement_vs_eta2', math.nan):.5g}."
        for row in models if row.get("model") in {"development_mean_baseline", "eta2", "eta2_plus_entropy", "all_geometry"}
    ]
    lines = [
        "# Weak Directional Cues and Magnitude Uncertainty", "",
        "## Methods", "",
        "For each held-out truth pair and reaction, M=|Δv| is evaluated over the 19×19 held-out anchor-weighted sample-pair support after excluding both truth samples and renormalizing exactly as in the frozen PL2 holdout. The tie threshold is 1×10⁻¹². A1 reuses the frozen BIO-0 strong-anchor weights; A2 reuses the frozen A20 dual-anchor weights. No flux vectors or candidate identities are regenerated.", "",
        "The primary uncertainty statistic is V=Σwᵢ(Mᵢ−μ)². Absolute reduction is V₀−Vλ; relative reduction is 1−Vλ/V₀ and is omitted when V₀≤100ε·max(E[M²],E[M]²). The threshold and excluded counts are recorded per case. Secondary statistics are weighted SD, inverse-CDF weighted IQR, central 90% inverse-CDF width, and discrete weighted CRPS. Weighted quantiles use the first sorted support point whose cumulative normalized mass reaches q, matching PL1.", "",
        "For sign state s∈{−1,0,+1}, let πₛ be baseline mass and let μₛ and σ²ₛ be conditional magnitude mean and variance. A cue c changes only mixture masses: π'ₛ=cₛπₛ/Σₜcₜπₜ, where cₛ=exp(λ d s) and d is cue direction. Thus V'=Σₛπ'ₛσ²ₛ+Σₛπ'ₛ(μₛ−μ')². Conditional magnitude distributions do not change. The support remains unchanged whenever all baseline weights are positive; variance and quantile widths can still shrink or expand as direction-specific mixture mass shifts.", "",
        "Accuracy is the unchanged PL2 expected-absolute-magnitude prediction gain. CRPS improvement is baseline CRPS minus cue CRPS. Neither a mean shift nor a directional entropy change is treated as uncertainty reduction; all three endpoints are analyzed separately.", "",
        "CRPS uses CRPS(F,y)=E|X−y|−½E|X−X′| for the weighted discrete predictive distribution. Conditional sign-group expectations and cross-group absolute-distance terms are computed once per truth case and reused across both cue directions and all λ values; cue CRPS then uses only a three-state mixture calculation. The synthetic benchmark precedes any empirical run.", "",
        "## Prospective geometry prediction", "",
        "Correct-cue and incorrect-cue relative variance reduction are separate prediction targets; both use the identical complete-case reaction population and frozen Development/Confirmation split for A1 and A2. Prespecified ridge models compare η²_dir alone with η²_dir+H_dir, η²_dir+dominant-direction mass, η²_dir+non-tie coverage, η²_dir+supported sign states, and the full baseline geometry set. Scaling is fit on Development only. Confirmation usefulness and prediction gains never enter as prospective predictors.", "",
        "## Explanatory association", "",
        "Observed usefulness associations are descriptive and are kept separate from prospective model fitting. EXPLANATORY_USEFULNESS_ASSOCIATION.tsv reports rank association of η²_dir with uncertainty reduction before and after adjusting for observed directionally-useful outcome. It does not establish prospective prediction or experimental calibration.", "",
        "## Results", "",
        f"Matched distinct truth pairs: {summary['matched_truth_pairs']}; common reactions: {summary['reaction_population']}; case rows per anchor: {summary['case_rows_per_anchor']}. The complete-case prospective population contains {summary['prospective_finite_common_reactions']} reactions ({summary['prospective_split_counts']}).",
        "", "### Variance change frequencies", "", *frequency_lines, "",
        "### Prediction accuracy versus uncertainty", "", *outcome_lines, "",
        "CRPS improvement is a separate proper-score endpoint and does not define variance reduction or frozen prediction-error gain.", "",
        "### Prospective Confirmation scores", "", *model_lines, "",
        "Case frequencies stratified by reconstruction method and RNA context are in CASE_STRATUM_OUTCOME_COUNTS.tsv. Descriptor associations, continuous summaries, relative-variance exclusions, non-finite policy, and mouse-cluster intervals are in the corresponding tables.", "",
        "## Title verdict", "",
        f"**{title_verdict}.** {verdict_reason}", "",
        "## Limitations", "",
        "These results describe finite weighted synthetic-truth ensembles, not experimentally validated uncertainty estimates. The crossed bootstrap resamples only five mouse identities per cohort and preserves repeated evaluations/reactions within sampled mouse pairs; it does not capture uncertainty from new reconstruction methods, new cohorts, or the finite candidate-generation process. Prospective models have a single frozen reaction holdout and the all-geometry combination is descriptive absent an independently replicated confirmation set.", "",
    ]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--preflight", action="store_true", help="verify required input paths and hashes only")
    group.add_argument("--benchmark", action="store_true", help="run synthetic performance benchmark only")
    group.add_argument("--run", action="store_true", help="run empirical analysis after the hash gate")
    args = parser.parse_args()
    if args.preflight:
        result = preflight()
    elif args.benchmark:
        result = benchmark()
    else:
        result = run_empirical()
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
