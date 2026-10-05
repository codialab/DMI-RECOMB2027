#!/usr/bin/env python3
"""Produce the compact DMI-BRIDGE-PL2B sign-only utility bundle."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import lzma
import math
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np

try:
    from scripts import dmi_bridge_pl2_sign_only_core_v1 as core
    from scripts import dmi_bridge_pl2a_prepare_v1 as pl2a
    from scripts import dmi_bridge_pl2b_batch_v1 as batch
except ModuleNotFoundError:  # direct script execution
    import dmi_bridge_pl2_sign_only_core_v1 as core
    import dmi_bridge_pl2a_prepare_v1 as pl2a
    import dmi_bridge_pl2b_batch_v1 as batch


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_REL = Path("outputs/dmi_bridge_pl2b_sign_only_utility_v1")
PL2A_REL = Path("outputs/dmi_bridge_pl2a_sign_only_prepare_v1")
PL2A_MANIFEST_SHA256 = "8ae3eb60eb23614e379e38e3aa3cb72ec16c6bdb69aa9c74cc94499e31f1c312"
FOUNDATION_PATCH_REL = Path("dmi_bridge_pl2b_compact_production_foundation_v1.patch")
FOUNDATION_PATCH_SHA256 = "bbfffd7e0a45848f15bad164609112c066cc72ba03320766ef76ccc179b8246b"
CONTRACT_REL = Path("docs/DMI_BRIDGE_PL2B_PRODUCTION_CONTRACT.md")
COMPLETE_STATUS = "PL2B_SIGN_ONLY_UTILITY_COMPLETE"
NOOP_STATUS = "NO_OP_EXISTING_IDENTICAL_PL2B"
TRUTH_PAIR_SCHEMA = "bridge.pl2b.truth_pair_identity.v1"
PART_SCHEMA = "bridge.pl2b.case_outcomes.partitioned_tsv.v1"
WORK_SCHEMA = "bridge.pl2b.production_work.v1"
PART_COUNT = 128
MAX_PART_BYTES = 64 * 1024 * 1024
REACTION_BLOCK_SIZE = 64
EXPECTED_CASE_ROWS = 3_494_480
EXPECTED_COUNTS = {
    "truth_selection_rows": 480,
    "truth_selection_zero_post_holdout_mass_rows": 18,
    "truth_pair_aliases_total": 1200,
    "truth_pair_aliases_evaluable": 1110,
    "truth_pair_aliases_non_evaluable": 90,
    "reaction_registry": 4180,
    "pl1_predictor_join_rows": 1_672_000,
}
LAMBDA_PREFIX = {0.25: "lambda_0_25", 0.5: "lambda_0_5", 1.0: "lambda_1_0"}

IDENTITY_FIELDS = [
    "truth_pair_id", "evaluation_id", "algorithm", "rna_context_key",
    "ct2a_mouse", "gl261_mouse", "alias_count", "quantile_aliases",
    "reaction_id", "subsystem", "same_subsystem_as_HEX1", "proximal_set_member",
]
BASELINE_FIELDS = [
    "truth_delta_b", "truth_magnitude", "truth_direction", "truth_status",
    "baseline_magnitude_estimate", "baseline_signed_estimate", "baseline_abs_error",
    "baseline_joint_ess", "baseline_negative_mass", "baseline_tie_mass",
    "baseline_positive_mass",
]
LAMBDA_SUFFIXES = [
    "correct_magnitude_estimate", "wrong_magnitude_estimate",
    "correct_absolute_error", "wrong_absolute_error",
    "correct_absolute_error_gain", "wrong_absolute_error_gain",
    "correct_gain_status", "wrong_gain_status",
    "correct_posterior_ess", "wrong_posterior_ess",
    "random_sign_expected_gain", "random_sign_gain_status",
    "break_even_reliability",
]
CASE_FIELDS = IDENTITY_FIELDS + BASELINE_FIELDS + [
    f"{prefix}_{suffix}"
    for prefix in LAMBDA_PREFIX.values()
    for suffix in LAMBDA_SUFFIXES
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def json_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True) + "\n").encode("ascii")


def atomic_json(path: Path, value: object) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_bytes(json_bytes(value))
    os.replace(tmp, path)


def stable_hash(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


def format_value(value: object) -> str:
    return pl2a.format_value(value)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field)) for field in fields})


def resolve_regular(root: Path, relative: Path) -> Path:
    path = root / relative
    resolved = path.resolve(strict=True)
    try:
        resolved.relative_to(root.resolve())
    except ValueError as exc:
        raise RuntimeError(f"source resolves outside repository: {relative}") from exc
    if path.is_symlink() or not resolved.is_file():
        raise RuntimeError(f"source is not a regular repository file: {relative}")
    return resolved


def _implementation_hashes(root: Path = ROOT) -> dict[str, str]:
    relative = [
        Path("scripts/dmi_bridge_pl2_sign_only_core_v1.py"),
        Path("scripts/dmi_bridge_pl2b_batch_v1.py"),
        Path("scripts/dmi_bridge_pl2b_production_v1.py"),
        CONTRACT_REL,
        FOUNDATION_PATCH_REL,
    ]
    return {str(path): sha256_file(resolve_regular(root, path)) for path in relative}


def verify_admission(root: Path = ROOT) -> dict[str, object]:
    patch = resolve_regular(root, FOUNDATION_PATCH_REL)
    if sha256_file(patch) != FOUNDATION_PATCH_SHA256:
        raise RuntimeError("PL2B foundation patch hash mismatch")
    manifest_path = resolve_regular(root, PL2A_REL / "BRIDGEPL2A_MANIFEST.json")
    if sha256_file(manifest_path) != PL2A_MANIFEST_SHA256:
        raise RuntimeError("PL2A manifest hash mismatch")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "PL2A_READY_FOR_PL2B_SIGN_ONLY_UTILITY":
        raise RuntimeError("PL2A status does not admit PL2B")
    if manifest.get("counts") != EXPECTED_COUNTS:
        raise RuntimeError("PL2A manifest counts differ from the frozen contract")
    artifact_hashes = manifest.get("artifact_sha256", {})
    for name, expected in artifact_hashes.items():
        path = resolve_regular(root, PL2A_REL / name)
        if sha256_file(path) != expected:
            raise RuntimeError(f"PL2A artifact hash mismatch: {name}")
    observed_sources = pl2a.verify_sources()
    live_manifest = json.loads((root / pl2a.PL1_MANIFEST_REL).read_text(encoding="utf-8"))
    pl1_artifacts = pl2a.verify_pl1_artifacts(live_manifest)
    return {
        "pl2a_manifest_sha256": PL2A_MANIFEST_SHA256,
        "pl2a_status": manifest["status"],
        "pl2a_counts": manifest["counts"],
        "pl2a_artifact_sha256": artifact_hashes,
        "upstream_source_sha256": observed_sources,
        "pl1_consumed_artifacts": pl1_artifacts,
        "foundation_patch_sha256": FOUNDATION_PATCH_SHA256,
    }


def truth_pair_id(evaluation_id: str, ct2a_candidate_id: str, gl261_candidate_id: str) -> str:
    return stable_hash({
        "schema": TRUTH_PAIR_SCHEMA,
        "evaluation_id": evaluation_id,
        "ct2a_candidate_id": ct2a_candidate_id,
        "gl261_candidate_id": gl261_candidate_id,
    })


def collapse_truth_pairs(rows: list[dict[str, str]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    if len(rows) != 1200:
        raise RuntimeError("PL2A alias registry does not contain 1200 rows")
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (row["evaluation_id"], row["ct2a_candidate_id"], row["gl261_candidate_id"])
        grouped[key].append(row)
    distinct: list[dict[str, object]] = []
    alias_map: list[dict[str, object]] = []
    invariant_fields = [
        "algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse",
        "ct2a_candidate_id", "gl261_candidate_id", "ct2a_sample_index",
        "gl261_sample_index", "ct2a_post_holdout_mass", "gl261_post_holdout_mass",
        "ct2a_post_holdout_ess", "gl261_post_holdout_ess", "ct2a_holdout_status",
        "gl261_holdout_status", "pair_evaluable", "pair_status", "non_evaluable_reason",
    ]
    for key, aliases in grouped.items():
        aliases = sorted(aliases, key=lambda row: float(row["quantile_alias"]))
        representative = aliases[0]
        for row in aliases[1:]:
            if any(row[field] != representative[field] for field in invariant_fields):
                raise RuntimeError("quantile aliases disagree within a canonical truth pair")
        pair_id = truth_pair_id(*key)
        alias_names = [row["quantile_alias"] for row in aliases]
        record: dict[str, object] = {
            "truth_pair_id": pair_id,
            "evaluation_id": representative["evaluation_id"],
            "algorithm": representative["algorithm"],
            "rna_context_key": representative["rna_context_key"],
            "ct2a_mouse": representative["ct2a_mouse"],
            "gl261_mouse": representative["gl261_mouse"],
            "ct2a_candidate_id": representative["ct2a_candidate_id"],
            "gl261_candidate_id": representative["gl261_candidate_id"],
            "ct2a_sample_index": int(representative["ct2a_sample_index"]),
            "gl261_sample_index": int(representative["gl261_sample_index"]),
            "ct2a_post_holdout_mass": float(representative["ct2a_post_holdout_mass"]),
            "gl261_post_holdout_mass": float(representative["gl261_post_holdout_mass"]),
            "ct2a_post_holdout_ess": float(representative["ct2a_post_holdout_ess"]),
            "gl261_post_holdout_ess": float(representative["gl261_post_holdout_ess"]),
            "ct2a_holdout_status": representative["ct2a_holdout_status"],
            "gl261_holdout_status": representative["gl261_holdout_status"],
            "pair_evaluable": representative["pair_evaluable"] == "true",
            "pair_status": representative["pair_status"],
            "non_evaluable_reason": representative["non_evaluable_reason"],
            "alias_count": len(aliases),
            "quantile_aliases": ",".join(alias_names),
        }
        distinct.append(record)
        for row in aliases:
            alias_map.append({
                "truth_pair_id": pair_id,
                "alias_count": len(aliases),
                "quantile_aliases": ",".join(alias_names),
                **row,
            })
    distinct.sort(key=lambda row: str(row["truth_pair_id"]))
    alias_map.sort(key=lambda row: (str(row["truth_pair_id"]), float(str(row["quantile_alias"]))))
    evaluable = [row for row in distinct if bool(row["pair_evaluable"])]
    non_evaluable = [row for row in distinct if not bool(row["pair_evaluable"])]
    if len(distinct) != 885 or len(evaluable) != 836 or len(non_evaluable) != 49:
        raise RuntimeError("distinct truth-pair accounting mismatch")
    if Counter(int(row["alias_count"]) for row in evaluable) != Counter({1: 642, 2: 114, 3: 80}):
        raise RuntimeError("evaluable truth-pair alias multiplicity mismatch")
    if Counter(int(row["alias_count"]) for row in non_evaluable) != Counter({1: 24, 2: 9, 3: 16}):
        raise RuntimeError("non-evaluable truth-pair alias multiplicity mismatch")
    if len(alias_map) != 1200:
        raise RuntimeError("alias map row count mismatch")
    return distinct, alias_map


def _close(value: float, expected: float, *, tol: float = 1e-12) -> float:
    deviation = abs(value - expected)
    if not math.isfinite(deviation) or deviation > tol:
        raise RuntimeError(f"frozen holdout statistic mismatch: {value} versus {expected}")
    return deviation


def build_holdout_specs(
    distinct: list[dict[str, object]],
    candidates: list[dict[str, str]],
    data: dict,
    mats: np.ndarray,
    cache_info: dict,
) -> tuple[list[dict[str, object]], dict[str, float]]:
    candidate_by_id: dict[str, dict[str, str]] = {}
    for row in candidates:
        cid = pl2a.candidate_identity(row)
        if cid in candidate_by_id:
            raise RuntimeError("candidate identity is not unique")
        candidate_by_id[cid] = row
    cache_index = {
        (algorithm, tumor, ensemble): index
        for index, (algorithm, tumor, ensemble) in enumerate(
            zip(cache_info["alg"], cache_info["tumor"], cache_info["eh"])
        )
    }
    candidate_by_stratum_sample = {
        (row["algorithm"], row["tumor"], row["ensemble_hash"], int(row["sample_index"])): row
        for row in candidates
    }
    hex_index = cache_info["rxns"].index(pl2a.STRONG_ANCHOR)
    specs: list[dict[str, object]] = []
    max_mass_deviation = 0.0
    max_ess_deviation = 0.0
    for pair in distinct:
        if not pair["pair_evaluable"]:
            continue
        condition_specs = {}
        for tumor, prefix in (("CT2A", "ct2a"), ("GL261", "gl261")):
            cid = str(pair[f"{prefix}_candidate_id"])
            candidate = candidate_by_id.get(cid)
            if candidate is None:
                raise RuntimeError("selected candidate identity is absent from the frozen table")
            if candidate["algorithm"] != pair["algorithm"] or candidate["tumor"] != tumor or candidate["rna_context_key"] != pair["rna_context_key"]:
                raise RuntimeError("selected candidate identity disagrees with truth-pair context")
            sample_index = int(pair[f"{prefix}_sample_index"])
            if int(candidate["sample_index"]) != sample_index:
                raise RuntimeError("selected candidate sample index mismatch")
            mouse = str(pair[f"{prefix}_mouse"])
            weight_key = (
                mouse, tumor, str(pair["algorithm"]), candidate["ensemble_hash"],
                str(pair["rna_context_key"]),
            )
            raw_rows = data["weights"].get(weight_key, [])
            if len(raw_rows) != 20 or len({sample for sample, _ in raw_rows}) != 20:
                raise RuntimeError("frozen condition weight pool is not exactly 20 candidates")
            cache_key = (str(pair["algorithm"]), tumor, candidate["ensemble_hash"])
            if cache_key not in cache_index:
                raise RuntimeError("selected candidate stratum is absent from the flux cache")
            cache_position = cache_index[cache_key]
            raw_by_sample = dict(raw_rows)
            ordered_candidates = [
                candidate_by_stratum_sample[(cache_key[0], cache_key[1], cache_key[2], sample)]
                for sample in raw_by_sample
            ]
            ordered_candidates.sort(key=lambda row: (
                float(mats[cache_position, int(row["sample_index"]), hex_index]),
                pl2a.candidate_identity(row),
            ))
            samples = np.asarray([int(row["sample_index"]) for row in ordered_candidates], dtype=int)
            weights = core.normalize_weights([raw_by_sample[int(row["sample_index"])] for row in ordered_candidates])
            locations = np.flatnonzero(samples == sample_index)
            if locations.size != 1:
                raise RuntimeError("selected candidate is not unique in its inference pool")
            selected_position = int(locations[0])
            post_mass = float(1.0 - weights[selected_position])
            if post_mass <= 0.0:
                raise RuntimeError("evaluable pair has zero post-holdout mass")
            held_samples = np.delete(samples, selected_position)
            held_weights = np.delete(weights, selected_position) / post_mass
            post_ess = float(1.0 / np.sum(held_weights * held_weights))
            max_mass_deviation = max(
                max_mass_deviation,
                _close(post_mass, float(pair[f"{prefix}_post_holdout_mass"])),
            )
            max_ess_deviation = max(
                max_ess_deviation,
                _close(post_ess, float(pair[f"{prefix}_post_holdout_ess"])),
            )
            condition_specs[prefix] = {
                "cache_index": cache_position,
                "truth_sample": sample_index,
                "held_samples": held_samples,
                "held_weights": held_weights,
            }
        specs.append({"pair": pair, **condition_specs})
    specs.sort(key=lambda item: str(item["pair"]["truth_pair_id"]))
    if len(specs) != 836:
        raise RuntimeError("evaluable holdout specification count mismatch")
    return specs, {
        "maximum_post_holdout_mass_deviation": max_mass_deviation,
        "maximum_post_holdout_ess_deviation": max_ess_deviation,
    }


def _part_name(index: int) -> str:
    return f"BRIDGEPL2B_CASE_OUTCOMES.part-{index:03d}.tsv.xz"


def _part_ranges(total_pairs: int) -> list[tuple[int, int]]:
    return [
        ((total_pairs * index) // PART_COUNT, (total_pairs * (index + 1)) // PART_COUNT)
        for index in range(PART_COUNT)
    ]


def _max_deviation(current: float, *values: tuple[float, float]) -> float:
    for left, right in values:
        if math.isfinite(left) and math.isfinite(right):
            current = max(current, abs(left - right))
    return current


def _scalar_parity(scalar: dict[str, object], block_result: dict[str, object], position: int) -> tuple[float, float, float]:
    maximum = 0.0
    break_even_maximum = 0.0
    break_even_gain_scale_maximum = 0.0
    baseline_keys = [
        "truth_magnitude", "baseline_magnitude_estimate", "baseline_signed_estimate",
        "baseline_abs_error", "baseline_joint_ess", "baseline_negative_mass",
        "baseline_tie_mass", "baseline_positive_mass",
    ]
    for key in baseline_keys:
        maximum = _max_deviation(maximum, (float(scalar[key]), float(block_result[key][position])))
    if scalar["status"] == "TRUTH_TIE":
        if int(block_result["truth_direction"][position]) != 0:
            raise RuntimeError("scalar/batch truth-tie mismatch")
        return maximum, break_even_maximum, break_even_gain_scale_maximum
    scalar_lambdas = {float(row["lambda_total"]): row for row in scalar["lambda_results"]}
    for lam, result in block_result["lambda_results"].items():
        row = scalar_lambdas[float(lam)]
        pairs = [
            (row["correct"]["magnitude_mean"], result["correct_magnitude_estimate"][position]),
            (row["wrong"]["magnitude_mean"], result["wrong_magnitude_estimate"][position]),
            (row["correct"]["signed_mean"], result["correct_signed_estimate"][position]),
            (row["wrong"]["signed_mean"], result["wrong_signed_estimate"][position]),
            (row["correct"]["joint_ess"], result["correct_posterior_ess"][position]),
            (row["wrong"]["joint_ess"], result["wrong_posterior_ess"][position]),
            (row["correct_abs_error"], result["correct_absolute_error"][position]),
            (row["wrong_abs_error"], result["wrong_absolute_error"][position]),
            (row["correct_abs_error_gain"], result["correct_absolute_error_gain"][position]),
            (row["wrong_abs_error_gain"], result["wrong_absolute_error_gain"][position]),
            (row["random_sign_expected_abs_error_gain"], result["random_sign_expected_gain"][position]),
        ]
        maximum = _max_deviation(maximum, *[(float(a), float(b)) for a, b in pairs])
        left = float(row["break_even_reliability"])
        right = float(result["break_even_reliability"][position])
        if math.isfinite(left) and math.isfinite(right):
            difference = abs(left - right)
            break_even_maximum = max(break_even_maximum, difference)
            scalar_slope = float(row["correct_abs_error_gain"]) - float(row["wrong_abs_error_gain"])
            batch_slope = float(result["correct_absolute_error_gain"][position]) - float(result["wrong_absolute_error_gain"][position])
            gain_scale_deviation = difference * max(abs(scalar_slope), abs(batch_slope))
            break_even_gain_scale_maximum = max(break_even_gain_scale_maximum, gain_scale_deviation)
        elif math.isfinite(left) != math.isfinite(right):
            raise RuntimeError("scalar/batch break-even existence mismatch")
    return maximum, break_even_maximum, break_even_gain_scale_maximum


def _case_row(pair: dict[str, object], reaction: dict[str, str], result: dict[str, object], position: int) -> dict[str, object]:
    direction = int(result["truth_direction"][position])
    row: dict[str, object] = {
        "truth_pair_id": pair["truth_pair_id"],
        "evaluation_id": pair["evaluation_id"],
        "algorithm": pair["algorithm"],
        "rna_context_key": pair["rna_context_key"],
        "ct2a_mouse": pair["ct2a_mouse"],
        "gl261_mouse": pair["gl261_mouse"],
        "alias_count": pair["alias_count"],
        "quantile_aliases": pair["quantile_aliases"],
        "reaction_id": reaction["reaction_id"],
        "subsystem": reaction["subsystem"],
        "same_subsystem_as_HEX1": reaction["same_subsystem_as_HEX1"],
        "proximal_set_member": reaction["proximal_set_member"],
        "truth_delta_b": result["truth_delta_b"][position],
        "truth_magnitude": result["truth_magnitude"][position],
        "truth_direction": "POSITIVE" if direction > 0 else "NEGATIVE" if direction < 0 else "TIE",
        "truth_status": "NON_TIE" if direction else "TRUTH_TIE",
        "baseline_magnitude_estimate": result["baseline_magnitude_estimate"][position],
        "baseline_signed_estimate": result["baseline_signed_estimate"][position],
        "baseline_abs_error": result["baseline_abs_error"][position],
        "baseline_joint_ess": result["baseline_joint_ess"][position],
        "baseline_negative_mass": result["baseline_negative_mass"][position],
        "baseline_tie_mass": result["baseline_tie_mass"][position],
        "baseline_positive_mass": result["baseline_positive_mass"][position],
    }
    for lam, prefix in LAMBDA_PREFIX.items():
        values = result["lambda_results"][lam]
        if direction == 0:
            for suffix in LAMBDA_SUFFIXES:
                row[f"{prefix}_{suffix}"] = "TRUTH_TIE" if suffix.endswith("gain_status") else float("nan")
            continue
        correct_gain = float(values["correct_absolute_error_gain"][position])
        wrong_gain = float(values["wrong_absolute_error_gain"][position])
        random_gain = float(values["random_sign_expected_gain"][position])
        row.update({
            f"{prefix}_correct_magnitude_estimate": values["correct_magnitude_estimate"][position],
            f"{prefix}_wrong_magnitude_estimate": values["wrong_magnitude_estimate"][position],
            f"{prefix}_correct_absolute_error": values["correct_absolute_error"][position],
            f"{prefix}_wrong_absolute_error": values["wrong_absolute_error"][position],
            f"{prefix}_correct_absolute_error_gain": correct_gain,
            f"{prefix}_wrong_absolute_error_gain": wrong_gain,
            f"{prefix}_correct_gain_status": core.classify_gain(correct_gain),
            f"{prefix}_wrong_gain_status": core.classify_gain(wrong_gain),
            f"{prefix}_correct_posterior_ess": values["correct_posterior_ess"][position],
            f"{prefix}_wrong_posterior_ess": values["wrong_posterior_ess"][position],
            f"{prefix}_random_sign_expected_gain": random_gain,
            f"{prefix}_random_sign_gain_status": core.classify_gain(random_gain),
            f"{prefix}_break_even_reliability": values["break_even_reliability"][position],
        })
    return row


def write_case_part(
    path: Path,
    specs: list[dict[str, object]],
    reactions: list[dict[str, str]],
    reaction_columns: np.ndarray,
    mats: np.ndarray,
) -> dict[str, object]:
    tmp = path.with_name(path.name + ".tmp")
    tmp.unlink(missing_ok=True)
    rows_written = 0
    lambda_zero_max = 0.0
    scalar_batch_max = 0.0
    scalar_batch_break_even_max = 0.0
    scalar_batch_break_even_gain_scale_max = 0.0
    with lzma.open(tmp, "wt", encoding="utf-8", newline="", preset=0, format=lzma.FORMAT_XZ) as handle:
        writer = csv.DictWriter(handle, fieldnames=CASE_FIELDS, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        for spec in specs:
            pair = spec["pair"]
            ct = spec["ct2a"]
            gl = spec["gl261"]
            for start in range(0, len(reactions), REACTION_BLOCK_SIZE):
                stop = min(start + REACTION_BLOCK_SIZE, len(reactions))
                cols = reaction_columns[start:stop]
                result = batch.evaluate_truth_block(
                    ct2a_values=mats[int(ct["cache_index"])][ct["held_samples"]][:, cols],
                    gl261_values=mats[int(gl["cache_index"])][gl["held_samples"]][:, cols],
                    ct2a_weights=ct["held_weights"],
                    gl261_weights=gl["held_weights"],
                    truth_ct2a=mats[int(ct["cache_index"]), int(ct["truth_sample"]), cols],
                    truth_gl261=mats[int(gl["cache_index"]), int(gl["truth_sample"]), cols],
                )
                non_tie = np.asarray(result["truth_direction"]) != 0
                zero = result["lambda_results"][0.0]
                if np.any(non_tie):
                    for key, baseline_key in (
                        ("correct_magnitude_estimate", "baseline_magnitude_estimate"),
                        ("wrong_magnitude_estimate", "baseline_magnitude_estimate"),
                        ("correct_signed_estimate", "baseline_signed_estimate"),
                        ("wrong_signed_estimate", "baseline_signed_estimate"),
                        ("correct_posterior_ess", "baseline_joint_ess"),
                        ("wrong_posterior_ess", "baseline_joint_ess"),
                    ):
                        deviation = np.max(np.abs(np.asarray(zero[key])[non_tie] - np.asarray(result[baseline_key])[non_tie]))
                        lambda_zero_max = max(lambda_zero_max, float(deviation))
                if start == 0:
                    delta, weights = core.joint_delta_distribution(
                        mats[int(ct["cache_index"])][ct["held_samples"], cols[0]],
                        mats[int(gl["cache_index"])][gl["held_samples"], cols[0]],
                        ct["held_weights"], gl["held_weights"],
                    )
                    scalar = core.evaluate_truth_case(
                        delta_b=delta,
                        baseline_weights=weights,
                        truth_delta_b=float(result["truth_delta_b"][0]),
                    )
                    parity, break_even_parity, break_even_gain_scale = _scalar_parity(scalar, result, 0)
                    scalar_batch_max = max(scalar_batch_max, parity)
                    scalar_batch_break_even_max = max(scalar_batch_break_even_max, break_even_parity)
                    scalar_batch_break_even_gain_scale_max = max(
                        scalar_batch_break_even_gain_scale_max, break_even_gain_scale
                    )
                for local, reaction in enumerate(reactions[start:stop]):
                    row = _case_row(pair, reaction, result, local)
                    writer.writerow({field: format_value(row[field]) for field in CASE_FIELDS})
                    rows_written += 1
    os.replace(tmp, path)
    if path.stat().st_size >= MAX_PART_BYTES:
        raise RuntimeError(f"PL2B case part exceeds the 64 MiB safety ceiling: {path.name}")
    if lambda_zero_max > 1e-12:
        raise RuntimeError("lambda-zero recovery exceeded tolerance")
    if scalar_batch_max > 5e-12:
        raise RuntimeError("scalar/batch sentinel parity exceeded tolerance")
    if scalar_batch_break_even_gain_scale_max > 5e-12:
        raise RuntimeError("scalar/batch break-even gain-scale parity exceeded tolerance")
    return {
        "filename": path.name,
        "rows": rows_written,
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "lambda_zero_maximum_deviation": lambda_zero_max,
        "scalar_batch_maximum_deviation": scalar_batch_max,
        "scalar_batch_break_even_maximum_deviation": scalar_batch_break_even_max,
        "scalar_batch_break_even_gain_scale_maximum_deviation": scalar_batch_break_even_gain_scale_max,
    }


def _validate_completed_part(work: Path, record: dict[str, object], expected_rows: int) -> None:
    path = work / str(record["filename"])
    if not path.is_file() or path.stat().st_size != int(record["bytes"]) or sha256_file(path) != record["sha256"]:
        raise RuntimeError("completed PL2B work part failed hash/size validation")
    if int(record["rows"]) != expected_rows:
        raise RuntimeError("completed PL2B work part row count mismatch")


def compute_parts(
    *,
    output: Path,
    specs: list[dict[str, object]],
    reactions: list[dict[str, str]],
    reaction_columns: np.ndarray,
    mats: np.ndarray,
    fingerprint: dict[str, object],
) -> tuple[Path, list[dict[str, object]]]:
    work = output.with_name(output.name + ".work_v1")
    manifest_path = work / "WORK_MANIFEST.json"
    if not work.exists():
        work.mkdir(parents=True)
        manifest = {"schema": WORK_SCHEMA, "fingerprint": fingerprint, "completed_parts": {}}
        atomic_json(manifest_path, manifest)
    if not manifest_path.is_file():
        raise RuntimeError("PL2B work directory is incomplete")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != WORK_SCHEMA or manifest.get("fingerprint") != fingerprint:
        raise RuntimeError("PL2B work directory fingerprint mismatch")
    completed = manifest.get("completed_parts", {})
    ranges = _part_ranges(len(specs))
    records: list[dict[str, object]] = []
    for index, (start, stop) in enumerate(ranges):
        name = _part_name(index)
        expected_rows = (stop - start) * len(reactions)
        if name in completed:
            record = completed[name]
            _validate_completed_part(work, record, expected_rows)
        else:
            record = write_case_part(work / name, specs[start:stop], reactions, reaction_columns, mats)
            if int(record["rows"]) != expected_rows:
                raise RuntimeError("new PL2B work part row count mismatch")
            completed[name] = record
            manifest["completed_parts"] = completed
            atomic_json(manifest_path, manifest)
        records.append(record)
    return work, records


def inspect_case_parts(
    directory: Path,
    *,
    expected_manifest: dict[str, object] | None = None,
    collect_summary: bool = False,
) -> tuple[dict[str, object], dict[str, object] | None]:
    expected_names = [_part_name(index) for index in range(PART_COUNT)]
    actual_names = sorted(path.name for path in directory.glob("BRIDGEPL2B_CASE_OUTCOMES.part-*.tsv.xz"))
    if actual_names != expected_names:
        raise RuntimeError("PL2B case parts are missing, extra, or not canonically named")
    logical = hashlib.sha256()
    records = []
    total = 0
    previous_key: tuple[str, str] | None = None
    gains = np.full((len(LAMBDA_PREFIX), EXPECTED_CASE_ROWS, 2), np.nan) if collect_summary else None
    context_codes = np.empty(EXPECTED_CASE_ROWS, dtype=np.uint8) if collect_summary else None
    context_to_code: dict[str, int] = {}
    non_tie = np.zeros(EXPECTED_CASE_ROWS, dtype=bool) if collect_summary else None
    gain_columns = []
    for prefix in LAMBDA_PREFIX.values():
        gain_columns.extend([f"{prefix}_correct_absolute_error_gain", f"{prefix}_wrong_absolute_error_gain"])
    for part_index, name in enumerate(expected_names):
        path = directory / name
        if path.stat().st_size >= MAX_PART_BYTES:
            raise RuntimeError("PL2B case part exceeds the fixed safety ceiling")
        count = 0
        with lzma.open(path, "rb", format=lzma.FORMAT_XZ) as handle:
            header_line = handle.readline()
            if not header_line:
                raise RuntimeError("empty PL2B case part")
            header = header_line.rstrip(b"\n").decode("utf-8").split("\t")
            if header != CASE_FIELDS:
                raise RuntimeError("PL2B case schema mismatch")
            if part_index == 0:
                logical.update(header_line)
            positions = {field: header.index(field) for field in [
                "truth_pair_id", "reaction_id", "algorithm", "rna_context_key", "truth_status", *gain_columns,
            ]}
            for line in handle:
                logical.update(line)
                fields = line.rstrip(b"\n").decode("utf-8").split("\t")
                if len(fields) != len(header):
                    raise RuntimeError("PL2B case row width mismatch")
                key = (fields[positions["truth_pair_id"]], fields[positions["reaction_id"]])
                if previous_key is not None and key <= previous_key:
                    raise RuntimeError("PL2B case rows are not strictly sorted")
                previous_key = key
                if collect_summary:
                    row_index = total + count
                    label = fields[positions["algorithm"]] + "\x1f" + fields[positions["rna_context_key"]]
                    if label not in context_to_code:
                        if len(context_to_code) >= 255:
                            raise RuntimeError("too many PL2B summary contexts")
                        context_to_code[label] = len(context_to_code)
                    context_codes[row_index] = context_to_code[label]
                    non_tie[row_index] = fields[positions["truth_status"]] == "NON_TIE"
                    for lam_index, prefix in enumerate(LAMBDA_PREFIX.values()):
                        gains[lam_index, row_index, 0] = float(fields[positions[f"{prefix}_correct_absolute_error_gain"]])
                        gains[lam_index, row_index, 1] = float(fields[positions[f"{prefix}_wrong_absolute_error_gain"]])
                count += 1
        record = {"filename": name, "rows": count, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
        records.append(record)
        total += count
    if total != EXPECTED_CASE_ROWS:
        raise RuntimeError(f"PL2B case row total is {total}, expected {EXPECTED_CASE_ROWS}")
    manifest = {
        "schema": PART_SCHEMA,
        "logical_artifact": "BRIDGEPL2B_CASE_OUTCOMES",
        "format": "tsv.xz.parts",
        "partition_scheme": "128_contiguous_canonical_truth_pair_ranges",
        "sort_key": ["truth_pair_id", "reaction_id"],
        "number_of_parts": PART_COUNT,
        "column_names": CASE_FIELDS,
        "parts": records,
        "total_data_rows": total,
        "logical_content_sha256": logical.hexdigest(),
        "maximum_part_bytes": max(record["bytes"] for record in records),
        "part_size_ceiling_bytes": MAX_PART_BYTES,
    }
    if expected_manifest is not None:
        for key in (
            "schema", "logical_artifact", "format", "partition_scheme", "sort_key",
            "number_of_parts", "column_names", "parts", "total_data_rows",
            "logical_content_sha256", "maximum_part_bytes", "part_size_ceiling_bytes",
        ):
            if expected_manifest.get(key) != manifest[key]:
                raise RuntimeError(f"PL2B case-parts manifest mismatch: {key}")
    summary_data = None if not collect_summary else {
        "gains": gains,
        "context_codes": context_codes,
        "context_to_code": context_to_code,
        "non_tie": non_tie,
    }
    return manifest, summary_data


def build_context_summary(summary_data: dict[str, object]) -> list[dict[str, object]]:
    gains = np.asarray(summary_data["gains"])
    context_codes = np.asarray(summary_data["context_codes"])
    context_to_code = dict(summary_data["context_to_code"])
    non_tie = np.asarray(summary_data["non_tie"], dtype=bool)
    contexts = sorted(context_to_code)
    scopes: list[tuple[str, str, str, np.ndarray]] = [
        ("OVERALL", "ALL", "ALL", np.ones(non_tie.size, dtype=bool))
    ]
    for label in contexts:
        algorithm, rna = label.split("\x1f", 1)
        scopes.append(("CONTEXT", algorithm, rna, context_codes == context_to_code[label]))
    rows: list[dict[str, object]] = []
    lambdas = list(LAMBDA_PREFIX)
    for scope, algorithm, rna, scope_mask in scopes:
        nt_mask = scope_mask & non_tie
        total_count = int(np.count_nonzero(scope_mask))
        non_tie_count = int(np.count_nonzero(nt_mask))
        tie_count = total_count - non_tie_count
        for lam_index, lam in enumerate(lambdas):
            correct = gains[lam_index, :, 0]
            wrong = gains[lam_index, :, 1]
            for q in core.RELIABILITY_GRID:
                expected = q * correct + (1.0 - q) * wrong
                values = expected[nt_mask]
                if values.size != non_tie_count or not np.all(np.isfinite(values)):
                    raise RuntimeError("non-tie reliability gains are incomplete or non-finite")
                improved = int(np.count_nonzero(values > core.GAIN_TOL))
                harmed = int(np.count_nonzero(values < -core.GAIN_TOL))
                tied = non_tie_count - improved - harmed
                quantiles = np.quantile(values, [0.10, 0.50, 0.90], method="linear") if values.size else np.full(3, np.nan)
                rows.append({
                    "summary_scope": scope,
                    "algorithm": algorithm,
                    "rna_context_key": rna,
                    "lambda_total": lam,
                    "lambda_role": "PRIMARY" if lam == core.PRIMARY_LAMBDA else "SENSITIVITY",
                    "reliability_q": q,
                    "total_distinct_reaction_truth_pair_cases": total_count,
                    "non_tie_count": non_tie_count,
                    "truth_tie_count": tie_count,
                    "improved_count": improved,
                    "improved_fraction_of_non_ties": improved / non_tie_count if non_tie_count else float("nan"),
                    "tied_count": tied,
                    "tied_fraction_of_non_ties": tied / non_tie_count if non_tie_count else float("nan"),
                    "harmed_count": harmed,
                    "harmed_fraction_of_non_ties": harmed / non_tie_count if non_tie_count else float("nan"),
                    "gain_median": quantiles[1],
                    "gain_q10": quantiles[0],
                    "gain_q90": quantiles[2],
                    "fraction_denominator": "NON_TIE_CASES",
                    "quantile_alias_weighting": "DISTINCT_TRUTH_PAIR_ONCE",
                })
    return rows


def _registry_fields() -> list[str]:
    return [
        "truth_pair_id", "evaluation_id", "algorithm", "rna_context_key", "ct2a_mouse",
        "gl261_mouse", "ct2a_candidate_id", "gl261_candidate_id", "ct2a_sample_index",
        "gl261_sample_index", "ct2a_post_holdout_mass", "gl261_post_holdout_mass",
        "ct2a_post_holdout_ess", "gl261_post_holdout_ess", "ct2a_holdout_status",
        "gl261_holdout_status", "pair_evaluable", "pair_status", "non_evaluable_reason",
        "alias_count", "quantile_aliases",
    ]


def _alias_fields(source_rows: list[dict[str, str]]) -> list[str]:
    return ["truth_pair_id", "alias_count", "quantile_aliases", *source_rows[0].keys()]


def _summary_fields() -> list[str]:
    return [
        "summary_scope", "algorithm", "rna_context_key", "lambda_total", "lambda_role",
        "reliability_q", "total_distinct_reaction_truth_pair_cases", "non_tie_count",
        "truth_tie_count", "improved_count", "improved_fraction_of_non_ties",
        "tied_count", "tied_fraction_of_non_ties", "harmed_count",
        "harmed_fraction_of_non_ties", "gain_median", "gain_q10", "gain_q90",
        "fraction_denominator", "quantile_alias_weighting",
    ]


def validate_existing_output(output: Path, fingerprint: dict[str, object]) -> dict[str, object]:
    manifest_path = output / "BRIDGEPL2B_MANIFEST.json"
    if not manifest_path.is_file():
        raise RuntimeError("refusing existing PL2B directory without a manifest")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != COMPLETE_STATUS or manifest.get("fingerprint") != fingerprint:
        raise RuntimeError("existing PL2B output identity/status mismatch")
    for name, expected in manifest.get("artifact_sha256", {}).items():
        path = output / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"existing PL2B artifact mismatch: {name}")
    parts_manifest = json.loads((output / "BRIDGEPL2B_CASE_OUTCOMES.parts.json").read_text(encoding="utf-8"))
    inspect_case_parts(output, expected_manifest=parts_manifest, collect_summary=False)
    return manifest


def produce(output: Path) -> dict[str, object]:
    admission = verify_admission()
    implementation_hashes = _implementation_hashes()
    fingerprint = {
        "schema": "bridge.pl2b.production_fingerprint.v1",
        "admission": admission,
        "implementation_sha256": implementation_hashes,
        "fixed_parameters": {
            "tie_tolerance": core.TIE_TOL,
            "gain_tolerance": core.GAIN_TOL,
            "lambda_grid": list(core.LAMBDA_GRID),
            "reliability_grid": list(core.RELIABILITY_GRID),
            "part_count": PART_COUNT,
            "reaction_block_size": REACTION_BLOCK_SIZE,
        },
    }
    if output.exists():
        manifest = validate_existing_output(output, fingerprint)
        return {"status": NOOP_STATUS, "manifest": manifest}

    alias_source = read_tsv(ROOT / PL2A_REL / "BRIDGEPL2A_TRUTH_PAIR_REGISTRY.tsv")
    distinct, alias_map = collapse_truth_pairs(alias_source)
    pl1_manifest, candidates, strong, data, mats, cache_info = pl2a.load_and_validate_inputs()
    del strong, pl1_manifest
    specs, holdout_qc = build_holdout_specs(distinct, candidates, data, mats, cache_info)
    reactions = read_tsv(ROOT / PL2A_REL / "BRIDGEPL2A_REACTION_REGISTRY.tsv")
    if len(reactions) != 4180 or [row["reaction_id"] for row in reactions] != sorted(row["reaction_id"] for row in reactions):
        raise RuntimeError("PL2A reaction registry identity/order mismatch")
    reaction_index = {reaction: index for index, reaction in enumerate(cache_info["rxns"])}
    if set(row["reaction_id"] for row in reactions) - set(reaction_index):
        raise RuntimeError("reaction registry is not covered by the frozen flux cache")
    reaction_columns = np.asarray([reaction_index[row["reaction_id"]] for row in reactions], dtype=int)
    work, work_records = compute_parts(
        output=output, specs=specs, reactions=reactions, reaction_columns=reaction_columns,
        mats=mats, fingerprint=fingerprint,
    )
    del mats

    with tempfile.TemporaryDirectory(prefix="dmi_bridge_pl2b_publish_", dir=str(output.parent)) as temp_name:
        temp = Path(temp_name)
        for record in work_records:
            shutil.copy2(work / str(record["filename"]), temp / str(record["filename"]))
        parts_manifest, summary_data = inspect_case_parts(temp, collect_summary=True)
        context_summary = build_context_summary(summary_data)
        del summary_data

        write_tsv(temp / "BRIDGEPL2B_DISTINCT_TRUTH_PAIR_REGISTRY.tsv", _registry_fields(), distinct)
        write_tsv(temp / "BRIDGEPL2B_ALIAS_MAP.tsv", _alias_fields(alias_source), alias_map)
        write_tsv(temp / "BRIDGEPL2B_CONTEXT_SUMMARY.tsv", _summary_fields(), context_summary)
        (temp / "BRIDGEPL2B_CASE_OUTCOMES.parts.json").write_bytes(json_bytes(parts_manifest))

        source_audit = {
            "schema": "bridge.pl2b.source_audit.v1",
            "status": "PL2B_SOURCE_AUDIT_COMPLETE",
            **admission,
            "implementation_sha256": implementation_hashes,
            "data_boundary": "all consumed data resolved to regular files inside the repository",
            "reaction_B_access_after_admission": True,
        }
        (temp / "BRIDGEPL2B_SOURCE_AUDIT.json").write_bytes(json_bytes(source_audit))
        analysis_contract = core.analysis_contract()
        analysis_contract.update({
            "production_contract_sha256": implementation_hashes[str(CONTRACT_REL)],
            "case_storage": {
                "representation": "correct_and_wrong_endpoints_only",
                "reliability_rows_materialized": False,
                "parts": PART_COUNT,
                "sort_key": ["truth_pair_id", "reaction_id"],
            },
            "summary_fraction_denominator": "non_tie_cases",
            "quantile_alias_summaries_produced": False,
            "pl1_predictors_duplicated_in_case_rows": False,
            "pl1_join_fields": ["evaluation_id", "reaction_id"],
        })
        (temp / "BRIDGEPL2B_ANALYSIS_CONTRACT.json").write_bytes(json_bytes(analysis_contract))

        lambda_zero_max = max(float(record["lambda_zero_maximum_deviation"]) for record in work_records)
        scalar_batch_max = max(float(record["scalar_batch_maximum_deviation"]) for record in work_records)
        scalar_batch_break_even_max = max(float(record["scalar_batch_break_even_maximum_deviation"]) for record in work_records)
        scalar_batch_break_even_gain_scale_max = max(
            float(record["scalar_batch_break_even_gain_scale_maximum_deviation"])
            for record in work_records
        )
        total_ties = next(
            int(row["truth_tie_count"])
            for row in context_summary
            if row["summary_scope"] == "OVERALL" and row["lambda_total"] == 0.25 and row["reliability_q"] == 0.5
        )
        qc = {
            "schema": "bridge.pl2b.qc.v1",
            "status": "PASS",
            "counts": {
                "quantile_aliases": 1200,
                "distinct_truth_pairs": 885,
                "distinct_evaluable_truth_pairs": 836,
                "distinct_non_evaluable_truth_pairs": 49,
                "reactions": 4180,
                "case_rows": parts_manifest["total_data_rows"],
                "truth_tie_case_rows": total_ties,
                "non_tie_case_rows": EXPECTED_CASE_ROWS - total_ties,
                "context_summary_rows": len(context_summary),
            },
            "alias_multiplicity": {
                "evaluable": {"1": 642, "2": 114, "3": 80},
                "non_evaluable": {"1": 24, "2": 9, "3": 16},
            },
            "holdout": {**holdout_qc, "tolerance": 1e-12, "status": "PASS"},
            "cartesian_support": "separate_19_candidate_condition_pools_crossed_as_19x19",
            "lambda_zero": {"maximum_deviation": lambda_zero_max, "tolerance": 1e-12, "status": "PASS"},
            "scalar_batch_sentinel": {
                "endpoint_and_baseline_maximum_deviation": scalar_batch_max,
                "endpoint_and_baseline_tolerance": 5e-12,
                "break_even_maximum_deviation": scalar_batch_break_even_max,
                "break_even_gain_scale_maximum_deviation": scalar_batch_break_even_gain_scale_max,
                "break_even_gain_scale_tolerance": 5e-12,
                "break_even_note": "raw q deviation is reported; the hard gate scales it by the endpoint-gain slope because near-zero slopes are condition-sensitive",
                "status": "PASS",
            },
            "storage": {
                "number_of_parts": PART_COUNT,
                "maximum_part_bytes": parts_manifest["maximum_part_bytes"],
                "part_size_ceiling_bytes": MAX_PART_BYTES,
                "logical_content_sha256": parts_manifest["logical_content_sha256"],
                "status": "PASS",
            },
            "forbidden_work": {
                "solver_invoked": False,
                "optimization_invoked": False,
                "fva_invoked": False,
                "sampling_invoked": False,
                "reconstruction_invoked": False,
                "new_strong_anchor_weights_generated": False,
                "gene_score_concordance_used": False,
                "historical_tanh_outcomes_used": False,
                "historical_bridge3_outcomes_used": False,
                "outcome_driven_reaction_selection": False,
                "outcome_driven_lambda_tuning": False,
                "outcome_driven_pl1_threshold_fitting": False,
                "new_pl1_predictors_constructed": False,
            },
        }
        (temp / "BRIDGEPL2B_QC.json").write_bytes(json_bytes(qc))

        artifact_names = [
            "BRIDGEPL2B_SOURCE_AUDIT.json", "BRIDGEPL2B_ANALYSIS_CONTRACT.json",
            "BRIDGEPL2B_DISTINCT_TRUTH_PAIR_REGISTRY.tsv", "BRIDGEPL2B_ALIAS_MAP.tsv",
            "BRIDGEPL2B_CASE_OUTCOMES.parts.json", "BRIDGEPL2B_CONTEXT_SUMMARY.tsv",
            "BRIDGEPL2B_QC.json", *[_part_name(index) for index in range(PART_COUNT)],
        ]
        artifact_hashes = {name: sha256_file(temp / name) for name in artifact_names}
        manifest = {
            "schema": "bridge.pl2b.sign_only_utility.manifest.v1",
            "status": COMPLETE_STATUS,
            "fingerprint": fingerprint,
            "counts": qc["counts"],
            "artifact_sha256": artifact_hashes,
            "case_logical_content_sha256": parts_manifest["logical_content_sha256"],
            "case_parts_manifest_sha256": artifact_hashes["BRIDGEPL2B_CASE_OUTCOMES.parts.json"],
            "primary_endpoint": "baseline_abs_error_minus_posterior_abs_error_for_E_abs_delta_B",
            "primary_lambda": core.PRIMARY_LAMBDA,
            "primary_benchmark_unit": "distinct_truth_pair_once",
            "reliability_storage": "correct_and_wrong_endpoints_only",
            "pl1_pl2a_outputs_modified": False,
            **qc["forbidden_work"],
        }
        (temp / "BRIDGEPL2B_MANIFEST.json").write_bytes(json_bytes(manifest))
        output.parent.mkdir(parents=True, exist_ok=True)
        os.replace(temp, output)
    validate_existing_output(output, fingerprint)
    return {"status": COMPLETE_STATUS, "manifest": manifest}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / OUTPUT_REL)
    args = parser.parse_args()
    result = produce(args.output.resolve())
    print(json.dumps({
        "status": result["status"],
        "counts": result["manifest"]["counts"],
        "case_logical_content_sha256": result["manifest"]["case_logical_content_sha256"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
