#!/usr/bin/env python3
"""Resumable production adapter for DMI-BRIDGE-PL1.

This module changes execution only.  Scientific quantities remain defined by
``target_landscape_metrics`` and the validated block accelerator used by
``dmi_bridge_pl1_predictability_landscape_v1``.  The work directory is a
checkpoint cache, never a published PL1 result.
"""
from __future__ import annotations

import argparse
import csv
import json
import lzma
import os
import shutil
from pathlib import Path

import numpy as np

from scripts import dmi_bridge_pl1_predictability_landscape_v1 as base
from scripts import dmi_bridge_pl1_storage_v1 as storage
from scripts.dmi_bridge_pl1_predictability_batch_v1 import target_landscape_metrics_block
from scripts.dmi_bridge_pl1_predictability_core_v1 import (
    COMPLETE_STATUS,
    STRONG_ANCHOR_REACTION,
    analysis_contract,
    target_landscape_metrics,
)

WORK_SCHEMA = "bridge.pl1.predictability_landscape.resume_work.v1"
RESUME_SCHEMA = "bridge.pl1.predictability_landscape.resume_execution.v1"
CONTEXT_COUNT = 16
EVALS_PER_CONTEXT = 25
FEATURE_ROWS_PER_CONTEXT = 104_500
WORK_SUFFIX = ".resume_work_v1"


def _atomic_json(path: Path, obj: dict) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_bytes(json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=True).encode("ascii") + b"\n")
    os.replace(tmp, path)


def _code_hashes() -> dict[str, str]:
    paths = (
        Path(__file__).resolve(),
        Path(base.__file__).resolve(),
        base.ROOT / "scripts/dmi_bridge_pl1_storage_v1.py",
        base.ROOT / "scripts/dmi_bridge_pl1_predictability_batch_v1.py",
        base.ROOT / "scripts/dmi_bridge_pl1_predictability_core_v1.py",
    )
    return {str(p.relative_to(base.ROOT)): base.sha256(p) for p in paths}


def _metric_layout() -> tuple[list[str], list[str], dict[str, int]]:
    probe = target_landscape_metrics(
        anchor_ct2a=[0], anchor_gl261=[0], target_ct2a=[0], target_gl261=[0],
        strong_weights_ct2a=[1], strong_weights_gl261=[1],
    )
    metric_fields = [
        k for k, v in probe.items()
        if isinstance(v, (int, float, np.number)) and not isinstance(v, bool)
    ]
    status_fields = [
        "target_delta_degenerate",
        "target_abs_magnitude_degenerate",
        "anchor_delta_degenerate",
        "sign_magnitude_status",
        "anchor_target_coupling_status",
    ]
    return metric_fields, status_fields, {m: i for i, m in enumerate(metric_fields)}


def _context_specs(by_key: dict, wm: dict, cache_index: dict) -> tuple[list[dict], list[dict]]:
    specs: list[dict] = []
    registry: list[dict] = []
    stratum_keys = sorted(by_key, key=lambda k: (k[0], k[1], by_key[k][0]["rna_context_key"]))
    eval_index = 0
    for key in stratum_keys:
        method, tumor, ensemble = key
        if tumor != "CT2A":
            continue
        rows = by_key[key]
        rna = rows[0]["rna_context_key"]
        if any(r["rna_context_key"] != rna for r in rows):
            raise RuntimeError("candidate stratum has mixed RNA identities")
        other_keys = [
            k for k in stratum_keys
            if k[0] == method and k[1] == "GL261" and by_key[k][0]["rna_context_key"] == rna
        ]
        if len(other_keys) != 1:
            raise RuntimeError("paired tumor stratum is missing")
        gkey = other_keys[0]
        grows = by_key[gkey]
        ge = grows[0]["ensemble_hash"]
        cmice = sorted({k[0] for k in wm if k[1] == "CT2A" and k[2] == method and k[3] == ensemble and k[4] == rna})
        gmice = sorted({k[0] for k in wm if k[1] == "GL261" and k[2] == method and k[3] == ge and k[4] == rna})
        if len(cmice) != 5 or len(gmice) != 5:
            raise RuntimeError("expected five mice per tumor")
        cpos = len(specs)
        specs.append({
            "context_index": cpos,
            "method": method,
            "rna_context_key": rna,
            "ct2a_key": key,
            "gl261_key": gkey,
            "ct2a_cache_index": cache_index[key],
            "gl261_cache_index": cache_index[gkey],
            "ct2a_ensemble": ensemble,
            "gl261_ensemble": ge,
            "ct2a_mice": cmice,
            "gl261_mice": gmice,
            "eval_start": eval_index,
            "eval_stop": eval_index + EVALS_PER_CONTEXT,
        })
        for cmouse in cmice:
            for gmouse in gmice:
                registry.append({
                    "evaluation_id": f"E{eval_index:04d}",
                    "algorithm": method,
                    "rna_context_key": rna,
                    "ct2a_mouse": cmouse,
                    "gl261_mouse": gmouse,
                    "cartesian_support": 400,
                    "joint_weight_sum": 1.0,
                    "status": "ADMITTED",
                })
                eval_index += 1
    if len(specs) != CONTEXT_COUNT or len(registry) != 400:
        raise RuntimeError(f"unexpected PL1 context/evaluation count: {len(specs)}/{len(registry)}")
    return specs, registry


def _work_paths(out: Path) -> dict[str, Path]:
    work = out.with_name(out.name + WORK_SUFFIX)
    return {
        "dir": work,
        "manifest": work / "WORK_MANIFEST.json",
        "registry": work / "BRIDGEPL1_EVALUATION_REGISTRY.tsv",
        "metrics": work / "metrics.npy.xz",
        "legacy_metrics": work / "metrics.npy",
        "bools": work / "bool_metrics.npy",
        "statuses": work / "status_metrics.npy",
        "shards": work / "feature_shards",
    }


def _initial_work_manifest(source_hashes: dict, code_hashes: dict, metric_fields: list[str]) -> dict:
    return {
        "schema": WORK_SCHEMA,
        "source_sha256": source_hashes,
        "code_sha256": code_hashes,
        "metric_fields": metric_fields,
        "metrics_storage": "metrics.npy.xz",
        "completed_contexts": {},
        "counts": {
            "contexts": CONTEXT_COUNT,
            "evaluations": 400,
            "primary_reactions": 4180,
            "feature_rows_per_context": FEATURE_ROWS_PER_CONTEXT,
        },
    }


def _open_or_create_work(
    *, out: Path, source_hashes: dict, code_hashes: dict, registry: list[dict], metric_fields: list[str]
):
    paths = _work_paths(out)
    work = paths["dir"]
    shape = (400, 4180, len(metric_fields))
    if not work.exists():
        work.mkdir(parents=True)
        paths["shards"].mkdir()
        base.write_tsv(paths["registry"], list(registry[0]), registry)
        metric_store = storage.NpyXZWork(paths["metrics"], paths["legacy_metrics"], shape=shape, dtype="f8")
        metrics = metric_store.array
        bools = np.lib.format.open_memmap(paths["bools"], mode="w+", dtype="u1", shape=(400, 4180, 3))
        statuses = np.lib.format.open_memmap(paths["statuses"], mode="w+", dtype="u1", shape=(400, 4180, 2))
        metrics.flush(); bools.flush(); statuses.flush()
        metric_store.checkpoint()
        paths["legacy_metrics"].unlink(missing_ok=True)
        _atomic_json(paths["manifest"], _initial_work_manifest(source_hashes, code_hashes, metric_fields))
    if not all(paths[k].is_file() for k in ("manifest", "registry", "bools", "statuses")) or not (paths["metrics"].is_file() or paths["legacy_metrics"].is_file()):
        raise RuntimeError("PL1 resume work directory is incomplete")
    paths["shards"].mkdir(exist_ok=True)
    manifest = json.loads(paths["manifest"].read_text())
    if manifest.get("schema") != WORK_SCHEMA:
        raise RuntimeError("PL1 resume work schema mismatch")
    if manifest.get("source_sha256") != source_hashes:
        raise RuntimeError("PL1 resume work source hashes do not match current frozen inputs")
    if manifest.get("code_sha256") != code_hashes:
        raise RuntimeError("PL1 resume work code hashes do not match current implementation")
    if manifest.get("metric_fields") != metric_fields:
        raise RuntimeError("PL1 resume work metric layout mismatch")
    current_registry = paths["registry"].read_bytes()
    tmp_registry = work / ".registry_check.tsv"
    base.write_tsv(tmp_registry, list(registry[0]), registry)
    expected_registry = tmp_registry.read_bytes()
    tmp_registry.unlink()
    if current_registry != expected_registry:
        raise RuntimeError("PL1 resume registry mismatch")
    metric_store = storage.NpyXZWork(paths["metrics"], paths["legacy_metrics"])
    metrics = metric_store.array
    if paths["legacy_metrics"].is_file() and paths["metrics"].is_file():
        paths["legacy_metrics"].unlink()
    bools = np.load(paths["bools"], mmap_mode="r+")
    statuses = np.load(paths["statuses"], mmap_mode="r+")
    if metrics.shape != shape or bools.shape != (400, 4180, 3) or statuses.shape != (400, 4180, 2):
        raise RuntimeError("PL1 resume array shape mismatch")
    return paths, manifest, metric_store, metrics, bools, statuses


def _feature_fields(metric_fields: list[str], status_fields: list[str]) -> list[str]:
    return [
        "evaluation_id", "algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse",
        "reaction_id", "subsystem", "same_subsystem_as_HEX1", "proximal_set_member",
    ] + metric_fields + status_fields


def _write_feature_shard(
    *, path: Path, include_header: bool, eval_start: int, eval_stop: int, registry: list[dict],
    primary_rxns: list[str], panel: dict[str, str], proximal: set[str], metric_fields: list[str],
    status_fields: list[str], metric_pos: dict[str, int], all_metrics: np.ndarray,
    bool_metrics: np.ndarray, status_metrics: np.ndarray,
) -> int:
    tmp = path.with_name(path.name + ".tmp")
    tmp.unlink(missing_ok=True)
    fields = _feature_fields(metric_fields, status_fields)
    hex_subsystem = panel.get(STRONG_ANCHOR_REACTION)
    rows_written = 0
    with lzma.open(tmp, "wt", encoding="utf-8", newline="", preset=base.XZ_PRESET) as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        if include_header:
            writer.writeheader()
        for e in range(eval_start, eval_stop):
            reg = registry[e]
            for rpos, rid in enumerate(primary_rxns):
                subsystem = panel.get(rid, "UNMAPPED_PATHWAY")
                row = {
                    "evaluation_id": reg["evaluation_id"],
                    "algorithm": reg["algorithm"],
                    "rna_context_key": reg["rna_context_key"],
                    "ct2a_mouse": reg["ct2a_mouse"],
                    "gl261_mouse": reg["gl261_mouse"],
                    "reaction_id": rid,
                    "subsystem": subsystem,
                    "same_subsystem_as_HEX1": subsystem == hex_subsystem,
                    "proximal_set_member": subsystem in proximal,
                }
                row.update({m: all_metrics[e, rpos, metric_pos[m]] for m in metric_fields})
                row.update({m: bool_metrics[e, rpos, j] for j, m in enumerate(status_fields[:3])})
                row.update({
                    m: ("OK" if status_metrics[e, rpos, j] else "DEGENERATE")
                    for j, m in enumerate(status_fields[3:])
                })
                writer.writerow({k: base.f(row[k]) for k in fields})
                rows_written += 1
    os.replace(tmp, path)
    return rows_written


def _compute_context(
    *, spec: dict, wm: dict, mats: np.ndarray, primary_rxns: list[str], rxn_index: dict[str, int],
    hex_idx: int, metric_fields: list[str], status_fields: list[str], metric_pos: dict[str, int],
    all_metrics: np.ndarray, bool_metrics: np.ndarray, status_metrics: np.ndarray,
) -> None:
    ci, gi = spec["ct2a_cache_index"], spec["gl261_cache_index"]
    method, rna = spec["method"], spec["rna_context_key"]
    anchor_c = mats[ci, :, hex_idx]
    anchor_g = mats[gi, :, hex_idx]
    local_eval = spec["eval_start"]
    for cmouse in spec["ct2a_mice"]:
        wc = np.asarray([
            w for _, w in sorted(wm[(cmouse, "CT2A", method, spec["ct2a_ensemble"], rna)])
        ], dtype=float)
        for gmouse in spec["gl261_mice"]:
            wg = np.asarray([
                w for _, w in sorted(wm[(gmouse, "GL261", method, spec["gl261_ensemble"], rna)])
            ], dtype=float)
            for start in range(0, len(primary_rxns), base.REACTION_BLOCK_SIZE):
                stop = min(start + base.REACTION_BLOCK_SIZE, len(primary_rxns))
                block_rxns = primary_rxns[start:stop]
                cols = [rxn_index[r] for r in block_rxns]
                target_c = mats[ci][:, cols]
                target_g = mats[gi][:, cols]
                block = target_landscape_metrics_block(
                    anchor_ct2a=anchor_c,
                    anchor_gl261=anchor_g,
                    target_ct2a=target_c,
                    target_gl261=target_g,
                    strong_weights_ct2a=wc,
                    strong_weights_gl261=wg,
                )
                for name in metric_fields:
                    all_metrics[local_eval, start:stop, metric_pos[name]] = np.asarray(block[name], dtype=float)
                for j, name in enumerate(status_fields[:3]):
                    bool_metrics[local_eval, start:stop, j] = np.asarray(block[name], dtype=bool)
                for j, name in enumerate(status_fields[3:]):
                    status_metrics[local_eval, start:stop, j] = np.asarray(
                        block[name] == np.asarray("OK"), dtype=np.uint8
                    )
                if cmouse == spec["ct2a_mice"][0] and gmouse == spec["gl261_mice"][0]:
                    base._assert_scalar_parity(
                        batch=block,
                        local_pos=0,
                        anchor_ct2a=anchor_c,
                        anchor_gl261=anchor_g,
                        target_ct2a=target_c[:, 0],
                        target_gl261=target_g[:, 0],
                        strong_weights_ct2a=wc,
                        strong_weights_gl261=wg,
                    )
            local_eval += 1
    if local_eval != spec["eval_stop"]:
        raise RuntimeError("context evaluation span mismatch")
    all_metrics.flush(); bool_metrics.flush(); status_metrics.flush()


def _mapped_summary_arrays(all_metrics: np.ndarray, eval_slice: slice, metric_pos: dict[str, int]):
    out = {}
    for metric in base.REACTION_STATS:
        x = np.asarray(all_metrics[eval_slice, :, metric_pos[metric]], dtype=float)
        finite = np.isfinite(x)
        count = np.sum(finite, axis=0)
        with np.errstate(invalid="ignore"), np.testing.suppress_warnings() as sup:
            sup.filter(RuntimeWarning)
            q10, med, q90 = np.nanquantile(x, (0.10, 0.50, 0.90), axis=0, method="linear")
        out[metric] = (med, q10, q90, count, x.shape[0] - count)
    return out


def _concat_files(paths: list[Path], dest: Path) -> None:
    tmp = dest.with_name(dest.name + ".tmp")
    tmp.unlink(missing_ok=True)
    with tmp.open("wb") as out:
        for p in paths:
            with p.open("rb") as fh:
                shutil.copyfileobj(fh, out, length=1 << 20)
    os.replace(tmp, dest)


def _iter_feature_shard_rows(shard_paths: list[Path], fields: list[str]):
    for index, path in enumerate(shard_paths):
        with lzma.open(path, "rt", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(
                fh,
                fieldnames=None if index == 0 else fields,
                delimiter="\t",
            )
            for row in reader:
                yield row


def _validate_completed_shards(paths: dict, manifest: dict) -> None:
    completed = manifest.get("completed_contexts", {})
    for i in range(CONTEXT_COUNT):
        rec = completed.get(str(i))
        if not rec:
            raise RuntimeError(f"PL1 context {i} is not checkpoint-complete")
        p = paths["shards"] / rec["file"]
        if not p.is_file() or base.sha256(p) != rec["sha256"]:
            raise RuntimeError(f"PL1 context {i} feature shard is missing or changed")
        if int(rec["rows"]) != FEATURE_ROWS_PER_CONTEXT:
            raise RuntimeError(f"PL1 context {i} feature row count mismatch")


def _existing_complete_noop(out: Path, source_hashes: dict) -> dict | None:
    if not out.exists():
        return None
    if not out.is_dir():
        raise RuntimeError("PL1 output path exists and is not a directory")
    manifest_path = out / "BRIDGEPL1_MANIFEST.json"
    if not manifest_path.is_file():
        raise RuntimeError("refusing to overwrite existing non-complete PL1 directory")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("status") != COMPLETE_STATUS or manifest.get("source_sha256") != source_hashes:
        raise RuntimeError("refusing to overwrite differing PL1 outputs")
    artifact_hashes = manifest.get("artifact_sha256", {})
    if not artifact_hashes:
        raise RuntimeError("existing PL1 output lacks resumable artifact hash audit")
    for name, expected in artifact_hashes.items():
        p = out / name
        if not p.is_file() or base.sha256(p) != expected:
            raise RuntimeError("existing PL1 output artifact hash mismatch")
    storage.discover_feature_artifact(out)
    return manifest


def build(out: Path = base.OUT) -> dict:
    candidates, by_key, wm, rxns, mats, cache_index, hashes = base.load_inputs()
    noop = _existing_complete_noop(out, hashes)
    if noop is not None:
        return noop

    panel_rows = base.read_tsv(base.ROOT / "outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_REACTION_PANEL.tsv")
    panel = {r["reaction_id"]: r["subsystem"] for r in panel_rows}
    proximal = {"Glycolysis/gluconeogenesis", "Pyruvate metabolism", "Citric acid cycle", "Glutamate metabolism"}
    hex_idx = int(np.flatnonzero(rxns == STRONG_ANCHOR_REACTION)[0])
    primary_rxns = sorted(r for r in rxns.tolist() if r != STRONG_ANCHOR_REACTION)
    rxn_index = {r: i for i, r in enumerate(rxns.tolist())}
    metric_fields, status_fields, metric_pos = _metric_layout()
    specs, registry = _context_specs(by_key, wm, cache_index)
    code_hashes = _code_hashes()
    paths, work_manifest, metric_store, all_metrics, bool_metrics, status_metrics = _open_or_create_work(
        out=out,
        source_hashes=hashes,
        code_hashes=code_hashes,
        registry=registry,
        metric_fields=metric_fields,
    )

    for spec in specs:
        idx = spec["context_index"]
        rec = work_manifest.get("completed_contexts", {}).get(str(idx))
        if rec:
            shard = paths["shards"] / rec["file"]
            if shard.is_file() and base.sha256(shard) == rec.get("sha256") and int(rec.get("rows", -1)) == FEATURE_ROWS_PER_CONTEXT:
                continue
            raise RuntimeError(f"checkpoint for context {idx} is inconsistent")
        _compute_context(
            spec=spec, wm=wm, mats=mats, primary_rxns=primary_rxns, rxn_index=rxn_index,
            hex_idx=hex_idx, metric_fields=metric_fields, status_fields=status_fields,
            metric_pos=metric_pos, all_metrics=all_metrics, bool_metrics=bool_metrics,
            status_metrics=status_metrics,
        )
        shard_name = f"features_context_{idx:02d}.tsv.xz"
        shard = paths["shards"] / shard_name
        rows = _write_feature_shard(
            path=shard,
            include_header=(idx == 0),
            eval_start=spec["eval_start"], eval_stop=spec["eval_stop"], registry=registry,
            primary_rxns=primary_rxns, panel=panel, proximal=proximal,
            metric_fields=metric_fields, status_fields=status_fields, metric_pos=metric_pos,
            all_metrics=all_metrics, bool_metrics=bool_metrics, status_metrics=status_metrics,
        )
        if rows != FEATURE_ROWS_PER_CONTEXT:
            raise RuntimeError("feature shard row count mismatch")
        work_manifest.setdefault("completed_contexts", {})[str(idx)] = {
            "file": shard_name,
            "rows": rows,
            "sha256": base.sha256(shard),
            "evaluation_start": spec["eval_start"],
            "evaluation_stop": spec["eval_stop"],
        }
        metric_store.checkpoint()
        _atomic_json(paths["manifest"], work_manifest)

    _validate_completed_shards(paths, work_manifest)
    if base.source_hashes() != hashes:
        raise RuntimeError("source changed during PL1 resumable computation")

    tmp = out.with_name(out.name + ".publish_tmp")
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    try:
        base.write_tsv(tmp / "BRIDGEPL1_EVALUATION_REGISTRY.tsv", list(registry[0]), registry)
        shard_paths = [paths["shards"] / work_manifest["completed_contexts"][str(i)]["file"] for i in range(CONTEXT_COUNT)]
        feature_fields = _feature_fields(metric_fields, status_fields)
        feature_manifest = storage.write_partitioned_tsv(
            tmp,
            feature_fields,
            _iter_feature_shard_rows(shard_paths, feature_fields),
            total_rows=CONTEXT_COUNT * FEATURE_ROWS_PER_CONTEXT,
        )

        cfields = ["algorithm", "rna_context_key", "reaction_id"] + [
            f"{m}_{s}" for m in base.REACTION_STATS for s in ("median", "q10", "q90", "finite_count", "nonfinite_count")
        ] + list(base.DEGENERACY_COUNT_FIELDS)
        cfh, cwriter = base.open_xz_tsv(tmp / "BRIDGEPL1_CONTEXT_REACTION_SUMMARY.tsv.xz", cfields)
        try:
            for spec in specs:
                summaries = _mapped_summary_arrays(
                    all_metrics, slice(spec["eval_start"], spec["eval_stop"]), metric_pos
                )
                degeneracy = base._degeneracy_count_arrays(
                    np.asarray(bool_metrics[spec["eval_start"]:spec["eval_stop"], :, :]),
                    np.asarray(status_metrics[spec["eval_start"]:spec["eval_stop"], :, :]),
                )
                for rpos, rid in enumerate(primary_rxns):
                    row = base._summary_row_from_arrays(
                        {"algorithm": spec["method"], "rna_context_key": spec["rna_context_key"], "reaction_id": rid},
                        summaries,
                        rpos,
                    )
                    for name in base.DEGENERACY_COUNT_FIELDS:
                        row[name] = int(degeneracy[name][rpos])
                    cwriter.writerow({k: base.f(row[k]) for k in cfields})
        finally:
            cfh.close()

        reaction_summaries = _mapped_summary_arrays(all_metrics, slice(0, 400), metric_pos)
        reaction_degeneracy = base._degeneracy_count_arrays(np.asarray(bool_metrics), np.asarray(status_metrics))
        rrows = []
        for rpos, rid in enumerate(primary_rxns):
            row = base._summary_row_from_arrays({"reaction_id": rid}, reaction_summaries, rpos)
            for name in base.DEGENERACY_COUNT_FIELDS:
                row[name] = int(reaction_degeneracy[name][rpos])
            rrows.append(row)
        rfields = ["reaction_id"] + [
            f"{m}_{s}" for m in base.REACTION_STATS for s in ("median", "q10", "q90", "finite_count", "nonfinite_count")
        ] + list(base.DEGENERACY_COUNT_FIELDS)
        base.write_tsv(tmp / "BRIDGEPL1_REACTION_SUMMARY.tsv", rfields, rrows)

        pathway_vectors: dict[str, list[list[float]]] = {}
        for rpos, rid in enumerate(primary_rxns):
            subsystem = panel.get(rid, "UNMAPPED_PATHWAY")
            pathway_vectors.setdefault(subsystem, []).append([rrows[rpos][f"{m}_median"] for m in base.REACTION_STATS])
        prows = [
            base._stats_row({"subsystem": subsystem}, np.asarray(vectors, dtype=float))
            for subsystem, vectors in sorted(pathway_vectors.items())
        ]
        pfields = ["subsystem"] + [
            f"{m}_{s}" for m in base.REACTION_STATS for s in ("median", "q10", "q90", "finite_count", "nonfinite_count")
        ]
        base.write_tsv(tmp / "BRIDGEPL1_PATHWAY_SUMMARY.tsv", pfields, prows)

        stat_pos = {m: i for i, m in enumerate(base.REACTION_STATS)}
        sensitivity = []
        hex_subsystem = panel.get(STRONG_ANCHOR_REACTION)
        for scope, predicate in (
            ("primary", lambda same, prox: True),
            ("exclude_same_subsystem", lambda same, prox: not same),
            ("exclude_proximal_set", lambda same, prox: not prox),
        ):
            positions = []
            for i, rid in enumerate(primary_rxns):
                subsystem = panel.get(rid, "UNMAPPED_PATHWAY")
                if predicate(subsystem == hex_subsystem, subsystem in proximal):
                    positions.append(i)
            selected = np.asarray(
                [[rrows[i][f"{metric}_median"] for metric in base.REACTION_STATS] for i in positions], dtype=float
            )
            row = {"scope": scope, "reaction_count": len(positions), "feature_count": 400 * len(positions)}
            for metric in base.REACTION_STATS:
                s = base.qstats(selected[:, stat_pos[metric]].tolist())
                for k, v in s.items():
                    row[f"{metric}_{k}"] = v
            sensitivity.append(row)
        base.write_tsv(tmp / "BRIDGEPL1_SENSITIVITY_SUMMARY.tsv", list(sensitivity[0]), sensitivity)

        row_counts = {
            "evaluation_registry": 400,
            "reaction_evaluation_features": 1_672_000,
            "context_reaction_summary": 66_880,
            "reaction_summary": 4_180,
            "pathway_summary": len(prows),
        }
        source_audit = {
            "schema": "bridge.pl1.source_audit.v1",
            "status": COMPLETE_STATUS,
            "source_sha256": hashes,
            "counts": row_counts,
        }
        (tmp / "BRIDGEPL1_SOURCE_AUDIT.json").write_bytes(
            json.dumps(source_audit, sort_keys=True, indent=2, ensure_ascii=True).encode("ascii") + b"\n"
        )
        if base.source_hashes() != hashes:
            raise RuntimeError("source changed before PL1 publication")
        artifact_names = [
            "BRIDGEPL1_SOURCE_AUDIT.json",
            "BRIDGEPL1_EVALUATION_REGISTRY.tsv",
            storage.FEATURE_MANIFEST,
            *[p["filename"] for p in feature_manifest["parts"]],
            "BRIDGEPL1_CONTEXT_REACTION_SUMMARY.tsv.xz",
            "BRIDGEPL1_REACTION_SUMMARY.tsv",
            "BRIDGEPL1_PATHWAY_SUMMARY.tsv",
            "BRIDGEPL1_SENSITIVITY_SUMMARY.tsv",
        ]
        artifact_hashes = {name: base.sha256(tmp / name) for name in artifact_names}
        manifest = {
            "schema": base.CORE_SCHEMA,
            "status": COMPLETE_STATUS,
            "pl2_status": "NOT_RUN",
            "row_counts": row_counts,
            "source_sha256": hashes,
            "analysis_contract": analysis_contract(),
            "reaction_inventory": {
                "cache_reactions": 4181,
                "primary_reactions": 4180,
                "hex1_index": hex_idx,
                "panel_reactions": 3893,
                "unmapped_reactions": 288,
            },
            "execution": {
                "profile": RESUME_SCHEMA,
                "reaction_block_size": base.REACTION_BLOCK_SIZE,
                "xz_preset": base.XZ_PRESET,
                "feature_storage": "partitioned_tsv_v1",
                "feature_parts": storage.FEATURE_PARTS,
                "feature_rows_per_stream": FEATURE_ROWS_PER_CONTEXT,
                "scalar_authority": "target_landscape_metrics",
                "scalar_parity_rtol": base.SCALAR_PARITY_RTOL,
                "scalar_parity_atol": base.SCALAR_PARITY_ATOL,
                "checkpoint_work_schema": WORK_SCHEMA,
                "code_sha256": code_hashes,
            },
            "artifact_sha256": artifact_hashes,
            "publication": "PL1_PREDICTABILITY_LANDSCAPE_COMPLETE",
        }
        (tmp / "BRIDGEPL1_MANIFEST.json").write_bytes(
            json.dumps(manifest, sort_keys=True, indent=2, ensure_ascii=True).encode("ascii") + b"\n"
        )
        if out.exists():
            raise RuntimeError("refusing to overwrite existing PL1 outputs")
        metric_store.close()
        os.replace(tmp, out)
        return manifest
    except Exception:
        shutil.rmtree(tmp, ignore_errors=True)
        raise


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=base.OUT)
    ap.add_argument(
        "--cleanup-work",
        action="store_true",
        help="remove the checkpoint work directory only after a complete output has been validated",
    )
    args = ap.parse_args()
    manifest = build(args.output)
    print(json.dumps(manifest, sort_keys=True))
    if args.cleanup_work:
        work = _work_paths(args.output)["dir"]
        if args.output.is_dir() and (args.output / "BRIDGEPL1_MANIFEST.json").is_file():
            shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
