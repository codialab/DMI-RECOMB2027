#!/usr/bin/env python3
"""DMI-BRIDGE-3A configuration registry and operator-qualification gate.

This module implements BRIDGE C0/C1 only:
  C0. qualify measurement-specific strong/weak operators on one frozen panel;
  C1. freeze the complete 3^3 registry and the six primary incremental-weak
      contrasts.

It deliberately performs no posterior evaluation and imports no solver stack.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Any, Iterable

OBSERVABLES = ("M", "L", "G")
STATES = ("0", "W", "S")
SCHEMA_VERSION = "bridge3a.operator_audit.v1"
REGISTRY_VERSION = "bridge3a.registry.v1"

# Ordering is M=Vmax, L=Vlac, G=Vglx.
PRIMARY_MIXED_CONFIGS = ("SW0", "S0W", "WS0", "0SW", "W0S", "0WS")
PRIMARY_CONTRASTS = {
    "SW0": "S00",  # weak L added to strong M
    "S0W": "S00",  # weak G added to strong M
    "WS0": "0S0",  # weak M added to strong L
    "0SW": "0S0",  # weak G added to strong L
    "W0S": "00S",  # weak M added to strong G
    "0WS": "00S",  # weak L added to strong G
}

EXPECTED_INFORMATION_LEVEL = {
    "W": "weak_direction_only",
    "S": "strong_quantitative",
}
ALLOWED_OPERATOR_STATUS = {"QUALIFIED", "UNAVAILABLE", "BLOCKED"}


@dataclass(frozen=True)
class ConfigRow:
    config: str
    M: str
    L: str
    G: str
    n_absent: int
    n_weak: int
    n_strong: int
    config_class: str
    primary_mixed: bool
    baseline_config: str
    incremental_weak_observable: str
    strong_observable: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def classify(states: Iterable[str]) -> str:
    states = tuple(states)
    nw = states.count("W")
    ns = states.count("S")
    if ns == 0 and nw == 0:
        return "none"
    if ns == 0:
        return "weak_only" if nw == 1 else "multi_weak"
    if nw == 0:
        return "strong_only" if ns == 1 else "multi_strong"
    if ns == 1 and nw == 1 and states.count("0") == 1:
        return "primary_mixed"
    return "mixed_other"


def generate_registry() -> list[ConfigRow]:
    rows: list[ConfigRow] = []
    for state_tuple in product(STATES, repeat=3):
        config = "".join(state_tuple)
        primary = config in PRIMARY_CONTRASTS
        baseline = PRIMARY_CONTRASTS.get(config, "")
        weak_obs = ""
        strong_obs = ""
        if primary:
            weak_obs = OBSERVABLES[state_tuple.index("W")]
            strong_obs = OBSERVABLES[state_tuple.index("S")]
        rows.append(
            ConfigRow(
                config=config,
                M=state_tuple[0],
                L=state_tuple[1],
                G=state_tuple[2],
                n_absent=state_tuple.count("0"),
                n_weak=state_tuple.count("W"),
                n_strong=state_tuple.count("S"),
                config_class=classify(state_tuple),
                primary_mixed=primary,
                baseline_config=baseline,
                incremental_weak_observable=weak_obs,
                strong_observable=strong_obs,
            )
        )
    rows.sort(key=lambda row: row.config)
    assert len(rows) == 27
    assert {r.config for r in rows if r.primary_mixed} == set(PRIMARY_MIXED_CONFIGS)
    return rows


def _require_nonempty(value: Any, label: str) -> str:
    if value is None or not str(value).strip():
        raise ValueError(f"{label} must be non-empty")
    return str(value)


def validate_operator_audit(audit: dict[str, Any]) -> dict[str, Any]:
    if audit.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"schema_version must equal {SCHEMA_VERSION}")
    panel_id = _require_nonempty(audit.get("candidate_panel_id"), "candidate_panel_id")
    operators = audit.get("operators")
    if not isinstance(operators, dict):
        raise ValueError("operators must be an object")

    normalized: dict[str, dict[str, dict[str, Any]]] = {}
    for obs in OBSERVABLES:
        obs_entry = operators.get(obs)
        if not isinstance(obs_entry, dict):
            raise ValueError(f"operators.{obs} must be an object")
        normalized[obs] = {}
        for state in ("W", "S"):
            entry = obs_entry.get(state)
            if not isinstance(entry, dict):
                raise ValueError(f"operators.{obs}.{state} must be an object")
            status = _require_nonempty(entry.get("status"), f"operators.{obs}.{state}.status")
            if status not in ALLOWED_OPERATOR_STATUS:
                raise ValueError(f"invalid operator status for {obs}/{state}: {status}")

            item = dict(entry)
            item["status"] = status
            if status == "QUALIFIED":
                level = _require_nonempty(item.get("information_level"), f"operators.{obs}.{state}.information_level")
                if level != EXPECTED_INFORMATION_LEVEL[state]:
                    raise ValueError(
                        f"{obs}/{state} information_level={level!r}; expected {EXPECTED_INFORMATION_LEVEL[state]!r}"
                    )
                operator_panel = _require_nonempty(item.get("candidate_panel_id"), f"operators.{obs}.{state}.candidate_panel_id")
                _require_nonempty(item.get("candidate_coordinate"), f"operators.{obs}.{state}.candidate_coordinate")
                _require_nonempty(item.get("measurement_coordinate"), f"operators.{obs}.{state}.measurement_coordinate")
                _require_nonempty(item.get("operator_kind"), f"operators.{obs}.{state}.operator_kind")
                if bool(item.get("uses_truth", False)):
                    raise ValueError(f"qualified {obs}/{state} operator must be outcome-blind (uses_truth=false)")
                if bool(item.get("solver_required", False)):
                    raise ValueError(f"qualified {obs}/{state} operator cannot require a solver in BRIDGE-3A")
                item["panel_matches_registry"] = operator_panel == panel_id
            else:
                item["panel_matches_registry"] = False
            normalized[obs][state] = item

    result = dict(audit)
    result["candidate_panel_id"] = panel_id
    result["operators"] = normalized
    return result


def qualify_registry(audit: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    audit = validate_operator_audit(audit)
    rows = generate_registry()
    qualified_rows: list[dict[str, Any]] = []

    for row in rows:
        blockers: list[str] = []
        for obs, state in zip(OBSERVABLES, (row.M, row.L, row.G)):
            if state == "0":
                continue
            operator = audit["operators"][obs][state]
            if operator["status"] != "QUALIFIED":
                blockers.append(f"{obs}{state}:{operator['status']}")
            elif not operator["panel_matches_registry"]:
                blockers.append(f"{obs}{state}:PANEL_MISMATCH")
        qualified_rows.append(
            {
                **row.__dict__,
                "qualification_status": "EXECUTABLE" if not blockers else "BLOCKED",
                "blockers": ";".join(blockers),
            }
        )

    by_config = {r["config"]: r for r in qualified_rows}
    contrasts: list[dict[str, Any]] = []
    for mixed in PRIMARY_MIXED_CONFIGS:
        base = PRIMARY_CONTRASTS[mixed]
        mixed_row = by_config[mixed]
        base_row = by_config[base]
        pair_ready = (
            mixed_row["qualification_status"] == "EXECUTABLE"
            and base_row["qualification_status"] == "EXECUTABLE"
        )
        contrasts.append(
            {
                "contrast_id": f"{mixed}_vs_{base}",
                "mixed_config": mixed,
                "baseline_config": base,
                "strong_observable": mixed_row["strong_observable"],
                "incremental_weak_observable": mixed_row["incremental_weak_observable"],
                "contrast_status": "READY" if pair_ready else "BLOCKED",
                "mixed_blockers": mixed_row["blockers"],
                "baseline_blockers": base_row["blockers"],
            }
        )

    executable = sum(r["qualification_status"] == "EXECUTABLE" for r in qualified_rows)
    primary_ready = sum(r["contrast_status"] == "READY" for r in contrasts)
    manifest = {
        "registry_version": REGISTRY_VERSION,
        "operator_audit_schema": SCHEMA_VERSION,
        "candidate_panel_id": audit["candidate_panel_id"],
        "registry_configs": 27,
        "executable_configs": executable,
        "primary_mixed_contrasts": 6,
        "primary_mixed_contrasts_ready": primary_ready,
        "status": "READY_FOR_MIXED_SCREEN" if primary_ready == 6 else "BLOCKED_OPERATOR_QUALIFICATION",
        "stage": "BRIDGE_C0_C1",
        "scientific_results_generated": False,
        "solver_invoked": False,
        "optimization_invoked": False,
        "fva_invoked": False,
        "sampling_invoked": False,
        "reconstruction_invoked": False,
        "new_weight_generation_invoked": False,
        "gene_flux_concordance_invoked": False,
        "bridge2a_numeric_results_used": False,
    }
    return qualified_rows, contrasts, manifest


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("cannot write empty TSV")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(audit_path: Path, outdir: Path) -> dict[str, Any]:
    audit = json.loads(audit_path.read_text())
    registry, contrasts, manifest = qualify_registry(audit)
    outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(outdir / "BRIDGE3A_CONFIG_REGISTRY.tsv", registry)
    write_tsv(outdir / "BRIDGE3A_PRIMARY_CONTRASTS.tsv", contrasts)
    audit_copy = dict(audit)
    audit_copy["source_audit_sha256"] = sha256_file(audit_path)
    (outdir / "BRIDGE3A_OPERATOR_AUDIT.json").write_text(json.dumps(audit_copy, indent=2, sort_keys=True) + "\n")
    manifest = dict(manifest)
    manifest["operator_audit_sha256"] = sha256_file(audit_path)
    (outdir / "BRIDGE3A_MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--operator-audit", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = run(args.operator_audit.resolve(strict=True), args.outdir.resolve())
    print(json.dumps(manifest, sort_keys=True))
    return 0 if manifest["status"] == "READY_FOR_MIXED_SCREEN" else 2


if __name__ == "__main__":
    raise SystemExit(main())
