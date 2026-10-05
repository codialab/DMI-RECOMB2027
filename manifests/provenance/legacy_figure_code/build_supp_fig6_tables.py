#!/usr/bin/env python3
"""Build compact S6 tables from manifest-verified paired Figure 4 geometry."""
from __future__ import annotations

import hashlib
import json
import lzma
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FIG4 = ROOT / "figures/fig4/data"
OUT = ROOT / "figures/supp_fig6"
DATA = OUT / "data"
MANIFEST = FIG4 / "fig4_data_manifest.json"
METHODS = ("iMAT", "GIMME", "CORDA", "RIPTiDe")
KEYS = ["evaluation_id", "reaction_id"]
METRICS = {
    "eta2": {"column": "direction_explained_magnitude_variance", "label": "Sign–magnitude η²", "lower": 0.0, "upper": 1.0},
    "non_tie_coverage": {"column": "non_tie_coverage", "label": "Non-tie coverage (1 − p_tie)", "lower": 0.0, "upper": 1.0},
    "supported_sign_states": {"column": "supported_sign_state_count", "label": "Supported sign-state count", "lower": 1, "upper": 3},
}
TOL = 5e-12
ROWS = 1_671_600


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_parquet_xz(path: Path, columns: list[str]) -> pd.DataFrame:
    with tempfile.TemporaryDirectory(prefix="supp-s6-") as tmp:
        parquet = Path(tmp) / "geometry.parquet"
        with lzma.open(path, "rb", format=lzma.FORMAT_XZ) as src, parquet.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=1024 * 1024)
        return pd.read_parquet(parquet, columns=columns)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    source_manifest = json.loads(MANIFEST.read_text())
    if (source_manifest.get("status") != "PASS"
            or source_manifest.get("row_counts", {}).get("geometry_paired") != ROWS
            or source_manifest.get("unique_key_counts", {}).get("geometry_evaluation_reaction") != ROWS
            or source_manifest.get("duplicate_key_checks", {}).get("geometry_evaluation_reaction") != "PASS"
            or source_manifest.get("a1_a2_matching", {}).get("geometry_rows_matched") != ROWS):
        raise RuntimeError("Figure 4 geometry manifest is not the expected complete paired population")

    columns = KEYS + ["rna_context_key", "ct2a_mouse", "gl261_mouse"]
    for spec in METRICS.values():
        col = spec["column"]
        columns.extend([f"A1_{col}", f"A2_{col}"])
    parts = []
    source_hashes = {}
    per_method = {}
    for method in METHODS:
        rel = Path("figures/fig4/data/fig4_geometry_paired") / f"algorithm={method}" / "part.parquet.xz"
        path = ROOT / rel
        if path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(ROOT):
            raise RuntimeError(f"Missing, symlinked, or external geometry source: {rel}")
        digest = sha256(path)
        if digest != source_manifest.get("generated_sha256", {}).get(rel.as_posix()):
            raise RuntimeError(f"Figure 4 manifest checksum mismatch: {rel}")
        source_hashes[rel.as_posix()] = digest
        frame = read_parquet_xz(path, columns)
        if len(frame) != 417_900:
            raise RuntimeError(f"Unexpected method partition rows: {method}")
        frame.insert(2, "algorithm", method)
        if frame[KEYS].isna().any().any() or frame.duplicated(KEYS).any():
            raise RuntimeError(f"Missing or duplicate evaluation/reaction keys for {method}")
        per_method[method] = int(len(frame))
        parts.append(frame)

    paired = pd.concat(parts, ignore_index=True)
    if len(paired) != ROWS or paired.duplicated(["algorithm", *KEYS]).any():
        raise RuntimeError("Unexpected combined geometry population or duplicate keys")
    if paired[KEYS + ["algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse"]].isna().any().any():
        raise RuntimeError("Missing pairing or context metadata")
    if paired.evaluation_id.nunique() != 400 or paired.reaction_id.nunique() != 4_179:
        raise RuntimeError("Unexpected evaluation/reaction population")
    eval_metadata = paired.groupby("evaluation_id")[["algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse"]].nunique()
    if not eval_metadata.eq(1).all().all():
        raise RuntimeError("Evaluation metadata are not consistent within evaluation_id")

    out = paired[KEYS + ["algorithm", "rna_context_key", "ct2a_mouse", "gl261_mouse"]].copy()
    metric_summary = {}
    for key, spec in METRICS.items():
        col = spec["column"]
        a1 = pd.to_numeric(paired[f"A1_{col}"], errors="coerce").to_numpy(float)
        a2 = pd.to_numeric(paired[f"A2_{col}"], errors="coerce").to_numpy(float)
        f1, f2 = np.isfinite(a1), np.isfinite(a2)
        finite_pair = f1 & f2
        for arm, values, finite in (("A1", a1, f1), ("A2", a2, f2)):
            if np.any(values[finite] < spec["lower"] - TOL) or np.any(values[finite] > spec["upper"] + TOL):
                raise RuntimeError(f"{key} {arm} value is outside its valid range")
        if key == "supported_sign_states":
            for arm, values, finite in (("A1", a1, f1), ("A2", a2, f2)):
                if np.any(np.abs(values[finite] - np.rint(values[finite])) > TOL):
                    raise RuntimeError(f"{arm} supported sign-state count is not integer-valued")
        delta = a2 - a1
        out[f"A1_{key}"] = a1
        out[f"A2_{key}"] = a2
        out[f"delta_A2_minus_A1_{key}"] = delta
        metric_summary[key] = {
            "definition": {
                "eta2": "weighted fraction of variance in |Δv_B| explained by sign category {-1, 0, +1}",
                "non_tie_coverage": "1 - p_tie, where p_tie is the frozen probability mass with |Δv_B| within the production tie tolerance",
                "supported_sign_states": "number of {-1, 0, +1} sign categories with positive weighted probability mass",
            }[key],
            "valid_range": [spec["lower"], spec["upper"]],
            "input_rows": ROWS,
            "finite_A1": int(f1.sum()),
            "finite_A2": int(f2.sum()),
            "finite_paired": int(finite_pair.sum()),
            "excluded_from_paired_comparison": int((~finite_pair).sum()),
            "a1_range": [float(np.min(a1[f1])), float(np.max(a1[f1]))] if f1.any() else None,
            "a2_range": [float(np.min(a2[f2])), float(np.max(a2[f2]))] if f2.any() else None,
            "paired_delta_range": [float(np.min(delta[finite_pair])), float(np.max(delta[finite_pair]))] if finite_pair.any() else None,
            "per_method_finite_paired": {m: int(np.sum((paired.algorithm.to_numpy() == m) & finite_pair)) for m in METHODS},
        }

    state_a1 = out["A1_supported_sign_states"]
    state_a2 = out["A2_supported_sign_states"]
    state_finite = np.isfinite(state_a1) & np.isfinite(state_a2)
    if set(np.unique(np.concatenate([state_a1[np.isfinite(state_a1)], state_a2[np.isfinite(state_a2)]]))) - {1.0, 2.0, 3.0}:
        raise RuntimeError("Unexpected supported sign-state category")
    transitions = pd.crosstab(state_a1[state_finite].astype(int), state_a2[state_finite].astype(int), normalize="all")
    if transitions.size != 9:
        transitions = transitions.reindex(index=[1, 2, 3], columns=[1, 2, 3], fill_value=0)
    transition_counts = pd.crosstab(state_a1[state_finite].astype(int), state_a2[state_finite].astype(int)).reindex(index=[1, 2, 3], columns=[1, 2, 3], fill_value=0)
    transition_table = transition_counts.stack().rename("count").reset_index()
    transition_table.columns = ["A1_supported_sign_states", "A2_supported_sign_states", "count"]
    transition_table["fraction_of_paired_rows"] = transition_table["count"] / int(state_finite.sum())
    transition_table.to_csv(DATA / "supp_fig6_sign_state_transitions.tsv", sep="\t", index=False)
    out.to_csv(DATA / "supp_fig6_geometry.tsv.gz", sep="\t", index=False, compression="gzip", float_format="%.17g")

    summary = {
        "status": "PASS",
        "source_manifest": MANIFEST.relative_to(ROOT).as_posix(),
        "source_manifest_sha256": sha256(MANIFEST),
        "source_manifest_status": source_manifest["status"],
        "source_artifacts": source_hashes,
        "anchor_definitions": {"A1": "1-Strong-Anchor (HEX1)", "A2": "2-Strong-Anchors (HEX1 + LDH_L)"},
        "pairing_keys": KEYS,
        "methods": list(METHODS),
        "population": {"rows": ROWS, "evaluations": 400, "reactions": 4_179, "per_method_rows": per_method},
        "metrics": metric_summary,
        "sign_state_transition": {"paired_finite_rows": int(state_finite.sum()), "categories": [1, 2, 3], "output": "figures/supp_fig6/data/supp_fig6_sign_state_transitions.tsv", "aggregation": "pooled paired reaction × evaluation rows; fractions use paired finite rows"},
        "excluded_descriptors": {"joint_ess": "overlaps S1 effective-support diagnostics", "width80": "not selected; contract says not to revive solely to fill the figure", "directional_entropy3": "covered by S2 and S3", "dominant_sign_mass": "covered by S2 and S3"},
        "outputs": {},
    }
    for name in ("supp_fig6_geometry.tsv.gz", "supp_fig6_sign_state_transitions.tsv"):
        path = DATA / name
        summary["outputs"][name] = {"rows": int(ROWS if name.endswith("geometry.tsv.gz") else len(transition_table)), "sha256": sha256(path)}
    (DATA / "supp_fig6_build_manifest.json").write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(f"S6 tables validated: {ROWS:,} geometry pairs; eta2 {metric_summary['eta2']['finite_paired']:,} finite pairs; state transitions {int(state_finite.sum()):,} pairs")


if __name__ == "__main__":
    main()
