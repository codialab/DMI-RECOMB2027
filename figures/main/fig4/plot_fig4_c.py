#!/usr/bin/env python3
"""Render the current three-method Figure 4C from frozen candidate support."""
from __future__ import annotations

import hashlib
import json
import lzma
import os
import shutil
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(os.environ["RECOMB_OUTPUT_DIR"])
DATA = ROOT / "data/figure_inputs/fig4"
GROUP_ID = "FCG_ed8b1873c8dbc03fe5eb"
METHODS = ["CORDA", "GIMME", "iMAT"]
ARMS = ["A1", "A2-L"]

manifest_path = DATA / "fig4_data_manifest.json"
groups_path = DATA / "fig4_panelC_candidate_groups.parquet.xz"
metrics_path = DATA / "fig4_panelC_candidate_metrics.parquet.xz"
support_path = DATA / "fig4_panelC_candidate_support.parquet.xz"
manifest = json.loads(manifest_path.read_text())
assert manifest["status"] == "PASS"
assert manifest["candidate_groups_exported"] == 100
assert manifest["candidate_ranking"]["utility_used_for_ranking"] is False

def read_parquet_xz(path: Path) -> pd.DataFrame:
    with tempfile.TemporaryDirectory(prefix="fig4c-read-") as directory:
        parquet_path = Path(directory) / "table.parquet"
        with lzma.open(path, "rb", format=lzma.FORMAT_XZ) as source, parquet_path.open("wb") as target:
            shutil.copyfileobj(source, target, length=1024 * 1024)
        return pd.read_parquet(parquet_path, engine="pyarrow")

groups = read_parquet_xz(groups_path)
metrics = read_parquet_xz(metrics_path)
support = read_parquet_xz(support_path)
assert len(groups) == 100 and groups.candidate_group_id.is_unique
selected = groups.loc[groups.candidate_group_id.eq(GROUP_ID)]
assert len(selected) == 1
selection = selected.iloc[0]
assert selection.reaction_id == "FACOAL204"
assert selection.rna_context_key == "training_samples=setx1,setx3"
assert selection.ct2a_mouse == "C1" and selection.gl261_mouse == "G1"
assert int(selection.geometry_rank) == 23

plot_df = support.loc[
    support.candidate_group_id.eq(GROUP_ID)
    & support.reaction_id.eq("FACOAL204")
    & support.method.isin(METHODS)
].copy()
assert len(plot_df) == 3 * 2 * 400
assert set(plot_df.method) == set(METHODS) and set(plot_df.arm) == set(ARMS)
assert not plot_df.duplicated(["method", "arm", "support_id"]).any()
assert np.isfinite(plot_df[["delta_v_B", "weight"]].to_numpy(float)).all()
assert (plot_df.weight >= 0).all()
assert np.isclose(plot_df.groupby(["method", "arm"]).weight.sum().to_numpy(float), 1.0, rtol=0, atol=5e-12).all()
assert np.array_equal(np.where(plot_df.delta_v_B > manifest["truth_atol"], 1,
                               np.where(plot_df.delta_v_B < -manifest["truth_atol"], -1, 0)),
                      plot_df.sign.to_numpy(int))

# Geometry table and support table must refer to the same selected group/method rows.
metric_rows = metrics.loc[metrics.candidate_group_id.eq(GROUP_ID)]
assert len(metric_rows) == 4 and set(metric_rows.algorithm) == {"CORDA", "GIMME", "iMAT", "RIPTiDe"}
assert metric_rows.reaction_id.eq("FACOAL204").all()
for method in METHODS:
    metric_row = metric_rows.loc[metric_rows.algorithm.eq(method)].iloc[0]
    for arm, prefix in [("A1", "A1"), ("A2-L", "A2")]:
        block = plot_df.loc[plot_df.method.eq(method) & plot_df.arm.eq(arm)]
        signs = block.sign.to_numpy(int)
        weights = block.weight.to_numpy(float)
        mass = np.array([weights[signs < 0].sum(), weights[signs == 0].sum(), weights[signs > 0].sum()])
        positive = mass[mass > 0]
        entropy = 0.0 if len(positive) <= 1 else float(-(positive * np.log(positive)).sum() / np.log(3.0))
        assert np.isclose(entropy, float(metric_row[f"{prefix}_H_dir"]), rtol=0,
                          atol=float(manifest["support_entropy_atol"]))

OUT.mkdir(parents=True, exist_ok=True)
fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.35), sharey=True, constrained_layout=True)
colors = {"A1": "#3568a8", "A2-L": "#d06a3b"}
labels = {"A1": "1-Strong-Anchor (HEX1)", "A2-L": "2-Strong-Anchors (HEX1 + LDH_L)"}
for ax, method in zip(axes, METHODS):
    method_df = plot_df.loc[plot_df.method.eq(method)]
    all_delta = method_df.delta_v_B.to_numpy(float)
    lo, hi = float(all_delta.min()), float(all_delta.max())
    if hi <= lo:
        pad = max(abs(lo) * .05, 1e-12)
        lo, hi = lo - pad, hi + pad
    else:
        pad = .03 * (hi - lo)
        lo, hi = lo - pad, hi + pad
    bins = np.linspace(lo, hi, 31)
    for arm in ARMS:
        block = method_df.loc[method_df.arm.eq(arm)]
        values = block.delta_v_B.to_numpy(float)
        weights = block.weight.to_numpy(float)
        density, edges = np.histogram(values, bins=bins, weights=weights, density=True)
        assert np.isclose(float(np.sum(density * np.diff(edges))), 1.0, rtol=0, atol=5e-12)
        ax.stairs(density, edges, linewidth=2, color=colors[arm], label=labels[arm])
    ax.axvline(0, color="black", linestyle="--", linewidth=1)
    ax.set_title(method, fontsize=12, fontweight="bold")
    ax.set_xlabel(r"$\Delta v_B$ (native GEM flux units)")
    ax.grid(axis="y", color="#dddddd", linewidth=.6)
    ax.spines[["top", "right"]].set_visible(False)
axes[0].set_ylabel("Weighted probability density")
axes[-1].legend(frameon=False, fontsize=8, loc="upper right")
fig.suptitle("Figure 4C | Anchor-dependent $P(\\Delta v_B)$ for FACOAL204\nRNA context setx1 + setx3; CT2A C1 vs GL261 G1", fontsize=12)
fig.savefig(OUT / "fig4C_FACOAL204_three_methods.svg", bbox_inches="tight")
fig.savefig(OUT / "fig4C_FACOAL204_three_methods.png", dpi=300, bbox_inches="tight")

provenance = {
    "panel": "Figure 4C",
    "selection": {
        "candidate_group_id": GROUP_ID,
        "reaction_id": "FACOAL204",
        "rna_context_key": "training_samples=setx1,setx3",
        "ct2a_mouse": "C1",
        "gl261_mouse": "G1",
        "truth_selection": str(selection.truth_selection),
        "geometry_rank_in_exported_top_100": int(selection.geometry_rank),
        "range_A2_minus_A1_H_dir_shift": float(selection.range_A2_minus_A1_H_dir_shift),
        "selected_methods": METHODS,
        "excluded_from_main_panel": {"RIPTiDe": "Current manuscript assigns RIPTiDe to S7; older four-method notebook composition retained as provenance-only."},
        "selection_basis": "Existing manuscript-selected group, read from the frozen candidate group table. It is not reselected by this plotter. The frozen top-100 screen ranked by the manifest rule; the technical audit reports this group has the largest cross-method A2-L-minus-A1 entropy-shift range within the exported 100.",
        "top_100_ranking_rule": manifest["candidate_ranking"]["rule"],
    },
    "population": "400 weighted support states per selected method and anchor; A1 and A2-L plotted with frozen product weights.",
    "source_files": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [manifest_path, groups_path, metrics_path, support_path]},
    "plot_script": "figures/main/fig4/plot_fig4_c.py",
    "plot_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "binning": "30 shared bins within each method across A1 and A2-L, matching the original distribution construction.",
}
(OUT / "fig4C_selection_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
