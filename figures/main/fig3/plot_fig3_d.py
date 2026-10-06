#!/usr/bin/env python3
"""Render the pair-weighted Figure 3D endpoint panel from frozen result tables."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(os.environ["RECOMB_OUTPUT_DIR"])
SUMMARY = ROOT / "tables/manuscript/fig3_pair_gain_summary.tsv"
SEPARATION = ROOT / "tables/manuscript/fig3_full_endpoint_separation.tsv"

summary = pd.read_csv(SUMMARY, sep="\t", float_precision="round_trip")
separation = pd.read_csv(SEPARATION, sep="\t", float_precision="round_trip")

rows = []
for setting, label in [("A1", "A1 (HEX1)"), ("A2-L", "A2 (HEX1 + LDH_L)")]:
    block = summary.loc[summary.anchor_setting.eq(setting)]
    correct = block.loc[block.gain_measure.eq("correct")]
    wrong = block.loc[block.gain_measure.eq("wrong")]
    assert len(correct) == len(wrong) == 1
    correct_mean = float(correct.iloc[0].mean_gain)
    wrong_mean = float(wrong.iloc[0].mean_gain)
    n_correct, n_wrong = int(correct.iloc[0].n_reaction_evaluations), int(wrong.iloc[0].n_reaction_evaluations)
    assert n_correct == n_wrong == 660_237
    derived = correct_mean - wrong_mean
    endpoint = separation.loc[separation.anchor_setting.eq(setting)]
    assert len(endpoint) == 1
    assert np.isclose(derived, float(endpoint.iloc[0]["mean"]), rtol=0, atol=5e-12)
    rows.append({
        "anchor_setting": setting,
        "display_label": label,
        "n_reaction_evaluations": n_correct,
        "mean_wrong_cue_gain": wrong_mean,
        "mean_correct_cue_gain": correct_mean,
        "correct_minus_wrong_separation": derived,
    })

panel = pd.DataFrame(rows)
assert np.isfinite(panel[["mean_wrong_cue_gain", "mean_correct_cue_gain", "correct_minus_wrong_separation"]].to_numpy()).all()
OUT.mkdir(parents=True, exist_ok=True)
panel.to_csv(OUT / "fig3D_pair_weighted_summary.tsv", sep="\t", index=False, float_format="%.17g")

colors = {"wrong": "#777777", "correct": "#2878B5"}
fig, ax = plt.subplots(figsize=(9.2, 3.7), constrained_layout=True)
y_positions = {"A1": 1.2, "A2-L": .4}
for row in rows:
    y = y_positions[row["anchor_setting"]]
    wrong, correct = row["mean_wrong_cue_gain"], row["mean_correct_cue_gain"]
    ax.plot([wrong, correct], [y, y], color="#999999", linewidth=1.8, zorder=1)
    ax.scatter(wrong, y, s=65, color=colors["wrong"], edgecolor="white", linewidth=.7,
               zorder=3, label="Wrong-cue mean" if row["anchor_setting"] == "A1" else None)
    ax.scatter(correct, y, s=65, color=colors["correct"], edgecolor="white", linewidth=.7,
               zorder=3, label="Correct-cue mean" if row["anchor_setting"] == "A1" else None)
    ax.annotate(f"{wrong:.6g}", (wrong, y), xytext=(0, 9), textcoords="offset points",
                ha="center", va="bottom", fontsize=8, color=colors["wrong"])
    ax.annotate(f"{correct:.6g}", (correct, y), xytext=(0, 9), textcoords="offset points",
                ha="center", va="bottom", fontsize=8, color=colors["correct"])
    ax.text(.99, y + .22, f"Correct − wrong = {row['correct_minus_wrong_separation']:.6g}",
            transform=ax.get_yaxis_transform(), ha="right", va="bottom", fontsize=8.5)
ax.axvline(0, color="#222222", linewidth=.8, linestyle=(0, (3, 3)), zorder=0)
ax.set_yticks([y_positions["A1"], y_positions["A2-L"]], [r["display_label"] for r in rows])
ax.tick_params(axis="y", length=0)
ax.set_ylim(0, 1.6)
ax.set_xlabel("Mean cue gain (native GEM flux units)")
fig.suptitle("Figure 3D | Pair-weighted directional-cue gains\nEqual-weight means across eligible reaction–evaluation pairs (n = 660,237 per cue and anchor)", fontsize=12)
ax.grid(axis="x", color="#dddddd", linewidth=.6)
ax.spines[["top", "right", "left"]].set_visible(False)
ax.legend(frameon=True, facecolor="white", edgecolor="white", framealpha=.95,
          loc="center", bbox_to_anchor=(.5, .5), ncol=2)
fig.savefig(OUT / "fig3D_pair_weighted.svg", bbox_inches="tight")
fig.savefig(OUT / "fig3D_pair_weighted.png", dpi=300, bbox_inches="tight")

provenance = {
    "panel": "Figure 3D",
    "source_tables": {
        "fig3_pair_gain_summary.tsv": hashlib.sha256(SUMMARY.read_bytes()).hexdigest(),
        "fig3_full_endpoint_separation.tsv": hashlib.sha256(SEPARATION.read_bytes()).hexdigest(),
    },
    "plot_script": "figures/main/fig3/plot_fig3_d.py",
    "plot_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "population": "Eligible reaction–evaluation pairs, cases averaged within pairs; each pair receives equal weight.",
    "endpoints": panel.to_dict(orient="records"),
    "units": "native GEM flux units",
}
(OUT / "fig3D_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
