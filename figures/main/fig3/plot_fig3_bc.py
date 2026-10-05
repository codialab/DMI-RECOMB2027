#!/usr/bin/env python3
"""Selected source plotting cells; see manifests/provenance/NOTEBOOK_EXTRACTIONS.json."""
from pathlib import Path
import os
import hashlib, io, json, lzma, shutil, tempfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
ROOT = Path(__file__).resolve().parents[3]
OUT = Path(os.environ["RECOMB_OUTPUT_DIR"])
OUT.mkdir(parents=True, exist_ok=True)
DATA = ROOT / "data/figure_inputs/fig3"
with lzma.open(DATA / "fig3_gain_reaction_evaluation.parquet.xz", "rb") as f:
    reaction_eval = pd.read_parquet(io.BytesIO(f.read()))
anchor_order = ["A1", "A2-L"]  # Historical input schema; display labels use A2.
anchor_colors = {"A1": "#2878B5", "A2-L": "#D9534F"}
display_labels = {"A1": "A1 (HEX1)", "A2-L": "A2 (HEX1 + LDH_L)"}
def export_plot(fig, stem):
    fig.savefig(OUT / (stem + ".svg"), bbox_inches="tight")
    fig.savefig(OUT / (stem + ".png"), dpi=300, bbox_inches="tight")


# Source notebook cell 20
fig, axex = plt.subplots(figsize=(14, 5), nrows=1, ncols=2, sharey=True)


thresholds = np.logspace(-15, 0, 250)

ax = axex[0]
for anchor in anchor_order:
    v = reaction_eval.loc[
        reaction_eval["anchor_setting"].eq(anchor),
        "correct_gain_mean"
    ].to_numpy(float)

    frac = np.array([
        np.mean(v > t)
        for t in thresholds
    ])

    ax.plot(
        thresholds,
        frac,
        linewidth=1.8,
        color=anchor_colors[anchor],
        label=display_labels[anchor],
    )

ax.set_xscale("log")
ax.set_xlabel(r"Threshold $\tau$")
ax.set_ylabel(r"Fraction with $g_{\mathrm{correct}} > \tau$", fontsize=12)
ax.set_title("B Benefit from a correct directional cue", x=0.3, y=1.02, fontsize=14, fontweight="bold")

ax.legend(frameon=False)
ax.spines[["top", "right"]].set_visible(False)

# export_plot(fig, "fig3_panelB_correct_gain_exceedance")





thresholds = np.logspace(-15, 0, 250)

ax = axex[1]
for anchor in anchor_order:
    v = reaction_eval.loc[
        reaction_eval["anchor_setting"].eq(anchor),
        "directional_advantage"
    ].to_numpy(float)

    frac = np.array([
        np.mean(v > t)
        for t in thresholds
    ])

    ax.plot(
        thresholds,
        frac,
        linewidth=1.8,
        color=anchor_colors[anchor],
        label=display_labels[anchor],
    )

ax.set_xscale("log")
ax.set_xlabel(r"Threshold $\tau$")
ax.set_ylabel(r"Fraction with $(\Delta g_{\mathrm{correct}} - \Delta g_{\mathrm{wrong}}) > \tau$", fontsize=12)
ax.set_title("C Advantage of the correct direction", x=0.3, y=1.02, fontsize=14, fontweight="bold")

ax.legend(frameon=False)
ax.spines[["top", "right"]].set_visible(False)

export_plot(fig, "fig3_panelBC")
