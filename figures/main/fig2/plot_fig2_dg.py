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
fig2_dir = OUT
data_path = ROOT / "data/figure_inputs/fig2/figure2_DG_points.csv.xz"


# Source notebook cell 3
# Panel order and labels
panel_order = ["D", "E", "F", "G"]

panel_titles = {
    "D": "Development\n1-Strong-Anchor",
    "E": "Development\n2-Strong-Anchors",
    "F": "Held-out\n1-Strong-Anchor",
    "G": "Held-out\n2-Strong-Anchors",
}

# Column names from the frozen point-table contract
x_col = "x_direction_explained_magnitude_variance"
y_col = "y_mean_g_info"

# Human-facing axis labels
x_label = "Fraction of magnitude variance explained by direction"
y_label = r"Mean information advantage, $g_{\mathrm{info}}$"

# Easy-to-edit visual parameters
figure_width = 7.2
figure_height = 6.4
point_size = 12
point_alpha = 0.45
panel_label_size = 12
title_size = 10
axis_label_size = 10
tick_label_size = 8
rho_text_size = 9

# Shared-axis padding as a fraction of the full data range
x_padding_fraction = 0.04
y_padding_fraction = 0.06

# Output filenames
combined_svg = fig2_dir / "figure2_DG_g_info_scatter.svg"
combined_png = fig2_dir / "figure2_DG_g_info_scatter.png"

# Optional separate panel files. Set True when useful for manual assembly.
save_individual_panels = False

# Source notebook cell 5
df = pd.read_csv(data_path, compression="xz", float_precision="round_trip")

expected_columns = [
    "schema_version",
    "panel",
    "panel_order",
    "anchor_setting",
    "reaction_split",
    "reaction_id",
    "pathway",
    x_col,
    y_col,
    "n_eta2_evaluations",
    "n_g_info_evaluations",
    "n_directional_usefulness_evaluations",
    "n_non_tie_cases",
]

missing_columns = []
for col in expected_columns:
    if col not in df.columns:
        missing_columns.append(col)

if missing_columns:
    raise ValueError(f"Missing expected columns: {missing_columns}")

if df.duplicated(["panel", "reaction_id"]).any():
    duplicate_rows = df.loc[df.duplicated(["panel", "reaction_id"], keep=False), ["panel", "reaction_id"]]
    raise ValueError(f"Duplicate (panel, reaction_id) rows found:\n{duplicate_rows.head(20)}")

if not np.isfinite(df[x_col].to_numpy(dtype=float)).all():
    raise ValueError(f"Non-finite values found in {x_col}")

if not np.isfinite(df[y_col].to_numpy(dtype=float)).all():
    raise ValueError(f"Non-finite values found in {y_col}")

expected_panel_metadata = {
    "D": ("A1", "development"),
    "E": ("A2-L", "development"),
    "F": ("A1", "held_out"),
    "G": ("A2-L", "held_out"),
}

for panel in panel_order:
    panel_df = df.loc[df["panel"] == panel]
    if panel_df.empty:
        raise ValueError(f"Panel {panel} has no rows")

    expected_anchor, expected_split = expected_panel_metadata[panel]
    observed_anchors = sorted(panel_df["anchor_setting"].dropna().astype(str).unique().tolist())
    observed_splits = sorted(panel_df["reaction_split"].dropna().astype(str).unique().tolist())

    if observed_anchors != [expected_anchor]:
        raise ValueError(f"Panel {panel}: expected anchor {expected_anchor}, observed {observed_anchors}")
    if observed_splits != [expected_split]:
        raise ValueError(f"Panel {panel}: expected split {expected_split}, observed {observed_splits}")

print(f"Loaded {len(df):,} plotted reaction rows")
df.head()

# Source notebook cell 18
usefulness_y_col = "y_directional_usefulness_fraction"
usefulness_qc_rows = []

for panel in panel_order:
    panel_df = df.loc[df["panel"] == panel].copy()
    finite_rows = np.isfinite(panel_df[x_col]) & np.isfinite(panel_df[usefulness_y_col])
    plotted_df = panel_df.loc[finite_rows].copy()
    rho_result = spearmanr(plotted_df[x_col], plotted_df[usefulness_y_col])

    usefulness_qc_rows.append({
        "panel": panel,
        "anchor_setting": panel_df["anchor_setting"].iloc[0],
        "split": panel_df["reaction_split"].iloc[0],
        "n_plotted_reactions": len(plotted_df),
        "spearman_rho": float(rho_result.statistic),
        "median_n_eta2_evaluations": float(plotted_df["n_eta2_evaluations"].median()),
        "median_n_directional_usefulness_evaluations": float(plotted_df["n_directional_usefulness_evaluations"].median()),
    })

usefulness_qc_table = pd.DataFrame(usefulness_qc_rows)
expected_usefulness = {
    "D": (2454, 0.8571160663367106),
    "E": (2485, 0.8339440056039364),
    "F": (838, 0.8351435823183291),
    "G": (843, 0.8493657186460156),
}
for _, row in usefulness_qc_table.iterrows():
    expected_n, expected_rho = expected_usefulness[row["panel"]]
    if int(row["n_plotted_reactions"]) != expected_n or not np.isclose(row["spearman_rho"], expected_rho, rtol=0, atol=1e-6):
        raise ValueError(f"A23 usefulness parity mismatch for panel {row['panel']}: n={row['n_plotted_reactions']} rho={row['spearman_rho']}")

usefulness_qc_table[[
    "panel", "anchor_setting", "split", "n_plotted_reactions", "spearman_rho",
    "median_n_eta2_evaluations", "median_n_directional_usefulness_evaluations",
]]

# Source notebook cell 20
usefulness_x_min = float(df[x_col].min())
usefulness_x_max = float(df[x_col].max())
usefulness_x_range = usefulness_x_max - usefulness_x_min
if usefulness_x_range == 0:
    usefulness_x_range = 1.0
usefulness_shared_xlim = (
    usefulness_x_min - x_padding_fraction * usefulness_x_range,
    usefulness_x_max + x_padding_fraction * usefulness_x_range,
)
usefulness_shared_ylim = (0.0, 1.0)
usefulness_y_label = "Benchmark directional-usefulness\nfrequency"
usefulness_svg = fig2_dir / "figure2_DG_usefulness_scatter.svg"
usefulness_png = fig2_dir / "figure2_DG_usefulness_scatter.png"

fig_usefulness, axes_usefulness = plt.subplots(
    1, 4, figsize=(figure_width*2, figure_height/2), sharex=True, sharey=True,
)
axes_usefulness_flat = axes_usefulness.ravel()

for i, panel in enumerate(panel_order):
    ax = axes_usefulness_flat[i]
    panel_df = df.loc[df["panel"] == panel].copy()
    panel_df = panel_df.loc[np.isfinite(panel_df[x_col]) & np.isfinite(panel_df[usefulness_y_col])].copy()
    rho_result = spearmanr(panel_df[x_col], panel_df[usefulness_y_col])

    ax.scatter(panel_df[x_col], panel_df[usefulness_y_col], s=point_size, alpha=point_alpha, linewidths=0)
    ax.set_xlim(usefulness_shared_xlim)
    ax.set_ylim(usefulness_shared_ylim)
    ax.tick_params(labelsize=tick_label_size)
    ax.text(-0.16, 1.07, panel, transform=ax.transAxes, fontsize=panel_label_size*1.2, fontweight="bold", va="top")
    ax.set_title(panel_titles[panel], fontsize=title_size*1.2)
    ax.text(0.04, 0.96, rf"$\rho$ = {rho_result.statistic:.3f}" + "\n" + rf"$n_{{\mathrm{{reactions}}}}$ = {len(panel_df):,}",
            transform=ax.transAxes, fontsize=rho_text_size*1.2, va="top", ha="left")

fig_usefulness.supxlabel(x_label, fontsize=axis_label_size*1.2)
fig_usefulness.supylabel(usefulness_y_label, fontsize=axis_label_size*1.2)
fig_usefulness.tight_layout()
fig_usefulness.savefig(usefulness_svg, bbox_inches="tight")
fig_usefulness.savefig(usefulness_png, dpi=300, bbox_inches="tight")
plt.show()
print("Saved:", usefulness_svg)
print("Saved:", usefulness_png)
