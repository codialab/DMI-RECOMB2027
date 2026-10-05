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
DATA = ROOT / "data/figure_inputs/fig4"
PANEL = OUT
manifest = json.loads((DATA / "fig4_data_manifest.json").read_text())



# Source notebook cell 11
# Figure 4C — representative weighted P(Δv_B) distributions
# Selected outcome-blind candidate:
# FACOAL204 | RNA context setx1,setx3 | C1 vs G1
# candidate_group_id = FCG_ed8b1873c8dbc03fe5eb

import lzma
import shutil
import tempfile
from pathlib import Path

candidate_group_id = 'FCG_ed8b1873c8dbc03fe5eb'
reaction_id = 'FACOAL204'
methods_to_plot = ['CORDA', 'GIMME', 'iMAT', 'RIPTiDe']

# ------------------------------------------------------------------
# Load the XZ-wrapped Parquet support table
# Skip this block if panelC_support is already loaded in the notebook.
# ------------------------------------------------------------------

support_path = DATA / 'fig4_panelC_candidate_support.parquet.xz'

with tempfile.TemporaryDirectory(prefix='fig4C-support-') as tmpdir:
    parquet_path = Path(tmpdir) / 'support.parquet'

    with lzma.open(support_path, 'rb') as src, parquet_path.open('wb') as dst:
        shutil.copyfileobj(src, dst, length=1024 * 1024)

    panelC_support = pd.read_parquet(parquet_path)

# ------------------------------------------------------------------
# Select the frozen Figure 4C candidate
# ------------------------------------------------------------------

plot_df = panelC_support.loc[
    panelC_support['candidate_group_id'].eq(candidate_group_id)
    & panelC_support['reaction_id'].eq(reaction_id)
    & panelC_support['method'].isin(methods_to_plot)
].copy()

# Expected: 3 methods × 2 anchor settings × 400 support points
assert len(plot_df) == 4 * 2 * 400, len(plot_df)

counts = (
    plot_df
    .groupby(['method', 'arm'])
    .size()
)

assert set(counts.values) == {400}
assert set(plot_df['method']) == set(methods_to_plot)
assert set(plot_df['arm']) == {'A1', 'A2-L'}

# ------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------

fig, axes = plt.subplots(
    2, 2,
    figsize=(8, 8),
    constrained_layout=True,
    sharey=False
)

n_bins = 30

for ii, method in enumerate(methods_to_plot):
    iii, jjj = ii // 2, ii % 2
    ax = axes[iii, jjj]

    method_df = plot_df.loc[plot_df['method'].eq(method)]

    # Common binning for A1 and A2-L within this method
    all_delta = method_df['delta_v_B'].to_numpy(dtype=float)

    lo = np.nanmin(all_delta)
    hi = np.nanmax(all_delta)

    if not np.isfinite(lo) or not np.isfinite(hi):
        raise RuntimeError(f'Non-finite Δv_B range for {method}')

    if hi <= lo:
        pad = max(abs(lo) * 0.05, 1e-12)
        lo -= pad
        hi += pad
    else:
        pad = 0.03 * (hi - lo)
        lo -= pad
        hi += pad

    bins = np.linspace(lo, hi, n_bins + 1)

    hdir_values = {}

    for arm, color, label in [
        ('A1', '#3568a8', '1-Strong-Anchor\n(HEX1)'),
        ('A2-L', '#d06a3b', '2-Strong-Anchors\n(HEX1 + LDH_L)')
    ]:

        arm_df = method_df.loc[method_df['arm'].eq(arm)]

        delta = arm_df['delta_v_B'].to_numpy(dtype=float)
        weight = arm_df['weight'].to_numpy(dtype=float)

        # Normalize defensively
        weight = weight / weight.sum()

        # Weighted density
        density, edges = np.histogram(
            delta,
            bins=bins,
            weights=weight,
            density=True
        )

        ax.stairs(
            density,
            edges,
            linewidth=2.0,
            color=color,
            label=label
        )

        # Reconstruct direction entropy directly from frozen sign labels
        signs = arm_df['sign'].to_numpy(dtype=int)

        sign_mass = np.array([
            weight[signs < 0].sum(),
            weight[signs == 0].sum(),
            weight[signs > 0].sum()
        ])

        nonzero_mass = sign_mass[sign_mass > 0]

        if len(nonzero_mass) <= 1:
            hdir = 0.0
        else:
            hdir = -np.sum(
                nonzero_mass * np.log(nonzero_mass)
            ) / np.log(3.0)

        hdir_values[arm] = hdir

    # Δv_B = 0 reference
    ax.axvline(
        0,
        color='black',
        linestyle='--',
        linewidth=1.1
    )

    ax.set_title(method, fontsize=12, fontweight="bold")
    # if iii == 1:
    #     ax.set_xlabel(r'$\Delta v_B$')
    if jjj == 0:
        ax.set_ylabel('Probability density')

    # Compact geometry annotation
    ax.text(
        0.96,
        0.96,
        rf'$H_{{dir}}$: '
        f'{hdir_values["A1"]:.3f} (1-Anchor)\n'
        r' $\rightarrow$ '
        f'{hdir_values["A2-L"]:.3f} (2-Anchors)',
        transform=ax.transAxes,
        va='top',
        ha='right',
        fontsize=10
    )

# One legend only
axes[1,1].legend(
    frameon=False,
    fontsize=9,
    loc='lower right'
)

axes[1,1].set_xlim(-0.5e-12, 1.5e-12)

fig.suptitle(
    r'Anchor-dependent reshaping of $P(\Delta v_B)$ of reaction FACOAL204'
    '\n'
    r'RNA context | CT2A: GSM4577664, GSM4577666'
    '\n'
    r'GL261: GSM4577667, GSM4577669',
    y=1.08,
    fontsize=12,
)
fig.text(
    0.0, 1.03,
    r'C',
    ha='left',
    fontsize=14,
    fontweight='bold'
)
fig.text(
    0.5, 0.00,
    r'$\Delta v_B$ [native GEM flux units]',
    ha='center',
    fontsize=12
)
# fig.tight_layout()

fig.savefig(
    PANEL / 'fig4C_FACOAL204_distributions.svg',
    bbox_inches='tight'
)

fig.savefig(
    PANEL / 'fig4C_FACOAL204_distributions.pdf',
    bbox_inches='tight'
)

plt.show()
