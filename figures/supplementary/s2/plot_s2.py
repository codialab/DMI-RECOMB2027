#!/usr/bin/env python3
"""Render audited Supplementary Figure S2 from included frozen paired shifts.

Notebook provenance: figures/supp_fig2/supp_fig2.ipynb, cells 1–2;
see manifests/provenance/NOTEBOOK_EXTRACTIONS.json.
"""
from pathlib import Path
import json
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = ROOT / 'data/figure_inputs/s2'
OUTPUT_DIR = Path(os.environ['RECOMB_OUTPUT_DIR'])
summary = json.loads((DATA_DIR / 'supp_fig2_table_summary.json').read_text())
assert summary['status'] == 'PASS'
METHODS = ['iMAT', 'GIMME', 'CORDA', 'RIPTiDe']
METRICS = [('H_dir', 'Directional entropy shift (A2 − A1)'),
           ('dominant_direction_mass', 'Dominant-direction mass shift (A2 − A1)')]
COLORS = dict(zip(METHODS, plt.get_cmap('tab10').colors[:4]))
plt.rcParams.update({'font.size': 10, 'axes.titlesize': 11, 'axes.labelsize': 10, 'legend.fontsize': 8})
fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.25), constrained_layout=True)
for ax, (metric, xlabel) in zip(axes, METRICS):
    info = summary['metrics'][metric]
    table = pd.read_csv(ROOT / info['output'], sep='\t', compression='gzip')
    assert len(table) == info['output_rows'] == 1_671_600
    for method in METHODS:
        values = table.loc[table.algorithm.eq(method), 'delta_A2_minus_A1'].to_numpy(float)
        assert len(values) == info['panel_n_by_method'][method] == 417_900
        assert np.isfinite(values).all()
        x = np.sort(values); y = np.arange(1, len(x) + 1, dtype=float) / len(x)
        ax.step(x, y, where='post', color=COLORS[method], linewidth=1.45, label=f'{method} (n={len(x):,})')
    ax.axvline(0, color='black', linewidth=1.35, linestyle=(0, (4, 3)), zorder=10)
    ax.set_title('Directional entropy' if metric == 'H_dir' else 'Dominant-direction mass')
    ax.set_xlabel(xlabel + ' (dimensionless)'); ax.set_ylabel('Cumulative fraction')
    ax.set_ylim(0, 1); ax.grid(axis='y', color='#d9d9d9', linewidth=.6, alpha=.8)
    ax.legend(frameon=False, loc='lower right')
fig.suptitle('Paired shifts in residual flux-distribution geometry', fontsize=13, y=1.04)
fig.text(.5, -.025, 'A1 = 1-Strong-Anchor (glucose uptake); A2 = 2-Strong-Anchors (glucose uptake + LDH_L)', ha='center', fontsize=9)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
fig.savefig(OUTPUT_DIR / 'supp_fig2.svg', bbox_inches='tight')
fig.savefig(OUTPUT_DIR / 'supp_fig2.png', dpi=300, bbox_inches='tight')
