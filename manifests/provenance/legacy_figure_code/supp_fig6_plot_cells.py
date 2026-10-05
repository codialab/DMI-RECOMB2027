#!/usr/bin/env python3
"""Selected source plotting cells; see manifests/provenance/NOTEBOOK_EXTRACTIONS.json."""
# Provenance only. Do not execute as a current manuscript figure.


# Source notebook cell 1
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ROOT = Path.cwd().resolve()
while ROOT != ROOT.parent and not (ROOT / 'figures/supp_fig6/data/supp_fig6_build_manifest.json').is_file():
    ROOT = ROOT.parent
if not (ROOT / 'figures/supp_fig6/data/supp_fig6_build_manifest.json').is_file():
    raise FileNotFoundError('Run this notebook from the repository or one of its subdirectories.')
DATA = ROOT / 'figures/supp_fig6/data'
OUTPUTS = ROOT / 'figures/supp_fig6/outputs'
OUTPUTS.mkdir(parents=True, exist_ok=True)
manifest = json.loads((DATA / 'supp_fig6_build_manifest.json').read_text())
assert manifest['status'] == 'PASS'
table_path = DATA / 'supp_fig6_geometry.tsv.gz'
assert hashlib.sha256(table_path.read_bytes()).hexdigest() == manifest['outputs']['supp_fig6_geometry.tsv.gz']['sha256']
paired = pd.read_csv(table_path, sep='\t')
assert len(paired) == manifest['population']['rows'] == 1_671_600
assert not paired.duplicated(['algorithm', 'evaluation_id', 'reaction_id']).any()
print(f"Loaded {len(paired):,} paired reaction–evaluation observations across {paired.algorithm.nunique()} methods.")

# Source notebook cell 2
methods = ['iMAT', 'GIMME', 'CORDA', 'RIPTiDe']
colors = {'iMAT': '#4477AA', 'GIMME': '#EE6677', 'CORDA': '#228833', 'RIPTiDe': '#AA3377'}
panel_specs = [
    ('eta2', r'Sign–magnitude coupling ($\eta^2$)', (-0.02, 1.02)),
    ('non_tie_coverage', 'Non-tie coverage (1 − p_tie)', (-1.02, 1.02)),
]
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.5), constrained_layout=True, gridspec_kw={'width_ratios': [1.15, 1.15, 0.95]})
for ax, (metric, title, xlim) in zip(axes[:2], panel_specs):
    col = f'delta_A2_minus_A1_{metric}'
    for method in methods:
        values = paired.loc[paired.algorithm.eq(method), col].to_numpy(float)
        values = np.sort(values[np.isfinite(values)])
        y = np.arange(1, len(values) + 1) / len(values)
        ax.plot(values, y, lw=1.7, color=colors[method], label=f'{method} (n={len(values):,})')
    ax.axvline(0, color='#444444', lw=1, ls='--')
    ax.set_xlim(*xlim); ax.set_ylim(0, 1)
    ax.set_xlabel('Paired change (A2 − A1)')
    ax.set_ylabel('Cumulative fraction')
    ax.set_title(title)
    ax.grid(axis='y', color='#dddddd', lw=0.6)
    ax.legend(frameon=False, fontsize=7, loc='lower right')

trans = pd.read_csv(DATA / 'supp_fig6_sign_state_transitions.tsv', sep='\t')
mat = trans.pivot(index='A1_supported_sign_states', columns='A2_supported_sign_states', values='fraction_of_paired_rows').reindex(index=[1,2,3], columns=[1,2,3]).fillna(0)
ax = axes[2]
im = ax.imshow(mat.to_numpy(), origin='lower', cmap='Blues', norm=Normalize(vmin=0, vmax=max(mat.to_numpy().max(), 1e-12)), aspect='equal')
for i in range(3):
    for j in range(3):
        frac = mat.iloc[i, j]
        count = int(trans.loc[trans.A1_supported_sign_states.eq(i+1) & trans.A2_supported_sign_states.eq(j+1), 'count'].sum())
        ax.text(j, i, f'{frac:.1%}\n({count:,})', ha='center', va='center', fontsize=7, color='white' if frac > mat.to_numpy().max()*0.55 else '#222222')
ax.set_xticks([0,1,2], ['1','2','3']); ax.set_yticks([0,1,2], ['1','2','3'])
ax.set_xlabel('A2 supported sign states'); ax.set_ylabel('A1 supported sign states')
ax.set_title(f'Sign-state support transitions\n(n={int(trans["count"].sum()):,})')
fig.colorbar(im, ax=ax, shrink=0.78, pad=0.04, label='Fraction of paired rows')
fig.suptitle('Supplementary Figure S6 | Extended residual geometry', fontsize=13, fontweight='bold')
fig.savefig(OUTPUTS / 'supp_fig6.svg', bbox_inches='tight')
fig.savefig(OUTPUTS / 'supp_fig6.png', dpi=300, bbox_inches='tight')
plt.show()
