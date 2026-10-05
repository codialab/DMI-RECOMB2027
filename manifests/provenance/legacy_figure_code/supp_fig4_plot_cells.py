#!/usr/bin/env python3
"""Selected source plotting cells; see manifests/provenance/NOTEBOOK_EXTRACTIONS.json."""
# Provenance only. Do not execute as a current manuscript figure.


# Source notebook cell 1
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents] if (p / 'figures/supp_fig4/data/supp_fig4_validation_summary.json').is_file())
FIGURE_DIR = ROOT / 'figures/supp_fig4'
DATA_DIR = FIGURE_DIR / 'data'
OUTPUT_DIR = FIGURE_DIR / 'outputs'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
validation = json.loads((DATA_DIR / 'supp_fig4_validation_summary.json').read_text())
manifest = json.loads((DATA_DIR / 'supp_fig4_build_manifest.json').read_text())
assert validation['status'] == manifest['status'] == 'PASS'
correlations = pd.read_csv(DATA_DIR / 'supp_fig4_correlations.tsv', sep='\t')
contexts = pd.read_csv(DATA_DIR / 'supp_fig4_contexts.tsv', sep='\t')
bootstrap = pd.read_csv(DATA_DIR / 'supp_fig4_confirmation_bootstrap.tsv.gz', sep='\t')
assert len(correlations) == 2 and len(contexts) == 16 and len(bootstrap) == 5000
assert set(correlations['cohort']) == {'DEVELOPMENT', 'CONFIRMATION'}
assert np.isfinite(bootstrap['rho']).all()


# Source notebook cell 2
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9})
fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.55), constrained_layout=True, gridspec_kw={'width_ratios': [0.88, 1.15, 2.0]})
colors = {'DEVELOPMENT': '#3B6FB6', 'CONFIRMATION': '#D95F02'}

# A: Keep the populations as distinct estimates and denominators.
ax = axes[0]
rows = correlations.set_index('cohort').loc[['DEVELOPMENT', 'CONFIRMATION']].reset_index()
for y, row in enumerate(rows.itertuples(index=False)):
    ax.scatter(row.rho, y, s=55, color=colors[row.cohort], zorder=3)
    ax.text(row.rho + 0.025, y, f'{row.rho:.3f}  (n={row.finite_n:,})', va='center', fontsize=8)
ax.set_yticks([0, 1], ['Development', 'Held-out confirmation'])
ax.set_xlim(0, 1.12)
ax.set_ylim(-0.6, 1.6)
ax.set_xlabel('Spearman correlation, ρ')
ax.set_title('A  Separate populations')
ax.grid(axis='x', color='#dddddd', linewidth=0.6)
ax.spines[['top', 'right', 'left']].set_visible(False)
ax.tick_params(axis='y', length=0)

# B: Display the frozen reaction bootstrap and percentile stability interval.
ax = axes[1]
boot = validation['bootstrap']
ax.hist(bootstrap['rho'], bins=36, color='#7D9FCA', edgecolor='white', linewidth=0.45)
ax.axvspan(boot['q025'], boot['q975'], color='#E9B44C', alpha=0.23, label=f"95% interval [{boot['q025']:.3f}, {boot['q975']:.3f}]")
observed_rho = float(rows.loc[rows.cohort.eq('CONFIRMATION'), 'rho'].iloc[0])
ax.axvline(observed_rho, color='#B33A3A', linewidth=1.6, label=f'Observed ρ={observed_rho:.3f}')
ax.axvline(boot['q025'], color='#9C6B16', linestyle='--', linewidth=1)
ax.axvline(boot['q975'], color='#9C6B16', linestyle='--', linewidth=1)
ax.set_xlabel('Confirmation Spearman correlation, ρ')
ax.set_ylabel('Bootstrap replicates')
ax.set_title('B  Reaction bootstrap (5,000)')
ax.legend(frameon=False, fontsize=7, loc='upper left')
ax.grid(axis='y', color='#dddddd', linewidth=0.6)

# C: Preserve all method × RNA-context results and show each exact finite n.
ax = axes[2]
methods = ['GIMME', 'iMAT', 'CORDA', 'RIPTiDe']
method_colors = dict(zip(methods, ['#3B6FB6', '#D95F02', '#4C956C', '#845EC2']))
context_labels = ['setx1 + setx2', 'setx1 + setx2 + setx3', 'setx1 + setx3', 'setx2 + setx3']
contexts['method_order'] = contexts['algorithm'].map({m: i for i, m in enumerate(methods)})
contexts['context_order'] = contexts['context_label'].map({c: i for i, c in enumerate(context_labels)})
offsets = dict(zip(methods, [-0.27, -0.09, 0.09, 0.27]))
for method in methods:
    part = contexts.loc[contexts.algorithm.eq(method)].sort_values('context_order')
    yy = part.context_order.to_numpy(float) + offsets[method]
    xx = part.rho.to_numpy(float)
    ax.scatter(xx, yy, s=34, color=method_colors[method], label=method, zorder=3)
    for x, y, n in zip(xx, yy, part.finite_n):
        ax.annotate(f'{x:.2f} · {int(n)}', (x, y), xytext=(5, 0), textcoords='offset points', va='center', fontsize=6.6, color=method_colors[method])
ax.axvline(0, color='#666666', linewidth=0.8, linestyle=(0, (3, 3)))
ax.set_yticks(range(4), context_labels)
ax.invert_yaxis()
ax.set_xlim(0, 1.28)
ax.set_xlabel('Spearman correlation, ρ  ·  finite n')
ax.set_title('C  Confirmation by method and RNA context')
ax.legend(frameon=False, ncol=4, loc='lower center', bbox_to_anchor=(0.5, -0.24), fontsize=7)
ax.grid(axis='x', color='#dddddd', linewidth=0.6)
ax.spines[['top', 'right', 'left']].set_visible(False)
ax.tick_params(axis='y', length=0)

fig.suptitle('Held-out robustness of the geometry-to-utility relationship', fontsize=13)
fig.text(0.5, -0.045, 'A1 = 1-Strong-Anchor (HEX1). Predictor: sign-magnitude η²; response: directional-usefulness fraction. Bootstrap interval is computational reaction-level stability.', ha='center', fontsize=8)
fig.savefig(OUTPUT_DIR / 'supp_fig4.svg', bbox_inches='tight')
fig.savefig(OUTPUT_DIR / 'supp_fig4.png', dpi=300, bbox_inches='tight')
plt.show()
print('Wrote SVG and 300 dpi PNG to', OUTPUT_DIR)
