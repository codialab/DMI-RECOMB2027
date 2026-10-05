#!/usr/bin/env python3
"""S1 source plotting cell; original display behavior preserved for review."""
import matplotlib
matplotlib.use("Agg")
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import Normalize

import os
ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / 'data/figure_inputs/s1'
OUT = Path(os.environ['RECOMB_OUTPUT_DIR'])
OUT.mkdir(parents=True, exist_ok=True)
weights = pd.read_csv(DATA / 'candidate_weights.tsv', sep='\t')
strata = pd.read_csv(DATA / 'conditional_stratum_summary.tsv', sep='\t')
global_ess = pd.read_csv(DATA / 'global_ess_summary.tsv', sep='\t')
contrast = pd.read_csv(DATA / 'contrast_product_support.tsv', sep='\t')
holdout = pd.read_csv(DATA / 'post_holdout_pair_support.tsv', sep='\t')
assert len(strata) == 320 and len(weights) == 6400
assert holdout.truth_pair_id.nunique() == 836
assert len(global_ess) == 20 and len(contrast) == 800
methods = ['iMAT', 'GIMME', 'CORDA', 'RIPTiDe']
contexts = ['setx1,setx2', 'setx1,setx3', 'setx2,setx3', 'setx1,setx2,setx3']
context_order = ['training_samples=setx1,setx2', 'training_samples=setx1,setx3', 'training_samples=setx2,setx3', 'training_samples=setx1,setx2,setx3']
context_short = {c: f'RNA context {i+1}' for i, c in enumerate(context_order)}
colors = {'iMAT':'#4477AA', 'GIMME':'#EE6677', 'CORDA':'#228833', 'RIPTiDe':'#AA3377'}
anchor_titles = {'A1':'A1 · 1-Strong-Anchor (HEX1)', 'A2':'A2 · 2-Strong-Anchors (HEX1 + LDH_L)'}
fig = plt.figure(figsize=(15.5, 18), constrained_layout=True)
gs = fig.add_gridspec(3, 4, height_ratios=[1.0, 1.35, 0.85])
fig.suptitle('Supplementary Figure S1 | Candidate weights and effective support', fontsize=17, fontweight='bold')

# A: every conditional stratum, split by tumor and anchor.
for r, tumor in enumerate(['CT2A', 'GL261']):
    for c, anchor in enumerate(['A1', 'A2']):
        ax = fig.add_subplot(gs[0, c + 2*r])
        sub = strata[(strata.tumor == tumor) & (strata.anchor == anchor)]
        xlabels = [str(j) for _m in methods for j, _ctx in enumerate(contexts, start=1)]
        # Stable display order matches the four frozen RNA contexts.
        xpos = {(m, 'training_samples='+ctx): i for i, (m, ctx) in enumerate((m, ctx) for m in methods for ctx in contexts)}
        for m in methods:
            for ctx in contexts:
                key = 'training_samples='+ctx
                block = sub[(sub.method == m) & (sub.rna_context_key == key)].sort_values('mouse_id')
                x = xpos[(m, key)]
                jitter = np.linspace(-0.13, 0.13, len(block))
                ax.scatter(x + jitter, block.conditional_ess, s=15, color=colors[m], alpha=.78, linewidths=0, zorder=3)
        ax.axhline(20, color='#555555', linestyle=':', linewidth=1)
        ax.set_xticks(range(16), xlabels, fontsize=7)
        for boundary in [3.5, 7.5, 11.5]: ax.axvline(boundary, color='#aaaaaa', linewidth=.6, alpha=.7)
        ax.set_ylim(0, 20.5)
        ax.set_ylabel('Conditional ESS (20 candidates)')
        ax.set_title(f'{tumor} · {anchor_titles[anchor]} · 80 strata', fontsize=9)
        ax.grid(axis='y', color='#dddddd', linewidth=.5)
        ax.spines[['top','right']].set_visible(False)
        ax.set_xlabel('Context index (1–4); method shown by color', fontsize=7)
        if r == 0 and c == 0:
            ax.legend(handles=[Patch(facecolor=colors[m], label=m) for m in methods], loc='upper left', ncol=2, fontsize=6, frameon=False)

# B: all per-candidate conditional weights, rows are the 320 strata.
for r, tumor in enumerate(['CT2A', 'GL261']):
    for c, anchor in enumerate(['A1', 'A2']):
        ax = fig.add_subplot(gs[1, c + 2*r])
        sub = weights[(weights.tumor == tumor) & (weights.anchor == anchor)].copy()
        sub['context_display'] = sub.rna_context_key.map(context_short)
        sub = sub.sort_values(['method', 'rna_context_key', 'mouse_id', 'candidate_index'], kind='stable')
        matrix = sub.pivot(index=['method','rna_context_key','mouse_id'], columns='candidate_index', values='conditional_weight')
        matrix = matrix.reindex(pd.MultiIndex.from_frame(sub[['method','rna_context_key','mouse_id']].drop_duplicates().sort_values(['method','rna_context_key','mouse_id'])), columns=range(20))
        im = ax.imshow(matrix.to_numpy(), aspect='auto', interpolation='nearest', cmap='magma', norm=Normalize(0, 1))
        ticks, labels = [], []
        for i, (m, ctx) in enumerate((m, 'training_samples='+ctx) for m in methods for ctx in contexts):
            ticks.append(i*5 + 2)
            labels.append(f'{m} · {context_short[ctx]}')
        ax.set_yticks(ticks, labels, fontsize=6)
        ax.set_xticks([0,4,9,14,19], ['1','5','10','15','20'], fontsize=7)
        ax.set_xlabel('Candidate index (ordered as frozen)')
        ax.set_title(f'{tumor} · {anchor_titles[anchor]} · 80 strata', fontsize=9)
        ax.set_ylabel('Method · RNA context (five mouse strata each)')
        cb = fig.colorbar(im, ax=ax, fraction=.04, pad=.02)
        cb.set_label('Conditional candidate weight', fontsize=7)
        cb.ax.tick_params(labelsize=6)

# C: three non-interchangeable support diagnostics.
ax = fig.add_subplot(gs[2, 0])
for i, (anchor, color) in enumerate([('A1','#4477AA'),('A2','#CC6677')]):
    d = global_ess[global_ess.anchor == anchor]
    jitter = np.linspace(-.08,.08,len(d))
    ax.scatter(np.full(len(d),i)+jitter,d.global_ess,s=20,color=color,alpha=.8)
ax.axhline(20,color='#555555',linestyle=':')
ax.set_xticks([0,1], ['A1','A2'])
ax.set_ylim(19.5,20.5)
ax.set_title('Global ESS · 320 states · target ≈20', fontsize=9)
ax.set_ylabel('Global ESS')
ax.grid(axis='y',color='#dddddd',linewidth=.5)
ax.spines[['top','right']].set_visible(False)

ax = fig.add_subplot(gs[2, 1])
for i, (anchor, color) in enumerate([('A1','#4477AA'),('A2','#CC6677')]):
    d = contrast[contrast.anchor == anchor].contrast_product_ess
    parts = ax.violinplot(d, positions=[i], widths=.72, showmedians=True, showextrema=False)
    for body in parts['bodies']:
        body.set_facecolor(color); body.set_edgecolor(color); body.set_alpha(.45)
    parts['cmedians'].set_color('#222222')
ax.set_xticks([0,1], ['A1','A2'])
ax.set_title('Product-weight ESS · 400 pairs/evaluation', fontsize=9)
ax.set_ylabel('Product-weight ESS')
ax.grid(axis='y',color='#dddddd',linewidth=.5)
ax.spines[['top','right']].set_visible(False)

ax = fig.add_subplot(gs[2, 2:])
for i, (anchor, color) in enumerate([('A1','#4477AA'),('A2','#CC6677')]):
    d = holdout[holdout.anchor == anchor]
    for j, cond in enumerate(['CT2A','GL261']):
        v = d.loc[d.condition == cond, 'post_holdout_ess'].to_numpy(float)
        x = i*2+j
        jitter = np.linspace(-.18,.18,len(v))
        ax.scatter(np.full(len(v),x)+jitter,v,s=4,color=color,alpha=.22,rasterized=True)
        ax.plot([x-.18,x+.18],[np.median(v)]*2,color='#222222',linewidth=1.2)
ax.set_xticks([0,1,2,3], ['A1 · CT2A','A1 · GL261','A2 · CT2A','A2 · GL261'], rotation=18, ha='right')
ax.set_title('Post-holdout ESS · fixed 836 matched pairs', fontsize=9)
ax.set_ylabel('Condition-specific post-holdout ESS')
ax.grid(axis='y',color='#dddddd',linewidth=.5)
ax.spines[['top','right']].set_visible(False)

fig.text(.01, .005, 'Weights and ESS are finite-weight diagnostics, not counts of independent observations. Flux units do not apply to these dimensionless quantities.', fontsize=8)
for ext in ['svg','png']:
    output_path = OUT / f'supp_fig1.{ext}'
    fig.savefig(output_path, dpi=320, bbox_inches='tight')
    if ext == 'svg':
        output_path.write_text('\n'.join(line.rstrip() for line in output_path.read_text().splitlines()) + '\n')
plt.show()
