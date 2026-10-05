#!/usr/bin/env python3
"""Selected source plotting cells; see manifests/provenance/NOTEBOOK_EXTRACTIONS.json."""
# Provenance only. Do not execute as a current manuscript figure.


# Source notebook cell 1
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path.cwd()
FIGDIR = HERE if (HERE / 'data').is_dir() else HERE / 'figures/supp_fig5'
DATA = FIGDIR / 'data'
OUT = FIGDIR / 'outputs'
context = pd.read_csv(DATA / 'supp_fig5_context_gains.tsv', sep='\t')
outcomes = pd.read_csv(DATA / 'supp_fig5_outcome_composition.tsv', sep='\t')
coverage = pd.read_csv(DATA / 'supp_fig5_truth_tie_coverage.tsv', sep='\t')
manifest = json.loads((DATA / 'supp_fig5_build_manifest.json').read_text())
assert manifest['status'] == 'PASS'
assert len(context) == 64 and len(outcomes) == 4 and len(coverage) == 2
assert context[['mean_gain','q05','q95']].notna().all().all()
assert np.isfinite(context[['mean_gain','q05','q95']].to_numpy()).all()
assert np.allclose(outcomes[['improved_fraction','gain_tie_fraction','harmed_fraction']].sum(axis=1), 1.0, atol=1e-12, rtol=0)
assert (context['n_reaction_evaluations'] > 0).all()

plt.rcParams.update({
    'font.family':'DejaVu Sans', 'font.size':9, 'axes.titlesize':10,
    'axes.labelsize':9, 'xtick.labelsize':8, 'ytick.labelsize':8,
    'svg.fonttype':'none', 'axes.spines.top':False, 'axes.spines.right':False,
})
colors={'A1':'#2878B5','A2':'#D9534F'}
markers={'wrong':'v','correct':'^'}
contexts=sorted(context['rna_context_key'].unique())
algorithms=['CORDA','GIMME','iMAT','RIPTiDe']
fig=plt.figure(figsize=(14.2,10.4))
gs=fig.add_gridspec(2,4,height_ratios=[1.65,1.0],left=.11,right=.99,top=.83,bottom=.22,wspace=.10,hspace=.46)
fig.suptitle('Supplementary Figure S5 | Directional-cue utility controls',fontsize=15,fontweight='bold',y=.99)
fig.text(.11,.945,'A  Method- and context-stratified gain',fontsize=12,fontweight='bold',va='top')
fig.text(.11,.908,'Whiskers show descriptive 5th–95th percentiles; each method has its own x-axis scale. Ranges are not confidence intervals.',fontsize=8.2,va='top')
axes=[fig.add_subplot(gs[0,i]) for i in range(4)]

for ax, algorithm in zip(axes, algorithms):
    block=context.loc[context.algorithm.eq(algorithm)]
    short_contexts=sorted(block.rna_context_key.unique())
    # Within each RNA context, offset the four anchor/cue summaries for legibility.
    offsets={('A1','wrong'):-0.24,('A1','correct'):-0.08,('A2','wrong'):0.08,('A2','correct'):0.24}
    for yi, rna_context in enumerate(short_contexts):
        for (anchor,cue),offset in offsets.items():
            row=block.loc[block.anchor_setting.eq(anchor)&block.cue_direction.eq(cue)&block.rna_context_key.eq(rna_context)]
            if len(row)!=1: raise ValueError(f'Expected one summary for {algorithm}, {rna_context}, {anchor}, {cue}')
            row=row.iloc[0]
            low=row['mean_gain']-row['q05']; high=row['q95']-row['mean_gain']
            ax.errorbar(row['mean_gain'],yi+offset,xerr=np.array([[max(0,low)],[max(0,high)]]),
                        fmt=markers[cue],markersize=4.5,color=colors[anchor],
                        ecolor=colors[anchor],elinewidth=1,capsize=1.5,alpha=.9,zorder=3)
    ax.axvline(0,color='#555555',linewidth=.8,linestyle='--',zorder=1)
    ax.set_title(algorithm)
    ax.set_yticks(range(len(short_contexts)))
    ax.set_yticklabels([c.replace('training_samples=','').replace(',',' + ') for c in short_contexts])
    ax.invert_yaxis()
    ax.grid(axis='x',color='#dddddd',linewidth=.6)
    ax.set_xlabel('Mean gain (mmol gDW$^{-1}$ h$^{-1}$)')
    ax.ticklabel_format(axis='x',style='sci',scilimits=(-3,3),useMathText=True)
    ax.xaxis.get_offset_text().set_fontsize(7)
    if ax is not axes[0]:
        ax.tick_params(axis='y',labelleft=False)
axes[0].set_ylabel('RNA training context')
legend=[
 Line2D([0],[0],marker='v',color='none',markerfacecolor=colors['A1'],markeredgecolor=colors['A1'],label='A1 — wrong cue'),
 Line2D([0],[0],marker='^',color='none',markerfacecolor=colors['A1'],markeredgecolor=colors['A1'],label='A1 — correct cue'),
 Line2D([0],[0],marker='v',color='none',markerfacecolor=colors['A2'],markeredgecolor=colors['A2'],label='A2 — wrong cue'),
 Line2D([0],[0],marker='^',color='none',markerfacecolor=colors['A2'],markeredgecolor=colors['A2'],label='A2 — correct cue'),
]
fig.legend(handles=legend,loc='upper right',bbox_to_anchor=(.99,.955),ncol=2,frameon=False,fontsize=7.6,handletextpad=.4,columnspacing=1.1)

ax=fig.add_subplot(gs[1,:])
order=[('A1',0.0),('A1',1.0),('A2',0.0),('A2',1.0)]
labels=[]
segments=[('improved_fraction','Improved','#4C956C'),('gain_tie_fraction','No gain change','#B8B8B8'),('harmed_fraction','Harmed','#C75B5B')]
for yi,(anchor,q) in enumerate(order):
    row=outcomes.loc[outcomes.anchor_setting.eq(anchor)&outcomes.reliability_q.eq(q)]
    if len(row)!=1: raise ValueError(f'Missing outcome row for {anchor}, q={q}')
    row=row.iloc[0]
    left=0.0
    for col,label,color in segments:
        value=float(row[col])
        ax.barh(yi,value,left=left,height=.62,color=color,edgecolor='white',linewidth=.6)
        if value>=.075:
            ax.text(left+value/2,yi,f'{value:.1%}',ha='center',va='center',fontsize=8,color='white' if color!='#B8B8B8' else '#222222',fontweight='bold')
        left+=value
    labels.append(f"{anchor} · {row['cue_label']}  (n={int(row['non_tie_denominator']):,})")
ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels); ax.invert_yaxis()
ax.set_xlim(0,1); ax.set_xticks(np.linspace(0,1,6)); ax.set_xticklabels([f'{x:.0%}' for x in np.linspace(0,1,6)])
ax.set_xlabel('Share of non-tie reaction × truth cases')
ax.set_title('B  Frozen gain outcomes at wrong- and correct-cue endpoints (lambda = 0.25)',loc='left',fontweight='bold',pad=10)
ax.grid(axis='x',color='#dddddd',linewidth=.6,zorder=0)
fig.legend(handles=[plt.Rectangle((0,0),1,1,color=color,label=label) for _,label,color in segments],
          loc='lower center',bbox_to_anchor=(.55,.145),ncol=3,frameon=False)
truth_note='Truth ties excluded from these bars: '+ '  |  '.join(
    f"{r.anchor_setting}: {int(r.truth_tie_count):,}/{int(r.total_cases):,} ({r.truth_tie_fraction:.1%})"
    for r in coverage.itertuples())
fig.text(.11,.085,truth_note,ha='left',va='top',fontsize=8)

OUT.mkdir(parents=True,exist_ok=True)
fig.savefig(OUT/'supp_fig5.svg',bbox_inches='tight')
svg_path=OUT/'supp_fig5.svg'
svg_path.write_text('\n'.join(line.rstrip() for line in svg_path.read_text().splitlines())+'\n')
fig.savefig(OUT/'supp_fig5.png',dpi=300,bbox_inches='tight')
plt.show()
print('Panel A matched reaction-evaluation rows:',manifest['gain_validation']['anchor_rows'])
print('Panel B non-tie denominators:',outcomes.set_index(['anchor_setting','reliability_q'])['non_tie_denominator'].to_dict())
print('Truth tie counts:',coverage.set_index('anchor_setting')[['truth_tie_count','total_cases']].to_dict(orient='index'))
print('Saved:',OUT/'supp_fig5.svg',OUT/'supp_fig5.png')

