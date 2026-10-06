#!/usr/bin/env python3
"""Render audited Supplementary Figure S3 from included paired geometry.

Notebook provenance: figures/supp_fig3/supp_fig3.ipynb, cells 1–2;
see manifests/provenance/NOTEBOOK_EXTRACTIONS.json.
"""
from pathlib import Path
import hashlib, json, lzma, os, shutil, tempfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FormatStrFormatter
ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / 'data/figure_inputs/s3'
OUTPUTS = Path(os.environ['RECOMB_OUTPUT_DIR'])
manifest = json.loads((DATA / 'supp_fig3_build_manifest.json').read_text())
table_path = ROOT / manifest['output_table']
assert hashlib.sha256(table_path.read_bytes()).hexdigest() == manifest['output_sha256'], 'S3 plotting table checksum mismatch'
with tempfile.TemporaryDirectory(prefix='supp-s3-parquet-read-') as temporary:
    parquet_path = Path(temporary) / 'paired_geometry.parquet'
    with lzma.open(table_path, 'rb', format=lzma.FORMAT_XZ) as source, parquet_path.open('wb') as destination:
        shutil.copyfileobj(source, destination, length=1024 * 1024)
    paired = pd.read_parquet(parquet_path)
assert len(paired) == manifest['output_rows'] == 1_671_600
assert not paired.duplicated(['evaluation_id', 'reaction_id']).any()
assert set(paired['algorithm']) == {'CORDA', 'GIMME', 'iMAT', 'RIPTiDe'}
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 130})
metrics = [
    {'key': 'H_dir', 'a1': 'A1_H_dir', 'a2': 'A2_H_dir', 'title': r'Directional entropy $H_{\mathrm{dir}}$'},
    {'key': 'dominant_direction_mass', 'a1': 'A1_dominant_direction_mass', 'a2': 'A2_dominant_direction_mass', 'title': 'Dominant-direction mass'},
]
fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.6), constrained_layout=True)
last_hexbin = None; extent = (-1e-12, 1 + 1e-12, -1e-12, 1 + 1e-12)
for ax, metric in zip(axes, metrics):
    x = paired[metric['a1']].to_numpy(dtype=float); y = paired[metric['a2']].to_numpy(dtype=float)
    finite = np.isfinite(x) & np.isfinite(y); x, y = x[finite], y[finite]
    n = len(x); assert n == manifest['metric_checks'][metric['key']]['finite_pairs']
    hb = ax.hexbin(x, y, gridsize=55, bins='log', mincnt=1, cmap='YlOrRd', extent=extent, linewidths=0); last_hexbin = hb
    ax.plot([0, 1], [0, 1], linestyle='--', linewidth=1.15, color='#303030', zorder=3, label=r'Identity ($y=x$)')
    ax.set_xlim(extent[0], extent[1]); ax.set_ylim(extent[2], extent[3]); ax.set_aspect('equal', adjustable='box')
    ticks = np.linspace(0, 1, 6); ax.xaxis.set_major_locator(FixedLocator(ticks)); ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f')); ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.set_xlabel('A1 — 1-Strong-Anchor (glucose uptake)'); ax.set_ylabel('A2 — 2-Strong-Anchors (glucose uptake + LDH_L)'); ax.set_title(metric['title'])
    delta = y-x; tol = manifest['pairwise_direction_tolerance']; above=np.mean(delta>tol); equal=np.mean(np.abs(delta)<=tol); below=np.mean(delta < -tol)
    ax.text(.035,.965,f'n = {n:,}\nA2 > A1: {above:.1%}\nA2 ≈ A1: {equal:.1%}\nA2 < A1: {below:.1%}',transform=ax.transAxes,va='top',ha='left',fontsize=8,bbox={'facecolor':'white','edgecolor':'none','alpha':.82,'pad':3})
    ax.legend(loc='lower right',frameon=False,fontsize=8)
cbar=fig.colorbar(last_hexbin,ax=axes,shrink=.88,pad=.025); cbar.set_label('Matched reaction–evaluation pairs per hexagon (log scale)')
fig.suptitle('Supplementary Figure S3 | A1 versus A2 paired geometry',fontsize=12)
OUTPUTS.mkdir(parents=True,exist_ok=True); fig.savefig(OUTPUTS/'supp_fig3.svg',bbox_inches='tight'); fig.savefig(OUTPUTS/'supp_fig3.png',dpi=300,bbox_inches='tight')
