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



# Source notebook cell 3
geometry_files = sorted((DATA / 'fig4_geometry_paired').rglob('*.parquet.xz'))
geometry_frames = []
for path in geometry_files:
    with tempfile.TemporaryDirectory(prefix='fig4-panelB-parquet-') as temporary:
        parquet_path = Path(temporary) / 'geometry.parquet'
        with lzma.open(path, 'rb', format=lzma.FORMAT_XZ) as source, parquet_path.open('wb') as destination:
            shutil.copyfileobj(source, destination, length=1024 * 1024)
        geometry_frames.append(pd.read_parquet(parquet_path).assign(algorithm=path.parent.name.split('=', 1)[1]))
geometry = pd.concat(geometry_frames, ignore_index=True)
assert not geometry.duplicated(['evaluation_id', 'reaction_id']).any()
assert len(geometry) == manifest['row_counts']['geometry_paired']
assert set(geometry.algorithm.unique()) == {'CORDA', 'GIMME', 'iMAT', 'RIPTiDe'}
assert geometry.paired_geometry_core_finite.notna().all()
required_panelB_columns = {
    'A1_H_dir','A2_H_dir',
    'A1_dominant_direction_mass','A2_dominant_direction_mass',
    'A1_delta_v_B_width80','A2_delta_v_B_width80',
}
missing_panelB_columns = sorted(required_panelB_columns - set(geometry.columns))
assert not missing_panelB_columns, f'Geometry table is missing Panel B columns: {missing_panelB_columns}. Rebuild FIG4-DATA-B with schema v1.1.0.'
metric_config = [
    {'key':'H_dir','a1':'A1_H_dir','a2':'A2_H_dir','label':'Direction entropy','scale':'linear'},
    {'key':'dominant_direction_mass','a1':'A1_dominant_direction_mass','a2':'A2_dominant_direction_mass','label':'Probability mass in dominant direction','scale':'linear'},
    {'key':'delta_v_B_width80','a1':'A1_delta_v_B_width80','a2':'A2_delta_v_B_width80','label':r'Central 80% width of $P(\Delta v_B)$','scale':'linear'},
]
paired_data = {}
for metric in metric_config:
    a1 = pd.to_numeric(geometry[metric['a1']], errors='coerce').to_numpy(float)
    a2 = pd.to_numeric(geometry[metric['a2']], errors='coerce').to_numpy(float)
    mask = np.isfinite(a1) & np.isfinite(a2)
    paired_data[metric['key']] = (a1[mask], a2[mask], a2[mask]-a1[mask], geometry.loc[mask,'algorithm'].to_numpy())
[(m['key'], len(paired_data[m['key']][2])) for m in metric_config]

# Source notebook cell 6
fig, axes = plt.subplots(
    1, 2,
    figsize=(10, 4.2),
    constrained_layout=True
)

max_curve_points = 2000

for ii, (ax, metric) in enumerate(zip(axes, metric_config[:2])):
    a1, a2, shift, methods = paired_data[metric['key']]

    a1 = np.asarray(a1, dtype=float)
    a2 = np.asarray(a2, dtype=float)

    a1 = a1[np.isfinite(a1)]
    a2 = a2[np.isfinite(a2)]

    # Exact sorting
    a1_sorted = np.sort(a1)
    a2_sorted = np.sort(a2)

    # Downsample the ECDF by rank for plotting only
    n1 = len(a1_sorted)
    n2 = len(a2_sorted)

    idx1 = np.unique(
        np.linspace(
            0, n1 - 1,
            min(max_curve_points, n1),
            dtype=int
        )
    )

    idx2 = np.unique(
        np.linspace(
            0, n2 - 1,
            min(max_curve_points, n2),
            dtype=int
        )
    )

    x1 = a1_sorted[idx1]
    y1 = (idx1 + 1) / n1

    x2 = a2_sorted[idx2]
    y2 = (idx2 + 1) / n2

    ax.step(
        x1, y1,
        where='post',
        linewidth=2.0,
        color='#3568a8',
        label='1-Strong-Anchor\n(HEX1)'
    )

    ax.step(
        x2, y2,
        where='post',
        linewidth=2.0,
        color='#d06a3b',
        label='2-Strong-Anchors\n(HEX1 + LDH_L)'
    )

    if metric['label'] == 'Direction entropy':
        ax.set_xlabel(r'Direction entropy $H_{dir}$', fontsize=12)
    else:
        ax.set_xlabel(metric['label'], fontsize=12)
    if ii == 0:
        ax.set_ylabel('Cumulative fraction')
    ax.set_ylim(0, 1)

    if ii == 1:        
        ax.legend(frameon=False)

fig.suptitle(
    r'Flux-distribution P$(\Delta v_B)$ descriptors after anchoring',
    y=1.04,
    fontsize=13,
)
fig.text(
    0.0, 1.00,
    r'B',
    ha='left',
    fontsize=14,
    fontweight='bold'
)
fig.tight_layout(pad=0.5, h_pad=0.5, w_pad=0.5)

fig.savefig(
    PANEL/'fig4B_metric_distributions.svg',
    bbox_inches='tight'
)

fig.savefig(
    PANEL/'fig4B_metric_distributions.pdf',
    bbox_inches='tight'
)

plt.show()
