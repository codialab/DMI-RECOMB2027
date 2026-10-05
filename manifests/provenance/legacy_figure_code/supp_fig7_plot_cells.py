#!/usr/bin/env python3
"""Selected source plotting cells; see manifests/provenance/NOTEBOOK_EXTRACTIONS.json."""
# Provenance only. Do not execute as a current manuscript figure.


# Source notebook cell 1
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def find_root(start):
    for candidate in (start, *start.parents):
        if (candidate / 'figures/fig4/data/fig4_data_manifest.json').is_file():
            return candidate
    raise FileNotFoundError('Could not locate repository root from notebook working directory')

ROOT = find_root(Path.cwd().resolve())
OUT = ROOT / 'figures/supp_fig7'
DATA = OUT / 'data'
OUTPUTS = OUT / 'outputs'
OUTPUTS.mkdir(parents=True, exist_ok=True)
manifest = json.loads((DATA / 'supp_fig7_manifest.json').read_text())
if manifest.get('status') != 'PASS':
    raise RuntimeError('S7 table build did not pass')
examples = pd.read_csv(DATA / 'supp_fig7_example_selection.tsv', sep='\t')
support = pd.read_csv(DATA / 'supp_fig7_weighted_distributions.tsv.gz', sep='\t')
examples = examples.sort_values(['example_order', 'method'], kind='stable')
support = support.sort_values(['example_order', 'method', 'arm', 'support_id'], kind='stable')
print(examples[['example_order', 'example_label', 'reaction_id', 'method', 'rna_context_key', 'ct2a_mouse', 'gl261_mouse', 'evaluation_id', 'geometry_rank']].to_string(index=False))
print("Support rows: {:,}; example/method/arm distributions: {}".format(len(support), support.groupby(['candidate_group_id', 'method', 'arm']).ngroups))

# Source notebook cell 2
methods = ['CORDA', 'GIMME', 'iMAT', 'RIPTiDe']
arms = [('A1', '#3568a8', 'A1 — 1-Strong-Anchor (HEX1)'),
        ('A2-L', '#d06a3b', 'A2 — 2-Strong-Anchors (HEX1 + LDH_L)')]
example_rows = (examples[['example_order', 'example_label', 'candidate_group_id', 'reaction_id',
                        'rna_context_key', 'ct2a_mouse', 'gl261_mouse', 'geometry_rank',
                        'median_absolute_A2_minus_A1_H_dir_shift']]
               .drop_duplicates('candidate_group_id').sort_values('example_order'))
fig, axes = plt.subplots(3, 4, figsize=(15.5, 10.2), sharey=False, constrained_layout=True)
max_mass_error = 0.0
for row_index, example in enumerate(example_rows.itertuples(index=False)):
    group = support.loc[support.candidate_group_id.eq(example.candidate_group_id)]
    for col_index, method in enumerate(methods):
        ax = axes[row_index, col_index]
        method_data = group.loc[group.method.eq(method)]
        all_delta = method_data.delta_v_B.to_numpy(dtype=float)
        lo, hi = float(np.min(all_delta)), float(np.max(all_delta))
        if hi <= lo:
            pad = max(abs(lo) * 0.05, 1e-12)
            lo, hi = lo - pad, hi + pad
        else:
            pad = 0.03 * (hi - lo)
            lo, hi = lo - pad, hi + pad
        bins = np.linspace(lo, hi, 31)
        for arm, color, label in arms:
            arm_data = method_data.loc[method_data.arm.eq(arm)]
            delta = arm_data.delta_v_B.to_numpy(dtype=float)
            weight = arm_data.weight.to_numpy(dtype=float)
            if len(delta) != 400 or not np.isfinite(delta).all() or not np.isfinite(weight).all():
                raise RuntimeError(f'Invalid support in {example.candidate_group_id}/{method}/{arm}')
            weight = weight / weight.sum()
            density, edges = np.histogram(delta, bins=bins, weights=weight, density=True)
            mass = float(np.sum(density * np.diff(edges)))
            max_mass_error = max(max_mass_error, abs(mass - 1.0))
            if not np.isclose(mass, 1.0, rtol=0, atol=1e-10):
                raise RuntimeError(f'Weighted histogram mass is {mass} for {example.candidate_group_id}/{method}/{arm}')
            ax.stairs(density, edges, linewidth=1.8, color=color, label=label)
        ax.axvline(0, color='black', linestyle='--', linewidth=0.9)
        ax.set_title(method, fontsize=10, fontweight='bold')
        ax.grid(axis='y', alpha=0.18, linewidth=0.5)
        if row_index == 2:
            ax.set_xlabel(r'$\Delta v_B$ [mmol gDW$^{-1}$ h$^{-1}$]', fontsize=8)
        if col_index == 0:
            ax.set_ylabel('Probability density', fontsize=8)
        if row_index == 0 and col_index == 0:
            ax.legend(frameon=False, fontsize=8, loc='best')
    row_label = (f'{example.example_label}: {example.reaction_id}\n'
                 f'{example.rna_context_key.replace("training_samples=", "")}; CT2A {example.ct2a_mouse} vs GL261 {example.gl261_mouse}; rank {example.geometry_rank}\n'
                 f'median |A2−A1 H_dir| = {example.median_absolute_A2_minus_A1_H_dir_shift:.3g}')
    axes[row_index, 0].text(-0.32, 0.5, row_label, transform=axes[row_index, 0].transAxes,
                           ha='right', va='center', fontsize=8, wrap=True)

fig.suptitle('Representative weighted residual distributions: A1 vs A2', fontsize=14)
fig.text(0.5, -0.025, 'Each curve uses 400 frozen Cartesian states; weights are normalized product weights.',
         ha='center', fontsize=9)
svg_path = OUTPUTS / 'supp_fig7.svg'
png_path = OUTPUTS / 'supp_fig7.png'
fig.savefig(svg_path, bbox_inches='tight')
fig.savefig(png_path, dpi=300, bbox_inches='tight')
svg_text = svg_path.read_text()
svg_path.write_text('\n'.join(line.rstrip() for line in svg_text.splitlines()) + '\n')
plt.show()

def sha256_file(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()

manifest['render_validation'] = {'maximum_weighted_histogram_mass_error': float(max_mass_error), 'histogram_mass_tolerance': 1e-10}
manifest['rendered_output_hashes'] = {path.relative_to(ROOT).as_posix(): sha256_file(path) for path in (svg_path, png_path)}
(DATA / 'supp_fig7_manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + '\n')
print(f'Maximum weighted histogram mass error: {max_mass_error:.3g}')
print(f'Saved {svg_path.relative_to(ROOT)} and {png_path.relative_to(ROOT)}')
