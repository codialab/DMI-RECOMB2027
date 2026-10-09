# S6 validation report

## Inputs and provenance

- The figure uses the four frozen, checksum-verified paired-geometry partitions indexed by `manifests/FIGURE_INPUTS.json` and described by `data/figure_inputs/fig4/fig4_data_manifest.json`.
- No reconstruction, optimization, sampling, candidate generation, external data, or new S6 input table is used.
- Plot implementation and source notebook adaptation are identified in `figures/supplementary/s6/PROVENANCE.json`.

## Population and metric checks

- **Matched key population:** 1,671,600 unique method/evaluation/reaction keys = 400 evaluations × 4,179 reactions; each method contributes 417,900 keys.
- **S6A η²:** values are read under A1 and A2 and paired changes are computed as A2 − A1. Both arms must be finite. The 573,882 finite paired keys are CORDA 212,164; GIMME 141,365; iMAT 207,930; RIPTiDe 12,423. Degenerate values are not imputed.
- **S6B supported sign states:** all 1,671,600 keys are included. A state count is the number of negative, tie/near-zero, and positive sign masses with positive probability mass. Transition counts are 1→1 1,293,850; 1→2 1,055; 2→2 376,470; 3→3 225; all other cells are zero.
- **Retained non-tie coverage diagnostic:** all 1,671,600 keys are available in the same frozen inputs. Mean paired change is −0.00003375; 776 keys (0.046%) differ by more than 10⁻¹². This metric remains documented but is not a graphical panel.

## Interpretation and limits

The paired η² changes are heterogeneous: 24.8% decreased, 26.5% were unchanged, and 48.7% increased; median change is zero and mean change is +0.001429 over the finite η² population. Supported-state count is unchanged for 99.9369% of keys. These are descriptive reaction–evaluation summaries, not independent biological replicates or evidence of a uniform method effect. They do not imply expansion or contraction of a feasible flux region.

S6 contains only the two panels specified by the current Supplementary Materials. The eight paired descriptor summaries belong to Supplementary Table S7. S6 panel labels, metric names, and denominators match that editorial specification.

## Rendering verification

`python reproduce.py figures --only fig1,s6 --output-root reproduced/figures/sync_20261009_final` completed both jobs successfully. The regenerated S6 PNG is pixel-identical to the current source notebook PNG. SVG output is vector-equivalent; byte differences are expected from Matplotlib creation timestamps and generated element IDs. The regenerated Figure 1 PNG retains panels A–D, the caption-aligned panel C examples, and corrected panel B sentence. Output hashes are recorded in `manifests/ARTWORK.json` and `manifests/checksums/SHA256SUMS.txt`.
