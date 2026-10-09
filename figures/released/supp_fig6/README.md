# Supplementary Figure S6

S6 compares the frozen A1 glucose-uptake anchor with A2, which adds the LDH_L-based lactate soft anchor. It has two panels, matching the current Supplementary Materials: method-stratified paired changes in direction-explained magnitude variance (η²), and the pooled transition matrix for supported sign-state counts.

## Reproduction

From the repository root, run:

```bash
python reproduce.py validate
python reproduce.py figures --only s6
```

The plotting script reads the four checksum-verified Figure 4 paired-geometry partitions already included in `data/figure_inputs/fig4/fig4_geometry_paired/`. It writes SVG and 300 dpi PNG outputs under `reproduced/figures/s6/`. The standalone source is `figures/supplementary/s6/plot_s6.py`; source notebook and adaptation provenance are recorded in `figures/supplementary/s6/PROVENANCE.json`. No S6-specific numerical inputs are added or rebuilt.

## Panels and frozen populations

- **S6A:** Method-stratified ECDFs of A2 − A1 η² across the 573,882 finite paired reaction–evaluation keys. The legend gives the per-method counts: CORDA 212,164; GIMME 141,365; iMAT 207,930; and RIPTiDe 12,423. Paired changes are calculated from the frozen A1 and A2 values. η² retains its production definition and nonfinite values for degenerate magnitude distributions are excluded from this panel.
- **S6B:** Pooled A1-to-A2 transition matrix for the number of supported sign states among the negative, tie/near-zero, and positive probability masses. It uses all 1,671,600 matched reaction–evaluation keys (400 evaluations × 4,179 reactions); cells show counts and fractions.

The reaction–evaluation rows are descriptive computational units, not independent biological replicates. Non-tie coverage remains in the frozen input and analysis: its mean A2 − A1 change is −0.00003375 across all keys, and 776 keys (0.046%) change by more than 10⁻¹². It is not plotted in S6. The eight paired descriptor summaries are reported in Supplementary Table S7, not as additional S6 panels. Direction entropy and dominant-direction mass remain in S2/S3, and product-weight ESS remains in S1.

The figure does not alter the frozen evaluation population, anchor operators, reaction split, or numerical inputs. It makes no claim that the feasible flux region expanded or contracted.
