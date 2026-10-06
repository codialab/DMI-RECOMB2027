# Supplementary Figure S6

## Objective

S6 extends the residual-geometry diagnostics with sign–magnitude coupling, non-tie coverage, and the number of supported sign states. It compares the frozen production A1 setting (1-Strong-Anchor, glucose uptake) with A2 (2-Strong-Anchors, glucose uptake + LDH_L).

## Authoritative inputs and build

The builder reads the four manifest-verified paired geometry partitions in `figures/fig4/data/fig4_geometry_paired/` and verifies their hashes against `figures/fig4/data/fig4_data_manifest.json`. That Figure 4 table pairs A1 and A2 on `(evaluation_id, reaction_id)`. Upstream geometry is from frozen PL1 (A1) and A21 (A2) artifacts; no reconstruction, optimization, sampling, candidate generation, or production pipeline is rerun.

From the repository root:

```bash
python figures/supp_fig6/build_supp_fig6_tables.py
jupyter nbconvert --to notebook --execute --inplace figures/supp_fig6/supp_fig6.ipynb
```

The builder writes `data/supp_fig6_geometry.tsv.gz`, `data/supp_fig6_sign_state_transitions.tsv`, and `data/supp_fig6_build_manifest.json`. The notebook writes `outputs/supp_fig6.svg` and `outputs/supp_fig6.png` (300 dpi).

## Panels and aggregation

- **A:** Method-stratified ECDFs of paired change in sign–magnitude η². Production defines η² as the weighted fraction of variance in `|Δv_B|` explained by sign category {-1, 0, +1}. It lies in [0, 1]; it is nonfinite when magnitude variance is degenerate. Only rows finite under both anchors enter paired differences.
- **B:** Method-stratified ECDFs of paired change in non-tie coverage, defined as `1 − p_tie`, where ties use the production tolerance of 1e-12. Values lie in [0, 1].
- **C:** Pooled A1-to-A2 transition matrix for supported sign-state count, the number of sign categories with positive probability mass. Categories are 1, 2, and 3. Cells show counts and fractions of paired rows.

The starting geometry population has 1,671,600 paired reaction–evaluation observations: 400 evaluations × 4,179 reactions, with 417,900 observations per method. η² has 573,882 finite paired rows (CORDA 212,164; GIMME 141,365; iMAT 207,930; RIPTiDe 12,423). Non-tie coverage and sign-state count each use all 1,671,600 rows. These descriptive rows are not independent biological replicates.

## Observed patterns

The median paired change is zero for all three descriptors. Mean η² change is 0.00143 across finite pairs; mean non-tie coverage change is −0.0000337 across the full population. Sign-state count is unchanged in 99.94% of paired rows; the remaining observed transitions are from one to two supported states. These summaries describe the frozen population and do not imply a uniform effect across methods or reactions.

## Descriptor selection and limitations

Joint ESS was excluded because S1 already presents candidate and effective-support diagnostics; it describes effective weighted support, not independent sample size. `width80` was excluded as uninformative in the prior figure planning and is not revived to fill a panel. Directional entropy and dominant-direction mass were excluded because S2/S3 and the main geometry figure already cover them. S6 therefore focuses on complementary dimensions of residual geometry.
