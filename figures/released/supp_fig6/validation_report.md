# S6 validation report

## Inputs and provenance

- Figure 4 manifest: `figures/fig4/data/fig4_data_manifest.json`, status `PASS`; SHA-256 recorded in `data/supp_fig6_build_manifest.json`.
- Frozen inputs: the four `figures/fig4/data/fig4_geometry_paired/algorithm=*/part.parquet.xz` partitions. All four SHA-256 values are checked against the Figure 4 manifest and recorded in the S6 build manifest.
- The paired table's source lineage is frozen PL1 A1 geometry and A21 A2 geometry. The builder uses only the persisted paired table and does not rerun production analyses.

## Dimensions and keys

- Input and compact geometry output: 1,671,600 rows; 400 evaluations × 4,179 reactions.
- Method counts: 417,900 each for CORDA, GIMME, iMAT, and RIPTiDe.
- Pairing key: `(evaluation_id, reaction_id)`; no missing or duplicate keys. Output uniqueness is checked on `(algorithm, evaluation_id, reaction_id)`.
- Evaluation metadata (`algorithm`, RNA context, CT2A mouse, GL261 mouse) are internally consistent.

## Metric checks

- **η²:** production definition is weighted sign-category explained variance for `|Δv_B|`; valid range [0, 1]. Finite A1: 580,507; finite A2: 650,226; paired finite: 573,882; paired comparison excludes 1,097,718 rows with at least one nonfinite arm. Per-method paired counts: CORDA 212,164; GIMME 141,365; iMAT 207,930; RIPTiDe 12,423. Degenerate values remain missing; they are not imputed.
- **Non-tie coverage:** `1 − p_tie`, with frozen production tie tolerance 1e-12; valid range [0, 1]. All 1,671,600 pairs are finite. Values within numerical tolerance of the range boundary are retained without clipping.
- **Supported sign states:** count of sign categories {-1, 0, +1} with positive mass; valid categories {1, 2, 3}. All 1,671,600 pairs are finite and integer-valued. Transition counts are: 1→1 1,293,850; 1→2 1,055; 2→2 376,470; 3→3 225; all other transitions 0.
- Paired differences are computed explicitly as A2 − A1. No stored delta column is substituted.

## Exclusions and nonduplication

- η² paired comparisons exclude 1,097,718 rows because at least one anchor value is nonfinite from a degenerate magnitude distribution.
- No rows are excluded from non-tie coverage or sign-state counts.
- Joint ESS (S1 overlap), width80 (not to be revived for figure fill), directional entropy (S2/S3 overlap), and dominant-direction mass (S2/S3 overlap) are omitted.
- No truth-pair, utility, alternative-anchor, exploratory, or newly generated production data enter S6.

## Scientific observations

Across the frozen population, median paired changes are zero. Mean η² change is 0.001429 over its 573,882 finite pairs; mean non-tie coverage change is −0.00003375 over all pairs. Supported-state count is unchanged in 99.9369% of rows. These are descriptive reaction–evaluation summaries, not independent biological replicates or evidence of a uniform method effect.

## Issues

- **BLOCKING:** none identified in table construction or metric checks.
- **IMPORTANT:** η² is unavailable for many pairs due to the production degeneracy rule, so its panel denominator is substantially smaller and method-specific; panel counts are shown.
- **COSMETIC:** none identified from data validation. Visual inspection is required after notebook rendering.
