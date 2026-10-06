# Supplementary Figure S3 validation report

## Inputs and provenance

S3 uses the frozen paired geometry table produced by FIG4-DATA-B, not exploratory alternative-anchor outputs.

| Input | Role and provenance |
|---|---|
| figures/fig4/data/fig4_data_manifest.json | Figure 4 data manifest; stage FIG4-DATA-B, version 1.1.0, status PASS; recorded geometry population 1,671,600. |
| figures/fig4/data/fig4_geometry_paired/algorithm=CORDA/part.parquet.xz | Frozen paired geometry partition; SHA-256 0fd4d308c6d03d5ef144fb5e8c48e497e115af7bc27cdddd1010ddef4e3ba457. |
| figures/fig4/data/fig4_geometry_paired/algorithm=GIMME/part.parquet.xz | Frozen paired geometry partition; SHA-256 f492947266440ed378540d6c0d6a40aabb8952fafef5cc728dbc8c9aa70f17ad. |
| figures/fig4/data/fig4_geometry_paired/algorithm=iMAT/part.parquet.xz | Frozen paired geometry partition; SHA-256 a44c629ffdc5ddcb1956d6132a8ec15d7dbb358e86b7888b512e315a73d4a070. |
| figures/fig4/data/fig4_geometry_paired/algorithm=RIPTiDe/part.parquet.xz | Frozen paired geometry partition; SHA-256 642fc972db1aa17854fbd631afb909fc07e9be54c7d5fdc51d23fbb60ffc2e42. |

The Figure 4 manifest maps A1 geometry to the frozen PL1 output and A2 geometry to the frozen A21 output. Their source statuses are PL1_PREDICTABILITY_LANDSCAPE_COMPLETE and BRIDGE_A21_DUAL_ANCHOR_GEOMETRY_FROZEN. In publication-facing terminology, A1 is 1-Strong-Anchor (glucose uptake) and A2 is 2-Strong-Anchors (glucose uptake + LDH_L).

## Dimensions and pairing

- Input geometry rows: **1,671,600**.
- Output paired-table rows: **1,671,600**.
- Methods: **4**, with **417,900** rows per method.
- Unique evaluations: **400**; unique reactions: **4,179**.
- Exact S3 pair key: (evaluation_id, reaction_id).
- Pair metadata retained: method, RNA context, CT2A mouse, and GL261 mouse.
- Duplicate pair keys: **0**.
- Matched A1 rows / matched A2 rows: **1,671,600 / 1,671,600**.
- Unmatched A1 / unmatched A2: **0 / 0**, as verified by the Figure 4 source builder's one-to-one outer merge and all-keys-matched check; the S3 builder also verifies the frozen table's pair-key count and uniqueness.
- Metadata consistency across each pair: **PASS**.

## Panel-specific metric checks

Both panels use A1 as x and A2 as y. The panel-specific sample size is **n = 1,671,600** for each; all input rows have finite values in both arms, and **0** rows are excluded.

| Metric | A1 observed range | A2 observed range | Values outside [0,1] beyond 1e-12 | Values just outside exact [0,1] due to floating-point roundoff |
|---|---:|---:|---:|---:|
| Directional entropy \(H_{\mathrm{dir}}\) | \([-2.02\times10^{-16}, 0.758945]\) | \([-2.02\times10^{-16}, 0.845123]\) | 0 | 8,445 across both arms |
| Dominant-direction mass | \([0.500000000000952, 1.0000000000000002]\) | \([0.500015314047320, 1.0000000000000002]\) | 0 | 156,943 across both arms |

The entropy definition in the production code is \(-\sum_s p_s\log(p_s)/\log(3)\) over nonzero negative, tie, and positive sign masses. Dominant-direction mass is the maximum of those three probabilities. All exact-bound excursions are at machine precision; the maximum violation beyond the unit interval is below 1e-12. The plotting limits include that roundoff, and no values are clipped.

Fractions relative to the identity line use a tolerance of \(10^{-12}\):

| Metric | A2 > A1 | A2 approximately A1 | A2 < A1 |
|---|---:|---:|---:|
| Directional entropy \(H_{\mathrm{dir}}\) | 17.19% | 78.87% | 3.94% |
| Dominant-direction mass | 3.91% | 79.15% | 16.94% |

These pooled descriptive fractions show the direction of changed pairs while retaining the many unchanged observations. They do not imply a uniform shift across reactions or methods.

## Exclusions and population checks

- Nonfinite A1/A2 pairs: **0** excluded for either metric.
- Duplicate pairs: **0**.
- Methods match the expected four-method panel; paired metadata are consistent, and A1/A2 metric columns follow the Figure 4 manifest mapping.
- Development/confirmation mixing: **none**; the inputs are the Figure 4 frozen paired geometry partitions only.
- Exploratory alternative-anchor data: **none**.
- The common population is the full reaction × evaluation geometry population, not the 836 truth-pair utility population.

## Scientific checks

- A1/A2 labels follow the frozen definitions: glucose uptake and glucose uptake + LDH_L.
- Each observation is paired by the production key; no unpaired comparison or aggregation across reactions/evaluations was used.
- The identity line means equal descriptor values under A1 and A2. Values above it indicate a larger A2 value.
- Entropy and probability mass are dimensionless. Flux units do not apply.
- The density representation is appropriate for 1.67 million observations and preserves observation-level pairing in the underlying plotted table.
- The pooled method display is intentional: S2 presents the method-stratified paired-difference analysis.

## Current Figure 4 assembly and open issues

The direct pooled scatter code remains in figures/fig4/fig4_panelB_distributions_updated.ipynb, but it is a legacy supplementary candidate, not part of the current assembled main Figure 4. The notebook's assembly cell combines Figure 4A, the absolute distributions for Figure 4B, and the FACOAL204 weighted-distribution example for Figure 4C. The current draft manuscript caption likewise defines Figure 4B as absolute distributions and Figure 4C as the FACOAL204 example, while referring to the paired A1-versus-A2 comparison as Figure S6. The authoritative S1–S3 contract and dedicated S3 task assign this paired geometry comparison to S3.

- **BLOCKING:** None.
- **IMPORTANT:** The manuscript draft's Figure S6 cross-reference conflicts with the authoritative S3 assignment and should be reconciled by the manuscript/controller owner.
- **COSMETIC:** None identified.
