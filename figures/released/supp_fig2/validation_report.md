# S2 Validation Report

## Inputs

Frozen production inputs are the four Figure 4 A1/A2 paired geometry partitions listed below. `figures/fig4/data/fig4_data_manifest.json` reports `status: PASS`, identifies its frozen upstream source lineage, records source and generated artifact digests, and records 1,671,600 matched geometry rows with a passing unique-key check. The S2 builder verifies the manifest SHA-256 for each partition before use.

| Input artifact | SHA-256 | Status |
|---|---|---|
| `figures/fig4/data/fig4_geometry_paired/algorithm=iMAT/part.parquet.xz` | `a44c629ffdc5ddcb1956d6132a8ec15d7dbb358e86b7888b512e315a73d4a070` | Frozen Figure 4 geometry partition; hash verified |
| `figures/fig4/data/fig4_geometry_paired/algorithm=GIMME/part.parquet.xz` | `f492947266440ed378540d6c0d6a40aabb8952fafef5cc728dbc8c9aa70f17ad` | Frozen Figure 4 geometry partition; hash verified |
| `figures/fig4/data/fig4_geometry_paired/algorithm=CORDA/part.parquet.xz` | `0fd4d308c6d03d5ef144fb5e8c48e497e115af7bc27cdddd1010ddef4e3ba457` | Frozen Figure 4 geometry partition; hash verified |
| `figures/fig4/data/fig4_geometry_paired/algorithm=RIPTiDe/part.parquet.xz` | `642fc972db1aa17854fbd631afb909fc07e9be54c7d5fdc51d23fbb60ffc2e42` | Frozen Figure 4 geometry partition; hash verified |

The Figure 4 geometry builder uses frozen PL1 feature outputs for A1 and frozen A21 feature outputs for A2. The Figure 4 build contract matches on `evaluation_id` and `reaction_id`, requires one-to-one matching, checks matching evaluation metadata, and removes the quantitative anchor reaction from the non-anchor reaction inventory. The S2 population is that resulting paired geometry table. S2 does not use the 836 truth-pair population, utility tables, or exploratory alternative-anchor outputs.

Anchor definitions: A1 = 1-Strong-Anchor (glucose uptake); A2 = 2-Strong-Anchors (glucose uptake + LDH_L). Plotted descriptors and their shifts are dimensionless.

## Dimensions and pairing integrity

| Check | Result |
|---|---:|
| Figure 4 paired geometry rows in manifest | 1,671,600 |
| Unique `evaluation_id` × `reaction_id` keys | 1,671,600 |
| Methods | iMAT, GIMME, CORDA, RIPTiDe |
| Rows per method | 417,900 |
| S2 output rows per metric | 1,671,600 |
| Duplicate keys in S2 method tables | 0 |
| Unmatched A1/A2 pairs | 0 (upstream one-to-one Figure 4 pairing check: PASS) |
| Missing pairing keys | 0 |
| Missing/nonfinite values for either S2 metric | 0 |

Each paired shift was checked against `A2_metric - A1_metric` at absolute tolerance (5\times10^{-12}). No sorting, sampling, or finite-value filtering changed the population.

## Panel-specific sample sizes and observed shifts

| Metric | Method | Original paired (n) | Nonfinite excluded | Plotted (n) | Mean A2 − A1 | Median A2 − A1 |
|---|---|---:|---:|---:|---:|---:|
| Directional entropy (H_{dir}) | iMAT | 417,900 | 0 | 417,900 | 0.02224 | 0 |
| Directional entropy (H_{dir}) | GIMME | 417,900 | 0 | 417,900 | 0.03469 | 0 |
| Directional entropy (H_{dir}) | CORDA | 417,900 | 0 | 417,900 | 0.04872 | 0 |
| Directional entropy (H_{dir}) | RIPTiDe | 417,900 | 0 | 417,900 | 0.00068 | 0 |
| Dominant-direction mass | iMAT | 417,900 | 0 | 417,900 | −0.01069 | 0 |
| Dominant-direction mass | GIMME | 417,900 | 0 | 417,900 | −0.01765 | 0 |
| Dominant-direction mass | CORDA | 417,900 | 0 | 417,900 | −0.02485 | 0 |
| Dominant-direction mass | RIPTiDe | 417,900 | 0 | 417,900 | −0.00002 | 0 |

The mean shifts have the expected signs for all four methods: entropy increases and dominant-direction mass decreases. Medians are zero for all methods, indicating the shift is not uniform across paired observations.

## Integrity and scientific checks

- Anchor values were drawn from the manifest-verified A1/A2 paired geometry table; publication-facing terminology is used in the plot and documentation.
- Both entropy and dominant-direction mass lie in the valid probability/entropy range ([0,1]) within (5\times10^{-12}) numerical tolerance. The source includes values within floating-point tolerance of zero and one; values are retained as stored.
- Stored deltas equal A2 minus A1 within tolerance for every plotted row.
- The expected four methods are present; each metric includes the same paired reaction–evaluation population.
- No development/confirmation populations are mixed. No exploratory alternative-anchor data entered the build.
- ECDF curves include all finite shifts and zero is marked with a dashed vertical line. (n) is annotated per method and panel.
- No exclusion rule was applied; excluded observations: zero for every method and panel.

## Outputs

- Plotting tables and machine-readable build summary: `data/`.
- Plotting notebook: `supp_fig2.ipynb`.
- Rendered outputs: `outputs/supp_fig2.svg` and `outputs/supp_fig2.png`.
- Rebuild command: `python figures/supp_fig2/build_supp_fig2_tables.py`, followed by execution of the notebook from the repository root.

## Open issues

- **BLOCKING:** None.
- **IMPORTANT:** None.
- **COSMETIC:** None.
