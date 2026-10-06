# Supplementary Figure S2

## Purpose

S2 shows paired changes in residual flux-distribution geometry from A1 to A2. Each panel is an empirical cumulative distribution of **A2 − A1** shifts, with one curve per reconstruction method. The zero line marks no shift. This method-stratified view extends the pooled Figure 4 result and shows whether its direction is shared across methods.

The anchor definitions are A1 = **1-Strong-Anchor (glucose uptake)** and A2 = **2-Strong-Anchors (glucose uptake + LDH_L)**. Both plotted descriptors are dimensionless.

## Inputs and provenance

The table builder reads only the manifest-verified Figure 4 paired geometry partitions:

- `figures/fig4/data/fig4_geometry_paired/algorithm=iMAT/part.parquet.xz`
- `figures/fig4/data/fig4_geometry_paired/algorithm=GIMME/part.parquet.xz`
- `figures/fig4/data/fig4_geometry_paired/algorithm=CORDA/part.parquet.xz`
- `figures/fig4/data/fig4_geometry_paired/algorithm=RIPTiDe/part.parquet.xz`
- Provenance and digest source: `figures/fig4/data/fig4_data_manifest.json`

The Figure 4 data manifest reports `PASS`, records the source lineage and generated hashes, and reports 1,671,600 uniquely keyed A1/A2 geometry rows. The builder checks each input partition's SHA-256 against that manifest before reading it. Figure 4's paired geometry builder establishes the A1/A2 match using `evaluation_id` and `reaction_id` with a one-to-one join and records a passing duplicate-key check. Its source feature tables come from the frozen PL1 A1 and A21 A2 geometry outputs.

No truth-pair utility data, alternative-anchor tables, or exploratory flux cache enter S2. The 836 matched truth pairs do not define this figure's population.

## Tables and plotting

Build and validate the compact tables from the repository root:

```bash
python figures/supp_fig2/build_supp_fig2_tables.py
```

The builder validates source hashes, method membership, row/key counts, finite paired values, metric ranges, and the stored subtraction against A2 minus A1. It writes:

- `data/supp_fig2_H_dir_paired_shifts.tsv.gz`
- `data/supp_fig2_dominant_direction_mass_paired_shifts.tsv.gz`
- `data/supp_fig2_table_summary.json`

Execute `supp_fig2.ipynb` from the repository root to create:

- `outputs/supp_fig2.svg`
- `outputs/supp_fig2.png` (300 dpi)

## Panels and sample sizes

- **Directional entropy:** ECDF of paired changes in (H_{dir}).
- **Dominant-direction mass:** ECDF of paired changes in probability mass assigned to the dominant direction.
- Both panels show iMAT, GIMME, CORDA, and RIPTiDe separately, include a dashed zero reference, and report (n) in each method legend entry.
- Each method contributes 417,900 paired reaction–evaluation observations to each panel; total (n=1,671,600) per panel. There are no nonfinite-value exclusions.

## Limitations

The rows describe paired reaction–evaluation geometry observations, not independent biological replicates. The descriptors summarize finite weighted candidate distributions. Mean shifts are modest relative to the full range, and each method has a median shift of zero; interpret the ECDFs as distributions of paired changes rather than a uniform shift in every reaction.
