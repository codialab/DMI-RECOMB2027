# Supplementary Figure S7

## Objective

Show weighted residual `P(Δv_B)` distributions for the existing Main Figure 4C example and two examples selected descriptively from the same frozen Figure 4C candidate pool. The additional examples illustrate the largest and smallest frozen median absolute A1-to-A2 directional-entropy shifts in that pool after excluding FACOAL204. A zero directional-entropy shift does not imply that the full flux-magnitude distributions are unchanged. The examples are not statistically representative.

## Authoritative inputs

- `figures/fig4/data/fig4_data_manifest.json` and its checksum-verified `fig4_panelC_candidate_groups.parquet.xz`, `fig4_panelC_candidate_metrics.parquet.xz`, and `fig4_panelC_candidate_support.parquet.xz`.
- The Figure 4C notebook and builder document the distribution implementation and source provenance. BIO0 A1 weights, A20 A2 weights, the candidate flux cache, and the Stage-12 candidate inventory are checked against hashes in the Figure 4 manifest; the A20 weights are also checked against the A20 manifest.

## Rebuild

From the repository root:

```bash
python figures/supp_fig7/build_supp_fig7_tables.py
jupyter nbconvert --to notebook --execute --inplace figures/supp_fig7/supp_fig7.ipynb
```

The builder writes compact weighted-distribution, example-selection, and support-summary tables under `data/`. The notebook renders `outputs/supp_fig7.svg` and `outputs/supp_fig7.png`.

## Panels and aggregation

The figure is a 3-row by 4-column grid: examples by reconstruction method. Each panel overlays A1 (glucose uptake) and A2 (glucose uptake + LDH_L) weighted densities using common bins within that example/method. Axis ranges vary between panels. Each method/anchor distribution contains 400 Cartesian support states. The production support's `delta_v_B` is CT2A flux minus GL261 flux, in mmol gDW⁻¹ h⁻¹; the frozen `weight` column is the product weight over the two mouse-specific candidate supports.

FACOAL204 is retained as the existing Main Figure 4C example. The other examples are chosen from the 100 frozen candidate groups by descending and ascending `median_absolute_A2_minus_A1_H_dir_shift`; ties resolve by ascending `geometry_rank`, then `candidate_group_id`. FACOAL204 is excluded from these two additional choices. Selection is based on the frozen summary table before distribution support is loaded or plotted.

## Known limitations

The repository identifies FACOAL204 as an outcome-blind Main Figure 4C candidate and provides a general candidate ranking, but it does not record a frozen rule specific to why FACOAL204 was selected from the candidate pool. The two S7 additions are descriptive post hoc examples chosen by an explicit S7 rule. They do not establish statistical representativeness. Effective support can be concentrated; the builder reports it in the support summary for QC, not as an independent sample count.
