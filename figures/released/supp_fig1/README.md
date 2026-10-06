# Supplementary Figure S1

## Purpose

S1 documents how finite candidate-state weights behave across reconstruction contexts. It separates global 320-candidate ESS targeting, conditional 20-candidate stratum ESS, and product-weight support for contrast distributions. These are diagnostics of finite weighted supports, not counts of statistically independent observations.

## Frozen inputs

The figure uses the A1 `strong_anchor_baseline` weights from `outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_FOUR_ARM_WEIGHTS.tsv.xz` and the A2 `A2-L` weights from `outputs/dmi_bridge_a20_dual_anchor_qualification_v1/BRIDGEA20_DUAL_ANCHOR_WEIGHTS.tsv.xz`. Publication-facing labels call these A1 (glucose uptake) and A2 (glucose uptake + LDH_L). The A1 artifact hash is checked against the frozen Figure 4 data manifest; A20, PL1, PL2A, PL2B, and A22 artifacts are verified against their own manifests. Truth support is restricted to the fixed 836 matched A1/A2 pairs.

No exploratory alternative-anchor outputs, reconstruction, optimization, candidate generation, or sampling are used.

## Rebuild

From the repository root:

```bash
python figures/supp_fig1/build_supp_fig1_tables.py
jupyter nbconvert --to notebook --execute --inplace figures/supp_fig1/supp_fig1.ipynb
```

The table builder writes compact plotting tables and `data/validation_summary.json`. The notebook writes `outputs/supp_fig1.svg` and `outputs/supp_fig1.png`.

## Panels

- **A:** Conditional ESS across all 320 method × RNA context × tumor × mouse strata, split by tumor and anchor. Each stratum contains 20 candidate states. RNA contexts are numbered in the frozen source order; method colors and group order are shown in the panel.
- **B:** Heatmaps of the normalized 20 candidate weights for every stratum. Every stratum is shown, including concentrated supports.
- **C:** Global ESS across each 320-candidate tumor/mouse pool; product-weight ESS over each 400-state contrast support; post-holdout, condition-specific ESS for the fixed 836 matched truth pairs.

## Sample sizes

There are 6,400 candidate-weight rows and 320 conditional strata (160 per anchor setting). Global ESS is summarized in 20 tumor × mouse × anchor pools. Product-weight contrast support has 800 evaluation × anchor rows. Post-holdout diagnostics cover 836 truth pairs × two conditions × two anchors (3,344 rows).

## Limitations

The conditional weight heatmap preserves actual candidate-level weights but displays candidates in their frozen index order. The global ESS target is approximately 20; the exact frozen weights yield values numerically equal to 20. Product-weight ESS is the effective support of the Cartesian product distribution and is distinct from both global and conditional ESS. These values do not establish statistical independence.
