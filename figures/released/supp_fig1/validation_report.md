# S1 validation report

## Inputs

All inputs are repository-contained frozen production outputs. Source manifests and hashes are recorded in `data/validation_summary.json`.

| Input | Role and provenance |
|---|---|
| `outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_FOUR_ARM_WEIGHTS.tsv.xz` | A1 `strong_anchor_baseline`; its exact SHA256 is cross-checked against `figures/fig4/data/fig4_data_manifest.json`, which records it as a Figure 4 source. Other BIO0 arms are not used. |
| `outputs/dmi_bridge_a20_dual_anchor_qualification_v1/BRIDGEA20_DUAL_ANCHOR_WEIGHTS.tsv.xz` | Qualified frozen A2-L weights, displayed as A2 per the master contract; verified against `BRIDGEA20_MANIFEST.json`. |
| `outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_EVALUATION_REGISTRY.tsv` | Frozen 400-evaluation population for contrast-product support; verified against the PL1 manifest. |
| `outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_TRUTH_SELECTION.tsv` and `BRIDGEPL2A_TRUTH_PAIR_REGISTRY.tsv` | Frozen holdout candidate support and truth-pair population; verified against the PL2A manifest. |
| `outputs/dmi_bridge_pl2b_sign_only_utility_v1/BRIDGEPL2B_DISTINCT_TRUTH_PAIR_REGISTRY.tsv` | A1 post-holdout support paired to the fixed truth population; verified against the PL2B manifest. |
| `outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_MATCHED_TRUTH_PAIR_REGISTRY.tsv` and `BRIDGEA22_TRUTH_SUPPORT_AUDIT.tsv` | A2 post-holdout support for the same 836 pairs; verified against the A22 manifest. |

No exploratory alternative-anchor evidence is used. No reconstruction, optimization, sampling, or candidate generation was rerun.

## Dimensions and sample sizes

- BIO0 four-arm input: 12,800 rows; selected A1 strong-anchor baseline: 3,200 rows.
- A2 frozen weights: 3,200 rows.
- Figure candidate-weight table: 6,400 rows; 320 unique strata; 20 candidates per stratum.
- Strata: 4 methods × 4 RNA contexts × 2 tumors × 5 mice × 2 anchors = 320.
- Global ESS: 20 pools, each covering 320 candidates.
- Contrast product support: 800 rows (400 evaluations × two anchor settings), each representing 400 candidate pairs.
- Fixed matched truth population: 836 unique truth pairs.
- Post-holdout pair support: 3,344 rows (836 pairs × two conditions × two anchors).
- PL2A candidate-selection diagnostics: 480 rows (retained in the plotting tables; they do not define the fixed 836-pair comparison denominator).

## Integrity and scientific checks

- Expected methods (`iMAT`, `GIMME`, `CORDA`, `RIPTiDe`), four RNA contexts, two tumors, A1 and A2: PASS.
- Candidate count per conditional stratum (20), stratum count (320), global pool count (320): PASS.
- Duplicate candidate keys, missing/nonfinite weights, negative weights: 0.
- Global weight normalization and conditional within-stratum renormalization: PASS (absolute tolerance 1e-10).
- Recomputed conditional ESS uses `1 / sum(w²)` after within-stratum normalization. Range: 1.00–16.24.
- Global ESS uses the frozen across-context/method weights over 320 candidates. Range: 20.000000–20.000000, matching the approximately 20 target.
- Product-weight contrast ESS uses the Cartesian product of two condition-specific 20-candidate distributions. Range: 1.00–201.03 across 800 evaluation/anchor rows.
- Post-holdout support pair matching: PASS; exactly 836 matched pairs, each with CT2A and GL261 support for both anchors.
- Post-holdout conditional ESS range: 1.00–10.30 (minor floating-point deviation below one is within numerical tolerance).
- Highly concentrated cases remain present: 153 strata have conditional ESS ≤ 2; 181 have maximum candidate weight ≥ 0.5; no stratum has singleton positive support.
- Publication-facing terminology: A1 = 1-Strong-Anchor (glucose uptake); A2 = 2-Strong-Anchors (glucose uptake + LDH_L). ESS and probabilities are dimensionless; these panels do not display flux quantities.

## Exclusions

- A1 `strong_anchor_baseline` is selected from the 12,800-row BIO0 table; 9,600 rows from the three other decomposition arms are excluded because they are not the A1 baseline weights.
- PL2A includes 1,200 initial truth pairs, including 90 non-evaluable due to zero post-holdout mass. These do not enter the S1 fixed comparison population. The matched registry contains 836 pairs; no additional pairs are dropped from that population.
- No rows are removed for concentration. The plotted strata and weights include all 320 conditional strata.

## Open issues

- **BLOCKING:** None.
- **IMPORTANT:** Effective support is a finite-weight diagnostic and must not be interpreted as an independent sample count. The figure and caption state this explicitly.
- **COSMETIC:** None identified after visual inspection.
