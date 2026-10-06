# S7 validation report

## Inputs and provenance

All inputs are repository-contained frozen Figure 4C production artifacts. The builder verifies generated-table hashes against `figures/fig4/data/fig4_data_manifest.json`; verifies A1 weights, A2 weights, candidate inventory, and flux-cache hashes against the same manifest; and also verifies A2 weights against the A20 manifest. Exact hashes are recorded in `data/supp_fig7_manifest.json`.

The authoritative candidate-selection inputs are `fig4_panelC_candidate_groups.parquet.xz` (100 groups) and its Figure 4 manifest. The authoritative weighted distributions are in `fig4_panelC_candidate_support.parquet.xz` (candidate-state delta, product weight, sign, arm, method, evaluation, and support ID). `fig4_panelC_candidate_metrics.parquet.xz` provides the frozen method/evaluation/context mapping and A1/A2 geometry summaries. The source distribution implementation is `figures/fig4/build_fig4_tables.py::mats_delta`, reproduced in the existing Figure 4C notebook.

No reconstruction, optimization, sampling, candidate generation, or production pipeline was rerun.

## Selection provenance

FACOAL204 is retained because it is already used in Main Figure 4C. The notebook calls it an “outcome-blind candidate,” and the frozen group table gives it geometry rank 23. The repository records a general geometry ranking rule but no specific frozen rule that explains why FACOAL204 was selected for the main panel. S7 does not infer one.

The two additions are descriptive post hoc illustrations selected before loading support distributions or plotting: excluding FACOAL204, the builder chooses the maximum `median_absolute_A2_minus_A1_H_dir_shift`, then the minimum of that metric among remaining groups. Ties use ascending `geometry_rank`, then `candidate_group_id`. There is no claim of statistical representativeness. A zero directional-entropy shift does not imply that the full flux-magnitude distributions are unchanged.

## Dimensions and integrity

The selected examples are FACOAL204 (rank 23; median absolute entropy shift 0.5411), AMPDA (rank 63; maximum among the other 99 candidates, 0.2261), and RE3597C (rank 1; zero median absolute entropy shift). AMPDA and RE3597C are chosen mechanically by the predeclared summary rule. The builder records the exact selected reaction, context, CT2A/GL261 mouse IDs, method-specific evaluation and truth-pair IDs, source anchor arm, selection reason, metric value, and geometry rank in `data/supp_fig7_example_selection.tsv`.

The input candidate pool has 100 groups. The selected output contains 12 method/evaluation records and 9,600 support rows, arranged as 24 example × method × arm distributions with 400 unique support states each. Effective support ranges from 1.84 to 152.18 across those distributions.

It requires three distinct groups, four methods per group, two anchor settings, 400 unique support states per example/method/anchor, matching A1/A2 Cartesian support IDs, finite `delta_v_B` and weights, nonnegative weights, normalized weight mass within `1e-10`, and sign labels consistent with the frozen `1e-12` tie tolerance. It checks support evaluation/reaction/context identities against the frozen metric table. Duplicate support keys: 0; missing/nonfinite delta or weight: 0; negative weights: 0; maximum absolute weight-sum error: `2.22e-16`. The notebook checks weighted histogram mass (maximum error `2.22e-16`) and emits both SVG and PNG. The manifest records exact source and output hashes.

The exact row counts, weight-sum error, effective-support range, hashes, and output checksums are recorded in `data/supp_fig7_manifest.json`. Any exclusion from the frozen group pool is counted there with its reason.

## Scientific definitions

- A1 is 1-Strong-Anchor (glucose uptake).
- A2 is 2-Strong-Anchors (glucose uptake + LDH_L); the production source arm is `A2-L`.
- The production `delta_v_B` convention is CT2A flux minus GL261 flux, as implemented in the Figure 4C `mats_delta` function.
- Flux units are mmol gDW⁻¹ h⁻¹. Within each example/method panel, A1 and A2 use the same histogram bins and the frozen support weights.

## Issues

- **BLOCKING:** None if all automated provenance and integrity checks pass.
- **IMPORTANT:** The main-panel FACOAL204 selection rationale is not specified beyond the “outcome-blind candidate” note; S7 additions are explicitly descriptive post hoc illustrations. Weighted support effective sample size can be concentrated and does not count independent observations. X-axis limits vary by example/method panel, with A1/A2 comparisons sharing limits and bins within each panel.
- **COSMETIC:** None identified after visual inspection of the rendered outputs.
