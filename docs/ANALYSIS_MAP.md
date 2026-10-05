# Manuscript analysis and figure map

Authority: the current main and supplementary DOCX drafts at source commit
`16379a169f45ec556274953eb66ae12293e1e985`, interpreted with the manuscript technical
audit. Draft supplement panels are not finalized. Existing folder numbering is
supporting evidence and does not override the drafts.

[`FIGURE_REGISTRY.json`](../manifests/FIGURE_REGISTRY.json) records executable
readiness for every figure/panel group. [`FIGURE_INPUTS.json`](../manifests/FIGURE_INPUTS.json)
lists exact prospective numerical paths and hashes. Those files are currently
withheld for rights review. `reproduce.py` reads only destination files.

## Main figures

| Figure/panel | Plotting code or artwork | Frozen inputs / upstream code | Expected output and manuscript metric |
|---|---|---|---|
| 1A | Conceptual artwork only; source `figures/fig1/recomb_figure1_final_download.svg`, withheld for artwork/rights review | No numerical dependency | DMI observables mapped to a few GEM coordinates; most fluxes remain unmeasured |
| 1B | Same conceptual SVG | No numerical dependency; illustrative densities | Anchoring may change a distribution while preserving a projected feasible range |
| 1C | Same conceptual SVG | No numerical dependency | Separates feasibility, feasible ranges and sampled distributions; motivates two questions |
| 2A–B | Conceptual artwork; source `recomb_figure2_ABC_rough_design_b.svg` requires editorial/rights review | No numerical dependency | Strong anchors plus a separate weak sign-only cue; direction–magnitude coupling |
| 2C | Same conceptual source | No quantitative plot dependency; context design audited in Stage-11 configuration | Four methods and four RNA specifications; development/holdout distinction |
| 2D | `figures/main/fig2/plot_fig2_dg.py` | `data/figure_inputs/fig2/figure2_DG_points.csv.xz`; upstream export builder in `figures/main/fig2/upstream/`; PL1 + PL2C + A23 | `reproduced/figures/fig2DG/figure2_DG_usefulness_scatter.{svg,png}`; A1 development usefulness frequency, n=2,454, rho≈0.857 |
| 2E | Same plotting script and input | A21 + A22 + A23, same frozen development split | Same four-panel output; A2 development, n=2,485, rho≈0.834 |
| 2F | Same plotting script and input | PL2D original one-shot confirmation | Same output; A1 held-out, n=838, rho≈0.835 |
| 2G | Same plotting script and input | A21 + A22 + A23, frozen held-out split | Same output; A2 robustness, n=843, rho≈0.849; not independent confirmation |
| 3A | Conceptual source `figure3_panelA_directional_information_gain_v2b.svg`; artwork/rights pending | No empirical data dependency | Baseline and sign-informed magnitude errors; illustrative gain definition |
| 3B | `figures/main/fig3/plot_fig3_bc.py`, extracted source cell 20 | `data/figure_inputs/fig3/fig3_gain_reaction_evaluation.parquet.xz`; upstream `make_fig3_gain_tables.py`, PL2B + A22 + A23 | `reproduced/figures/fig3BC/fig3_panelBC.{svg,png}`; correct-cue pair-mean gain exceedance, 660,237 pairs per anchor |
| 3C | Same script and table | Same fixed population | Same output; full correct-minus-wrong pair-mean gain exceedance, not half-separation information advantage |
| 3D | No authoritative corrected plotting source identified | `fig3_gain_summary.parquet.xz` and `fig3_directional_advantage_summary.parquet.xz` in the same prospective input folder | Final endpoint means are traceable; legacy D artwork contains superseded case-pooled values and is excluded |
| 4A | Conceptual source `figure4_panelA_draft_b.svg`; artwork/rights pending | No empirical density dependency; labels describe the A23 matched design | Schematic transition from HEX1 to HEX1 + LDH_L |
| 4B | `figures/main/fig4/plot_fig4_b.py`, source cells 3 and 6 | Four `data/figure_inputs/fig4/fig4_geometry_paired/algorithm=<method>/part.parquet.xz` partitions plus `fig4_data_manifest.json`; updated builder, PL1 + A21 | `reproduced/figures/fig4B/fig4B_metric_distributions.{svg,pdf}`; absolute directional-entropy and dominant-direction-mass ECDFs, 1,671,600 matched rows |
| 4C | `figures/main/fig4/plot_fig4_c_source.py`, source cell 11; blocked for main-assembly/selection review | `fig4_panelC_candidate_support.parquet.xz`, candidate groups/metrics and source manifest; updated Figure 4 builder | Source four-method drawing `fig4C_FACOAL204_distributions.{svg,pdf}`; final main draft requires three methods and moves RIPTiDe to S7 |

The Figure 2 source notebook contains both information-advantage and usefulness
variants. Only the usefulness plotting cells matching the draft are selected.
The display label was changed to “Benchmark directional-usefulness frequency” to
avoid implying calibrated probability. The endpoint, aggregation and parity checks
are unchanged.

Figure 3B/C preserve the original plotting grid `np.logspace(-15, 0, 250)`. This
is a display grid, not evidence of a prespecified nonzero scientific cutoff.
The source technical audit establishes native/model flux units. Final artwork
must not claim an unrecorded absolute DMI-to-GEM conversion.

Figure 4C is FACOAL204, context x1+x3, C1 versus G1, group
`FCG_ed8b1873c8dbc03fe5eb`. The original notebook retains four methods. Its physical
unit label is adapted to native GEM flux units; weights, histogram definition,
supported sign labels and entropy calculation remain unchanged. Rendering it as
final main artwork is blocked until author review reconciles composition and
selection intent.

## Supplementary figures

These entries follow the **current draft**, not the legacy folder labels. No
additional panel letters or scientific comparisons have been invented.

| Draft figure | Code/input evidence | Required scientific message and unresolved coverage |
|---|---|---|
| S1 | `figures/supplementary/s1/plot_s1.py`; prospective `data/figure_inputs/s1/` compact weight, stratum, global/product/post-holdout support tables; source builder retained under `upstream/` | Candidate weights and distinct ESS definitions on the fixed 836-pair population. Scope matches; source heatmap row/label order needs review. Rendering blocked. |
| S2 | PL2D context results and M1 confirmation-context evidence; legacy S4 plotting cells in provenance contain context plots | Draft requests context-specific A1 held-out associations. Legacy S2 is a paired-geometry shift plot and is excluded from current production. Matching final composition/code is missing. |
| S3 | PL2C/PL2D supportive results and M1 facts identify alternative descriptors and the continuous information-advantage endpoint | Draft requests descriptor comparisons and separate continuous information advantage. Legacy S3 is a paired geometry hexbin. Final selections, plotting code and panel denominators need author review. |
| S4 | Legacy `supp_fig4` builder/plot cells retained as provenance; PL2D bootstrap/supportive results provide upstream evidence | Bootstrap and separate development/confirmation/context plots exist. Draft also requires coverage-adjusted and reaction-exclusion analyses; current workflow does not cover all requested panels. |
| S5 | Legacy `supp_fig5` builder/plot cells retained as provenance; final Figure 3 pair tables and PL2B/A22 control tables | Existing context gains/outcome composition partially match. Draft requests negative tails and reliability curves with explicit case/pair weighting. Missing portions must not be replaced by an exploratory scan. |
| S6 | Legacy `supp_fig6` builder/plot cells retained as provenance; Figure 4 geometry partitions, PL1/A21 descriptors and A23 comparison evidence | Legacy workflow displays three geometry diagnostics. Draft requests eight descriptors with their own finite populations and unduplicated evaluation-level ESS. Partial coverage only. |
| S7 | Legacy `supp_fig7` builder/plot cells retained as provenance; Figure 4 candidate/support tables | Draft requests the FACOAL204 selection audit and exact degenerate RIPTiDe counterpart. Legacy workflow adds two post hoc examples and says no specific frozen main-example selection criterion was found. Author reconciliation required. |

Legacy code cells in `manifests/provenance/legacy_figure_code/` contain no notebook
outputs and are not executable publication entry points. Their source identities
and extracted cells are recorded in `NOTEBOOK_EXTRACTIONS.json`.

## Scientific populations and upstream boundary

The publication names the glucose anchor HEX1. However,
`scripts/dmi_bridge_a20_dual_anchor_qualification_v1.py::load_panel` constructs
the candidate Vmax observable as `max(-EX_glc__D_e, 0)` and ranks that observable.
The A20 contract describes Vmax→HEX1. This implementation/interpretation difference
requires author verification; no operator substitution or equivalence assumption
was made during migration.

The common geometry population is 400 evaluations × 4,179 reactions. The utility
population is the A23 matched 836 truth-pair intersection; q10/q50/q90 aliases and
held-out-vector handling remain governed by PL2A and A22. Truth ties and numerical
gain ties are different concepts. Figure 3 averages distinct non-tie truths within
each reaction–evaluation pair, then weights those pairs equally.

The original A1 holdout was frozen before confirmation. A2 reuses the split for
robustness. Reaction bootstraps quantify computational stability, not biological
replication. The original 2,455-row development snapshot and the final 2,454-point
Figure 2D scatter have different frozen eligibility definitions and are not merged.

Upstream dependency hashes are in `FROZEN_LINEAGE.json`. The candidate-vector cache
has a historical path inside a six-arm exploratory directory; final Bridge scripts
explicitly pin that individual cache. This does not justify migrating the exploratory
directory. Bulk inputs, missing foundation patches and external model/data rights
define the current expensive-upstream boundary.

## Supplementary tables

| Draft table | Frozen evidence and reproduction boundary |
|---|---|
| S1: population/denominator accounting | A23 matched audit, Figure 3 metadata, PL2 split/confirmation manifests and Figure 4 manifest. Counts are documented; original numerical artifacts are withheld. |
| S2: matched gain summaries under two aggregation conventions | Final pair-level Figure 3 summaries plus `fig3_gain_definition_audit.parquet.xz` as a provenance-only comparison of case pooling. That audit is relevant to this explicit table; it must not replace main Figure 3's pair-weighted endpoints. |
| S3: original A1 held-out context associations | PL2D context results and M1 confirmation-context table. No matching final table-release workflow identified; upstream files remain withheld. |
| S4: primary/supportive A1 held-out associations | PL2D primary/supportive result records and frozen exclusion/coverage definitions. Author must finalize the table's frozen artifact locator and release rights. |
| S5: full paired descriptor summaries | A23 geometry comparison table and PL1/A21 descriptors. Preserve metric-specific populations; do not duplicate evaluation-level ESS over reactions. Original tables remain in the source. |

The lightweight table command exports only traced Figure 2/3 summaries; it does
not claim to reconstruct all five supplementary tables.
