# Manuscript reproducibility matrix

Source authority: clean `project_MRI_brain_tumor` checkout at `225279129`.
“Reproducible” means regenerated from the frozen files included here by the
root-level CLI, with population and numerical checks passing. A retained current
render is evidence/artwork, not by itself a successful rebuild.

## Main figures

| Manuscript element | Included source/input | Reproduction status | Population / endpoint and qualification |
|---|---|---|---|
| Fig. 1A–C | Current editable SVG and PNG in `figures/released/fig1/` | Artwork retained; conceptual and composition review remains | No numerical endpoint. Editable component retained; final assembly identity is not asserted. |
| Fig. 2A–C | Current composite and rough conceptual SVG in `figures/released/fig2/` | Artwork retained; conceptual panels not data-regenerated | No numerical endpoint. |
| Fig. 2D–G | Frozen `figure2_DG_points.csv.xz`; `plot_fig2_dg.py` | **Reproducible** (`fig2DG`) | Plotted reaction counts 2,454 / 2,485 / 838 / 843; Spearman ρ = 0.8571160 / 0.8339440 / 0.8351436 / 0.8493657. Point table and aggregation are source-frozen. |
| Fig. 3A | Current editable SVG in `figures/released/fig3/` | Artwork retained; conceptual content | No numerical endpoint. |
| Fig. 3B–C | Frozen reaction–evaluation table; `plot_fig3_bc.py` | **Reproducible** (`fig3BC`) | 660,237 eligible reaction–evaluation observations per anchor, 3,625 reactions, 370 evaluations; pair-weighted correct/wrong gain summaries are in `tables/manuscript/fig3_pair_gain_summary.tsv`. |
| Fig. 3D | Historical panel-D SVG retained only for provenance | **Superseded** | The source audit identifies the old artwork as case-pooled; it must not be presented as the final pair-weighted panel. Corrected standalone assembly is not available. |
| Fig. 4A | Current draft SVG in `figures/released/fig4/` | Artwork retained; conceptual content | No numerical endpoint. |
| Fig. 4B | Four frozen method partitions; `plot_fig4_b.py` | **Reproducible** (`fig4B`) | 400 evaluations × 4,179 reaction coordinates; four methods, 417,900 rows per method, unique evaluation/reaction pairs. |
| Fig. 4C | Candidate group/metric/support tables are included | Data retained; assembly/selection wording unresolved | Source notebook uses four methods while the manuscript draft specifies three; no method or selection rule was changed here. |

## Supplementary figures and tables

| Element | Included frozen evidence/artwork | Reproduction status | Scope and outstanding qualification |
|---|---|---|---|
| Fig. S1 | Candidate weights, conditional/global ESS, contrast-product and post-holdout support tables; current SVG/PNG | **Reproducible** (`s1`; standalone script) | 6,400 candidate rows, 320 strata, 20 candidates per stratum; 836 matched truth pairs for post-holdout support. ESS values are descriptive, not an exclusion rule. |
| Fig. S2 | Current source render and source audit retained | **Reproducible** (`s2`; notebook-cell plot extraction) | Method-stratified A2−A1 shifts; source audit reports 1,671,600 paired rows, 417,900 per method. |
| Fig. S3 | Current source render and source audit retained | **Reproducible** (`s3`; notebook-cell plot extraction) | Pooled A1/A2 geometry relationship; same 1,671,600 paired rows. Tiny serialization differences are recorded in the current source audit. |
| Fig. S4 | Current source render and its validation report | Artwork retained; full manuscript scope not established | Existing source covers bootstrap/development/context panels; source audit notes additional draft coverage and reaction-exclusion requests. |
| Fig. S5 | Current source render and its validation report | Artwork retained; full manuscript scope not established | Existing source covers context gains/outcome composition; draft also requests negative tails and reliability curves. |
| Fig. S6 | Current source render and its validation report | Artwork retained; full descriptor set not established | Existing source plots three geometry diagnostics; the draft requests eight descriptors with metric-specific denominators and unduplicated evaluation ESS. |
| Fig. S7 | Current source render and selection manifest | Artwork retained; example lineage unresolved | Two post-hoc examples are recorded; the exact frozen main-example selection rule remains unresolved. |
| Table S2 / aggregation audit | `fig3_gain_definition_audit.parquet.xz` under `data/processed/aggregation_audit/` | Provenance-only | Case-pooled definitions are superseded for primary Figure 3 results; primary summaries use equal-weight matched-pair aggregation. |

The root CLI regenerates the quantitative figure panels for which current standalone
plotting code and an authoritative input contract are available (Fig. 2D–G, 3B–C,
4B, and S1–S3). Current S1–S7 artwork, notebook/source provenance, and source audits
are retained under `figures/released/` and `manifests/provenance/source/` as
appropriate; unresolved items are not relabeled as reproducible. `status`,
`validate`, `figures`, and `tables` never run upstream reconstruction, Gurobi,
flux sampling, or candidate generation.
