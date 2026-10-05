# RECOMB Supplementary Figure Wave 1 — Controller/QC Audit

**Audit date:** 2026-10-04  
**Authority:** `docs/FIGURE_CONTRACT_S1_S2_S3.md` and its S1–S3 task files  
**Scope:** Frozen worker commits only; no reconstruction, optimization, sampling, candidate generation, or pipeline rerun.

## Summary and verdict

The three supplied commits share the same parent, `01c2389f99be95f66089f7196a8eb38b6db437d8`, and each commit changes files only under its assigned `figures/supp_figN/` directory. All have a table builder, executed notebook, README, validation report, SVG, PNG, and compact plotting tables. No source artifact, manifest, manuscript, or shared production code was changed.

The frozen-table checks support the reported S1–S3 counts and principal results. In particular, S2 and S3 have identical full `(method, evaluation_id, reaction_id)` key sets (1,671,600 observations); their paired metric values agree within at most `2.3e-16`, consistent with serialization-level floating-point differences. Both use the Figure 4 paired geometry source and the same A1/A2 definitions. The figures serve distinct roles and neither changes the production analysis.

**Verdict: PASS WITH IMPORTANT/COSMETIC FOLLOW-UP — safe to merge.** The worker commits can be merged unchanged. The IMPORTANT items below concern manuscript integration and wording in the S1 validation report; neither changes the computed figures or invalidates the frozen tables.

## Findings

| Class | Figure | File / location | Issue and evidence | Recommended action |
|---|---|---|---|---|
| IMPORTANT | S1 | `figures/supp_fig1/validation_report.md`, “Integrity and scientific checks” | The report says “Highly concentrated cases remain present” and then reports 153 strata with `ESS ≤ 2`. The figure and README do not use this as a cutoff, and the report records it as a validation count, but it does not explicitly say the value 2 is descriptive rather than a new scientific threshold. | Clarify in the report that `ESS ≤ 2` is a descriptive validation count only, not a prespecified scientific classification. Keep it out of publication-facing interpretation unless the contract is explicitly revised. |
| IMPORTANT | S3 / manuscript integration | `figures/supp_fig3/validation_report.md`, “Current Figure 4 assembly and open issues” | S3 reports that the manuscript draft calls the paired scatter Figure S6, while the current plan assigns it to S3. The authoritative contract and task assign this content to S3. The report’s claim could not be independently confirmed in the available manuscript text: the two repository-root DOCX files do not contain a Figure S6 or paired A1/A2 reference, and no manuscript source text with that cross-reference was found. | During manuscript integration, locate and update the stale paired-scatter citation from Figure S6 to Figure S3, then check all supplement cross-references and caption numbering. |
| COSMETIC | S1 | `figures/supp_fig1/outputs/supp_fig1.png` | At the PNG’s native 1472 × 1668 size the panels are clear, all strata are represented in the heatmaps, and global, conditional, contrast-product, and post-holdout support have separate panels. The dense method/context labels and candidate heatmaps will become small at manuscript reduction. The x-axis identifies RNA contexts only by index 1–4; the raw context keys are listed in the README/validation artifacts. | At final layout, check legibility at actual print width and provide a compact reader-facing mapping for context indices if the caption or surrounding text does not already define them. |

No BLOCKING findings were identified.

## Figure-specific verification

### S1 — candidate weights and effective support

- **Structure:** 6,400 candidate rows form 320 unique conditional strata: 160 for A1 and 160 for A2. Each stratum has 20 candidates. The dimensions are four methods (iMAT, GIMME, CORDA, RIPTiDe), four RNA contexts, two tumors (CT2A, GL261), five mice per tumor, and two anchor settings.
- **Anchor definitions:** A1 is 1-Strong-Anchor (HEX1); A2 is 2-Strong-Anchors (HEX1 + LDH_L). S1 uses the frozen A1 `strong_anchor_baseline` weights and qualified A2 weights. Source hashes are recorded in the S1 validation summary and checked against the applicable frozen manifests / Figure 4 source manifest.
- **Weights and conditional ESS:** Candidate weights are finite and nonnegative, candidate keys are unique, and conditional weights normalize within strata (maximum observed normalization error about `1.3e-15`). Recomputing `1 / sum(w²)` from the 20 conditional weights gives a range of 1.00–16.24 and matches the saved stratum summaries within `2.9e-14`. The descriptive count `ESS ≤ 2` recomputes to 153 strata.
- **Separate support quantities:** The global pool contains 320 states per tumor/mouse/anchor and has ESS approximately 20 in all 20 pools. Product-weight contrast support is separately represented for 400 evaluations × two anchors (800 rows), with ESS range approximately 1.00–201.03. Post-holdout support is based on exactly 836 matched truth pairs; the 3,344 pair-support rows cover two tumors and two anchors. Neither 885 nor 1,110 replaces this fixed population.
- **Concentrated support:** All strata remain represented; there is no exclusion based on ESS or maximum candidate weight. The heatmaps show the candidate weights for each stratum. The figure footer states that effective support is a finite-weight diagnostic, not a count of independent observations.
- **Visual check:** Four tumor/anchor conditional-ESS panels, four all-strata weight heatmaps, and separate global, product-weight, and post-holdout diagnostics are visible. The global target is distinguished from the conditional ESS and product-weight ESS. S1 is information-dense; final print-scale legibility deserves the cosmetic check above.

### S2 — method-stratified paired shifts

- **Source and provenance:** The builder reads the four Figure 4 paired geometry partitions and verifies their SHA-256 values against `figures/fig4/data/fig4_data_manifest.json` (manifest status PASS). The manifest identifies frozen PL1 A1 geometry and A21 A2 geometry. No exploratory alternative-anchor input or truth-pair restriction is used.
- **Population and matching:** Each metric table contains 1,671,600 rows, 417,900 per method, with unique `(algorithm, evaluation_id, reaction_id)` keys and no missing pairing keys. Across the tables, there are 400 evaluations × 4,179 reactions. The authoritative A1/A2 pairing is `(evaluation_id, reaction_id)`; Figure 4’s source builder performs one-to-one matching and checks the paired metadata.
- **Recalculation:** Both tables have finite values, no nonfinite exclusions, and stored differences equal `A2 − A1` within `5e-12` (the deterministic recheck found no mismatches). For directional entropy, mean shifts are positive for all methods (iMAT 0.02224, GIMME 0.03469, CORDA 0.04872, RIPTiDe 0.00068); for dominant-direction mass, means are negative for all methods (−0.01069, −0.01765, −0.02485, and approximately −0.00002, respectively). Every method’s median shift is zero.
- **Presentation:** The two method-stratified ECDFs show A2 − A1, a visible zero reference, and `n=417,900` for each method in each panel. The curves are distinguishable at native 2048 × 1024 PNG size. Both metrics are dimensionless, as labeled.

### S3 — pooled paired A1-versus-A2 values

- **Source and provenance:** S3 uses the same four Figure 4 paired geometry partitions and manifest hashes as S2. The source manifest is FIG4-DATA-B v1.1.0, status PASS. The compact S3 table records the same four source hashes.
- **Population and pairing:** The compact table has 1,671,600 unique `(evaluation_id, reaction_id)` rows, 400 evaluations, 4,179 reactions, and 417,900 rows per method. A1/A2 metadata and metric columns are retained together. Both panels use all rows; no finite-value exclusions, truth-pair restrictions, or alternative-anchor data enter.
- **Metric definitions and ranges:** The builder and README define directional entropy as Shannon entropy over negative, tie, and positive sign masses, normalized by `log(3)`, and dominant-direction mass as the largest of those three probabilities. Both are probability-like dimensionless descriptors expected in `[0,1]`. Observed entropy values range from about `−2.1e-16` to 0.8452; dominant-direction mass ranges from 0.5 to `1.0000000000000002`. There are no excursions larger than `1e-12`; values are not clipped and the plot limits retain machine-precision roundoff.
- **Identity summaries:** With the implemented `1e-12` equality tolerance, entropy is above/on/below the identity for about 17.19% / 78.87% / 3.94% of pairs. Dominant-direction mass is above/on/below for about 3.91% / 79.15% / 16.94%. The conclusions retain the many unchanged observations and do not imply every pair shifts in one direction.
- **Presentation:** Both panels pool methods as specified, use A1 on x and A2 on y, display log-scaled hexbin density, show a dashed identity line and equal 0–1 axes, and report the same exact `n=1,671,600`. The density does not hide the overall identity relationship at native 2048 × 1024 PNG size.

## Cross-figure consistency, terminology, and roles

- **S2 ↔ S3:** Their builder code reads the same hash-verified Figure 4 geometry partitions and manifest. A deterministic comparison of the compact tables found identical complete `(method, evaluation_id, reaction_id)` key sets with 1,671,600 entries. A1/A2 metric columns agree to within `2.3e-16`; S2 differences agree with the S3 paired values to floating-point precision. The tiny numerical differences do not alter sample membership, ranges, mean directions, or reported summaries.
- **Terms:** Figure labels and README definitions use A1 = 1-Strong-Anchor (HEX1), A2 = 2-Strong-Anchors (HEX1 + LDH_L), and the exact method names. No prohibited “target reaction” wording appears in S1–S3 figure labels or README/report prose. The raw S1 RNA-context provenance keys include `setx1`/`setx3`; these are not used as plot labels. `A2-L` appears in the legacy Figure 4 QC notebook and provenance, outside the audited supplementary labels and plots.
- **Flux units:** S1’s weights and ESS and S2/S3’s entropy and probability-like geometry measures are dimensionless. No flux-valued quantity is plotted in these figures, so flux units are not required on their axes.
- **Distinct purpose:** S1 describes finite candidate weighting and effective support; S2 shows method-stratified A2 − A1 geometry shifts; S3 shows the pooled absolute paired A1/A2 relationship. The current Figure 4 assembly described in the repository uses absolute metric distributions and a selected weighted `P(Δv_B)` example, so S2/S3 add paired information rather than duplicating those assembled panels.

## Reproducibility and artifact hygiene

| Figure | Builder | Executed notebook | README / report | Compact tables | SVG / PNG | Frozen input checks |
|---|---|---|---|---|---|---|
| S1 | Present | 1 code cell executed; output present | Present | Candidate, stratum, global, contrast, and post-holdout tables | Present | Source manifests and hashes recorded; A1 hash tied to Figure 4 manifest |
| S2 | Present | 2 code cells executed; outputs present | Present | Two compressed shift tables and JSON summary | Present | All four Figure 4 partition hashes verified |
| S3 | Present | 2 code cells executed; outputs present | Present | Paired geometry table and JSON build manifest | Present | All four Figure 4 partition hashes verified |

The committed S2 tables are about 11.9 MB and 11.4 MB compressed; S3’s paired table is about 7.0 MB. These sizes are reasonable for 1.67 million paired rows and are not individually close to common Git/GitHub large-file limits. The notebooks contain executed output; none is disproportionately large relative to the supplied plotting tables and renderings. No duplicated copy of the frozen Figure 4 source partitions was added.

The worker commits contain only their owned figure directories. This audit verified the specified refs and their parent commits from the repository. It did not inspect the separate worker worktree directories, which are outside the repository data boundary; cleanliness here refers to the committed snapshots and the clean main worktree at audit start.

## Summary table

| Figure | Commit | Scientific checks | Reproducibility | Visual QC | Merge status |
|---|---|---|---|---|---|
| S1 | `8164b43add876c7f690dce887fe6ae4a816173c0` | PASS; counts, weights, ESS, fixed 836-pair support verified | PASS; builder, executed notebook, provenance, compact data, report, SVG/PNG present | PASS with cosmetic final-size and context-index clarity follow-up | Safe to merge unchanged |
| S2 | `cb6d966cb1e1aeffa4cc493b22126e27992d6c77` | PASS; 1,671,600 matched pairs, actual shift directions and summaries verified | PASS; common Figure 4 manifest and source hashes, builder, executed notebook, compact tables | PASS; ECDFs and zero reference legible | Safe to merge unchanged |
| S3 | `41e7fc43e3e29b79f243ca17ecd43afa2e89b105` | PASS; same paired population, metric ranges, roundoff, and summaries verified | PASS; common Figure 4 manifest and source hashes, builder, executed notebook, compact table | PASS; density, identity line, equal limits, and n legible | Safe to merge unchanged |

## Checks performed

- Compared each commit’s complete diff against the common base and inspected committed file inventories and sizes.
- Read the authoritative contract/tasks, worker builders, READMEs, validation reports, machine-readable summaries, and Figure 4 manifest.
- Recomputed S1 conditional and global ESS from committed weights and verified candidate counts, normalization, contrast support, and fixed holdout rows.
- Read both frozen S2 shift tables and the S3 paired table; checked keys, duplicate counts, finite values, A2 − A1 arithmetic, S2/S3 population equality, metric ranges, and reported summaries.
- Checked notebook execution counts and visually reviewed the committed PNGs. No production workflow, table builder, notebook, reconstruction, optimization, sampling, or candidate-generation pipeline was run.
- Searched repository manuscript text and inspected text from the two repository-root DOCX files for the S6 reference; no independent primary-manuscript occurrence was found.
