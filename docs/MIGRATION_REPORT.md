# Migration report

## Reproducibility verdict

**NOT YET REPRODUCIBLE.** The authorized code, provenance, documentation and
lightweight entry-point curation is committed. Complete manuscript reproduction
is blocked by pending numerical-data rights review, incomplete final artwork and
unresolved scientific/figure lineages. No frozen numerical input or ambiguous
artwork has been copied, and no figure has been regenerated. This is not a complete
publication release.

## Source state

Source: `/home/pty/work/project_MRI_brain_tumor/`.
Commit: `16379a169f45ec556274953eb66ae12293e1e985`, clean `main` checkout when
inspected. The source is strictly read-only regardless of filesystem permissions.
No source script, tests, build, solver or sampling workflow was run. No source
working tree, index, branch/ref, cache, generated output or intentional timestamp
write was performed. Source Git inspection disables optional locks. No external
dataset was accessed.

## Destination state

Initial commit: `63edf6a9f3e0f3ba84c17dd7a828db5b6d7170de`, clean `main` checkout
containing `DATA_PROVENANCE.md`. That original file is preserved unchanged.
Destination branch: `publication-migration`. All persistent writes, logs, temporary
work and caches are under `/home/pty/work/recomb2027-dmi/`.

Logical commits group inventory/rights review, production implementations,
main-figure entry points, supplementary provenance, reproducibility documentation,
licenses/environment/citation, and final audit. A frozen-data commit is intentionally
absent because rights remain pending. No history rewrite or remote push occurred.
The final destination SHA is obtained with `git rev-parse HEAD` and reported with
this report; a commit cannot contain its own hash.

## Files copied

The file-level record is `manifests/provenance/SELECTION.tsv`. The current transfer
comparison contains **101 mappings**: 88 byte-identical files, four path-adapted
builders and nine source plotting-code extractions.

- 34 native Bridge implementation/dependency modules. BIO0/qualification helpers
  retain provenance-only roles; dropped concordance is not an entry point.
- 21 native Stage-11 reconstruction, qualification and retained-vector modules.
- 11 existing synthetic scientific-core test files.
- Three main-figure data builders and one S1 builder, with only packaging/path edits.
- Five selected main/S1 plotting scripts extracted without executed notebook outputs.
- Four legacy supplementary builders and four extracted plotting-cell files,
  retained only as provenance for manuscript discrepancies.
- 18 original contracts, audit documents, environment/configuration records.

New documentation, licenses, citation metadata, registries, checksums and the
lightweight `reproduce.py` entry point are separate from source transfers.
Numerical data copied: **zero**. Artwork copied: **zero**. Fifteen frozen stage
identities and source artifact hashes are recorded without releasing their data.

## Files intentionally excluded

Alternate secondary-anchor branches, AKGDm and six-arm explorations, failed pilots,
superseded work directories, abandoned figure candidates, old Figure 3 composites,
obsolete Figure 4 builder, dropped RNA–flux concordance and alternative scoring
pipelines remain in the source. Standalone historical reconstruction methods,
raw inputs, third-party GEMs/papers, bulk ensembles/models, solver logs, caches and
scratch notebooks are excluded.

Original library modules retain dormant optional/historical branches when removing
them would require rewriting validated code or change implementation fingerprints.
They are not selected by publication entry points or promoted into production.

The individual cache under a historical exploratory path is a hash-pinned input
to final Bridge code. That lineage is recorded as a narrow upstream exception;
neither the cache nor its surrounding exploratory directory is copied.

The case-pooled Figure 3 audit is relevant only as the draft Table S2 aggregation
comparison. It remains withheld for rights review and must not replace final
pair-weighted main-figure endpoints. Decompressed or redundant dataset copies
are not committed.

## Code modifications

### Path, packaging and display changes

Original Bridge, Stage-11 and test files retain their source bytes. Native locations
are preserved to retain imports, root discovery and hash-pinned implementation
identities. This is the documented modest variation from the illustrative `src/`
layout; no duplicate subsystem implementation was created.

The four relocated builders adjust root discovery and output locations. Main
builder outputs go under `reproduced/derived/`; S1 also locates the prospective
curated Figure 4 manifest. Exact diffs are in `PACKAGING_DIFFS.patch`.

Notebook code is extracted by recorded cell index and hash. Wrappers provide
repository-relative frozen inputs, ignored output directories and noninteractive
plotting. Figure 2 selects the existing usefulness plotting cells and labels the
axis as a benchmark frequency. Figure 4C's unsupported physical-unit label becomes
native GEM flux units, following the technical audit. Original source calculations,
histograms, weights, thresholds and plot grids are preserved. Trailing whitespace
in extracted code is normalized only after AST equality is verified.

`reproduce.py` adds checksum/input validation, repository-boundary checks,
per-job identity/log/output manifests, valid-output skipping and refusal to overwrite
incompatible outputs. It performs no reconstruction, optimization or flux sampling.
All selected rendering temporary files/caches are directed inside the repository.

Three narrow `.gitattributes` whitespace exceptions preserve an original package
EOF, original audit Markdown formatting and valid unified-patch context. They do
not remove any file from checksum verification. Newly authored code/metadata pass
whitespace checks.

### Scientific changes

**Zero.** No equation, estimator, threshold, tolerance, random seed, reaction axis,
sampling rule, truth-pair definition, population, aggregation, objective or model
constraint was changed. Existing hash/qualification gates are retained. Scientific
ambiguities below are recorded rather than repaired by changing code.

## Unresolved issues for author review

| ID | Issue | Required follow-up |
|---|---|---|
| R1 | Redistribution rights pending | Review 50 exact numerical artifacts, four conceptual SVGs and three external-source groups in `RIGHTS_REVIEW.tsv`. Confirm original ownership/derivative rights before copying. |
| S1 | Glucose operator/anchor interpretation | A20 constructs candidate Vmax from `max(-EX_glc__D_e, 0)` while manuscript/contracts describe Vmax→HEX1. Establish the intended interpretation/equivalence; migration does not substitute HEX1 or change weights. |
| S2 | Evaluation/draw wording | The frozen registry has 4 methods × 4 RNA contexts × 5 CT2A subjects × 5 GL261 subjects = 400 evaluations. Reconcile the proposed two-mice/five-draws description; do not invent an independent draw axis. |
| F1 | Conceptual/final artwork | Figure 1 A-C, 2 A-C, 3A and 4A have identified conceptual sources but ownership/final editorial identity remains pending. |
| F2 | Figure 3D | Existing D SVG/raster labels are superseded case-pooled means. Final pair-level summaries are traceable, but corrected authoritative artwork/source is missing. |
| F3 | Figure 4C | Source draws four methods; main draft uses three and moves RIPTiDe to S7. Confirm final composition and the two-stage example-selection intent. |
| F4 | S1 heatmap order | Source matrix rows use alphabetical method order while labels use iMAT/GIMME/CORDA/RIPTiDe. Review/fix the label/data alignment without changing the underlying population. Original behavior is preserved; rendering is blocked. |
| F5 | S2/S3 numbering and scientific content | Legacy S2/S3 show paired geometry; draft requests context associations/alternative descriptors. Final plotting workflows are missing. |
| F6 | S4/S5 incomplete panels | Existing S4 lacks the draft's full exclusion/coverage display; existing S5 lacks all requested negative-tail/reliability displays. Do not promote exploratory replacements. |
| F7 | S6 descriptor coverage | Existing S6 shows three diagnostics; draft requests eight with metric-specific finite populations and unduplicated evaluation ESS. |
| F8 | S7 selection lineage | Legacy manifest adds two post hoc examples and says no specific frozen main-example selection criterion was found. Reconcile with the draft's selection-audit claim and retain the exact degenerate counterpart. |
| U1 | Missing foundation artifacts | Twelve hash-gated patch files are absent from the source checkout. Restore exact artifacts for upstream reproduction; do not delete/reset their gates. |
| U2 | External/model/archive boundary | Model/raw/medium/projection inputs and bulk artifacts are not redistributed. Establish licenses and any release archive/download process before promising a full upstream rerun. |
| U3 | Environment and citation | Optional solver/modeling dependencies were not installed/imported/solved here. Final title, DOI, ORCIDs, conference citation and archive identifiers remain TODOs. |

The main/supplementary manuscript snapshots and technical audit remain authority.
No legacy folder name or prior PASS report resolves these disagreements.

## Validation

- Source transfers: 88 byte-identical files verified; four builder adaptations
  checked against an explicit path-only transformation; all nine extracted source
  cell bodies checked against their recorded transformations and AST-preserving
  whitespace normalization. Audit record: `SOURCE_COPY_AUDIT.tsv`.
- Python: 84 relevant files compiled; 34 lightweight Bridge modules imported
  without executing their main functions. Optional Stage-11 modeling imports were
  not attempted. Existing synthetic scientific-core suite: **108 passed**.
- Prepared CLI: status, pending-input validation/figures/tables, unknown-job and
  output-boundary failures checked. Missing/withheld inputs cannot yield a success
  claim or invoke an upstream solve.
- Production path search: no dependency on the old absolute repository path in
  production-facing Python. Historical provenance documents may retain it.
- Terminology: public definitions use A2; legacy labels remain only in original
  schemas/code/provenance or explanatory historical notes.
- Figure dependency audit: all four main figures and draft S1-S7 mapped; incomplete
  sources/data/compositions are explicit. Supplementary Tables S1-S5 also mapped.
- Frozen data verification/rendering: **not performed on destination data**, because
  no numerical input is cleared/copied. Fifty source input hashes are recorded;
  that is not copied-artifact or full-render validation.
- Included release-file checksums, CFF YAML/required-field checks, full migration
  diff whitespace checks and final Git/source-state checks are recorded in
  `manifests/provenance/VALIDATION.json`.

## Final repository audit

| Question | Result |
|---|---|
| Code for every main figure identifiable? | Partial: quantitative B/C/D-G source mappings exist; corrected Figure 3D/final artwork gaps are explicit. |
| Frozen quantitative main data identifiable? | Yes in the source with hashes; no released copies while rights remain pending. |
| A1/A2 unambiguous? | Publication labels are defined; glucose operator interpretation remains a scientific author issue. |
| Development versus confirmation clear? | Yes; A2 is post-freeze robustness, not independent confirmation. |
| All manuscript plots regenerable without upstream reconstruction? | Not yet; prepared quantitative commands await data rights and final panel reconciliation. |
| Outputs traceable to source? | Yes, via source commit, native paths, manifests and SHA-256 records. |
| Exploratory/superseded analyses absent from production workflow? | Yes; source/provenance exceptions and dormant native branches are explicitly separated from publication entry points. |
| External materials/licensing boundaries explicit? | Yes; ambiguous products withheld. |
| Paths portable? | Prepared entry points and relocated builders use this root; full upstream is explicitly incomplete. |
| Relevant files left only in source explained? | Identified main/supplementary/table dependencies and upstream omissions are recorded; release completeness is not claimed. |

## Next author actions

Confirm redistribution decisions, reconcile the glucose operator/evaluation
terminology, finalize artwork and supplementary panel/table choices, and restore
missing upstream identity artifacts where a full upstream reproduction is desired.
Then include only cleared frozen inputs, verify their source-copy checksums,
execute selected plotting workflows and repeat the full figure/release audit.
