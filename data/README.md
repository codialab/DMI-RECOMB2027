# Data availability and redistribution policy

No pending research dataset has been copied. The exact candidate figure-input
paths, SHA-256 identities, sizes and rights decisions are listed in
[`RIGHTS_REVIEW.tsv`](../manifests/provenance/RIGHTS_REVIEW.tsv).
[`FIGURE_INPUTS.json`](../manifests/FIGURE_INPUTS.json) records prospective locations
for 23 main-figure/S1 inputs and one provenance-only Table S2 aggregation audit. A
prospective path does not mean the file is included or approved.

## Included material

The repository includes code, configurations as provenance, source commit/file
identities, and documentation. It includes no raw biological datasets, full parent
model, reconstructed model files, bulk flux ensembles, solver caches, or ambiguous
legacy artwork. The original numerical results remain in the read-only source.

## Derived products awaiting author review

- Figure 2 reaction-level point table: one reaction per panel after frozen
  within-evaluation and reaction aggregation.
- Figure 3 equal-weight reaction–evaluation gains, cue-direction records, and
  final endpoint summaries, with compressed authoritative representations retained.
- Figure 4 paired geometry partitions, selected-candidate/support tables and
  the frozen data manifest. The redundant paired utility partitions are unnecessary
  for the selected main plotting entry points and are not slated for Git.
- S1 compact candidate-weight and support-diagnostic tables. Other supplement
  tables also require scientific-lineage reconciliation before inclusion.

These are computational derivatives of candidate GEM flux vectors and rank-based
DMI weighting. Traceable computational origin does not by itself establish all
redistribution rights. Author review must confirm ownership and derivative rights.

## External sources and upstream boundary

The source records identify Khalsa et al. GSE151414 processed TPM as the principal
RNA input, Simões et al. 2025 DMI data at Dryad DOI
`10.5061/dryad.905qfttwb`, and the Khodaee et al. 2020 iMM1865 GEM. Consult those
original deposits/publications and their licenses to obtain raw inputs. No external
data source was accessed during migration.

Local Dryad metadata declares the DMI deposit CC0-1.0. Redistribution permission for
the parent GEM and Khalsa source files was not established from inspected local
records. The distinct RNA and DMI cohorts have no subject-level pairing crosswalk.

Rebuilding upstream requires the Stage-11 frozen medium, expression/projection
inputs, per-job seeds and qualification/provenance gates as well as the data above.
Some hash-pinned foundation patch files are absent from the source checkout.
Restoring a missing dependency is an author follow-up; deleting a gate or substituting
a solver is not an acceptable reproduction shortcut.

## Size and storage decisions

The 50 reviewed figure-related files total about 306.8 MB compressed. The selected
main/S1 subset is smaller and preserves `.xz` containers. Do not commit decompressed
Parquet copies. Per-method Figure 4 geometry partitions are about 2.3–29.8 MB each;
the larger redundant utility partitions (about 5.9–53.3 MB each) are omitted.

Full candidate ensembles, reconstructed models, large case partitions and caches
belong in an author-approved release archive with checksums if required. No archive,
Git LFS upload, external publication or remote push was performed. Their availability,
licenses and download instructions must be established before promising a complete
upstream rerun. Figures should depend on cleared frozen inputs, not bulk archives.
