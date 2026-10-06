# Data availability and redistribution policy

Included compact project-derived figure inputs and their SHA-256 identities,
sizes and release decisions are listed in
[`RIGHTS_REVIEW.tsv`](../manifests/provenance/RIGHTS_REVIEW.tsv).
[`FIGURE_INPUTS.json` records each source and destination identity. The selected
release set includes main-figure inputs plus the audited S1–S3 diagnostic tables;
all payloads have source SHA-256 and destination SHA-256 records. Two JSON
manifests have repository-relative paths adapted for this checkout, with the
source identity retained. `MANUSCRIPT_TABLES.json` records the generated result
summary identities and their frozen source tables.

## Included material

The repository includes code, configurations as provenance, source commit/file
identities, current artwork, and compact frozen figure inputs. It includes no raw
biological datasets, full parent model, reconstructed model files, bulk flux
ensembles, or solver caches. Superseded or assembly-ambiguous artwork is labeled
provenance-only in `manifests/ARTWORK.json`.

## Remaining upstream requirements

- Raw RNA and DMI source datasets, parent GEM, frozen medium/projection inputs,
  per-job seeds, historical identity artifacts, and missing foundation patches
  are not bundled. A complete upstream reconstruction/solver rerun is outside
  the included-data reproduction workflow.
- Figure 4 paired geometry is included to support Figure 4B and S2/S3 plotting;
  redundant paired utility partitions are omitted.
- S4–S7 current renders are retained as source artwork, while matching standalone
  data-to-plot pipelines and all manuscript-requested panels remain incomplete.

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
