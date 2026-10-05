# Data and implementation provenance

Source repository: `/home/pty/work/project_MRI_brain_tumor/`.
Source commit: **`16379a169f45ec556274953eb66ae12293e1e985`**.
Migration date: **2026-10-05**.
Destination initial commit: `63edf6a9f3e0f3ba84c17dd7a828db5b6d7170de`.
Destination branch: `publication-migration`.

The source repository is strictly read-only even where OS permissions might allow
writes. No source working-tree, index, branch/ref, cache, generated-output or
intentional timestamp change is made. All persistent migration files, runtime
caches and logs are under the destination. Source Git inspection uses optional
locks disabled; nothing is moved or deleted from the source. No external dataset
was accessed, and no source script or expensive scientific computation was run.

## Authority and manuscript snapshots

The main draft SHA-256 is
`6b837c1229b9f495aa1dcfab4b1512c7ef4c1df85e3aef3ab07b174ac68e467b`;
the supplement draft SHA-256 is
`ae62ba34b22d6bf1d81a7c39535430abb091940ba8063e4b0aa51471143c5fcc`.
Their source locators are `manuscript_audit_inputs/RECOMB_with_A2_Manuscript_Draft.docx`
and `manuscript_audit_inputs/RECOMB_with_A2_Supplementary_Materials_Draft.docx`.
They were inspected in memory; manuscript files and embedded artwork were not
copied. The source technical audit and selected frozen contracts are preserved
verbatim under `manifests/provenance/source/`.

Source evidence may be internally inconsistent. It is not used to silently resolve
the manuscript's figure numbering, panel selection or operator interpretation.
The panel and table dependency audit is in [ANALYSIS_MAP.md](ANALYSIS_MAP.md).

## Artifact lineage

[`FROZEN_LINEAGE.json`](../manifests/provenance/FROZEN_LINEAGE.json) captures the
identities of 15 frozen stages: A1 weight qualification, A1 geometry/truth/cue
utility, development/freeze/amendment/confirmation, A2 qualification/geometry/cue
utility, matched synthesis, and manuscript evidence consolidation.

The final A1 geometry lineage uses the checksum-bound resumable PL1 output, not an
arbitrary newer runtime variant. A20/A21/A22 establish the qualified A2 lineage;
A23 freezes the 836-pair common population. Figure 3 final means weight eligible
reaction–evaluation pairs equally. Historical case pooling is retained only as a
prospective provenance comparison required by draft Table S2.

The common candidate axis has 4,181 coordinates; excluding HEX1 and LDH_L yields
4,179. This is not the parent GEM reaction count or the distinct Stage-4 projection
dimension. Original code and schemas retain these separate definitions.

[`EVALUATION_LAYOUT.json`](../manifests/provenance/EVALUATION_LAYOUT.json) records
the exact registry identity and the 4 × 4 × 5 × 5 evaluation layout. The user's
proposed “two mice / five draws” wording requires reconciliation with that frozen
registry; the migration does not invent an independent reconstruction-draw axis.

## Copies and path/import changes

[`SELECTION.tsv`](../manifests/provenance/SELECTION.tsv) is the file/directory
classification and mapping record. Original Bridge scripts and Stage-11 modules
retain their native locations and code bytes, preserving imports and implementation
fingerprints. Test copies and source contracts retain their source identities.

Three main-figure data builders and the S1 builder are relocated. Their repository
root discovery/output paths are adapted; the S1 Figure 4 manifest path points to
the prospective curated input. All diffs are recorded in
[`PACKAGING_DIFFS.patch`](../manifests/provenance/PACKAGING_DIFFS.patch).
Those upstream builders are not executed during migration.

[`NOTEBOOK_EXTRACTIONS.json`](../manifests/provenance/NOTEBOOK_EXTRACTIONS.json)
records notebook SHA-256, selected source cell indices/hashes, destination scripts
and textual adaptations. Selected code is retained without executed notebook
outputs or abandoned figure variants. Main Figure 2 uses the usefulness variant;
Figure 3 uses the existing B/C plotting cell; Figure 4B and the source Figure 4C
drawing use distinct selected cells. S1's source display-order issue is preserved
for review rather than silently corrected.

Figure 2's axis wording is changed to benchmark frequency. The source Figure 4C
physical-unit label is changed to native GEM flux units, following the technical
audit. These are display/documentation changes. Scientific computations, endpoints,
parameters, populations and random seeds remain unchanged. The new entry point
performs integrity checks, manages output paths and executes selected plotting
scripts; it does not invoke upstream scientific production.

## Numerical inputs and licensing boundary

No ambiguous frozen numerical artifact is copied. Fifty exact numerical
source identities, four conceptual artwork files and three external-source groups are tracked in
[`RIGHTS_REVIEW.tsv`](../manifests/provenance/RIGHTS_REVIEW.tsv). The selected
prospective main/S1/table-audit subset is in `FIGURE_INPUTS.json` with inclusion
false and rights pending. Raw/model dependencies remain omitted.

The source technical audit reports original Figure 3 pair-table SHA-256
`7a83a31d29317f7b0d906fff73691428d98427d93d2fdebd4ef91294f8a46676`;
the inventory independently records the physical source-file hash. Current A23
manifest identity is `b8e57074f3df23fe7708145bc01f51a19f946545b1c16e61628c2678dbcabd18`.
Do not confuse that manifest hash with its implementation hash or source-audit hash.

[`UPSTREAM_DEPENDENCIES.json`](../manifests/provenance/UPSTREAM_DEPENDENCIES.json)
lists 12 narrowly identified upstream inputs, including the individual pinned cache
under a historical exploratory path. Its containing six-arm directory is excluded.
`MISSING_FOUNDATION_ARTIFACTS.json` identifies 12 absent hash-gated patch files.
No checksum gate is removed or reset.

Third-party datasets, GEMs and software retain their original licenses. Local Dryad
metadata declares CC0-1.0 for the DMI deposit. Other raw/model redistribution rights
and derivative permissions require author confirmation. Nothing was uploaded to an
archive or remote, and no release availability was invented.

## Integrity and limitations

[`SOURCE_COPY_AUDIT.tsv`](../manifests/provenance/SOURCE_COPY_AUDIT.tsv) compares
source identities with destination copies/extractions. Included-file SHA-256
values are in `manifests/checksums/SHA256SUMS.txt`; pending source inputs are not
presented as verified release data. The original root `DATA_PROVENANCE.md` is
preserved as an initial historical note; this document describes current status.

Migration validation covers static compilation, lightweight imports, existing
synthetic scientific-core tests, portable paths, terminology and dependency mapping.
It does not establish fresh model qualification, solver feasibility, redistribution
permission or complete manuscript rendering. See the machine-readable validation
record and [MIGRATION_REPORT.md](MIGRATION_REPORT.md).
