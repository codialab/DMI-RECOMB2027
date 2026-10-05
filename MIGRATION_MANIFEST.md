# Migration manifest

Source: `/home/pty/work/project_MRI_brain_tumor/`, commit
`16379a169f45ec556274953eb66ae12293e1e985`.
Destination initial commit: `63edf6a9f3e0f3ba84c17dd7a828db5b6d7170de`.
Migration branch: `publication-migration`. Migration date: 2026-10-05.

The source is strictly read-only, including its Git state, ignored files, caches,
and generated outputs. All persistent migration writes occur in the destination.
Files are copied, never moved. No expensive scientific workflow is authorized.

## Classification contract

Every important source path receives exactly one of these classifications:
`INCLUDE_PRODUCTION`, `INCLUDE_PROVENANCE_ONLY`, `EXCLUDE_EXPLORATORY`,
`EXCLUDE_SUPERSEDED`, or `EXCLUDE_LARGE_OR_EXTERNAL`.

The final file-level decisions, destinations, reasons, authority status and copy
status are recorded in [SELECTION.tsv](manifests/provenance/SELECTION.tsv).
Directory exclusions are recorded in that same table; they are not bulk-copy rules.
Redistribution decisions are separate from scientific relevance and are recorded
in [RIGHTS_REVIEW.tsv](manifests/provenance/RIGHTS_REVIEW.tsv).

## Authority and dependencies

The current manuscript drafts and manuscript technical audit govern figure
identity, numbering, panel composition and scientific message. A frozen manifest
establishes artifact identity, not automatic agreement with the manuscript.
Historical filenames are retained for provenance; public terminology is A2.

The Bridge scripts retain their original `scripts/` paths because imports and
frozen implementation fingerprints name those paths. Stage-11 retains its native
package hierarchy to preserve relative imports, cache identities and root discovery.
This is a deliberate variation from the illustrative `src/` layout, avoiding an
unnecessary reorganization of validated scientific modules.

Figure workflows may be moved under `figures/main/` or `figures/supplementary/`.
Original source paths and any packaging edits are recorded individually. Scientific
defaults, equations, populations, seeds, tolerances and estimators are preserved.

## Pending author decisions

Ambiguous research products are withheld until redistribution is confirmed.
Legacy S2/S3 numbering conflicts with the current supplement draft; those artifacts
are not promoted into the manuscript workflow. Figure 3D contains superseded
case-pooled labels. See [RIGHTS_REVIEW.md](docs/RIGHTS_REVIEW.md) and the final
migration report for author-review items.
