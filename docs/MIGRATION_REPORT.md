# Migration report

Migration is in progress. The final audit will replace this status with a supported
reproducibility verdict.

## Source state

Source commit: `16379a169f45ec556274953eb66ae12293e1e985`, clean `main` checkout.
The source is strictly read-only. Inspection uses `GIT_OPTIONAL_LOCKS=0` for source
Git operations. No source scripts, tests, reconstruction, optimization or sampling
are executed.

## Destination state

Initial commit: `63edf6a9f3e0f3ba84c17dd7a828db5b6d7170de`, clean `main` checkout
containing `DATA_PROVENANCE.md`. That historical document is preserved.
Work proceeds on `publication-migration`; no remote push is authorized.

## Unresolved issues found during inventory

- The supplement draft assigns S2 to context associations and S3 to alternative
  descriptors. Existing S2/S3 workflows instead show paired geometry. Preserve the
  disagreement for author review.
- Legacy S4, S5, S6 and S7 workflows also require panel-by-panel comparison; their
  presence and prior audits do not prove draft-caption coverage.
- Figure 3D's legacy SVG reports case-pooled endpoints, superseded by the final
  equal-weight reaction–evaluation means. It must not be shipped as final artwork.
- Figure 4C's draft caption omits the first stage of the top-100 screen. The current
  builder/manifest records a first-stage geometry ranking followed by the selected
  example's cross-method entropy-shift range. Preserve that discrepancy.
- The current Figure 4 notebook displays four methods, while the main draft shows
  three and assigns the degenerate RIPTiDe counterpart to S7. Final assembly needs
  author review.
- Gain and flux axes require native/model flux units; no DMI-to-GEM physical-unit
  conversion is established by the frozen rank operators.
- Frozen derived-table redistribution requires the review recorded in
  `RIGHTS_REVIEW.tsv`. No pending table has been copied.
- Several production adapters require root-level foundation patch files not found
  during inspection. Their hash gates must not be bypassed to obtain a successful run.

## Code modifications

Scientific changes: none. Packaging changes and source-to-destination comparisons
will be recorded before the final audit.
