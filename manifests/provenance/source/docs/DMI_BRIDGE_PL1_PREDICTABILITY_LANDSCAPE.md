# DMI-BRIDGE-PL1 — network-wide predictability landscape

PL1 is the first stage of the revised DMI-BRIDGE roadmap. It asks a question that must
be answered before another weak-prior performance experiment: **under the qualified
strong HEX1 anchor, which reaction-B geometries could in principle allow directional
information about B to reduce uncertainty about the magnitude of the CT2A–GL261 B
contrast?**

PL1 is deliberately outcome-blind. It uses the existing qualified strong-anchor flux
ensembles and scans **all admissible finite cache reactions except HEX1**. It does not
preselect LDH/TCA reactions, does not exclude same-pathway or proximal reactions from
the primary analysis, and does not use truth, gene-score concordance, a weak prior,
solver/FVA, new sampling, reconstruction, or new weight fitting. Pathway/proximity is
annotation; exclusion is sensitivity-only.

The primary evaluation unit is an exact reconstruction-method × RNA-context stratum and
one CT2A-mouse × GL261-mouse contrast. Strong-anchor weights are normalized within that
exact stratum. The finite reference uses uniform weights on the **same candidate
support**. Consequently, `anchor_*_contraction` is weight-induced concentration on a
fixed finite support. It is not FVA contraction and must not be described as shrinkage
of the exact feasible set.

For every reaction B, PL1 records: strong-anchor versus uniform SD and q10–q90 width;
positive/tie/negative mass; directional entropy; joint and sign-specific ESS; magnitude
concentration within positive and negative branches; weighted eta-squared and
correlation between sign(B contrast) and |B contrast|; and weighted coupling between
the residual HEX1 contrast and B contrast. Degenerate cases are retained with explicit
status rather than dropped.

No PL1 metric is a bridgeability score. PL1 must not rank or select reaction B and must
not fit a threshold using later performance. Its purpose is to freeze a network-wide,
outcome-blind geometric feature table. Only after PL1 is frozen should PL2 apply a
**sign-only** weak-prior operator with correct, wrong, randomized, and reliability-
varying directions and ask whether the pre-existing PL1 geometry predicts improvement
in truth-referenced |ΔB| error.

## Production execution acceleration

The production adapter may evaluate reactions in bounded vectorized blocks and use a
low deterministic XZ compression preset. These are execution-only changes: the scalar
`target_landscape_metrics()` function remains the scientific authority, all PL1
features and weighted inverse-CDF quantile definitions are unchanged, and deterministic
sentinel reactions are checked against the scalar authority before publication. PL1
must fail closed if batch/scalar parity exceeds the pinned numerical tolerance.

Summary tables distinguish per-metric `finite_count` / `nonfinite_count` from explicit
PL1 degeneracy counts (`target_delta`, `target_abs_magnitude`, `anchor_delta`,
`sign_magnitude`, and `anchor_target_coupling`). A finite metric value is not treated
as evidence that the corresponding reaction geometry is non-degenerate.

## Resumable production execution

The scientific PL1 contract is unchanged by the resumable production adapter
`scripts/dmi_bridge_pl1_predictability_resume_v1.py`.  It exists solely to make
production robust to bounded execution windows.  The adapter checkpoints one of
16 fixed method/RNA contexts only after all 25 mouse-pair evaluations for that
context have passed scalar/batch parity and its 104,500-row feature shard has
been closed and hashed.  A terminated invocation can therefore be rerun with the
same command; completed contexts are hash-verified and skipped, while an
incomplete context is recomputed from its beginning.

Checkpoint arrays and feature shards live in a sibling `*.resume_work_v1`
directory and are not publication artifacts. After all contexts are complete,
the feature rows are emitted as four deterministic files
`BRIDGEPL1_REACTION_EVALUATION_FEATURES.part-000.tsv.xz` through
`part-003.tsv.xz`, with a header in every part. The adjacent
`BRIDGEPL1_REACTION_EVALUATION_FEATURES.parts.json` manifest defines canonical
numeric ordering, row counts, schema, hashes, and the logical full-table digest.
Readers treat these parts as one logical feature artifact. Historical monolithic
`BRIDGEPL1_REACTION_EVALUATION_FEATURES.tsv.xz` files remain read-only fallback
inputs; writers never create them, and conflicting old/new forms fail closed.
The resumable metrics checkpoint is canonically `metrics.npy.xz`. Because XZ
cannot be memory-mapped directly, resumable execution uses a controlled external
temporary NPY mmap and checkpoints it back to deterministic XZ after each
completed context. Legacy `metrics.npy` is accepted only as a verified fallback.
Global summaries are computed one metric at a time from the checkpointed memory
maps so that the previous approximately 0.5-GB advanced-index copy is not
materialized.

No truth, weak prior, gene score, solver, FVA, sampling, reconstruction, target
selection, threshold, regime label, or scientific metric is introduced or
changed by resumability.  The final manifest records the checkpoint execution
profile and hashes all published artifacts.  PL2 remains `NOT_RUN`.
