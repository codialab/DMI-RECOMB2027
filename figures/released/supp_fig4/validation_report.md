# S4 validation report

## Inputs and provenance

The builder uses the production PL2C development bundle, PL2D-H development snapshot, and PL2D confirmation bundle. It checks the pinned PL2C manifest SHA-256 `6abcf149cd6adb294c0931585cec28f992b453451a5b5f03e10b3ffb2f9b9c39`, PL2D-H manifest SHA-256 `a1dd76d6b40e741e07efe64c34bc39eaad289a46bf6c870917a4f59cab688427`, and PL2D manifest SHA-256 `9f9eb730473e7ad209c1d5ec3340808ab614087b8ae1121d929eff61be12b699`, then checks every artifact declared by those manifests. Exact paths and checksums are in `data/supp_fig4_validation_summary.json` and `data/supp_fig4_build_manifest.json`.

The builder independently aggregates the frozen PL2C evaluation × reaction data by reaction and recomputes development Spearman rho with the frozen SciPy implementation. It checks the PL2D reaction table against the frozen holdout registry and recomputes confirmation rho using round-trip float parsing. It also cross-checks the frozen M1 development/confirmation table.

## Dimensions and panel denominators

- PL2C development split: 3,135 reactions; 1,159,950 evaluation × reaction rows; 2,455 finite reaction-level predictor/response pairs; Spearman ρ = 0.8570601823.
- PL2D held-out confirmation: 1,045 reaction rows; 838 finite pairs; Spearman ρ = 0.8351435823.
- Confirmation contexts: 16 method × RNA-context rows, all evaluable and positive; finite context counts range from 34 to 638.
- Confirmation bootstrap: 5,000 finite reaction-level replicates.

Context correlations range from ρ = 0.2756103148 for RIPTiDe in `setx2 + setx3` (n=39) to ρ = 0.9255128993 for RIPTiDe in `setx1 + setx2 + setx3` (n=37). The low-support context remains in the figure and prevents implying uniform performance.

Development and confirmation reaction identifiers are disjoint. No confirmation observations enter the development estimate.

## Integrity and exclusions

Reaction identifiers are unique in the frozen split and reaction-level outputs; PL2C evaluation × reaction keys are unique. Spearman correlations use pairwise finite predictor/response values. This excludes 680 development reactions and 207 confirmation reactions from their respective correlation estimates because at least one value is nonfinite. No bootstrap replicate is excluded. Context results retain all 16 frozen rows and pass the minimum support gate of 30 finite pairs.

Round-trip float parsing is used for serialized PL2D reaction values so the stored numerical values and rank tie structure reproduce the frozen confirmation correlation exactly. Default pandas float parsing shifted rho slightly and is not used.

All 16 rows in the PL2D context table match the corresponding context records embedded in the frozen primary result by method, RNA context, finite count, and rho.

The notebook's Python cells were executed in order directly with Python, producing both image files. `nbconvert --execute` could not start a Jupyter kernel because the sandbox denies local socket creation; the direct execution exercised the same notebook code without kernel startup. Both PNG and SVG were rendered and visually inspected.

## Bootstrap definition

The table contains the frozen 5,000 PL2D bootstrap replicates (NumPy Generator PCG64, seed 20260929), each resampling 1,045 reaction indices and filtering to finite pairs within the replicate. The percentile interval uses NumPy quantile with `method="linear"`: median 0.8354193419; 95% interval [0.8015629119, 0.8631275423]. This is computational reaction-level stability, not biological replication.

## Scientific validation and issues

Panel definitions follow the S4 contract. Development and confirmation are shown separately; the 16 context correlations include the weakest context; no exploratory alternative-anchor results enter the figure. The displayed S4 panels do not reproduce S1–S3 diagnostics.

- **BLOCKING:** None.
- **IMPORTANT:** The reusable M1 `verify_sources()` helper could not run because its expected foundation patch file is absent from this worktree. The S4 builder instead directly verified the pinned PL2C, PL2D-H, and PL2D manifests and all artifacts listed in them. This does not block S4 provenance checks.
- **COSMETIC:** None identified.
