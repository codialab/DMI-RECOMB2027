# DMI-BRIDGE-A2.1 — residual geometry after the qualified dual strong anchor

Status: pre-outcome foundation contract.

## Purpose

A2.1 is the first scientific stage after A2.0. It asks whether the network-wide
residual geometry that made sign-only information useful under the original
one-anchor baseline persists after a second **qualified quantitative strong
measurement** has already been incorporated.

A2.0 ended at `BRIDGE_A20_PARTIAL_QUALIFICATION`:

- A2-L (`Vmax + Vlac` -> `HEX1 + LDH_L`) is `QUALIFIED` on all 10 mice;
- A2-G (`Vmax + Vglx`) is `NOT_QUALIFIED` because the archived quantitative Vglx
  operator used PYRt2m while the current `Vglx -> AKGDm` mapping is only a
  sensitivity association without an independently qualified current-panel
  quantitative Vglx operator.

Therefore **A2.1 uses A2-L only**. A2-G must remain present in provenance as a
failed qualification arm, but it must not be approximated, repaired, substituted,
or silently replaced by AKGDm, PYRt2m, PDHm, or another reaction.

A2.1 is outcome-blind. It does not apply a weak prior, use truth, read PL2 utility
outcomes, use gene scores, fit new weights, run optimization/FVA/sampling, or rank
reactions.

## Frozen A2.0 admission

The adapter must fail closed unless it verifies at minimum:

- `outputs/dmi_bridge_a20_dual_anchor_qualification_v1/BRIDGEA20_MANIFEST.json`
  SHA256 `01f8004a0bf9d5358d7ce7281e827cc1a1fbe00a6750a5c063a1edad2eedf380`;
- manifest status `BRIDGE_A20_PARTIAL_QUALIFICATION`;
- A2-L status `QUALIFIED` and A2-G status `NOT_QUALIFIED`;
- `BRIDGEA20_DUAL_ANCHOR_WEIGHTS.tsv.xz`
  SHA256 `b153e2fde45dca9573c71c548783e704313f056d0980f64af60180ad90ebc1f4`;
- exactly 3,200 A2-L weight rows, 10 mouse summaries, and no A2-G weight rows;
- current panel ID
  `9c27c8e75aa128a4a4cb15547f23557845f37d35a7323e875af5a1711a527f42`;
- A2.0 records `scientific_outcomes_computed=false` and
  `a2_g_outcome_accessed=false`.

The adapter must also verify the A2.0 source audit and all bound current-panel,
flux-cache, and mapping identities rather than trusting filename freshness.

## Reuse the PL1 scientific definitions

A2.1 is deliberately a **PL1-style rerun under different frozen weights**, not a
new geometry framework.

Use the existing PL1 scalar definitions as the scientific authority. Preserve the
same CT2A-minus-GL261 contrast orientation, tie tolerance, weighted inverse-CDF
quantiles, sign-state definitions, conditional ESS calculations, eta-squared,
directional entropy, and uniform-on-identical-support reference.

For A1-versus-A2 comparability, preserve the original eight PL2C predictor
semantics:

- `sign_magnitude_eta2`;
- `sign_magnitude_correlation`;
- `directional_entropy3`;
- `dominant_sign_mass`;
- `n_supported_sign_states`;
- `anchor_width80_contraction`;
- `joint_ess`;
- `abs_delta_anchor_target_correlation`.

The last field continues to mean **absolute HEX1-contrast versus absolute
reaction-B-contrast coupling**, because HEX1 is the common strong anchor present in
both A1 and A2. Do not silently redefine this field as a composite two-anchor
quantity.

A2.1 may additionally report separately named lactate-anchor coupling diagnostics:

- `delta_lactate_target_correlation`;
- `abs_delta_lactate_target_correlation`;
- lactate contrast mean/SD and degeneracy/status fields.

The lactate anchor observable is exactly the qualified A2-L candidate observable
`lac_production = max(-LDH_L, 0)`.

Do **not** create a composite dual-anchor bridgeability/coupling score in A2.1.

## Weight and support semantics

A2.0 freezes 320 candidate weights per mouse over the tumor-matched current panel.
For each exact reconstruction-method x RNA-context stratum, take the corresponding
20 candidates and normalize the frozen A2-L weights **within that exact stratum**,
matching PL1's conditioning semantics. Do not refit the A2-L kernel or retune ESS.

For every CT2A-mouse x GL261-mouse pair within a matched method/RNA context, form
the 20 x 20 Cartesian contrast support. The uniform reference uses uniform weights
on those exact same 20 + 20 candidate sets. Thus `anchor_*_contraction` remains
weight-induced concentration on fixed finite support; it is not FVA or exact
feasible-set contraction.

## Reaction universe

The current cache contains 4,181 reactions. A2.1 excludes both reactions used as
qualified strong anchors from candidate weak target B:

- `HEX1`;
- `LDH_L`.

If both are uniquely present and all other admission gates pass, the expected
primary reaction universe is therefore **4,179 reactions**.

Do not preselect LDH/TCA reactions or favorable pathways. Same-pathway/proximal
status remains annotation and sensitivity-only. Add explicit annotations for
relation to HEX1 and LDH_L while retaining the historical proximal-set annotation
for comparability.

Expected full accounting if the frozen 16 method/RNA contexts and 25 mouse pairs
per context reproduce:

- evaluation registry: `16 x 25 = 400` rows;
- reaction-evaluation features: `400 x 4,179 = 1,671,600` rows;
- context-reaction summary: `16 x 4,179 = 66,864` rows;
- reaction summary: `4,179` rows;
- pathway summary: dynamic from the frozen subsystem annotations.

These are admission expectations, not numbers to force if an identity gate fails.

## A1/A2 comparability boundary

A2.1 itself does not perform the A1-versus-A2 synthesis. It only freezes the A2-L
geometry table.

Later A2.3 may compare A1 and A2 on the **4,179-reaction intersection** by excluding
LDH_L from the historical A1 table. A2.1 must not rewrite or filter historical PL1
outputs to manufacture that intersection.

The primary mechanistic predictor already confirmed in the original BRIDGE branch,
`sign_magnitude_eta2` (fraction of magnitude variance explained by direction), is
computed here with the same definition but under A2-L weights. This is not itself
a test of utility; utility remains A2.2.

## Required implementation behavior

Prefer reuse of the existing PL1 vectorized/scalar-validated and resumable execution
machinery rather than a scientifically different reimplementation. Any new batch
path must reproduce the scalar A2.1 authority on deterministic sentinels at pinned
numerical tolerances.

Because the feature table is large, production should be resumable by the same
fixed 16 method/RNA contexts: publish a context checkpoint only after all 25 mouse
pairs and all admitted reactions for that context pass validation. Incomplete
contexts are recomputed from their beginning. Publication artifacts must be
partitioned/deterministic rather than relying on an in-memory monolith.

## Required outputs

Write only under a new versioned directory, preferably:

`outputs/dmi_bridge_a21_dual_anchor_geometry_v1/`

The final implementation should publish at least:

- `BRIDGEA21_SOURCE_AUDIT.json`;
- `BRIDGEA21_ANALYSIS_CONTRACT.json`;
- `BRIDGEA21_EVALUATION_REGISTRY.tsv`;
- partitioned `BRIDGEA21_REACTION_EVALUATION_FEATURES` plus parts manifest;
- `BRIDGEA21_CONTEXT_REACTION_SUMMARY.tsv.xz`;
- `BRIDGEA21_REACTION_SUMMARY.tsv.xz`;
- `BRIDGEA21_PATHWAY_SUMMARY.tsv`;
- `BRIDGEA21_SENSITIVITY_SUMMARY.tsv`;
- `BRIDGEA21_QC.json`;
- `BRIDGEA21_MANIFEST.json`.

An identical completed rerun must be a no-op; conflicting pre-existing outputs
must fail closed.

## Prohibited work and terminal state

A2.1 must not compute or inspect A2.2 correct/wrong/random/reliability outcomes,
truth-referenced error, geometry-to-utility correlations, gene-score concordance,
reaction rankings, new thresholds, or manuscript inferential statistics.

A2.1 may terminate as

`BRIDGE_A21_DUAL_ANCHOR_GEOMETRY_FROZEN`

only after all source, support, scalar/batch parity, row-accounting, determinism,
and no-forbidden-work gates pass.

Stop after A2.1. Do not automatically launch A2.2.
