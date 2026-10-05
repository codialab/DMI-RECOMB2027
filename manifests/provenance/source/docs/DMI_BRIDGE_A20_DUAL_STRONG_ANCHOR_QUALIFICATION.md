# DMI-BRIDGE-A2.0 — dual-strong-anchor qualification and freeze

Status: pre-outcome foundation contract.

## Purpose

DMI-BRIDGE-A2 is a post-M1 robustness/generalization extension. It must not modify,
replace, reopen, or retrospectively redefine the frozen one-strong-anchor
DMI-BRIDGE PL1/PL2/M1 evidence.

A2 asks whether the same residual-geometry -> weak-direction-utility mechanism
persists when quantitative DMI information is already supplied through two strong
anchors rather than one.

The two A2 strong-anchor arms are frozen before A2 scientific outcomes are read:

- A2-L: Vmax + Vlac, interpreted through the already-established mappings
  Vmax -> HEX1 and Vlac -> LDH_L.
- A2-G: Vmax + Vglx, interpreted through the already-established mappings
  Vmax -> HEX1 and Vglx -> AKGDm.

Neither arm may be dropped, promoted, tuned, or relabelled based on A2 outcomes.

## Frozen predecessor boundary

Require the repository to contain the frozen manuscript-evidence state

`BRIDGE_M1_MANUSCRIPT_EVIDENCE_FROZEN`

under

`outputs/dmi_bridge_m1_manuscript_evidence_v1/`.

A2 is additive. No Stage-10/11/12, CUP, SCOPE, SHIFT, SA, PL1, PL2A/B/C/D,
or M1 artifact may be edited.

The existing PL2D held-out set has already been opened. A2 therefore must not be
described as a new untouched confirmation experiment. It is a post-freeze
generalization/sensitivity analysis using frozen rules for comparability.

## Canonical source families

The repository adapter must fail closed unless it can verify the literal source
identities and the producing-code semantics.

Known canonical sources include:

1. Canonical DMI mouse input

`pre-simulation_constraint_analysis/results/run_20260824T164247Z_corrected_analysis_A/mouse_level_dmi_input.csv`

expected SHA256:

`ae3cbad6d422405780d0dc44e0e31f95770ebc204d1f80f2ccdada12bb119082`

Expected population: 10 mice, 5 CT2A and 5 GL261.

2. Current 640-candidate panel

`12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/candidate_fraction_table.tsv.xz`

expected SHA256:

`421c2c1d196fbd088a4e92e2ecc998ad578aab814a025fff8e1c843597fed1da`

Expected population: 640 unique candidates =
4 reconstruction methods x 2 tumors x 4 RNA contexts x 20 samples.

3. Historical Stage-12 single-anchor provenance

`12_recomb_method_comparison/outputs_stage12b/dmi_single_anchor_weights.tsv.xz`

expected SHA256:

`b3357587b70329fca7d77b99e4b6c6487feccb54ac7f251d8a0a8437353a8030`

This file is provenance/qualification evidence, not an authority for stale
candidate identities. The current panel is authoritative for A2 candidate
membership.

4. The current Vmax baseline and its qualification receipts created by the
frozen BRIDGE/BIO qualification workflow. The adapter must discover them by
their binding manifest/contract, never by filename freshness, and must pin their
literal hashes in the A2 manifest.

## Coordinate semantics

A2 must use the archived Stage-12 quantitative-DMI rank-kernel semantics, not a
newly invented likelihood.

For every observed DMI coordinate:

- candidate values are converted within each tumor's 320 current candidates to
  an average-rank coordinate `(rank - 0.5) / 320`;
- mouse values are converted within each tumor's five DMI mice to
  `(rank - 1) / 4`;
- exact ties use average ranks;
- candidate identity and row order are frozen before weights are calculated.

For Vmax, the repository adapter must reproduce the already-qualified current
Vmax baseline at numerical tolerance before A2 may proceed.

For Vlac, require the established candidate observable
`lac_production = max(-LDH_L, 0)` and reproduce the historical Stage-12
single-anchor semantics on the identity-overlap subset before rebuilding the
coordinate on the exact current panel.

For Vglx, do not substitute PDHm or a newly calibrated coefficient. Require the
established Vglx -> AKGDm association and recover the exact Stage-12 Vglx
candidate-coordinate semantics from the producing code/provenance. If the
producing semantics cannot be reproduced on the current panel, A2-G is
`NOT_QUALIFIED` and must not be approximated.

## Dual-anchor kernel

The adapter must verify that the archived multi-observable Stage-12 operator is
the equal-weight sum of squared rank-coordinate differences. If verified, for
anchor set A = {a1, a2}:

`D_i = (q_i,a1 - q_mouse,a1)^2 + (q_i,a2 - q_mouse,a2)^2`.

Weights are

`w_i proportional to exp(-(D_i - min(D)) / T)`,

with T solved by the archived deterministic temperature procedure to the same
target ESS used by the one-anchor baseline (20), independently for each mouse.
The verified archived procedure is:

- initialize `lo = 1e-12` and `hi = 1.0`;
- while ESS at `hi` is below target, multiply `hi` by 10, stopping at `1e6`;
- require the target ESS to be bracketed;
- run exactly 80 geometric-mean bisection steps;
- after each step, move `lo` upward when ESS is below target and `hi` downward
  otherwise;
- return the normalized weights evaluated at the final `hi`, not at the final
  geometric mean.

This ESS-matched construction is deliberate: A2 tests the effect of an
additional quantitative coordinate while holding overall finite-panel weight
concentration comparable to A1.

If the archived producer uses materially different multi-observable semantics,
stop with `A20_ARCHIVED_OPERATOR_MISMATCH`; do not silently adopt this formula.

## A2.0 scope

A2.0 is qualification only. It may compute and freeze:

- source hashes and candidate identities;
- current-panel Vmax, Vlac, and Vglx rank coordinates;
- mouse rank coordinates;
- A2-L and A2-G distances, temperatures, and weights;
- per-mouse achieved ESS and numerical reproduction residuals;
- anchor-fit diagnostics;
- exact excluded/failed rows and reasons;
- deterministic manifest and artifact hashes.

A2.0 must not compute or inspect:

- weak-prior correct/wrong/random q outcomes;
- truth-referenced gain or harm;
- PL2B/PL2C/PL2D outcome metrics;
- geometry -> utility correlations;
- reaction ranking by A2 benefit;
- manuscript claim selection.

## Admission

A2-L is qualified only if all of the following hold:

- all current candidate identities are exact and unique;
- Vmax reproduction passes;
- Vlac semantics pass;
- all 10 mouse baselines are finite and normalized;
- target ESS 20 is reproduced within numerical tolerance;
- no candidate support is silently dropped.

A2-G uses the same requirements with Vglx/AKGDm semantics.

One arm may fail qualification without invalidating the other, but failed arms
must remain recorded in the manifest and may not be replaced with a different
second anchor.

## Required outputs

Write only under:

`outputs/dmi_bridge_a20_dual_anchor_qualification_v1/`

with at least:

- `BRIDGEA20_MANIFEST.json`
- `BRIDGEA20_SOURCE_AUDIT.tsv`
- `BRIDGEA20_ANCHOR_REGISTRY.tsv`
- `BRIDGEA20_MOUSE_COORDINATES.tsv`
- `BRIDGEA20_DUAL_ANCHOR_WEIGHTS.tsv.xz`
- `BRIDGEA20_WEIGHT_SUMMARY.tsv`
- `BRIDGEA20_STATUS.json`

Recommended terminal states:

- `BRIDGE_A20_DUAL_ANCHOR_BASELINE_FROZEN`
- `BRIDGE_A20_PARTIAL_QUALIFICATION`
- `BRIDGE_A20_QUALIFICATION_FAILED`

Rerunning an identical frozen qualification must be a no-op.

## Planned downstream stages

Only after A2.0 is frozen:

- A2.1 — dual-anchor residual-geometry landscape.
  Recompute the PL1-style network-wide geometry under A2-L and A2-G, excluding
  the strong anchors themselves from candidate weak target B.

- A2.2 — sign-only weak-prior utility.
  Reuse the frozen PL2 sign-only operator and truth construction on the A2
  baselines. Correct, wrong, randomized, and reliability-varying weak directions
  retain their frozen definitions. Do not retune q, thresholds, or reaction
  subsets.

- A2.3 — A1 vs A2 synthesis.
  Ask whether the geometry-utility relationship and the bridgeable regime persist
  after a second quantitative anchor. Report A2 as a post-freeze robustness
  extension, not a new independent confirmation set.

Stop after A2.3 unless a new manuscript question is explicitly authorized.
