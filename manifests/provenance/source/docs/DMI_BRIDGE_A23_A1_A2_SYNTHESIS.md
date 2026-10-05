# DMI-BRIDGE-A2.3 — matched A1 versus A2 synthesis

Status: pre-synthesis foundation contract.

## Purpose

A2.3 is the terminal post-freeze synthesis stage for the dual-strong-anchor
extension. It asks whether the geometry -> weak-direction-utility relationship
established in the frozen A1 BRIDGE branch persists, attenuates, or changes after
the qualified second quantitative anchor `Vlac -> LDH_L` is added.

A2.3 performs no new metabolic simulation and no new weak-prior production. It
consumes only frozen A1 and A2 artifacts.

A2 is not an untouched confirmation set. The PL2D holdout has already been
opened. A2.3 is a matched post-freeze robustness/generalization analysis.

## Required frozen predecessors

A2.3 must fail closed unless it verifies:

- A2.0 terminal state `BRIDGE_A20_PARTIAL_QUALIFICATION`;
- A2-L qualified and A2-G not qualified;
- A2.1 terminal state `BRIDGE_A21_DUAL_ANCHOR_GEOMETRY_FROZEN`;
- A2.2 terminal state `BRIDGE_A22_SIGN_ONLY_UTILITY_FROZEN`;
- A2.2 QC PASS and exact frozen A2.2 logical case identity;
- the frozen A1 PL2A/PL2B/PL2C/PL2D artifacts and the manuscript-evidence M1
  state used by the final one-anchor BRIDGE branch.

A2.3 must discover A1 artifacts through their binding manifests/contracts, never
through filename freshness.

## Matched comparison population

The primary A1-versus-A2 comparison is restricted before synthesis outcomes are
read.

Use exactly:

- the 836 distinct truth pairs that were A1-evaluable and were reused by A2.2;
- the 4,179 reaction targets common to A1 and A2, equal to the original A1
  non-HEX1 target universe with `LDH_L` removed because it is a strong anchor in
  A2-L;
- exact frozen truth-pair and reaction identities.

Do not add the 49 A1-non-evaluable pairs merely because A2-L supports them.

Do not compare A1 `LDH_L` weak-target outcomes with A2, because `LDH_L` is a
strong anchor in A2-L.

The expected fully matched case-key population is therefore
`836 * 4179 = 3,493,644` before any estimator-specific tie or finite-value
filtering required by the frozen A1 contract.

## No new tuning

A2.3 must not tune or redefine:

- lambda;
- reliability q;
- gain tolerance;
- truth construction;
- reaction subset;
- geometry predictor;
- development/confirmation split;
- pathway exclusions;
- ESS thresholds;
- bridgeability thresholds.

The primary weak-prior arm remains the frozen A1 primary arm, expected to be
`lambda_total = 0.25` and correct direction / reliability `q = 1`, but the
repository adapter must verify the literal A1 contract before computing
synthesis results.

The wrong-direction and analytic `q = 0.5` random-sign arms remain controls.
Sensitivity lambdas/reliabilities retain their frozen roles.

## A1-comparable geometry

For the common reaction/evaluation population, preserve the original PL1/PL2C
meanings of the A1-comparable geometry features:

- `sign_magnitude_eta2`;
- `sign_magnitude_correlation`;
- `directional_entropy3`;
- `dominant_sign_mass`;
- `n_supported_sign_states`;
- `anchor_width80_contraction`;
- `joint_ess`;
- `abs_delta_anchor_target_correlation`.

In both A1 and A2, `abs_delta_anchor_target_correlation` continues to refer to
HEX1-target coupling. Do not replace it with a two-anchor composite.

A2-only lactate-target coupling features may be reported as supportive
descriptors, but they are not substitutes for the frozen A1 primary predictor
and must not be used to retune the A1 development/confirmation rule.

## Frozen PL2C/PL2D estimator contract

Before computing synthesis outcomes, the adapter must recover and bind the exact
frozen PL2C/PL2D definitions for:

- utility endpoint;
- case/evaluation/reaction aggregation level;
- truth-tie handling;
- finite/evaluable filters;
- primary geometry predictor;
- development and confirmation reaction identities;
- `rho_dev` and `rho_confirm` computation;
- any bootstrap, interval, permutation, or other uncertainty procedure;
- any proximal/pathway sensitivity definition.

If any of these cannot be recovered exactly, stop with an explicit provenance
blocker. Do not recreate them from memory.

A2.3 applies the same estimator to A2-L. It does not design a new predictor.

## Required synthesis questions

Once the frozen A1 estimator is verified, A2.3 must answer, on matched identities:

1. Geometry shift:
   how the A1-comparable residual-geometry descriptors change from A1 to A2-L.

2. Utility shift:
   how the frozen primary correct-direction utility changes from A1 to A2-L.

3. Directional specificity:
   how the contrast between correct and wrong direction changes from A1 to A2-L.
   The analytic q=0.5 arm remains the random-sign reference.

4. Geometry -> utility persistence:
   apply the exact frozen PL2C/PL2D estimator to A2-L and report the A2 analogues
   of the frozen development/confirmation quantities without changing the split
   or predictor.

5. Matched interpretation:
   distinguish:
   - already determined / concentrated residual geometry;
   - remaining directional ambiguity;
   - weak-direction utility;
   - finite-support/ESS limitations.

A2.3 must not turn these questions into a new reaction-ranking or classifier
search.

## Statistical scope

Use the exact frozen A1 uncertainty/statistical procedure where a matched A2
analogue is required. Do not introduce new p-values, thresholds, multiplicity
procedures, or outcome-selected subgroups unless a separately documented
sensitivity analysis is explicitly authorized.

Any direct A1-minus-A2 descriptive quantity must be paired on exact frozen
identities. The adapter must state the unit of analysis and denominator for every
summary.

## Required outputs

Write only under:

`outputs/dmi_bridge_a23_a1_a2_synthesis_v1/`

with at least:

- `BRIDGEA23_SOURCE_AUDIT.json`
- `BRIDGEA23_ANALYSIS_CONTRACT.json`
- `BRIDGEA23_MATCHED_POPULATION_AUDIT.tsv`
- `BRIDGEA23_GEOMETRY_COMPARISON.tsv`
- `BRIDGEA23_UTILITY_COMPARISON.tsv`
- `BRIDGEA23_GEOMETRY_UTILITY_PERSISTENCE.tsv`
- `BRIDGEA23_CONTROL_SUMMARY.tsv`
- `BRIDGEA23_SENSITIVITY_SUMMARY.tsv`
- `BRIDGEA23_QC.json`
- `BRIDGEA23_MANIFEST.json`

The exact row counts for estimator-dependent outputs must be frozen only after
the read-only audit reconstructs the A1 PL2C/PL2D contract.

Preferred terminal state:

`BRIDGE_A23_A1_A2_SYNTHESIS_FROZEN`

An identical completed rerun must be a no-op. Conflicting existing outputs must
fail closed.

## Forbidden work

A2.3 must not:

- generate new A1 or A2 candidate weights;
- select new truths;
- run PL2B/A2.2 weak-prior production;
- run optimization, FVA, reconstruction, or sampling;
- use A2-G;
- replace the frozen dev/confirm split;
- tune lambda, q, thresholds, predictors, or reaction subsets;
- select reactions or pathways because of A2 outcomes;
- call A2 an untouched confirmation experiment.

Stop after A2.3. Any further biological application or manuscript expansion
requires separate authorization.
