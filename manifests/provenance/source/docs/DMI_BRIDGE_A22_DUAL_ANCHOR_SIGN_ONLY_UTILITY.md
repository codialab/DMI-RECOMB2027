# DMI-BRIDGE-A2.2 — matched sign-only weak-prior utility after two strong anchors

Status: pre-outcome foundation contract.

## Purpose

DMI-BRIDGE-A2.2 tests whether the already-frozen PL2 sign-only weak-direction
operator retains utility after the qualified dual-strong-anchor A2-L baseline is
present.  A2.2 is a post-M1 robustness/generalization experiment.  It is not a
new untouched confirmation set because the original PL2D holdout has already
been opened.

A2.2 consumes only the qualified A2-L arm:

- Vmax -> HEX1;
- Vlac -> LDH_L using `lac_production = max(-LDH_L, 0)`;
- A2-G remains `NOT_QUALIFIED` and is not repaired, substituted, or used.

A2.2 must stop before any A1-versus-A2 geometry/utility synthesis.  That belongs
to A2.3.

## Frozen predecessor gates

Require:

1. A2.0 terminal status `BRIDGE_A20_PARTIAL_QUALIFICATION`, with A2-L
   `QUALIFIED`, A2-G `NOT_QUALIFIED`, and no A2-G outcome access.
2. A2.1 terminal status `BRIDGE_A21_DUAL_ANCHOR_GEOMETRY_FROZEN` and
   `scientific_outcomes_computed = false`.
3. The frozen original PL2A truth-selection bundle and its exact identity
   contract.
4. The frozen original PL2 sign-only numerical core and vectorized batch
   semantics.

All bound files and logical-content identities must be hash verified before any
reaction-B truth value is read.

## Preserve the original held-out truth identities

A2.2 must not select new truth candidates from A2-L weights.

The original PL2A truth construction is frozen and remains the benchmark:

- 480 quantile-selection aliases;
- 1,200 truth-pair aliases;
- 885 distinct selected truth pairs under
  `(evaluation_id, ct2a_candidate_id, gl261_candidate_id)`;
- 836 distinct truth pairs that were evaluable in the original A1 PL2B
  benchmark;
- 49 distinct truth pairs that were A1-non-evaluable.

The **primary A2.2 matched population is exactly the 836 original A1-evaluable
distinct truth pairs**.  This restriction is frozen before A2.2 outcomes so that
A1 and A2 later compare the same held-out benchmark population.

The 49 original A1-non-evaluable pairs remain in the A2.2 support audit but are
not promoted into primary A2.2 outcome production even if A2-L happens to make
them evaluable.  No replacement truth pair, pseudomass, or newly selected truth
candidate is allowed.

A2.2 must preserve the original canonical `truth_pair_id` construction and all
1,200 aliases for audit.  Quantile aliases never duplicate primary scientific
weight.

## Recompute post-holdout support under A2-L

For each of the 836 matched truth pairs:

1. locate the exact frozen A2-L 20-candidate weight vector for each condition,
   mouse, reconstruction method, and RNA context;
2. verify the selected truth candidate identity/sample index against the frozen
   PL2A identity;
3. remove that candidate from its own 20-candidate inference pool;
4. renormalize the remaining 19 A2-L weights separately in CT2A and GL261;
5. record post-holdout mass, positive-weight count, and ESS for each condition.

A2.2 has no new ESS threshold.  Low conditioned ESS is a property of the
qualified dual-anchor baseline, not an exclusion criterion.  A pair is admitted
whenever both post-holdout distributions retain positive finite mass.

For strict A1/A2 matched production, all 836 original A1-evaluable pairs must
remain A2-L evaluable.  If any fail this support gate, stop before reaction-B
outcomes with a support-mismatch status; do not silently shrink or replace the
primary population.

The 49 A1-non-evaluable pairs may have their A2-L support status recorded
outcome-blind, but no reaction-B utility outcomes are computed for them in A2.2
v1.

## Common reaction universe

The qualified cache has 4,181 reactions.  A2-L uses two strong anchors, so
neither may be weak target B:

- exclude `HEX1`;
- exclude `LDH_L`.

Expected A2.2 target count: 4,179 reactions.

This is the exact common A1/A2 reaction universe for later matched synthesis:
A1 PL2B contains one additional target (`LDH_L`) that must not enter A1-versus-A2
comparisons.

Do not preselect favorable pathways or remove same-subsystem/proximal reactions
from the primary scan.  Pathway exclusions remain sensitivity-only.

## Frozen weak-direction operator

Reuse the original PL2 operator without scientific modification.

For reaction B and a CT2A/GL261 candidate pair:

`delta_B = B_CT2A - B_GL261`.

For supplied direction `d in {-1,+1}`:

- score `+1` if `d * delta_B > 1e-12`;
- score `-1` if `d * delta_B < -1e-12`;
- score `0` otherwise.

Posterior candidate-pair weights are proportional to

`p0(pair) * exp(lambda_total * score)`.

Freeze exactly:

- `TIE_TOL = 1e-12`;
- `GAIN_TOL = 1e-12`;
- lambda grid `{0, 0.25, 0.5, 1.0}`;
- primary lambda `0.25`;
- reliability grid `{0, 0.25, 0.5, 0.75, 1}`.

The primary estimate remains `E[|delta_B|]`; the primary error is absolute error
to the held-out `|delta_B_truth|`; positive gain means lower absolute error than
the frozen A2-L baseline.

For every non-tie truth case compute the correct direction and its exact wrong
opposite.  The reliability curve is analytic:

`g(q) = q*g_correct + (1-q)*g_wrong`.

`q=0.5` is the exact expected random-sign null.  Do not introduce stochastic
Monte Carlo random-sign controls or retune the reliability grid.

Truth ties remain `TRUTH_TIE` and receive no fabricated correct/wrong direction.
Lambda zero is QC-only and must recover the A2-L held-out baseline.

## Outcome storage and expected accounting

If all 836 matched truth pairs pass A2-L post-holdout support, expected primary
case rows are:

`836 * 4179 = 3,493,644`.

Compute each distinct truth pair exactly once per reaction.  Do not materialize
five reliability-q rows per case.  Store correct/wrong endpoint gains, exact
random-sign expected gain, and break-even reliability as in original PL2B.

A2.2 case rows must preserve the original `truth_pair_id`, `evaluation_id`,
algorithm, RNA context, mouse pair, reaction ID, and pathway annotations needed
for later exact A1/A2 joins.

A2.2 must not duplicate A2.1 geometry features into each outcome row.  Bind the
A2.1 logical feature identity and join later by `(evaluation_id, reaction_id)`.

## Permitted summaries

A2.2 may publish the same fixed descriptive utility summaries as PL2B:

- primary lambda 0.25 correct/wrong/random-sign expected gains;
- improved/tied/harmed fractions;
- gain quantiles;
- frozen reliability-grid summaries derived analytically;
- lambda 0.5 and 1.0 sensitivities;
- exact method × RNA-context summaries;
- truth-tie and finite-support accounting.

Primary summaries weight each distinct truth pair once.  Alias-weighted summaries
are sensitivity-only if retained.

A2.2 must not relate utility to A2.1 geometry, fit a bridgeability score, select
predictors, optimize a threshold, rank reactions, or perform A1-versus-A2
synthesis.  Those are A2.3 questions.

## Required outputs

Write only under:

`outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/`

with at least:

- `BRIDGEA22_SOURCE_AUDIT.json`
- `BRIDGEA22_ANALYSIS_CONTRACT.json`
- `BRIDGEA22_TRUTH_SUPPORT_AUDIT.tsv`
- `BRIDGEA22_MATCHED_TRUTH_PAIR_REGISTRY.tsv`
- `BRIDGEA22_ALIAS_MAP.tsv`
- `BRIDGEA22_REACTION_REGISTRY.tsv`
- partitioned `BRIDGEA22_CASE_OUTCOMES` plus a parts manifest
- `BRIDGEA22_CONTEXT_SUMMARY.tsv`
- `BRIDGEA22_QC.json`
- `BRIDGEA22_MANIFEST.json`

Use deterministic partitioning with explicit row counts, byte sizes, per-part
SHA256 values, schema, and a logical-content SHA256.  An identical completed
rerun must be a no-op; conflicting existing outputs fail closed.

Desired terminal state:

`BRIDGE_A22_SIGN_ONLY_UTILITY_FROZEN`.

A support failure before outcome access must be recorded separately and must not
be repaired by changing truth identities or support rules.

## Forbidden work

A2.2 must not:

- select new held-out truths;
- use the 49 A1-non-evaluable truth pairs as primary outcome cases;
- refit or retune A2-L weights;
- impose a new ESS cutoff;
- use A2-G or substitute AKGDm/PYRt2m/PDHm as a second strong anchor;
- run solver/optimization, FVA, sampling, reconstruction, or new strong-anchor
  fitting;
- use gene scores;
- compute geometry-to-utility associations;
- inspect outcomes to choose reactions, pathways, lambda, q, thresholds, or
  predictor subsets;
- perform A1-versus-A2 synthesis;
- claim untouched confirmation.

Stop after A2.2.  A2.3 requires separate authorization.
