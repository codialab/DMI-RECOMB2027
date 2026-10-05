# DMI-BRIDGE-PL2C — frozen PL1 geometry to PL2B utility development contract

PL2C is the first stage allowed to relate the frozen, outcome-blind PL1 geometry
features to the frozen PL2B sign-only utility outcomes.  It is a **development**
stage, not the final confirmation stage.  Before any geometry-to-utility result is
computed, PL2C freezes a reaction-level development/confirmation split so that a
later PL2D can test the resulting rule or hypothesis on reactions that were not
used to choose it.

PL2C does not run a metabolic solver, perform FVA or sampling, reconstruct a
model, refit the HEX1 anchor, generate new flux weights, or change PL1/PL2A/PL2B.
It consumes only already-qualified artifacts.

## Frozen admission

PL2C must fail closed unless all upstream identities verify.  In particular:

- PL2A manifest SHA256 is
  `8ae3eb60eb23614e379e38e3aa3cb72ec16c6bdb69aa9c74cc94499e31f1c312`.
- PL2B manifest SHA256 is
  `7d9ca97a3f4ea5aa8b0b983e09ee83fc3c82f326d1b1cbf9c6137cfadf25539d`.
- PL2B status is `PL2B_SIGN_ONLY_UTILITY_COMPLETE`.
- PL2B case logical-content SHA256 is
  `c3670992a0ad33698af6c7620429bed48e0ecaa98adb3c87e648eefb1ba70363`.
- The frozen PL1 reaction-evaluation feature logical-content SHA256 remains
  `46aa9f0410f89d03b020e01fefa517348faae1d55d9bfe8387ac04b39bb9ac93`.
- PL2B still contains 885 distinct truth pairs, 836 evaluable pairs, 49
  non-evaluable pairs, 4,180 reactions, and exactly 3,494,480 compact case rows.

Every source artifact consumed by PL2C must be verified against the hashes already
recorded in the qualified upstream manifests before analysis.

## Freeze the confirmation reactions before outcome analysis

The primary network-wide unit of generalization is reaction B.  PL2C therefore
freezes a reaction-level holdout using only the 4,180 frozen non-HEX1 reaction
identities.  No PL2B gain, truth, baseline-error, or other outcome value may enter
this split.

For each reaction ID compute

`SHA256("DMI-BRIDGE-PL2C-REACTION-HOLDOUT-V1\\0" + reaction_id)`.

Sort reactions by `(hash, reaction_id)`.  The first 1,045 reactions are the
`CONFIRMATION_HOLDOUT`; the remaining 3,135 are `DEVELOPMENT`.  The complete split
registry and its digest are published before any geometry-to-utility summaries.
Split assignment must be invariant to input ordering.

PL2C may read PL2B case files sequentially for storage reasons, but confirmation
rows must be rejected immediately by reaction ID and their outcome fields must not
be converted, summarized, fitted, ranked, or emitted into a PL2C development
artifact.  The only permitted PL2C information about confirmation reactions is
outcome-blind identity/annotation needed to audit the split.  PL2D is the first
stage permitted to analyze confirmation outcomes.

## Primary scientific question

The question is not simply whether a correct sign arm has positive gain.  Random
sign perturbation can also improve some cases.  PL2C therefore separates:

1. **correct-sign utility versus the strong-anchor baseline**, and
2. **direction-specific information value beyond the random-sign control**.

For the primary PL2B strength `lambda = 0.25`, define for every non-tie case:

`g_correct = baseline_abs_error - correct_abs_error`

`g_wrong = baseline_abs_error - wrong_abs_error`

`g_random = 0.5 * (g_correct + g_wrong)`

and the primary direction-specific response

`information_advantage = g_correct - g_random`

which is identically

`0.5 * (g_correct - g_wrong)`.

Positive `information_advantage` means that knowing the correct direction adds
value beyond a 50/50 random sign prior.  This is distinct from `g_correct > 0`,
which asks whether the correct weak prior improves over the strong-anchor baseline.

Use the already-frozen PL2B gain tolerance `1e-12`.  A case is
`DIRECTIONALLY_USEFUL` only when both `g_correct > 1e-12` and
`information_advantage > 1e-12`.  Truth-tie cases remain accounted for but have
no directional response and must not be assigned a correct/wrong direction.

## Align PL2B outcomes to the PL1 predictor unit

PL1 has one predictor row per `(evaluation_id, reaction_id)`.  PL2B can contain
more than one distinct held-out truth pair for a given evaluation.  PL2C must not
repeat the PL1 predictor row once per truth pair as though those were independent
predictor observations.

For each DEVELOPMENT `(evaluation_id, reaction_id)`:

- retain each distinct PL2B truth pair exactly once;
- exclude `TRUTH_TIE` pairs from directional-response denominators but report their
  count;
- aggregate non-tie pair responses with equal weight across distinct truth pairs;
- record `n_distinct_pairs`, `n_non_tie_pairs`, and `n_truth_tie_pairs`;
- record `non_tie_pair_fraction` so directional-response coverage is visible rather than silently conditioned away;
- compute the mean and median `g_correct`;
- compute the mean and median `information_advantage`;
- compute the fraction with positive correct gain;
- compute the fraction with positive information advantage;
- compute the fraction classified `DIRECTIONALLY_USEFUL`.

Rows with zero non-tie truth pairs are retained with an explicit
`NO_DIRECTIONAL_TRUTH` status and are excluded from directional association
statistics.  Alias multiplicity must never weight these aggregates.

The primary PL2C response is the scale-free
`directionally_useful_fraction`.  `non_tie_pair_fraction` is a required coverage
diagnostic and must be reported alongside utility so the high PL2B truth-tie rate
cannot be hidden by conditioning on non-ties.  The primary continuous secondary response is
`mean_information_advantage`.  `mean_correct_gain` is a secondary baseline-utility
response.  Do not normalize by truth magnitude or fit a denominator threshold in
PL2C v1.

## Frozen PL1 feature panel

PL2C v1 uses the following pre-specified PL1 features only.  They were frozen by
PL1 before PL2 outcomes and are selected for mechanistic relevance to sign-only
information, not by PL2B association strength:

- `sign_magnitude_eta2`
- `sign_magnitude_correlation`
- `directional_entropy3`
- `dominant_sign_mass`
- `n_supported_sign_states`
- `anchor_width80_contraction`
- `joint_ess`
- `abs_delta_anchor_target_correlation`

Do not add/drop features after inspecting PL2C outcomes.  PL1 status fields and
pathway flags may be used for missingness/degeneracy and sensitivity summaries,
but not as a route to outcome-driven feature selection.

Non-finite PL1 values remain non-finite; do not silently impute them to zero.
Every association reports its finite denominator and the corresponding PL1
status/degeneracy counts.

## Development analyses permitted in PL2C

PL2C v1 is deliberately descriptive and rule-development oriented.  It may:

- compute deterministic Spearman associations between each frozen feature and
  each frozen response on DEVELOPMENT reactions;
- compute the same univariate associations within each algorithm x RNA context
  and summarize their direction/stability;
- construct outcome-blind quartile bins from all DEVELOPMENT PL1 predictor values
  using NumPy linear quantiles and report response summaries by bin; for the discrete
  `n_supported_sign_states` feature, use its exact integer levels instead of forcing
  artificial quartiles;
- construct a fixed 4 x 4 DEVELOPMENT map of `sign_magnitude_eta2` quartile by
  `directional_entropy3` quartile, with explicit non-finite/degenerate accounting;
- repeat descriptive summaries after excluding proximal/same-pathway reactions
  only as a labeled sensitivity analysis.

PL2C v1 must not yet publish a reaction ranking, optimize a threshold, tune a
lambda, select a feature subset, fit an unrestricted machine-learning model, or
inspect the confirmation-holdout outcomes.  If the development landscape shows a
coherent relationship, a subsequent rule-freeze step can define the exact PL2D
confirmation hypothesis before the holdout is opened.

No inferential biological p-values may treat reaction rows, CT2A x GL261 mouse
crosses, or repeated evaluation rows as independent biological replicates.  These
are computational benchmark associations.

## Required PL2C artifacts

Publish under
`outputs/dmi_bridge_pl2c_geometry_utility_development_v1/`:

- `BRIDGEPL2C_SOURCE_AUDIT.json`
- `BRIDGEPL2C_ANALYSIS_CONTRACT.json`
- `BRIDGEPL2C_REACTION_SPLIT.tsv`
- `BRIDGEPL2C_DEVELOPMENT_EVALUATION_REACTION.tsv.xz`
- `BRIDGEPL2C_FEATURE_ASSOCIATIONS.tsv`
- `BRIDGEPL2C_CONTEXT_ASSOCIATIONS.tsv`
- `BRIDGEPL2C_FEATURE_QUARTILES.tsv`
- `BRIDGEPL2C_ETA2_ENTROPY_MAP.tsv`
- `BRIDGEPL2C_QC.json`
- `BRIDGEPL2C_MANIFEST.json`

The manifest must state explicitly that confirmation outcomes were not analyzed
and that PL2D has not started.

## Terminal gate

PL2C may end as

`PL2C_DEVELOPMENT_GEOMETRY_UTILITY_COMPLETE`

only when source hashes, split identity, development-only filtering, PL1 join
coverage, distinct-pair aggregation, response identities, deterministic summaries,
and no-forbidden-work gates all pass.  The terminal manifest must also record

`pl2d_confirmation_status = UNTOUCHED`.

Stop after PL2C.  Do not automatically open the 1,045 confirmation reactions.
