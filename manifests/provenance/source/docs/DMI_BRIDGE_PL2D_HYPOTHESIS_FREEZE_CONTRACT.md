# DMI-BRIDGE-PL2D-H — confirmation-hypothesis freeze

## Purpose

PL2D-H converts the completed PL2C development result into a one-shot confirmation contract **before any of the 1,045 confirmation-reaction outcomes are opened**. It performs no PL2D confirmation analysis and no new metabolic simulation.

The confirmation unit is the held-out **reaction**, because the PL2C split was made at reaction level. Evaluation/context rows are aggregated within reaction before the primary confirmation statistic is calculated.

## Frozen source

PL2D-H may read only the canonical PL2C development bundle. It must verify the exact PL2C manifest and all PL2C artifact hashes before use.

The PL2C manifest must have:

- status `PL2C_DEVELOPMENT_GEOMETRY_UTILITY_COMPLETE`;
- `pl2d_confirmation_status = UNTOUCHED`;
- `confirmation_outcomes_analyzed = false`;
- exactly 3,135 development reactions and 1,045 confirmation reactions.

PL2D-H may read the reaction-split registry to construct an identity-only confirmation registry. It must not read PL1 predictor values or PL2B outcomes for confirmation reactions.

## Why the primary feature is eta2

PL2C showed several strong geometry associations with `directionally_useful_fraction`. For confirmation, only one is primary:

`sign_magnitude_eta2`

This is the frozen sign-magnitude coupling measure and directly represents whether the sign state of reaction B organizes its magnitude under the strong HEX1 anchor. It also avoids making the numerically degenerate entropy quartile map a confirmation rule.

The other development-supported geometry features remain prespecified supportive checks:

- `directional_entropy3`: expected positive;
- `dominant_sign_mass`: expected negative;
- `n_supported_sign_states`: expected positive.

They cannot rescue a failed primary hypothesis.

## Frozen reaction-level aggregation

For each reaction, using DEVELOPMENT rows only during hypothesis freeze and CONFIRMATION_HOLDOUT rows only during PL2D confirmation:

1. Predictor: mean of all finite `sign_magnitude_eta2` values across evaluation rows for that reaction. Predictor aggregation is outcome-blind and must not be conditioned on response availability.
2. Primary response: mean of finite `directionally_useful_fraction` values across rows with a defined directional response for that reaction.
3. Coverage diagnostic: mean `non_tie_pair_fraction` across all finite evaluation rows for that reaction.
4. Secondary magnitude response: mean finite `mean_information_advantage` across rows with a defined directional response.

A reaction with no finite eta2 or no defined directional response is retained in accounting but is not a finite pair for the primary Spearman statistic.

### Authoritative numerical reference

The frozen development values below were computed from the canonical PL2C table using **pandas grouped-mean semantics**, not `math.fsum` or another algebraically equivalent summation. This distinction is binding because 1e-16-scale changes in grouped floating-point means can alter exact ties and therefore Spearman ranks.

The PL2D-H/PL2D reference implementation must therefore:

1. parse the relevant numeric columns with `pandas.to_numeric(errors="coerce")`;
2. convert `+/-inf` to missing values;
3. compute reaction means with `DataFrame.groupby("reaction_id", sort=True).mean()`, which skips missing values;
4. compute Spearman rho with `scipy.stats.spearmanr` on the finite paired vectors (SciPy average-rank handling for exact ties);
5. compute the one-control partial Spearman from the three corresponding `scipy.stats.spearmanr` coefficients using the standard closed-form formula; and
6. record the Python, pandas, NumPy, and SciPy versions used by the freeze and later confirmation adapters.

No alternative summation, tie tolerance, rank method, or correlation implementation may be substituted after seeing the development or confirmation values. A failure to reproduce the frozen development snapshot is fail-closed.

The primary statistic is Spearman rho across reactions between the reaction-level eta2 mean and reaction-level directional-usefulness mean.

## Frozen development snapshot

Using the canonical PL2C development table and the aggregation above:

- finite primary reaction pairs: 2,455;
- primary reaction-level Spearman rho: `0.8570601823118907`;
- reaction-level partial Spearman controlling `non_tie_pair_fraction`: `0.8356721975298953`;
- eta2 versus reaction-level mean information advantage: `0.3730769973633655`;
- the eta2/usefulness association is positive in all 16 algorithm × RNA contexts.

Supportive reaction-level development associations with directional usefulness are:

- `directional_entropy3`: `+0.8102865680476654` (2,716 finite reactions);
- `dominant_sign_mass`: `-0.7966029638193884` (2,716 finite reactions);
- `n_supported_sign_states`: `+0.7964386893537818` (2,716 finite reactions).

These values are discovery evidence, not confirmation results.

## One-shot PL2D primary confirmation rule

PL2D is confirmatory only for the 1,045 reactions already frozen as `CONFIRMATION_HOLDOUT`.

Primary support requirements:

- at least 700 holdout reactions with both finite reaction-level eta2 and a defined reaction-level directional-usefulness response;
- all 16 algorithm × RNA contexts must be evaluable for the context-stability audit, with at least 100 finite reaction pairs per context.

Primary effect requirements:

1. overall reaction-level Spearman rho between eta2 and directional usefulness is at least `+0.50`;
2. the association direction is positive in at least 12 of the 16 algorithm × RNA contexts.

The `+0.50` threshold is frozen before holdout opening as a substantive retained-effect criterion; it is below the development rho of 0.857 and is not to be changed after seeing confirmation data.

A deterministic reaction-bootstrap stability interval is supportive, not a biological confidence interval and not an additional success gate:

- 5,000 bootstrap replicates;
- resample confirmation reaction IDs with replacement;
- seed `20260929`;
- percentile 95% interval;
- recompute reaction-level Spearman rho in every replicate.

The bootstrap is a computational stability diagnostic only because GEM reactions are not independent biological replicates.

## Confirmation statuses

- `PL2D_INSUFFICIENT_SUPPORT`: a primary support requirement is not met.
- `PL2D_NOT_CONFIRMED`: support is adequate but either primary effect requirement fails.
- `PL2D_CONFIRMED`: all primary support and primary effect requirements pass.

Secondary/supportive analyses cannot rescue `PL2D_NOT_CONFIRMED`.

## Prespecified supportive analyses

After the primary status is fixed, PL2D may report:

- reaction-level eta2 versus mean information advantage, expected positive;
- reaction-level partial Spearman between eta2 and directional usefulness controlling non-tie coverage, expected positive;
- reaction-level directional entropy versus usefulness, expected positive;
- reaction-level dominant sign mass versus usefulness, expected negative;
- reaction-level supported-sign-state count versus usefulness, expected positive;
- pathway-proximal and same-subsystem exclusions using the same frozen definitions as PL2C;
- context-specific rho values and finite support counts.

These are supportive only. No new feature, threshold, bin boundary, composite score, model, or subgroup may be selected from confirmation outcomes.

## Information firewall

PL2D-H itself must not read or derive for confirmation reactions:

- PL1 predictor values;
- PL2B truth deltas, truth-tie status, gains, or utility responses;
- PL2C-like aggregate responses;
- any confirmation geometry-to-utility association.

It may publish only an identity-only 1,045-reaction confirmation registry copied from the already-frozen PL2C split registry. Its canonical TSV SHA256 is frozen as `02a38afe3cb25b8c85a7774c1f20118ad9601ba74f0df0fbff021ba2494bd2e5`.

The later PL2D confirmation adapter must verify this hypothesis-freeze contract hash before opening confirmation predictors/outcomes.

## Forbidden work

Neither PL2D-H nor the later one-shot confirmation may tune the eta2 threshold, change the `+0.50` rho criterion, change support minima, change the 12/16 context criterion, select lambda, fit an unrestricted ML model, create a new bridgeability score, rank/select reactions by confirmation outcome, rerun flux sampling/FVA/reconstruction/optimization, or modify PL1/PL2A/PL2B/PL2C outputs.

PL2D confirmation is a computational holdout confirmation, not independent biological validation.
