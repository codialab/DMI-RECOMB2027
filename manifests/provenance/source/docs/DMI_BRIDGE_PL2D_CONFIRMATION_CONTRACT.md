# DMI-BRIDGE-PL2D — one-shot confirmation contract

## Purpose

PL2D is the one-shot computational confirmation of the geometry-to-utility hypothesis frozen in PL2D-H and amended, before holdout opening, by the PL2D-H support-gate amendment. PL2D opens only the 1,045 reactions already assigned to `CONFIRMATION_HOLDOUT`. It performs no new metabolic simulation, fitting, feature selection, threshold tuning, or reaction selection.

The primary confirmation unit is the held-out reaction. Evaluation rows and distinct truth-pair cases are intermediate computational units only.

## Immutable prerequisite chain

Before reading any confirmation predictor or outcome value, the confirmation adapter must verify all of the following exact identities:

- PL2D-H manifest SHA256: `a1dd76d6b40e741e07efe64c34bc39eaad289a46bf6c870917a4f59cab688427`;
- PL2D-H confirmation registry SHA256: `02a38afe3cb25b8c85a7774c1f20118ad9601ba74f0df0fbff021ba2494bd2e5`;
- PL2D-H amendment manifest SHA256: `fb7df9e208da2e6c9841d43a3b572038e59943a4221535705eef9f5dc0ed6a3a`;
- amended confirmation contract SHA256: `08baccb983fa7690d1e10dd965b95d42f0223a798e797dbf0b0c1a830491cfb4`;
- PL2C manifest SHA256: `6abcf149cd6adb294c0931585cec28f992b453451a5b5f03e10b3ffb2f9b9c39`;
- PL2B manifest SHA256: `7d9ca97a3f4ea5aa8b0b983e09ee83fc3c82f326d1b1cbf9c6137cfadf25539d`;
- PL2B case logical-content SHA256: `c3670992a0ad33698af6c7620429bed48e0ecaa98adb3c87e648eefb1ba70363`;
- PL2A manifest SHA256: `8ae3eb60eb23614e379e38e3aa3cb72ec16c6bdb69aa9c74cc94499e31f1c312`;
- PL1 feature logical-content SHA256: `46aa9f0410f89d03b020e01fefa517348faae1d55d9bfe8387ac04b39bb9ac93`.

The original PL2D-H bundle and amendment bundle are immutable inputs. PL2D must not rewrite either one.

The effective confirmation contract is the original PL2D-H contract plus exactly the append-only support amendment `minimum_finite_reactions_per_context: 100 -> 30`. Any other semantic contract difference is fail-closed.

## Frozen confirmation population

The confirmation registry contains exactly 1,045 reactions at frozen split ranks 0 through 1044. No development reaction may contribute to a PL2D scientific output.

Expected raw accounting is fixed before holdout access:

- PL1 confirmation predictor rows: `370 × 1,045 = 386,650` evaluation × reaction rows;
- PL2B confirmation case rows: `836 × 1,045 = 873,620` distinct-truth-pair × reaction cases;
- PL2B development case rows skipped before outcome parsing: `836 × 3,135 = 2,620,860`;
- reaction-level confirmation rows: exactly 1,045.

When PL2B case parts are streamed, `reaction_id` must be inspected first. Rows not in the confirmation registry must be skipped before truth status or gain fields are converted. This makes the PL2D dataset confirmation-only even though the upstream PL2B parts contain both splits.

## Frozen response reconstruction

PL2D must reproduce the exact PL2C response definition at primary lambda `0.25`.

For every confirmation `(evaluation_id, reaction_id)`, include every evaluable distinct `truth_pair_id` exactly once. Quantile aliases must never be used as weights.

For a non-tie truth-pair case:

- `g_correct` is `lambda_0_25_correct_absolute_error_gain`;
- `g_wrong` is `lambda_0_25_wrong_absolute_error_gain`;
- `g_random` is `lambda_0_25_random_sign_expected_gain`;
- require `g_random = 0.5 * (g_correct + g_wrong)` within the frozen deterministic tolerance;
- `information_advantage = g_correct - g_random`;
- `directionally_useful = (g_correct > 1e-12) and (information_advantage > 1e-12)`.

Truth ties remain truth ties and do not receive a fabricated direction.

Aggregate distinct truth pairs with equal weight to the evaluation × reaction unit, producing the same PL2C fields including `directionally_useful_fraction`, `mean_information_advantage`, and `non_tie_pair_fraction`. Evaluation × reaction rows with no non-tie truth are retained with undefined directional responses.

## Frozen predictor reconstruction

Join confirmation response aggregates to the frozen PL1 feature table on exactly `(evaluation_id, reaction_id)`.

The primary predictor is `sign_magnitude_eta2`. The supportive predictors are:

- `directional_entropy3`;
- `dominant_sign_mass`;
- `n_supported_sign_states`.

Carry the frozen pathway annotations `proximal_set_member` and `same_subsystem_as_HEX1` for supportive sensitivity analyses.

Predictor aggregation must not condition on response availability.

## Authoritative reaction-level numerical reference

Use exactly the numerical reference frozen in PL2D-H:

1. `pandas.to_numeric(errors="coerce")` for numeric columns;
2. convert `+/-inf` to missing;
3. reaction means from `DataFrame.groupby("reaction_id", sort=True).mean()` with missing values skipped;
4. Spearman rho from `scipy.stats.spearmanr` on finite paired vectors;
5. one-control partial Spearman from the three corresponding SciPy Spearman coefficients using the frozen closed-form formula.

The adapter must record Python, pandas, NumPy, and SciPy versions.

For each held-out reaction define:

- primary predictor: mean finite `sign_magnitude_eta2` across all evaluation rows;
- primary response: mean finite `directionally_useful_fraction` across evaluation rows with defined directional response;
- coverage diagnostic: mean finite `non_tie_pair_fraction`;
- secondary magnitude response: mean finite `mean_information_advantage`;
- supportive predictors: corresponding reaction-level finite means.

## One-shot primary confirmation rule

The effective amended primary support requirements are:

- at least 700 reactions with both finite reaction-level eta2 and defined reaction-level directional usefulness;
- exactly 16 algorithm × RNA context results;
- every context has at least 30 finite reaction pairs;
- every context-specific Spearman rho is finite (this is the frozen meaning of all 16 contexts being evaluable);
- overall primary Spearman rho is finite.

The primary effect requirements are:

1. overall reaction-level Spearman rho between eta2 and directional usefulness is at least `+0.50`;
2. at least 12 of 16 algorithm × RNA context associations are strictly positive.

Classify exactly:

- `PL2D_INSUFFICIENT_SUPPORT` if any primary support requirement fails;
- `PL2D_NOT_CONFIRMED` if support is adequate but either effect requirement fails;
- `PL2D_CONFIRMED` only if every support and effect requirement passes.

The primary status must be computed and frozen before supportive results are interpreted. Supportive analyses cannot rescue or downgrade the primary status.

## Supportive stability bootstrap

The reaction bootstrap is supportive only and is not a biological confidence interval or a confirmation gate.

Freeze:

- population: all 1,045 confirmation reaction rows;
- resampling unit: reaction row / reaction ID;
- replicates: 5,000;
- random generator: NumPy `Generator(PCG64(seed))`;
- seed: `20260929`;
- each replicate samples 1,045 reaction indices with replacement, filters to finite predictor-response pairs inside that replicate, and recomputes SciPy Spearman rho;
- interval: 2.5th and 97.5th percentiles of finite replicate rho values using NumPy `quantile(..., method="linear")`.

Record the number of finite and non-finite replicates. The bootstrap interval is descriptive computational stability only.

## Prespecified supportive analyses

Only after the primary status is fixed, report:

- eta2 versus reaction-level mean information advantage, expected positive;
- partial Spearman eta2 versus usefulness controlling mean non-tie coverage, expected positive;
- directional entropy versus usefulness, expected positive;
- dominant sign mass versus usefulness, expected negative;
- supported sign-state count versus usefulness, expected positive;
- all 16 context-specific eta2/usefulness rhos and finite support counts;
- eta2/usefulness after excluding proximal reactions;
- eta2/usefulness after excluding same-subsystem-as-HEX1 reactions;
- eta2/usefulness after excluding either flag.

No supportive result changes the primary status.

Do not calculate a new composite score, fit an unrestricted model, search thresholds, create a favorable subgroup, rank reactions, or inspect alternate lambdas.

## Required publication

Publish atomically under `outputs/dmi_bridge_pl2d_confirmation_v1/`:

- `BRIDGEPL2D_SOURCE_AUDIT.json`;
- `BRIDGEPL2D_CONFIRMATION_EVALUATION_REACTION.tsv.xz` — exactly 386,650 rows;
- `BRIDGEPL2D_CONFIRMATION_REACTION.tsv` — exactly 1,045 rows;
- `BRIDGEPL2D_CONTEXT_RESULTS.tsv` — exactly 16 context rows;
- `BRIDGEPL2D_PRIMARY_RESULT.json`;
- `BRIDGEPL2D_SUPPORTIVE_RESULTS.json`;
- `BRIDGEPL2D_BOOTSTRAP.tsv.xz` — exactly 5,000 replicate rows;
- `BRIDGEPL2D_QC.json`;
- `BRIDGEPL2D_MANIFEST.json`.

The manifest status must equal the resulting primary status: `PL2D_CONFIRMED`, `PL2D_NOT_CONFIRMED`, or `PL2D_INSUFFICIENT_SUPPORT`.

The manifest must explicitly record that PL2D opened the previously reserved confirmation set under the frozen original contract plus the qualified support amendment, and that no post-opening threshold, feature, lambda, population, or status-rule changes occurred.

A rerun against an existing valid bundle must return a validated no-op. Any source or output mutation must fail closed.

## Forbidden work

PL2D must not run solver/optimization/FVA, perform new sampling or reconstruction, refit anchor weights, modify PL1/PL2A/PL2B/PL2C/PL2D-H artifacts, tune lambda, select features, select reactions, alter the 700/30/0.50/12-of-16 gates, create a new bridgeability score, or characterize the computational reaction holdout as independent biological validation.

PL2D is the terminal confirmatory analysis for this frozen DMI-BRIDGE geometry hypothesis. Any later exploratory analysis must be clearly separated from this one-shot confirmation result.
