# DMI-BRIDGE-PL2B — compact sign-only utility production contract

PL2B is the first stage permitted to access reaction-B truth values.  It consumes
the qualified PL2A preparation bundle and computes the truth-referenced sign-only
utility experiment.  It must not change PL1 or PL2A scientific identities.

## Frozen admission gate

PL2B is admitted only when all of the following hold:

- `outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_MANIFEST.json` has
  SHA256 `8ae3eb60eb23614e379e38e3aa3cb72ec16c6bdb69aa9c74cc94499e31f1c312`.
- PL2A status is `PL2A_READY_FOR_PL2B_SIGN_ONLY_UTILITY`.
- PL2A counts remain 4,180 reactions, 480 selection rows, 1,200 aliases,
  1,110 evaluable aliases, 90 non-evaluable aliases, and 18 zero-post-holdout
  selection aliases.
- Every PL2A artifact listed in the PL2A manifest verifies against its recorded
  SHA256 before any B value is read.
- The candidate table, strong-anchor weights, flux cache, and PL1 predictor
  artifacts verify against the source identities already frozen by PL2A.

PL2B must fail closed on any mismatch.

## Canonical truth-pair identity

Define a distinct held-out truth pair by the exact tuple:

`(evaluation_id, ct2a_candidate_id, gl261_candidate_id)`.

The frozen PL2A alias registry must deterministically collapse to:

- 885 distinct truth pairs total;
- 836 distinct evaluable truth pairs;
- 49 distinct non-evaluable truth pairs;
- evaluable alias multiplicities: 642 pairs with one alias, 114 with two, and
  80 with three;
- non-evaluable alias multiplicities: 24 pairs with one alias, 9 with two, and
  16 with three.

All 1,200 aliases remain preserved in a PL2B alias map.  Scientific outcomes are
computed once per distinct evaluable truth pair, giving exactly
`836 * 4180 = 3,494,480` reaction/truth-pair case rows before truth-tie handling.
No replacement truth, pseudomass, alternate arm, or alias duplication is allowed.

## Held-out inference distribution

For each distinct evaluable pair, use the exact frozen candidate identities and
strong-anchor weights.  Remove the selected CT2A candidate from the CT2A pool
and the selected GL261 candidate from the GL261 pool, then renormalize each
condition separately.  The recomputed post-holdout masses and ESS values must
match PL2A within deterministic floating-point tolerance before any reaction-B
outcome is accepted.

For each reaction B, the truth is the selected-candidate contrast

`delta_B_truth = B_CT2A_truth - B_GL261_truth`.

The inference support is the Cartesian product of the two held-out condition
pools.  The strong-HEX1 baseline uses only the renormalized frozen strong-anchor
weights.

## Outcome contract

Use `TIE_TOL = 1e-12`, `GAIN_TOL = 1e-12`, lambda grid
`{0, 0.25, 0.5, 1.0}`, and reliability grid
`{0, 0.25, 0.5, 0.75, 1}`.  Lambda `0.25` is primary; `0.5` and `1.0` are
sensitivity strengths; lambda `0` is QC only.

For every non-tie truth case, compute the correct and exactly wrong binary
direction arms.  The primary estimate is `E[|delta_B|]`; primary error is its
absolute error to `|delta_B_truth|`; gain is baseline absolute error minus
posterior absolute error.  Classify gain as `IMPROVED`, `TIED`, or `HARMED`
using the frozen `1e-12` tolerance.

For truth ties, retain the case as `TRUTH_TIE` and do not fabricate correct or
wrong directions.

Lambda-zero must recover the baseline for both directions to tolerance and is
reported only in QC.  Do not materialize lambda-zero scientific rows.

Reliability is derived exactly from the endpoint gains:

`g(q) = q*g_correct + (1-q)*g_wrong`.

Do not materialize the five-q reliability grid at case level.  Store endpoint
gains plus break-even reliability; derive q-specific summaries analytically.

## Compact case-level storage

Publish one deterministic partitioned case table with exactly 3,494,480 rows,
sorted by `(truth_pair_id, reaction_id)`.  Use enough fixed parts to keep every
compressed artifact below the repository file-size ceiling; the part inventory,
row counts, sizes, SHA256 values, schema, and a logical-content SHA256 must be
manifested and validated.

Each row must contain, at minimum:

- truth-pair/evaluation identity, algorithm, RNA context, mouse pair, alias
  multiplicity and alias list, reaction identity and pathway flags;
- truth delta, truth magnitude, truth direction/status;
- baseline magnitude/signed estimate, absolute error, joint ESS, and
  negative/tie/positive masses;
- for lambda 0.25, 0.5, and 1.0: correct/wrong magnitude estimates, absolute
  errors, absolute-error gains, gain statuses, posterior ESS values, exact
  random-sign expected gain, and break-even reliability.

Do not duplicate frozen PL1 predictors into every case row.  Preserve their
cryptographic source identity and join them by `(evaluation_id, reaction_id)`
for downstream PL2C geometry-to-utility analysis.

## PL2B summaries

PL2B may publish fixed descriptive utility summaries only.  Primary overall
summaries weight each distinct truth pair once.  Preserve method x RNA context
and report non-tie/tie counts, gain quantiles, and improved/tied/harmed fractions
for the primary lambda across the frozen reliability grid.  Lambda sensitivities
are secondary.  Quantile-alias-weighted summaries are sensitivity-only.

Do not fit a bridgeability score, choose predictor thresholds, optimize lambda,
rank/select reactions from PL2 outcomes, or assign inferential p-values treating
mouse crosses/reactions as independent replicates.  The PL1 predictor-to-utility
relationship remains a separate downstream analysis over the frozen PL2B output.

## Required artifacts

Publish under `outputs/dmi_bridge_pl2b_sign_only_utility_v1/`:

- `BRIDGEPL2B_SOURCE_AUDIT.json`
- `BRIDGEPL2B_ANALYSIS_CONTRACT.json`
- `BRIDGEPL2B_DISTINCT_TRUTH_PAIR_REGISTRY.tsv`
- `BRIDGEPL2B_ALIAS_MAP.tsv`
- partitioned `BRIDGEPL2B_CASE_OUTCOMES` data plus its part manifest
- `BRIDGEPL2B_CONTEXT_SUMMARY.tsv`
- `BRIDGEPL2B_QC.json`
- `BRIDGEPL2B_MANIFEST.json`

The terminal status is `PL2B_SIGN_ONLY_UTILITY_COMPLETE` only if all admission,
identity, row-count, deterministic-storage, lambda-zero, and no-forbidden-work
gates pass.

## Forbidden work

PL2B must not invoke solver/optimization, FVA, sampling, reconstruction, new
strong-anchor fitting, gene-score concordance, historical tanh outcomes,
historical BRIDGE-3 outcome values, PL2 outcome-driven reaction selection,
outcome-driven threshold fitting, or new PL1 predictor construction.
