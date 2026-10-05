# DMI-BRIDGE-PL2 — truth-referenced network-wide sign-only utility

PL2 follows the frozen, outcome-blind PL1 landscape.  It asks the operational
question that PL1 alone cannot answer: **if the only new information about a
reaction B is the direction of its CT2A–GL261 difference, does that direction
improve prediction of the magnitude |ΔB| beyond the qualified strong HEX1
anchor, and under what pre-existing PL1 geometry and direction reliability?**

PL2 does not change PL1, rank reactions from PL1, or use PL2 outcomes to choose
reaction B.  All 4,180 non-HEX1 cache reactions remain in the network-wide
registry.  Same-subsystem and the pre-existing proximal pathway set remain in
the primary analysis; their exclusion is sensitivity-only.

## Separation from historical BRIDGE experiments

The historical BRIDGE-2A posterior used a magnitude-dependent bounded tanh
operator.  It remains a separate mechanism experiment and its numeric outcomes
are not PL2 inputs.  Historical BRIDGE-3C/3D used a sign-only operator on a
small L/G supported subset.  PL2 may reuse the already-qualified mathematical
sign-score/budget semantics as provenance, but it is a new network-wide
truth-referenced experiment and does not use historical C2/C3 outcome values to
select reactions, thresholds, lambda, or regimes.

## Frozen weak operator

For a candidate CT2A/GL261 pair and reaction B,

`delta_B = B_CT2A - B_GL261`.

For supplied direction `d in {-1,+1}`,

- `s = +1` if `d * delta_B > 1e-12`,
- `s = -1` if `d * delta_B < -1e-12`,
- `s = 0` for a numerical tie.

The posterior is proportional to

`p_lambda(pair) = p0(pair) * exp(lambda_total * s)`.

The primary weak budget is `lambda_total = 0.25`, fixed before PL2 outcome
access.  `0.5` and `1.0` are sensitivity strengths; `0` is the exact-recovery
QC baseline.  Lambda must not be chosen from PL2 outcomes.

## Primary endpoint: magnitude of the condition difference

PL2 explicitly targets the user's scientific question about the **magnitude of
the difference**, not exact signed flux recovery.  For a held-out truth pair,

`m_truth = |delta_B_truth|`

and the inference estimate is

`m_hat = E_p[|delta_B|]`.

Primary absolute error is `|m_hat - m_truth|`; positive error gain means the
sign-only posterior reduces this error relative to the strong-HEX1 baseline.
Signed-delta estimates may be reported secondarily, but they must not replace
the primary magnitude endpoint.

## Outcome-blind held-out truth selection

Truth candidate identities must be selected without using reaction-B values or
PL2 gains.  Within each exact condition × mouse × reconstruction-method × RNA
context, sort candidates by the strong-anchor `HEX1` coordinate (with exact
identity fields as deterministic tie breakers), use the frozen strong-anchor
weights, and select weighted quantiles `0.10`, `0.50`, and `0.90`.  Collapse
duplicate candidate selections while retaining quantile aliases.  Across CT2A
and GL261, pair only selections sharing a quantile alias.

For each truth pair, remove the selected candidate from its own condition's
inference pool and renormalize the remaining strong-anchor weights before
forming the Cartesian product.  A reaction-B truth tie is retained with an
explicit `TRUTH_TIE` status and is not assigned a fabricated correct/wrong
binary direction.

This selection is intentionally based on A=HEX1 rather than B so that the same
truth identities can benchmark the entire network without selecting easy or
large-effect B reactions.

### Canonical PL2B truth-pair unit and quantile aliases

The frozen PL2A registry contains 1,200 quantile aliases, but alias collapse
means these correspond to 885 distinct selected truth pairs under the canonical
identity `(evaluation_id, ct2a_candidate_id, gl261_candidate_id)`.  Of these,
836 distinct truth pairs are evaluable and 49 are non-evaluable under the
frozen normalized-weight representation.

PL2B must compute each distinct evaluable truth pair exactly once.  All 1,200
quantile aliases remain preserved in an alias map for audit and quantile-specific
sensitivity summaries.  Primary overall summaries give each distinct truth pair
equal weight; they must not count the same selected truth pair two or three
times merely because it carries multiple quantile aliases.  Alias-weighted
overall summaries are sensitivity-only and must be labeled as such.

## Correct, wrong, randomized, and reliability-varying direction

For every non-tie truth case, compute both the correct truth direction and its
exact opposite under the same lambda.  These two arms are sufficient to derive
a reliability curve without stochastic Monte Carlo.

If a directional source is correct with probability `q`, its expected gain is

`g(q) = q * g_correct + (1-q) * g_wrong`.

The frozen reliability grid is `q = {0, 0.25, 0.5, 0.75, 1}`.  Thus `q=1` is
the correct arm, `q=0` is the wrong arm, and `q=0.5` is the exact expected
random-sign null.  A per-case break-even reliability may be reported when the
linear curve has a unique crossing in `[0,1]`.

PL2B stores the correct and wrong endpoint gains and derives the reliability
grid analytically.  It must not materialize five reliability rows per
case/lambda.  The frozen gain classification tolerance is `1e-12`: gains above
that value are `IMPROVED`, gains below `-1e-12` are `HARMED`, and values within
the tolerance are `TIED`.

## PL1 predictors and the "when does B become predictable?" question

PL1 remains immutable and supplies outcome-blind predictors only.  The primary
prospective predictor is `sign_magnitude_eta2`, because it directly measures
the fraction of baseline |ΔB| variance explained by the sign category when
finite.  Secondary frozen predictors include directional entropy / dominant
sign mass, sign-specific width contraction, anchor-induced concentration,
residual HEX1–B coupling, and joint/sign-specific ESS.

Do not fit a new bridgeability score or threshold from PL2 outcomes.  Report
PL2 utility stratified by exact reconstruction/RNA context and relate utility
to the frozen continuous PL1 predictors.  Degenerate PL1 values remain named
categories.  Any PL1-only descriptive bin boundaries must be written and
hashed before PL2 outcome computation.

## Dependence and interpretation

The 25 CT2A×GL261 mouse crosses inside a reconstruction/RNA context are
repeated/dependent computational benchmark units, not 25 biological
replicates.  Reactions are also network coordinates, not independent subjects.
Primary summaries therefore preserve method/RNA context and report descriptive
fractions improved/tied/harmed and error-gain distributions rather than treating
millions of reaction rows as independent observations.

PL2 is solver-free.  It must not call optimization, FVA, new sampling,
reconstruction, or new strong-anchor weight fitting.  Gene-score concordance
is not part of the primary PL2 test.

## Required staging

Before numerical PL2 outcome generation, freeze a PL2A preparation bundle that
contains: source/hash audit; the complete 4,180-reaction registry; held-out
HEX1-based truth selections and truth-pair registry; exact PL1 predictor join
contract; primary/sensitivity lambda values; reliability grid; pathway
sensitivity definitions; and all expected row-count/accounting identities.
PL2B production may start only after this preparation bundle is validated.

## PL2A holdout support terminology

PL2A preserves all 1,200 truth-pair aliases.  Under the frozen floating-point
normalized-weight representation, 1,110 aliases are evaluable and 90 are
reported as:

`NON_EVALUABLE_ZERO_POST_HOLDOUT_MASS under the frozen normalized-weight representation`

This status does not assert that raw-weight support is absent.  No truth
candidate is replaced, no pseudomass is introduced, no alternate weight arm is
substituted, and no reaction-B information is used to decide evaluability.
