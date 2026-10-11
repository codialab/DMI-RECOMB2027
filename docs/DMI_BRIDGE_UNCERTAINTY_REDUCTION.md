# DMI-BRIDGE weak-cue uncertainty analysis

This additive analysis asks when the frozen weak directional cue changes the
weighted uncertainty of held-out flux magnitude, and whether that change is
associated with accuracy. It does not modify frozen analysis code, input
identities, manuscript files, or published outputs.

## Reproduction

Run from the repository root:

```bash
python scripts/dmi_bridge_uncertainty_reduction_v1.py --preflight
python scripts/dmi_bridge_uncertainty_reduction_v1.py --benchmark
python scripts/dmi_bridge_uncertainty_reduction_v1.py --run
```

`--preflight` checks the pinned SHA-256 for every required input. The empirical
runner refuses to start if any asset is missing, mismatched, a symlink, or
resolves outside this repository. When an authorized frozen bundle is
available, place each file under
`.cache/recomb2027-dmi/frozen/<repository-relative-path>`; the exact relative
paths and hashes are enumerated in the runner's `EXPECTED` table and the
repository `FROZEN_LINEAGE.json`. Do not regenerate candidates, flux vectors,
weights, anchors, or truth pairs to satisfy this gate.

The synthetic benchmark creates no biological results. On the implementation
host, the latest 64-case synthetic run with 361 support points took 0.0601
seconds in the fastest of two runs (1,065 cases/second; peak process RSS
37,784 KiB). A simple linear extrapolation for 2 × 836 × 4,179 cases is 6,563
seconds (1.82 hours),
before file I/O, compression, aggregation, and plotting. Rerun the benchmark
on the execution host before authorizing a full empirical run.

## Frozen definitions and estimands

For each reaction B and held-out mouse pair, the support is the existing
19-by-19 Cartesian product after removing the held-out truth sample in each
cohort and renormalizing the frozen A1 BIO-0 strong-anchor weights or A2 A20
dual-anchor weights. The response is `M = |v_CT2A − v_GL261|`. Truth ties use
the frozen `1e-12` threshold and are excluded from cue comparisons. The same
836 distinct truth pairs and 4,179 non-anchor reactions are required for A1
and A2.

The primary statistic is weighted population variance
`V = Σᵢ wᵢ(Mᵢ − μ)²`. Absolute reduction is `V₀ − Vλ`; relative reduction is
`1 − Vλ/V₀`. Relative reduction is omitted when
`V₀ ≤ 100 ε max(E[M²], E[M]²)` in float64; the absolute change remains
available. A change is called negligible when its magnitude is at most
`100 ε` times the largest baseline/updated second-moment or squared-mean scale.
Weighted SD, inverse-CDF weighted IQR, central 90% inverse-CDF width, and
weighted discrete CRPS are also calculated. The inverse-CDF quantile is the
first support value whose cumulative normalized weight reaches q, matching
the existing PL1 convention.

For sign state s in {−1, 0, +1}, write baseline mass, conditional magnitude
mean, and conditional variance as πₛ, μₛ, and σ²ₛ. A cue with direction d and
strength λ multiplies mass by `cₛ = exp(λ d s)`, giving
`π'ₛ = cₛπₛ / Σₜ cₜπₜ`. The resulting variance is

```text
Vλ = Σₛ π'ₛ σ²ₛ + Σₛ π'ₛ (μₛ − μλ)²,
μλ = Σₛ π'ₛ μₛ.
```

The cue preserves each within-sign distribution. It changes the within-group
contribution only through mixture proportions in the outer sum, and changes
the between-group contribution through both proportions and the new overall
mean. Positive exponential weights leave supported values unchanged, so a
variance or quantile reduction is a concentration change, not contraction of
the feasible range. The numerical core verifies the baseline CRPS mixture
identity and unchanged weighted support.

Accuracy is the existing PL2 expected-absolute-magnitude prediction gain.
CRPS improvement is separately baseline CRPS minus cue CRPS. A shifted mean,
lower directional entropy, or lower prediction error is not itself counted as
magnitude uncertainty reduction.

Uncertainty statistics are computed at λ=0, 0.25, 0.5, and 1.0. Existing
prediction-gain estimates at those strengths are retained independently, and
the λ=0 no-op is checked against the baseline distribution.

CRPS is computed from
`E|X−y| − 0.5 E|X−X′|`. Within each held-out case, conditional sign-group
truth distances and all nine cross-group absolute-distance expectations are
calculated once. Their three-state mixture is reused for both cue directions
and all four λ values; this avoids repeated pairwise support calculations.

## Aggregation, prediction, and inference

Case-level counts retain repeated reaction/evaluation/truth-pair rows and
must not be interpreted as independent biological observations. The
reaction–evaluation table averages distinct truth-pair cases equally within
each unit. Reaction-level outcomes then average eligible evaluation units
equally. A1/A2 analyses and all prospective prediction targets use matched
populations.

Prospective geometry models predict correct-cue and incorrect-cue relative
variance reduction separately. Both use identical finite reactions and the
existing frozen Development/Confirmation split. Predictor sets are fixed:
η²_dir alone; η²_dir with H_dir; η²_dir with dominant-direction mass; η²_dir
with baseline non-tie coverage (`p_negative + p_positive`); η²_dir with
supported sign-state count; and all five baseline geometry descriptors.
Only Development outcomes train standardized Ridge models. Confirmation
truth-derived gains and usefulness are evaluation outcomes only, never
predictors. Confirmation MAE/RMSE/R² are compared with the Development-mean
baseline and η²_dir-only model.

Observed usefulness associations are explanatory only. The output reports
raw rank association and η²_dir association adjusted for observed usefulness,
for uncertainty reduction and prediction gain separately. They are not
prospective estimates.

The report assigns title alternative A only when correct cues reduce absolute
variance in a majority of eligible cases under both anchors and η²_dir alone
beats the Development-mean Confirmation MAE baseline for the correct-cue
target under both anchors. It assigns B when reductions occur but are not a
majority under both anchors, while the full prespecified baseline geometry
set beats that Confirmation baseline for the correct-cue target under both
anchors. Otherwise it assigns C. The rule is intentionally conservative and
is reported alongside the component counts and scores; prediction-error gains
do not count as evidence for variance reduction.

The planned uncertainty intervals resample the five CT2A and five GL261 mouse
identities independently with replacement, cross their multiplicities, and
reuse each draw for A1 and A2. This preserves repeated evaluations and
reactions within sampled mouse pairs, while not representing new cohorts,
methods, or candidate-generation uncertainty. The table also includes
leave-one-mouse-out estimates for each cohort. All intervals must be read with
the five-mouse-per-cohort limitation in mind.

## Outputs after a successful empirical run

The runner writes to `outputs/dmi_bridge_uncertainty_reduction_v1/` and resumes
only from input-hash-matched, checkpointed case partitions. Outputs include partitioned case metrics,
reaction–evaluation summaries, reaction-level plotting tables, the
illustrative support/weight table, prospective model scores, descriptive and
explanatory association tables, case outcome counts, mouse-cluster bootstrap
intervals, publication-style PNG/SVG figures, a scientific report, QC, and an
input/output hash manifest. Accuracy gain, uncertainty reduction, and CRPS
remain separate columns and report sections.

## Current data gate

At implementation time, 165 of 331 pinned required artifacts were available
in the repository, including the complete published A21/A22 bundles and PL2C
split. The candidate table, flux cache, A1/A2 weights, and PL1 and PL2A/PL2B
source bundles and partitioned outputs accounted for the 166 missing files.
No empirical cases were
run, no substitute inputs were used, and no simulated results are reported as
findings. The report must not assign a title verdict until the frozen empirical
analysis completes and the prospective Confirmation results are available.

## Focused validation

```bash
pytest -q tests/test_dmi_bridge_uncertainty_reduction_core_v1.py
python -m py_compile scripts/dmi_bridge_uncertainty_reduction_core_v1.py scripts/dmi_bridge_uncertainty_reduction_v1.py
```

The tests cover tie handling, λ=0 identity, weighted quantiles, scale-aware
relative-variance exclusions, variance mixture calculations, CRPS parity,
support invariance, exact matched populations, split-safe separate prediction
targets, and exclusion of observed usefulness from prospective predictors.
