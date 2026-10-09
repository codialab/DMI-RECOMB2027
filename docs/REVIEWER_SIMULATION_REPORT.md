# Reviewer-requested sensitivity analyses

Status: Experiments 1–4 have been implemented and evaluated as post-hoc
robustness analyses. Experiments 1–2 use released compact tables; Experiments
3–4 use the authorized read-only scientific checkout. No manuscript or
supplementary material, production pipeline, frozen input, candidate vector, or
canonical result was modified.

Run with `python analyses/reviewer_checks/run_reviewer_checks.py`. The scripts,
synthetic tests, compressed tables, and exact input hashes are listed below and
under `reproduced/reviewer_checks/`.

## Inputs, definitions, and parity

The analysis used these repository-contained inputs (SHA-256 values are recorded
again in `reproduced/reviewer_checks/manifest.json`):

| Artifact | SHA-256 |
|---|---|
| `data/figure_inputs/fig2/figure2_DG_points.csv.xz` | `aedc195237389ccd131cb983d7948483fc84ab1d0df1d88b5f1c1601577cf161` |
| `data/figure_inputs/fig3/fig3_gain_reaction_evaluation.parquet.xz` | `7a83a31d29317f7b0d906fff73691428d98427d93d2fdebd4ef91294f8a46676` |
| `data/figure_inputs/fig4/fig4_geometry_paired/algorithm=CORDA/part.parquet.xz` | `0fd4d308c6d03d5ef144fb5e8c48e497e115af7bc27cdddd1010ddef4e3ba457` |
| `data/figure_inputs/fig4/fig4_geometry_paired/algorithm=GIMME/part.parquet.xz` | `f492947266440ed378540d6c0d6a40aabb8952fafef5cc728dbc8c9aa70f17ad` |
| `data/figure_inputs/fig4/fig4_geometry_paired/algorithm=RIPTiDe/part.parquet.xz` | `642fc972db1aa17854fbd631afb909fc07e9be54c7d5fdc51d23fbb60ffc2e42` |
| `data/figure_inputs/fig4/fig4_geometry_paired/algorithm=iMAT/part.parquet.xz` | `a44c629ffdc5ddcb1956d6132a8ec15d7dbb358e86b7888b512e315a73d4a070` |
| `tables/manuscript/fig2_panel_statistics.tsv` | `2488bdea4236a9897f33d21023b54e6b5257ef982c80537a1177d15b607e29ff` |

Frozen Figure 2 panels supplied the original panel-specific reaction eligibility
and usefulness response. Figure 4 supplied evaluation-level geometry and pathway
labels; eta2 means over the same 370 eligible evaluations per anchor reproduce
the Figure 2 eta2 coordinates to a maximum absolute difference below
`1.2e-16`. The A1 confirmation panel retains its original **838 finite reaction
rows**. Other Figure 2 row counts are A1 development 2,454, A2 development 2,485,
and A2 held-out 843. They were not forced to a common denominator.

**Parity passed.** Correlations recomputed from Figure 2 points using the same
`float_precision="round_trip"` parser as the original plotter reproduce the
checked-in statistics to numerical precision, with panel sizes 2,454 / 2,485 /
838 / 843. The geometry eta2 join matches the Figure 2 x coordinate to a maximum
absolute difference below `1.2e-16`. Using pandas' default float parser perturbs
response ties; the analysis explicitly uses the original round-trip parser.
`parity_checks.tsv.gz` records the panel association and eta2 checks.
The script asserts these checks before starting new analyses: exact D/E/F/G counts
and anchor/split labels, disjoint original development/confirmation reactions,
published association parity, reaction-level eta2 evaluation-count parity,
complete reaction/pathway joins without label conflicts, and GroupKFold test
coverage with no pathway in both train and test. Each asserted check passes.

The released Figure 3 reaction/evaluation table is pre-aggregated over truth
cases and contains no `truth_pair_id`; therefore exact pair identities cannot be
re-exported or independently audited from this compact input. The stated 836
matched-pair denominator is retained from the frozen A23 design and is not
expanded into invented IDs. `eligible_evaluations.tsv.gz` preserves the 370
evaluation IDs per anchor used to compute the released reaction/evaluation
aggregates.

The original endpoint remains the fraction of non-tie truth cases satisfying
both `g_correct > 1e-12` and `g_info > 1e-12`, with primary lambda 0.25 and truth
and numerical-gain tolerances of `1e-12`. Descriptors and usefulness are averaged
according to the released reaction/evaluation summaries. A1 development,
confirmation, and A2 rows retain their own frozen eligibility. A1–A2 paired
comparisons use reaction intersections whose outcomes derive from the 836 matched
truth-pair identities per anchor; the truth pairs are not treated as biological
replicates. No HEX1 coordinate was substituted for the implemented glucose
observable.

## Experiment 1: matched descriptors and prediction

All descriptors and usefulness were finite for every original Figure 2 panel
row. Thus the matched descriptor denominators are 2,454 (A1 development), 2,485
(A2 development), 838 (A1 confirmation), and 843 (A2 held-out robustness). The
table below gives reaction-level Spearman rho; paired reaction-bootstrap 95%
descriptive intervals for each rho are in
`matched_descriptor_comparisons.tsv.gz`.

| Panel | eta2 | Direction entropy H_dir | Dominant-direction mass D_dir | Supported sign states K_dir | eta2 partial H_dir | eta2 partial H_dir + coverage |
|---|---:|---:|---:|---:|---:|---:|
| A1 development, n=2,454 | 0.8571 | 0.7894 | -0.7754 | 0.7701 | 0.5948 | 0.5154 |
| A2 development, n=2,485 | 0.8339 | 0.8562 | -0.8471 | 0.8196 | 0.4846 | 0.3984 |
| A1 confirmation, n=838 | 0.8351 | 0.8225 | -0.8086 | 0.8020 | 0.4925 | 0.4120 |
| A2 held-out robustness, n=843 | 0.8494 | 0.8661 | -0.8562 | 0.8461 | 0.4829 | 0.3864 |

For eta2 minus H_dir, matched correlation differences (95% reaction-bootstrap
intervals) are +0.0677 `[0.0501, 0.0863]` in A1 development, -0.0223
`[-0.0392, -0.0060]` in A2 development, +0.0126 `[-0.0197, 0.0419]` in A1
confirmation, and -0.0167 `[-0.0404, 0.0069]` in A2 held-out robustness.
Therefore eta2 has a clear matched association advantage over entropy in A1
development, but that ranking does not hold across anchors and held-out sets.
Partial rank association remains positive after entropy and coverage adjustment,
though it is smaller than the marginal eta2 association.

Dominant-direction mass is negatively associated in every population, so its
predictive strength is compared sign-aligned (`-D_dir`) and by absolute rho;
the original signed `D_dir` correlations remain in the table above. The paired
bootstrap difference below is `|rho(eta2)| - |rho(D_dir)|`:

| Panel | Scope | n | rho eta2 | rho D_dir (signed) | Difference in absolute rho [95% interval] |
|---|---|---:|---:|---:|---:|
| A1 development | Full | 2,454 | 0.8571 | -0.7754 | +0.0818 [0.0635, 0.1006] |
| A1 development | Mapped only | 2,272 | 0.8599 | -0.7931 | +0.0668 [0.0489, 0.0850] |
| A2 development robustness | Full | 2,485 | 0.8339 | -0.8471 | -0.0131 [-0.0295, 0.0033] |
| A2 development robustness | Mapped only | 2,303 | 0.8283 | -0.8551 | -0.0268 [-0.0436, -0.0100] |
| A1 original confirmation | Full | 838 | 0.8351 | -0.8086 | +0.0265 [-0.0056, 0.0598] |
| A1 original confirmation | Mapped only | 779 | 0.8342 | -0.8184 | +0.0158 [-0.0190, 0.0523] |
| A2 confirmation-split robustness | Full | 843 | 0.8494 | -0.8562 | -0.0068 [-0.0310, 0.0163] |
| A2 confirmation-split robustness | Mapped only | 784 | 0.8453 | -0.8579 | -0.0126 [-0.0363, 0.0119] |

Full and mapped-only bootstrap samples use their exact finite reaction
intersections. These intervals describe reaction-level computational stability,
not biological confidence.

Standardized Ridge (alpha=1) models used identical five-fold pathway GroupKFold
splits within each development panel and population scope, with scaling fit only
on each training fold. Alongside M0–M5, two post-hoc comparisons use the same
finite reactions and folds within each panel/scope: cubic entropy (`H_dir`,
`H_dir²`, `H_dir³`) with versus without eta2, and entropy plus non-tie coverage
with versus without eta2. A1 specification used development only. A1 confirmation
predictions used the fixed development fit; A2 fits are post-freeze robustness.

| Population | Model | MAE | Spearman(prediction, usefulness) | MAE improvement over M1 |
|---|---|---:|---:|---:|
| A1 development OOF, n=2,454 | M1 H_dir | 0.06251 | 0.7745 | 0 |
| A1 development OOF | M2 eta2 | 0.05183 | 0.8525 | 0.01068 |
| A1 development OOF | M3 H_dir + eta2 | 0.05017 | 0.8599 | 0.01234 |
| A1 development OOF | M5 H_dir + D_dir + eta2 | 0.04885 | 0.8648 | 0.01366 |
| A1 confirmation fixed fit, n=838 | M1 H_dir | 0.05608 | 0.8225 | 0 |
| A1 confirmation fixed fit | M3 H_dir + eta2 | 0.04854 | 0.8528 | 0.00754 |
| A1 confirmation fixed fit | M5 H_dir + D_dir + eta2 | 0.04674 | 0.8596 | 0.00933 |
| A2 development OOF, n=2,485 | M1 H_dir | 0.07173 | 0.8518 | 0 |
| A2 development OOF | M2 eta2 | 0.07749 | 0.8286 | -0.00576 |
| A2 development OOF | M3 H_dir + eta2 | 0.06445 | 0.8882 | 0.00728 |
| A2 development OOF | M5 H_dir + D_dir + eta2 | 0.06156 | 0.8934 | 0.01017 |
| A2 held-out fixed fit, n=843 | M1 H_dir | 0.06639 | 0.8661 | 0 |
| A2 held-out fixed fit | M3 H_dir + eta2 | 0.06061 | 0.8988 | 0.00578 |
| A2 held-out fixed fit | M5 H_dir + D_dir + eta2 | 0.05742 | 0.9051 | 0.00897 |

These predictive comparisons support incremental eta2 information in the A1
development model and in the combined A2 model, but not universal eta2
superiority: eta2 alone performs worse than entropy alone in A2 development.
Confirmation performance is post-hoc support for this model exercise, not a new
independent confirmation.

The additional model pairs report the improvement from adding eta2
(`baseline MAE − plus-eta2 MAE`; `plus-eta2 rho − baseline rho`). Each model pair
uses an identical eligible population and shared pathway folds:

| Population | Scope | n | Cubic entropy ΔMAE / Δrho | Entropy + coverage ΔMAE / Δrho |
|---|---|---:|---:|---:|
| A1 development | Full | 2,454 | 0.00997 / +0.07745 | 0.01117 / +0.09189 |
| A1 development | Mapped only | 2,272 | 0.00911 / +0.06473 | 0.01038 / +0.07901 |
| A2 development robustness | Full | 2,485 | 0.00519 / +0.03467 | 0.00565 / +0.02425 |
| A2 development robustness | Mapped only | 2,303 | 0.00471 / +0.02837 | 0.00504 / +0.02036 |
| A1 original confirmation, fixed development fit | Full | 838 | 0.00542 / +0.03513 | 0.00673 / +0.05870 |
| A1 original confirmation, fixed development fit | Mapped only | 779 | 0.00492 / +0.03017 | 0.00591 / +0.04962 |
| A2 confirmation-split robustness | Full | 843 | 0.00496 / +0.03255 | 0.00513 / +0.02773 |
| A2 confirmation-split robustness | Mapped only | 784 | 0.00499 / +0.03023 | 0.00491 / +0.02653 |

Adding eta2 improves both reported OOF metrics in these fixed-alpha model pairs
for the development populations and the fixed-fit confirmation populations.
These are post-hoc model results, not a newly frozen confirmatory hypothesis.
Full model-level MAE, Spearman, baseline, and sample-size records are in
`predictive_models.tsv.gz`.

## Experiment 2: pathway dependence

The pathway labels contain **49 mapped metabolic subsystems plus
`UNMAPPED_PATHWAY`** where all mapped subsystems occur. A1/A2 development each
have 49 mapped groups plus the unmapped group (50 groups total); each confirmation
split has 48 mapped groups plus unmapped (49 total). Unmapped reactions are
retained in FULL results and excluded only in MAPPED_ONLY sensitivity. Unmapped
counts are 182 in each development panel and 59 in each confirmation panel.
Panel denominators are 2,454/2,485/838/843 FULL and 2,272/2,303/779/784
MAPPED_ONLY. Confirmation mapped pathway sizes range 1–144 (median 6.5); 17
mapped pathways have fewer than five reactions and six are singletons.
Within-pathway associations are reported for pathways with at least five finite
reactions; smaller groups remain listed and flagged in
`pathway_within_associations.tsv.gz`.

The eta2 association remains positive under both reaction-level pathway-block
bootstrap (2,000 replicates) and equal-pathway summaries:

| Panel | Reaction-level rho | Pathway-block interval | Equal-pathway rho | Equal-pathway bootstrap interval |
|---|---:|---:|---:|---:|
| A1 development | 0.8571 | [0.8117, 0.8771] | 0.8841 | [0.7872, 0.9285] |
| A2 development | 0.8339 | [0.7799, 0.8720] | 0.8351 | [0.6743, 0.9205] |
| A1 confirmation | 0.8351 | [0.7612, 0.8733] | 0.8216 | [0.6701, 0.9069] |
| A2 held-out robustness | 0.8494 | [0.8020, 0.8776] | 0.8283 | [0.6850, 0.9175] |

Pathway grouping therefore does not remove the positive eta2 association in these
released populations. It does not guarantee statistical independence; broad
intervals, especially on confirmation panels with small pathway sizes, reflect
the limited number and uneven size of groups. After excluding
`UNMAPPED_PATHWAY`, eta2 rho remains positive: 0.8599 (A1 development), 0.8283
(A2 development), 0.8342 (A1 confirmation), and 0.8453 (A2 confirmation-split
robustness). The mapped-only equal-pathway summaries are 0.8795, 0.8290, 0.8225,
and 0.8286; all corresponding 2,000-replicate block-bootstrap intervals remain
positive. Grouped predictive scores are in `predictive_models.tsv.gz`, and exact
fold identities are in `fold_assignments.tsv.gz`.

## Superseded Experiments 3–4 results (initial implementation)

The numerical findings in this historical section were superseded after review
identified pooled truth-pair aggregation, confirmation-fold leakage, and a
changing eligible truth-pair population in ESS comparisons. Do not cite these
values as the current results. They are retained to preserve the analysis
history; corrected results appear in the next section.

These analyses use the authorized, read-only scientific checkout. The four
primary source files passed their expected SHA-256 checks. The carbon-normalized
flux cache is used for candidate-vector geometry and frozen prediction outcomes;
A1 glucose weighting uses raw panel `max(-EX_glc__D_e, 0)`, and A2 lactate
weighting uses raw panel `max(-LDH_L, 0)`. Those operators reproduce the frozen
candidate observables. HEX1 is used only for synthetic-truth selection and
geometry diagnostics. No reconstruction, sampling, frozen source data, or
production output was changed.

The source panel contains 4,181 reactions. The PL1 descriptor inventory is
4,180 (excluding HEX1); weak-cue outcomes use the distinct 4,179-reaction
population (excluding HEX1 and LDH_L). The 400 evaluations are 16 method/RNA
settings crossed with 25 CT2A–GL261 mouse pairs. They are computational cases,
not independent biological replicates. Reaction-level associations and grouped
cross-validation are post-hoc summaries; correlated reactions and pathways
limit independent-reaction interpretations.

### Experiment 3A: exclude the selected whole truth vectors from geometry

The recomputed 20×20 descriptors reproduced the frozen PL1/A21 reaction-summary
medians and finite counts on the shared 4,179-reaction population. The maximum
absolute median difference was `1.64e-15`; a four-method deterministic 19×19
pilot also passed the frozen PL1 descriptor implementation to floating-point
precision. The full comparison retained 836 distinct evaluable truth pairs in
both A1 and A2 and removed each selected whole candidate vector before
renormalizing the corresponding 19-vector weights.

For every reaction, 20×20 and 19×19 geometry summaries and their differences
are reported over the same distinct truth-pair population. Median absolute
changes were 0.0098 for η², 0.0056 for direction entropy, and 0.0034 for
dominant-direction mass. The η²–usefulness comparisons use the same finite-η²
reaction set for both geometries within each arm: 3,289 reactions in A1 and
3,326 in A2. These 19×19 descriptors depend on which truth vectors were
selected and are therefore **truth-exclusion sensitivity measures**, not
independently available pre-cue predictors.

| Arm | Geometry | Partial rank η² given entropy | Grouped-CV MAE improvement, all | Confirmation-holdout improvement |
|---|---|---:|---:|---:|
| A1 | Original 20×20 | 0.395 | 0.00409 | 0.00297 |
| A1 | Truth-excluded 19×19 | 0.452 | 0.00443 | 0.00317 |
| A2 | Original 20×20 | 0.288 | 0.00254 | 0.00287 |
| A2 | Truth-excluded 19×19 | 0.288 | 0.00324 | 0.00357 |

η² retains positive incremental association beyond entropy for both geometries
and arms. Grouped-CV MAE improves in these matched finite-reaction populations,
although the size varies by arm and split. This supports robustness to removing
the selected whole truth vectors from geometry calculation; it does not make the
truth-dependent descriptors available before truth selection.

### Experiment 3B: alternative synthetic-truth selections

The frozen q10/q50/q90 truth selection was compared with weighted HEX1
q20/q50/q80 selection and four seeded weighted random selections (seeds
20271028–20271031; three whole-vector selections per candidate pool, with
replacement). Candidate IDs were selected once and reused in paired A1–A2
outcomes. Duplicate selections are recorded. Each scenario has 1,200 selection
rows, of which 1,110 are evaluable. After collapsing repeated identities to
distinct evaluation/truth-pair cases, support is 836 for the frozen selection,
803 for q20/q50/q80, and 721–742 across the four random selections. The reduced
alternative counts are reported rather than treating duplicate pairings as
additional cases.

Across these six scenarios, the partial rank association of the frozen 20×20
η² descriptor with usefulness, conditional on entropy, ranges from 0.836–0.848
in A1 and 0.871–0.884 in A2. Adding η² reduces all-reaction grouped-CV MAE by
0.00021–0.00032 in A1 and 0.01294–0.01329 in A2. On the confirmation holdout,
A1’s MAE change is slightly unfavorable (about −0.0003) across scenarios,
whereas A2’s improvement remains positive (0.0076–0.0092). Thus the association
is stable under these alternative selections, while A1’s incremental predictive
gain is small and does not consistently improve the holdout score.

These scenarios resample the same 20 frozen vectors and weights. They alter the
synthetic benchmark population; they do not provide independent biological
validation.

### Experiment 4: global ESS calibration sensitivity

A1 and A2 weights were recalibrated with the original rank-based operators at
the tumor–mouse pool level of 320 candidates. Recalibration at target ESS 20
exactly reproduced the frozen weights for all ten mice (maximum absolute weight
difference 0 for both arms). Raw candidate-panel observables were used for both
anchor operators. λ remained 0.25 and the original truth identities were
retained.

At ESS 10, conditional support remained in 137/160 A1 and 156/160 A2
method/context strata; some strata had zero conditional support and were not
redistributed. The common paired truth-pair support is 428/836 at ESS 10 and
836/836 at ESS 20 and ESS 40. Conditional ESS and product-weight ESS are
provided for all settings.

| Global ESS | Mean usefulness A1 | Mean usefulness A2 | Correct-minus-wrong gain A1 | Correct-minus-wrong gain A2 |
|---:|---:|---:|---:|---:|
| 10 | 0.1122 | 0.2209 | 0.00200 | 0.00235 |
| 20 | 0.1692 | 0.2500 | 0.03309 | 0.05454 |
| 40 | 0.2506 | 0.2924 | 0.15118 | 0.14630 |

Wrong-cue mean gains remain negative at each target. A2 has higher usefulness
at all three targets, and a larger correct-minus-wrong separation at ESS 10 and
20. The separation reverses slightly at ESS 40, where A1 is larger. Thus the
larger A2 gain is not qualitatively stable across the full calibration range.
The η²–usefulness Spearman association remains positive at all targets
(0.823–0.864 across arms); partial rank η² given entropy also remains positive
(0.476–0.586). These comparisons are calibration sensitivities, not estimates
of the isolated causal effect of adding a measurement.

### Updated answers

1. **Does η² retain incremental information when truth vectors are excluded
   from geometry?** Yes in the matched finite-reaction populations for both
   arms. Treat the 19×19 values as truth-exclusion sensitivity descriptors.
2. **Are associations sensitive to alternative truth selection?** The partial
   η² association remains positive across the q20/q50/q80 and four seeded
   selections. A1’s incremental MAE improvement is very small and is not
   positive on the confirmation holdout; A2’s is positive in these scenarios.
3. **Is the larger A2 gain robust to weighting calibration?** Not across the
   entire ESS range. A2 usefulness is higher at ESS 10, 20, and 40, but its
   correct-minus-wrong separation is lower than A1 at ESS 40.
4. **Do these results provide independent biological confirmation?** No. They
   reuse the frozen candidates, mice, and synthetic truth-generation framework.
   The 400 evaluations and repeated reaction/truth-pair rows are not independent
   biological replicates.

## Corrected Experiments 3–4: current findings

All primary endpoints now follow the frozen two-stage rule: distinct truth-pair
cases are averaged within each evaluation after excluding ties, then defined
evaluations receive equal weight in the reaction summary. Evaluation counts and
finite denominators are retained. Directional usefulness still requires both
correct-cue gain and information advantage to exceed `1e-12`; λ is 0.25. The
synthetic unequal-case-count test confirms this differs from direct pooling.
Baseline case outcomes use round-trip parsing of frozen PL2B/A22 values. This
avoids near-tolerance flips from recomputing binary endpoints through a different
floating-point path.

Predictive models now use fixed α=1 Ridge, with scaling and fitting restricted
to development reactions. Five-fold pathway GroupKFold is used for development
out-of-fold scores. The complete development fit is then applied unchanged to
confirmation reactions. Confirmation has already informed the broader
manuscript revision, so these are post-hoc robustness comparisons, not
prospective independent validation.

### Experiment 3A

Original 20×20 geometry parity passes with maximum median discrepancy
`1.64e-15`. The 836 frozen matched distinct truth pairs produce 1,546,230
evaluation×reaction rows per arm; round-trip frozen-case endpoint parity passes
for A1/PL2B and A2/A22. The weak-cue universe remains 4,179 reactions. The
19×19 descriptor excludes each selected whole vector and renormalizes the
remaining weights. It depends on truth identity and remains a
**truth-exclusion sensitivity descriptor**, not a predictor available before
truth selection.

η²–usefulness Spearman associations remain positive for both geometries:
development/confirmation 0.857/0.835 (A1, 20×20), 0.887/0.856 (A1, 19×19),
0.834/0.849 (A2, 20×20), and 0.827/0.839 (A2, 19×19). Partial rank
associations controlling for entropy are also positive (A1 0.414/0.343 and
0.485/0.406; A2 0.283/0.310 and 0.256/0.301 for development/confirmation).
Adding η² to the entropy model reduces MAE in development grouped CV and in
confirmation from the development fit: 0.0044/0.0027 and 0.0055/0.0040 for A1
(20×20/19×19), and 0.0032/0.0029 and 0.0028/0.0029 for A2. Thus the
truth-exclusion association and incremental score persist, though the A2
partial association weakens slightly under exclusion.

### Experiment 3B

The six scenarios retain their own distinct-pair populations: 836 frozen
q10/q50/q90 pairs, 803 q20/q50/q80 pairs, and 742, 727, 739, and 721 pairs for
random seeds 20271028–20271031. Each has 1,200 selection rows and 1,110
evaluable selection rows; duplicate vector selections are recorded and are not
counted as additional distinct truth pairs. Each scenario retains 370 distinct
contributing evaluations. The same selected candidate IDs are used in paired
A1–A2 comparisons; every scenario has the same 3,291 shared finite reactions
(2,453 development and 838 confirmation). The exported tables include exact
scenario-level pair and evaluation counts.

η² associations and partial associations controlling for entropy are positive
for every scenario. Across the four alternatives and the frozen baseline, η²
Spearman ranges 0.830–0.870 (A1) and 0.915–0.922 (A2) in development, with
partial rank ranges 0.841–0.848 and 0.893–0.898, respectively. Incremental
development-CV MAE improvement is small for A1 (0.00077–0.00089) and larger for
A2 (0.01618–0.01661). On confirmation, η² worsens A1 MAE by 0.00070–0.00078
and improves A2 MAE by 0.01165–0.01225. Therefore the η² association is
selection-robust, but incremental prediction is not equally robust across
anchors: the A1 confirmation benefit disappears, while A2’s remains positive.

These selections resample the same frozen vectors and weights and change the
synthetic benchmark population. They do not add biological validation.

### Experiment 4

ESS20 weight parity passes exactly for all ten mouse pools. The frozen original
truth identities and raw candidate-panel observables are retained; glucose uses
`max(-EX_glc__D_e, 0)` and lactate uses `max(-LDH_L, 0)`. Calibration remains at
the global 320-candidate tumor–mouse pool. The ESS20 endpoint table matches the
corrected PL2B/A22 outcomes for both arms, including usefulness and pair/tie
counts. An initial vector-recomputed endpoint differed for one evaluation×
reaction per arm near the threshold; all ESS20 results below use the frozen
case decisions.

Each arm has 3,625 finite reaction-level endpoint summaries. Paired A1–A2
predictive models use the identical finite reaction intersections at each target
(3,275, 3,291, and 3,331 reactions for ESS10, ESS20, and ESS40 respectively).

| Population | ESS | Truth pairs | Evaluations | A1 correct−wrong | A2 correct−wrong | A1 usefulness | A2 usefulness |
|---|---:|---:|---:|---:|---:|---:|---:|
| ESS-specific | 10 | 428 | 164 | 0.00142 | 0.00209 | 0.0919 | 0.1940 |
| ESS-specific | 20 | 836 | 370 | 0.02018 | 0.04975 | 0.1386 | 0.2156 |
| ESS-specific | 40 | 836 | 370 | 0.10129 | 0.11020 | 0.2117 | 0.2541 |
| Fixed common | 10 | 428 | 164 | 0.00142 | 0.00209 | 0.0919 | 0.1940 |
| Fixed common | 20 | 428 | 164 | 0.03820 | 0.09422 | 0.1585 | 0.2349 |
| Fixed common | 40 | 428 | 164 | 0.19237 | 0.21144 | 0.2303 | 0.2641 |

The fixed population is the exact intersection of truth identities supported at
ESS10/20/40 in both arms (428 pairs). A2’s correct-minus-wrong gain separation
is greater at all targets in both population definitions. The previously
reported ESS40 A1-over-A2 reversal does **not** persist after equal-evaluation
aggregation and frozen endpoint restoration. ESS10 has conditional support in
137/160 A1 and 156/160 A2 context strata; conditional ESS and product-weight
ESS are in the output tables. Across ESS settings, η² associations remain
positive; the corrected development-fitted confirmation MAE improvement ranges
0.0050–0.0104 (A1) and 0.0053–0.0105 (A2) across the reported ESS/population
comparisons. Paired models use identical finite reactions (3,275, 3,291, and
3,331 by ESS target).

These contrasts describe calibration sensitivity, not the isolated causal
effect of adding lactate measurement. Reactions remain correlated features, and
reaction-level rows are not biological replicates.

### Assessment and manuscript use

1. η² retains a positive association after truth-vector exclusion.
2. Development-fitted η² improves confirmation MAE in 3A and across ESS
   calibration in Experiment 4. In 3B it improves A2 but worsens A1 confirmation
   MAE, so avoid a universal incremental-prediction claim.
3. Positive η² associations are robust across all alternative truth selections;
   prediction gains are anchor-dependent.
4. The A2–A1 gain-separation advantage is stable across ESS targets after
   correcting the aggregation and controlling truth-pair eligibility.
5. The ESS40 reversal disappears on both ESS-specific and fixed-common
   populations.
6. The defensible findings are methodological sensitivity results about
   synthetic candidate benchmarks. None establishes biological validity.
7. A manuscript may report the positive association and the A2-specific
   robustness, with the A1 3B failure, calibration dependence, post-hoc status,
   and truth-dependent geometry clearly stated. Do not claim independent
   biological validation.

## Outputs, tests, and suggested manuscript follow-up

- `analyses/reviewer_checks/run_reviewer_checks.py` implements the released-table
  analysis; `tests/test_reviewer_checks.py` contains nine synthetic tests.
- Compressed tables: `parity_checks.tsv.gz`, `reaction_level_analysis.tsv.gz`,
  `eligible_evaluations.tsv.gz`, `matched_descriptor_comparisons.tsv.gz`,
  `paired_anchor_comparisons.tsv.gz`, `correlation_strength_comparisons.tsv.gz`,
  `predictive_models.tsv.gz`, `pathway_within_associations.tsv.gz`,
  `pathway_blocked_robustness.tsv.gz`, `fold_assignments.tsv.gz`, and
  `development_oof_predictions.tsv.gz`.
- `manifest.json` and the existing tables retain Experiments 1–2 provenance.
  Corrected results are in `reproduced/reviewer_checks/corrected/`; these files
  are gitignored. The compact tracked reviewer bundle is
  `analyses/reviewer_checks/exports/candidate_sensitivity_review_bundle.zip`.
  No candidate vectors or raw case rows are included.
- The candidate-sensitivity suite and complete test suite passed. The existing Experiments 1–2
  numerical results and test records above are preserved. No manuscript or
  supplementary file was edited, and no reconstruction, OptGP sampling, or
  GapSplit rerun occurred.

Suggested manuscript/supplement follow-up: (i) retain the exact Figure 2
round-trip float parsing convention in any derived correlation table; (ii)
report matched eta2/entropy comparisons with panel-specific denominators; (iii)
present the truth-exclusion analysis explicitly as a post-hoc sensitivity
measure; and (iv) report that the corrected ESS40 reversal disappears on both
ESS-specific and fixed-common populations, alongside the loss of conditional
support at ESS10. These sensitivity analyses do not warrant claims of
independent biological validation. This report does not edit manuscript or
supplementary files.
