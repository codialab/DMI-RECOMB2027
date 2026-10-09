# Reviewer-requested sensitivity analyses

Status: completed for Experiments 1–2 from released compact tables; Experiments
3–4 deferred because their required upstream artifacts are not included in this
repository. All results below are post-hoc robustness analyses. No manuscript or
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

## Experiments 3–4 and explicit answers

Experiments 3 and 4 were not run. This checkout contains compact summaries, not
the original whole candidate vectors, exact candidate IDs/truth-pair selection
registries, candidate-level raw observables, per-mouse ranks, or the 320-candidate
weight inputs. Relevant upstream `outputs/` artifacts and bulk `external_data/`
inputs are omitted from the public reproducibility checkout. No `--source-root`
was supplied or accessed. The manifest records these missing artifact classes;
providing the original vector cache plus frozen identities/weights would enable
the 19×19 truth-excluded descriptor pilot and alternative truth selections.
Providing original observables, ranks, and baseline weighting inputs at the
320-candidate tumor–mouse pool level would enable ESS 10/20/40 calibration on
the fixed original truth pairs. Until then, truth-vector inclusion sensitivity,
truth-selection sensitivity, and calibration robustness remain unknown.

1. **Does eta2 add predictive value beyond linear entropy?** Yes in A1 and A2
   development grouped CV; adding eta2 reduces MAE by 0.01234 and 0.00728,
   respectively. A2 remains post-freeze robustness, and marginal rho rankings
   are not uniform across anchors.
2. **Does that improvement remain with nonlinear entropy?** Yes in these
   post-hoc fixed-alpha comparisons: eta2 reduces MAE by 0.00997 in A1 and
   0.00519 in A2 development, with positive Spearman changes.
3. **Does the improvement remain after accounting for non-tie coverage?** Yes
   in development grouped CV: adding eta2 to H_dir plus coverage reduces MAE by
   0.01117 in A1 and 0.00565 in A2.
4. **Does the positive association survive exclusion of unmapped pathways?**
   Yes. Eta2 reaction and equal-pathway associations and block intervals stay
   positive in all four mapped-only populations.
5. **Are conclusions consistent across A1 and A2?** Incremental model
   performance is directionally consistent, but descriptor rankings are not:
   entropy rho exceeds eta2 rho in A2 development and held-out robustness.
   A2 is never independent confirmation.
6. **Which conclusions are supported only by post-hoc analyses?** The nonlinear
   entropy, coverage-adjusted, mapped-only, all A2, and model-based confirmation
   comparisons are post-hoc and do not establish independent biological
   confirmation.
7. **Is eta2 sensitive to including truth vectors in descriptor construction?**
   Unknown; Experiment 3A could not run without original vectors and truth IDs.
8. **Is the result sensitive to synthetic-truth selection?** Unknown; Experiment
   3B could not run without those same frozen candidate identities and vectors.
9. **Is the A2 gain difference robust to weighting calibration?** Unknown;
   Experiment 4 could not run without original calibration observables, ranks,
   and weights. Current A2 results do not isolate a causal effect of adding a
   second measurement.

The results permit stating that eta2 is positively associated with the frozen
usefulness response under reaction and pathway-level analyses, and that eta2
adds predictive value beyond linear entropy, cubic entropy, and entropy plus
coverage in these post-hoc grouped model comparisons. They do **not**
establish universal descriptor superiority, biological replication across the
400 evaluations, truth-selection independence, robustness to truth-vector
exclusion or weighting recalibration, or independent A2 confirmation. The
Figure 2 point-table/statistics-table association parity passed under the
original float parsing convention.

## Outputs, tests, and suggested manuscript follow-up

- `analyses/reviewer_checks/run_reviewer_checks.py` implements the released-table
  analysis; `tests/test_reviewer_checks.py` contains nine synthetic tests.
- Compressed tables: `parity_checks.tsv.gz`, `reaction_level_analysis.tsv.gz`,
  `eligible_evaluations.tsv.gz`, `matched_descriptor_comparisons.tsv.gz`,
  `paired_anchor_comparisons.tsv.gz`, `correlation_strength_comparisons.tsv.gz`,
  `predictive_models.tsv.gz`, `pathway_within_associations.tsv.gz`,
  `pathway_blocked_robustness.tsv.gz`, `fold_assignments.tsv.gz`, and
  `development_oof_predictions.tsv.gz`.
- `manifest.json` records input hashes, seeds, denominators, deferred input
  requirements, and post-hoc interpretation limits. No figure was added because
  the tables directly show the comparisons and intervals.
- Tests: **9 passed**; `python reproduce.py validate` passed all eight scientific
  checks. Full numerical analysis completed from released tables;
  no solver, reconstruction, OptGP sampling, or GapSplit rerun occurred.

Suggested manuscript/supplement follow-up: (i) retain the exact Figure 2
round-trip float parsing convention in any derived correlation table; (ii)
report the matched eta2/entropy model comparison and its panel-specific limits;
(iii) add a pathway-blocked robustness table with reaction and pathway
denominators and explicitly non-biological interpretation; and (iv) label
Experiments 3–4 as unavailable until the required original source artifacts can
be analyzed. This report does not edit manuscript or supplementary files.
