# Geometry usefulness regression v2: scientific summary

This is a post hoc computational follow-up on the frozen Development/Confirmation split. The Confirmation reactions have already been examined in v1; these results are not a new independent confirmation experiment. Bootstrap intervals are conditional on the fitted models and do not account for model or cutoff selection.

## 1. Usefulness-frequency prediction and uncertainty
- **A1:** M3−M4 MAE difference=0.0067 (95% reaction-bootstrap CI 0.0040 to 0.0094; subsystem-bootstrap CI 0.0031 to 0.0098). Relative MAE reduction=12.03% (reaction CI 7.38% to 16.36%). Spearman M4−M3=0.0580 (reaction CI 0.0315 to 0.0827).
- **A2:** M3−M4 MAE difference=0.0051 (95% reaction-bootstrap CI 0.0025 to 0.0076; subsystem-bootstrap CI -0.0002 to 0.0100). Relative MAE reduction=7.68% (reaction CI 3.86% to 11.28%). Spearman M4−M3=0.0275 (reaction CI 0.0167 to 0.0384).

Reaction-level resampling treats reactions as exchangeable and does not fully represent metabolic dependence. Subsystem resampling is a sensitivity analysis over 49 represented clusters per anchor; the clusters are highly imbalanced (finite Confirmation sizes: median 7, maximum 144), so those intervals should be read cautiously.

## 2. Top-k usefulness prioritization
- **A1:** top 10% M3 0.2709 vs M4 0.2697; top 20% M3 0.2528 vs M4 0.2521. Paired reaction-bootstrap changes (M4−M3) were top10=-0.0011 (95% CI -0.0150 to 0.0151) and top20=-0.0006 (CI -0.0112 to 0.0116; subsystem CI -0.0154 to 0.0206). Thus the evaluated top-k metrics do not show a material gain. These previously inspected Confirmation rankings are descriptive and do not establish independent prioritization performance.
- **A2:** top 10% M3 0.4109 vs M4 0.4178; top 20% M3 0.3829 vs M4 0.4013. Paired reaction-bootstrap changes (M4−M3) were top10=0.0069 (95% CI -0.0085 to 0.0179) and top20=0.0184 (CI 0.0048 to 0.0282; subsystem CI -0.0040 to 0.0414). Thus the top-20% usefulness lift is positive under reaction resampling but its subsystem interval includes zero. These previously inspected Confirmation rankings are descriptive and do not establish independent prioritization performance.

## 3. Mean correct-cue gain prediction
- **A1:** Gain-target M4 confirmation MAE=0.0174, RMSE=0.0463, Spearman=0.1625; M0 MAE=0.0173. M4 versus M0 MAE improvement=-0.0001; M4 versus M3 changes were MAE improvement -0.0002, RMSE improvement 0.0000, and Spearman change 0.0515. Across models, M4 does not improve gain MAE over M0 under either anchor; the gain target therefore has weak overall predictive support in this cohort despite some rank associations and top-k contrasts.
- **A2:** Gain-target M4 confirmation MAE=0.0445, RMSE=0.1265, Spearman=0.2743; M0 MAE=0.0412. M4 versus M0 MAE improvement=-0.0033; M4 versus M3 changes were MAE improvement -0.0033, RMSE improvement 0.0008, and Spearman change 0.3063. Across models, M4 does not improve gain MAE over M0 under either anchor; the gain target therefore has weak overall predictive support in this cohort despite some rank associations and top-k contrasts.

The gain target reuses v1's signed `mean_correct_gain`: within each evaluation, correct-cue gain is averaged over non-tie evaluable truth pairs; the reaction target is the arithmetic mean of finite evaluation-level values. It is not clipped or reweighted by the number of pairs per evaluation.

Gain-targeted and usefulness-frequency-targeted rankings are compared on identical per-anchor reaction populations. For M4, the actual top-k comparison is:
- **A1:** top 10% gain ranking: mean gain 0.0161, positive fraction 0.7976, usefulness 0.1008; frequency ranking: 0.0000, 0.9286, 0.2697; top 20% gain ranking: mean gain 0.0177, positive fraction 0.8214, usefulness 0.1101; frequency ranking: 0.0005, 0.9464, 0.2521.
- **A2:** top 10% gain ranking: mean gain 0.0782, positive fraction 0.9529, usefulness 0.3142; frequency ranking: 0.0038, 0.8941, 0.4178; top 20% gain ranking: mean gain 0.0640, positive fraction 0.9053, usefulness 0.2889; frequency ranking: 0.0021, 0.8876, 0.4013.

The tables also compare each ranking with random selection and the entropy-only ranking. The endpoints remain distinct: usefulness frequency measures how often outcomes are beneficial, while mean correct-cue gain is a signed magnitude. Positive-gain fraction is reported separately from mean gain.

## 4. Calibration
- **A1:** M1 bias overall/top10/top20=-0.0008/-0.0329/-0.0179; M3 bias overall/top10/top20=-0.0010/-0.0326/-0.0161; M4 bias overall/top10/top20=-0.0035/-0.0481/-0.0321. Bias is observed minus predicted usefulness; calibration is descriptive and does not establish external generalization.
- **A2:** M1 bias overall/top10/top20=0.0015/-0.0476/-0.0314; M3 bias overall/top10/top20=0.0015/-0.0586/-0.0332; M4 bias overall/top10/top20=-0.0029/-0.0617/-0.0301. Bias is observed minus predicted usefulness; calibration is descriptive and does not establish external generalization.

## 5. Anchor comparison and manuscript implication
The M3/M4 usefulness-frequency accuracy gains occur under both anchors, while top-k usefulness changes are small for A1 and more favorable for A2 at 20%; subsystem resampling makes that A2 top-20% interval uncertain. Gain-targeted M4 ranking selects higher realized mean gain than matched frequency-targeted M4 ranking under both anchors, especially A2, while its selected usefulness frequency is lower under A2. The mismatch shows why practical prioritization cannot be inferred from regression accuracy alone.

For the manuscript, these results support an anchor-specific predictive association between η²_dir and beneficial-outcome frequency conditional on the fitted linear models. They do not establish reliable prediction of signed mean magnitude gain or prove the theoretical direction–magnitude coupling mechanism. Keep these claims separate; η²_dir is not universally superior across targets and prioritization metrics.

All modeling used Development outcomes only for fitting, fold-local standardization, and alpha selection. No Confirmation outcome was used to fit or select gain models. No GEM reconstruction, flux sampling, truth generation, or case-level aggregation was rerun.
