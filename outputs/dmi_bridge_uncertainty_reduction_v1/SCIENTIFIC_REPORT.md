# Weak Directional Cues and Magnitude Uncertainty

## Methods

For each held-out truth pair and reaction, M=|Δv| is evaluated over the 19×19 held-out anchor-weighted sample-pair support after excluding both truth samples and renormalizing exactly as in the frozen PL2 holdout. The tie threshold is 1×10⁻¹². A1 reuses the frozen BIO-0 strong-anchor weights; A2 reuses the frozen A20 dual-anchor weights. No flux vectors or candidate identities are regenerated.

The primary uncertainty statistic is V=Σwᵢ(Mᵢ−μ)². Absolute reduction is V₀−Vλ; relative reduction is 1−Vλ/V₀ and is omitted when V₀≤100ε·max(E[M²],E[M]²). The threshold and excluded counts are recorded per case. Secondary statistics are weighted SD, inverse-CDF weighted IQR, central 90% inverse-CDF width, and discrete weighted CRPS. Weighted quantiles use the first sorted support point whose cumulative normalized mass reaches q, matching PL1.

For sign state s∈{−1,0,+1}, let πₛ be baseline mass and let μₛ and σ²ₛ be conditional magnitude mean and variance. A cue c changes only mixture masses: π'ₛ=cₛπₛ/Σₜcₜπₜ, where cₛ=exp(λ d s) and d is cue direction. Thus V'=Σₛπ'ₛσ²ₛ+Σₛπ'ₛ(μₛ−μ')². Conditional magnitude distributions do not change. The support remains unchanged whenever all baseline weights are positive; variance and quantile widths can still shrink or expand as direction-specific mixture mass shifts.

Accuracy is the unchanged PL2 expected-absolute-magnitude prediction gain. CRPS improvement is baseline CRPS minus cue CRPS. Neither a mean shift nor a directional entropy change is treated as uncertainty reduction; all three endpoints are analyzed separately.

CRPS uses CRPS(F,y)=E|X−y|−½E|X−X′| for the weighted discrete predictive distribution. Conditional sign-group expectations and cross-group absolute-distance terms are computed once per truth case and reused across both cue directions and all λ values; cue CRPS then uses only a three-state mixture calculation. The synthetic benchmark precedes any empirical run.

## Prospective geometry prediction

Correct-cue and incorrect-cue relative variance reduction are separate prediction targets; both use the identical complete-case reaction population and frozen Development/Confirmation split for A1 and A2. Prespecified ridge models compare η²_dir alone with η²_dir+H_dir, η²_dir+dominant-direction mass, η²_dir+non-tie coverage, η²_dir+supported sign states, and the full baseline geometry set. Scaling is fit on Development only. Confirmation usefulness and prediction gains never enter as prospective predictors.

## Explanatory association

Observed usefulness associations are descriptive and are kept separate from prospective model fitting. EXPLANATORY_USEFULNESS_ASSOCIATION.tsv reports rank association of η²_dir with uncertainty reduction before and after adjusting for observed directionally-useful outcome. It does not establish prospective prediction or experimental calibration.

## Results

Matched distinct truth pairs: 836; common reactions: 4179; case rows per anchor: 3493644. The complete-case prospective population contains 3291 reactions ({'CONFIRMATION_HOLDOUT': 838, 'DEVELOPMENT': 2453}).

### Variance change frequencies

- A1 correct: decreased 315824/1440313 (21.9%), increased 154149, negligible 970340; near-zero relative exclusions 41423.
- A1 incorrect: decreased 169914/1440313 (11.8%), increased 302379, negligible 968020; near-zero relative exclusions 41423.
- A2 correct: decreased 456548/1440313 (31.7%), increased 222700, negligible 761065; near-zero relative exclusions 5200.
- A2 incorrect: decreased 242051/1440313 (16.8%), increased 438205, negligible 760057; near-zero relative exclusions 5200.

### Prediction accuracy versus uncertainty

- A1 correct: both improve=139733; overconfidence risk=127087; accuracy improves as dispersion increases=104349; both worsen=37693; negligible/tied=1031451.
- A1 incorrect: both improve=40212; overconfidence risk=116860; accuracy improves as dispersion increases=121491; both worsen=132380; negligible/tied=1029370.
- A2 correct: both improve=227129; overconfidence risk=192825; accuracy improves as dispersion increases=133189; both worsen=76904; negligible/tied=810266.
- A2 incorrect: both improve=80704; overconfidence risk=148283; accuracy improves as dispersion increases=182603; both worsen=220396; negligible/tied=808327.

CRPS improvement is a separate proper-score endpoint and does not define variance reduction or frozen prediction-error gain.

### Prospective Confirmation scores

- A1 correct_relative_variance_reduction development_mean_baseline: Confirmation MAE=0.01446, RMSE=0.02114, R²=-0.009703, MAE improvement vs η²-only=nan.
- A1 correct_relative_variance_reduction eta2: Confirmation MAE=0.01487, RMSE=0.021308, R²=-0.02579, MAE improvement vs η²-only=0.
- A1 correct_relative_variance_reduction eta2_plus_entropy: Confirmation MAE=0.014976, RMSE=0.021583, R²=-0.05247, MAE improvement vs η²-only=-0.00010584.
- A1 correct_relative_variance_reduction all_geometry: Confirmation MAE=0.014782, RMSE=0.020935, R²=0.009858, MAE improvement vs η²-only=8.7561e-05.
- A1 incorrect_relative_variance_reduction development_mean_baseline: Confirmation MAE=0.018973, RMSE=0.025287, R²=-0.00913, MAE improvement vs η²-only=nan.
- A1 incorrect_relative_variance_reduction eta2: Confirmation MAE=0.017571, RMSE=0.024452, R²=0.05637, MAE improvement vs η²-only=0.
- A1 incorrect_relative_variance_reduction eta2_plus_entropy: Confirmation MAE=0.017622, RMSE=0.024587, R²=0.04594, MAE improvement vs η²-only=-5.0179e-05.
- A1 incorrect_relative_variance_reduction all_geometry: Confirmation MAE=0.015938, RMSE=0.023126, R²=0.1559, MAE improvement vs η²-only=0.001633.
- A2 correct_relative_variance_reduction development_mean_baseline: Confirmation MAE=0.01086, RMSE=0.014968, R²=-6.049e-06, MAE improvement vs η²-only=nan.
- A2 correct_relative_variance_reduction eta2: Confirmation MAE=0.0089776, RMSE=0.013899, R²=0.1377, MAE improvement vs η²-only=0.
- A2 correct_relative_variance_reduction eta2_plus_entropy: Confirmation MAE=0.0088756, RMSE=0.013976, R²=0.1282, MAE improvement vs η²-only=0.00010198.
- A2 correct_relative_variance_reduction all_geometry: Confirmation MAE=0.0088835, RMSE=0.01373, R²=0.1586, MAE improvement vs η²-only=9.4117e-05.
- A2 incorrect_relative_variance_reduction development_mean_baseline: Confirmation MAE=0.013623, RMSE=0.018371, R²=-0.0004033, MAE improvement vs η²-only=nan.
- A2 incorrect_relative_variance_reduction eta2: Confirmation MAE=0.010228, RMSE=0.015877, R²=0.2528, MAE improvement vs η²-only=0.
- A2 incorrect_relative_variance_reduction eta2_plus_entropy: Confirmation MAE=0.010253, RMSE=0.015977, R²=0.2433, MAE improvement vs η²-only=-2.5619e-05.
- A2 incorrect_relative_variance_reduction all_geometry: Confirmation MAE=0.010031, RMSE=0.015263, R²=0.3094, MAE improvement vs η²-only=0.00019645.

Case frequencies stratified by reconstruction method and RNA context are in CASE_STRATUM_OUTCOME_COUNTS.tsv. Descriptor associations, continuous summaries, relative-variance exclusions, non-finite policy, and mouse-cluster intervals are in the corresponding tables.

## Title verdict

**C. Not supported.** The empirical criteria for common correct-cue variance reduction and replicated Confirmation prediction were not both met; prediction-error improvement remains a separate possible conclusion.

## Limitations

These results describe finite weighted synthetic-truth ensembles, not experimentally validated uncertainty estimates. The crossed bootstrap resamples only five mouse identities per cohort and preserves repeated evaluations/reactions within sampled mouse pairs; it does not capture uncertainty from new reconstruction methods, new cohorts, or the finite candidate-generation process. Prospective models have a single frozen reaction holdout and the all-geometry combination is descriptive absent an independently replicated confirmation set.
