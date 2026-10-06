# Supplementary Figure S4 — Held-out robustness

## Scientific objective

S4 shows whether the frozen relationship between residual geometry and the usefulness of weak directional information persists in the held-out reaction population. It keeps development, confirmation, bootstrap stability, and method/context summaries distinct.

## Authoritative inputs

The table builder reads the frozen PL2C development bundle, PL2D-H development snapshot, and PL2D confirmation bundle under `outputs/`. It verifies the PL2C, PL2D-H, and PL2D manifest hashes and every artifact hash listed in those manifests before using data. The build manifest records exact input and generated-table SHA-256 values. The frozen M1 development/confirmation table is a cross-check.

PL2C assigns 3,135 reactions to development; 2,455 have finite reaction-level primary predictor and response values. PL2D assigns 1,045 reactions to confirmation; 838 have finite pairs. The two sets are disjoint. These are reaction-level units, not independent biological replicates.

## Tables and panels

Run the builder from the repository root:

```bash
python figures/supp_fig4/build_supp_fig4_tables.py
```

It creates compact correlation, context, bootstrap, validation, and build-manifest tables in `data/`. The notebook plots:

- **A:** development and held-out confirmation correlations separately, with finite pair counts.
- **B:** the frozen 5,000-replicate reaction bootstrap, observed confirmation correlation, and 2.5th/97.5th percentile interval.
- **C:** every confirmation method × RNA-context correlation, with finite pair counts. The weakest context is retained.

Execute the notebook from the repository root:

```bash
jupyter nbconvert --to notebook --execute --inplace figures/supp_fig4/supp_fig4.ipynb
```

The notebook writes `outputs/supp_fig4.svg` and `outputs/supp_fig4.png` at 300 dpi. The predictor is mean sign-magnitude η² and the response is mean directional-usefulness fraction. S4 uses A1 = 1-Strong-Anchor (glucose uptake).

## Limitations

The confirmation bootstrap resamples reaction rows and describes computational stability. It is not a biological confidence interval. Correlations summarize finite reaction pairs, so denominators vary by population and context. Positive context results do not imply uniform performance; the weakest method/context result is shown explicitly.
