# DMI-BRIDGE-PL2D-H — support-gate amendment v1

## Purpose

This document is an **append-only amendment** to the qualified PL2D-H hypothesis freeze. It corrects one support-adequacy gate before any PL2D confirmation predictor values or outcomes are opened.

The original PL2D-H bundle remains immutable and auditable. The later PL2D confirmation must verify both the original PL2D-H manifest and the canonical amendment manifest. The amendment does not replace, rewrite, or silently republish the original freeze.

## Why an amendment is required

The original freeze required all 16 algorithm × RNA contexts to contain at least 100 finite confirmation reaction pairs. After publication of the freeze, review of the **already-open DEVELOPMENT snapshot only** showed that the four RIPTiDe contexts contain just 129, 159, 134, and 162 finite reactions among 3,135 development reactions.

The confirmation split contains 1,045 reactions, exactly one third of the 3,135-reaction development set. Proportional support expected from the sparsest development context is therefore:

`129 × 1045 / 3135 = 43` finite confirmation reactions.

Thus the original 100-per-context requirement exceeds even the development-calibrated proportional support of the sparsest context by more than two-fold and could mechanically force `PL2D_INSUFFICIENT_SUPPORT` even if the association direction generalizes.

This defect was identified while:

- `pl2d_confirmation_status = NOT_RUN`;
- confirmation PL1 predictor values had not been accessed;
- confirmation PL2B outcomes had not been accessed; and
- no confirmation geometry-to-utility association had been computed.

Therefore the support gate can still be amended without conditioning on holdout results.

## Frozen amendment

Exactly one field changes:

- original `minimum_finite_reactions_per_context = 100`;
- amended `minimum_finite_reactions_per_context = 30`.

The value 30 is frozen from development support only. The sparsest context has 129 finite development reactions, implying 43 under exact 1:3 proportional scaling to the holdout. A prespecified 70% retention margin gives `floor(43 × 0.70) = 30`.

The value is a **support-adequacy floor**, not an effect threshold and not a power claim. It ensures that context-direction auditing is not declared evaluable from only a handful of reactions while avoiding a gate that is incompatible with the support density already observed in development.

## Unchanged primary confirmation contract

Every other primary and supportive definition remains exactly as frozen in the qualified PL2D-H bundle:

- primary predictor: reaction-level mean `sign_magnitude_eta2`;
- primary response: reaction-level mean `directionally_useful_fraction`;
- primary unit: held-out reaction;
- overall minimum finite reactions: 700;
- all 16 algorithm × RNA contexts must be evaluable;
- overall Spearman rho threshold: `+0.50`;
- at least 12 of 16 context associations must be positive;
- numerical reference: pandas grouped means plus `scipy.stats.spearmanr`;
- 5,000-replicate reaction bootstrap remains supportive and nongating;
- all previously frozen supportive analyses remain supportive only.

The primary status rule after amendment is:

- `PL2D_INSUFFICIENT_SUPPORT`: fewer than 700 finite overall reaction pairs, or any of the 16 contexts has fewer than 30 finite reaction pairs, or the overall rho is non-finite;
- `PL2D_NOT_CONFIRMED`: support is adequate but overall rho is below `+0.50` or fewer than 12 of 16 contexts have positive rho;
- `PL2D_CONFIRMED`: all amended support requirements and all original effect requirements pass.

Supportive analyses cannot rescue a failed primary result.

## Frozen development support used for calibration

The amendment may use only the finite reaction counts already published in the canonical PL2D-H DEVELOPMENT snapshot:

| Algorithm | RNA context | Finite development reactions |
| --- | --- | ---: |
| CORDA | setx1,setx2 | 1847 |
| CORDA | setx1,setx2,setx3 | 1815 |
| CORDA | setx1,setx3 | 1854 |
| CORDA | setx2,setx3 | 1841 |
| GIMME | setx1,setx2 | 1235 |
| GIMME | setx1,setx2,setx3 | 1264 |
| GIMME | setx1,setx3 | 1237 |
| GIMME | setx2,setx3 | 1283 |
| RIPTiDe | setx1,setx2 | 129 |
| RIPTiDe | setx1,setx2,setx3 | 159 |
| RIPTiDe | setx1,setx3 | 134 |
| RIPTiDe | setx2,setx3 | 162 |
| iMAT | setx1,setx2 | 1729 |
| iMAT | setx1,setx2,setx3 | 1737 |
| iMAT | setx1,setx3 | 1701 |
| iMAT | setx2,setx3 | 1726 |

No confirmation support count may be inspected to choose or revise the value 30.

## Provenance and publication

The amendment publication must:

1. verify the original PL2D-H manifest SHA256 `a1dd76d6b40e741e07efe64c34bc39eaad289a46bf6c870917a4f59cab688427`;
2. verify all five original PL2D-H artifact hashes;
3. verify terminal status `PL2D_HYPOTHESIS_FROZEN_READY_FOR_CONFIRMATION` and `pl2d_confirmation_status = NOT_RUN`;
4. verify all original confirmation-access firewall flags are false;
5. reproduce the 16 frozen development-context support counts above from `BRIDGEPL2DH_DEVELOPMENT_SNAPSHOT.json`;
6. derive the amended floor of 30 using only the frozen development counts and split sizes;
7. prove that every original primary/supportive contract field is unchanged except the per-context support floor and the corresponding status-rule text; and
8. publish a separate amendment bundle atomically.

The original PL2D-H output directory must not be modified.

## Information firewall

The amendment stage must not read or derive for confirmation reactions:

- PL1 predictor values;
- PL2B outcomes, truth ties, gains, or utility responses;
- finite confirmation support counts;
- any confirmation geometry-to-utility association.

It may read only the qualified original PL2D-H artifacts needed to verify and amend the contract. It must not open PL1 or PL2B sources.

## Interpretation

This amendment changes **support adequacy only**. It does not alter the hypothesis, predictor, endpoint, effect-size threshold, direction criterion, feature set, aggregation, lambda, or confirmation population. It was frozen before holdout opening and therefore must not be characterized as a response to PL2D results.
