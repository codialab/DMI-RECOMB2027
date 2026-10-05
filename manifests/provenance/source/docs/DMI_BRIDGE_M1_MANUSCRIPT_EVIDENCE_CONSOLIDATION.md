# DMI-BRIDGE-M1 — manuscript evidence consolidation contract

## Purpose

DMI-BRIDGE-M1 is a deterministic, read-only consolidation stage after the qualified
PL2D confirmation. It creates a manuscript-facing evidence ledger and figure-ready
fact tables from already frozen PL2B/PL2C/PL2D evidence.

M1 performs **no new scientific analysis**. It must not fit, tune, select, rank, or
rescue any result. It must not alter any upstream artifact.

Terminal status:

`BRIDGE_M1_MANUSCRIPT_EVIDENCE_FROZEN`

## Frozen source identities

At minimum, M1 must verify these qualified sources before producing any artifact:

- PL2B manifest SHA256:
  `7d9ca97a3f4ea5aa8b0b983e09ee83fc3c82f326d1b1cbf9c6137cfadf25539d`
- PL2B context summary SHA256:
  `d59b0eec6956081a760a52359ca489e461d4b8dbeac784a223c78dc3d20a9a01`
- PL2C manifest SHA256:
  `6abcf149cd6adb294c0931585cec28f992b453451a5b5f03e10b3ffb2f9b9c39`
- PL2D-H manifest SHA256:
  `a1dd76d6b40e741e07efe64c34bc39eaad289a46bf6c870917a4f59cab688427`
- PL2D-H development snapshot SHA256:
  `ebd60bbdd480ec83fe233e068a398776b4bd4ff964cc4bc6bf8c615c66cc28b8`
- PL2D-H support-amendment manifest SHA256:
  `fb7df9e208da2e6c9841d43a3b572038e59943a4221535705eef9f5dc0ed6a3a`
- PL2D manifest SHA256:
  `9f9eb730473e7ad209c1d5ec3340808ab614087b8ae1121d929eff61be12b699`
- PL2D primary-result SHA256:
  `a01ca592b819ab2e6b4463bcbb91ba3bc113bd80dbde2c7aa96b505521bfa1b8`
- PL2D supportive-result SHA256:
  `fbd1e8398a248eeb043e9e985c09dbdcd5d75f067168d8d35b8d18cbec2735b0`
- PL2D context-result SHA256:
  `c063d9809ee0999ce030f7b557fc8064b53348ccda551f3175dd71ef477c744b`
- PL2D QC SHA256:
  `108befc356110ced2ec19ff0fcec6fc01ce38d3d6dd5dc584f2b9f11915936ea`

M1 must also verify all artifact hashes declared by each consumed manifest.

When exact reproduction of rank-based statistics requires re-reading `.17g` TSV
numeric values with pandas, use `float_precision="round_trip"`. Do not recompute
frozen statistics with a different parser and silently replace their qualified values.

## Frozen scientific hierarchy

### Claim C1 — sign-only utility is not universal

At primary lambda = 0.25, the overall PL2B reliability summaries must be reported
without hiding tied or harmed cases.

Required rows:

| q | Improved | Tied | Harmed |
|---|---:|---:|---:|
| 0.00 | 0.1123 | 0.7145 | 0.1732 |
| 0.50 | 0.1541 | 0.7235 | 0.1223 |
| 1.00 | 0.1695 | 0.7160 | 0.1145 |

These displayed fractions are manuscript-rounded values. The adapter must retain
and export the full-precision source values from PL2B and verify that manuscript
rounding reproduces the table above.

The q = 0.50 row is the random-sign control, **not** the strong-anchor baseline.
A positive correct-direction gain alone is not sufficient to claim useful directional
information.

### Claim C2 — geometry predicts whether direction is useful

The central BRIDGE claim is about the **occurrence/fraction of useful directional
information**, not exact recovery of flux and not the magnitude of gain.

Development reaction-level result frozen before confirmation:

- finite reactions: 2,455
- Spearman rho, mean sign-magnitude eta2 vs mean directional-usefulness fraction:
  `0.8570601823118907`
- positive algorithm × RNA context directions: 16/16

Independent reaction-holdout confirmation:

- terminal status: `PL2D_CONFIRMED`
- finite reactions: 838
- Spearman rho: `0.8351435823183292`
- positive contexts: 16/16
- minimum finite context support: 34
- supportive reaction-bootstrap median: `0.8354193419051692`
- supportive 95% computational stability interval:
  `[0.8015629118742916, 0.8631275423376382]`

The bootstrap interval must be described only as computational stability across
reaction benchmark units, not as independent biological replication.

### Claim C3 — utility occurrence is more predictable than gain magnitude

In the confirmation set:

- eta2 vs directional-usefulness fraction rho: `0.8351435823183292`
- eta2 vs mean information advantage rho: `0.34542078458252307`

The manuscript-safe interpretation is that strong-anchor sign-magnitude geometry
predicts **where correct directional information is useful more strongly than how
large the numerical information advantage will be**.

Do not claim that PL1 geometry accurately predicts gain magnitude.

### Claim C4 — confirmation is stable across contexts and pathway sensitivity

Supportive confirmation results:

- partial eta2/usefulness rho controlling mean non-tie coverage:
  `0.8056871576315723`
- directional entropy vs usefulness rho:
  `0.8383836292673364`
- dominant sign mass vs usefulness rho:
  `-0.8257657485170178`
- supported sign-state count vs usefulness rho:
  `0.8207118788803361`
- exclude proximal reactions: `0.8386984413243371`
- exclude same-subsystem-as-HEX1 reactions: `0.8378450478802439`
- exclude either flag: `0.8386984413243371`

These are supportive and must not be presented as additional confirmation gates.

### Limitation L1 — directional truth is unavailable for many cases

The manuscript must retain the truth-tie/coverage limitation.

PL2D confirmation case accounting:

- total confirmation cases: 873,620
- non-tie cases: 359,099
- truth-tie cases: 514,521
- evaluation × reaction rows: 386,650
- rows with defined directional response: 164,690
- rows with no directional truth: 221,960
- mean non-tie-pair fraction: `0.425937` (manuscript-rounded)

Directional usefulness is therefore a conditional benchmark response. Do not imply
that every reaction/context admits a meaningful binary directional truth.

## Frozen manuscript wording boundaries

Allowed central interpretation:

> Strong-anchor flux geometry identifies where additional weak directional
> information is likely to be useful.

Allowed more specific interpretation:

> Sign-magnitude coupling under the HEX1-conditioned ensemble predicts the fraction
> of held-out condition contrasts for which a correct weak directional prior improves
> |Delta B| beyond both the strong-anchor baseline and the random-sign control.

Prohibited overclaims include:

- "the weak prior generally recovers reaction B flux";
- "directional priors improve most reactions";
- "eta2 predicts the exact gain magnitude";
- "the bootstrap is biological replication";
- "the 370 evaluation rows or reaction rows are independent biological replicates";
- any causal claim that geometry causes weak-prior utility;
- any claim that the confirmation threshold was selected after opening the holdout.

## Figure-ready panel contract

M1 must create a figure-data plan with these logical panels. It does not render the
final publication figure.

- **Panel A — BRIDGE concept.** Strong anchor A (HEX1) induces a flux-contrast
  geometry; candidate weak directional information at B is useful only in some
  geometries. No quantitative claim beyond the frozen definitions.
- **Panel B — PL2B utility benchmark.** Primary lambda = 0.25; stacked
  improved/tied/harmed fractions for q = 0, 0.5, 1.0. q = 0.5 must be labelled
  random-sign control.
- **Panel C — development and confirmation.** Show reaction-level eta2/usefulness
  association for DEVELOPMENT and CONFIRMATION as separate datasets. The holdout
  must be visually identified as prespecified confirmation. Do not pool the two
  datasets into a single rho.
- **Panel D — context confirmation.** Show all 16 confirmation context rhos with
  finite reaction counts. No context may be omitted because it is weak.
- **Panel E — specificity/robustness.** Contrast primary usefulness rho with the
  weaker information-advantage rho and show the prespecified pathway sensitivities.

The eta2 × entropy quartile map from PL2C may be retained as a development-only
supplementary visualization, but must not be presented as the confirmation rule.

## Required M1 outputs

Publish under:

`outputs/dmi_bridge_m1_manuscript_evidence_v1/`

Required artifacts:

- `BRIDGEM1_SOURCE_AUDIT.json`
- `BRIDGEM1_CLAIM_LEDGER.tsv`
- `BRIDGEM1_RESULTS_FACTS.tsv`
- `BRIDGEM1_FIGURE_PANEL_PLAN.tsv`
- `BRIDGEM1_PL2B_UTILITY.tsv`
- `BRIDGEM1_DEVELOPMENT_CONFIRMATION.tsv`
- `BRIDGEM1_CONFIRMATION_CONTEXTS.tsv`
- `BRIDGEM1_LIMITATIONS.tsv`
- `BRIDGEM1_MANIFEST.json`

The source audit must pin every consumed upstream file by SHA256. The adapter must
fail closed on source mutation.

The claim ledger must include, at minimum:

- claim ID;
- claim level (`PRIMARY`, `SUPPORTIVE`, `LIMITATION`);
- manuscript-safe statement;
- source stage/artifact;
- exact supporting values;
- interpretation boundary/prohibited overclaim.

## Forbidden work

M1 must not:

- run solver, FVA, sampling, or reconstruction;
- compute new flux weights;
- refit or tune any prior;
- change lambda, reliability, or reaction populations;
- fit a classifier or bridgeability score;
- search a new threshold, subgroup, pathway, or feature;
- rerun PL2C/PL2D with modified rules;
- calculate new inferential p-values;
- pool development and confirmation to report a new primary rho;
- rank reactions for manuscript inclusion;
- omit weak RIPTiDe contexts from the context panel;
- overwrite any PL1/PL2A/PL2B/PL2C/PL2D artifact.

## Completion gate

M1 may finish only if all source identities, exact frozen values, row counts,
manuscript rounding, claim wording boundaries, and deterministic-output checks pass.

Terminal status:

`BRIDGE_M1_MANUSCRIPT_EVIDENCE_FROZEN`

M1 is an evidence-consolidation stage. Final prose editing and final vector figure
rendering are downstream authoring tasks, not reasons to reopen the scientific
analysis.
