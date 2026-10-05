# DMI-BRIDGE-M2 — A2 robustness manuscript-evidence consolidation

## Purpose

DMI-BRIDGE-M2 is a deterministic, read-only evidence-consolidation stage after
the frozen A2.3 matched A1/A2 synthesis.

M2 is an additive successor to the already frozen M1 manuscript-evidence bundle.
It must not modify, replace, reinterpret, or republish M1 as if M1 had originally
contained A2. M1 remains the frozen one-strong-anchor evidence bundle.

M2 adds a compact manuscript-facing ledger for the post-freeze A2-L robustness
branch:

- A2.0: dual-anchor qualification;
- A2.1: residual geometry under `Vmax + Vlac -> HEX1 + LDH_L`;
- A2.2: matched sign-only weak-prior utility under A2-L;
- A2.3: matched A1-versus-A2 synthesis.

M2 performs **no new scientific analysis**. It may verify hashes, copy frozen
facts, apply manuscript rounding, assemble deterministic tables, and state
interpretation boundaries. It must not refit, resplit, rerank, retest, or search.

Preferred terminal status:

`BRIDGE_M2_A2_ROBUSTNESS_EVIDENCE_FROZEN`

## Required frozen source chain

M2 must fail closed unless it verifies the exact live frozen identities of:

- M1: `BRIDGE_M1_MANUSCRIPT_EVIDENCE_FROZEN`;
- A2.0: `BRIDGE_A20_PARTIAL_QUALIFICATION`;
- A2.1: `BRIDGE_A21_DUAL_ANCHOR_GEOMETRY_FROZEN`;
- A2.2: `BRIDGE_A22_SIGN_ONLY_UTILITY_FROZEN`;
- A2.3: `BRIDGE_A23_A1_A2_SYNTHESIS_FROZEN`.

A2-G must remain `NOT_QUALIFIED` and unused. M2 must not turn A2-G into a
manuscript result.

The repository adapter must verify all artifact hashes declared by each consumed
manifest and must bind literal paths plus SHA256 values in the M2 source audit.

The current frozen A2.3 snapshot expected by this foundation has:

- manifest SHA256:
  `b8e57074f3df23fe7708145bc01f51a19f946545b1c16e61628c2678dbcabd18`
- QC SHA256:
  `ed502de0e91e793d4c0ed3c7a0237355ae30a0bf69d81ccde2e1a818e358ee20`
- utility comparison SHA256:
  `a59b1870d50efc87c75d3be32c5e890f0b63413e38c87a47adb65dc99f00f354`
- geometry comparison SHA256:
  `e9e0561248a231d5f9480365888a76cd8c99f6d06dd04ab4d92141601fa21f2f`
- geometry-utility persistence SHA256:
  `33af7c3b59a32872aee0070af10f35b72f974537bb6e1745dcc2c5fbd42afc03`
- control summary SHA256:
  `e70e9361e68cfd3e8d72eb27095a0781eb564626913951b3bbe5631d83f45cb9`
- sensitivity summary SHA256:
  `2d33a5358ad38741660585ffbb52f68636907e825467cedaab46ad7df5a6497f`

If the live A2.3 manifest differs, do not silently repin. Establish whether the
difference is a qualified serialization/provenance successor or a scientific
change and stop if identity cannot be reconciled.

## Frozen matched population

The primary A2 robustness evidence is based on the frozen matched population:

- truth pairs: 836;
- reactions: 4,179;
- matched case keys: 3,493,644;
- `LDH_L` excluded from the A1 weak-target universe because it is a strong
  anchor in A2-L;
- 49 original A1-non-evaluable truth pairs remain excluded from primary outcome
  production.

No population may be changed for manuscript convenience.

## Frozen A2 robustness facts

### R1 — matched correct-direction utility is larger under A2-L

Across all matched non-tie case keys:

- A1 mean correct-direction gain:
  `0.014234484268356776`
- A2-L mean correct-direction gain:
  `0.022801865806147577`
- paired A2-minus-A1 difference:
  `0.0085673815377913642`

The corresponding correct-minus-wrong directional-specificity means are:

- A1: `0.03309389154965113`
- A2-L: `0.054535148930311576`
- paired A2-minus-A1 difference:
  `0.021441257380664738`

This is a descriptive matched post-freeze comparison. Do not call it an
independent confirmation or a newly prespecified hypothesis test.

### R2 — the effect remains direction-specific

Across the same matched non-tie cases:

- wrong-direction mean gain:
  - A1: `-0.018859407281294222`
  - A2-L: `-0.031733283124166602`
- analytic random-sign (`q=0.5`) mean gain:
  - A1: `-0.0023124615064686166`
  - A2-L: `-0.0044657086590093476`

The manuscript-safe interpretation is that the increased correct-direction gain
is accompanied by stronger penalty for wrong direction; the A2 result is not a
generic regularization benefit.

`q=0.5` remains a random-sign control carrying no directional information on
average.

### R3 — the frozen geometry-to-utility relationship persists

Using the exact frozen PL2C/PL2D estimator on the matched target universe:

Development split:

- matched A1 rho: `0.8571160663367106`
- A2-L rho: `0.8339440056039364`
- finite reactions: A1 `2454`, A2 `2485`

Confirmation split:

- A1 rho: `0.8351435823183291`
- A2-L rho: `0.8493657186460156`
- finite reactions: A1 `838`, A2 `843`

The A2-L confirmation bootstrap is supportive computational stability only:

- replicates: 5,000;
- median rho: `0.8488457925736193`;
- 2.5th percentile: `0.823810132333918`;
- 97.5th percentile: `0.8714723352450245`.

Do not claim that the small A1/A2 rho differences are statistically significant
or that A2 improves the predictor unless a separately authorized test exists.

### R4 — residual geometry changes in mixed ways

On the exact 1,671,600 paired evaluation × reaction geometry keys:

- directional entropy mean increases under A2-L;
- dominant-direction mass mean decreases under A2-L;
- joint ESS mean increases under A2-L;
- the width-contraction descriptor changes substantially;
- finite support for `sign_magnitude_eta2` remains incomplete.

For `sign_magnitude_eta2`:

- paired candidate rows: `1,671,600`;
- finite paired rows: `573,882`.

M2 must export exact source values for all eight A1-comparable geometry
descriptors rather than summarize A2-L as universally "more constrained" or
"more ambiguous."

## Required limitations

### L1 — A2 is post-freeze robustness, not independent confirmation

PL2D confirmation was opened before A2. A2.0–A2.3 therefore provide matched
robustness/generalization evidence only.

### L2 — only the lactate second-anchor arm qualified

A2-L (`HEX1 + LDH_L`) qualified. A2-G did not pass the quantitative-operator
provenance gate. Do not generalize the A2 result to arbitrary second strong
anchors.

### L3 — directional truth remains conditional

Across all matched A2.3 cases:

- all case keys: `3,493,644`;
- non-tie pairs: `1,440,313`;
- truth ties: `2,053,331`.

Gain summaries for correct/wrong/random direction are therefore conditional on
non-tie directional truth.

### L4 — finite geometry support remains incomplete

Do not interpret missing/undefined eta2 as zero. M2 must retain exact finite
denominators.

## Manuscript wording boundaries

Allowed:

> The geometry-to-utility relationship remained strong after adding a second
> qualified quantitative anchor.

Allowed:

> In the matched post-freeze robustness analysis, correct directional
> information produced a larger mean gain under the dual-anchor A2-L baseline,
> while wrong direction became more harmful.

Allowed:

> These results indicate that weak directional information can remain useful
> even when an additional quantitative metabolic constraint is available.

Prohibited:

- "A2 independently confirms the BRIDGE hypothesis";
- "two strong anchors always make the weak prior more useful";
- "adding LDH_L causes better prediction";
- "the A2 confirmation rho is significantly higher than A1";
- "the random-sign control improves inference";
- generalization from A2-L to A2-G or arbitrary second anchors;
- describing finite computational reaction units as biological replicates.

## Figure-ready evidence contract

M2 must create a compact figure-data plan but must not render or redesign the
already frozen Figure 1.

Recommended A2 robustness logical panels:

- **A — matched design.** A1: HEX1 strong anchor; A2-L: HEX1 + LDH_L strong
  anchors; identical 836 truth pairs and 4,179 weak-target reactions.
- **B — utility and specificity.** A1 versus A2-L correct, wrong, and analytic
  random-sign mean gains, plus correct-minus-wrong specificity.
- **C — geometry-to-utility persistence.** Matched development and confirmation
  rho values for A1 and A2-L, shown without significance stars.
- **D — residual geometry descriptors.** Paired shifts in the eight original
  A1-comparable geometry features, with finite denominators explicit.
- **E — limitations.** A2-L only; A2-G provenance failure; post-freeze
  robustness status; truth-tie and finite-support limitations.

The eventual manuscript may place these panels in a main or supplementary
figure. M2 must not decide placement by outcome favorability.

## Required M2 outputs

Publish only under:

`outputs/dmi_bridge_m2_a2_robustness_evidence_v1/`

Required artifacts:

- `BRIDGEM2_SOURCE_AUDIT.json`
- `BRIDGEM2_CLAIM_LEDGER.tsv`
- `BRIDGEM2_RESULTS_FACTS.tsv`
- `BRIDGEM2_A1_A2_UTILITY.tsv`
- `BRIDGEM2_A1_A2_GEOMETRY.tsv`
- `BRIDGEM2_A1_A2_PERSISTENCE.tsv`
- `BRIDGEM2_CONTROLS.tsv`
- `BRIDGEM2_LIMITATIONS.tsv`
- `BRIDGEM2_FIGURE_PANEL_PLAN.tsv`
- `BRIDGEM2_MANIFEST.json`

All rows must be deterministic and source-traceable.

## Forbidden work

M2 must not:

- run solver, FVA, reconstruction, or sampling;
- generate or refit strong-anchor weights;
- rerun A2.1, A2.2, or A2.3;
- select new truth pairs or reactions;
- change lambda, q, tie tolerance, gain tolerance, or the dev/confirmation split;
- compute new hypothesis tests or p-values;
- search new subgroups, thresholds, predictors, or pathways;
- rank reactions for manuscript inclusion;
- modify M1 or any predecessor;
- choose manuscript claims based on whether A2 appears favorable.

## Completion gate

M2 may freeze only if every consumed manifest/artifact hash verifies, all frozen
A2.3 facts reproduce exactly from the source tables, manuscript rounding is
deterministic, wording boundaries are encoded in the claim ledger, and an
identical completed rerun is a no-op.

Terminal status:

`BRIDGE_M2_A2_ROBUSTNESS_EVIDENCE_FROZEN`
