# Anchor operator audit: S1

Audit date: 2026-10-06. Source commit:
`16379a169f45ec556274953eb66ae12293e1e985`. Starting destination commit:
`ae0ce8e9885c70f686708e1ee5ac35e3fa231ee6` (`publication-migration`).
Inspection was read-only; no solver, reconstruction, sampling or operator rerun
was performed. No rights-pending data were copied.

## Verdict

**SCIENTIFIC DISCREPANCY REQUIRES AUTHOR DECISION.** The implemented quantitative
weighting coordinate for glucose is `max(-EX_glc__D_e, 0)`, not HEX1. HEX1 is
used downstream for truth selection, contrast-coupling descriptors and exclusions.
The inspected lineage establishes neither an implemented Vmax-to-HEX1 conversion
nor equality or rank equivalence between glucose exchange and HEX1. Their biological
relationship does not establish mathematical equivalence.

Calling A1 a direct quantitative constraint on HEX1, or A2 direct constraints on
HEX1 + LDH_L, therefore does not accurately describe the implemented weighting
operator. The code does implement one- and two-observable soft conditioning;
neither imposes new mouse-specific flux bounds in the Bridge analysis. An author
must decide whether to describe the actual operator or pursue a separately
authorized scientific change. This audit changes no scientific definition or code.

## Production call chain and exact locations

Paths below are relative to the source root. Migrated `scripts/` and Stage-11
modules have identical source bytes. Other cited upstream files remain source-only.

| Stage | Variables/reactions and action | Code locations |
|---|---|---|
| Measurement provenance | Established loader reads fitted denoised DMI maps and subject metadata. `tumor_filtered` requires finite/nonnegative Vmax/Vlac/Vglx, Vmax ≤ 2.5, and branch agreement within `max(0.05, 0.05*Vmax)`. | `scripts/load_brain_glioma_GL261_CT2A_example.py:167–252` |
| Mouse observable | Median across retained tumor voxels, separately for Vmax, Vlac and Vglx; identifiers derive from the subject manifest. | `pre-simulation_constraint_analysis/analysis_a_capacity.py:557–584`, `mouse_dmi_table` |
| Frozen mouse table | Corrected Analysis A reads its configured existing mouse table and writes the canonical `mouse_level_dmi_input.csv`; Bridge consumes this frozen table, not a fresh fit. Historical capacity constraints are a separate analysis. | `pre-simulation_constraint_analysis/run_corrected_analysis_a.py:374–375` |
| Candidate generation | Qualified Stage-11 jobs are `rna_only_primary`, `condition_type=rna_only`, with no DMI mouse ID. Configuration disables RNA/DMI hard constraints and anchor protection. Primary requirement policy has an empty measurement-required list. | `11_recomb_flux_expert_generation/provenance/master_config.yaml:11–14`; `src/models.py:94–131`; `src/jobs.py:293–331` |
| Candidate panel | `collect_candidates` reads raw full-parent vector columns including EX_glc__D_e and LDH_L and attaches ensemble, projection, RNA and sample identities. The separate DMI-fraction ensemble is not a fifth production Bridge method. | `12_recomb_method_comparison/scripts/run_dmi_fraction_ensemble.py:260–305` |
| Raw candidate observables | Vmax proxy is `max(-EX_glc__D_e,0)`; lactate proxy is `max(-LDH_L,0)`. Minus sign selects uptake/production orientation; clipping retains nonnegative proxy magnitudes. Tumor-wide candidate midranks are formed. | `outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/build_current_candidates.py:6–12` |
| Original current A1 weights | Per-tumor mouse Vmax ranks are compared with the candidate glucose-uptake ranks. Squared distance enters the ESS-20 kernel. Saved `Vmax_coordinate` is the candidate rank; misleading historical `q_Vmax` in the weight table is the saved weight. | `13_recomb_identifiability/scripts/run_sa_u1_u2_followup.py:285–296,334–337` |
| A1 qualification/preservation | BIO0R independently reproduces the exchange proxy and rank kernel; BIO0 copies the saved baseline into `strong_anchor_baseline` without a HEX1 remapping. PL1 loads that arm. | `scripts/dmi_bridge_bio0r_repo_qualification_v1.py:271–330`; `dmi_bridge_bio0_repo_audit_v1.py:89–123`; `dmi_bridge_pl1_predictability_landscape_v1.py:169–216` |
| A20 qualification | `load_panel` joins cached HEX1/LDH_L/AKGDm values, but computes candidate Vmax from the raw panel exchange column. Lactate weights use raw panel LDH_L, not carbon-normalized cached LDH_L. Candidate ranks cover 320 vectors/tumor; mouse ranks cover five subjects/tumor. | `scripts/dmi_bridge_a20_dual_anchor_qualification_v1.py:119–170` |
| A2 weight production | `qualify` reproduces A1 using `vdist=(q_Vmax-mouse_q_Vmax)^2`; A2 adds `(q_Vlac-mouse_q_Vlac)^2`. HEX1 does not appear in either distance. Both calibrate temperature over all 320 tumor-matched candidates to ESS 20. | Same file `:235–277`; `scripts/dmi_bridge_a20_dual_anchor_core_v1.py:43–120` |
| Conditional geometry | Within each method/RNA/tumor/mouse stratum, restrict the global weights to its 20 candidates and normalize. Contrast weights are products on the CT2A–GL261 Cartesian support. Cached HEX1 supplies the common diagnostic contrast; its value does not generate the strong weights. | `scripts/dmi_bridge_pl1_predictability_resume_v1.py:270–299`; `dmi_bridge_pl1_predictability_core_v1.py:111–125,199–219`; `dmi_bridge_a21_dual_anchor_geometry_v1.py:148–163` |
| Truth and holdout | A1 sorts the cached HEX1 coordinate, uses A1 weights to select q10/q50/q90 whole vectors, and removes selected vectors before inference. A2 reuses the same truth identities, substitutes A20 weights, removes each held-out vector and renormalizes. No additional HEX1 measurement matching occurs. | `scripts/dmi_bridge_pl2a_prepare_v1.py:224–313`; `dmi_bridge_a22_dual_anchor_sign_only_utility_v1.py:114–181,217–238` |
| Matched endpoints | A23 joins the frozen PL2/A21/A22 outputs. HEX1 coupling is preserved as the common named descriptor; the matching stage does not reinterpret the glucose weighting coordinate. | `scripts/dmi_bridge_a23_a1_a2_synthesis_v1.py:167,248`; A23 frozen contract |

## Operator definition and interpretation

For candidate `i`, let `g_i=max(-v_EX_glc__D_e,i,0)` and
`l_i=max(-v_LDH_L,i,0)` from the raw panel. Their average-rank coordinates are
`q_g,i=(rank(g_i)-0.5)/320` and `q_l,i=(rank(l_i)-0.5)/320` within tumor.
Measured mouse coordinates are `(rank(Vmax_m)-1)/4` and `(rank(Vlac_m)-1)/4`.

The implemented distances are:

```text
A1: D_i = (q_g,i - q_Vmax,m)^2
A2: D_i = (q_g,i - q_Vmax,m)^2 + (q_l,i - q_Vlac,m)^2
w_i ∝ exp(-(D_i - min(D))/temperature_m), normalized over 320 candidates
```

Temperature is calibrated separately for each mouse and setting to global ESS 20.
Restriction/renormalization to 20 candidates is a later step; conditional ESS is
not fixed at 20. Thus the exchange observable is both a ranking input and the
coordinate that actually determines glucose soft conditioning. It is not an
unused measurement label. Changing it to HEX1 could change ranks, temperatures,
weights, truth selection and downstream endpoints; no such change is authorized.

The descriptor cache uses total-carbon-normalized vectors, whereas the quantitative
proxies originate in raw full-parent vectors. Stage-12 normalization is implemented
in `run_dmi_fraction_ensemble.py:679–734`; cache reconstruction identifies these
normalized matrices in `reconstruct_flux_cache.py:263–291`. This further prevents
identifying the raw exchange weighting coordinate with the cached HEX1 diagnostic
merely by comparing names. No candidate-flux equivalence calculation was run here.

The canonical mouse table matches the exact A20 pin. Its configured predecessor
has identical subject/QC metadata but is not byte-identical and has very small
rate differences (recorded as aggregate maxima in the evidence index). These are
compatible with CSV read/write precision; causation was not established by rerunning
the producer. The canonical pinned table remains authoritative and was not replaced
with the earlier table. Neither table implements a Vmax-to-HEX1 mapping.

## Mapping evidence and conflicts

`analysis/brain_glioma_gem/config/dmi_to_gem_mapping.tsv:2–3` explicitly distinguishes:

- `vmax_glucose_exchange`: EX_glc__D_e, active by default; scale not validated.
- `vmax_hexokinase`: HEX1, inactive alternative; other fates/transport may differ.

These are declarations of alternative mappings, not an executed exchange-to-HEX1
transformation. A20 reads this mapping file at lines 285–286 only to verify the
inactive Vglx/AKGDm association. It never uses the HEX1 mapping row for glucose
weights. The BIO0 registry names EX_glc__D_e explicitly for `Vmax_strong`.

This explains why `load_panel` uses negative exchange uptake: it deliberately
reproduces the pre-existing, qualified Vmax proxy and its candidate ranks, with
hash-gated regression checks. The earlier mapping favors an interface-level glucose
observable and marks intracellular HEX1 as an alternative. That provenance explains
the implemented choice; it does not supply a scientific equivalence proof.

The frozen A20 anchor registry simultaneously records
`reaction_associations=HEX1;LDH_L` and
`candidate_observables=max(-EX_glc__D_e,0);max(-LDH_L,0)`.
The association label is not proof of coordinate equivalence. Contracts describe
Vmax→HEX1, but qualification verifies reproduction of the archived exchange-based
operator, not a separately established HEX1 operator.

Earlier capacity/FBA analyses contain aggregate intracellular glucose ceilings,
including multiple reactions, and historical HEX1-containing hard constraints.
They must not be imported as proof of this Bridge operator: the final candidate
jobs are RNA-only and Bridge reweights frozen vectors without changing bounds.

The preserved `RECOMB_MANUSCRIPT_TECHNICAL_AUDIT.md` section E2 asserts association
with a HEX1 *candidate coordinate*. That interpretation is not supported by the
traced weighting calculation. This audit supersedes E2's operator interpretation;
the original audit is preserved as provenance. Its distinction between rank
weights and absolute unit conversion remains supported.

## Frozen evidence

All exact file identities are in
[`S1_S2_EVIDENCE.json`](../manifests/provenance/S1_S2_EVIDENCE.json), alongside
[`FROZEN_LINEAGE.json`](../manifests/provenance/FROZEN_LINEAGE.json).

| Artifact | Production role |
|---|---|
| Canonical corrected-Analysis-A mouse input, SHA-256 `ae3cbad6d422405780d0dc44e0e31f95770ebc204d1f80f2ccdada12bb119082` | Frozen measured-observable source pinned by A20 |
| `candidate_fraction_table.tsv.xz`, SHA-256 `421c2c1d196fbd088a4e92e2ecc998ad578aab814a025fff8e1c843597fed1da` | Raw candidate proxy and stable identity lineage |
| `current_vmax_base_weights.tsv.xz`, SHA-256 `91121bf1251354505cde5f5ae783dc8d14bc9a22ecdf8daa06be1e5d82e3385d` | A1 exchange-rank baseline; manifest names archived corrected operator |
| `BRIDGEBIO0_FOUR_ARM_WEIGHTS.tsv.xz`, SHA-256 `843e17b7012a6fa212d1e9ac9e1dc6e309e854d94ab5c4b3f9e316b589768f4f` | PL1/PL2A select only `strong_anchor_baseline` |
| BIO0R manifest `d155b63f59db3479b1dfd2d87c0b9ca470358c0166496abf288a649f04890cea` | Qualified raw/rank/weight reproduction, not an equivalence proof |
| A20 registry `18ae24892e4989fac21dbccd33a53ce78b8aaaaebe72ab01b198a86e3ab3da69` | Separates reaction association from implemented observable |
| A20 dual weights `b153e2fde45dca9573c71c548783e704313f056d0980f64af60180ad90ebc1f4` | Actual A2 weights; only historical A2-L arm qualified |
| Cache `f9644da4ebd875cf5f0b34cc58b20faf0a0457d7c7a44635e8a010d8896c3dbb` | 32×20×4181 normalized diagnostic vectors; contains HEX1 |
| PL2A, A21, A22, A23 frozen manifests/contracts | HEX1 truth/coupling roles, held-out vector logic and final matched population |

## Manuscript implication

Supplementary Methods S2's statement that the one-anchor distance uses HEX1 alone
is contradicted by the production distance. Main Methods, captions and shorthand
definitions require author review; the discrepancy cannot be resolved merely by
noting that the weights are dimensionless. Safe factual wording for review is:
“A1 matches measured Vmax ranks to candidate glucose-uptake ranks; A2 additionally
matches measured Vlac ranks to the lactate-producing LDH_L proxy. HEX1 defines the
common truth-selection and coupling coordinate.” Acceptance of this interpretation
is an author decision. Numerical results are unchanged and not recalculated.

Affected locations and proposed corrections are in
[`AUTHOR_CORRECTIONS_REQUIRED.md`](AUTHOR_CORRECTIONS_REQUIRED.md).
