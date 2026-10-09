# Author corrections required after S1/S2 investigation

Date: 2026-10-06. Source commit `16379a169f45ec556274953eb66ae12293e1e985`.
Findings are supported by [the anchor audit](ANCHOR_OPERATOR_AUDIT.md),
[the evaluation audit](EVALUATION_DESIGN_AUDIT.md) and
[the evidence index](../manifests/provenance/S1_S2_EVIDENCE.json).
Scientific definitions, plotting labels, original contracts and manuscript DOCX
files were not changed. Only unambiguous S2 documentation was clarified.

| Issue | Authoritative finding | Affected locations | Proposed correction | Scientific impact |
|---|---|---|---|---|
| S1: literal anchor definition | A1 distance uses raw `max(-EX_glc__D_e,0)` ranks; A2 adds raw `max(-LDH_L,0)` ranks. No direct HEX1-distance term or conversion exists in the traced weighting stage. | README Analysis settings; both `analyses/*/README.md` openings; main draft abstract/Methods/results/captions; Supplementary Methods S2 and S8 | Author to decide whether to describe A1 as Vmax-to-glucose-uptake soft conditioning and A2 as additionally Vlac-to-LDH_L conditioning. Retain A1/A2 identifiers; do not imply equivalence with HEX1. | Interpretation correction is required for claims of literal HEX1 weighting. Substitution of HEX1 would change the scientific operator and require separate authorization/validation. |
| S1: HEX1 diagnostic roles | HEX1 genuinely orders deterministic truth selections and defines common contrast-coupling/subsystem diagnostics under exchange-based weights. | PL2 truth-selection descriptions, eight-descriptor panel, exclusion labels and reaction counts; Supplementary Methods S3/S4/S7, Table S5 and Figure S6 | Preserve these exact HEX1 roles; explicitly distinguish them from the coordinate generating the glucose weights. Do not globally replace HEX1 with EX_glc__D_e. | Truth identities, reaction exclusions, 4,179 population and coupling endpoints remain unchanged. Whether the weak-cue population should exclude the weighting exchange coordinate is a separate author question; no reaction-set change is made here. |
| S1: older authority overstatement | A20 registry separates HEX1 association from exchange observable; the preserved technical audit E2 overstates a HEX1 candidate-coordinate match. | `manifests/provenance/source/RECOMB_MANUSCRIPT_TECHNICAL_AUDIT.md` E2 and its concluding no-conflict statement; frozen PL1/PL2/A20–A23/M1/M2 contract prose | Add an author-approved erratum referencing the executed observable. Preserve historical contracts/audit bytes and hashes. | Qualification proves reproduction of the archived operator, not biological or mathematical equivalence. |
| S1: figure labels | Existing display legends call A1 HEX1 and A2 HEX1 + LDH_L. | Figure 2 D–G titles; Figure 3 B/C legends and D caption; Figure 4 A schematic and B/C legends; S1 titles/manifest; provenance-only supplementary S4–S7 labels | After author decision, replace literal reaction-anchor legends with A1/A2 and define the actual weighting observables in captions. Retain explicit HEX1 labels on genuine HEX1 diagnostic axes. | Display/interpretation change only if frozen data/operator are retained; no current script edit or rendering performed. |
| S2: 400 evaluations | Exact full product: 4 methods × 4 RNA contexts × 5 CT2A mice × 5 GL261 mice = 400. The two fives are DMI subjects. | README design paragraph, publication provenance, migration issue S2; main Methods paragraph 28; prospective design/overview captions and Table S1 | Use “25 cross-tumor mouse comparisons per method/RNA context; 400 computational evaluations sharing ten subjects and fixed vectors.” Remove “two mice × five draws” wherever used. Current main draft paragraph 28 is already correct. | Terminology only; no count, grouping or analysis change. README/provenance edits applied in this commit. |
| S2: candidates versus contrasts | 32 RNA-only ensembles × 20 vectors = 640. All frozen jobs use ensemble replicate 0. Each evaluation has 20×20=400 pre-holdout contrasts; held-out utility has 19×19=361. | Methods candidate-panel paragraph; Supplementary S1/S2/S4 and Table S1; S1 product-ESS title; Figure 4/S7 density-support descriptions | Distinguish 400 evaluations from 400 candidate contrasts/evaluation and 836 matched truths. Do not interpret candidate states or crosses as independent animals or five reconstruction repeats. | Terminology/support interpretation only; existing counts preserved. |
| S2: actual randomness | Stored vectors have seeded upstream sampler provenance; weighted quantile truth selections are deterministic; q=0.5 random-sign control is analytic; reaction bootstrap is separate. | Supplementary sampler placeholders, truth “sampling” language, control/uncertainty descriptions; provenance contract terminology | Say “deterministic weighted quantile truth selection” and name upstream sampler/seed and bootstrap separately. Do not invent a common generated-draw count or a new evaluation axis. | No random seeds or sampling logic change. |

## Manuscript snapshot locators

Paragraph numbers are one-based OOXML `w:p` document-order positions, including
table-cell paragraphs, not Word page or displayed paragraph numbers. The evidence
index records all 54 paragraphs matching anchor/design terms without copying the
manuscript text. These are review locators: many diagnostic uses of HEX1 and most
current design statements are correct and should be retained.

- Main draft (`RECOMB_with_A2_Manuscript_Draft.docx`, SHA-256
  `6b837c1229b9f495aa1dcfab4b1512c7ef4c1df85e3aef3ab07b174ac68e467b`):
  S1 claims p15, p30, p57–58, p60, p64–65, p71 and p73; retain genuine HEX1 truth
  selection p38 and reaction registry/exclusion p48. Design/support review p15,
  p28, p30, p33, p52, p68, p71 and p79. **p28 already states the correct 400 design.**
- Supplement (`RECOMB_with_A2_Supplementary_Materials_Draft.docx`, SHA-256
  `ae62ba34b22d6bf1d81a7c39535430abb091940ba8063e4b0aa51471143c5fcc`):
  S1 weighting/mapping claims p17, p22, p177, p191 and S1 caption p369;
  p22 explicitly says the one-anchor distance uses HEX1 alone and requires
  correction. Preserve HEX1 descriptor/truth/exclusion roles at p45, p47, p50–51,
  p90, p93, p119, p181, p189, p301, p356 and p384, with a distinction from weighting.
  DMI observable provenance p10/p16 requires no new coordinate substitution.
  S2 design/support locators: p7, p13–14, p16, p20, p27, p30–31, p34–35, p51,
  p63, p69, p71–72, p78, p111, p119, p177, p314, p362, p365, p369, p384 and p387.
  The final example's C1/G1 identity is a particular mouse comparison, not an
  assertion that only two mice were used in the study.

## Plot-label locations requiring S1 review

| Migrated source | Locations |
|---|---|
| `figures/main/fig2/plot_fig2_dg.py` | `panel_titles`, lines 23–28 use generic one-/two-anchor counts; these titles can stay, but their caption must define the actual observables |
| `figures/main/fig3/plot_fig3_bc.py` | `display_labels`, line 20 |
| `figures/main/fig4/plot_fig4_b.py` | Legends, lines 108 and 116 |
| `figures/main/fig4/plot_fig4_c_source.py` | Legends, lines 116–117 |
| `figures/supplementary/s1/plot_s1.py` | `anchor_titles`, line 30; product-support label line 107 correctly denotes 400 pairs/evaluation |
| `figures/supplementary/s1/upstream/build_supp_fig1_tables.py` | Manifest anchor labels, line 249 |
| `manifests/provenance/legacy_figure_code/build_supp_fig5_tables.py` | Display labels, lines 238 and 253 |
| `manifests/provenance/legacy_figure_code/build_supp_fig7_tables.py` | ARMS line 33 and manifest labels lines 231–232 |

The Figure 1 four-panel SVG has since been adopted for the current manuscript,
with its PNG export and provenance recorded in the publication package. The
separate Figure 2 anchor schematics, Figure 3A, and Figure 4A artwork remain
subject to their own review. Their source history is documented in
`ANALYSIS_MAP.md`. Generic A1/A2 identifiers and genuine HEX1 coupling labels
require no automatic replacement.

## Documentation changes applied

README and `docs/DATA_PROVENANCE.md` now give the verified design rather than
leaving S2 pending. `docs/MIGRATION_REPORT.md` marks S2 resolved and records S1's
substantive operator/interpretation mismatch. These audits and identity/census
metadata were added, and included-file checksums refreshed. Original source
documents, scientific code, figure labels and scientific input registries remain
unchanged. Release completion and rights review remain pending.

## Existing-document occurrence inventory

The following complete file-level inventory covers matching anchor/design terms
in migrated Markdown, excluding this new audit and the migration report already
listed above. Line locators are recorded in `S1_S2_EVIDENCE.json`; they describe
this follow-up working state. A match is a review candidate, not a claim that every
occurrence needs correction. Frozen source documents remain unchanged.

- `README.md`
- `analyses/A1_single_anchor/README.md`
- `analyses/A2_dual_anchor/README.md`
- `data/README.md`
- `docs/ANALYSIS_MAP.md`
- `docs/DATA_PROVENANCE.md`
- `manifests/provenance/source/RECOMB_MANUSCRIPT_TECHNICAL_AUDIT.md`
- `manifests/provenance/source/SUPPLEMENTARY_FIGURES_WAVE1_AUDIT.md`
- `manifests/provenance/source/docs/DMI_BRIDGE_A20_DUAL_STRONG_ANCHOR_QUALIFICATION.md`
- `manifests/provenance/source/docs/DMI_BRIDGE_A21_DUAL_ANCHOR_GEOMETRY.md`
- `manifests/provenance/source/docs/DMI_BRIDGE_A22_DUAL_ANCHOR_SIGN_ONLY_UTILITY.md`
- `manifests/provenance/source/docs/DMI_BRIDGE_A23_A1_A2_SYNTHESIS.md`
- `manifests/provenance/source/docs/DMI_BRIDGE_M1_MANUSCRIPT_EVIDENCE_CONSOLIDATION.md`
- `manifests/provenance/source/docs/DMI_BRIDGE_M2_A2_ROBUSTNESS_EVIDENCE_CONSOLIDATION.md`
- `manifests/provenance/source/docs/DMI_BRIDGE_PL1_PREDICTABILITY_LANDSCAPE.md`
- `manifests/provenance/source/docs/DMI_BRIDGE_PL2B_PRODUCTION_CONTRACT.md`
- `manifests/provenance/source/docs/DMI_BRIDGE_PL2C_GEOMETRY_UTILITY_DEVELOPMENT_CONTRACT.md`
- `manifests/provenance/source/docs/DMI_BRIDGE_PL2D_CONFIRMATION_CONTRACT.md`
- `manifests/provenance/source/docs/DMI_BRIDGE_PL2D_HYPOTHESIS_FREEZE_CONTRACT.md`
- `manifests/provenance/source/docs/DMI_BRIDGE_PL2D_H_SUPPORT_GATE_AMENDMENT_V1.md`
- `manifests/provenance/source/docs/DMI_BRIDGE_PL2_SIGN_ONLY_UTILITY.md`
