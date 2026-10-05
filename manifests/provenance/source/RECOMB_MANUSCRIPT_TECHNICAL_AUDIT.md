# RECOMB manuscript technical verification audit

Read-only audit of the repository state at commit `39844bab8129e75ea2a02a7727070bfe2c2a77fe` (2026-10-04) and the two DOCX snapshots supplied under `manuscript_audit_inputs/`. DOCX files were extracted as text for comparison; neither snapshot nor scientific code was modified. No reconstruction, optimization, solver, or flux-sampling run was performed. Parquet checks and summaries cited below were read-only; figures were not regenerated.

| ID | Topic | Verdict | Manuscript action |
|---|---|---|---|
| A1–A3 | Transcriptomic source, samples, preprocessing | VERIFIED_WITH_CORRECTION | REPLACE / COMPLETE PLACEHOLDER |
| B1–B2 | Parent GEM and flux reaction universe | VERIFIED_WITH_CORRECTION | REPLACE |
| C | Production contextualization settings | VERIFIED_WITH_CORRECTION | COMPLETE PLACEHOLDER |
| D | Frozen panel and sampling provenance | VERIFIED_WITH_CORRECTION | REPLACE / COMPLETE PLACEHOLDER |
| E1–E4 | DMI provenance, mapping, flux units | VERIFIED_WITH_CORRECTION | REPLACE / AUTHOR INPUT REQUIRED |
| F1–F5 | Figure 3 population, values, thresholds, old assets | VERIFIED_WITH_CORRECTION | REPLACE |
| G1–G3 | Held-out reaction exclusions | VERIFIED_WITH_CORRECTION | COMPLETE PLACEHOLDER |
| H | Figure 4C candidate provenance | VERIFIED_WITH_CORRECTION | REPLACE |
| I1–I5 | Supplementary figures and denominators | VERIFIED_WITH_CORRECTION | COMPLETE PLACEHOLDER |
| J | Frozen identities | VERIFIED | KEEP AS WRITTEN |
| K–L | Release and author-side information | OUTSIDE_REPOSITORY_SCOPE | AUTHOR INPUT REQUIRED |
| Literal scan | Full two-document comparison | VERIFIED_WITH_CORRECTION | REPLACE |

## A. Transcriptomic source and sample mapping

### A1. Dataset identity

**Claim/question.** Do final production inputs use Khalsa et al. GSE151414, with CT2A and GL261, three samples each?

**Verdict.** VERIFIED.

**Evidence.** `external_data/brain_glioma_GL261_CT2A/khalsa_GSE151414/metadata/sample_manifest.tsv` maps GSM4577664–7666 to CT2A and GSM4577667–7669 to GL261. The downloaded processed TPM source is recorded in `2_independent_transcriptomic_concordance/inputs/SOURCE_TPM_PROVENANCE.json` and `external_data/brain_glioma_GL261_CT2A/manifests/download_manifest.tsv`. Stage-11 `provenance/master_config.yaml`, sample/context manifests, and production ensemble references bind the final four-method panel to the two tumor inputs and four sample subsets. The manifest SHA256 recorded for `sample_manifest.tsv` is `88a1a614ec3d41fa86695b060102b9942250557176f055d86352a20f0b035fd0`.

**Correct technical statement.** The production panel uses the Khalsa et al. GSE151414 processed TPM input, with three CT2A and three GL261 samples.

**Manuscript action.** REPLACE: remove the verification request and state the confirmed identity.

### A2. Exact sample IDs and context labels

**Claim/question.** What are the six accessions, internal sample mapping, and relation to x/setx labels?

**Verdict.** VERIFIED_WITH_CORRECTION.

**Evidence.** The sample manifest gives the exact mapping below. Production RNA subsets are named by the `seta`/`setb` sample IDs. The Bridge/Figure 2–4 exports label the same positional subsets as `setx1`, etc.; those are analysis-context labels, not GEO accessions. The repository’s figure context keys preserve the subset positions, while the tumor identity selects whether the positions mean SetA/CT2A or SetB/GL261.

| Tumor | Position | Production ID | GEO accession | Pair/triple subsets |
|---|---:|---|---|---|
| CT2A | 1 | seta1 | GSM4577664 | x1+x2=`setx1,setx2`; x1+x3=`setx1,setx3`; x2+x3=`setx2,setx3`; all=`setx1,setx2,setx3` |
| CT2A | 2 | seta2 | GSM4577665 | same positional scheme |
| CT2A | 3 | seta3 | GSM4577666 | same positional scheme |
| GL261 | 1 | setb1 | GSM4577667 | same positional scheme |
| GL261 | 2 | setb2 | GSM4577668 | same positional scheme |
| GL261 | 3 | setb3 | GSM4577669 | same positional scheme |

The exact individual sample IDs represented in each subset are therefore seta1+seta2, seta1+seta3, seta2+seta3, or all three for CT2A; analogous setb combinations for GL261. `setx1,setx3` in Figure 4C is the positional x1+x3 context.

**Correct technical statement.** Keep the manuscript’s x1+x2 / x1+x3 / x2+x3 / all-three wording, and define x1–x3 as tumor-specific positions. Add the six GSM mappings; clarify `setx*` are Bridge context aliases.

**Manuscript action.** REPLACE.

### A3. Expression preprocessing

**Claim/question.** Aggregation, normalization, filtering, identifier handling, duplicate/missing genes, transformations, and GPR propagation.

**Verdict.** VERIFIED_WITH_CORRECTION.

**Evidence.** Stage-11 `master_config.yaml` explicitly gives `pooled_statistic: median`, detection threshold 0.1, and 0.25/0.75 quantiles. Thus combined samples use the median per gene, not the arithmetic mean. `SOURCE_TPM_PROVENANCE.json` identifies the input as already processed TPM; the downloaded gzip and uncompressed source hashes are recorded there (compressed SHA256 `19d3d8261c38907a153e2908b3b5bc5b907269a10a9f6a5598c0362765b5b5da`; decompressed data SHA256 begins `96e502`). Stage-11 expression processing maps gene identifiers to model genes, applies its configured TPM detection/quantile rules, and maps expression to reaction evidence through each algorithm’s GPR handling. Method implementations differ: CORDA constructs confidence categories; GIMME/iMAT use method-specific expression thresholds; RIPTiDe uses its expression-based reaction scoring and GPR processing. The repository does not establish one common post-aggregation log transformation or a universal duplicate/missing-gene rule shared by all algorithms. Do not state such a common rule without method-specific evidence. Source tables preserve gene-level TPM; missing model genes do not acquire measured TPM merely by identifier conversion.

**Correct technical statement.** “Within each selected sample subset, gene TPM values were summarized by the configured median (`pooled_statistic: median`); a 0.1 detection threshold and 0.25/0.75 quantile thresholds were used in the production expression-processing configuration. Reaction evidence was derived through method-specific GPR and contextualization procedures (CORDA confidence categories; GIMME/iMAT thresholds; RIPTiDe scoring). The processed source is TPM.” Complete any method-specific duplicate/missing-ID or transformation description only where explicit code/config establishes it; do not claim a single shared aggregation or GPR rule.

**Manuscript action.** REPLACE.

## B. Parent GEM and reaction universe

### B1. iMM1865 provenance

**Claim/question.** Parent model file, release, curation, and hash.

**Verdict.** VERIFIED_WITH_CORRECTION.

**Evidence.** The distributed local model is `metabolic_network/iMM1865/iMM1865.xml.xz`, SHA256 `3152e7854ff37dbf512483b6f1b15312958923def6f590d4e9180dab0ca3bc06`. It is identified in repository documentation as iMM1865 and cited to Khodaee et al. 2020. The companion exchange-bound file is `metabolic_network/iMM1865/iMM1865_exchange_bounds.csv` (its digest is present in the model provenance records). Production model records document local processing/projection, but no separate upstream release tag/version is established beyond the file identity and cited publication. The full parent count in Stage-11 QC is 10,612 reactions.

**Correct technical statement.** Identify the exact compressed model file and SHA256; cite Khodaee et al. (2020). Describe the production transformations as the documented boundary/medium treatment and reaction projection, not as an unrecorded upstream release.

**Manuscript action.** REPLACE.

### B2. Cache and reaction count semantics

**Claim/question.** Whether 4,181 is parent, projected universe, or cache dimension; anchor and weak-cue counts.

**Verdict.** VERIFIED_WITH_CORRECTION.

**Evidence.** `outputs/dmi_bridge_bio0r_qualification_v1/BRIDGEBIO0R_IDENTITY_AUDIT.json` records cache shape `[32,20,4181]` and reaction axis length 4,181; this is the frozen flux-cache axis, not the parent-model size. In that axis HEX1 and LDH_L each occur once. Excluding both yields 4,179 weak-cue reactions. The cache is a common aligned analysis axis; it is not a claim that all 4,181 reactions survive in every method-specific reconstructed GEM. Parent: 10,612; Stage-4 projection documents 4,747 in a separate projection context and must not be substituted for the Stage-11 parent or Bridge cache.

**Correct technical statement.** “The analysis cache stores 4,181 aligned reaction coordinates; HEX1 and LDH_L each occur once, leaving 4,179 reactions in the common weak-cue population. The iMM1865 parent has 10,612 reactions.” Projection, absent-reaction and blocked-reaction treatment is method/cache alignment documented in Stage-11 projection manifests; the 4,181 number is not the unreduced parent count.

**Manuscript action.** REPLACE.

## C. Contextualization settings

**Claim/question.** Production settings recoverable for iMAT, GIMME, CORDA, RIPTiDe and shared physiological constraints.

**Verdict.** VERIFIED_WITH_CORRECTION.

**Evidence.** Main authority is `11_recomb_flux_expert_generation/provenance/master_config.yaml`, together with `src/imat_canonical.py`, `src/gimme_repair.py`, `src/corda_canonical.py`, `src/riptide_native.py`, `src/qualification.py`, and method job manifests. Values explicitly configured: iMAT epsilon 0.01, feasibility/activity tolerance 1e-8, MIP gap 0.01, biomass floor 0.90 of parent Bmax; sole core is `BIOMASS_reaction`. GIMME threshold 0.1, biomass objective `BIOMASS_reaction`, retained objective fraction 0.90, numerical epsilon 1e-4. CORDA uses expression confidence categories and package 0.5.1; empty CORDA parameter map intentionally invokes constructor defaults `n=3`, `penalty_factor=100`, `support=5` (source comment/config documents this). RIPTiDe objective fraction 0.8, pruning true, task fraction 0.01, set_bounds true, GPR true, empty tasks/exclusions; native RIPTiDe GapSplit fitting requests 50 method-native samples. Common production objective/floor and medium/boundary settings are in the Stage-11 configuration and qualification records; they are not to be replaced by generic package documentation defaults. Oxygen/medium bounds are specified by the frozen medium/boundary files. The repository does not encode a single universal blocked-reaction rescue rule; reactions are retained/projected according to each method’s output and aligned to the cache. Reversibility is represented by model bounds. No independent new qualification was run for this audit.

**Correct technical statement.** Insert the above method-specific settings and state the common production biomass floor (0.90 of the relevant qualified capacity) and frozen boundary-safe medium. Mark other values “not recorded” where no production field exists; do not infer defaults. Current iMAT contract must be followed: one canonical `imat.run()`, biomass sole core, final topology union of raw retained indices and biomass core, with hard parent-Bmax 0.90 fit floor.

**Manuscript action.** COMPLETE PLACEHOLDER.

## D. Flux sampling and candidate states

**Claim/question.** Count and provenance of 32 ensembles, 20 vectors each; upstream draw details and downstream reuse.

**Verdict.** VERIFIED_WITH_CORRECTION.

**Evidence.** `11_recomb_flux_expert_generation/qc/PRODUCTION_FLUX_ENSEMBLE_AUDIT.tsv` records 32 qualified ensemble artifacts with 20 vectors each (640 total). `src/ensembles.py` uses COBRApy `OptGPSampler` for GIMME/iMAT/CORDA and native RIPTiDe GapSplit for RIPTiDe: OptGP was not used for every method. Requirements pin COBRApy 0.31.1 and RIPTiDe 3.4.81 (also Gurobi 13.0.1, Troppo 0.7, CORDA 0.5.1). Configured thinning is 100, process/chain count 1, seeds are in `manifests/ensemble_jobs.tsv`; 20 are the retained whole-vector rows. RIPTiDe configuration separately records native fit sample count 50. Solver is Gurobi where invoked; configured common feasibility tolerance is recorded in Stage-11 configuration. The repository does not provide a defensible single “all methods generated N then retained 20” draw total, explicit uniform burn-in, or common convergence diagnostic. Do not invent one. Frozen downstream PL1/PL2/A1/A2 stages reuse candidate vectors: their manifests/source audits mark solver and sampling not invoked; they reweight/describe existing states.

**Correct technical statement.** “The frozen candidate panel comprises 32 ensembles × 20 stored whole-flux vectors = 640. GIMME, iMAT and CORDA ensembles use COBRApy OptGPSampler with configured thinning 100, one process, and per-job seeds; RIPTiDe uses native GapSplit. The 20 stored candidates are not an effective sample size or proof of convergence. Downstream descriptor/utility/A1/A2 synthesis reused this panel without new flux sampling.” Omit unspecified warm-up/generated-draw/convergence claims.

**Manuscript action.** REPLACE.

## E. DMI provenance, anchors, and physical units

### E1. Repository DMI provenance

**Claim/question.** Animals, values, units, study, acquisition, model, cohort relation.

**Verdict.** VERIFIED_WITH_CORRECTION.

**Evidence.** `external_data/brain_glioma_GL261_CT2A/simoes_2025_DMI/metadata/README.md` identifies Simões et al., eLife 2025 and Dryad DOI `10.5061/dryad.905qfttwb`, describes DGE-DMI in CT2A and GL261, and states kinetic maps include Vmax, Vglx, and Vlac in mM/min. `metadata/subject_manifest.tsv` lists C1–C5 as CT2A and G1–G5 as GL261 (five each); repository MATLAB/source/data documentation describes tumor-region data and kinetic fitting outputs. RNA samples have GSE accession IDs while DMI subjects have C/G IDs from a distinct source study; no subject crosswalk exists, so repository evidence supports distinct cohorts, not paired animals. Full acquisition details may exist in archived source materials, but the manuscript’s requested final acquisition citation, ethics approval and author/legal declarations are not all encoded as project decisions.

**Correct technical statement.** State five DMI subjects per tumor and cite/identify the Simões DMI data deposit; Vmax, Vlac, Vglx are reported in mM/min. Describe RNA and DMI as distinct cohorts with no animal-level pairing. Exact acquisition/ethics details still require author verification against the source publication/protocol.

**Manuscript action.** COMPLETE PLACEHOLDER / AUTHOR INPUT REQUIRED.

### E2. Vmax-to-HEX1 operator

**Claim/question.** Whether Vmax is quantitatively converted/matched to HEX1 or used as rank information.

**Verdict.** VERIFIED_WITH_CORRECTION.

**Evidence.** Frozen A1/A2 operator registry and source code in the Bridge `a20`/`a21`/`a23` records bind the Vmax rank coordinate to HEX1. The operator compares rank coordinates and applies squared rank-distance exponential weights with per-mouse temperature calibration to global ESS 20; it does not equate absolute Vmax with an absolute HEX1 flux. Rank transforms make this comparison dimensionless; no tissue mass/volume or minute-to-hour conversion is in the operator.

**Correct technical statement.** “Vmax is associated with the HEX1 candidate coordinate through a dimensionless rank-based soft-weighting operator; it is not an absolute flux equality or calibrated unit conversion.”

**Manuscript action.** REPLACE.

### E3. Lactate orientation

**Claim/question.** A2-L lactate coordinate and sign inversion.

**Verdict.** VERIFIED.

**Evidence.** Frozen anchor operator/registry defines the lactate candidate observable as `max(-v_LDH_L, 0)`. The minus sign follows the model’s LDH_L reaction orientation. It is used as an anchor coordinate, then ranked; it is not a direct absolute DMI/GEM unit conversion.

**Correct technical statement.** Retain `max(−v_LDH_L, 0)` and explain it is the nonnegative lactate-producing orientation of the model reaction.

**Manuscript action.** KEEP AS WRITTEN.

### E4. Flux units and gain scale

**Claim/question.** Whether gains can be labeled mmol gDW⁻¹ h⁻¹ and whether DMI conversion occurred.

**Verdict.** CONTRADICTED.

**Evidence.** Candidate vectors and all gain tables are in the stored native/model flux coordinate scale. DMI input is mM/min; frozen rank operators use only normalized rank coordinates. Inspection of A1/A2 operator code, manifests, and source tables found no conversion formula, tissue-volume/dry-weight factor, or minutes-to-hours conversion applied to candidate fluxes or gains. The DMI rates inform ranks only. Thus reported gains are numerical model-flux-unit gains, not calibrated physical flux rates. Figure3 source tables do not carry a physical-unit calibration.

**Correct technical statement.** Replace all `mmol gDW⁻¹ h⁻¹` gain/axis labels with “native GEM flux units” or “model flux units” and state that the DMI-to-GEM operator is rank based and dimensionless. Absolute gains have the candidate model’s stored scale; they are not physically calibrated.

**Manuscript action.** REPLACE.

## F. Figure 3 aggregation and values

### F1. Population and aggregation

**Claim/question.** 660,237 pairs per anchor and 1,440,313 non-tie cases; definitions and weighting.

**Verdict.** VERIFIED.

**Evidence.** Authoritative final source is `figures/fig3/tables/fig3_gain_reaction_evaluation.parquet.xz`, SHA256 `7a83a31d29317f7b0d906fff73691428d98427d93d2fdebd4ef91294f8a46676`; metadata `fig3_gain_distribution_metadata.json` SHA256 `d6ebfaa08aec0255104463bd327c55a74f24a1d074589d563782cb5953ab14a7`. A case is one eligible reaction × truth-pair/evaluation record for a distinct matched non-tie synthetic truth. A reaction–evaluation pair is one reaction and one mouse/method/RNA evaluation; distinct non-tie truths within it are first arithmetically averaged, then each eligible pair has equal weight. Each anchor has 660,237 pairs. Matched population: 836 truth pairs × 4,179 reactions = 3,493,644 cases, of which 1,440,313 non-tie and 2,053,331 tie cases.

**Correct technical statement.** Keep stated counts and define case and pair as above; Figure 3 pair means do not pool truth cases with unequal counts.

**Manuscript action.** KEEP AS WRITTEN, after replacing unit claims.

### F2–F3. Means and endpoint separation

**Claim/question.** Pair-level means and correct-minus-wrong separation.

**Verdict.** VERIFIED.

**Evidence.** `fig3_gain_summary.parquet.xz` (SHA256 `0529bab84c6ad4cfc4a3c9e0d07ad418afcab3be52d2254183400dcd679af482`) yields A1 correct 0.010387 and wrong −0.013747; A2-L correct 0.022294 and wrong −0.030611. `fig3_directional_advantage_summary.parquet.xz` (SHA256 `d2ec6a86e344877f878b63e9ad03c5de0823f821b7527b8e4191f1ab5d39afab`) stores unrounded separations A1 `0.024133718031918716`, A2-L `0.052905717474604916`; these round to 0.02413 and 0.05291. Values are equal-weight pair-level summaries.

**Correct technical statement.** The supplied values are accurate after rounding; identify them as model-flux-unit gains, not physical units.

**Manuscript action.** KEEP AS WRITTEN with unit replacement.

### F4. Threshold curves

**Claim/question.** Final zero fractions, exported thresholds, and prespecified nonzero threshold.

**Verdict.** VERIFIED_WITH_CORRECTION.

**Evidence.** Read-only calculation over `fig3_gain_reaction_evaluation.parquet.xz`, grouped by `anchor_setting`, gives correct-cue pair mean gain above zero: A1 `0.4149918893`, A2-L `0.4389287483`; correct-minus-wrong specificity above zero: A1 `0.3575261611`, A2-L `0.4101845246`. The final source table contains pair rows, not a frozen threshold/fraction export. The exploratory notebook `figures/fig3/fig3_gain_exploration.ipynb` evaluates `np.logspace(-15,0,250)`; this is an exploratory grid, not evidence of a prespecified nonzero threshold. No prespecified nonzero threshold was found.

**Correct technical statement.** Report zero exceedance fractions only if desired, with the endpoint definition and pair-level denominator. Do not claim a prespecified nonzero threshold. If retaining a nonzero cutoff, first document an independent prespecification; do not select one retrospectively.

**Manuscript action.** REPLACE / REMOVE CLAIM.

### F5. Old and superseded Figure 3 assets

**Claim/question.** Whether exploratory pooled-case numbers/labels remain in assets.

**Verdict.** VERIFIED_WITH_CORRECTION.

**Evidence.** `figures/fig3/candidates/fig3_panelC_exploratory_histogram_A1.svg` embeds `n=660,237; mean=0.02413; median=0; positive=0.358`; the `mean=0.02413` is endpoint separation, not that panel’s displayed per-pair gain distribution. Exploratory candidate plots and `fig3_gain_definition_audit.parquet.xz` retain pooled-case alternatives. In the case-pooled audit, A1 correct mean 0.014234, wrong −0.018859, separation 0.03309389; A2-L correct 0.02280187, wrong −0.03173328, separation 0.05453515. These are superseded for current pair-weighted Figure 3. The older raster files `figures/fig3/fig3D.png`, `fig3ABCD.png`, and exploratory candidate SVG/PNGs should not be treated as final artwork without label/source review.

**Correct technical statement.** Use only the pair-level source table and replace/remove any pooled-case values where a pair denominator is claimed. Audit all candidate artwork labels before submission.

**Manuscript action.** REPLACE.

## G. Held-out sensitivity and exclusions

### G1–G3. Proximal/pathway sets and held-out correlations

**Claim/question.** Frozen proximal set, pathway convention, union, and denominators/correlations.

**Verdict.** VERIFIED_WITH_CORRECTION.

**Evidence.** Frozen registry is `outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_REACTION_REGISTRY.tsv` (SHA256 `30dad204d329bf7fec9cc0156a0944c7785f7676c31f2fc3a96f599b8f98095c`), fields `reaction_id`, `subsystem`, `same_subsystem_as_HEX1`, `proximal_set_member`. Exactly 80 rows are marked `proximal_set_member=true`: ABTArm, ACITL, ACONT, ACONTm, ACTLMO, ACTNMO, ACYP, AKGDm, ALATA_L, ALATA_Lm, ALCD21_D, ALCD22_D, ALCD22_L, ALDD2x, ALDD2xm, ALDD2y, ALR2, ALR3, ASPT, CBPPer, CITL, CSm, DPGM, DPGase, ENO, FBA, FUM, FUMm, G3PD2m, G6PPer, GAPD, GDHm, GLNS, GLUDC, GLUDym, GLUNm, GLYOX, GLYOXm, GTHOm, GTHOr, HMR_3855, HMR_7748, HMR_7749, ICDHxm, ICDHym, ICDHyp, ICDHyr, LALDO, LALDO2, LCADi, LCADi_D, LCADm, LCARS, LDH_L, LDH_Lm, LGTHL, MDH, MDHm, MDHp, ME1m, ME2m, MGSA2, PDHm, PFK, PGI, PGK, PGM, PGMT, PPDOy, PYK, SUCD1m, SUCOAS1m, SUCOASm, TPI, r0173, r0202, r0354, r0355, r0509, r0818. The source code `scripts/dmi_bridge_pl2a_prepare_v1.py` defines proximal membership as any reaction annotated to one of four subsystems: Glycolysis/gluconeogenesis, Pyruvate metabolism, Citric acid cycle, or Glutamate metabolism. The registry has 27 rows marked `same_subsystem_as_HEX1=true`; this requires exact equality with HEX1’s subsystem annotation. Multi-annotation handling follows the single `subsystem` field in the pinned reaction panel; unmapped annotations are not shared membership. Union is logical OR of the flags. `outputs/dmi_bridge_pl2d_confirmation_v1/BRIDGEPL2D_SUPPORTIVE_RESULTS.json` and `outputs/dmi_bridge_m1_manuscript_evidence_v1/BRIDGEM1_DEVELOPMENT_CONFIRMATION.tsv` provide original A1 held-out evidence: primary 838 finite reactions, rho 0.8351436; coverage-adjusted partial Spearman rho 0.8056872, n=838; proximal exclusion rho 0.8386984, n=816; same-subsystem exclusion rho 0.8378450, n=829; union rho 0.8386984, n=816. The frozen proximal set’s exact IDs and cardinality are the rows marked `proximal_set_member=true` in the named registry; finite denominators above are after exclusion and finite-value intersection. The manuscript should copy that exact registry list/count into the supplement, not derive a new set.

**Correct technical statement.** Correlations and counts above are the original A1 one-shot held-out confirmation, not development or A2-L. Same pathway uses frozen registry membership; union uses OR. Cite registry and report its marked reaction IDs/size in the final supplement.

**Manuscript action.** COMPLETE PLACEHOLDER.

## H. Figure 4C geometry-only selection

**Claim/question.** FACOAL204 example identity, candidate screen/ranking, cross-method behavior and leakage control.

**Verdict.** VERIFIED_WITH_CORRECTION.

**Evidence.** `figures/fig4/data/fig4_data_manifest.json` records 100 geometry candidates and `utility_used_for_ranking:false`; `figures/fig4/build_fig4_tables_updated.py` ranks candidate groups by descending `max(cross_method_range_A1_H_dir,cross_method_range_A2_H_dir)`, then geometry completeness, paired-shift finite count, minimum ESS, and candidate ID ascending. The final notebook hardcodes candidate `FCG_ed8b1873c8dbc03fe5eb`: FACOAL204, `setx1,setx3`, C1 vs G1. In the exported 100-row candidate table this row is geometry rank 23, and has `range_A2_minus_A1_H_dir_shift=1.085471`, the largest such shift in the 100-row screen. Thus it is accurate to say the example has the largest cross-method range of A2-L-minus-A1 entropy shift among those 100, but not that this was the primary sorting rule used to produce the top-100 screen. No utility field is used by ranking code; manifest says false.

Exact metric comparison from candidate metrics: GIMME H_dir approximately 0→0.630930 (more mixed); iMAT approximately 0→0.627654 (more mixed); CORDA 0.628923→0.174381 (more one-sided, dominant mass 0.533193→0.952336); RIPTiDe H_dir=0 and ΔvB width=0 under both anchors, with exact degenerate support at zero. The source rows are keyed by candidate group ID and method/evaluation IDs in `fig4_panelC_candidate_metrics.parquet.xz`; truth/evaluation IDs are included there. Ties are resolved by the listed deterministic secondary keys, ending in ascending candidate_group_id. Utility outcomes were absent from the ranking expression.

**Correct technical statement.** Clarify two-stage provenance: “From the 100 geometry-screen candidates ranked by the frozen geometry score, FACOAL204 (setx1,setx3; C1 vs G1) has the largest cross-method range in the A2-L-minus-A1 entropy shift. Ranking used geometry/support fields only; no utility outcome entered selection.” Method behavior in the passage is supported. Figure S7 may present the zero-degenerate RIPTiDe example as a limited-change counterpart.

**Manuscript action.** REPLACE.

## I. Supplementary figures, denominators, and controls

### I1. Figure S1 populations

**Claim/question.** ESS, within-stratum, product-weight, holdout, and matched-pair denominators.

**Verdict.** VERIFIED_WITH_CORRECTION.

**Evidence.** `outputs/dmi_bridge_a20_dual_anchor_qualification_v1/` and A23 manifests define global calibration over 320 candidates per tumor/mouse, ESS target 20; 160 method×RNA×tumor×mouse conditional strata per anchor, 20 candidate support each (3,200 weight rows per anchor); 20×20 product contrast support per evaluation and 400 evaluations; 1,671,600 evaluation-reaction keys per anchor; holdout utility has 836 matched truth pairs. Each global ESS is computed once per tumor×mouse×anchor weight system, not once per reaction. Conditional ESS is per stratum/mouse/anchor; product ESS is per evaluation. Small effective supports and degenerate sign/magnitude support are real finite-support diagnostics and should be shown from frozen source tables, not silently removed.

**Correct technical statement.** Label all three ESS levels separately; do not replicate a global ESS reaction-wise. State 836 matched truth pairs and the relevant support level for each panel.

**Manuscript action.** COMPLETE PLACEHOLDER.

### I2. Figure S2 context-specific confirmation

**Claim/question.** Four method × four RNA-context A1 held-out correlations, development vs confirmation vs A2-L.

**Verdict.** VERIFIED.

**Evidence.** `outputs/dmi_bridge_m1_manuscript_evidence_v1/BRIDGEM1_CONFIRMATION_CONTEXTS.tsv` contains original A1 one-shot held-out values (rho, finite n):

| Method | x1+x2 | x1+x3 | x2+x3 | x1+x2+x3 |
|---|---:|---:|---:|---:|
| CORDA | 0.861571 (638) | 0.835509 (636) | 0.803326 (625) | 0.917413 (622) |
| GIMME | 0.872555 (417) | 0.862306 (428) | 0.752206 (428) | 0.829166 (427) |
| RIPTiDe | 0.818001 (34) | 0.863043 (34) | 0.275610 (39) | 0.925513 (37) |
| iMAT | 0.879262 (593) | 0.771652 (592) | 0.867976 (588) | 0.835907 (593) |

Overall original A1 confirmation is n=838, rho=0.835144. Matched A2-L held-out is descriptive/post-freeze robustness (overall n=843, rho=0.849366), not a new independent confirmation. Development matched A1 n=2,454/rho=0.857116; original broad A1 analysis has n=2,455 because it includes the extra reaction outside the matched 4,179 universe. These are distinct populations and labels must preserve that distinction.

**Correct technical statement.** Table above is A1 held-out confirmation. Label A2-L as post-freeze descriptive analysis and development separately.

**Manuscript action.** KEEP AS WRITTEN, with explicit split/anchor labels.

### I3. Figure S3 descriptors

**Claim/question.** Primary association, information advantage definition and finite counts.

**Verdict.** VERIFIED.

**Evidence.** Frozen A1 held-out supportive results: direction-explained magnitude variance vs usefulness frequency rho=0.835144, n=838; vs information advantage rho=0.345421, n=838. Information advantage is half the correct-minus-wrong gain difference, i.e. the difference from the average of the two endpoints. Undefined descriptors (notably zero/near-zero weighted magnitude variance for direction-explained variance) remain missing and are excluded pairwise from finite correlations; they are not zero-imputed. Additional descriptor rows/counts are in `BRIDGEPL2D_SUPPORTIVE_RESULTS.json` and frozen S3 tables.

**Correct technical statement.** Retain values and definition; state finite counts and do not replace undefined descriptor values with zero.

**Manuscript action.** KEEP AS WRITTEN.

### I4. Figure S4 bootstrap and exclusions

**Claim/question.** Bootstrap procedure, interval, original/adjusted/exclusion statistics and interpretation.

**Verdict.** VERIFIED.

**Evidence.** `BRIDGEM1_MANUSCRIPT_EVIDENCE` records 5,000 reaction-level bootstrap replicates, seed 20260929, percentile interval 0.801562911874–0.863127542338 around original rho 0.835143582318 (n=838). Coverage-adjusted partial Spearman is 0.8056871576 (n=838); proximal exclusion 0.8386984413 (n=816); same subsystem/pathway exclusion 0.8378450479 (n=829); union exclusion 0.8386984413 (n=816). Bootstrap resamples finite reaction-level pairs, not animals or biological cohorts; interval is computational/resampling stability.

**Correct technical statement.** Retain these values with the stated finite n and describe the interval as reaction-level resampling stability, not biological sampling uncertainty.

**Manuscript action.** KEEP AS WRITTEN.

### I5. Figure S5 gain tails and reliability

**Claim/question.** Gain-level unit, λ/q values, classification level, tie tolerance, random-sign interpretation.

**Verdict.** VERIFIED_WITH_CORRECTION.

**Evidence.** Frozen utility tables and scripts define gains per non-tie truth case, then Figure 3 aggregation averages cases within reaction–evaluation pairs. Reliability/negative-tail exports therefore must identify whether they summarize case-level gains or pair-averaged gains; do not transfer one denominator label to the other. λ=0 is baseline recovery; λ=0.25 is primary. The frozen reliability summary evaluates q=0, 0.5, and 1 at λ=0.25; q=0.5 is the random-sign control and carries no directional information. Information advantage for a cue setting is half correct-minus-wrong gain; relative to itself at q=0.5 it is zero by symmetry. A finite-support/operator-induced baseline-relative change at q=0.5 is not information supplied by random data. Classification tables (where included) classify the stated averaged gain at the table’s level; truth ties have no correct/wrong direction and are excluded from directional classification. Tied classification uses the frozen numerical gain tolerance, not truth-tie status.

**Correct technical statement.** Figure S5 captions must state unit of aggregation for each panel, λ and q definitions, and tolerance for improved/tied/harmed. Retain the interpretation that the analysis does not estimate reliability of a real biological source. Replace physical units per E4.

**Manuscript action.** COMPLETE PLACEHOLDER.

## J. Frozen artifact identities

**Claim/question.** What each listed digest identifies and whether current frozen manifests reproduce it.

**Verdict.** VERIFIED.

**Evidence.** The A23 synthesis manifest `outputs/dmi_bridge_a23_a1_a2_synthesis_v1/BRIDGEA23_MANIFEST.json` binds the qualified upstream sources and source audit (`BRIDGEA23_SOURCE_AUDIT.json`, PASS). The recorded logical content identities are reproduced by their source audit/manifest fields: A1 descriptor `46aa9f0410f89d03b020e01fefa517348faae1d55d9bfe8387ac04b39bb9ac93` identifies the A1 descriptor logical table content; A2-L descriptor `50bc3e47239e4b22d0931246f4fa0d5d126a5b69226d47cb342cb8a5b3ee99d2` identifies the A2-L descriptor content; original A1 case `c3670992a0ad33698af6c7620429bed48e0ecaa98adb3c87e648eefb1ba70363` identifies the original case table/content; A2-L case `70321eb7d322f06a3c1e29c048d08f9b4c55d291763b5194f6d626325ed786f4` identifies the A2-L case content. These are logical-content identities as defined by the generating audit, not archive file SHA256s. Current A23 source audit confirms the identity bindings.

**Correct technical statement.** Preserve the hashes with their logical dataset meaning and label them “logical-content identity”; do not describe them as archive checksums.

**Manuscript action.** KEEP AS WRITTEN, but bind to the future public release manifest when it exists.

## K–L. Reproducibility and author-side matters

**Claim/question.** Current commit, environment, commands, public release/DOI, controlled access, citations, ethics, declarations.

**Verdict.** OUTSIDE_REPOSITORY_SCOPE.

**Evidence.** Current Git HEAD is `39844bab8129e75ea2a02a7727070bfe2c2a77fe`. Stage-11 `requirements.txt` and configuration document package/solver versions; project README and stage scripts provide commands and artifact/table locations. The repository includes GSE151414 documentation and the DMI Dryad DOI above. No public repository URL, archived release DOI, future tag, or final controlled-access/legal decision is established by current repository state. Final DMI acquisition/kinetic-model citation formatting, animal ethics identifiers, conflict-of-interest and author declarations require author confirmation. Existing local input-data files do not establish permission to redistribute them.

**Correct technical statement.** Populate reproducibility details that are already versioned (commit, requirements, commands, files) only after authors select a release. Leave public URL/archive DOI, future tag, controlled-access decision and declarations as author/publishing actions. Repository evidence does not establish the requested legal or declaration decisions.

**Manuscript action.** AUTHOR INPUT REQUIRED.

## Literal manuscript consistency findings

Both complete DOCX text streams were extracted and searched, including abstract, methods, results, discussion, captions, supplementary methods and tables. The following literal conflicts or unresolved phrases were found; repeated occurrences of an affected unit/claim should be changed consistently.

| Document | Location | Current text/value | Issue | Verified replacement | Evidence |
|---|---|---|---|---|---|
| Main | Abstract | “0.01039 to 0.02229 (in flux units [mmol DW-1 h-1])” | Unsupported physical calibration | “in native GEM flux units” | Frozen candidate vectors/source tables; no conversion in anchor/gain pipeline |
| Main | Methods, Figure 3 caption and Results | “mmol gDW⁻¹ h⁻¹” / “under the accepted flux conversion” | No accepted conversion exists in repository | “native GEM flux units” | E4 evidence above |
| Supplement | S1 preprocessing | asks whether sample combination was averaged | Production uses configured median, not arithmetic mean | State median across selected sample TPM values | Stage-11 master config `pooled_statistic: median` |
| Main and Supplement | Figure 4C caption/prose | says candidate was selected “by the largest … shift” | The top-100 ranking rule differs; chosen row maximizes shift within exported 100 | State two-stage selection and disclose the primary top-100 sort rule | Figure4 build code and manifest |
| Supplement | S1 | calls x1/x2/x3 tumor-specific samples without GSM mapping | Incomplete traceability | Add six accessions and explain positional setx aliases | GSE sample manifest |
| Supplement | S1 model | “4,181-reaction flux cache” adjacent to parent-model statement | Could be read as parent count | Explicitly distinguish 10,612 parent reactions from 4,181 cache coordinates | Stage-11 QC and Bridge identity audit |
| Main | Abstract/discussion and Figure 3 text | physical-unit gain wording repeated | Repeated instance of unsupported calibration | Replace all with model-flux-unit wording | Source tables and operator contain ranks, no conversion |
| Supplement | DMI references | `[CITE DMI dataset and kinetic model source]`, `[Add the verified DMI acquisition and kinetic-model references.]` | Citation details are not fully established by project provenance | Cite identified Simões/Dryad source; authors verify acquisition and kinetic model citations | Local DMI metadata README and author-side gap |
| Main and Supplement | Availability sections | public URL, archive DOI, controlled-access prompts | Future/external author decisions | Complete at release; do not invent values | Current repository metadata |

No other contradictory occurrence was found for the checked pair-level gain numbers, 660,237 pair denominator, 1,440,313 non-tie count, 836 matched truth pairs, 4,179 weak-cue count, λ=0.25 primary setting, held-out A1 overall result (0.835/838), A2-L descriptive result, HEX1/LDH_L mapping, or FACOAL204/context/mouse identity. One contextual-count distinction must remain explicit: matched development A1 is n=2,454, whereas the original wider-universe A1 analysis is n=2,455; these are not contradictory when each is labeled by its reaction universe. The manuscript’s direction-explained magnitude variance definition, information advantage definition, λ and q interpretation, and development/held-out labels otherwise agree with frozen analysis descriptions.

## Targeted highlighted-passage dispositions

The following maps each supplied unresolved passage to the technical finding above; each passage has one primary disposition.

| Passage | Verdict | Manuscript action | Disposition |
|---|---|---|---|
| MAIN-Q1 | VERIFIED | REPLACE | Confirmed iMM1865 and GSE151414; remove query. |
| MAIN-Q2 | VERIFIED_WITH_CORRECTION | REPLACE | Add exact IDs and median aggregation; define positional x labels. |
| MAIN-Q3 | CONTRADICTED | REPLACE | Physical unit unsupported; correct artwork labels and use model-flux units. |
| MAIN-Q4 | CONTRADICTED | REPLACE | Remove physical calibration assertion. |
| MAIN-Q5 | VERIFIED_WITH_CORRECTION | REPLACE | Zero fractions available; no prespecified nonzero cutoff established. |
| MAIN-Q6 | VERIFIED_WITH_CORRECTION | REPLACE | Behavior correct; revise selection description to two-stage ranking. |
| MAIN-Q7 | OUTSIDE_REPOSITORY_SCOPE | AUTHOR INPUT REQUIRED | Supply release URL/DOI and access decision after publication decisions. |
| SUPP-Q1 | VERIFIED_WITH_CORRECTION | REPLACE | Add six GSM accessions and clarify aliases. |
| SUPP-Q2 | VERIFIED_WITH_CORRECTION | REPLACE | Median and configured thresholds; method-specific GPR; some ID edge rules are not universal/recoverable. |
| SUPP-Q3 | VERIFIED_WITH_CORRECTION | REPLACE | File/hash established; distinguish parent count from cache dimension. |
| SUPP-Q4 | VERIFIED_WITH_CORRECTION | COMPLETE PLACEHOLDER | Settings listed in section C; report missing values as unrecovered. |
| SUPP-Q5 | VERIFIED_WITH_CORRECTION | REPLACE | Sampler differs for RIPTiDe; generated draw/burn-in/convergence specifics incomplete. |
| SUPP-Q6 | VERIFIED_WITH_CORRECTION | AUTHOR INPUT REQUIRED | Deposit and cohort facts established; final protocol/citation/ethics details need source/author confirmation. |
| SUPP-Q7 | CONTRADICTED | REPLACE | No physical conversion; retain rank-operator explanation and model units. |
| SUPP-Q8 | VERIFIED | KEEP AS WRITTEN | Bind to named pair-level table and SHA256 in F1. |
| SUPP-Q9 | VERIFIED_WITH_CORRECTION | COMPLETE PLACEHOLDER | Frozen registry fields and statistics given in G; include exact flagged set listing from registry. |
| SUPP-Q10 | OUTSIDE_REPOSITORY_SCOPE | AUTHOR INPUT REQUIRED | Versioned metadata exists; public archive and distribution choices remain. |
| SUPP-Q11 | VERIFIED_WITH_CORRECTION | COMPLETE PLACEHOLDER | Qualify aggregation levels and replace unsupported physical units. |
| SUPP-Q12 | VERIFIED_WITH_CORRECTION | REPLACE | Explain top-100 then maximum-shift selection. |
| SUPP-Q13 | VERIFIED_WITH_CORRECTION | COMPLETE PLACEHOLDER | Population levels given in I1. |
| SUPP-Q14 | VERIFIED | KEEP AS WRITTEN | Values and n in I2; label split and anchor. |
| SUPP-Q15 | VERIFIED | KEEP AS WRITTEN | Definition and values in I3. |
| SUPP-Q16 | VERIFIED | KEEP AS WRITTEN | Values/interval interpretation in I4. |
| SUPP-Q17 | VERIFIED_WITH_CORRECTION | COMPLETE PLACEHOLDER | Counts and q/λ meanings are supported; exact table classification label/tolerance must accompany any percentages. |
| SUPP-Q18 | VERIFIED | KEEP AS WRITTEN | Logical content identities bind as described in J; public release binding remains future. |
| SUPP-Q19 | VERIFIED_WITH_CORRECTION | REPLACE | Priorities divide into repository-resolved and author-side items above. |
| SUPP-Q20 | OUTSIDE_REPOSITORY_SCOPE | AUTHOR INPUT REQUIRED | Local deposit citation found; final acquisition/kinetic-model references need author verification. |

## 1. Verified manuscript corrections

- Replace every calibrated physical flux unit with “native GEM flux units” unless authors can provide an independently documented conversion with model-scale assumptions.
- State median aggregation of TPM samples and provide GSM/sample mapping.
- Distinguish the 10,612-reaction parent from the 4,181-coordinate cache and 4,179 weak-cue subset.
- Describe the sampler separately for OptGPSampler methods and native RIPTiDe GapSplit; distinguish stored vectors from sampler internals.
- Correct Figure 4C’s selection provenance to describe the top-100 rule and the selected candidate’s maximum entropy-shift range within that set.
- Use the pair-level Figure 3 source table, pair means, and counts; do not reuse case-pooled exploratory outputs as pair-level results.
- State there was no prespecified nonzero threshold found; zero-exceedance fractions are available from read-only calculation.

## 2. Remaining author-side questions

- Final DMI acquisition and kinetic-model citations, exact protocol/analysis details and ethics identifiers.
- Final author contribution, conflict-of-interest and other declarations.
- Public repository URL, archived DOI/release tag, redistribution and controlled-access decision.
- Confirm final citation formatting against the source DMI paper/materials before submission.

## 3. Remaining computational questions

- None requires new expensive computation for the conclusions reported here. Exact sampler warm-up/generated draw counts and a common convergence diagnostic are not recoverable from current frozen records; resolving those requires provenance records or a new explicitly authorized sampling run. This audit did not run one.
- No proximal-set recomputation was performed; its exact 80 IDs and rule are reported under G.

## 4. Stale/superseded artifacts discovered

- Exploratory Figure 3 candidate plots and `fig3_gain_definition_audit.parquet.xz` retain case-pooled means/separations; these differ from final pair-level values. In particular the exploratory A1 separation is 0.03309389 and A2-L 0.05453515.
- `figures/fig3/candidates/fig3_panelC_exploratory_histogram_A1.svg` embeds a misleading `mean=0.02413` with `n=660,237`; 0.02413 is the pair-level endpoint separation, not the depicted gain mean. Do not use it as final artwork.
- `figures/fig3/fig3D.png`, `fig3ABCD.png` and other candidate PNG/SVGs require label/source review before submission; current authoritative numeric data are the final tables in `figures/fig3/tables/`.
- Historical Stage-11 E-Flux and superseded-generation artifacts exist in provenance/QC directories. Current primary panel is four methods; do not include E-Flux based on historical files.
- Old `figures/fig4/build_fig4_tables.py` and earlier draft figure assets coexist with the updated builder and current Figure 4 candidate tables. Use the current updated script/manifest and verify artwork identity.

AUDIT COMPLETE: 6 verified, 15 verified with correction, 3 contradicted, 0 not determinable, 3 outside repository scope, 0 requiring recomputation.
