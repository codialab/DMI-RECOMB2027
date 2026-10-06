# Evaluation design audit: S2

Audit date: 2026-10-06. Source commit:
`16379a169f45ec556274953eb66ae12293e1e985`. Metadata were inspected read-only;
no reconstruction, optimization, sampling, solver or numerical-result rerun occurred.

## Verdict and arithmetic

**RESOLVED — authoritative terminology established.** The design is:

```text
4 methods × 4 RNA contexts × 5 CT2A mice × 5 GL261 mice
= 16 method/RNA contexts × 25 cross-tumor mouse pairs
= 400 computational evaluations.
```

Each factor of five denotes distinct biological DMI subjects, not random draws.
One evaluation uses one CT2A mouse and one GL261 mouse at one fixed method/RNA
context. The full study uses ten DMI mice. There is no independent five-draw
dimension in the frozen evaluation or ensemble-job registry. “Two mice × five
draws” is incorrect as a description of the study design; “two mice per evaluation”
is accurate only when the two separate five-subject cohorts are also explained.

## Authoritative identity construction

The registry is source
`outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_EVALUATION_REGISTRY.tsv`,
SHA-256 `aa468dcec11dea764682bd3d7eb14c03c003063e74d92080d0bfaabb366290a3`.
Its 400 rows have 400 unique tuples
`(algorithm, rna_context_key, ct2a_mouse, gl261_mouse)` and exactly equal the
Cartesian product of their four identity axes. This census explicitly checked
missing/duplicated combinations, not just the total row count. Every method/RNA
context contains all 25 mouse crosses, once each. There is no draw field.

| Axis | Frozen identities |
|---|---|
| Method | iMAT, GIMME, CORDA, RIPTiDe |
| RNA context | `training_samples=setx1,setx2`; `setx1,setx3`; `setx2,setx3`; `setx1,setx2,setx3` (all carry the `training_samples=` prefix) |
| CT2A DMI mouse | C1, C2, C3, C4, C5 |
| GL261 DMI mouse | G1, G2, G3, G4, G5 |

`scripts/dmi_bridge_pl1_predictability_resume_v1.py::_context_specs`, lines 77–132,
pairs one CT2A and one GL261 ensemble sharing method and RNA-context key, verifies
five mice on each side, then nests the CT2A and GL261 mouse loops. The original
PL1 adapter has the same construction at `dmi_bridge_pl1_predictability_landscape_v1.py:346–414`.
PL2A reads that exact registry (`dmi_bridge_pl2a_prepare_v1.py:281–282`); A21/A22
preserve these evaluation identities for the matched analyses.

For one evaluation, mouse-specific weights are restricted and renormalized over
the fixed 20 vectors of each condition's ensemble. Contrasts are CT2A minus GL261;
joint weights are the products of the two normalized condition weights
(`dmi_bridge_pl1_predictability_core_v1.py:111–125`). The resulting **20×20 = 400
candidate contrasts per evaluation** must not be confused with **400 evaluations
across the study**. Neither count supplies 400 independent animals.

## Biological subjects and RNA contexts

Source `external_data/brain_glioma_GL261_CT2A/simoes_2025_DMI/metadata/subject_manifest.tsv`
has ten rows, C1–C5 mapped to CT2A and G1–G5 to GL261; SHA-256
`bf631c2549b137100dd1a4cee3cdb986d19cd35e9091bbfba64cbdbeb5a43f8f`.
The accompanying local README describes these as two mouse cohorts and maps the
source subject axes to C/G identifiers. The canonical DMI loader and median table
producer carry these identities into `mouse_id`, used by the weighting operators.
No measurement values or raw maps were copied or printed in the evidence census.

RNA contexts are training sample subsets, not four random draws: the three
two-sample subsets and the complete three-sample subset. Tumor-specific sample
names are `seta1–seta3` and `setb1–setb3`, mapped to the common `setx` context key.
The mapping is implemented by `_canonical_rna_sample_context` and
`_add_rna_context_key` in source
`12_recomb_method_comparison/src/recomb_compare/stage11.py:21–69`, which replace
the tumor-specific `seta`/`setb` prefixes while preserving sample-subset identity.
The DMI subjects are distinct from RNA sample IDs; no animal-level crosswalk is
established. Cross-tumor evaluation pairs are all computational crosses, not
experimentally matched pairs of animals. These shared inputs preclude treating
the 25 crosses as independent biological replication.

## Ensembles, vectors and actual random operations

| Object | Verified count/meaning | Evidence |
|---|---|---|
| Reconstruction ensembles | 4 methods × 2 tumor contexts × 4 RNA subsets = **32** | Stage-11 `manifests/ensemble_jobs.tsv`, 32 unique ensemble identities, SHA-256 `3357f216e0efed3ebf532bf7bf3954945bbc03b15749331b519ba801fde6e91f` |
| Ensemble replicate | Frozen jobs all have `ensemble_replicate=0`, `condition_type=rna_only`, blank DMI mouse ID, `expert_variant=rna_only_primary` | Same job registry; `src/jobs.py:348–357`, `src/cache_identity.py:320–332` |
| Stored whole vectors | **20 per ensemble**, sample indices 0–19; **32×20=640 total** | Candidate identity gates in PL1 `:182–202`; Stage-11 production QC SHA-256 `59dc08dfbd8d0193dc017f1bb1ffb904dd9d45bd9dd9ef93b37d14a43a8095ed` has 32 ensemble + 32 projection records, each 20 rows |
| Mouse-weight strata | 10 DMI mice × 4 methods × 4 contexts = **160**, each reusing 20 candidate vectors; 3,200 weight rows/setting | PL1 weight gate `:203–216`, A20 qualification; these are not new reconstructions |
| Upstream stochastic sampling | Per-job seeds generate the stored vectors; 32 distinct seeds, configured thinning 100 and process count 1 | Stage-11 jobs/config; `src/ensembles.py:81–102` uses OptGP for GIMME/iMAT/CORDA; `:106–148` uses native RIPTiDe GapSplit |
| Truth selections | q10/q50/q90 are deterministic weighted quantile selections ordered by cached HEX1, with stable identity tie breakers | `scripts/dmi_bridge_pl2a_prepare_v1.py:224–279`; no random truth-sampling draw axis |
| Truth aliases/pairs | 160×3 = **480** condition selection aliases; 400×3 = **1,200** same-quantile pair aliases; **885** distinct pairs, **836** primary evaluable matched pairs | PL2A pair registry; A22 `:114–122,398`; alias collapse/holdout is not independent biological replication |
| Holdout support | Remove selected whole vector from each own condition pool: 19×19 = **361** contrasts for an evaluable held-out pair | A22 `:145–185`; 20×20 geometry support and 19×19 utility support are separate stages |
| Random-sign null | Exact expectation at cue correctness q=0.5, not Monte Carlo draws | `scripts/dmi_bridge_a22_dual_anchor_sign_only_core_v1.py:150`; PL2 scalar contract |
| Bootstrap | Reaction-level computational stability resampling; PL2D records 5,000 replicates | `scripts/dmi_bridge_pl2d_confirmation_v1.py:449–455`; not an evaluation identity factor |

Twenty retained vectors are neither ESS 20 nor evidence of convergence. Global
weight calibration targets ESS 20 over 320 tumor-matched candidates. Conditional
and joint ESS arise after restriction/product weighting. This audit makes no claim
about a common generated-draw count or burn-in across different upstream samplers.
No sampling was run to establish these metadata findings.

## Manuscript-ready terminology

“The frozen candidate panel comprised 32 method-by-tumor-by-RNA-context ensembles,
each with 20 stored whole-flux vectors (640 vectors total). For each of the 16
method-by-RNA-context settings, we compared every CT2A DMI mouse with every GL261
DMI mouse, using five subjects per tumor. The 25 mouse comparisons per setting
yielded 400 computational evaluations. Each evaluation reused the fixed 20-vector
supports with mouse-specific weights; the evaluations share animals and candidate
vectors and are not independent biological replicates. Weighted q10/q50/q90 truth
selections are deterministic benchmark selections, not additional random draws.”

The current main draft paragraph 28 already states the correct design. The earlier
“two mice/five draws” description should be replaced wherever used; it is not a
new scientific dimension missing from the frozen pipeline.

## Edits and evidence boundary

README and publication provenance now state the verified full Cartesian design
and absence of an independent reconstruction-draw factor. The migration report
marks S2 resolved; the historical migration metadata remain preserved. No code,
caption, frozen scientific artifact or manuscript DOCX was modified. Full source
identities, registry census and paragraph locators are in
[`S1_S2_EVIDENCE.json`](../manifests/provenance/S1_S2_EVIDENCE.json).
