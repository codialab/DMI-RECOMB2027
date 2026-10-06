# Privileged reaction exclusion audit

Date: 2026-10-06. Source: `/home/pty/work/project_MRI_brain_tumor`, commit
`16379a169f45ec556274953eb66ae12293e1e985`. Starting destination:
`ce621c74ed4a5a325ebb7521d1247d9c2791b94f`, branch `publication-migration`.

**Verdict: CURRENT POPULATION DEFENSIBLE**, specifically with respect to the
three coordinates investigated here. The actual common population already omits
HEX1, glucose exchange and LDH_L. No production population change is needed to
implement policy C on the frozen axis. This does not settle the separate S1
anchor-operator/manuscript issue or authorize a publication release.

## 1. Finding: glucose exchange is not a weak-cue reaction

`EX_glc__D_e` is **absent from the 4,181-coordinate frozen flux cache**. It is
therefore absent from PL1, PL2, A21/A22/A23 and the final quantitative Figure
2/3/4 reaction tables. Subtracting it from the current 4,179-reaction population
removes zero reactions; the proposed count of 4,178 does not describe this cache.

The raw candidate panel used for weighting has a glucose-exchange column even
though the downstream analysis cache does not. These are different objects:

| Role | Coordinate / observable | Production use |
| --- | --- | --- |
| Glucose conditioning | Raw `max(-EX_glc__D_e, 0)` | Vmax-to-candidate rank-distance weights in A1 and A2 |
| Additional conditioning in A2 | Raw `max(-LDH_L, 0)` | Additional Vlac-to-candidate rank-distance term |
| Deterministic truth definition | Cached HEX1 | Orders candidates for weighted q10/q50/q90 truth selection |
| Common coupling descriptor | Cached HEX1 contrast | Defines HEX1-to-B contrast correlation under the respective frozen weights |
| Ordinary weak cue | B in the admitted common reaction registry | Its truth contrast sets the correct/wrong sign cue; its residual distribution is evaluated |

See [the anchor operator audit](ANCHOR_OPERATOR_AUDIT.md) for the complete
weighting lineage and the absence of an established EX-to-HEX1 equivalence.
No biological or mathematical equivalence is assumed here. The raw uptake
observable and the normalized cached flux coordinate are not interchangeable.

## 2. Where glucose exchange disappeared

The source cache reconstruction's `recover_rxns` function is at
`outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstruct_flux_cache.py:134`.
Its actual selection rule, at lines 143–168, is:

1. Load the 4,747 ordered Stage 4 reaction identities from
   `dmi_sa_r1_cup_reuse_20260926/selected_model/ensembles/stage4_projection_4747/metadata.json`.
2. Collect the 68 subsystem identities in
   `13_recomb_identifiability/exports/dmi_sa_offline_feasibility/validation/rna_effects_used.csv`.
3. Use the historical runner's
   `13_recomb_identifiability/exports/dmi_sa_offline_feasibility/reaction_metadata/reaction_metadata.tsv`
   to retain reactions annotated to one of those subsystems, preserving Stage 4 order.

Stage 4 **does** contain `EX_glc__D_e`. The runner's annotation assigns it only
to `Exchange/demand reaction`, which is not one of the 68 RNA subsystems. Thus it
does not enter the 4,181-coordinate cache. The audit independently recovered the
entire ordered axis from these identity/annotation inputs and obtained an exact
match to the stored `rxns` array. Only identities and annotations were needed;
no expression values or cached flux vectors were loaded for this recovery.

The reconstructor explicitly distinguishes this runner-faithful rule from the
alternative `subsystem_memberships.csv.xz` rule (lines 170–203), which would
produce a different population. We did not substitute that alternative. This
audit does not establish that all exchange reactions are absent or justify
extending results to the full metabolic-model axis.

The cache SHA-256 is
`f9644da4ebd875cf5f0b34cc58b20faf0a0457d7c7a44635e8a010d8896c3dbb`.
All inspected artifact hashes and membership counts are in
[the evidence manifest](../manifests/provenance/PRIVILEGED_REACTION_EXCLUSION_EVIDENCE.json).

## 3. Downstream lineage and direct membership checks

Paths in the artifact column are relative to the source repository. Code links
point to the preserved destination copies; line references below identify the
production decisions, not code executed by this audit.

| Stage | Production decision / artifact | Direct audit result |
| --- | --- | --- |
| PL1 | [Resume adapter](../scripts/dmi_bridge_pl1_predictability_resume_v1.py), lines 399–400, excludes HEX1 from `rxns`; `outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_REACTION_SUMMARY.tsv` | 4,180 unique B identities; EX has no row |
| PL2A | [Preparation](../scripts/dmi_bridge_pl2a_prepare_v1.py), lines 338–344; `outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_REACTION_REGISTRY.tsv` | Same 4,180 identities; EX has no row |
| PL2B | Uses the frozen PL2A registry/truths; the preserved PL2B contract admits that population | No EX B can be produced by that admitted registry; case shards were not reprocessed |
| PL2C | [Split core](../scripts/dmi_bridge_pl2c_geometry_utility_core_v1.py), lines 59–80; `outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_REACTION_SPLIT.tsv` | 4,180 identities; EX has neither development nor confirmation assignment |
| PL2D | `outputs/dmi_bridge_pl2d_confirmation_v1/BRIDGEPL2D_CONFIRMATION_REACTION.tsv` | 1,045 held-out reaction rows; no EX row |
| A21 | [Geometry adapter](../scripts/dmi_bridge_a21_dual_anchor_geometry_v1.py), lines 320–329, subtracts HEX1 and LDH_L; `outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_REACTION_SUMMARY.tsv.xz` | 4,179 identities, exactly cache minus those two; no EX row |
| A22 | [Common-target core](../scripts/dmi_bridge_a22_dual_anchor_sign_only_core_v1.py), lines 66–78; `outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_REACTION_REGISTRY.tsv` | 4,179 identities, exactly the A21 set; no EX row |
| A23 | [Synthesis](../scripts/dmi_bridge_a23_a1_a2_synthesis_v1.py), `identity`, lines 172–187, requires A22 registry = PL2C registry minus LDH_L | 3,134 development + 1,045 confirmation reactions; no EX contribution; frozen summary artifacts verified against A23 manifest |
| Figure 2 D–G | `figures/fig2/figure2_DG_points.csv.xz` | D/E/F/G have 2,454/2,485/838/843 finite plotted reactions respectively; zero EX rows |
| Figure 3 B/C | `figures/fig3/tables/fig3_gain_reaction_evaluation.parquet.xz`; `fig3_gain_cue_direction.parquet.xz` | 660,237 eligible evaluation–reaction rows per setting in the first table; zero EX rows in both tables |
| Figure 4 B | All four `figures/fig4/data/fig4_geometry_paired/algorithm=*/part.parquet.xz` partitions | Each has 100 × 4,179 rows; union has 400 × 4,179 = 1,671,600; no EX rows |
| Figure 4 C | `figures/fig4/data/fig4_panelC_candidate_groups.parquet.xz`; `fig4_panelC_candidate_metrics.parquet.xz` | 100 candidate example groups and 400 method-specific rows; no EX rows; policy C retains every example candidate |

The different finite Figure 2 counts are endpoint support, not different
admitted reaction universes. Aggregated Figure 3 summaries and A23 comparisons
inherit the populations above; an exchange endpoint is not hidden in a summary.
This audit did not scan every large PL1/PL2/A21/A22 case shard. Their exclusion
is established by the exact cache/registry identity and production iteration;
the final figure tables were checked directly.

## 4. Individual glucose-exchange results

| Requested result | A1 | A2 |
| --- | --- | --- |
| B geometry descriptors | Not evaluated; no B row | Not evaluated; no B row |
| B information-gain endpoints | Not evaluated; no B row | Not evaluated; no B row |
| B directional-usefulness endpoints | Not evaluated; no B row | Not evaluated; no B row |
| Development/confirmation assignment | None | None |
| Reaction-specific contribution to the frozen association/mapping | No included observations | No included observations |

These are **N/A endpoints**, not measured zero gain or zero usefulness.
The exchange observable affects the upstream weights of other B reactions;
its absence as B does not mean that it has no conditioning role. Adding an EX
endpoint would require a changed candidate axis and fresh scientific analysis,
which is outside this audit.

## 5. Why HEX1 is not an ordinary weak cue under this protocol

[PL2A `select_truth`](../scripts/dmi_bridge_pl2a_prepare_v1.py), lines 224–245,
sorts candidates by cached HEX1 (with candidate identity breaking ties), then
uses cumulative A1 baseline weights to choose q10/q50/q90 candidates. Selection
is deterministic. A22's contract in
[the dual-anchor core](../scripts/dmi_bridge_a22_dual_anchor_sign_only_core_v1.py),
line 140, requires reuse of those PL2A truth identities, with no A2 reselection.

If B were HEX1, the coordinate whose contrast is evaluated would also be the
coordinate used to select the synthetic truths. Candidate holdout removes the
selected vector from inference; it does not make that truth-selection rule
independent of HEX1. Consequently, a HEX1 endpoint could be calculated, but it
would be a privileged diagnostic, not an ordinary target under an outcome-blind
truth definition. Inclusion on equal terms would require an explicit author
decision about a scientifically different truth protocol or claim.

HEX1 also supplies the common contrast-coupling coordinate: the PL1 resume
adapter, lines 270–288, and A21 adapter, lines 149–162, pass HEX1 contrasts as
`anchor_ct2a/anchor_gl261`. A21 preserves that descriptor for A1/A2 comparability;
A23 verifies the preservation at lines 167 and 248. If B = HEX1, the descriptor
becomes self-correlation: mathematically one for a nondegenerate distribution,
and undefined under the implementation's degeneracy tolerance. See
[the scalar core](../scripts/dmi_bridge_pl1_predictability_core_v1.py),
`weighted_correlation`, lines 87–108, and its use at line 253. That is not an
ordinary between-coordinate coupling measurement. The audit did not calculate
HEX1 weak-cue results or alter this code.

## 6. Exclusion policies on the actual frozen axis

Let U be the stored 4,181-coordinate axis. It contains HEX1 and LDH_L and does
not contain EX_glc__D_e.

| Policy | Intended exclusion | Actual count | Consequence |
| --- | --- | ---: | --- |
| A: current common population | HEX1 + LDH_L | 4,179 | Excludes truth-selection/common-coupling coordinate and additional A2 conditioning coordinate; EX already absent upstream |
| B: conditioning coordinates only | EX_glc__D_e + LDH_L | **4,180** | Restores HEX1; does not address direct truth-definition dependency or self-coupling |
| C: all three privileged coordinates | HEX1 + EX_glc__D_e + LDH_L | **4,179** | Exactly the same set as A; no additional deletion |

As a conceptual policy, C most clearly distinguishes both conditioning
coordinates from the truth-selection/coupling coordinate and excludes all
three from ordinary B endpoints. On the frozen population it is already
implemented as a set, although the final subtraction only names two members.
Policy A is therefore defensible for this frozen lineage. Policy B would not
be a harmless label correction: it adds a privileged reaction. The asserted
4,179/4,179/4,178 comparison would only hold for an axis containing all three,
which is not the authoritative axis. There is no justification to delete a
different reaction to force a count of 4,178.

Exclusion alone does not establish that the strong weighting coordinate equals
HEX1. S1 still requires author review of literal HEX1 conditioning claims.

## 7. Lightweight sensitivity: C versus frozen common A population

The [audit calculation](../analyses/audit/privileged_reaction_exclusion.py)
reads authorized frozen inputs, verifies their identities, filters rows by C,
and writes only aggregate audit outputs. No frozen artifact was overwritten;
no weights, truths, splits, candidate states or scientific defaults changed.
No source production module was imported or run.

```bash
PYTHONDONTWRITEBYTECODE=1 python analyses/audit/privileged_reaction_exclusion.py \
  --source-root /home/pty/work/project_MRI_brain_tumor
```

This command requires explicit access to the read-only source inputs; those
rights-pending tables remain omitted from the publication repository. Use
`--replace-audit-outputs` only to replace this script's two destination audit
summaries on a deliberate rerun. It cannot replace source production artifacts.

[The sensitivity table](../tables/audit/privileged_reaction_exclusion_sensitivity.tsv)
reports aggregation units, source paths, retained denominators, signed and
absolute changes, and whether a statistic was recomputed from filtered rows or
preserved by the proven identity of its input set. Every reported policy C
absolute change is **0**, and every reported removed-row count is **0**.

### Headline comparisons

Values below are rounded for display; the TSV retains full computed precision.
The current and C values are identical.

| Quantity | A1 current = C | A2 current = C | Absolute policy C change |
| --- | ---: | ---: | ---: |
| Figure 2 development Spearman rho (D/E) | 0.857116 | 0.833944 | 0 |
| Figure 2 confirmation Spearman rho (F/G) | 0.835144 | 0.849366 | 0 |
| Figure 3 mean correct-cue gain | 0.010386903 | 0.022294330 | 0 |
| Figure 3 mean wrong-cue gain | -0.013746815 | -0.030611388 | 0 |
| Figure 3 mean correct-minus-wrong separation | 0.024133718 | 0.052905717 | 0 |
| Figure 3 mean information advantage | 0.012066859 | 0.026452859 | 0 |
| Figure 2 development mean usefulness fraction | 0.146601285 | 0.223616366 | 0 |
| Figure 2 confirmation mean usefulness fraction | 0.142007574 | 0.218401555 | 0 |
| Figure 4 mean direction entropy | 0.020117155 | 0.046696558 | 0 |
| Figure 4 mean dominant-direction mass | 0.989652161 | 0.976350807 | 0 |
| Figure 4 mean central 80% width | 0.616939113 | 0.992334026 | 0 |
| Joint ESS, equal mean over 400 evaluations | 5.094012412 | 27.186039153 | 0 |

**Aggregation distinctions:** Figure 3 first averages distinct non-tie truths
within each evaluation–reaction and then gives equal weight to the 660,237
eligible evaluation–reaction pairs per setting. Figure 2 usefulness is the
truth-level predicate `correct_gain > 1e-12 AND information_advantage > 1e-12`,
averaged within evaluation and then within reaction; the displayed audit mean
averages the finite plotted reaction means. It is not inferred by thresholding
Figure 3 pair-mean gains. A23 case-pooled means use a different unit and are
identified separately in the TSV; they are not substituted for Figure 3 means.
ESS is summarized once per evaluation, not as 1,671,600 independent values.

The table also covers Figure 3 medians, selected quantiles and positive-fraction
ECDF diagnostics; six Figure 4 descriptors with paired means/medians; all eight
A23 geometry-comparison features; A23 case-pooled controls/separation and its
development/confirmation/context associations. Figure 3/4 distribution arrays
are unchanged, so every ECDF threshold and plotting downsample is unchanged,
not just the reported moments. Figure 4's example candidate set and rankings
remain unchanged without reranking or regenerating the examples. Conceptual
panels have no numerical reaction aggregate to refilter.

**Older A1 development result:** M1's original PL2 population contains 4,180
reactions and includes LDH_L in development. Its 2,455 finite-reaction rho is
`0.8570601823118907`. It is not the frozen 4,179-reaction baseline for this
comparison. The final matched Figure 2 D baseline already excludes LDH_L and
has 2,454 finite reactions and rho approximately `0.8571160663367106`.
That historical difference is not a newly induced policy C effect. The older
held-out population contains no LDH_L and remains identical; M1 confirmation
rho, 16/16 positive contexts and existing bootstrap summary are preserved by
identity of their inputs. The bootstrap was not rerun. We do not claim that
all older 4,180-reaction development aggregates equal the matched aggregates.

Figure 2 recomputed rho agrees with the frozen reference within the existing
exporter's 1e-6 parity tolerance for last-bit/rank-tie effects. The comparison
of current versus C is exact because the arrays are identical.

## 8. Manuscript consequences and scope of the verdict

No interpretation of the frozen reaction aggregates changes under C: the
positive development/confirmation associations, correct-versus-wrong cue
separation, and A1/A2 geometry differences are unchanged. There is no EX-specific
weak-cue contamination to remove from these tables. This finding does not
validate synthetic truths as biological truth or resolve other method issues.

Author review should distinguish the following prospective statements. This
pass changes none of the manuscript text, existing documentation definitions,
final figure labels or production files.

| Affected statement/location | Supported clarification | Scientific impact |
| --- | --- | --- |
| Reaction-population Methods; Supplementary Methods S3/S4/S7; Table S5 and eight-descriptor descriptions | The 4,181-axis is a historical RNA-subsystem-filtered cache; the 4,179 common B set omits HEX1/LDH_L and already lacks EX_glc__D_e | Clarifies analyzed universe; no numeric change |
| Prior `AUTHOR_CORRECTIONS_REQUIRED.md`, S1 HEX1 diagnostic-role row | Its open question about excluding exchange is answered for this cache: exchange is already absent upstream | Follow-up resolved here; preserve that historical audit unchanged |
| README/analysis introductions; Figure 2/3/4 and supplementary anchor labels | Keep the S1 distinction between actual glucose-exchange weighting and HEX1 truth/coupling roles | S1 operator/interpretation decision remains open |
| Any proposed “policy C leaves 4,178” or “policy B leaves 4,179” sentence | For the authoritative axis C leaves 4,179 and B leaves 4,180 | Corrects arithmetic/premise; does not redefine science |
| Any prospective ordinary HEX1 weak-cue endpoint | HEX1 truth-dependent and self-coupling endpoints require separate interpretation or a revised truth protocol | Scientific author decision if inclusion is desired |
| M1 versus final Figure 2 development rho | Distinguish original 4,180/2,455-finite and matched 4,179/2,454-finite populations | Prevents aggregation/population conflation |

Remaining follow-up: author-approved S1 wording; clear disclosure of the
historical reaction-axis restriction; redistribution approval for frozen data.
No new scientific operator is recommended or implemented. The publication
release remains subject to its existing migration verdict and unresolved issues.

## 9. Validation

The evidence manifest records identical source HEAD, empty tracked-tree/index
diff and identical index-listing hash before/after the calculation. Frozen
Figure 2/3/4 inputs and cache were verified against their recorded SHA-256
identities; A23 summary hashes were checked against its admitted manifest.
The aggregate comparison asserts identical retained arrays/sets and zero C
changes. Additional checks, including Python syntax, restricted destination
diff scope and `git diff --check`, are recorded in
`manifests/provenance/PRIVILEGED_REACTION_EXCLUSION_VALIDATION.json`.

All persistent writes are destination audit material and its checksum manifest.
No rights-pending input table, source artifact or artwork was copied. The new
aggregate sensitivity table is an author-review audit product, not a finalized
publication data release or a resolution of third-party redistribution rights.
