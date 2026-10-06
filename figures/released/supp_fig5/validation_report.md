# Supplementary Figure S5 validation report

## Inputs and provenance

| Input | Use | Status / checksum |
|---|---|---|
| `outputs/dmi_bridge_pl2b_sign_only_utility_v1/BRIDGEPL2B_MANIFEST.json` and `BRIDGEPL2B_CONTEXT_SUMMARY.tsv` | A1 frozen source, endpoint composition, truth-tie counts | `PL2B_SIGN_ONLY_UTILITY_COMPLETE`; manifest SHA-256 `7d9ca97a3f4ea5aa8b0b983e09ee83fc3c82f326d1b1cbf9c6137cfadf25539d`; context table matches manifest artifact checksum |
| `outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_MANIFEST.json` and `BRIDGEA22_CONTEXT_SUMMARY.tsv` | A2 frozen matched source, endpoint composition, truth-tie counts | `BRIDGE_A22_SIGN_ONLY_UTILITY_FROZEN`; manifest SHA-256 `15b237151005538e5077a352868f873ef00466f4f5ea49307b781a8a8cd0e7cb`; context table matches manifest artifact checksum |
| `figures/fig3/tables/fig3_gain_cue_direction.parquet.xz` | Matched correct/wrong gain distributions | SHA-256 `3c571ad4b97bb451d165c00fd464ccb51b0b5de8c8e69a49d3b423d78eeee30d`; matches Figure 3 metadata |
| `figures/fig3/tables/fig3_gain_distribution_metadata.json` | Gain definitions, population and aggregation provenance | SHA-256 `d6ebfaa08aec0255104463bd327c55a74f24a1d074589d563782cb5953ab14a7` |
| `outputs/dmi_bridge_pl1_predictability_landscape_v1/BRIDGEPL1_EVALUATION_REGISTRY.tsv` | Method and context labels for Panel A | SHA-256 `aa468dcec11dea764682bd3d7eb14c03c003063e74d92080d0bfaabb366290a3`; matches PL2B input fingerprint |
| `outputs/dmi_bridge_a23_a1_a2_synthesis_v1/BRIDGEA23_MANIFEST.json` | Fixed A1/A2 matched population provenance | SHA-256 `b8e57074f3df23fe7708145bc01f51a19f946545b1c16e61628c2678dbcabd18` |

The builder also verifies the hashes of all five compact Figure 3 gain tables against both its frozen checks and Figure 3 metadata. No exploratory, alternative-anchor, or superseded utility table enters the panels.

## Dimensions and aggregation

- Figure 3 cue input: **2,640,948** rows, comprising 660,237 reaction-evaluation observations × two cue directions × two anchors.
- Panel A input per anchor: **660,237** reaction-evaluation cue rows per direction; **370** represented evaluation IDs and **3,625** represented reaction IDs. The shared eligible truth-pair population is **836**; the common reaction universe is **4,179**.
- Panel A output: **64** unique summaries, keyed by anchor × method × RNA context × cue direction (four methods × four contexts × two anchors × two cue directions). Each summary reports its exact reaction-evaluation count, arithmetic mean, median, and descriptive 5th/95th percentiles.
- Panel B output: **4** rows: two cue endpoints for each source anchor at λ = 0.25. Fractions use the source's declared non-tie reaction × truth case denominator.
- Truth-tie coverage output: **2** rows, one per anchor source.
- All plotted gains are finite; output tables have unique declared keys. No additional rows were excluded by the S5 builder.

## Population and exclusions

The matched A1/A2 gain tables contain **3,493,644** cases per arm before truth-tie exclusion: 836 truth pairs × 4,179 reactions. Each arm has **2,053,331** truth ties and **1,440,313** valid non-tie cases for the matched gain table. For A1, the Figure 3 builder filters the extra LDH_L target reaction: **836** cases, one per matched truth pair. A2's 4,179-reaction universe already excludes LDH_L as an anchor reaction.

The full-source composition summaries retain their source denominators. PL2B A1 contains **3,494,480** cases over 4,180 reactions, with **1,441,149** non-ties and **2,053,331** truth ties. A22 A2 contains **3,493,644** cases over 4,179 reactions, with **1,440,313** non-ties and **2,053,331** truth ties. Panel B labels the separate denominators and does not present these rows as a matched A1/A2 contrast.

No nonfinite gain values, duplicate reaction-evaluation/cue keys, missing evaluation mappings, or invalid categories were found. Gain-tie outcomes are distinct from truth ties: improved + gain-tied + harmed counts sum to the non-tie denominator; truth ties are excluded and reported separately.

## Definition and quantitative checks

- Figure 3 metadata defines `correct_gain` and `wrong_gain` and gives the production aggregation: average valid matched synthetic-truth gains within reaction × evaluation. The production identity for `g_info` and the reaction-evaluation directional advantage is preserved in the source metadata; S5 does not recompute either definition.
- Independent checks on the stored Figure 3 reaction-evaluation table found maximum absolute identity errors of **3.55 × 10⁻¹⁵** for `g_info_mean = 0.5 × (correct_gain_mean − wrong_gain_mean)` and **0** for `directional_advantage = correct_gain_mean − wrong_gain_mean`.
- The source utility contracts define gain classification using a tolerance of (10^{-12}), with negative/positive outcomes classified as harmed/improved and values within tolerance as gain-tied.
- At λ = 0.25, the selected endpoints are q = 0 (wrong cue) and q = 1 (correct cue); both are present exactly once per source.
- Matched truth-level mean gains (mmol gDW⁻¹ h⁻¹): A1 correct **0.014234**, wrong **−0.018859**, difference **0.033094**; A2 correct **0.022802**, wrong **−0.031733**, difference **0.054535**. These are reference summaries, not Panel A's reaction-evaluation aggregation.
- Panel B outcomes: A1 improved/tied/harmed changes from **11.23% / 71.45% / 17.32%** at q = 0 to **16.95% / 71.60% / 11.45%** at q = 1. A2 changes from **18.29% / 56.10% / 25.61%** to **25.02% / 56.25% / 18.73%**.
- Truth-tie coverage is **2,053,331 / 3,494,480 (58.76%)** for A1 and **2,053,331 / 3,493,644 (58.77%)** for A2.

## Render and reproducibility checks

- Table builder: `python figures/supp_fig5/build_supp_fig5_tables.py` — **PASS**, source and derived-table hash gates and row checks succeeded.
- Notebook: executed with `jupyter nbconvert --to notebook --execute --inplace figures/supp_fig5/supp_fig5.ipynb` — **PASS**. Kernel launch required host-side execution because the sandbox denied Jupyter's local kernel socket; the approved invocation read S5 tables and wrote only S5 artifacts.
- Outputs: `outputs/supp_fig5.png` (300 dpi) and `outputs/supp_fig5.svg`; both inspected. Panel spacing, labels, scales, legends, and denominator annotation are legible.
- Production reconstruction, optimization, sampling, candidate generation, and simulation pipelines were not run.

## Issues

- **BLOCKING:** None.
- **IMPORTANT:** Panel A's four method facets use method-specific x-axis scales because the frozen gain magnitudes span different ranges; the plot explicitly discloses this. The source-specific A1/A2 denominators in Panel B differ by the additional A1 LDH_L reaction and are labeled.
- **COSMETIC:** None identified.
