# Supplementary Figure S5 — Directional-cue utility controls

## Scientific objective

S5 expands the utility analysis by showing how correct- and wrong-direction cues behave across reconstruction methods and RNA contexts, and by reporting the frozen utility outcome and tie structure at the wrong- and correct-cue endpoints. It complements Figure 3B–3D with context stratification and explicit outcome coverage; it does not reproduce those panels unchanged.

## Authoritative inputs

- `outputs/dmi_bridge_pl2b_sign_only_utility_v1/`: frozen A1 sign-only utility production (`BRIDGEPL2B_SIGN_ONLY_UTILITY_COMPLETE`); manifest SHA-256 `7d9ca97a3f4ea5aa8b0b983e09ee83fc3c82f326d1b1cbf9c6137cfadf25539d`.
- `outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/`: frozen matched A2 utility production (`BRIDGE_A22_SIGN_ONLY_UTILITY_FROZEN`); manifest SHA-256 `15b237151005538e5077a352868f873ef00466f4f5ea49307b781a8a8cd0e7cb`.
- `figures/fig3/tables/fig3_gain_cue_direction.parquet.xz` and its Figure 3 metadata: frozen correct/wrong cue gains on the exact 836 truth-pair, 4,179-reaction common population. The S5 builder verifies the table hashes and source manifest hashes recorded in the metadata.
- `outputs/dmi_bridge_pl1_predictability_landscape_v1/BRIDGEPL1_EVALUATION_REGISTRY.tsv`: verified against the PL2B fingerprint; maps evaluation IDs to the four current methods and RNA contexts.
- `outputs/dmi_bridge_a23_a1_a2_synthesis_v1/BRIDGEA23_MANIFEST.json`: frozen matched-population provenance, SHA-256 `b8e57074f3df23fe7708145bc01f51a19f946545b1c16e61628c2678dbcabd18`.

The compact Figure 3 tables and source context summaries remain authoritative. S5 performs deterministic aggregation and plotting only; it does not rerun utility or production calculations.

## Build and panels

From the repository root:

```bash
python figures/supp_fig5/build_supp_fig5_tables.py
jupyter nbconvert --to notebook --execute --inplace figures/supp_fig5/supp_fig5.ipynb
```

- **Panel A — Method- and context-stratified gain:** for each anchor × method × RNA context × cue direction, plot the arithmetic mean across reaction-evaluation gain means, with descriptive 5th–95th percentiles. These ranges are not confidence intervals. Each method facet has its own x-axis scale, which is disclosed in the plot.
- **Panel B — Overall utility and tie structure:** show improved, gain-tied, and harmed fractions among non-tie reaction × truth cases at λ = 0.25 for q = 0 (wrong cue) and q = 1 (correct cue). Truth ties are excluded from these bars and reported separately below them.

A1 is 1-Strong-Anchor (glucose uptake); A2 is 2-Strong-Anchors (glucose uptake + LDH_L). The Figure 3 paired gain tables use 836 matched truth pairs and a shared 4,179-reaction universe. PL2B's full A1 utility summary has 4,180 reactions because it includes the additional LDH_L reaction; its counts are retained only for the separate source-composition display and are not described as a paired A1/A2 population.

## Generated files

- `build_supp_fig5_tables.py` — checksum-gated deterministic table builder.
- `data/supp_fig5_context_gains.tsv` — 64 method/context/anchor/cue summaries.
- `data/supp_fig5_outcome_composition.tsv` — four frozen q-endpoint rows.
- `data/supp_fig5_truth_tie_coverage.tsv` — separate source-specific truth-tie counts.
- `data/supp_fig5_build_manifest.json` — input/output hashes, dimensions, and validation summary.
- `supp_fig5.ipynb` — plotting and visual checks.
- `outputs/supp_fig5.svg`, `outputs/supp_fig5.png` — editable vector figure and 300 dpi preview.

## Principal observations

Across the 836-pair common population, the frozen truth-level mean correct-cue gains are 0.014234 for A1 and 0.022802 for A2; wrong-cue gains are −0.018859 and −0.031733, respectively. On the primary source summaries, moving from q = 0 to q = 1 increases the improved fraction (A1 11.23% to 16.95%; A2 18.29% to 25.02%) and reduces the harmed fraction (A1 17.32% to 11.45%; A2 25.61% to 18.73%). Gain-tie fractions remain about 71.5% for A1 and 56.2% for A2.

## Limitations

Panel A's observational unit is a reaction × evaluation after aggregation over valid matched synthetic truths; these observations are descriptive and not independent replicates. Panel B uses source-specific full production denominators, which differ by the additional A1 LDH_L reaction. Method-specific x-axis scales differ because the frozen gain ranges differ substantially across methods.
