# A2 dual-anchor workflow

A2 uses HEX1 + LDH_L as strong anchors. It is the historical internal A2-L setting.
The lactate-producing coordinate is `max(-v_LDH_L, 0)` before ranking, according
to the GEM reaction orientation. No absolute DMI/model unit equality is assumed.

A20 qualifies/fixes the dual-anchor weighting, A21 computes geometry from existing
vectors, A22 evaluates the separate weak sign-only cue, and A23 fixes the matched
A1/A2 synthesis population. Implementation modules are `scripts/dmi_bridge_a2*`.

Primary λ is 0.25; only frozen manuscript-retained λ = 0, 0.5 and 1.0 sensitivities
are relevant. Alternative secondary anchors are excluded. A2 is post-freeze
robustness on the original reaction split, not independent held-out confirmation.

Frozen A1/A2 geometry pairs cover 400 evaluations × 4,179 reactions = 1,671,600
rows. Matched utility uses 836 truth pairs × 4,179 reactions = 3,493,644 cases,
including truth ties. The non-tie cases and equal-weight eligible reaction–evaluation
pairs used in Figure 3 are separate denominators.

Historical table schemas retain their stored labels and hashes. User-facing
figures/documentation use A2. Redistribution and final-figure readiness are recorded
in `manifests/` and `docs/MIGRATION_REPORT.md`.
