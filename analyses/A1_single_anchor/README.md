# A1 single-anchor workflow

A1 uses HEX1 as the single strong anchor through the frozen rank-based weighting
operator. Downstream geometry describes existing candidate vectors; it does not
generate new flux samples. PL1 computes geometry, PL2A fixes reaction/truth registries,
PL2B evaluates weak sign-only cues, PL2C examines development reactions, and the
PL2D-H freeze/amendment records precede the original PL2D held-out confirmation.

Implementation modules are in `scripts/dmi_bridge_pl*`. Frozen identities are in
`manifests/provenance/FROZEN_LINEAGE.json`; original contracts are preserved under
`manifests/provenance/source/docs/`. These identities document provenance without
implying that all upstream datasets are included.

The draft's development usefulness scatter has 2,454 finite reactions. The original
development snapshot reports 2,455 for its own aggregation/support definition.
These are distinct frozen populations; do not silently replace one denominator
with the other. The one-shot confirmation plot uses 838 finite reactions.

Truth construction includes held-out candidate exclusion, q10/q50/q90 selection
and alias handling as fixed by PL2A. The final A1/A2 comparison uses the 836-pair
intersection specified by A23. The source intermediate split originally contains
4,180 reactions; the common A1/A2 weak-cue population excludes both anchors and
contains 4,179. Preserve both definitions and label their uses.

Run `python reproduce.py status` for frozen-input readiness. Do not run original
upstream adapters expecting absent artifacts or bypass their checksum gates.
