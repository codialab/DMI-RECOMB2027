# RECOMB 2027 DMI reproducibility repository

This study examines when weak qualitative directional information can further
reduce uncertainty after strong quantitative flux constraints, using residual
flux distributions in genome-scale metabolic models.

**Release status: NOT YET REPRODUCIBLE.** Production code and its provenance have
been curated. Frozen numerical inputs are withheld pending author redistribution
review, and several figure lineages require author reconciliation. The current
state is documented in [the migration report](docs/MIGRATION_REPORT.md).

## Analysis settings

- **A1:** HEX1 is the single strong anchor.
- **A2:** HEX1 + LDH_L are the two strong anchors.
- A separate weak directional cue has primary strength **λ = 0.25**. Frozen
  manuscript sensitivities use λ = 0, 0.5 and 1.0 where retained.

The production operators use dimensionless rank-based soft weighting. The DMI
measurements do not establish absolute flux equalities or physical-unit conversion.
Gains are in native GEM/model flux units. An operator-lineage issue remains: A20
ranks glucose uptake `max(-EX_glc__D_e, 0)` as its candidate Vmax observable while
the manuscript/contracts describe the glucose anchor through HEX1. Their intended
equivalence requires author verification; the migration preserves the operator.
The historical label `A2-L` denotes A2
in internal filenames and frozen table schemas; those schemas are preserved.

The frozen design uses iMAT, GIMME, CORDA and RIPTiDe; four RNA contexts; two tumor
contexts; five DMI subjects per tumor; and a full 5 × 5 mouse-comparison grid.
An individual evaluation compares one CT2A mouse with one GL261 mouse: 4 methods ×
4 contexts × 25 mouse pairs = **400 evaluations**. This differs from describing
the five subjects as independent reconstruction draws. There are 32 reconstruction
ensembles with 20 stored vectors each, and 4,179 common non-anchor reaction
coordinates. Finite candidate/support counts are not biological replication.

The final utility comparison uses **836 matched truth pairs**, with 1,440,313
distinct non-tie cases and 660,237 eligible reaction–evaluation pairs per anchor.
Cases are averaged within pairs before pairs receive equal weight. Original A1
held-out confirmation is distinct from development; A2 uses the same split as
post-freeze robustness and is not independent confirmation.

## Repository structure

| Location | Contents |
|---|---|
| `scripts/` | Original Bridge geometry, directional update, evaluation and evidence modules; paths retained for import/hash compatibility |
| `11_recomb_flux_expert_generation/src/` | Original reconstruction and candidate-vector modules; audit/upstream code, not a turnkey production release |
| `analyses/` | A1/A2 workflow descriptions and scientific boundaries |
| `data/` | Data policy; prospective frozen inputs are listed but currently withheld |
| `figures/main/` | Selected source plotting cells and upstream table builders |
| `figures/supplementary/` | S1 support-diagnostic code; manuscript discrepancies are recorded |
| `manifests/` | File mappings, source identities, checksums, rights decisions and figure registry |
| `docs/` | Reproducibility, panel map, provenance, author-review items and migration audit |
| `environment/` | Minimal plotting dependencies and separate optional upstream dependencies |
| `tests/` | Existing synthetic scientific-core regression checks |

The native code hierarchy is intentionally retained instead of reorganizing
scientific modules into new subsystem packages. Empty template directories are
not created, and code is not duplicated to satisfy a layout.

## Reproducing figures

Install the [plotting environment](environment/README.md), then run from this root:

```bash
python reproduce.py status
python reproduce.py checksums
python reproduce.py validate
python reproduce.py figures --only fig2DG,fig3BC,fig4B
python reproduce.py tables
```

`validate`, `figures` and `tables` currently report pending numerical inputs and
return a nonzero status. These commands are concrete prepared entry points, not
a claim that the absent inputs have been released. Once approved inputs are
included and their manifest decisions updated, the selected plotting cells can
regenerate quantitative panels without reconstruction, optimization or sampling.

Outputs go under ignored `reproduced/` directories. Existing outputs are reused
only when input/code/environment identities and output hashes agree; incompatible
outputs are not overwritten. Each plotting job records logs and a status manifest.
All persistent writes, temporary work and plot caches are directed into this
repository by the entry point. The source repository is never needed by it.

Full main/supplementary artwork is not yet reproducible. Conceptual panels are
marked explicitly, Figure 3D's existing artwork is superseded, and supplementary
numbering/composition conflicts remain. See [ANALYSIS_MAP.md](docs/ANALYSIS_MAP.md).

## Upstream analysis and data availability

The reproduction hierarchy is external inputs → production computation → frozen
processed tables → figure inputs → figures and manuscript statistics. Frozen
figure inputs should make expensive upstream regeneration unnecessary.

The preserved upstream modules require external datasets, the parent iMM1865 model,
frozen medium/projection inputs, and historical identity artifacts. Some foundation
patch files are missing from the source checkout. Their gates have not been
bypassed; a complete upstream rerun is not currently supported or scientifically
revalidated. No expensive computation was performed during migration.

[data/README.md](data/README.md) explains omitted datasets, archive boundaries,
compressed-artifact policy and redistribution review. Original frozen inputs are
identified by source SHA-256 even when withheld.

## Citation and license

Use [CITATION.cff](CITATION.cff). Final title, DOI, ORCIDs and conference citation
remain TODOs; author names are taken from the current manuscript draft.

Source code is licensed under [MIT](LICENSE). Cleared original research products,
figures and documentation are licensed under [CC BY 4.0](LICENSE-DATA).
Third-party datasets, metabolic models, software, reconstruction methods and
externally sourced materials retain their original licenses and are not relicensed
by this repository. Ambiguous artifacts are withheld pending author review.
