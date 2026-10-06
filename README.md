# RECOMB 2027 DMI reproducibility repository

This study examines when weak qualitative directional information can further
reduce uncertainty after strong quantitative flux constraints, using residual
flux distributions in genome-scale metabolic models.

**Release status: frozen-input reproduction is available for Figure 2D–G, Figure
3B–C, and Figure 4B.** The current source is RECOMB 2027 DMI commit
`225279129`; included compact inputs have source SHA-256 records and pass the
population and endpoint checks. Other figure assemblies and supplement lineages
remain explicitly marked in the [reproducibility matrix](docs/REPRODUCIBILITY_MATRIX.md).
The migration history and remaining scientific decisions are in the
[migration report](docs/MIGRATION_REPORT.md).

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
The [anchor operator audit](docs/ANCHOR_OPERATOR_AUDIT.md) records the traced
distinction and required author decision; the labels above remain under review.
The historical label `A2-L` denotes A2
in internal filenames and frozen table schemas; those schemas are preserved.

The frozen design uses iMAT, GIMME, CORDA and RIPTiDe; four RNA contexts; two tumor
contexts; five DMI subjects per tumor; and a full 5 × 5 mouse-comparison grid.
An individual evaluation compares one CT2A mouse with one GL261 mouse: 4 methods ×
4 contexts × 5 CT2A mice × 5 GL261 mice = **400 evaluations**. Each of the 16
method/RNA contexts contains all 25 mouse pairs exactly once. The two sets of five
are biological DMI subjects; there is no independent five-draw factor. These
comparisons share animals and vectors and are not independent biological replicates.
See the [evaluation design audit](docs/EVALUATION_DESIGN_AUDIT.md). There are 32 reconstruction
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
| `data/` | Frozen plotting inputs, derived summaries, and data-boundary policy |
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

`validate`, `figures`, and `tables` run from included files only. `validate`
checks frozen panel populations and key endpoints. The plotting command renders
Figure 2D–G, Figure 3B–C, Figure 4B, and the audited supplementary Figures S1–S3; `tables` exports Figure 2 correlation
statistics and frozen pair-weighted Figure 3 summaries. See
[`tables/manuscript/`](tables/manuscript/) for a checked-in machine-readable
snapshot and generation notes. Plot jobs record code, input, and environment
identities under ignored `reproduced/` output directories.

Outputs go under ignored `reproduced/` directories. Existing outputs are reused
only when input/code/environment identities and output hashes agree; incompatible
outputs are not overwritten. Each plotting job records logs and a status manifest.
All persistent writes, temporary work and plot caches are directed into this
repository by the entry point. The source repository is never needed by it.

The numeric panels above are reproducible; this does not establish that each
saved composite is the final manuscript assembly. Current source artwork and
editable components are retained under [`figures/released/`](figures/released/)
with identity and status recorded in `manifests/ARTWORK.json`. Figure 3D is
provenance-only because the case-pooled artwork is superseded. Supplementary
figures S1–S3 have validated standalone reproduction paths. S4–S7 have current
source artwork retained, but lack a validated release reproduction path or have unresolved manuscript mapping; see the
[panel matrix](docs/REPRODUCIBILITY_MATRIX.md).

## Upstream analysis and data availability

The reproduction hierarchy is external inputs → production computation → frozen
processed tables → figure inputs → figures and manuscript statistics. Frozen
figure inputs should make expensive upstream regeneration unnecessary.

The preserved upstream modules require external datasets, the parent iMM1865 model,
frozen medium/projection inputs, and historical identity artifacts. Some foundation
patch files are missing from the source checkout. Their gates have not been
bypassed; a complete upstream rerun is not currently supported or scientifically
revalidated. No expensive computation was performed during migration.

[data/README.md](data/README.md) explains omitted upstream datasets, archive
boundaries, compressed-artifact policy and included derived inputs. No full raw
datasets, parent model, credentials, or bulk upstream caches are included.

## Citation and license

Use [CITATION.cff](CITATION.cff). Final title, DOI, ORCIDs and conference citation
remain TODOs; author names are taken from the current manuscript draft.

Source code is licensed under [MIT](LICENSE). Cleared original research products,
figures and documentation are licensed under [CC BY 4.0](LICENSE-DATA).
Third-party datasets, metabolic models, software, reconstruction methods and
externally sourced materials retain their original licenses and are not relicensed
by this repository. Ambiguous artifacts are withheld pending author review.
