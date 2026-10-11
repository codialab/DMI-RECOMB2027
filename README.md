# RECOMB 2027 DMI reproducibility repository

This study predicts when weak qualitative directional information is useful for
estimating flux magnitudes after strong quantitative flux anchoring, using
residual flux distributions in genome-scale metabolic models. The central
contribution is identifying geometries in which a weak directional cue helps
estimate magnitude; the results do not support a general claim that such cues
reduce flux-magnitude uncertainty.

**Release status: frozen-input reproduction is available for Figure 2D–G, Figure
3B–D, Figure 4B–C, and S1–S3 and S6.** The source analysis is RECOMB 2027 DMI
commit `22527912996cecedd17bb23c0b7ba93e889ad4ad`; included compact inputs have
source SHA-256 records and pass the population and endpoint checks. Other figure assemblies and supplement lineages
remain explicitly marked in the [reproducibility matrix](docs/REPRODUCIBILITY_MATRIX.md).
The migration history and remaining scientific decisions are in the
[migration report](docs/MIGRATION_REPORT.md).

## Analysis settings

- **A1:** glucose-uptake anchor, implemented as `max(-EX_glc__D_e, 0)`.
- **A2:** glucose uptake plus the `LDH_L` lactate proxy, implemented as
  `max(-EX_glc__D_e, 0)` and `max(-LDH_L, 0)`.
- **HEX1:** not a quantitative anchor; it is used for candidate-state ordering
  and geometry analysis.
- A separate weak directional cue has primary strength **λ = 0.25**. Frozen
  manuscript sensitivities use λ = 0, 0.5 and 1.0 where retained.

A1/A2 use measured Vmax and, for A2, Vlac ranks to softly weight candidate
states. The weak cue supplies a reaction's direction only; its usefulness for
flux-magnitude estimation is evaluated after the strong-anchor weighting. It
does not impose a flux magnitude or replace A1/A2.

The production operators use dimensionless rank-based soft weighting. The DMI
measurements do not establish absolute flux equalities or physical-unit conversion.
Gains are in native GEM/model flux units. The [anchor operator audit](docs/ANCHOR_OPERATOR_AUDIT.md)
traces the implemented observables and HEX1's diagnostic roles. The historical
label `A2-L` denotes A2 in internal filenames and frozen table schemas; those
schemas are preserved.

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
| `figures/main/` | Selected source plotting cells and standalone main-panel plotters |
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
python reproduce.py figures --only fig1,fig2DG,fig3BC,fig3D,fig4B,fig4C,s1,s2,s3,s6 --output-root reproduced/reviewer-v1-final-figures
python reproduce.py tables --output-root reproduced/reviewer-v1-tables
```

`validate`, `figures`, and `tables` run from included files only. `validate`
independently checks frozen panel populations and key endpoints. The plotting
command regenerates supported numeric panels from frozen inputs; `tables`
exports Figure 2 correlation statistics and frozen pair-weighted Figure 3
summaries. See
[`tables/manuscript/`](tables/manuscript/) for a checked-in machine-readable
snapshot and generation notes. Plot jobs record code, input, and environment
identities under ignored `reproduced/` output directories.

New outputs go under ignored `reproduced/` directories. Existing release
snapshots are retained and never overwritten. Existing outputs are reused
only when input/code/environment identities and output hashes agree; incompatible
outputs are not overwritten. Each plotting job records logs and a status manifest.
All persistent writes, temporary work and plot caches are directed into this
repository by the entry point. The source repository is never needed by it.

The numeric panels above are reproducible; this does not establish that each
saved composite is the final manuscript assembly. Current source artwork and
editable components are retained under [`figures/released/`](figures/released/)
with identity and status recorded in `manifests/ARTWORK.json`. Corrected Figure 3D
and current three-method Figure 4C assets are generated from bundled frozen tables.
Earlier case-pooled Figure 3D artwork and four-method Figure 4C composite remain
provenance-only. Supplementary
figures S1–S3 and S6 have validated standalone reproduction paths. S4, S5, and
S7 have current source artwork retained, but lack complete supported release
reproduction paths or have unresolved manuscript mapping; see the
[panel matrix](docs/REPRODUCIBILITY_MATRIX.md).

## What each reproduction level means

| Level | Status | What is available |
|---|---|---|
| Regenerate supported figures and tables | Fully reproducible for the panels listed in the matrix | The CLI uses included, checksum-pinned plotting inputs. |
| Independently validate included numerical results | Reproducible validation workflow; partial manuscript-result coverage | It checks populations, counts, selected endpoints and frozen table relationships. This is not a fresh simulation. |
| Repeat full simulations from frozen ensembles | Not supported | The repository does not include the complete production ensembles and executable run inputs. |
| Reconstruct metabolic models from original biological inputs | Provenance-only | Raw RNA/DMI inputs, the parent GEM and required historical artifacts are not bundled. |

See the [reproducibility matrix](docs/REPRODUCIBILITY_MATRIX.md) for figure and
table coverage. Development supports the fitted usefulness relationship, while
held-out Confirmation is reported separately. A2 uses the same split for its
post-freeze robustness analysis and is not an independent confirmation. Weak
directional cues are evaluated after A1/A2 strong-anchor weighting and do not
replace those quantitative anchors.

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

Use [CITATION.cff](CITATION.cff). Final manuscript title, DOI, ORCIDs and
conference citation require author-supplied publication metadata.

Source code is licensed under [MIT](LICENSE). Cleared original research products,
figures and documentation are licensed under [CC BY 4.0](LICENSE-DATA).
Third-party datasets, metabolic models, software, reconstruction methods and
externally sourced materials retain their original licenses and are not relicensed
by this repository. Ambiguous artifacts are withheld pending author review.
