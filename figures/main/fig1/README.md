# Figure 1 artwork

The current manuscript uses four panels (1A–D). The editable artwork is
`figures/released/fig1/recomb_figure1_editable_2.svg`; the publication PNG is
`figures/released/fig1/recomb_figure_1.png`. The SVG was adopted from source
revision `3e2fffc4d7d18469cc9e5affe9c274178413c212`; its panel B copy is corrected
to “A directional cue does not directly convey magnitude information.”

Panel C agrees with the canonical Figure 1 caption: both four-state examples have
equal initial weights and sign probabilities; the contrast sets are
(-7, -3, +3, +7) and (-9, -7, +1, +3), with η² values 0 and 0.90 and updated
magnitude estimates 5 and approximately 4.27. Panel D depicts 4 methods × 4 RNA
contexts × 25 cross-tumor mouse comparisons = 400 evaluations across 4,179
non-anchor reactions.

This conceptual figure has no numerical data dependency. The publication
workflow packages the editable SVG and renders its PNG with Inkscape at 2,400
pixels wide, matching the updated source render:

```bash
python reproduce.py figures --only fig1
```

For a direct Inkscape export, use the same source and `--export-width=2400`,
`--export-background='#ffffff'`, and `--export-background-opacity=255` options.

Source identity, export details and current hashes are recorded in
`manifests/ARTWORK.json` and `manifests/checksums/SHA256SUMS.txt`.
