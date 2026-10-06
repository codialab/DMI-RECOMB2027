# Supplementary Figure S3

## Purpose

S3 shows the absolute A1-versus-A2 relationship for two descriptors of the residual sign distribution. Each plotted pair represents the same reaction and evaluation under both anchoring settings. This view complements S2's paired differences and preserves reaction/evaluation heterogeneity.

## Frozen inputs and pairing

The table builder reads the four method partitions in figures/fig4/data/fig4_geometry_paired/ and verifies each partition checksum against figures/fig4/data/fig4_data_manifest.json. That manifest identifies the source as FIG4-DATA-B, version 1.1.0, status PASS. Its geometry combines the frozen A1 PL1 and A2 A21 descriptors. The S3 builder does not rerun geometry calculations or reconstruct the pairing.

The authoritative pair key is (evaluation_id, reaction_id). Method, RNA context, CT2A mouse, and GL261 mouse are retained as pair metadata. Figure 4's builder created the table with a one-to-one A1/A2 outer merge, required every key to match, and checked equality of those metadata fields. S3 verifies the resulting unique pair keys and metadata.

## Panels

- **Directional entropy \(H_{\mathrm{dir}}\):** A1 on x and A2 on y. The production definition is Shannon entropy across negative, tie, and positive sign masses, normalized by \(\log(3)\); valid values are in \([0,1]\).
- **Dominant-direction mass:** A1 on x and A2 on y. This is the largest of the negative, tie, and positive sign probabilities; valid values are in \([0,1]\).

The four methods are pooled for both panels as specified for S3. The hexbin colors show matched reaction–evaluation pair counts on a log scale. The dashed \(y=x\) line marks equal A1 and A2 values. Each panel reports its exact \(n\) and the fractions above, within tolerance of, and below the identity line. Both metrics are dimensionless.

## Build and outputs

From the repository root, run:

    python figures/supp_fig3/build_supp_fig3_tables.py
    jupyter nbconvert --to notebook --execute --inplace figures/supp_fig3/supp_fig3.ipynb

The builder writes the compact paired table and its provenance/check summary to data/. The notebook writes:

- outputs/supp_fig3.svg
- outputs/supp_fig3.png

The frozen population has 1,671,600 pairs: 400 evaluations × 4,179 reactions, with 417,900 rows per method. Both panels use all 1,671,600 finite pairs; no observations are excluded.

## Limitations

The density display summarizes rather than individually labels the observations. The pooled view does not show method-specific shifts; S2 supplies the method-stratified paired-difference view. Reaction–evaluation pairs are the defined descriptive units and should not be interpreted as independent statistical observations.
