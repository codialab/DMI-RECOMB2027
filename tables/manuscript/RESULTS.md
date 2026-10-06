# Frozen manuscript result tables

Generated with `python reproduce.py tables` from included frozen Figure 2/3
source tables. The generation function is `tables()` in `reproduce.py`; table
outputs are TSV with source column names preserved.

- `fig2_panel_statistics.tsv`: all plotted reactions per panel and Spearman
  correlation between direction-explained magnitude variance and directional
  usefulness fraction. Population is the frozen point table, one row per
  plotted reaction within a panel.
- `fig3_pair_gain_summary.tsv`: frozen pair-weighted gain descriptors, retained
  by anchor and cue direction. Cases are averaged within matched truth pairs;
  pairs receive equal weight. Zero and missing counts are retained as stored.
- `fig3D_pair_weighted_summary.tsv`: Figure 3D panel values regenerated from the two frozen summaries; includes both cue means, their correct-minus-wrong separation, and pair counts for A1 and A2-L. The standalone plotter verifies the subtraction against the frozen endpoint table.
- `fig3_full_endpoint_separation.tsv`: frozen full correct-minus-wrong endpoint
  summary. This is the pair-level primary aggregation, not the superseded
  case-pooled audit table.

The generation is deterministic apart from library-level floating-point
serialization. `python reproduce.py validate` checks source hashes, Figure 2
panel populations/correlations, Figure 3 pair populations and arithmetic, and
Figure 4 paired keys. No upstream science is rerun.
