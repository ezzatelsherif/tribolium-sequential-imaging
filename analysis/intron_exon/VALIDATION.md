# Publication-package verification

The package was checked during preparation on 11 September 2026. The scientific
kernel equations, input data, parameter grid, fit restrictions, and recorded
model-fit JSON outputs were retained from the completed analysis.

The following checks passed without rerunning the long central fits or the
199-draw fit-assessment optimizations:

- Reconstructed all 360 background-corrected numbered slots from the three
  included source workbooks and matched the supplied corrected-data table.
- Reconstructed all empirical point estimates using the included threshold,
  interpolation, and pulse-window definitions.
- Reconstructed all 12 stored fitted profiles, with maximum prediction
  difference zero in the checked environment, and reproduced their objectives.
- Recomputed gene-level counts from all 1,944 recorded mapped sweep scenarios
  and independently reran one scenario per gene.
- Recomputed every conditional lack-of-fit p value from the retained null draws.
- Verified that all six supplied target sequences match the archived genomic
  intervals exactly, except for the specified excluded mlpt terminal 17 A's.
- Ran a 20-draw empirical-bootstrap smoke check in a separate temporary output
  directory. The published 5,000-draw intervals were retained unchanged.
- Regenerated analysis figures and fit-summary/normalized-curve tables.
- Regenerated final-fit mesh diagnostics at 0.025, 0.0125, and 0.00625 stage
  units. The diagnostics match the original final-fit table within 1e-12.
- Compiled all Python modules successfully.
- Ran the portable manuscript-figure script without the original Illustrator
  file, reproducing the replacement numerical panels and S3/S4 figure files in
  a separate temporary output directory. The optional FFmpeg movie execution
  was verified during preparation of the accompanying manuscript figures.

The publication entrypoint adds empirical regeneration and numerical checks to
the full execution sequence. Figure-only regeneration now includes the sweep
figures. Final-fit table exports are regenerated alongside fit plots. An explicit
`--fresh` option allows optimization without archived parameter warm starts.
Machine-specific fallback paths and the earlier document-packaging script were
removed. These are portability and packaging changes, not changes to the
scientific model or the selected results.

The preserved adaptive-quadrature tables document the original independent
numerical check. They can be regenerated with
`python reproduce_analysis.py --stage numerics`.

These checks establish internal numerical and data consistency. They do not
verify specimen identity, stage timing, biological assumptions, global optimizer
optimality in every resample, or predictive generalization.
