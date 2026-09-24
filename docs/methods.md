# Numerical processing of expression profiles

The original viewer implementations are credited in [CITATION.md](../CITATION.md).
Their numerical transformations are implemented in `seqimaging/data.py` and
`seqimaging/figures.py`. Source selection is explicit in `config/figures.json`.

## Spatial line profiles

For each selected embryo and channel, the original coordinate and fluorescence
are independently mapped to their observed 0-1 ranges. Fluorescence is then
smoothed by a Savitzky-Golay filter of polynomial order 3. The window is the
original viewer's odd integer rule, approximately 10% of the AP bins:

```python
window = int(max(order + 2, (n * fraction) // 2 * 2 + 1))
window = min(window, n if n % 2 else n - 1)
```

The line-plot fraction is 0.1; the block-display fraction is 0.15. Smoothing
uses SciPy's default interpolated boundary treatment. Small overshoots after
smoothing are retained in the numerical table and line curves. Display blocks
clip their input to 0-1. Each spatial curve is one embryo, without averaging
across biological replicates.

## Temporal summaries

The measurements are posterior ROI A minus anterior background ROI B. Channels
are summarized separately within each cohort and assigned stage. The source
viewer uses cached `BASE` values and treats cached zeros as missing. For the
included source data, no complete raw A-B observation equals a discarded zero;
the original raw/cached completeness discrepancy is recorded in
[data_notes.md](data_notes.md). The current workflow applies the confirmed
background correction from `data/temporal_corrections.csv` and plots A-B
values. Original source fields and cached values remain in the exported table.

For each channel, let `mean[s]` be the stage mean and let `lo` and `hi` be the
minimum and maximum of those means over all 24 stages. The display values are
`(mean[s] - lo) / (hi - lo)`. SEM is the sample standard deviation (`ddof=1`)
divided by the square root of the available measurement count, then divided
by `hi - lo`. SEM describes variability conditional on assigned stage. It
does not include staging uncertainty.

PCHIP interpolation through all 24 stage means supplies a 400-point display
curve. The SEM envelope is interpolated separately. Temporal comparison panels
display this same full-series interpolation cropped at the indicated stage;
they do not renormalize or refit interpolation separately for each crop.
Stages occupy indices 1-24 with equal spacing, without conversion to minutes.
The comparison eve curve pools the six cohorts, counting the shared run/odd
eve measurements once, as detailed in the data notes.

Fig. 3/S4 uses its own paired specimen filter and maximum-only normalization.
See the [intron/exon protocol](../analysis/intron_exon/analysis_protocol.md)
for the empirical half-height estimates, paired bootstrap, mapped observation
kernels, kinetic sweep, joint fits, and restrictions on interpretation.

## Color blocks

Blocks are a nonlinear visualization of expression. They do not supply
additional measurements and are not inputs to the intron/exon model.

Normalized spatial curves are smoothed with fraction 0.15. Temporal blocks use
linear interpolation of stage means. Each is clipped to 0-1 and resampled at
800 equally spaced positions. The local neighborhood excludes the current
position and has half-width 75 samples for eve, 65 for run/odd, and 50 for the
other genes. At boundaries the neighborhood is shortened.

For each position, use the neighborhood minimum as a baseline, except that
when the intensity exceeds that minimum by more than 0.0025, use the mean of
the ten lowest neighboring values. The local contrast is the absolute
intensity-minus-baseline difference. The contrast weight is 1 for intensities
above 0.75; otherwise it is the contrast divided by 0.175, capped at 1.

Above intensity 0.1, the displayed value is `0.3 * intensity + 0.7 * weight`.
At or below 0.1 it remains the unenhanced intensity. The cad-specific rule in
the viewer uses coefficient 0.9 where its intensity exceeds 0.05. Values are
clipped to 0-1 and mapped from white to the gene's color. These settings follow
the viewer for the included developmental stages, all at or after stage 3.3.

## Files and reproducibility

Input tables are read without modification. Figure products and derived tables
are written under the selected output directory. `settings.json` records the
figure configuration. The original S1 Code package has a separate release
manifest; model reruns take place in a copy so its recorded results can still
be checked against that manifest.
