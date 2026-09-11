# Prespecified constrained fitting design

Written before new fits were evaluated, 2026-09-10.

The primary shared initiation input is one compact asymmetric raised-cosine
pulse for Kr and two for mlpt and svb. Each pulse has a center c, rise width R,
fall width F and nonnegative amplitude. It is sin-squared between c-R and c,
cos-squared between c and c+F, and zero elsewhere. First centers are in [0,12],
second centers in [12,26], and each width in [0.25,12] eve-en substage units.
Pulse counts and center windows follow the previously defined first and second
experimental pulse regions; they are restrictions, not independently known
transcriptional timings. A split Gaussian with separate rise/fall standard
deviations uses the same bounds as a single declared shape sensitivity. We will
report it regardless of whether its fit is better or worse.

The previously validated full mapped-region age kernels are retained unchanged.
One traversal duration T, target-intron processing mean, and postprocessing
decay mean applies to all pulses within a gene. Uniform signal per mapped target
base is assumed. T and processing means are bounded by [0.005,6], decay by
[0.01,6], and relative channel gain by [0.2,5]. An instantaneous-traversal
comparison fixes T=0 with all other input, kernel, and fitting assumptions
matched. These are effective stage-coordinate durations, not measured rates in
minutes. Stage indices remain equally spaced for this analysis; no absolute
timing calibration is inferred or imported.

Both channels are fit jointly to the 24 stage means, formed from complete paired
specimen slots. Negative background-subtracted observations are retained. Each
channel is divided by its observed full-series mean maximum, with no minimum
subtraction. Unweighted least squares is used. Input amplitudes and two
nonnegative channel residual offsets are profiled by NNLS. All other parameters
are bounded nonlinear parameters. The primary mapped model has 10 adjustable
quantities for Kr and 14 for mlpt/svb, against 48 channel-stage means (one fewer
for the collapsed comparison). Kernel unit-area normalization changes amplitude
and gain conventions only; it does not change the allowed predicted profiles.

Optimization uses differential evolution plus multiple local starts, including
the matched collapsed solution as starting points for mapped models. Numerical
convolution integrates age kernels analytically against a piecewise-linear
approximation to the pulse. Mesh refinement is checked on fitted predictions
and objective values. Initial history covers the complete compact pulse support,
and at least nine rising Gaussian standard deviations before the earliest
allowed center. No extra pulses or stage-specific initiation amplitudes will be
added in response to residuals. Fits establish compatibility, not unique
biochemical parameter identification. Boundary hits and structured residuals
will be reported explicitly.
