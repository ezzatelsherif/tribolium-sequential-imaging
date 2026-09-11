# Constrained fits: implementation and numerical checks

The analysis replaces stage-specific initiation amplitudes with one asymmetric
raised-cosine pulse for Kr and two for mlpt and svb. The full mapped-region age kernels are used. Ten fitted quantities are
available for Kr and fourteen for mlpt/svb, including input amplitudes and
channel offsets. A split Gaussian is a declared pulse-family sensitivity,
reported alongside the primary family rather than selected after seeing fits.
See fitdesign.md for the design declared before fitting.

For convolution, the analytic cumulative zeroth and first age-kernel moments
are used to compute exact kernel weights for linear interpolation of the input
pulse on a uniform grid. These weights are applied by FFT. This integrates
kernel transitions and short exponential waiting times, even when they are
shorter than the grid step. The pulse, not the age kernel, is discretized.
Negative FFT roundoff does not clip observed data. Tiny negative nodal weights
from cancellation in the remote kernel tail are clipped to zero.

Central optimization starts with differential evolution and multiple bounded
least-squares optimizations at grid step 0.05, followed by local refinement at
0.025 substage units. Input pulse amplitudes and two nonnegative residual
backgrounds are solved by nonnegative least squares at each nonlinear parameter
evaluation. The objective is unweighted joint squared residual error on 48
channel-stage means. The two channels are scaled by fixed observed mean maxima.
No noise clipping, minimum subtraction, or stage-dependent kinetic rates are
used. Mean input amplitudes do not have a biochemical interpretation because
the observation kernels are normalized to unit area for conditioning and the
relative channel gain remains free.

The final 0.025 predictions were compared with predictions at 0.0125 and
0.00625 without altering their biological or input parameters. Across all
twelve gene/pulse-family/model combinations the largest normalized prediction
change from 0.025 to 0.00625 was below 0.00049, or 0.049% of the normalized
channel maximum. The central 0.025 fits were retained. The additional script
validate_pulse_numerics.py independently compares FFT-based predictions against
adaptive quadrature of the original age kernel times the analytic input pulse.
Across 72 gene/family/model/time checks, its maximum unit-amplitude response
error was 0.000509 at step 0.025, 0.000127 at 0.0125, and 0.000032 at 0.00625,
consistent with second-order convergence of the input interpolation.
Numerical agreement verifies the implementation, not the biological assumptions.

One fall width reaches the lower bound of 0.25 in every fit (first pulse for Kr
and svb, second pulse for mlpt). Thus sharp cessation of initiation is not
resolved within that prescribed bound. All primary mapped fits put traversal T
at its 0.005 lower bound; this does not establish instantaneous biological
traversal. Matched T=0 comparisons produce essentially indistinguishable
predicted profiles. The Gaussian sensitivity permits different kinetic
trade-offs, especially for mlpt, without demonstrating a uniquely identified
traversal duration.

The restricted initiation family creates visibly larger residuals than the
stage-spaced free-input fits. Those residuals and the paired-resampling
lack-of-fit assessment must be retained when interpreting compatibility. This
reduction in flexibility is an explicit safeguard against following noise,
not an independent measurement of initiation and not a guarantee against
overfitting.

API: PulseFit(gene, mode='mapped', shape='cosine', dt=0.025).fit(y,
starts=[parameter_dictionary], global_search=False) performs a local refit.
Use global_search=True for differential evolution plus local searches.
PulseFit.predict(fit, at) predicts an existing fit, and load() returns complete
paired specimen arrays by gene and stage. The fixed central scale should be
retained when refitting resampled stage means. The command-line fit pipeline
explicitly uses a coarse global search followed by the 0.025 refinement.
