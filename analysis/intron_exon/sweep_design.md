# Hypothetical-pulse sweep design

Written before computing sweep outcomes, 2026-09-10.

## Scope and fixed assumptions

Use the previously implemented analytic mapped biochemical observation kernels, unchanged. Their input is a common transcription-initiation pulse, not a measured intronic curve. Target signal is uniform per mapped target base, not a reconstruction of manufactured oligonucleotide efficiencies. Mapped target-span boundaries proxy transcription-unit boundaries. The kernels include progressive target transcription, exponential target-intron processing after its 3-prime splice site, and exponential exonic RNA decay starting after both transcript completion and target-intron processing. They omit other introns, excised-intron persistence, termination and export delays, and nonlinear detection. Gene geometries are the existing Tcas5.2 Kr, mlpt, and svb mappings.

## Primary design

The hypothetical initiation pulse has support duration D=1, onset zero, amplitude one, and peak at r in {0.25,0.50,0.75}. For 0 <= t <= r, u(t)=sin^2(pi*t/(2*r)); for r <= t <= 1, u(t)=cos^2(pi*(t-r)/(2*(1-r))); otherwise u(t)=0. These three shapes vary rise/fall asymmetry while holding duration and integrated production constant.

Total traversal duration T/D, mean target-intron processing wait/D, and mean completed/processed exonic RNA decay wait/D each take values {0.01,0.03,0.1,0.3,1,3}. All Cartesian-product combinations are used. This gives 216 kinetic combinations x 3 shapes x 3 genes = 1,944 primary scenarios. The ranges are exploratory, not calibrated biological priors. Equal weights define this discrete design and are not probabilities of biological mechanisms. No ranges or weights will be selected in response to the direction of the resulting offsets. Rate values are durations, not rates or half-lives.

## Measurements and classifications

For each channel separately, normalize to its own maximum above the model's zero baseline. Record the first upward and final downward crossings of 50% maximum, with 40% and 60% as prespecified threshold sensitivities, and the time of the global peak. Record additional crossings if present rather than silently selecting a local pulse. Signed offsets are exon minus intron in units of D. Positive means intronic profile reaches the chosen feature earlier. The primary practical classification uses tolerance 0.01 D: intronic lead >+0.01, small separation in [-0.01,+0.01], and exonic lead <-0.01. This is a declared descriptive tolerance, not an experimental noise threshold or significance test. Also retain raw signed offsets so near-zero positive values are visible. Distinguish these profile statistics from molecular first appearance. Failed or invalid crossings must be counted and retained.

## Matched comparison

Compute instantaneous-traversal kernels at T=0 for each processing wait, decay wait, and pulse shape. This gives 108 unique comparison scenarios because collapsed kernels do not depend on gene geometry. Report these separately from the 1,944 mapped scenarios and never inflate their weighting by counting repeated T/gene combinations.

## Numerical computation and checks

Evaluate convolution using exact analytic kernel integrals in age bins and midpoint initiation samples, accelerated with FFT. Initial step is dt=0.001 D. Include an age horizon T + 12*max(processing wait,decay wait), extending the signal through completion of the input pulse; this is far below the smallest crossing threshold at the tail. Validate at dt=0.0005 and 0.00025 on all 8 kinetic-grid corners for each of 3 genes and 3 shapes, plus up to 12 scenarios closest to each classification boundary (+/-0.01 D) for each primary phase. Record changes in primary crossing offsets, classifications, peak offset and peak-normalized signal. If crossing classifications are unstable, refine those scenarios further and explicitly disclose refinements. Peak times of exceptionally flat profiles may be less well conditioned than half-height crossings. Record tail heights and mass conservation errors. Check at least a few independent convolution values using numerical quadrature of the unchanged point kernel and the analytic pulse. Any failure is reported, not discarded.

## Outputs and interpretation

Save all scenario rows, stratified counts by gene and pulse shape, counts across parameter combinations, matched comparison summaries, and numerical checks. Main figures show the fraction of the specified scenarios with positive practical separation within defined parameter slices and make the denominator explicit. The sweep establishes conditional behavior of this defined model family; it neither estimates biological parameter frequencies nor establishes the cause of observed offsets.

## Numerical-check clarification after initial run

Because raw sign-only counts are also reported, the nearest-to-zero scenarios were added to convergence checks in addition to the already specified +/-0.01 D classification boundaries. This changed no model, parameter range, shape, weight, or statistic. The final check set contains 132 distinct mapped scenarios, each evaluated at both finer time steps.
