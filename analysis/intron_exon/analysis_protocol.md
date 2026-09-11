# Full biochemical sweep and constrained pulse fitting

Protocol recorded before inspection of the new sweep and fitting results, 10 September 2026.

## Purpose

Evaluate where the full mapped-target model produces intronic profile leads, and test descriptive compatibility with the experimental intron/exon profiles using a restricted shared initiation input. The experimental offsets remain independent empirical evidence. The analysis does not infer protein exposure or claim uniquely measured biochemical rates.

## Fixed biological model

Reuse the analytical age-response kernels from the completed position-aware analysis. Transcription progresses at constant speed across the modeled span. Signal accumulates as mapped target bases are transcribed, with uniform contribution per base. The targeted intron is removed after an exponential processing wait from its 3-prime splice site. Exonic signal persists through transcription and processing, then decays exponentially after both target-intron removal and modeled transcript completion. Other introns, excised-intron persistence, export, individual oligo efficiency, and nonlinear HCR response are not separately represented.

All analyses use the previously mapped Tcas5.2 target intervals, not independently measured individual manufactured probe-pair coordinates. Transcript boundaries are the previous target-based proxies. Mapping and assembly qualifications remain applicable.

## Hypothetical-pulse sweep

- One compact asymmetric raised-cosine initiation pulse, total support duration D=1, with peaks at 0.25, 0.5, or 0.75 of its support.
- Traversal time, mean processing wait, and mean decay wait each take values 0.01, 0.03, 0.1, 0.3, 1, or 3 times D.
- The full factorial design contains 216 kinetic combinations per shape per gene, or 1,944 mapped scenarios overall.
- These are deliberately broad hypothetical scenarios. They are not an experimentally calibrated distribution of biological rates.
- Primary outcomes are exon-minus-intron rising and falling half-height offsets, each relative to its channel's own peak. Positive means the intronic crossing occurs earlier. Levels 0.4 and 0.6 assess threshold dependence.
- Classify offsets exceeding +0.01D as positive, below -0.01D as reversed, and the remainder as near zero. This is a declared display tolerance, not the experimental noise floor. Preserve continuous offsets and report sign-only fractions as well.
- Report outcomes per gene and shape, including reversals and any undefined crossings. Fractions describe equal weighting of this particular design, not biological probabilities.
- Strict molecular onset is a separate quantity. The mapped upstream exonic intervals preclude interpreting a profile half-height lead as proof that intronic targets first become detectable earlier.
- A zero-traversal comparison, if calculated, is reported separately and does not enter mapped scenario denominators.
- Check numerical convergence at extreme and near-classification-boundary scenarios.

## Constrained experimental fits

- Use the existing background-corrected paired specimen measurements at 24 eve-en substages. Pairing is inferred from numbered workbook slots. Distinct stages contain distinct embryos.
- Fit both observed channel-stage means jointly after scaling each channel by its observed mean maximum, with no subtraction of the observed minimum. Retain negative individual measurements.
- Primary initiation input is the sum of one asymmetric raised-cosine pulse for Kr and two for mlpt and svb. Each pulse has amplitude, peak time, rise duration, and fall duration.
- Use one kinetic parameter set and one relative channel gain per gene, shared across its pulses. Fit nonnegative residual channel offsets, as in the previous primary analysis.
- Maximum adjustable parameter counts are 10 for Kr and 14 for mlpt/svb in the mapped model, versus 9 and 13 in the zero-traversal model. Nonnegativity and boundary constraints can reduce effective flexibility.
- Kinetic duration search bounds match the previous analysis: traversal and processing 0.005 to 6 substages, decay 0.01 to 6; relative gain 0.2 to 5. Pulse peaks are ordered into [0,12] and, where applicable, [12,26]; rise/fall durations are 0.25 to 12 substages.
- Fit the full mapped model and the matched zero-traversal model. One declared alternative input family, asymmetric Gaussian pulses, checks dependence on the imposed shape. The raised-cosine family remains primary regardless of which fits better.
- Use multiple optimization starts and inspect numerical convergence. Do not add pulses or individual-stage adjustments to remove residual structure.
- Report fit errors, residuals relative to sampling uncertainty, boundary hits, and the number of parameters. A descriptive fit is not predictive validation or proof of the mechanism.

## Time and uncertainty

Hypothetical sweep times are ratios to pulse duration. Experimental fits retain the previous equally spaced eve-en substage coordinate. Unequal actual substage durations remain a limitation; no biochemical durations in minutes or hours are claimed. This execution does not undertake the separately discussed historical timing calibration.

Retain the existing 5,000-draw within-stage paired bootstrap estimates for experimental offsets. Those intervals condition on staging, the inferred pairing, and the selected pulse windows. They do not include stage-assignment errors, batch effects, or multiplicity correction. New fit assessment, if performed, must likewise preserve paired channels within each stage and report its assumptions explicitly.

