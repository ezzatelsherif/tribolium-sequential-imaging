# Full-model hypothetical-pulse sweep

The sweep supports a widespread later exonic falling half-height, while rising-phase ordering depends strongly on gene geometry. It does not support a universal intronic rising lead.

## Prespecified design

There are 1,944 mapped-model scenarios: three genes, three compact raised-cosine initiation pulses with peak at 25%, 50%, or 75% of pulse duration D, and 216 Cartesian-product combinations of total traversal, mean intron-processing wait, and mean RNA decay wait. Each duration takes the six values 0.01, 0.03, 0.1, 0.3, 1, and 3 times D. The protocol was written before results were computed. All scenarios are retained and equally weighted as a description of this discrete design. These are not biological probabilities or independently calibrated parameter ranges.

Mapped target regions and the existing biochemical kernels are unchanged. Signal contributions are uniform per target base; individual manufactured probe-pair efficiencies remain unknown. The model follows progressive target appearance, removal of the targeted intron after its splice site and an exponential processing wait, and exonic RNA decay after both modeled transcript completion and target-intron processing. These are expected signals, not individual stochastic trajectories.

## Primary result: half-height offsets

Positive offset means the exonic profile reaches the feature later than the intronic profile. Each channel is normalized to its own maximum above zero. The descriptive tolerance is 0.01 D and is not an experimental noise threshold.

| Gene | Phase | Intronic lead >0.01 D | Small separation <=0.01 D in magnitude | Exonic lead >0.01 D | Raw positive offsets |
|---|---|---:|---:|---:|---:|
| Kr | Rising | 616/648 (95.1%) | 32/648 | 0/648 | 648/648 |
| Kr | Falling | 621/648 (95.8%) | 27/648 | 0/648 | 648/648 |
| mlpt | Rising | 251/648 (38.7%) | 143/648 | 254/648 (39.2%) | 322/648 |
| mlpt | Falling | 542/648 (83.6%) | 104/648 | 2/648 (0.3%) | 618/648 |
| svb | Rising | 613/648 (94.6%) | 35/648 | 0/648 | 648/648 |
| svb | Falling | 621/648 (95.8%) | 27/648 | 0/648 | 648/648 |

The mlpt result is an important counterexample to a general rising-lead claim. About 54% of its mapped exonic target bases are upstream of the targeted intron, compared with approximately 8-9% for Kr and svb. In this defined model, finite traversal can therefore shift substantial exonic signal forward. At total traversal durations of one or three pulse durations, none of the sampled mlpt combinations has an intronic rising half-height lead greater than 0.01 D. This comparison concerns modeled target geometry, not experimentally determined oligonucleotide weights.

The two mlpt falling reversals greater than 0.01 D occur at traversal=0.3 D, processing mean=0.1 D, decay mean=0.01 D, with pulse peaks at 0.25 D or 0.5 D. Their offsets are -0.01818 D and -0.01599 D. Thus even the falling-phase ordering is not a mathematical invariant of every combination.

All molecular first-appearance claims must remain separate: each mapped gene includes exonic targets beginning at coordinate zero, upstream of the intronic targets. The modeled exonic signal can therefore first appear earlier even when its half-height is reached later. The sweep counts half-height ordering, not molecular onset.

## Threshold and comparison checks

Changing half-height to 40% or 60% leaves Kr and svb positive raw offsets in all sampled scenarios and leaves their falling counts unchanged at 621/648 beyond the tolerance. The mlpt falling intronic-lead count ranges from 538 to 543/648, while the rising count ranges from 247 to 313/648. Rising conclusions are therefore both geometry-dependent and more threshold-sensitive.

The separately tabulated instantaneous-traversal comparison has 108 unique processing/decay/shape combinations and is independent of gene geometry. All have positive raw rising and falling offsets. With the 0.01 D tolerance, 70/108 have a rising intronic lead and 76/108 a falling lead. These comparison scenarios are not included in the 1,944 primary counts.

## Numerical verification

Convolution uses analytic age-kernel bin integrals and midpoint pulse evaluation at dt=0.001 D. Checks used dt=0.0005 and 0.00025 D in 132 distinct mapped cases: every kinetic-grid corner for every gene and pulse shape, plus cases nearest the -0.01, 0, and +0.01 D boundaries. Maximum changes from the primary grid were 7.34e-7 D for rising offsets and 2.23e-6 D for falling offsets. No checked primary classifications or raw signs changed. Maximum change in normalized signal was 3.05e-6; maximum peak-offset change was 3.29e-5 D.

Independent quadrature of the point kernels and analytic pulses, checked at rise, peak, and decline for three diverse parameter/geometry combinations, agreed within 1.30e-6 of the channel maximum. There were no failed or missing crossings and no multiple threshold crossings in the 1,944 primary scenarios. Maximum mass truncation/discretization error was 4.30e-5 relative to expected full signal area. The tail extends well beyond every threshold crossing.

## Suggested interpretation

Across the specified parameter design, the full mapped model commonly produces later exonic decline, including in regimes where rising-phase order varies or reverses. This is consistent with using falling-phase separation as the clearer qualitative expectation. The result supports a plausible model behavior and should accompany direct experimental uncertainty estimates. It does not show that most biological parameter values have this behavior, uniquely identify mature-RNA decay as its cause, or establish a universal intronic onset lead.

## Files

`sweep.py` reproduces the sweep and all summaries and checks. `sweep_results/sweep_results.csv` contains every mapped scenario and measurement; `summary_by_gene.csv`, `summary_by_gene_shape.csv`, and parameter-stratified summaries retain exact counts. `collapsed_results.csv` is the separate matched comparison. `numerical_convergence.csv` and `independent_quadrature.csv` provide numerical validation. `sweep_phase_map.png/svg` shows intronic-lead counts out of 18 per traversal/decay cell; `sweep_signed_offsets.png/svg` retains the signed rise/fall offsets for all scenarios.
