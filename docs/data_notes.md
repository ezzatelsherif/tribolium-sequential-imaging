# Data provenance and interpretation notes

These notes record the relationship between the supplied measurements,
original plotting code, and `Manuscript_ver44.docx`. The repository preserves
the source measurements and reproduces the recorded plotting conventions.
It does not silently correct measurements or change the manuscript's results.

## Spatial normalization

The original spatial viewer subtracts each gene profile's minimum and divides
by its range, separately for each embryo:

`normalized = (intensity - minimum) / (maximum - minimum)`

This is implemented before Savitzky-Golay smoothing. The manuscript's spatial
methods currently describe division by the maximum alone. Those operations
have the same affine curve shape but different baselines and amplitudes.
The repository follows the source code. The manuscript wording should be
reconciled with this actual computation before submission.

This does not change the separate Fig. 3/S4 intron/exon convention: those
paired temporal profiles are divided by their maximum stage means without
subtracting a stage-mean minimum.

## eve reference in temporal-to-spatial comparisons

The original viewer concatenates eve observations from all seven temporal
workbooks at each substage before calculating the reference curve in Fig. 6,
S1, and S2. `run.xlsx` and `odd.xlsx` contain identical eve measurements, so this
procedure includes that shared cohort twice. It is a weighted display
reference, not a set of independent additional eve embryos. Counts in the
pooled output table therefore count workbook entries, not unique embryos.

The six gene-specific temporal curve shapes and this pooled eve curve match
the supplied comparison figures. The default configuration retains that
pooling for reproducibility. The Fig. 5 cohort summaries are exported
separately, with each gene's own co-measured eve. The Python function
`comparison_temporal(measurements, eve_reference="runt")` can instead calculate
the shared pair-rule cohort's eve reference for a sensitivity comparison.
That choice changes the comparison reference and should be reported if used.

The manuscript correctly identifies shared run/odd/eve measurements, but its
Fig. 6 interpretation should explicitly distinguish that cohort from the
pooled eve reference actually plotted. Do not interpret the pooled reference
as additional independent biological replication or as a single co-measured
pair-rule cohort.

## One missing anterior-background measurement

In `KR.xlsx`, `Sheet1`, stage **8.1**, **eve specimen slot 5** has posterior
intensity **1380** and a blank anterior-background measurement. Its cached
`BASE` value is nevertheless **1380**, consistent with a worksheet treating
the blank as zero. The original viewer reads that cached value.

The raw correction exported by this repository is missing for that slot, and
`complete` is false. The separate `workbook_display` column retains 1380 to
reproduce the original viewer. No background has been measured or imputed by
this release. All complete raw A-B values agree with the cached values.
The verifier explicitly checks and reports this one exception.

This affects the late-stage eve summary for the Kr cohort and the pooled eve
reference; it does not remove a Kr intron/exon pair. The background should be
resolved from the measurement record, or the choice to retain/exclude this
cached value should be reported. Do not overwrite the source workbook to hide
the discrepancy. The raw-correction summaries provided alongside the cached
display summaries allow its effect to be examined.

## Spatial replicate source versions

The stage 4.2 EMB_3 file in the Drive collection differs from the version in
Theophilus Frimpong's intron/exon repository. The GitHub version matches all
five curves of the additional stage 4.2 embryo in Fig. S5 and is included here.
Its upstream blob is `2c483e35fe2e7b54065afcb5222046083a57fb5c`, at commit
`7121c65d65b4d8979282d9ebb1edb97a5af5b304`. `data/manifest.csv` links to this
immutable version. The additional stage 4.3 embryo matches the supplied Drive
EMB_3 file. Original filenames are retained.

The input selection was established by comparison with the figure curves,
not by assuming that file order or embryo numbering implied figure order.

## Fig. 6 panel selection

In the supplied `Fig 6 - Comparing time to space - ver09.pdf`, panels B and C
contain the two selected gap genes without an eve curve. Their manuscript
legend still names eve. The supplemental counterparts include eve. The
repository's configuration follows the supplied figure curves; the final
legend and panel selection should be made consistent before submission.

## Scope of reproduction

The numerical workflow starts after fluorescence extraction. It does not
reproduce image acquisition, DAPI-based registration, manual ROI placement,
or complete microscopy composites. The original registration implementation
was not present in the supplied revision code package. The registration
description in the manuscript alone is insufficient to recover its original
implementation and settings. A complete reproduction from raw images would
also require those materials.

The manuscript's current data-availability paragraph cites the earlier
intron/exon-only commit `4d3951c`. A submission using this expanded release
should cite its new immutable commit and describe the additional spatial and
temporal datasets.
