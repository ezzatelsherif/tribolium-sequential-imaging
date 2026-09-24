# Data provenance and interpretation notes

The original source workbooks and spatial CSVs are preserved byte for byte.
Confirmed measurement corrections are recorded separately and applied by the
reader. The current settings accompany `Manuscript_ver46.docx`.

## Spatial normalization

Spatial profiles are normalized separately for each channel and embryo:

`normalized = (intensity - minimum) / (maximum - minimum)`

This precedes Savitzky-Golay smoothing, which can produce small excursions
beyond 0-1. Manuscript versions 45 and later describe this computation.
The separate Fig. 3/S4 intron/exon analysis divides paired temporal profiles
by their maximum stage means without subtracting a stage-mean minimum.

## Confirmed anterior-background correction

In `KR.xlsx`, `Sheet1`, stage **8.1**, **eve specimen slot 5**, the original
posterior value is **1380**, the anterior-background cell is blank, and the
cached `BASE` value is **1380**. Theophilus Frimpong rechecked the embryo and
reported a mean background of **328**, in a reply supplied by Ezzat El-Sherif
on 2026-09-24. The corrected fluorescence is therefore **1052**.

`data/temporal_corrections.csv` records that measurement and its provenance.
The default reader applies it to `anterior_background`, recomputes `corrected`,
and sets `complete` to true. The original values remain in the `source_*`
columns, while `workbook_display` retains the historical cached 1380.
`correction_id` identifies the affected row. No background was imputed.
`temporal_measurements(apply_corrections=False)` returns the historical data.

Current plots use the corrected A-B measurements. The correction affects eve
in the Kr cohort and the pooled eve reference, including normalization based
on the full 24-stage series. It does not change any Kr intron/exon measurement,
pair count, fit, offset estimate, or simulation.

## eve reference in temporal-to-spatial comparisons

Fig. 6, S1, and S2 use a common eve reference pooled from six temporal cohorts:
the five gap-gene cohorts and the shared run/odd cohort. The run and odd
workbooks contain the same eve measurements. The default `pooled_unique`
setting verifies that correspondence and includes the run copy once, omitting
the repeated odd copy. It does not deduplicate unrelated observations merely
because their numerical values happen to match. Available specimen values
are pooled within substage before averaging and full-series normalization.

The original viewer concatenated all seven workbooks, so the shared pair-rule
cohort entered twice. Frimpong confirmed this duplication in the same reply.
The `pooled` option retains that historical weighting for comparison. Use
`comparison_temporal(temporal_measurements(apply_corrections=False),
"pooled", "workbook_display")` to recover the original curves exactly.

Fig. 5 retains separate summaries for each gene and its co-measured eve.
The common comparison reference is not a single co-measured pair-rule cohort.
Fig. S1/S2 add spatial embryos, not new temporal cohorts.

## Effect of the corrections

`python -m seqimaging.correction_review` reproduces the comparison tables and
plot in `outputs/correction_review/`. Relative to the original displayed pool,
the combined background correction and single counting of the shared cohort
change normalized eve stage means by at most **0.035492** on the 0-1 scale
(stage 5.1; RMS change **0.015406**). All eight local peak substages remain
the same. Pooled counts are now 25-30 available specimens per substage.

Using the run cohort alone instead changes stage means by as much as
**0.202694** and changes three local peak substages. This alternative is
provided for comparison, not used as the default. Unchanged peak substages
do not establish that every threshold crossing or phase estimate is identical.
The other temporal gene profiles and all spatial measurements are unchanged.

The regenerated numerical panels supersede the temporal curves in the earlier
assembled figures. The original Illustrator/microscopy composites have not
been edited by this numerical workflow.

## Spatial replicate source versions

The stage 4.2 EMB_3 file in the Drive collection differs from the version in
Frimpong's intron/exon repository. The GitHub version matches all five curves
of the additional stage 4.2 embryo in Fig. S5 and is included here. Its upstream
blob is `2c483e35fe2e7b54065afcb5222046083a57fb5c`, at commit
`7121c65d65b4d8979282d9ebb1edb97a5af5b304`. `data/manifest.csv` records this
immutable provenance. The additional stage 4.3 embryo matches the supplied
Drive EMB_3 file. Original filenames are retained. Source selection was
established by comparing figure curves, not by assuming file or embryo order.

## Fig. 6 panel selection

Fig. 6B and C intentionally omit eve at the reviewer's request, as confirmed
by El-Sherif. Their supplemental counterparts include eve. The configuration
and manuscript legends retain those selections.

## Scope of reproduction

The workflow starts after fluorescence extraction. It does not reproduce
image acquisition, DAPI-based registration, manual ROI placement, or complete
microscopy composites. The original registration implementation was not in
the supplied revision code package. Reproduction from raw images would also
require that implementation, its settings, and the image inputs.
