# Source measurements

Only source datasets used in the manuscript are included here. File contents
are preserved, including original headers, unused columns within a selected
source file, and missing values. Confirmed corrections are supplied separately in
`temporal_corrections.csv` and applied by the reader. The figure workflow selects the channels
listed in `config/figures.json`; it does not treat every column as a separate
published result. Original provenance links are recorded for traceability,
but reading the repository does not require access to those links.

## Spatial CSVs

`spatial/exonic/` contains 12 embryos: EMB_1, EMB_2, and EMB_3 at each of stages
3.3, 4.1, 4.3, and 5.2. These supply Fig. 6, Fig. S4, and Fig. S5, respectively.
Stage 4.1 EMB_1 also supplies the Fig. 2C profile example.

`spatial/intron_exon/` contains five embryos. EMB_1 at stages 3.3, 4.2, and 4.3
supplies Fig. 4B. EMB_1 and EMB_3 at stages 4.2 and 4.3 supply Fig. S3. Its
"replicate 2" is the original **EMB_3** file, not EMB_2. Stage 4.2 EMB_1 also
supplies Fig. 3F spatial overlays.

Embryo identifiers are local to the assay and stage. EMB_1 at two different
stages does not mean the same embryo followed through time. Columns within
one spatial file refer to channels measured in that embryo.

| Field | Meaning |
|---|---|
| `Distance_(microns)` | Original AP-bin coordinate, preserved as supplied; plotted after mapping its observed range to 0-1 |
| `hb`, `kr`, `mlpt`, `gt`, `svb`, `eve`, `runt`, `odd`, `en`, `cad` | Gene fluorescence channel; `runt` is the file label for run |
| `krint`, `mlptint`, `svbint` | Intronic fluorescence channels |
| `eve-en` | Original combined/reference-channel label in intron/exon files |

Intensities are fluorescence measurements in arbitrary units. The repository
does not infer a new physical calibration from the distance-column header.
Entirely empty rows are ignored on reading. Missing values in a requested
channel cause an error; no per-channel position deletion or imputation is used.

The selected stage 4.2 EMB_3 source has 49 missing `svbint` values at coordinates
500-548 and four entirely empty trailing rows. Fig. S3 uses kr, krint, mlpt,
mlptint, and cad, whose measurements are complete. That unused svbint column
is retained in the original CSV and is excluded from the figure workflow.

## Temporal workbooks

`temporal/` contains `HB.xlsx`, `KR.xlsx`, `MLPT.xlsx`, `GT.xlsx`, `SVB.xlsx`,
`run.xlsx`, and `odd.xlsx`, covering 24 ordered eve-en substages from 1.1 to 8.3.
Stages are labels, not decimal-valued elapsed time. Plots use equally spaced
stage indices 1-24.

| Worksheet or field | Meaning |
|---|---|
| `Sheet1` | Original posterior (A) and anterior-background (B) measurements |
| `Stage` | Assigned eve-en substage |
| `eve1A` / `eve1B`, etc. | Posterior and background for one numbered eve specimen slot |
| `mrna1A` / `mrna1B`, etc. | Corresponding gene's exonic/mRNA channel |
| `intr1A` / `intr1B`, etc. | Intronic channel where present |
| `BASE` | Cached background-corrected values read by the original plotting viewer |

Numbered slots 1-5 identify paired measurements **within one workbook row**.
They are not persistent embryo IDs across stages or across gap-gene cohorts.
The run and odd workbooks share the same eve measurements. Gap genes were
measured in separate cohorts. Original source blanks remain visible, and negative corrected intensities
remain negative. One recovered background value is supplied separately in
`temporal_corrections.csv`, with its measurement provenance.

`python reproduce.py profiles` exports `temporal_measurements.csv` with raw
posterior/background, `corrected` (A-B), `workbook_display`, and `complete`
columns. `source_anterior_background`, `source_corrected`, and `source_complete`
preserve the original workbook values, and `correction_id` links any corrected
row to the correction table. `workbook_display` always means the original
cached display value, not the corrected plotting value.
`temporal_stage_summary.csv` reports counts, means, sample standard deviations,
SEM, and normalization parameters for both definitions. Current figures use
`corrected`. The [data notes](../docs/data_notes.md) document the recovered
background of 328 and the shared run/odd eve cohort, which is counted once.

## Intron/exon model inputs

The original S1 Code inputs remain under `analysis/intron_exon/`: paired
workbooks, all 360 numbered corrected slots, target sequences, reference spans,
and mapped coordinates. Their [data dictionary](../analysis/intron_exon/data/README.md)
defines coordinate conventions and target-region qualifications. Hypothetical
simulation outputs are labeled separately from observed embryo measurements.
