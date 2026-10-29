# S1 Data workbook

`S1_Data_Landmark_Results.xlsx` consolidates the three result tables into one
supplementary workbook. It contains the pooled summary, separate peak and
half-height summaries, complete landmark inventory, original peak audit, source
summary, reader notes and column definitions.

The numerical reconstruction workflow regenerates the original CSV tables. The
workbook is the fixed publication layout for those results. Reproducing the
analysis does not require its spreadsheet-authoring software or regenerating
that layout. The original CSV filenames remain stable in the code for
provenance, although they are no longer separate supplementary files.

## Verify all workbook results

`verify_results_workbook.py` uses only the Python standard library. Supply the
publication workbook and the three CSV result tables. The CSV arguments can
point to freshly reconstructed outputs instead of the supplied reference files.

```sh
python verify_results_workbook.py \
  --workbook reference/S1_Data_Landmark_Results.xlsx \
  --peak-audit data/source_tables/peak_inventory.csv \
  --landmark-inventory outputs/S2_Data_Landmark_Comparisons.csv \
  --summary-source outputs/S3_Data_Landmark_Summary.csv
```

Use `--report verification.json` to save the report. The script checks:

- Every source-table cell, including blanks, text, Booleans and numerical values.
- All 859 landmark-inventory rows, 1,167 original peak-audit rows and 45
  summary-source rows.
- Each summary count and denominator against the eligible landmark inventory.
- Every displayed count and percentage against its source summary.
- The saved summary formula values and absence of Excel error cells.

The verifier does not execute Excel formulas. The supplied workbook includes
working summary formulas and saved results; its formula recalculation was
separately checked by changing and restoring a source count. The verifier reads
the saved numerical results and independently recomputes the relevant counts
and percentages from the CSVs.

## Counting units

The landmark inventory has one row for an event pair, landmark and spatial
stage, with all three spatial embryos recorded on that row. Its 233 eligible
comparisons define the Fig. 7B denominator. The original peak audit has one row
per spatial embryo and retains its earlier eligibility fields. Do not combine
the two inventories' denominators. The 15 pooled summary rows and 30 stratified
rows are two views of the same eligible comparisons, not additional data.

Source-table numerical values are retained as binary64 numbers. Display formats
shorten decimal places without changing the stored values. The source summary
uses percentages on a 0 to 100 scale, while formatted summary formulas use
fractions on a 0 to 1 scale.
