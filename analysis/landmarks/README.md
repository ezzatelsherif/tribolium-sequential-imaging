# Reconstructing the temporal to spatial landmark analysis

This directory reproduces the complete peak and half-height comparison inventory,
the pooled and stratified counts, and the final six-panel Figure 7. Figure 7A
introduces the three landmark types, and Figure 7B combines their ordering
results in one summary. The numerical inputs, original analysis code and saved
biological-bootstrap arrays are all included in this repository. No Drive
access or external data download is needed.

## S1 Data workbook

The final submission uses one supplementary results workbook with these sheets:

| Sheet | Contents and counting unit |
|---|---|
| `ReadMe` | Navigation, provenance and the distinction between counting units |
| `Overall summary` | Pooled landmark counts and percentages underlying Fig. 7B |
| `By landmark` | Separate peak and combined half-height summaries |
| `Landmark inventory` | All 859 candidate comparisons, with 233 eligible comparisons and the reasons for 626 exclusions |
| `Original peak audit` | The original 1,167 embryo-by-event records, retained separately because their counting unit and eligibility differ |
| `Summary source` | The complete 45-row summary source: 15 pooled rows and 30 stratified rows |
| `Data dictionary` | Definitions, units and source-table mappings |

The original peak audit and current landmark inventory must not be combined
into one denominator. The pooled and stratified summary rows describe the same
observations from two views and must not be added together. CSV filenames
containing `S2_Data` and `S3_Data` are retained as stable numerical source names;
these are now sheets in S1 Data, not separate supplements.

## Reproduce the submitted results

Use Python 3.12 and run from this directory:

```bash
python -m pip install -r requirements.txt
python reproduce.py
```

This rebuilds the complete comparison table from the supplied temporal feature arrays and spatial landmark tables, then writes:

- `outputs/S2_Data_Landmark_Comparisons.csv`
- `outputs/S3_Data_Landmark_Summary.csv`
- `outputs/S1_Data_Landmark_Results.xlsx`
- `validation/replay.json` and `validation/workbook_replay.json`

The command independently verifies all 70,406 source cells and every displayed
summary in the publication workbook against the reconstructed CSV tables, then
copies that verified workbook into `outputs/`. The workbook retains its fixed
publication layout and working formulas. The verifier uses only the Python
standard library; reproducing the numerical results does not require the
original spreadsheet-authoring software. See [WORKBOOK_README.md](WORKBOOK_README.md)
for the full check and its scope.

The saved temporal arrays contain the original observation followed by 20,000 biological-bootstrap draws. The command calculates the direct temporal differences from those joint draws. It does not create new embryo resamples. `reproduce.py` checks both output tables against the supplied references, checks the original comparison table, and verifies the sizes of all saved arrays.

An isolated replay containing only the scripts, source tables, temporal feature cache and reference files reproduced all three tables byte-for-byte. The ver03 pooling changes the presentation and summary denominators; all 859 original measurements, eligibility decisions, exclusions and category assignments match ver02. All 109 temporal feature arrays and all 226 temporal contrast arrays matched the reviewed analysis exactly, including missing values. The validation record is in `validation/replay.json`.

## Reproduce the updated Figure 7

Install the figure dependencies and the Inkscape command-line application, with `inkscape` available on the command path. Then run:

```bash
python -m pip install -r requirements-figure.txt
python reproduce.py --figure
```

This command first reconstructs the numerical analysis, copies the regenerated comparison and summary tables into the figure inputs, and builds the final Figure 7 as editable single-layer SVG, vector PDF, 600 dpi TIFF and a PNG preview under `figure7/`. The four existing case panels retain their underlying curves, gene colors, peak markers and intervals, with their current C-F labels unchanged. Their plots are rearranged into a landscape layout. Panel A illustrates the peak, rising half-height and falling half-height on an expression profile. Panel B pools the three landmark types in one summary, followed by the existing examples as panels C-F. The illustration explains the definitions and is not an additional biological measurement.

To rebuild the artwork directly from the included frozen input tables, run:

```bash
python figure7/build_fig7_ver12.py
```

An isolated figure build reproduced the frozen SVG byte-for-byte and verified all three pooled aggregate count vectors. The builder checks preservation of the original plotted data and verifies that the PDF contains no raster images. The recorded result is in `validation/figure_replay.json`. The figure builder's temporary render checks and preview images are generated when the command runs; they are not duplicated in this archive.

## What is counted

One entry is a comparison between two assigned gene-expression events at one spatial stage for one landmark type. Each eligible entry is evaluated in the three spatial embryos at that stage. Figure 7B pools peaks, rising half-heights and falling half-heights. One gene/event pair at one stage can therefore contribute more than one landmark when each is eligible. The pooled count is the sum of the eligible peak and half-height counts, not an average of their percentages.

The current inventory contains 859 candidate entries covering all 21 gene pairs, with 233 eligible entries: 116 peaks and 117 half-heights. The pooled denominators are 33 pair-rule, 41 gap-gene and 159 gap / pair-rule comparisons. The 626 excluded entries remain in the S1 Data `Landmark inventory` sheet with their exclusion reasons. Pair-rule pairs contain two of eve, run and odd; gap-gene pairs contain two of Kr, gt, mlpt and svb; gap / pair-rule pairs contain one gene from each set.

Percentages use each bar's complete eligible denominator. No category is removed to increase a preservation or reversal percentage. Repeated comparisons share genes, events, temporal estimates and spatial embryos. The counts are descriptive summaries of the analyzed relationships, not independent biological replicates or formal tests of a difference between gene classes.

## Landmark definitions

Expression peaks are the maxima assigned to the original expression-event windows. The spatial position of a split eve family is represented by the midpoint of its two daughter peaks; the alternative daughter-span outcome is retained in the S1 Data `Landmark inventory` sheet. Pair-rule comparisons use corresponding cycles. Mixed gene pairs retain the assigned event inventory and, when eve is involved, the target gene's own matched temporal eve cohort.

A half-height is the profile minimum plus half the difference between that minimum and the assigned episode maximum. The temporal rising half-height is matched to the anterior spatial half-height, and the temporal falling half-height to the posterior spatial half-height. Borders describe the outer expression episode. The two gt peaks remain separate peak events but share one outer episode for border measurements. A half-height is a reproducible profile landmark, not a direct measurement of the first activation or final cessation of transcription.

The temporal profiles use stage-mean measurements at 24 ordered substages, full-series normalization, PCHIP interpolation to 800 positions, and a Savitzky-Golay filter with window length 81 and polynomial degree 3. Spatial profiles retain their recorded gene-, stage- and embryo-specific settings. The complete implementation and configuration are included under `../peak_order/engine/`.

## Eligibility and outcomes

The spatial stages are 3.3, 4.1, 4.3 and 5.2, corresponding to ordinal indices 8, 9, 11 and 13. A comparison requires identified, finite landmarks in all three spatial embryos within the verified posterior analysis region. Source event-identity exclusions and the existing late Kr/svb domain exclusions remain in force.

The same temporal availability criterion applies to every gene group and landmark: the upper marginal 95% bootstrap limit for each landmark must precede the current spatial-stage bin, expressed as `stage_index - 0.5`. This criterion prevents an unfinished temporal event from entering an earlier spatial comparison.

For each eligible entry, the five plotted categories are:

1. `preserved_all3`: the temporal difference interval excludes zero and all three spatial embryos preserve that order.
2. `reversed_all3`: the temporal difference interval excludes zero and all three spatial embryos show the opposite order.
3. `spatial_consistent_temporal_unresolved`: the temporal interval includes zero but all three spatial differences have the same nonzero sign.
4. `mixed`: temporal order is resolved but spatial ordering differs among embryos.
5. `other_unresolved`: the remaining eligible temporal or spatial ambiguities.

The third category is labeled **Unresolved in time, separated in space**. Its temporal difference interval includes zero, while the two spatial landmarks occupy distinct positions in the same relative direction in all three embryos. The label describes these criteria; it does not establish simultaneous temporal events or a minimum biological separation distance. Numerical zero uses a tolerance of 1e-12 embryo length. Missing landmarks remain missing, rather than becoming zero. Intervals are pointwise percentiles conditional on successful landmark recovery; recovery fractions are retained for every entry. Temporal coordinates are ordered-substage intervals, not elapsed minutes. The bootstrap samples embryos within their assigned stages and preserves co-measured channels together.

## Files and source lineage

- `build_dataset.py`: reconstructs the complete landmark inventory and temporal contrasts from the numerical inputs.
- `summarize_landmarks.py`: applies the five-category partition and calculates the inventory and both pooled and stratified summary tables. The 15 rows with `panel=all_landmarks` supply Fig. 7B. The remaining 30 rows retain separate peak and combined half-height summaries for supplementary reporting; they repeat the same observations by landmark stratum and must not be added to the pooled rows.
- `data/source_tables/`: the assigned peak and border positions, event identities, posterior masks and comparison source tables.
- `data/temporal_features_with_observed.npz`: 109 arrays, each with an observed estimate followed by 20,000 draws.
- `data/temporal_contrast_draws.npz`: 226 arrays of 20,000 direct B-minus-A temporal differences.
- `data/landmark_distances.csv`: the original complete comparison table. Its signed differences support the ordering categories.
- `data/dataset_manifest.json`: the source-table hashes and original code commit.
- `reference/`: the reviewed tables used for reproducibility checks.
- `DATA_DICTIONARY.md`: column definitions and units.
- `../../data/`: original temporal workbooks and spatial profile measurements, shared with the existing repository workflows.
- `../peak_order/`: original analysis code, curated event definitions, source hashes, requirements, fonts and numerical references.
- `figure7/`: the final figure builder, input tables, original editable artwork, source notes and frozen final SVG reference.
- `requirements-figure.txt`: Python dependencies for figure reconstruction. Inkscape is a separate system application.
- `validation/`: the independently checked numerical and figure reconstruction records.

Repository: https://github.com/ezzatelsherif/tribolium-sequential-imaging

Original peak-analysis commit: `a09fafea149c097516646d0ad8bec3d5b0ad8169`.

The current border arrays use the same paired-channel draw state as their associated peaks and eve references. The same cycle/trough extraction rule is used for matched eve borders across cohorts. The accompanying audit tables retain the checks linking these arrays to the earlier recorded landmark tables.

## Optional reconstruction from measured profiles

The routine command above is sufficient to reconstruct the submitted summary. To regenerate the original biological resampling from the measured profiles, install the full pinned upstream requirements and run:

```bash
python -m pip install -r ../peak_order/requirements.txt
python ../peak_order/run_analysis.py --skip-figures --output-dir full_reconstruction
```

This separate, longer run uses the original 20,000-draw implementation and its recorded seed. It checks the original peak-analysis reference tables. Its historical worksheet missing-value convention is retained, matching the saved draws used here. A separate confirmed-correction option remains available in the upstream code and source records.

For regeneration of the extended landmark cache, place copies of `build_dataset.py` and `summarize_landmarks.py` in a new directory that does not contain the saved `data/temporal_features_with_observed.npz`, then run its dataset builder with `--source-work /absolute/path/to/full_reconstruction/work --repository /absolute/path/to/repository`. This extracts the peak and border arrays from the newly generated paired-channel means using the recorded algorithms. Run the copied summary script afterwards. This optional fresh extraction was not repeated when preparing this package; the supplied cache-based reconstruction was executed and verified as described above.

The final Fig. 7 layout and editable artwork are also supplied with the submission's figure source files. This repository contains everything needed to reconstruct the numerical analysis, aggregate tables and updated Figure 7 after installing the stated dependencies.

This extension does not assign a new license to the original scientific data or code. The bundled Liberation Sans font files retain their existing license.


## Publication and reproducibility version

Cite the exact repository commit used for the final submission. The original
peak workflow remains unchanged; the new `landmarks` command adds the current
analysis and figure without replacing historical outputs. The command writes
to a fresh output directory so the reviewed source files are not overwritten.

The default command replays saved biological draws. This is sufficient to
reconstruct the submitted inventory, counts and figure. The optional full
regeneration instructions above describe the separate, longer embryo-resampling
workflow and subsequent extraction of the extended landmark cache.
