# Temporal and spatial peak order

This package reproduces the peak-order analysis accompanying the second revision of **High-multiplexity imaging reveals distinct modes of translating time into space during embryonic patterning**. It starts from the measured spatial-profile CSVs and temporal Excel workbooks, regenerates the 20,000 biological-bootstrap draws, reapplies the documented smoothing, identifies peaks, and produces Figure 7 and the complete comparison table.

The question is whether the order of two expression peaks in developmental time agrees with their anterior-to-posterior order in space. Gene pairs are screened without assigning an expected outcome from their traditional gap or pair-rule classification. A split eve domain is represented by its two observed daughter maxima and their midpoint. The midpoint is a domain centre, not an observed or reconstructed presplitting peak.

## Reproduce the results

Tested with Python 3.12.14. Use a new virtual environment and the pinned dependencies. Run these commands from the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r analysis/peak_order/requirements.txt
python analysis/peak_order/run_analysis.py --output-dir peak_order_results
```

On Windows, replace the activation command with `.venv\Scripts\activate`. The commands use only local inputs after dependency installation. No Google Drive, handover notebook, private file path, ChatGPT runtime, Illustrator installation, or saved bootstrap cache is required.

`--data-root /path/to/data` can select a different checkout’s repository data folder. Input files are byte-identical to the reviewed source measurements; their hashes are recorded in `inputs/scientific_input_manifest.json`.

The default command performs the full **20,000-draw** analysis with seed **20261019**. It refuses to overwrite a nonempty output directory. Intermediate files are written under `peak_order_results/work`, not into the repository source files. Logs for every stage are under `peak_order_results/logs`.

For a quicker installation check, use a smaller bootstrap count and a separate output directory:

```bash
python analysis/peak_order/run_analysis.py --bootstraps 100 --output-dir peak_order_smoke
```

The smaller run checks execution only; its intervals are not the reported results. `--skip-figures` runs the numerical analysis without figure rendering. Reference comparisons are enabled for the default full run.

## Main outputs

| File | Content |
| --- | --- |
| `Figure_7.pdf` | Vector PDF of the main peak-order figure |
| `Figure_7.svg` | Editable single-layer SVG with live text and vector paths |
| `Figure_7.png` | Figure preview |
| `Peak_Order_Companion.pdf` | Seven-page Fig. S6 companion: temporal cohorts and all four spatial stages, with three embryos separate |
| `S1_Data_Peak_Order_Comparisons.csv` | Complete distinct-gene comparison inventory, including exclusions and qualifications |
| `Gene_Pair_Coverage.csv` | Coverage of all 21 distinct gene pairs |
| `Posterior_Domain_Mask.csv` | Posterior family identities, ranks, and verified analysis limits |
| `Fig7_peak_time_differences_ver08.csv` | Six direct temporal contrasts displayed in Figure 7 |
| `Fig7_peak_markers_ver08.csv` | Every Figure 7 peak or domain-centre marker and temporal interval |
| `run_manifest.json` | Input hashes, runtime, stage timings, and comparisons against reviewed source tables |

The full reconstruction completed in approximately 95 seconds in the validation environment; runtimes depend on hardware. All eight numerical reference tables passed (including 1,167 comparison rows) at the stated tolerance. The recorded run is in `validation/historical_run.json`.

The companion matches the seven-page submitted Fig. S6 in numerical content and text, and applies the current Figure 7 gene palette. The submitted historical Fig. S6 used the older palette; the regenerated companion is the colour-consistent replacement.

The historical version numbers in `engine/` and `work/` preserve the provenance of the reviewed analysis. They do not require earlier versions to be installed. The original filenames are retained so the implementation can be compared with the audit history. The named files in the output directory are the publication-facing results.

## Data and pipeline

1. **Temporal data.** The repository’s `data/temporal/*.xlsx` contains the measured intensity and background for five embryo slots at each of 24 ordered substages. The analysis uses measured intensity minus background from `Sheet1`; missing slots remain missing, not zero. Cached `BASE` values are checked against the raw difference. The run and odd workbooks contain identical eve measurements, so those three channels form one jointly resampled cohort. Target/eve channels measured together stay paired. GT, MLPT, KR and SVB are independently sampled cohorts. The HB workbook is retained because it was part of the source cohort inventory; HB is not a gene in the peak-order screen.
2. **Biological resampling.** Embryos are resampled within each assigned stage, with all available channels for the same embryo resampled together. The observed stage means occupy draw zero; draws 1 through 20,000 are bootstrap means. There is no resampling of spatial pixels or treatment of expression cycles as biological replicates.
3. **Spatial data.** The primary analysis uses embryos 1, 2 and 3 separately at stages 3.3, 4.1, 4.3 and 5.2. Only the 12 primary spatial CSVs are used; the earlier unverified EMB4 audit is omitted. The source CSVs are read from `data/spatial/exonic/`. Whole-embryo anterior-to-posterior coordinates are scaled to 0–1 before smoothing. No cross-embryo peak alignment or spatial averaging is used in Figure 7.
4. **Smoothing.** The shared helper and `original_smoothing_config.json` preserve the original temporal PCHIP interpolation to 800 positions followed by Savitzky–Golay filtering with window 81 and polynomial degree 3. Spatial profiles use the original gene/stage/embryo fractional window before cropping; the full native sample count determines the filter width. The original source has no svb-specific spatial setting, so svb uses the disclosed default extension, fraction 0.10 and degree 3. The package reruns its sensitivity check at fractions 0.05, 0.10, 0.15 and 0.20.
5. **Episode identity and peak extraction.** `inputs/reconstructed_points.json` contains the original coarse episode configuration and its source hash. These windows identify corresponding expression episodes; measured maxima are re-extracted from the curves. Old hand-entered landmarks are not used as the final peak coordinates. Temporal eve/run/odd searches within a cycle use the same window, extending half a sampled substage beyond the first and last sample of that cycle. Unresolved maxima remain flagged.
6. **Posterior restriction.** The three most posterior eve domain families define the permitted region, counting a split family once. At stages 4.1 and 4.3 the families are cycles 4, 3 and 2; at 5.2 they are 5, 4 and 3. Only two source-identified families are available at 3.3. The unresolved third-family boundary is recorded rather than invented. Developing and unresolved families retain their posterior rank. The primary split rule requires two daughter maxima with at least 10% episode-range prominence and 4% embryo-length separation; a 15% prominence sensitivity is retained.
7. **Comparison.** Paired temporal peak differences are calculated within each bootstrap draw, then summarized by their pointwise 2.5th and 97.5th percentiles. Independently sampled cohorts are compared using their independently resampled draws and assigned stage coordinates. Spatial differences are calculated within each embryo. The complete inventory retains excluded, uncertain and centre-sensitive cases.

## Confirmed late-stage background correction

The historical peak-order run used the original workbook values and therefore omitted the KR-cohort eve observation at stage 8.1, slot 5, whose background was originally missing. The main repository profile workflow separately applies the recovered background of 328 (posterior 1380; corrected intensity 1052), documented in `data/temporal_corrections.csv`. The default command deliberately reconstructs the reviewed peak-order analysis, including its original missing-value pattern. It does not silently reinterpret that run as corrected.

To recompute the complete analysis with the confirmed correction, use a separate output folder:

```bash
python analysis/peak_order/run_analysis.py --apply-confirmed-corrections --output-dir peak_order_corrected
```

To compare both completed runs, use:

```bash
python analysis/peak_order/compare_correction_runs.py peak_order_results peak_order_corrected --output-dir peak_order_correction_audit
```

This sensitivity run retains the original algorithm and seed but applies the recovered observation before resampling. Adding a fifth complete paired observation at stage 8.1 also advances the original shared random-number stream differently for subsequent cohorts. Bootstrap-interval differences can therefore include Monte Carlo variation even where biological inputs are unchanged. Historical reference-table equality is disabled only for this explicitly corrected mode. The correction source is verified against the original workbook, and the applied value is recorded in `work/all_gene_direct_analysis/bootstrap_audit.json`. `temporal_raw_tidy.csv` retains original workbook-derived values; stage summaries and bootstrap curves reflect the explicitly selected correction mode.

The full 20,000-draw correction sensitivity preserved every observed peak coordinate, all six Figure 7 direct contrasts and their intervals, and every outcome/eligibility classification in the 1,167-row comparison inventory. Six of the 74 broader gap/mixed contrasts had interval endpoints move by at most 0.028786 ordinal substage steps; one svb marginal peak-time upper endpoint also shifted by that amount. These changes occur in cohorts following KR in the shared random-number stream. The 66 pair-rule contrasts were unchanged. Before substage index 15, raw smoothed profiles were identical; full-series rescaling changed the KR-cohort eve display by at most 0.007529 on its 0–1 scale. The machine-readable comparison and both run manifests are in `validation/`.

## What the uncertainty means

The temporal intervals quantify biological sampling under the assigned stages, the chosen smoother, and successful identification of the specified peaks. They do not include unknown embryo ages, mis-staging, unequal substage durations, a systematic offset between independently staged cohorts, or smoothing-choice uncertainty. Coordinates are **ordered substage steps, not minutes**. These are fixed-embryo reconstructions, not repeated observations of one embryo.

Intervals are conditional on valid peak recovery, so the recovery fraction accompanies each result. The odd cycle-3 peak is recovered in 86.99% of the 20,000 draws; its interval is conditional on those draws. The intervals are pointwise and exploratory, without correction for testing many related comparisons. The 1,167 rows in the combined inventory share genes, events and embryos, so they are not 1,167 independent biological tests. Agreement across the three observed spatial embryos is replicate consistency, not strong population-level significance.

## Interpretation of the displayed examples

- The run-before-odd relationship in cycle 2 is a preserved actual-peak example in all three embryos at stages 3.3, 4.1 and 4.3.
- svb follows its same-cohort eve cycle-3 peak temporally, but its spatial maximum is anterior to eve in all three embryos at 4.1 and 4.3. It is anterior to both daughters when eve is split.
- The smoothed GT/MLPT comparison has temporal order MLPT, GT1, GT2 and spatial order GT1, MLPT, GT2 in all three stage-4.3 embryos. The temporal cohorts are independent, and the first separation is conditional on smoothing and stage alignment. In the unsmoothed stage means, GT1 and MLPT both peak at stage 3.3.
- The cycle-3 eve-centre, run, odd spatial order in stage-4.3 embryo 1 depends on the chosen eve centre because run lies between the daughters. The complete lobe-span and split-rule checks retain that qualification.
- The Kr/eve cycle-3 temporal difference includes zero. Their temporal order remains unresolved, even though their spatial positions are distinct.

These are event- and stage-specific statements. Peak order alone does not establish whole-profile identity, complete domain containment, or a universal distinction between gap and pair-rule systems. No Damerau–Levenshtein score is used in this package.

## Implementation provenance and validation

`source_code_provenance.json` identifies the reviewed source modules and their original hashes. The historical numerical functions are preserved. Packaging changes comprise orchestration, an explicitly separate confirmed-correction sensitivity mode, skipping unavailable historical raw-versus-smoothed diagnostic tables, portable bundled fonts, and using PyMuPDF for preview rendering. The table-only assembler omits an obsolete narrative report. The current figure export uses the submitted gene palette with the requested darker gt yellow and individually colour-matched gene names.

`reference/` contains the reviewed numerical tables used by the full-run regression check. Paths that merely identify an output directory are ignored; table values, flags and dimensions must match to a tolerance of 1e-12. The full run rebuilds inputs, means, smoothing and features before making that comparison. It does not copy reference tables into the results.

This package does not assign a new license to the project code or scientific measurements. The bundled Liberation Sans files retain their SIL Open Font License in `fonts/LICENSE_LIBERATION`.
