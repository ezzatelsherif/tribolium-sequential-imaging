# Tribolium sequential imaging

Data and analysis code for **High-multiplexity imaging reveals distinct modes of
translating time into space during embryonic patterning**, Garcia-Guillen,
Ahmadi, Frimpong and colleagues, El-Sherif laboratory.

This repository reproduces the numerical expression profiles, intron/exon
analyses, biochemical simulations, and the final second-revision temporal-to-spatial
landmark analysis, including peaks and rising/falling half-height borders. Its inputs
are the extracted fluorescence measurements used in the paper: 17 spatial
embryo CSVs, seven temporal workbooks, and the mapped probe-target information
in the intron/exon analysis package. Additional embryos and developmental
stages from the larger experimental collection are not included.

The workflow starts from measured fluorescence tables. Raw microscopy,
registration software, manual ImageJ ROIs, and the assembled microscopy figure
layouts are outside this release. Everything needed for the numerical analyses
below is included; no Drive access or external data download is needed.

## Quick start

Use **Python 3.12**. Clone this repository, open a terminal in its root directory,
and create a dedicated environment. The repository is public.

```sh
git clone https://github.com/ezzatelsherif/tribolium-sequential-imaging.git
cd tribolium-sequential-imaging
python -m venv .venv
```

Activate the environment on macOS or Linux:

```sh
source .venv/bin/activate
```

On Windows PowerShell, use `.\.venv\Scripts\Activate.ps1`. If environment
activation is restricted, use `.\.venv\Scripts\python.exe` in place of `python`
in the commands below. A short checkout path such as `C:\seq` avoids Windows
path-length problems.

```sh
python -m pip install -r requirements.txt
python reproduce.py verify
python reproduce.py profiles
python reproduce.py model --stage manuscript-figures
python reproduce.py landmarks
```

`verify` checks source checksums, historical manuscript curve matches, the
confirmed background correction, unique-cohort eve pooling, complete intron/exon
pair counts, target sequences, saved fits, and simulation summaries.

`profiles` writes CSV tables and PDF, SVG, and PNG plots to `outputs/profiles/`.
The model command writes Fig. 3 numerical components and supplemental model
figures (S1 and S2 in the final revised manuscript) to
`outputs/model/figures/`. These commands regenerate numerical panels with clear
labels; they do not reconstruct the complete microscopy composites.

`landmarks` reconstructs the complete 859-candidate landmark inventory and the
pooled and stratified results from the saved 20,000-draw biological bootstrap.
There are 233 eligible comparisons. The combined landmark results and original
peak audit are provided in **S1 Data**, a workbook with separate, clearly labeled
sheets. The CSV filenames containing `S2_Data` and `S3_Data` are retained as stable
internal filenames; they are no longer separate submission supplements.

To reproduce the final six-panel Fig. 7, install the additional plotting
dependencies and the Inkscape command-line application, then use a fresh output
directory:

```sh
python -m pip install -r analysis/landmarks/requirements-figure.txt
python reproduce.py landmarks --figure --output outputs/landmarks_figure
```

This rebuilds the editable single-layer SVG, vector PDF, 600 dpi TIFF and PNG
preview. The SVG is checked against the exact submitted artwork. See the
[landmark documentation](analysis/landmarks/README.md) for the workbook map,
source measurements, numerical checks, uncertainty definitions and full
regeneration from measured profiles.

`peaks` remains available as the historical peak-only workflow. It rebuilds the
original 20,000-draw analysis, its earlier Fig. 7 layout, the Fig. S6 companion
and the original peak audit. It does not generate the current six-panel Fig. 7.
The original input policies and correction sensitivity are documented in the
[peak-order documentation](analysis/peak_order/README.md).

## Recompute the analyses

Run these commands from the repository root. Model computations operate on a
copy in `outputs/model/`, leaving the archived analysis package intact. Use
`--output PATH` to keep a separate run. An existing output copy is reused, so
choose a new path when starting from a different repository version.

| Command | Result |
|---|---|
| `python reproduce.py verify` | Source integrity, 121 historical curve-shape checks and correction checks, and archived model verification |
| `python reproduce.py profiles` | Spatial and temporal profiles, stage summaries, and comparison panels |
| `python reproduce.py landmarks` | Current landmark inventory and pooled/stratified counts from saved biological draws, with reference checks |
| `python reproduce.py landmarks --figure` | Current six-panel Fig. 7 in SVG, PDF, TIFF and PNG, plus numerical checks |
| `python reproduce.py peaks` | Historical full 20,000-draw peak-order analysis, earlier Fig. 7, Fig. S6 companion, original peak audit, and regression checks |
| `python reproduce.py peaks --apply-confirmed-corrections --output outputs/peak_order_corrected` | Separate full sensitivity with the recovered late KR-cohort eve background |
| `python -m seqimaging.correction_review` | Effects of the recovered background and removal of repeated eve measurements |
| `python reproduce.py model --stage empirical` | Background correction and 5,000 paired within-stage bootstrap resamples |
| `python reproduce.py model --stage sweep` | All 1,944 hypothetical kinetic scenarios and numerical convergence checks |
| `python reproduce.py model --stage fits` | Joint fits for both initiation-pulse families and traversal modes |
| `python reproduce.py model --stage numerics` | Final-fit mesh checks and independent adaptive quadrature |
| `python reproduce.py model --stage assessment` | Traversal profiles, 199 conditional bootstrap refits per gene, optimizer checks |
| `python reproduce.py model --stage all` | All computational analysis stages in dependency order |
| `python reproduce.py model --stage manuscript-figures` | Fig. 3 numerical components and supplemental model figures (S1 and S2 in the final revised manuscript) from results in the output copy |
| `python reproduce.py model --stage movie` | Movie S1; also requires FFmpeg with the `libx264` encoder |

Verification and profile plotting are suitable starting points. Full fitting
and the bootstrap refits are substantially more expensive. The fit command
uses archived parameters as additional starts alongside seeded searches. Add
`--fresh` to `fits` or `all` to omit these warm starts. Different numerical
platforms or fresh optimization can reach different local minima. Random seeds,
model restrictions, and numerical qualifications are documented in the
[intron/exon package](analysis/intron_exon/README.md).

## Find the data and methods

| Location | Contents |
|---|---|
| [data/README.md](data/README.md) | Data dictionary and specimen identifiers |
| [data/manifest.csv](data/manifest.csv) | Exact source files, provenance links, sizes, checksums, and figure uses |
| [data/temporal_corrections.csv](data/temporal_corrections.csv) | Confirmed background correction and its provenance |
| [config/figures.json](config/figures.json) | Embryo assignments, plotted channels, and display settings |
| [docs/figure_map.md](docs/figure_map.md) | Manuscript panels, source files, commands, and output filenames |
| [docs/data_notes.md](docs/data_notes.md) | Normalization, unique-cohort eve pooling, recovered background, and source-version notes |
| [docs/methods.md](docs/methods.md) | Fluorescence processing and block-display calculations |
| [docs/validation.md](docs/validation.md) | What was checked and the limits of those checks |
| [seqimaging/](seqimaging/) | Profile readers, summaries, plotting, and verification |
| [analysis/intron_exon/](analysis/intron_exon/) | Original S1 Code package, source measurements, saved results, and executable notebook |
| [analysis/peak_order/](analysis/peak_order/) | Original peak-order code, curated event identities, reviewed tables, and full-run validation |
| [analysis/landmarks/](analysis/landmarks/) | Final peak and half-height analysis, saved bootstrap arrays, full inventory, pooled statistics, S1 Data workbook, and current Fig. 7 reconstruction |

The seven temporal workbooks retain the original measurements and cached
worksheet values. Three are also present in the unchanged S1 Code package;
verification checks that these copies are identical. The regenerated temporal
CSV retains the original values alongside the confirmed correction. Current
profile plots use corrected posterior-minus-background values and count the shared
run/odd eve cohort once. Historical cached values remain available for
comparison. The default peak-order run separately preserves the historical
missing-value pattern. See the [data notes](docs/data_notes.md) and the
[peak-order package](analysis/peak_order/README.md) for provenance and impact.

The original profile source-data selection and figure assignments were verified
against the first revised manuscript (Drive file
`1iWvbs7Pprpb97Rsep7eCD2Y3o_wvY7dX`, `Manuscript_ver51.docx`). The
underlying data and panel assignments are those of analysis commit
`19370cf7e16210569157a784997fbce0f1933fa0`. Between manuscript versions
46 and 51, the supplemental figure labels changed: former S3/S4/S5 became
S1/S2/S3, and former S1/S2 became S4/S5. The numerical workflows retain
their original output filenames; [the figure map](docs/figure_map.md)
identifies the corresponding final panels. The second revision adds Fig. 7,
Fig. S6 and S1 Data, audited against `Manuscript_ver57.docx` and
`S1_Text_ver03.docx` on 20 October 2026. These additions reuse 12 of the existing
spatial files and the seven existing temporal workbooks; they introduce no
additional biological specimens. The final landmark extension adds half-height
borders and the pooled Fig. 7 summary using the same measured specimens. Its
provenance and numerical reference checks are recorded in `analysis/landmarks/`.

## Attribution and citation

Spatial and temporal measurements and the original profile viewers were
provided by Theophilus Frimpong. This release supplies a command-line workflow
following those viewers' numerical transformations. See
[CITATION.md](CITATION.md) for manuscript attribution and upstream code links.
The intron/exon package from commit `4d3951c` is preserved unchanged.

Cite the manuscript and the exact repository commit used for an analysis.
No software DOI or open-source license has yet been assigned.
