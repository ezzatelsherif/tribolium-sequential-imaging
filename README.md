# Tribolium sequential imaging

Data and analysis code for **High-multiplexity imaging reveals distinct modes of
translating time into space during embryonic patterning**, Garcia-Guillen,
Ahmadi, Frimpong and colleagues, El-Sherif laboratory.

This repository reproduces the numerical expression profiles, intron/exon
analyses, and biochemical simulations in the revised manuscript. Its inputs
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
and create a dedicated environment. While the repository is private, cloning
requires a GitHub account with access.

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
```

`verify` checks source checksums, historical manuscript curve matches, the
confirmed background correction, unique-cohort eve pooling, complete intron/exon
pair counts, target sequences, saved fits, and simulation summaries.

`profiles` writes CSV tables and PDF, SVG, and PNG plots to `outputs/profiles/`.
The last command writes Fig. 3 numerical components and Figs. S3 and S4 to
`outputs/model/figures/`. These commands regenerate numerical panels with clear
labels; they do not reconstruct the complete microscopy composites.

## Recompute the analyses

Run these commands from the repository root. Model computations operate on a
copy in `outputs/model/`, leaving the archived analysis package intact. Use
`--output PATH` to keep a separate run. An existing output copy is reused, so
choose a new path when starting from a different repository version.

| Command | Result |
|---|---|
| `python reproduce.py verify` | Source integrity, 121 historical curve-shape checks and correction checks, and archived model verification |
| `python reproduce.py profiles` | Spatial and temporal profiles, stage summaries, and comparison panels |
| `python -m seqimaging.correction_review` | Effects of the recovered background and removal of repeated eve measurements |
| `python reproduce.py model --stage empirical` | Background correction and 5,000 paired within-stage bootstrap resamples |
| `python reproduce.py model --stage sweep` | All 1,944 hypothetical kinetic scenarios and numerical convergence checks |
| `python reproduce.py model --stage fits` | Joint fits for both initiation-pulse families and traversal modes |
| `python reproduce.py model --stage numerics` | Final-fit mesh checks and independent adaptive quadrature |
| `python reproduce.py model --stage assessment` | Traversal profiles, 199 conditional bootstrap refits per gene, optimizer checks |
| `python reproduce.py model --stage all` | All computational analysis stages in dependency order |
| `python reproduce.py model --stage manuscript-figures` | Fig. 3 numerical components and Figs. S3 and S4 from results in the output copy |
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

The seven temporal workbooks retain the original measurements and cached
worksheet values. Three are also present in the unchanged S1 Code package;
verification checks that these copies are identical. The regenerated temporal
CSV retains the original values alongside the confirmed correction. Current
plots use corrected posterior-minus-background values and count the shared
run/odd eve cohort once. Historical cached values remain available for
comparison. See the [data notes](docs/data_notes.md) for provenance and impact.

The source-data selection and figure assignments correspond to
`Manuscript_ver46.docx` in the September 14 revision folder. Regenerated
numerical panels incorporate the corrections; earlier assembled microscopy
figures require the corresponding temporal curves to be replaced.

## Attribution and citation

Spatial and temporal measurements and the original profile viewers were
provided by Theophilus Frimpong. This release supplies a command-line workflow
following those viewers' numerical transformations. See
[CITATION.md](CITATION.md) for manuscript attribution and upstream code links.
The intron/exon package from commit `4d3951c` is preserved unchanged.

Cite the manuscript and the exact repository commit used for an analysis.
No software DOI or open-source license has yet been assigned.
