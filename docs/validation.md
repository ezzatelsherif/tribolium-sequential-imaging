# Reproduction checks

The release was checked on 2026-09-24 with Python 3.12.14 on Linux, using the
versions in the root `requirements.txt`: NumPy 2.3.5, SciPy 1.17.0, pandas 2.2.3,
Matplotlib 3.10.8, openpyxl 3.1.5, and PyMuPDF 1.26.6. Windows and macOS commands
are documented, but this release was not executed on those systems.

## Source and figure checks

`python reproduce.py verify` checks all 24 source-file SHA256 hashes and
requires every included spatial CSV to have a figure assignment. The three
temporal workbooks used by both workflows are checked for byte identity.
Raw corrections, complete intron/exon pair counts, shared run/odd eve entries,
and the single missing-background exception are checked explicitly.

The file `validation/figure_traces.json` contains 121 sampled curve shapes
from the supplied manuscript figure PDFs: 114 spatial curves spanning all
17 included embryos, and seven temporal comparison curves (Kr, mlpt, gt, svb,
run, odd, and pooled eve). Each reference records the figure and vector-path
index used. Repeated identical curves are represented once per source/channel.

The reference and regenerated curves are mapped to their own 0-1 x and y
ranges and compared at 120 positions. The maximum observed root-mean-square
error was **0.000420**, below the **0.002** acceptance tolerance. The tolerance
allows for PDF coordinate rounding and resampling. This comparison establishes
source identity, relative curve shape, and the smoothing choice. It cannot
distinguish fluorescence baselines related by an affine transformation,
establish a physical axis calibration, or validate the underlying biology.
The min-max normalization convention is established from the source viewer.

| Reference name in the trace file | Supplied figure source |
|---|---|
| `Figure_2.pdf` | `Fig 2 - Methodology Intro ver01 copy.pdf` |
| `Figure_3.pdf` | `Figure_3.pdf` in the revision figures folder |
| `Figure_6.pdf` | `Fig 6 - Comparing time to space - ver09.pdf` |
| `Figure_S1.pdf` | `Fig S1 - Comparing time to space - ver03 - replicate 2 copy.pdf` |
| `Figure_S2.pdf` | `Fig S2 - Comparing time to space - ver03 - replicate 3 copy.pdf` |
| `Fig 4 - Comparing_time_with_Space_INTRONIC_ver04.pdf` | Same filename in the revision figures folder |
| `Fig_S5_Spatial_intronic_exonic_replicates_ver03.pdf` | Same filename in the revision figures folder |

Only numerical line traces are stored here; microscopy images and PDF
composites are not required to rerun the checks. The Fig. 5 cohort summaries
follow the workbook measurements and manuscript normalization, and are not
included in the independent PDF curve-shape check. Block displays follow the
original viewer's algorithm but are not tested by pixel comparison against
the assembled figure.

## Intron/exon package checks

The original package verifier confirmed:

- All 99 files in its release manifest match their recorded hashes.
- The six probe-design target sequences match the reference intervals, with
  the documented terminal-A exclusion for mlpt.
- Background corrections match all 360 source slots, with 116 complete pairs
  for Kr, 118 for mlpt, and 116 for svb.
- All 12 stored fits reconstruct their archived predictions, with maximum
  prediction difference zero in this environment.
- The sweep summaries and three independently recomputed scenarios agree.
- Conditional lack-of-fit summaries agree with the recorded null draws.

The 5,000-resample empirical bootstrap and the complete 1,944-scenario sweep
were rerun from the included inputs. Corrected pairs, bootstrap estimates,
sensitivity results, and scenario tables agreed with their archived tables
using relative tolerance `1e-11` and absolute tolerance `1e-12`.
The model consistency verifier also passed on the rerun copy using
`--skip-manifest`, which permits regenerated file bytes to differ.

Final-fit mesh checks and the independent adaptive-quadrature calculations
were rerun. Numerical manuscript panels, the profile plots and CSV exports,
and Movie S1 were regenerated successfully. Representative spatial and
temporal comparison plots were inspected visually.

All saved optimizations and recorded conditional bootstrap fits are included
and verified for internal consistency. The full fresh global-optimization
search and all 199 bootstrap refits per gene were **not** repeated as part of
this packaging check. Commands for those computations are in the main README.

## Limits

Reproducibility checks do not resolve the source and manuscript qualifications
in [data_notes.md](data_notes.md). In particular, the display's cached
missing-background value, duplicated eve weighting, spatial-normalization
wording, and Fig. 6B/C legend mismatch remain explicit. The original source
measurements have not been altered to remove these qualifications.
