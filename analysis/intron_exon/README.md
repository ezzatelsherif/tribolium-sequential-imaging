# Intron–exon signal analysis for sequential imaging

Reproducible source code, input data, and recorded results for the intron–exon
analysis of *Tribolium castaneum* **Kr**, **mlpt**, and **svb**. This package is
**S1 Code** accompanying the revised manuscript. The methods are described in
the manuscript and **S1 Text**. No repository URL or DOI has been assigned here.

The analysis has three separate parts:

1. Experimental rising and falling profile offsets, with a paired within-stage
   bootstrap using 5,000 resamples.
2. A hypothetical-pulse sweep of 1,944 parameter combinations using mapped probe
   design target regions and transcription, processing, and persistence dynamics.
3. Joint fits to both experimental channels using a constrained shared
   transcriptional initiation input, with numerical and descriptive fit checks.

## Start here

Open `analysis.ipynb` in Jupyter to inspect the recorded results. The files can
also be used directly from a terminal. Extract to a short path, for example
`C:\model`, to avoid Windows path-length problems.

The recorded environment is Python 3.12.14 with the package versions pinned in
`requirements.txt`; `environment.json` records it in machine-readable form.
Create a Python environment, activate it, and run from this directory:

```text
python -m pip install -r requirements.txt
python reproduce_analysis.py --stage verify
```

The verification checks release hashes, source workbook corrections, all 12
stored fits, sweep summary counts, three rerun scenarios, and the conditional
fit-assessment summaries. It does not rerun the long optimizations or establish
biological validity. Jupyter is optional and is not required for command-line
reproduction; install `jupyterlab` separately if a notebook interface is desired.

## Reproduction

| Command | Operation |
|---|---|
| `python reproduce_analysis.py --stage verify` | Verify the unmodified release and saved results. |
| `python reproduce_analysis.py --stage figures` | Regenerate the analysis plots from saved results, without refitting. |
| `python reproduce_analysis.py --stage empirical` | Recompute paired offsets and corrected measurements from the three workbooks. |
| `python reproduce_analysis.py --stage sweep` | Recompute the full hypothetical sweep, summaries, numerical checks, and sweep plots. |
| `python reproduce_analysis.py --stage fits` | Refit both pulse families and both traversal modes for every gene. |
| `python reproduce_analysis.py --stage numerics` | Regenerate final-fit mesh and independent adaptive-quadrature checks. |
| `python reproduce_analysis.py --stage assessment` | Refit traversal profiles and 199 model-centered residual-bootstrap draws per gene, including optimizer checks. |
| `python reproduce_analysis.py --stage all` | Recompute all analysis stages in dependency order, then verify internal consistency. |

Running analysis stages overwrites derived files **in this extracted copy**.
Original measurements in `empirical/temporal/` and target regions in
`data/geometry.json` are never written by these stages. For a fresh independent
optimization without archived parameter warm starts, add `--fresh` to `fits`
or `all`. The default uses stored parameters as additional starts, as well as
the specified seeded global searches; this is useful for recovering the archived
minima. A fresh optimization may converge to a different local solution.

The sweep and full-fit/assessment stages can take substantially longer than
verification or plotting. All central random seeds are specified in the source.
Numerical minima and figure file bytes can vary by operating system and software
versions. This package does not promise bit-for-bit regeneration of optimizer
metadata, timestamps, or SVG identifiers. After an intentional rerun, use
`python verify_saved_results.py --skip-manifest` to check scientific consistency
without comparing changed files against the release checksums.

## File map

| Path | Contents |
|---|---|
| `data/geometry.json` | Fixed target coordinates and reference/annotation qualifications. |
| `data/target_sequences.fasta` | Six supplied sequences used for probe design, in DNA alphabet. |
| `data/reference_spans.fasta` | Corresponding Tcas5.2 modeled spans, oriented in transcription direction. |
| `data/reference_provenance.tsv`, `data/target_alignments.tsv` | Accessions, coordinates, sequence hashes, and recorded alignment evidence. |
| `data/README.md` | Input provenance and definitions, including all coordinate conventions. |
| `empirical/temporal/` | Source intensity workbooks, with target-region A and background-region B measurements. |
| `data/background_corrected_pairs.csv` | All 360 numbered slots after A-minus-B correction, including incomplete pairs. |
| `empirical/original_stage_lag_analysis.py` | Bootstrap analysis and source-workbook parser. |
| `empirical/stage_lag_*.csv` | Recorded empirical intervals and threshold/interpolation sensitivity. |
| `position_kernels.py` | Analytic single-transcript age-response kernels and integrals. |
| `sweep.py`, `sweep_results/` | Hypothetical scenarios, classifications, summaries, and numerical checks. |
| `pulse_forward.py`, `pulse_fit.py`, `fitresults/` | Convolution, constrained joint fitting, stored fits, and final-fit numerical checks. |
| `assess_constrained.py`, `assessment/` | Conditional model-centered residual bootstrap, traversal profiles, and optimizer checks. |
| `analysis_protocol.md`, `sweep_design.md`, `fitdesign.md` | Recorded analysis design and fixed restrictions. |
| `fitting_methods_notes.md`, `assessment/README.md` | Implementation details and assessment interpretation. |
| `provenance/` | Historical coarse-grid/refinement logs, not new independent analyses. |
| `manifest.json` | SHA256 checksums of this release. |
| `references.ris` | Five added scientific/methodological references as importable RIS metadata. |

Standalone analysis figures are supplied as PNG and editable SVG and regenerated
by `make_sweep_figures.py`, `make_fit_figures.py`, and `make_empirical_figure.py`.
`make_manuscript_figures.py` supplies the numerical panels arranged for the
revised manuscript, S3 Fig, S4 Fig, and the replacement S1 Movie. Install its
additional PDF-layout dependency with:

```text
python -m pip install -r requirements-figures.txt
python reproduce_analysis.py --stage manuscript-figures
```

These outputs are written to `figures/`. The script needs no original microscopy
or Illustrator files in this mode. To reassemble the entire Fig 3 composite,
provide the PDF-compatible original as
`python make_manuscript_figures.py --original-figure PATH_TO_ORIGINAL.ai`.
That original is read without alteration and is not distributed in S1 Code.

For the movie, install FFmpeg with the `libx264` encoder available on the system
PATH, then run `python reproduce_analysis.py --stage movie`. The manuscript
figure and movie stages are optional presentation steps and are not invoked by
the core `--stage all` command. `--output-dir PATH` on the figure script directs
its output elsewhere. The PDF-layout dependency used here is PyMuPDF 1.26.6;
the original numerical-analysis environment remains recorded separately.

## Interpretation limits

- Times in the hypothetical sweep are ratios to pulse duration. The fit and
  empirical analyses use consecutive equally spaced eve–en substages, without
  absolute time calibration. Fitted waits are not measured biochemical lifetimes.
- Target intervals describe design sequences, not manufactured HCR oligo sites
  or measured fluorescence efficiencies. Signal is uniform per modeled target
  base. Target boundaries proxy transcription-unit boundaries.
- The model includes nascent exonic signal, processing of the targeted intron,
  and subsequent exonic RNA persistence. Other introns, residual excised-intron
  signal, export, and nonlinear detection are not separately modeled.
- Scenario percentages describe this specified parameter grid, not probabilities
  of biological mechanisms. Half-height ordering is distinct from molecular onset.
- Pairing follows numbered intron/exon slots in each workbook stage. Independent
  specimen identifiers were not available. Complete-pair totals are 116, 118,
  and 116 for Kr, mlpt, and svb. Stages contain different embryos.
- Empirical percentile intervals condition on assigned stages, slot pairing, and
  pulse windows. They omit stage-assignment and batch uncertainty and are not
  corrected for multiple comparisons.
- Fits minimize unweighted joint error in 48 channel-stage means, after fixed
  channel-maximum scaling. Restricting the initiation input to one/two smooth
  pulses reduces flexibility but does not provide held-out predictive validation.
- The primary constrained model retains a clear mismatch for mlpt. Biochemical
  parameter trade-offs, pulse-shape dependence, and boundary hits must remain
  visible when discussing fit quality. The conditional fit-assessment p values
  are approximate; 0.005 is the minimum resolution with 199 draws.

## Citation and reuse

See `CITATION.md` for attribution. A software DOI and public repository URL have
not been assigned in this package. No open-source license has been selected;
the package must not be represented as carrying an MIT, BSD, GPL, or other
license that its authors have not specified. Source measurements and original
code are provided for inspection and reproducibility of the accompanying study.

`references.ris` supplies the five added references from the accompanying
manuscript in a format importable by Sciwheel and other reference managers. It
does not invent Sciwheel database identifiers or represent newly imported
records as already linked Word citation fields.
