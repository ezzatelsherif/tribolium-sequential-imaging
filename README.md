# Tribolium sequential imaging

Code and reproducibility materials for the sequential imaging paper on gene
expression dynamics in *Tribolium castaneum* from the El-Sherif laboratory.
This repository is organized to hold the paper's analysis and imaging code.

## Current contents

`analysis/intron_exon/` contains the finalized **S1 Code** package accompanying
the revised manuscript. It includes source measurements, mapped probe target
regions, code, notebooks, recorded results, and checksums for:

- Experimental intronic–exonic profile offsets with bootstrap uncertainty.
- The full biochemical model sweep using hypothetical transcription pulses.
- Constrained joint fits to the measured intronic and exonic profiles.
- Numerical checks and scripts for the corresponding figures and movie.

The existing package is preserved unchanged, including its relative paths and
release manifest. Its [README](analysis/intron_exon/README.md) documents the
methods, inputs, limitations, and complete reproduction commands. See its
[CITATION.md](analysis/intron_exon/CITATION.md) for attribution and
[references.ris](analysis/intron_exon/references.ris) for methodological references.

## Run the intron–exon analysis

Use Python 3.12 in a dedicated environment, then run:

```sh
cd analysis/intron_exon
python -m pip install -r requirements.txt
python reproduce_analysis.py --stage verify
```

The verification command checks the saved package and results without rerunning
the long optimizations. See the package README for the full analysis and optional
figure-generation dependencies. On Windows, clone to a short path such as
`C:\seq` to avoid path-length problems.

## Adding the paper's other code

Add other analyses in named subdirectories under `analysis/`, and imaging or
registration workflows under `imaging/`. Keep each workflow's dependencies,
input provenance, run instructions, and relationship to the manuscript together
in its own README. Other paper workflows have not yet been added to this initial
package. Original manuscript and microscopy/Illustrator figure files are not
included here.

No software DOI or open-source license has been assigned to this repository.
