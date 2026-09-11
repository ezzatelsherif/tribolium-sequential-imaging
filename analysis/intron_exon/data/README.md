# Input data and probe-target provenance

## Experimental measurements

`../empirical/temporal/KR.xlsx`, `MLPT.xlsx`, and `SVB.xlsx` are the source
workbooks supplied for the Figure 3 intron–exon profiles. `Sheet1` contains
24 stages, labeled 1.1 through 8.3, with five numbered specimen slots for each
channel. For each channel and slot, the corrected value is region **A minus B**.
The parser preserves negative corrected measurements. A slot is included in
paired calculations only when both corrected channel measurements are finite.

`background_corrected_pairs.csv` contains all 360 slots with columns `gene`,
`stage`, `stage_index`, `specimen_label`, `exon`, `intron`, and `complete_pair`.
The 24 integer indices order the stages; they are not elapsed minutes.
There are 116 complete Kr pairs, 118 mlpt pairs, and 116 svb pairs. Pair identity
is inferred from matched numbered slots, not independently verified identifiers.
Different stages use different embryos.

## Supplied probe-design targets

`target_sequences.fasta` contains only the six target sequences used in this
analysis. Whitespace was removed and RNA U was converted to DNA T. The duplicate
Kr intron entry in the source document was identical and included once.
These are **sequences used for probe design**, not individual manufactured HCR
oligonucleotide sequences. Their lengths are Kr exon 2,158 nt, Kr intron
2,879 nt, mlpt exon 616 nt, mlpt intron 569 nt, svb exon 3,421 nt, and svb intron
4,368 nt. The original source-document SHA256 is
`4e2904fa329da5b4dc942a177296ee33faff4844817a2c1c6fb299dd9a1f7025`.

`reference_spans.fasta` contains the Tcas5.2 genomic span used for each gene,
oriented in transcription direction. Sequence zero is the aligned exonic
target's 5-prime boundary. `geometry.json` uses **zero-based, half-open** target
intervals relative to that origin. Chromosomal coordinates in the provenance
and alignment tables are **one-based, inclusive**. On the minus strand, the
oriented span is reverse complemented; this is the case for mlpt.

| Gene | Reference accession | Chromosomal span | Strand | Modeled span |
|---|---|---|---|---|
| Kr | NC_007425.3 | 6,495,452–6,502,496 | + | 7,045 bp |
| mlpt | NC_007424.3 | 15,958,600–15,959,870 | − | 1,271 bp |
| svb | NC_007419.2 | 5,009,353–5,018,302 | + | 8,950 bp |

`target_alignments.tsv` retains the six recorded Tcas5.2 alignments, including
query coverage, matches, strand, and CIGAR. `reference_provenance.tsv` gives
versioned accessions, original public NCBI retrieval URLs, and SHA256 hashes of
the retrieved GenBank files and oriented modeled sequences. No network access
is needed to verify target identity against these archived spans.

For each intronic target, the selected interval reproduces the complete supplied
sequence exactly. Concatenating each gene's exonic intervals likewise reproduces
the exonic target, except that the supplied mlpt exon target has 17 terminal A
residues not included in the genomic mapping. These residues are retained in
the supplied query FASTA but excluded from signal weights and coordinates.

## Mapping qualifications

- Model start and completion positions are target-based proxies, not measured
  transcription start and termination sites.
- The Kr modeled Tcas5.2 span contains 984 N bases within the intron locus;
  the supplied target bases themselves match exactly. These Ns contribute to
  the modeled span length, so reference-gap uncertainty is retained.
- The supplied svb exonic sequence matches the genomic exon chain in
  `geometry.json`, but that chain differs from the RefSeq transcript models
  checked in the preceding mapping work. Canonical splice junctions do not
  establish which isoform is expressed in these embryos.
- Alignments were performed within the downloaded candidate loci. They do not
  establish genome-wide probe specificity.
- Target-length fractions and uniform-per-base modeling weights are not
  measurements of manufactured probe counts or fluorescence contribution.
