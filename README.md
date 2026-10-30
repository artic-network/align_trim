# align_trim

Stand alone version of ARTIC's fieldbioinformatics align_trim.py

## Installation  

From conda
```bash
conda install bioconda::align_trim 
```
from pypi
```bash
pip install align_trim
```
from source
```bash
git clone https://github.com/artic-network/align_trim.git
cd align_trim
uv sync
uv run align_trim --help
```

## Command Line Interface

### Basic Usage

```bash
align_trim [OPTIONS] BEDFILE
```

The tool reads alignment data from either a SAM/BAM file or stdin and outputs trimmed alignments to stdout in SAM format by default.

### Required Arguments

- `BEDFILE`: BED file containing the amplicon primer scheme in [v3](https://doi.org/10.5281/zenodo.16366659) format. 

### Optional Arguments

#### Input/Output Options

- `--samfile`, `-s` : Sorted SAM/BAM file containing the aligned reads, if this is not provided (or '-') then 'align_trim' will read from stdin.
- `--output`, `-o` : Output file path. Format determined by extension (.sam/.bam). If not provided or '-', writes SAM to stdout

#### Processing Options

- `--normalise`, `-n` : Normalise to target depth N per amplicon using a greedy per-read algorithm. Each read (or read pair) is kept only if it brings the amplicon depth closer to the target. Depth is counted in fragments: a read pair counts once where its mates overlap, so N means N fragments for both single-end and paired-end data. Use 0 for no normalisation (default: 0)
- `--min-mapq`, `-m` : Minimum mapping quality to keep an aligned read (default: 20)
- `--primer-match-threshold`, `-p` : Add this many bases of padding to the 5' end of primer coordinates to allow fuzzy matching for reads with barcodes/adapters (default: 35)

#### Primer and Read Handling

- `--no-trim-primers` : Do not trim primers from reads (by default, primers are trimmed)
- `--allow-incorrect-pairs` : Allow reads to be assigned to amplicons even if primers are not correctly paired
- `--require-full-length` : Require all reads to start and stop in primer sites (do not use with rapid barcoding)

#### Output and Reporting

- `--report`, `-r` : Output detailed report TSV to specified filepath
- `--amp-depth-report`, `-a` : Output mean fragment depth for each amplicon as TSV to specified filepath
- `--genome-coverage-report`, `-g` : Output per-position genome coverage TSV(s) using the given prefix. Summary statistics (% genome covered at >=1x, >=10x, >=20x, >=100x) are printed to stderr. See [Genome Coverage Report](#genome-coverage-report) for details
- `--qc-prefix` : Output amplicon QC reports as `PREFIX.<name>.tsv`, with a one line summary per report printed to stderr. See [QC Reports](#qc-reports) for details
- `--no-read-groups` : Do not divide reads into pool-based read groups in SAM/BAM output

#### General Options

- `--verbose`, `-v` : Enable debug mode with detailed logging to stderr
- `--version` : Show version information
- `--help` : Show help message

### Examples

#### Basic trimming with primer removal
```bash
align_trim primers.bed --samfile input.bam --output trimmed.bam
```

#### Normalize coverage and generate reports
```bash
align_trim primers.bed --samfile input.bam --normalise 100 \
  --report alignment_report.tsv --amp-depth-report depth_report.tsv \
  --genome-coverage-report sample1 \
  --output normalized.bam
```

#### Process from stdin with verbose output
```bash
samtools view -h input.bam | align_trim primers.bed --verbose > trimmed.sam 2> verbose.out.txt
```

#### Strict full-length read filtering
```bash
align_trim primers.bed --samfile input.bam --require-full-length \
  --min-mapq 30 --output filtered.bam
```

#### Allow mismatched primer pairs with custom threshold
```bash
align_trim primers.bed --samfile input.bam --allow-incorrect-pairs \
  --primer-match-threshold 50 --output relaxed.bam
```

### Output Formats

The tool supports multiple output formats based on file extension:
- `.sam` - SAM format (text)
- `.bam` - BAM format (binary, compressed)
- No extension or `-` - SAM format to stdout

### Report Files

When using `--report`, a tab-separated file is generated with the following columns:
- `chrom`: Reference chromosome/contig
- `QueryName`: Read name
- `ReferenceStart`/`ReferenceEnd`: Alignment coordinates
- `PrimerPair`: Primer pair assignment
- `Primer1`/`Primer2`: Individual primer information
- `CorrectlyPaired`: Boolean indicating proper primer pairing
- Additional alignment metrics

The `--amp-depth-report` generates a summary of depth per amplicon with columns `chrom`, `amplicon` and `mean_fragment_depth`. Depth is counted in fragments (a read pair counts once where its mates overlap), the same unit used by `--normalise`. The mean is taken over the full amplicon including primer sites.

### Genome Coverage Report

When using `--genome-coverage-report PREFIX`, per-position genome coverage is written as tab-separated files with columns `chrom`, `pos` (1-based), and `depth`.

**Without `--normalise`**, a single file is produced:
- `PREFIX.pre-normalisation.coverage.tsv` — coverage of all reads passing filtering and trimming

**With `--normalise`**, two files are produced:
- `PREFIX.pre-normalisation.coverage.tsv` — coverage of all reads passing filtering and trimming (before normalisation subsampling)
- `PREFIX.post-normalisation.coverage.tsv` — coverage of reads retained after normalisation

In both cases, a coverage summary is printed to stderr showing the percentage of genome positions covered at >=1x, >=10x, >=20x, >=100x, and >=1000x.

Genome coverage counts reads, matching `samtools depth`, so where paired mates overlap a position is counted twice. This differs from the fragment units used by `--normalise` and `--amp-depth-report`.

### QC Reports

`--qc-prefix PREFIX` writes the following tab-separated reports. They do not change the trimmed output.

#### `PREFIX.read_fate.tsv`

One row per amplicon (including amplicons with no reads) plus `unassigned` rows, counting what happened to each template (a single-end read or a read pair). This separates amplicons that genuinely dropped out (`total` of 0) from amplicons whose reads arrived but were filtered.

Columns: `chrom`, `amplicon`, `pool`, `total`, then one column per fate:
- `kept`: written to the output
- `normalised_out`: rejected by `--normalise`
- `not_full_length`: rejected by `--require-full-length`
- `fully_masked`: no aligned bases left after primer trimming
- `softmask_error`: primer trimming failed (e.g. an indel inside the primer)
- `incorrectly_paired`: primers from different amplicons (see `primer_pairs.tsv`)
- `no_primer`: no primer found for one end. Attributed to the amplicon of the primer that was found, if any
- `low_mapq`, `no_mapping_data`, `unmapped`, `supplementary`, `secondary`, `missing_mate`: filtered before primer assignment, so always counted in the `unassigned` row for the chromosome (`*` if unknown)

Low mapping quality reads are not attributed to an amplicon, so amplicons inside repeated regions (e.g. poxvirus inverted terminal repeats, where reads map equally well to both copies) will show 0 reads here. These reads are also absent from the depth reports.

#### `PREFIX.primer_pairs.tsv`

Counts of each left/right primer combination, one row per combination seen, with columns `chrom_left`, `left_primer`, `left_amplicon`, `left_pool`, `chrom_right`, `right_primer`, `right_amplicon`, `right_pool`, `category` and `count`. Pairs are counted before `--allow-incorrect-pairs` is applied, so counts are the same with or without it. Adjacency is defined by amplicon number within a chromosome.

| category                    | primers                             | meaning                                                                       |
| --------------------------- | ----------------------------------- | ----------------------------------------------------------------------------- |
| `correct`                   | LEFT_n / RIGHT_n                    | expected product                                                              |
| `short_overlap`             | LEFT_n+1 / RIGHT_n                  | short product from the overlap between neighbouring amplicons                 |
| `long_same_pool`            | LEFT_n / RIGHT_n+k, same pool       | skip product, the intervening amplicon(s) failed or there is a large deletion |
| `long_cross_pool`           | LEFT_n / RIGHT_n+k, different pools | suggests pools were mixed or contaminated                                     |
| `cross_chrom`               | primers on different chromosomes    | inter-segment product or chimera                                              |
| `other`                     | anything else                       |                                                                               |
| `left_left` / `right_right` | two primers of the same direction   | mispriming by a same-direction primer                                         |

The first six categories come from primer coordinates. `left_left` and `right_right` come from primer sequence: when a read end is not within `--primer-match-threshold` of its assigned primer's 5' end, the end is searched for an exact match to the 3' 12 bases of a primer in the wrong orientation. These rows are counted in addition to the coordinate-based row for the same template. They are a lower bound, as sequencing errors in the primer prevent an exact match.

#### `PREFIX.split_alignments.tsv` and `PREFIX.split_breakpoints.tsv`

Reads with split alignments, taken from the `SA` tag on primary alignments (so it works when supplementary records have been removed). Each read passing the mapping filters is checked on its own, so paired mates are counted separately. Split reads are classed by where the breakpoints fall:
- `amplicon_junction`: every breakpoint is within `--primer-match-threshold` of an amplicon end, the signature of whole amplicons joined during library preparation (chimeras)
- `internal_breakpoint`: a breakpoint falls inside an amplicon, more likely a real structural difference from the reference (a large deletion, a rearrangement, or a change in repeat length)
- `low_confidence`: a part of the read has mapping quality below `--min-mapq`, as is common for repeated regions

`orientation` is `inverted` if parts align to different strands, `collinear` if parts are in reference order, otherwise `non_collinear`.

`split_alignments.tsv` counts reads by class, orientation and the amplicons of the first two parts (`chrom_a`, `amplicon_a`, `chrom_b`, `amplicon_b`, `reads`). `split_breakpoints.tsv` counts reads per breakpoint pair, binned to 10 bp (`chrom_a`, `bin_a`, `chrom_b`, `bin_b`, `class`, `orientation`, `reads`). A breakpoint pair recurring across many reads indicates a real structural variant, rather than chimeras which join amplicons at random. The chimera rate printed to stderr only counts `amplicon_junction` reads.

The BED file does not record whether a genome is circular. On a circular genome, reads crossing the origin are split at the end of the reference and will appear as a recurring `internal_breakpoint` between the end and start of the chromosome.

#### `PREFIX.overlap_positions.tsv` and `PREFIX.overlap_summary.tsv`

The region shared by neighbouring amplicons n and n+1, excluding primers, is sequenced by two independent PCR reactions. Comparing the alleles seen in each amplicon there acts as a technical replicate, so disagreement points to low template input, contamination of one pool, or an amplicon-specific artefact.

Only reads written to the output from correctly paired primers are counted, so with `--normalise` this is after normalisation. Each fragment is counted once per position.

`overlap_positions.tsv` has one row per overlap position with per-amplicon counts (`left_*`, `right_*` for `A`, `C`, `G`, `T`, `N`, `DEL`, `INS`; insertions are counted at the preceding position), the depth (`A`+`C`+`G`+`T`+`DEL`) and major allele for each side, `af_diff` (the largest difference in allele frequency between the two amplicons) and `evaluable` (both sides have depth >= 10).

`overlap_summary.tsv` has one row per overlap with `start`/`end` (1-based, inclusive), `evaluable_positions`, the median depth of each side, `max_af_diff` and its position, `discordant_positions` (evaluable with `af_diff` >= 0.2) and `major_allele_mismatches` (evaluable positions where the amplicons disagree on the major allele, the case most likely to change a consensus call).

At low depth (e.g. after normalising to a low target) random sequencing errors can produce some discordant positions. `major_allele_mismatches` is the more robust signal.

### Behaviour changes

- `--normalise` now counts paired data in fragments, so paired-end output at the same `-n` is roughly double the size of previous versions.
- The `--amp-depth-report` column is renamed from `mean_depth` to `mean_fragment_depth`. For paired data without `--normalise` it previously counted only one mate of each pair.
- Supplementary and secondary alignments are no longer paired in place of the real mate for paired-end data, and secondary alignments are no longer written for single-end data.
- Reads where primer trimming would produce an invalid CIGAR (e.g. an indel inside the primer) are now dropped. Previously they were written with their start position moved but the CIGAR unchanged.
- Reads aligned with `=`/`X` CIGAR operations are now trimmed. Previously they were dropped.
- A deletion left directly before the right primer softmask is now removed, as was already done for the left primer. Previously these reads ended in e.g. `1D22S`, reporting a spurious deletion at the last base before the primer.
