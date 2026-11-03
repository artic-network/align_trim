# Usage

```bash
align_trim [OPTIONS] BEDFILE
```

The tool reads alignment data from either a SAM/BAM file or stdin and outputs trimmed alignments to stdout in SAM format by default.

## Required arguments

- `BEDFILE`: BED file containing the amplicon primer scheme in [v3](https://doi.org/10.5281/zenodo.16366659) format.

## Optional arguments

### Input/output

- `--samfile`, `-s` : Sorted SAM/BAM file containing the aligned reads, if this is not provided (or '-') then 'align_trim' will read from stdin.
- `--output`, `-o` : Output file path. Format determined by extension (.sam/.bam). If not provided or '-', writes SAM to stdout

### Processing

- `--normalise`, `-n` : Normalise to target depth N per amplicon using a greedy per-read algorithm. Each read (or read pair) is kept only if it brings the amplicon depth closer to the target. Depth is counted in fragments: a read pair counts once where its mates overlap, so N means N fragments for both single-end and paired-end data. Use 0 for no normalisation (default: 0)
- `--min-mapq`, `-m` : Minimum mapping quality to keep an aligned read (default: 20)
- `--primer-match-threshold`, `-p` : Maximum distance in bases between a read end and a primer 5' end for the read to be assigned to that primer. Aligners soft-clip barcodes and adapters, so genuine amplicon reads end within a few bases of the primer. The same value is used by the `--qc-prefix` reports to decide whether a read end is at a primer (default: 5)

### Primer and read handling

- `--no-trim-primers` : Do not trim primers from reads (by default, primers are trimmed)
- `--allow-incorrect-pairs` : Allow reads to be assigned to amplicons even if primers are not correctly paired. Reads lying inside two neighbouring amplicons (left primer of amplicon n+1, right primer of amplicon n) are still removed, as they cannot be assigned to either
- `--require-full-length` : Require all reads to start and stop in primer sites (do not use with rapid barcoding)

### Output and reporting

- `--report`, `-r` : Output detailed report TSV to specified filepath. See [outputs](outputs.md#per-read-report)
- `--amp-depth-report`, `-a` : Output mean fragment depth for each amplicon as TSV to specified filepath. See [outputs](outputs.md#amplicon-depth-report)
- `--genome-coverage-report`, `-g` : Output per-position genome coverage TSV(s) using the given prefix. Summary statistics (% genome covered at >=1x, >=10x, >=20x, >=100x, >=1000x) are printed to stderr. See [outputs](outputs.md#genome-coverage-report)
- `--qc-prefix` : Output amplicon QC reports as `PREFIX.<name>.tsv`, with a one line summary per report printed to stderr. See [QC reports](qc-reports.md)
- `--rejected-output FILE` : Also write the reads which were removed to FILE (.bam or .sam), as they were read in (untrimmed) and tagged with the reason. See [outputs](outputs.md#rejected-reads)
- `--no-read-groups` : Do not divide reads into pool-based read groups in SAM/BAM output

### General

- `--verbose`, `-v` : Debug mode: print per-read decisions to stderr, including each read skipped because primer soft masking failed. Disables the live progress bar
- `--quiet`, `-q` : Only print the coverage summary, QC summary and warnings to stderr (no run details or progress)
- `--version` : Show version information
- `--help` : Show help message

## Terminal output

The trimmed SAM/BAM is the only thing written to stdout. Everything else goes to stderr, grouped under headings in this order (a section only appears when it has something to show):

1. `Run`: version, input, scheme and output, progress while reads are processed, then one completion line (records, read pairs, time and rate)
2. `Coverage`: with `--genome-coverage-report`
3. `QC summary`: with `--qc-prefix`
4. `Warnings`: e.g. a majority base difference between overlapping amplicons, or how many fragments were skipped because primer soft masking failed

In a terminal, headings are coloured rules and a live progress bar shows the percentage done and time remaining for a BAM file, or a count and rate when reading from stdin. The live bar is not shown with `--verbose`, or when SAM is written to the same terminal.

When stderr is redirected to a file (e.g. a pipeline log) the output is plain text with no colour or progress bar: headings are written as `== Name ==` lines, and a `progress:` line is written every 30 seconds. Lines keep their `QC:`, `WARNING:` and `Pre-normalisation:`/`Post-normalisation:` prefixes, so a section can be pulled out with e.g. `grep '^QC:' sample.log`.

Set `NO_COLOR=1` to turn off colour in a terminal, or `FORCE_COLOR=1` to keep it when stderr is redirected.

```
== Run ==
align_trim 2.0.0
  input   sample.sorted.bam (paired-end)
  scheme  primer.bed (48 amplicons, 1 reference)
  output  trimmed.bam
progress: 1,163,264 records (73%), 38.5k records/s
processed 1,601,770 records (760,958 read pairs) in 40.2s, 39.8k records/s
== Coverage ==
Pre-normalisation: Genome coverage summary
...
== QC summary ==
QC: 800518 fragments, 41236 kept (5.15%); ...
...
== Warnings ==
WARNING: majority base differs between overlapping amplicons ..., which can indicate contamination; see sample1.overlaps.tsv
```

## Examples

### Basic trimming with primer removal
```bash
align_trim primer.bed --samfile input.bam --output trimmed.bam
```

### Normalise coverage and generate reports
```bash
align_trim primer.bed --samfile input.bam --normalise 100 \
  --report alignment_report.tsv --amp-depth-report depth_report.tsv \
  --genome-coverage-report sample1 \
  --output normalized.bam
```

### Amplicon QC reports
```bash
align_trim primer.bed --samfile input.bam --qc-prefix sample1 --output trimmed.bam
```

### Keep the rejected reads
```bash
align_trim primer.bed --samfile input.bam --output trimmed.bam \
  --rejected-output rejected.bam
```

### Process from stdin with verbose output
```bash
samtools view -h input.bam | align_trim primer.bed --verbose > trimmed.sam 2> verbose.out.txt
```

### Strict full-length read filtering
```bash
align_trim primer.bed --samfile input.bam --require-full-length \
  --min-mapq 30 --output filtered.bam
```

### Allow mismatched primer pairs with custom threshold
```bash
align_trim primer.bed --samfile input.bam --allow-incorrect-pairs \
  --primer-match-threshold 50 --output relaxed.bam
```

## Output formats

The output format is set by the file extension:
- `.sam` - SAM format (text)
- `.bam` - BAM format (binary, compressed)
- No extension or `-` - SAM format to stdout
