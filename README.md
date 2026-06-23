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

- `--normalise`, `-n` : Normalise to target depth N per amplicon using a greedy per-read algorithm. Each read is kept only if it brings the amplicon depth closer to the target. Use 0 for no normalisation (default: 0)
- `--min-mapq`, `-m` : Minimum mapping quality to keep an aligned read (default: 20)
- `--primer-match-threshold`, `-p` : Add this many bases of padding to the 5' end of primer coordinates to allow fuzzy matching for reads with barcodes/adapters (default: 35)

#### Primer and Read Handling

- `--no-trim-primers` : Do not trim primers from reads (by default, primers are trimmed)
- `--allow-incorrect-pairs` : Allow reads to be assigned to amplicons even if primers are not correctly paired
- `--require-full-length` : Require all reads to start and stop in primer sites (do not use with rapid barcoding)

#### Output and Reporting

- `--report`, `-r` : Output detailed report TSV to specified filepath
- `--amp-depth-report`, `-a` : Output mean depth for each amplicon as TSV to specified filepath
- `--genome-coverage-report`, `-g` : Output per-position genome coverage TSV(s) using the given prefix. Summary statistics (% genome covered at >=1x, >=10x, >=20x, >=100x) are printed to stderr. See [Genome Coverage Report](#genome-coverage-report) for details
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

The `--amp-depth-report` generates a summary of coverage depth per amplicon.

### Genome Coverage Report

When using `--genome-coverage-report PREFIX`, per-position genome coverage is written as tab-separated files with columns `chrom`, `pos` (1-based), and `depth`.

**Without `--normalise`**, a single file is produced:
- `PREFIX.pre-normalisation.coverage.tsv` — coverage of all reads passing filtering and trimming

**With `--normalise`**, two files are produced:
- `PREFIX.pre-normalisation.coverage.tsv` — coverage of all reads passing filtering and trimming (before normalisation subsampling)
- `PREFIX.post-normalisation.coverage.tsv` — coverage of reads retained after normalisation

In both cases, a coverage summary is printed to stderr showing the percentage of genome positions covered at >=1x, >=10x, >=20x, >=100x, and >=1000x.
