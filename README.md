# align_trim

Stand alone version of ARTIC's fieldbioinformatics align_trim.py. It assigns aligned reads from a tiling amplicon scheme to their amplicons, soft-clips the primers and can normalise depth per amplicon.

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

## Quick start

```bash
align_trim primer.bed --samfile input.bam --output trimmed.bam
```

- `primer.bed` is the primer scheme in [primer.bed v3](https://doi.org/10.5281/zenodo.16366659) format.
- `input.bam` is a sorted BAM of reads aligned to the scheme's reference. Leave out `--samfile` to read SAM/BAM from stdin.
- The output format follows the extension (`.bam` or `.sam`). Leave out `--output` to write SAM to stdout.

Normalise to 200x per amplicon and write the amplicon QC reports:

```bash
align_trim primer.bed --samfile input.bam --output trimmed.bam \
  --normalise 200 --qc-prefix sample1
```

## Common options

| option | meaning |
| --- | --- |
| `--normalise N`, `-n` | keep reads until each amplicon reaches a depth of N fragments (a read pair counts once). 0 turns it off (default) |
| `--min-mapq`, `-m` | minimum mapping quality (default: 20) |
| `--primer-match-threshold`, `-p` | how far, in bases, a read end can be from a primer and still match it (default: 5) |
| `--allow-incorrect-pairs` | keep reads whose two ends match primers from different amplicons (reads inside two neighbouring amplicons are still removed) |
| `--require-full-length` | only keep reads spanning a whole amplicon (not for rapid barcoding or other fragmented libraries) |
| `--no-trim-primers` | assign reads but leave primers untrimmed |
| `--qc-prefix PREFIX` | write the amplicon QC reports as `PREFIX.<name>.tsv` |

Run `align_trim --help` for every option. Progress, coverage and QC summaries and warnings are printed to stderr; see [terminal output](https://github.com/artic-network/align_trim/blob/main/docs/usage.md#terminal-output).

## Documentation

- [Usage](https://github.com/artic-network/align_trim/blob/main/docs/usage.md): every option, examples and output formats
- [Report outputs](https://github.com/artic-network/align_trim/blob/main/docs/outputs.md): `--report`, `--amp-depth-report` and `--genome-coverage-report`
- [QC reports](https://github.com/artic-network/align_trim/blob/main/docs/qc-reports.md): the `--qc-prefix` reports and how to read them
- [Changelog](https://github.com/artic-network/align_trim/blob/main/CHANGELOG.md): changes in each release, including behaviour changes in 2.0.0
