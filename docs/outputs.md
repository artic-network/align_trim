# Report outputs

Besides the trimmed SAM/BAM, `align_trim` can write the reports below. For the amplicon QC reports written by `--qc-prefix`, see [QC reports](qc-reports.md).

## Per-read report

When using `--report`, a tab-separated file is generated with the following columns:
- `chrom`: Reference chromosome/contig
- `QueryName`: Read name
- `ReferenceStart`/`ReferenceEnd`: Alignment coordinates
- `PrimerPair`: Primer pair assignment
- `Primer1`/`Primer2`: Individual primer information
- `CorrectlyPaired`: Boolean indicating proper primer pairing
- Additional alignment metrics

## Amplicon depth report

`--amp-depth-report` writes the depth of each amplicon with columns `chrom`, `amplicon` and `mean_fragment_depth`. Depth is counted in fragments (a read pair counts once where its mates overlap), the same unit used by `--normalise`. The mean is taken over the full amplicon including primer sites.

## Genome coverage report

When using `--genome-coverage-report PREFIX`, per-position genome coverage is written as tab-separated files with columns `chrom`, `pos` (1-based), and `depth`.

**Without `--normalise`**, a single file is produced:
- `PREFIX.pre-normalisation.coverage.tsv` — coverage of all reads passing filtering and trimming

**With `--normalise`**, two files are produced:
- `PREFIX.pre-normalisation.coverage.tsv` — coverage of all reads passing filtering and trimming (before normalisation subsampling)
- `PREFIX.post-normalisation.coverage.tsv` — coverage of reads retained after normalisation

In both cases, a coverage summary is printed to stderr under the `Coverage` heading, showing the percentage of genome positions covered at >=1x, >=10x, >=20x, >=100x, and >=1000x. In a terminal it is shown as a table with pre- and post-normalisation side by side; in a log it is one prefixed line per value (see [terminal output](usage.md#terminal-output)).

Genome coverage counts reads, matching `samtools depth`, so where paired mates overlap a position is counted twice. This differs from the fragment units used by `--normalise` and `--amp-depth-report`.
