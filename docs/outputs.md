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

## Rejected reads

`--rejected-output FILE` writes the reads that are removed to a separate SAM/BAM file (format from the `.bam`/`.sam` extension; it cannot be stdout). The trimmed output is unchanged.

- Every primary record is in exactly one of the two files: read pairs are written together, so a read name is never split between them.
- Records are written as they were read in, untrimmed, with their input tags and the input header (plus an align_trim `@PG` line).
- Each record has a `ZR:Z:<reason>` tag. The reasons are the same as the `removed_<reason>` columns of the [QC summary](qc-reports.md#prefixsummarytsv) (`unmapped`, `not_proper_pair`, `missing_mate`, `low_mapq`, `no_alignment`, `no_primer_match`, `ambiguous_amplicon`, `mismatched_primers`, `trim_failed`, `primer_only`, `not_full_length`, `by_normalise`), and the counts match.
- Supplementary and secondary alignments are not written to either file.
- With `--normalise`, most rejected reads are usually `by_normalise`, which makes the file large and the run slower (about a quarter longer for a 1.6M-record paired-end run at `-n 200`).
- Single-end output stays in input order. Paired-end output is written as each pair completes, with unpaired mates at the end, so the header is `SO:unsorted`: run `samtools sort` before indexing.

To pull out the reads removed for one reason:

```bash
samtools view -b -d ZR:ambiguous_amplicon rejected.bam > ambiguous_amplicon.bam
```
