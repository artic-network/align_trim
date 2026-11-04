# Changelog

All notable changes to align_trim are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0]

This release changes which reads are kept and how paired-end data is normalised, so trimmed output will differ from 1.2.0 for the same command. See **Changed** before upgrading a pipeline.

### Added

- `--qc-prefix PREFIX` writes five amplicon QC reports: a run summary, per-amplicon fragment fates, primer pair combinations, split (supplementary) alignments and concordance between overlapping amplicons. A one-line summary of each is printed to stderr. The trimmed output is unchanged by this option. See [QC reports](docs/qc-reports.md).
- A warning on stderr when the majority base differs between two overlapping amplicons, which can indicate contamination.
- stderr is grouped under `Run`, `Coverage`, `QC summary` and `Warnings` headings. A completion line reports the records and read pairs processed, the time taken and the rate. See [terminal output](docs/usage.md#terminal-output).
- In a terminal: coloured headings, a live progress bar (percentage and time remaining for a BAM file, count and rate for stdin) and the coverage summary as a pre/post-normalisation table. Logs get plain `== Name ==` headings and a `progress:` line every 30 seconds. `NO_COLOR` and `FORCE_COLOR` are respected.
- `--quiet`/`-q` prints only the coverage summary, QC summary and warnings to stderr.
- Fragments lying inside two neighbouring amplicons (left primer of amplicon n+1, right primer of amplicon n) are removed as `ambiguous_amplicon` and never counted against either amplicon. `overlaps.tsv` gives the count for each overlap.
- `--rejected-output FILE` writes the reads that were removed, untrimmed, to a separate SAM/BAM file, each tagged `ZR:Z:<reason>` with the same reasons as the QC summary. See [rejected reads](docs/outputs.md#rejected-reads).

### Changed

- With `--allow-incorrect-pairs`, fragments lying inside two neighbouring amplicons are now removed. Previously they were kept, trimmed and normalised as amplicon n+1.
- **Breaking:** the `--primer-match-threshold` default is now 5 (previously 35). Reads with an end more than 5 bases outside their primer's 5' end, i.e. containing genomic sequence beyond the amplicon, are no longer assigned to that primer. On the paired-end test data this removes about 0.4% of read pairs; single-end ONT test output is unchanged.
- **Breaking:** `--normalise` counts paired-end data in fragments (a read pair counts once where its mates overlap), so paired-end output at the same `-n` is roughly double the size of previous versions.
- **Breaking:** the `--amp-depth-report` column `mean_depth` is renamed `mean_fragment_depth`.
- Reads skipped because primer soft masking failed are no longer printed one per line by default. A single warning gives the count, and `--verbose` prints each read. `--verbose` also disables the live progress bar.
- The README covers installation and a quick start; the full option reference and report descriptions are in [`docs/`](docs/).

### Fixed

- A read ending exactly `--primer-match-threshold` bases past its right primer's 5' end was given the next amplicon's right primer and removed as mismatched, and a read ending at the last base of a reference matched no right primer. Right ends are now matched within the threshold, as left ends already were.
- `--amp-depth-report` for paired-end data without `--normalise` counted only one mate of each pair.
- Supplementary and secondary alignments could be paired in place of the real mate for paired-end data, and secondary alignments were written for single-end data.
- Reads where primer trimming would produce an invalid CIGAR (e.g. an indel inside the primer) were written with their start position moved but the CIGAR unchanged. They are now dropped.
- Reads aligned with `=`/`X` CIGAR operations were dropped instead of trimmed.
- A deletion left directly before the right primer softmask was kept, so reads ended in e.g. `1D22S` and reported a spurious deletion at the last base before the primer. It is now removed, as was already done for the left primer.
- For paired-end data the right primer was looked up at the end of whichever mate started later, rather than the rightmost mate end, so a pair with one mate contained in the other could be given the wrong right primer or removed by `--require-full-length`.

## [1.2.0] - 2026-06-23

### Added

- `--genome-coverage-report PREFIX` writes per-position genome coverage before and after normalisation, with a coverage summary on stderr ([#23](https://github.com/artic-network/align_trim/pull/23)).

## [1.1.0] - 2026-04-09

### Fixed

- Normalisation logic ([#22](https://github.com/artic-network/align_trim/pull/22)).

## [1.0.2] - 2025-08-28

### Added

- Check for paired reads ([#19](https://github.com/artic-network/align_trim/pull/19)).

## [1.0.1] - 2025-08-22

### Changed

- Unpinned numpy and pysam versions for wider compatibility.

## [1.0.0] - 2025-08-22

Initial release of align_trim as a standalone package. It was originally developed for primer trimming in [artic-network/fieldbioinformatics](https://github.com/artic-network/fieldbioinformatics).

[2.0.0]: https://github.com/artic-network/align_trim/compare/1.2.0...HEAD
[1.2.0]: https://github.com/artic-network/align_trim/compare/1.1.0...1.2.0
[1.1.0]: https://github.com/artic-network/align_trim/compare/1.0.2...1.1.0
[1.0.2]: https://github.com/artic-network/align_trim/compare/1.0.1...1.0.2
[1.0.1]: https://github.com/artic-network/align_trim/compare/1.0.0...1.0.1
[1.0.0]: https://github.com/artic-network/align_trim/releases/tag/1.0.0
