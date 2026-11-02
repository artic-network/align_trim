# Behaviour changes

Changes to the trimmed output and existing reports compared with align_trim 1.2.0.

- The `--primer-match-threshold` default is now 5 (previously 35). Reads with an end more than 5 bases outside its primer's 5' end, i.e. containing genomic sequence beyond the amplicon, are no longer assigned to that primer. On the paired-end test data this removes about 0.4% of read pairs; single-end ONT test output is unchanged.
- `--normalise` now counts paired data in fragments, so paired-end output at the same `-n` is roughly double the size of previous versions.
- The `--amp-depth-report` column is renamed from `mean_depth` to `mean_fragment_depth`. For paired data without `--normalise` it previously counted only one mate of each pair.
- Supplementary and secondary alignments are no longer paired in place of the real mate for paired-end data, and secondary alignments are no longer written for single-end data.
- Reads where primer trimming would produce an invalid CIGAR (e.g. an indel inside the primer) are now dropped. Previously they were written with their start position moved but the CIGAR unchanged.
- Reads aligned with `=`/`X` CIGAR operations are now trimmed. Previously they were dropped.
- stderr is now grouped under `Run`, `Coverage`, `QC summary` and `Warnings` headings (`== Name ==` lines in a log) with a completion line after the reads are processed, and a progress bar in a terminal. Existing coverage and `QC:` lines are unchanged. See [terminal output](usage.md#terminal-output).
- Reads skipped because primer soft masking failed are no longer printed one per line by default. A single warning gives the count; `--verbose` prints each read.
- A deletion left directly before the right primer softmask is now removed, as was already done for the left primer. Previously these reads ended in e.g. `1D22S`, reporting a spurious deletion at the last base before the primer.
