# QC reports

`--qc-prefix PREFIX` writes five tab-separated reports. They do not change the trimmed output.

Conventions used in every report:
- **Counts are in fragments.** A fragment is one sequenced molecule: a single-end read, or a read pair (counted once). Supplementary and secondary alignments are extra alignments of a fragment already counted, so are not counted again.
- Coordinates are 1-based and inclusive.
- `NA` means a value could not be calculated (e.g. too little depth). It never stands in for 0.
- On every row with a `kept` column, `kept` plus all `removed_*` columns equals `fragments`.

A short summary of each report is printed to stderr under the `QC summary` heading, prefixed `QC:` (see [terminal output](usage.md#terminal-output)).

## `PREFIX.summary.tsv`

One row for the whole input, so files from many samples can be stacked. Columns are `fragments`, `kept` and one `removed_<reason>` column per reason, in the order the checks are applied:

| column                       | fragments removed because                                                                                                                    |
| ---------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| `removed_unmapped`           | the read, or either read of a pair, is unmapped                                                                                              |
| `removed_not_proper_pair`    | the aligner did not flag the reads as a proper pair (e.g. different chromosomes, wrong orientation or too far apart)                         |
| `removed_missing_mate`       | only one read of the pair was in the input                                                                                                   |
| `removed_low_mapq`           | mapping quality below `--min-mapq`                                                                                                           |
| `removed_no_alignment`       | flagged as mapped but with no aligned bases                                                                                                  |
| `removed_no_primer_match`    | a fragment end is not within `--primer-m  ch-threshold` of any primer                                                                        |
| `removed_partial_in_overlap` | a partial fragment lying inside the region shared by two neighbouring amplicons, so it cannot be assigned to either (see `primer_pairs.tsv`) |
| `removed_mismatched_primers` | the two ends match primers from different amplicons (see `primer_pairs.tsv`)                                                                 |
| `removed_trim_failed`        | primer trimming could not be applied (e.g. an indel inside the primer)                                                                       |
| `removed_primer_only`        | nothing is left after trimming primers                                                                                                       |
| `removed_not_full_length`    | rejected by `--require-full-length`                                                                                                          |
| `removed_by_normalise`       | passed every check but was not needed to reach the `--normalise` depth                                                                       |

The two mismatched primer reasons only remove fragments without `--allow-incorrect-pairs`.

## `PREFIX.amplicons.tsv`

One row per amplicon in the scheme, including amplicons with no fragments, so a genuine dropout (`fragments` of 0) can be told apart from an amplicon whose fragments were sequenced but removed.

Columns: `chrom`, `amplicon`, `pool`, `start`, `end` (the amplicon including primers), `fragments`, `kept`, and the removal reasons that apply once a fragment is assigned to an amplicon (`removed_no_primer_match` to `removed_by_normalise` above). A fragment with only one end matching a primer is counted against that primer's amplicon. The earlier reasons are only counted in `summary.tsv`.

Fragments with low mapping quality are not assigned to an amplicon, so amplicons inside repeated regions (e.g. poxvirus inverted terminal repeats, where reads map equally well to both copies) will show 0 fragments here. These fragments are also absent from the depth reports.

## `PREFIX.primer_pairs.tsv`

Counts of each combination of primers matched at the fragment's left (lower coordinate) and right ends, one row per combination seen. Columns: `left_chrom`, `left_primer`, `left_amplicon`, `right_chrom`, `right_primer`, `right_amplicon`, `type` and `fragments`. Counts are made before `--allow-incorrect-pairs` is applied, so they are the same with or without it.

| type                                     | primers                             | meaning                                                                                                                                                                                                                                     |
| ---------------------------------------- | ----------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `expected`                               | LEFT_n / RIGHT_n                    | the intended product                                                                                                                                                                                                                        |
| `overlap_product`                        | LEFT_n+1 / RIGHT_n                  | a short product of only the region shared by neighbouring amplicons, with both ends at those primers. The two primers are in different pools, so this suggests the pools were mixed. Can be overcounted in fragmented libraries (see below) |
| `partial_in_overlap`                     | LEFT_n+1 / RIGHT_n                  | a fragment inside the shared region whose ends are not at those primers: a truncated read or a piece of a fragmented (e.g. tagmented) library that only lies in the overlap by chance                                                       |
| `long_same_pool`                         | LEFT_n / RIGHT_n+k, same pool       | a long product spanning a primer site that should have stopped it: the primer failed (e.g. a mutation in its binding site) or there is a large deletion                                                                                     |
| `long_cross_pool`                        | LEFT_n / RIGHT_n+k, different pools | suggests the pools were mixed or contaminated                                                                                                                                                                                               |
| `cross_chrom`                            | primers on different chromosomes    | an inter-segment product or chimera                                                                                                                                                                                                         |
| `other`                                  | anything else                       | e.g. a rearrangement                                                                                                                                                                                                                        |
| `two_left_primers` / `two_right_primers` | two primers of the same direction   | mispriming by a primer of the wrong orientation                                                                                                                                                                                             |

An end is at a primer when it is within `--primer-match-threshold` of the primer's 5' end. Long products are not checked this way, as their span cannot come from any single amplicon wherever their ends fall.

**Fragmented libraries.** When a library prep cuts amplicons at random positions (e.g. Illumina Nextera, other tagmentation preps, or ONT rapid barcoding), some fragments lie entirely inside an overlap region. Most of these are counted as `partial_in_overlap`, but any whose ends happen to fall within `--primer-match-threshold` of both overlap primers are counted as `overlap_product`, even though no cross-pool product was made. How often this happens depends on the overlap length relative to the mean fragment size: the shorter the fragments are compared with the overlap, the more fragments fit inside it. ONT rapid barcoding fragments are usually long enough that this is negligible for typical overlaps (under 100 bp), while short-fragment Illumina tagmentation preps are more affected. Where many fragments are short relative to the overlap, treat `overlap_product` counts as an upper bound rather than as evidence of mixed pools on their own. A high `partial_in_overlap` count is expected in that case and is not a fault.

`two_left_primers` and `two_right_primers` are found by primer sequence: when a fragment end is not at its matched primer, the end is searched for an exact match to the 3' 12 bases of a primer in the wrong orientation. These rows are counted in addition to the coordinate-based row for the same fragment, and are a lower bound, as sequencing errors in the primer prevent an exact match.

## `PREFIX.split_reads.tsv`

Joins in split alignments, taken from the `SA` tag of primary alignments (so this works when supplementary records have been removed). One row per join, grouped into 10 bp windows, with the lower coordinate side first:

| column                                          | meaning                                                                                                                                                                                                                                                                                                                                                                                                                          |
| ----------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `chrom`, `pos`, `amplicon`                      | one side of the join: the start of the 10 bp window containing the aligned base next to the join on that side, and the amplicon of that part of the read                                                                                                                                                                                                                                                                         |
| `joined_chrom`, `joined_pos`, `joined_amplicon` | the other side of the join                                                                                                                                                                                                                                                                                                                                                                                                       |
| `type`                                          | `amplicon_chimera`: both sides are at amplicon ends (within `--primer-match-threshold`), the signature of whole amplicons joined during library preparation. `within_amplicon`: a side falls inside an amplicon, more likely a real difference from the reference (a large deletion, a rearrangement or a change in repeat length). `low_mapq`: a part has mapping quality below `--min-mapq`, as is common for repeated regions |
| `orientation`                                   | `collinear`: both parts on the same strand and in reference order (deletion-like). `non_collinear`: same strand but out of order or on different chromosomes. `inverted`: parts on opposite strands                                                                                                                                                                                                                              |
| `fragments`                                     | fragments containing this join (a join seen in both reads of a pair is counted once)                                                                                                                                                                                                                                                                                                                                             |

A join recurring across many fragments indicates a real structural variant, whereas chimeras join amplicons at random. The chimera rate printed to stderr counts fragments with an `amplicon_chimera` join.

The BED file does not record whether a genome is circular. On a circular genome, reads crossing the origin are split at the end of the reference and will appear as a recurring `within_amplicon` join between the end and start of the chromosome.

## `PREFIX.overlaps.tsv`

The region shared by neighbouring amplicons n and n+1, excluding primers, is sequenced by two independent PCR reactions. Comparing the bases seen in each amplicon there acts as a technical replicate, so disagreement points to contamination of one amplicon or pool, low template input, or an amplicon-specific artefact. Only fragments written to the output with matching primers are counted, so with `--normalise` this is after normalisation.

| column                                                | meaning                                                                                                                                                                                                                                                                             |
| ----------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `chrom`, `amplicon`, `next_amplicon`                  | the neighbouring amplicons n and n+1                                                                                                                                                                                                                                                |
| `start`, `end`, `length`                              | the shared region, excluding primers                                                                                                                                                                                                                                                |
| `compared_positions`                                  | positions with a depth of at least 10 in both amplicons                                                                                                                                                                                                                             |
| `amplicon_median_depth`, `next_amplicon_median_depth` | median depth over the region from each amplicon's fragments                                                                                                                                                                                                                         |
| `max_freq_difference`, `max_freq_difference_pos`      | the largest difference between the amplicons in the frequency of any base (A, C, G, T, deletion or insertion) at a compared position, from 0 to 1, and where it is                                                                                                                  |
| `differing_positions`                                 | compared positions where some base's frequency differs by at least 0.2                                                                                                                                                                                                              |
| `majority_base_differences`                           | compared positions where the most common base differs between the amplicons. This compares base frequencies in the reads, not variant or consensus calls, but it is the case most likely to change a consensus made from this output                                                |
| `majority_base_difference_details`                    | each majority base difference as `pos:frequencies/next_frequencies`, separated by `;`. Frequencies are the relative frequency of each base (A, C, G, T or DEL) in that amplicon's reads, most common first, with bases below 0.01 left out, e.g. `3623:C=0.95,T=0.05/T=0.88,C=0.12` |

The last five columns are `NA` when no positions could be compared. At low depth (e.g. after normalising to a low target) random sequencing errors can produce some differing positions, so `majority_base_differences` is the more robust signal.

When any overlap has a majority base difference, a `WARNING: majority base differs between overlapping amplicons ...` line naming those overlaps is printed to stderr under the `Warnings` heading. It can indicate contamination, but also low template input or an amplicon-specific artefact, and does not say which of the two amplicons is affected. The frequencies in `majority_base_difference_details` help tell these apart, e.g. a near 50/50 split in one amplicon against a clean call in the other.
