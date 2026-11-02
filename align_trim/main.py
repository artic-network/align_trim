import argparse
import csv
import itertools
import sys
from collections import defaultdict
from copy import copy
from importlib.metadata import version
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pysam
from primalbedtools.amplicons import Amplicon, create_amplicons
from primalbedtools.bedfiles import BedLine, merge_primers
from primalbedtools.scheme import Scheme

from align_trim.console import BOLD, bam_fraction, console
from align_trim.qc import DEFAULT_PRIMER_MATCH_THRESHOLD, QC

# consumesReference lookup for if a CIGAR operation consumes the reference sequence
consumesReference = [True, False, True, True, False, False, False, True, True]

# consumesQuery lookup for if a CIGAR operation consumes the query sequence
consumesQuery = [True, True, False, False, True, False, False, True, True]

# CIGAR operations which align a query base to a reference base (M, =, X)
matchOps = (0, 7, 8)


def find_primer_with_lookup(lookup, pos, direction, chrom) -> Optional[BedLine]:
    pos_amps = lookup[chrom][:, pos]  # Search both pools for amplicons at this position
    closest_dist = float("inf")
    closest_p = None
    if direction == "+":
        # Loops over pool O(N)
        for amp in pos_amps:
            if amp is None:
                continue
            dist = abs(amp.coverage_start - pos)
            if dist < closest_dist:
                closest_p = amp.left[0]
                closest_dist = dist
    elif direction == "-":
        for amp in pos_amps:
            if amp is None:
                continue
            dist = abs(amp.coverage_end - pos)
            if dist < closest_dist:
                closest_p = amp.right[0]
                closest_dist = dist
    else:
        pass
    return closest_p


def find_primer(primers: list[BedLine], pos, direction, chrom, threshold=35):
    """Given a reference position and a direction of travel, walk out and find the nearest primer site.

    Parameters
    ----------
    bed : list
        A list of dictionaries, where each dictionary contains a row of bedfile data
    pos : int
        The position in the reference sequence to start from
    direction : string
        The direction to search along the reference sequence

    Returns
    -------
    tuple[int, int, dict] | bool
        A tuple containing the distance to the primer, the relative position of the primer, and the primer site, or False if no primer found
    """
    from operator import itemgetter

    if direction == "+":
        primer_distances = [
            (abs(bl.start - pos), bl.start - pos, bl)
            for bl in primers
            if (pos >= (bl.start - threshold)) and chrom == bl.chrom
        ]

    else:
        primer_distances = [
            (abs(bl.end - pos), bl.end - pos, bl)
            for bl in primers
            if (pos <= (bl.end + threshold)) and chrom == bl.chrom
        ]

    if not primer_distances:
        return False

    closest = min(
        primer_distances,
        key=itemgetter(0),
    )

    return closest


def trim(segment, primer_pos, end, verbose=False):
    """Soft mask an alignment to fit within primer start/end sites.

    Parameters
    ----------
    segment : pysam.AlignedSegment
        The aligned segment to mask
    primer_pos : int
        The position in the reference to soft mask up to (equates to the start/end position of the primer in the reference)
    end : bool
        If True, the segment is being masked from the end (i.e. for the reverse primer)
    verbose : bool
        If True, will print soft masking info during trimming
    """
    if verbose:
        console.debug(
            f"{segment.query_name}: Trimming {'end' if end else 'start'} of read to primer position {primer_pos}"
        )
    # get a copy of the cigar tuples to work with
    cigar = copy(segment.cigartuples)

    # get the segment position in the reference (depends on if start or end of the segment is being processed)
    if not end:
        pos = segment.pos
    else:
        pos = segment.reference_end

    # process the CIGAR to determine how much softmasking is required
    eaten = 0
    while 1:
        # chomp CIGAR operations from the start/end of the CIGAR
        try:
            if end:
                flag, length = cigar.pop()
            else:
                flag, length = cigar.pop(0)
            if verbose:
                console.debug(f"{segment.query_name}: Chomped a {flag}, {length}")
        except IndexError:
            if verbose:
                console.debug(
                    f"{segment.query_name}: Ran out of cigar during soft masking - completely masked read will be ignored"
                )
            break

        # if the CIGAR operation consumes the reference sequence, increment/decrement the position by the CIGAR operation length
        if consumesReference[flag]:
            if not end:
                pos += length
            else:
                pos -= length

        # if the CIGAR operation consumes the query sequence, increment the number of CIGAR operations eaten by the CIGAR operation length
        if consumesQuery[flag]:
            eaten += length

        # stop processing the CIGAR if we've gone far enough to mask the primer
        if not end and pos >= primer_pos and flag in matchOps:
            break
        if end and pos <= primer_pos and flag in matchOps:
            break

    # calculate how many extra matches are needed in the CIGAR
    extra = abs(pos - primer_pos)
    if verbose:
        console.debug(f"{segment.query_name}: extra {extra}")
    if extra:
        if verbose:
            console.debug(f"{segment.query_name}: Inserted a 0, {extra}")
        if end:
            cigar.append((0, extra))
        else:
            cigar.insert(0, (0, extra))
        eaten -= extra

    # softmask the left primer
    if not end:
        # update the position of the leftmost mapping base
        segment.pos = pos - extra
        if verbose:
            console.debug(f"{segment.query_name}: New pos - {segment.pos}")

        # if proposed softmask leads straight into a deletion, shuffle leftmost mapping base along and ignore the deletion
        if cigar[0][0] == 2:
            if verbose:
                console.debug(
                    f"{segment.query_name}: softmask created a leading deletion in the CIGAR, shuffling the alignment"
                )
            while 1:
                if cigar[0][0] != 2:
                    break
                _, length = cigar.pop(0)
                segment.pos += length

        # now add the leading softmask
        cigar.insert(0, (4, eaten))

    # softmask the right primer
    else:
        # if proposed softmask leads straight out of a deletion, drop the trailing deletion (the reference end moves left with it)
        while cigar and cigar[-1][0] == 2:
            if verbose:
                console.debug(
                    f"{segment.query_name}: softmask created a trailing deletion in the CIGAR, removing it"
                )
            cigar.pop()

        cigar.append((4, eaten))

    # check the new CIGAR and replace the old one, segment.pos may already have moved so the segment must be dropped
    if cigar[0][1] <= 0 or cigar[-1][1] <= 0:
        raise ValueError(
            "invalid cigar operation created - possibly due to INDEL in primer"
        )

    segment.cigartuples = cigar
    return


def fragment_spans(segments) -> list[tuple[int, int]]:
    """Reference intervals covered by a template (one read or a read pair), with overlapping mates merged so each position is counted once.

    Parameters
    ----------
    segments : tuple[pysam.AlignedSegment, ...]
        The segment(s) making up the template

    Returns
    -------
    list[tuple[int, int]]
        Disjoint (start, end) intervals in reference coordinates
    """
    spans = sorted((seg.reference_start, seg.reference_end) for seg in segments)
    merged = [spans[0]]
    for start, end in spans[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def add_fragment_depth(depths, segments, offset):
    """Add +1 per position covered by a template to an amplicon depth array starting at reference position `offset`."""
    for start, end in fragment_spans(segments):
        depths[max(0, start - offset) : max(0, end - offset)] += 1


def _skip(qc, fate, chrom=None, amplicon=None):
    """Record why a fragment (read or read pair) was removed (if QC reports were requested) and return False."""
    if qc is not None:
        qc.fate(fate, chrom, amplicon)
    return False


def handle_segments(
    segment: Union[
        pysam.AlignedSegment, tuple[pysam.AlignedSegment, pysam.AlignedSegment]
    ],
    lookup: dict,
    args: argparse.Namespace,
    min_mapq: int,
    outfile_writer: pysam.AlignmentFile,
    amp_depths: dict,
    report_writer: csv.DictWriter = False,  # type: ignore
    genome_coverage: Optional[dict] = None,
    qc=None,
):
    """Handle the alignment segment(s) including filtering, soft masking, and reporting.

    Args:
        segment (pysam.AlignedSegment | tuple): The alignment segment to process, can be a single segment or a tuple of paired segments
        bed (dict): The primer scheme
        reportfh (typing.IO): The report file handle
        args (argparse.Namespace): The command line arguments
        qc (align_trim.qc.QC | None): QC report accumulator, None if --qc-prefix was not given

    Returns:
        tuple [int, pysam.AlignedSegment | bool, bool] | bool: A tuple containing the amplicon number, the alignment segment and whether the primers were correctly paired, or False if the segment is to be skipped
    """
    paired = isinstance(segment, tuple)
    if paired:
        segment1, segment2 = segment
        if not segment1 or not segment2:
            segment = segment1 if segment1 else segment2
            if args.verbose:
                console.debug(
                    f"{segment.query_name}: Pair skipped as at least one segment in pair does not exist"
                )
            return _skip(qc, "missing_mate", segment.reference_name)
        segments = (segment1, segment2)
    else:
        segments = (segment,)
    chrom = segments[0].reference_name

    # filter out unmapped, supplementary and secondary alignment segments
    if not paired:
        if segment.is_unmapped:
            if args.verbose:
                console.debug(f"{segment.query_name}: skipped as unmapped")
            return _skip(qc, "unmapped", chrom)
    else:
        if segment1.is_unmapped or segment2.is_unmapped:
            if args.verbose:
                console.debug(f"{segment1.query_name}: skipped as unmapped")
            return _skip(qc, "unmapped", chrom)

    # supplementary/secondary records are extra alignments of a read counted via its primary record, so no fate is recorded
    if not paired:
        if segment.is_supplementary:
            if args.verbose:
                console.debug(f"{segment.query_name}: skipped as supplementary")
            return False
        # read_pair_generator already drops secondary alignments for paired data
        if segment.is_secondary:
            if args.verbose:
                console.debug(f"{segment.query_name}: skipped as secondary")
            return False
    else:
        if segment1.is_supplementary or segment2.is_supplementary:
            if args.verbose:
                console.debug(f"{segment1.query_name}: skipped as supplementary")
            return False

    if not paired:
        if segment.mapping_quality < min_mapq:
            if args.verbose:
                console.debug(
                    f"{segment.query_name}: skipped as mapping quality below threshold"
                )
            return _skip(qc, "low_mapq", chrom)
    else:
        if segment1.mapping_quality < min_mapq or segment2.mapping_quality < min_mapq:
            if args.verbose:
                console.debug(
                    f"{segment1.query_name}: skipped as mapping quality below threshold"
                )
            return _skip(qc, "low_mapq", chrom)

    if not paired:
        if segment.reference_end is None:
            if args.verbose:
                console.debug(f"{segment.query_name}: skipped as no mapping data")
            return _skip(qc, "no_alignment", chrom)
    else:
        if segment1.reference_end is None or segment2.reference_end is None:
            if args.verbose:
                console.debug(f"{segment1.query_name}: skipped as no mapping data")
            return _skip(qc, "no_alignment", chrom)

    if qc is not None:
        qc.split(segments)

    if not paired:
        # locate the nearest primers to this alignment segment
        p1 = find_primer_with_lookup(
            lookup=lookup,
            pos=segment.reference_start,
            direction="+",
            chrom=segment.reference_name,
        )

        p2 = find_primer_with_lookup(
            lookup=lookup,
            pos=segment.reference_end,
            direction="-",
            chrom=segment.reference_name,
        )
    else:
        # locate the nearest primers to this alignment segment pair
        if segment1.reference_start < segment2.reference_start:
            # if segment1 starts before segment2, then segment1 is the left segment relative to the reference
            p1 = find_primer_with_lookup(
                lookup=lookup,
                pos=segment1.reference_start,
                direction="+",
                chrom=segment1.reference_name,
            )
            p2 = find_primer_with_lookup(
                lookup=lookup,
                pos=segment2.reference_end,
                direction="-",
                chrom=segment2.reference_name,
            )
        else:
            # otherwise then segment2 is the left segment relative to the reference
            p1 = find_primer_with_lookup(
                lookup=lookup,
                pos=segment2.reference_start,
                direction="+",
                chrom=segment2.reference_name,
            )
            p2 = find_primer_with_lookup(
                lookup=lookup,
                pos=segment1.reference_end,
                direction="-",
                chrom=segment1.reference_name,
            )

    left_seg, right_seg = (
        sorted(segments, key=lambda s: s.reference_start)
        if paired
        else (segment, segment)
    )
    if qc is not None:
        qc.same_direction(left_seg, right_seg, p1, p2)

    if not p1 or not p2:
        if paired:
            segment = segment1 if segment1 else segment2
        if args.verbose:
            console.debug(
                f"{segment.query_name}: skipped as no primer found for segment"
            )
        found = p1 or p2
        return _skip(
            qc,
            "no_primer_match",
            found.chrom if found else chrom,
            found.amplicon_number if found else None,
        )

    # check if primers are correctly paired and then assign read group
    correctly_paired = p1.amplicon_number == p2.amplicon_number
    pair_type = None
    if qc is not None:
        pair_type = qc.pair(
            p1, p2, start=left_seg.reference_start, end=right_seg.reference_end
        )

    if not paired:
        if not args.no_read_groups:
            if correctly_paired:
                segment.set_tag("RG", str(p1.pool))
            else:
                segment.set_tag("RG", "unmatched")
    else:
        if not args.no_read_groups:
            if correctly_paired:
                segment1.set_tag("RG", str(p1.pool))
                segment2.set_tag("RG", str(p2.pool))
            else:
                segment1.set_tag("RG", "unmatched")
                segment2.set_tag("RG", "unmatched")

    # get the amplicon number
    amplicon = p1.amplicon_number
    chrom = p1.chrom

    if args.report:
        # update the report with this alignment segment + primer details
        report_segment = segment if not paired else segment1
        report = {
            "chrom": report_segment.reference_name,
            "QueryName": report_segment.query_name,
            "ReferenceStart": report_segment.reference_start,
            "ReferenceEnd": report_segment.reference_end,
            "PrimerPair": f"{p1.primername}_{p2.primername}",
            "Primer1": p1.primername,
            "Primer1Start": p1.start,
            "Primer2": p2.primername,
            "Primer2Start": p2.start,
            "IsSecondary": report_segment.is_secondary,
            "IsSupplementary": report_segment.is_supplementary,
            "Start": p1.start,
            "End": p2.end,
            "CorrectlyPaired": correctly_paired,
        }
        report_writer.writerow(report)

    if not args.allow_incorrect_pairs and not correctly_paired:
        segment = segment if not paired else segment1
        if args.verbose:
            console.debug(f"{segment.query_name}: skipped as not correctly paired")
        return _skip(
            qc,
            "partial_in_overlap"
            if pair_type == "partial_in_overlap"
            else "mismatched_primers",
            chrom,
            amplicon,
        )

    # get the primer positions
    if not args.no_trim_primers:
        p1_position = p1.end
        p2_position = p2.start
    else:
        p1_position = p1.start
        p2_position = p2.end

    # softmask the alignment if left primer start/end inside alignment
    if not paired:
        if segment.reference_start < p1_position:
            try:
                trim(segment, p1_position, False, args.verbose)
                if args.verbose:
                    console.debug(
                        f"{segment.query_name}: ref start {segment.reference_start} >= primer_position {p1_position}"
                    )
            except Exception as e:
                console.read_warning(
                    f"{segment.query_name}: problem soft masking left primer (error: {e}), skipping"
                )
                return _skip(qc, "trim_failed", chrom, amplicon)

        # softmask the alignment if right primer start/end inside alignment
        if segment.reference_end > p2_position:  # type: ignore
            try:
                trim(segment, p2_position, True, args.verbose)
                if args.verbose:
                    console.debug(
                        f"{segment.query_name}: ref start {segment.reference_start} >= primer_position {p2_position}"
                    )
            except Exception as e:
                console.read_warning(
                    f"{segment.query_name}: problem soft masking right primer (error: {e}), skipping"
                )
                return _skip(qc, "trim_failed", chrom, amplicon)

        # check the the alignment still contains bases matching the reference
        if not any(op in matchOps for op, _ in segment.cigartuples):  # type: ignore
            if args.verbose:
                console.debug(
                    f"{segment.query_name}:  dropped as does not match reference post masking"
                )
            return _skip(qc, "primer_only", chrom, amplicon)

        # Check require-full-length
        if args.require_full_length:
            if segment.reference_start > p1.end or segment.reference_end < p2.start:  # type: ignore
                if args.verbose:
                    console.debug(
                        f"{segment.query_name}: ref_start {segment.reference_start} > p1.end {p1.end} or ref_end {segment.reference_end} < p2.start {p2.start}, does not span a full amplicon, skipping"
                    )
                return _skip(qc, "not_full_length", chrom, amplicon)

        if genome_coverage is not None:
            genome_coverage[segment.reference_name][
                segment.reference_start : segment.reference_end
            ] += 1

        # If not normalising, write the segment to the output file and add it to amplicon depth numpy array
        if not args.normalise:
            outfile_writer.write(segment)
            add_fragment_depth(amp_depths[chrom][amplicon], segments, p1.start)
            if qc is not None:
                qc.kept(segments, chrom, amplicon, correctly_paired)

            return (amplicon, False, correctly_paired)

        return (amplicon, segment, correctly_paired)

    else:
        for segment_of_pair in (segment1, segment2):
            if segment_of_pair.reference_start < p1_position:
                try:
                    trim(
                        segment=segment_of_pair,
                        primer_pos=p1_position,
                        end=False,
                        verbose=args.verbose,
                    )
                    if args.verbose:
                        console.debug(
                            f"{segment_of_pair.query_name}: ref start {segment_of_pair.reference_start} >= primer_position {p1_position}"
                        )
                except Exception as e:
                    console.read_warning(
                        f"{segment_of_pair.query_name}: Problem soft masking left primer (error: {e}), skipping"
                    )
                    return _skip(qc, "trim_failed", chrom, amplicon)

            if segment_of_pair.reference_end > p2_position:  # type: ignore
                try:
                    trim(
                        segment=segment_of_pair,
                        primer_pos=p2_position,
                        end=True,
                        verbose=args.verbose,
                    )
                    if args.verbose:
                        console.debug(
                            f"{segment_of_pair.query_name}: ref_end {segment_of_pair.reference_end} >= primer_position {p2_position}"
                        )
                except Exception as e:
                    console.read_warning(
                        f"{segment_of_pair.query_name}: Problem soft masking right primer (error: {e}), skipping"
                    )
                    return _skip(qc, "trim_failed", chrom, amplicon)

        # check the the alignment still contains bases matching the reference
        if not all(
            any(op in matchOps for op, _ in seg.cigartuples)  # type: ignore
            for seg in segments
        ):
            if args.verbose:
                console.debug(
                    f"{segment1.query_name}: Paired segment dropped as does not match reference post masking"
                )
            return _skip(qc, "primer_only", chrom, amplicon)

        if args.require_full_length:
            if segment1.reference_start < segment2.reference_start:
                if (
                    segment1.reference_start > p1.end  # type: ignore
                    or segment2.reference_end < p2.start  # type: ignore
                ):
                    if args.verbose:
                        console.debug(
                            f"{segment1.query_name}: ref_start {segment1.reference_start} > p1.end {p1.end} or ref_end {segment2.reference_end} < p2.start {p2.start}, does not span a full amplicon, skipping"
                        )
                    return _skip(qc, "not_full_length", chrom, amplicon)
            else:
                if (
                    segment2.reference_start > p1.end
                    or segment1.reference_end < p2.start  # type: ignore
                ):
                    if args.verbose:
                        console.debug(
                            f"{segment1.query_name}: ref_end {segment1.reference_end} < p2.start {p2.start} or ref_start {segment2.reference_start} > p1.end {p1.end}, does not span a full amplicon, skipping"
                        )
                    return _skip(qc, "not_full_length", chrom, amplicon)

        if genome_coverage is not None:
            for seg in (segment1, segment2):
                genome_coverage[seg.reference_name][
                    seg.reference_start : seg.reference_end
                ] += 1

        # If not normalising, write the segments to the output file and add the fragment to amplicon depth numpy array
        if not args.normalise:
            outfile_writer.write(segment1)
            outfile_writer.write(segment2)
            add_fragment_depth(amp_depths[chrom][amplicon], segments, p1.start)
            if qc is not None:
                qc.kept(segments, chrom, amplicon, correctly_paired)

            return (amplicon, False, correctly_paired)

    return (amplicon, segment, correctly_paired)


def read_pair_generator(bam, region_string=None, qc=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found. Only proper pairs are yielded, pairs which are
    not (and mates which never find a partner) are recorded in the QC read fate report if requested.
    """

    def chrom_of(*reads):
        mapped = [r for r in reads if r is not None and not r.is_unmapped]
        return mapped[0].reference_name if mapped else None

    read_dict = defaultdict(lambda: [None, None])
    for read in bam:
        # supplementary/secondary records share the qname and would be paired in place of the real mate
        if read.is_supplementary or read.is_secondary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                pair = (read, read_dict[qname][1])
            else:
                pair = (read_dict[qname][0], read)
            del read_dict[qname]
            if all(r is not None and r.is_proper_pair for r in pair):
                yield pair
            elif any(r is None or r.is_unmapped for r in pair):
                _skip(qc, "unmapped", chrom_of(*pair))
            else:
                _skip(qc, "not_proper_pair", chrom_of(*pair))

    # mates whose partner never appeared (e.g. filtered out upstream)
    for read1, read2 in read_dict.values():
        read = read1 or read2
        if read.is_unmapped or read.mate_is_unmapped:  # type: ignore
            _skip(qc, "unmapped", chrom_of(read))
        else:
            _skip(qc, "missing_mate", chrom_of(read))


def create_primer_lookup(
    ref_len_tuple, amplicons: list[Amplicon], padding=DEFAULT_PRIMER_MATCH_THRESHOLD
):
    """
    Create a lookup table for efficient primer position queries across reference genomes.

    Each chromosome gets its own 2D lookup array where:
    - Rows represent non-overlapping "pools"* of amplicons at their corresponding positions.
    - Columns represent genomic positions
    - Values are Amplicon objects or None

    The function automatically determines the minimum number of rows needed to ensure
    no amplicons overlap within the same row when accounting for padding.

    * Amplicons are placed in the first available row where they don't overlap, not their pool index.

    Parameters
    ----------
    ref_len_tuple : list[tuple[str, int]]
        List of tuples containing (chromosome_name, chromosome_length) pairs
        from the reference genome
    amplicons : list[Amplicon]
        List of Amplicon objects containing primer scheme information
    padding : int, optional
        Number of bases to extend amplicon boundaries on both sides, so reads ending
        slightly outside a primer 5' end are still assigned to it (default: 5)

    Returns
    -------
    dict[str, np.ndarray]
        Dictionary mapping chromosome names to 2D numpy arrays of shape (N, chrom_len+1)
        where N is the minimum number of rows needed to prevent amplicon overlap.
        Array elements are either Amplicon objects or None.


    """
    lookups = {}
    for chrom, chromlen in ref_len_tuple:
        lookup_array = np.empty_like(None, shape=(1, chromlen + 1))
        for amp in amplicons:
            added = False
            if amp.chrom == chrom:
                # If amplicon clashes with any in same pool add new row
                amp_slice = lookup_array[
                    :,
                    max(amp.amplicon_start - padding, 0) : min(
                        amp.amplicon_end + padding, chromlen
                    ),
                ]
                for i, row in enumerate(amp_slice):  # Check each row for collision
                    if row[row != None].size == 0:
                        lookup_array[
                            i,
                            max(amp.amplicon_start - padding, 0) : min(
                                amp.amplicon_end + padding, chromlen
                            ),
                        ] = amp
                        added = True
                # If not added, create new row, add the amplicon to that then add back to original array
                if not added:
                    new_row = np.empty_like(None, shape=(1, chromlen + 1))
                    new_row[
                        0,
                        max(amp.amplicon_start - padding, 0) : min(
                            amp.amplicon_end + padding, chromlen
                        ),
                    ] = amp
                    lookup_array = np.vstack((lookup_array, new_row))

        lookups[chrom] = lookup_array
    return lookups


COVERAGE_THRESHOLDS = [1, 10, 20, 100, 1000]


def write_genome_coverage(filepath, genome_coverage, label=None):
    """Write per-position coverage and return (label, total positions, [(threshold, count)])."""
    total_len = sum(len(arr) for arr in genome_coverage.values())
    with open(filepath, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["chrom", "pos", "depth"], delimiter="\t"
        )
        writer.writeheader()
        for chrom, depths in genome_coverage.items():
            for i, d in enumerate(depths):
                writer.writerow({"chrom": chrom, "pos": i + 1, "depth": int(d)})
    all_depths = np.concatenate(list(genome_coverage.values()))
    counts = [(t, int(np.sum(all_depths >= t))) for t in COVERAGE_THRESHOLDS]
    return label, total_len, counts


def print_coverage_summaries(summaries):
    """Print genome coverage summaries: a table on a terminal, prefixed lines in a log."""
    console.heading("Coverage")
    if not console.live:
        for label, total_len, counts in summaries:
            prefix = f"{label}: " if label else ""
            console.line(f"{prefix}Genome coverage summary")
            console.line(f"{prefix}  Total positions: {total_len}")
            for t, count in counts:
                pct = (count / total_len) * 100 if total_len > 0 else 0
                console.line(f"{prefix}  >= {t}x: {count}/{total_len} ({pct:.2f}%)")
        return

    col = 24
    header = "".join(f"{(label or '').lower():>{col}}" for label, _, _ in summaries)
    console.line(console.style(f"{'':<11}{header}", BOLD))
    console.line(
        f"{'positions':<11}"
        + "".join(f"{total_len:>{col}}" for _, total_len, _ in summaries)
    )
    for i, t in enumerate(COVERAGE_THRESHOLDS):
        cells = []
        for _, total_len, counts in summaries:
            count = counts[i][1]
            pct = (count / total_len) * 100 if total_len > 0 else 0
            cell = f"{count}  {pct:6.2f}%"
            cells.append(
                " " * (col - len(cell))
                + f"{count}  "
                + console.style(f"{pct:6.2f}%", BOLD)
            )
        console.line(f"{f'>= {t}x':<11}" + "".join(cells))


def go(args):
    """Filter and soft mask an alignment file so that the alignment boundaries match the primer start and end sites.

    Based on the most likely primer position, based on the alignment coordinates.
    """
    console.reset(quiet=getattr(args, "quiet", False), verbose=args.verbose)

    # guard for negative normalise
    if args.normalise is not None and args.normalise < 0:
        console.error("normalise must be >= 0, exiting.")
        sys.exit(1)

    # prepare the report outfile
    if args.report:
        reportfh = open(args.report, "w")
        report_headers = [
            "chrom",
            "QueryName",
            "ReferenceStart",
            "ReferenceEnd",
            "PrimerPair",
            "Primer1",
            "Primer1Start",
            "Primer2",
            "Primer2Start",
            "IsSecondary",
            "IsSupplementary",
            "Start",
            "End",
            "CorrectlyPaired",
        ]
        report_writer = csv.DictWriter(
            reportfh, fieldnames=report_headers, delimiter="\t"
        )
        report_writer.writeheader()

    # open the primer scheme and get the pools
    scheme = Scheme.from_file(args.bedfile)
    # keep the unmerged primers (including alts) for primer sequence lookups
    unmerged_bedlines = list(scheme.bedlines)

    # Merge the primers
    scheme.bedlines = merge_primers(scheme.bedlines)

    amplicon_list = create_amplicons(scheme.bedlines)
    amplicons = {}
    for amplicon in amplicon_list:
        amplicon.length = amplicon.amplicon_end - amplicon.amplicon_start  # type: ignore
        amplicons.setdefault(amplicon.chrom, {})[amplicon.amplicon_number] = amplicon

    pools = set([bl.pool for bl in scheme.bedlines])

    pools_str = {str(x) for x in pools}
    pools_str.add("unmatched")

    # open the input samfile and process read groups
    if args.samfile and args.samfile != "-":
        infile = pysam.AlignmentFile(args.samfile, "rb")
    else:
        infile = pysam.AlignmentFile("-", "rb")

    first_segment = next(infile, None)
    if not first_segment:
        console.error("No segments found in the input file, exiting.")
        sys.exit(1)

    # check if the first segment is paired, then chain the saved first segment with the infile iterator so nothing is lost
    paired = first_segment.is_paired
    chained_iterator = itertools.chain([first_segment], infile)

    bam_header = infile.header.copy().to_dict()
    if not args.no_read_groups:
        bam_header["RG"] = []
        for pool in sorted(pools_str):  # set order can be non deterministic
            read_group = {}
            read_group["ID"] = pool
            bam_header["RG"].append(read_group)

    cli_cmd = " ".join(sys.argv)
    bam_header["PG"].append(
        {
            "PN": "align_trim",
            "ID": "align_trim",
            "VN": version("align_trim"),
            "CL": cli_cmd,
        }
    )

    # prepare the alignment outfile
    if args.output and args.output != "-":
        if args.output.name.endswith(".bam"):
            outfile = pysam.AlignmentFile(args.output, "wb", header=bam_header)
        elif args.output.name.endswith(".sam"):
            outfile = pysam.AlignmentFile(args.output, "wh", header=bam_header)
        else:
            console.error(
                "Output file path must end with either .bam or .sam, exiting."
            )
            sys.exit(1)

    else:
        outfile = pysam.AlignmentFile("-", "wh", header=bam_header)

    if not console.quiet:
        console.heading("Run")
        console.line(f"align_trim {version('align_trim')}")
        reads = "paired-end" if paired else "single-end"
        console.line(f"  input   {args.samfile or 'stdin'} ({reads})")
        console.line(
            f"  scheme  {args.bedfile} ({len(amplicon_list)} amplicons, "
            f"{len(amplicons)} reference{'s' if len(amplicons) != 1 else ''})"
        )
        console.line(f"  output  {args.output or 'stdout'}")

    # a live bar would be drawn in among SAM records written to the same terminal
    sam_on_terminal = (
        not args.output or str(args.output) == "-"
    ) and sys.stdout.isatty()
    progress = console.progress(
        fraction=bam_fraction(infile, args.samfile),
        enabled=not console.quiet,
        live=not sam_on_terminal,
    )
    chained_iterator = progress.track(chained_iterator)
    n_pairs = 0

    # Initialise the amplicon depth dict
    amp_depths = {}
    for amp in amplicon_list:
        amp_depths.setdefault(amp.chrom, {})
        amp_depths[amp.chrom].setdefault(
            amp.amplicon_number,
            np.zeros(amp.length, dtype=int),  # type: ignore
        )

    # Initialise the mean depths dictionary, this will get stomped over if normalisation is requested
    mean_amp_depths = {}
    for chrom in amplicons:
        for amplicon in amplicons[chrom]:
            mean_amp_depths[(chrom, amplicon)] = 0

    # Create a lookup table for primer location
    ref_lengths = [(r, infile.get_reference_length(r)) for r in infile.references]
    primer_lookup = create_primer_lookup(
        ref_len_tuple=ref_lengths,
        amplicons=amplicon_list,
        padding=args.primer_match_threshold,
    )

    qc = QC.from_args(args, amplicon_list, primer_lookup, unmerged_bedlines)

    # Initialise genome-wide coverage arrays if requested
    genome_coverage_pre = None
    genome_coverage_post = None
    if args.genome_coverage_report:
        genome_coverage_pre = {}
        for ref_name, ref_len in ref_lengths:
            genome_coverage_pre[ref_name] = np.zeros(ref_len, dtype=int)
        if args.normalise:
            genome_coverage_post = {}
            for ref_name, ref_len in ref_lengths:
                genome_coverage_post[ref_name] = np.zeros(ref_len, dtype=int)

    # Per-amplicon normalisation state: running depth array and current MAD from target
    if args.normalise:
        norm_state = {}
        for amp in amplicon_list:
            norm_state[(amp.chrom, amp.amplicon_number)] = {
                "depth": np.zeros(amp.length, dtype=int),  # type: ignore
                "distance": float(args.normalise),
            }

    if paired:
        read_pairs = read_pair_generator(chained_iterator, qc=qc)

        for segments in read_pairs:
            n_pairs += 1
            if args.report:
                trimming_tuple = handle_segments(
                    segment=segments,  # type: ignore
                    lookup=primer_lookup,
                    args=args,
                    report_writer=report_writer,  # type: ignore
                    min_mapq=args.min_mapq,
                    outfile_writer=outfile,
                    amp_depths=amp_depths,
                    genome_coverage=genome_coverage_pre,
                    qc=qc,
                )
            else:
                trimming_tuple = handle_segments(
                    segment=segments,  # type: ignore
                    lookup=primer_lookup,
                    args=args,
                    min_mapq=args.min_mapq,
                    outfile_writer=outfile,
                    amp_depths=amp_depths,
                    genome_coverage=genome_coverage_pre,
                    qc=qc,
                )

            if not trimming_tuple:
                continue

            # unpack the trimming tuple since segment passed trimming
            amplicon, trimmed_pair, correctly_paired = trimming_tuple

            # If we aren't normalising the segments will have already been written to the outfile
            if not args.normalise and not trimmed_pair:
                continue

            if args.normalise and trimmed_pair:
                chrom = trimmed_pair[0].reference_name  # type: ignore
                state = norm_state[(chrom, amplicon)]
                p_start = amplicons[chrom][amplicon].amplicon_start
                test_depths = np.copy(state["depth"])
                # count the fragment once where the mates overlap
                add_fragment_depth(test_depths, trimmed_pair, p_start)
                test_distance = np.mean(np.abs(test_depths - args.normalise))
                if test_distance < state["distance"]:
                    state["depth"] = test_depths
                    state["distance"] = test_distance
                    outfile.write(trimmed_pair[0])  # type: ignore
                    outfile.write(trimmed_pair[1])  # type: ignore
                    if genome_coverage_post is not None:
                        for seg in trimmed_pair:  # type: ignore
                            genome_coverage_post[seg.reference_name][
                                seg.reference_start : seg.reference_end
                            ] += 1
                    if qc is not None:
                        qc.kept(trimmed_pair, chrom, amplicon, correctly_paired)
                elif qc is not None:
                    qc.fate("by_normalise", chrom, amplicon)

        if args.normalise:
            mean_amp_depths = {k: np.mean(v["depth"]) for k, v in norm_state.items()}
        else:
            mean_amp_depths = {}
            for chrom, chrom_amps in amp_depths.items():
                for amplicon, depths in chrom_amps.items():
                    mean_amp_depths[(chrom, amplicon)] = np.mean(depths)

        # write mean amplicon depths to file
        if args.amp_depth_report:
            with open(args.amp_depth_report, "w") as amp_depth_report_fh:
                writer = csv.DictWriter(
                    amp_depth_report_fh,
                    fieldnames=["chrom", "amplicon", "mean_fragment_depth"],
                    delimiter="\t",
                )
                writer.writeheader()
                for (chrom, amplicon), depth in mean_amp_depths.items():
                    writer.writerow(
                        {
                            "chrom": chrom,
                            "amplicon": amplicon,
                            "mean_fragment_depth": depth,
                        }
                    )

    else:
        # iterate over the alignment segments in the input SAM file
        for segment in chained_iterator:
            if args.report:
                trimming_tuple = handle_segments(
                    segment=segment,
                    args=args,
                    report_writer=report_writer,  # type: ignore
                    min_mapq=args.min_mapq,
                    lookup=primer_lookup,
                    outfile_writer=outfile,
                    amp_depths=amp_depths,
                    genome_coverage=genome_coverage_pre,
                    qc=qc,
                )

            else:
                trimming_tuple = handle_segments(
                    segment=segment,
                    args=args,
                    min_mapq=args.min_mapq,
                    lookup=primer_lookup,
                    outfile_writer=outfile,
                    amp_depths=amp_depths,
                    genome_coverage=genome_coverage_pre,
                    qc=qc,
                )

            if not trimming_tuple:
                continue

            # unpack the trimming tuple since segment passed trimming
            amplicon, trimmed_segment, correctly_paired = trimming_tuple

            # If we aren't normalising the segments will have already been written to the outfile
            if not args.normalise and not trimmed_segment:
                continue

            if args.normalise and trimmed_segment:
                chrom = trimmed_segment.reference_name  # type: ignore
                state = norm_state[(chrom, amplicon)]
                p_start = amplicons[chrom][amplicon].amplicon_start
                test_depths = np.copy(state["depth"])
                relative_start = max(0, trimmed_segment.reference_start - p_start)  # type: ignore
                relative_end = trimmed_segment.reference_end - p_start  # type: ignore
                test_depths[relative_start:relative_end] += 1
                test_distance = np.mean(np.abs(test_depths - args.normalise))
                if test_distance < state["distance"]:
                    state["depth"] = test_depths
                    state["distance"] = test_distance
                    outfile.write(trimmed_segment)  # type: ignore
                    if genome_coverage_post is not None:
                        genome_coverage_post[trimmed_segment.reference_name][  # type: ignore
                            trimmed_segment.reference_start : trimmed_segment.reference_end  # type: ignore
                        ] += 1
                    if qc is not None:
                        qc.kept((trimmed_segment,), chrom, amplicon, correctly_paired)
                elif qc is not None:
                    qc.fate("by_normalise", chrom, amplicon)

        # normalise if requested
        if args.normalise:
            mean_amp_depths = {k: np.mean(v["depth"]) for k, v in norm_state.items()}

        else:
            mean_amp_depths = {}
            for chrom, chrom_amps in amp_depths.items():
                for amplicon, depths in chrom_amps.items():
                    mean_amp_depths[(chrom, amplicon)] = np.mean(depths)

        # write mean amplicon depths to file
        if args.amp_depth_report:
            with open(args.amp_depth_report, "w") as amp_depth_report_fh:
                writer = csv.DictWriter(
                    amp_depth_report_fh,
                    fieldnames=["chrom", "amplicon", "mean_fragment_depth"],
                    delimiter="\t",
                )
                writer.writeheader()

                for (chrom, amplicon), depth in mean_amp_depths.items():
                    writer.writerow(
                        {
                            "chrom": chrom,
                            "amplicon": amplicon,
                            "mean_fragment_depth": depth,
                        }
                    )

    progress.done(f"{n_pairs:,} read pairs" if paired else None)

    # Write genome coverage reports
    if args.genome_coverage_report:
        pre_path = f"{args.genome_coverage_report}.pre-normalisation.coverage.tsv"
        summaries = [
            write_genome_coverage(
                pre_path, genome_coverage_pre, label="Pre-normalisation"
            )
        ]
        if args.normalise and genome_coverage_post is not None:
            post_path = f"{args.genome_coverage_report}.post-normalisation.coverage.tsv"
            summaries.append(
                write_genome_coverage(
                    post_path, genome_coverage_post, label="Post-normalisation"
                )
            )
        print_coverage_summaries(summaries)

    if qc is not None:
        qc.write(args.qc_prefix)

    console.flush_warnings()

    # close up the file handles
    infile.close()
    outfile.close()
    if args.report:
        reportfh.close()  # type: ignore


def main():
    parser = argparse.ArgumentParser(
        description="Trim alignments from an amplicon scheme. Bam (input) can be provided by --samfile or stdin"
    )
    parser.add_argument(
        "bedfile",
        help="BED file containing the amplicon scheme",
        type=Path,
        metavar="BEDFILE",
    )
    parser.add_argument(
        "--samfile",
        "-s",
        help="Sorted SAM/BAM file containing the aligned reads, if this is not provided (or '-') then 'align_trim' will read from stdin.",
        required=False,
    )
    parser.add_argument(
        "--normalise",
        "-n",
        type=int,
        help="Subsample to N coverage per amplicon, counted in fragments (a read pair counts once where its mates overlap). Use 0 for no normalisation. (default: %(default)s)",
        default=0,
    )
    parser.add_argument(
        "--min-mapq",
        "-m",
        type=int,
        default=20,
        help="Minimum mapping quality to keep an aligned read (default: %(default)s)",
    )
    parser.add_argument(
        "--primer-match-threshold",
        "-p",
        type=int,
        default=DEFAULT_PRIMER_MATCH_THRESHOLD,
        help="Max distance in bases between a read end and a primer 5' end for the read to be assigned to that primer. Also used by the --qc-prefix reports to decide whether a read end is at a primer. (default: %(default)s)",
    )
    parser.add_argument(
        "--report", "-r", type=Path, help="Output report TSV to filepath"
    )
    parser.add_argument(
        "--amp-depth-report",
        "-a",
        type=Path,
        help="Output amplicon depth TSV (mean fragment depth per amplicon) to filepath",
    )
    parser.add_argument(
        "--genome-coverage-report",
        "-g",
        type=str,
        default=None,
        metavar="PREFIX",
        help="Output per-position genome coverage TSV(s) using PREFIX. Produces PREFIX.pre-normalisation.coverage.tsv (always) and PREFIX.post-normalisation.coverage.tsv (when --normalise is used). Summary stats are printed to stderr.",
    )
    parser.add_argument(
        "--qc-prefix",
        type=str,
        default=None,
        metavar="PREFIX",
        help="Output amplicon QC reports as PREFIX.<name>.tsv (summary, amplicons, primer_pairs, split_reads, overlaps). Summary lines are printed to stderr.",
    )
    parser.add_argument(
        "--no-trim-primers",
        action="store_true",
        help="Do not trim primers from reads",
    )
    parser.add_argument(
        "--no-read-groups",
        dest="no_read_groups",
        help="Do not divide reads into groups in samfile output",
        action="store_true",
    )
    parser.add_argument(
        "--allow-incorrect-pairs",
        action="store_true",
        help="Allow reads to be assigned to amplicons even if the primers are not correctly paired, i.e. primer1 and primer2 are not from the same amplicon.",
    )
    parser.add_argument(
        "--require-full-length",
        action="store_true",
        help="Requires all reads to start and stop in a primer site, do not use this option if you are using rapid barcoding since the reads will not be full length.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=None,
        metavar="OUTPUT",
        help="Location to write the output samfile to, the output type will be determined by the file extension. If no <OUTPUT> or '-' provided, will write plaintext samfile to stdout",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Debug mode: print per-read decisions to stderr (disables the progress bar)",
    )
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Only print the coverage summary, QC summary and warnings to stderr (no run details or progress)",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {version('align_trim')}",
        help="Show the version of align_trim",
    )

    args = parser.parse_args()

    go(args)


if __name__ == "__main__":
    main()
