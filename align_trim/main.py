from copy import copy
import csv
import pysam
import sys
import numpy as np
import random
import argparse
from collections import defaultdict
from typing import Optional
from pathlib import Path
import itertools
from typing import Union

from importlib.metadata import version

from primalbedtools.scheme import Scheme
from primalbedtools.bedfiles import BedLine, merge_primers
from primalbedtools.amplicons import Amplicon, create_amplicons

RANDOM_SEED = 42


# consumesReference lookup for if a CIGAR operation consumes the reference sequence
consumesReference = [True, False, True, True, False, False, False, True]

# consumesQuery lookup for if a CIGAR operation consumes the query sequence
consumesQuery = [True, True, False, False, True, False, False, True]


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
        print(
            f"{segment.query_name}: Trimming {'end' if end else 'start'} of read to primer position {primer_pos}",
            file=sys.stderr,
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
                print(
                    f"{segment.query_name}: Chomped a {flag}, {length}",
                    file=sys.stderr,
                )
        except IndexError:
            if verbose:
                print(
                    f"{segment.query_name}: Ran out of cigar during soft masking - completely masked read will be ignored",
                    file=sys.stderr,
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
        if not end and pos >= primer_pos and flag == 0:
            break
        if end and pos <= primer_pos and flag == 0:
            break

    # calculate how many extra matches are needed in the CIGAR
    extra = abs(pos - primer_pos)
    if verbose:
        print(f"{segment.query_name}: extra {extra}", file=sys.stderr)
    if extra:
        if verbose:
            print(
                f"{segment.query_name}: Inserted a 0, {extra}",
                file=sys.stderr,
            )
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
            print(
                f"{segment.query_name}: New pos - {segment.pos}",
                file=sys.stderr,
            )

        # if proposed softmask leads straight into a deletion, shuffle leftmost mapping base along and ignore the deletion
        if cigar[0][0] == 2:
            if verbose:
                print(
                    f"{segment.query_name}: softmask created a leading deletion in the CIGAR, shuffling the alignment",
                    file=sys.stderr,
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
        cigar.append((4, eaten))

    # check the new CIGAR and replace the old one
    if cigar[0][1] <= 0 or cigar[-1][1] <= 0:
        if verbose:
            print(
                f"{segment.query_name}: invalid cigar operation created - possibly due to INDEL in primer",
                file=sys.stderr,
            )
        return

    segment.cigartuples = cigar
    return


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
):
    """Handle the alignment segment(s) including filtering, soft masking, and reporting.

    Args:
        segment (pysam.AlignedSegment | tuple): The alignment segment to process, can be a single segment or a tuple of paired segments
        bed (dict): The primer scheme
        reportfh (typing.IO): The report file handle
        args (argparse.Namespace): The command line arguments

    Returns:
        tuple [int, pysam.AlignedSegment | bool] | bool: A tuple containing the amplicon number and the alignment segment, or False if the segment is to be skipped
    """
    paired = isinstance(segment, tuple)
    if paired:
        segment1, segment2 = segment
        if not segment1 or not segment2:
            segment = segment1 if segment1 else segment2
            if args.verbose:
                print(
                    f"{segment.query_name}: Pair skipped as at least one segment in pair does not exist",
                    file=sys.stderr,
                )
            return False

    # filter out unmapped and supplementary alignment segments
    if not paired:
        if segment.is_unmapped:
            if args.verbose:
                print(
                    f"{segment.query_name}: skipped as unmapped",
                    file=sys.stderr,
                )
            return False
    else:
        if segment1.is_unmapped or segment2.is_unmapped:
            if args.verbose:
                print(
                    f"{segment1.query_name}: skipped as unmapped",
                    file=sys.stderr,
                )
            return False

    if not paired:
        if segment.is_supplementary:
            if args.verbose:
                print(
                    f"{segment.query_name}: skipped as supplementary",
                    file=sys.stderr,
                )
            return False
    else:
        if segment1.is_supplementary or segment2.is_supplementary:
            if args.verbose:
                print(
                    f"{segment1.query_name}: skipped as supplementary",
                    file=sys.stderr,
                )
            return False

    if not paired:
        if segment.mapping_quality < min_mapq:
            if args.verbose:
                print(
                    f"{segment.query_name}: skipped as mapping quality below threshold",
                    file=sys.stderr,
                )
            return False
    else:
        if segment1.mapping_quality < min_mapq or segment2.mapping_quality < min_mapq:
            if args.verbose:
                print(
                    f"{segment1.query_name}: skipped as mapping quality below threshold",
                    file=sys.stderr,
                )
            return False

    if not paired:
        if segment.reference_end is None:
            if args.verbose:
                print(
                    f"{segment.query_name}: skipped as no mapping data",
                    file=sys.stderr,
                )
            return False
    else:
        if segment1.reference_end is None or segment2.reference_end is None:
            if args.verbose:
                print(
                    f"{segment1.query_name}: skipped as no mapping data",
                    file=sys.stderr,
                )
            return False
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

    if not p1 or not p2:
        segment = segment1 if segment1 else segment2
        if args.verbose:
            print(
                f"{segment.query_name}: skipped as no primer found for segment",
                file=sys.stderr,
            )
        return False

    # check if primers are correctly paired and then assign read group
    correctly_paired = p1.amplicon_number == p2.amplicon_number

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
            print(
                f"{segment.query_name}: skipped as not correctly paired",
                file=sys.stderr,
            )
        return False

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
                    print(
                        f"{segment.query_name}: ref start {segment.reference_start} >= primer_position {p1_position}",
                        file=sys.stderr,
                    )
            except Exception as e:
                print(
                    f"{segment.query_name}: problem soft masking left primer (error: {e}), skipping",
                    file=sys.stderr,
                )
                return False

        # softmask the alignment if right primer start/end inside alignment
        if segment.reference_end > p2_position:
            try:
                trim(segment, p2_position, True, args.verbose)
                if args.verbose:
                    print(
                        f"{segment.query_name}: ref start {segment.reference_start} >= primer_position {p2_position}",
                        file=sys.stderr,
                    )
            except Exception as e:
                print(
                    f"{segment.query_name}: problem soft masking right primer (error: {e}), skipping",
                    file=sys.stderr,
                )
                return False

        # check the the alignment still contains bases matching the reference
        if "M" not in segment.cigarstring:  # type: ignore
            if args.verbose:
                print(
                    f"{segment.query_name}:  dropped as does not match reference post masking",
                    file=sys.stderr,
                )
            return False

        # Check require-full-length
        if args.require_full_length:
            if segment.reference_start > p1.end or segment.reference_end < p2.start:
                if args.verbose:
                    print(
                        f"{segment.query_name}: ref_start {segment.reference_start} > p1.end {p1.end} or ref_end {segment.reference_end} < p2.start {p2.start}, does not span a full amplicon, skipping",
                        file=sys.stderr,
                    )
                return False

        # If not normalising, write the segment to the output file and add it to amplicon depth numpy array
        if not args.normalise:
            outfile_writer.write(segment)
            segment_amp_relative_start = segment.reference_start - p1.start
            segment_amp_relative_end = segment.reference_end - p1.start
            if segment_amp_relative_start < 0:
                segment_amp_relative_start = 0

            amp_depths[segment.reference_name][amplicon][
                segment_amp_relative_start:segment_amp_relative_end
            ] += 1

            return (amplicon, False)

        return (amplicon, segment)

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
                        print(
                            f"{segment_of_pair.query_name}: ref start {segment_of_pair.reference_start} >= primer_position {p1_position}",
                            file=sys.stderr,
                        )
                except Exception as e:
                    print(
                        f"{segment_of_pair.query_name}: Problem soft masking left primer (error: {e}), skipping",
                        file=sys.stderr,
                    )
                    return False

            if segment_of_pair.reference_end > p2_position:  # type: ignore
                try:
                    trim(
                        segment=segment_of_pair,
                        primer_pos=p2_position,
                        end=True,
                        verbose=args.verbose,
                    )
                    if args.verbose:
                        print(
                            f"{segment_of_pair.query_name}: ref_end {segment_of_pair.reference_end} >= primer_position {p2_position}",
                            file=sys.stderr,
                        )
                except Exception as e:
                    print(
                        f"{segment_of_pair.query_name}: Problem soft masking right primer (error: {e}), skipping",
                        file=sys.stderr,
                    )
                    return False

        # check the the alignment still contains bases matching the reference
        if "M" not in segment1.cigarstring or "M" not in segment2.cigarstring:  # type: ignore
            if args.verbose:
                print(
                    f"{segment1.query_name}: Paired segment dropped as does not match reference post masking",
                    file=sys.stderr,
                )
            return False

        if args.require_full_length:
            if segment1.reference_start < segment2.reference_start:
                if (
                    segment1.reference_start > p1.end
                    or segment2.reference_end < p2.start
                ):
                    if args.verbose:
                        print(
                            f"{segment1.query_name}: ref_start {segment1.reference_start} > p1.end {p1.end} or ref_end {segment2.reference_end} < p2.start {p2.start}, does not span a full amplicon, skipping",
                            file=sys.stderr,
                        )
                    return False
            else:
                if (
                    segment2.reference_start > p1.end
                    or segment1.reference_end < p2.start
                ):
                    if args.verbose:
                        print(
                            f"{segment1.query_name}: ref_end {segment1.reference_end} < p2.start {p2.start} or ref_start {segment2.reference_start} > p1.end {p1.end}, does not span a full amplicon, skipping",
                            file=sys.stderr,
                        )
                    return False

        # If not normalising, write the segments to the output file and add them to amplicon depth numpy array
        if not args.normalise:
            outfile_writer.write(segment1)
            outfile_writer.write(segment2)
            for segment_in_pair in (segment1, segment2):
                segment_amp_relative_start = segment_in_pair.reference_start - p1.start
                segment_amp_relative_end = segment_in_pair.reference_end - p1.start
                if segment_amp_relative_start < 0:
                    segment_amp_relative_start = 0
            amp_depths[segment1.reference_name][amplicon][
                segment_amp_relative_start:segment_amp_relative_end
            ] += 1

            return (amplicon, False)

    return (amplicon, segment)


def normalise(
    trimmed_segments: dict,
    normalise: int,
    primers: list[BedLine],
    outfile: pysam.AlignmentFile,
    verbose: bool = False,
):
    """Normalise the depth of the trimmed segments to a given value. Perform per-amplicon normalisation using numpy vector maths to determine whether the segment in question would take the depth closer to the desired depth across the amplicon.

    Args:
        trimmed_segments (dict): Dict containing amplicon number as key and list of pysam.AlignedSegment as value, if paired segments are used, the value will be a list of tuples containing the two segments.
        normalise (int): Desired normalised depth
        bed (list): Primer scheme as a list of BedLine objects
        outfile (pysam.AlignmentFile): Output file handle to write the normalised segments to
        verbose (bool): If True, will print normalisation info during processing

    Raises:
        ValueError: Amplicon assigned to segment not found in primer scheme file

    Returns:
        dict: A dictionary containing the mean depth for each amplicon post normalisation
    """
    amplicons = {}

    for amplicon in create_amplicons(primers):
        amplicons.setdefault(amplicon.chrom, {})
        amplicons[amplicon.chrom].setdefault(
            amplicon.amplicon_number,
            {
                "length": amplicon.amplicon_end - amplicon.amplicon_start,
                "p_start": amplicon.amplicon_start,
            },
        )

    # mean_depths = {x: {} for x in amplicons}
    mean_depths = {}
    for chrom in amplicons:
        for amplicon in amplicons[chrom]:
            mean_depths[(chrom, amplicon)] = 0

    for chrom, amplicon_dict in trimmed_segments.items():
        for amplicon, segments in amplicon_dict.items():
            if amplicon not in amplicons[chrom]:
                raise ValueError(f"Amplicon {amplicon} not found in primer scheme file")

            desired_depth = np.full_like(
                (amplicons[chrom][amplicon]["length"],), normalise, dtype=int
            )

            amplicon_depth = np.zeros(
                (amplicons[chrom][amplicon]["length"],), dtype=int
            )

            if not segments:
                if verbose:
                    print(
                        f"No segments assigned to amplicon {amplicon}, skipping",
                        file=sys.stderr,
                    )
                continue

            random.Random(RANDOM_SEED).shuffle(segments)

            distance = np.mean(np.abs(amplicon_depth - desired_depth))

            for segment in segments:
                paired = isinstance(segment, tuple)

                if paired:
                    test_depths = np.copy(amplicon_depth)
                    segment1, segment2 = segment
                    for segment in (segment1, segment2):
                        relative_start = (
                            segment.reference_start
                            - amplicons[chrom][amplicon]["p_start"]
                        )

                        if relative_start < 0:
                            relative_start = 0

                        relative_end = (
                            segment.reference_end
                            - amplicons[chrom][amplicon]["p_start"]
                        )

                        test_depths[relative_start:relative_end] += 1

                    test_distance = np.mean(np.abs(test_depths - desired_depth))

                    if test_distance < distance:
                        amplicon_depth = test_depths
                        distance = test_distance
                        # write the segments to the output file
                        outfile.write(segment1)
                        outfile.write(segment2)
                else:
                    test_depths = np.copy(amplicon_depth)

                    relative_start = (
                        segment.reference_start - amplicons[chrom][amplicon]["p_start"]
                    )

                    if relative_start < 0:
                        relative_start = 0

                    relative_end = (
                        segment.reference_end - amplicons[chrom][amplicon]["p_start"]
                    )

                    test_depths[relative_start:relative_end] += 1

                    test_distance = np.mean(np.abs(test_depths - desired_depth))

                    if test_distance < distance:
                        amplicon_depth = test_depths
                        distance = test_distance
                        outfile.write(segment)

            mean_depths[(chrom, amplicon)] = np.mean(amplicon_depth)

    return mean_depths


def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam:
        if not read.is_proper_pair:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def create_primer_lookup(ref_len_tuple, pools, amplicons: list[Amplicon], padding=35):
    """
    Returns a dict of chroms, each containing a (N, chrom_len) shaped array
    """
    lookups = {}
    for chrom, chromlen in ref_len_tuple:
        a = np.empty_like(None, shape=(1, chromlen + 1))
        for amp in amplicons:
            added = False
            if amp.chrom == chrom:
                # If amplicon clashes with any in same pool add new row
                amp_slice = a[
                    :,
                    max(amp.amplicon_start - padding, 0) : min(
                        amp.amplicon_end + padding, chromlen
                    ),
                ]
                for i, row in enumerate(amp_slice):  # Check each row for collision
                    if row[row != None].size == 0:
                        a[
                            i,
                            max(amp.amplicon_start - padding, 0) : min(
                                amp.amplicon_end + padding, chromlen
                            ),
                        ] = amp
                        added = True
                # If not added, add new row to array and add to that.
                if not added:
                    b = np.empty_like(None, shape=(1, chromlen + 1))
                    b[
                        0,
                        max(amp.amplicon_start - padding, 0) : min(
                            amp.amplicon_end + padding, chromlen
                        ),
                    ] = amp
                    a = np.vstack((a, b))

        lookups[chrom] = a
    return lookups


def go(args):
    """Filter and soft mask an alignment file so that the alignment boundaries match the primer start and end sites.

    Based on the most likely primer position, based on the alignment coordinates.
    """
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

    # Merge the primers
    scheme.bedlines = merge_primers(scheme.bedlines)

    amplicon_list = create_amplicons(scheme.bedlines)
    amplicons = {}
    for amplicon in amplicon_list:
        amplicon.length = amplicon.amplicon_end - amplicon.amplicon_start
        amplicons.setdefault(amplicon.chrom, {})[amplicon.amplicon_number] = amplicon

    pools = set([bl.pool for bl in scheme.bedlines])
    chroms = set([bl.chrom for bl in scheme.bedlines])

    pools_str = {str(x) for x in pools}
    pools_str.add("unmatched")

    # open the input samfile and process read groups
    if args.samfile and args.samfile != "-":
        infile = pysam.AlignmentFile(args.samfile, "rb")
    else:
        infile = pysam.AlignmentFile("-", "rb")

    first_segment = next(infile, None)
    if not first_segment:
        print("No segments found in the input file, exiting.", file=sys.stderr)
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
            print(
                "Output file path must end with either .bam or .sam, exiting.",
                file=sys.stderr,
            )
            sys.exit(1)

    else:
        outfile = pysam.AlignmentFile("-", "wh", header=bam_header)

    # Initialise the amplicon depth dict
    amp_depths = {}
    for amp in amplicon_list:
        amp_depths.setdefault(amp.chrom, {})
        amp_depths[amp.chrom].setdefault(
            amp.amplicon_number, np.zeros(amp.length, dtype=int)
        )

    # Initialise the mean depths dictionary, this will get stomped over if normalisation is requested
    mean_amp_depths = {}
    for chrom in amplicons:
        for amplicon in amplicons[chrom]:
            mean_amp_depths[(chrom, amplicon)] = 0

    # Create a lookup table for primer location
    ref_lengths = [(r, infile.get_reference_length(r)) for r in infile.references]
    primer_lookup = create_primer_lookup(
        ref_len_tuple=ref_lengths, pools=pools, amplicons=amplicon_list, padding=35
    )

    trimmed_segments = {x: {} for x in chroms}

    if paired:
        read_pairs = read_pair_generator(chained_iterator)

        for segments in read_pairs:
            if args.report:
                trimming_tuple = handle_segments(
                    segment=segments,  # type: ignore
                    lookup=primer_lookup,
                    args=args,
                    report_writer=report_writer,  # type: ignore
                    min_mapq=args.min_mapq,
                    outfile_writer=outfile,
                    amp_depths=amp_depths,
                )
            else:
                trimming_tuple = handle_segments(
                    segment=segments,  # type: ignore
                    lookup=primer_lookup,
                    args=args,
                    min_mapq=args.min_mapq,
                    outfile_writer=outfile,
                    amp_depths=amp_depths,
                )

            if not trimming_tuple:
                continue

            # unpack the trimming tuple since segment passed trimming
            amplicon, trimmed_pair = trimming_tuple

            # If we aren't normalising the segments will have already been written to the outfile
            if not args.normalise and not trimmed_pair:
                continue

            trimmed_segments[trimmed_pair[0].reference_name].setdefault(amplicon, [])  # type: ignore

            if trimmed_segments:
                trimmed_segments[trimmed_pair[0].reference_name][amplicon].append(  # type: ignore
                    trimmed_pair
                )

        # normalise if requested and write normalised segments to outfile
        if args.normalise:
            mean_amp_depths = normalise(
                trimmed_segments=trimmed_segments,
                normalise=args.normalise,
                primers=scheme.bedlines,
                outfile=outfile,
                verbose=args.verbose,
            )
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
                    fieldnames=["chrom", "amplicon", "mean_depth"],
                    delimiter="\t",
                )
                writer.writeheader()
                for (chrom, amplicon), depth in mean_amp_depths.items():
                    writer.writerow(
                        {"chrom": chrom, "amplicon": amplicon, "mean_depth": depth}
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
                )

            else:
                trimming_tuple = handle_segments(
                    segment=segment,
                    args=args,
                    min_mapq=args.min_mapq,
                    lookup=primer_lookup,
                    outfile_writer=outfile,
                    amp_depths=amp_depths,
                )

            if not trimming_tuple:
                continue

            # unpack the trimming tuple since segment passed trimming
            amplicon, trimmed_segment = trimming_tuple

            # If we aren't normalising the segments will have already been written to the outfile
            if not args.normalise and not trimmed_segment:
                continue

            trimmed_segments[trimmed_segment.reference_name].setdefault(amplicon, [])  # type: ignore

            if trimmed_segment and args.normalise:
                trimmed_segments[trimmed_segment.reference_name][amplicon].append(  # type: ignore
                    trimmed_segment
                )

        # normalise if requested
        if args.normalise:
            mean_amp_depths = normalise(
                trimmed_segments=trimmed_segments,
                normalise=args.normalise,
                primers=scheme.bedlines,
                outfile=outfile,
                verbose=args.verbose,
            )

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
                    fieldnames=["chrom", "amplicon", "mean_depth"],
                    delimiter="\t",
                )
                writer.writeheader()

                for (chrom, amplicon), depth in mean_amp_depths.items():
                    writer.writerow(
                        {"chrom": chrom, "amplicon": amplicon, "mean_depth": depth}
                    )

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
        help="Subsample to N coverage per amplicon. Use 0 for no normalisation. (default: %(default)s)",
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
        default=35,
        help="Fuzzy match primer positions within this threshold",
    )
    parser.add_argument(
        "--report", "-r", type=Path, help="Output report TSV to filepath"
    )
    parser.add_argument(
        "--amp-depth-report",
        "-a",
        type=Path,
        help="Output amplicon depth TSV to filepath",
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
    parser.add_argument("--verbose", "-v", action="store_true", help="Debug mode")
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
