#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2024 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Script: compute_coverage.py
#
# Description:
#     compute_coverage.py calculates the binned coverage from a BAM file. It
#     can output coverage as either BEDGRAPH or BIGWIG files, or both. Users
#     can choose to normalize the coverage by fragment length and total reads,
#     or apply a scaling factor to adjust for, e.g., spike-in controls or to
#     model physical quantities of chromatin with a siQ-ChIP 'alpha' scaling
#     factor. The bin size and fragment length can be customized. The script
#     supports parallel processing, in which coverage calculations for each
#     chromosome are distributed across multiple CPU threads.
#
# Usage:
#     python compute_coverage.py \
#         [--verbose] --infile <str> --outfile <str> -outtype <option> \
#         --bin_siz <posint> [--scl_fct <flt>] [--norm] [--usr_frg <posint>]
#
# Arguments:
#      -v, --verbose  Run script in 'verbose mode'.
#      -t, --threads  Number of threads for parallel processing (default: 1)
#      -i, --infile   Path to the BAM infile.
#      -o, --outfile  Path and base name for the outfile (without extension).
#     -ot, --outtype  Output format: 'bedgraph', 'bigwig', or 'both' (default:
#                     bigwig).
#     -bs, --bin_siz  Bin size for coverage calculation in base pairs (default:
#                     30).
#     -sf, --scl_fct  Optional scaling factor to apply to the coverage (cannot
#                     be used with --norm).
#     -no, --norm     Normalize coverage by fragment length and total reads in
#                     the style of Dickson et al., Sci Rep 2023 (cannot be used
#                     with --scl_fct).
#     -uf, --usr_frg  Optional fixed fragment length to use instead of the BAM
#                     infile's query (read) or template lengths.
#
# Example:
#     ```
#     python compute_coverage.py \
#         --verbose \
#         --threads 4 \
#         --infile sample.bam \
#         --outfile sample_coverage \
#         --outttyp both \
#         --bin_siz 50 \
#         --norm
#     ```
#
# Output:
#     BEDGRAPH or BIGWIG file (or both) containing binned coverage data.
#
# License:
#     Distributed under terms of the MIT license.

#  Import libraries, etc.
import argparse
import gzip
import os
import pyBigWig
import pysam
import sys

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor


#  Run script in interactive/test mode (True) or command-line mode (False)
interactive = False

#  Define dictionary to map Roman numerals to integers for sorting
roman_to_int = {
    "I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6, "VII": 7, 
    "VIII": 8, "IX": 9, "X": 10, "XI": 11, "XII": 12, "XIII": 13, 
    "XIV": 14, "XV": 15, "XVI": 16
}


#  Define functions
def set_interactive():
    """Set up paths and parameters for interactive mode."""
    #  Set general paths
    dir_bas = "/home/kalavatt/tsukiyamalab/Kris"  # WARNING: Change as needed #
    dir_rep = f"{dir_bas}/202X_protocol_ChIP"
    dir_dat = f"{dir_rep}/data"
    dir_pro = f"{dir_dat}/processed"

    aln = "bowtie2"
    atp = "global"
    flg = 2
    mpq = 1
    det = f"flag-{flg}_mapq-{mpq}"
    cov = "spike"  # "alpha", "norm", "raw", or "spike"

    dir_aln = f"{dir_pro}/align_{aln}_{atp}/{det}"
    dir_bam = f"{dir_aln}/sc"
    fil_bam = "IP_WT_Q_Hho1_6337.sc.bam"  # WARNING: Change as needed #
    pth_bam = f"{dir_bam}/{fil_bam}"
    
    dir_cov = f"{dir_pro}/compute_coverage/{det}/{cov}"
    dir_trk = f"{dir_cov}/tracks"
    fil_trk = os.path.splitext(fil_bam)[0]
    pth_trk = f"{dir_trk}/{fil_trk}"
    
    #  Check if paths exist (optional for debugging/development purposes)
    assert os.path.exists(dir_bam), f"Indirectory {dir_bam} doesn't exist."
    assert os.path.exists(dir_trk), f"Outdirectory {dir_trk} doesn't exist."
    assert os.path.exists(pth_bam), f"Infile {os.path.basename(pth_bam)} doesn't exist."
    # assert os.path.exists(pth_trk), f"Outfile {os.path.basename(pth_trk)} doesn't exist."

    #  Set argument values
    verbose = True
    threads = 8
    infile = pth_bam
    outfile = pth_trk
    outtype = "bigwig"  # "both"
    bin_siz = 1
    scl_fct = 0.70371  # None
    norm = False
    usr_frg = None

    #  Return the arguments wrapped in argparse.Namespace
    return argparse.Namespace(
        verbose=verbose,
        threads=threads,
        infile=infile,
        outfile=outfile,
        outtype=outtype,
        bin_siz=bin_siz,
        scl_fct=scl_fct,
        norm=norm,
        usr_frg=usr_frg
    )


def sort_chrom_roman_arabic(chrom):
    """
    Map Roman numerals to integer values, otherwise return as is.
    """
    return roman_to_int.get(chrom, chrom)


def parse_bam(path_bam, usr_frg=None):
    """
    Reads a BAM file and returns a dictionary of lists, where each key is a
    chromosome (chrom), and each value is a list of tuples representing
    fragments. Each tuple contains: (start, end, len_frag), where 'start'
    and 'end' define the fragment's coordinates, and 'len_frag' is the
    fragment's length.

    For paired-end alignments, the function processes fragments based on
    properly paired reads, using FLAGs 99 or 1123 (forward strand) and
    163 or 1187 (reverse strand). For single-end alignments, the function
    estimates fragments using, by default, the read length and considers FLAGs
    0 or 1024 (forward strand) and 16 or 1040 (reverse strand). If 'usr_frg', a
    user-defined fragment length, is provided, then this value is used instead
    of 'template_length' for paired-end alignments or 'query_length' (i.e.,
    read length) for single-end alignments.
    """
    entries_bam = defaultdict(list)
    with pysam.AlignmentFile(path_bam, 'rb') as bam_file:
        for read in bam_file.fetch():
            #  Run paired-end logic: Work with properly paired alignments with
            #  FLAGs 99 or 1123, or 163 or 1187
            if read.flag in {99, 163, 1123, 1187}:
                chrom = bam_file.get_reference_name(read.reference_id)
                start = read.reference_start
                end = start + (usr_frg if usr_frg else read.template_length)
                len_frag = usr_frg if usr_frg else read.template_length
                entries_bam[chrom].append((start, end, len_frag))

            #  Run single-end logic: Work with alignments with FLAGs 0 or 1024
            #  (forward strand), or 16 or 1040 (reverse strand)
            elif read.flag in {0, 16, 1024, 1040}:
                chrom = bam_file.get_reference_name(read.reference_id)
                start = read.reference_start
                end = start + (usr_frg if usr_frg else read.query_length)
                len_frag = usr_frg if usr_frg else read.query_length
                entries_bam[chrom].append((start, end, len_frag))

    return entries_bam


def calculate_covg_chrom(
    chrom, entries_bam, frags_total, bin_siz, norm, scl_fct=None
):
    """
    Calculate the binned coverage for a set of BED entries, either 'normalized'
    or 'traditional'. Optionally apply a scaling factor.
    """
    #  Check for invalid bin size
    if bin_siz <= 0:
        raise ValueError("--bin_siz must be greater than 0.")

    #  Ensure scaling factor and normalization are not used together
    if scl_fct is not None and norm:
        raise ValueError("--scl_fct cannot be used together with --norm.")

    #  Proceed with coverage calculation
    coverage = defaultdict(float)
    for start, end, len_frag in entries_bam:
        contribution = 1 / len_frag if norm else 1
        for pos in range(start, end + 1):
            coverage[(chrom, pos)] += contribution

    #  Normalize by aligned reads total if specified
    if norm:
        for key in coverage:
            coverage[key] /= frags_total

    #  Apply scaling factor if provided
    if scl_fct is not None:
        if scl_fct <= 0:
            raise ValueError("--scl_fct must be greater than 0.")
        for key in coverage:
            coverage[key] *= scl_fct

    #  Bin the coverage data
    coverage_bin = defaultdict(float)
    for (chrom, pos), value in coverage.items():
        bin_strt = (pos // bin_siz) * bin_siz
        coverage_bin[(chrom, bin_strt)] += value / bin_siz

    return coverage_bin


def calculate_covg_task(data):
    return calculate_covg_chrom(*data)


def write_bigwig(coverage, outfile, bin_siz, chrom_siz):
    """
    Write the binned coverage data to a BIGWIG file.
    """
    #  Create a dictionary of chromosome sizes
    dict_size = {chrom: size for chrom, size in chrom_siz.items()}

    #  Open a BIGWIG file for writing
    bw = pyBigWig.open(outfile, 'w')

    #  Add the chromosome sizes to the BIGWIG header
    bw.addHeader(list(dict_size.items()))

    #  Initialize empty lists to store chromosomes, start positions, end
    #  positions, and coverage values
    chroms = []
    starts = []
    ends = []
    values = []

    #  Iterate over the coverage data, sorted by chromosome (with
    #  Roman-to-Arabic numeral conversion) and bin start position
    for (chrom, bin_start), value in sorted(
        coverage.items(), key=lambda x: (
            sort_chrom_roman_arabic(x[0][0]), x[0][1]
        )
    ):
        #  Append the current chromosome, start position, end position
        #  (start + bin size), and value to respective lists
        chroms.append(chrom)
        starts.append(bin_start)
        ends.append(bin_start + bin_siz)
        values.append(value)

    #  Write the chromosome, start, end, and value lists to the BIGWIG file
    bw.addEntries(chroms, starts, ends=ends, values=values)

    #  Close the BIGWIG file
    bw.close()


def write_bedgraph(coverage, outfile, bin_siz):
    """
    Write the binned coverage data to a BEDGRAPH file.
    """
    #  Open a gzip-compressed BEDGRAPH file for writing
    with gzip.open(outfile, 'wt') as outfile:
        #  Initialize lists to store chromosomes, start positions, end
        #  positions, and coverage values
        chroms = []
        starts = []
        ends = []
        values = []

        #  Iterate over the coverage data, sorted by chromosome (using
        #  Roman-to-Arabic numeral conversion) and bin start position
        for (chrom, bin_start), value in sorted(
            coverage.items(), key=lambda x: (
                sort_chrom_roman_arabic(x[0][0]), x[0][1]
            )
        ):
            #  Append the current chromosome, start position, end position
            #  (start + bin size), and value to respective lists
            chroms.append(chrom)
            starts.append(bin_start)
            ends.append(bin_start + bin_siz)
            values.append(value)

        #  Write the contents of chroms, starts, ends, and values to the
        #  BEDGRAPH file
        for chrom, start, end, value in zip(chroms, starts, ends, values):
            # Format each line as per BEDGRAPH requirements
            output_line = f"{chrom}\t{start}\t{end}\t{value:.24f}\n"
            outfile.write(output_line)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Calculate binned coverage, normalized or not, from a BAM file."
        ),
        argument_default=argparse.SUPPRESS
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Increase output verbosity for debugging."
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help=(
            "Number of threads for parallel processing (default: "
            "%(default)s). Each thread processes a different chromosome's "
            "coverage in parallel."
        )
    )
    parser.add_argument(
        "-i",
        "--infile",
        required=True,
        help="Path to the input BAM file."
    )
    parser.add_argument(
        "-o",
        "--outfile",
        required=True,
        help=(
            "Path to the output file without extension; extensions will be "
            "automatically determined by --outtype."
        )
    )
    parser.add_argument(
        "-ot",
        "--outtype",
        choices=['bedgraph', 'bigwig', 'both'],
        default='bigwig',
        help=(
            "Specify the output format: 'bedgraph', 'bigwig', or 'both' "
            "(default: %(default)s)."
        )
    )
    parser.add_argument(
        "-bs",
        "--bin_siz",
        type=int,
        default=10,
        help=(
            "Bin size for coverage calculation in base pairs (default: "
            "%(default)s)."
        )
    )
    parser.add_argument(
        "-sf",
        "--scl_fct",
        type=float,
        default=None,
        help=(
            "Optional scaling factor to apply to the coverage (default: "
            "%(default)s). Cannot be used with --norm."
        )
    )
    parser.add_argument(
        "-no",
        "--norm",
        action="store_true",
        help=(
            "Normalize coverage by fragment length and total reads, "
            "generating 'normalized coverage' as described in Dickson et "
            "al., Sci Rep 2023. If not specified, then 'traditional coverage' "
            "(unadjusted or scaling factor-adjusted aligned read depth) is "
            "calculated. Cannot be used with --scl_fct."
        )
    )
    parser.add_argument(
        "-uf",
        "--usr_frg",
        type=int,
        default=None,
        help=(
            "Optional user-specified fixed fragment length to use instead of "
            "read lengths (single-end alignments) or template lengths "
            "(paired-end alignments; default: %(default)s)."
        )
    )

    return parser.parse_args()


def main():
    #  Use command-line arguments or interactive setup based on `interactive`
    if interactive:
        args = set_interactive()
    else:
        args = parse_args()

    #  Check if --outfile has a disallowed extension
    filename, ext = os.path.splitext(args.outfile)
    disallowed = {"bigwig", "bw", "bedgraph", "bg"}
    if ext.lower().lstrip('.') in disallowed:
        print(
            f"Error: --outfile '{args.outfile}' has a disallowed extension "
            f"'{ext}'. Provide only the base name (e.g., "
            f"'{os.path.basename(filename)}') without a disallowed extension.",
            file=sys.stderr
        )
        sys.exit(1)

    #  Print argument assignments in 'verbose mode'
    if args.verbose:
        print(f"--threads  {args.threads}")
        print(f"--infile   {args.infile}")
        print(f"--outfile  {args.outfile}")
        print(f"--outtype  {args.outtype}")
        print(f"--bin_siz  {args.bin_siz}")
        print(f"--norm     {args.norm}")
        print(f"--scl_fct  {args.scl_fct if args.scl_fct else None}")
        print(f"--usr_frg  {args.usr_frg if args.usr_frg else None}")

    #  Parse and process BAM file
    entries_bam = parse_bam(args.infile, args.usr_frg)
    frags_total = sum(len(entries) for entries in entries_bam.values())
    covg_comb = defaultdict(float)

    #  Collect chromosome sizes from BAM file
    with pysam.AlignmentFile(args.infile, 'rb') as bam_file:
        chrom_siz = {seq['SN']: seq['LN'] for seq in bam_file.header['SQ']}

    #  Prepare data for parallel processing
    task_data = [
        (chrom, entries, frags_total, args.bin_siz, args.norm, args.scl_fct)
        for chrom, entries in entries_bam.items()
    ]

    #  Execute parallel processing
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        results = executor.map(calculate_covg_task, task_data)

    #  Combine results from all processes
    for result in results:
        for key, value in result.items():
            covg_comb[key] += value

    #  Write outfile based on user-specified format
    if args.outtype in ['bedgraph', 'both']:
        #  Write output to BEDGRAPH file
        write_bedgraph(covg_comb, f"{args.outfile}.bdg.gz", args.bin_siz)

    if args.outtype in ['bigwig', 'both']:
        #  Write output to BIGWIG file
        write_bigwig(covg_comb, f"{args.outfile}.bw", args.bin_siz, chrom_siz)


if __name__ == "__main__":
    main()
