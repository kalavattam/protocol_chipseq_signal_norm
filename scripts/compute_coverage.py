#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2024-2025 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Script: compute_coverage.py
#
# Description:
#     'compute_coverage.py' calculates the binned coverage from a BAM file. It
#     can output coverage as either BEDGRAPH or BIGWIG files, or both. Users
#     can choose the type of coverage calculation ('raw', 'frag', or
#     'normalized'), optionally applying a scaling factor for modeling
#     purposes. The bin size and fragment length can be customized. The script
#     supports parallel processing, in which coverage calculations for each
#     chromosome are distributed across multiple CPU threads (although this
#     needs to be further optimized).
#
# Usage:
#     python compute_coverage.py \
#         [--verbose] --infile <str> --outfile <str> --typ_out <option> \
#         --siz_bin <posint> [--scl_fct <flt>] [--typ_cvg <option>]
#         [--usr_frg <posint>]
#
# Arguments:
#      -v, --verbose  Run script in 'verbose mode'.
#      -t, --threads  Number of threads for parallel processing (default: 1)
#      -i, --infile   Path to the BAM infile.
#      -o, --outfile  Path and base name for the outfile without extension.
#     -to, --typ_out  Output format: 'bedgraph', 'bdg', 'bg', 'bigwig', 'bw',
#                     or 'both' (default: 'bdg'). Also, the user can supply
#                     'bed' to write a BED file of alignments processed by
#                     function parse_bam(). Note: When '--typ_out bed' is used,
#                     the script will output alignment data processed by
#                     parse_bam() in BED-like format; this does not calculate
#                     coverage.
#     -sb, --siz_bin  Bin size for coverage calculation in base pairs (default:
#                     10).
#     -tv, --typ_cvg  Specify coverage calculation type (default: 'norm').
#                     Options:
#                       - 'raw', 'unadj', 'unadjusted': Compute unadjusted
#                         ("raw") coverage.
#                       - 'frag', 'len_frag': Normalize coverage by fragment
#                         length.
#                       - 'norm', 'normalized': Normalize coverage by fragment
#                         length and total fragments such that coverage sums to
#                         unity.
#     -sf, --scl_fct  Optional scaling factor to apply to the coverage. Can be
#                     used with any '--typ_cvg' value.
#     -uf, --usr_frg  Optional fixed fragment length to use instead of the BAM
#                     infile's query (read) or template (fragment) lengths.
#
# Example:
#     ```
#     python compute_coverage.py \
#         --verbose \
#         --threads 4 \
#         --infile sample.bam \
#         --outfile sample_coverage \
#         --typ_out bdg \
#         --siz_bin 5 \
#         --typ_cvg norm
#     ```
#
# Output:
#     BEDGRAPH or BIGWIG file (or both) containing binned coverage data.
#
# License:
#     Distributed under the MIT license.

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
    """Set paths and parameters for interactive mode."""
    #  Set general paths
    dir_bas = "/home/kalavatt/tsukiyamalab/Kris"  # WARNING: Change as needed #
    dir_rep = f"{dir_bas}/protocol_chipseq_signal_norm"
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
    assert os.path.exists(dir_bam), f"No indir: {dir_bam}."
    assert os.path.exists(dir_trk), f"No outdir: {dir_trk}."
    assert os.path.exists(pth_bam), f"No infile: {os.path.basename(pth_bam)}."

    #  Set argument values
    verbose = True
    threads = 8
    infile = pth_bam
    outfile = pth_trk
    typ_out = "bigwig"  # "both"
    siz_bin = 1
    scl_fct = 0.70371  # None
    norm = False
    usr_frg = None

    #  Return the arguments wrapped in argparse.Namespace
    return argparse.Namespace(
        verbose=verbose,
        threads=threads,
        infile=infile,
        outfile=outfile,
        typ_out=typ_out,
        siz_bin=siz_bin,
        scl_fct=scl_fct,
        norm=norm,
        usr_frg=usr_frg
    )


def sort_chrom_roman_arabic(chrom):
    """
    Map Roman numerals to integer values, otherwise return as is.
    """
    return roman_to_int.get(chrom, chrom)


def validate_outfile(value):
    """
    Validates that --outfile does not contain a disallowed extension.
    """
    disallowed = {"bedgraph", "bdg", "bg", "bigwig", "bw"}
    filename, ext = os.path.splitext(value)
    if ext.lower().lstrip('.') in disallowed:
        raise argparse.ArgumentTypeError(
            f"Invalid extension '{ext}' in --outfile. Provide only the base "
            f"name (e.g., '{os.path.basename(filename)}')."
        )
    return value

    
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
    try:
        with pysam.AlignmentFile(path_bam, 'rb') as bam_file:
            for read in bam_file.fetch():
                #  Run paired-end logic: Work with properly paired alignments
                #  with FLAGs 99 or 1123, or 163 or 1187
                if read.flag in {99, 163, 1123, 1187}:
                    chrom = bam_file.get_reference_name(read.reference_id)
                    start = read.reference_start
                    end = start + (
                        usr_frg if usr_frg else read.template_length
                    )
                    len_frag = usr_frg if usr_frg else read.template_length
                    entries_bam[chrom].append((start, end, len_frag))

                #  Run single-end logic: Work with alignments with FLAGs 0 or
                #  1024 (forward strand), or 16 or 1040 (reverse strand)
                elif read.flag in {0, 16, 1024, 1040}:
                    chrom = bam_file.get_reference_name(read.reference_id)
                    start = read.reference_start
                    end = start + (usr_frg if usr_frg else read.query_length)
                    len_frag = usr_frg if usr_frg else read.query_length
                    entries_bam[chrom].append((start, end, len_frag))
    except FileNotFoundError:
        print(
            f"Error: BAM file '{path_bam}' not found.",
            file=sys.stderr
        )
        sys.exit(1)
    except ValueError as e:
        print(
            f"Error processing BAM file '{path_bam}': {e}",
            file=sys.stderr
        )
        sys.exit(1)
    except Exception as e:
        print(
            f"Unexpected error with BAM file '{path_bam}': {e}",
            file=sys.stderr
        )
        sys.exit(1)

    return entries_bam


def calculate_covg_chrom(
    chrom, entries_bam, frags_total, siz_bin, is_norm, is_len, scl_fct=None
):
    """
    Calculate the binned coverage for a set of BED entries, with options to
    normalize by fragment length and/or total fragments. Optionally applies a
    scaling factor.
    """
    #  Check for invalid bin size
    if siz_bin <= 0:
        raise ValueError("'--siz_bin' must be greater than 0.")

    #  Proceed with coverage calculation
    coverage = defaultdict(float)
    for start, end, len_frag in entries_bam:
        #  Normalize by fragment length if specified
        contribution = (1 / len_frag) if (is_norm or is_len) else 1

        #  Add the contribution across the span of the fragment
        for pos in range(start, end + 1):
            coverage[(chrom, pos)] += contribution

    #  Normalize by total fragments if 'is_norm' is true
    if is_norm:
        for key in coverage:
            coverage[key] /= frags_total

    #  Apply scaling factor if specified
    if scl_fct is not None:
        if scl_fct <= 0:
            raise ValueError("'--scl_fct' must be greater than 0.")
        for key in coverage:
            coverage[key] *= scl_fct

    #  Bin the coverage data
    coverage_bin = defaultdict(float)
    for (chrom, pos), value in coverage.items():
        bin_strt = (pos // siz_bin) * siz_bin
        coverage_bin[(chrom, bin_strt)] += value

    return coverage_bin


def calculate_covg_task(data):
    return calculate_covg_chrom(*data)


def write_bigwig(coverage, outfile, siz_bin, chrom_siz):
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
        ends.append(bin_start + siz_bin)
        values.append(value)

    #  Write the chromosome, start, end, and value lists to the BIGWIG file
    bw.addEntries(chroms, starts, ends=ends, values=values)

    #  Close the BIGWIG file
    bw.close()


# def write_chrom_bigwig(data):
#     """
#     # TODO: Description of function.
#     """
#     chrom, chrom_coverage, siz_bin, outfile = data
#     with pyBigWig.open(outfile, 'w') as bw:
#         chrom_coverage.sort(key=lambda x: x[0])  # Sort by start position
#         starts, ends, values = zip(*chrom_coverage)
#         bw.addEntries(
#             [chrom] * len(starts), starts, ends=ends, values=values
#         )


# def write_bigwig_parallel(coverage, outfile, siz_bin, chrom_siz, threads):
#     """
#     Parallelized BIGWIG writing by chromosome.
#     """
#     #  Prepare data for each chromosome
#     tasks = []
#     for chrom in sorted(chrom_siz.keys(), key=sort_chrom_roman_arabic):
#         chrom_coverage = [
#             (start, start + siz_bin, coverage[(chrom, start)])
#             for (chrom_name, start) in coverage.keys()
#             if chrom_name == chrom
#         ]
#         tasks.append(
#             (chrom, chrom_coverage, siz_bin, f"{outfile}.{chrom}.tmp.bw")
#         )
#
#     #  Process chromosomes in parallel
#     with ProcessPoolExecutor(max_workers=threads) as executor:
#         executor.map(write_chrom_bigwig, tasks)
#
#     #  Merge temporary BIGWIG files into one (requires pyBigWig or similar
#     #  tool)
#     with pyBigWig.open(outfile, 'w') as bw:
#         bw.addHeader(list(chrom_siz.items()))
#         for task in tasks:
#             tmp_file = task[3]
#             with pyBigWig.open(tmp_file, 'r') as tmp_bw:
#                 for chrom in tmp_bw.chroms():
#                     for entry in tmp_bw.entries(chrom):
#                         bw.addEntries(
#                             [chrom], entry[0], ends=entry[1], values=entry[2]
#                         )
#             os.remove(tmp_file)


def write_bedgraph(coverage, outfile, siz_bin):
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
            ends.append(bin_start + siz_bin)
            values.append(value)

        #  Write the contents of chroms, starts, ends, and values to the
        #  BEDGRAPH file
        for chrom, start, end, value in zip(chroms, starts, ends, values):
            # Format each line as per BEDGRAPH requirements
            output_line = f"{chrom}\t{start}\t{end}\t{value:.24f}\n"
            outfile.write(output_line)


def parse_args():
    """
    Parse command line arguments.

    Args:
         -v, --verbose (flag): Increase output verbosity.
         -t, --threads  (int): Number of threads for parallel processing.
         -i, --infile   (str): Path to input BAM file.
         -o, --outfile  (str): Path to output file without extension.
        -to, --typ_out  (str): Specify output format. Options: 'bedgraph',
                               'bdg', 'bg', 'bigwig', 'bw', 'both'.
        -sb, --siz_bin  (int): Bin size for coverage calculation in base pairs.
        -tc, --typ_cvg  (str): Specify coverage calculation type. Options:
                               'raw', 'unadj', 'unadjusted', 'frag',
                               'len_frag', 'norm', 'normalized'.
        -sf, --scl_fct  (flt): Optional scaling factor to apply to coverage.
        -uf, --usr_frg  (int): Optional fixed fragment length.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Calculate binned coverage from a BAM file with optional "
            "normalization. Alternatively, output processed alignments in a "
            "BED-like format using the parse_bam() function."
        ),
        argument_default=argparse.SUPPRESS
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Increase output verbosity."
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
        type=validate_outfile,
        required=True,
        help=(
            "Path to the output file without extension; extensions will be "
            "automatically determined by '--typ_out'."
        )
    )
    parser.add_argument(
        "-to", "--typ_out",
        choices=['bedgraph', 'bdg', 'bg', 'bigwig', 'bw', 'both', 'bed'],
        default='bdg',
        help=(
            "Specify the output format. Options: 'bedgraph', 'bdg', 'bg', "
            "'bigwig', 'bw', 'both', 'bed'. Use 'bedgraph', 'bdg', 'bg' to "
            "output binned coverage in BEDGRAPH format. Use 'bigwig', 'bw' to "
            "output binned coverage in BIGWIG format. Use 'both' to output "
            "binned coverage in both BEDGRAPH and BIGWIG formats. Use 'bed' "
            "to output alignments processed by function parse_bam() in "
            "BED-like format (note: these are not coverage values; setting "
            "'--typ_out bed' results in arguments such as '--siz_bin', "
            "'--typ_cvg', '--scl_fct' '--usr_frg' are ignored). Default: "
            "%(default)s."
        )
    )
    parser.add_argument(
        "-sb", "--siz_bin",
        type=int,
        default=10,
        help=(
            "Bin size for coverage calculation in base pairs (default: "
            "%(default)s)."
        )
    )
    parser.add_argument(
        "-tv", "--typ_cvg",
        choices=[
            'raw', 'unadj', 'unadjusted', 'frag', 'len_frag', 'norm',
            'normalized'
        ],
        default='norm',
        help=(
            "Specify coverage calculation type. Options: 'raw', 'unadj', "
            "'unadjusted', 'frag', 'len_frag', 'norm', 'normalized'. Use "
            "'raw', 'unadj', 'unadjusted' to calculate unadjusted coverage. "
            "Use 'frag', 'len_frag' to normalize coverage by fragment length. "
            "Use 'norm', 'normalized' to compute 'normalized coverage' (PMID: "
            "37160995), i.e., coverage normalized for both fragment length "
            "and total fragments such that coverage sums to unity. Default: "
            "%(default)s."
        )
    )
    parser.add_argument(
        "-sf",
        "--scl_fct",
        type=float,
        default=None,
        help=(
            "Optional scaling factor to apply to the coverage (default: "
            "%(default)s)."
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

    #  Display help and exit if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    return parser.parse_args()


def main():
    #  Use command line arguments or interactive setup based on 'interactive'
    args = set_interactive() if interactive else parse_args()

    #  Determine whether to normalize based on 'typ_cvg'
    is_len = args.typ_cvg in {'frag', 'len_frag'}
    is_norm = args.typ_cvg in {'norm', 'normalized'}

    #  Print verbose output
    if args.verbose:
        print("## Arguments for compute_coverage.py ##")
        for arg, value in vars(args).items():
            print(f"--{arg:<6}: {value}")

    #  Parse and process BAM file
    entries_bam = parse_bam(args.infile, args.usr_frg)
    if args.typ_out == 'bed':
        #  Write alignments processed by parse_bam() in a BED-like format
        with gzip.open(f"{args.outfile}.bed.gz", "wt") as bed_file:
            for chrom, fragments in entries_bam.items():
                for start, end, length in fragments:
                    bed_file.write(f"{chrom}\t{start}\t{end}\t{length}\n")
        return  # Omit coverage calculations, etc.

    frags_total = sum(len(entries) for entries in entries_bam.values())
    if is_norm and frags_total == 0:
        raise ValueError(
            "Normalization requires non-zero total fragments. Check the BAM "
            "file."
        )

    #  Collect chromosome sizes from BAM file
    with pysam.AlignmentFile(args.infile, 'rb') as bam_file:
        chrom_siz = {seq['SN']: seq['LN'] for seq in bam_file.header['SQ']}

    #  Prepare and execute parallel tasks
    task_data = [
        (
            chrom, entries, frags_total, args.siz_bin, is_norm, is_len,
            args.scl_fct
        )
        for chrom, entries in entries_bam.items()
    ]

    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        results = executor.map(calculate_covg_task, task_data)

    #  Combine results from all processes
    covg_comb = defaultdict(float)
    for result in results:
        for key, value in result.items():
            covg_comb[key] += value

    #  Write outfile based on user-specified format
    if args.typ_out in ['bedgraph', 'bdg', 'bg', 'both']:
        #  Write output to BEDGRAPH file
        write_bedgraph(covg_comb, f"{args.outfile}.bdg.gz", args.siz_bin)

    if args.typ_out in ['bigwig', 'bw', 'both']:
        #  Write output to BIGWIG file
        write_bigwig(covg_comb, f"{args.outfile}.bw", args.siz_bin, chrom_siz)

    # if args.typ_out in ['bigwig', 'both']:
    #     #  Write output to BIGWIG file using the parallel implementation
    #     write_bigwig_parallel(
    #         covg_comb, f"{args.outfile}.bw", args.siz_bin, chrom_siz,
    #         args.threads
    #     )


if __name__ == "__main__":
    main()
