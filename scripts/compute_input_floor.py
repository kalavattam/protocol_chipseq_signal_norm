#!/usr/bin/env python3

import argparse
import gzip
import os
import pysam
import sys


def count_alignments_bam(pth_bam):
    """
    Counts alignments in a BAM file based on specific FLAGs:
        - For paired-end alignments, counts alignments with flags 99, 1123,
          163, or 1187.
        - For single-end alignments, counts alignments with flags 0, 16, 1024,
          or 1040.

    Args:
        pth_bam (str): Path to BAM file.

    Returns:
        int: Total count of valid alignments (both paired- and single-end).
    """
    #  Define FLAG sets
    flg_pe = {99, 1123, 163, 1187}  # Proper paired-end alignments
    flg_se = {0, 16, 1024, 1040}    # Single-end alignments

    n_in = 0

    try:
        with pysam.AlignmentFile(pth_bam, 'rb') as fil_bam:
            for read in fil_bam.fetch():
                if read.flag in flg_pe or read.flag in flg_se:
                    n_in += 1
    except (FileNotFoundError, OSError, ValueError) as e:
        print(
            f"Error: Cannot process BAM file '{pth_bam}': {e}",
            file=sys.stderr
        )
        sys.exit(1)

    return n_in


def count_alignments_bed(pth_bed):
    """
    Count the number of alignments (lines) in a BED file.

    Args:
        pth_bed (str): Path to BED file (can be gzipped).

    Returns:
        int: Total count of alignments (lines).
    """
    try:
        #  Handle gzip-compressed files
        open_func = gzip.open if pth_bed.endswith(".gz") else open
        with open_func(pth_bed, "rt") as f:
            return sum(1 for _ in f)
    except (FileNotFoundError, PermissionError) as e:
        print(
            f"Error: Cannot open BED file '{pth_bed}': {e}",
            file=sys.stderr
        )
        sys.exit(1)
    except Exception as e:
        print(
            f"Unexpected error with BED file '{pth_bed}': {e}",
            file=sys.stderr
        )
        sys.exit(1)


def compute_factor(infile, siz_bin, siz_gen, mode):
    """
    Compute depth factor for input normalization.

    Args:
        infile  (str): Path to input BAM or BED file.
        siz_bin (int): Size of bins in BEDGRAPH.
        siz_gen (int): Effective genome size of model organism.
        mode    (str): 'frag' for fragment-length normalized signal, 'norm' for
                       fragment length- and unity-normalized signal.

    Returns:
        float: Computed depth factor.
    """
    if mode not in {"frag", "norm"}:
        print(
            f"Error: Invalid mode '{mode}'. Use 'frag' or 'norm'.",
            file=sys.stderr
        )
        sys.exit(1)

    #  Determine count of valid alignments
    if infile.endswith(".bam"):
        n_in = count_alignments_bam(infile)
    elif infile.endswith(".bed"):
        n_in = count_alignments_bed(infile)
    else:
        print(
            f"Error: Unsupported file type: {infile}. Provide a BAM or BED "
            "file.",
            file=sys.stderr
        )
        sys.exit(1)

    #  Compute depth factor
    fct_dep = ((n_in * siz_bin) / siz_gen) / (1 - (siz_bin / siz_gen))

    #  Prevent division by zero in 'norm' mode
    if mode == "norm":
        if n_in == 0:
            print(
                "Error: No valid alignments found. Cannot compute depth "
                "factor for 'norm' mode.",
                file=sys.stderr
            )
            sys.exit(1)
        return fct_dep / n_in  # Apply unity normalization

    return fct_dep


def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Compute depth factor used to avoid extreme or erroneous"
            "divisions in input normalization."
        ),
        argument_default=argparse.SUPPRESS
    )
    parser.add_argument(
        "-i", "--infile",
        type=str,
        required=True,
        help="Path to input BAM or BED file."
    )
    parser.add_argument(
        "-sb", "--siz_bin",
        type=int,
        default=10,
        help="BEDGRAPH bin size in base pairs (default: %(default)s)."
    )
    parser.add_argument(
        "-sg", "--siz_gen",
        type=int,
        default=12157105,
        help=(
            "Effective genome size of the model organism (default: "
            "%(default)s)."
        )
    )
    parser.add_argument(
        "-m", "--mode",
        choices=["frag", "norm"],
        default="norm",
        help=(
            "Normalization mode (default: '%(default)s'). Use 'frag' to "
            "compute the input depth factor for fragment-length normalized "
            "signal (intermediate data in the siQ-ChIP workflow). Use 'norm' "
            "to compute the input depth factor for 'normalized' signal, which "
            "is both fragment length- and unity-normalized."
        )
    )
    parser.add_argument(
        "-r", "--rnd",
        default=24,
        type=int,
        help=(
            "Number of decimal places for rounding depth factor (default: "
            "%(default)s)."
        )
    )

    #  Display help and exit if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    return parser.parse_args()


def main():
    """Main script execution."""
    args = parse_args()

    #  Validate input file(s) exist(s) before processing
    if not os.path.isfile(args.infile):
        print(
            f"Error: Input file not found: {args.infile}",
            file=sys.stderr
        )
        sys.exit(1)

    #  Check that bin and genome-size integers are positive
    if args.siz_bin <= 0:
        print(
            "Error: siz_bin must be greater than zero.",
            file=sys.stderr
        )
        sys.exit(1)
    if args.siz_gen <= 0:
        print(
            "Error: siz_gen must be greater than zero.",
            file=sys.stderr
        )
        sys.exit(1)

    #  Compute depth factor
    fct_dep = compute_factor(
        args.infile, args.siz_bin, args.siz_gen, args.mode
    )

    #  Print result
    print(f"{fct_dep:.{args.rnd}f}")


if __name__ == "__main__":
    main()
