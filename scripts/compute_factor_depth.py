#!/usr/bin/env python3

import argparse
import pysam
import sys


def count_alignments(pth_bam):
    """
    Counts alignments in a BAM file based on specific FLAGs.

    - For paired-end reads, counts alignments with flags 99, 1123, 163, 1187.
    - For single-end reads, counts alignments with flags 0, 16, 1024, 1040.

    Args:
        pth_bam (str): Path to BAM file.

    Returns:
        int: Total count of valid alignments (both paired- and single-end).
    """
    # Define FLAG sets
    flg_pe = {99, 1123, 163, 1187}  # Proper paired-end alignments
    flg_se = {0, 16, 1024, 1040}    # Single-end alignments

    # Initialize counter
    n_in = 0

    try:
        with pysam.AlignmentFile(pth_bam, 'rb') as fil_bam:
            for read in fil_bam.fetch():
                if read.flag in flg_pe or read.flag in flg_se:
                    n_in += 1
    except FileNotFoundError:
        print(
            f"Error: BAM file not found: '{pth_bam}'.",
            file=sys.stderr
        )
        sys.exit(1)
    except ValueError as e:
        print(
            f"Error processing BAM file '{pth_bam}': {e}",
            file=sys.stderr
        )
        sys.exit(1)
    except Exception as e:
        print(
            f"Unexpected error with BAM file '{pth_bam}': {e}",
            file=sys.stderr
        )
        sys.exit(1)

    return n_in


def compute_factor_depth(fil_bam, siz_bin, siz_gen):
    """
    Compute depth factor for input normalization.

    Args:
        fil_bam (str): Path to input BAM file.
        siz_bin  (int): Size of bins in BEDGRAPH.
        siz_gen  (int): Effective genome size of model organism.

    Returns:
        float: Computed depth factor.
    """
    #  Count valid alignments with 'count_alignments'
    n_in = count_alignments(fil_bam)

    #  Compute depth factor
    fct_dep = ((n_in * siz_bin) / siz_gen) / (1 - (siz_bin / siz_gen))

    return fct_dep


def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Compute depth factor for input normalization.",
        argument_default=argparse.SUPPRESS
    )
    parser.add_argument(
        "-b",
        "--fil_bam",
        required=True,
        type=str,
        help="Path to input BAM file."
    )
    parser.add_argument(
        "-sb",
        "--siz_bin",
        required=True,
        type=int,
        help="Bin size in base pairs (BEDGRAPH bin size)."
    )
    parser.add_argument(
        "-sg",
        "--siz_gen",
        required=True,
        type=int,
        help="Effective genome size of the model organism."
    )
    parser.add_argument(
        '-r', '--rnd',
        type=int,
        default=20,
        help=(
            'Number of decimal places for rounding depth factor (default: '
            '%(default)s.'
        )
    )

    return parser.parse_args()


def main():
    """Main script execution."""
    args = parse_args()

    #  Compute depth factor
    fct_dep = compute_factor_depth(args.fil_bam, args.siz_bin, args.siz_gen)

    #  Print result
    print(f"{fct_dep:.{args.rnd}f}")


if __name__ == "__main__":
    main()
