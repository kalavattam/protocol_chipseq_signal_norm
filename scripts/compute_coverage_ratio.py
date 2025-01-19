#!/usr/bin/env python3

"""
Notes on normalization and log2 ratio handling in function
compute_coverage_ratio()

1. Allowing 'scl_fct' in log2 transformation
    - 'scl_fct' (scaling factor) is applied after log2 transformation:
      'scl_fct * log2(sig_ip / max(sig_in, dep_min) + 1)'
    - This is set up so that scaling is applied uniformly, maintaining
      consistency in how ratios are transformed.

2. Avoiding log2(0) errors using a pseudocount
    - A pseudocount ('+ 1') is added inside the log2 function to prevent
      undefined log2(0) cases:
      'log2(scl_fct * (sig_ip / max(sig_in, dep_min)) + 1)'
    - Without the pseudocount, bins with zero signal would cause log errors.

3. Handling missing bins with zero imputation
    - If a bin is missing in 'fil_in' but is present in 'fil_ip', the input
      value set to '0.0', thereby defaulting 'dep_min' and ensuring meaningful
      ratios.
    - If a bin is missing in 'fil_ip' but is present in 'fil_in', the IP value
      is set to '0.0'.
    - The calculation method ensures that missing bins do not cause division
      errors or discontinuities in the track.

4. Modularizing computation with 'compute_ratio()'
    - To eliminate redundant code, we introduced 'compute_ratio()', which
      applies the following:
        + Standard ratio computation
        + Log2 transformation with pseudocount of '1'
        + Zero imputation for missing bins

5. Output consistency
    - The output BEDGRAPH file will contain either:
        + A standard ratio, optionally scaled by '--scl_fct'.
        + A log2-transformed ratio, optionally scaled by '--scl_fct'.
    - The '--log2' flag determines whether log2 transformation is applied.
    - By default (when '--log2' is not used), the script computes the standard
      ratio.

6. Choosing the appropriate ratio calculation
    - To compute siQ- or spike-in scaled coverage:
        + Supply '--scl_fct' with the appropriate scaling factor obtained using
          the siQ-ChIP or spike-in method.
        + Do not use the '--log2' flag (i.e., omit it from the command).
        + This outputs scaled ratio coverage values.
    - To evaluate IP enrichment (log2 ratio of fil_ip/fil_in):
        + Leave '--scl_fct' at its default value (1.0).
        + Enable the '--log2' flag.
        + This outputs log2(sig_ip/sig_in) values, which is useful for
          visualizing enrichment.
"""

import argparse
import gzip
import math
import os
import sys


def open_file(pth_fil, mode="rt"):
    """
    Open a file normally or as a gzip depending on extension.
    """
    if pth_fil.endswith(".gz"):
        return gzip.open(pth_fil, mode)
    else:
        return open(pth_fil, mode)


# def get_chr_first(fil_ip):
#     """
#     Determine first chromosome in BEDGRAPH file.
#
#     Args:
#         fil_ip: Path to IP BEDGRAPH file (can be gzipped).
#
#     Returns:
#         str: Chromosome name found in the first line of the file.
#     """
#     with open_file(fil_ip) as f:
#         for line in f:
#             if line.strip():
#                 return line.split()[0]  # Extract the first field
#     raise ValueError("Error: The IP BEDGRAPH file appears to be empty.")


def check_bin_size(fil_ip, fil_in):
    """
    Verify that the bin sizes in IP and input BEDGRAPH files are identical.

    Args:
        fil_ip (str): Path to IP BEDGRAPH file.
        fil_in (str): Path to input BEDGRAPH file.

    Raises:
        ValueError: If the bin sizes differ between the two files.
    """
    with open_file(fil_ip) as ip_f, open_file(fil_in) as in_f:
        for _ in range(5):  # Check the first 5 bins
            lin_ip = ip_f.readline().strip()
            lin_in = in_f.readline().strip()

            if not lin_ip or not lin_in:
                continue

            fld_ip = lin_ip.split()
            fld_in = lin_in.split()

            # Extract bin sizes
            bin_ip = int(fld_ip[2]) - int(fld_ip[1])
            bin_in = int(fld_in[2]) - int(fld_in[1])

            if bin_ip != bin_in:
                raise ValueError(
                    f"Error: Mismatched bin sizes detected.\n"
                    f"  - {fil_ip}: {bin_ip}-bp bins\n"
                    f"  - {fil_in}: {bin_in}-bp bins\n"
                    f"Ensure both files are binned identically."
                )


def compute_ratio(sig_ip, sig_in, scl_fct, dep_min, log2):
    """
    Compute the ratio or log2 ratio of IP and input BEDGRAPH signals.

    Args:
        sig_ip  (flt): Signal value from IP BEDGRAPH file.
        sig_in  (flt): Signal value from input BEDGRAPH file.
        scl_fct (flt): Scaling factor applied to ratio or log2 ratio
                       computation.
        dep_min (flt): Minimum input depth to avoid extreme/erroneous division.
        log2    (bol): If 'True', compute 'scl_fct * log2(ratio + 1)';
                       otherwise, compute standard ratio.

    Returns:
        float: Computed signal ratio or log2-transformed ratio.
    """
    ratio = sig_ip / max(sig_in, dep_min)  # Prevent div-by-zero with 'dep_min'
    if log2:
        return scl_fct * math.log2(ratio + 1)
    else:
        return scl_fct * ratio


def compute_coverage_ratio(
    fil_ip, fil_in, fil_out, scl_fct, dep_min, rnd, log2
):
    """
    Compute the ratio or log2 ratio of IP and input BEDGRAPH tracks, applying
    a scaling factor. Handles cases where bins are missing in one of the files.

    Args:
        fil_ip  (str): Path to IP BEDGRAPH file (can be gzipped).
        fil_in  (str): Path to input BEDGRAPH file (can be gzipped).
        fil_out (str): Output file path (if '.gz' extension, gzip compression).
        scl_fct (flt): Scaling factor to apply to ratio, standard or log2.
        dep_min (flt): Minimum input depth to avoid extreme/erroneous division.
        rnd     (int): Decimal places for rounding binned signal values.
        log2    (bol): If 'True', compute 'scl_fct * log2(ratio + 1)';
                       otherwise, compute standard ratio.

    Writes an output BEDGRAPH file with column 4 as...
        - Default: 'scl_fct * (sig_ip / max(sig_in, dep_min))'
        - If log2: 'scl_fct * log2(sig_ip / max(sig_in, dep_min) + 1)'
    """
    #  Check bin size consistency before proceeding
    check_bin_size(fil_ip, fil_in)

    #  Determine if output should be gzipped
    gz_out = fil_out.endswith(".gz")

    #  Open input and output files safely
    with open_file(fil_ip) as ip_f, open_file(fil_in) as in_f, (
        gzip.open(fil_out, "wt") if gz_out else open(fil_out, "w")
    ) as f_out:
        
        lin_ip = ip_f.readline().strip()
        lin_in = in_f.readline().strip()

        while lin_ip or lin_in:  # Continue until both files are exhausted
            fld_ip = lin_ip.split() if lin_ip else None
            fld_in = lin_in.split() if lin_in else None
            
            #  Parse BEDGRAPH format, handling missing bins
            chr_ip, start_ip, end_ip, sig_ip = (
                (fld_ip[0], int(fld_ip[1]), int(fld_ip[2]), float(fld_ip[3]))
                if fld_ip else (None, None, None, 0.0)  # IP zero imputation
            )

            chr_in, start_in, end_in, sig_in = (
                (fld_in[0], int(fld_in[1]), int(fld_in[2]), float(fld_in[3]))
                if fld_in else (None, None, None, 0.0)  # Input zero imputation
            )

            #  Determine which bin to process
            if chr_ip == chr_in and start_ip == start_in:
                #  Both bins exist
                sig_mrg = compute_ratio(sig_ip, sig_in, scl_fct, dep_min, log2)
                f_out.write(
                    f"{chr_ip}\t{start_ip}\t{end_ip}\t{sig_mrg:.{rnd}f}\n"
                )

                #  Read next entries from both files
                lin_ip = ip_f.readline().strip()
                lin_in = in_f.readline().strip()

            elif chr_ip is None or (
                chr_in and (chr_ip > chr_in or start_ip > start_in)
            ):
                #  IP file is behind: Treat IP as 0 and move input forward
                sig_mrg = compute_ratio(0.0, sig_in, scl_fct, dep_min, log2)
                f_out.write(
                    f"{chr_in}\t{start_in}\t{end_in}\t{sig_mrg:.{rnd}f}\n"
                )

                lin_in = in_f.readline().strip()

            else:
                #  Input file is behind: Treat input as 0, thereby defaulting
                #  to 'dep_min', and move IP forward
                sig_mrg = compute_ratio(sig_ip, 0.0, scl_fct, dep_min, log2)
                f_out.write(
                    f"{chr_ip}\t{start_ip}\t{end_ip}\t{sig_mrg:.{rnd}f}\n"
                )

                lin_ip = ip_f.readline().strip()


# def compute_coverage_ratio(fil_ip, fil_in, fil_out, scl_fct, dep_min, rnd):
#     """
#     Compute the ratio of IP and input BEDGRAPH tracks, applying a scaling
#     factor.
#    
#     Args:
#         fil_ip:  Path to IP BEDGRAPH file (can be gzipped).
#         fil_in:  Path to input BEDGRAPH file (can be gzipped).
#         fil_out: Output file path (if '.gz' extension, gzip compression).
#         scl_fct: Scaling factor to apply to the ratio.
#         dep_min: Minimum expected input depth to avoid division by zero.
#         rnd:     Number of decimal places for rounding binned signal values.
#
#     Writes an output BEDGRAPH file with column 4 as
#     'scl_fct * (IP / max(sig_in, dep_min))' rounded to 'rnd' decimal places.
#     """
#
#     # Determine the first chromosome dynamically
#     first = get_chr_first(fil_ip)
#
#     #  Determine if output should be gzipped
#     gz_out = fil_out.endswith(".gz")
#
#     #  Open input and output files safely
#     with open_file(fil_ip) as ip_f, open_file(fil_in) as in_f, (
#         gzip.open(fil_out, "wt") if gz_out else open(fil_out, "w")
#     ) as f_out:
#
#         lin_ip = ip_f.readline().strip()
#         lin_in = in_f.readline().strip()
#
#         while lin_ip and lin_in:
#             fld_ip = lin_ip.split()
#             fld_in = lin_in.split()
#
#             #  Parse BEDGRAPH format
#             chr_ip, start_ip, end_ip, sig_ip = (
#                 fld_ip[0], int(fld_ip[1]), int(fld_ip[2]), float(fld_ip[3])
#             )
#
#             chr_in, start_in, end_in, sig_in = (
#                 fld_in[0], int(fld_in[1]), int(fld_in[2]), float(fld_in[3])
#             )
#
#             #  Ensure we are on the correct chromosome
#             if chr_ip != chr_in or chr_ip != first:
#                 lin_ip = ip_f.readline().strip() if chr_ip < chr_in else lin_ip
#                 lin_in = in_f.readline().strip() if chr_in < chr_ip else lin_in
#                 continue
#
#             #  Overlap condition: Intervals must intersect
#             if start_ip <= end_in and start_in <= end_ip:
#                 sig_mrg = scl_fct * (sig_ip / max(sig_in, dep_min))
#                 f_out.write(
#                     f"{chr_ip}\t{start_ip}\t{end_ip}\t{sig_mrg:.{rnd}f}\n"
#                 )
#
#                 #  Read next IP entry
#                 lin_ip = ip_f.readline().strip()
#             elif start_ip > end_in:
#                 #  Move input file forward
#                 lin_in = in_f.readline().strip()
#             else:
#                 #  Move IP file forward
#                 lin_ip = ip_f.readline().strip()


def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Merge IP and input signal tracks with siQ scaling.",
        argument_default=argparse.SUPPRESS
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        default=False,
        help="Increase output verbosity."
    )
    parser.add_argument(
        "-fp", "--fil_ip",
        required=True,
        type=str,
        help=(
            "Path to the IP (immunoprecipitation) BEDGRAPH file (can be "
            "gzipped)."
        )
    )
    parser.add_argument(
        "-fn", "--fil_in",
        required=True,
        type=str,
        help="Path to the input BEDGRAPH file (can be gzipped)."
    )
    parser.add_argument(
        "-fo", "--fil_out",
        required=True,
        type=str,
        help=(
            "Path to the output BEDGRAPH file (gzip compression applied if "
            "'.gz' extension is given)."
        )
    )
    parser.add_argument(
        "-sf", "--scl_fct",
        type=float,
        default=1.0,
        help=(
            "Multiplicative scaling factor applied to data ratio (default: "
            "%(default)s)."
        )
    )
    parser.add_argument(
        "-dm", "--dep_min",
        required=True,
        type=float,
        help=(
            "Minimum expected input depth (used to prevent extreme or "
            "erroneous division operations)."
        )
    )
    parser.add_argument(
        "-l2", "--log2",
        action="store_true",
        default=False,
        help=(
            "Compute log2(sig_ip/sig_in) with a pseudocount of 1 to prevent "
            "log(0) errors. The scaling factor ('--scl_fct') is applied after "
            "the log2 transformation: 'scl_fct * log2(sig_ip/sig_in + 1)'."
        )
    )
    parser.add_argument(
        "-r", "--rnd",
        type=int,
        default=24,
        help=(
            "Number of decimal places for rounding binned signal ratio values "
            "(default: %(default)s)."
        )
    )

    #  Display help and exit if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    return parser.parse_args()


def main():
    """Main script execution."""
    args = parse_args()

    #  Validate infile existence
    if not os.path.exists(args.fil_ip):
        raise FileNotFoundError(f"Error: IP file not found: {args.fil_ip}")
    if not os.path.exists(args.fil_in):
        raise FileNotFoundError(f"Error: Input file not found: {args.fil_in}")

    #  Validate output directory existence, writability
    out_dir = os.path.dirname(args.fil_out)
    if not os.path.isdir(out_dir):
        raise FileNotFoundError(
            f"Error: Output directory does not exist: {out_dir}"
        )
    if not os.access(out_dir, os.W_OK):
        raise PermissionError(
            f"Error: No write permission for output directory: {out_dir}"
        )

    #  Validate numeric arguments
    if args.scl_fct <= 0.0:
        raise ValueError("Error: Scaling factor (--scl_fct) must be > 0.")
    if args.dep_min <= 0.0:
        raise ValueError("Error: Expected depth (--dep_min) must be > 0.")
    if args.rnd <= 0:
        raise ValueError(
            "Error: Number of decimal places for rounding (--rnd) must be > 0."
        )

    #  Print verbose output
    if args.verbose:
        print("#############################################")
        print("## Arguments for compute_coverage_ratio.py ##")
        print("#############################################")
        print("")
        print(f"--fil_ip  {args.fil_ip}")
        print(f"--fil_in  {args.fil_in}")
        print(f"--fil_out {args.fil_out}")
        print(f"--scl_fct {args.scl_fct}")
        print(f"--dep_min {args.dep_min}")
        print(f"--log2    {args.log2}")
        print(f"--rnd     {args.rnd}")
        print("")
        print("")

    compute_coverage_ratio(
        args.fil_ip, args.fil_in, args.fil_out, args.scl_fct, args.dep_min,
        args.rnd, args.log2
    )


if __name__ == "__main__":
    main()
