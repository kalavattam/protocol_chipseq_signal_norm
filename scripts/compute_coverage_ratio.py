#!/usr/bin/env python3

import argparse
import gzip
import os


def open_file(pth_fil, mode="rt"):
    """
    Open a file normally or as a gzip depending on extension.
    """
    if pth_fil.endswith(".gz"):
        return gzip.open(pth_fil, mode)
    else:
        return open(pth_fil, mode)


def get_chr_first(fil_ip):
    """
    Determine first chromosome in BEDGRAPH file.

    Arguments:
        - fil_ip: Path to IP BEDGRAPH file (can be gzipped).

    Returns:
        - str: Chromosome name found in the first line of the file.
    """
    with open_file(fil_ip) as f:
        for line in f:
            if line.strip():
                return line.split()[0]  # Extract the first field
    raise ValueError("Error: The IP BEDGRAPH file appears to be empty.")


def compute_ratio(fil_ip, fil_in, fil_out, fctr, dep_exp, rnd):
    """
    Compute the ratio of IP and input BEDGRAPH tracks, applying a scaling
    factor.
    
    Arguments:
        - fil_ip:  Path to IP BEDGRAPH file (can be gzipped).
        - fil_in:  Path to input BEDGRAPH file (can be gzipped).
        - fil_out: Output file path (if '.gz' extension, gzip compression).
        - fctr:    Scaling factor to apply to the ratio.
        - dep_exp: Minimum expected depth to avoid division by zero.
        - rnd:     Number of decimal places for rounding binned signal values.

    Writes an output BEDGRAPH file with column 4 as 'fctr * (IP / max(sig_in,
    dep_exp))' rounded to 'rnd' decimal places.
    """

    # Determine the first chromosome dynamically
    first = get_chr_first(fil_ip)

    #  Determine if output should be gzipped
    gz_out = fil_out.endswith(".gz")

    #  Open input and output files safely
    with open_file(fil_ip) as ip_f, open_file(fil_in) as in_f, \
         (gzip.open(fil_out, "wt") if gz_out else open(fil_out, "w")) as f_out:
        
        lin_ip = ip_f.readline().strip()
        lin_in = in_f.readline().strip()

        while lin_ip and lin_in:
            fld_ip = lin_ip.split()
            fld_in = lin_in.split()
            
            #  Parse BEDGRAPH format
            chr_ip, start_ip, end_ip, sig_ip = (
                fld_ip[0], int(fld_ip[1]), int(fld_ip[2]), float(fld_ip[3])
            )

            chr_in, start_in, end_in, sig_in = (
                fld_in[0], int(fld_in[1]), int(fld_in[2]), float(fld_in[3])
            )
            
            #  Ensure we are on the correct chromosome
            if chr_ip != chr_in or chr_ip != first:
                lin_ip = ip_f.readline().strip() if chr_ip < chr_in else lin_ip
                lin_in = in_f.readline().strip() if chr_in < chr_ip else lin_in
                continue

            #  Overlap condition: Intervals must intersect
            if start_ip <= end_in and start_in <= end_ip:
                sig_mrg = fctr * (sig_ip / max(sig_in, dep_exp))
                f_out.write(
                    f"{chr_ip}\t{start_ip}\t{end_ip}\t{sig_mrg:.{rnd}f}\n"
                )
                
                #  Read next IP entry
                lin_ip = ip_f.readline().strip()
            elif start_ip > end_in:
                #  Move input file forward
                lin_in = in_f.readline().strip()
            else:
                #  Move IP file forward
                lin_ip = ip_f.readline().strip()


def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Merge IP and input signal tracks with siQ scaling.",
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
        "-fp",
        "--fil_ip",
        required=True,
        type=str,
        help=(
            "Path to the IP (immunoprecipitation) BEDGRAPH file (can be "
            "gzipped)."
        )
    )
    parser.add_argument(
        "-fn",
        "--fil_in",
        required=True,
        type=str,
        help="Path to the input BEDGRAPH file (can be gzipped)."
    )
    parser.add_argument(
        "-fo",
        "--fil_out",
        required=True,
        type=str,
        help=(
            "Path to the output BEDGRAPH file (gzip compression applied if "
            "'.gz' extension is given)."
        )
    )
    parser.add_argument(
        "-sf",
        "--scl_fct",
        type=float,
        default=1.0,
        help=(
            "Scaling factor multiplicatively applied to data ratio (default: "
            "%(default)s)."
        )
    )
    parser.add_argument(
        "-de",
        "--dep_exp",
        required=True,
        type=float,
        help="Expected input depth to prevent division errors."
    )
    parser.add_argument(
        '-r', '--rnd',
        type=int,
        default=20,
        help='Number of decimal places for rounding binned signal values.'
    )

    return parser.parse_args()


def main():
    """Main script execution."""
    args = parse_arguments()

    #  Validate file existence
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
    if args.dep_exp <= 0.0:
        raise ValueError("Error: Expected depth (--dep_exp) must be > 0.")

    #  Print verbose output
    if args.verbose:
        print("#############################################")
        print("## Arguments for compute_coverage_ratio.py ##")
        print("#############################################")
        print("")
        print(f"--fil_ip  {args.fil_ip}")
        print(f"--fil_in  {args.fil_in}")
        print(f"--fil_out {args.fil_out}")
        print(f"--scl_fct {args.fctr}")
        print(f"--dep_exp {args.dep_exp}")
        print(f"--rnd     {args.rnd}")
        print("")
        print("")

    compute_ratio(
        args.fil_ip, args.fil_in, args.fil_out, args.scl_fct, args.dep_exp,
        args.rnd
    )


if __name__ == "__main__":
    main()
