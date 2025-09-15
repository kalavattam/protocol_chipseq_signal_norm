#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2024 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Script: parse_metadata_siq_chip.py
#
# Description:
#     Parses a metadata table in CSV or TSV format to extract siQ-ChIP
#     measurement values (e.g., volume, mass, concentration, length) based on
#     characteristics inferred from a BAM file (e.g., state, factor, strain).
#
# Usage:
#     python parse_metadata_siq_chip.py \
#         --text table.tsv --bam sample.bam [--shell] [--verbose]
#
# Arguments:
#      -v, --verbose  Enables verbose output for debugging or information.
#     -tm, --text     CSV or TSV siQ-ChIP metadata infile.
#     -eq, --eqn      Equation to compute. Options: '5', '5nd', '6', or '6nd'.
#      -b, --bam      BAM infile used to identify metadata row.
#     -sh, --shell    Outputs values in a shell-parseable format (for export).
#
# Example:
#     python parse_metadata_siq_chip.py \
#         --text "${HOME}/path/to/table.tsv" \
#         --bam "${HOME}/path/to/sample.bam" \
#         --shell
#
# Output:
#     When --shell is specified, outputs key-value pairs in shell-exportable
#     format. Otherwise, prints measurement values to standard output.
#
# License:
#     Distributed under terms of the MIT license.

#  Import libraries
import argparse
import csv
import os
import re
import sys


#  Run script in interactive mode (true) or command-line mode (false)
interactive = False

#  Define a dictionary of alternative names for each column
col_nam_map = {
    'vol_in': [
        'volume_in',
        'vol_in',
        'v_in',
        'input sample volume'
    ],
    'vol_all': [
        'volume_all',
        'volume_ip',
        'vol_all',
        'vol_ip',
        'v_all',
        'v_ip',
        'total volume before removal of input'
    ],
    'mass_in': [
        'mass_in',
        'm_in',
        'input DNA mass (ng)',
        'input DNA mass'
    ],
    'mass_ip': [
        'mass_ip',
        'm_ip',
        'IP DNA mass (ng)',
        'IP DNA mass'
    ],
    'conc_in': [
        'conc_in',
        'con_in',
        'concentration_in',
        'input DNA concentration (ng/µL)',
        'input DNA concentration'
    ],
    'conc_ip': [
        'conc_ip',
        'con_ip',
        'concentration_ip',
        'IP DNA concentration (ng/µL)',
        'IP DNA concentration'
    ],
    'len_in': [
        'length_in',
        'len_in',
        'input average fragment length (from TapeStation)',
        'input average fragment length (from Bioanalyzer)',
        'input average fragment length'
    ],
    'len_ip': [
        'length_ip',
        'len_ip',
        'IP average fragment length (from TapeStation)',
        'IP average fragment length (from Bioanalyzer)',
        'IP average fragment length'
    ],
    'dep_in': [
        'depth_in',
        'dep_in',
        'input sequencing depth',
        'input depth'
    ],
    'dep_ip': [
        'depth_ip',
        'dep_ip',
        'IP sequencing depth',
        'IP depth'
    ]
}


def set_interactive():
    """Set parameters for interactive mode."""
    verbose = True
    tbl_met = ...  # TODO
    eqn = '6nd'
    bam = ...  # TODO
    shell = True

    #  Return the arguments wrapped in argparse.Namespace
    return argparse.Namespace(
        verbose=verbose,
        tbl_met=tbl_met,
        eqn=eqn,
        bam=bam,
        shell=shell
    )


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# ###############################
# ##  Start: Work in progress  ##
# ###############################

# OLD
# def parse_bam_filename(filename):
#     """
#     Parse a BAM filename into its components based on the naming scheme.
# 
#     The expected filename format is as follows:
#         assay_genotype_state_treatment_factor_strain/replicate.
# 
#     Required filename components:
#         - assay: Must be 'IP' or 'in' and is always present. It must be
#                  followed by an underscore.
#         - factor: A required component preceded by an underscore.
#         - strain/replicate: A required component preceded by an underscore;
#                             it marks the end of the pattern.
# 
#     Optional filename components:
#         - genotype: If present, must be preceded by an underscore.
#         - state: An optional component with preferred values (e.g., 'G1',
#                  'G2M', 'log', or 'Q'), but can also be flexible. If present,
#                  it must be preceded by an underscore.
#         - treatment: An optional component preceded by an underscore if
#                      present.
# 
#     Args:
#         filename (str): Path to the BAM file.
# 
#     Returns:
#         dict: Parsed components as a dictionary with keys 'assay',
#               'genotype', 'state', 'treatment', 'factor', and 'strain'.
#     """
#     nam_bas = os.path.basename(filename)
#     bam_id = nam_bas.split('.')[0]  # Remove extension
# 
#     #  Define regex pattern to match the filename structure
#     pattern = (
#         r'^(?P<assay>IP|in)'                    # assay: mandatory (IP or in)
#         r'(?:_(?P<genotype>[^_]+))?'            # genotype: optional
#         r'(?:_(?P<state>G1|G2M|log|Q|[^_]+))?'  # state: optional (flexible)
#         r'(?:_(?P<treatment>[^_]+))?'           # treatment: optional
#         r'_(?P<factor>[^_]+)'                   # factor: mandatory
#         r'_(?P<strain>[^_]+)$'                  # strain/replicate: mandatory
#     )
# 
#     match = re.match(pattern, bam_id)
# 
#     if not match:
#         raise ValueError(
#             f"Filename '{filename}' does not match the expected pattern."
#         )
# 
#     #  Extract matched groups as a dictionary
#     comp = match.groupdict()
#     return comp

# "Parsed BAM filename components:\n",
# f"  - eqn={args.eqn}\n",
# f"  - assay={assay}\n",
# f"  - genotype={genotype}\n",
# f"  - state={state}\n",
# f"  - treatment={treatment}\n",
# f"  - factor={factor}\n",
# f"  - strain/replicate={strain}"

# NEW
def parse_bam_filename(filename, verbose=False):
    """
    Parse a BAM filename into its components.

    Expected tokens (underscored):
    assay _ [genotype] _ ... (state/factor in any order) ... [treatment?] _ strain
    """
    nam_bas = os.path.basename(filename)

    #  Strip BAM/CRAM/SAM and an extra trailing tag like ".sc"
    stem = nam_bas
    for ext in ('.bam', '.cram', '.sam'):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
            break
    stem = re.sub(r'\.[^.]+$', '', stem)

    toks = stem.split('_')
    if len(toks) < 3:
        raise ValueError(f"Filename '{filename}' is too short to parse.")

    assay = toks[0]
    if assay not in ('IP', 'in'):
        raise ValueError(
            f"Assay must be 'IP' or 'in' (got '{assay}') in '{filename}'."
        )

    strain = toks[-1]
    mid = toks[1:-1]  # everything between assay and strain

    # Known states; extend if needed.
    STATES = {'G1', 'G2M', 'log', 'Q'}

    def is_state(tok: str) -> bool:
        return tok in STATES

    # Accept common histone marks, "*-flag" proteins, or whitelisted proteins.
    def is_factor(tok: str) -> bool:
        return bool(re.match(
            r'^(?:'
            r'H[1-4](?:[A-Za-z0-9]+)?'  # H2B, H3, H4K20...
            r'|[A-Za-z0-9]+-flag'       # Isw1-flag, Swr1-flag, Htz1-flag
            r'|Hmo1|Hho1'
            r')$',
            tok
        ))

    # Genotype heuristics:
    # - 'WT' is a genotype
    # - token ending with '-KO' (case-insensitive) is a genotype
    # - explicitly NOT tokens ending with '-flag' (those are factors)
    def is_genotype(tok: str) -> bool:
        if tok == 'WT':
            return True
        if tok.lower().endswith('-ko'):
            return True
        if tok.lower().endswith('-flag'):
            return False
        return False

    genotype = 'N/A'
    state = 'N/A'
    factor = None
    treatment = 'N/A'

    if verbose:
        eprint(f"[parse] tokens={toks}  mid={mid}")

    # 1) Pull off genotype if present at the front of mid.
    cand = mid[:]
    if cand:
        if is_genotype(cand[0]):
            genotype = cand.pop(0)
            if verbose:
                eprint(f"[parse] genotype={genotype}, remaining={cand}")

    # 2) Find state and factor anywhere in remaining tokens (prefer first hits).
    state_idx = next((i for i, t in enumerate(cand) if is_state(t)), None)
    factor_idx = next((i for i, t in enumerate(cand) if is_factor(t)), None)

    # If ambiguous ordering caused wrong picks, try to favor a legit pair.
    if state_idx is not None and factor_idx is not None:
        # both found—good
        state = cand[state_idx]
        factor = cand[factor_idx]
        # 3) Anything left (other than those indices) could be treatment
        rest = [t for i, t in enumerate(cand) if i not in {state_idx, factor_idx}]
        if rest:
            treatment = rest[0]
    else:
        # Fallback heuristics similar to your original logic.
        if len(cand) == 1:
            t = cand[0]
            if is_state(t):
                state = t
            elif is_factor(t):
                factor = t
        elif len(cand) >= 2:
            a, b = cand[0], cand[1]
            rest = cand[2:]
            if is_state(a) and is_factor(b):
                state, factor = a, b
            elif is_factor(a) and is_state(b):
                factor, state = a, b
            else:
                if is_state(a):
                    state = a
                if is_state(b) and state == 'N/A':
                    state = b
                if is_factor(a) and factor is None:
                    factor = a
                if is_factor(b) and factor is None:
                    factor = b
            if rest:
                treatment = rest[0]

    if verbose:
        eprint(
            "[parse] resolved:",
            f"assay={assay} genotype={genotype} state={state} "
            f"factor={factor} treatment={treatment} strain={strain}"
        )

    if factor is None:
        raise ValueError(
            f"Could not determine factor from filename '{filename}'. "
            f"Parsed: assay={assay}, genotype={genotype}, state={state}, "
            f"strain={strain}"
        )

    return {
        'assay': assay,
        'genotype': genotype,
        'state': state,
        'treatment': treatment,
        'factor': factor,
        'strain': strain,
    }


# OLD
# def load_file(file_path):
#     """
#     Load CSV/TSV file as a list of dictionaries, standardizing column names.
#     """
#     if file_path.endswith('.csv'):
#         delm = ','
#     elif file_path.endswith(('.tsv', '.txt')):
#         delm = '\t'
#     else:
#         raise ValueError("Input file must be a CSV, TSV, or TXT file.")
# 
#     with open(file_path, newline='', encoding='utf-8') as file:
#         reader = csv.DictReader(file, delimiter=delm)
# 
#         # Standardize column names
#         standardized_columns = {
#             col: std_col for std_col, alt_nam in col_nam_map.items()
#             for col in reader.fieldnames if col in alt_nam
#         }
# 
#         # Rename columns and store as a list of dictionaries
#         data = []
#         for row in reader:
#             standardized_row = {
#                 standardized_columns.get(k, k): v for k, v in row.items()
#             }
#             data.append(standardized_row)
# 
#     return data


# NEW
def load_file(file_path, verbose=False):
    if file_path.endswith('.csv'):
        delm = ','
    elif file_path.endswith(('.tsv', '.txt')):
        delm = '\t'
    else:
        raise ValueError(
            "Input file must be a CSV, TSV, or tab-separated TXT file."
        )

    with open(file_path, newline='', encoding='utf-8') as file:
        reader = csv.DictReader(file, delimiter=delm)
        header = reader.fieldnames or []
        standardized_columns = {
            col: std_col
            for std_col, alt_nam in col_nam_map.items()
            for col in header if col in alt_nam
        }
        if verbose:
            eprint("[load] header:", header)
            eprint("[load] standardization map:", standardized_columns)

        data = []
        for row in reader:
            standardized_row = {
                standardized_columns.get(k, k): v for k, v in row.items()
            }
            data.append(standardized_row)

    if verbose:
        eprint(f"[load] rows loaded: {len(data)}")
    return data


# OLD
# def find_matching_row(data, factor, strain, genotype='N/A', state='N/A'):
#     """
#     Find a row in the data list matching factor, strain, and optional
#     genotype/state.
#     """
#     for row in data:
#         if (
#             row.get('factor', '').lower() == factor.lower()
#             and row.get('strain', '').lower() == strain.lower()
#             and (
#                 genotype == 'N/A' 
#                 or row.get('genotype', '').lower() == genotype.lower()
#             )
#             and (
#                 state == 'N/A' 
#                 or row.get('state', '').lower() == state.lower()
#             )
#         ):
#             return row
# 
#     raise ValueError(
#         f"No matching row found for state '{state}', factor '{factor}', "
#         f"strain '{strain}', genotype '{genotype}'."
#     )


# NEW
def find_matching_row(data, factor, strain, genotype='N/A', state='N/A',
                      verbose=False):
    """
    Find a metadata row by (factor or factor_abbrev), strain, and optional
    genotype/state.
    """
    def norm(s):
        return (s or "").strip().lower()

    nf, ns, ng, nst = map(norm, (factor, strain, genotype, state))
    if verbose:
        eprint(f"[match] want: factor={nf} strain={ns} "
               f"genotype={ng or 'n/a'} state={nst or 'n/a'}")

    def row_matches(row, fac, st):
        rfac = norm(row.get('factor'))
        rfab = norm(row.get('factor_abbrev'))
        rstr = norm(row.get('strain'))
        rgen = norm(row.get('genotype'))
        rstate = norm(row.get('state'))

        fac_ok = (fac == rfac) or (fac == rfab)
        gen_ok = (ng == 'n/a') or (ng == rgen)
        st_ok = (nst == 'n/a') or (st == rstate)

        if verbose and (rstr == ns):
            eprint(f"[match] row strain={rstr} fac=({rfac}|{rfab}) "
                   f"gen={rgen} state={rstate} → "
                   f"fac_ok={fac_ok} gen_ok={gen_ok} st_ok={st_ok}")
        return fac_ok and (ns == rstr) and gen_ok and st_ok

    # Pass 1: direct factor/state
    for row in data:
        if row_matches(row, nf, nst):
            if verbose:
                eprint("[match] PASS1 matched")
            return row

    # Pass 2: swap factor/state if filename order fooled upstream logic
    for row in data:
        if row_matches(row, nst, nf):
            if verbose:
                eprint("[match] PASS2 matched with swapped factor/state")
            return row

    raise ValueError(
        f"No matching row found for state '{state or 'N/A'}', "
        f"factor '{factor}', strain '{strain}', genotype '{genotype or 'N/A'}'."
    )


# ##############################
# ##  Stop: Work in progress  ##
# ##############################

def output_for_shell(**kwargs):
    """Print the extracted values in a shell-parseable format."""
    for key, value in kwargs.items():
        print(f"export {key}={value}")


def parse_args():
    """
    Parse command line arguments.

    Args:
        ... # TODO
    """
    parser = argparse.ArgumentParser(
        description=(
            "Parse BAM file and retrieve matching row from TSV/CSV table of "
            "siQ-ChIP metadata."
        )
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output."
    )
    parser.add_argument(
        "-tm", "--tbl_met",
        help="TSV or CSV table of siQ-ChIP metadata.",
        required=True
    )
    parser.add_argument(
        '-eq', '--eqn',
        type=str,
        required=True,
        choices=['5', '5nd', '6', '6nd'],
        default='6nd',
        help=(
            "Equation to compute the alpha scaling factor (PMID: 37160995; "
            "default: %(default)s). Options: '5' applies Equation 5 for use "
            "with fragment length-normalized coverage, '5nd' uses Equation 5 "
            "without depth terms for use with 'normalized coverage', '6' "
            "applies Equation 6 for use with fragment length-normalized "
            "coverage, and '6nd' uses Equation 6 without depth terms for use "
            "with 'normalized coverage'."
        )
    )
    parser.add_argument(
        "-b", "--bam",
        help="Input BAM file.",
        required=True
    )
    parser.add_argument(
        "-sh", "--shell",
        action="store_true",
        help="Output values in shell export format."
    )

    #  Display help and exit if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)
    
    return parser.parse_args()


def main():
    #  Use command line arguments or interactive setup
    if interactive:
        args = set_interactive()
    else:
        args = parse_args()

    #  Validate the alpha equation assignment
    if args.eqn not in {'5', '5nd', '6', '6nd'}:
        raise ValueError(f"Invalid equation selection: {args.eqn}")

    #  Parse the BAM filename
    comp = parse_bam_filename(args.bam, verbose=args.verbose)
    assay = comp['assay']
    genotype = comp.get('genotype', 'N/A')
    state = comp.get('state', 'N/A')
    treatment = comp.get('treatment', 'N/A')
    factor = comp['factor']
    strain = comp['strain']

    if args.verbose:
        eprint(
            "Parsed BAM filename components:\n"
            f"  - eqn={args.eqn}\n"
            f"  - assay={assay}\n"
            f"  - genotype={genotype}\n"
            f"  - state={state}\n"
            f"  - treatment={treatment}\n"
            f"  - factor={factor}\n"
            f"  - strain/replicate={strain}"
        )

    #  Load the TSV/CSV file
    data = load_file(args.tbl_met, verbose=args.verbose)

    #  Find the matching row
    try:
        row = find_matching_row(
            data, factor, strain, genotype, state, verbose=args.verbose
        )
        if args.verbose:
            eprint(f"Matching row found:\n{row}")
    except ValueError as e:
        eprint(e)
        sys.exit(1)

    #  Output based on --shell flag
    if args.shell:
        output_for_shell(
            eqn=args.eqn,
            vol_in=row['vol_in'],
            vol_all=row['vol_all'],
            mass_in=row['mass_in'],
            mass_ip=row['mass_ip'],
            conc_in=row.get('conc_in', 'N/A'),
            conc_ip=row.get('conc_ip', 'N/A'),
            len_in=row['len_in'],
            len_ip=row['len_ip'],
            dep_in=row.get('dep_in', 'N/A'),
            dep_ip=row.get('dep_ip', 'N/A')
        )
    else:
        print(f"vol_in: {row['vol_in']}")
        print(f"vol_all: {row['vol_all']}")
        print(f"mass_in: {row['mass_in']}")
        print(f"mass_ip: {row['mass_ip']}")
        print(f"conc_in: {row.get('conc_in', 'N/A')}")
        print(f"conc_ip: {row.get('conc_ip', 'N/A')}")
        print(f"len_in: {row['len_in']}")
        print(f"len_ip: {row['len_ip']}")
        print(f"dep_in: {row.get('dep_in', 'N/A')}")
        print(f"dep_ip: {row.get('dep_ip', 'N/A')}")
        # TODO: Change 'N/A' to 'NA' here and above, and in script using output


if __name__ == "__main__":
    main()
