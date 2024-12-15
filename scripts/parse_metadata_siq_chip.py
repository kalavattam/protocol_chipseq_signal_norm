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
#     -tx, --text     CSV or TSV siQ-ChIP metadata infile.
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
import os
import pandas as pd
import re
import sys


#  Run script in interactive/test mode (True) or command-line mode (False)
interactive = False

#  Define a dictionary of alternative names for each column
column_name_map = {
    'volume_in': [
        'volume_in',
        'vol_in',
        'v_in',
        'input sample volume'
    ],
    'volume_ip': [
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
    'length_in': [
        'length_in',
        'len_in',
        'input average fragment length (from TapeStation)',
        'input average fragment length (from Bioanalyzer)',
        'input average fragment length'
    ],
    'length_ip': [
        'length_ip',
        'len_ip',
        'IP average fragment length (from TapeStation)',
        'IP average fragment length (from Bioanalyzer)',
        'IP average fragment length'
    ],
    'depth_in': [
        'depth_in',
        'dep_in',
        'input sequencing depth',
        'input depth'
    ],
    'depth_ip': [
        'depth_ip',
        'dep_ip',
        'IP sequencing depth',
        'IP depth'
    ]
}


def set_interactive():
    """Set parameters for interactive mode."""
    verbose = True
    text = ...  # TODO
    bam = ...  # TODO
    shell = True

    #  Return the arguments wrapped in argparse.Namespace
    return argparse.Namespace(
        verbose=verbose,
        text=text,
        bam=bam,
        shell=shell
    )


def standardize_columns(df):
    """Standardize dataframe column names to match expected names."""
    reverse_mapping = {
        alt: std for std, alts in column_name_map.items() for alt in alts
    }
    return df.rename(columns=lambda col: reverse_mapping.get(col, col))


def parse_bam_filename(filename):
    """
    Parse a BAM filename into its components based on the naming scheme.

    The expected filename format is as follows:
        assay_genotype_state_treatment_factor_strain/replicate.
    
    Required filename components:
        - assay: Must be 'IP' or 'in' and is always present. It must be
                 followed by an underscore.
        - factor: A required component preceded by an underscore.
        - strain/replicate: A required component preceded by an underscore; it
                            marks the end of the pattern.

    Optional filename components:
        - genotype: If present, must be preceded by an underscore.
        - state: An optional component with preferred values (e.g., 'G1',
                 'G2M', 'log', or 'Q'), but can also be flexible. If present,
                 it must be preceded by an underscore.
        - treatment: An optional component preceded by an underscore if
                     present.

    Args:
        filename (str): Path to the BAM file.

    Returns:
        dict: Parsed components as a dictionary with keys 'assay', 'genotype',
              'state', 'treatment', 'factor', and 'strain'.
    """
    base_name = os.path.basename(filename)
    bam_id = base_name.split('.')[0]  # Remove extension

    #  Define regex pattern to match the filename structure
    pattern = (
        r'^(?P<assay>IP|in)'                    # assay: mandatory (IP or in)
        r'(?:_(?P<genotype>[^_]+))?'            # genotype: optional
        r'(?:_(?P<state>G1|G2M|log|Q|[^_]+))?'  # state: optional (flexible)
        r'(?:_(?P<treatment>[^_]+))?'           # treatment: optional
        r'_(?P<factor>[^_]+)'                   # factor: mandatory
        r'_(?P<strain>[^_]+)$'                  # strain/replicate: mandatory
    )

    match = re.match(pattern, bam_id)

    if not match:
        raise ValueError(
            f"Filename '{filename}' does not match the expected pattern."
        )

    #  Extract matched groups as a dictionary
    components = match.groupdict()
    return components


def load_file(file_path):
    """Load dataframe based on file extension: CSV, TSV, or TXT."""
    if file_path.endswith('.csv'):
        delimiter = ','
    elif file_path.endswith(('.tsv', '.txt')):
        delimiter = '\t'
    else:
        raise ValueError("Input file must be a CSV, TSV, or TXT file.")
    
    return pd.read_csv(file_path, delimiter=delimiter)


def find_matching_row(df, factor, strain, genotype='N/A', state='N/A'):
    #  Ensure 'state', 'factor', and 'strain' columns are strings
    df['genotype'] = df['genotype'].astype(str)
    df['state'] = df['state'].astype(str)
    df['factor'] = df['factor'].astype(str)
    df['strain'] = df['strain'].astype(str)

    #  Handle optional values (e.g., state and genotype can be 'N/A')
    conditions = (df['factor'].str.lower() == factor.lower()) & \
                 (df['strain'].str.lower() == strain.lower())

    if genotype != 'N/A':
        conditions &= (df['genotype'].str.lower() == genotype.lower())
    
    if state != 'N/A':
        conditions &= (df['state'].str.lower() == state.lower())

    matched_row = df[conditions]

    if matched_row.empty:
        raise ValueError(
            f"No matching row found for state '{state}', factor '{factor}', "
            f"strain '{strain}', genotype '{genotype}'."
        )

    #  Return the first match
    return matched_row.iloc[0]


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
        description="Parse BAM file and retrieve matching row from TSV/CSV."
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output."
    )
    parser.add_argument(
        "-tx", "--text",
        help="Input TEXT (TSV/CSV) file",
        required=True
    )
    parser.add_argument(
        "-b", "--bam",
        help="Input BAM file",
        required=True
    )
    parser.add_argument(
        "-sh", "--shell",
        action="store_true",
        help="Output values in shell export format"
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

    #  Parse the BAM filename
    parsed_components = parse_bam_filename(args.bam)
    assay = parsed_components['assay']
    genotype = parsed_components.get('genotype', 'N/A')
    state = parsed_components.get('state', 'N/A')
    treatment = parsed_components.get('treatment', 'N/A')
    factor = parsed_components['factor']
    strain = parsed_components['strain']

    if args.verbose:
        print(
            "Parsed BAM filename components:\n",
            f"  - assay={assay}\n",
            f"  - genotype={genotype}\n",
            f"  - state={state}\n",
            f"  - treatment={treatment}\n",
            f"  - factor={factor}\n",
            f"  - strain/replicate={strain}"
        )

    #  Load the TSV/CSV file
    df = load_file(args.text)

    #  Standardize the column names
    df = standardize_columns(df)

    #  Find the matching row, adjusting for missing components
    try:
        row = find_matching_row(df, factor, strain, genotype, state)
        if args.verbose:
            print(f"Matching row found:\n{row}")
    except ValueError as e:
        print(e, file=sys.stderr)
        sys.exit(1)

    #  Output based on --shell flag
    if args.shell:
        output_for_shell(
            volume_in=row['volume_in'],
            volume_ip=row['volume_ip'],
            mass_in=row['mass_in'],
            mass_ip=row['mass_ip'],
            conc_in=row.get('conc_in', 'N/A'),
            conc_ip=row.get('conc_ip', 'N/A'),
            length_in=row['length_in'],
            length_ip=row['length_ip'],
            depth_in=row.get('depth_in', 'N/A'),
            depth_ip=row.get('depth_ip', 'N/A')
        )
    else:
        print(f"volume_in: {row['volume_in']}")
        print(f"volume_ip: {row['volume_ip']}")
        print(f"mass_in: {row['mass_in']}")
        print(f"mass_ip: {row['mass_ip']}")
        print(f"conc_in: {row.get('conc_in', 'N/A')}")
        print(f"conc_ip: {row.get('conc_ip', 'N/A')}")
        print(f"length_in: {row['length_in']}")
        print(f"length_ip: {row['length_ip']}")
        print(f"depth_in: {row.get('depth_in', 'N/A')}")
        print(f"depth_ip: {row.get('depth_ip', 'N/A')}")


if __name__ == "__main__":
    main()
