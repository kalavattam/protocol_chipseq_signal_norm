#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2024 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Script: parse_metadata_siq_chip.py
#
# Description:
#     #TODO
#
# Usage:
#     python parse_metadata_siq_chip.py \
#         #TODO
#
# Arguments:
#     #TODO
#
# Example:
#     python parse_metadata_siq_chip.py \
#         #TODO
#
# Output:
#     #TODO
#
# License:
#     Distributed under terms of the MIT license.

#  Import libraries
import argparse
import os
import pandas as pd


#  Run script in interactive/test mode (True) or command-line mode (False)
interactive = False

#  Define a dictionary of alternative names for each column
column_name_map = {
    'volume_in': [
        'volume_in',
        'vol_in',
        'input sample volume'
    ],
    'volume_ip': [
        'volume_all',
        'volume_ip',
        'vol_all',
        'vol_ip',
        'total volume before removal of input'
    ],
    'mass_in': [
        'mass_in',
        'input DNA mass (ng)',
        'input DNA mass'
    ],
    'mass_ip': [
        'mass_ip',
        'IP DNA mass (ng)',
        'IP DNA mass'
    ],
    'conc_in': [
        'conc_in',
        'concentration_in',
        'input DNA concentration (ng/µL)',
        'input DNA concentration'
    ],
    'conc_ip': [
        'conc_ip',
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
        'input sequencing depth',
        'input depth'
    ],
    'depth_ip': [
        'depth_ip',
        'IP sequencing depth',
        'IP depth'
    ]
}


def standardize_columns(df):
    """Standardize dataframe column names to match expected names."""
    reverse_mapping = {
        alt: std for std, alts in column_name_map.items() for alt in alts
    }
    return df.rename(columns=lambda col: reverse_mapping.get(col, col))


def parse_bam_filename(filename):
    """Strip file extensions and split by '_'."""
    base_name = os.path.basename(filename)
    bam_id = base_name.split('.')[0]  # Remove extension
    parts = bam_id.split('_')

    if len(parts) < 4:
        raise ValueError(
            f"Filename '{filename}' does not follow the expected structure."
        )
    
    #  Extract relevant parts: exp, state, factor, strain
    exp, state, factor, strain = parts[0], parts[1], parts[2], parts[3]
    return exp, state, factor, strain


def load_file(file_path):
    """Load as CSV or TSV based on file extension."""
    if file_path.endswith('.csv'):
        delimiter = ','
    elif file_path.endswith('.tsv'):
        delimiter = '\t'
    else:
        raise ValueError("Input file must be a CSV or TSV file.")
    
    return pd.read_csv(file_path, delimiter=delimiter)


def find_matching_row(df, state, factor, strain):
    #  Ensure 'state', 'factor', and 'strain' columns are strings
    df['state'] = df['state'].astype(str)
    df['factor'] = df['factor'].astype(str)
    df['strain'] = df['strain'].astype(str)

    #  Find the matching row based on state, factor, and strain
    matched_row = df[
        (df['state'].str.lower() == state.lower()) &
        (df['factor'].str.lower() == factor.lower()) &
        (df['strain'].str.lower() == strain.lower())
    ]

    if matched_row.empty:
        raise ValueError(
            f"No matching row found for state '{state}', factor '{factor}', "
            f"strain '{strain}'."
        )

    #  Return the first match
    return matched_row.iloc[0]


def output_for_shell(**kwargs):
    """Print the extracted values in a shell-parseable format."""
    for key, value in kwargs.items():
        print(f"export {key}={value}")


def parse_args():
    """
    Parse command-line arguments or use hardcoded values for interactive mode.
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
        "--text",
        help="Input TEXT (TSV/CSV) file",
        required=not interactive
    )
    parser.add_argument(
        "--bam",
        help="Input BAM file",
        required=not interactive
    )
    parser.add_argument(
        "--shell",
        action="store_true",
        help="Output values in shell export format"
    )
    
    return parser.parse_args() if not interactive else argparse.Namespace(
        verbose=True,
        text="/home/kalavatt/projects-etc/202X_protocol_ChIP/data/measurements_siq_chip.tsv",
        # bam="/home/kalavatt/projects-etc/202X_protocol_ChIP/03_bam/bowtie2/global/flag-2_mapq-1/sc/IP_G1_Hho1_6336.sc.bam",
        bam="/home/kalavatt/tsukiyamalab/Kris/202X_protocol_ChIP/03_bam/bowtie2/global/flag-2_mapq-1/sc/IP_Q_Esa5_7041.sc.bam",
        shell=True
    )


def main():
    #  Parse arguments
    args = parse_args()

    #  Parse the BAM filename
    exp, state, factor, strain = parse_bam_filename(args.bam)
    if args.verbose:
        print(
            "Parsed BAM filename:\n",
            f"  - exp={exp}\n",
            f"  - state={state}\n",
            f"  - factor={factor}\n",
            f"  - strain={strain}"
        )

    #  Load the TSV/CSV file
    df = load_file(args.text)

    #  Standardize the column names
    df = standardize_columns(df)

    #  Find the matching row
    row = find_matching_row(df, state, factor, strain)
    if args.verbose:
        print(f"Matching row found:\n{row}")

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
