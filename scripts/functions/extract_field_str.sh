#!/bin/bash
# -*- coding: utf-8 -*-
# 
# Script: extract_field_str.sh
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5-series models) was used in development.
# 
# Distributed under the MIT license.

function extract_field_str() {
    local tbl="${1:-}"       # Path to the TSV file
    local fld="${2:-}"       # 1-based index of the column to extract
    local hdr="${3:-false}"  # Skip header (true/false)
    local hdr_lc             # Lowercase-converted skip-header Boolean
    local num_fld            # No. tab-delimited fields in first data line
    local show_help          # Help message

    show_help=$(cat << EOM
Usage:
  extract_field_str [-h|--hlp|--help] tbl fld [hdr]

Description:
  Extract a specific field (column) from a tab-separated value (TSV) file (e.g., a 'table' or 'tbl') and return it as a single, comma-separated string.

  Validates inputs, ensuring file exists, is readable, has more than just a header, and is properly formatted as a TSV file.

  Checks that the specified field index is within the valid range of columns.

Positional arguments:
  1  tbl  <str>  Path to the TSV file.
  2  fld  <int>  The 1-based index of the column to extract.
  3  hdr  <bol>  'true' to skip header, 'false' to include it (default: ${hdr}).

Returns:
  0 and a comma-separated string containing the values from the specified column; otherwise, 1 with an error message.

Examples:
  1. Extract the second column from a valid table
  '''bash
  tbl="example_table.tsv"
  extract_field_str "\${tbl}" 2
  '''

  '''txt
  "value1,value2,value3"
  '''

  2. Confirm invalid field index errors
  '''bash
  extract_field_str "\${tbl}" 0
  '''

  '''txt
  Error: Field index must be a positive integer.
  '''

  3. Confirm nonexistent file errors
  '''bash
  extract_field_str "nonexistent_file.tsv" 1
  '''

  '''txt
  Error: Table file does not exist: 'nonexistent_file.tsv'.
  '''

  4. Confirm out-of-range field index errors
  '''bash
  tbl="example_table.tsv"
  extract_field_str "\${tbl}" 10
  '''

  '''txt
  Error: Field index is out of range: '10'. The table has 3 fields.
  '''
EOM
    )

    #  Parse and check function arguments
    if [[ "${tbl}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${tbl}" ]]; then
        echo "Error: Positional argument 1, 'tbl', is missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if [[ -z "${fld}" ]]; then
        echo "Error: Positional argument 2, 'fld', is missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Check that 'tbl' exists, is readable, and is not empty
    if [[ ! -e "${tbl}" ]]; then
        echo "Error: Table file does not exist: '${tbl}'." >&2
        return 1
    elif [[ ! -r "${tbl}" ]]; then
        echo "Error: Table file is not readable: '${tbl}'." >&2
        return 1
    elif [[ ! -s "${tbl}" ]]; then
        echo "Error: Table file is empty: '${tbl}'." >&2
        return 1
    fi

    #  Validate 'fld' is a positive integer
    if [[ ! "${fld}" =~ ^[1-9][0-9]*$ ]]; then
        echo "Error: Field index must be a positive integer." >&2
        return 1
    fi

    #  Convert 'hdr' to lowercase letters so Boolean matching accepts values
    #+ like 'T' / 'tRuE' / 'False' / 'FALSE' / etc.
    hdr_lc=$(printf '%s' "${hdr}" | tr '[:upper:]' '[:lower:]')

    #  Check that 'hdr' is a properly formatted Boolean
    case "${hdr_lc}" in
        t|true)  hdr=true  ;;
        f|false) hdr=false ;;
        *)
            echo \
                "Error: Input for 'hdr', positional argument 3, is" \
                "'${hdr}'. Expected 't', 'true', 'f', or 'false'." >&2
            return 1
            ;;
    esac

    #  Check that the table has enough data rows: 1 or more lines, or header
    #+ and one or more lines
    if ! \
        awk \
            -v hdr="${hdr}" \
            'END {
                min_rows = (hdr == "true" ? 2 : 1)
                if (NR < min_rows) exit 1
            }' \
                "${tbl}"
    then
        echo \
            "Error: Table file does not contain enough data rows:" \
            "'${tbl}'." >&2
        return 1
    fi

    #  Check that 'tbl' is tab-separated by checking the first data line
    if ! \
        awk \
            -F '\t' \
            -v hdr="${hdr}" \
            'NR == (hdr == "true" ? 2 : 1) { exit (NF > 1 ? 0 : 1) }' \
                "${tbl}"
    then
        echo \
            "Error: Table file does not appear to be tab-separated:" \
            "'${tbl}'." >&2
        return 1
    fi

    #  Determine the number of fields in 'tbl' from the first data line
    num_fld="$(
        awk \
            -F '\t' \
            -v hdr="${hdr}" \
            'NR == (hdr == "true" ? 2 : 1) { print NF; exit }' \
                "${tbl}"
    )"

    #  Check that 'fld' is within the valid range of fields
    if (( fld > num_fld )); then
        echo \
            "Error: Field index is out of range: '${fld}'. The table has" \
            "${num_fld} fields." >&2
        return 1
    fi

    #  Extract the specified field and return it as a comma-separated string
    awk \
        -F '\t' \
        -v fld="${fld}" \
        -v hdr="${hdr}" \
        '(hdr == "true" && NR == 1) {
            next
        } {
            val = val ? val "," $fld : $fld
        } END {
            print val
        }' \
            "${tbl}"
}
