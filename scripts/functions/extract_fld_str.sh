#!/bin/bash

extract_fld_str() {
    local fld="${1}"
    local tbl="${2}"
    local show_help
    local num_fld

    show_help=$(cat << EOM
---------------
extract_fld_str
---------------

Description:
  Extract a specific column (field) from a tab-separated value (TSV) file
  (e.g., a 'table' or 'tbl') and return it as a single, comma-separated string.
  The function validates inputs and ensures the file exists, is readable, has
  more than just a header, and is properly formatted as a TSV file. It also
  checks that the specified field index is within the valid range of columns.

Positional parameters:
  1, fld (int): The 1-based index of the column to extract.
  2, tbl (str): Path to the TSV file.

Returns:
  0 and a comma-separated string containing the values from the specified 
  column. Otherwise, returns 1 with an appropriate error message.

Usage:
  extract_fld_str "\${fld}" "\${tbl}"

Examples:
  \`\`\`
  #  Extract the second column from a valid table
  tbl="example_table.tsv"
  extract_fld_str 2 "\${tbl}"
  "value1,value2,value3"

  #  Error: Invalid field index
  extract_fld_str 0 "\${tbl}"
  Error: Field index must be a positive integer.

  #  Error: File does not exist
  extract_fld_str 1 "nonexistent_file.tsv"
  Error: Table file does not exist, is not readable, or is empty: 'nonexistent_file.tsv'.

  #  Error: Field index out of range
  tbl="example_table.tsv"
  extract_fld_str 10 "\${tbl}"
  Error: Field index is out of range: '10'. The table has 3 fields.
  \`\`\`
EOM
    )

    #  Parse and check function parameter
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    #  Validate 'fld' is a positive integer
    if [[ ! "${fld}" =~ ^[1-9][0-9]*$ ]]; then
        echo "Error: Field index must be a positive integer." >&2
        return 1
    fi

    #  Validate 'tbl' exists, is readable, and is not empty
    if [[ ! -r "${tbl}" || ! -s "${tbl}" ]]; then
        echo \
            "Error: Table file does not exist, is not readable, or is empty:" \
            "'${tbl}'." >&2
        return 1
    fi

    #  Validate 'tbl' has more than one row (header + data)
    if (( $(awk 'END { print NR }' "${tbl}") <= 1 )); then
        echo \
            "Error: Table file does not have more than a header: '${tbl}'." >&2
        return 1
    fi

    #  Validate 'tbl' is tab-separated by checking the first data line
    if ! awk -F '\t' 'NR == 2 { exit (NF > 1 ? 0 : 1) }' "${tbl}"; then
        echo \
            "Error: Table file does not appear to be tab-separated:" \
            "'${tbl}'." >&2
        return 1
    fi

    #  Determine the number of fields in 'tbl' from the first data line
    num_fld=$(awk -F '\t' 'NR == 2 { print NF; exit }' "${tbl}")

    #  Check that 'fld' is within the valid range of fields
    if (( fld > num_fld )); then
        echo \
            "Error: Field index is out of range: '${fld}'. The table has" \
            "${num_fld} fields." >&2
        return 1
    fi

    #  Extract the specified field and return it as a comma-separated string
    awk \
        -v fld="${fld}" \
        'NR > 1 { val = val ? val "," $fld : $fld } END { print val }' \
        "${tbl}"
}
