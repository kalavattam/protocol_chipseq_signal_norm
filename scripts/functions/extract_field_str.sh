#!/bin/bash

function extract_field_str() {
    local tbl="${1}"         # Path to the TSV file
    local fld="${2}"         # 1-based index of the column to extract
    local hdr="${3:-false}"  # Skip header (true/false)
    local show_help
    local num_fld

    show_help=$(cat << EOM
-----------------
extract_field_str
-----------------

Description:
  Extract a specific column (field) from a tab-separated value (TSV) file
  (e.g., a 'table' or 'tbl') and return it as a single, comma-separated string.
  The function validates inputs and ensures the file exists, is readable, has
  more than just a header, and is properly formatted as a TSV file. It also
  checks that the specified field index is within the valid range of columns.

Positional parameters:
  1, tbl (str): Path to the TSV file.
  2, fld (int): The 1-based index of the column to extract.
  3, hdr (bol): Optional: 'true' to skip header (default: 'false').

Returns:
  0 and a comma-separated string containing the values from the specified 
  column. Otherwise, returns 1 with an appropriate error message.

Usage:
  extract_field_str "\${tbl}" "\${fld}" ["\${hdr}"]

Examples:
  \`\`\`
  #  Extract the second column from a valid table
  tbl="example_table.tsv"
  extract_field_str "\${tbl}" 2
  "value1,value2,value3"

  #  Error: Invalid field index
  extract_field_str "\${tbl}" 0
  Error: Field index must be a positive integer.

  #  Error: File does not exist
  extract_field_str "nonexistent_file.tsv" 1
  Error: Table file does not exist, is not readable, or is empty: 'nonexistent_file.tsv'.

  #  Error: Field index out of range
  tbl="example_table.tsv"
  extract_field_str "\${tbl}" 10
  Error: Field index is out of range: '10'. The table has 3 fields.
  \`\`\`
EOM
    )

    #  Parse and check function parameter
    if [[ -z "${tbl}" || "${tbl}" == "-h" || "${tbl}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    if [[ -z "${fld}" ]]; then
        echo "Error: Missing input for 'fld', positional parameter 2." >&2
        echo "" >&2
        echo "${show_help}"
        return 1
    fi

    #  Validate 'tbl' exists, is readable, and is not empty
    if [[ ! -r "${tbl}" || ! -s "${tbl}" ]]; then
        echo \
            "Error: Table file does not exist, is not readable, or is empty:" \
            "'${tbl}'." >&2
        return 1
    fi

    #  Validate 'fld' is a positive integer
    if [[ ! "${fld}" =~ ^[1-9][0-9]*$ ]]; then
        echo "Error: Field index must be a positive integer." >&2
        return 1
    fi

    #  Validate 'hdr' is properly formatted Boolean
    case "${hdr}" in
        t|true)  hdr=true  ;;
        f|false) hdr=false ;;
        *)
            echo \
                "Error: Input for 'hdr', positional parameter 3, is" \
                "'${hdr}'. Expected 't', 'true', 'f', or 'false'." >&2
            return 1
            ;;
    esac

    #  Validate that the table has enough data rows: 1 or more lines, or header
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

    #  Validate 'tbl' is tab-separated by checking the first data line
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
