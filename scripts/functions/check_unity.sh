#!/bin/bash

#  Helper function to validate that a value is a positive float
function val_pos_flt() {
    local val="${1}"  # Value assigned to range variable
    local nam="${2}"  # Name of range variable

    if [[
           ! "${val}" =~ ^[0-9]*\.?[0-9]+$ \
        || $(echo "${val} <= 0" | bc -l) -eq 1
    ]]; then
        echo "Error: '${nam}' must be a positive float: '${val}'." >&2
        return 1
    fi
}


function check_unity() {
    local fil_in="${1}"        # Input BEDGRAPH file (can be gzipped)
    local rng_gt="${2:-0.98}"  # Lower bound for unity check
    local rng_lt="${3:-1.02}"  # Upper bound for unity check
    local reader               # Determines whether to use 'cat' or 'zcat'
    local show_help            # Help message/documentation

    show_help=$(cat << EOM
-----------
check_unity
-----------

Description:
  This function checks whether the sum of column 4 in a BEDGRAPH file is 
  approximately 1 (i.e., unity). The input file can be plain text or gzipped.
  The acceptable range for summation can be set with optional parameters.

Positional parameters:
  1, fil_in (str): Path to the BEDGRAPH file (can be '.gz' compressed).
  2, rng_gt (flt): Lower bound for unity check (default: '${rng_gt}').
  3, rng_lt (flt): Upper bound for unity check (default: '${rng_lt}').

Returns:
  0 if the sum is within the defined range, 1 and error message otherwise.

Usage:
  check_unity fil_in [rng_gt] [rng_lt]

Examples:
  \`\`\`
  check_unity "coverage_track.bdg"
  check_unity "coverage_track.bdg.gz"
  check_unity "coverage_track.bdg.gz" 0.99 1.01
  \`\`\`
EOM
    )

    #  Display help message if no arguments or -h/--help is given
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    #  Validate input file
    if [[ ! -f "${fil_in}" ]]; then
        echo "Error: File not found or not accessible: '${fil_in}'." >&2
        return 1
    fi

    #  Validate range values are positive floats
    val_pos_flt "${rng_gt}" "rng_gt" || return 1
    val_pos_flt "${rng_lt}" "rng_lt" || return 1

    #  Determine whether to read file with 'gunzip -cd' or 'cat'
    if [[ "${fil_in}" == *.gz ]]; then
        reader="gunzip -cd"
    else
        reader="cat"
    fi

    #  Process the file and check its sum
    ${reader} "${fil_in}" \
        | awk -v rng_gt="${rng_gt}" -v rng_lt="${rng_lt}" '{
            sum += $4
        } END {
            if (sum >= rng_gt && sum <= rng_lt) {
                print "File sums to approximately unity:", sum
            } else {
                print "File does not sum to unity:", sum
            }
        }'
}
