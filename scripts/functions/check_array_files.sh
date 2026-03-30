#!/bin/bash

#  Check that all supplied file paths exist
function check_array_files() {
    local desc="${1:-}"
    local file

    #  Check that 'desc' input is not empty
    if [[ -z "${desc}" ]]; then
        echo "Error: Positional argument 1, 'desc', is required." >&2
        return 1
    fi

    shift

    #  Check that files were supplied
    if (( $# < 1 )); then
        echo "Error: No files supplied to validate for '${desc}'." >&2
        return 1
    fi
    
    #  Check each supplied file path
    for file in "$@"; do
        if [[ ! -f "${file}" ]]; then
            echo "Error: '${desc}' file does not exist: '${file}'." >&2
            return 1
        fi
    done

    return 0
}
