#!/bin/bash

#  Function to validate existence of files in arrays
function check_array_files() {
    local desc="${1}"
    shift

    #  Check that files are supplied
    if [[ "$#" -eq 0 ]]; then
        echo "Error: No files supplied to validate for ${desc}." >&2
        return 1
    fi
    
    #  Iterate through the supplied files and check their existence
    for file in "$@"; do
        if [[ ! -f "${file}" ]]; then
            echo "Error: ${desc} file does not exist: '${file}'." >&2
            return 1
        fi
    done

    #  If all files exist, return success
    return 0
}
