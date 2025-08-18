#!/bin/bash

#  Function to populate array with 'NA' if it is empty
function populate_array_empty() {
    local arr_nam="${1}"
    local target="${2}"

    #  Use indirect reference to access the array
    eval "local arr_siz=\${#${arr_nam}[@]}"

    # shellcheck disable=SC2154
    if [[ "${arr_siz}" -eq 0 ]]; then
        for ((i = 0; i < target; i++)); do
            eval "${arr_nam}+=( \"NA\" )"
        done
    fi
}
