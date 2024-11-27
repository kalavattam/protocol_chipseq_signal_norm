#!/bin/bash

#  Function to validate that two arrays have matching lengths
function check_arrays_lengths() {
    local arr_nam_1="${1}"
    local -n arr_1="${2}"
    local arr_nam_2="${3}"
    local -n arr_2="${4}"

    if [[ "${#arr_1[@]}" -ne "${#arr_2[@]}" ]]; then
        echo \
            "Error: '${arr_nam_1}' must match the number of '${arr_nam_2}'." \
            "Got ${#arr_1[@]} for '${arr_nam_1}' but ${#arr_2[@]} for" \
            "'${arr_nam_2}'." >&2
        return 1
    fi
}
