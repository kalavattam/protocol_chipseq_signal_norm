#!/bin/bash

#  Function to populate array with '#N/A' if it is empty
function populate_array_empty() {
    local -n arr="${1}"
    local target="${2}"

    if [[ -z "${arr[*]}" ]]; then
        for ((i = 0; i < target; i++)); do arr+=( "#N/A" ); done
    fi
}
