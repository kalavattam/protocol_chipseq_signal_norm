#!/bin/bash

#  Function to validate TSV table infile is not empty
function check_table() {
    local table="${1}"

    if [[ $(wc -l < "${table}") -le 1 ]]; then
        echo \
            "Error: Table file '${table}' is empty or contains only a" \
            "header." >&2
        return 1
    fi
}
