#!/bin/bash

#  Function to vaildate existence of column in table
function check_table_column() {
    local table="${1}"
    local column="${2}"

    if \
        ! awk -F '\t' -v col="${column}" '
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    if ($i == col) exit 0
                } exit 1
            }
        ' "${table}"
    then
        echo \
            "Error: Column '${column}' not found in table header:" \
            "$(awk 'NR == 1' "${table}" | tr '\t' ' ')"
        return 1
    fi
}
