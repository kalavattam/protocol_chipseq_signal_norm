#!/bin/bash

#  Debug array contents (usable with Bash >=3.2)
function debug_array_contents() {
    local arr_nam
    local -a arr

    for arr_nam in "$@"; do
        #  Skip invalid variable names
        if ! [[ "${arr_nam}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ ]]; then
            continue
        fi

        #  Check array exists (if unset / not an array, behavior can get messy)
        if ! eval 'declare -p '"${arr_nam}"' >/dev/null 2>&1'; then
            continue
        fi

        #  Skip non-indexed variables
        if ! eval '
            declare -p '"${arr_nam}"' 2>/dev/null | grep -q "^declare \-a "
        '; then
            continue
        fi

        #  Access the array indirectly using eval
        eval "arr=( \"\${${arr_nam}[@]}\" )"
        if [[ -n "${arr[*]}" ]]; then
            echo "  - ${arr_nam}=( ${arr[*]} )"
        fi
    done
}
