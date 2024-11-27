#!/bin/bash

#  Function to debug array contents
function debug_array_contents() {
    for arr_nam in "$@"; do
        #  Access the array indirectly using eval
        eval "arr=( \"\${${arr_nam}[@]}\" )"
        if [[ -n "${arr[*]}" ]]; then
            echo "  - ${arr_nam}: ( ${arr[*]} )"
        fi
    done
}
