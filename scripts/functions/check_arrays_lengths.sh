#!/bin/bash

#  Function to validate that two indexed arrays have matching lengths (usable
#+ with Bash ≥3.2)
function check_arrays_lengths() {
    local arr_nam_1="${1}"
    local arr_nam_2="${2}"
    local arr_siz_1 arr_siz_2

    #  Ensure array names are valid
    if [[
           ! "${arr_nam_1}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ \
        || ! "${arr_nam_2}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$
    ]]; then
        echo "Error: Invalid array name '${arr_nam_1}' and/or '${arr_nam_2}'." >&2
        return 1
    fi

    #  Ensure both names refer to arrays
    if ! \
        eval 'declare -p '"${arr_nam_1}"' &>/dev/null'
    then
        echo "Error: '${arr_nam_1}' is unset." >&2
        return 1
    fi

    if ! \
        eval 'declare -p '"${arr_nam_2}"' &>/dev/null'
    then
        echo "Error: '${arr_nam_2}' is unset." >&2
        return 1
    fi

    if ! \
        eval 'declare -p '"${arr_nam_1}"' 2>/dev/null | grep -q "^declare \-a "'
    then
        echo "Error: '${arr_nam_1}' is not an indexed array." >&2
        return 1
    fi

    if ! \
        eval 'declare -p '"${arr_nam_2}"' 2>/dev/null | grep -q "^declare \-a "'
    then
        echo "Error: '${arr_nam_2}' is not an indexed array." >&2
        return 1
    fi

    #  Get array sizes using indirect references
    eval "arr_siz_1=\${#${arr_nam_1}[@]}"
    eval "arr_siz_2=\${#${arr_nam_2}[@]}"

    #  Check that array sizes match
    # shellcheck disable=SC2154
    if [[ "${arr_siz_1}" -ne "${arr_siz_2}" ]]; then
        echo \
            "Error: '${arr_nam_1}' must match the number of '${arr_nam_2}'." \
            "Got '${arr_siz_1}' for '${arr_nam_1}' but '${arr_siz_2}' for" \
            "'${arr_nam_2}'." >&2
        return 1
    fi

    return 0
}
