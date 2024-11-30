#!/bin/bash

#  Function to validate matching values between file pairs, i.e.,
#+ first file/numerator and second file/denominator
function check_match() {
    local var_dsc="${1}"  # Description of variable being compared
    local val_num="${2}"  # Numerator value
    local val_den="${3}"  # Denominator value
    local out_nam="${4}"  # Name of output variable to assign if matched

    if [[ "${val_num}" == "${val_den}" ]]; then
        #  Assign the matched value to the output variable if specified
        if [[ -n "${out_nam}" ]]; then eval "${out_nam}=\"${val_num}\""; fi
    else
        echo \
            "Error: ${var_dsc} does not match between numerator" \
            "(${val_num}) and denominator (${val_den})." >&2
        return 1
    fi
}
