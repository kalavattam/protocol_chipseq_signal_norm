#!/bin/bash

#  Function to note that scaling factor(s) or normalization will be multiplied
#+ with those in the table.
function check_table_scaling_factor() {
    local type="${1,,}"   # Expected values: 'string' or 'boolean', etc.
    local table="${2}"    # Path to table file
    local scl_fct="${3}"  # Scaling factor(s) (string) or norm. flag (boolean)
    local name="${4}"     # Option name being checked (e.g., 'scl_fct', 'typ_cvg')
    local msg

    #  Validate positional parameters
    if [[ -z "${type}" ]]; then
        echo "Error: Positional parameter 1, 'type', is required." >&2
        return 1
    fi

    if [[ ! "${type}" =~ ^(str|string|bol|bool|boolean|flg|flag)$ ]]; then
        echo \
            "Error: Invalid positional parameter 1, 'type': '${type}'." \
            "Expected 'str', 'string', 'bol', 'bool' or 'boolean'." >&2
        return 1
    fi

    if [[ -z "${table}" ]]; then
        echo "Error: Positional parameter 2, 'table', is required." >&2
        return 1
    fi

    if [[ ! -f "${table}" ]]; then
        echo \
            "Error: Positional parameter 2, 'table', does not exist:" \
            "'${table}'." >&2
        return 1
    fi

    if [[ -z "${scl_fct}" ]]; then
        echo "Error: Positional parameter 3, 'scl_fct', is required." >&2
        return 1
    fi

    if [[ "${type}" == "string" && "${scl_fct}" == "NA" ]]; then
        echo \
            "Error: Positional parameter 3, 'scl_fct', is assigned 'NA'," \
            "which is invalid for type (positional parameter 1) 'string'." >&2
        return 1
    fi  #TODO: Not sure if this is needed...

    if [[
        "${type}" == "boolean" && "${scl_fct}" != true && "${scl_fct}" != false
    ]]; then
        echo \
            "Error: Positional parameter 3, 'scl_fct', must be 'true' or" \
            "'false' for type (positional parameter 1) 'boolean'." >&2
        return 1
    fi

    if [[ -z "${name}" ]]; then
        echo "Error: Positional parameter 4, 'name', is required." >&2
        return 1
    fi

    #  Determine appropriate message
    if [[ "${name}" == "typ_cvg" ]]; then
        msg="Note: --${name} scaling factors will be multiplied with those"
        msg+=" from --table (--tbl_col) if present; otherwise, --${name}"
        msg+=" scaling factors will be applied directly to raw coverage."
    elif [[ "${name}" == "scl_fct" && -n "${scl_fct}" ]]; then
        msg="Note: --${name} will override scaling factors from --table"
        msg+=" (--tbl_col)."
    fi
    
    #  Output the message if applicable
    if [[ -n "${msg}" ]]; then echo "${msg}"; fi
}
