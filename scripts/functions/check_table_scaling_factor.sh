#!/bin/bash

#  Function to warn that scaling factor(s) or normalization override(s) those
#+ in table
function check_table_scaling_factor() {
    local type="${1,,}"   # Expected values: 'string' or 'boolean', etc.
    local table="${2}"    # Path to table file
    local scl_fct="${3}"  # Scaling factor(s) (string) or norm. flag (boolean)
    local name="${4}"     # Option name being checked (e.g., 'scl_fct', 'norm')

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

    if [[ "${type}" == "string" && "${scl_fct}" == "#N/A" ]]; then
        echo \
            "Error: Positional parameter 3, 'scl_fct', is assigned '#N/A'," \
            "which is invalid for type (positional parameter 2) 'string'." >&2
        return 1
    fi  #TODO: Not sure if this is needed...

    if [[
        "${type}" == "boolean" && "${scl_fct}" != true && "${scl_fct}" != false
    ]]; then
        echo \
            "Error: Positional parameter 3, 'scl_fct', must be 'true' or" \
            "'false' for type (positional parameter 2) 'boolean'." >&2
        return 1
    fi

    if [[ -z "${name}" ]]; then
        echo "Error: Positional parameter 4, 'name', is required." >&2
        return 1
    fi

    #  Define warning message
    local msg="Warning: --${name} will override scaling factors from --table."
    
    #  Check the type and issue a warning if necessary
    case "${type}" in
        str|string)
            if [[ -n "${table}" && -n "${scl_fct}" ]]; then
                echo "${msg}" >&2
            fi
            ;;
        bol|bool|boolean|flg|flag)
            if [[ -n "${table}" ]] && ${scl_fct}; then
                echo "${msg}" >&2
            fi
            ;;
        *)
            echo \
                "Error: Invalid positional parameter 1, 'type': '${type}'." \
                "Expected 'str', 'string', 'bol', 'bool' or 'boolean'." >&2
            return 1
            ;;
    esac
}
