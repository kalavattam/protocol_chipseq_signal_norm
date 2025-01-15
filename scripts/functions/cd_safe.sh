#!/bin/bash

#  Function to safely change directories
function cd_safe() {
    local dir="${1}"
    local msg=${2:-false}

    #  Validate input
    if [[ -z "${dir}" ]]; then
        echo "Error: No directory specified." >&2
        return 1
    fi

    #  Validate second parameter as Boolean
    case "${msg}" in
        true|false) : ;;
        *)
            echo \
                "Error: Positional parameter must be Boolean, but got" \
                "'${msg}'." >&2
            return 1
            ;;
    esac

    #  Check that directory exists
    if [[ ! -d "${dir}" ]]; then
        echo "Error: Directory does not exist: '${dir}'." >&2
        return 1
    fi

    #  Attempt to change to directory
    if ! cd "${dir}"; then
        echo "Error: Failed to change directory to '${dir}'." >&2
        return 1
    fi

    #  Print optional success message
    if ${msg}; then echo "Changed directory to '${dir}'."; fi
}
