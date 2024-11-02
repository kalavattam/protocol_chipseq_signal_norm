#!/bin/bash

function check_int_pos() {
    local value="${1}"
    local name="${2}"
    local show_help

    show_help=$(cat << EOM
-------------
check_int_pos
-------------

Check that a value is an integer greater than or equal to (gte) 1.

Positional parameters:
  1, value (int): The value to check (required).
  2, name  (str): Name of argument or variable associated with the value
                  (optional).

Returns:
  0 if the check passes, otherwise calls return 1.

Dependencies:
  - Bash

Example:
  Check that the value assigned to the argument --threads is a positive integer
  greater than or equal to 1:
  \`\`\`
  ❯ check_int_pos 4 "threads"  # Returns 0
  
  ❯ check_int_pos a "threads"  # Returns 1 and the below message
  Error: --threads was assigned 'a' but must be a positive integer greater than or equal to 1.

  ❯ check_int_pos -2  # Returns 1 and the below message
  Error: -2 is not a positive integer greater than or equal to 1.
  \`\`\`
EOM
    )

    #  Parse and check function parameters
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    if [[ -z "${value}" ]]; then
        echo "Error: Positional parameter 1, value, is missing." >&2
        echo "" >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Perform the check and return an error message if it fails
    if ! [[ "${value}" =~ ^[1-9][0-9]*$ ]]; then
        if [[ -n "${name}" ]]; then
            echo \
                "Error: --${name} was assigned '${value}' but must be an" \
                "integer greater than or equal to 1." >&2
        else
            echo \
                "Error: ${value} is not an integer greater than or equal to" \
                "1." >&2
        fi
        return 1
    fi
}
