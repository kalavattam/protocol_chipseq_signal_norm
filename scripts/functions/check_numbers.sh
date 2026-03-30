#!/bin/bash

function check_flt_pos() {
    local value="${1}"
    local name="${2}"
    local show_help

    show_help=$(cat << EOM
-------------
check_flt_pos
-------------

Description:
  check_flt_pos checks that a value is a float greater than (gt) 0.

Positional parameters:
  1, value (int): The value to check (required).
  2, name  (str): Name of argument or variable associated with the value
                  (optional).

Returns:
  0 if the check passes, otherwise calls return 1.

Dependencies:
  - Bash or Zsh

Example:
  Check that the value assigned to the argument --scaling_factor is a positive
  float:
  \`\`\`
  ❯ check_flt_pos 4.5 "scaling_factor"  # Returns 0

  ❯ check_flt_pos a "scaling_factor"  # Returns 1 and the below message
  Error: --scaling_factor was assigned 'a' but must be a positive float.

  ❯ check_flt_pos -2.5  # Returns 1 and the below message
  Error: -2.5 is not a positive float.
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

    #  Check if the value is a positive float or integer
    if [[
            "${value}" =~ ^[+]?[0-9]+$ \
        || (
               "${value}" =~ ^[+]?[0-9]*\.[0-9]+$ \
            && "${value}" != "." \
            && "${value:0:1}" != "-" \
        )
    ]]; then
        return 0
    else
        if [[ -n "${name}" ]]; then
            echo \
                "Error: --${name} was assigned '${value}' but must be a" \
                "positive float." >&2
        else
            echo \
                "Error: ${value} is not a positive float." >&2
        fi
        return 1
    fi
}
