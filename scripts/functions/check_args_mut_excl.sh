#!/bin/bash

function check_mut_excl_args() {
    local nam_1="${1}"
    local val_1="${2}"
    local nam_2="${3}"
    local val_2="${4}"
    local show_help

    show_help=$(cat << EOM
Description:
  Checks that two mutually exclusive arguments are not specified at the same
  time. If both are specified or neither is specified, an error is returned.

Positional parameters:
  1, nam_1 (str): Name of the first argument (e.g., "norm").
  2, val_1 (str): Value of the first argument.
  3, nam_2 (str): Name of the second argument (e.g., "raw").
  4, val_2 (str): Value of the second argument.

Returns:
  0 if the arguments are mutually exclusive and valid, 1 and an error message
  if both arguments are specified or neither is specified.

Usage:
  check_mut_excl_args
      "arg_1_name" "\${arg_1_value}"
      "arg_2_name" "\${arg_2_value}"

Dependencies:
  - Bash or Zsh

Example:
  \`\`\`
  #TODO
  \`\`\`
EOM
    )

    #  Print help message if no arguments are passed or if help is requested
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    #  Validate mutually exclusive arguments
    if [[ -n "${val_1}" && -n "${val_2}" ]]; then
        echo \
            "Error: Only one of --${nam_1} or --${nam_2} can be specified at" \
            "a time." >&2
        return 1
    elif [[ -z "${val_1}" && -z "${val_2}" ]]; then
        echo "Error: One of --${nam_1} or --${nam_2} must be specified." >&2
        return 1
    fi
}
