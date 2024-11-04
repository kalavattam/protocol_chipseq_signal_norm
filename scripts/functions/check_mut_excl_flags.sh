#!/bin/bash

function check_mut_excl_flags() {
    local flag_1="${1}"
    local name_1="${2}"
    local flag_2="${3}"
    local name_2="${4}"
    local show_help

    show_help=$(cat << EOM
--------------------
check_mut_excl_flags
--------------------

Description:
  check_mut_excl_flags checks if two mutually exclusive flags are erroneously
  specified at the same time. If both flags are true, it prints an error
  message and returns exit code 1.

Positional parameters:
  1, flag_1 (bool): The first flag to check (required).
  2, name_1  (str): The name of the first flag (required).
  3, flag_2 (bool): The second flag to check (required).
  4, name_2  (str): The name of the second flag (required).

Returns:
  0 if the flags are mutually exclusive; otherwise, if both flags are
  specified, then returns and error message along with exit code 1.

Dependencies:
  - Bash or Zsh

Example:
  \`\`\`
  #  Check that either --flag_1 or --flag_2 is specified, but not both.
  ❯ flag_1=true
  ❯ flag_2=false
  ❯ check_mut_excl_flags "\${flag_1}" "flag_1" "\${flag_2}" "flag_2"
  #  Returns 0.

  #  What happens if --flag_1 are --flag_2 are both specified?
  ❯ flag_1=true
  ❯ flag_2=true
  ❯ check_mut_excl_flags "\${flag_1}" "flag_1" "\${flag_2}" "flag_2"
  Error: Only one of --flag_1 or --flag_2} can be specified at a time.
  \`\`\`
EOM
    )

    #  Parse and check function parameters
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    #  Perform the check for two-flag mutual exclusivity
    if ${flag_1} && ${flag_2}; then
        echo \
            "Error: Only one of --${name_1} or --${name_2} can be specified" \
            "at a time." >&2
        return 1
    fi
}
