#!/bin/bash

function check_mut_excl_flags() {
    local flg_1="${1}"
    local nam_1="${2}"
    local flg_2="${3}"
    local nam_2="${4}"
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
  1, flg_1 (bool): The first flag to check.
  2, nam_1  (str): The name of the first flag.
  3, flg_2 (bool): The second flag to check.
  4, nam_2  (str): The name of the second flag.

Returns:
  0 if the flags are mutually exclusive; otherwise, if both flags are
  specified, then returns and error message along with exit code 1.

Dependencies:
  - Bash or Zsh

Example:
  \`\`\`
  #  Check that either --flg_1 or --flg_2 is specified, but not both.
  flg_1=true
  flg_2=false
  check_mut_excl_flags "\${flg_1}" "flg_1" "\${flg_2}" "flg_2"  # Returns 0

  #  What happens if --flg_1 are --flg_2 are both specified?
  flg_1=true
  flg_2=true
  check_mut_excl_flags "\${flg_1}" "flg_1" "\${flg_2}" "flg_2"
  # Error: Only one of --flg_1 or --flg_2} can be specified at a time.
  \`\`\`
EOM
    )

    #  Parse and check function parameters
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    #  Perform the check for two-flag mutual exclusivity
    if ${flg_1} && ${flg_2}; then
        echo \
            "Error: Only one of --${nam_1} or --${nam_2} can be specified at" \
            "a time." >&2
        return 1
    fi
}
