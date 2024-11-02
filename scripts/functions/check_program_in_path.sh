#!/bin/bash

function check_program_in_path() {
    local program="${1}"
    local show_help

    show_help=$(cat << EOM
---------------------
check_program_in_path
---------------------

Description:
  check_program_in_path checks if a given program is available in the system's
  PATH. If the program is not found, it prints an error message and returns
  exit code 1.

Positional parameter:
  1, program (str): The name of the program to check (required).

Returns:
  0 if the program is found in the PATH; otherwise, prints an error message and
  returns exit code 1.

Dependencies:
  - Bash or Zsh

Example:
  \`\`\`
  #  Check if the program "samtools" is in the PATH.
  check_program_in_path "samtools"
  # If "samtools" is not in the PATH, an error message is printed.
  \`\`\`
EOM
    )

    #  Parse and check function parameters
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    if ! command -v "${program}" &> /dev/null; then
        echo \
            "Error: ${program} is not in PATH. Please install ${program} or" \
            "add it to PATH." >&2
        return 1
    fi
}
