#!/bin/bash

function check_installed_env() {
    local env_name="${1}"
    local show_help

    show_help=$(cat << EOM
-------------------
check_installed_env
-------------------

Description:
  check_installed_env checks that a specific Conda or Mamba environment is
  installed. If the environment is not found, it prints an error message and
  returns exit code 1.

Positional parameters:
  1, env_name (str): The name of the Conda/Mamba environment to check
                     (required).

Returns:
  0 if the environment is installed; otherwise, prints an error message and 
  returns exit code 1.

Dependencies:
  - Bash
  - Conda and/or Mamba

Example:
  Check that the environment "env_R" is installed.
  \`\`\`
  ❯ check_installed_env "env_R"
  Error: Environment "env_R" is not installed.
  \`\`\`
EOM
    )

    #  Parse and check function parameters
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    #  Check that the specific Conda or Mamba environment is installed
    if ! conda env list | grep -q "^${env_name} "; then
        return 1
    fi
}
