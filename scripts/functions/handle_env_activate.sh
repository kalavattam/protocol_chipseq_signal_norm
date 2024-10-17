#!/bin/bash

function handle_env_activate() {
    local env_name="${1}"
    local show_help

    show_help=$(cat << EOM
-------------------
handle_env_activate
-------------------

Description:
  handle_env_activate activates a specified Conda or Mamba environment. It
  attempts to activate the environment using 'mamba'. If that fails, it tries
  using 'conda'. If both 'mamba' and 'conda' fail, it attempts to activate the
  environment using the 'source activate' command.

Positional parameters:
  1, env_name (str): The name of the environment to activate (required).

Returns:
  0 if the environment is activated successfully, 1 if the environment fails to
  activate.

Dependencies:
  - Bash or Zsh
  - Mamba (preferred) or Conda

Example:
  \`\`\`
  # Activate the environment named "my_env"
  handle_env_activate "my_env"
  # Output: Environment "my_env" activated successfully.
  \`\`\`
EOM
    )

    #  Parse and check function parameters
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    if [[ -z "${env_name}" ]]; then
        echo "Error: Positional parameter 1, env_name, is required." >&2
        echo "" >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! mamba activate "${env_name}" &> /dev/null; then
        if ! conda activate "${env_name}" &> /dev/null; then
            # shellcheck disable=1091
            if ! source activate "${env_name}" &> /dev/null; then
                echo \
                    "Error: Failed to activate environment" \
                    "\"${env_name}\"." >&2
                return 1
            fi
        fi
    fi
    
    echo "Environment \"${env_name}\" activated successfully."
    return 0
}
