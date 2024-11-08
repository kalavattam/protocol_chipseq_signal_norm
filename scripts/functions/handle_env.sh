#!/bin/bash

#  Function to deactivate a Conda/Mamba environment; requirements: Bash or Zsh,
#+ and Conda
#+ 
#+ handle_env_deactivate first checks if deactivation is necessary; if so, it
#+ sources an appropriate Conda shell hook and uses 'conda deactivate' for
#+ environment deactivation
function handle_env_deactivate() {
    local shl_cur

    #  If no Conda environment is active or the active environment is "base",
    #  then deactivation is not needed; return successfully (0)
    if [[
        -z "${CONDA_DEFAULT_ENV}"
        || "${CONDA_DEFAULT_ENV}" == "base"
    ]]; then
        return 0
    fi

    #  Determine the current shell (e.g., bash, zsh)
    shl_cur=$(basename "${SHELL}")

    #  Source the appropriate Conda shell hook based on the detected shell
    case "${shl_cur}" in
        bash)
            eval "$(conda shell.bash hook)"
            ;;
        zsh)
            eval "$(conda shell.zsh hook)"
            ;;
        *)
            #  If an unsupported shell is detected, print a warning
            echo \
                "Warning: Unsupported shell '${shl_cur}'; environment" \
                "deactivation may not work as expected." >&2
            ;;
    esac

    #  Attempt to deactivate using Conda if available
    if command -v conda &> /dev/null; then
        conda deactivate ||
            {
                echo "Error: Failed to deactivate environment using conda." >&2
                return 1
            }
    fi
}


function handle_env_activate() {
    local env_nam="${1}"
    local shl_cur
    local show_help

    show_help=$(cat << EOM
-------------------
handle_env_activate
-------------------

Description:
  handle_env_activate activates a specified Conda or Mamba environment. First,
  it sources the appropriate Conda shell hook for the current shell (e.g., Bash
  or Zsh). It then attempts to activate the environment using 'mamba'. If that
  fails, it tries using 'conda'. If both 'mamba' and 'conda' fail, it attempts
  to activate the environment using the 'source activate' command.

Positional parameter:
  1, env_nam (str): The name of the environment to activate (required).

Returns:
  0 if the environment is activated successfully, 1 if the environment fails to
  activate.

Dependencies:
  - Bash or Zsh
  - Mamba or Conda

Example:
  \`\`\`
  #  Activate an environment named "my_env"
  ❯ handle_env_activate "my_env"
  Environment "my_env" activated successfully.

  #  Activate an environment that does not exist
  ❯ handle_env_activate "monkey"

  \`\`\`
EOM
    )

    #  Parse and check function parameter
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    if [[ -z "${env_nam}" ]]; then
        echo "Error: Positional parameter 1, env_nam, is required." >&2
        echo "" >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Determine the current shell (e.g., bash, zsh) and source the appropriate
    #+ Conda shell hook
    shl_cur=$(basename "${SHELL}")
    case "${shl_cur}" in
        bash)
            eval "$(conda shell.bash hook)"
            ;;
        zsh)
            eval "$(conda shell.zsh hook)"
            ;;
        *)
            #  Warn if the shell is unsupported; environment activation may
            #+ still proceed but may not work as expected
            echo \
                "Warning: Unsupported shell '${shl_cur}'; environment" \
                "activation may not work as expected." >&2
            ;;
    esac

    #  Check that the specified environment exists using conda/mamba
    #+ environment list
    if \
           ! (conda env list | grep -qw "${env_nam}") \
        && ! (mamba env list | grep -qw "${env_nam}");
    then
        echo "Error: Environment \"${env_nam}\" does not exist." >&2
        return 1
    fi

    #  Attempt to activate the specified environment using Mamba, Conda, or a
    #+ fallback call to 'source activate'
    if ! mamba activate "${env_nam}" &> /dev/null; then
        if ! conda activate "${env_nam}" &> /dev/null; then
            # shellcheck disable=SC1091
            if ! source activate "${env_nam}" &> /dev/null; then
                echo \
                    "Error: Failed to activate environment" \
                    "\"${env_nam}\"." >&2
                return 1
            fi
        fi
    fi
    
    echo "Environment \"${env_nam}\" activated successfully."
    return 0
}


#  Function to ensure the specified Conda/Mamba environment is activated; if a
#+ different environment is currently active, it switches to the specified one
#+ 
#+ Requirements:
#+   - Bash or Zsh
#+   - Mamba (for handle_env_activate), Conda (for handle_env_activate and
#+     handle_env_deactivate})
function handle_env() {
    local env_nam="${1}"  # The name of the environment to activate

    #  Check if no Conda environment is currently active or if the active
    #+ environment is "base"; in either case, activate the specified
    #+ environment
    if [[
        -z "${CONDA_DEFAULT_ENV}"
        || "${CONDA_DEFAULT_ENV}" == "base"
    ]]; then
        handle_env_activate "${env_nam}"
    #  If a different environment is active, deactivate it first, then activate
    #+ the specified one
    elif [[ "${CONDA_DEFAULT_ENV}" != "${env_nam}" ]]; then
        handle_env_deactivate
        handle_env_activate "${env_nam}"
    fi
}
