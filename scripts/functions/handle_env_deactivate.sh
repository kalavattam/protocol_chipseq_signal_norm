#!/bin/bash

#  Function to deactivate a Conda/Mamba environment; requirements: Bash or Zsh,
#+ Mamba or Conda
function handle_env_deactivate() {
    if [[ "${CONDA_DEFAULT_ENV}" != "base" ]]; then
        if ! mamba deactivate &> /dev/null; then
            if ! conda deactivate &> /dev/null; then
                # shellcheck disable=1091
                if ! source deactivate &> /dev/null; then
                    echo "Error: Failed to deactivate environment." >&2
                    return 1
                fi
            fi
        fi
    fi
}
