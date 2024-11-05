#!/bin/bash

#  Function to ensure the specified environment is activated, switching from
#+ the current one if needed; requirements: Bash or Zsh, Mamba or Conda
function handle_env() {
    local env_nam="${1}"

    # Check the current environment and switch if necessary
    if [[ -z "${CONDA_DEFAULT_ENV}" || "${CONDA_DEFAULT_ENV}" == "base" ]]; then
        handle_env_activate "${env_nam}"
    elif [[ "${CONDA_DEFAULT_ENV}" != "${env_nam}" ]]; then
        handle_env_deactivate
        handle_env_activate "${env_nam}"
    fi
}