#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: stub_package_manager.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5.5).
#
# Distributed under the MIT license.


# Require Bash >= 4.4 before doing any work.
if [[ -z "${BASH_VERSION:-}" ]] || (( BASH_VERSINFO[0] < 4 )); then
    echo "error(fake_package_manager):" \
        "Bash >= 4.4 is required." >&2
    exit 1
fi

# Run in safe mode, exiting on errors, unset variables, and pipe failures.
set -euo pipefail

# Identify the fake command selected through its symlink.
cmd_nam="${0##*/}"


# Record every package-manager call in the caller-provided fixture log.
function log_invocation() {
    {
        printf '%s' "${cmd_nam}"
        printf ' %q' "$@"
        printf '\n'
    } >> "${INSTALL_ENVS_FAKE_LOG}"
}


# Print environments recorded in the caller-provided fixture state file.
function print_environments() {
    local dir_envs=""
    local env_nam=""

    [[ -f "${INSTALL_ENVS_FAKE_ENVS}" ]] || : > "${INSTALL_ENVS_FAKE_ENVS}"

    dir_envs="$(dirname "${INSTALL_ENVS_FAKE_ENVS}")"

    printf '# conda environments:\n'
    while IFS= read -r env_nam; do
        [[ -n "${env_nam}" ]] || continue
        printf '%s  %s/%s\n' "${env_nam}" "${dir_envs}" "${env_nam}"
    done < "${INSTALL_ENVS_FAKE_ENVS}"
}


# Record a newly created environment by the name declared in its YAML file.
function create_environment() {
    local env_nam=""
    local yml=""

    if [[ "${INSTALL_ENVS_FAKE_CREATE_FAIL:-false}" == "true" ]]; then
        echo "error(fake_package_manager):" \
            "forced environment-create failure." >&2
        return 1
    fi

    while (( $# > 0 )); do
        case "${1}" in
            -f)
                yml="${2:-}"
                shift 2
                ;;

            *) shift ;;
        esac
    done

    env_nam="$(sed -n 's/^name:[[:space:]]*//p' "${yml}" | head -n 1)"
    [[ -n "${env_nam}" ]] || {
        echo "error(fake_package_manager):" \
            "YAML environment name is missing." >&2
        return 1
    }

    printf '%s\n' "${env_nam}" >> "${INSTALL_ENVS_FAKE_ENVS}"
}


# Emulate the package-manager operations used by the installer contract.
function main() {
    case "${cmd_nam}:${1:-}:${2:-}" in
        mamba:--version:|conda:--version:)
            printf '%s 1.5.9\n' "${cmd_nam}"
            ;;

        conda:env:list) print_environments ;;

        mamba:env:create)
            create_environment "$@"
            ;;

        mamba:install:*)
            if [[ "${INSTALL_ENVS_FAKE_INSTALL_FAIL:-false}" == "true" ]]
            then
                echo "error(fake_package_manager):" \
                    "forced install failure." >&2
                return 1
            fi
            ;;

        mamba:run:*)
            if [[ "${INSTALL_ENVS_FAKE_PIP_FAIL:-false}" == "true" ]]; then
                echo "error(fake_package_manager):" \
                    "forced pip failure." >&2
                return 1
            fi
            ;;

        *)
            echo "error(fake_package_manager):" \
                "unsupported call: ${cmd_nam} $*." >&2
            return 1
            ;;
    esac
}


log_invocation "$@"
main "$@"
