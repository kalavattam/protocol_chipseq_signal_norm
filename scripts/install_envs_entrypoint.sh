#!/bin/sh
# -*- coding: utf-8 -*-
#
# Script: install_envs_entrypoint.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.4) was used in development.
#
# Distributed under the MIT license.


show_help() {
    cat << EOM >&2
Usage:
  install_envs_entrypoint.sh [-h|--hlp|--help]
  install_envs_entrypoint.sh [keyword arguments for install_envs.sh]

Description:
  POSIX entrypoint for repository environment setup.

  This script is intended to be invokable from various user shell environments (e.g., terminals used with 'zsh', 'csh', 'fish', etc.) and systems before a modern Bash environment has been set up.

  If Bash >= 5 is already available in PATH, this script hands off directly to 'install_envs.sh'.

  If Bash >= 5 is not available, this script prints guidance for the next step. If Conda or Mamba is available, it explains how to install Bash >= 5 in the current environment. If Conda or Mamba is not available, it explains that Miniforge should be installed and initialized first.

Returns:
  0 on success; non-zero on error.

Notes:
  - This script is POSIX-sh compatible.
  - Repository-specific installation logic is in 'install_envs.sh'.
EOM
}


err() { printf '%s\n' "error(install_envs_entrypoint.sh): $*" >&2; }


# warn() { printf '%s\n' "warning(install_envs_entrypoint.sh): $*" >&2; }


# info() { printf '%s\n' "info(install_envs_entrypoint.sh): $*" >&2; }


find_bash_5() {
    #  Prefer an already-available 'bash' if it is Bash >= 5
    if command -v bash >/dev/null 2>&1; then
        if bash -c '
            [ -n "${BASH_VERSION:-}" ] || exit 1
            [ "${BASH_VERSINFO[0]:-0}" -ge 5 ] || exit 1
        ' >/dev/null 2>&1; then
            command -v bash
            return 0
        fi
    fi

    #  Common fallback locations for user-installed Bash on macOS/Linux
    for pth_bash in \
        /opt/homebrew/bin/bash \
        /usr/local/bin/bash \
        /bin/bash \
        /usr/bin/bash
    do
        if [ -x "${pth_bash}" ]; then
            # shellcheck disable=SC2016
            if "${pth_bash}" -c '
                [ -n "${BASH_VERSION:-}" ] || exit 1
                [ "${BASH_VERSINFO[0]:-0}" -ge 5 ] || exit 1
            ' >/dev/null 2>&1; then
                printf '%s\n' "${pth_bash}"
                return 0
            fi
        fi
    done

    return 1
}


print_guidance_bash() {
    cat << EOM >&2
Bash >= 5 was not found in PATH.

If Bash >= 5 is already installed, then make sure it is available in PATH, for example by activating the appropriate Conda/Mamba environment or adjusting shell initialization.

If Bash >= 5 is not installed, then install it in the current Conda/Mamba context, e.g.,

  mamba install -c conda-forge bash

or

  conda install -c conda-forge bash

Then verify it, e.g.,

  bash --version

After that, rerun

  sh install_envs_entrypoint.sh [arguments]

or run 'install_envs.sh' directly with the newer Bash, e.g.,

  /path/to/bash install_envs.sh [arguments]

EOM
}


print_guidance_miniforge() {
    cat << EOM >&2
The package managers Conda and Mamba were not found in PATH.

If Conda and/or Mamba is already installed, make sure at least one is initialized and available in PATH.

If neither Conda nor Mamba is installed, Miniforge (https://github.com/conda-forge/miniforge/) is the recommended starting point.

  1. Determine operating system and architecture:

       uname -s
       uname -m

  2. Download the appropriate Miniforge installer. Installation instructions are available here:

       https://github.com/conda-forge/miniforge/

  3. Run the installer and allow Conda initialization.

  4. Restart the terminal, or source the updated shell configuration file.

  5. Verify installation:

       conda --version
       mamba --version

(For more information, see the Tsukiyama Lab Bio-protocol manuscript [PMID 40364978].)

After Conda and/or Mamba has been installed and initialized, install Bash >= 5, e.g.,

  mamba install -c conda-forge bash

or

  conda install -c conda-forge bash

Then rerun

  sh install_envs_entrypoint.sh [keyword arguments]

EOM
}


main() {
    argv="$*"

    #  Parse arguments lightly so that help requests and clearly unexpected
    #+ options can be handled even before Bash >= 5 is available
    while [ "$#" -gt 0 ]; do
        case "${1}" in
            -h|--hlp|--help)
                show_help
                exit 0
                ;;
            -dr|--dry_run|--dry-run|-y|--yes)
                shift 1
                ;;
            -en|--env_nam|--env-nam)
                if [ "$#" -lt 2 ] || [ -z "${2}" ] || [ "${2#-}" != "${2}" ]; then
                    err "option '${1}' requires a value."
                    echo >&2
                    show_help
                    exit 1
                fi
                shift 2
                ;;
            -*)
                echo "## Unknown keyword argument passed: '${1}' ##" >&2
                echo >&2
                show_help
                exit 1
                ;;
            *)
                echo "## Unknown positional argument passed: '${1}' ##" >&2
                echo >&2
                show_help
                exit 1
                ;;
        esac
    done

    dir_scr=$(
        CDPATH=
        cd -- "$(dirname -- "$0")" 2>/dev/null || exit 1
        pwd
    ) || {
        err "failed to determine script directory."
        exit 1
    }

    pth_main="${dir_scr}/install_envs.sh"

    if [ ! -f "${pth_main}" ]; then
        err "required handoff script does not exist: '${pth_main}'."
        exit 1
    fi

    if [ ! -r "${pth_main}" ]; then
        err "required handoff script is not readable: '${pth_main}'."
        exit 1
    fi

    if pth_bash="$(find_bash_5)"; then
        # shellcheck disable=SC2086
        exec "${pth_bash}" "${pth_main}" ${argv}
    fi

    if ! (
           command -v conda >/dev/null 2>&1 \
        || command -v mamba >/dev/null 2>&1
    ); then
        err \
            "Bash >= 5 was not found in PATH, and neither 'conda' nor" \
            "'mamba' is currently available in PATH."
        echo >&2
        print_guidance_miniforge
        exit 1
    fi

    err "Bash >= 5 was not found in PATH."
    echo >&2
    print_guidance_bash
    exit 1
}


main "$@"
