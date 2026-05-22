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


show_help_entrypoint() {
    cat << EOM >&2
Usage:
  install_envs_entrypoint.sh [-h|--hlp|--help]
  install_envs_entrypoint.sh --env_nam <spec> [--dry_run] [--if_exis <spec>] [--yes]

Description:
  POSIX entrypoint for repository environment setup. The examples use Bash, but this entrypoint can also be invoked with 'sh' or run directly. It will attempt to locate an appropriate Bash for the handoff to 'install_envs.sh'.

Recommended first-time command:
  '''bash
  bash scripts/install_envs_entrypoint.sh --env_nam env_protocol --yes
  '''

Environment names:
  env_protocol  Main workflow environment.
  env_analyze   Analysis environment.
  env_siqchip   siQ-ChIP environment.

Keyword arguments:
  -h, --hlp, --help  <flag>
    Display this help message and exit.

  -dr, --dry, --dry_run  <flag>
    Print the resolved installation command without creating an environment.

  -en, --env, --env_nam  <spec>
    Environment to create: 'env_protocol', 'env_analyze', or 'env_siqchip'.

  -ie, --if_ex, --if_exis, --if_exists  <spec>
    What to do if the requested environment already exists: 'fail' or 'reuse' (default: 'fail').

  -y, --yes  <flag>
    Automatically answer yes to package-manager prompts.

Notes:
  - Mamba is preferred and will be used by 'install_envs.sh' when available.
  - Old Conda/Anaconda installations may be slow. Miniforge is recommended for a current Conda/Mamba setup.
  - Invoke this script with 'bash', 'sh', or as an executable script. Do not force another interpreter such as 'zsh'.
  - This entrypoint requires Mamba or Conda and a Bash >= 4.4 handoff interpreter.
  - Full package details are available with:
    '''bash
    bash scripts/install_envs.sh --help
    '''

Returns:
  0 on success; non-zero on error.
EOM
}


err() { printf '%s\n' "error(install_envs_entrypoint.sh): $*" >&2; }


# warn() { printf '%s\n' "warning(install_envs_entrypoint.sh): $*" >&2; }


# info() { printf '%s\n' "info(install_envs_entrypoint.sh): $*" >&2; }


find_bash_required() {
    #  Prefer an already-available 'bash' if it is Bash >= 4.4
    if command -v bash >/dev/null 2>&1; then
        if bash -c '
            [ -n "${BASH_VERSION:-}" ] || exit 1
            major="${BASH_VERSINFO[0]:-0}"
            minor="${BASH_VERSINFO[1]:-0}"
            [ "${major}" -gt 4 ] || { [ "${major}" -eq 4 ] && [ "${minor}" -ge 4 ]; }
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
                major="${BASH_VERSINFO[0]:-0}"
                minor="${BASH_VERSINFO[1]:-0}"
                [ "${major}" -gt 4 ] || { [ "${major}" -eq 4 ] && [ "${minor}" -ge 4 ]; }
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
Bash >= 4.4 was not found in PATH.

If Bash >= 4.4 is already installed, then make sure it is available in PATH, for example by activating the appropriate environment or adjusting shell initialization.

Conda or Mamba must also be available in PATH before this entrypoint can hand off to 'install_envs.sh'.

If Bash >= 4.4 is not installed, then install it in the current Conda/Mamba context, e.g.,

  mamba install -c conda-forge bash

or

  conda install -c conda-forge bash

Then verify it, e.g.,

  bash --version

After that, rerun

  sh scripts/install_envs_entrypoint.sh [arguments]

or run 'install_envs.sh' directly with the newer Bash, e.g.,

  /path/to/bash scripts/install_envs.sh [arguments]

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

After Conda and/or Mamba have been installed and initialized, install Bash >= 4.4, e.g.,

  mamba install -c conda-forge bash

or

  conda install -c conda-forge bash

Then rerun

  sh scripts/install_envs_entrypoint.sh [keyword arguments]

EOM
}


check_args_light() {
    while [ "$#" -gt 0 ]; do
        case "${1}" in
            -h|--hlp|--help)
                show_help_entrypoint
                exit 0
                ;;
            -dr|--dry|--dry_run|--dry-run|-y|--yes)
                shift 1
                ;;
            -en|--env|--env_nam|--env-nam)
                if \
                       [ "$#" -lt 2 ] \
                    || [ -z "${2}" ]  \
                    || [ "${2#-}" != "${2}" ]
                then
                    err "option '${1}' requires a value."
                    echo >&2
                    show_help_entrypoint
                    exit 1
                fi
                shift 2
                ;;
            -ie|--if_ex|--if-ex|--if_exis|--if-exis|--if_exists|--if-exists)
                if \
                       [ "$#" -lt 2 ] \
                    || [ -z "${2}" ]  \
                    || [ "${2#-}" != "${2}" ]
                then
                    err "option '${1}' requires a value."
                    echo >&2
                    show_help_entrypoint
                    exit 1
                fi

                case "${2}" in
                    fail|reuse) : ;;
                    *)
                        err \
                            "invalid value for '${1}': '${2}'. Must be" \
                            "'fail' or 'reuse'."
                        echo >&2
                        show_help_entrypoint
                        exit 1
                        ;;
                esac

                shift 2
                ;;
            -*)
                echo "## Unknown keyword argument passed: '${1}' ##" >&2
                echo >&2
                show_help_entrypoint
                exit 1
                ;;
            *)
                echo "## Unknown positional argument passed: '${1}' ##" >&2
                echo >&2
                show_help_entrypoint
                exit 1
                ;;
        esac
    done
}


main() {
    if [ "$#" -eq 0 ]; then
        show_help_entrypoint
        exit 0
    fi

    check_args_light "$@"

    dir_nam=$(dirname "$0")
    dir_scr=$(
        CDPATH=
        cd "${dir_nam}" 2>/dev/null || exit 1
        pwd
    ) || {
        err "failed to determine script directory."
        exit 1
    }
    unset dir_nam

    pth_main="${dir_scr}/install_envs.sh"

    if [ ! -f "${pth_main}" ]; then
        err "required handoff script does not exist: '${pth_main}'."
        exit 1
    fi

    if [ ! -r "${pth_main}" ]; then
        err "required handoff script is not readable: '${pth_main}'."
        exit 1
    fi

    if pth_bash="$(find_bash_required)"; then
        if (
               command -v conda >/dev/null 2>&1 \
            || command -v mamba >/dev/null 2>&1
        ); then
            exec "${pth_bash}" "${pth_main}" "$@"
        fi

        err \
            "Bash >= 4.4 was found in PATH, but neither 'conda' nor 'mamba'" \
            "is currently available in PATH."
        echo >&2
        print_guidance_miniforge
        exit 1
    fi

    if ! (
           command -v conda >/dev/null 2>&1 \
        || command -v mamba >/dev/null 2>&1
    ); then
        err \
            "Bash >= 4.4 was not found in PATH, and neither 'conda' nor" \
            "'mamba' is currently available in PATH."
        echo >&2
        print_guidance_miniforge
        exit 1
    fi

    err "Bash >= 4.4 was not found in PATH."
    echo >&2
    print_guidance_bash
    exit 1
}


main "$@"
