#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: handle_env.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# activate_env
# _current_errexit_nounset
# _restore_errexit_nounset
# _handle_env_deactivate
# _handle_env_activate_success
# _handle_env_activate
# handle_env


#  Require Bash >= 4.4 before defining functions
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be sourced or run under Bash >= 4.4." >&2

    if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
        return 1
    else
        exit 1
    fi
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2

    if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
        return 1
    else
        exit 1
    fi
fi

#  Source required helper functions if needed
# shellcheck disable=SC1091
{
    _dir_src_env="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    source "${_dir_src_env}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_env}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_env}" \
        format_outputs check_source || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_env
}


#  Legacy public wrapper; prefer 'handle_env'
function activate_env() { handle_env "$@"; }
#TODO: remove once all call sites have been updated to use 'handle_env'


#  Report current 'set -e' and 'set -u' states as, respectively, 'had_e' and
#+ 'had_u' on stdout, each 0 (off) or 1 (on); does not change shell options
function _current_errexit_nounset() {
    local had_e=0
    local had_u=0

    case "$-" in
        *e*) had_e=1 ;;
    esac

    case "$-" in
        *u*) had_u=1 ;;
    esac

    printf '%s %s\n' "${had_e}" "${had_u}"
}


#  Restore prior 'set -e', 'set -u' states based on flags '0' (off), '1' (on)
function _restore_errexit_nounset() {
    local had_e="${1:-0}"
    local had_u="${2:-0}"

    if [[ "${had_e}" -eq 1 ]]; then set -e; fi
    if [[ "${had_u}" -eq 1 ]]; then set -u; fi
}


function _handle_env_deactivate() {
    local shl_cur="${1:-}"
    local cur_env="${CONDA_DEFAULT_ENV:-}"
    local had_e had_u
    local rc=0
    local show_help

    show_help=$(cat << EOM
Usage:
  _handle_env_deactivate [-h|--hlp|--help] [shl_cur]

Description:
  Safely deactivate a Conda/Mamba environment.

  If deactivation is needed, temporarily disables 'set -e' and 'set -u' so Conda shell-hook and deactivation code can run without aborting on non-zero check steps or unset-variable references. Restores the prior shell-option state before returning.

Positional argument:
  1  shl_cur  <str>  Optional shell-name hint for Conda hook initialization. This helper is intended for Bash >= 5 workflows. If supplied, 'bash' is the primary supported case and 'zsh' is handled only as a limited fallback for hook setup.

Returns:
  0 if no deactivation is needed or if deactivation succeeds; otherwise, 1.

Dependencies:
  - Bash >= 5
  - Conda

#TODO:
  - Add usage example(s).
  - If 'CONDA_DEFAULT_ENV' is set but 'conda' is somehow unavailable in 'PATH', the function currently restores flags and returns success unless hook init failed. That is a kind of subtle edge case; maybe it’s acceptable in practice for the workflow?
EOM
    )

    #  Parse help request
    if [[ "${shl_cur}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    #  If no env is active, or only 'base' is active, nothing to do
    if [[ -z "${cur_env}" || "${cur_env}" == "base" ]]; then
        return 0
    fi

    read -r had_e had_u < <(_current_errexit_nounset)
    set +e +u

    if [[ -z "${shl_cur}" ]]; then
        shl_cur=$(basename "${SHELL:-/bin/bash}")
    fi

    case "${shl_cur}" in
        bash)
            # shellcheck disable=SC1091
            source "$(conda info --base)/etc/profile.d/conda.sh" || rc=$?
            ;;

        zsh)
            #  Limited fallback branch for Conda hook initialization; this
            #+ helper is otherwise intended for Bash >= 5 workflows
            eval "$(conda shell.zsh hook)" || rc=$?
            ;;

        *)
            echo_warn_func "${FUNCNAME[0]}" \
                "unsupported shell '${shl_cur}' for Conda hook" \
                "initialization; only 'bash' and 'zsh' are handled" \
                "explicitly. Proceeding to try 'conda deactivate' anyway."
            ;;
    esac

    #  If hook-sourcing failed, bail still with 'set -e' / 'set -u' off, but
    #+ restore flags before returning
    if [[ "${rc}" -ne 0 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to initialize Conda shell hook."
        _restore_errexit_nounset "${had_e}" "${had_u}"
        return 1
    fi

    if command -v conda &> /dev/null; then
        if ! conda deactivate; then
            echo_err_func "${FUNCNAME[0]}" \
                "failed to deactivate environment using conda."
            rc=1
        fi
    fi

    _restore_errexit_nounset "${had_e}" "${had_u}"
    return "${rc}"
}


#  Print activation success message and restore prior strict-mode flags
function _handle_env_activate_success() {
    local env_nam="${1:-}"
    local had_e="${2:-0}"
    local had_u="${3:-0}"

    echo "success($(basename "${0}")::_handle_env_activate):" \
        "environment '${env_nam}' activated successfully." >&2
    _restore_errexit_nounset "${had_e}" "${had_u}"
}


function _handle_env_activate() {
    local env_nam="${1:-}"
    local shl_cur show_help
    local had_e had_u
    local has_env=false
    local rc=0

    show_help=$(cat << EOM
Usage:
  _handle_env_activate [-h|--hlp|--help] env_nam

Description:
  Safely activate a specified Conda or Mamba environment.

  In Bash >= 5 workflows, this helper first initializes the appropriate Conda shell hook for the current shell context, then attempts environment activation, preferring 'mamba' when available and falling back to 'conda' and finally 'source activate'.

  Temporarily disables 'set -e' and 'set -u' while sourcing Conda shell hooks and attempting activation, as those shell-initialization commands may reference unset variables or return non-zero statuses during checking/fallback steps. The original 'set' state is restored before returning.

Positional argument:
  1  env_nam  <str>  The name of the environment to activate (required).

Returns:
  0 if the requested environment is activated successfully; otherwise 1.

Dependencies:
  - Bash >= 5
  - Conda
  - Mamba (optional, preferred)

Examples:
  1. Activate an extant environment named "my_env"
  '''bash
  _handle_env_activate "my_env"
  '''

  '''txt
  Environment 'my_env' activated successfully.
  '''

  2. Try to activate an environment that does not exist
  '''bash
  _handle_env_activate "monkey"
  '''

  '''txt
  Error: Environment 'monkey' does not exist.
  '''
EOM
    )

    #  Parse and check function arguments
    if [[ "${env_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${env_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'env_nam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Conda/Mamba shell-hook and activation logic may not be strict-mode-safe,
    #+ so capture current '-e' / '-u' state and temporarily disable both
    read -r had_e had_u < <(_current_errexit_nounset)
    set +e +u

    #  Source the appropriate Conda shell hook
    shl_cur=$(basename "${SHELL:-/bin/bash}")
    case "${shl_cur}" in
        bash)
            # shellcheck disable=SC1091
            source "$(conda info --base)/etc/profile.d/conda.sh" || rc=$?
            ;;

        zsh)
            #  Limited fallback branch for Conda hook initialization; this
            #+ helper is otherwise intended for Bash >= 5 workflows
            eval "$(conda shell.zsh hook)" || rc=$?
            ;;

        *)
            echo_warn_func "${FUNCNAME[0]}" \
                "unsupported shell '${shl_cur}'; environment activation may" \
                "not work as expected."
            ;;
    esac

    #  If sourcing the hook failed, bail still with '-e' / '-u' off, but
    #+ ensure flags are restored before returning
    if [[ "${rc}" -ne 0 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to initialize Conda shell hook."
        _restore_errexit_nounset "${had_e}" "${had_u}"
        return 1
    fi

    #  Optional: check that the env exists; now safe because '-u' / '-e' are
    #+ off
    if command -v conda &> /dev/null; then
        if conda env list | grep -Fqw -- "${env_nam}"; then
            has_env=true
        fi
    fi

    if [[ "${has_env}" == false ]] && command -v mamba &> /dev/null; then
        if mamba env list | grep -Fqw -- "${env_nam}"; then
            has_env=true
        fi
    fi

    if [[ "${has_env}" == false ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "environment '${env_nam}' does not exist."
        _restore_errexit_nounset "${had_e}" "${had_u}"
        return 1
    fi

    #  Try activation with 'mamba', then 'conda', then 'source activate'
    if command -v mamba &> /dev/null; then
        #TODO: when Mamba is installed/working, 'mamba activate' can still be
        #+     unreliable depending on shell-hook/setup details; revisit this
        if mamba activate "${env_nam}" &> /dev/null; then
            _handle_env_activate_success \
                "${env_nam}" "${had_e}" "${had_u}"
            return 0
        fi
    fi

    if command -v conda &> /dev/null; then
        if conda activate "${env_nam}" &> /dev/null; then
            _handle_env_activate_success \
                "${env_nam}" "${had_e}" "${had_u}"
            return 0
        fi
    fi

    # shellcheck disable=SC1091
    if source activate "${env_nam}" &> /dev/null; then
        _handle_env_activate_success \
            "${env_nam}" "${had_e}" "${had_u}"
        return 0
    fi

    echo_err_func "${FUNCNAME[0]}" \
        "failed to activate environment '${env_nam}'."
    _restore_errexit_nounset "${had_e}" "${had_u}"
    return 1
}


function handle_env() {
    local env_nam="${1:-}"  # Name of environment to activate
    local cur_env="${CONDA_DEFAULT_ENV:-}"
    local show_help

    show_help=$(cat << EOM
Usage:
  handle_env [-h|--hlp|--help] env_nam

Description:
  Ensure that the requested Conda/Mamba environment is active.

  If no environment is currently active, or if only 'base' is active, activates 'env_nam'. If a different environment is active, first deactivates it and then activates 'env_nam'. If the requested environment is already active, nothing changes.

  Delegates environment switching to internal helper functions that temporarily relax 'set -e' and 'set -u' around Conda/Mamba shell-hook logic, then restore the prior shell-option state before returning.

Positional argument:
  1  env_nam  <str>  Name of environment to ensure is active (required).

Returns:
  0 if the requested environment is already active or is activated successfully; otherwise, 1.

Dependencies:
  - Bash >= 5
  - Conda
  - Mamba (optional; tried first in '_handle_env_activate' when available)

Example:
  '''bash
  handle_env "env_protocol"
  '''
EOM
    )

    #  Parse and check function arguments
    if [[ "${env_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${env_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'env_nam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  If no env or env is base, activate requested one
    if [[ -z "${cur_env}" || "${cur_env}" == "base" ]]; then
        _handle_env_activate "${env_nam}" || return 1
    elif [[ "${cur_env}" != "${env_nam}" ]]; then
        #  If different env is active, deactivate then activate
        _handle_env_deactivate || return 1
        _handle_env_activate "${env_nam}" || return 1
    fi

    #  Otherwise, requested environment is already active, so do nothing
    return 0
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
