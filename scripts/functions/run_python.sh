#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: run_python.sh
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# _resolve_dir_rep_run_py
# to_module
# run_py


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

#  Set default Python runner mode
PY_INVOKE="${PY_INVOKE:-module}"  # "file" or "module" (default)

#  Set internal fallback repo root for module conversion and PYTHONPATH setup
_run_py_dir_rep="$(
    cd "$(dirname "${BASH_SOURCE[0]}")/../.." > /dev/null 2>&1 && pwd
)"


# shellcheck disable=SC2120
function _resolve_dir_rep_run_py() {
    local dir_rep_lcl="${dir_rep:-${_run_py_dir_rep:-}}"
    local show_help

    show_help=$(cat << EOM
Usage:
  _resolve_dir_rep_run_py [-h|--hlp|--help]

Description:
  Resolve the repository root used by 'run_python.sh'.

  If global variable 'dir_rep' is set and non-empty, that value is used.
  Otherwise, an internal default derived relative to this script is used.

Returns:
  Prints the resolved repository-root path.
EOM
    )

    if [[ "${1:-}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    fi

    if [[ -z "${dir_rep_lcl}" ]]; then
        echo \
            "error(run_python.sh::_resolve_dir_rep_run_py): failed to" \
            "resolve repository root." >&2
        return 1
    fi

    printf '%s\n' "${dir_rep_lcl}"
}


function to_module() {
    local p="${1:-}"
    local dir_rep_lcl
    local show_help

    show_help=$(cat << EOM
Usage:
  to_module [-h|--hlp|--help] p

Description:
  Convert a repo-relative or repo-anchored Python file path to a dotted module path suitable for 'python -m'.

Positional argument:
  1  p  <str>  Python file path to convert.

Returns:
  Prints the corresponding dotted module path.

Examples:
  '''bash
  to_module scripts/python/calculate_scaling_factor_spike.py
  to_module "\${dir_rep}/scripts/python/utils/helpers.py"
  '''

Notes:
  - The repository root is resolved from global variable 'dir_rep' when set; otherwise, an internal default relative to this script is used.
  - If 'p' begins with that resolved repository-root path, the prefix is removed. A trailing '.py' suffix is also removed.
EOM
    )

    if [[ "${p}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${p}" ]]; then
        echo \
            "error(run_python.sh::to_module): positional argument 1, 'p', is" \
            "missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    # shellcheck disable=SC2119
    dir_rep_lcl="$(_resolve_dir_rep_run_py)" || return 1

    p="${p#"${dir_rep_lcl%/}"/}"  # Strip repo root
    p="${p%.py}"                  # Drop '.py'

    printf '%s\n' "${p//\//.}"    # Convert '/' to '.'
}


function run_py() {
    local entry="${1:-}"
    local dir_rep_lcl
    local show_help

    show_help=$(cat << EOM
Usage:
  run_py [-h|--hlp|--help] entry [args...]

Description:
  Run a Python entry point using either module invocation ('python -m <module>') or file-path invocation ('python <file>'), depending on global variable 'PY_INVOKE'.

  If 'PY_INVOKE=module' and 'entry' looks like a file path (contains '/' or ends in '.py'), 'to_module' is used to convert it to a dotted module path before execution.

Positional arguments:
  1   entry   <str>  Python module name or file path to execute.
  2+  args    <str>  Additional arguments passed through unchanged.

Expected global variables:
  Read:
    PY_INVOKE  Invocation mode: 'module' or 'file'.
    dir_rep    Optional repository root used to help convert file paths to modules; if unset, an internal default relative to this script is used.

Returns:
  Executes the requested Python entry point and returns its exit status.

Examples:
  '''bash
  run_py scripts/python/calculate_scaling_factor_spike.py --coef fractional
  run_py scripts.python.calculate_scaling_factor_spike --coef fractional
  '''
EOM
    )

    if [[ "${entry}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${entry}" ]]; then
        echo \
            "error(run_python.sh::run_py): positional argument 1, 'entry', is" \
            "missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    shift || true

    case "${PY_INVOKE}" in
        module|file) : ;;
        *)
            echo \
                "error(run_python.sh::run_py): global variable 'PY_INVOKE'" \
                "must be 'module' or 'file': '${PY_INVOKE}'." >&2
            return 1
            ;;
    esac

    if [[ "${PY_INVOKE}" == "module" ]]; then
        local mod="${entry}"

        if [[ "${entry}" == *.py || "${entry}" == */* ]]; then
            mod="$(to_module "${entry}")" || return 1
        fi

        # shellcheck disable=SC2119
        dir_rep_lcl="$(_resolve_dir_rep_run_py)" || return 1
        PYTHONPATH="${dir_rep_lcl}:${PYTHONPATH:-}" python -m "${mod}" "$@"
    else
        python "${entry}" "$@"
    fi
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script is intended to be sourced before use; do not run it as," \
        "e.g., 'bash ${BASH_SOURCE[0]}'." >&2
    exit 1
fi
