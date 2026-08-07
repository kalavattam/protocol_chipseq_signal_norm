#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: run_python.sh
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6) were
# used in design, development, and documentation, with all output reviewed,
# edited, and approved by the author.
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

#  Use installed console commands locally; bundled Slurm runs opt into module
#+ dispatch with an explicit Python interpreter and bundled source root.
PY_INVOKE="${PY_INVOKE:-console}"
PROTOCOL_PYTHON="${PROTOCOL_PYTHON:-python}"

#  Set internal fallback repo root for module conversion and PYTHONPATH setup
_run_py_dir_rep="$(
    cd "$(dirname "${BASH_SOURCE[0]}")/../../.." > /dev/null 2>&1 && pwd
)"


# shellcheck disable=SC2120
function _resolve_dir_rep_run_py() {
    local dir_rep_lcl="${dir_rep:-${_run_py_dir_rep:-}}"
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _resolve_dir_rep_run_py
    [--help]

  Resolve the repository root used by 'run_python.sh'.

  If global variable 'dir_rep' is set and non-empty, that value is used.
  Otherwise, an internal default derived relative to this script is used.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

Returns
-------
  Prints the resolved repository-root path.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - dirname

Examples
--------
  1. Print the repository root derived by 'run_python.sh'.
    '''bash
    _resolve_dir_rep_run_py
    '''

  2. Display the helper's contract without resolving a path.
    '''bash
    _resolve_dir_rep_run_py --help
    '''
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
Usage
-----
  to_module
    [--help] p

  Convert a console-command stem or package source path to its explicit module path.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  p : str
    Console-command stem or Python package source path to convert.

Expected globals
----------------
  dir_rep : dir
    Optional repository root used when removing a repository prefix; when unset, the function derives its default relative to this script.

Returns
-------
  Prints the corresponding dotted module path.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - dirname

  - The repository root is resolved from global variable 'dir_rep' when set; otherwise, an internal default relative to this script is used.
  - If 'p' begins with that resolved repository-root path, the prefix is removed. A trailing '.py' suffix is also removed.

Examples
--------
  1. Convert a repository-relative package source path.
    '''bash
    to_module src/protocol_chipseq_signal_norm/cli/compute_signal.py
    '''

  2. Convert an installed console-command stem.
    '''bash
    to_module compute_signal
    '''
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

    p="${p#"${dir_rep_lcl%/}"/}"
    p="${p#src/}"
    p="${p%.py}"
    if [[ "${p}" != */* && "${p}" != *.* ]]; then
        p="protocol_chipseq_signal_norm/cli/${p}"
    fi

    printf '%s\n' "${p//\//.}"
}


function run_py() {
    local entry="${1:-}"
    local dir_rep_lcl
    local show_help

    show_help=$(cat << EOM
Usage
-----
  run_py
    [--help] entry [args...]

  Run an installed console command, or invoke the same CLI from a bundled package source tree when 'PY_INVOKE=module'.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1   entry : structured string
    Installed console-command stem or package CLI module/path.

  2+  args : str
    Additional arguments passed through unchanged.

Expected globals
----------------
  PY_INVOKE : {'console', 'module'}
    Python invocation mode. The default 'console' resolves the installed command. Use 'module' only for the provenance-isolated Slurm bundle.

  PROTOCOL_PYTHON : file
    Python interpreter used in module mode (default: 'python').

  dir_rep : dir
    Optional repository root used when converting file paths to module names; when unset, the function derives its default relative to this script.

Returns
-------
  Executes the requested Python entry point and returns its exit status.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - dirname
    - python >= 3.11

Examples
--------
  1. Invoke the installed compute-signal console command.
    '''bash
    PY_INVOKE=console
    run_py compute_signal --help
    '''

  2. Invoke the bundled package with an explicit Python interpreter.
    '''bash
    if PY_INVOKE=module PROTOCOL_PYTHON=/opt/env/bin/python3 \
        run_py compute_pseudo --help
    then
        printf '%s\n' 'bundled module resolved'
    fi
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
        console|module) : ;;
        *)
            echo \
                "error(run_python.sh::run_py): global variable 'PY_INVOKE'" \
                "must be 'console' or 'module': '${PY_INVOKE}'." >&2
            return 1
            ;;
    esac

    if [[ "${PY_INVOKE}" == "console" ]]; then
        "${entry}" "$@"
    else
        local mod
        mod="$(to_module "${entry}")" || return 1
        # shellcheck disable=SC2119
        dir_rep_lcl="$(_resolve_dir_rep_run_py)" || return 1
        PYTHONPATH="${dir_rep_lcl}/src${PYTHONPATH:+:${PYTHONPATH}}" \
            "${PROTOCOL_PYTHON}" -m "${mod}" "$@"
    fi
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script is intended to be sourced before use; do not run it as," \
        "e.g., 'bash ${BASH_SOURCE[0]}'." >&2
    exit 1
fi
