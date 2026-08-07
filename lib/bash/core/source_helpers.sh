#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: source_helpers.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in design, development,
# and documentation, with all output reviewed, edited, and approved by the
# author.
#
# Distributed under the MIT license.


# _source_helper_err
# _source_helper_resolve
# source_once
# source_helpers


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

#  Initialize global registry of sourced helper files
#+
#+ Key: canonical absolute file path
#+ Value:
#+   loading = currently being sourced
#+   loaded  = successfully sourced
if ! declare -p __SOURCED_HELPERS >/dev/null 2>&1; then
    declare -gA __SOURCED_HELPERS=()
fi

#  Resolve the directory containing this source helper
__DIR_SOURCE_HELPERS="$(
    cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
)"


function _source_helper_err() {
    echo "error($(basename "${0}")::${FUNCNAME[1]}): $*" >&2
}


function _source_helper_resolve() {
    local helper="${1:-}"
    local dir_fnc="${2:-${__DIR_SOURCE_HELPERS}}"
    local candidate=""
    local path=""
    local -a arr_candidates=()
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _source_helper_resolve
    [--help] helper [dir_fnc]

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  helper : str
    Absolute path, responsibility-qualified path, or unique helper basename.

  2  dir_fnc : dir
    Base Bash-library directory for relative helper names.

Returns
-------
  Returns 0 after printing help or resolving one helper; otherwise 1.

Output
------
  Writes the resolved canonical helper path to stdout. Writes help and errors to stderr so command substitutions cannot mistake documentation for a path.

Notes
-----
  Runtime requirements:
    - basename
    - bash >= 4.4
    - dirname

Examples
--------
  1. Resolve a unique helper basename from the default library root.
    '''bash
    _source_helper_resolve check_args "\${__DIR_SOURCE_HELPERS}"
    '''

  2. Resolve one responsibility-qualified helper path.
    '''bash
    _source_helper_resolve \\
        workflows/process_tables "\${__DIR_SOURCE_HELPERS}/.."
    '''
EOM
)

    if [[ "${helper}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ -z "${helper}" ]]; then
        _source_helper_err "helper name/path is missing."
        return 1
    fi

    #  Allow an absolute path, a responsibility-qualified relative path, or a
    #+ unique basename found under the Bash-library responsibility roots.
    if [[ "${helper}" == /* ]]; then
        path="${helper}"
        [[ "${path}" == *.sh ]] || path="${path}.sh"
    elif [[ "${helper}" == */* ]]; then
        path="${dir_fnc}/${helper}"
        [[ "${path}" == *.sh ]] || path="${path}.sh"
    else
        path="${dir_fnc}/${helper}"
        [[ "${path}" == *.sh ]] || path="${path}.sh"

        if [[ ! -f "${path}" ]]; then
            for candidate in core dispatch workflows help; do
                path="${dir_fnc}/${candidate}/${helper}"
                [[ "${path}" == *.sh ]] || path="${path}.sh"
                [[ -f "${path}" ]] && arr_candidates+=( "${path}" )
            done

            if (( ${#arr_candidates[@]} > 1 )); then
                _source_helper_err \
                    "helper basename is ambiguous under '${dir_fnc}':" \
                    "'${helper}'."
                return 1
            elif (( ${#arr_candidates[@]} == 1 )); then
                path="${arr_candidates[0]}"
            fi
        fi
    fi

    #  Canonicalize without requiring GNU realpath behavior
    if [[ ! -f "${path}" ]]; then
        _source_helper_err "helper file not found: '${path}'."
        return 1
    fi

    (
        cd "$(dirname "${path}")" > /dev/null 2>&1 \
            && printf '%s/%s\n' "${PWD}" "$(basename "${path}")"
    )
}


function source_once() {
    local helper="${1:-}"
    local dir_fnc="${2:-${__DIR_SOURCE_HELPERS}}"
    local path=""

    if [[ "${helper}" =~ ^(-h|--h[e]?lp)$ ]]; then
        cat >&2 << EOM
Usage
-----
  source_once
    [--help] helper [dir_fnc]

  Source a Bash helper file at most once per shell process.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  helper : str
    Helper name or path; examples: 'check_args', 'help/help_find_files', '/abs/path/helper.sh'.

  2  dir_fnc : dir
    Base function directory for relative helper names.

Returns
-------
  0 if the helper is already loaded or is sourced successfully; otherwise 1.

Notes
-----
  Runtime requirements:
    - basename
    - bash >= 4.4
    - dirname

  - Uses a global associative array registry named '__SOURCED_HELPERS'.
  - Detects circular sourcing while a file is already being loaded and returns 0 instead of re-sourcing it.

Examples
--------
  1. Load one shared validation helper by aggregate name.
    '''bash
    source_once check_args "\${__DIR_SOURCE_HELPERS}"
    '''

  2. Repeat a load request and let the registry suppress re-sourcing.
    '''bash
    source_once check_inputs "\${__DIR_SOURCE_HELPERS}"
    source_once check_inputs "\${__DIR_SOURCE_HELPERS}"
    '''
EOM
        return 0
    fi

    path="$(_source_helper_resolve "${helper}" "${dir_fnc}")" || return 1

    case "${__SOURCED_HELPERS[${path}]:-}" in
        loaded)
            return 0
            ;;

        loading)
            #  Break recursive source cycles. The first source chain continues
            #+ loading the file; the recursive request is treated as already
            #+ satisfied for import-control purposes.
            return 0
            ;;
    esac

    __SOURCED_HELPERS["${path}"]="loading"

    # shellcheck disable=SC1090
    if ! source "${path}"; then
        unset '__SOURCED_HELPERS[${path}]'
        _source_helper_err "failed to source '${path}'."
        return 1
    fi

    __SOURCED_HELPERS["${path}"]="loaded"
    return 0
}


function source_helpers() {
    local dir_fnc="${1:-}"
    local helper=""

    if [[ "${dir_fnc}" =~ ^(-h|--h[e]?lp)$ ]]; then
        cat >&2 << EOM
Usage
-----
  source_helpers
    [--help] dir_fnc helper [helper ...]

  Source multiple Bash helper files at most once each.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  dir_fnc : dir
    Base function directory.

  2+  helper : str
    Helper names or relative paths.

Returns
-------
  Returns 0 when every requested helper is sourced successfully; otherwise 1.

Notes
-----
  Runtime requirements:
    - basename
    - bash >= 4.4
    - dirname

Examples
--------
  1. Load argument and input validation helpers together.
    '''bash
    source_helpers "\${dir_fnc}" check_args check_inputs
    '''

  2. Load one nested help module relative to the function directory.
    '''bash
    source_helpers "\${dir_fnc}" help/help_find_files
    '''
EOM
        return 0
    elif [[ -z "${dir_fnc}" ]]; then
        _source_helper_err "positional argument 1, 'dir_fnc', is missing."
        return 1
    elif (( $# < 2 )); then
        _source_helper_err "at least one helper name must be supplied."
        return 1
    fi

    shift

    for helper in "$@"; do
        source_once "${helper}" "${dir_fnc}" || return 1
    done
}
