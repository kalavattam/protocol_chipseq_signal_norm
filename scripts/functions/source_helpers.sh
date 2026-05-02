#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: source_helpers.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
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
    local path=""

    if [[ -z "${helper}" ]]; then
        _source_helper_err "helper name/path is missing."
        return 1
    fi

    #  Allow either:
    #+   - aggregate helper name: check_args
    #+   - relative helper path: help/help_find_files
    #+   - explicit shell filename: check_args.sh
    #+   - absolute path: /path/to/check_args.sh
    if [[ "${helper}" == /* ]]; then
        path="${helper}"
    else
        path="${dir_fnc}/${helper}"
    fi

    if [[ "${path}" != *.sh ]]; then
        path="${path}.sh"
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
        cat >&2 <<'EOM'
Usage:
  source_once helper [dir_fnc]

Description:
  Source a Bash helper file at most once per shell process.

Positional arguments:
  1  helper   <str>  Helper name or path; examples: 'check_args', 'help/help_find_files', '/abs/path/helper.sh'.
  2  dir_fnc  <str>  Optional base function directory for relative helper names.

Returns:
  0 if the helper is already loaded or is sourced successfully; otherwise 1.

Notes:
  - Uses a global associative array registry named '__SOURCED_HELPERS'.
  - Detects circular sourcing while a file is already being loaded and returns 0 instead of re-sourcing it.
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
        cat >&2 <<'EOM'
Usage:
  source_helpers dir_fnc helper1 [helper2 ...]

Description:
  Source multiple Bash helper files at most once each.

Positional arguments:
  1   dir_fnc  <str>  Base function directory.
  2+  helper   <mlt>  Helper names or relative paths.

Example:
  source_helpers "${dir_fnc}" \
      format_outputs \
      check_args \
      check_inputs \
      help/help_find_files
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
