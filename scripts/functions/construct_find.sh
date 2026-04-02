#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: construct_find.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
# 
# Distributed under the MIT license.


function construct_find() {
    local dir_fnd=""
    local pattern=""
    local follow=false
    local depth=""
    local include=""
    local exclude=""
    local arr_nam=""
    local show_help
    local patt
    local -a arr_inc
    local -a arr_exc

    show_help=$(cat << EOM
Description:
  Build a 'find' command as an argument array for safe downstream execution.

Usage:
  construct_find
    [--help] --arr_nam <str> --dir_fnd <str> --pattern <str> [--follow] [--depth <int>] [--include <str_csv>] [--exclude <str_csv>]

Keyword parameters:
  -an, --arr_nam  <str>  Name of the output array variable to populate.
  -df, --dir_fnd  <str>  Directory in which to search for files (required).
  -pa, --pattern  <str>  Primary filename pattern to match (required).
  -fl, --follow   <flg>  Follow symbolic links during the search.
  -de, --depth    <int>  Maximum search depth.
  -in, --include  <str>  Comma-separated list of additional filename patterns to include; each is applied with logical AND.
  -ex, --exclude  <str>  Comma-separated list of filename patterns to exclude; each is applied with logical AND NOT.

Returns:
  Populates the named array variable with a 'find' argument array.

Notes:
  - This function builds an argument array and is intended to be consumed as "\${arr_cmd[@]}" rather than by reconstructing a command string.
  - Logical OR is not implemented; '--include' patterns are applied as repeated '-name' predicates, i.e., logical AND.
EOM
    )

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -an|--arr[_-]nam)           arr_nam="${2:-}"; shift $(( $# >= 2 ? 2 : 1 )) ;;
            -df|--dir[_-]fnd)           dir_fnd="${2:-}"; shift $(( $# >= 2 ? 2 : 1 )) ;;
            -pa|--pattern)              pattern="${2:-}"; shift $(( $# >= 2 ? 2 : 1 )) ;;
            -fl|--follow|-sl|--symlink) follow=true;      shift 1 ;;
            -de|--dpth|--depth)         depth="${2:-}";   shift $(( $# >= 2 ? 2 : 1 )) ;;
            -in|--incld|--include)      include="${2:-}"; shift $(( $# >= 2 ? 2 : 1 )) ;;
            -ex|--excld|--exclude)      exclude="${2:-}"; shift $(( $# >= 2 ? 2 : 1 )) ;;
            *)
                echo "## Unknown parameter passed: '${1}' ##" >&2
                echo >&2
                echo "${show_help}" >&2
                return 1
                ;;
        esac
    done

    if [[ -z "${arr_nam}" ]]; then
        echo "Error: '--arr_nam' is required." >&2
        return 1
    elif [[ ! "${arr_nam}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ ]]; then
        echo \
            "Error: '--arr_nam' must be a valid shell variable name; got" \
            "'${arr_nam}'." >&2
        return 1
    fi

    if [[ -z "${dir_fnd}" ]]; then
        echo "Error: '--dir_fnd' is required." >&2
        return 1
    fi

    if [[ ! -d "${dir_fnd}" ]]; then
        echo \
            "Error: Directory associated with '--dir_fnd' does not exist:" \
            "'${dir_fnd}'." >&2
        return 1
    fi

    if [[ -z "${pattern}" ]]; then
        echo "Error: '--pattern' is required." >&2
        return 1
    fi

    if [[ -n "${depth}" && ! "${depth}" =~ ^[1-9][0-9]*$ ]]; then
        echo \
            "Error: '--depth' must be a positive integer >= 1; got" \
            "'${depth}'." >&2
        return 1
    fi

    eval "unset ${arr_nam}"
    eval "${arr_nam}=()"

    eval "${arr_nam}+=( find )"

    if [[ "${follow}" == "true" ]]; then
        eval "${arr_nam}+=( -L )"
    fi

    eval "${arr_nam}+=( \"\${dir_fnd}\" )"

    if [[ -n "${depth}" ]]; then
        eval "${arr_nam}+=( -maxdepth \"\${depth}\" )"
    fi

    eval "${arr_nam}+=( -type f -name \"\${pattern}\" )"

    if [[ -n "${include}" ]]; then
        IFS=',' read -r -a arr_inc <<< "${include}"
        for patt in "${arr_inc[@]}"; do
            [[ -n "${patt}" ]] || continue
            eval "${arr_nam}+=( -name \"\${patt}\" )"
        done
    fi

    if [[ -n "${exclude}" ]]; then
        IFS=',' read -r -a arr_exc <<< "${exclude}"
        for patt in "${arr_exc[@]}"; do
            [[ -n "${patt}" ]] || continue
            eval "${arr_nam}+=( ! -name \"\${patt}\" )"
        done
    fi
}
