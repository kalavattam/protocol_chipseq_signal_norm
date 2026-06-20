#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: fixture_helpers.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development.
#
# Distributed under the MIT license.


#  Print a fixture-generator error and stop
function die() {
    local src="${BASH_SOURCE[1]:-${0}}"

    echo "error($(basename "${src}")):" "$@" >&2
    exit 1
}


#  Print a fixture-generator note
function note() {
    local src="${BASH_SOURCE[1]:-${0}}"

    echo "note($(basename "${src}")):" "$@" >&2
}


#  Print a fixture-generator success message
function succeed() {
    local src="${BASH_SOURCE[1]:-${0}}"

    echo "success($(basename "${src}")):" "$@" >&2
}


#  Require one command to be available
function require_cmd() {
    local cmd="${1:-}"
    shift || true

    if [[ -z "${cmd}" ]]; then
        die "'require_cmd' requires a command name."
    fi

    if ! \
        command -v "${cmd}" > /dev/null 2>&1
    then
        die "'${cmd}' must be available" "$@"
    fi
}


#  Require one or more commands to be available
function require_cmds() {
    local context="${1:-}"
    shift || true

    local cmd=""

    if (( $# == 0 )); then
        die "'require_cmds' requires at least one command name."
    fi

    for cmd in "$@"; do
        require_cmd "${cmd}" "${context}"
    done
}


#  Require a named Conda/Mamba environment to be active
function require_env() {
    local env_nam="${1:-}"
    shift || true

    if [[ -z "${env_nam}" ]]; then
        die "'require_env' requires an environment name."
    fi

    if [[ "${CONDA_DEFAULT_ENV:-}" != "${env_nam}" ]]; then
        die \
            "activate '${env_nam}' before generating fixtures;" \
            "current environment: '${CONDA_DEFAULT_ENV:-none}'" \
            "$@"
    fi
}


#  Create one or more directories
function mkdirs() {
    if (( $# == 0 )); then
        die "'mkdirs' requires at least one directory path."
    fi

    mkdir -p "$@"
}


#  Remove a generated fixture file inside a fixture root
function rm_file() {
    local dir_fix="${1:-}"
    local file="${2:-}"

    if [[ -z "${dir_fix}" ]]; then
        die "refusing to remove a file without a fixture root."
    elif [[ -z "${file}" ]]; then
        die "refusing to remove an empty file path."
    elif [[ "${file}" != "${dir_fix}/"* ]]; then
        die "refusing to remove path outside fixture directory: '${file}'."
    elif [[ -d "${file}" ]]; then
        die "refusing to remove directory: '${file}'."
    fi

    rm -f -- "${file}"
}


#  Remove generated fixture files inside a fixture root
function rm_files() {
    local dir_fix="${1:-}"
    shift || true

    local file=""

    if (( $# == 0 )); then
        die "rm_files requires at least one file path."
    fi

    for file in "$@"; do
        rm_file "${dir_fix}" "${file}"
    done
}


#  Register a cleanup function for shell exit
function register_cleanup() {
    local fnc="${1:-}"

    if [[ -z "${fnc}" ]]; then
        die "register_cleanup requires a function name."
    elif ! \
        declare -F "${fnc}" > /dev/null
    then
        die "cleanup function is not defined: '${fnc}'."
    fi

    # shellcheck disable=SC2064
    trap "${fnc}" EXIT
}


#  Write one tab-delimited line/row
function write_line() {
    local field

    if (( $# == 0 )); then
        die \
            "'write_line', 'write_tsv_row', and 'write_sam_line' require at" \
            "least one field."
    fi

    printf '%s' "${1}"
    shift

    for field in "${@}"; do
        printf '\t%s' "${field}"
    done

    printf '\n'
}


#  Write one tab-delimited TSV row (alias for 'write_line')
function write_tsv_row() { write_line "${@}"; }


#  Write one tab-delimited SAM line (alias for 'write_line')
function write_sam_line() { write_line "${@}"; }


#  Write gzip output with deterministic headers
function gzip_n() {
    local fil_in="${1:-}"
    local fil_out="${2:-}"

    if [[ -z "${fil_in}" ]]; then
        die "'gzip_n' requires an input path."
    elif [[ -z "${fil_out}" ]]; then
        die "'gzip_n' requires an output path."
    fi

    gzip -n -c "${fil_in}" > "${fil_out}"
}
