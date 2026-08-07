#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: clean_artifacts.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in design, development,
# and documentation, with all output reviewed, edited, and approved by the
# author.
#
# Distributed under the MIT license.


#TODO: refactor for lifecycle-function formatting

#  Require Bash >= 4.4 before doing any work
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be run under Bash >= 4.4." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

#  Run in safe mode, exiting on errors, unset variables, and pipe failures
set -euo pipefail

#  Resolve paths relative to 'tests'
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_rep="$(cd "${dir_scr}/../.." > /dev/null 2>&1 && pwd)"


#  Print help text
function show_help_clean() {
    cat << EOM
Usage
-----
  clean_artifacts.sh
    [--help] [--dry_run] (--fixtures | --outputs | --all)

  Remove ignored generated test artifacts from explicitly scoped paths. Tracked fixture README files and other tracked files are never removed.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -dr, --dry, --dry_run : flag
    Run script in dry-run mode. Print the cleanup plan without removing files.

  -fx, --fix, --fixtures : flag
    Remove ignored generated files under tests/fixtures/*.

  -o, --out, --outputs : flag
    Remove ignored disposable logs, temporary work, Python/pytest caches, and package-build output. Retained checkpoint and migration evidence is never selected.

  -a, --all : flag
    Remove both generated fixtures and test outputs.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - git

  Destructive cleanup requires the exact environment value 'CONFIRM_TEST_CLEANUP=1'. True-like alternatives are rejected.

Returns
-------
  0 on success; non-zero on error.

Examples
--------
  1. Preview cleanup of ignored test outputs only.
    '''bash
    bash tests/support/clean_artifacts.sh --dry_run --outputs
    '''

  2. Preview cleanup of ignored generated fixture files only.
    '''bash
    bash tests/support/clean_artifacts.sh --dry_run --fixtures
    '''
EOM
}


#  Print an error
function error() {
    echo "error($(basename "${BASH_SOURCE[0]}")):" "$@" >&2
}


#  Print an error and stop
function die() {
    echo "error($(basename "${BASH_SOURCE[0]}")):" "$@" >&2
    exit 1
}


#  Print a cleanup note
function note() {
    echo "note($(basename "${BASH_SOURCE[0]}")):" "$@" >&2
}


#  Require execution from the expected repository checkout
function require_git_repo() {
    local dir_git=""

    if ! \
        command -v git > /dev/null 2>&1
    then
        die "'git' must be available to clean generated test artifacts."
    fi

    if ! \
        dir_git="$(git -C "${dir_rep}" rev-parse --show-toplevel 2> /dev/null)"
    then
        die "could not resolve a Git repository from '${dir_rep}'."
    fi

    if [[ "${dir_git}" != "${dir_rep}" ]]; then
        die \
            "resolved Git root '${dir_git}' does not match expected root" \
            "'${dir_rep}'."
    fi
}


#  Require fixture roots to contain only their hand-written source files
function check_fixture_safe() {
    local dir_fix=""
    local file=""
    local readme=""
    local recipe=""
    local -a arr_tracked=()

    for dir_fix in "${arr_fix[@]}"; do
        readme="${dir_fix}/README.md"
        recipe="${dir_fix}/make.sh"
        if [[ ! -f "${readme}" || ! -f "${recipe}" ]]; then
            die \
                "fixture root '${dir_fix}' must contain README.md and make.sh."
        fi
        mapfile -t arr_tracked < <(
            git -C "${dir_rep}" ls-files -- "${dir_fix}"
        )

        for file in "${arr_tracked[@]}"; do
            if [[ "${file}" != "${readme}" && "${file}" != "${recipe}" ]]; then
                die \
                    "refusing fixture cleanup; unexpected tracked path:" \
                    "'${file}'."
            fi
        done
    done
}


#  Require test output directories to contain no tracked files
function check_output_safe() {
    local -a arr_tracked=()
    local dir_out=""
    local file=""

    for dir_out in "${arr_out[@]}"; do
        mapfile -t arr_tracked < <(
            git -C "${dir_rep}" ls-files -- "${dir_out}"
        )

        if (( ${#arr_tracked[@]} > 0 )); then
            error "refusing output cleanup; tracked paths found:"
            for file in "${arr_tracked[@]}"; do
                printf '  %s\n' "${file}" >&2
            done
            exit 1
        fi
    done
}


#  Print selected cleanup targets
function print_targets() {
    local target=""

    note "selected cleanup targets:"
    for target in "${arr_tgt[@]}"; do
        printf '  %s\n' "${target}" >&2
    done
}


#  Remove or report ignored files under the selected scoped paths
function clean_paths() {
    local -a arr_cmd=(
        git -C "${dir_rep}" clean -dX
    )

    if [[ "${dry_run}" == "true" ]]; then
        arr_cmd+=( -n )
        note "dry-run mode: reporting ignored generated artifacts."
    else
        arr_cmd+=( -f )
        note "removing ignored generated artifacts."
    fi

    #  This helper is called only for fixture cleanup. Passing ignored output
    #+ children to 'git clean' can cause Git to select their ignored parent and
    #+ retained sibling evidence, so disposable outputs are handled literally
    #+ by 'clean_output_paths' instead.
    arr_cmd+=( -- "${arr_fix[@]}" )
    "${arr_cmd[@]}"
}


#  Remove or report each literal disposable output root
function clean_output_paths() {
    local dir_out=""
    local path_out=""

    if [[ "${dry_run}" == "true" ]]; then
        note "dry-run mode: reporting disposable output roots."
    else
        note "removing literal disposable output roots."
    fi

    for dir_out in "${arr_out[@]}"; do
        path_out="${dir_rep}/${dir_out}"
        if [[ ! -e "${path_out}" ]]; then
            continue
        elif [[ -L "${path_out}" || ! -d "${path_out}" ]]; then
            die "output target must be a nonsymlink directory: '${dir_out}'."
        fi

        if [[ "${dry_run}" == "true" ]]; then
            printf 'Would remove %s\n' "${dir_out}"
        else
            rm -rf -- "${path_out}"
        fi
    done
}


#  Define repository-relative cleanup targets ================================
arr_fix=(
    tests/fixtures/align_fastqs
    tests/fixtures/calculate_scaling_factor
    tests/fixtures/compute_signal
    tests/fixtures/download_fastqs
    tests/fixtures/filter_alignments
    tests/fixtures/trim_fastqs
)

arr_out=(
    .pytest_cache
    artifacts/tests/logs
    artifacts/tests/package
    artifacts/tests/pycache
    artifacts/tests/pytest_cache
    artifacts/tests/tmp
    dev/__pycache__
    dev/audit/__pycache__
    dev/tools/__pycache__
    tests/__pycache__
    tests/unit/dev_audit/__pycache__
)


#  Parse arguments ============================================================
cln_fix=false
cln_out=false
dry_run=false

while [[ "$#" -gt 0 ]]; do
    case "${1}" in
        -h|--hlp|--help)
            show_help_clean >&2
            exit 0
            ;;

        -dr|--dry|--dry[_-]run)
            dry_run=true
            shift 1
            ;;

        -fx|--fix|--fixtures)
            cln_fix=true
            shift 1
            ;;

        -o|--out|--outputs)
            cln_out=true
            shift 1
            ;;

        -a|--all)
            cln_fix=true
            cln_out=true
            shift 1
            ;;
        *)
            error "unknown option/parameter passed: '${1}'."
            echo >&2
            show_help_clean >&2
            exit 1
            ;;
    esac
done

if [[ "${cln_fix}" == "false" && "${cln_out}" == "false" ]]; then
    error \
        "one of '--fixtures', '--outputs', or '--all' must be specified."
    echo >&2
    show_help_clean >&2
    exit 1
fi


#  Validate and run scoped cleanup ============================================
require_git_repo

arr_tgt=()
if [[ "${cln_fix}" == "true" ]]; then
    check_fixture_safe
    arr_tgt+=( "${arr_fix[@]}" )
fi

if [[ "${cln_out}" == "true" ]]; then
    check_output_safe
    arr_tgt+=( "${arr_out[@]}" )
fi

print_targets
if [[ \
    "${dry_run}" == "false" && \
    "${CONFIRM_TEST_CLEANUP:-}" != "1"
]]; then
    die "destructive cleanup requires exact CONFIRM_TEST_CLEANUP=1."
fi
if [[ "${cln_fix}" == "true" ]]; then
    clean_paths
fi
if [[ "${cln_out}" == "true" ]]; then
    clean_output_paths
fi
