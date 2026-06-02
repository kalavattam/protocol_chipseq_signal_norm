#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: clean_test_outputs.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


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

#  Resolve paths relative to 'tests/scripts'
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_rep="$(cd "${dir_scr}/../.." > /dev/null 2>&1 && pwd)"


#  Print help text
function show_help_clean() {
    cat << EOM
Usage:
  clean_test_outputs.sh [-h|--hlp|--help]
  clean_test_outputs.sh [--dry_run] (--fixtures | --outputs | --all)

Description:
  Remove ignored generated test artifacts from explicitly scoped paths. Tracked fixture README files and other tracked files are never removed.

Arguments:
  -h, --hlp, --help  <flag>
    Display this help message and exit.

  -dr, --dry, --dry_run  <flag>
    Print the cleanup plan without removing files.

  -fx, --fix, --fixtures  <flag>
    Remove ignored generated files under tests/*/fixtures.

  -o, --out, --outputs  <flag>
    Remove ignored generated files under tests/output.

  -a, --all  <flag>
    Remove both generated fixtures and smoke-test outputs.

Dependencies:
  External programs:
    git
      Required for repository-root validation, tracked-file checks, and scoped cleanup of ignored generated artifacts using 'git clean'.

  Sourced function scripts:
    None.

Returns:
  0 on success; non-zero on error.
EOM
}


#  Print an error
function error() {
    echo "error($(basename "${BASH_SOURCE[0]}")):" "$@" >&2
}


#  Print a warning
function warn() {
    echo "warning($(basename "${BASH_SOURCE[0]}")):" "$@" >&2
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


#  Require fixture roots to contain only their tracked README files
function check_fixture_safe() {
    local dir_fix=""
    local readme=""
    local -a arr_tracked=()

    for dir_fix in "${arr_fix[@]}"; do
        readme="${dir_fix}/README.md"
        mapfile -t arr_tracked < <(
            git -C "${dir_rep}" ls-files -- "${dir_fix}"
        )

        if (( ${#arr_tracked[@]} != 1 )); then
            die \
                "fixture root '${dir_fix}' must track only '${readme}';" \
                "found ${#arr_tracked[@]} tracked files."
        elif [[ "${arr_tracked[0]}" != "${readme}" ]]; then
            die \
                "refusing fixture cleanup; unexpected tracked path:" \
                "'${arr_tracked[0]}'."
        fi
    done
}


#  Require the smoke-test output directory to contain no tracked files
function check_output_safe() {
    local -a arr_tracked=()
    local file=""

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

    arr_cmd+=( -- "${arr_tgt[@]}" )
    "${arr_cmd[@]}"
}


#  Define repository-relative cleanup targets ================================
arr_fix=(
    tests/align_fastqs/fixtures
    tests/calculate_scaling_factor/fixtures
    tests/compute_signal/fixtures
    tests/download_fastqs/fixtures
    tests/filter_alignments/fixtures
    tests/trim_fastqs/fixtures
)

dir_out="tests/output"


#  Parse arguments ============================================================
cln_fix=false
cln_out=false
dry_run=false
shw_hlp_cln=false

while [[ "$#" -gt 0 ]]; do
    case "${1}" in
        -h|--hlp|--help)
            show_help_clean >&2
            exit 0
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

        -dr|--dry|--dry[_-]run)
            dry_run=true
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
    if [[ "${dry_run}" == "true" ]]; then
        warn \
            "'--dry_run' was specified without a cleanup target; showing" \
            "what '--dry_run --all' would do."
        cln_fix=true
        cln_out=true
        shw_hlp_cln=true
    else
        error \
            "one of '--fixtures', '--outputs', or '--all' must be specified."
        echo >&2
        show_help_clean >&2
        exit 1
    fi
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
    arr_tgt+=( "${dir_out}" )
fi

print_targets
clean_paths

if [[ "${shw_hlp_cln}" == "true" ]]; then
    echo >&2
    show_help_clean >&2
fi
