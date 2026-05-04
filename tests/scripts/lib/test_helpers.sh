#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_helpers.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


#  Require Bash >= 4.4 before doing any work
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(test_helpers.sh):" \
        "this test suite requires Bash >= 4.4." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error(test_helpers.sh):" \
        "this test suite requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi


#  Set paths used by smoke-test scripts =======================================
# shellcheck disable=SC2034
{
    TEST_DIR_LIB="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
    TEST_DIR_SCR="$(cd "${TEST_DIR_LIB}/.." > /dev/null 2>&1 && pwd)"
    TEST_DIR="$(cd "${TEST_DIR_SCR}/.." > /dev/null 2>&1 && pwd)"
    ROOT_REPO="$(cd "${TEST_DIR}/.." > /dev/null 2>&1 && pwd)"
    TEST_DIR_OUT="${ROOT_REPO}/tests/output"
    TEST_DIR_TMP="${TEST_DIR_OUT}/tmp"
    TEST_DIR_LOG="${TEST_DIR_OUT}/logs"
    TEST_BASH="${BASH}"
}

mkdir -p "${TEST_DIR_TMP}" "${TEST_DIR_LOG}"


#  Initialize result counters =================================================
TEST_PASS=0
TEST_FAIL=0
TEST_WARN=0
TEST_SKIP=0
TEST_NAME="${TEST_NAME:-$(basename "${0}")}"


#  Print a repository-relative path when possible
function rec_relpath() {
    local path="${1:-}"

    if [[ "${path}" == "${ROOT_REPO}"/* ]]; then
        printf '%s\n' "${path#"${ROOT_REPO}/"}"
    else
        printf '%s\n' "${path}"
    fi
}


#  Print a smoke-test rec_section heading
function rec_section() {
    printf '\n-- %s --\n' "${1:-${TEST_NAME}}"
}


#  Record a passing assertion
function rec_pass() {
    TEST_PASS=$(( TEST_PASS + 1 ))
    printf 'PASS: %s\n' "$*"
}


#  Record a failing assertion
function rec_fail() {
    TEST_FAIL=$(( TEST_FAIL + 1 ))
    printf 'FAIL: %s\n' "$*" >&2
}


#  Record a non-fatal warning
function rec_warn() {
    TEST_WARN=$(( TEST_WARN + 1 ))
    printf 'WARN: %s\n' "$*" >&2
}


#  Record a skipped check
function rec_skip() {
    TEST_SKIP=$(( TEST_SKIP + 1 ))
    printf 'SKIP: %s\n' "$*"
}


#  Check whether a command is available
function check_cmd_exists() {
    command -v "${1:-}" > /dev/null 2>&1
}


#  Find a usable Python command
function find_python() {
    if check_cmd_exists python; then
        command -v python
    elif check_cmd_exists python3; then
        command -v python3
    else
        return 1
    fi
}


#  Check whether Python is at least version 3.10
function check_python_ge_310() {
    local py="${1:-}"

    "${py}" - << PY
import sys
raise SystemExit(0 if sys.version_info >= (3, 10) else 1)
PY
}


#  Run a command and capture stdout/stderr in a log file
function run_capture() {
    local name="${1:-command}"
    local outfile="${2:-}"

    shift 2

    if [[ -z "${outfile}" ]]; then
        outfile="${TEST_DIR_LOG}/${name//[^A-Za-z0-9_.-]/_}.log"
    fi

    mkdir -p "$(dirname "${outfile}")"
    "$@" > "${outfile}" 2>&1
}


#  Assert that a log file contains a pattern
function assert_grep_pattern() {
    local file="${1:-}"
    local pattern="${2:-}"
    local label="${3:-${pattern}}"

    if \
        grep -q -- "${pattern}" "${file}"
    then
        rec_pass "${label}"
    else
        rec_fail "${label}; see $(rec_relpath "${file}")"
    fi
}


#  Warn when a help-output rec_section is missing
function warn_grep_help() {
    local file="${1:-}"
    local pattern="${2:-}"
    local label="${3:-${pattern}}"

    if \
        grep -q -- "${pattern}" "${file}"
    then
        rec_pass "${label}"
    else
        rec_warn "${label}; see $(rec_relpath "${file}")"
    fi
}


#  Print and return the final status for one smoke-test group
function finish() {
    printf \
        '\nSummary for %s: rec_pass=%d fail=%d warn=%d skip=%d\n' \
        "${TEST_NAME}" "${TEST_PASS}" "${TEST_FAIL}" "${TEST_WARN}" \
        "${TEST_SKIP}"

    if (( TEST_FAIL > 0 )); then
        return 1
    fi

    return 0
}
