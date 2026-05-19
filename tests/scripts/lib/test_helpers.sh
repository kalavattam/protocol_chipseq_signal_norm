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


#  Print a smoke-test section heading
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


#  Check whether GNU Parallel smoke tests were explicitly requested
function is_parallel_enabled() {
    [[ "${RUN_PARALLEL:-0}" == "1" ]]
}


#  Check whether Atria smoke tests were explicitly requested
function is_atria_enabled() {
    [[ "${RUN_ATRIA:-0}" == "1" ]]
}


#  Require Atria and compression helpers in the requested project environment
function require_atria_env() {
    local env_nam="${1:-env_protocol}"
    local log_lcl="${2:-${TEST_DIR_LOG}/atria_project_env.log}"
    local rc=0

    mkdir -p "$(dirname "${log_lcl}")"

    {
        echo "Current shell:"
        echo "SHELL=${SHELL:-UNSET}"
        echo "BASH=${BASH:-UNSET}"
        echo "PATH=${PATH:-UNSET}"
        echo "CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-UNSET}"
        echo "CONDA_PREFIX=${CONDA_PREFIX:-UNSET}"
        echo

        echo "Current shell command checks:"
        for cmd in atria pigz pbzip2 gzip; do
            printf '%s: ' "${cmd}"
            command -v "${cmd}" || true
        done
        echo

        if [[ "${CONDA_DEFAULT_ENV:-}" == "${env_nam}" ]]; then
            echo "Project-env command checks:"
            for cmd in atria pigz pbzip2 gzip; do
                command -v "${cmd}"
            done
        elif check_cmd_exists conda; then
            echo "Project-env command checks via Conda activation:"
            ENV_NAM="${env_nam}" bash -lc '
                eval "$(conda shell.bash hook)"
                conda activate "${ENV_NAM}"
                echo "CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-UNSET}"
                echo "CONDA_PREFIX=${CONDA_PREFIX:-UNSET}"
                echo "PATH=${PATH:-UNSET}"
                echo
                for cmd in atria pigz pbzip2 gzip; do
                    command -v "${cmd}"
                done
            '
        else
            echo "conda is unavailable; cannot inspect '${env_nam}'."
            rc=1
        fi
    } > "${log_lcl}" 2>&1 || rc=1

    if (( rc == 0 )); then
        rec_pass "Atria trim dependencies are available in project environment"
        return 0
    fi

    rec_fail \
        "Atria trim dependencies unavailable in project environment; see" \
        "$(rec_relpath "${log_lcl}")"
    return 1
}


#  Require GNU Parallel in the requested project environment
function require_parallel_env() {
    local env_nam="${1:-env_protocol}"
    local log_lcl="${2:-${TEST_DIR_LOG}/parallel_project_env.log}"
    local rc=0

    mkdir -p "$(dirname "${log_lcl}")"

    {
        echo "Current shell:"
        echo "SHELL=${SHELL:-UNSET}"
        echo "BASH=${BASH:-UNSET}"
        echo "PATH=${PATH:-UNSET}"
        echo "CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-UNSET}"
        echo "CONDA_PREFIX=${CONDA_PREFIX:-UNSET}"
        echo "MAMBA_EXE=${MAMBA_EXE:-UNSET}"
        echo

        echo "Current shell command -v parallel:"
        command -v parallel || true
        echo

        if [[ "${CONDA_DEFAULT_ENV:-}" == "${env_nam}" ]]; then
            echo "Project-env command -v parallel:"
            command -v parallel
        elif check_cmd_exists conda; then
            echo "Project-env command -v parallel via Conda activation:"
            ENV_NAM="${env_nam}" bash -lc '
                eval "$(conda shell.bash hook)"
                conda activate "${ENV_NAM}"
                echo "CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-UNSET}"
                echo "CONDA_PREFIX=${CONDA_PREFIX:-UNSET}"
                echo "PATH=${PATH:-UNSET}"
                echo
                command -v parallel
            '
        else
            echo "conda is unavailable; cannot inspect '${env_nam}'."
            rc=1
        fi
    } > "${log_lcl}" 2>&1 || rc=1

    if (( rc == 0 )); then
        rec_pass "GNU Parallel is available in project environment"
        return 0
    fi

    rec_fail \
        "GNU Parallel unavailable in project environment; see" \
        "$(rec_relpath "${log_lcl}")"
    return 1
}


#  Require wget and gzip in the requested project environment
function require_download_env() {
    local env_nam="${1:-env_protocol}"
    local log_lcl="${2:-${TEST_DIR_LOG}/download_fastqs_project_env.log}"
    local rc=0

    mkdir -p "$(dirname "${log_lcl}")"

    {
        echo "Current shell:"
        echo "SHELL=${SHELL:-UNSET}"
        echo "BASH=${BASH:-UNSET}"
        echo "PATH=${PATH:-UNSET}"
        echo "CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-UNSET}"
        echo "CONDA_PREFIX=${CONDA_PREFIX:-UNSET}"
        echo

        echo "Current shell command checks:"
        for cmd in wget gzip; do
            printf '%s: ' "${cmd}"
            command -v "${cmd}" || true
        done
        echo

        if [[ "${CONDA_DEFAULT_ENV:-}" == "${env_nam}" ]]; then
            echo "Project-env command checks:"
            for cmd in wget gzip; do
                command -v "${cmd}"
            done
        elif check_cmd_exists conda; then
            echo "Project-env command checks via Conda activation:"
            ENV_NAM="${env_nam}" bash -lc '
                eval "$(conda shell.bash hook)"
                conda activate "${ENV_NAM}"
                echo "CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-UNSET}"
                echo "CONDA_PREFIX=${CONDA_PREFIX:-UNSET}"
                echo "PATH=${PATH:-UNSET}"
                echo
                for cmd in wget gzip; do
                    command -v "${cmd}"
                done
            '
        else
            echo "conda is unavailable; cannot inspect '${env_nam}'."
            rc=1
        fi
    } > "${log_lcl}" 2>&1 || rc=1

    if (( rc == 0 )); then
        rec_pass "download dependencies are available in project environment"
        return 0
    fi

    rec_fail \
        "download dependencies unavailable in project environment; see" \
        "$(rec_relpath "${log_lcl}")"
    return 1
}


#  Resolve the active project environment or require a named fallback
function require_env_project() {
    local env_ref="${1:-env_nam}"
    local env_fallback="${2:-env_protocol}"

    if [[ -n "${CONDA_DEFAULT_ENV:-}" && "${CONDA_DEFAULT_ENV}" != "base" ]]; then
        printf -v "${env_ref}" '%s' "${CONDA_DEFAULT_ENV}"
        return 0
    fi

    if ! \
        check_cmd_exists conda
    then
        rec_fail \
            "setup requires active project env or conda to activate" \
            "'${env_fallback}'"
        return 1
    fi

    if ! \
        conda env list 2> /dev/null \
            | awk -v env="${env_fallback}" '
                $1 == env { found = 1 } END { exit !found }
            '
    then
        rec_fail "setup requires Conda/Mamba environment '${env_fallback}'"
        return 1
    fi

    printf -v "${env_ref}" '%s' "${env_fallback}"
}


#  Require one or more fixture files to exist and be non-empty
function require_files_exist() {
    local file=""
    local rc=0

    for file in "$@"; do
        if [[ ! -s "${file}" ]]; then
            rec_fail "missing fixture $(rec_relpath "${file}")"
            rc=1
        fi
    done

    return "${rc}"
}


#  Assert one predictable output file exists, without requiring non-emptiness
function assert_file_exists() {
    local file="${1:-}"
    local label="${2:-file exists}"

    if [[ -f "${file}" ]]; then
        rec_pass "${label}"
    else
        rec_fail "${label}; missing $(rec_relpath "${file}")"
    fi
}


#  Assert exactly one path was found by a caller-specific search
function assert_one_path_found() {
    local arr_ref="${1:-}"
    local label="${2:-path}"
    local dir_search="${3:-}"
    local out_ref="${4:-}"
    local path_found=""

    local -n arr_lcl="${arr_ref}"

    if (( ${#arr_lcl[@]} == 1 )); then
        path_found="${arr_lcl[0]}"
        printf -v "${out_ref}" '%s' "${path_found}"
        rec_pass "one ${label} found"
    else
        printf -v "${out_ref}" ''
        rec_fail \
            "expected exactly one ${label} in" \
            "$(rec_relpath "${dir_search}"), found ${#arr_lcl[@]}"
    fi
}


#  Run samtools from the active shell or the resolved project environment
function run_samtools() {
    local env_lcl="${env_nam:-env_protocol}"

    if check_cmd_exists samtools; then
        samtools "$@"
    else
        conda run -n "${env_lcl}" samtools "$@"
    fi
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


#  Find an available loopback TCP port for local HTTP smoke tests
function find_free_port() {
    local py="${1:-}"

    "${py}" - << PY
import socket
with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
    sock.bind(("127.0.0.1", 0))
    print(sock.getsockname()[1])
PY
}


#  Wait until a local HTTP server is reachable
function wait_for_local_http() {
    local py="${1:-}"
    local url="${2:-}"
    local tries="${3:-50}"
    local i=0

    for (( i = 1; i <= tries; i++ )); do
        if \
            URL_LCL="${url}" "${py}" - << PY
import os
import urllib.request

url = os.environ["URL_LCL"]
try:
    with urllib.request.urlopen(url, timeout=0.5) as response:
        raise SystemExit(0 if response.status == 200 else 1)
except Exception:
    raise SystemExit(1)
PY
        then
            return 0
        fi

        sleep 0.1
    done

    return 1
}


#  Stop a local HTTP server by PID
function cleanup_http_server() {
    local pid="${1:-}"

    if [[ -n "${pid}" ]] && kill -0 "${pid}" > /dev/null 2>&1; then
        kill "${pid}" > /dev/null 2>&1 || true
        wait "${pid}" > /dev/null 2>&1 || true
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


#  Assert that a log or output file does not contain a pattern
function assert_no_grep_pattern() {
    local file="${1:-}"
    local pattern="${2:-}"
    local label="${3:-${pattern}}"

    if \
        grep -q -- "${pattern}" "${file}"
    then
        rec_fail "${label}; see $(rec_relpath "${file}")"
    else
        rec_pass "${label}"
    fi
}


#  Assert that a generated output file exists and is non-empty
function assert_file_nonempty() {
    local file="${1:-}"
    local label="${2:-output}"

    if [[ -s "${file}" ]]; then
        rec_pass "${label} exists and is non-empty"
    else
        rec_fail "${label} missing or empty: $(rec_relpath "${file}")"
    fi
}


#  Assert that a CRAM index exists using common Samtools naming variants
function assert_cram_index() {
    local cram="${1:-}"
    local label="${2:-CRAM index}"
    local idx_1="${cram}.crai"
    local idx_2="${cram%.cram}.crai"

    if [[ -s "${idx_1}" ]]; then
        rec_pass "${label} exists and is non-empty at $(rec_relpath "${idx_1}")"
    elif [[ -s "${idx_2}" ]]; then
        rec_pass "${label} exists and is non-empty at $(rec_relpath "${idx_2}")"
    else
        rec_fail \
            "${label} missing or empty; checked $(rec_relpath "${idx_1}") and" \
            "$(rec_relpath "${idx_2}")"
    fi
}


#  Assert a reference-backed CRAM read count for one contig
function assert_cram_count() {
    local cram="${1:-}"
    local ref_fa="${2:-}"
    local contig="${3:-}"
    local expected="${4:-}"
    local label="${5:-CRAM count}"
    local count_file="${6:-}"

    run_capture \
        "${label}" "${count_file}" \
        run_samtools view -T "${ref_fa}" -c "${cram}" "${contig}"

    assert_grep_pattern "${count_file}" "^${expected}$" "${label}"
}


#  Assert a downloaded gzip FASTQ matches its source and expected content
function assert_downloaded_fastq() {
    local source_fastq="${1:-}"
    local outfile="${2:-}"
    local label="${3:-downloaded FASTQ}"
    local read_pattern="${4:-}"
    local view_fastq="${5:-}"
    local count_reads="${6:-}"

    assert_file_nonempty "${outfile}" "${label}"

    if [[ -s "${outfile}" ]]; then
        if cmp -s "${source_fastq}" "${outfile}"; then
            rec_pass "${label} matches source fixture byte-for-byte"
        else
            rec_fail "${label} differs from source fixture"
        fi

        if gzip -t "${outfile}"; then
            rec_pass "${label} passes gzip integrity"
        else
            rec_fail "${label} fails gzip integrity"
        fi

        if gzip -cd "${outfile}" > "${view_fastq}"; then
            rec_pass "${label} can be decompressed"
        else
            rec_fail "${label} cannot be decompressed"
        fi

        if [[ -s "${view_fastq}" ]]; then
            assert_grep_pattern "${view_fastq}" "${read_pattern}" \
                "${label} contains expected read name"

            awk 'NR % 4 == 1 { n++ } END { print n + 0 }' \
                "${view_fastq}" > "${count_reads}"
            assert_grep_pattern "${count_reads}" '^1$' \
                "${label} contains one read"
        fi
    fi
}


#  Assert a custom FASTQ path exists and is represented as a symlink
function assert_custom_symlink() {
    local symlink="${1:-}"
    local label="${2:-custom symlink}"

    assert_file_nonempty "${symlink}" "${label} target"

    if [[ -L "${symlink}" ]]; then
        rec_pass "${label} path is a symlink"
    else
        rec_fail "${label} path is not a symlink"
    fi
}


#  Warn when a help-output section is missing
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
        '\nSummary for %s: pass=%d fail=%d warn=%d skip=%d\n' \
        "${TEST_NAME}" "${TEST_PASS}" "${TEST_FAIL}" "${TEST_WARN}" \
        "${TEST_SKIP}"

    if (( TEST_FAIL > 0 )); then
        return 1
    fi

    return 0
}
