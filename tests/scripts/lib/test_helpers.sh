#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_helpers.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development.
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
    TEST_DIR_LIB="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"
    TEST_DIR_SCR="$(cd "${TEST_DIR_LIB}/.." > /dev/null 2>&1 && pwd)"
    TEST_DIR="$(cd "${TEST_DIR_SCR}/.." > /dev/null 2>&1 && pwd)"
    ROOT_REPO="$(cd "${TEST_DIR}/.." > /dev/null 2>&1 && pwd)"
    TEST_DIR_OUT="${ROOT_REPO}/tests/outputs"
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
function print_relpath() {
    local path="${1:-}"

    if [[ "${path}" == "${ROOT_REPO}"/* ]]; then
        printf '%s\n' "${path#"${ROOT_REPO}/"}"
    else
        printf '%s\n' "${path}"
    fi
}


#  Print a smoke-test section heading
function print_section() {
    printf '\n-- %s --\n' "${1:-${TEST_NAME}}"
}


#  Record a passing assertion
function record_pass() {
    TEST_PASS=$(( TEST_PASS + 1 ))
    printf 'PASS: %s\n' "$*"
}


#  Record a failing assertion
function record_fail() {
    TEST_FAIL=$(( TEST_FAIL + 1 ))
    printf 'FAIL: %s\n' "$*" >&2
}


#  Record a non-fatal warning
function record_warn() {
    TEST_WARN=$(( TEST_WARN + 1 ))
    printf 'WARN: %s\n' "$*" >&2
}


#  Record a skipped check
function record_skip() {
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


#  Check whether download smoke tests were explicitly requested
function is_download_enabled() {
    [[ "${RUN_DOWNLOAD:-0}" == "1" ]]
}


#TODO: is_slurm_enabled and is_slurm_wait
# function is_slurm_enabled() {
#     [[ "${RUN_SLURM:-0}" == "1" ]]
# }
#
#
# function is_slurm_wait() {
#     [[ "${SLURM_WAIT:-0}" == "1" || "${WAIT_SLURM:-0}" == "1" ]]
# }


#  Log command availability in the current shell and project environment
function check_env_cmds() {
    local env_nam="${1:-env_protocol}"
    shift || true

    local -a arr_cmd=( "$@" )
    local cmd=""

    echo "Current shell:"
    echo "SHELL=${SHELL:-UNSET}"
    echo "BASH=${BASH:-UNSET}"
    echo "PATH=${PATH:-UNSET}"
    echo "CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-UNSET}"
    echo "CONDA_PREFIX=${CONDA_PREFIX:-UNSET}"
    echo "MAMBA_EXE=${MAMBA_EXE:-UNSET}"
    echo

    echo "Current shell command checks:"
    for cmd in "${arr_cmd[@]}"; do
        printf '%s: ' "${cmd}"
        command -v "${cmd}" || true
    done
    echo

    if [[ "${CONDA_DEFAULT_ENV:-}" == "${env_nam}" ]]; then
        echo "Project-env command checks:"
        for cmd in "${arr_cmd[@]}"; do
            printf '%s: ' "${cmd}"
            command -v "${cmd}"
        done
    elif \
        check_cmd_exists conda
    then
        echo "Project-env command checks via Conda activation:"
        ENV_NAM="${env_nam}" \
        CMD_LCL="${arr_cmd[*]}" \
        bash -lc '
            eval "$(conda shell.bash hook)"
            conda activate "${ENV_NAM}"
            echo "CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-UNSET}"
            echo "CONDA_PREFIX=${CONDA_PREFIX:-UNSET}"
            echo "PATH=${PATH:-UNSET}"
            echo

            read -r -a arr_cmd <<< "${CMD_LCL}"
            for cmd in "${arr_cmd[@]}"; do
                printf "%s: " "${cmd}"
                command -v "${cmd}"
            done
        '
    else
        echo "conda is unavailable; cannot inspect '${env_nam}'."
        return 1
    fi
}


#  Require Atria and compression helpers in the requested project environment
function require_env_atria() {
    local env_nam="${1:-env_protocol}"
    local log_lcl="${2:-${TEST_DIR_LOG}/atria_project_env.log}"
    local rc=0

    mkdir -p "$(dirname "${log_lcl}")"

    check_env_cmds \
        "${env_nam}" \
        atria pigz pbzip2 gzip \
        > "${log_lcl}" 2>&1 || rc=1

    if (( rc == 0 )); then
        record_pass \
            "Atria trim dependencies are available in project environment"
        return 0
    fi

    record_fail \
        "Atria trim dependencies unavailable in project environment; see" \
        "$(print_relpath "${log_lcl}")"
    return 1
}


#  Require GNU Parallel in the requested project environment
function require_env_parallel() {
    local env_nam="${1:-env_protocol}"
    local log_lcl="${2:-${TEST_DIR_LOG}/parallel_project_env.log}"
    local rc=0

    mkdir -p "$(dirname "${log_lcl}")"

    check_env_cmds \
        "${env_nam}" \
        parallel \
        > "${log_lcl}" 2>&1 || rc=1

    if (( rc == 0 )); then
        record_pass "GNU Parallel is available in project environment"
        return 0
    fi

    record_fail \
        "GNU Parallel unavailable in project environment; see" \
        "$(print_relpath "${log_lcl}")"
    return 1
}


#  Require wget and gzip in the requested project environment
function require_env_download() {
    local env_nam="${1:-env_protocol}"
    local log_lcl="${2:-${TEST_DIR_LOG}/download_fastqs_project_env.log}"
    local rc=0

    mkdir -p "$(dirname "${log_lcl}")"

    check_env_cmds \
        "${env_nam}" \
        wget gzip \
        > "${log_lcl}" 2>&1 || rc=1

    if (( rc == 0 )); then
        record_pass \
            "download dependencies are available in project environment"
        return 0
    fi

    record_fail \
        "download dependencies unavailable in project environment; see" \
        "$(print_relpath "${log_lcl}")"
    return 1
}


#  Resolve the active project environment or require a named fallback
function require_env_project() {
    local env_ref="${1:-env_nam}"
    local env_fb="${2:-env_protocol}"

    if [[
        -n "${CONDA_DEFAULT_ENV:-}" && "${CONDA_DEFAULT_ENV}" != "base"
    ]]; then
        printf -v "${env_ref}" '%s' "${CONDA_DEFAULT_ENV}"
        return 0
    fi

    if ! \
        check_cmd_exists conda
    then
        record_fail \
            "setup requires active project env or conda to activate" \
            "'${env_fb}'"
        return 1
    fi

    if ! \
        conda env list 2> /dev/null \
            | awk -v env="${env_fb}" '
                $1 == env { found = 1 } END { exit !found }
            '
    then
        record_fail "setup requires Conda/Mamba environment '${env_fb}'"
        return 1
    fi

    printf -v "${env_ref}" '%s' "${env_fb}"
}


#  Require one or more fixture files to exist and be non-empty
function require_files_nonempty() {
    local file=""
    local rc=0

    for file in "$@"; do
        if [[ ! -s "${file}" ]]; then
            record_fail "missing fixture $(print_relpath "${file}")"
            rc=1
        fi
    done

    return "${rc}"
}


#  Assert one predictable output file exists, without requiring non-emptiness
function assert_file_exists() {
    local file="${1:-}"
    shift || true

    local lbl="${*:-file exists}"

    if [[ -f "${file}" ]]; then
        record_pass "${lbl}"
    else
        record_fail "${lbl}; missing $(print_relpath "${file}")"
    fi
}


#  Assert exactly one path was found by a caller-specific search
function assert_path_found() {
    local arr_ref="${1:-}"
    local lbl="${2:-path}"
    local dir_srch="${3:-}"
    local out_ref="${4:-}"
    local pth_fnd=""

    local -n arr_lcl="${arr_ref}"

    if (( ${#arr_lcl[@]} == 1 )); then
        pth_fnd="${arr_lcl[0]}"
        printf -v "${out_ref}" '%s' "${pth_fnd}"
        record_pass "one ${lbl} found"
    else
        printf -v "${out_ref}" ''
        record_fail \
            "expected exactly one ${lbl} in $(print_relpath "${dir_srch}")," \
            "found ${#arr_lcl[@]}"
    fi
}


#  Run samtools from the active shell or the resolved project environment
function run_samtools() {
    local env_lcl="${env_nam:-env_protocol}"

    if \
        check_cmd_exists samtools
    then
        samtools "$@"
    else
        conda run -n "${env_lcl}" samtools "$@"
    fi
}


#  Build and index a BAM input fixture from a committed SAM fixture
function build_filter_alignments_bam_fixture() {
    local in_sam="${1:-}"
    local out_bam="${2:-}"
    local log_lcl="${3:-}"
    shift 3 || true

    local lbl="${*:-BAM input fixture}"

    if \
        run_capture \
            "prepare ${lbl}" "${log_lcl}" \
            run_samtools view -bS -o "${out_bam}" "${in_sam}" \
        && run_capture \
            "index ${lbl}" "${log_lcl}.index" \
            run_samtools index "${out_bam}"
    then
        record_pass "${lbl} is prepared"
        return 0
    fi

    record_fail "failed to prepare ${lbl}; see $(print_relpath "${log_lcl}")"
    return 1
}


#  Build and index a CRAM input fixture from committed SAM/reference fixtures
function build_filter_alignments_cram_fixture() {
    local in_sam="${1:-}"
    local ref_fa="${2:-}"
    local out_cram="${3:-}"
    local log_lcl="${4:-}"
    shift 4 || true

    local lbl="${*:-filter CRAM input fixture}"

    if \
        run_capture \
            "prepare ${lbl}" "${log_lcl}" \
            run_samtools view -C -T "${ref_fa}" -o "${out_cram}" "${in_sam}" \
        && run_capture \
            "index ${lbl}" "${log_lcl}.index" \
            run_samtools index "${out_cram}"
    then
        record_pass "${lbl} is prepared"
        return 0
    fi

    record_fail "failed to prepare ${lbl}; see $(print_relpath "${log_lcl}")"
    return 1
}


#  Find a usable Python command
function find_python() {
    if \
        check_cmd_exists python
    then
        command -v python
    elif \
        check_cmd_exists python3
    then
        command -v python3
    else
        return 1
    fi
}


#  Find a Python command that can bind loopback sockets
function find_python_loopback() {
    local py=""
    local -a arr_py=()

    if check_cmd_exists python; then
        arr_py+=( "$(command -v python)" )
    fi

    if check_cmd_exists python3; then
        arr_py+=( "$(command -v python3)" )
    fi

    if [[ -x /usr/bin/python3 ]]; then
        arr_py+=( /usr/bin/python3 )
    fi

    for py in "${arr_py[@]}"; do
        if \
            "${py}" - << PY > /dev/null 2>&1
import socket
with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
    sock.bind(("127.0.0.1", 0))
PY
        then
            printf '%s\n' "${py}"
            return 0
        fi
    done

    return 1
}


#  Find an available loopback TCP port for local HTTP smoke tests
function find_port_free() {
    local py="${1:-}"

    "${py}" - << PY
import socket
with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
    sock.bind(("127.0.0.1", 0))
    print(sock.getsockname()[1])
PY
}


#  Wait until a local HTTP server is reachable
function wait_http_local() {
    local py="${1:-}"
    local url="${2:-}"
    local try="${3:-50}"
    local i=0

    for (( i = 1; i <= try; i++ )); do
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
function cleanup_server_http() {
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
    local nam="${1:-command}"
    local out="${2:-}"

    shift 2

    if [[ -z "${out}" ]]; then
        out="${TEST_DIR_LOG}/${nam//[^A-Za-z0-9_.-]/_}.log"
    fi

    mkdir -p "$(dirname "${out}")"
    "$@" > "${out}" 2>&1
}


#  Run one filter-alignments retain-mode wrapper case
function run_case_filter() {
    local retain="${1:-}"
    local nam_cas="${2:-}"
    local log_lcl="${3:-}"
    local arr_cmd_ref="${4:-}"
    local nam_job_pfx="${5:-}"
    local arr_arg_ref="${6:-}"
    shift 6 || true

    local lbl="${*:-filter-alignments retain=${retain} ${nam_cas}}"

    if [[ -z "${arr_cmd_ref}" ]]; then
        record_fail "${lbl}; missing command array reference"
        return 1
    fi

    if [[ -z "${nam_job_pfx}" ]]; then
        record_fail "${lbl}; missing job-name prefix"
        return 1
    fi

    if [[ -z "${arr_arg_ref}" ]]; then
        record_fail "${lbl}; missing argument array reference"
        return 1
    fi

    # shellcheck disable=SC2178
    local -n arr_cmd="${arr_cmd_ref}"
    local -n arr_arg="${arr_arg_ref}"

    if (( ${#arr_cmd[@]} == 0 )); then
        record_fail "${lbl}; command array is empty"
        return 1
    fi

    if \
        run_capture \
            "${lbl}" "${log_lcl}" \
            "${arr_cmd[@]}" \
                --retain "${retain}" \
                --nam_job "${nam_job_pfx}_${nam_cas}" \
                "${arr_arg[@]}"
    then
        record_pass "${lbl} exits 0"
    else
        record_fail "${lbl} failed; see $(print_relpath "${log_lcl}")"
    fi
}


#  Run one compute-signal wrapper case
function run_case_compute_signal() {
    local wrap="${1:-}"
    local fmt_in="${2:-}"
    local nam_cas="${3:-}"
    local mode="${4:-}"
    local in_lcl="${5:-}"
    local out_spc="${6:-}"
    local log_lcl="${7:-}"
    local dir_out_lcl="${8:-}"
    local dir_err_lcl="${9:-}"
    local ref_fa_lcl="${10:-}"

    shift 10 || true

    local fmt_up="${fmt_in^^}"
    local script="${ROOT_REPO}/scripts/${wrap}_compute_signal.sh"
    local nam_job=""
    local -a arr_cmd=()

    case "${wrap}" in
        submit)
            nam_job="test_compute_${fmt_in}_${nam_cas}"
            # shellcheck disable=SC2154
            arr_cmd=(
                "${TEST_BASH}" "${script}"
                    --env_nam "${env_nam}"
                    --dir_scr "${ROOT_REPO}/scripts"
                    --threads 1
                    --mode "${mode}"
                    --csv_infile "${in_lcl}"
                    --csv_outfile "${out_spc}"
            )
            ;;
        execute)
            nam_job="test_execute_compute_${fmt_in}_${nam_cas}"
            arr_cmd=(
                "${TEST_BASH}" "${script}"
                    --threads 1
                    --mode "${mode}"
                    --csv_infile "${in_lcl}"
                    --dir_out "${dir_out_lcl}"
                    --typ_out "${out_spc}"
            )
            ;;
        *)
            record_fail "unknown compute-signal wrapper '${wrap}'"
            return 1
            ;;
    esac

    if [[ -n "${ref_fa_lcl}" ]]; then
        arr_cmd+=( --ref_fa "${ref_fa_lcl}" )
    fi

    arr_cmd+=(
        --csv_usr_frg NA
    )

    if [[ "${wrap}" == "execute" ]]; then
        arr_cmd+=( --dp 3 )
    fi

    arr_cmd+=(
        --err_out "${dir_err_lcl}"
        --nam_job "${nam_job}"
    )

    if [[ "${wrap}" == "execute" ]]; then
        arr_cmd+=( --max_job 1 )
    fi

    if \
        run_capture \
            "${wrap} compute-signal ${fmt_up} ${nam_cas}" \
            "${log_lcl}" \
            "${arr_cmd[@]}" \
            "$@"
    then
        record_pass "${wrap}_compute_signal.sh ${mode} ${nam_cas} exits 0"
    else
        record_fail \
            "${wrap}_compute_signal.sh ${mode} ${nam_cas} failed; see" \
            "$(print_relpath "${log_lcl}")"
    fi
}


#  Run one compute-signal ratio wrapper case
function run_case_compute_signal_ratio() {
    local wrap="${1:-}"
    local cas_nam="${2:-}"
    local out_spec="${3:-}"
    local log_lcl="${4:-}"
    local method="${5:-unadj}"
    local fil_A="${6:-}"
    local fil_B="${7:-}"
    local dir_out_lcl="${8:-}"
    local dir_err_lcl="${9:-}"

    shift 9 || true

    local script="${ROOT_REPO}/scripts/${wrap}_compute_signal.sh"
    local -a arr_cmd=()

    case "${wrap}" in
        submit)
            # shellcheck disable=SC2154
            arr_cmd=(
                "${TEST_BASH}" "${script}"
                --env_nam "${env_nam}"
                --dir_scr "${ROOT_REPO}/scripts"
                --threads 1
                --mode ratio
                --method "${method}"
                --csv_fil_A "${fil_A}"
                --csv_fil_B "${fil_B}"
                --csv_outfile "${out_spec}"
                --err_out "${dir_err_lcl}"
                --nam_job "test_submit_compute_ratio_${cas_nam}"
            )
            ;;
        execute)
            arr_cmd=(
                "${TEST_BASH}" "${script}"
                --threads 1
                --mode ratio
                --method "${method}"
                --csv_fil_A "${fil_A}"
                --csv_fil_B "${fil_B}"
                --dir_out "${dir_out_lcl}"
                --typ_out bdg
                --prefix "${out_spec}"
                --eps 0
                --dp 3
                --err_out "${dir_err_lcl}"
                --nam_job "test_execute_compute_ratio_${cas_nam}"
                --max_job 1
            )
            ;;
        *)
            record_fail "unknown compute-signal ratio wrapper '${wrap}'"
            return 1
            ;;
    esac

    if \
        run_capture \
            "${wrap} compute-signal ratio ${cas_nam}" \
            "${log_lcl}" \
            "${arr_cmd[@]}" \
            "$@"
    then
        record_pass "${wrap}_compute_signal.sh ratio ${cas_nam} exits 0"
    else
        record_fail \
            "${wrap}_compute_signal.sh ratio ${cas_nam} failed; see" \
            "$(print_relpath "${log_lcl}")"
    fi
}


#  Assert that a log file contains a pattern
function assert_pattern_found() {
    local file="${1:-}"
    local patn="${2:-}"
    shift 2 || true

    local lbl="${*:-${patn}}"

    if \
        grep -q -- "${patn}" "${file}"
    then
        record_pass "${lbl}"
    else
        record_fail "${lbl}; see $(print_relpath "${file}")"
    fi
}


#  Assert that a log or output file does not contain a pattern
function assert_pattern_absent() {
    local file="${1:-}"
    local patn="${2:-}"
    shift 2 || true

    local lbl="${*:-${patn}}"

    if \
        grep -q -- "${patn}" "${file}"
    then
        record_fail "${lbl}; see $(print_relpath "${file}")"
    else
        record_pass "${lbl}"
    fi
}


#  Assert that a generated output file exists and is non-empty
function assert_file_nonempty() {
    local file="${1:-}"
    shift || true

    local lbl="${*:-output}"

    if [[ -s "${file}" ]]; then
        record_pass "${lbl} exists and is non-empty"
    else
        record_fail "${lbl} missing or empty: $(print_relpath "${file}")"
    fi
}


#  Assert that a file contains exactly one expected line
function assert_file_exact_line() {
    local file="${1:-}"
    local expected="${2:-}"
    shift 2 || true

    local lbl="${*:-output line}"
    local observed
    local n_line

    if [[ ! -s "${file}" ]]; then
        record_fail "${lbl} missing or empty: $(print_relpath "${file}")"
        return 1
    fi

    n_line="$(wc -l < "${file}")"
    if [[ "${n_line}" -eq 1 ]]; then
        record_pass "${lbl} has one row"
    else
        record_fail "${lbl} has unexpected row count; see $(print_relpath "${file}")"
    fi

    observed="$(cat "${file}")"
    if [[ "${observed}" == "${expected}" ]]; then
        record_pass "${lbl} has expected row"
    else
        record_fail "${lbl} differs from expected; see $(print_relpath "${file}")"
    fi
}


#  Assert whether a scaling-factor output contains an expected header
function assert_scaling_factor_header() {
    local file="${1:-}"
    local header="${2:-}"
    local expect="${3:-true}"
    shift 3 || true

    local lbl="${*:-scaling-factor header}"

    if [[ "${expect}" == "true" ]]; then
        assert_pattern_found "${file}" "${header}" "${lbl} has core header"
    else
        assert_pattern_absent "${file}" "${header}" "${lbl} omits core header"
    fi
}


#  Run one execute-level calculate-scaling-factor two-row case
function run_case_scaling_factor_execute() {
    local cas="${1:-}"
    local mode="${2:-}"
    local arr_cmd_nam="${3:-}"
    local dir_out_lcl="${4:-}"
    local dir_log_lcl="${5:-}"
    local out_suffix="${6:-}"
    local job_prefix="${7:-}"
    local header="${8:-}"
    local row_0="${9:-}"
    local row_1="${10:-}"
    local tail_0="${11:-}"
    local tail_1="${12:-}"
    shift 12 || true

    local -n arr_cmd_ref="${arr_cmd_nam}"

    local expect_header=true
    local fil_out="${dir_out_lcl}/scaling.${cas}.${out_suffix}.tsv"
    local prt_0="${fil_out}.part.000000"
    local prt_1="${fil_out}.part.000001"
    local nam_job="${job_prefix}_${cas}"
    local log="${dir_log_lcl}/execute_${out_suffix}_${cas}.log"
    local -a arr_case=(
        "${arr_cmd_ref[@]}"
        --fil_out "${fil_out}"
        --nam_job "${nam_job}"
    )

    for arg in "$@"; do
        if [[ "${arg}" =~ ^--no[_-]header$ ]]; then
            expect_header=false
        fi
    done

    arr_case+=( "$@" )

    if \
        run_capture \
            "execute calculate-scaling-factor ${mode} ${cas}" \
            "${log}" \
            "${arr_case[@]}"
    then
        record_pass "execute_calculate_scaling_factor.sh ${mode} ${cas} exits 0"
    else
        record_fail \
            "execute_calculate_scaling_factor.sh ${mode} ${cas} failed; see" \
            "$(print_relpath "${log}")"
    fi

    assert_file_nonempty \
        "${fil_out}" \
        "execute scaling-factor ${mode} ${cas} final TSV"

    assert_file_nonempty \
        "${prt_0}" \
        "execute scaling-factor ${mode} ${cas} first retained part"

    assert_file_nonempty \
        "${prt_1}" \
        "execute scaling-factor ${mode} ${cas} second retained part"

    assert_scaling_factor_header \
        "${fil_out}" \
        "${header}" \
        "${expect_header}" \
        "execute scaling-factor ${mode} ${cas} final TSV"

    assert_pattern_found \
        "${fil_out}" \
        "^${row_0}"$'\t'"${tail_0}"'$' \
        "execute scaling-factor ${mode} ${cas} final TSV has first row"

    assert_pattern_found \
        "${fil_out}" \
        "^${row_1}"$'\t'"${tail_1}"'$' \
        "execute scaling-factor ${mode} ${cas} final TSV has second row"
}


#  Run one submit-level calculate-scaling-factor one-part case
function run_case_scaling_factor_submit_part() {
    local cas="${1:-}"
    local mode="${2:-}"
    local arr_cmd_nam="${3:-}"
    local fil_out="${4:-}"
    local idx_out="${5:-}"
    local fil_log="${6:-}"
    local row_exp="${7:-}"
    local header="${8:-}"
    shift 8 || true

    local -n arr_cmd_ref="${arr_cmd_nam}"

    local fil_prt="${fil_out}.part.$(printf '%06d' "${idx_out}")"

    if \
        run_capture \
            "submit calculate-scaling-factor ${mode} ${cas}" \
            "${fil_log}" \
            "${arr_cmd_ref[@]}" \
            "$@"
    then
        record_pass "submit_calculate_scaling_factor.sh ${mode} ${cas} exits 0"
    else
        record_fail \
            "submit_calculate_scaling_factor.sh ${mode} ${cas} failed; see" \
            "$(print_relpath "${fil_log}")"
    fi

    assert_file_nonempty \
        "${fil_prt}" \
        "submit scaling-factor ${mode} ${cas} indexed part"

    if [[ -n "${header}" ]]; then
        assert_pattern_absent \
            "${fil_prt}" \
            "${header}" \
            "submit scaling-factor ${mode} ${cas} part stays data-only"
    fi

    if [[ ! -e "${fil_out}" ]]; then
        record_pass "submit_calculate_scaling_factor.sh ${mode} ${cas} writes no final TSV"
    else
        record_fail "submit_calculate_scaling_factor.sh ${mode} ${cas} wrote final TSV"
    fi

    assert_file_exact_line \
        "${fil_prt}" \
        "${row_exp}" \
        "submit scaling-factor ${mode} ${cas} part"
}


#  Assert that a CRAM index exists using common Samtools naming variants
function assert_cram_index() {
    local cram="${1:-}"
    shift || true

    local lbl="${*:-CRAM index}"
    local idx_1="${cram}.crai"
    local idx_2="${cram%.cram}.crai"

    if [[ -s "${idx_1}" ]]; then
        record_pass \
            "${lbl} exists and is non-empty at $(print_relpath "${idx_1}")"
    elif [[ -s "${idx_2}" ]]; then
        record_pass \
            "${lbl} exists and is non-empty at $(print_relpath "${idx_2}")"
    else
        record_fail \
            "${lbl} missing or empty; checked $(print_relpath "${idx_1}")" \
            "and $(print_relpath "${idx_2}")"
    fi
}


#  Assert a reference-backed CRAM read count for one or more regions
function assert_cram_count() {
    local cram="${1:-}"
    local ref_fa="${2:-}"
    local exp="${3:-}"
    local fil_cnt="${4:-}"
    local arr_ref="${5:-}"
    shift 5 || true

    local lbl="${*:-CRAM count}"

    if [[ -z "${arr_ref}" ]]; then
        record_fail "${lbl}; missing region array reference"
        return 1
    fi

    local -n arr_rgn="${arr_ref}"

    if (( ${#arr_rgn[@]} == 0 )); then
        record_fail "${lbl}; region array is empty"
        return 1
    fi

    run_capture \
        "${lbl}" "${fil_cnt}" \
        run_samtools view -T "${ref_fa}" -c "${cram}" "${arr_rgn[@]}"

    assert_pattern_found \
        "${fil_cnt}" \
        "^${exp}$" \
        "${lbl}"
}


#  Assert filter_alignments @PG provenance in a BAM or CRAM header
function assert_filter_alignments_pg_header() {
    local infile="${1:-}"
    local ref_fa="${2:-}"
    local pg_id="${3:-}"
    local retain="${4:-}"
    local out_ext="${5:-}"
    local header_file="${6:-}"
    shift 6 || true

    local lbl="${*:-filter_alignments @PG header}"
    local -a ref_arg=()

    if [[ "${infile,,}" == *.cram ]]; then
        ref_arg=( -T "${ref_fa}" )
    fi

    if \
        run_capture \
            "header ${lbl}" "${header_file}" \
            run_samtools view -H "${ref_arg[@]}" "${infile}"
    then
        record_pass "${lbl} header can be read"
    else
        record_fail \
            "${lbl} header cannot be read; see" \
            "$(print_relpath "${header_file}")"
        return 1
    fi

    if [[ -s "${header_file}" ]]; then
        assert_pattern_found \
            "${header_file}" \
            $'^@PG\tID:'"${pg_id}"$'\tPN:filter_alignments\tCL:'"${pg_id} retain=${retain}" \
            "${lbl} contains ${pg_id} @PG record"
        assert_pattern_found \
            "${header_file}" \
            "out_ext=${out_ext}" \
            "${lbl} records out_ext=${out_ext}"
    fi
}


#  Assert a gzip FASTQ exists and has expected read content
function assert_fastq_gzip() {
    local fil_out="${1:-}"
    local rd_patn="${2:-}"
    local fq_vw="${3:-}"
    local rd_cnt="${4:-}"
    local fq_src="${5:-}"
    shift 5 || true

    local lbl="${*:-gzip FASTQ}"

    if [[ -z "${fil_out}" ]]; then
        return 0
    fi

    assert_file_nonempty \
        "${fil_out}" \
        "${lbl}"

    if [[ -s "${fil_out}" && -n "${fq_src}" ]]; then
        if \
            cmp -s "${fq_src}" "${fil_out}"
        then
            record_pass "${lbl} matches source fixture byte-for-byte"
        else
            record_fail "${lbl} differs from source fixture"
        fi
    fi

    if \
        gzip -t "${fil_out}"
    then
        record_pass "${lbl} passes gzip integrity"
    else
        record_fail "${lbl} fails gzip integrity"
    fi

    if \
        gzip -cd "${fil_out}" > "${fq_vw}"
    then
        record_pass "${lbl} can be decompressed"
    else
        record_fail "${lbl} cannot be decompressed"
    fi

    if [[ -s "${fq_vw}" ]]; then
        assert_pattern_found \
            "${fq_vw}" \
            "${rd_patn}" \
            "${lbl} contains expected read name"

        awk 'NR % 4 == 1 { n++ } END { print n + 0 }' "${fq_vw}" > "${rd_cnt}"

        assert_pattern_found \
            "${rd_cnt}" \
            '^1$' \
            "${lbl} contains one read"
    fi
}


#  Assert a custom FASTQ path exists and is represented as a symlink
function assert_custom_symlink() {
    local sym="${1:-}"
    shift || true

    local lbl="${*:-custom symlink}"

    assert_file_nonempty \
        "${sym}" \
        "${lbl}" \
        "target"

    if [[ -L "${sym}" ]]; then
        record_pass "${lbl} path is a symlink"
    else
        record_fail "${lbl} path is not a symlink"
    fi
}


#  Warn when a help-output section is missing
function warn_help_pattern_missing() {
    local file="${1:-}"
    local patn="${2:-}"
    shift 2 || true

    local lbl="${*:-${patn}}"

    if \
        grep -q -- "${patn}" "${file}"
    then
        record_pass "${lbl}"
    else
        record_warn "${lbl}; see $(print_relpath "${file}")"
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
