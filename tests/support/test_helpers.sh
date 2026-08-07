#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_helpers.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in design, development,
# and documentation, with all output reviewed, edited, and approved by the
# author.
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


#  Set paths used by test scripts =============================================
# shellcheck disable=SC2034
{
    TEST_DIR_LIB="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"
    TEST_DIR="$(cd "${TEST_DIR_LIB}/.." > /dev/null 2>&1 && pwd)"
    TEST_DIR_SCR="${TEST_DIR}"
    ROOT_REPO="$(cd "${TEST_DIR}/.." > /dev/null 2>&1 && pwd)"
    TEST_DIR_OUT="${TEST_ARTIFACT_ROOT:-${ROOT_REPO}/artifacts/tests}"
    if [[ "${TEST_DIR_OUT}" != /* || "${TEST_DIR_OUT}" == "/" ]]; then
        echo "error(test_helpers.sh):" \
            "TEST_ARTIFACT_ROOT must be an absolute non-root path." >&2
        exit 1
    fi
    TEST_OUT_PARENT="$(
        cd "$(dirname "${TEST_DIR_OUT}")" > /dev/null 2>&1 && pwd -P
    )" || {
        echo "error(test_helpers.sh):" \
            "TEST_ARTIFACT_ROOT parent does not exist." >&2
        exit 1
    }
    TEST_DIR_OUT="${TEST_OUT_PARENT}/$(basename "${TEST_DIR_OUT}")"
    case "${TEST_DIR_OUT}/" in
        "${ROOT_REPO}/"*)
            if [[ "${TEST_DIR_OUT}" != "${ROOT_REPO}/artifacts/tests" ]]; then
                echo "error(test_helpers.sh):" \
                    "TEST_ARTIFACT_ROOT must be outside the repository." >&2
                exit 1
            fi
            ;;
    esac
    TEST_DIR_TMP="${TEST_DIR_OUT}/tmp"
    TEST_DIR_LOG="${TEST_DIR_OUT}/logs"
    TEST_ENV_PREFIX="${TEST_MANAGED_PREFIX:-}"
    if [[ -z "${TEST_ENV_PREFIX}" ]]; then
        if [[ "${CONDA_DEFAULT_ENV:-}" == "env_protocol" && \
            -n "${CONDA_PREFIX:-}" ]]
        then
            TEST_ENV_PREFIX="${CONDA_PREFIX}"
        else
            CONDA_EXECUTABLE="${CONDA_EXE:-}"
            if [[ -z "${CONDA_EXECUTABLE}" || \
                ! -x "${CONDA_EXECUTABLE}" ]]
            then
                CONDA_EXECUTABLE="$(command -v conda || true)"
            fi
            if [[ -z "${CONDA_EXECUTABLE}" ]]; then
                echo "error(test_helpers.sh):" \
                    "Conda is required to resolve env_protocol." >&2
                exit 1
            fi
            TEST_ENV_PREFIX="$(${CONDA_EXECUTABLE} run -n env_protocol \
                /bin/sh -c 'printf "%s\n" "$CONDA_PREFIX"')"
        fi
    fi
    TEST_BASH="${TEST_ENV_PREFIX}/bin/bash"
    TEST_PYTHON="${TEST_ENV_PREFIX}/bin/python"
    if [[ ! -x "${TEST_BASH}" || ! -x "${TEST_PYTHON}" ]]; then
        echo "error(test_helpers.sh):" \
            "env_protocol managed executables are unavailable." >&2
        exit 1
    fi
    TEST_BASH_VERSION="$(
        "${TEST_BASH}" -c 'printf "%s\n" "${BASH_VERSION}"'
    )"
}

export TEST_ARTIFACT_ROOT="${TEST_DIR_OUT}"
export TEST_MANAGED_PREFIX="${TEST_ENV_PREFIX}"
export TEST_MANAGED_BASH="${TEST_BASH}"
export TEST_MANAGED_PYTHON="${TEST_PYTHON}"
export CONDA_DEFAULT_ENV=env_protocol
export CONDA_PREFIX="${TEST_ENV_PREFIX}"
export PATH="${TEST_ENV_PREFIX}/bin:${PATH}"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache"
export PYTEST_ADDOPTS="${PYTEST_ADDOPTS:+${PYTEST_ADDOPTS} }-o cache_dir=${TEST_DIR_OUT}/pytest_cache"
mkdir -p "${TEST_DIR_TMP}" "${TEST_DIR_LOG}"


#  Source canonical production validation helpers used by the test harness
# shellcheck source=lib/bash/core/format_outputs.sh
# shellcheck source=lib/bash/core/check_args.sh
{
    source "${ROOT_REPO}/lib/bash/core/format_outputs.sh"
    source "${ROOT_REPO}/lib/bash/core/check_args.sh"
}


#  Initialize result counters =================================================
TEST_PASS=0
TEST_FAIL=0
TEST_WARN=0
TEST_SKIP=0
TEST_NAME="${TEST_NAME:-$(basename "${0}")}"


#  Print a repository-relative path when possible
function print_relpath() {
    local path="${1:-}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  print_relpath
    [--help] path

  Print a repository-relative path when the path is under the repository root.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  path : path
    Absolute or relative path to display.

Returns
-------
  Prints a display path to stdout and returns 0.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Display a repository-owned test path relative to the repository root.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    print_relpath "${ROOT_REPO}/tests/support/test_helpers.sh"
    '''

  2. Preserve an external temporary path that is outside the repository.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    print_relpath /private/tmp/test.log
    '''
EOM
    )

    if [[ "${path}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ "${path}" == "${ROOT_REPO}"/* ]]; then
        printf '%s\n' "${path#"${ROOT_REPO}/"}"
    else
        printf '%s\n' "${path}"
    fi
}


#  Print a test section heading
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


#  Normalize one optional Boolean test gate
function invalid_test_gate_name() {
    echo_err_func "normalize_test_gate" \
        "positional argument 1, 'name', must be a valid shell variable" \
        "name: '${1:-}'."
}


function normalize_test_gate() {
    local name="${1:-}" value patn='^[A-Za-z_][A-Za-z0-9_]*$'
    [[ "${name}" =~ ${patn} ]] || { invalid_test_gate_name "${name}"; return 1; }
    value="${!name-}"
    [[ -n "${value}" ]] || { printf '%s\n' false; return 0; }
    normalize_bool "${value}" "${name}"
}


#  Check whether GNU Parallel tests were explicitly requested
function is_parallel_enabled() {
    local value

    value="$(normalize_test_gate RUN_PARALLEL)" || exit 1
    [[ "${value}" == "true" ]]
}


#  Check whether Atria tests were explicitly requested
function is_atria_enabled() {
    local value

    value="$(normalize_test_gate RUN_ATRIA)" || exit 1
    [[ "${value}" == "true" ]]
}


#  Check whether download tests were explicitly requested
function is_download_enabled() {
    local value

    value="$(normalize_test_gate RUN_DOWNLOAD)" || exit 1
    [[ "${value}" == "true" ]]
}


#  Check whether Slurm tests were explicitly requested
function is_slurm_enabled() {
    local value

    value="$(normalize_test_gate RUN_SLURM)" || exit 1
    [[ "${value}" == "true" ]]
}


#  Check whether Slurm output polling was explicitly requested
function is_slurm_wait_enabled() {
    local value

    value="$(normalize_test_gate WAIT_SLURM)" || exit 1
    [[ "${value}" == "true" ]]
}


#  Log command availability in the current shell and project environment
function check_env_cmds() {
    local env_nam="${1:-env_protocol}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  check_env_cmds
    [--help] [env_nam] cmd [cmd...]

  Print command-resolution diagnostics for the current shell and the requested Conda environment.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  env_nam : str
    Conda environment to activate. Environment name to inspect (default: env_protocol).

  2+ cmd : str
    Commands to resolve in the current shell and project environment.

Returns
-------
  Prints diagnostic information to stdout. Returns 0 when commands can be inspected; returns 1 if Conda is unavailable for environment inspection.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment
    - bash >= 4.4
    - conda (when the requested environment is not active)

Examples
--------
  1. Inspect one command in the active project environment.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    check_env_cmds "${CONDA_DEFAULT_ENV:-env_protocol}" bash
    '''

  2. Inspect the download tools in the named fallback environment.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    check_env_cmds env_protocol wget gzip
    '''
EOM
    )

    if [[ "${env_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  require_env_atria
    [--help] [env_nam] [log_lcl]

  Require Atria and compression helpers in the requested project environment.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  env_nam : str
    Conda environment to activate. Environment name to inspect (default: env_protocol).

  2  log_lcl : file
    Log file that receives command-resolution diagnostics.

Returns
-------
  Records a pass/fail assertion and returns 0 when dependencies are available; otherwise records a failure and returns 1.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - atria
    - bash >= 4.4
    - conda (when the requested environment is not active)
    - dirname
    - gzip
    - mkdir
    - pbzip2
    - pigz

Examples
--------
  1. Require the Atria toolchain using the default environment and log path.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    require_env_atria
    '''

  2. Capture the expected failure for an environment without the Atria tools.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    if ! \
        require_env_atria missing_env "${tmp}/atria.log"
    then
        printf '%s\n' 'Atria environment unavailable as expected'
    fi
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${env_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  require_env_parallel
    [--help] [env_nam] [log_lcl]

  Require GNU Parallel in the requested project environment.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  env_nam : str
    Conda environment to activate. Environment name to inspect (default: env_protocol).

  2  log_lcl : file
    Log file that receives command-resolution diagnostics.

Returns
-------
  Records a pass/fail assertion and returns 0 when GNU Parallel is available; otherwise records a failure and returns 1.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing parallel
    - bash >= 4.4
    - conda (when the requested environment is not active)
    - dirname
    - mkdir
    - parallel

Examples
--------
  1. Require GNU Parallel using the default environment and log path.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    require_env_parallel
    '''

  2. Capture the expected failure for an environment without GNU Parallel.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    if ! \
        require_env_parallel missing_env "${tmp}/parallel.log"
    then
        printf '%s\n' 'GNU Parallel environment unavailable as expected'
    fi
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${env_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  require_env_download
    [--help] [env_nam] [log_lcl]

  Require download workflow tools in the requested project environment.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  env_nam : str
    Conda environment to activate. Environment name to inspect (default: env_protocol).

  2  log_lcl : file
    Log file that receives command-resolution diagnostics.

Returns
-------
  Records a pass/fail assertion and returns 0 when dependencies are available; otherwise records a failure and returns 1.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - bash >= 4.4
    - conda (when the requested environment is not active)
    - dirname
    - gzip
    - mkdir
    - wget

Examples
--------
  1. Require wget and gzip using the default environment and log path.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    require_env_download
    '''

  2. Capture the expected failure for an environment without download tools.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    if ! \
        require_env_download missing_env "${tmp}/download.log"
    then
        printf '%s\n' 'download environment unavailable as expected'
    fi
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${env_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  require_env_project
    [--help] [env_ref] [env_fb]

  Resolve the active project environment or require a named fallback environment to exist.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  env_ref : str
    Caller variable name that receives the resolved environment name.

  2  env_fb : str
    Fallback environment name to require when no project environment is active.

Returns
-------
  Writes the resolved environment name to the caller variable. Returns 0 on success; records a failure and returns 1 when no usable environment exists.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment
    - awk (when no non-base project environment is active)
    - bash >= 4.4
    - conda (when no non-base project environment is active)

Examples
--------
  1. Reuse the active non-base environment.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    CONDA_DEFAULT_ENV=env_protocol
    require_env_project env_nam
    printf '%s\n' "${env_nam}"
    '''

  2. Resolve the named fallback when no project environment is active.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    unset CONDA_DEFAULT_ENV
    require_env_project env_nam env_protocol
    printf '%s\n' "${env_nam}"
    '''
EOM
    )

    if [[ "${env_ref}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
        record_fail "setup requires Conda environment '${env_fb}'"
        return 1
    fi

    printf -v "${env_ref}" '%s' "${env_fb}"
}


#  Require one or more fixture files to exist and be non-empty
function require_files_nonempty() {
    local file=""
    local rc=0
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  require_files_nonempty
    [--help] file [file...]

  Require one or more fixture files to exist and be non-empty.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1+ file : file
    Fixture paths to check.

Returns
-------
  Records one failure per missing or empty file. Returns 0 when all files are non-empty; otherwise returns 1.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Accept one non-empty fixture file.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    printf '%s\n' fixture > "${tmp}/present.txt"
    require_files_nonempty "${tmp}/present.txt"
    rm -r -- "${tmp}"
    '''

  2. Reject a set containing an empty and a missing fixture.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    : > "${tmp}/empty.txt"
    if ! \
        require_files_nonempty "${tmp}/empty.txt" "${tmp}/missing.txt"
    then
        printf '%s\n' 'fixture validation failed as expected'
    fi
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${1:-}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  assert_file_exists
    [--help] file [label...]

  Assert that one predictable output file exists.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  file : file
    File path to check.

  2+ label : str
    Optional assertion label.

Returns
-------
  Records a pass when the file exists; otherwise records a failure.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Record a pass for an existing empty status file.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    : > "${tmp}/status.txt"
    assert_file_exists "${tmp}/status.txt"
    rm -r -- "${tmp}"
    '''

  2. Record a failure with a caller-supplied label for a missing file.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    assert_file_exists "${tmp}/missing.txt" "expected status file"
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${file}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  assert_path_found
    [--help] arr_ref lbl dir_srch out_ref

  Assert that a caller-specific search found exactly one path.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  arr_ref : str
    Name of the indexed array containing candidate paths.

  2  lbl : str
    Human-readable path label.

  3  dir_srch : dir
    Directory searched by the caller.

  4  out_ref : str
    Caller variable name that receives the found path.

Returns
-------
  Records a pass and writes the path when exactly one candidate exists; otherwise records a failure and clears the output variable.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Accept one discovered BAM and write it to the caller variable.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    matches=( "${tmp}/sample.bam" )
    assert_path_found matches BAM "${tmp}" found
    printf '%s\n' "${found}"
    rm -r -- "${tmp}"
    '''

  2. Reject multiple candidates and clear the caller variable.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    matches=( /private/tmp/a.bam /private/tmp/b.bam )
    found=stale
    assert_path_found matches "BAM output" /private/tmp found
    printf 'found=%s\n' "${found}"
    '''
EOM
    )

    if [[ "${arr_ref}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  run_samtools
    [--help] samtools_arg [samtools_arg...]

  Run Samtools from the active shell when available, or through the resolved project environment otherwise.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1+ samtools_arg : str
    Arguments passed through to 'samtools'.

Returns
-------
  Returns the exit status from the selected Samtools invocation.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing samtools (when 'samtools' is unavailable)
    - bash >= 4.4
    - conda (when 'samtools' is unavailable)
    - samtools

Examples
--------
  1. Print the Samtools version from the active shell or project environment.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    run_samtools --version
    '''

  2. Count reads in a committed BAM fixture.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    run_samtools view -c tests/fixtures/compute_signal/bam/se/tiny_se.bam
    '''
EOM
    )

    if [[ "${1:-}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if \
        check_cmd_exists samtools
    then
        samtools "$@"
    else
        conda run -n "${env_lcl}" samtools "$@"
    fi
}


#  Build and index a BAM input fixture from a committed SAM fixture
function build_filter_alignments_fixture_bam() {
    local in_sam="${1:-}"
    local out_bam="${2:-}"
    local log_lcl="${3:-}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  build_filter_alignments_fixture_bam
    [--help] in_sam out_bam log_lcl [label...]

  Build and index a BAM input fixture from a committed SAM fixture.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  in_sam : file
    Input SAM fixture.

  2  out_bam : file
    Output BAM fixture to create.

  3  log_lcl : file
    Log file for the conversion command.

  4+ label : str
    Optional assertion label.

Returns
-------
  Records a pass and returns 0 when conversion and indexing succeed; otherwise records a failure and returns 1.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing samtools (when 'samtools' is unavailable)
    - bash >= 4.4
    - conda (when 'samtools' is unavailable)
    - dirname
    - mkdir
    - samtools

Examples
--------
  1. Build and index a BAM from the committed filter-alignments SAM fixture.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    build_filter_alignments_fixture_bam tests/fixtures/filter_alignments/sam/filter_sc_sp.sam "${tmp}/filter.bam" "${tmp}/build.log"
    rm -r -- "${tmp}"
    '''

  2. Record a conversion failure for a missing SAM fixture.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    if ! \
        build_filter_alignments_fixture_bam \
            "${tmp}/missing.sam" \
            "${tmp}/filter.bam" \
            "${tmp}/build.log" \
            "missing-SAM BAM fixture"
    then
        printf '%s\n' 'BAM construction failed as expected'
    fi
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${in_sam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
function build_filter_alignments_fixture_cram() {
    local in_sam="${1:-}"
    local ref_fa="${2:-}"
    local out_cram="${3:-}"
    local log_lcl="${4:-}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  build_filter_alignments_fixture_cram
    [--help] in_sam ref_fa out_cram log_lcl [label...]

  Build and index a CRAM input fixture from committed SAM/reference fixtures.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  in_sam : file
    Input SAM fixture.

  2  ref_fa : file
    Reference FASTA file. Reference FASTA used for CRAM writing.

  3  out_cram : file
    Output CRAM fixture to create.

  4  log_lcl : file
    Log file for the conversion command.

  5+ label : str
    Optional assertion label.

Returns
-------
  Records a pass and returns 0 when conversion and indexing succeed; otherwise records a failure and returns 1.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing samtools (when 'samtools' is unavailable)
    - bash >= 4.4
    - conda (when 'samtools' is unavailable)
    - dirname
    - mkdir
    - Reference FASTA for CRAM writing
    - samtools

Examples
--------
  1. Build and index a CRAM with the committed reference FASTA.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    build_filter_alignments_fixture_cram tests/fixtures/filter_alignments/sam/filter_sc_sp.sam tests/fixtures/filter_alignments/reference/filter_sc_sp.fa "${tmp}/filter.cram" "${tmp}/build.log"
    rm -r -- "${tmp}"
    '''

  2. Record a CRAM conversion failure for a missing reference.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    if ! \
        build_filter_alignments_fixture_cram \
            tests/fixtures/filter_alignments/sam/filter_sc_sp.sam \
            "${tmp}/missing.fa" \
            "${tmp}/filter.cram" \
            "${tmp}/build.log" \
            "missing-reference CRAM fixture"
    then
        printf '%s\n' 'CRAM construction failed as expected'
    fi
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${in_sam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  find_python
    [--help]

  Find a usable Python command on PATH.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

Returns
-------
  Prints a Python executable path to stdout and returns 0 when found; otherwise returns 1.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - python >= 3.11 or python3 >= 3.11

Examples
--------
  1. Resolve the preferred Python command from the current PATH.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    find_python
    '''

  2. Handle a PATH with neither 'python' nor 'python3'.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    if ! \
        PATH=/nonexistent find_python
    then
        printf '%s\n' 'Python unavailable as expected'
    fi
    '''
EOM
    )

    if [[ "${1:-}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  find_python_loopback
    [--help]

  Find a Python command that can bind a local loopback socket.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

Returns
-------
  Prints a Python executable path to stdout and returns 0 when found; otherwise returns 1.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - python >= 3.11 or python3 >= 3.11, each with loopback socket support

Examples
--------
  1. Resolve a Python interpreter that can bind a loopback socket.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    find_python_loopback
    '''

  2. Allow the helper to use its explicit '/usr/bin/python3' fallback when PATH has no Python commands.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    if PATH=/nonexistent find_python_loopback; then
        printf '%s\n' 'system Python accepts loopback sockets'
    fi
    '''
EOM
    )

    if [[ "${1:-}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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


#  Find an available loopback TCP port for local HTTP tests
function find_port_free() {
    local py="${1:-}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  find_port_free
    [--help] py

  Find an available loopback TCP port using a Python interpreter.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  py : file
    Python executable to run.

Returns
-------
  Prints an available port number to stdout. Returns the Python command's exit status.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - Caller-supplied Python interpreter >= 3.11

Examples
--------
  1. Allocate a free loopback port with a discovered Python interpreter.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    py="$(find_python_loopback)"
    find_port_free "${py}"
    '''

  2. Propagate failure from an unusable interpreter path.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    if ! \
        find_port_free /bin/false
    then
        printf '%s\n' 'port discovery failed as expected'
    fi
    '''
EOM
    )

    if [[ "${py}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  wait_http_local
    [--help] py url [try]

  Wait until a local HTTP server returns status 200.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  py : file
    Python executable used for URL probing.

  2  url : str
    Local HTTP URL to probe.

  3  try : int
    Maximum number of polling attempts (default: 50).

Returns
-------
  Returns 0 when the URL becomes reachable; otherwise returns 1 after all attempts are exhausted.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - Caller-supplied Python interpreter >= 3.11 with 'urllib' support
    - sleep

Examples
--------
  1. Wait for a bounded loopback HTTP server and then clean it up.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    py="$(find_python_loopback)"
    port="$(find_port_free "${py}")"
    "${py}" -m http.server "${port}" --bind 127.0.0.1 > /dev/null 2>&1 &
    pid=$!
    trap 'cleanup_server_http "${pid}"' EXIT
    wait_http_local "${py}" "http://127.0.0.1:${port}/" 20
    cleanup_server_http "${pid}"
    trap - EXIT
    '''

  2. Report timeout when no server is listening on a newly released port.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    py="$(find_python_loopback)"
    port="$(find_port_free "${py}")"
    if ! \
        wait_http_local "${py}" "http://127.0.0.1:${port}/" 2
    then
        printf '%s\n' 'HTTP readiness timed out as expected'
    fi
    '''
EOM
    )

    if [[ "${py}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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


#  Check whether Python is at least version 3.11
function check_python_ge_311() {
    local py="${1:-}"

    "${py}" - << PY
import sys
raise SystemExit(0 if sys.version_info >= (3, 11) else 1)
PY
}


#  Run a command and capture stdout/stderr in a log file
function run_capture() {
    local nam="${1:-command}"
    local out="${2:-}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  run_capture
    [--help] nam out cmd [arg...]

  Run a command and capture stdout/stderr in a log file.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  nam : str
    Case name used when deriving a default log path.

  2  out : file
    Log file path. If empty, a sanitized path under TEST_DIR_LOG is used.

  3+ cmd : str
    Command and arguments to execute.

Returns
-------
  Returns the command exit status after writing stdout and stderr to the log.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - Caller-supplied command executable
    - dirname
    - mkdir

Examples
--------
  1. Capture successful stdout in an explicit log file.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    run_capture greeting "${tmp}/greeting.log" printf '%s\n' hello
    grep -q '^hello$' "${tmp}/greeting.log"
    rm -r -- "${tmp}"
    '''

  2. Capture stderr and retain the command's nonzero status.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    if ! \
        run_capture \
            rejection \
            "${tmp}/rejection.log" \
            sh -c 'printf "%s\n" rejected >&2; exit 3'
    then
        grep -q '^rejected$' "${tmp}/rejection.log"
    fi
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  run_case_filter
    [--help] retain nam_cas log_lcl arr_cmd_ref nam_job_pfx arr_arg_ref [label...]

  Run one filter-alignments retain-mode wrapper case.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  retain : str
    Species chromosomes to retain. Retain mode passed to the wrapper.

  2  nam_cas : str
    Case name suffix used for the job name.

  3  log_lcl : file
    Log file for the wrapper command.

  4  arr_cmd_ref : str
    Name of the command array to execute.

  5  nam_job_pfx : str
    Job name. Job-name prefix.

  6  arr_arg_ref : str
    Name of the array containing additional wrapper arguments.

  7+ label : str
    Optional assertion label.

Returns
-------
  Records a pass/fail assertion for the wrapper command. Returns 1 for missing array metadata or empty command arrays.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - dirname
    - mkdir

Examples
--------
  1. Run an execute-wrapper retain=sc case against a deterministic BAM fixture.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    require_env_project env_nam
    tmp="$(mktemp -d)"
    mkdir -p "${tmp}/input" "${tmp}/out" "${tmp}/logs"
    build_filter_alignments_fixture_bam tests/fixtures/filter_alignments/sam/filter_sc_sp.sam "${tmp}/input/filter.bam" "${tmp}/prepare.log"
    cmd=( "${TEST_BASH}" bin/execute_filter_alignments.sh --threads 1 --csv_fil_in "${tmp}/input/filter.bam" --dir_out "${tmp}/out" --dir_eo "${tmp}/logs" --max_job 1 )
    args=()
    run_case_filter sc canonical "${tmp}/filter-sc.log" cmd test_filter args
    rm -r -- "${tmp}"
    '''

  2. Run retain=sp while explicitly retaining supported optional contigs.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    require_env_project env_nam
    tmp="$(mktemp -d)"
    mkdir -p "${tmp}/input" "${tmp}/out" "${tmp}/logs"
    build_filter_alignments_fixture_bam tests/fixtures/filter_alignments/sam/filter_sc_sp.sam "${tmp}/input/filter.bam" "${tmp}/prepare.log"
    cmd=( "${TEST_BASH}" bin/execute_filter_alignments.sh --threads 1 --csv_fil_in "${tmp}/input/filter.bam" --dir_out "${tmp}/out" --dir_eo "${tmp}/logs" --max_job 1 )
    args=( --tg --mtr --mito )
    run_case_filter sp optional "${tmp}/filter-sp.log" cmd test_filter args && assert_file_nonempty "${tmp}/out/filter.sp.bam"
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${retain}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
                --env_nam "${env_nam}" \
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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  run_case_compute_signal
    [--help] wrap fmt_in nam_cas mode in_lcl out_spc log_lcl dir_out_lcl dir_err_lcl ref_fa_lcl [arg...]

  Run one compute-signal wrapper case for BAM or CRAM input.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  wrap : {'submit', 'execute'}
    Wrapper family to run.

  2  fmt_in : str
    Input format label used in case names.

  3  nam_cas : str
    Case name suffix.

  4  mode : str
    Workflow mode. Compute-signal mode to pass to the wrapper.

  5  in_lcl : file
    Input alignment file or serialized input string.

  6  out_spc : str
    Output file or output type spec, depending on wrapper family.

  7  log_lcl : file
    Log file for the wrapper command.

  8  dir_out_lcl : dir
    Output directory for execute-mode cases.

  9  dir_err_lcl : dir
    Log directory passed to the wrapper.

  10  ref_fa_lcl : file
    Reference FASTA file path.

Returns
-------
  Records a pass/fail assertion for the wrapper command.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - dirname
    - mkdir

Examples
--------
  1. Run the submit wrapper for signal mode with a BAM fixture.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    require_env_project env_nam
    tmp="$(mktemp -d)"
    mkdir -p "${tmp}/out" "${tmp}/logs"
    run_case_compute_signal submit bam se_signal signal tests/fixtures/compute_signal/bam/se/tiny_se.bam "${tmp}/out/signal.bdg" "${tmp}/submit.log" "${tmp}/out" "${tmp}/logs" '' --chr_sizes tests/fixtures/compute_signal/reference/tiny.fa.fai --method unadj --siz_bin 10 --engine window --csv_scl_fct NA --dp 3
    rm -r -- "${tmp}"
    '''

  2. Run the execute wrapper for coordinate mode with a CRAM fixture.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    require_env_project env_nam
    tmp="$(mktemp -d)"
    mkdir -p "${tmp}/out" "${tmp}/logs"
    run_case_compute_signal execute cram se_coord coord tests/fixtures/compute_signal/cram/se/tiny_se.cram bed "${tmp}/execute.log" "${tmp}/out" "${tmp}/logs" tests/fixtures/compute_signal/reference/tiny.fa
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${wrap}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    shift 10 || true

    local fmt_up="${fmt_in^^}"
    local script="${ROOT_REPO}/bin/${wrap}_compute_signal.sh"
    local nam_job=""
    local -a arr_cmd=()

    case "${wrap}" in
        submit)
            nam_job="test_compute_${fmt_in}_${nam_cas}"
            # shellcheck disable=SC2154
            arr_cmd=(
                "${TEST_BASH}" "${script}"
                    --env_nam "${env_nam}"
                    --dir_scr "${ROOT_REPO}/bin"
                    --threads 1
                    --mode "${mode}"
                    --csv_fil_in "${in_lcl}"
                    --csv_fil_out "${out_spc}"
            )
            ;;

        execute)
            nam_job="test_execute_compute_${fmt_in}_${nam_cas}"
            arr_cmd=(
                "${TEST_BASH}" "${script}"
                    --env_nam "${env_nam}"
                    --threads 1
                    --mode "${mode}"
                    --csv_fil_in "${in_lcl}"
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
        --dir_eo "${dir_err_lcl}"
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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  run_case_compute_signal_ratio
    [--help] wrap cas_nam out_spec log_lcl [method] fil_A fil_B dir_out_lcl dir_err_lcl [arg...]

  Run one compute-signal ratio wrapper case.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  wrap : {'submit', 'execute'}
    Wrapper family to run.

  2  cas_nam : str
    Case name suffix.

  3  out_spec : str
    Output file or prefix, depending on wrapper family.

  4  log_lcl : file
    Log file for the wrapper command.

  5  method : str
    Workflow method. Ratio method to request (default: unadj).

  6  fil_A : file
    First bedGraph input file, file A.

  7  fil_B : file
    Second bedGraph input file, file B.

  8  dir_out_lcl : dir
    Output directory for execute-mode cases.

  9  dir_err_lcl : dir
    Log directory passed to the wrapper.

Returns
-------
  Records a pass/fail assertion for the wrapper command.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - dirname
    - mkdir

Examples
--------
  1. Run an unadjusted submit-level ratio case with explicit output.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    require_env_project env_nam
    tmp="$(mktemp -d)"
    mkdir -p "${tmp}/out" "${tmp}/logs"
    run_case_compute_signal_ratio submit unadj "${tmp}/out/ratio.bdg" "${tmp}/submit.log" unadj tests/fixtures/compute_signal/bedgraph/ratio_A.bdg tests/fixtures/compute_signal/bedgraph/ratio_B.bdg "${tmp}/out" "${tmp}/logs" --csv_scl_fct NA --csv_dep_min NA --csv_pseudo NA --eps 0 --dp 3
    rm -r -- "${tmp}"
    '''

  2. Run an execute-level unadjusted ratio case with a track sidecar.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    require_env_project env_nam
    tmp="$(mktemp -d)"
    mkdir -p "${tmp}/out" "${tmp}/logs"
    run_case_compute_signal_ratio execute track ratio_track "${tmp}/execute.log" unadj tests/fixtures/compute_signal/bedgraph/ratio_A.bdg tests/fixtures/compute_signal/bedgraph/ratio_B.bdg "${tmp}/out" "${tmp}/logs" --track
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${wrap}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    shift 9 || true

    local script="${ROOT_REPO}/bin/${wrap}_compute_signal.sh"
    local -a arr_cmd=()

    case "${wrap}" in
        submit)
            # shellcheck disable=SC2154
            arr_cmd=(
                "${TEST_BASH}" "${script}"
                    --env_nam "${env_nam}"
                    --dir_scr "${ROOT_REPO}/bin"
                    --threads 1
                    --mode ratio
                    --method "${method}"
                    --csv_fil_A "${fil_A}"
                    --csv_fil_B "${fil_B}"
                    --csv_fil_out "${out_spec}"
                    --dir_eo "${dir_err_lcl}"
                    --nam_job "test_submit_compute_ratio_${cas_nam}"
            )
            ;;

        execute)
            arr_cmd=(
                "${TEST_BASH}" "${script}"
                    --env_nam "${env_nam}"
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
                    --dir_eo "${dir_err_lcl}"
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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  run_case_scaling_factor_execute
    [--help] cas mode arr_cmd_nam dir_out_lcl dir_log_lcl out_suffix job_prefix header row_0 row_1 tail_0 tail_1 [arg...]

  Run one execute-level calculate-scaling-factor two-row case and assert the generated final and part files.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  cas : str
    Case name.

  2  mode : str
    Workflow mode. Scaling-factor mode label.

  3  arr_cmd_nam : str
    Name of the base command array.

  4  dir_out_lcl : dir
    Output directory for generated TSV files.

  5  dir_log_lcl : dir
    Log directory for command output.

  6  out_suffix : str
    Output file path. Output filename suffix.

  7  job_prefix : str
    Job name. Job-name prefix.

  8  header : str
    Expected header pattern.

  9  row_0 : str
    Expected prefix for the first data row.

  10  row_1 : str
    Expected prefix for the second data row.

  11  tail_0 : str
    Expected suffix for the first data row.

  12  tail_1 : str
    Expected suffix for the second data row.

Returns
-------
  Records command, file-existence, header, and row-content assertions.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - dirname
    - grep
    - mkdir

Examples
--------
  1. Run a two-row spike-in case that retains the final header.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    require_env_project env_nam
    tmp="$(mktemp -d)"
    dir_fx=tests/fixtures/calculate_scaling_factor/bam/se
    mip_0="${dir_fx}/IP_WT_log_Brn1_rep1.sc.bam"
    mip_1="${dir_fx}/IP_WT_log_Brn1_rep2.sc.bam"
    min_0="${dir_fx}/in_WT_log_Brn1_rep1.sc.bam"
    min_1="${dir_fx}/in_WT_log_Brn1_rep2.sc.bam"
    sip_0="${dir_fx}/IP_WT_log_Brn1_rep1.sp.bam"
    sip_1="${dir_fx}/IP_WT_log_Brn1_rep2.sp.bam"
    sin_0="${dir_fx}/in_WT_log_Brn1_rep1.sp.bam"
    sin_1="${dir_fx}/in_WT_log_Brn1_rep2.sp.bam"
    cmd=( "${TEST_BASH}" bin/execute_calculate_scaling_factor.sh --threads 1 --mode spike --aln_typ auto --csv_mip "${mip_0},${mip_1}" --csv_min "${min_0},${min_1}" --csv_sip "${sip_0},${sip_1}" --csv_sin "${sin_0},${sin_1}" --dir_eo "${tmp}/logs" --max_job 1 )
    row_0="${mip_0}"$'\t'"${sip_0}"$'\t'"${min_0}"$'\t'"${sin_0}"
    row_1="${mip_1}"$'\t'"${sip_1}"$'\t'"${min_1}"$'\t'"${sin_1}"
    run_case_scaling_factor_execute default spike cmd "${tmp}/out" "${tmp}/logs" spike test_spike $'^main_ip\tspike_ip\tmain_in\tspike_in\tspike\tcoef\tnum_mp\tnum_sp\tnum_mn\tnum_sn$' "${row_0}" "${row_1}" $'2\tchiprx_alpha_ratio\t3\t1\t2\t2' $'0.5\tchiprx_alpha_ratio\t2\t2\t3\t1'
    rm -r -- "${tmp}"
    '''

  2. Run a two-row siQ-ChIP case without adding a final-table header.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    require_env_project env_nam
    tmp="$(mktemp -d)"
    dir_fx=tests/fixtures/calculate_scaling_factor/bam/pe
    mip_0="${dir_fx}/IP_WT_G1_Hho1_6336.sc.bam"
    mip_1="${dir_fx}/IP_WT_G1_Hho1_6337.sc.bam"
    min_0="${dir_fx}/in_WT_G1_Hho1_6336.sc.bam"
    min_1="${dir_fx}/in_WT_G1_Hho1_6337.sc.bam"
    cmd=( "${TEST_BASH}" bin/execute_calculate_scaling_factor.sh --threads 1 --mode siq --aln_typ auto --csv_mip "${mip_0},${mip_1}" --csv_min "${min_0},${min_1}" --tbl_met tests/fixtures/calculate_scaling_factor/metadata/measurements_siqchip.tsv --cfg_met data/raw/docs/parse_metadata_siqchip.yml --eqn 6nd --dir_eo "${tmp}/logs" --max_job 1 )
    row_0="${mip_0}"$'\t'"${min_0}"
    row_1="${mip_1}"$'\t'"${min_1}"
    run_case_scaling_factor_execute pe_bam siQ cmd "${tmp}/out" "${tmp}/logs" siq test_siq $'^fil_ip\tfil_in\tsiq\teqn\tmass_ip\tmass_in\tvol_all\tvol_in\tdep_ip\tdep_in\tlen_ip\tlen_in$' "${row_0}" "${row_1}" $'0.001912211397724232486706\t6nd\t2.7\t72.5\t300\t20\t3\t2\t626\t450' $'0.002902612244746139314594\t6nd\t5\t81.1\t300\t20\t2\t3\t663\t437' --no_header
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${cas}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
        --env_nam "${env_nam}"
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
        record_pass \
            "execute_calculate_scaling_factor.sh ${mode} ${cas} exits 0"
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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  run_case_scaling_factor_submit_part
    [--help] cas mode arr_cmd_nam fil_out idx_out fil_log row_exp header [arg...]

  Run one submit-level calculate-scaling-factor one-part case and assert the generated part file.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  cas : str
    Case name.

  2  mode : str
    Workflow mode. Scaling-factor mode label.

  3  arr_cmd_nam : str
    Name of the base command array.

  4  fil_out : file
    Output file path. Final TSV path used to derive the part-file path.

  5  idx_out : int
    Zero-based part index expected from the submit wrapper.

  6  fil_log : file
    Log file for command output.

  7  row_exp : str
    Expected complete part-file row.

  8  header : str
    Header pattern that should be absent from the part file.

Returns
-------
  Records command, part-file, final-file, and row-content assertions.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - cat
    - dirname
    - grep
    - mkdir
    - wc

Examples
--------
  1. Run one spike-in submit part and require a data-only row.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    require_env_project env_nam
    tmp="$(mktemp -d)"
    dir_fx=tests/fixtures/calculate_scaling_factor/bam/se
    mip="${dir_fx}/IP_WT_log_Brn1_rep1.sc.bam"
    min="${dir_fx}/in_WT_log_Brn1_rep1.sc.bam"
    sip="${dir_fx}/IP_WT_log_Brn1_rep1.sp.bam"
    sin="${dir_fx}/in_WT_log_Brn1_rep1.sp.bam"
    out="${tmp}/scaling.tsv"
    cmd=( "${TEST_BASH}" bin/submit_calculate_scaling_factor.sh --env_nam "${env_nam}" --dir_scr scripts --threads 1 --mode spike --aln_typ auto --csv_mip "${mip}" --csv_min "${min}" --csv_sip "${sip}" --csv_sin "${sin}" --fil_out "${out}" --idx_out 0 --dir_eo "${tmp}/logs" --nam_job test_spike )
    row="${mip}"$'\t'"${sip}"$'\t'"${min}"$'\t'"${sin}"$'\t2\tchiprx_alpha_ratio\t3\t1\t2\t2'
    run_case_scaling_factor_submit_part se_default spike cmd "${out}" 0 "${tmp}/submit.log" "${row}" $'^main_ip\tspike_ip\tmain_in\tspike_in\tspike\tcoef'
    rm -r -- "${tmp}"
    '''

  2. Run one siQ-ChIP submit part using gzip-compressed metadata.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    require_env_project env_nam
    tmp="$(mktemp -d)"
    dir_fx=tests/fixtures/calculate_scaling_factor/bam/pe
    mip="${dir_fx}/IP_WT_G1_Hho1_6336.sc.bam"
    min="${dir_fx}/in_WT_G1_Hho1_6336.sc.bam"
    out="${tmp}/scaling.tsv"
    cmd=( "${TEST_BASH}" bin/submit_calculate_scaling_factor.sh --env_nam "${env_nam}" --dir_scr scripts --threads 1 --mode siq --aln_typ auto --csv_mip "${mip}" --csv_min "${min}" --fil_out "${out}" --idx_out 19 --tbl_met tests/fixtures/calculate_scaling_factor/metadata/measurements_siqchip.tsv.gz --cfg_met data/raw/docs/parse_metadata_siqchip.yml --eqn 6nd --dir_eo "${tmp}/logs" --nam_job test_siq )
    row="${mip}"$'\t'"${min}"$'\t0.001912211397724232486706\t6nd\t2.7\t72.5\t300\t20\t3\t2\t626\t450'
    run_case_scaling_factor_submit_part pe_gzip siQ cmd "${out}" 19 "${tmp}/submit.log" "${row}" $'^fil_ip\tfil_in\tsiq\teqn\tmass_ip\tmass_in' && assert_file_nonempty "${out}.part.000019"
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${cas}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    shift 8 || true

    local -n arr_cmd_ref="${arr_cmd_nam}"

    local idx_fmt=""
    local fil_prt=""

    printf -v idx_fmt '%06d' "${idx_out}"
    fil_prt="${fil_out}.part.${idx_fmt}"

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
        record_pass \
            "submit_calculate_scaling_factor.sh ${mode} ${cas} writes no" \
            "final TSV"
    else
        record_fail \
            "submit_calculate_scaling_factor.sh ${mode} ${cas} wrote final TSV"
    fi

    assert_file_exact_line \
        "${fil_prt}" \
        "${row_exp}" \
        "submit scaling-factor ${mode} ${cas} part"
}


#  Assert that a log file contains a pattern
function assert_pattern_found() {
    local file="${1:-}"
    local patn="${2:-}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  assert_pattern_found
    [--help] file patn [label...]

  Assert that a file contains a grep-compatible pattern.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  file : file
    File to inspect.

  2  patn : str
    Pattern passed to 'grep -q --'.

  3+ label : str
    Optional assertion label.

Returns
-------
  Records a pass when the pattern is found; otherwise records a failure.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - grep

Examples
--------
  1. Record a pass when a generated log contains its success marker.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    printf '%s\n' 'workflow complete' > "${tmp}/run.log"
    assert_pattern_found "${tmp}/run.log" '^workflow complete$'
    rm -r -- "${tmp}"
    '''

  2. Record a labeled failure when an expected marker is absent.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    printf '%s\n' 'workflow started' > "${tmp}/run.log"
    assert_pattern_found "${tmp}/run.log" '^workflow complete$' "completion marker"
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${file}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  assert_pattern_absent
    [--help] file patn [label...]

  Assert that a file does not contain a grep-compatible pattern.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  file : file
    File to inspect.

  2  patn : str
    Pattern passed to 'grep -q --'.

  3+ label : str
    Optional assertion label.

Returns
-------
  Records a pass when the pattern is absent; otherwise records a failure.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - grep

Examples
--------
  1. Record a pass when a data-only part file has no header.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    printf '%s\n' $'sample_A\t1.0' > "${tmp}/part.tsv"
    assert_pattern_absent "${tmp}/part.tsv" '^sample_name'
    rm -r -- "${tmp}"
    '''

  2. Record a labeled failure when a forbidden pattern is present.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    printf '%s\n' $'sample_name\tscaling_factor' > "${tmp}/part.tsv"
    assert_pattern_absent "${tmp}/part.tsv" '^sample_name' "part stays data-only"
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${file}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  assert_file_nonempty
    [--help] file [label...]

  Assert that a generated output file exists and is non-empty.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  file : file
    File to inspect.

  2+ label : str
    Optional assertion label.

Returns
-------
  Records a pass when the file is non-empty; otherwise records a failure.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Record a pass for a generated non-empty bedGraph.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    printf '%s\n' $'I\t0\t10\t1' > "${tmp}/signal.bdg"
    assert_file_nonempty "${tmp}/signal.bdg"
    rm -r -- "${tmp}"
    '''

  2. Record a labeled failure for an empty output file.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    : > "${tmp}/signal.bdg"
    assert_file_nonempty "${tmp}/signal.bdg" "signal output"
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${file}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  assert_file_exact_line
    [--help] file expected [label...]

  Assert that a file contains exactly one line matching an expected string.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  file : file
    File to inspect.

  2  expected : str
    Exact expected line.

  3+ label : str
    Optional assertion label.

Returns
-------
  Records row-count and content assertions. Returns 1 immediately when the file is missing or empty.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - cat
    - wc

Examples
--------
  1. Accept a one-row submit part with the expected content.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    printf '%s\n' $'sample_A\t1.0' > "${tmp}/part.tsv"
    assert_file_exact_line "${tmp}/part.tsv" $'sample_A\t1.0'
    rm -r -- "${tmp}"
    '''

  2. Record labeled row-count and content failures for a two-row file.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    printf '%s\n' first second > "${tmp}/part.tsv"
    assert_file_exact_line "${tmp}/part.tsv" first "submit part"
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${file}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
        record_fail \
            "${lbl} has unexpected row count; see $(print_relpath "${file}")"
    fi

    observed="$(cat "${file}")"
    if [[ "${observed}" == "${expected}" ]]; then
        record_pass "${lbl} has expected row"
    else
        record_fail \
            "${lbl} differs from expected; see $(print_relpath "${file}")"
    fi
}


#  Assert that two files have exactly identical contents
function assert_files_equal() {
    local observed="${1:-}"
    local expected="${2:-}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  assert_files_equal
    [--help] observed expected [label...]

  Assert that two files have exactly identical contents.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  observed : file
    Generated or observed file.

  2  expected : file
    Expected fixture file.

  3+ label : str
    Optional assertion label.

Returns
-------
  Records a pass when files match byte-for-byte; otherwise records a failure.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - cmp

Examples
--------
  1. Accept two byte-identical generated files.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    printf '%s\n' expected > "${tmp}/expected.txt"
    cp "${tmp}/expected.txt" "${tmp}/observed.txt"
    assert_files_equal "${tmp}/observed.txt" "${tmp}/expected.txt"
    rm -r -- "${tmp}"
    '''

  2. Record a labeled failure for files with different contents.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    printf '%s\n' observed > "${tmp}/observed.txt"
    printf '%s\n' expected > "${tmp}/expected.txt"
    assert_files_equal "${tmp}/observed.txt" "${tmp}/expected.txt" "normalized output"
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${observed}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    shift 2 || true

    local lbl="${*:-file content}"

    if \
        cmp -s "${observed}" "${expected}"
    then
        record_pass "${lbl}"
    else
        record_fail "${lbl}; see $(print_relpath "${observed}")"
    fi
}


#  Assert whether a scaling-factor output contains an expected header
function assert_scaling_factor_header() {
    local file="${1:-}"
    local header="${2:-}"
    local expect="${3:-true}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  assert_scaling_factor_header
    [--help] file header [expect] [label...]

  Assert whether a scaling-factor output contains an expected header.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  file : file
    Scaling-factor TSV file.

  2  header : str
    Header pattern to search for.

  3  expect : bool
    Whether the header should be present (default: true).

  4+ label : str
    Optional assertion label.

Returns
-------
  Delegates to pattern-present or pattern-absent assertions and records their result.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - grep

  Accepted true-like values: true, t, yes, y, 1. Accepted false-like values: false, f, no, n, 0. Matching is case-insensitive; surrounding whitespace and empty values are invalid. Successful normalization emits 'true' or 'false'.

Examples
--------
  1. Require the core header in a standard scaling-factor TSV.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    printf '%s\n' $'sample_name\tscaling_factor' > "${tmp}/scaling.tsv"
    assert_scaling_factor_header "${tmp}/scaling.tsv" '^sample_name' true
    rm -r -- "${tmp}"
    '''

  2. Require a no-header output to omit the core header.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    printf '%s\n' $'sample_A\t1.0' > "${tmp}/scaling.tsv"
    assert_scaling_factor_header "${tmp}/scaling.tsv" '^sample_name' false "no-header TSV"
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${file}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    expect="$(normalize_bool "${expect}" "expect")" || return 1

    shift 3 || true

    local lbl="${*:-scaling-factor header}"

    if [[ "${expect}" == "true" ]]; then
        assert_pattern_found \
            "${file}" \
            "${header}" \
            "${lbl} has core header"
    else
        assert_pattern_absent \
            "${file}" \
            "${header}" \
            "${lbl} omits core header"
    fi
}


#  Assert that a CRAM index exists using common Samtools naming variants
function assert_cram_index() {
    local cram="${1:-}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  assert_cram_index
    [--help] cram [label...]

  Assert that a CRAM index exists using common Samtools naming variants.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  cram : file
    CRAM path whose index should exist.

  2+ label : str
    Optional assertion label.

Returns
-------
  Records a pass when either '<cram>.crai' or '<stem>.crai' exists and is non-empty; otherwise records a failure.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Accept the '<cram>.crai' naming convention emitted by Samtools.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    : > "${tmp}/sample.cram"
    printf '%s\n' index > "${tmp}/sample.cram.crai"
    assert_cram_index "${tmp}/sample.cram"
    rm -r -- "${tmp}"
    '''

  2. Accept the alternate '<stem>.crai' naming convention.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    : > "${tmp}/sample.cram"
    printf '%s\n' index > "${tmp}/sample.crai"
    assert_cram_index "${tmp}/sample.cram" "alternate CRAM index"
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${cram}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  assert_cram_count
    [--help] cram ref_fa exp fil_cnt arr_ref [label...]

  Assert a reference-backed CRAM read count for one or more regions.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  cram : file
    CRAM file to count.

  2  ref_fa : file
    Reference FASTA file for CRAM decoding.

  3  exp : int
    Expected read count.

  4  fil_cnt : file
    File that receives the observed count.

  5  arr_ref : str
    Name of the indexed array containing regions to count.

  6+ label : str
    Optional assertion label.

Returns
-------
  Records failures for missing region metadata and records a pattern assertion for the observed count.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing samtools (when 'samtools' is unavailable)
    - bash >= 4.4
    - conda (when 'samtools' is unavailable)
    - dirname
    - grep
    - mkdir
    - Reference FASTA for CRAM decoding
    - samtools

Examples
--------
  1. Assert the read count for one chromosome in a committed CRAM fixture.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    regions=( I )
    assert_cram_count tests/fixtures/compute_signal/cram/se/tiny_se.cram tests/fixtures/compute_signal/reference/tiny.fa 2 "${tmp}/count.txt" regions
    rm -r -- "${tmp}"
    '''

  2. Assert a combined count across two requested regions.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    regions=( I:1-10 I:21-30 )
    assert_cram_count tests/fixtures/compute_signal/cram/se/tiny_se.cram tests/fixtures/compute_signal/reference/tiny.fa 2 "${tmp}/count.txt" regions "two-region CRAM count"
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${cram}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local fil_in="${1:-}"
    local ref_fa="${2:-}"
    local pg_id="${3:-}"
    local retain="${4:-}"
    local out_ext="${5:-}"
    local header_file="${6:-}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  assert_filter_alignments_pg_header
    [--help] fil_in ref_fa pg_id retain out_ext header_file [label...]

  Assert filter_alignments @PG provenance in a BAM or CRAM header.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fil_in : file
    Input file path. BAM or CRAM file to inspect.

  2  ref_fa : file
    Reference FASTA file. Reference FASTA required for CRAM header reads.

  3  pg_id : str
    Expected @PG ID.

  4  retain : str
    Species chromosomes to retain. Expected retain value in the @PG command line.

  5  out_ext : str
    Final output extension. Expected output extension recorded in the @PG command line.

  6  header_file : file
    File that receives the extracted header.

  7+ label : str
    Optional assertion label.

Returns
-------
  Records header-read, @PG-presence, and out_ext assertions.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing samtools (when 'samtools' is unavailable)
    - bash >= 4.4
    - conda (when 'samtools' is unavailable)
    - dirname
    - grep
    - mkdir
    - Reference FASTA (when inspecting CRAM)
    - samtools

Examples
--------
  1. Inspect provenance in a BAM output without a reference argument.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    awk 'BEGIN { pg = "@PG\tID:filter_alignments\tPN:filter_alignments\tCL:filter_alignments retain=sc out_ext=bam" } /^@/ { print; next } !seen { print pg; seen = 1 } { print }' tests/fixtures/filter_alignments/sam/filter_sc_sp.sam > "${tmp}/provenance.sam"
    run_samtools view -bS -o "${tmp}/filtered.bam" "${tmp}/provenance.sam"
    assert_filter_alignments_pg_header "${tmp}/filtered.bam" '' filter_alignments sc bam "${tmp}/header.sam"
    rm -r -- "${tmp}"
    '''

  2. Inspect provenance in a CRAM output with its decoding reference.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    awk 'BEGIN { pg = "@PG\tID:filter_alignments\tPN:filter_alignments\tCL:filter_alignments retain=sp out_ext=cram" } /^@/ { print; next } !seen { print pg; seen = 1 } { print }' tests/fixtures/filter_alignments/sam/filter_sc_sp.sam > "${tmp}/provenance.sam"
    run_samtools view -C -T tests/fixtures/filter_alignments/reference/filter_sc_sp.fa -o "${tmp}/filtered.cram" "${tmp}/provenance.sam"
    assert_filter_alignments_pg_header "${tmp}/filtered.cram" tests/fixtures/filter_alignments/reference/filter_sc_sp.fa filter_alignments sp cram "${tmp}/header.sam" "CRAM provenance"
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${fil_in}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    shift 6 || true

    local lbl="${*:-filter_alignments @PG header}"
    local patn_pg=""
    local -a arr_ref_arg=()

    patn_pg=$'^@PG\tID:'"${pg_id}"$'\tPN:filter_alignments\tCL:'
    patn_pg+="${pg_id} retain=${retain}"

    if [[ "${fil_in,,}" == *.cram ]]; then
        arr_ref_arg=( -T "${ref_fa}" )
    fi

    if \
        run_capture \
            "header ${lbl}" "${header_file}" \
            run_samtools view -H "${arr_ref_arg[@]}" "${fil_in}"
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
            "${patn_pg}" \
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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  assert_fastq_gzip
    [--help] fil_out rd_patn fq_vw rd_cnt fq_src [label...]

  Assert a gzip FASTQ exists, passes integrity checks, and has expected read content.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fil_out : file
    Output file path. Gzip FASTQ output to inspect.

  2  rd_patn : str
    Expected read-name pattern.

  3  fq_vw : file
    Temporary decompressed FASTQ view.

  4  rd_cnt : file
    Temporary file that receives the read count.

  5  fq_src : file
    Optional source fixture for byte-for-byte comparison.

  6+ label : str
    Optional assertion label.

Returns
-------
  Records file, gzip integrity, decompression, read-name, and read-count assertions. Returns 0 immediately if 'fil_out' is empty.

Notes
-----
  Runtime requirements:
    - awk
    - bash >= 4.4
    - cmp (when a source fixture is supplied)
    - grep
    - gzip

Examples
--------
  1. Validate a gzip FASTQ and compare it with its source fixture byte-for-byte.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    cp tests/fixtures/download_fastqs/source/se/tiny_download_se.fastq.gz "${tmp}/read.fastq.gz"
    assert_fastq_gzip "${tmp}/read.fastq.gz" '^@tiny_download_se_read_1$' "${tmp}/read.fastq" "${tmp}/count.txt" tests/fixtures/download_fastqs/source/se/tiny_download_se.fastq.gz
    rm -r -- "${tmp}"
    '''

  2. Validate read content without requesting a byte-for-byte source comparison.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    cp tests/fixtures/trim_fastqs/fastq/se/tiny_se.fastq.gz "${tmp}/trimmed.fastq.gz"
    assert_fastq_gzip "${tmp}/trimmed.fastq.gz" '^@tiny_trim_se_read_1$' "${tmp}/trimmed.fastq" "${tmp}/count.txt" '' "trimmed FASTQ"
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${fil_out}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  assert_custom_symlink
    [--help] sym [label...]

  Assert that a custom FASTQ path exists and is represented as a symlink.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  sym : path
    Symlink path to inspect.

  2+ label : str
    Optional assertion label.

Returns
-------
  Records file-existence and symlink assertions.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Accept a custom FASTQ symlink whose target is non-empty.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    printf '%s\n' '@read' > "${tmp}/source.fastq"
    ln -s "${tmp}/source.fastq" "${tmp}/custom.fastq"
    assert_custom_symlink "${tmp}/custom.fastq"
    rm -r -- "${tmp}"
    '''

  2. Record a labeled failure when the custom path is a regular file.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    printf '%s\n' '@read' > "${tmp}/custom.fastq"
    assert_custom_symlink "${tmp}/custom.fastq" "custom read path"
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${sym}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  warn_help_pattern_missing
    [--help] file patn [label...]

  Warn when a rendered help-output section is missing.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  file : file
    Help-output capture to inspect.

  2  patn : str
    Pattern passed to 'grep -q --'.

  3+ label : str
    Optional warning label.

Returns
-------
  Records a pass when the pattern is present; otherwise records a warning.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - grep

Examples
--------
  1. Record a pass when rendered help contains the requested section.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    printf '%s\n' Usage ----- > "${tmp}/help.txt"
    warn_help_pattern_missing "${tmp}/help.txt" '^Usage$'
    rm -r -- "${tmp}"
    '''

  2. Record a labeled warning when rendered help omits Examples.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    tmp="$(mktemp -d)"
    printf '%s\n' Usage ----- > "${tmp}/help.txt"
    warn_help_pattern_missing "${tmp}/help.txt" '^Examples$' "full help has Examples"
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${file}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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


#  Print and return the final status for one test group
function finish() {
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  finish
    [--help]

  Print and return the final status for one test group.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

Returns
-------
  Prints the pass/fail/warn/skip summary. Returns 1 when any failures were recorded; otherwise returns 0.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Finish a smoke group with passes and skips but no failures.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    TEST_NAME='example smoke group'
    TEST_PASS=2
    TEST_SKIP=1
    finish
    '''

  2. Propagate a nonzero final status when failures were recorded.
    '''bash
    # shellcheck source=tests/support/test_helpers.sh
    source tests/support/test_helpers.sh
    TEST_NAME='failing smoke group'
    TEST_FAIL=1
    if ! \
        finish
    then
        printf '%s\n' 'smoke group failed as expected'
    fi
    '''
EOM
    )

    if [[ "${1:-}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    printf \
        '\nSummary for %s: pass=%d fail=%d warn=%d skip=%d\n' \
        "${TEST_NAME}" "${TEST_PASS}" "${TEST_FAIL}" "${TEST_WARN}" \
        "${TEST_SKIP}"

    if (( TEST_FAIL > 0 )); then
        return 1
    fi

    return 0
}
