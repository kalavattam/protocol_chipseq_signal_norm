#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: manage_parallel.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in development and documentation.
#
# Distributed under the MIT license.


# determine_cores
# print_parallel_info
# reset_max_job
# set_params_parallel


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

#  Source required helper functions if needed
{
    _dir_src_parl="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    # shellcheck source=lib/bash/core/source_helpers.sh
    source "${_dir_src_parl}/../core/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_parl}/../core/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_parl}/.." \
        check_args check_inputs check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_parl
}

#MAYBE: make function "private"?
# shellcheck disable=SC2120
function determine_cores() {
    local cores=""   # Number of CPU cores available on system
    local show_help  # Help message

    show_help=$(cat << EOM
Usage
-----
  determine_cores
    [--help]

  Determine the number of CPU cores available on the current system.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

Returns
-------
  Prints a positive integer core count to stdout; otherwise returns 1.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - nproc, sysctl, or getconf

  - Uses 'nproc' when available and valid.
  - Otherwise falls back to 'sysctl -n hw.ncpu' when available and valid.
  - Otherwise falls back to 'getconf _NPROCESSORS_ONLN' when available and valid.

Examples
--------
  1. Print the detected positive core count.
    '''bash
    determine_cores
    '''

  2. Confirm that an unexpected positional argument is rejected.
    '''bash
    if ! determine_cores unexpected; then
        printf '%s\n' 'unexpected argument rejected as expected'
    fi
    '''
EOM
    )

    if [[ "${1:-}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif (( $# > 0 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "this function does not accept positional arguments."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if \
        command -v nproc > /dev/null 2>&1
    then
        cores="$(nproc 2> /dev/null || true)"
    fi

    if ! \
        [[ "${cores}" =~ ^[1-9][0-9]*$ ]] \
        && command -v sysctl > /dev/null 2>&1
    then
        cores="$(sysctl -n hw.ncpu 2> /dev/null || true)"
    fi

    if ! \
        [[ "${cores}" =~ ^[1-9][0-9]*$ ]] \
        && command -v getconf > /dev/null 2>&1
    then
        cores="$(getconf _NPROCESSORS_ONLN 2> /dev/null || true)"
    fi

    if ! [[ "${cores}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "unable to determine a valid CPU core count using 'nproc'," \
            "'sysctl -n hw.ncpu', or 'getconf _NPROCESSORS_ONLN'."
        return 1
    fi

    echo "${cores}"
}


function print_parallel_info() {
    local slurm="${1:-}"    # Boolean-like flag for Slurm: 'true' or 'false'
    local max_job="${2:-}"  # No. concurrent jobs: Slurm
    local par_job="${3:-}"  # No. concurrent jobs: GNU Parallel or serial
    local threads="${4:-}"  # No. threads per job
    local show_help         # Help message

    show_help=$(cat << EOM
Usage
-----
  print_parallel_info
    [--help] slurm max_job par_job threads [arr1 arr2 ...]

  Print parsed array contents and resolved parallelization settings.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  slurm : bool
    Submit jobs to the Slurm scheduler. Whether Slurm mode is active.

  2  max_job : int
    Maximum number of jobs to run concurrently. Maximum concurrent jobs for Slurm mode.

  3  par_job : int
    Maximum concurrent jobs for GNU Parallel or serial mode.

  4  threads : int
    Number of threads to use. Threads per job.

  5+ arr_nam : str
    Optional names of indexed arrays to print via 'debug_arr_contents'.

Returns
-------
  0 after printing the summary; otherwise 1.

Notes
-----
  Runtime requirements:
    bash >= 4.4

  Accepted true-like values: true, t, yes, y, 1. Accepted false-like values: false, f, no, n, 0. Matching is case-insensitive; surrounding whitespace and empty values are invalid. Successful normalization emits 'true' or 'false'.

Examples
--------
  1. Report a serial local configuration using four threads per job.
    '''bash
    print_parallel_info false 0 1 4
    '''

  2. Report a Slurm configuration and include one parsed job array.
    '''bash
    declare -a arr_jobs=(sample_A sample_B)
    print_parallel_info true 8 1 4 arr_jobs
    '''
EOM
    )

    if [[ "${slurm}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${slurm}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'slurm', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${max_job}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'max_job', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${par_job}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'par_job', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${threads}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 4, 'threads', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    slurm="$(normalize_bool "${slurm}" "slurm")" || return 1

    if ! [[ "${max_job}" =~ ^[0-9]+$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'max_job', must be a non-negative integer:" \
            "'${max_job}'."
        return 1
    fi

    if ! [[ "${par_job}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'par_job', must be a positive integer:" \
            "'${par_job}'."
        return 1
    fi

    if ! [[ "${threads}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 4, 'threads', must be a positive integer:" \
            "'${threads}'."
        return 1
    fi

    shift 4

    echo "########################################################"
    echo "## Parsed vector(s) and arguments for parallelization ##"
    echo "########################################################"
    echo

    if (( $# > 0 )); then debug_arr_contents "$@"; fi

    if [[ "${slurm}" == "true" ]]; then
        echo "  - Max concurrent jobs (Slurm): ${max_job}"
    elif [[ "${par_job}" -gt 1 ]]; then
        echo "  - Max concurrent jobs (GNU Parallel): ${par_job}"
    else
        echo "  - Jobs running in serial mode: ${par_job}"
    fi

    echo "  - Number of threads to use. Threads per job: ${threads}"
    echo
    echo
}


function reset_max_job() {
    local max_job="${1:-}"
    local num_fil="${2:-}"
    local show_help

    show_help=$(cat << EOM
Usage
-----
  reset_max_job
    [--help] max_job num_fil

  Cap the requested Slurm maximum concurrent job count at the number of input files.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  max_job : int
    Maximum number of jobs to run concurrently. Requested maximum number of concurrent jobs.

  2  num_fil : int
    Number of input files/jobs available.

Returns
-------
  Prints the resolved maximum concurrent job count to stdout; otherwise 1.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Cap eight requested jobs at three available inputs.
    '''bash
    reset_max_job 8 3
    '''

  2. Confirm that a negative requested job count is rejected.
    '''bash
    if ! reset_max_job -1 5; then
        printf '%s\n' 'negative maximum rejected as expected'
    fi
    '''
EOM
    )

    if [[ "${max_job}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${max_job}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'max_job', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${num_fil}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'num_fil', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! [[ "${max_job}" =~ ^[0-9]+$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'max_job', must be a non-negative" \
            "integer: '${max_job}'."
        return 1
    fi

    if ! [[ "${num_fil}" =~ ^[0-9]+$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'num_fil', must be a non-negative" \
            "integer: '${num_fil}'."
        return 1
    fi

    if [[ "${max_job}" -gt "${num_fil}" ]]; then
        echo "${num_fil}"
    else
        echo "${max_job}"
    fi
}


function set_params_parallel() {
    local threads="${1:-}"  # Total thread budget (not per-job threads)
    local max_job="${2:-}"  # Maximum number of jobs allowed
    local n_cores           # Total available CPU cores on system
    local par_job           # Resolved GNU Parallel job count
    local show_help

    show_help=$(cat << EOM
Usage
-----
  set_params_parallel
    [--help] threads max_job

  For the non-Slurm path, determine a safe combination of threads per job and GNU Parallel job count from a requested total thread budget and maximum job count.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  threads : int
    Number of threads to use. Requested total local CPU/thread budget.

  2  max_job : int
    Maximum number of jobs to run concurrently.

Returns
-------
  Prints '<threads_per_job>;<par_job>' to stdout; otherwise 1.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - nproc, sysctl, or getconf

  - Intended for the non-Slurm / GNU Parallel execution path.
  - Interprets 'threads' as a total local CPU/thread budget, not as per-job threads.

Examples
--------
  1. Request serial execution with the full local thread budget.
    '''bash
    set_params_parallel 8 1
    '''

  2. Request four concurrent jobs from a sixteen-thread budget.
    '''bash
    if params="\$(set_params_parallel 16 4)"; then
        printf 'threads_per_job;parallel_jobs=%s\n' "\${params}"
    fi
    '''
EOM
    )

    #  Handle help requests and missing required positional arguments
    if [[ "${threads}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${threads}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'threads', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${max_job}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'max_job', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Check integer inputs
    if ! [[ "${threads}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'threads', must be a positive integer:" \
            "'${threads}'."
        return 1
    fi

    if ! [[ "${max_job}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'max_job', must be a positive integer:" \
            "'${max_job}'."
        return 1
    fi

    #  Get system CPU core count
    n_cores=$(determine_cores) || return 1

    #  Do not allow requested total threads to exceed available cores
    if [[ "${threads}" -gt "${n_cores}" ]]; then
        echo_warn_func "${FUNCNAME[0]}" \
            "requested threads ('${threads}') exceed available cores" \
            "('${n_cores}'). Using '${n_cores}'."
        threads="${n_cores}"
    fi

    #  If only one job is allowed, run serially
    if [[ "${max_job}" -le 1 ]]; then
        par_job=1
        echo "${threads};${par_job}"
        return 0
    fi

    #  Start from requested max parallel job count, but do not exceed total
    #+ available cores
    par_job="${max_job}"
    if [[ "${par_job}" -gt "${n_cores}" ]]; then
        echo_warn_func "${FUNCNAME[0]}" \
            "requested parallel jobs ('${par_job}') exceed available cores" \
            "('${n_cores}'). Using '${n_cores}'."
        par_job="${n_cores}"
    fi

    #  Convert total thread budget into per-job thread count, enforcing a
    #+ minimum of 1 thread per job
    if [[ "${threads}" -lt "${par_job}" ]]; then
        echo_warn_func "${FUNCNAME[0]}" \
            "requested total threads ('${threads}') are fewer than parallel" \
            "jobs ('${par_job}'). Using 1 thread per job."
        threads=1
    else
        threads=$(( threads / par_job ))
    fi

    echo "${threads};${par_job}"  # Return adjusted threads-per-job, job count
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
