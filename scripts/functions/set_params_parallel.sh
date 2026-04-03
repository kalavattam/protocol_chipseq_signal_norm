#!/bin/bash
# -*- coding: utf-8 -*-
# 
# Script: set_params_parallel.sh
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5-series models) was used in development.
# 
# Distributed under the MIT license.

#  Determine number of available CPU cores
function determine_cores() {
    local cores  # Number of CPU cores available on system

    if command -v nproc &> /dev/null; then
        cores=$(nproc)
    elif command -v sysctl &> /dev/null; then
        cores=$(sysctl -n hw.ncpu 2> /dev/null)
    else
        echo "Error: Unable to determine the number of CPU cores." >&2
        return 1
    fi

    #  Check that output is a valid integer
    if [[ ! "${cores}" =~ ^[0-9]+$ ]]; then
        echo "Error: Failed to retrieve a valid core count." >&2
        return 1
    fi

    echo "${cores}"
}


function set_params_parallel() {
    local threads="${1:-}"  # Total thread budget (not per-job threads)
    local max_job="${2:-}"  # Maximum number of jobs allowed
    local n_cores           # Total available CPU cores on system
    local par_job           # Resolved GNU Parallel job count
    local show_help

    show_help=$(cat << EOM
Usage:
  set_params_parallel [-h|--hlp|--help] threads max_job

Description:
  For the non-SLURM path, determine a safe combination of threads per job and GNU Parallel job count from a requested total thread budget and maximum job count.

Positional arguments:
  1  threads  <int>  Requested total local CPU/thread budget.
  2  max_job  <int>  Maximum number of concurrent jobs allowed.

Returns:
  '<threads_per_job>;<par_job>' to stdout; otherwise, 1 and error message to stderr.

Notes:
  - Intended for the non-SLURM / GNU Parallel execution path.
  - Interprets 'threads' as a total local CPU/thread budget, not as per-job threads.

#TODO:
  - Add usage example(s).
EOM
    )

    #  Handle help requests and missing required positional arguments
    if [[ "${threads}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif (( $# < 1 )); then
        echo "Error: Positional argument 1, 'threads', is missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif (( $# < 2 )); then
        echo "Error: Positional argument 2, 'max_job', is missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Check integer inputs
    if [[ ! "${threads}" =~ ^[0-9]+$ || "${threads}" -lt 1 ]]; then
        echo "Error: 'threads' must be a positive integer >= 1." >&2
        return 1
    fi

    if [[ ! "${max_job}" =~ ^[0-9]+$ || "${max_job}" -lt 1 ]]; then
        echo "Error: 'max_job' must be a positive integer >= 1." >&2
        return 1
    fi

    #  Get system CPU core count
    n_cores=$(determine_cores) || return 1

    #  Do not allow requested total threads to exceed available cores
    if [[ "${threads}" -gt "${n_cores}" ]]; then
        echo \
            "Warning: Requested threads ('${threads}') exceed available" \
            "cores ('${n_cores}'). Using '${n_cores}'." >&2
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
        echo \
            "Warning: Requested parallel jobs ('${par_job}') exceed" \
            "available cores ('${n_cores}'). Using '${n_cores}'." >&2
        par_job="${n_cores}"
    fi

    #  Convert total thread budget into per-job thread count, enforcing a
    #+ minimum of 1 thread per job
    if [[ "${threads}" -lt "${par_job}" ]]; then
        echo \
            "Warning: Requested total threads ('${threads}') are fewer than" \
            "parallel jobs ('${par_job}'). Using 1 thread per job." >&2
        threads=1
    else
        threads=$(( threads / par_job ))
    fi

    echo "${threads};${par_job}"  # Return adjusted threads-per-job, job count
}
