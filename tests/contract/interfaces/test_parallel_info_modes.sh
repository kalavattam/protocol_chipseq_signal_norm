#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_parallel_info_modes.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="print_parallel_info validates only the active mode's job count"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Run 'print_parallel_info' in a subshell, returning its status and output. The
# function is sourced per invocation rather than once at file scope, so that a
# failing call cannot leave shared state behind for the next case.
function run_parallel_info() {
    # shellcheck disable=SC2016  # Expand in the child shell, not this one.
    "${TEST_BASH}" -c '
        source "${1}/lib/bash/core/source_helpers.sh"
        source_helpers "${1}/lib/bash" \
            core/check_args \
            core/check_inputs \
            core/format_outputs \
            dispatch/manage_parallel
        shift 1
        print_parallel_info "$@"
    ' _ "${ROOT_REPO}" "$@" 2>&1
}


# Assert that a mode reports successfully and names its own job count.
function check_mode_accepts() {
    local label="${1}"   # Case name recorded by the harness.
    local expect="${2}"  # Substring the successful report must contain.
    shift 2

    local out=""      # Captured combined output.
    local status=0    # Captured exit status.

    set +e
    out="$(run_parallel_info "$@")"
    status=$?
    set -e

    if [[ "${status}" -ne 0 ]]; then
        record_fail "${label}: exited ${status}; expected success"
        return 0
    fi

    if [[ "${out}" == *"${expect}"* ]]; then
        record_pass "${label}"
    else
        record_fail "${label}: output did not contain '${expect}'"
    fi
}


# Assert that a genuinely invalid value is still rejected.
function check_mode_rejects() {
    local label="${1}"   # Case name recorded by the harness.
    local expect="${2}"  # Substring the rejection message must contain.
    shift 2

    local out=""      # Captured combined output.
    local status=0    # Captured exit status.

    set +e
    out="$(run_parallel_info "$@")"
    status=$?
    set -e

    if [[ "${status}" -eq 0 ]]; then
        record_fail "${label}: succeeded; expected rejection"
        return 0
    fi

    if [[ "${out}" == *"${expect}"* ]]; then
        record_pass "${label}"
    else
        record_fail "${label}: error did not mention '${expect}'"
    fi
}


print_section "${TEST_NAME}"

# Each dispatch mode leaves the other mode's job count unresolved, and the
# callers in 'bin/execute_*.sh' pass a sentinel for it. Accepting that sentinel
# is the contract: 'print_parallel_info' reads 'max_job' only under Slurm and
# 'par_job' only outside it, so validating both unconditionally turned
# '--verbose' into a hard failure in every mode but serial.
check_mode_accepts \
    "Slurm mode accepts an unresolved 'par_job'" \
    "Max concurrent jobs (Slurm): 4" \
    true 4 UNSET 12

check_mode_accepts \
    "GNU Parallel mode accepts an unresolved 'max_job'" \
    "Max concurrent jobs (GNU Parallel): 4" \
    false UNSET 4 12

check_mode_accepts \
    "serial mode reports without a sentinel" \
    "Jobs running in serial mode: 1" \
    false 1 1 12

# The relaxation must not become an absence of validation, so the argument the
# active mode does read is still checked.
check_mode_rejects \
    "Slurm mode rejects a non-numeric 'max_job'" \
    "'max_job'" \
    true bogus UNSET 12

check_mode_rejects \
    "local mode rejects a non-numeric 'par_job'" \
    "'par_job'" \
    false UNSET bogus 12

check_mode_rejects \
    "local mode rejects a zero 'par_job'" \
    "'par_job'" \
    false UNSET 0 12

check_mode_rejects \
    "every mode rejects a non-numeric 'threads'" \
    "'threads'" \
    true 4 UNSET bogus

# shellcheck disable=SC2119  # 'finish' receives no script arguments.
finish
