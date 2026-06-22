#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_relativize_scaling_factors.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="relativize scaling factors"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"


function run_relativize_success() {
    local label="${1:-success}"
    local log="${2:-}"
    local tbl="${3:-}"
    shift 3

    if \
        PYTHONPATH="${ROOT_REPO}" \
        run_capture \
            "relativize scaling factors ${label}" \
            "${log}" \
            "${py_cmd[@]}" "${scr_rel}" --infile "${tbl}" "$@"
    then
        record_pass "relativize scaling factors ${label} exits 0"
    else
        record_fail \
            "relativize scaling factors ${label} failed; see" \
            "$(print_relpath "${log}")"
    fi
}


function run_relativize_failure() {
    local label="${1:-failure}"
    local log="${2:-}"
    local tbl="${3:-}"
    local patn="${4:-Error}"
    shift 4

    if ! \
        PYTHONPATH="${ROOT_REPO}" \
        run_capture \
            "relativize scaling factors ${label}" \
            "${log}" \
            "${py_cmd[@]}" "${scr_rel}" --infile "${tbl}" "$@"
    then
        assert_pattern_found \
            "${log}" \
            "${patn}" \
            "relativize scaling factors ${label}"
    else
        record_fail "relativize scaling factors ${label} unexpectedly passed"
    fi
}


print_section "${TEST_NAME}"

scr_rel="${ROOT_REPO}/scripts/relativize_scaling_factors.py"

dir_tmp="${TEST_DIR_TMP}/relativize_scaling_factors"
dir_log="${TEST_DIR_LOG}/relativize_scaling_factors"

tbl_spike="${dir_tmp}/spike.tsv"
tbl_siq="${dir_tmp}/siq.tsv"
tbl_no_scl="${dir_tmp}/missing_scaling.tsv"
tbl_no_sample="${dir_tmp}/missing_sample.tsv"

exp_spike="${dir_tmp}/spike.expected.tsv"
exp_spike_input="${dir_tmp}/spike_input.expected.tsv"
exp_siq="${dir_tmp}/siq.expected.tsv"

out_spike_no_input="${dir_log}/spike_no_input.tsv"
out_spike_with_input="${dir_log}/spike_with_input.tsv"
out_siq_preferred="${dir_log}/siq_preferred.tsv"

log_missing_scaling_column="${dir_log}/missing_scaling_column.err"
log_missing_sample_column="${dir_log}/missing_sample_column.err"

rm -rf "${dir_tmp}"
mkdir -p "${dir_tmp}" "${dir_log}"

require_files_nonempty "${scr_rel}" || {
    finish
    exit $?
}

if [[ -n "${PYTHON:-}" ]]; then
    py_cmd=( "${PYTHON}" )
elif py="$(find_python)"; then
    py_cmd=( "${py}" )
else
    record_fail "python is available for relativize_scaling_factors.py"
    finish
    exit $?
fi

{
    printf 'sample\tspike\tother\n'
    printf 'IP_a\t2\ta\n'
    printf 'IP_b\t4\tb\n'
    printf 'in_a\t10\tc\n'
} > "${tbl_spike}"

{
    printf 'sample\tspike\tscaled\tother\n'
    printf 'IP_a\t2\t0.5\ta\n'
    printf 'IP_b\t4\t1\tb\n'
    printf 'in_a\t10\t1\tc\n'
} > "${exp_spike}"

{
    printf 'sample\tspike\tscaled\tother\n'
    printf 'IP_a\t2\t0.5\ta\n'
    printf 'IP_b\t4\t1\tb\n'
    printf 'in_a\t10\t2.5\tc\n'
} > "${exp_spike_input}"

{
    printf 'sample\tsiq\tspike\n'
    printf 'IP_a\t3\t100\n'
    printf 'IP_b\t6\t200\n'
} > "${tbl_siq}"

{
    printf 'sample\tsiq\tscaled\tspike\n'
    printf 'IP_a\t3\t0.5\t100\n'
    printf 'IP_b\t6\t1\t200\n'
} > "${exp_siq}"

{
    printf 'sample\tcoef\n'
    printf 'IP_a\t1\n'
} > "${tbl_no_scl}"

{
    printf 'spike\n'
    printf '1\n'
} > "${tbl_no_sample}"


run_relativize_success \
    spike_no_input \
    "${out_spike_no_input}" \
    "${tbl_spike}"

assert_files_equal \
    "${out_spike_no_input}" \
    "${exp_spike}" \
    "relativize scaling factors spike_no_input has expected TSV"

run_relativize_success \
    spike_with_input \
    "${out_spike_with_input}" \
    "${tbl_spike}" \
    --input

assert_files_equal \
    "${out_spike_with_input}" \
    "${exp_spike_input}" \
    "relativize scaling factors spike_with_input has expected TSV"

run_relativize_success \
    siq_preferred \
    "${out_siq_preferred}" \
    "${tbl_siq}"

assert_files_equal \
    "${out_siq_preferred}" \
    "${exp_siq}" \
    "relativize scaling factors siq_preferred has expected TSV"

run_relativize_failure \
    missing_scaling_column \
    "${log_missing_scaling_column}" \
    "${tbl_no_scl}" \
    "No supported scaling-factor column"

run_relativize_failure \
    missing_sample_column \
    "${log_missing_sample_column}" \
    "${tbl_no_sample}" \
    "missing required column 'sample'"

finish
