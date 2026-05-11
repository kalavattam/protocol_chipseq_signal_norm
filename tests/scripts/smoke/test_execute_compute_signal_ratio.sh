#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_compute_signal_ratio.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute compute-signal ratio"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

#  Define fixture and output paths for the execute-to-submit ratio path
dir_fx="${ROOT_REPO}/tests/compute_signal/fixtures/bedgraph"
fil_A="${dir_fx}/ratio_A.bdg"
fil_B="${dir_fx}/ratio_B.bdg"

tmp="${TEST_DIR_TMP}/execute_compute_signal_ratio"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/compute_signal"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_exist "${fil_A}" "${fil_B}" || {
    finish
    exit $?
}


#  Run a local serial execute-wrapper ratio case through submit and Python
function run_case_ratio() {
    local nam_case="${1:-}"
    local log_lcl="${2:-}"
    local pfx_lcl="${3:-exec}"

    shift 3

    # shellcheck disable=SC2154
    if \
        run_capture \
            "execute compute-signal ratio ${nam_case}" "${log_lcl}" \
            "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_compute_signal.sh" \
                --threads 1 \
                --mode ratio \
                --method unadj \
                --csv_fil_A "${fil_A}" \
                --csv_fil_B "${fil_B}" \
                --dir_out "${dir_out}" \
                --typ_out bdg \
                --prefix "${pfx_lcl}" \
                --eps 0 \
                --dp 3 \
                --err_out "${dir_err}" \
                --nam_job "test_execute_compute_ratio_${nam_case}" \
                --max_job 1 \
                "$@"
    then
        rec_pass "execute_compute_signal.sh ratio ${nam_case} exits 0"
    else
        rec_fail \
            "execute_compute_signal.sh ratio ${nam_case} failed; see" \
            "$(rec_relpath "${log_lcl}")"
    fi
}


#  Baseline unadjusted ratio with three-decimal rounding
outfile="${dir_out}/exec_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_unadj.log"

run_case_ratio "unadj" "${log}" "exec"

assert_file_nonempty "${outfile}" "execute ratio output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t2$' \
        "execute ratio output has I:0-10 = 2"
    assert_grep_pattern "${outfile}" $'^I\t10\t20\t0$' \
        "execute ratio output has I:10-20 = 0"
    assert_grep_pattern "${outfile}" $'^I\t40\t50\t4$' \
        "execute ratio output has I:40-50 = 4"
    assert_grep_pattern "${outfile}" $'^I\t60\t70\t0.333$' \
        "execute ratio output has I:60-70 = 0.333"
    assert_grep_pattern "${outfile}" $'^I\t70\t80\t1$' \
        "execute ratio output has I:70-80 = 1"
fi


#  Denominator floor: B=0.04 is floored to 0.1, so 1 / 0.1 = 10
outfile="${dir_out}/exec_dep_min_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_dep_min.log"

run_case_ratio "dep_min" "${log}" "exec_dep_min" --csv_dep_min 0.1

assert_file_nonempty "${outfile}" "execute dep_min ratio output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t50\t60\t10$' \
        "execute dep_min ratio output has I:50-60 = 10"
fi


#  Pseudocounts: (0 + 1) / (2 + 1) = 0.333 at three decimals
outfile="${dir_out}/exec_pseudo_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_pseudo.log"

run_case_ratio "pseudo" "${log}" "exec_pseudo" --csv_pseudo 1:1

assert_file_nonempty "${outfile}" "execute pseudo ratio output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t10\t20\t0.333$' \
        "execute pseudo ratio output has I:10-20 = 0.333"
fi


#  Drop non-finite rows while preserving finite ratio rows
outfile="${dir_out}/exec_drp_nan_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_drp_nan.log"

run_case_ratio "drp_nan" "${log}" "exec_drp_nan" --drp_nan

assert_file_nonempty "${outfile}" "execute drp_nan ratio output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t2$' \
        "execute drp_nan ratio output retains I:0-10 = 2"
    assert_grep_pattern "${outfile}" $'^I\t60\t70\t0.333$' \
        "execute drp_nan ratio output retains I:60-70 = 0.333"
    assert_no_grep_pattern "${outfile}" $'^I\t20\t30\t' \
        "execute drp_nan ratio output omits I:20-30"
    assert_no_grep_pattern "${outfile}" $'^I\t30\t40\t' \
        "execute drp_nan ratio output omits I:30-40"
fi


#  Zero-zero skipping before scaling removes the A=0, B=0 bin
outfile="${dir_out}/exec_skip_00_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_skip_00.log"

run_case_ratio "skip_00" "${log}" "exec_skip_00" --skip_00 pre_scale

assert_file_nonempty "${outfile}" "execute skip_00 ratio output"

if [[ -s "${outfile}" ]]; then
    assert_no_grep_pattern "${outfile}" $'^I\t30\t40\t' \
        "execute skip_00 ratio output omits I:30-40"
fi


#  Track sidecar should be generated and should omit non-finite rows
outfile="${dir_out}/exec_track_ratio_A.bdg"
trackfile="${dir_out}/exec_track_ratio_A.track.bdg"
log="${dir_log}/execute_compute_signal_ratio_track.log"

run_case_ratio "track" "${log}" "exec_track" --track

assert_file_nonempty "${outfile}" "execute track main ratio output"
assert_file_nonempty "${trackfile}" "execute track sidecar output"

if [[ -s "${trackfile}" ]]; then
    assert_grep_pattern "${trackfile}" $'^I\t0\t10\t2$' \
        "execute track sidecar retains I:0-10 = 2"
    assert_grep_pattern "${trackfile}" $'^I\t60\t70\t0.333$' \
        "execute track sidecar retains I:60-70 = 0.333"
    assert_no_grep_pattern "${trackfile}" $'^I\t20\t30\t' \
        "execute track sidecar omits I:20-30"
    assert_no_grep_pattern "${trackfile}" $'^I\t30\t40\t' \
        "execute track sidecar omits I:30-40"
fi


finish
