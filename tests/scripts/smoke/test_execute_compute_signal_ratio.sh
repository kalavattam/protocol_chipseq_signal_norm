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


#  Run the execute wrapper through submit_compute_signal.sh into Python
outfile="${dir_out}/exec_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_unadj.log"

# shellcheck disable=SC2154
if \
    run_capture \
        "execute compute-signal ratio unadj" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_compute_signal.sh" \
            --threads 1 \
            --mode ratio \
            --method unadj \
            --csv_fil_A "${fil_A}" \
            --csv_fil_B "${fil_B}" \
            --dir_out "${dir_out}" \
            --typ_out bdg \
            --prefix exec \
            --eps 0 \
            --dp 3 \
            --err_out "${dir_err}" \
            --nam_job "test_execute_compute_ratio" \
            --max_job 1
then
    rec_pass "execute_compute_signal.sh ratio exits 0"
else
    rec_fail \
        "execute_compute_signal.sh ratio failed; see" \
        "$(rec_relpath "${log}")"
fi

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

finish
