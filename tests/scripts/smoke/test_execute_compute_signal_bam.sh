#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_compute_signal_bam.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute compute-signal BAM"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

#  Define fixture and output paths for the execute-to-submit BAM path
dir_fx="${ROOT_REPO}/tests/compute_signal/fixtures/bam"
infile_se="${dir_fx}/tiny_se.bam"
infile_pe="${dir_fx}/tiny_pe.bam"

tmp="${TEST_DIR_TMP}/execute_compute_signal_bam"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/compute_signal"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_exist "${infile_se}" "${infile_pe}" || {
    finish
    exit $?
}


#  Run the execute wrapper through submit_compute_signal.sh into Python
function run_case_bam() {
    local nam_case="${1:-}"
    local mode="${2:-}"
    local infile_lcl="${3:-}"
    local typ_out="${4:-}"
    local log_lcl="${5:-}"

    shift 5

    # shellcheck disable=SC2154
    if \
        run_capture \
            "execute compute-signal BAM ${nam_case}" "${log_lcl}" \
            "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_compute_signal.sh" \
                --threads 1 \
                --mode "${mode}" \
                --csv_infile "${infile_lcl}" \
                --dir_out "${dir_out}" \
                --typ_out "${typ_out}" \
                --csv_usr_frg NA \
                --dp 3 \
                --err_out "${dir_err}" \
                --nam_job "test_execute_compute_bam_${nam_case}" \
                --max_job 1 \
                "$@"
    then
        rec_pass "execute_compute_signal.sh ${mode} ${nam_case} exits 0"
    else
        rec_fail \
            "execute_compute_signal.sh ${mode} ${nam_case} failed; see" \
            "$(rec_relpath "${log_lcl}")"
    fi
}


#  Signal mode: two 10-bp SE alignments produce chromosome-I bedGraph bins
outfile="${dir_out}/tiny_se.bdg"
log="${dir_log}/execute_compute_signal_bam_se_signal.log"

run_case_bam \
    "se_signal" "signal" "${infile_se}" "bdg" "${log}" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA

assert_file_nonempty "${outfile}" "execute SE signal bedGraph output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t10$' \
        "execute SE signal output has chromosome-I bin I:0-10"
    assert_grep_pattern "${outfile}" $'^I\t20\t30\t10$' \
        "execute SE signal output has chromosome-I bin I:20-30"
fi


#  Coord mode: the same SE alignments emit BED-like processed fragments
outfile="${dir_out}/tiny_se.bed"
log="${dir_log}/execute_compute_signal_bam_se_coord.log"

run_case_bam \
    "se_coord" "coord" "${infile_se}" "bed" "${log}"

assert_file_nonempty "${outfile}" "execute SE coord BED output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t10$' \
        "execute SE coord output has chromosome-I fragment I:0-10"
    assert_grep_pattern "${outfile}" $'^I\t20\t30\t10$' \
        "execute SE coord output has chromosome-I fragment I:20-30"
fi


#  Signal mode: two PE fragments cover bins from I:10-60
outfile="${dir_out}/tiny_pe.bdg"
log="${dir_log}/execute_compute_signal_bam_pe_signal.log"

run_case_bam \
    "pe_signal" "signal" "${infile_pe}" "bdg" "${log}" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA

assert_file_nonempty "${outfile}" "execute PE signal bedGraph output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t10\t20\t10$' \
        "execute PE signal output has chromosome-I bin I:10-20"
    assert_grep_pattern "${outfile}" $'^I\t40\t50\t10$' \
        "execute PE signal output has chromosome-I bin I:40-50"
fi


#  Coord mode: PE output emits one BED-like row per leftmost proper pair
outfile="${dir_out}/tiny_pe.bed"
log="${dir_log}/execute_compute_signal_bam_pe_coord.log"

run_case_bam \
    "pe_coord" "coord" "${infile_pe}" "bed" "${log}"

assert_file_nonempty "${outfile}" "execute PE coord BED output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t10\t40\t30$' \
        "execute PE coord output has chromosome-I fragment I:10-40"
    assert_grep_pattern "${outfile}" $'^I\t40\t60\t20$' \
        "execute PE coord output has chromosome-I fragment I:40-60"
fi


#  Scaling factor and prefix propagation: raw 10-bp SE bins scaled by 2
outfile="${dir_out}/scaled.tiny_se.bdg"
log="${dir_log}/execute_compute_signal_bam_se_signal_scaled.log"

run_case_bam \
    "se_signal_scaled" "signal" "${infile_se}" "bdg" "${log}" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct 2 \
    --prefix scaled

assert_file_nonempty "${outfile}" "execute scaled SE signal output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t20$' \
        "execute scaled SE signal output has I:0-10 = 20"
    assert_grep_pattern "${outfile}" $'^I\t20\t30\t20$' \
        "execute scaled SE signal output has I:20-30 = 20"
fi


#  Fixed SE fragment length extends reads to 20 bp in signal mode
outfile="${dir_out}/usr_frg_signal.tiny_se.bdg"
log="${dir_log}/execute_compute_signal_bam_se_signal_usr_frg.log"

run_case_bam \
    "se_signal_usr_frg" "signal" "${infile_se}" "bdg" "${log}" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA \
    --csv_usr_frg 20 \
    --prefix usr_frg_signal

assert_file_nonempty "${outfile}" "execute usr_frg SE signal output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t10$' \
        "execute usr_frg SE signal output has I:0-10 = 10"
    assert_grep_pattern "${outfile}" $'^I\t10\t20\t20$' \
        "execute usr_frg SE signal output has I:10-20 = 20"
    assert_grep_pattern "${outfile}" $'^I\t20\t30\t10$' \
        "execute usr_frg SE signal output has I:20-30 = 10"
fi


#  Fixed SE fragment length is reflected in coord-mode BED intervals
outfile="${dir_out}/usr_frg_coord.tiny_se.bed"
log="${dir_log}/execute_compute_signal_bam_se_coord_usr_frg.log"

run_case_bam \
    "se_coord_usr_frg" "coord" "${infile_se}" "bed" "${log}" \
    --csv_usr_frg 20 \
    --prefix usr_frg_coord

assert_file_nonempty "${outfile}" "execute usr_frg SE coord output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t20\t20$' \
        "execute usr_frg SE coord output has chromosome-I fragment I:0-20"
    assert_grep_pattern "${outfile}" $'^I\t10\t30\t20$' \
        "execute usr_frg SE coord output has chromosome-I fragment I:10-30"
fi

finish
