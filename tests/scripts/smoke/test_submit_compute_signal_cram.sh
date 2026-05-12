#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_compute_signal_cram.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="submit compute-signal CRAM"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

#  Define fixture and output paths for local serial CRAM-backed smoke tests
dir_fx="${ROOT_REPO}/tests/compute_signal/fixtures"
infile_se="${dir_fx}/cram/tiny_se.cram"
infile_pe="${dir_fx}/cram/tiny_pe.cram"
ref_fa="${dir_fx}/reference/tiny.fa"

tmp="${TEST_DIR_TMP}/submit_compute_signal_cram"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/compute_signal"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_exist "${infile_se}" "${infile_pe}" "${ref_fa}" || {
    finish
    exit $?
}


#  Run a local serial submit-wrapper CRAM case into compute_signal.py
function run_case_cram() {
    local nam_case="${1:-}"
    local mode="${2:-}"
    local infile_lcl="${3:-}"
    local outfile_lcl="${4:-}"
    local log_lcl="${5:-}"

    shift 5

    # shellcheck disable=SC2154
    if \
        run_capture \
            "submit compute-signal CRAM ${nam_case}" "${log_lcl}" \
            "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_compute_signal.sh" \
                --env_nam "${env_nam}" \
                --dir_scr "${ROOT_REPO}/scripts" \
                --threads 1 \
                --mode "${mode}" \
                --csv_infile "${infile_lcl}" \
                --csv_outfile "${outfile_lcl}" \
                --ref_fa "${ref_fa}" \
                --csv_usr_frg NA \
                --err_out "${dir_err}" \
                --nam_job "test_compute_cram_${nam_case}" \
                "$@"
    then
        rec_pass "submit_compute_signal.sh ${mode} ${nam_case} exits 0"
    else
        rec_fail \
            "submit_compute_signal.sh ${mode} ${nam_case} failed; see" \
            "$(rec_relpath "${log_lcl}")"
    fi
}


#  Signal mode: two 10-bp SE alignments produce chromosome-I bedGraph bins
outfile="${dir_out}/tiny_se_signal_unadj.bdg"
log="${dir_log}/submit_compute_signal_cram_se_signal.log"

run_case_cram \
    "se_signal" "signal" "${infile_se}" "${outfile}" "${log}" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA \
    --dp 3

assert_file_nonempty "${outfile}" "CRAM SE signal bedGraph output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t10$' \
        "CRAM SE signal output has chromosome-I bin I:0-10"
    assert_grep_pattern "${outfile}" $'^I\t20\t30\t10$' \
        "CRAM SE signal output has chromosome-I bin I:20-30"
fi


#  Coord mode: the same SE alignments emit BED-like processed fragments
outfile="${dir_out}/tiny_se_coord.bed"
log="${dir_log}/submit_compute_signal_cram_se_coord.log"

run_case_cram \
    "se_coord" "coord" "${infile_se}" "${outfile}" "${log}" \
    --dp 3

assert_file_nonempty "${outfile}" "CRAM SE coord BED output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t10$' \
        "CRAM SE coord output has chromosome-I fragment I:0-10"
    assert_grep_pattern "${outfile}" $'^I\t20\t30\t10$' \
        "CRAM SE coord output has chromosome-I fragment I:20-30"
fi


#  Signal mode: two PE fragments cover bins from I:10-60
outfile="${dir_out}/tiny_pe_signal_unadj.bdg"
log="${dir_log}/submit_compute_signal_cram_pe_signal.log"

run_case_cram \
    "pe_signal" "signal" "${infile_pe}" "${outfile}" "${log}" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA \
    --dp 3

assert_file_nonempty "${outfile}" "CRAM PE signal bedGraph output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t10\t20\t10$' \
        "CRAM PE signal output has chromosome-I bin I:10-20"
    assert_grep_pattern "${outfile}" $'^I\t40\t50\t10$' \
        "CRAM PE signal output has chromosome-I bin I:40-50"
fi


#  Coord mode: PE output emits one BED-like row per leftmost proper pair
outfile="${dir_out}/tiny_pe_coord.bed"
log="${dir_log}/submit_compute_signal_cram_pe_coord.log"

run_case_cram \
    "pe_coord" "coord" "${infile_pe}" "${outfile}" "${log}" \
    --dp 3

assert_file_nonempty "${outfile}" "CRAM PE coord BED output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t10\t40\t30$' \
        "CRAM PE coord output has chromosome-I fragment I:10-40"
    assert_grep_pattern "${outfile}" $'^I\t40\t60\t20$' \
        "CRAM PE coord output has chromosome-I fragment I:40-60"
fi


#  CRAM input without a reference FASTA should fail clearly in the wrapper
outfile="${dir_out}/tiny_se_missing_ref.bdg"
log="${dir_log}/submit_compute_signal_cram_missing_ref.log"

if \
    run_capture \
        "submit compute-signal CRAM missing_ref" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_compute_signal.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --mode signal \
            --csv_infile "${infile_se}" \
            --csv_outfile "${outfile}" \
            --csv_usr_frg NA \
            --err_out "${dir_err}" \
            --nam_job "test_compute_cram_missing_ref" \
            --method unadj \
            --siz_bin 10 \
            --csv_scl_fct NA \
            --dp 3
then
    rec_fail "submit_compute_signal.sh CRAM without --ref_fa unexpectedly passed"
else
    rec_pass "submit_compute_signal.sh CRAM without --ref_fa fails"
fi

assert_grep_pattern "${log}" \
    "'--ref_fa' is required when '--csv_infile' contains CRAM" \
    "submit CRAM missing-ref error mentions --ref_fa"

if [[ ! -s "${outfile}" ]]; then
    rec_pass "submit CRAM missing-ref output is absent or empty"
else
    rec_fail "submit CRAM missing-ref output was unexpectedly written"
fi

finish
