#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_compute_signal_cram.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="execute compute-signal CRAM"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

#  Define fixture and output paths for the execute-to-submit CRAM path
dir_fx="${ROOT_REPO}/tests/compute_signal/fixtures"
in_se="${dir_fx}/cram/se/tiny_se.cram"
in_pe="${dir_fx}/cram/pe/tiny_pe.cram"
ref_fa="${dir_fx}/reference/tiny.fa"

tmp="${TEST_DIR_TMP}/execute_compute_signal_cram"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/compute_signal"

fil_out_se_signal="${dir_out}/tiny_se.bdg"
fil_out_se_coord="${dir_out}/tiny_se.bed"
fil_out_pe_signal="${dir_out}/tiny_pe.bdg"
fil_out_pe_coord="${dir_out}/tiny_pe.bed"
fil_out_missing_ref="${dir_out}/tiny_se_missing_ref.bdg"

log_se_signal="${dir_log}/execute_compute_signal_cram_se_signal.log"
log_se_coord="${dir_log}/execute_compute_signal_cram_se_coord.log"
log_pe_signal="${dir_log}/execute_compute_signal_cram_pe_signal.log"
log_pe_coord="${dir_log}/execute_compute_signal_cram_pe_coord.log"
log_missing_ref="${dir_log}/execute_compute_signal_cram_missing_ref.log"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${in_se}" \
    "${in_pe}" \
    "${ref_fa}" || {
    finish
    exit $?
}


#  Signal mode: two 10-bp SE alignments produce chromosome-I bedGraph bins
run_case_compute_signal \
    execute \
    cram \
    "se_signal" \
    "signal" \
    "${in_se}" \
    "bdg" \
    "${log_se_signal}" \
    "${dir_out}" \
    "${dir_err}" \
    "${ref_fa}" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA

assert_file_nonempty \
    "${fil_out_se_signal}" \
    "execute CRAM SE signal bedGraph output"

if [[ -s "${fil_out_se_signal}" ]]; then
    assert_pattern_found \
        "${fil_out_se_signal}" \
        $'^I\t0\t10\t10$' \
        "execute CRAM SE signal output has chromosome-I bin I:0-10"

    assert_pattern_found \
        "${fil_out_se_signal}" \
        $'^I\t20\t30\t10$' \
        "execute CRAM SE signal output has chromosome-I bin I:20-30"
fi


#  Coord mode: the same SE alignments emit BED-like processed fragments
run_case_compute_signal \
    execute \
    cram \
    "se_coord" \
    "coord" \
    "${in_se}" \
    "bed" \
    "${log_se_coord}" \
    "${dir_out}" \
    "${dir_err}" \
    "${ref_fa}"

assert_file_nonempty \
    "${fil_out_se_coord}" \
    "execute CRAM SE coord BED output"

if [[ -s "${fil_out_se_coord}" ]]; then
    assert_pattern_found \
        "${fil_out_se_coord}" \
        $'^I\t0\t10\t10$' \
        "execute CRAM SE coord output has chromosome-I fragment I:0-10"

    assert_pattern_found \
        "${fil_out_se_coord}" \
        $'^I\t20\t30\t10$' \
        "execute CRAM SE coord output has chromosome-I fragment I:20-30"
fi


#  Signal mode: two PE fragments cover bins from I:10-60
run_case_compute_signal \
    execute \
    cram \
    "pe_signal" \
    "signal" \
    "${in_pe}" \
    "bdg" \
    "${log_pe_signal}" \
    "${dir_out}" \
    "${dir_err}" \
    "${ref_fa}" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA

assert_file_nonempty \
    "${fil_out_pe_signal}" \
    "execute CRAM PE signal bedGraph output"

if [[ -s "${fil_out_pe_signal}" ]]; then
    assert_pattern_found \
        "${fil_out_pe_signal}" \
        $'^I\t10\t20\t10$' \
        "execute CRAM PE signal output has chromosome-I bin I:10-20"

    assert_pattern_found \
        "${fil_out_pe_signal}" \
        $'^I\t40\t50\t10$' \
        "execute CRAM PE signal output has chromosome-I bin I:40-50"
fi


#  Coord mode: PE output emits one BED-like row per leftmost proper pair
run_case_compute_signal \
    execute \
    cram \
    "pe_coord" \
    "coord" \
    "${in_pe}" \
    "bed" \
    "${log_pe_coord}" \
    "${dir_out}" \
    "${dir_err}" \
    "${ref_fa}"

assert_file_nonempty \
    "${fil_out_pe_coord}" \
    "execute CRAM PE coord BED output"

if [[ -s "${fil_out_pe_coord}" ]]; then
    assert_pattern_found \
        "${fil_out_pe_coord}" \
        $'^I\t10\t40\t30$' \
        "execute CRAM PE coord output has chromosome-I fragment I:10-40"

    assert_pattern_found \
        "${fil_out_pe_coord}" \
        $'^I\t40\t60\t20$' \
        "execute CRAM PE coord output has chromosome-I fragment I:40-60"
fi


#  CRAM input without a reference FASTA should fail clearly in the wrapper
if \
    run_capture \
        "execute compute-signal CRAM missing_ref" \
        "${log_missing_ref}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_compute_signal.sh" \
            --threads 1 \
            --mode signal \
            --csv_infile "${in_se}" \
            --dir_out "${dir_out}" \
            --typ_out bdg \
            --csv_usr_frg NA \
            --dp 3 \
            --err_out "${dir_err}" \
            --nam_job "test_execute_compute_cram_missing_ref" \
            --max_job 1 \
            --method unadj \
            --siz_bin 10 \
            --csv_scl_fct NA
then
    record_fail \
        "execute_compute_signal.sh CRAM without --ref_fa unexpectedly passed"
else
    record_pass "execute_compute_signal.sh CRAM without --ref_fa fails"
fi

assert_pattern_found \
    "${log_missing_ref}" \
    "'--ref_fa' is required when '--csv_infile' contains CRAM" \
    "execute CRAM missing-ref error mentions --ref_fa"

if [[ ! -s "${fil_out_missing_ref}" ]]; then
    record_pass "execute CRAM missing-ref output is absent or empty"
else
    record_fail "execute CRAM missing-ref output was unexpectedly written"
fi

finish
