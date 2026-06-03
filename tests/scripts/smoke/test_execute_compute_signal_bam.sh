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

print_section "${TEST_NAME}"

#  Define fixture and output paths for the execute-to-submit BAM path
dir_fx="${ROOT_REPO}/tests/compute_signal/fixtures"
in_se="${dir_fx}/bam/se/tiny_se.bam"
in_pe="${dir_fx}/bam/pe/tiny_pe.bam"

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

require_files_nonempty \
    "${in_se}" \
    "${in_pe}" || {
    finish
    exit $?
}


#  Signal mode: two 10-bp SE alignments produce chromosome-I bedGraph bins
fil_out="${dir_out}/tiny_se.bdg"
log="${dir_log}/execute_compute_signal_bam_se_signal.log"

run_case_compute_signal \
    execute \
    bam \
    "se_signal" \
    "signal" \
    "${in_se}" \
    "bdg" \
    "${log}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA

assert_file_nonempty \
    "${fil_out}" \
    "execute SE signal bedGraph output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t0\t10\t10$' \
        "execute SE signal output has chromosome-I bin I:0-10"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t20\t30\t10$' \
        "execute SE signal output has chromosome-I bin I:20-30"
fi


#  Coord mode: the same SE alignments emit BED-like processed fragments
fil_out="${dir_out}/tiny_se.bed"
log="${dir_log}/execute_compute_signal_bam_se_coord.log"

run_case_compute_signal \
    execute \
    bam \
    "se_coord" \
    "coord" \
    "${in_se}" \
    "bed" \
    "${log}" \
    "${dir_out}" \
    "${dir_err}" \
    ""

assert_file_nonempty \
    "${fil_out}" \
    "execute SE coord BED output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t0\t10\t10$' \
        "execute SE coord output has chromosome-I fragment I:0-10"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t20\t30\t10$' \
        "execute SE coord output has chromosome-I fragment I:20-30"
fi


#  Signal mode: two PE fragments cover bins from I:10-60
fil_out="${dir_out}/tiny_pe.bdg"
log="${dir_log}/execute_compute_signal_bam_pe_signal.log"

run_case_compute_signal \
    execute \
    bam \
    "pe_signal" \
    "signal" \
    "${in_pe}" \
    "bdg" \
    "${log}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA

assert_file_nonempty \
    "${fil_out}" \
    "execute PE signal bedGraph output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t10\t20\t10$' \
        "execute PE signal output has chromosome-I bin I:10-20"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t40\t50\t10$' \
        "execute PE signal output has chromosome-I bin I:40-50"
fi


#  Coord mode: PE output emits one BED-like row per leftmost proper pair
fil_out="${dir_out}/tiny_pe.bed"
log="${dir_log}/execute_compute_signal_bam_pe_coord.log"

run_case_compute_signal \
    execute \
    bam \
    "pe_coord" \
    "coord" \
    "${in_pe}" \
    "bed" \
    "${log}" \
    "${dir_out}" \
    "${dir_err}" \
    ""

assert_file_nonempty \
    "${fil_out}" \
    "execute PE coord BED output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t10\t40\t30$' \
        "execute PE coord output has chromosome-I fragment I:10-40"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t40\t60\t20$' \
        "execute PE coord output has chromosome-I fragment I:40-60"
fi


#  Scaling factor and prefix propagation: raw 10-bp SE bins scaled by 2
fil_out="${dir_out}/scaled.tiny_se.bdg"
log="${dir_log}/execute_compute_signal_bam_se_signal_scaled.log"

run_case_compute_signal \
    execute \
    bam \
    "se_signal_scaled" \
    "signal" \
    "${in_se}" \
    "bdg" \
    "${log}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct 2 \
    --prefix scaled

assert_file_nonempty \
    "${fil_out}" \
    "execute scaled SE signal output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t0\t10\t20$' \
        "execute scaled SE signal output has I:0-10 = 20"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t20\t30\t20$' \
        "execute scaled SE signal output has I:20-30 = 20"
fi


#  Fixed SE fragment length extends reads to 20 bp in signal mode
fil_out="${dir_out}/usr_frg_signal.tiny_se.bdg"
log="${dir_log}/execute_compute_signal_bam_se_signal_usr_frg.log"

run_case_compute_signal \
    execute \
    bam \
    "se_signal_usr_frg" \
    "signal" \
    "${in_se}" \
    "bdg" \
    "${log}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA \
    --csv_usr_frg 20 \
    --prefix usr_frg_signal

assert_file_nonempty \
    "${fil_out}" \
    "execute usr_frg SE signal output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t0\t10\t10$' \
        "execute usr_frg SE signal output has I:0-10 = 10"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t10\t20\t20$' \
        "execute usr_frg SE signal output has I:10-20 = 20"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t20\t30\t10$' \
        "execute usr_frg SE signal output has I:20-30 = 10"
fi


#  Fixed SE fragment length is reflected in coord-mode BED intervals
fil_out="${dir_out}/usr_frg_coord.tiny_se.bed"
log="${dir_log}/execute_compute_signal_bam_se_coord_usr_frg.log"

run_case_compute_signal \
    execute \
    bam \
    "se_coord_usr_frg" \
    "coord" \
    "${in_se}" \
    "bed" \
    "${log}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --csv_usr_frg 20 \
    --prefix usr_frg_coord

assert_file_nonempty \
    "${fil_out}" \
    "execute usr_frg SE coord output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t0\t20\t20$' \
        "execute usr_frg SE coord output has chromosome-I fragment I:0-20"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t10\t30\t20$' \
        "execute usr_frg SE coord output has chromosome-I fragment I:10-30"
fi

finish
