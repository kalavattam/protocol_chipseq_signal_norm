#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_write_header.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="write scaling-factor header"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"


function run_header_success() {
    local label="${1:-success}"
    local log="${2:-}"
    shift 2

    if \
        run_capture \
            "write scaling-factor header ${label}" \
            "${log}" \
            "${TEST_BASH}" "${scr_hdr}" "$@"
    then
        record_pass "write_header.sh ${label} exits 0"
    else
        record_fail \
            "write_header.sh ${label} failed; see" \
            "$(print_relpath "${log}")"
    fi
}


function run_header_failure() {
    local label="${1:-failure}"
    local log="${2:-}"
    local patn="${3:-Error}"
    shift 3

    if \
        run_capture \
            "write scaling-factor header ${label}" \
            "${log}" \
            "${TEST_BASH}" "${scr_hdr}" "$@"
    then
        record_fail "write_header.sh ${label} unexpectedly succeeded"
    else
        assert_pattern_found \
            "${log}" \
            "${patn}" \
            "write_header.sh ${label}"
    fi
}


function assert_header_once() {
    local file="${1:-}"
    local header="${2:-}"
    local label="${3:-}"
    local num_hdr=""

    assert_file_nonempty \
        "${file}" \
        "write_header.sh ${label} output"

    assert_pattern_found \
        "${file}" \
        "^${header}$" \
        "write_header.sh ${label} writes expected header"

    num_hdr="$(grep -c "^${header}$" "${file}" || true)"
    if [[ "${num_hdr}" -eq 1 ]]; then
        record_pass "write_header.sh ${label} writes one header"
    else
        record_fail "write_header.sh ${label} header count is not one"
    fi

    if [[ "$(sed -n '1p' "${file}")" == "${header}" ]]; then
        record_pass "write_header.sh ${label} puts header first"
    else
        record_fail "write_header.sh ${label} did not put header first"
    fi
}


print_section "${TEST_NAME}"

scr_hdr="${ROOT_REPO}/scripts/write_header.sh"

dir_tmp="${TEST_DIR_TMP}/write_header"
dir_log="${TEST_DIR_LOG}/write_header"
dir_in="${dir_tmp}/in"
dir_out="${dir_tmp}/out"

fil_spk_new="${dir_out}/spike.header_only.tsv"
fil_spk_vrb="${dir_out}/spike.verbose.tsv"
fil_spk_hdr="${dir_out}/spike.headered.tsv"

fil_siq_def="${dir_out}/siq.default.tsv"
fil_siq_dat="${dir_in}/siq.data.tsv"
fil_siq_cpy="${dir_out}/siq.copy.tsv"
fil_siq_inp="${dir_out}/siq.in_place.tsv"

fil_dry="${dir_out}/dry.tsv"
fil_dry_hdr="${dir_out}/dry.headered.tsv"
fil_dry_dat="${dir_in}/dry.data.tsv"
fil_dry_out="${dir_out}/dry.copy.tsv"

fil_old_amb="${dir_out}/old_ambiguous.tsv"
fil_bogus="${dir_out}/bogus.tsv"
fil_both="${dir_out}/both.tsv"
fil_unknown="${dir_out}/unknown.tsv"
fil_missing_dir="${dir_tmp}/missing/spike.tsv"

log_header_only="${dir_log}/header_only.log"
log_help="${dir_log}/help.log"
log_default="${dir_log}/default.log"
log_verbose="${dir_log}/verbose.log"
log_siq_copy="${dir_log}/siq_copy.log"
log_siq_in_place="${dir_log}/siq_in_place.log"
log_spike_idempotent="${dir_log}/spike_idempotent.log"
log_dry_header_only="${dir_log}/dry_header_only.log"
log_dry_in_place="${dir_log}/dry_in_place.log"
log_dry_copy="${dir_log}/dry_copy.log"
log_invalid_mode="${dir_log}/invalid_mode.log"
log_missing_output_mode="${dir_log}/missing_output_mode.log"
log_in_place_without_input="${dir_log}/in_place_without_input.log"
log_output_and_in_place="${dir_log}/output_and_in_place.log"
log_same_input_output="${dir_log}/same_input_output.log"
log_old_fil_out_existing="${dir_log}/old_fil_out_existing.log"
log_missing_mode_value="${dir_log}/missing_mode_value.log"
log_unknown_option="${dir_log}/unknown_option.log"
log_missing_output_dir="${dir_log}/missing_output_dir.log"

rm -rf "${dir_tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_log}"

require_files_nonempty "${scr_hdr}" || {
    finish
    exit $?
}

printf -v hdr_spk \
    '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "main_ip" "spike_ip" "main_in" "spike_in" "spike" "coef" \
    "num_mp" "num_sp" "num_mn" "num_sn"

printf -v hdr_siq \
    '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "fil_ip" "fil_in" "siq" "eqn" "mass_ip" "mass_in" \
    "vol_all" "vol_in" "dep_ip" "dep_in" "len_ip" "len_in"

printf -v row_spk \
    '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "/path/IP.sc.bam" "/path/IP.sp.bam" \
    "/path/in.sc.bam" "/path/in.sp.bam" \
    "chiprx_alpha_ratio" "2" "10" "5" "20" "10"

printf -v row_siq \
    '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "/path/IP.sc.bam" "/path/in.sc.bam" "0.5" "6nd" \
    "2.7" "72.5" "300" "20" "10" "20" "100" "100"


#  Header-only utility mode remains supported when no input is supplied
run_header_success \
    header_only \
    "${log_header_only}" \
    --mode spike \
    --fil_out "${fil_spk_new}"

assert_header_once \
    "${fil_spk_new}" \
    "${hdr_spk}" \
    "header_only"

if [[ "$(wc -l < "${fil_spk_new}")" -eq 1 ]]; then
    record_pass "write_header.sh header_only output has one line"
else
    record_fail "write_header.sh header_only output has extra lines"
fi

run_header_success \
    help \
    "${log_help}" \
    --help

assert_pattern_found \
    "${log_help}" \
    '^Usage:' \
    "write_header.sh help prints usage"

run_header_success \
    default \
    "${log_default}" \
    --fil_out "${fil_siq_def}"

assert_header_once \
    "${fil_siq_def}" \
    "${hdr_siq}" \
    "default"

run_header_success \
    verbose \
    "${log_verbose}" \
    --verbose \
    --mode spike \
    --fil_out "${fil_spk_vrb}"

assert_header_once \
    "${fil_spk_vrb}" \
    "${hdr_spk}" \
    "verbose"

assert_pattern_found \
    "${log_verbose}" \
    "^${hdr_spk}$" \
    "write_header.sh verbose prints spike header"


#  Explicit input/output and in-place modes
printf '%s\n' "${row_siq}" > "${fil_siq_dat}"

run_header_success \
    siq_copy \
    "${log_siq_copy}" \
    --mode siq \
    --fil_in "${fil_siq_dat}" \
    --fil_out "${fil_siq_cpy}"

assert_header_once \
    "${fil_siq_cpy}" \
    "${hdr_siq}" \
    "siq_copy"

assert_pattern_found \
    "${fil_siq_cpy}" \
    "^${row_siq}$" \
    "write_header.sh siq_copy keeps data row"

if [[ "$(sed -n '1p' "${fil_siq_dat}")" == "${row_siq}" ]]; then
    record_pass "write_header.sh siq_copy leaves input unchanged"
else
    record_fail "write_header.sh siq_copy modified input"
fi

printf '%s\n' "${row_siq}" > "${fil_siq_inp}"

run_header_success \
    siq_in_place \
    "${log_siq_in_place}" \
    --mode siq \
    --fil_in "${fil_siq_inp}" \
    --in_place

assert_header_once \
    "${fil_siq_inp}" \
    "${hdr_siq}" \
    "siq_in_place"

assert_pattern_found \
    "${fil_siq_inp}" \
    "^${row_siq}$" \
    "write_header.sh siq_in_place keeps data row"

{
    printf '%s\n' "${hdr_spk}"
    printf '%s\n' "${row_spk}"
} > "${fil_spk_hdr}"

run_header_success \
    spike_idempotent \
    "${log_spike_idempotent}" \
    --mode spike \
    --fil_in "${fil_spk_hdr}" \
    --in_place

assert_header_once \
    "${fil_spk_hdr}" \
    "${hdr_spk}" \
    "spike_idempotent"

if [[ "$(wc -l < "${fil_spk_hdr}")" -eq 2 ]]; then
    record_pass "write_header.sh spike_idempotent output has two lines"
else
    record_fail "write_header.sh spike_idempotent output changed line count"
fi


#  Dry-run modes report planned work without writing or modifying files
run_header_success \
    dry_header_only \
    "${log_dry_header_only}" \
    --dry_run \
    --mode spike \
    --fil_out "${fil_dry}"

if [[ ! -e "${fil_dry}" ]]; then
    record_pass "write_header.sh dry_header_only writes no output"
else
    record_fail "write_header.sh dry_header_only wrote output"
fi

assert_pattern_found \
    "${log_dry_header_only}" \
    'would create' \
    "write_header.sh dry_header_only reports planned creation"

{
    printf '%s\n' "${hdr_spk}"
    printf '%s\n' "${row_spk}"
} > "${fil_dry_hdr}"
run_header_success \
    dry_in_place \
    "${log_dry_in_place}" \
    --dry_run \
    --mode spike \
    --fil_in "${fil_dry_hdr}" \
    --in_place

assert_pattern_found \
    "${log_dry_in_place}" \
    'header already present' \
    "write_header.sh dry_in_place detects existing header"

if [[ "$(wc -l < "${fil_dry_hdr}")" -eq 2 ]]; then
    record_pass "write_header.sh dry_in_place preserves line count"
else
    record_fail "write_header.sh dry_in_place changed line count"
fi

printf '%s\n' "${row_siq}" > "${fil_dry_dat}"
run_header_success \
    dry_copy \
    "${log_dry_copy}" \
    --dry_run \
    --mode siq \
    --fil_in "${fil_dry_dat}" \
    --fil_out "${fil_dry_out}"

assert_pattern_found \
    "${log_dry_copy}" \
    'would prepend header' \
    "write_header.sh dry_copy reports headered copy"

if [[ ! -e "${fil_dry_out}" ]]; then
    record_pass "write_header.sh dry_copy writes no output"
else
    record_fail "write_header.sh dry_copy wrote output"
fi

if [[ "$(sed -n '1p' "${fil_dry_dat}")" == "${row_siq}" ]]; then
    record_pass "write_header.sh dry_copy preserves input"
else
    record_fail "write_header.sh dry_copy modified input"
fi


#  Invalid mode and invalid argument combinations should fail clearly
run_header_failure \
    invalid_mode \
    "${log_invalid_mode}" \
    "must be.*'siq' or 'spike'" \
    --mode bogus \
    --fil_out "${fil_bogus}"

run_header_failure \
    missing_output_mode \
    "${log_missing_output_mode}" \
    "requires either '--fil_out' or '--in_place'" \
    --mode spike \
    --fil_in "${fil_spk_hdr}"

run_header_failure \
    in_place_without_input \
    "${log_in_place_without_input}" \
    "'--in_place' requires '--fil_in'" \
    --mode spike \
    --in_place

run_header_failure \
    output_and_in_place \
    "${log_output_and_in_place}" \
    "use either '--fil_out' or '--in_place'" \
    --mode spike \
    --fil_in "${fil_spk_hdr}" \
    --fil_out "${fil_both}" \
    --in_place

run_header_failure \
    same_input_output \
    "${log_same_input_output}" \
    "same path; use '--in_place'" \
    --mode spike \
    --fil_in "${fil_spk_hdr}" \
    --fil_out "${fil_spk_hdr}"

printf '%s\n' "${row_spk}" > "${fil_old_amb}"
run_header_failure \
    old_fil_out_existing \
    "${log_old_fil_out_existing}" \
    "already exists in header-only mode" \
    --mode spike \
    --fil_out "${fil_old_amb}"

run_header_failure \
    missing_mode_value \
    "${log_missing_mode_value}" \
    "option '--mode' requires a value" \
    --mode

assert_pattern_found \
    "${log_missing_mode_value}" \
    '^Usage:' \
    "write_header.sh missing mode value prints usage"

run_header_failure \
    unknown_option \
    "${log_unknown_option}" \
    "unknown option/parameter passed: '--bogus'" \
    --bogus \
    --fil_out "${fil_unknown}"

run_header_failure \
    missing_output_dir \
    "${log_missing_output_dir}" \
    "directory for 'dir_out' does not exist" \
    --mode spike \
    --fil_out "${fil_missing_dir}"

finish
