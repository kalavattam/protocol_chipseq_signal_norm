#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_write_header.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="write scaling-factor header"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

scr_hdr="${ROOT_REPO}/scripts/write_header.sh"
tmp="${TEST_DIR_TMP}/write_header"
dir_log="${TEST_DIR_LOG}/write_header"
mkdir -p "${tmp}/out" "${dir_log}"
rm -f "${tmp}/out/"*

hdr_spk=$'main_ip\tspike_ip\tmain_in\tspike_in\tspike\tcoef\tnum_mp\tnum_sp\tnum_mn\tnum_sn'
hdr_siq=$'fil_ip\tfil_in\tsiq\teqn\tmass_ip\tmass_in\tvol_all\tvol_in\tdep_ip\tdep_in\tlen_ip\tlen_in'

row_spk=$'/path/IP.sc.bam\t/path/IP.sp.bam\t/path/in.sc.bam\t/path/in.sp.bam\tchiprx_alpha_ratio\t2\t10\t5\t20\t10'
row_siq=$'/path/IP.sc.bam\t/path/in.sc.bam\t0.5\t6nd\t2.7\t72.5\t300\t20\t10\t20\t100\t100'


#  Assert write_header.sh creates or prepends one expected header
function assert_write_header() {
    local mode="${1:-}"
    local fil_out="${2:-}"
    local header="${3:-}"
    local label="${4:-}"
    local log="${5:-}"
    local num_hdr=""

    if \
        run_capture \
            "write ${mode} scaling-factor header" \
            "${log}" \
            "${TEST_BASH}" "${scr_hdr}" \
                --mode "${mode}" \
                --fil_out "${fil_out}"
    then
        record_pass "write_header.sh ${label} exits 0"
    else
        record_fail \
            "write_header.sh ${label} failed; see" \
            "$(print_relpath "${log}")"
        return
    fi

    assert_file_nonempty \
        "${fil_out}" \
        "write_header.sh ${label} output"

    assert_pattern_found \
        "${fil_out}" \
        "^${header}$" \
        "write_header.sh ${label} writes expected header"

    num_hdr="$(grep -c "^${header}$" "${fil_out}" || true)"
    if [[ "${num_hdr}" -eq 1 ]]; then
        record_pass "write_header.sh ${label} writes one header"
    else
        record_fail "write_header.sh ${label} header count is not one"
    fi
}


#  Create a header-only spike table when the output does not exist
fil_spk_new="${tmp}/out/spike.header_only.tsv"
assert_write_header \
    spike \
    "${fil_spk_new}" \
    "${hdr_spk}" \
    "spike header-only" \
    "${dir_log}/spike_header_only.log"

if [[ "$(wc -l < "${fil_spk_new}")" -eq 1 ]]; then
    record_pass "write_header.sh spike header-only output has one line"
else
    record_fail "write_header.sh spike header-only output has extra lines"
fi


#  Help should print usage and exit cleanly
log="${dir_log}/help.log"
if \
    run_capture \
        "write scaling-factor header help" \
        "${log}" \
        "${TEST_BASH}" "${scr_hdr}" --help
then
    record_pass "write_header.sh --help exits 0"
    assert_pattern_found \
        "${log}" \
        '^Usage:' \
        "write_header.sh --help prints usage"
else
    record_fail "write_header.sh --help failed; see $(print_relpath "${log}")"
fi


#  Default mode should write a siQ-ChIP header
fil_siq_def="${tmp}/out/siq.default.tsv"
log="${dir_log}/siq_def.log"
if \
    run_capture \
        "write default scaling-factor header" \
        "${log}" \
        "${TEST_BASH}" "${scr_hdr}" \
            --fil_out "${fil_siq_def}"
then
    record_pass "write_header.sh default mode exits 0"
else
    record_fail \
        "write_header.sh default mode failed; see $(print_relpath "${log}")"
fi

assert_pattern_found \
    "${fil_siq_def}" \
    "^${hdr_siq}$" \
    "write_header.sh default mode writes siq header"


#  Verbose mode should print the selected header while writing output
fil_spk_verbose="${tmp}/out/spike.verbose.tsv"
log="${dir_log}/spike_verbose.log"
if \
    run_capture \
        "write scaling-factor header verbose" \
        "${log}" \
        "${TEST_BASH}" "${scr_hdr}" \
            --verbose \
            --mode spike \
            --fil_out "${fil_spk_verbose}"
then
    record_pass "write_header.sh --verbose exits 0"
else
    record_fail \
        "write_header.sh --verbose failed; see $(print_relpath "${log}")"
fi

assert_pattern_found \
    "${fil_spk_verbose}" \
    "^${hdr_spk}$" \
    "write_header.sh --verbose writes spike header"

assert_pattern_found \
    "${log}" \
    "^${hdr_spk}$" \
    "write_header.sh --verbose prints spike header"


#  Prepend a siQ-ChIP header to an existing data-only table
fil_siq_dat="${tmp}/out/siq.data.tsv"
printf '%s\n' "${row_siq}" > "${fil_siq_dat}"

assert_write_header \
    siq \
    "${fil_siq_dat}" \
    "${hdr_siq}" \
    "siq prepend" \
    "${dir_log}/siq_prepend.log"

assert_pattern_found \
    "${fil_siq_dat}" \
    "^${row_siq}$" \
    "write_header.sh siq prepend keeps data row"

if [[ "$(sed -n '1p' "${fil_siq_dat}")" == "${hdr_siq}" ]]; then
    record_pass "write_header.sh siq prepend puts header first"
else
    record_fail "write_header.sh siq prepend did not put header first"
fi


#  Leave an already headered spike table unchanged instead of duplicating header
fil_spk_hdr="${tmp}/out/spike.headered.tsv"
{
    printf '%s\n' "${hdr_spk}"
    printf '%s\n' "${row_spk}"
} > "${fil_spk_hdr}"

assert_write_header \
    spike \
    "${fil_spk_hdr}" \
    "${hdr_spk}" \
    "spike idempotent" \
    "${dir_log}/spike_idempotent.log"

if [[ "$(wc -l < "${fil_spk_hdr}")" -eq 2 ]]; then
    record_pass "write_header.sh spike idempotent output has two lines"
else
    record_fail "write_header.sh spike idempotent output changed line count"
fi


#  Dry-run reports the planned action without creating an output file
fil_dry="${tmp}/out/dry.tsv"
log="${dir_log}/dry_run.log"
if \
    run_capture \
        "write scaling-factor header dry-run" \
        "${log}" \
        "${TEST_BASH}" "${scr_hdr}" \
            --dry_run \
            --mode spike \
            --fil_out "${fil_dry}"
then
    record_pass "write_header.sh --dry_run exits 0"
else
    record_fail "write_header.sh --dry_run failed; see $(print_relpath "${log}")"
fi

if [[ ! -e "${fil_dry}" ]]; then
    record_pass "write_header.sh --dry_run writes no output"
else
    record_fail "write_header.sh --dry_run wrote output"
fi

assert_pattern_found \
    "${log}" \
    'Dry run:' \
    "write_header.sh --dry_run reports planned action"


#  Dry-run should not modify an existing headered file
fil_dry_hdr="${tmp}/out/dry.headered.tsv"
log="${dir_log}/dry_run_headered.log"
{
    printf '%s\n' "${hdr_spk}"
    printf '%s\n' "${row_spk}"
} > "${fil_dry_hdr}"

if \
    run_capture \
        "write scaling-factor header dry-run already headered" \
        "${log}" \
        "${TEST_BASH}" "${scr_hdr}" \
            --dry_run \
            --mode spike \
            --fil_out "${fil_dry_hdr}"
then
    record_pass "write_header.sh --dry_run headered exits 0"
else
    record_fail \
        "write_header.sh --dry_run headered failed; see" \
        "$(print_relpath "${log}")"
fi

assert_pattern_found \
    "${log}" \
    'header already present' \
    "write_header.sh --dry_run detects existing header"

if [[ "$(wc -l < "${fil_dry_hdr}")" -eq 2 ]]; then
    record_pass "write_header.sh --dry_run headered preserves line count"
else
    record_fail "write_header.sh --dry_run headered changed line count"
fi


#  Dry-run should not modify an existing data-only file
fil_dry_dat="${tmp}/out/dry.data.tsv"
log="${dir_log}/dry_run_data.log"
printf '%s\n' "${row_siq}" > "${fil_dry_dat}"

if \
    run_capture \
        "write scaling-factor header dry-run data-only" \
        "${log}" \
        "${TEST_BASH}" "${scr_hdr}" \
            --dry_run \
            --mode siq \
            --fil_out "${fil_dry_dat}"
then
    record_pass "write_header.sh --dry_run data-only exits 0"
else
    record_fail \
        "write_header.sh --dry_run data-only failed; see" \
        "$(print_relpath "${log}")"
fi

assert_pattern_found \
    "${log}" \
    'would prepend header' \
    "write_header.sh --dry_run reports header prepend"

if [[ "$(sed -n '1p' "${fil_dry_dat}")" == "${row_siq}" ]]; then
    record_pass "write_header.sh --dry_run data-only preserves first row"
else
    record_fail "write_header.sh --dry_run data-only modified first row"
fi


#  Invalid mode fails clearly
log="${dir_log}/invalid_mode.log"
if \
    run_capture \
        "write scaling-factor header invalid mode" \
        "${log}" \
        "${TEST_BASH}" "${scr_hdr}" \
            --mode bogus \
            --fil_out "${tmp}/out/bogus.tsv"
then
    record_fail "write_header.sh invalid mode unexpectedly succeeded"
else
    assert_pattern_found \
        "${log}" \
        "must be.*'siq' or 'spike'" \
        "write_header.sh rejects invalid mode"
fi


#  Missing required arguments and invalid options should fail clearly
log="${dir_log}/missing_fil_out.log"
if \
    run_capture \
        "write scaling-factor header missing output" \
        "${log}" \
        "${TEST_BASH}" "${scr_hdr}" \
            --mode spike
then
    record_fail "write_header.sh missing --fil_out unexpectedly succeeded"
else
    assert_pattern_found \
        "${log}" \
        "'fil_out' is empty or unset" \
        "write_header.sh rejects missing fil_out"
fi

log="${dir_log}/missing_mode_value.log"
if \
    run_capture \
        "write scaling-factor header missing mode value" \
        "${log}" \
        "${TEST_BASH}" "${scr_hdr}" \
            --mode
then
    record_fail "write_header.sh missing --mode value unexpectedly succeeded"
else
    assert_pattern_found \
        "${log}" \
        "option '--mode' requires a value" \
        "write_header.sh rejects missing mode value"

    assert_pattern_found \
        "${log}" \
        '^Usage:' \
        "write_header.sh missing mode value prints usage"
fi

log="${dir_log}/unknown_option.log"
if \
    run_capture \
        "write scaling-factor header unknown option" \
        "${log}" \
        "${TEST_BASH}" "${scr_hdr}" \
            --bogus \
            --fil_out "${tmp}/out/unknown.tsv"
then
    record_fail "write_header.sh unknown option unexpectedly succeeded"
else
    assert_pattern_found \
        "${log}" \
        "unknown option/parameter passed: '--bogus'" \
        "write_header.sh rejects unknown option"
fi

log="${dir_log}/missing_output_dir.log"
if \
    run_capture \
        "write scaling-factor header missing output directory" \
        "${log}" \
        "${TEST_BASH}" "${scr_hdr}" \
            --mode spike \
            --fil_out "${tmp}/missing/spike.tsv"
then
    record_fail \
        "write_header.sh missing output directory unexpectedly succeeded"
else
    assert_pattern_found \
        "${log}" \
        "directory for 'dir_out' does not exist" \
        "write_header.sh rejects missing output directory"
fi

finish
