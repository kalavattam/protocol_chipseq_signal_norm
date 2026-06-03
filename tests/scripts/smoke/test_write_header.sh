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

finish
