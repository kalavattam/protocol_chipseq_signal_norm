#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_combine_parts_scaling_factor.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="combine scaling-factor parts"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

scr_cmb="${ROOT_REPO}/scripts/combine_parts_scaling_factor.sh"
dir_fix="${ROOT_REPO}/tests/calculate_scaling_factor/fixtures"
dir_prt="${dir_fix}/parts"
dir_log="${TEST_DIR_LOG}/combine_parts_scaling_factor"
tmp="${TEST_DIR_TMP}/combine_parts_scaling_factor"
dir_prt_dry="${tmp}/parts_dry_no_parts"
dir_prt_out="${tmp}/parts_output_input"
dir_prt_tmp="${tmp}/parts_no_parts"

spk_0="${dir_prt}/example_scaling_factors.spike.tsv.part.000000"
spk_2="${dir_prt}/example_scaling_factors.spike.tsv.part.000002"
siq_0="${dir_prt}/example_scaling_factors.siq.tsv.part.000000"
siq_2="${dir_prt}/example_scaling_factors.siq.tsv.part.000002"
fil_bad="${dir_prt}/malformed_scaling_factors.spike.tsv.part.000003"
fil_hdr="${dir_prt}/header_scaling_factors.spike.tsv.part.000004"
dup_idx_a="${dir_prt}/duplicate_index_A.spike.tsv.part.000005"
dup_idx_b="${dir_prt}/duplicate_index_B.spike.tsv.part.000005"


#  Assert one successful combination with reverse-ordered input files
function assert_combined_mode() {
    local mode="${1:-}"
    local csv_in="${2:-}"
    local fil_out="${3:-}"
    local row_1="${4:-}"
    local row_2="${5:-}"
    local fil_log="${6:-}"

    if \
        run_capture \
            "combine scaling-factor ${mode} parts" \
            "${fil_log}" \
            "${TEST_BASH}" "${scr_cmb}" \
                --mode "${mode}" \
                --csv_infile "${csv_in}" \
                --fil_out "${fil_out}"
    then
        record_pass "combine_parts_scaling_factor.sh ${mode} exits 0"
    else
        record_fail \
            "combine_parts_scaling_factor.sh ${mode} failed; see" \
            "$(print_relpath "${fil_log}")"
        return
    fi

    assert_file_nonempty \
        "${fil_out}" \
        "combined scaling-factor ${mode} TSV"

    assert_pattern_found \
        "${fil_out}" \
        "^${row_1}" \
        "combined scaling-factor ${mode} TSV has first numeric-index row"

    assert_pattern_found \
        "${fil_out}" \
        "^${row_2}" \
        "combined scaling-factor ${mode} TSV has second numeric-index row"

    if [[ "$(sed -n '1p' "${fil_out}")" == "${row_1}"$'\t'* ]]; then
        record_pass "combined scaling-factor ${mode} TSV sorts first row"
    else
        record_fail \
            "combined scaling-factor ${mode} TSV first row is unsorted"
    fi

    if [[ "$(sed -n '2p' "${fil_out}")" == "${row_2}"$'\t'* ]]; then
        record_pass "combined scaling-factor ${mode} TSV sorts second row"
    else
        record_fail \
            "combined scaling-factor ${mode} TSV second row is unsorted"
    fi
}


#  Assert one expected combine_parts_scaling_factor.sh failure
function assert_combine_fails() {
    local label="${1:-expected failure}"
    local fil_log="${2:-}"
    local patn="${3:-}"
    shift 3 || true

    if \
        run_capture \
            "combine scaling-factor ${label}" \
            "${fil_log}" \
            "${TEST_BASH}" "${scr_cmb}" \
                "$@"
    then
        record_fail "combine_parts_scaling_factor.sh ${label} unexpectedly pass"
    else
        assert_pattern_found \
            "${fil_log}" \
            "${patn}" \
            "combine_parts_scaling_factor.sh ${label}"
    fi
}


print_section "${TEST_NAME}"

mkdir -p "${tmp}/out" "${dir_log}"
rm -f "${tmp}/out/"*
rm -rf \
    "${dir_prt_dry}" \
    "${dir_prt_out}" \
    "${dir_prt_tmp}"

if ! \
    require_files_nonempty \
        "${spk_0}" \
        "${spk_2}" \
        "${siq_0}" \
        "${siq_2}" \
        "${fil_bad}" \
        "${fil_hdr}" \
        "${dup_idx_a}" \
        "${dup_idx_b}"
then
    finish
    exit $?
fi


#  Help should print usage and exit cleanly
fil_log="${dir_log}/help.log"
if \
    run_capture \
        "combine scaling-factor help" \
        "${fil_log}" \
        "${TEST_BASH}" "${scr_cmb}" --help
then
    record_pass "combine_parts_scaling_factor.sh --help exits 0"
    assert_pattern_found \
        "${fil_log}" \
        '^Usage:' \
        "combine_parts_scaling_factor.sh --help prints usage"
else
    record_fail \
        "combine_parts_scaling_factor.sh --help failed; see" \
        "$(print_relpath "${fil_log}")"
fi


#  Combine spike-in and siQ-ChIP parts passed in reverse order
fil_out_spk="${tmp}/out/scaling.spike.tsv"
assert_combined_mode \
    spike \
    "${spk_2},${spk_0}" \
    "${fil_out_spk}" \
    '/path/to/IP_WT_G1_Hho1_6336.sc.bam' \
    '/path/to/IP_WT_G1_Hho1_6337.sc.bam' \
    "${dir_log}/spike.log"

assert_pattern_found \
    "${fil_out_spk}" \
    $'\t2.081558847888101748679901\tchiprx_alpha_ratio\t13492920\t217340\t12851824\t452406$' \
    "combined scaling-factor spike TSV records first default coefficient"

assert_pattern_found \
    "${fil_out_spk}" \
    $'\t2.898996981538645822951139\tchiprx_alpha_ratio\t13655994\t116947\t12030091\t339029$' \
    "combined scaling-factor spike TSV records second default coefficient"

assert_pattern_absent \
    "${fil_out_spk}" \
    $'^main_ip\tspike_ip\tmain_in\tspike_in\tspike\tcoef' \
    "combined scaling-factor spike TSV stays data-only"

fil_out_siq="${tmp}/out/scaling.siq.tsv"
assert_combined_mode \
    siq \
    "${siq_2},${siq_0}" \
    "${fil_out_siq}" \
    '/path/to/IP_WT_G1_Hho1_6336.sc.bam' \
    '/path/to/IP_WT_G1_Hho1_6337.sc.bam' \
    "${dir_log}/siq.log"

assert_pattern_absent \
    "${fil_out_siq}" \
    $'^fil_ip\tfil_in\tsiq\teqn\tmass_ip\tmass_in' \
    "combined scaling-factor siq TSV stays data-only"


#  Combine one spike-in part into a one-row table
fil_out_one="${tmp}/out/scaling.one_part.tsv"
if \
    run_capture \
        "combine scaling-factor one part" \
        "${dir_log}/one_part.log" \
        "${TEST_BASH}" "${scr_cmb}" \
            --mode spike \
            --csv_infile "${spk_0}" \
            --fil_out "${fil_out_one}"
then
    record_pass "combine_parts_scaling_factor.sh one-part exits 0"
else
    record_fail \
        "combine_parts_scaling_factor.sh one-part failed; see" \
        "$(print_relpath "${dir_log}/one_part.log")"
fi

assert_file_nonempty \
    "${fil_out_one}" \
    "combined scaling-factor one-part TSV"

if [[ "$(wc -l < "${fil_out_one}")" -eq 1 ]]; then
    record_pass "combine_parts_scaling_factor.sh one-part writes one row"
else
    record_fail "combine_parts_scaling_factor.sh one-part row count changed"
fi

assert_pattern_found \
    "${fil_out_one}" \
    '^/path/to/IP_WT_G1_Hho1_6336.sc.bam' \
    "combine_parts_scaling_factor.sh one-part writes expected row"


#  Confirm part files are retained by default
assert_file_nonempty \
    "${spk_0}" \
    "default combination retains first spike-in part"

assert_file_nonempty \
    "${spk_2}" \
    "default combination retains second spike-in part"


#  Confirm dry-run validation writes no final table
fil_out_dry="${tmp}/out/scaling.dry.tsv"
fil_log="${dir_log}/dry_run.log"
if \
    run_capture \
        "combine scaling-factor dry-run" \
        "${fil_log}" \
        "${TEST_BASH}" "${scr_cmb}" \
            --dry_run \
            --mode spike \
            --csv_infile "${spk_2},${spk_0}" \
            --fil_out "${fil_out_dry}"
then
    record_pass "combine_parts_scaling_factor.sh --dry_run exits 0"
else
    record_fail \
        "combine_parts_scaling_factor.sh --dry_run failed; see" \
        "$(print_relpath "${fil_log}")"
fi

if [[ ! -e "${fil_out_dry}" ]]; then
    record_pass "combine_parts_scaling_factor.sh --dry_run writes no output"
else
    record_fail "combine_parts_scaling_factor.sh --dry_run wrote output"
fi


#  Confirm dry-run no_parts reports removal but retains copied inputs
mkdir -p "${dir_prt_dry}"

spk_dry_0="${dir_prt_dry}/scaling.spike.tsv.part.000000"
spk_dry_2="${dir_prt_dry}/scaling.spike.tsv.part.000002"
cp "${spk_0}" "${spk_dry_0}"
cp "${spk_2}" "${spk_dry_2}"

fil_log="${dir_log}/dry_run_no_parts.log"
if \
    run_capture \
        "combine scaling-factor dry-run no-parts" \
        "${fil_log}" \
        "${TEST_BASH}" "${scr_cmb}" \
            --dry_run \
            --mode spike \
            --csv_infile "${spk_dry_2},${spk_dry_0}" \
            --fil_out "${tmp}/out/scaling.dry_no_parts.tsv" \
            --no_parts
then
    record_pass "combine_parts_scaling_factor.sh --dry_run --no_parts exits 0"
else
    record_fail \
        "combine_parts_scaling_factor.sh --dry_run --no_parts failed; see" \
        "$(print_relpath "${fil_log}")"
fi

assert_pattern_found \
    "${fil_log}" \
    'would remove validated part files' \
    "combine_parts_scaling_factor.sh --dry_run --no_parts reports removal"

if [[ -e "${spk_dry_0}" && -e "${spk_dry_2}" ]]; then
    record_pass "combine_parts_scaling_factor.sh --dry_run --no_parts keeps inputs"
else
    record_fail "combine_parts_scaling_factor.sh --dry_run --no_parts removed inputs"
fi


#  Confirm invalid inputs fail clearly
assert_combine_fails \
    "rejects missing mode" \
    "${dir_log}/missing_mode.log" \
    "'mode' is empty or unset" \
    --csv_infile "${spk_0}" \
    --fil_out "${tmp}/out/missing_mode.tsv"

assert_combine_fails \
    "rejects invalid mode" \
    "${dir_log}/invalid_mode.log" \
    "'--mode' must be 'siq' or 'spike'" \
    --mode bogus \
    --csv_infile "${spk_0}" \
    --fil_out "${tmp}/out/invalid_mode.tsv"

assert_combine_fails \
    "rejects missing csv_infile" \
    "${dir_log}/missing_csv_infile.log" \
    "'csv_infile' is empty or unset" \
    --mode spike \
    --fil_out "${tmp}/out/missing_csv_infile.tsv"

assert_combine_fails \
    "rejects unknown option" \
    "${dir_log}/unknown_option.log" \
    "unknown option/parameter passed: '--bogus'" \
    --bogus \
    --mode spike \
    --csv_infile "${spk_0}" \
    --fil_out "${tmp}/out/unknown_option.tsv"

assert_combine_fails \
    "rejects duplicate input path" \
    "${dir_log}/duplicate_input_path.log" \
    "duplicate input part file" \
    --mode spike \
    --csv_infile "${spk_0},${spk_0}" \
    --fil_out "${tmp}/out/duplicate_path.tsv"

mkdir -p "${dir_prt_out}"

spk_out_in="${dir_prt_out}/scaling.out.spike.tsv.part.000006"
cp "${spk_0}" "${spk_out_in}"

assert_combine_fails \
    "rejects output that is also input" \
    "${dir_log}/output_is_input.log" \
    "'--fil_out' is also an input part file" \
    --mode spike \
    --csv_infile "${spk_out_in}" \
    --fil_out "${spk_out_in}" \
    --force

assert_combine_fails \
    "rejects duplicate numeric index" \
    "${dir_log}/duplicate_numeric_index.log" \
    "duplicate numeric part index '5'" \
    --mode spike \
    --csv_infile "${dup_idx_a},${dup_idx_b}" \
    --fil_out "${tmp}/out/duplicate_index.tsv"

assert_combine_fails \
    "rejects header row" \
    "${dir_log}/header_row.log" \
    "appears to contain a header" \
    --mode spike \
    --csv_infile "${fil_hdr}" \
    --fil_out "${tmp}/out/header_row.tsv"

assert_combine_fails \
    "rejects malformed fields" \
    "${dir_log}/malformed_fields.log" \
    "expected '10' for mode 'spike'" \
    --mode spike \
    --csv_infile "${fil_bad}" \
    --fil_out "${tmp}/out/malformed.tsv"

assert_combine_fails \
    "rejects existing output without force" \
    "${dir_log}/existing_output.log" \
    'output file already exists' \
    --mode spike \
    --csv_infile "${spk_2},${spk_0}" \
    --fil_out "${fil_out_spk}"


#  Confirm force replaces output and no_parts removes only validated inputs
mkdir -p "${dir_prt_tmp}"

spk_tmp_0="${dir_prt_tmp}/scaling.spike.tsv.part.000000"
spk_tmp_2="${dir_prt_tmp}/scaling.spike.tsv.part.000002"
cp "${spk_0}" "${spk_tmp_0}"
cp "${spk_2}" "${spk_tmp_2}"

fil_out_nop="${tmp}/out/scaling.no_parts.tsv"
fil_log="${dir_log}/no_parts.log"
if \
    run_capture \
        "combine scaling-factor no-parts" \
        "${fil_log}" \
        "${TEST_BASH}" "${scr_cmb}" \
            --mode spike \
            --csv_infile "${spk_tmp_2},${spk_tmp_0}" \
            --fil_out "${fil_out_nop}" \
            --no_parts
then
    record_pass "combine_parts_scaling_factor.sh --no_parts exits 0"
else
    record_fail \
        "combine_parts_scaling_factor.sh --no_parts failed; see" \
        "$(print_relpath "${fil_log}")"
fi

assert_file_nonempty \
    "${fil_out_nop}" \
    "combined scaling-factor --no_parts TSV"

if [[ ! -e "${spk_tmp_0}" && ! -e "${spk_tmp_2}" ]]; then
    record_pass "combine_parts_scaling_factor.sh --no_parts removes inputs"
else
    record_fail "combine_parts_scaling_factor.sh --no_parts retained inputs"
fi

fil_log="${dir_log}/force.log"
if \
    run_capture \
        "combine scaling-factor force" \
        "${fil_log}" \
        "${TEST_BASH}" "${scr_cmb}" \
            --mode spike \
            --csv_infile "${spk_2},${spk_0}" \
            --fil_out "${fil_out_spk}" \
            --force
then
    record_pass "combine_parts_scaling_factor.sh --force exits 0"
else
    record_fail \
        "combine_parts_scaling_factor.sh --force failed; see" \
        "$(print_relpath "${fil_log}")"
fi

assert_pattern_found \
    "${fil_out_spk}" \
    '^/path/to/IP_WT_G1_Hho1_6336.sc.bam' \
    "combine_parts_scaling_factor.sh --force rewrites expected first row"

assert_pattern_found \
    "${fil_out_spk}" \
    '^/path/to/IP_WT_G1_Hho1_6337.sc.bam' \
    "combine_parts_scaling_factor.sh --force rewrites expected second row"

finish
