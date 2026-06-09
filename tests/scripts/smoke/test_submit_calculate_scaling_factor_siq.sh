#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_calculate_scaling_factor_siq.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="submit calculate-scaling-factor siQ"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

#  Define fixture and output paths for the direct submit-worker smoke test
scr_sub="${ROOT_REPO}/scripts/submit_calculate_scaling_factor.sh"
cfg_met="${ROOT_REPO}/data/raw/docs/parse_metadata_siqchip.yml"
dir_fix="${ROOT_REPO}/tests/calculate_scaling_factor/fixtures"
dir_bam="${dir_fix}/bam"
dir_bam_se="${dir_bam}/se"
dir_bam_pe="${dir_bam}/pe"
dir_met="${dir_fix}/metadata"
tbl_met="${dir_met}/measurements_siqchip.tsv"

bam_pe_mip="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sc.bam"
bam_pe_min="${dir_bam_pe}/in_WT_G1_Hho1_6336.sc.bam"
bam_se_mip="${dir_bam_se}/IP_WT_log_Brn1_rep1.sc.bam"
bam_se_min="${dir_bam_se}/in_WT_log_Brn1_rep1.sc.bam"

tmp="${TEST_DIR_TMP}/submit_calculate_scaling_factor_siq"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/calculate_scaling_factor"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

if ! \
    require_files_nonempty \
        "${scr_sub}" \
        "${cfg_met}" \
        "${tbl_met}" \
        "${bam_pe_mip}" \
        "${bam_pe_mip}.bai" \
        "${bam_pe_min}" \
        "${bam_pe_min}.bai" \
        "${bam_se_mip}" \
        "${bam_se_mip}.bai" \
        "${bam_se_min}" \
        "${bam_se_min}.bai"
then
    finish
    exit $?
fi


#  PE BAM input should compute lengths from paired-end fragment sizes
fil_out="${dir_out}/scaling.submit.pe.siq.tsv"
fil_prt="${fil_out}.part.000009"
fil_log="${dir_log}/submit_siq_pe.log"

printf -v row_exp '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${bam_pe_mip}" \
    "${bam_pe_min}" \
    '0.002660098522167487697376' \
    '6nd' \
    '2.7' \
    '72.5' \
    '300' \
    '20' \
    '3' \
    '2' \
    '20' \
    '20'

# shellcheck disable=SC2154
if \
    run_capture \
        "submit calculate-scaling-factor siQ PE" \
        "${fil_log}" \
        "${TEST_BASH}" "${scr_sub}" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --mode siq \
            --csv_mip "${bam_pe_mip}" \
            --csv_min "${bam_pe_min}" \
            --aln_typ auto \
            --fil_out "${fil_out}" \
            --idx_out 9 \
            --tbl_met "${tbl_met}" \
            --cfg_met "${cfg_met}" \
            --eqn 6nd \
            --err_out "${dir_err}" \
            --nam_job test_submit_calculate_scaling_factor_siq_pe
then
    record_pass "submit_calculate_scaling_factor.sh siQ PE exits 0"
else
    record_fail \
        "submit_calculate_scaling_factor.sh siQ PE failed; see" \
        "$(print_relpath "${fil_log}")"
fi

assert_file_nonempty \
    "${fil_prt}" \
    "submit scaling-factor siQ PE indexed part"

assert_pattern_absent \
    "${fil_prt}" \
    $'^fil_ip\tfil_in\tsiq\teqn' \
    "submit scaling-factor siQ PE part stays data-only"

if [[ ! -e "${fil_out}" ]]; then
    record_pass "submit_calculate_scaling_factor.sh siQ writes no final TSV"
else
    record_fail "submit_calculate_scaling_factor.sh siQ unexpectedly wrote TSV"
fi

if [[ "$(cat "${fil_prt}")" == "${row_exp}" ]]; then
    record_pass "submit scaling-factor siQ PE part has expected row"
else
    record_fail \
        "submit scaling-factor siQ PE part has unexpected row; see" \
        "$(print_relpath "${fil_prt}")"
fi

assert_pattern_found \
    "${fil_log}" \
    'idx_out=9' \
    "submit scaling-factor siQ PE uses idx_out=9"

assert_pattern_found \
    "${fil_log}" \
    'typ_ip=pe' \
    "submit scaling-factor siQ PE auto-detects IP as PE"

assert_pattern_found \
    "${fil_log}" \
    'dep_ip=3' \
    "submit scaling-factor siQ PE counts IP fragments"

assert_pattern_found \
    "${fil_log}" \
    'dep_in=2' \
    "submit scaling-factor siQ PE counts input fragments"

assert_pattern_found \
    "${fil_log}" \
    'len_ip=20' \
    "submit scaling-factor siQ PE computes IP fragment length"

assert_pattern_found \
    "${fil_log}" \
    'len_in=20' \
    "submit scaling-factor siQ PE computes input fragment length"


#  SE BAM input should use an explicit default fragment length
fil_out="${dir_out}/scaling.submit.se.siq.tsv"
fil_prt="${fil_out}.part.000010"
fil_log="${dir_log}/submit_siq_se.log"

printf -v row_exp '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${bam_se_mip}" \
    "${bam_se_min}" \
    '0.004761904761904762334312' \
    '6nd' \
    '4' \
    '60' \
    '300' \
    '20' \
    '3' \
    '2' \
    '150' \
    '150'

if \
    run_capture \
        "submit calculate-scaling-factor siQ SE" \
        "${fil_log}" \
        "${TEST_BASH}" "${scr_sub}" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --mode siq \
            --csv_mip "${bam_se_mip}" \
            --csv_min "${bam_se_min}" \
            --aln_typ auto \
            --len_def 150 \
            --fil_out "${fil_out}" \
            --idx_out 10 \
            --tbl_met "${tbl_met}" \
            --cfg_met "${cfg_met}" \
            --eqn 6nd \
            --err_out "${dir_err}" \
            --nam_job test_submit_calculate_scaling_factor_siq_se
then
    record_pass "submit_calculate_scaling_factor.sh siQ SE exits 0"
else
    record_fail \
        "submit_calculate_scaling_factor.sh siQ SE failed; see" \
        "$(print_relpath "${fil_log}")"
fi

assert_file_nonempty \
    "${fil_prt}" \
    "submit scaling-factor siQ SE indexed part"

if [[ "$(cat "${fil_prt}")" == "${row_exp}" ]]; then
    record_pass "submit scaling-factor siQ SE part has expected row"
else
    record_fail \
        "submit scaling-factor siQ SE part has unexpected row; see" \
        "$(print_relpath "${fil_prt}")"
fi

assert_pattern_found \
    "${fil_log}" \
    'idx_out=10' \
    "submit scaling-factor siQ SE uses idx_out=10"

assert_pattern_found \
    "${fil_log}" \
    'typ_ip=se' \
    "submit scaling-factor siQ SE auto-detects IP as SE"

assert_pattern_found \
    "${fil_log}" \
    'dep_ip=3' \
    "submit scaling-factor siQ SE counts IP alignments"

assert_pattern_found \
    "${fil_log}" \
    'dep_in=2' \
    "submit scaling-factor siQ SE counts input alignments"

assert_pattern_found \
    "${fil_log}" \
    'len_ip=150' \
    "submit scaling-factor siQ SE uses default IP fragment length"

assert_pattern_found \
    "${fil_log}" \
    'len_in=150' \
    "submit scaling-factor siQ SE uses default input fragment length"


#  Explicit SE mode without length information should fail before processing
fil_log="${dir_log}/submit_siq_se_missing_length.log"
if \
    run_capture \
        "submit calculate-scaling-factor siQ SE missing length" \
        "${fil_log}" \
        "${TEST_BASH}" "${scr_sub}" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --mode siq \
            --csv_mip "${bam_se_mip}" \
            --csv_min "${bam_se_min}" \
            --aln_typ se \
            --fil_out "${dir_out}/scaling.submit.se_missing_length.siq.tsv" \
            --idx_out 11 \
            --tbl_met "${tbl_met}" \
            --cfg_met "${cfg_met}" \
            --eqn 6nd \
            --err_out "${dir_err}" \
            --nam_job test_submit_calculate_scaling_factor_siq_se_missing_length
then
    record_fail \
        "submit_calculate_scaling_factor.sh siQ SE without length" \
        "unexpectedly passed"
else
    assert_pattern_found \
        "${fil_log}" \
        "supply '--len_def' or both '--len_mip' and '--len_min'" \
        "submit_calculate_scaling_factor.sh rejects SE siQ without length"
fi

finish
