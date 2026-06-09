#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_calculate_scaling_factor_spike.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="submit calculate-scaling-factor spike"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

#  Define fixture and output paths for the direct submit-worker smoke test
scr_sub="${ROOT_REPO}/scripts/submit_calculate_scaling_factor.sh"
dir_fix="${ROOT_REPO}/tests/calculate_scaling_factor/fixtures"
dir_bam_pe="${dir_fix}/bam/pe"
dir_cram_pe="${dir_fix}/cram/pe"
ref_fa="${dir_fix}/reference/tiny.fa"

bam_pe_mip="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sc.bam"
bam_pe_min="${dir_bam_pe}/in_WT_G1_Hho1_6336.sc.bam"
bam_pe_sip="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sp.bam"
bam_pe_sin="${dir_bam_pe}/in_WT_G1_Hho1_6336.sp.bam"

cram_pe_mip="${dir_cram_pe}/IP_WT_G1_Hho1_6336.sc.cram"
cram_pe_min="${dir_cram_pe}/in_WT_G1_Hho1_6336.sc.cram"
cram_pe_sip="${dir_cram_pe}/IP_WT_G1_Hho1_6336.sp.cram"
cram_pe_sin="${dir_cram_pe}/in_WT_G1_Hho1_6336.sp.cram"

tmp="${TEST_DIR_TMP}/submit_calculate_scaling_factor_spike"
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
        "${bam_pe_mip}" \
        "${bam_pe_mip}.bai" \
        "${bam_pe_min}" \
        "${bam_pe_min}.bai" \
        "${bam_pe_sip}" \
        "${bam_pe_sip}.bai" \
        "${bam_pe_sin}" \
        "${bam_pe_sin}.bai" \
        "${cram_pe_mip}" \
        "${cram_pe_min}" \
        "${cram_pe_sip}" \
        "${cram_pe_sin}" \
        "${ref_fa}"
then
    finish
    exit $?
fi


#  Direct submit-worker execution should write one indexed, data-only part row
fil_out="${dir_out}/scaling.submit.spike.tsv"
fil_prt="${fil_out}.part.000007"
fil_log="${dir_log}/submit_spike_pe_default.log"

printf -v row_exp '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${bam_pe_mip}" \
    "${bam_pe_sip}" \
    "${bam_pe_min}" \
    "${bam_pe_sin}" \
    2 \
    chiprx_alpha_ratio \
    3 \
    1 \
    2 \
    2

# shellcheck disable=SC2154
if \
    run_capture \
        "submit calculate-scaling-factor spike PE default" \
        "${fil_log}" \
        "${TEST_BASH}" "${scr_sub}" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --mode spike \
            --csv_mip "${bam_pe_mip}" \
            --csv_min "${bam_pe_min}" \
            --csv_sip "${bam_pe_sip}" \
            --csv_sin "${bam_pe_sin}" \
            --aln_typ auto \
            --fil_out "${fil_out}" \
            --idx_out 7 \
            --err_out "${dir_err}" \
            --nam_job test_submit_calculate_scaling_factor_spike_pe_default
then
    record_pass "submit_calculate_scaling_factor.sh spike PE exits 0"
else
    record_fail \
        "submit_calculate_scaling_factor.sh spike PE failed; see" \
        "$(print_relpath "${fil_log}")"
fi

assert_file_nonempty \
    "${fil_prt}" \
    "submit scaling-factor spike PE indexed part"

assert_pattern_absent \
    "${fil_prt}" \
    $'^main_ip\tspike_ip\tmain_in\tspike_in\tspike\tcoef' \
    "submit scaling-factor spike PE part stays data-only"

if [[ ! -e "${fil_out}" ]]; then
    record_pass "submit_calculate_scaling_factor.sh writes no final TSV"
else
    record_fail "submit_calculate_scaling_factor.sh unexpectedly wrote final TSV"
fi

if [[ "$(wc -l < "${fil_prt}")" -eq 1 ]]; then
    record_pass "submit scaling-factor spike PE part has one row"
else
    record_fail "submit scaling-factor spike PE part has unexpected row count"
fi

if [[ "$(cat "${fil_prt}")" == "${row_exp}" ]]; then
    record_pass "submit scaling-factor spike PE part has expected row"
else
    record_fail \
        "submit scaling-factor spike PE part has unexpected row; see" \
        "$(print_relpath "${fil_prt}")"
fi

assert_pattern_found \
    "${fil_log}" \
    'method=chiprx_alpha_ratio' \
    "submit scaling-factor spike PE resolves default method"

assert_pattern_found \
    "${fil_log}" \
    'idx_out=7' \
    "submit scaling-factor spike PE uses idx_out=7"

assert_pattern_found \
    "${fil_log}" \
    'typ_mp=pe' \
    "submit scaling-factor spike PE auto-detects main IP as PE"

assert_pattern_found \
    "${fil_log}" \
    'num_mp=3' \
    "submit scaling-factor spike PE counts main IP fragments"

assert_pattern_found \
    "${fil_log}" \
    'num_sp=1' \
    "submit scaling-factor spike PE counts spike IP fragments"

assert_pattern_found \
    "${fil_log}" \
    'num_mn=2' \
    "submit scaling-factor spike PE counts main input fragments"

assert_pattern_found \
    "${fil_log}" \
    'num_sn=2' \
    "submit scaling-factor spike PE counts spike input fragments"


#  Direct PE CRAM input should work with a reference and non-default method
fil_out="${dir_out}/scaling.submit.cram.spike.tsv"
fil_prt="${fil_out}.part.000008"
fil_log="${dir_log}/submit_spike_pe_cram_alpha_ip.log"

printf -v row_exp '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${cram_pe_mip}" \
    "${cram_pe_sip}" \
    "${cram_pe_min}" \
    "${cram_pe_sin}" \
    1000000 \
    chiprx_alpha_ip \
    3 \
    1 \
    2 \
    2

if \
    run_capture \
        "submit calculate-scaling-factor spike PE CRAM alpha IP" \
        "${fil_log}" \
        "${TEST_BASH}" "${scr_sub}" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --mode spike \
            --method chiprx_ip \
            --csv_mip "${cram_pe_mip}" \
            --csv_min "${cram_pe_min}" \
            --csv_sip "${cram_pe_sip}" \
            --csv_sin "${cram_pe_sin}" \
            --aln_typ auto \
            --ref_fa "${ref_fa}" \
            --fil_out "${fil_out}" \
            --idx_out 8 \
            --err_out "${dir_err}" \
            --nam_job test_submit_calculate_scaling_factor_spike_pe_cram
then
    record_pass "submit_calculate_scaling_factor.sh spike PE CRAM exits 0"
else
    record_fail \
        "submit_calculate_scaling_factor.sh spike PE CRAM failed; see" \
        "$(print_relpath "${fil_log}")"
fi

assert_file_nonempty \
    "${fil_prt}" \
    "submit scaling-factor spike PE CRAM indexed part"

assert_pattern_absent \
    "${fil_prt}" \
    $'^main_ip\tspike_ip\tmain_in\tspike_in\tspike\tcoef' \
    "submit scaling-factor spike PE CRAM part stays data-only"

if [[ "$(cat "${fil_prt}")" == "${row_exp}" ]]; then
    record_pass "submit scaling-factor spike PE CRAM part has expected row"
else
    record_fail \
        "submit scaling-factor spike PE CRAM part has unexpected row; see" \
        "$(print_relpath "${fil_prt}")"
fi

assert_pattern_found \
    "${fil_log}" \
    'method=chiprx_alpha_ip' \
    "submit scaling-factor spike PE CRAM canonicalizes chiprx_ip alias"

assert_pattern_found \
    "${fil_log}" \
    "ref_fa=${ref_fa}" \
    "submit scaling-factor spike PE CRAM records ref_fa"

assert_pattern_found \
    "${fil_log}" \
    'idx_out=8' \
    "submit scaling-factor spike PE CRAM uses idx_out=8"

assert_pattern_found \
    "${fil_log}" \
    'typ_mp=pe' \
    "submit scaling-factor spike PE CRAM auto-detects main IP as PE"

assert_pattern_found \
    "${fil_log}" \
    'num_mp=3' \
    "submit scaling-factor spike PE CRAM counts main IP fragments"

assert_pattern_found \
    "${fil_log}" \
    'num_sp=1' \
    "submit scaling-factor spike PE CRAM counts spike IP fragments"

assert_pattern_found \
    "${fil_log}" \
    'num_mn=2' \
    "submit scaling-factor spike PE CRAM counts main input fragments"

assert_pattern_found \
    "${fil_log}" \
    'num_sn=2' \
    "submit scaling-factor spike PE CRAM counts spike input fragments"


#  CRAM input without a reference FASTA should fail before worker execution
fil_log="${dir_log}/submit_spike_cram_missing_ref.log"
if \
    run_capture \
        "submit calculate-scaling-factor spike CRAM missing reference" \
        "${fil_log}" \
        "${TEST_BASH}" "${scr_sub}" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --mode spike \
            --csv_mip "${cram_pe_mip}" \
            --csv_min "${cram_pe_min}" \
            --csv_sip "${cram_pe_sip}" \
            --csv_sin "${cram_pe_sin}" \
            --aln_typ auto \
            --fil_out "${dir_out}/scaling.cram_missing_ref.spike.tsv" \
            --idx_out 3 \
            --err_out "${dir_err}" \
            --nam_job test_submit_calculate_scaling_factor_spike_cram_missing_ref
then
    record_fail \
        "submit_calculate_scaling_factor.sh CRAM without ref unexpectedly pass"
else
    assert_pattern_found \
        "${fil_log}" \
        "'--ref_fa' is required" \
        "submit_calculate_scaling_factor.sh rejects CRAM without ref_fa"
fi

finish
