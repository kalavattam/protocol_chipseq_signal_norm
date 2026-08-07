#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_calculate_scaling_factor_spike.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="submit calculate-scaling-factor spike"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


#  Run one direct submit-worker spike case and assert its part row
function run_submit_spike_case() {
    local cas="${1:-}"
    local idx_out="${2:-}"
    local mip="${3:-}"
    local min="${4:-}"
    local sip="${5:-}"
    local sin="${6:-}"
    local method="${7:-}"
    local row_exp="${8:-}"
    local typ_exp="${9:-}"
    shift 9 || true

    local fil_out="${dir_out}/scaling.submit.${cas}.spike.tsv"
    local fil_log="${dir_log}/submit_spike_${cas}.log"
    local nam_job="test_submit_calculate_scaling_factor_spike_"
    local -a arr_cmd=()

    nam_job+="${cas}"

    arr_cmd=(
        "${TEST_BASH}" "${scr_sub}"
            --env_nam "${env_nam}"
            --dir_scr "${ROOT_REPO}/bin"
            --threads 1
            --mode spike
            --aln_typ auto
            --csv_mip "${mip}"
            --csv_min "${min}"
            --csv_sip "${sip}"
            --csv_sin "${sin}"
            --fil_out "${fil_out}"
            --idx_out "${idx_out}"
            --dir_eo "${dir_err}"
            --nam_job "${nam_job}"
    )

    if [[ -n "${method}" ]]; then
        arr_cmd+=( --method "${method}" )
    fi

    arr_cmd+=( "$@" )

    run_case_scaling_factor_submit_part \
        "${cas}" \
        spike \
        arr_cmd \
        "${fil_out}" \
        "${idx_out}" \
        "${fil_log}" \
        "${row_exp}" \
        $'^main_ip\tspike_ip\tmain_in\tspike_in\tspike\tcoef'

    assert_pattern_found \
        "${fil_log}" \
        "typ_mp=${typ_exp}" \
        "submit scaling-factor spike ${cas} auto-detects main IP as ${typ_exp}"
}


#  Define fixture and output paths for the direct submit-worker test
scr_sub="${ROOT_REPO}/bin/submit_calculate_scaling_factor.sh"
dir_fix="${ROOT_REPO}/tests/fixtures/calculate_scaling_factor"
dir_bam_se="${dir_fix}/bam/se"
dir_bam_pe="${dir_fix}/bam/pe"
dir_cram_pe="${dir_fix}/cram/pe"
ref_fa="${dir_fix}/reference/tiny.fa"

bam_se_mip="${dir_bam_se}/IP_WT_log_Brn1_rep1.sc.bam"
bam_se_min="${dir_bam_se}/in_WT_log_Brn1_rep1.sc.bam"
bam_se_sip="${dir_bam_se}/IP_WT_log_Brn1_rep1.sp.bam"
bam_se_sin="${dir_bam_se}/in_WT_log_Brn1_rep1.sp.bam"

bam_pe_mip="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sc.bam"
bam_pe_min="${dir_bam_pe}/in_WT_G1_Hho1_6336.sc.bam"
bam_pe_sip="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sp.bam"
bam_pe_sin="${dir_bam_pe}/in_WT_G1_Hho1_6336.sp.bam"

cram_pe_mip="${dir_cram_pe}/IP_WT_G1_Hho1_6336.sc.cram"
cram_pe_min="${dir_cram_pe}/in_WT_G1_Hho1_6336.sc.cram"
cram_pe_sip="${dir_cram_pe}/IP_WT_G1_Hho1_6336.sp.cram"
cram_pe_sin="${dir_cram_pe}/in_WT_G1_Hho1_6336.sp.cram"

crai_pe_mip="${cram_pe_mip}.crai"
crai_pe_min="${cram_pe_min}.crai"
crai_pe_sip="${cram_pe_sip}.crai"
crai_pe_sin="${cram_pe_sin}.crai"

bai_se_mip="${bam_se_mip}.bai"
bai_se_min="${bam_se_min}.bai"
bai_se_sip="${bam_se_sip}.bai"
bai_se_sin="${bam_se_sin}.bai"
bai_pe_mip="${bam_pe_mip}.bai"
bai_pe_min="${bam_pe_min}.bai"
bai_pe_sip="${bam_pe_sip}.bai"
bai_pe_sin="${bam_pe_sin}.bai"

tmp="${TEST_DIR_TMP}/submit_calculate_scaling_factor_spike"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/calculate_scaling_factor"

fil_out_pe_default="${dir_out}/scaling.submit.spike.tsv"
fil_out_pe_cram_alpha_ip="${dir_out}/scaling.submit.cram.spike.tsv"
fil_out_cram_missing_ref="${dir_out}/scaling.cram_missing_ref.spike.tsv"
fil_out_invalid_method="${dir_out}/scaling.invalid_method.spike.tsv"

fil_prt_pe_default="${fil_out_pe_default}.part.000007"
fil_prt_pe_cram_alpha_ip="${fil_out_pe_cram_alpha_ip}.part.000008"

fil_log_pe_default="${dir_log}/submit_spike_pe_default.log"
fil_log_pe_cram_alpha_ip="${dir_log}/submit_spike_pe_cram_alpha_ip.log"
fil_log_cram_missing_ref="${dir_log}/submit_spike_cram_missing_ref.log"
fil_log_invalid_method="${dir_log}/submit_spike_invalid_method.log"

nam_job_pe_default="test_submit_calculate_scaling_factor_spike_pe_default"
nam_job_pe_cram_alpha_ip="test_submit_calculate_scaling_factor_spike_pe_cram"
nam_job_cram_missing_ref="test_submit_calculate_scaling_factor_spike_cram_missing_ref"
nam_job_invalid_method="test_submit_calculate_scaling_factor_spike_invalid_method"

printf -v row_exp_pe_default '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
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

printf -v row_exp_pe_cram_alpha_ip \
    '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
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

printf -v row_exp_se_default '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${bam_se_mip}" \
    "${bam_se_sip}" \
    "${bam_se_min}" \
    "${bam_se_sin}" \
    2 \
    chiprx_alpha_ratio \
    3 \
    1 \
    2 \
    2

printf -v row_exp_se_fractional '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${bam_se_mip}" \
    "${bam_se_sip}" \
    "${bam_se_min}" \
    "${bam_se_sin}" \
    2 \
    fractional \
    3 \
    1 \
    2 \
    2

printf -v row_exp_se_alpha_ip '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${bam_se_mip}" \
    "${bam_se_sip}" \
    "${bam_se_min}" \
    "${bam_se_sin}" \
    1000000 \
    chiprx_alpha_ip \
    3 \
    1 \
    2 \
    2

printf -v row_exp_se_alpha_in '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${bam_se_mip}" \
    "${bam_se_sip}" \
    "${bam_se_min}" \
    "${bam_se_sin}" \
    500000 \
    chiprx_alpha_in \
    3 \
    1 \
    2 \
    2

printf -v row_exp_se_rxinput '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${bam_se_mip}" \
    "${bam_se_sip}" \
    "${bam_se_min}" \
    "${bam_se_sin}" \
    500000 \
    rxinput_alpha \
    3 \
    1 \
    2 \
    2

printf -v row_exp_se_depth_override \
    '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${bam_se_mip}" \
    "${bam_se_sip}" \
    "${bam_se_min}" \
    "${bam_se_sin}" \
    1 \
    fractional \
    10 \
    2 \
    20 \
    4

print_section "${TEST_NAME}"


rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

if ! \
    require_files_nonempty \
        "${bam_se_mip}" \
        "${bai_se_mip}" \
        "${bam_se_min}" \
        "${bai_se_min}" \
        "${bam_se_sip}" \
        "${bai_se_sip}" \
        "${bam_se_sin}" \
        "${bai_se_sin}" \
        "${bam_pe_mip}" \
        "${bai_pe_mip}" \
        "${bam_pe_min}" \
        "${bai_pe_min}" \
        "${bam_pe_sip}" \
        "${bai_pe_sip}" \
        "${bam_pe_sin}" \
        "${bai_pe_sin}" \
        "${cram_pe_mip}" \
        "${crai_pe_mip}" \
        "${cram_pe_min}" \
        "${crai_pe_min}" \
        "${cram_pe_sip}" \
        "${crai_pe_sip}" \
        "${cram_pe_sin}" \
        "${crai_pe_sin}" \
        "${ref_fa}"
then
    finish
    exit $?
fi

#  Direct submit-worker execution should write one indexed, data-only part row
# shellcheck disable=SC2154
if \
    run_capture \
        "submit calculate-scaling-factor spike PE default" \
        "${fil_log_pe_default}" \
        "${TEST_BASH}" "${scr_sub}" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --mode spike \
            --aln_typ auto \
            --csv_mip "${bam_pe_mip}" \
            --csv_min "${bam_pe_min}" \
            --csv_sip "${bam_pe_sip}" \
            --csv_sin "${bam_pe_sin}" \
            --fil_out "${fil_out_pe_default}" \
            --idx_out 7 \
            --dir_eo "${dir_err}" \
            --nam_job "${nam_job_pe_default}"
then
    record_pass "submit_calculate_scaling_factor.sh spike PE exits 0"
else
    record_fail \
        "submit_calculate_scaling_factor.sh spike PE failed; see" \
        "$(print_relpath "${fil_log_pe_default}")"
fi

assert_file_nonempty \
    "${fil_prt_pe_default}" \
    "submit scaling-factor spike PE indexed part"

assert_pattern_absent \
    "${fil_prt_pe_default}" \
    $'^main_ip\tspike_ip\tmain_in\tspike_in\tspike\tcoef' \
    "submit scaling-factor spike PE part stays data-only"

if [[ ! -e "${fil_out_pe_default}" ]]; then
    record_pass "submit_calculate_scaling_factor.sh writes no final TSV"
else
    record_fail "submit_calculate_scaling_factor.sh unexpectedly wrote final TSV"
fi

if [[ "$(wc -l < "${fil_prt_pe_default}")" -eq 1 ]]; then
    record_pass "submit scaling-factor spike PE part has one row"
else
    record_fail "submit scaling-factor spike PE part has unexpected row count"
fi

if [[ "$(cat "${fil_prt_pe_default}")" == "${row_exp_pe_default}" ]]; then
    record_pass "submit scaling-factor spike PE part has expected row"
else
    record_fail \
        "submit scaling-factor spike PE part has unexpected row; see" \
        "$(print_relpath "${fil_prt_pe_default}")"
fi

assert_pattern_found \
    "${fil_log_pe_default}" \
    'method=chiprx_alpha_ratio' \
    "submit scaling-factor spike PE resolves default method"

assert_pattern_found \
    "${fil_log_pe_default}" \
    'idx_out=7' \
    "submit scaling-factor spike PE uses idx_out=7"

assert_pattern_found \
    "${fil_log_pe_default}" \
    'typ_mp=pe' \
    "submit scaling-factor spike PE auto-detects main IP as PE"

assert_pattern_found \
    "${fil_log_pe_default}" \
    'num_mp=3' \
    "submit scaling-factor spike PE counts main IP fragments"

assert_pattern_found \
    "${fil_log_pe_default}" \
    'num_sp=1' \
    "submit scaling-factor spike PE counts spike IP fragments"

assert_pattern_found \
    "${fil_log_pe_default}" \
    'num_mn=2' \
    "submit scaling-factor spike PE counts main input fragments"

assert_pattern_found \
    "${fil_log_pe_default}" \
    'num_sn=2' \
    "submit scaling-factor spike PE counts spike input fragments"


#  Direct PE CRAM input should work with a reference and non-default method
if \
    run_capture \
        "submit calculate-scaling-factor spike PE CRAM alpha IP" \
        "${fil_log_pe_cram_alpha_ip}" \
        "${TEST_BASH}" "${scr_sub}" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --mode spike \
            --method chiprx_ip \
            --aln_typ auto \
            --ref_fa "${ref_fa}" \
            --csv_mip "${cram_pe_mip}" \
            --csv_min "${cram_pe_min}" \
            --csv_sip "${cram_pe_sip}" \
            --csv_sin "${cram_pe_sin}" \
            --fil_out "${fil_out_pe_cram_alpha_ip}" \
            --idx_out 8 \
            --dir_eo "${dir_err}" \
            --nam_job "${nam_job_pe_cram_alpha_ip}"
then
    record_pass "submit_calculate_scaling_factor.sh spike PE CRAM exits 0"
else
    record_fail \
        "submit_calculate_scaling_factor.sh spike PE CRAM failed; see" \
        "$(print_relpath "${fil_log_pe_cram_alpha_ip}")"
fi

assert_file_nonempty \
    "${fil_prt_pe_cram_alpha_ip}" \
    "submit scaling-factor spike PE CRAM indexed part"

assert_pattern_absent \
    "${fil_prt_pe_cram_alpha_ip}" \
    $'^main_ip\tspike_ip\tmain_in\tspike_in\tspike\tcoef' \
    "submit scaling-factor spike PE CRAM part stays data-only"

if [[ "$(cat "${fil_prt_pe_cram_alpha_ip}")" == \
    "${row_exp_pe_cram_alpha_ip}" ]]
then
    record_pass "submit scaling-factor spike PE CRAM part has expected row"
else
    record_fail \
        "submit scaling-factor spike PE CRAM part has unexpected row; see" \
        "$(print_relpath "${fil_prt_pe_cram_alpha_ip}")"
fi

assert_pattern_found \
    "${fil_log_pe_cram_alpha_ip}" \
    'method=chiprx_alpha_ip' \
    "submit scaling-factor spike PE CRAM canonicalizes chiprx_ip alias"

assert_pattern_found \
    "${fil_log_pe_cram_alpha_ip}" \
    "ref_fa=${ref_fa}" \
    "submit scaling-factor spike PE CRAM records ref_fa"

assert_pattern_found \
    "${fil_log_pe_cram_alpha_ip}" \
    'idx_out=8' \
    "submit scaling-factor spike PE CRAM uses idx_out=8"

assert_pattern_found \
    "${fil_log_pe_cram_alpha_ip}" \
    'typ_mp=pe' \
    "submit scaling-factor spike PE CRAM auto-detects main IP as PE"

assert_pattern_found \
    "${fil_log_pe_cram_alpha_ip}" \
    'num_mp=3' \
    "submit scaling-factor spike PE CRAM counts main IP fragments"

assert_pattern_found \
    "${fil_log_pe_cram_alpha_ip}" \
    'num_sp=1' \
    "submit scaling-factor spike PE CRAM counts spike IP fragments"

assert_pattern_found \
    "${fil_log_pe_cram_alpha_ip}" \
    'num_mn=2' \
    "submit scaling-factor spike PE CRAM counts main input fragments"

assert_pattern_found \
    "${fil_log_pe_cram_alpha_ip}" \
    'num_sn=2' \
    "submit scaling-factor spike PE CRAM counts spike input fragments"


#  CRAM input without a reference FASTA should fail before worker execution
if \
    run_capture \
        "submit calculate-scaling-factor spike CRAM missing reference" \
        "${fil_log_cram_missing_ref}" \
        "${TEST_BASH}" "${scr_sub}" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --mode spike \
            --aln_typ auto \
            --csv_mip "${cram_pe_mip}" \
            --csv_min "${cram_pe_min}" \
            --csv_sip "${cram_pe_sip}" \
            --csv_sin "${cram_pe_sin}" \
            --fil_out "${fil_out_cram_missing_ref}" \
            --idx_out 3 \
            --dir_eo "${dir_err}" \
            --nam_job "${nam_job_cram_missing_ref}"
then
    record_fail \
        "submit_calculate_scaling_factor.sh CRAM without ref unexpectedly pass"
else
    assert_pattern_found \
        "${fil_log_cram_missing_ref}" \
        "'--ref_fa' is required" \
        "submit_calculate_scaling_factor.sh rejects CRAM without ref_fa"
fi


#  Unknown spike-in methods should fail before worker execution
if \
    run_capture \
        "submit calculate-scaling-factor spike invalid method" \
        "${fil_log_invalid_method}" \
        "${TEST_BASH}" "${scr_sub}" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --mode spike \
            --method not_a_method \
            --aln_typ auto \
            --csv_mip "${bam_pe_mip}" \
            --csv_min "${bam_pe_min}" \
            --csv_sip "${bam_pe_sip}" \
            --csv_sin "${bam_pe_sin}" \
            --fil_out "${fil_out_invalid_method}" \
            --idx_out 4 \
            --dir_eo "${dir_err}" \
            --nam_job "${nam_job_invalid_method}"
then
    record_fail \
        "submit_calculate_scaling_factor.sh spike invalid method" \
        "unexpectedly passed"
else
    assert_pattern_found \
        "${fil_log_invalid_method}" \
        "invalid '--method' value" \
        "submit_calculate_scaling_factor.sh rejects invalid spike method"
fi


#  SE BAM input should be accepted and auto-detected by the submit wrapper
run_submit_spike_case \
    se_default \
    9 \
    "${bam_se_mip}" \
    "${bam_se_min}" \
    "${bam_se_sip}" \
    "${bam_se_sin}" \
    "" \
    "${row_exp_se_default}" \
    se


#  Canonical spike-in methods should produce expected submit-worker rows
run_submit_spike_case \
    se_fractional \
    10 \
    "${bam_se_mip}" \
    "${bam_se_min}" \
    "${bam_se_sip}" \
    "${bam_se_sin}" \
    fractional \
    "${row_exp_se_fractional}" \
    se

run_submit_spike_case \
    se_alpha_ip \
    11 \
    "${bam_se_mip}" \
    "${bam_se_min}" \
    "${bam_se_sip}" \
    "${bam_se_sin}" \
    chiprx_alpha_ip \
    "${row_exp_se_alpha_ip}" \
    se

run_submit_spike_case \
    se_alpha_in \
    12 \
    "${bam_se_mip}" \
    "${bam_se_min}" \
    "${bam_se_sip}" \
    "${bam_se_sin}" \
    chiprx_alpha_in \
    "${row_exp_se_alpha_in}" \
    se

run_submit_spike_case \
    se_rxinput \
    13 \
    "${bam_se_mip}" \
    "${bam_se_min}" \
    "${bam_se_sip}" \
    "${bam_se_sin}" \
    rxinput_alpha \
    "${row_exp_se_rxinput}" \
    se


#  Explicit depth overrides should replace alignment-derived counts
run_submit_spike_case \
    se_depth_override \
    14 \
    "${bam_se_mip}" \
    "${bam_se_min}" \
    "${bam_se_sip}" \
    "${bam_se_sin}" \
    fractional \
    "${row_exp_se_depth_override}" \
    se \
    --csv_dep_mip 10 \
    --csv_dep_min 20 \
    --csv_dep_sip 2 \
    --csv_dep_sin 4

finish
