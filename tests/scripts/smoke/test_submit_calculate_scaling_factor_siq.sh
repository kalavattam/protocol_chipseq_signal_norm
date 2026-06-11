#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_calculate_scaling_factor_siq.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development.
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
dir_cram="${dir_fix}/cram"
dir_bam_se="${dir_bam}/se"
dir_bam_pe="${dir_bam}/pe"
dir_cram_pe="${dir_cram}/pe"
dir_met="${dir_fix}/metadata"

tbl_met="${dir_met}/measurements_siqchip.tsv"
tbl_gz="${dir_met}/measurements_siqchip.tsv.gz"
tbl_als="${dir_met}/measurements_siqchip_aliases.tsv"
tbl_pfx="${dir_met}/measurements_siqchip_skip_prefixes.tsv"
tbl_col="${dir_met}/measurements_siqchip_alias_collision.tsv"
tbl_dup="${dir_met}/measurements_siqchip_duplicate_match.tsv"
tbl_trt="${dir_met}/measurements_siqchip_treatment.tsv"
cfg_rgx="${dir_fix}/config/parse_metadata_siqchip_regex.yml"
cfg_trt="${dir_fix}/config/parse_metadata_siqchip_match_treatment.yml"
ref_fa="${dir_fix}/reference/tiny.fa"

bam_pe_mip="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sc.bam"
bam_pe_min="${dir_bam_pe}/in_WT_G1_Hho1_6336.sc.bam"
hu_pe_mip="${dir_bam_pe}/IP_WT_G1_HU_Hho1_6336.sc.bam"
hu_pe_min="${dir_bam_pe}/in_WT_G1_HU_Hho1_6336.sc.bam"
cram_pe_mip="${dir_cram_pe}/IP_WT_G1_Hho1_6336.sc.cram"
cram_pe_min="${dir_cram_pe}/in_WT_G1_Hho1_6336.sc.cram"
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
        "${cfg_rgx}" \
        "${cfg_trt}" \
        "${tbl_met}" \
        "${tbl_gz}" \
        "${tbl_als}" \
        "${tbl_pfx}" \
        "${tbl_col}" \
        "${tbl_dup}" \
        "${tbl_trt}" \
        "${bam_pe_mip}" \
        "${bam_pe_mip}.bai" \
        "${bam_pe_min}" \
        "${bam_pe_min}.bai" \
        "${hu_pe_mip}" \
        "${hu_pe_mip}.bai" \
        "${hu_pe_min}" \
        "${hu_pe_min}.bai" \
        "${cram_pe_mip}" \
        "${cram_pe_min}" \
        "${ref_fa}" \
        "${bam_se_mip}" \
        "${bam_se_mip}.bai" \
        "${bam_se_min}" \
        "${bam_se_min}.bai"
then
    finish
    exit $?
fi


#  Run one PE siQ submit-worker case and assert the exact part row
function run_submit_pe_siq() {
    local cas="${1:-}"
    local mip="${2:-}"
    local min="${3:-}"
    local tbl="${4:-}"
    local cfg="${5:-}"
    local idx="${6:-}"
    local alpha="${7:-}"
    local mass_ip="${8:-}"
    local mass_in="${9:-}"
    shift 9 || true

    local fil_out="${dir_out}/scaling.submit.${cas}.siq.tsv"
    local fil_log="${dir_log}/submit_siq_${cas}.log"
    local row_exp
    local -a arr_cmd=(
        "${TEST_BASH}" "${scr_sub}"
            --env_nam "${env_nam}"
            --dir_scr "${ROOT_REPO}/scripts"
            --threads 1
            --mode siq
            --csv_mip "${mip}"
            --csv_min "${min}"
            --aln_typ auto
            --fil_out "${fil_out}"
            --idx_out "${idx}"
            --tbl_met "${tbl}"
            --cfg_met "${cfg}"
            --eqn 6nd
            --err_out "${dir_err}"
            --nam_job "test_submit_calculate_scaling_factor_siq_${cas}"
    )

    printf -v row_exp '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
        "${mip}" \
        "${min}" \
        "${alpha}" \
        '6nd' \
        "${mass_ip}" \
        "${mass_in}" \
        '300' \
        '20' \
        '3' \
        '2' \
        '20' \
        '20'

    run_case_scaling_factor_submit_part \
        "${cas}" \
        siQ \
        arr_cmd \
        "${fil_out}" \
        "${idx}" \
        "${fil_log}" \
        "${row_exp}" \
        $'^fil_ip\tfil_in\tsiq\teqn\tmass_ip\tmass_in\tvol_all\tvol_in\tdep_ip\tdep_in\tlen_ip\tlen_in$' \
        "$@"
}


#  Run one PE siQ submit-worker case that should fail during metadata parsing
function run_submit_pe_siq_failure() {
    local cas="${1:-}"
    local tbl="${2:-}"
    local patn="${3:-}"
    local idx="${4:-}"
    local fil_log="${dir_log}/submit_siq_${cas}.log"

    if \
        run_capture \
            "submit calculate-scaling-factor siQ ${cas}" \
            "${fil_log}" \
            "${TEST_BASH}" "${scr_sub}" \
                --env_nam "${env_nam}" \
                --dir_scr "${ROOT_REPO}/scripts" \
                --threads 1 \
                --mode siq \
                --csv_mip "${bam_pe_mip}" \
                --csv_min "${bam_pe_min}" \
                --aln_typ auto \
                --fil_out "${dir_out}/scaling.submit.${cas}.siq.tsv" \
                --idx_out "${idx}" \
                --tbl_met "${tbl}" \
                --cfg_met "${cfg_met}" \
                --eqn 6nd \
                --err_out "${dir_err}" \
                --nam_job "test_submit_calculate_scaling_factor_siq_${cas}"
    then
        record_fail \
            "submit_calculate_scaling_factor.sh siQ ${cas} unexpectedly" \
            "passed"
    else
        assert_pattern_found \
            "${fil_log}" \
            "${patn}" \
            "submit_calculate_scaling_factor.sh siQ ${cas} fails clearly"
    fi
}


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


#  PE CRAM input should compute when an explicit reference is supplied
fil_out="${dir_out}/scaling.submit.pe_cram.siq.tsv"
fil_prt="${fil_out}.part.000013"
fil_log="${dir_log}/submit_siq_pe_cram.log"

printf -v row_exp '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${cram_pe_mip}" \
    "${cram_pe_min}" \
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

if \
    run_capture \
        "submit calculate-scaling-factor siQ PE CRAM" \
        "${fil_log}" \
        "${TEST_BASH}" "${scr_sub}" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --mode siq \
            --csv_mip "${cram_pe_mip}" \
            --csv_min "${cram_pe_min}" \
            --aln_typ auto \
            --ref_fa "${ref_fa}" \
            --fil_out "${fil_out}" \
            --idx_out 13 \
            --tbl_met "${tbl_met}" \
            --cfg_met "${cfg_met}" \
            --eqn 6nd \
            --err_out "${dir_err}" \
            --nam_job test_submit_calculate_scaling_factor_siq_pe_cram
then
    record_pass "submit_calculate_scaling_factor.sh siQ PE CRAM exits 0"
else
    record_fail \
        "submit_calculate_scaling_factor.sh siQ PE CRAM failed; see" \
        "$(print_relpath "${fil_log}")"
fi

assert_file_nonempty \
    "${fil_prt}" \
    "submit scaling-factor siQ PE CRAM indexed part"

if [[ "$(cat "${fil_prt}")" == "${row_exp}" ]]; then
    record_pass "submit scaling-factor siQ PE CRAM part has expected row"
else
    record_fail \
        "submit scaling-factor siQ PE CRAM part has unexpected row; see" \
        "$(print_relpath "${fil_prt}")"
fi

assert_pattern_found \
    "${fil_log}" \
    "ref_fa=${ref_fa}" \
    "submit scaling-factor siQ PE CRAM records ref_fa"

assert_pattern_found \
    "${fil_log}" \
    'typ_ip=pe' \
    "submit scaling-factor siQ PE CRAM auto-detects IP as PE"


#  CRAM input should fail clearly when no reference is supplied
fil_log="${dir_log}/submit_siq_cram_missing_ref.log"
if \
    run_capture \
        "submit calculate-scaling-factor siQ CRAM missing reference" \
        "${fil_log}" \
        "${TEST_BASH}" "${scr_sub}" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --mode siq \
            --csv_mip "${cram_pe_mip}" \
            --csv_min "${cram_pe_min}" \
            --aln_typ auto \
            --fil_out "${dir_out}/scaling.submit.cram_missing_ref.siq.tsv" \
            --idx_out 14 \
            --tbl_met "${tbl_met}" \
            --cfg_met "${cfg_met}" \
            --eqn 6nd \
            --err_out "${dir_err}" \
            --nam_job test_submit_calculate_scaling_factor_siq_cram_missing_ref
then
    record_fail \
        "submit_calculate_scaling_factor.sh siQ CRAM without reference" \
        "unexpectedly passed"
else
    assert_pattern_found \
        "${fil_log}" \
        "'--ref_fa' is required" \
        "submit_calculate_scaling_factor.sh rejects CRAM without ref_fa"
fi


#  PE BAM input should also work with accepted metadata column aliases
fil_out="${dir_out}/scaling.submit.pe_aliases.siq.tsv"
fil_prt="${fil_out}.part.000012"
fil_log="${dir_log}/submit_siq_pe_aliases.log"

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

if \
    run_capture \
        "submit calculate-scaling-factor siQ PE metadata aliases" \
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
            --idx_out 12 \
            --tbl_met "${tbl_als}" \
            --cfg_met "${cfg_met}" \
            --eqn 6nd \
            --err_out "${dir_err}" \
            --nam_job test_submit_calculate_scaling_factor_siq_pe_aliases
then
    record_pass "submit_calculate_scaling_factor.sh siQ PE aliases exits 0"
else
    record_fail \
        "submit_calculate_scaling_factor.sh siQ PE aliases failed; see" \
        "$(print_relpath "${fil_log}")"
fi

assert_file_nonempty \
    "${fil_prt}" \
    "submit scaling-factor siQ PE aliases indexed part"

if [[ "$(cat "${fil_prt}")" == "${row_exp}" ]]; then
    record_pass "submit scaling-factor siQ PE aliases part has expected row"
else
    record_fail \
        "submit scaling-factor siQ PE aliases part has unexpected row; see" \
        "$(print_relpath "${fil_prt}")"
fi

assert_pattern_found \
    "${fil_log}" \
    'mass_ip=2.7' \
    "submit scaling-factor siQ PE aliases propagates metadata aliases"


#  Metadata parser variants should propagate through the submit wrapper
run_submit_pe_siq \
    pe_gzip_metadata \
    "${bam_pe_mip}" \
    "${bam_pe_min}" \
    "${tbl_gz}" \
    "${cfg_met}" \
    19 \
    '0.002660098522167487697376' \
    '2.7' \
    '72.5'

run_submit_pe_siq \
    pe_skip_prefixes \
    "${bam_pe_mip}" \
    "${bam_pe_min}" \
    "${tbl_pfx}" \
    "${cfg_met}" \
    20 \
    '0.002660098522167487697376' \
    '2.7' \
    '72.5'

run_submit_pe_siq \
    pe_regex_config \
    "${bam_pe_mip}" \
    "${bam_pe_min}" \
    "${tbl_met}" \
    "${cfg_rgx}" \
    21 \
    '0.002660098522167487697376' \
    '2.7' \
    '72.5'

run_submit_pe_siq \
    pe_treatment \
    "${hu_pe_mip}" \
    "${hu_pe_min}" \
    "${tbl_trt}" \
    "${cfg_trt}" \
    22 \
    '0.007963320463320463019063' \
    '9.9' \
    '88.8'

run_submit_pe_siq_failure \
    pe_alias_collision \
    "${tbl_col}" \
    "multiple columns that normalize to the same canonical name" \
    23

run_submit_pe_siq_failure \
    pe_duplicate_match \
    "${tbl_dup}" \
    "Multiple matching rows found" \
    24


#  Depth-dependent equation 5 should use alignment-inferred depths
fil_out="${dir_out}/scaling.submit.pe_eqn5.siq.tsv"
fil_prt="${fil_out}.part.000015"
fil_log="${dir_log}/submit_siq_pe_eqn5.log"

printf -v row_exp '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${bam_pe_mip}" \
    "${bam_pe_min}" \
    '0.001655172413793103407958' \
    '5' \
    '2.7' \
    '72.5' \
    '300' \
    '20' \
    '3' \
    '2' \
    '20' \
    '20'

if \
    run_capture \
        "submit calculate-scaling-factor siQ PE equation 5" \
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
            --idx_out 15 \
            --tbl_met "${tbl_met}" \
            --cfg_met "${cfg_met}" \
            --eqn 5 \
            --err_out "${dir_err}" \
            --nam_job test_submit_calculate_scaling_factor_siq_pe_eqn5
then
    record_pass "submit_calculate_scaling_factor.sh siQ PE equation 5 exits 0"
else
    record_fail \
        "submit_calculate_scaling_factor.sh siQ PE equation 5 failed; see" \
        "$(print_relpath "${fil_log}")"
fi

assert_file_nonempty \
    "${fil_prt}" \
    "submit scaling-factor siQ PE equation 5 indexed part"

if [[ "$(cat "${fil_prt}")" == "${row_exp}" ]]; then
    record_pass "submit scaling-factor siQ PE equation 5 part has expected row"
else
    record_fail \
        "submit scaling-factor siQ PE equation 5 part has unexpected row; see" \
        "$(print_relpath "${fil_prt}")"
fi


#  Equation 5 without depth terms should ignore alignment-inferred depths
fil_out="${dir_out}/scaling.submit.pe_eqn5nd.siq.tsv"
fil_prt="${fil_out}.part.000025"
fil_log="${dir_log}/submit_siq_pe_eqn5nd.log"

printf -v row_exp '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${bam_pe_mip}" \
    "${bam_pe_min}" \
    '0.002482758620689655328778' \
    '5nd' \
    '2.7' \
    '72.5' \
    '300' \
    '20' \
    '3' \
    '2' \
    '20' \
    '20'

if \
    run_capture \
        "submit calculate-scaling-factor siQ PE equation 5nd" \
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
            --idx_out 25 \
            --tbl_met "${tbl_met}" \
            --cfg_met "${cfg_met}" \
            --eqn 5nd \
            --err_out "${dir_err}" \
            --nam_job test_submit_calculate_scaling_factor_siq_pe_eqn5nd
then
    record_pass "submit_calculate_scaling_factor.sh siQ PE equation 5nd exits 0"
else
    record_fail \
        "submit_calculate_scaling_factor.sh siQ PE equation 5nd failed; see" \
        "$(print_relpath "${fil_log}")"
fi

assert_file_nonempty \
    "${fil_prt}" \
    "submit scaling-factor siQ PE equation 5nd indexed part"

assert_file_exact_line \
    "${fil_prt}" \
    "${row_exp}" \
    "submit scaling-factor siQ PE equation 5nd part"


#  Depth-dependent equation 6 should use alignment-inferred depths
fil_out="${dir_out}/scaling.submit.pe_eqn6.siq.tsv"
fil_prt="${fil_out}.part.000016"
fil_log="${dir_log}/submit_siq_pe_eqn6.log"

printf -v row_exp '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${bam_pe_mip}" \
    "${bam_pe_min}" \
    '0.001773399014778325203864' \
    '6' \
    '2.7' \
    '72.5' \
    '300' \
    '20' \
    '3' \
    '2' \
    '20' \
    '20'

if \
    run_capture \
        "submit calculate-scaling-factor siQ PE equation 6" \
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
            --idx_out 16 \
            --tbl_met "${tbl_met}" \
            --cfg_met "${cfg_met}" \
            --eqn 6 \
            --err_out "${dir_err}" \
            --nam_job test_submit_calculate_scaling_factor_siq_pe_eqn6
then
    record_pass "submit_calculate_scaling_factor.sh siQ PE equation 6 exits 0"
else
    record_fail \
        "submit_calculate_scaling_factor.sh siQ PE equation 6 failed; see" \
        "$(print_relpath "${fil_log}")"
fi

assert_file_nonempty \
    "${fil_prt}" \
    "submit scaling-factor siQ PE equation 6 indexed part"

if [[ "$(cat "${fil_prt}")" == "${row_exp}" ]]; then
    record_pass "submit scaling-factor siQ PE equation 6 part has expected row"
else
    record_fail \
        "submit scaling-factor siQ PE equation 6 part has unexpected row; see" \
        "$(print_relpath "${fil_prt}")"
fi


#  Explicit length overrides should take precedence over PE fragment sizes
fil_out="${dir_out}/scaling.submit.pe_len_override.siq.tsv"
fil_prt="${fil_out}.part.000017"
fil_log="${dir_log}/submit_siq_pe_len_override.log"

printf -v row_exp '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${bam_pe_mip}" \
    "${bam_pe_min}" \
    '0.005320197044334976262114' \
    '6nd' \
    '2.7' \
    '72.5' \
    '300' \
    '20' \
    '3' \
    '2' \
    '100' \
    '200'

if \
    run_capture \
        "submit calculate-scaling-factor siQ PE length overrides" \
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
            --idx_out 17 \
            --tbl_met "${tbl_met}" \
            --cfg_met "${cfg_met}" \
            --eqn 6nd \
            --len_mip 100 \
            --len_min 200 \
            --err_out "${dir_err}" \
            --nam_job test_submit_calculate_scaling_factor_siq_pe_len_override
then
    record_pass "submit_calculate_scaling_factor.sh siQ PE length overrides exits 0"
else
    record_fail \
        "submit_calculate_scaling_factor.sh siQ PE length overrides failed;" \
        "see $(print_relpath "${fil_log}")"
fi

assert_file_nonempty \
    "${fil_prt}" \
    "submit scaling-factor siQ PE length overrides indexed part"

if [[ "$(cat "${fil_prt}")" == "${row_exp}" ]]; then
    record_pass \
        "submit scaling-factor siQ PE length overrides part has expected row"
else
    record_fail \
        "submit scaling-factor siQ PE length overrides part has unexpected" \
        "row; see $(print_relpath "${fil_prt}")"
fi


#  Explicit depth overrides should feed depth-dependent equations
fil_out="${dir_out}/scaling.submit.pe_dep_override.siq.tsv"
fil_prt="${fil_out}.part.000018"
fil_log="${dir_log}/submit_siq_pe_dep_override.log"

printf -v row_exp '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${bam_pe_mip}" \
    "${bam_pe_min}" \
    '0.005320197044334976262114' \
    '6' \
    '2.7' \
    '72.5' \
    '300' \
    '20' \
    '10' \
    '20' \
    '20' \
    '20'

if \
    run_capture \
        "submit calculate-scaling-factor siQ PE depth overrides" \
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
            --idx_out 18 \
            --tbl_met "${tbl_met}" \
            --cfg_met "${cfg_met}" \
            --eqn 6 \
            --dep_mip 10 \
            --dep_min 20 \
            --err_out "${dir_err}" \
            --nam_job test_submit_calculate_scaling_factor_siq_pe_dep_override
then
    record_pass "submit_calculate_scaling_factor.sh siQ PE depth overrides exits 0"
else
    record_fail \
        "submit_calculate_scaling_factor.sh siQ PE depth overrides failed;" \
        "see $(print_relpath "${fil_log}")"
fi

assert_file_nonempty \
    "${fil_prt}" \
    "submit scaling-factor siQ PE depth overrides indexed part"

if [[ "$(cat "${fil_prt}")" == "${row_exp}" ]]; then
    record_pass \
        "submit scaling-factor siQ PE depth overrides part has expected row"
else
    record_fail \
        "submit scaling-factor siQ PE depth overrides part has unexpected" \
        "row; see $(print_relpath "${fil_prt}")"
fi


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
