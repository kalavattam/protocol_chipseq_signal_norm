#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_calculate_scaling_factor_siq_yaml.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Distributed under the MIT license.


TEST_NAME="calculate scaling-factor siQ YAML"

# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

scr_prs="${ROOT_REPO}/scripts/parse_metadata_siqchip.py"
cfg_met="${ROOT_REPO}/data/raw/docs/parse_metadata_siqchip.yml"
cfg_p857="${ROOT_REPO}/data/raw/docs/parse_metadata_siqchip_PRJNA857063.yml"
tbl_p857="${ROOT_REPO}/data/raw/docs/measurements_siqchip_PRJNA857063.tsv"

dir_fix="${ROOT_REPO}/tests/calculate_scaling_factor/fixtures"
dir_cfg="${dir_fix}/config"
dir_met="${dir_fix}/metadata"
dir_tmp="${TEST_DIR_TMP}/calculate_scaling_factor_siq_yaml"
dir_log="${TEST_DIR_LOG}/calculate_scaling_factor_siq_yaml"

tbl_met="${dir_met}/measurements_siqchip.tsv"
tbl_gz="${dir_met}/measurements_siqchip.tsv.gz"
tbl_dup="${dir_met}/measurements_siqchip_duplicate_match.tsv"
tbl_mis="${dir_met}/measurements_siqchip_missing_required.tsv"
tbl_lib="${dir_met}/measurements_siqchip_lib_volume.tsv"
tbl_lib_one="${dir_met}/measurements_siqchip_lib_volume_one_sided.tsv"
tbl_lib_zero="${dir_met}/measurements_siqchip_lib_volume_zero.tsv"
tbl_pre="${dir_met}/measurements_siqchip_precomputed.tsv"
cfg_map="${dir_cfg}/parse_metadata_siqchip_field_to_column.yml"

bam_hho1="${dir_tmp}/IP_WT_G1_Hho1_6336.sc.bam"
bam_hho1_bad="${dir_tmp}/IP_WT_G1_Hho1_9999.sc.bam"
bam_extra="${dir_tmp}/IP_WT_G1_HU_Hho1_6336.sc.bam"
bam_p857="${dir_tmp}/IP_HeLa_DMSO_H3K18ac_rep1.hs.bam"


function run_parse_success() {
    local label="${1:-success}"
    local log="${2:-}"
    local tbl="${3:-}"
    local bam="${4:-}"
    local cfg="${5:-${cfg_met}}"

    if \
        PYTHONDONTWRITEBYTECODE=1 \
        PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
        PYTHONPATH="${ROOT_REPO}" \
        run_capture \
            "calculate scaling-factor siQ YAML ${label}" \
            "${log}" \
            "${py_cmd[@]}" -m scripts.parse_metadata_siqchip \
                --alignment "${bam}" \
                --tbl_met "${tbl}" \
                --cfg "${cfg}" \
                --shell
    then
        record_pass "siQ deterministic parser ${label} exits 0"
    else
        record_fail \
            "siQ deterministic parser ${label} failed; see" \
            "$(print_relpath "${log}")"
    fi
}


function run_parse_failure() {
    local label="${1:-failure}"
    local log="${2:-}"
    local tbl="${3:-}"
    local bam="${4:-}"
    local cfg="${5:-${cfg_met}}"
    local patn="${6:-Error}"

    if \
        PYTHONDONTWRITEBYTECODE=1 \
        PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
        PYTHONPATH="${ROOT_REPO}" \
        run_capture \
            "calculate scaling-factor siQ YAML ${label}" \
            "${log}" \
            "${py_cmd[@]}" -m scripts.parse_metadata_siqchip \
                --alignment "${bam}" \
                --tbl_met "${tbl}" \
                --cfg "${cfg}" \
                --shell
    then
        record_fail "siQ deterministic parser ${label} unexpectedly passed"
    else
        assert_pattern_found \
            "${log}" \
            "${patn}" \
            "siQ deterministic parser ${label}"
    fi
}


rm -rf "${dir_tmp}"
mkdir -p "${dir_tmp}" "${dir_log}"
touch "${bam_hho1}" "${bam_hho1_bad}" "${bam_extra}" "${bam_p857}"

require_files_nonempty \
    "${scr_prs}" \
    "${cfg_met}" \
    "${cfg_p857}" \
    "${cfg_map}" \
    "${tbl_met}" \
    "${tbl_gz}" \
    "${tbl_dup}" \
    "${tbl_mis}" \
    "${tbl_lib}" \
    "${tbl_lib_one}" \
    "${tbl_lib_zero}" \
    "${tbl_pre}" \
    "${tbl_p857}" || {
        finish
        exit $?
    }

require_env_project env_nam || {
    finish
    exit $?
}

# shellcheck disable=SC2154
if [[ -n "${CONDA_DEFAULT_ENV:-}" && "${CONDA_DEFAULT_ENV}" != "base" ]]; then
    py_cmd=( "$(find_python)" )
else
    log="${dir_log}/resolve_env_python.log"
    if \
        run_capture \
            "resolve env python" \
            "${log}" \
            conda run -n "${env_nam}" python -c \
                'import sys; print(sys.executable)'
    then
        IFS= read -r py < "${log}"
        py_cmd=( "${py}" )
    else
        record_fail "failed to resolve python from '${env_nam}'"
        finish
        exit $?
    fi
fi


run_parse_success \
    canonical \
    "${dir_log}/canonical.log" \
    "${tbl_met}" \
    "${bam_hho1}"

assert_pattern_found "${dir_log}/canonical.log" '^export mass_ip=2.7$' \
    "deterministic parser exports required mass_ip"
assert_pattern_found "${dir_log}/canonical.log" '^export vol_in=20$' \
    "deterministic parser maps volume_in to vol_in"
assert_pattern_found "${dir_log}/canonical.log" '^export len_ip=626$' \
    "deterministic parser maps length_ip to len_ip"
assert_pattern_absent "${dir_log}/canonical.log" '^export eqn=' \
    "deterministic parser does not export eqn"

run_parse_success \
    gzip_metadata \
    "${dir_log}/gzip.log" \
    "${tbl_gz}" \
    "${bam_hho1}"
assert_pattern_found "${dir_log}/gzip.log" '^export mass_in=72.5$' \
    "deterministic parser reads gzipped metadata"

run_parse_success \
    field_to_column \
    "${dir_log}/field_to_column.log" \
    "${tbl_met}" \
    "${bam_hho1}" \
    "${cfg_map}"
assert_pattern_found "${dir_log}/field_to_column.log" '^export mass_ip=2.7$' \
    "deterministic parser maps filename id to strain column"

run_parse_success \
    precomputed_values \
    "${dir_log}/precomputed.log" \
    "${tbl_pre}" \
    "${bam_hho1}"
assert_pattern_found "${dir_log}/precomputed.log" '^export dep_ip=3333$' \
    "deterministic parser exports precomputed dep_ip"
assert_pattern_found "${dir_log}/precomputed.log" '^export dep_in=2222$' \
    "deterministic parser exports precomputed dep_in"
assert_pattern_found "${dir_log}/precomputed.log" '^export len_ip=456$' \
    "deterministic parser exports precomputed len_ip"
assert_pattern_found "${dir_log}/precomputed.log" '^export len_in=123$' \
    "deterministic parser exports precomputed len_in"

run_parse_success \
    lib_volume \
    "${dir_log}/lib_volume.log" \
    "${tbl_lib}" \
    "${bam_hho1}"
assert_pattern_found "${dir_log}/lib_volume.log" '^export lib_vol_ip=2$' \
    "deterministic parser exports paired lib_vol_ip"
assert_pattern_found "${dir_log}/lib_volume.log" '^export lib_vol_in=4$' \
    "deterministic parser exports paired lib_vol_in"

run_parse_success \
    prjna857063 \
    "${dir_log}/prjna857063.log" \
    "${tbl_p857}" \
    "${bam_p857}" \
    "${cfg_p857}"
assert_pattern_found "${dir_log}/prjna857063.log" '^export mass_ip=22.68$' \
    "deterministic parser resolves PRJNA857063 reference config"

run_parse_failure \
    token_count \
    "${dir_log}/token_count.log" \
    "${tbl_met}" \
    "${bam_extra}" \
    "${cfg_met}" \
    "split into 6 token"

run_parse_failure \
    zero_rows \
    "${dir_log}/zero_rows.log" \
    "${tbl_met}" \
    "${bam_hho1_bad}" \
    "${cfg_met}" \
    "No metadata row matched"

run_parse_failure \
    multiple_rows \
    "${dir_log}/multiple_rows.log" \
    "${tbl_dup}" \
    "${bam_hho1}" \
    "${cfg_met}" \
    "Multiple metadata rows matched"

run_parse_failure \
    missing_required \
    "${dir_log}/missing_required.log" \
    "${tbl_mis}" \
    "${bam_hho1}" \
    "${cfg_met}" \
    "missing required calculator input"

run_parse_failure \
    one_sided_lib_volume \
    "${dir_log}/one_sided_lib_volume.log" \
    "${tbl_lib_one}" \
    "${bam_hho1}" \
    "${cfg_met}" \
    "lib_vol_ip.*lib_vol_in.*provided together"

run_parse_failure \
    nonpositive_lib_volume \
    "${dir_log}/nonpositive_lib_volume.log" \
    "${tbl_lib_zero}" \
    "${bam_hho1}" \
    "${cfg_met}" \
    "lib_vol_.*must be > 0"

finish
