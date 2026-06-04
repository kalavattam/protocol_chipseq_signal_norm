#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_calculate_scaling_factor_siq_python.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="calculate scaling-factor siQ Python"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

scr_clc="${ROOT_REPO}/scripts/calculate_scaling_factor_siq_chip.py"
scr_prs="${ROOT_REPO}/scripts/parse_metadata_siq_chip.py"
cfg_met="${ROOT_REPO}/data/raw/docs/parse_metadata_siq_chip.yaml"
dir_fix="${ROOT_REPO}/tests/calculate_scaling_factor/fixtures"
dir_met="${dir_fix}/metadata"
dir_tmp="${TEST_DIR_TMP}/calculate_scaling_factor_siq_python"
dir_log="${TEST_DIR_LOG}/calculate_scaling_factor_siq_python"

tbl_met="${dir_met}/measurements_siqchip.tsv"
tbl_als="${dir_met}/measurements_siqchip_aliases.tsv"
tbl_mis="${dir_met}/measurements_siqchip_missing_required.tsv"
tbl_uns="${dir_met}/measurements_siqchip_unsupported_alias.tsv"

bam_ip_hho1_g1="${dir_tmp}/IP_WT_G1_Hho1_6336.sc.bam"
bam_in_hho1_g1="${dir_tmp}/in_WT_G1_Hho1_6336.sc.bam"
cram_ip_hho1_g1="${dir_tmp}/IP_WT_G1_Hho1_6336.sp.cram"
bam_ip_hho1_g2m="${dir_tmp}/IP_WT_G2M_Hho1_6336.sc.bam"
bam_ip_hmo1_g1="${dir_tmp}/IP_WT_G1_Hmo1_7750.sc.bam"
bam_bad="${dir_tmp}/bad.bam"
bam_bad_match="${dir_tmp}/IP_WT_Q_Hho1_9999.sc.bam"
bam_bad_assay="${dir_tmp}/FLAGIP_WT_G1_Hho1_6336.sc.bam"

calc_arg=(
    --mass_ip 2.7
    --mass_in 72.5
    --vol_all 300
    --vol_in 20
    --len_ip 626
    --len_in 450
    --rnd 24
)


#  Run a successful siQ Python case
function run_siq_success() {
    local label="${1:-success}"
    local fil_log="${2:-}"
    local patn="${3:-}"
    shift 3 || true

    if \
        PYTHONDONTWRITEBYTECODE=1 \
        PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
        PYTHONPATH="${ROOT_REPO}" \
        run_capture \
            "calculate scaling-factor siQ ${label}" \
            "${fil_log}" \
            "${py_cmd[@]}" "$@"
    then
        assert_pattern_found \
            "${fil_log}" \
            "${patn}" \
            "siQ Python ${label}"
    else
        record_fail \
            "siQ Python ${label} failed; see" \
            "$(print_relpath "${fil_log}")"
    fi
}


#  Run an expected siQ Python failure
function run_siq_failure() {
    local label="${1:-expected failure}"
    local fil_log="${2:-}"
    local patn="${3:-}"
    shift 3 || true

    if \
        PYTHONDONTWRITEBYTECODE=1 \
        PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
        PYTHONPATH="${ROOT_REPO}" \
        run_capture \
            "calculate scaling-factor siQ ${label}" \
            "${fil_log}" \
            "${py_cmd[@]}" "$@"
    then
        record_fail "siQ Python ${label} unexpectedly passed"
    else
        assert_pattern_found \
            "${fil_log}" \
            "${patn}" \
            "siQ Python ${label}"
    fi
}


#  Run a successful metadata parser case
function run_parse_success() {
    local label="${1:-parse success}"
    local fil_log="${2:-}"
    local tbl_lcl="${3:-}"
    local bam_lcl="${4:-}"
    shift 4 || true

    run_siq_success \
        "${label}" \
        "${fil_log}" \
        '^export eqn=6nd$' \
        -m scripts.parse_metadata_siq_chip \
            --bam "${bam_lcl}" \
            --tbl_met "${tbl_lcl}" \
            --cfg "${cfg_met}" \
            --eqn 6nd \
            --shell \
            "$@"
}


rm -rf "${dir_tmp}"
mkdir -p "${dir_tmp}" "${dir_log}"

touch \
    "${bam_ip_hho1_g1}" \
    "${bam_in_hho1_g1}" \
    "${cram_ip_hho1_g1}" \
    "${bam_ip_hho1_g2m}" \
    "${bam_ip_hmo1_g1}" \
    "${bam_bad}" \
    "${bam_bad_match}" \
    "${bam_bad_assay}"

require_files_nonempty \
    "${scr_clc}" \
    "${scr_prs}" \
    "${cfg_met}" \
    "${tbl_met}" \
    "${tbl_als}" \
    "${tbl_mis}" \
    "${tbl_uns}" || {
        finish
        exit $?
    }

require_env_project env_nam || {
    finish
    exit $?
}

if [[ -n "${CONDA_DEFAULT_ENV:-}" && "${CONDA_DEFAULT_ENV}" != "base" ]]; then
    if ! \
        py="$(find_python)"
    then
        record_fail "active project environment has no python/python3 on PATH"
        finish
        exit $?
    fi

    py_cmd=( "${py}" )
else
    log="${dir_log}/resolve_env_python.log"

    # shellcheck disable=SC2154
    if \
        run_capture \
            "resolve env python" \
            "${log}" \
            conda run -n "${env_nam}" python -c \
                'import sys; print(sys.executable)'
    then
        IFS= read -r py < "${log}"
    else
        record_fail \
            "failed to resolve python from '${env_nam}'; see" \
            "$(print_relpath "${log}")"
        finish
        exit $?
    fi

    if [[ -z "${py}" || ! -x "${py}" ]]; then
        record_fail \
            "resolved python is not executable; see $(print_relpath "${log}")"
        finish
        exit $?
    fi

    py_cmd=( "${py}" )
fi

if ! \
    check_python_ge_310 "${py_cmd[0]}"
then
    record_fail \
        "$("${py_cmd[@]}" --version 2>&1) is older than Python 3.10"
    finish
    exit $?
fi


#  Metadata parser happy paths should resolve production-YAML shell fields
run_parse_success \
    "parses canonical IP BAM basename" \
    "${dir_log}/parse_hho1_ip_bam.log" \
    "${tbl_met}" \
    "${bam_ip_hho1_g1}"

assert_pattern_found \
    "${dir_log}/parse_hho1_ip_bam.log" \
    '^export mass_ip=2.7$' \
    "siQ metadata parser extracts Hho1 G1 mass_ip"

assert_pattern_found \
    "${dir_log}/parse_hho1_ip_bam.log" \
    '^export mass_in=72.5$' \
    "siQ metadata parser extracts Hho1 G1 mass_in"

assert_pattern_found \
    "${dir_log}/parse_hho1_ip_bam.log" \
    '^export vol_all=300$' \
    "siQ metadata parser extracts Hho1 G1 vol_all"

assert_pattern_found \
    "${dir_log}/parse_hho1_ip_bam.log" \
    '^export vol_in=20$' \
    "siQ metadata parser extracts Hho1 G1 vol_in"

assert_pattern_found \
    "${dir_log}/parse_hho1_ip_bam.log" \
    '^export len_ip=626$' \
    "siQ metadata parser extracts Hho1 G1 len_ip"

assert_pattern_found \
    "${dir_log}/parse_hho1_ip_bam.log" \
    '^export len_in=450$' \
    "siQ metadata parser extracts Hho1 G1 len_in"

run_parse_success \
    "parses canonical input BAM basename" \
    "${dir_log}/parse_hho1_input_bam.log" \
    "${tbl_met}" \
    "${bam_in_hho1_g1}"

run_parse_success \
    "strips CRAM extension and trailing organism tag" \
    "${dir_log}/parse_hho1_cram.log" \
    "${tbl_met}" \
    "${cram_ip_hho1_g1}"

run_parse_success \
    "parses second state" \
    "${dir_log}/parse_hho1_g2m.log" \
    "${tbl_met}" \
    "${bam_ip_hho1_g2m}"

assert_pattern_found \
    "${dir_log}/parse_hho1_g2m.log" \
    '^export mass_ip=6.6$' \
    "siQ metadata parser extracts Hho1 G2M mass_ip"

run_parse_success \
    "parses second factor" \
    "${dir_log}/parse_hmo1_g1.log" \
    "${tbl_met}" \
    "${bam_ip_hmo1_g1}"

assert_pattern_found \
    "${dir_log}/parse_hmo1_g1.log" \
    '^export mass_ip=8.4$' \
    "siQ metadata parser extracts Hmo1 G1 mass_ip"

run_parse_success \
    "parses curated column aliases" \
    "${dir_log}/parse_aliases.log" \
    "${tbl_als}" \
    "${bam_ip_hho1_g1}"

assert_pattern_found \
    "${dir_log}/parse_aliases.log" \
    '^export len_ip=626$' \
    "siQ metadata parser normalizes length IP alias"


#  Metadata parser failures should be clear and targeted
run_siq_failure \
    "rejects stdin BAM placeholder" \
    "${dir_log}/parse_stdin_bam.log" \
    "'--bam -' not supported" \
    -m scripts.parse_metadata_siq_chip \
        --bam - \
        --tbl_met "${tbl_met}" \
        --cfg "${cfg_met}" \
        --shell

run_siq_failure \
    "rejects unparseable filename" \
    "${dir_log}/parse_bad_filename.log" \
    "Could not parse filename" \
    -m scripts.parse_metadata_siq_chip \
        --bam "${bam_bad}" \
        --tbl_met "${tbl_met}" \
        --cfg "${cfg_met}" \
        --shell

run_siq_failure \
    "rejects unmatched metadata row" \
    "${dir_log}/parse_no_match.log" \
    "No matching row found" \
    -m scripts.parse_metadata_siq_chip \
        --bam "${bam_bad_match}" \
        --tbl_met "${tbl_met}" \
        --cfg "${cfg_met}" \
        --shell

run_siq_failure \
    "rejects unsupported assay token" \
    "${dir_log}/parse_bad_assay.log" \
    "expected field 'assay'" \
    -m scripts.parse_metadata_siq_chip \
        --bam "${bam_bad_assay}" \
        --tbl_met "${tbl_met}" \
        --cfg "${cfg_met}" \
        --shell

run_siq_failure \
    "reports missing required field" \
    "${dir_log}/parse_missing_required.log" \
    "Missing required column.*len_ip" \
    -m scripts.parse_metadata_siq_chip \
        --bam "${bam_ip_hho1_g1}" \
        --tbl_met "${tbl_mis}" \
        --cfg "${cfg_met}" \
        --shell

assert_pattern_found \
    "${dir_log}/parse_missing_required.log" \
    "Accepted aliases include: len_ip:" \
    "siQ metadata parser lists accepted aliases for missing field"

run_siq_failure \
    "rejects unsupported required-column alias" \
    "${dir_log}/parse_unsupported_alias.log" \
    "Missing required column.*vol_in" \
    -m scripts.parse_metadata_siq_chip \
        --bam "${bam_ip_hho1_g1}" \
        --tbl_met "${tbl_uns}" \
        --cfg "${cfg_met}" \
        --shell

assert_pattern_found \
    "${dir_log}/parse_unsupported_alias.log" \
    "Accepted aliases include: vol_in:" \
    "siQ metadata parser lists accepted aliases for unsupported alias"


#  Direct siQ calculator coverage should lock equation behavior
run_siq_success \
    "computes equation 5" \
    "${dir_log}/calc_eq5.log" \
    '^0.002974551063126584012769$' \
    -m scripts.calculate_scaling_factor_siq_chip \
        --eqn 5 \
        "${calc_arg[@]}" \
        --dep_ip 3 \
        --dep_in 5

run_siq_success \
    "computes equation 5nd" \
    "${dir_log}/calc_eq5nd.log" \
    '^0.001784730637875950407661$' \
    -m scripts.calculate_scaling_factor_siq_chip \
        --eqn 5nd \
        "${calc_arg[@]}"

run_siq_success \
    "computes equation 6" \
    "${dir_log}/calc_eq6.log" \
    '^0.003187018996207053710829$' \
    -m scripts.calculate_scaling_factor_siq_chip \
        --eqn 6 \
        "${calc_arg[@]}" \
        --dep_ip 3 \
        --dep_in 5

run_siq_success \
    "computes equation 6nd" \
    "${dir_log}/calc_eq6nd.log" \
    '^0.001912211397724232486706$' \
    -m scripts.calculate_scaling_factor_siq_chip \
        --eqn 6nd \
        "${calc_arg[@]}"

run_siq_success \
    "strips trailing zeros after rounding" \
    "${dir_log}/calc_rounding.log" \
    '^0.002$' \
    -m scripts.calculate_scaling_factor_siq_chip \
        --eqn 6nd \
        "${calc_arg[@]}" \
        --rnd 3


#  Direct calculator failures should reject invalid inputs
run_siq_failure \
    "rejects missing depths for equation 5" \
    "${dir_log}/calc_eq5_missing_depth.log" \
    "require both '--dep_ip' and '--dep_in'" \
    -m scripts.calculate_scaling_factor_siq_chip \
        --eqn 5 \
        "${calc_arg[@]}"

run_siq_failure \
    "rejects invalid volume relation" \
    "${dir_log}/calc_bad_volume.log" \
    "vol_all.*must be greater than.*vol_in" \
    -m scripts.calculate_scaling_factor_siq_chip \
        --eqn 6nd \
        --mass_ip 2.7 \
        --mass_in 72.5 \
        --vol_all 20 \
        --vol_in 20 \
        --len_ip 626 \
        --len_in 450

run_siq_failure \
    "rejects nonpositive mass" \
    "${dir_log}/calc_bad_mass.log" \
    "mass_ip.*> 0" \
    -m scripts.calculate_scaling_factor_siq_chip \
        --eqn 6nd \
        --mass_ip 0 \
        --mass_in 72.5 \
        --vol_all 300 \
        --vol_in 20 \
        --len_ip 626 \
        --len_in 450

run_siq_failure \
    "rejects invalid equation" \
    "${dir_log}/calc_bad_eqn.log" \
    "invalid choice" \
    -m scripts.calculate_scaling_factor_siq_chip \
        --eqn bogus \
        "${calc_arg[@]}"

finish
