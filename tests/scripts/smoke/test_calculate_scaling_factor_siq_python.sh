#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_calculate_scaling_factor_siq_python.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="calculate scaling-factor siQ Python"

# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"


function run_siq_success() {
    local label="${1:-success}"
    local log="${2:-}"
    local patn="${3:-}"
    shift 3 || true

    if \
        PYTHONDONTWRITEBYTECODE=1 \
        PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
        PYTHONPATH="${ROOT_REPO}" \
        run_capture \
            "calculate scaling-factor siQ ${label}" \
            "${log}" \
            "${py_cmd[@]}" "$@"
    then
        assert_pattern_found "${log}" "${patn}" "siQ Python ${label}"
    else
        record_fail "siQ Python ${label} failed; see $(print_relpath "${log}")"
    fi
}


function run_siq_failure() {
    local label="${1:-expected failure}"
    local log="${2:-}"
    local patn="${3:-}"
    shift 3 || true

    if \
        PYTHONDONTWRITEBYTECODE=1 \
        PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
        PYTHONPATH="${ROOT_REPO}" \
        run_capture \
            "calculate scaling-factor siQ ${label}" \
            "${log}" \
            "${py_cmd[@]}" "$@"
    then
        record_fail "siQ Python ${label} unexpectedly passed"
    else
        assert_pattern_found "${log}" "${patn}" "siQ Python ${label}"
    fi
}


print_section "${TEST_NAME}"

scr_clc="${ROOT_REPO}/scripts/calculate_scaling_factor_siqchip.py"
dir_log="${TEST_DIR_LOG}/calculate_scaling_factor_siq_python"

log_env_py="${dir_log}/resolve_env_python.log"
log_eq5="${dir_log}/calc_eq5.log"
log_eq5nd="${dir_log}/calc_eq5nd.log"
log_eq6="${dir_log}/calc_eq6.log"
log_eq6nd="${dir_log}/calc_eq6nd.log"
log_eq6nd_lib_volume="${dir_log}/calc_eq6nd_lib_volume.log"
log_prjna857063_lib_volume="${dir_log}/calc_prjna857063_lib_volume.log"
log_rounding="${dir_log}/calc_rounding.log"
log_eq5_missing_depth="${dir_log}/calc_eq5_missing_depth.log"
log_bad_volume="${dir_log}/calc_bad_volume.log"
log_lib_volume_one_sided_input="${dir_log}/calc_lib_volume_one_sided_input.log"
log_lib_volume_zero="${dir_log}/calc_lib_volume_zero.log"
log_bad_eqn="${dir_log}/calc_bad_eqn.log"

calc_arg=(
    --mass_ip 2.7
    --mass_in 72.5
    --vol_all 300
    --vol_in 20
    --len_ip 626
    --len_in 450
    --rnd 24
)

mkdir -p "${dir_log}"

require_files_nonempty "${scr_clc}" || {
    finish
    exit $?
}

require_env_project env_nam || {
    finish
    exit $?
}

if [[ -n "${CONDA_DEFAULT_ENV:-}" && "${CONDA_DEFAULT_ENV}" != "base" ]]; then
    py_cmd=( "$(find_python)" )
else
    # shellcheck disable=SC2154
    if \
        run_capture \
            "resolve env python" \
            "${log_env_py}" \
            conda run -n "${env_nam}" python -c \
                'import sys; print(sys.executable)'
    then
        IFS= read -r py < "${log_env_py}"
        py_cmd=( "${py}" )
    else
        record_fail "failed to resolve python from '${env_nam}'"
        finish
        exit $?
    fi
fi

run_siq_success \
    "computes equation 5" \
    "${log_eq5}" \
    '^0.002974551063126584012769$' \
    -m scripts.calculate_scaling_factor_siqchip \
        --eqn 5 \
        "${calc_arg[@]}" \
        --dep_ip 3 \
        --dep_in 5

run_siq_success \
    "computes equation 5nd" \
    "${log_eq5nd}" \
    '^0.001784730637875950407661$' \
    -m scripts.calculate_scaling_factor_siqchip \
        --eqn 5nd \
        "${calc_arg[@]}"

run_siq_success \
    "computes equation 6" \
    "${log_eq6}" \
    '^0.003187018996207053710829$' \
    -m scripts.calculate_scaling_factor_siqchip \
        --eqn 6 \
        "${calc_arg[@]}" \
        --dep_ip 3 \
        --dep_in 5

run_siq_success \
    "computes equation 6nd" \
    "${log_eq6nd}" \
    '^0.001912211397724232486706$' \
    -m scripts.calculate_scaling_factor_siqchip \
        --eqn 6nd \
        "${calc_arg[@]}"

run_siq_success \
    "computes equation 6nd with library-volume correction" \
    "${log_eq6nd_lib_volume}" \
    '^0.003824422795448464973411$' \
    -m scripts.calculate_scaling_factor_siqchip \
        --eqn 6nd \
        "${calc_arg[@]}" \
        --lib_vol_in 4 \
        --lib_vol_ip 2

run_siq_success \
    "computes PRJNA857063-like library-volume correction" \
    "${log_prjna857063_lib_volume}" \
    '^0.087954632669594509652988$' \
    -m scripts.calculate_scaling_factor_siqchip \
        --eqn 6nd \
        --mass_ip 22.68 \
        --mass_in 98.7 \
        --vol_all 250 \
        --vol_in 50 \
        --len_ip 499 \
        --len_in 382 \
        --lib_vol_in 4 \
        --lib_vol_ip 2 \
        --rnd 24

run_siq_success \
    "strips trailing zeros after rounding" \
    "${log_rounding}" \
    '^0.002$' \
    -m scripts.calculate_scaling_factor_siqchip \
        --eqn 6nd \
        "${calc_arg[@]}" \
        --rnd 3

run_siq_failure \
    "rejects missing depths for equation 5" \
    "${log_eq5_missing_depth}" \
    "require both '--dep_ip' and '--dep_in'" \
    -m scripts.calculate_scaling_factor_siqchip \
        --eqn 5 \
        "${calc_arg[@]}"

run_siq_failure \
    "rejects invalid volume relation" \
    "${log_bad_volume}" \
    "vol_all.*must be greater than.*vol_in" \
    -m scripts.calculate_scaling_factor_siqchip \
        --eqn 6nd \
        --mass_ip 2.7 \
        --mass_in 72.5 \
        --vol_all 20 \
        --vol_in 20 \
        --len_ip 626 \
        --len_in 450

run_siq_failure \
    "rejects one-sided input library volume" \
    "${log_lib_volume_one_sided_input}" \
    "requires both '--lib_vol_ip' and '--lib_vol_in'" \
    -m scripts.calculate_scaling_factor_siqchip \
        --eqn 6nd \
        "${calc_arg[@]}" \
        --lib_vol_in 4

run_siq_failure \
    "rejects zero library volume" \
    "${log_lib_volume_zero}" \
    "lib_vol_in.*> 0" \
    -m scripts.calculate_scaling_factor_siqchip \
        --eqn 6nd \
        "${calc_arg[@]}" \
        --lib_vol_in 0 \
        --lib_vol_ip 2

run_siq_failure \
    "rejects invalid equation" \
    "${log_bad_eqn}" \
    "invalid choice" \
    -m scripts.calculate_scaling_factor_siqchip \
        --eqn bogus \
        "${calc_arg[@]}"

finish
