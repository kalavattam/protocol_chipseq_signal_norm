#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_calculate_scaling_factor_spike_python.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="calculate scaling-factor spike Python"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

scr_py="${ROOT_REPO}/scripts/calculate_scaling_factor_spike.py"
dir_log="${TEST_DIR_LOG}/calculate_scaling_factor_spike_python"

cnt_arg=(
    --main_ip 90
    --spike_ip 10
    --main_in 80
    --spike_in 20
)


#  Run a successful Python calculator case
function run_spike_success() {
    local label="${1:-success}"
    local fil_log="${2:-}"
    local patn="${3:-}"
    shift 3 || true

    if \
        PYTHONDONTWRITEBYTECODE=1 \
        PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
        PYTHONPATH="${ROOT_REPO}" \
        run_capture \
            "calculate scaling-factor spike ${label}" \
            "${fil_log}" \
            "${py_cmd[@]}" -m scripts.calculate_scaling_factor_spike \
                "${cnt_arg[@]}" \
                "$@"
    then
        assert_pattern_found \
            "${fil_log}" \
            "${patn}" \
            "calculate_scaling_factor_spike.py ${label}"
    else
        record_fail \
            "calculate_scaling_factor_spike.py ${label} failed; see" \
            "$(print_relpath "${fil_log}")"
    fi
}


#  Run an expected Python calculator failure
function run_spike_failure() {
    local label="${1:-expected failure}"
    local fil_log="${2:-}"
    local patn="${3:-}"
    shift 3 || true

    if \
        PYTHONDONTWRITEBYTECODE=1 \
        PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
        PYTHONPATH="${ROOT_REPO}" \
        run_capture \
            "calculate scaling-factor spike ${label}" \
            "${fil_log}" \
            "${py_cmd[@]}" -m scripts.calculate_scaling_factor_spike \
                "$@"
    then
        record_fail \
            "calculate_scaling_factor_spike.py ${label} unexpectedly pass"
    else
        assert_pattern_found \
            "${fil_log}" \
            "${patn}" \
            "calculate_scaling_factor_spike.py ${label}"
    fi
}



mkdir -p "${dir_log}"

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

if [[ -f "${scr_py}" ]]; then
    record_pass "calculate_scaling_factor_spike.py script exists"
else
    record_fail \
        "calculate_scaling_factor_spike.py script missing:" \
        "$(print_relpath "${scr_py}")"
    finish
    exit $?
fi


#  Plain output should compute each canonical coefficient
run_spike_success \
    "defaults to chiprx_alpha_ratio" \
    "${dir_log}/default.log" \
    '^2$'

run_spike_success \
    "computes fractional" \
    "${dir_log}/fractional.log" \
    '^2$' \
    --coef fractional

run_spike_success \
    "computes chiprx_alpha_ip" \
    "${dir_log}/chiprx_alpha_ip.log" \
    '^100000$' \
    --coef chiprx_alpha_ip

run_spike_success \
    "computes chiprx_alpha_in" \
    "${dir_log}/chiprx_alpha_in.log" \
    '^50000$' \
    --coef chiprx_alpha_in

run_spike_success \
    "computes chiprx_alpha_ratio" \
    "${dir_log}/chiprx_alpha_ratio.log" \
    '^2$' \
    --coef chiprx_alpha_ratio

run_spike_success \
    "computes rxinput_alpha" \
    "${dir_log}/rxinput_alpha.log" \
    '^20000$' \
    --coef rxinput_alpha


#  Structured output should report canonical names and rounded values
run_spike_success \
    "canonicalizes chiprx_ratio alias in TSV" \
    "${dir_log}/alias_tsv.log" \
    $'^chiprx_alpha_ratio\t2$' \
    --coef chiprx_ratio \
    --format tsv

run_spike_success \
    "emits all coefficients in TSV" \
    "${dir_log}/all_tsv.log" \
    $'^rxinput_alpha\t20000$' \
    --coef all \
    --format tsv

assert_pattern_found \
    "${dir_log}/all_tsv.log" \
    $'^coef\tvalue$' \
    "calculate_scaling_factor_spike.py all TSV has header"

assert_pattern_found \
    "${dir_log}/all_tsv.log" \
    $'^fractional\t2$' \
    "calculate_scaling_factor_spike.py all TSV has fractional"

assert_pattern_found \
    "${dir_log}/all_tsv.log" \
    $'^chiprx_alpha_ip\t100000$' \
    "calculate_scaling_factor_spike.py all TSV has chiprx_alpha_ip"

assert_pattern_found \
    "${dir_log}/all_tsv.log" \
    $'^chiprx_alpha_in\t50000$' \
    "calculate_scaling_factor_spike.py all TSV has chiprx_alpha_in"

assert_pattern_found \
    "${dir_log}/all_tsv.log" \
    $'^chiprx_alpha_ratio\t2$' \
    "calculate_scaling_factor_spike.py all TSV has chiprx_alpha_ratio"

run_spike_success \
    "emits all coefficients in JSON" \
    "${dir_log}/all_json.log" \
    '"rxinput_alpha": 20000.0' \
    --coef all \
    --format json

assert_pattern_found \
    "${dir_log}/all_json.log" \
    '"fractional": 2.0' \
    "calculate_scaling_factor_spike.py all JSON has fractional"

assert_pattern_found \
    "${dir_log}/all_json.log" \
    '"chiprx_alpha_ratio": 2.0' \
    "calculate_scaling_factor_spike.py all JSON has chiprx_alpha_ratio"

run_spike_success \
    "strips trailing zeros after rounding" \
    "${dir_log}/rounding.log" \
    '^2$' \
    --coef chiprx_alpha_ratio \
    --rnd 24


#  Invalid CLI combinations and invalid counts should fail clearly
run_spike_failure \
    "rejects invalid coefficient" \
    "${dir_log}/invalid_coef.log" \
    "Invalid --coef 'bogus'" \
    "${cnt_arg[@]}" \
    --coef bogus

run_spike_failure \
    "rejects invalid format" \
    "${dir_log}/invalid_format.log" \
    "invalid choice" \
    "${cnt_arg[@]}" \
    --format bogus

run_spike_failure \
    "rejects plain all coefficients" \
    "${dir_log}/plain_all.log" \
    "Format 'plain' requires a single coefficient" \
    "${cnt_arg[@]}" \
    --coef all \
    --format plain

run_spike_failure \
    "rejects zero spike_ip denominator" \
    "${dir_log}/zero_spike_ip.log" \
    "'spike_ip' is 0" \
    --main_ip 90 \
    --spike_ip 0 \
    --main_in 80 \
    --spike_in 20 \
    --coef chiprx_alpha_ratio

run_spike_failure \
    "rejects negative counts" \
    "${dir_log}/negative_count.log" \
    "main_ip.*>= 0" \
    --main_ip -1 \
    --spike_ip 10 \
    --main_in 80 \
    --spike_in 20

finish
