#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_calculate_scaling_factor_spike_python.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="calculate scaling-factor spike Python"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


#  Run a successful Python calculator case
function run_spike_success() {
    local label="${1:-success}"
    local fil_log="${2:-}"
    local patn="${3:-}"
    shift 3 || true

    if \
        PYTHONDONTWRITEBYTECODE=1 \
        PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
        PYTHONPATH="${ROOT_REPO}/src" \
        run_capture \
            "calculate scaling-factor spike ${label}" \
            "${fil_log}" \
            "${py_cmd[@]}" -m protocol_chipseq_signal_norm.cli.calculate_scaling_factor_spike \
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
        PYTHONPATH="${ROOT_REPO}/src" \
        run_capture \
            "calculate scaling-factor spike ${label}" \
            "${fil_log}" \
            "${py_cmd[@]}" -m protocol_chipseq_signal_norm.cli.calculate_scaling_factor_spike \
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


scr_py="${ROOT_REPO}/src/protocol_chipseq_signal_norm/cli/calculate_scaling_factor_spike.py"
dir_log="${TEST_DIR_LOG}/calculate_scaling_factor_spike_python"

log_env_py="${dir_log}/resolve_env_python.log"
log_default="${dir_log}/default.log"
log_fractional="${dir_log}/fractional.log"
log_chiprx_alpha_ip="${dir_log}/chiprx_alpha_ip.log"
log_chiprx_alpha_in="${dir_log}/chiprx_alpha_in.log"
log_chiprx_alpha_ratio="${dir_log}/chiprx_alpha_ratio.log"
log_rxinput_alpha="${dir_log}/rxinput_alpha.log"
log_alias_tsv="${dir_log}/alias_tsv.log"
log_all_tsv="${dir_log}/all_tsv.log"
log_all_json="${dir_log}/all_json.log"
log_rounding="${dir_log}/rounding.log"
log_invalid_coef="${dir_log}/invalid_coef.log"
log_invalid_format="${dir_log}/invalid_format.log"
log_plain_all="${dir_log}/plain_all.log"
log_zero_spike_ip="${dir_log}/zero_spike_ip.log"
log_negative_count="${dir_log}/negative_count.log"

cnt_arg=(
    --main_ip 90
    --spike_ip 10
    --main_in 80
    --spike_in 20
)


print_section "${TEST_NAME}"

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
    # shellcheck disable=SC2154
    if \
        run_capture \
            "resolve env python" \
            "${log_env_py}" \
            conda run -n "${env_nam}" python -c \
                'import sys; print(sys.executable)'
    then
        IFS= read -r py < "${log_env_py}"
    else
        record_fail \
            "failed to resolve python from '${env_nam}'; see" \
            "$(print_relpath "${log_env_py}")"
        finish
        exit $?
    fi

    if [[ -z "${py}" || ! -x "${py}" ]]; then
        record_fail \
            "resolved python is not executable; see" \
            "$(print_relpath "${log_env_py}")"
        finish
        exit $?
    fi

    py_cmd=( "${py}" )
fi

if ! \
    check_python_ge_311 "${py_cmd[0]}"
then
    record_fail \
        "$("${py_cmd[@]}" --version 2>&1) is older than Python 3.11"
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
    "${log_default}" \
    '^2$'

run_spike_success \
    "computes fractional" \
    "${log_fractional}" \
    '^2$' \
    --coef fractional

run_spike_success \
    "computes chiprx_alpha_ip" \
    "${log_chiprx_alpha_ip}" \
    '^100000$' \
    --coef chiprx_alpha_ip

run_spike_success \
    "computes chiprx_alpha_in" \
    "${log_chiprx_alpha_in}" \
    '^50000$' \
    --coef chiprx_alpha_in

run_spike_success \
    "computes chiprx_alpha_ratio" \
    "${log_chiprx_alpha_ratio}" \
    '^2$' \
    --coef chiprx_alpha_ratio

run_spike_success \
    "computes rxinput_alpha" \
    "${log_rxinput_alpha}" \
    '^20000$' \
    --coef rxinput_alpha


#  Structured output should report canonical names and rounded values
run_spike_success \
    "canonicalizes chiprx_ratio alias in TSV" \
    "${log_alias_tsv}" \
    $'^chiprx_alpha_ratio\t2$' \
    --coef chiprx_ratio \
    --format tsv

run_spike_success \
    "emits all coefficients in TSV" \
    "${log_all_tsv}" \
    $'^rxinput_alpha\t20000$' \
    --coef all \
    --format tsv

assert_pattern_found \
    "${log_all_tsv}" \
    $'^coef\tvalue$' \
    "calculate_scaling_factor_spike.py all TSV has header"

assert_pattern_found \
    "${log_all_tsv}" \
    $'^fractional\t2$' \
    "calculate_scaling_factor_spike.py all TSV has fractional"

assert_pattern_found \
    "${log_all_tsv}" \
    $'^chiprx_alpha_ip\t100000$' \
    "calculate_scaling_factor_spike.py all TSV has chiprx_alpha_ip"

assert_pattern_found \
    "${log_all_tsv}" \
    $'^chiprx_alpha_in\t50000$' \
    "calculate_scaling_factor_spike.py all TSV has chiprx_alpha_in"

assert_pattern_found \
    "${log_all_tsv}" \
    $'^chiprx_alpha_ratio\t2$' \
    "calculate_scaling_factor_spike.py all TSV has chiprx_alpha_ratio"

run_spike_success \
    "emits all coefficients in JSON" \
    "${log_all_json}" \
    '"rxinput_alpha": 20000.0' \
    --coef all \
    --format json

assert_pattern_found \
    "${log_all_json}" \
    '"fractional": 2.0' \
    "calculate_scaling_factor_spike.py all JSON has fractional"

assert_pattern_found \
    "${log_all_json}" \
    '"chiprx_alpha_ratio": 2.0' \
    "calculate_scaling_factor_spike.py all JSON has chiprx_alpha_ratio"

run_spike_success \
    "strips trailing zeros after rounding" \
    "${log_rounding}" \
    '^2$' \
    --coef chiprx_alpha_ratio \
    --dp 24


#  Invalid CLI combinations and invalid counts should fail clearly
run_spike_failure \
    "rejects invalid coefficient" \
    "${log_invalid_coef}" \
    "Invalid --coef 'bogus'" \
    "${cnt_arg[@]}" \
    --coef bogus

run_spike_failure \
    "rejects invalid format" \
    "${log_invalid_format}" \
    "invalid choice" \
    "${cnt_arg[@]}" \
    --format bogus

run_spike_failure \
    "rejects plain all coefficients" \
    "${log_plain_all}" \
    "Format 'plain' requires a single coefficient" \
    "${cnt_arg[@]}" \
    --coef all \
    --format plain

run_spike_failure \
    "rejects zero spike_ip denominator" \
    "${log_zero_spike_ip}" \
    "'spike_ip' is 0" \
    --main_ip 90 \
    --spike_ip 0 \
    --main_in 80 \
    --spike_in 20 \
    --coef chiprx_alpha_ratio

run_spike_failure \
    "rejects negative counts" \
    "${log_negative_count}" \
    "main_ip.*>= 0" \
    --main_ip -1 \
    --spike_ip 10 \
    --main_in 80 \
    --spike_in 20

finish
