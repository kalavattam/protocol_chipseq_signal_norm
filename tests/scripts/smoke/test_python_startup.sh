#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_python_startup.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="python startup"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

#  Resolve the project environment locally for dependency-backed Python checks
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
    log="${TEST_DIR_LOG}/python_startup/resolve_env_python.log"

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
    "${py_cmd[@]}" - << PY
import sys
raise SystemExit(0 if sys.version_info >= (3, 10) else 1)
PY
then
    record_fail \
        "$("${py_cmd[@]}" --version 2>&1) is older than Python 3.10"
    finish
    exit $?
fi

#  Run lightweight Python syntax checks without writing pycache in-place
while IFS= read -r file; do
    rel="$(print_relpath "${file}")"
    log="${TEST_DIR_LOG}/python_compile/${rel//\//__}.log"

    if \
        PYTHONDONTWRITEBYTECODE=1 \
        PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
        PYTHONPATH="${ROOT_REPO}" \
        run_capture \
            "python compile ${rel}" \
            "${log}" \
            "${py_cmd[@]}" -m py_compile "${file}"
    then
        record_pass "python syntax ${rel}"
    else
        record_fail "python syntax ${rel}; see $(print_relpath "${log}")"
    fi
done < <(
    find "${ROOT_REPO}/scripts" \
        -path "${ROOT_REPO}/scripts/blog" \
        -prune -o -type f -name '*.py' -print \
        | sort
)

#  Run --help with the project environment so dependency imports are available
hlp_scr=(
    "scripts/add_coeffs_namespaced.py"
    "scripts/calculate_scaling_factor_siq_chip.py"
    "scripts/calculate_scaling_factor_spike.py"
    "scripts/compute_input_floor.py"
    "scripts/compute_pseudo.py"
    "scripts/compute_signal.py"
    "scripts/compute_signal_ratio.py"
    "scripts/merge_bins_bdg.py"
    "scripts/parse_metadata_siq_chip.py"
    "scripts/relativize_scaling_factors.py"
    "scripts/sum_bdg.py"
)

for rel in "${hlp_scr[@]}"; do
    file="${ROOT_REPO}/${rel}"
    [[ -f "${file}" ]] || {
        record_skip "python --help ${rel}: file not present"
        continue
    }

    log="${TEST_DIR_LOG}/python_help/${rel//\//__}.log"
    if \
        PYTHONDONTWRITEBYTECODE=1 \
        PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
        PYTHONPATH="${ROOT_REPO}" \
        run_capture \
            "python help ${rel}" \
            "${log}" \
            "${py_cmd[@]}" "${file}" --help
    then
        record_pass "python --help ${rel}"
    else
        record_warn \
            "python --help ${rel} failed; see $(print_relpath "${log}")"
    fi
done

#  Confirm the direct spike-in calculator defaults to the ChIP-Rx ratio
log="${TEST_DIR_LOG}/python_startup/calculate_scaling_factor_spike_default.log"
if \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
    PYTHONPATH="${ROOT_REPO}" \
    run_capture \
        "python calculate-scaling-factor spike default" \
        "${log}" \
        "${py_cmd[@]}" -m scripts.calculate_scaling_factor_spike \
            --main_ip 100 \
            --spike_ip 10 \
            --main_in 100 \
            --spike_in 20
then
    assert_pattern_found \
        "${log}" \
        '^2$' \
        "calculate_scaling_factor_spike.py defaults to chiprx_alpha_ratio"
else
    record_fail \
        "calculate_scaling_factor_spike.py default calculation failed; see" \
        "$(print_relpath "${log}")"
fi

finish
