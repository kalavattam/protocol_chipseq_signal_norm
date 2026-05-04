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

rec_section "${TEST_NAME}"

#  Skip Python checks cleanly when Python is unavailable or too old
if ! \
    py="$(find_python)"
then
    rec_skip "python/python3 unavailable; skipping Python smoke checks"
    finish
    exit $?
fi

if ! \
    check_python_ge_310 "${py}"
then
    rec_skip \
        "$("${py}" --version 2>&1) is older than Python 3.10; skipping" \
        "Python smoke checks"
    finish
    exit $?
fi

#  Run lightweight Python syntax checks without writing pycache in-place
while IFS= read -r file; do
    rel="$(rec_relpath "${file}")"
    log="${TEST_DIR_LOG}/python_compile/${rel//\//__}.log"

    if \
        PYTHONDONTWRITEBYTECODE=1 \
        PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
        PYTHONPATH="${ROOT_REPO}" \
        run_capture \
            "python compile ${rel}" "${log}" "${py}" -m py_compile "${file}"
    then
        rec_pass "python syntax ${rel}"
    else
        rec_fail "python syntax ${rel}; see $(rec_relpath "${log}")"
    fi
done < <(
    find "${ROOT_REPO}/scripts" \
        -path "${ROOT_REPO}/scripts/blog" -prune -o \
        -type f -name '*.py' -print \
        | sort
)

#  Run --help only for Python scripts known to be safe in a bare environment
safe_help=(
    "scripts/add_coeffs_namespaced.py"
    "scripts/calculate_scaling_factor_siq_chip.py"
    "scripts/calculate_scaling_factor_spike.py"
    "scripts/compute_input_floor.py"
    "scripts/compute_pseudo.py"
    "scripts/compute_signal_ratio.py"
    "scripts/merge_bins_bdg.py"
    "scripts/sum_bdg.py"
)

#  Keep optional-dependency startup imports visible as skips
rec_skip \
    "python --help skipped for scripts/compute_signal.py because it imports" \
    "pysam at startup"
rec_skip \
    "python --help skipped for scripts/parse_metadata_siq_chip.py because it" \
    "imports yaml at startup"
rec_skip \
    "python --help skipped for scripts/relativize_scaling_factors.py because" \
    "it imports pandas at startup"

for rel in "${safe_help[@]}"; do
    file="${ROOT_REPO}/${rel}"
    [[ -f "${file}" ]] || {
        rec_skip "python --help ${rel}: file not present"
        continue
    }

    log="${TEST_DIR_LOG}/python_help/${rel//\//__}.log"
    if PYTHONDONTWRITEBYTECODE=1 \
        PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
        PYTHONPATH="${ROOT_REPO}" \
        run_capture "python help ${rel}" "${log}" "${py}" "${file}" --help
    then
        rec_pass "python --help ${rel}"
    else
        rec_warn "python --help ${rel} failed; see $(rec_relpath "${log}")"
    fi
done

finish
