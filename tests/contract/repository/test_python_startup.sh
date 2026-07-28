#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_python_startup.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="python startup"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


log_env_rslv="${TEST_DIR_LOG}/python_startup/resolve_env_python.log"
log_pol_py="${TEST_DIR_LOG}/python_startup/python_version_policy.log"
log_def_spk="${TEST_DIR_LOG}/python_startup/calculate_scaling_factor_spike_default.log"
log_ret_in="${TEST_DIR_LOG}/python_startup/compute_signal_retired_input.log"
log_ret_out="${TEST_DIR_LOG}/python_startup/compute_signal_retired_output.log"
log_ret_i="${TEST_DIR_LOG}/python_startup/compute_signal_retired_short_i.log"
ret_alias_in="--in""file"
ret_alias_out="--out""file"
ret_alias_short="-""i"

hlp_scr=(
    "src/protocol_chipseq_signal_norm/cli/add_coeffs_namespaced.py"
    "src/protocol_chipseq_signal_norm/cli/calculate_scaling_factor_siqchip.py"
    "src/protocol_chipseq_signal_norm/cli/calculate_scaling_factor_spike.py"
    "src/protocol_chipseq_signal_norm/cli/compute_input_floor.py"
    "src/protocol_chipseq_signal_norm/cli/compute_pseudo.py"
    "src/protocol_chipseq_signal_norm/cli/compute_signal.py"
    "src/protocol_chipseq_signal_norm/cli/compute_signal_ratio.py"
    "src/protocol_chipseq_signal_norm/cli/merge_bins_bdg.py"
    "src/protocol_chipseq_signal_norm/cli/parse_metadata_siqchip.py"
    "src/protocol_chipseq_signal_norm/cli/relativize_scaling_factors.py"
    "src/protocol_chipseq_signal_norm/cli/sum_bdg.py"
)


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
    # shellcheck disable=SC2154
    if \
        run_capture \
            "resolve env python" \
            "${log_env_rslv}" \
            conda run -n "${env_nam}" python -c \
                'import sys; print(sys.executable)'
    then
        IFS= read -r py < "${log_env_rslv}"
    else
        record_fail \
            "failed to resolve python from '${env_nam}'; see" \
            "$(print_relpath "${log_env_rslv}")"
        finish
        exit $?
    fi

    if [[ -z "${py}" || ! -x "${py}" ]]; then
        record_fail \
            "resolved python is not executable; see" \
            "$(print_relpath "${log_env_rslv}")"
        finish
        exit $?
    fi

    py_cmd=( "${py}" )
fi

if ! \
    "${py_cmd[@]}" - << PY
import sys
raise SystemExit(0 if sys.version_info >= (3, 11) else 1)
PY
then
    record_fail \
        "$("${py_cmd[@]}" --version 2>&1) is older than Python 3.11"
    finish
    exit $?
fi

if \
    PYTHONDONTWRITEBYTECODE=1 \
    run_capture \
        "python version policy" \
        "${log_pol_py}" \
        "${py_cmd[@]}" -m dev.audit.python_version_policy \
            --root "${ROOT_REPO}" \
            --strict
then
    record_pass "repository Python >= 3.11 policy"
else
    record_fail \
        "repository Python >= 3.11 policy; see" \
        "$(print_relpath "${log_pol_py}")"
fi

#  Run lightweight Python syntax checks without writing pycache in-place
while IFS= read -r file; do
    rel="$(print_relpath "${file}")"
    log="${TEST_DIR_LOG}/python_compile/${rel//\//__}.log"

    if \
        PYTHONDONTWRITEBYTECODE=1 \
        PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
        PYTHONPATH="${ROOT_REPO}/src" \
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
    find "${ROOT_REPO}/src/protocol_chipseq_signal_norm" \
        -type f -name '*.py' -print \
        | sort
)

#  Run --help with the project environment so dependency imports are available
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
        PYTHONPATH="${ROOT_REPO}/src" \
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

#  Confirm retired compute_signal.py input/output aliases are not accepted
if \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
    PYTHONPATH="${ROOT_REPO}/src" \
    run_capture \
        "python compute_signal rejects retired input alias" \
        "${log_ret_in}" \
        "${py_cmd[@]}" "${ROOT_REPO}/src/protocol_chipseq_signal_norm/cli/compute_signal.py" \
            "${ret_alias_in}" x \
            --fil_in canonical.bam \
            --fil_out y
then
    record_fail "compute_signal.py unexpectedly accepts retired input alias"
else
    assert_pattern_found \
        "${log_ret_in}" \
        "unrecognized arguments: ${ret_alias_in}" \
        "compute_signal.py rejects retired input alias"
fi

if \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
    PYTHONPATH="${ROOT_REPO}/src" \
    run_capture \
        "python compute_signal rejects retired output alias" \
        "${log_ret_out}" \
        "${py_cmd[@]}" "${ROOT_REPO}/src/protocol_chipseq_signal_norm/cli/compute_signal.py" \
            --fil_in x \
            --fil_out canonical.bdg \
            "${ret_alias_out}" y
then
    record_fail "compute_signal.py unexpectedly accepts retired output alias"
else
    assert_pattern_found \
        "${log_ret_out}" \
        "unrecognized arguments: ${ret_alias_out}" \
        "compute_signal.py rejects retired output alias"
fi

if \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
    PYTHONPATH="${ROOT_REPO}/src" \
    run_capture \
        "python compute_signal rejects retired short input alias" \
        "${log_ret_i}" \
        "${py_cmd[@]}" "${ROOT_REPO}/src/protocol_chipseq_signal_norm/cli/compute_signal.py" \
            "${ret_alias_short}" x \
            --fil_in canonical.bam \
            --fil_out y
then
    record_fail "compute_signal.py unexpectedly accepts retired short input alias"
else
    assert_pattern_found \
        "${log_ret_i}" \
        "unrecognized arguments: ${ret_alias_short}" \
        "compute_signal.py rejects retired short input alias"
fi

#  Confirm the direct spike-in calculator defaults to the ChIP-Rx ratio
if \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
    PYTHONPATH="${ROOT_REPO}/src" \
    run_capture \
        "python calculate-scaling-factor spike default" \
        "${log_def_spk}" \
        "${py_cmd[@]}" -m protocol_chipseq_signal_norm.cli.calculate_scaling_factor_spike \
            --main_ip 100 \
            --spike_ip 10 \
            --main_in 100 \
            --spike_in 20
then
    assert_pattern_found \
        "${log_def_spk}" \
        '^2$' \
        "calculate_scaling_factor_spike.py defaults to chiprx_alpha_ratio"
else
    record_fail \
        "calculate_scaling_factor_spike.py default calculation failed; see" \
        "$(print_relpath "${log_def_spk}")"
fi

finish
