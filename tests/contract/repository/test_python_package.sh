#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_python_package.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5, Opus 5.5).
#
# Distributed under the MIT license.


set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

TEST_NAME="python package"

# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


print_section "${TEST_NAME}"

if check_cmd_exists mamba; then
    manager=mamba
elif check_cmd_exists conda; then
    manager=conda
else
    record_fail "neither mamba nor conda is available"
    finish
    exit $?
fi

env_cmd=( "${manager}" run --no-capture-output -n env_protocol )
package_artifacts="${TEST_DIR_OUT}/package"
mkdir -p "${package_artifacts}/wheel"
package_source="$(mktemp -d "${package_artifacts}/source.XXXXXX")"
cp -R \
    "${ROOT_REPO}/LICENSE" \
    "${ROOT_REPO}/README.md" \
    "${ROOT_REPO}/pyproject.toml" \
    "${ROOT_REPO}/setup.cfg" \
    "${ROOT_REPO}/src" \
    "${package_source}/"  # fixture-perturbation: the build writes beside it

chk_pkg="${ROOT_REPO}/tests/support/check_python_package.py"

if \
    "${env_cmd[@]}" python "${chk_pkg}" packages "${ROOT_REPO}"
then
    record_pass "setuptools discovers only the intended src package"
else
    record_fail "setuptools package discovery differs from the intended package"
fi

if (
    cd /private/tmp
    env -u PYTHONPATH "${env_cmd[@]}" python "${chk_pkg}" import-path
)
then
    record_pass "editable package imports outside the checkout root"
else
    record_fail "editable package import outside the checkout root failed"
fi

if \
    "${env_cmd[@]}" python "${chk_pkg}" entry-points
then
    record_pass "console and explicit-module forms share each typed main"
else
    record_fail "console-entrypoint registration or main identity failed"
fi

commands=(
    add_coeffs_namespaced
    calculate_scaling_factor_siqchip
    calculate_scaling_factor_spike
    compute_input_floor
    compute_pseudo
    compute_pseudo_deeptools
    compute_signal
    compute_signal_ratio
    merge_bins_bdg
    parse_metadata_siqchip
    relativize_scaling_factors
    sum_bdg
)
for command_name in "${commands[@]}"; do
    if "${env_cmd[@]}" sh -c \
        'command -v "$1" >/dev/null && "$1" --help >/dev/null 2>&1' \
        _ "${command_name}"
    then
        record_pass "console command resolves and supports --help: ${command_name}"
    else
        record_fail "console command resolution/help failed: ${command_name}"
    fi
done

if "${env_cmd[@]}" python -m pip wheel \
    --no-deps \
    --no-build-isolation \
    --wheel-dir "${package_artifacts}/wheel" \
    "${package_source}" >/dev/null
then
    record_pass "package wheel builds under the selected artifact root"
else
    record_fail "package wheel build failed"
fi

if rg -n \
    '(^|[^A-Za-z0-9_])(from|import)[[:space:]]+scripts([./]|$)|python[[:space:]]+-m[[:space:]]+scripts\.' \
    "${ROOT_REPO}" \
    -g '!artifacts/**' \
    -g '!tmp/**'
then
    record_fail "retired scripts package invocation or import remains"
else
    record_pass "no scripts package invocation or import remains"
fi

finish
