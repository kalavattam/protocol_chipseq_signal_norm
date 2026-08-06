#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_python_package.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
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
    "${package_source}/"

if "${env_cmd[@]}" python - << PY
from setuptools import find_packages

expected = {
    "protocol_chipseq_signal_norm",
    "protocol_chipseq_signal_norm.cli",
    "protocol_chipseq_signal_norm.utilities",
}
actual = set(find_packages(where="src"))
raise SystemExit(actual != expected)
PY
then
    record_pass "setuptools discovers only the intended src package"
else
    record_fail "setuptools package discovery differs from the intended package"
fi

if (
    cd /private/tmp
    env -u PYTHONPATH "${env_cmd[@]}" python - << PY
from pathlib import Path

import protocol_chipseq_signal_norm

path = Path(protocol_chipseq_signal_norm.__file__).resolve()
raise SystemExit(path.name != "__init__.py" or "src" not in path.parts)
PY
)
then
    record_pass "editable package imports outside the checkout root"
else
    record_fail "editable package import outside the checkout root failed"
fi

if "${env_cmd[@]}" python - << PY
from importlib import import_module
from importlib.metadata import entry_points

expected = {
    "add_coeffs_namespaced",
    "calculate_scaling_factor_siqchip",
    "calculate_scaling_factor_spike",
    "compute_input_floor",
    "compute_pseudo",
    "compute_signal",
    "compute_signal_ratio",
    "merge_bins_bdg",
    "parse_metadata_siqchip",
    "relativize_scaling_factors",
    "sum_bdg",
}
registered = {
    item.name: item
    for item in entry_points(group="console_scripts")
    if item.name in expected
}
if set(registered) != expected:
    raise SystemExit(1)
for name, entry in registered.items():
    module = import_module(f"protocol_chipseq_signal_norm.cli.{name}")
    if entry.load() is not module.main:
        raise SystemExit(1)
PY
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
