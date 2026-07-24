#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_ruff_inventory.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="Ruff approved core and deferred formatter inventory"
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"
readonly OUT_DIR="${TEST_DIR_OUT}/current_inventory"
mkdir -p "${OUT_DIR}"
cd "${ROOT_REPO}"
print_section "${TEST_NAME}"

if command -v mamba >/dev/null 2>&1; then
    manager=mamba
else
    manager=conda
fi

set +e
"${manager}" run --no-capture-output -n env_protocol \
    ruff check --no-fix src tests dev \
    > "${OUT_DIR}/ruff_check.txt" 2>&1
check_status=$?
"${manager}" run --no-capture-output -n env_protocol \
    ruff format --check src tests dev \
    > "${OUT_DIR}/ruff_format_check.txt" 2>&1
format_status=$?
set -e

(( check_status == 0 ))
(( format_status <= 1 ))
printf 'Ruff approved core: check=%d; deferred format inventory=%d\n' \
    "${check_status}" "${format_status}"
record_pass "Ruff core passed and formatter inventory completed"
finish
