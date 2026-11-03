#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_markdown_policy.sh
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

TEST_NAME="Markdown deterministic policy"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


readonly OUT_DIR="${TEST_DIR_OUT}/current_inventory"
readonly REPORT="${OUT_DIR}/markdown_findings.json"
readonly PREVIEW="${OUT_DIR}/markdown_format_preview.txt"

mkdir -p "${OUT_DIR}"
cd "${ROOT_REPO}"

print_section "${TEST_NAME}"

python_cmd=( "${TEST_PYTHON}" )

"${python_cmd[@]}" -m pytest -q tests/unit/dev_audit/test_markdown_policy.py
"${python_cmd[@]}" -m dev.tools.markdown_format \
    --root . > "${PREVIEW}"

status=0
"${python_cmd[@]}" -m dev.audit.markdown_policy --root . --json \
    > "${REPORT}" \
    || status=$?
(( status == 0 ))

"${python_cmd[@]}" -m json.tool "${REPORT}" > /dev/null

if \
    ! "${python_cmd[@]}" \
        "${ROOT_REPO}/tests/support/check_findings.py" \
        "${REPORT}"
then
    record_fail "Markdown checker reported deterministic findings"
fi

printf 'Markdown deterministic policy: status=%d report=%s\n' \
    "${status}" "${REPORT#"${ROOT_REPO}/"}"
printf 'Markdown deterministic formatter preview: clean report=%s\n' \
    "${PREVIEW#"${ROOT_REPO}/"}"

record_pass "Markdown deterministic checks and fixtures passed"

finish
