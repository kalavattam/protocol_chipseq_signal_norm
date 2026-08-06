#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_markdown_policy.sh
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

TEST_NAME="Markdown deterministic policy"
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
set +e
"${python_cmd[@]}" -m dev.audit.markdown_policy --root . --json > "${REPORT}"
status=$?
set -e
(( status == 0 ))
"${python_cmd[@]}" -m json.tool "${REPORT}" >/dev/null
if "${python_cmd[@]}" - "${REPORT}" << PY
import json
import sys
from pathlib import Path

findings = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
deterministic = [
    item for item in findings
    if item["classification"] == "deterministic"
]
raise SystemExit(1 if deterministic else 0)
PY
then
    :
else
    record_fail "Markdown checker reported deterministic findings"
fi
printf 'Markdown deterministic policy: status=%d report=%s\n' \
    "${status}" "${REPORT#"${ROOT_REPO}/"}"
printf 'Markdown deterministic formatter preview: clean report=%s\n' \
    "${PREVIEW#"${ROOT_REPO}/"}"
record_pass "Markdown deterministic checks and fixtures passed"
finish
