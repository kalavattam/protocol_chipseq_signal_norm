#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_standards_registry.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="standards registry"
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"
cd "${ROOT_REPO}"
print_section "${TEST_NAME}"

python_cmd=( "${TEST_PYTHON}" )

"${python_cmd[@]}" -m pytest -q tests/unit/dev_audit/test_standards_registry.py
"${python_cmd[@]}" -m dev.audit.standards_registry --root . --report-only
record_pass "registry checker and fixtures completed"
finish
