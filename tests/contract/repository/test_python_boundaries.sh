#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_python_boundaries.sh
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

TEST_NAME="Python package boundaries"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


cd "${ROOT_REPO}"

print_section "${TEST_NAME}"

python "${ROOT_REPO}/tests/support/check_import_layers.py" "${ROOT_REPO}"
record_pass "package import directions avoid retired and inverted layers"

finish
