#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_python_style.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5, Opus 5.5).
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="python style"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


log_style="${TEST_DIR_LOG}/python_style/style.log"
mkdir -p "$(dirname "${log_style}")"

print_section "${TEST_NAME}"

# Resolve the import root before the prefix, so the prefix does not expand a
# variable that it also assigns.
dir_src="${ROOT_REPO}/src"
aud_sty="${ROOT_REPO}/dev/audit/python_style.py"

if \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPATH="${dir_src}" \
    ROOT_REPO="${ROOT_REPO}" \
    python3 "${aud_sty}" > "${log_style}"
then
    :
else
    record_fail "python style scanner crashed; see $(print_relpath "${log_style}")"
fi

while IFS=: read -r sev rel line msg; do
    case "${sev}" in
        PASS) record_pass "${rel}" ;;
        WARN) record_warn "${rel}:${line}:${msg}" ;;
        FAIL) record_fail "${rel}:${line}:${msg}" ;;
    esac
done < "${log_style}"

finish
