#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_parameter_docs_consistency.sh
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

TEST_NAME="parameter description consistency"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


log_style="${TEST_DIR_LOG}/parameter_docs_consistency/style.log"
mkdir -p "$(dirname "${log_style}")"

print_section "${TEST_NAME}"

aud_pdc="${ROOT_REPO}/dev/audit/parameter_docs_consistency.py"

if \
    PYTHONDONTWRITEBYTECODE=1 \
    ROOT_REPO="${ROOT_REPO}" \
    python3 "${aud_pdc}" > "${log_style}"
then
    while IFS=: read -r sev file line msg; do
        [[ -n "${sev}" ]] || continue
        case "${sev}" in
            FAIL)
                record_fail "${file}:${line}:${msg}"
                ;;
            WARN)
                record_warn "${file}:${line}:${msg}"
                ;;
            *)
                record_fail \
                    "unrecognized parameter-doc finding:" \
                    "${sev}:${file}:${line}:${msg}"
                ;;
        esac
    done < "${log_style}"
else
    record_fail "parameter description consistency scanner crashed"
fi

if [[ ! -s "${log_style}" ]]; then
    record_pass "no parameter description consistency findings"
fi

finish
