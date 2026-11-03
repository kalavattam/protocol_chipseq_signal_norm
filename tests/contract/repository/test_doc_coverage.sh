#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_doc_coverage.sh
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

TEST_NAME="documentation coverage"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


#  Run the repository documentation-coverage scanner
function scan_doc_coverage() {
    local found=0
    local line sev file lineno msg log

    log="${TEST_DIR_LOG}/doc_coverage.tsv"

    if ! (
        cd "${ROOT_REPO}" && python3 dev/audit/doc_coverage.py > "${log}"
    ); then
        record_fail "documentation coverage scanner crashed"
        return
    fi

    while IFS=$'\t' read -r sev file lineno msg; do
        [[ -n "${sev}" ]] || continue

        found=1
        case "${sev}" in
            FAIL)
                record_fail "${file}:${lineno}:${msg}"
                ;;
            WARN)
                record_warn "${file}:${lineno}:${msg}"
                ;;
            *)
                record_fail "unknown doc coverage severity '${sev}'"
                ;;
        esac
    done < "${log}"

    if (( found == 0 )); then
        record_pass "documentation coverage queue is empty"
    fi
}


print_section "${TEST_NAME}"

scan_doc_coverage

finish
