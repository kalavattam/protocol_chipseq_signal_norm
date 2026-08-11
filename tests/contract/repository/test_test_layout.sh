#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_test_layout.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="test organization"

# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


function check_test_scripts() {
    local file=""
    local findings=0

    while IFS= read -r file; do
        if ! grep -q '^set -euo pipefail$' "${file}"; then
            record_fail "test missing safe mode: $(print_relpath "${file}")"
            findings=1
        fi
        if ! grep -q '^TEST_NAME=' "${file}"; then
            record_fail "test missing TEST_NAME: $(print_relpath "${file}")"
            findings=1
        fi
        if ! grep -q 'tests/support/test_helpers\.sh' "${file}"; then
            record_fail \
                "test does not source the shared helper: $(print_relpath "${file}")"
            findings=1
        fi
        if ! grep -q 'print_section "${TEST_NAME}"' "${file}"; then
            record_fail \
                "test missing print_section: $(print_relpath "${file}")"
            findings=1
        fi
    done < <(
        find "${ROOT_REPO}/tests/contract" \
            "${ROOT_REPO}/tests/integration/local" \
            "${ROOT_REPO}/tests/integration/parallel" \
            "${ROOT_REPO}/tests/integration/slurm" \
            -type f -name 'test_*.sh' -print | LC_ALL=C sort
    )

    (( findings == 0 )) && record_pass "shell tests follow the shared contract"
}


function check_boundaries() {
    local finding=""
    local findings=0

    while IFS= read -r finding; do
        [[ -z "${finding}" ]] && continue
        record_fail "test-like executable outside tests/: ${finding}"
        findings=1
    done < <(
        find "${ROOT_REPO}/bin" "${ROOT_REPO}/dev" \
            "${ROOT_REPO}/docs" "${ROOT_REPO}/notebooks" \
            -type f \( -name 'test_*.py' -o -name 'test_*.sh' \
            -o -name '*test_runner*.sh' \) -print \
            | sed "s#^${ROOT_REPO}/##" | LC_ALL=C sort
    )

    [[ -f "${ROOT_REPO}/tests/integration/slurm/coordinator.py" ]] || {
        record_fail "Slurm coordinator is missing from test integration"
        findings=1
    }
    [[ ! -e "${ROOT_REPO}/bin/slurm_bundle.py" ]] || {
        record_fail "test-only Slurm coordinator remains under bin/"
        findings=1
    }
    [[ -f "${ROOT_REPO}/artifacts/README.md" ]] || {
        record_fail "artifact boundary is undocumented"
        findings=1
    }

    (( findings == 0 )) && record_pass "test and artifact boundaries are clean"
}


# Fixture-root conformance is owned by 'test_fixture_paradigm.sh'. The checks
# that used to live here discovered their corpus with 'find ... -name make.sh',
# which could only ever inspect roots that already had a recipe, so the roots
# most needing inspection were the ones it could not see. The production
# boundary below is a different claim about a different corpus and stays.
function check_production_boundary() {
    local findings=0

    if rg -n \
        '^[[:space:]]*(source|\.)[[:space:]].*tests/(fixtures|support)/' \
        "${ROOT_REPO}/bin" "${ROOT_REPO}/lib/bash" \
        > /dev/null ||
        rg -n '^(from|import)[[:space:]]+tests([.]|[[:space:]])' \
            "${ROOT_REPO}/src" > /dev/null
    then
        record_fail "production source depends on test fixtures or support"
        findings=1
    fi

    (( findings == 0 )) && \
        record_pass "production source is independent of tests"
}


print_section "${TEST_NAME}"
check_test_scripts
check_boundaries
check_production_boundary
finish
