#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_python_state_isolation.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="Python test-state isolation"
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


function maintained_python_state() {
    find \
        "${ROOT_REPO}/src" \
        "${ROOT_REPO}/dev" \
        "${ROOT_REPO}/tests" \
        \( -type d -name __pycache__ -o \
            -type d -name .pytest_cache -o \
            -type f -name '*.pyc' \) \
        -print
}


print_section "${TEST_NAME}"

if [[ "${PYTHONDONTWRITEBYTECODE:-}" == "1" ]]; then
    record_pass "bytecode writes are disabled for test children"
else
    record_fail "PYTHONDONTWRITEBYTECODE is not enabled"
fi

canonical_root="${ROOT_REPO}/artifacts/tests"
case "${TEST_ARTIFACT_ROOT:-}" in
    "${canonical_root}")
        record_pass "canonical repository artifact root is selected"
        ;;

    ""|"${ROOT_REPO}"|"${ROOT_REPO}"/*)
        record_fail "artifact root is an arbitrary maintained-repository path"
        ;;

    /*)
        record_pass "approved absolute external artifact root is selected"
        ;;

    *)
        record_fail "artifact root is not an approved absolute path"
        ;;
esac

if [[ "${PYTHONPYCACHEPREFIX:-}" == "${TEST_ARTIFACT_ROOT}/pycache" ]]; then
    record_pass "Python cache prefix is rooted in the exact artifact root"
else
    record_fail "Python cache prefix is not rooted in the exact artifact root"
fi

if [[ "${PYTEST_ADDOPTS:-}" == \
    *"cache_dir=${TEST_ARTIFACT_ROOT}/pytest_cache"* ]]
then
    record_pass "nested pytest cache is rooted in the exact artifact root"
else
    record_fail "nested pytest cache is not rooted in the exact artifact root"
fi

mapfile -t arr_generated_state < <(maintained_python_state)
if (( ${#arr_generated_state[@]} == 0 )); then
    record_pass "maintained roots contain no generated Python test state"
else
    printf 'Generated Python state:\n' >&2
    printf '  %s\n' "${arr_generated_state[@]}" >&2
    record_fail "maintained roots contain generated Python test state"
fi

finish
