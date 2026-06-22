#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_shell_syntax.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="shell syntax"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"


#  Run Bash syntax checks over one shell-script group
function check_group_shell_syntax() {
    local nam_grp="${1:-shell scripts}"
    local file=""
    local rel=""
    local log=""

    shift

    printf '\n%s:\n' "${nam_grp}"

    for file in "$@"; do
        rel="$(print_relpath "${file}")"
        log="${TEST_DIR_LOG}/shell_syntax/${rel//\//__}.log"

        if \
            run_capture \
                "bash -n ${rel}" \
                "${log}" \
                "${TEST_BASH}" -n "${file}"
        then
            record_pass "bash -n ${rel}"
        else
            record_fail "bash -n ${rel}; see $(print_relpath "${log}")"
        fi
    done
}


print_section "${TEST_NAME}"

#  Collect active production and installation shell scripts
arr_scr=()
while IFS= read -r file; do
    arr_scr+=( "${file}" )
done < <(
    {
        find \
            "${ROOT_REPO}/scripts" \
            -path "${ROOT_REPO}/scripts/blog" \
            -prune -o -type f -name '*.sh' -print
        find \
            "${ROOT_REPO}/install/scripts" \
            -type f -name '*.sh' -print
    } \
        | sort
)


#  Collect test harness and fixture-generation shell scripts
arr_test=()
while IFS= read -r file; do
    arr_test+=( "${file}" )
done < <(
    find "${ROOT_REPO}/tests" \
        -path "${ROOT_REPO}/tests/outputs" -prune -o \
        -type f -name '*.sh' -print \
        | sort
)


check_group_shell_syntax "Production shell scripts" "${arr_scr[@]}"
check_group_shell_syntax "Test harness shell scripts" "${arr_test[@]}"

finish
