#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_shell_syntax.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="shell syntax"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"


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


#  Collect active production shell scripts
arr_scr=()
while IFS= read -r file; do
    arr_scr+=( "${file}" )
done < <(
    find "${ROOT_REPO}/scripts" \
        -path "${ROOT_REPO}/scripts/blog" -prune -o \
        -type f -name '*.sh' -print \
        | sort
)


#  Collect smoke-test harness scripts
arr_test=(
    "${ROOT_REPO}/tests/scripts/run_smoke_tests.sh"
    "${ROOT_REPO}/tests/scripts/lib/test_helpers.sh"
)

while IFS= read -r file; do
    arr_test+=( "${file}" )
done < <(
    find "${ROOT_REPO}/tests/scripts/smoke" \
        -type f -name '*.sh' -print \
        | sort
)

arr_test+=(
    "${ROOT_REPO}/tests/compute_signal/scripts/make_fixtures.sh"
)


check_group_shell_syntax "Production shell scripts" "${arr_scr[@]}"
check_group_shell_syntax "Test harness shell scripts" "${arr_test[@]}"

finish
