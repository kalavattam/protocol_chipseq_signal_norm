#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_shell_syntax.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="shell syntax"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


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


#  Check the POSIX installation bootstrap with sh, not Bash
function check_posix_bootstrap_syntax() {
    local file="${ROOT_REPO}/install/scripts/install_envs_entrypoint.sh"
    local log="${TEST_DIR_LOG}/shell_syntax/install__scripts__install_envs_entrypoint.sh.log"

    if run_capture "sh -n install/scripts/install_envs_entrypoint.sh" "${log}" \
        sh -n "${file}"
    then
        record_pass "POSIX sh -n install/scripts/install_envs_entrypoint.sh"
    else
        record_fail "POSIX sh -n install/scripts/install_envs_entrypoint.sh"
    fi
}


print_section "${TEST_NAME}"

printf 'Authoritative Bash: %s (%s)\n' "${TEST_BASH}" "${TEST_BASH_VERSION}"

if ! "${TEST_BASH}" -c '
    (( BASH_VERSINFO[0] > 4 || (BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] >= 4) ))
'; then
    record_fail "authoritative Bash must satisfy >= 4.4"
else
    record_pass "authoritative Bash satisfies >= 4.4"
fi

#  Collect active production and installation shell scripts
arr_scr=()
while IFS= read -r file; do
    arr_scr+=( "${file}" )
done < <(
    {
        find \
            "${ROOT_REPO}/bin" \
            -type f -name '*.sh' -print
        find \
            "${ROOT_REPO}/lib/bash" \
            -type f -name '*.sh' -print
        find \
            "${ROOT_REPO}/install/scripts" \
            -type f -name '*.sh' \
            ! -name 'install_envs_entrypoint.sh' -print
    } \
        | sort
)


#  Collect test harness and fixture-generation shell scripts
arr_test=()
while IFS= read -r file; do
    arr_test+=( "${file}" )
done < <(
    find "${ROOT_REPO}/tests" \
        -path "${ROOT_REPO}/artifacts/tests" -prune -o \
        -type f -name '*.sh' -print \
        | sort
)


check_group_shell_syntax "Production shell scripts" "${arr_scr[@]}"
check_group_shell_syntax "Test harness shell scripts" "${arr_test[@]}"
check_posix_bootstrap_syntax

finish "$@"
