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

rec_section "${TEST_NAME}"

#  Run Bash syntax checks over active production shell scripts
while IFS= read -r file; do
    rel="$(rec_relpath "${file}")"
    log="${TEST_DIR_LOG}/shell_syntax/${rel//\//__}.log"

    if \
        run_capture "bash -n ${rel}" "${log}" "${TEST_BASH}" -n "${file}"
    then
        rec_pass "bash -n ${rel}"
    else
        rec_fail "bash -n ${rel}; see $(rec_relpath "${log}")"
    fi
done < <(
    find "${ROOT_REPO}/scripts" \
        -path "${ROOT_REPO}/scripts/blog" -prune -o \
        -type f -name '*.sh' -print \
        | sort
)

finish
