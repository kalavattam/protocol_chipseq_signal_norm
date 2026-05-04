#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_startup_sources.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="startup and helper sourcing"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

#  Source each function script in isolation to catch immediate load failures
while IFS= read -r file; do
    rel="$(rec_relpath "${file}")"
    log="${TEST_DIR_LOG}/source/${rel//\//__}.log"

    # shellcheck disable=SC2016
    if \
        run_capture \
            "source ${rel}" "${log}" "${TEST_BASH}" -c 'source "$1"' _ \
            "${file}"
    then
        rec_pass "source ${rel}"
    else
        rec_fail "source ${rel}; see $(rec_relpath "${log}")"
    fi
done < <(
    find "${ROOT_REPO}/scripts/functions" -type f -name '*.sh' -print | sort
)

finish
