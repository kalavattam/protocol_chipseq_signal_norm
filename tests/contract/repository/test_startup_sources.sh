#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_startup_sources.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in design, development,
# and documentation, with all output reviewed, edited, and approved by the
# author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="startup and helper sourcing"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


print_section "${TEST_NAME}"

#  Source each function script in isolation to catch immediate load failures
while IFS= read -r file; do
    rel="$(print_relpath "${file}")"
    log="${TEST_DIR_LOG}/source/${rel//\//__}.log"

    # shellcheck disable=SC2016
    if \
        run_capture \
            "source ${rel}" \
            "${log}" \
            "${TEST_BASH}" -c 'source "$1"' _ "${file}"
    then
        record_pass "source ${rel}"
    else
        record_fail "source ${rel}; see $(print_relpath "${log}")"
    fi
done < <(
    find "${ROOT_REPO}/lib/bash" -type f -name '*.sh' -print | sort
)

finish
