#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_shell_line_length.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="shell line length"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"


#  Check one shell script for long non-heredoc lines
function check_line_length_shell() {
    local file="${1:-}"
    local found=0
    local line=""
    local rel=""

    rel="$(print_relpath "${file}")"

    while IFS= read -r line; do
        found=1
        record_warn "${line}"
    done < <(
        awk -v rel="${rel}" '
            /^[[:space:]]*EOM[[:space:]]*$/ {
                in_hdoc = ""
                next
            }

            /^[[:space:]]*EOF[[:space:]]*$/ {
                in_hdoc = ""
                next
            }

            in_hdoc { next }

            /<<[[:space:]-]*["'\'']?EOM["'\'']?/ {
                in_hdoc = "EOM"
                next
            }

            /<<[[:space:]-]*["'\'']?EOF["'\'']?/ {
                in_hdoc = "EOF"
                next
            }

            length($0) > 80 {
                printf "%s:%d:%d:%s\n", rel, FNR, length($0), $0
                found = 1
            }

            END { exit found }
        ' "${file}"
    )

    if (( found == 0 )); then
        record_pass "no long non-help shell lines in $(print_relpath "${file}")"
    fi
}


#  Collect production, installation, test harness, and fixture scripts
arr_fil=()
while IFS= read -r file; do
    arr_fil+=( "${file}" )
done < <(
    {
        find "${ROOT_REPO}/scripts" \
            -path "${ROOT_REPO}/scripts/blog" -prune -o \
            -type f -name '*.sh' -print

        find "${ROOT_REPO}/install/scripts" \
            -type f -name '*.sh' -print

        find "${ROOT_REPO}/tests" \
            -path "${ROOT_REPO}/tests/outputs" -prune -o \
            -type f -name '*.sh' -print
    } \
        | sort
)


for file in "${arr_fil[@]}"; do
    check_line_length_shell "${file}"
done

record_skip \
    "shell line-length findings are advisory until existing scripts are" \
    "cleaned up"

finish
