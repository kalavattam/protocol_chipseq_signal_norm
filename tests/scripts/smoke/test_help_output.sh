#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_help_output.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="help output"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

#  Check short and detailed help rendering for top-level shell scripts
while IFS= read -r file; do
    rel="$(rec_relpath "${file}")"
    log="${TEST_DIR_LOG}/help/${rel//\//__}.help.log"

    if \
        run_capture \
            "help ${rel}" "${log}" \
            "${TEST_BASH}" "${file}" --help
    then
        rec_pass "--help exits 0 for ${rel}"

        warn_grep_help \
            "${log}" \
            '^Usage:' \
            "${rel} help has Usage:"

        warn_grep_help \
            "${log}" \
            '^Description:' \
            "${rel} help has Description:"

        warn_grep_help \
            "${log}" \
            '^\(Arguments\|Positional arguments\|Keyword arguments\):' \
            "${rel} help has an argument section"

        warn_grep_help \
            "${log}" \
            '^Notes:' \
            "${rel} help has Notes:"
    else
        rec_fail "--help failed for ${rel}; see $(rec_relpath "${log}")"
        continue
    fi

    if \
        grep -q -- '--details\|--detail' "${log}"
    then
        dlog="${TEST_DIR_LOG}/help/${rel//\//__}.details.log"
        if \
            run_capture \
                "details ${rel}" "${dlog}" \
                "${TEST_BASH}" "${file}" --details
        then
            rec_pass "--details exits 0 for ${rel}"

            warn_grep_help \
                "${dlog}" \
                '^Usage:' \
                "${rel} details has Usage:"

            warn_grep_help \
                "${dlog}" \
                '^Description:' \
                "${rel} details has Description:"

            warn_grep_help \
                "${dlog}" \
                '^\(Arguments\|Positional arguments\|Keyword arguments\):' \
                "${rel} details has an argument section"

            warn_grep_help \
                "${dlog}" \
                '^Dependencies:' \
                "${rel} details has Dependencies:"

            warn_grep_help \
                "${dlog}" \
                '^Notes:' \
                "${rel} details has Notes:"

            warn_grep_help \
                "${dlog}" \
                '^Examples:' \
                "${rel} details has Examples:"
        else
            rec_warn \
                "--details advertised but failed for ${rel}; see" \
                "$(rec_relpath "${dlog}")"
        fi
    else
        rec_skip "no detailed help flag advertised for ${rel}"
    fi
done < <(
    find "${ROOT_REPO}/scripts" -maxdepth 1 -type f -name '*.sh' -print | sort
)

finish
