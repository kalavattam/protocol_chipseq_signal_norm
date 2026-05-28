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

print_section "${TEST_NAME}"

#  Check short and detailed help rendering for top-level shell scripts
while IFS= read -r file; do
    rel="$(print_relpath "${file}")"
    log="${TEST_DIR_LOG}/help/${rel//\//__}.help.log"

    if \
        run_capture \
            "help ${rel}" \
            "${log}" \
            "${TEST_BASH}" "${file}" --help
    then
        record_pass "--help exits 0 for ${rel}"

        warn_help_pattern_missing \
            "${log}" \
            '^Usage:' \
            "${rel} help has Usage:"

        warn_help_pattern_missing \
            "${log}" \
            '^Description:' \
            "${rel} help has Description:"

        warn_help_pattern_missing \
            "${log}" \
            '^\(Arguments\|Positional arguments\|Keyword arguments\):' \
            "${rel} help has an argument section"

        warn_help_pattern_missing \
            "${log}" \
            '^Notes:' \
            "${rel} help has Notes:"
    else
        record_fail "--help failed for ${rel}; see $(print_relpath "${log}")"
        continue
    fi

    if \
        grep -q -- '--details\|--detail' "${log}"
    then
        dlog="${TEST_DIR_LOG}/help/${rel//\//__}.details.log"
        if \
            run_capture \
                "details ${rel}" \
                "${dlog}" \
                "${TEST_BASH}" "${file}" --details
        then
            record_pass "--details exits 0 for ${rel}"

            warn_help_pattern_missing \
                "${dlog}" \
                '^Usage:' \
                "${rel} details has Usage:"

            warn_help_pattern_missing \
                "${dlog}" \
                '^Description:' \
                "${rel} details has Description:"

            warn_help_pattern_missing \
                "${dlog}" \
                '^\(Arguments\|Positional arguments\|Keyword arguments\):' \
                "${rel} details has an argument section"

            warn_help_pattern_missing \
                "${dlog}" \
                '^Dependencies:' \
                "${rel} details has Dependencies:"

            warn_help_pattern_missing \
                "${dlog}" \
                '^Notes:' \
                "${rel} details has Notes:"

            warn_help_pattern_missing \
                "${dlog}" \
                '^Examples:' \
                "${rel} details has Examples:"
        else
            record_warn \
                "--details advertised but failed for ${rel}; see" \
                "$(print_relpath "${dlog}")"
        fi
    else
        record_skip "no detailed help flag advertised for ${rel}"
    fi
done < <(
    find "${ROOT_REPO}/scripts" -maxdepth 1 -type f -name '*.sh' -print | sort
)

finish
