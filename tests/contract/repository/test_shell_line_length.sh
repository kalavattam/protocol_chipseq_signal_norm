#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_shell_line_length.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="shell line length"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


#  Check one shell script for long non-heredoc lines
function check_line_length_shell() {
    local file="${1:-}"
    local found=0
    local line=""
    local rel=""
    local strict="false"

    rel="$(print_relpath "${file}")"

    case "${rel}" in
        bin/execute_download_fastqs.sh|\
        bin/submit_download_fastqs.sh|\
        lib/bash/help/help_execute_download_fastqs.sh|\
        lib/bash/help/help_submit_download_fastqs.sh|\
        tests/integration/local/download_fastqs/test_submit_download_fastqs_local.sh|\
        tests/integration/local/download_fastqs/test_execute_download_fastqs_se_local.sh|\
        tests/integration/local/download_fastqs/test_execute_download_fastqs_pe_local.sh|\
        tests/integration/local/download_fastqs/test_execute_download_fastqs_mixed_local.sh|\
        tests/integration/parallel/download_fastqs/test_execute_download_fastqs_parallel.sh)
            strict="true"
            ;;
    esac

    while IFS= read -r line; do
        found=1
        if [[ "${strict}" == "true" ]]; then
            record_fail "${line}"
        else
            record_warn "${line}"
        fi
    done < <(
        awk -v rel="${rel}" -v strict="${strict}" '
            function allowed_long_line(line) {
                if (line ~ /https?:\/\// || line ~ /ftp:\/\//) {
                    return 1
                }
                if (line ~ /[[:xdigit:]]{64}/) {
                    return 1
                }
                if (line ~ /__BASE_URL__/) {
                    return 1
                }
                if (line ~ /\$\047\\t/ || line ~ /seq_[0-9]+=/) {
                    return 1
                }
                if (line ~ /(pat_|pattern|regex|_RE=|_re=)/) {
                    return 1
                }
                if (line ~ /^[[:space:]]*[^[:space:]]*\|[^[:space:]]*\)[[:space:]]*$/) {
                    return 1
                }
                if (index(line, "$") && index(line, "^")) {
                    return 1
                }
                if (index(line, "$") && index(line, "\\t")) {
                    return 1
                }
                if (line ~ /^[[:space:]]*(nam_job|fil_out|log_|hdr_)/) {
                    return 1
                }
                if (line ~ /^[[:space:]]*(tail_|row_|expected_)/) {
                    return 1
                }
                if (line ~ /(dir_log|dir_err|TEST_DIR_LOG)/) {
                    if (line ~ /\.(log|txt)"/) {
                        return 1
                    }
                }
                if (line ~ /heading = .*Usage\|Parameters\|Returns/) {
                    return 1
                }
                return 0
            }

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

            length($0) > 80 && (strict == "true" || ! allowed_long_line($0)) {
                printf "%s:%d:%d:%s\n", rel, FNR, length($0), $0
                found = 1
            }

            END { exit found }
        ' "${file}"
    )

    if (( found == 0 )); then
        record_pass "no long non-help shell lines in ${rel}"
    fi
}


#  Report bounded adjacent case-arm spacing as advisory
function warn_case_arm_spacing() {
    local line=""

    while IFS= read -r line; do
        if [[ "${line}" == SHELL.CASE_ARM.SPACING:*warning* ]]; then
            record_pass "${line}"
        else
            record_warn "${line}"
        fi
    done < <(
        PYTHONDONTWRITEBYTECODE=1 \
            python3 -m dev.audit.help_style \
                --root "${ROOT_REPO}" \
                --case-spacing
    )
}


print_section "${TEST_NAME}"

#  Collect production, installation, test harness, and fixture scripts
arr_fil=()
if (( $# > 0 )); then
    arr_fil=( "$@" )
else
    while IFS= read -r file; do
        arr_fil+=( "${file}" )
    done < <(
        {
            find "${ROOT_REPO}/bin" \
                -type f -name '*.sh' -print

            find "${ROOT_REPO}/lib/bash" \
                -type f -name '*.sh' -print

            find "${ROOT_REPO}/install/scripts" \
                -type f -name '*.sh' -print

            find "${ROOT_REPO}/tests" \
                -path "${ROOT_REPO}/artifacts/tests" -prune -o \
                -type f -name '*.sh' -print
        } \
            | sort
    )
fi


for file in "${arr_fil[@]}"; do
    check_line_length_shell "${file}"
done

warn_case_arm_spacing

record_skip \
    "non-pilot shell line-length findings remain advisory until existing" \
    "scripts are cleaned up"

finish
