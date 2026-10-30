#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_find_files_opt_applicability.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="find-files option applicability"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define the search directory. Per 'HELP.PARAMETER.APPLICABILITY', '--shw_str'
# and '--fastqs' are ignored with a warning when no one-line string is printed.
tmp="${TEST_DIR_TMP}/find_files_opt_applicability"
dir_fq="${ROOT_REPO}/tests/fixtures/utilities/fastq/find"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}"

require_files_nonempty \
    "${dir_fq}/a_R1.fastq.gz" \
    "${dir_fq}/a_R2.fastq.gz" || {
    finish
    exit $?
}


# Run 'find_files' from outside the search directory, as it requires.
function run_find() {
    (
        cd "${tmp}"
        "${TEST_BASH}" "${ROOT_REPO}/bin/find_files.sh" \
            --dir_fnd "${dir_fq}" \
            --pattern "*.fastq.gz" \
            "$@" 2>&1
    )
}


# Each row is label|arguments|option warned about, or empty when none is.
rows=(
    "--shw_str alone|--shw_str|--shw_str"
    "--shw_str with --chk_con|--chk_con --shw_str|--shw_str"
    "--fastqs with --chk_con|--chk_con --fastqs|--fastqs"
    "--fastqs with --chk_exc|--chk_exc --fastqs|--fastqs"
    "--fastqs alone|--fastqs|"
    "--shw_str with --chk_exc|--chk_exc --shw_str|"
    "both with --chk_exc|--chk_exc --shw_str --fastqs|"
)

for row in "${rows[@]}"; do
    IFS='|' read -r lbl arg opt <<< "${row}"
    read -r -a args <<< "${arg}"

    rc=0
    out="$(run_find "${args[@]}")" || rc=$?

    if [[ -n "${opt}" ]]; then
        if [[
            "${rc}" -eq 0
                && "${out}" == *"warning("*"'${opt}' has no effect"*
        ]]; then
            record_pass "${lbl} warns and still runs"
        else
            record_fail "${lbl} did not warn or failed (exit ${rc})"
        fi
    elif [[ "${rc}" -eq 0 && "${out}" != *"has no effect"* ]]; then
        record_pass "${lbl} runs with no warning"
    else
        record_fail "${lbl} warned or failed (exit ${rc})"
    fi
done


finish
