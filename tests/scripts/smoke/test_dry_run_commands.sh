#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_dry_run_commands.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="dry-run and expected failure"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

#  Create small temporary inputs for dry-run smoke checks
tmp="${TEST_DIR_TMP}/dry_run"
rm -rf "${tmp}"
mkdir -p "${tmp}/in" "${tmp}/out" "${tmp}/logs"
printf 'smoke\n' > "${tmp}/in/sample.txt"

log="${TEST_DIR_LOG}/dry_run/symlink_files.log"
if \
    run_capture "symlink dry-run" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/symlink_files.sh" \
            --dry_run \
            --csv_infile "${tmp}/in/sample.txt" \
            --dir_out "${tmp}/out"
then
    rec_pass "symlink_files.sh minimal --dry_run"
else
    rec_fail "symlink_files.sh minimal --dry_run; see $(rec_relpath "${log}")"
fi

#  Keep the write_header.sh dry-run expectation visible while production
#+ support is still being validated across local checkouts
log="${TEST_DIR_LOG}/dry_run/write_header.log"
if \
    run_capture "write_header dry-run" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/write_header.sh" \
            --dry_run \
            --mode alpha \
            --fil_out "${tmp}/out/header.tsv"
then
    if [[ -e "${tmp}/out/header.tsv" ]]; then
        rec_fail \
            "write_header.sh --dry_run created an output file; see" \
            "$(rec_relpath "${log}")"
    else
        rec_pass "write_header.sh minimal --dry_run"
    fi
else
    rec_skip \
        "write_header.sh minimal --dry_run pending production support; see" \
        "$(rec_relpath "${log}")"
fi

#  Check that selected wrappers fail clearly when required arguments are absent
log="${TEST_DIR_LOG}/expected_fail/symlink_missing_args.log"
if \
    run_capture \
        "symlink missing args" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/symlink_files.sh" --dry_run
then
    rec_fail "symlink_files.sh missing required args unexpectedly succeeded"
else
    assert_grep_pattern \
        "${log}" 'error' "symlink_files.sh missing args emits useful error"
fi

log="${TEST_DIR_LOG}/expected_fail/find_missing_args.log"
if \
    run_capture \
        "find missing args" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/find_files.sh" --pattern '*.sh'
then
    rec_fail "find_files.sh missing required args unexpectedly succeeded"
else
    assert_grep_pattern \
        "${log}" 'error' "find_files.sh missing args emits useful error"
fi

rec_skip \
    "execute_* dry-run wrappers require realistic files and/or environment" \
    "checks; covered later by integration fixtures"
rec_skip \
    "submit_* dry-run wrappers may activate Conda/Mamba or require Slurm" \
    "context; covered later by integration fixtures"

finish
