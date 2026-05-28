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

print_section "${TEST_NAME}"

#  Create small temporary inputs for dry-run smoke checks
tmp="${TEST_DIR_TMP}/dry_run"
rm -rf "${tmp}"
mkdir -p "${tmp}/in" "${tmp}/out" "${tmp}/logs"
printf 'smoke\n' > "${tmp}/in/sample.txt"

log="${TEST_DIR_LOG}/dry_run/symlink_files.log"
if \
    run_capture \
        "symlink dry-run" \
        "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/symlink_files.sh" \
            --dry_run \
            --csv_infile "${tmp}/in/sample.txt" \
            --dir_out "${tmp}/out"
then
    record_pass "symlink_files.sh minimal --dry_run"
else
    record_fail \
        "symlink_files.sh minimal --dry_run; see $(print_relpath "${log}")"
fi

#  Check that write_header.sh dry-run reports work without writing output
log="${TEST_DIR_LOG}/dry_run/write_header.log"
if \
    run_capture \
        "write_header dry-run" \
        "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/write_header.sh" \
            --dry_run \
            --mode siq \
            --fil_out "${tmp}/out/header.tsv"
then
    if [[ -e "${tmp}/out/header.tsv" ]]; then
        record_fail \
            "write_header.sh --dry_run created an output file; see" \
            "$(print_relpath "${log}")"
    else
        assert_pattern_found \
            "${log}" \
            "Dry run: would create" \
            "write_header.sh minimal --dry_run"
    fi
else
    record_fail \
        "write_header.sh minimal --dry_run; see $(print_relpath "${log}")"
fi

#  Check that selected wrappers fail clearly when required arguments are absent
log="${TEST_DIR_LOG}/expected_fail/symlink_missing_args.log"
if \
    run_capture \
        "symlink missing args" \
        "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/symlink_files.sh" --dry_run
then
    record_fail "symlink_files.sh missing required args unexpectedly succeeded"
else
    assert_pattern_found \
        "${log}" \
        'error' \
        "symlink_files.sh missing args emits useful error"
fi

log="${TEST_DIR_LOG}/expected_fail/find_missing_args.log"
if \
    run_capture \
        "find missing args" \
        "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/find_files.sh" --pattern '*.sh'
then
    record_fail "find_files.sh missing required args unexpectedly succeeded"
else
    assert_pattern_found \
        "${log}" \
        'error' \
        "find_files.sh missing args emits useful error"
fi

#TODO: write tests for these after codifying and implementing updates
record_skip \
    "execute_* dry-run wrappers require realistic files and/or environment" \
    "checks; covered later by integration fixtures"

record_skip \
    "submit_* dry-run wrappers may activate Conda/Mamba or require Slurm" \
    "context; covered later by integration fixtures"

finish
