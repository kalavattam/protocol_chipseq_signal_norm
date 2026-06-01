#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_clean_test_outputs.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="clean test outputs dry-run"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

scr_cln="${ROOT_REPO}/tests/scripts/clean_test_outputs.sh"
dir_log="${TEST_DIR_LOG}/clean_test_outputs"
mkdir -p "${dir_log}"


#  Assert one explicitly scoped cleanup dry-run exits successfully
function assert_clean_dry() {
    local nam_cas="${1:-cleanup}"
    local log="${2:-}"

    shift 2

    if \
        run_capture \
            "clean_test_outputs ${nam_cas} dry-run" \
            "${log}" \
            "${TEST_BASH}" "${scr_cln}" "$@"
    then
        record_pass "clean_test_outputs.sh ${nam_cas} dry-run exits 0"
        assert_pattern_found \
            "${log}" \
            "dry-run mode: reporting ignored generated artifacts." \
            "clean_test_outputs.sh ${nam_cas} uses dry-run mode"
    else
        record_fail \
            "clean_test_outputs.sh ${nam_cas} dry-run failed; see" \
            "$(print_relpath "${log}")"
    fi
}


#  Check concise help output
log="${dir_log}/help.log"
if \
    run_capture \
        "clean_test_outputs help" \
        "${log}" \
        "${TEST_BASH}" "${scr_cln}" --help
then
    record_pass "clean_test_outputs.sh --help exits 0"
    assert_pattern_found \
        "${log}" \
        '^Usage:' \
        "clean_test_outputs.sh --help prints Usage"
else
    record_fail \
        "clean_test_outputs.sh --help failed; see $(print_relpath "${log}")"
fi


#  Check selector validation and unknown-option handling
log="${dir_log}/no_selector.log"
if \
    run_capture \
        "clean_test_outputs no selector" \
        "${log}" \
        "${TEST_BASH}" "${scr_cln}"
then
    record_fail "clean_test_outputs.sh without selector unexpectedly succeeded"
else
    assert_pattern_found \
        "${log}" \
        "one of '--fixtures', '--outputs', or '--all' must be specified" \
        "clean_test_outputs.sh without selector emits useful error"
fi

log="${dir_log}/unknown_option.log"
if \
    run_capture \
        "clean_test_outputs unknown option" \
        "${log}" \
        "${TEST_BASH}" "${scr_cln}" --bogus
then
    record_fail "clean_test_outputs.sh unknown option unexpectedly succeeded"
else
    assert_pattern_found \
        "${log}" \
        "unknown option/parameter passed: '--bogus'" \
        "clean_test_outputs.sh unknown option emits useful error"
fi


#  Check explicitly scoped cleanup plans without deleting files
assert_clean_dry \
    "all" \
    "${dir_log}/all_dry_run.log" \
    --all \
    --dry_run

assert_clean_dry \
    "fixtures" \
    "${dir_log}/fixtures_dry_run.log" \
    --fixtures \
    --dry_run

assert_clean_dry \
    "outputs" \
    "${dir_log}/outputs_dry_run.log" \
    --outputs \
    --dry_run

finish
