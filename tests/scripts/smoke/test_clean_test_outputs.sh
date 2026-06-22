#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_clean_test_outputs.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="clean test outputs dry-run"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"


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


print_section "${TEST_NAME}"

scr_cln="${ROOT_REPO}/tests/scripts/clean_test_outputs.sh"
dir_log="${TEST_DIR_LOG}/clean_test_outputs"

log_help="${dir_log}/help.log"
log_no_selector="${dir_log}/no_selector.log"
log_opt_unknown="${dir_log}/unknown_option.log"
log_dry_all="${dir_log}/all_dry_run.log"
log_dry_fix="${dir_log}/fixtures_dry_run.log"
log_dry_out="${dir_log}/outputs_dry_run.log"

mkdir -p "${dir_log}"

#  Check concise help output
if \
    run_capture \
        "clean_test_outputs help" \
        "${log_help}" \
        "${TEST_BASH}" "${scr_cln}" --help
then
    record_pass "clean_test_outputs.sh --help exits 0"
    assert_pattern_found \
        "${log_help}" \
        '^Usage:' \
        "clean_test_outputs.sh --help prints Usage"
else
    record_fail \
        "clean_test_outputs.sh --help failed; see" \
        "$(print_relpath "${log_help}")"
fi


#  Check selector validation and unknown-option handling
if \
    run_capture \
        "clean_test_outputs no selector" \
        "${log_no_selector}" \
        "${TEST_BASH}" "${scr_cln}"
then
    record_fail "clean_test_outputs.sh without selector unexpectedly succeeded"
else
    assert_pattern_found \
        "${log_no_selector}" \
        "one of '--fixtures', '--outputs', or '--all' must be specified" \
        "clean_test_outputs.sh without selector emits useful error"
fi

if \
    run_capture \
        "clean_test_outputs unknown option" \
        "${log_opt_unknown}" \
        "${TEST_BASH}" "${scr_cln}" --bogus
then
    record_fail "clean_test_outputs.sh unknown option unexpectedly succeeded"
else
    assert_pattern_found \
        "${log_opt_unknown}" \
        "unknown option/parameter passed: '--bogus'" \
        "clean_test_outputs.sh unknown option emits useful error"
fi


#  Check explicitly scoped cleanup plans without deleting files
assert_clean_dry \
    "all" \
    "${log_dry_all}" \
    --all \
    --dry_run

assert_clean_dry \
    "fixtures" \
    "${log_dry_fix}" \
    --fixtures \
    --dry_run

assert_clean_dry \
    "outputs" \
    "${log_dry_out}" \
    --outputs \
    --dry_run

finish
