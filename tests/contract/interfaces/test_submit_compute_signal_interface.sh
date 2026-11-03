#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_compute_signal_interface.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="submit_compute_signal interface"

# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


function run_alias_case() {
    local label="${1}"
    local alias="${2}"
    local accepted="${3}"
    local log="${dir_log}/${label}.log"

    if "${TEST_BASH}" "${script}" \
        --dir_scr "${ROOT_REPO}/bin" \
        "${alias}" NA > "${log}" 2>&1
    then
        record_fail "${label} unexpectedly completed validation"
        return
    fi

    if [[ "${accepted}" == "true" ]]; then
        assert_pattern_absent \
            "${log}" \
            "unknown option/parameter passed: '${alias}'" \
            "${alias} is accepted by parse_args"
    else
        assert_pattern_found \
            "${log}" \
            "unknown option/parameter passed: '${alias}'." \
            "${alias} is rejected with the canonical error"
    fi
}


script="${ROOT_REPO}/bin/submit_compute_signal.sh"
execute_script="${ROOT_REPO}/bin/execute_compute_signal.sh"
dir_log="${TEST_DIR_LOG}/submit_compute_signal_interface"
log_help="${dir_log}/help.log"
log_execute_help="${dir_log}/execute_help.log"


print_section "${TEST_NAME}"

mkdir -p "${dir_log}"

if "${TEST_BASH}" "${script}" --help > "${log_help}" 2>&1; then
    record_pass "submit_compute_signal.sh --help exits 0"
else
    record_fail "submit_compute_signal.sh --help failed"
fi

assert_pattern_found \
    "${log_help}" \
    '^    \[--siz_bin .*\[--csv_usr_frg <csv>\]$' \
    "Usage exposes only canonical --csv_usr_frg"
assert_pattern_found \
    "${log_help}" \
    '^  -cuf, --csv_usr_frg : list of int$' \
    "Parameters exposes -cuf and --csv_usr_frg"
assert_pattern_absent "${log_help}" '(^|[[:space:],])-uf([,[:space:]]|$)' \
    "help omits retired -uf"
assert_pattern_absent "${log_help}" '--usr[_-]frg' \
    "help omits retired usr_frg long aliases"
if "${TEST_BASH}" "${execute_script}" --help > "${log_execute_help}" 2>&1; then
    record_pass "execute_compute_signal.sh --help exits 0"
else
    record_fail "execute_compute_signal.sh --help failed"
fi

run_alias_case "short_cuf" "-cuf" true
run_alias_case "canonical_long" "--csv_usr_frg" true
run_alias_case "hidden_hyphen" "--csv-usr-frg" true
run_alias_case "retired_short" "-uf" false
run_alias_case "retired_long_underscore" "--usr_frg" false
run_alias_case "retired_long_hyphen" "--usr-frg" false


# Positional arity: renumbering a dispatch helper broke the wrapper before, as
# a guard updated without its call sites rejects every dispatch while unit
# tests still pass. Derive both numbers from source and require each call site
# to match its own guard.
rep_arity="${TEST_DIR_TMP}/submit_compute_signal_interface/arity.txt"
mkdir -p "$(dirname "${rep_arity}")"

"${TEST_PYTHON:-python3}" \
    "${ROOT_REPO}/tests/support/check_submit_compute_signal.py" \
    arity \
    "${script}" \
    > "${rep_arity}" 2>&1

if [[ -s "${rep_arity}" ]]; then
    assert_pattern_absent \
        "${rep_arity}" \
        "MISMATCH" \
        "every run_comp_* call site matches its own arity guard"

    assert_pattern_found \
        "${rep_arity}" \
        "CALLS 3" \
        "all three run_comp_* call sites were found and checked"
else
    record_fail "arity report was not produced"
fi


# Report forwarding. A bare '--report_n_frg' makes the CLI derive its own path:
# intended under the 'DERIVE' guard, a silent fault anywhere else. Pin both
# directions, and keep path appends quoted under a non-empty guard so an empty
# variable fails here rather than deriving.
rep_fwd="${TEST_DIR_TMP}/submit_compute_signal_interface/forwarding.txt"
mkdir -p "$(dirname "${rep_fwd}")"

"${TEST_PYTHON:-python3}" \
    "${ROOT_REPO}/tests/support/check_submit_compute_signal.py" \
    forwarding \
    "${script}" \
    > "${rep_fwd}" 2>&1

if [[ -s "${rep_fwd}" ]]; then
    assert_pattern_absent \
        "${rep_fwd}" \
        "BARE-UNGUARDED" \
        "a bare report append occurs only under the 'DERIVE' guard"

    assert_pattern_absent \
        "${rep_fwd}" \
        "UNGUARDED --report" \
        "each valued report append sits under a non-empty guard"

    assert_pattern_absent \
        "${rep_fwd}" \
        "UNEXPECTED" \
        "no report append carries an unrecognized argument form"

    assert_pattern_found \
        "${rep_fwd}" \
        "BARE 2" \
        "the derive branch bare-appends both report flags"

    assert_pattern_found \
        "${rep_fwd}" \
        "VALUED 2" \
        "the explicit-path branch still forwards both report flags by value"
else
    record_fail "report forwarding report was not produced"
fi


finish "$@"
