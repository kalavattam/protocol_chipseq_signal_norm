#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_boolean_contracts.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="Boolean-like contracts"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"

#  Source the canonical Boolean helper and representative public callers
# shellcheck source=lib/bash/core/format_outputs.sh
source "${ROOT_REPO}/lib/bash/core/format_outputs.sh"
# shellcheck source=lib/bash/core/check_args.sh
source "${ROOT_REPO}/lib/bash/core/check_args.sh"
# shellcheck source=lib/bash/core/check_inputs.sh
source "${ROOT_REPO}/lib/bash/core/check_inputs.sh"


#  Run normalize_bool without letting expected failures terminate this test
function run_normalize() {
    local value="${1-}"
    local name="${2:-test_bool}"

    set +e
    NORMALIZE_OUT="$(normalize_bool "${value}" "${name}" 2> "${log_err}")"
    NORMALIZE_STATUS=$?
    set -e
}


#  Require one token to normalize to one canonical value
function assert_normalized() {
    local token="${1-}"
    local expected="${2:-}"

    run_normalize "${token}"
    if [[ "${NORMALIZE_STATUS}" -eq 0 && "${NORMALIZE_OUT}" == "${expected}" ]]; then
        record_pass "normalize_bool '${token}' -> '${expected}'"
    else
        record_fail \
            "normalize_bool '${token}' expected '${expected}', got" \
            "status=${NORMALIZE_STATUS} output='${NORMALIZE_OUT}'"
    fi
}


#  Require one invalid token to fail with the precise diagnostic
function assert_rejected() {
    local label="${1:-invalid token}"
    local token="${2-}"

    run_normalize "${token}" strict_value
    if (( NORMALIZE_STATUS != 0 )); then
        record_pass "${label} returns nonzero"
    else
        record_fail "${label} unexpectedly normalized to '${NORMALIZE_OUT}'"
    fi

    if grep -Fq \
        "argument 'strict_value' must be Boolean-like" \
        "${log_err}"
    then
        record_pass "${label} emits precise invalid-value diagnostic"
    else
        record_fail "${label} lacks precise invalid-value diagnostic"
    fi
}


#  Check the maintained test-gate environment contract
function check_test_gate() {
    local token="${1-}"
    local expected="${2:-}"
    local output status

    set +e
    output="$(RUN_PARALLEL="${token}" normalize_test_gate RUN_PARALLEL 2> "${log_err}")"
    status=$?
    set -e

    if [[ "${status}" -eq 0 && "${output}" == "${expected}" ]]; then
        record_pass "Boolean environment token '${token}' -> '${expected}'"
    else
        record_fail \
            "Boolean environment token '${token}' expected '${expected}', got" \
            "status=${status} output='${output}'"
    fi
}


print_section "${TEST_NAME}"

dir_log="${TEST_DIR_LOG}/boolean_contracts"
mkdir -p "${dir_log}"
log_err="${dir_log}/normalize.err"
log_audit="${dir_log}/audit.log"

for token in true t yes y 1 TRUE T YES Y; do
    assert_normalized "${token}" true
done

for token in false f no n 0 FALSE F NO N; do
    assert_normalized "${token}" false
done

assert_rejected "leading-space token" " true"
assert_rejected "trailing-space token" "true "
assert_rejected "leading-tab token" $'\ttrue'
assert_rejected "trailing-tab token" $'true\t'
assert_rejected "leading-newline token" $'\ntrue'
assert_rejected "trailing-newline token" $'true\n'
assert_rejected "arbitrary token" sometimes
assert_rejected "empty required token" ""

stale=previous
set +e
stale="$(normalize_bool sometimes stale_value 2> "${log_err}")"
status=$?
set -e
if [[ "${status}" -ne 0 && -z "${stale}" ]]; then
    record_pass "failed normalization does not retain a stale value"
else
    record_fail \
        "failed normalization leaked stale value: status=${status}" \
        "value='${stale}'"
fi

for token in true YES 1 false No 0; do
    expected=true
    case "${token,,}" in
        false|no|0) expected=false ;;
    esac
    check_test_gate "${token}" "${expected}"
done

unset RUN_PARALLEL
if [[ "$(normalize_test_gate RUN_PARALLEL)" == "false" ]]; then
    record_pass "unset optional Boolean environment gate stays disabled"
else
    record_fail "unset optional Boolean environment gate is not disabled"
fi

if [[ "$(RUN_PARALLEL='' normalize_test_gate RUN_PARALLEL)" == "false" ]]; then
    record_pass "empty optional Boolean environment gate stays disabled"
else
    record_fail "empty optional Boolean environment gate is not disabled"
fi

set +e
RUN_PARALLEL=sometimes normalize_test_gate RUN_PARALLEL \
    > /dev/null 2> "${log_err}"
status=$?
set -e
if (( status != 0 )); then
    record_pass "invalid nonempty Boolean environment gate returns nonzero"
else
    record_fail "invalid nonempty Boolean environment gate was accepted"
fi

tmp_file="${TEST_DIR_TMP}/boolean_contracts_file"
: > "${tmp_file}"
if validate_file test_file "${tmp_file}" 0 YES; then
    record_pass "validate_file accepts a canonical true-like token"
else
    record_fail "validate_file rejected a canonical true-like token"
fi

set +e
validate_file test_file "${tmp_file}" 0 sometimes \
    > /dev/null 2> "${log_err}"
status=$?
set -e
if (( status != 0 )); then
    record_pass "validate_file rejects an invalid Boolean token"
else
    record_fail "validate_file accepted an invalid Boolean token"
fi

before_fail="${TEST_FAIL}"
set +e
assert_scaling_factor_header \
    "${tmp_file}" '^sample$' sometimes "invalid expect" \
    > /dev/null 2> "${log_err}"
status=$?
set -e
if [[ "${status}" -ne 0 && "${TEST_FAIL}" -eq "${before_fail}" ]]; then
    record_pass "invalid assertion Boolean fails before mutating test results"
else
    record_fail \
        "invalid assertion Boolean did not preserve results:" \
        "status=${status} before=${before_fail} after=${TEST_FAIL}"
fi

if python3 -m dev.audit.boolean_contracts \
    --root "${ROOT_REPO}" > "${log_audit}"
then
    record_pass "repository Boolean-contract audit"
else
    record_fail \
        "repository Boolean-contract audit; see" \
        "$(print_relpath "${log_audit}")"
fi

# shellcheck disable=SC2119  # finish intentionally receives no script arguments.
finish
