#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_compute_signal_interface.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
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
dir_log="${TEST_DIR_LOG}/submit_compute_signal_interface"
log_help="${dir_log}/help.log"


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

run_alias_case "short_cuf" "-cuf" true
run_alias_case "canonical_long" "--csv_usr_frg" true
run_alias_case "hidden_hyphen" "--csv-usr-frg" true
run_alias_case "retired_short" "-uf" false
run_alias_case "retired_long_underscore" "--usr_frg" false
run_alias_case "retired_long_hyphen" "--usr-frg" false

finish "$@"
