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


function run_chunk_alias_case() {
    local interface="${1}"
    local script_path="${2}"
    local alias="${3}"
    local accepted="${4}"
    local safe_alias="${alias//[^A-Za-z0-9]/_}"
    local log="${dir_log}/${interface}_${safe_alias}.log"
    local -a command=( "${TEST_BASH}" "${script_path}" )

    if [[ "${interface}" == "submit" ]]; then
        command+=( --dir_scr "${ROOT_REPO}/bin" )
    fi
    command+=( "${alias}" 17 )

    if "${command[@]}" > "${log}" 2>&1; then
        record_fail \
            "${interface} ${alias} unexpectedly completed validation"
        return
    fi

    if [[ "${accepted}" == "true" ]]; then
        assert_pattern_absent \
            "${log}" \
            "unknown option/parameter passed: '${alias}'" \
            "${interface} accepts ${alias}"
    else
        assert_pattern_found \
            "${log}" \
            "unknown option/parameter passed: '${alias}'." \
            "${interface} rejects ${alias}"
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
assert_pattern_found \
    "${log_help}" \
    '^  -ck, --chunk_size : int$' \
    "submit help exposes canonical chunk-size aliases"
assert_pattern_absent \
    "${log_help}" \
    '--chunk-size|--chnk[_-]size' \
    "submit help hides chunk-size compatibility aliases"

if "${TEST_BASH}" "${execute_script}" --help > "${log_execute_help}" 2>&1; then
    record_pass "execute_compute_signal.sh --help exits 0"
else
    record_fail "execute_compute_signal.sh --help failed"
fi
assert_pattern_found \
    "${log_execute_help}" \
    '^  -ck, --chunk_size : int$' \
    "execute help exposes canonical chunk-size aliases"
assert_pattern_absent \
    "${log_execute_help}" \
    '--chunk-size|--chnk[_-]size' \
    "execute help hides chunk-size compatibility aliases"

run_alias_case "short_cuf" "-cuf" true
run_alias_case "canonical_long" "--csv_usr_frg" true
run_alias_case "hidden_hyphen" "--csv-usr-frg" true
run_alias_case "retired_short" "-uf" false
run_alias_case "retired_long_underscore" "--usr_frg" false
run_alias_case "retired_long_hyphen" "--usr-frg" false

for interface in submit execute; do
    if [[ "${interface}" == "submit" ]]; then
        chunk_script="${script}"
    else
        chunk_script="${execute_script}"
    fi
    for alias in \
        -ck \
        --chunk_size \
        --chunk-size \
        --chnk_size \
        --chnk-size
    do
        run_chunk_alias_case \
            "${interface}" "${chunk_script}" "${alias}" true
    done
    for alias in \
        --chunk_s \
        --chnk_s \
        --chunck_size \
        --chnk__size \
        --chunk_sizes
    do
        run_chunk_alias_case \
            "${interface}" "${chunk_script}" "${alias}" false
    done
done

assert_pattern_found \
    "${execute_script}" \
    "--chunk_size \"\${chunk_size}\"" \
    "execute forwards canonical --chunk_size"
assert_pattern_found \
    "${script}" \
    "--chunk_size \"\${chunk_size}\"" \
    "submit forwards canonical --chunk_size"

finish "$@"
