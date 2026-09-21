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

"${TEST_PYTHON:-python3}" - "${script}" > "${rep_arity}" 2>&1 << 'PY'
import re
import sys

source = open(sys.argv[1], encoding="utf-8").read().splitlines()

guards = {}
current = None
for line in source:
    named = re.match(r"function (run_comp_\w+)\(\)", line)
    if named:
        current = named.group(1)
    guard = re.search(r"\[\[ \$# -ne (\d+) \]\]", line)
    if guard and current:
        guards.setdefault(current, int(guard.group(1)))

calls = []
for index, line in enumerate(source):
    called = re.match(r"\s*(run_comp_\w+) \\\s*$", line)
    if not called:
        continue
    count, cursor = 0, index + 1
    while cursor < len(source):
        count += 1
        if not source[cursor].rstrip().endswith("\\"):
            break
        cursor += 1
    calls.append((called.group(1), index + 1, count))

for name, line_no, count in calls:
    want = guards.get(name)
    status = "OK" if want == count else "MISMATCH"
    print(f"{status} {name} line {line_no}: passes {count}, guard expects {want}")

print(f"CALLS {len(calls)}")
PY

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


# Report forwarding: a bare '--report_n_frg' is valid in 'compute_signal.py',
# which derives its path, so a regression here fails silently rather than
# loudly. Require every append to carry a quoted value under a non-empty guard.
rep_fwd="${TEST_DIR_TMP}/submit_compute_signal_interface/forwarding.txt"
mkdir -p "$(dirname "${rep_fwd}")"

"${TEST_PYTHON:-python3}" - "${script}" > "${rep_fwd}" 2>&1 << 'PY'
import re
import sys

source = open(sys.argv[1], encoding="utf-8").read().splitlines()

appends = 0
for index, line in enumerate(source):
    append = re.match(
        r'\s*cmd\+=\( --report_(n_frg|n_ovlp)(.*)\)\s*$',
        line,
    )

    if not append:
        continue

    appends += 1
    label, tail = append.group(1), append.group(2).strip()
    var = f"report_{label}"
    line_no = index + 1

    if tail != f'"${{{var}}}"':
        print(f"BARE --report_{label} line {line_no}: appends {tail!r}")
        continue

    window = source[max(0, index - 2):index]
    guard = f'if [[ -n "${{{var}}}" ]]; then'

    if any(entry.strip() == guard for entry in window):
        print(f"OK --report_{label} line {line_no}: quoted under '-n' guard")
    else:
        print(f"UNGUARDED --report_{label} line {line_no}: no '-n' on {var}")

print(f"APPENDS {appends}")
PY

if [[ -s "${rep_fwd}" ]]; then
    assert_pattern_absent \
        "${rep_fwd}" \
        "BARE" \
        "submit never appends a valueless '--report_n_frg' / '--report_n_ovlp'"

    assert_pattern_absent \
        "${rep_fwd}" \
        "UNGUARDED" \
        "each report append sits under a non-empty guard on its own variable"

    assert_pattern_found \
        "${rep_fwd}" \
        "APPENDS 2" \
        "both report forwarding sites were found and checked"
else
    record_fail "report forwarding report was not produced"
fi


finish "$@"
