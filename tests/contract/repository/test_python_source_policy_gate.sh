#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_python_source_policy_gate.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="Python source-policy gate"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# The registered command scans the whole maintained tree, but only this cohort
# is migrated and therefore blocking. The remainder is a recorded inventory
# under S3-MIG-007, never a silent omission.
declare -a arr_blocking_roots=(
    "src/protocol_chipseq_signal_norm"
)
declare -a arr_blocking_files=(
    "dev/audit/help_examples.py"
    "dev/audit/help_option_order.py"
)

# 'SOURCE.PROSE.WRAP' is advisory repository-wide until S3-MIG-007 clears its
# residual, so it is inventoried rather than enforced even on this cohort.
advisory_rule="SOURCE.PROSE.WRAP"

log_dir="${TEST_DIR_LOG}/python_source_policy_gate"
blocking_log="${log_dir}/blocking.log"
inventory_log="${log_dir}/inventory.log"
mkdir -p "${log_dir}"

print_section "${TEST_NAME}"


function count_findings() {
    local log="${1}"
    local pattern="${2}"

    # 'PY.SOURCE.POLICY:' is the checker's own summary row, not a finding.
    grep -E "${pattern}" "${log}" 2> /dev/null \
        | grep -cvE '^PY\.SOURCE\.POLICY: ' \
        || true
}


function scanned_paths() {
    local log="${1}"

    sed -nE 's/^PY\.SOURCE\.POLICY: ([0-9]+) path\(s\).*/\1/p' "${log}" \
        | tail -1
}


# The checker accepts files only, so expand each cohort root explicitly. A
# directory argument raises FileNotFoundError and would scan nothing.
declare -a arr_cohort=()

while IFS= read -r -d '' relative; do
    arr_cohort+=( "${relative}" )
done < <(
    git -C "${ROOT_REPO}" ls-files -z -- \
        "${arr_blocking_roots[@]/%//**/*.py}" "${arr_blocking_files[@]}"
)

if (( ${#arr_cohort[@]} == 0 )); then
    record_fail "migrated Python cohort resolved to no files"
    finish
fi


# Enforce every already-migrated rule over the repaired cohort.
set +e
(
    cd "${ROOT_REPO}" || exit 2
    PYTHONDONTWRITEBYTECODE=1 python3 -m dev.audit.python_source_policy \
        --root "${ROOT_REPO}" "${arr_cohort[@]}"
) > "${blocking_log}" 2>&1
blocking_status="$?"
set -e

blocking_scanned="$(scanned_paths "${blocking_log}")"
blocking_total="$(count_findings "${blocking_log}" '^[A-Z][A-Z0-9._]+: ')"
blocking_advisory="$(count_findings "${blocking_log}" "^${advisory_rule}: ")"
blocking_enforced=$(( blocking_total - blocking_advisory ))

# Prove the checker actually read the cohort before trusting a zero count.
if [[ "${blocking_scanned}" == "${#arr_cohort[@]}" ]]; then
    record_pass \
        "source-policy checker scanned all ${#arr_cohort[@]} cohort file(s)"
else
    record_fail \
        "source-policy checker scanned '${blocking_scanned:-none}' of" \
        "${#arr_cohort[@]} cohort file(s)"
    sed -n '1,20p' "${blocking_log}" >&2
fi

if (( blocking_status > 1 )); then
    record_fail "Python source-policy checker failed to run over the cohort"
fi

if (( blocking_enforced == 0 )); then
    record_pass "migrated Python cohort has no enforced source-policy findings"
else
    record_fail \
        "migrated Python cohort has ${blocking_enforced} enforced finding(s)"
    grep -vE "^${advisory_rule}: |^PY\.SOURCE\.POLICY: " "${blocking_log}" \
        | sed -n '1,40p' >&2
fi


# Report the unmigrated remainder without failing, so the residual stays
# visible instead of disappearing behind a green cohort.
set +e
(
    cd "${ROOT_REPO}" || exit 2
    PYTHONDONTWRITEBYTECODE=1 python3 -m dev.audit.python_source_policy \
        --root "${ROOT_REPO}" --all-maintained
) > "${inventory_log}" 2>&1
inventory_status="$?"
set -e

if (( inventory_status > 1 )); then
    record_fail "Python source-policy inventory failed to run"
else
    total_findings="$(count_findings "${inventory_log}" '^[A-Z][A-Z0-9._]+: ')"
    advisory_findings="$(
        count_findings "${inventory_log}" "^${advisory_rule}: "
    )"
    residual=$(( total_findings - blocking_total ))
    residual_advisory=$(( advisory_findings - blocking_advisory ))

    record_pass "maintained Python inventory completed"
    record_warn \
        "S3-MIG-007 residual: ${residual} finding(s) outside the migrated" \
        "cohort, of which ${residual_advisory} are advisory" \
        "${advisory_rule}; the cohort itself carries" \
        "${blocking_advisory} advisory finding(s)"
fi

finish
