#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_json_source_form_gate.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="JSON source-form gate"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# The checker discovers every maintained authored JSON path, but only the
# configuration and schema populations were swept into canonical form. Fixture
# case files are rewritten by the fixture-paradigm conversion, so they are an
# inventoried remainder here rather than a silent omission.
declare -a arr_blocking_roots=(
    "dev/config"
    "dev/schemas"
)

log_dir="${TEST_DIR_LOG}/json_source_form_gate"
blocking_log="${log_dir}/blocking.log"
inventory_log="${log_dir}/inventory.log"
inventory_json="${log_dir}/inventory.json"
mkdir -p "${log_dir}"

print_section "${TEST_NAME}"


function count_findings() {
    local log="${1}"

    # 'JSON.SOURCE.FORM: N path(s)' is the checker's own summary row.
    grep -cE '^JSON\.SOURCE\.FORM: [^ ]+:[0-9]+: ' "${log}" 2> /dev/null \
        || true
}


function scanned_paths() {
    local log="${1}"

    sed -nE 's/^JSON\.SOURCE\.FORM: ([0-9]+) path\(s\).*/\1/p' "${log}" \
        | tail -1
}


# Resolve the swept cohort explicitly, so a zero-finding result can be checked
# against a known inspected count rather than trusted on its own. Filtering a
# directory pathspec in the loop avoids the wildcard pathspec forms, one of
# which silently matched only the nested directory while this gate was written.
declare -a arr_cohort=()

while IFS= read -r -d '' relative; do
    [[ "${relative}" == *.json ]] || continue
    arr_cohort+=( "${relative}" )
done < <(git -C "${ROOT_REPO}" ls-files -z -- "${arr_blocking_roots[@]}")

if (( ${#arr_cohort[@]} == 0 )); then
    record_fail "swept JSON cohort resolved to no files"
    finish
fi


# Enforce the canonical rendering over the swept cohort.
set +e
(
    cd "${ROOT_REPO}" || exit 2
    PYTHONDONTWRITEBYTECODE=1 python3 -m dev.audit.json_source_form \
        --root "${ROOT_REPO}" "${arr_cohort[@]}"
) > "${blocking_log}" 2>&1
blocking_status="$?"
set -e

blocking_scanned="$(scanned_paths "${blocking_log}")"
blocking_findings="$(count_findings "${blocking_log}")"

# Prove the checker actually read the cohort before trusting a zero count.
if [[ "${blocking_scanned}" == "${#arr_cohort[@]}" ]]; then
    record_pass \
        "JSON checker scanned all ${#arr_cohort[@]} cohort file(s)"
else
    record_fail \
        "JSON checker scanned '${blocking_scanned:-none}' of" \
        "${#arr_cohort[@]} cohort file(s)"
    sed -n '1,20p' "${blocking_log}" >&2
fi

if (( blocking_status > 1 )); then
    record_fail "JSON source-form checker failed to run over the cohort"
fi

if (( blocking_findings == 0 )); then
    record_pass "swept JSON cohort is in canonical form"
else
    record_fail \
        "swept JSON cohort has ${blocking_findings} finding(s)"
    grep -vE '^JSON\.SOURCE\.FORM: [0-9]+ path' "${blocking_log}" \
        | sed -n '1,40p' >&2
fi


# Report the unswept remainder without failing, so it stays visible instead of
# disappearing behind a green cohort.
set +e
(
    cd "${ROOT_REPO}" || exit 2
    PYTHONDONTWRITEBYTECODE=1 python3 -m dev.audit.json_source_form \
        --root "${ROOT_REPO}"
) > "${inventory_log}" 2>&1
inventory_status="$?"
(
    cd "${ROOT_REPO}" || exit 2
    PYTHONDONTWRITEBYTECODE=1 python3 -m dev.audit.json_source_form \
        --root "${ROOT_REPO}" --json
) > "${inventory_json}" 2>&1
set -e


# Derive the cohort a second time from the checker's own discovery, which uses
# a '*.json' pathspec rather than a directory pathspec. Agreeing counts from
# two different derivations is what makes the inspected-count assertion above
# mean something; comparing one list against itself would not.
discovered_cohort="$(
    grep -oE '"(dev/config|dev/schemas)/[^"]+\.json"' "${inventory_json}" \
        | sort -u \
        | grep -c . \
        || true
)"

if [[ "${discovered_cohort}" == "${#arr_cohort[@]}" ]]; then
    record_pass \
        "cohort of ${#arr_cohort[@]} file(s) agrees with checker discovery"
else
    record_fail \
        "cohort selection resolved ${#arr_cohort[@]} file(s) but checker" \
        "discovery found ${discovered_cohort} under the same roots"
fi

if (( inventory_status > 1 )); then
    record_fail "JSON source-form inventory failed to run"
else
    total_scanned="$(scanned_paths "${inventory_log}")"
    total_findings="$(count_findings "${inventory_log}")"
    residual=$(( total_findings - blocking_findings ))

    if [[ -z "${total_scanned}" ]] \
        || (( total_scanned < ${#arr_cohort[@]} )); then
        record_fail \
            "repository inventory scanned '${total_scanned:-none}' path(s)," \
            "fewer than the ${#arr_cohort[@]} cohort file(s)"
    else
        record_pass \
            "maintained JSON inventory scanned ${total_scanned} path(s)"
    fi

    record_warn \
        "fixture-paradigm residual: ${residual} finding(s) outside the swept" \
        "cohort, deferred to the fixture-directory conversion that rewrites" \
        "tests/fixtures/*/cases.json"
fi

finish
