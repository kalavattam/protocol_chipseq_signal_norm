#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_fixture_paradigm.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="fixture paradigm"

# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"

dir_fix="${ROOT_REPO}/tests/fixtures"
fil_ignore="${ROOT_REPO}/.gitignore"


# List fixture roots as the filesystem sees them.
function list_roots_filesystem() {
    find "${dir_fix}" -mindepth 1 -maxdepth 1 -type d -print \
        | sed 's#.*/##' \
        | LC_ALL=C sort
}


# List fixture roots as Git's tracked inventory sees them.
function list_roots_tracked() {
    git -C "${ROOT_REPO}" ls-files -- 'tests/fixtures/*' \
        | cut -d/ -f3 \
        | LC_ALL=C sort -u
}


# List fixture roots as the ignore stanzas declare them.
function list_roots_declared() {
    sed -n 's#^tests/fixtures/\([^/]*\)/\*\*$#\1#p' "${fil_ignore}" \
        | LC_ALL=C sort -u
}


# Report whether a path is ignored by Git.
function path_is_ignored() {
    git -C "${ROOT_REPO}" check-ignore -q "${1}"
}


# Agree two independent derivations of the fixture-root population.
#
# An inspected-count assertion on its own proves only that the checker read
# what it was handed. The filesystem and Git's tracked inventory answer
# different questions, so a root present in one and absent from the other is a
# real defect rather than a counting error: every root must carry a tracked
# README and recipe, and nothing may sit under 'tests/fixtures' without them.
function check_population() {
    local observed=""
    local expected=""

    if (( ${#} == 0 )); then
        record_fail "no fixture roots were discovered under tests/fixtures"
        return 1
    fi

    expected="$(printf '%s\n' "$@")"
    observed="$(list_roots_tracked)"

    if [[ "${observed}" != "${expected}" ]]; then
        record_fail \
            "fixture roots on disk disagree with Git's tracked inventory"
        diff <(printf '%s\n' "${expected}") <(printf '%s\n' "${observed}") \
            >&2 || true
    else
        record_pass \
            "fixture-root population agrees across two derivations"
    fi

    record_pass "inspecting ${#} fixture root(s)"
    return 0
}


# Pin the ignore stanzas against the converted roots.
#
# A stanza naming no directory is a dead ignore rule left behind by a rename,
# and a converted root with no stanza would ship its generated products.
function check_declared_stanzas() {
    local observed=""
    local expected=""

    expected="$(printf '%s\n' "$@" | LC_ALL=C sort)"
    observed="$(list_roots_declared)"

    if [[ "${observed}" != "${expected}" ]]; then
        record_fail \
            "ignore stanzas disagree with the roots that carry a recipe"
        diff <(printf '%s\n' "${expected}") <(printf '%s\n' "${observed}") \
            >&2 || true
        return 1
    fi

    record_pass \
        "${#} ignore stanza(s) match the roots that carry a recipe"
    return 0
}


# Property 1: generated products are ignored and provenance is not.
function check_ignore_stanza() {
    local root="${1}"
    local findings=0

    if ! path_is_ignored "tests/fixtures/${root}/probe.generated"; then
        record_fail "generated products are not ignored: ${root}"
        findings=1
    fi

    if path_is_ignored "tests/fixtures/${root}/README.md"; then
        record_fail "fixture README is ignored: ${root}"
        findings=1
    fi

    if path_is_ignored "tests/fixtures/${root}/make.sh"; then
        record_fail "fixture recipe is ignored: ${root}"
        findings=1
    fi

    return "${findings}"
}


# Property 2: nothing but the README and the recipe is tracked.
function check_tracked_contents() {
    local root="${1}"
    local observed=""

    observed="$(
        git -C "${ROOT_REPO}" ls-files -- "tests/fixtures/${root}" \
            | sed "s#^tests/fixtures/${root}/##" \
            | LC_ALL=C sort \
            | paste -sd, -
    )"

    if [[ "${observed}" != "README.md,make.sh" ]]; then
        record_fail \
            "fixture root tracks more than its README and recipe:" \
            "${root} tracks '${observed}'"
        return 1
    fi

    return 0
}


# Property 3: the runner regenerates the root from a sentinel output.
function check_runner_wiring() {
    local root="${1}"
    local runner="${ROOT_REPO}/tests/run_tests.sh"

    if ! grep -q "ensure_fixture ${root} " "${runner}" \
        && ! grep -q "ensure_fixture ${root}\$" "${runner}"; then
        record_fail "fixture root is not wired into run_tests.sh: ${root}"
        return 1
    fi

    return 0
}


# Properties 3, 6, and 7: the recipe generates, and says so in the shared
# vocabulary.
#
# The validation check is the point of this function. A recipe that runs a
# checker over its own output has inverted the separation the standard states:
# 'make.sh' generates, and tests and checkers validate.
function check_recipe_shape() {
    local root="${1}"
    local recipe="${dir_fix}/${root}/make.sh"
    local findings=0
    local body=""

    if [[ ! -r "${recipe}" ]]; then
        record_fail "fixture recipe is missing or unreadable: ${root}"
        return 1
    fi

    if ! grep -q '^set -euo pipefail$' "${recipe}"; then
        record_fail "fixture recipe omits safe mode: ${root}"
        findings=1
    fi

    if ! grep -q 'BASH_VERSINFO' "${recipe}"; then
        record_fail "fixture recipe omits the Bash version guard: ${root}"
        findings=1
    fi

    if ! grep -q 'tests/support/fixture_helpers\.sh' "${recipe}"; then
        record_fail "fixture recipe omits the shared helper: ${root}"
        findings=1
    fi

    if ! grep -q '^succeed ' "${recipe}"; then
        record_fail "fixture recipe does not end in 'succeed': ${root}"
        findings=1
    fi

    # Read only the recipe's own statements: everything after safe mode is set,
    # with heredoc bodies removed. The Bash version guard above that line
    # necessarily writes its own diagnostics, because it runs before the shared
    # helpers can be sourced, and a fixture that is itself a shell script must
    # not be mistaken for the recipe that writes it. The opening delimiter is
    # captured by name and closed only on that same name, so a fixture carrying
    # its own heredoc does not end the recipe's.
    body="$(
        awk '
            /^set -euo pipefail$/ { body = 1; next }
            ! body { next }
            hdoc != "" {
                line = $0
                sub(/^[[:space:]]+/, "", line)
                if (line == hdoc) { hdoc = "" }
                next
            }
            match($0, /<<[[:space:]-]*["'"'"']?[A-Za-z_][A-Za-z0-9_]*/) {
                word = substr($0, RSTART, RLENGTH)
                sub(/^<<[[:space:]-]*["'"'"']?/, "", word)
                hdoc = word
                next
            }
            { print }
        ' "${recipe}"
    )"

    if grep -qE '(^|[^_[:alnum:]])(python3?|conda)[[:space:]]' \
        <<< "${body}"; then
        record_fail \
            "fixture recipe validates instead of generating:" \
            "${root} invokes a checker or interpreter"
        findings=1
    fi

    # Diagnostics go through the shared vocabulary. An 'echo' writing file
    # content is not a vocabulary fault, so this looks only at output sent to
    # standard error, which is what 'die' and 'note' exist to produce.
    if grep -qE '^[[:space:]]*echo[[:space:]].*>&2' <<< "${body}"; then
        record_fail \
            "fixture recipe writes diagnostics with bare 'echo' rather than" \
            "'die' or 'note': ${root}"
        findings=1
    fi

    return "${findings}"
}


# Property 6: the README carries the shared sections
function check_readme() {
    local root="${1}"
    local readme="${dir_fix}/${root}/README.md"
    local section=""
    local findings=0

    if [[ ! -s "${readme}" ]]; then
        record_fail "fixture README is missing or empty: ${root}"
        return 1
    fi

    for section in '## Files' '## Current and deferred test coverage'; do
        if ! grep -qF "${section}" "${readme}"; then
            record_fail \
                "fixture README omits the '${section}' section: ${root}"
            findings=1
        fi
    done

    return "${findings}"
}


# Property 5, the mechanism: an explicit stale sweep before writing.
#
# A recipe that rewrites every output unconditionally is idempotent without
# this, and two recipes relied on that. It is not enough. The sweep is what
# keeps regeneration correct after a later revision stops writing an output the
# previous revision wrote, which no amount of unconditional overwriting
# can undo. It also holds for the recipes this contract will not run, where the
# empirical check below is skipped.
function check_stale_sweep() {
    local root="${1}"

    if ! grep -qE '^[[:space:]]*rm_files?[[:space:]]' \
        "${dir_fix}/${root}/make.sh"; then
        record_fail "fixture recipe has no stale-output sweep: ${root}"
        return 1
    fi

    return 0
}


# Property 4 and 5, proved rather than argued: regenerate twice and compare.
#
# Recipes needing a managed environment or an external tool are not run here;
# the environment is the integration suite's concern, and running Samtools
# inside a repository contract would make this test's result depend on it.
# Those roots are reported as skipped rather than silently passed.
function run_recipe() {
    local root="${1}"
    local log="${2}"

    PATH="${TEST_ENV_PREFIX}/bin:${PATH}" \
        "${TEST_BASH}" "${dir_fix}/${root}/make.sh" > "${log}" 2>&1
}


# Report whether a failed recipe stopped on a missing prerequisite.
#
# 'require_cmd' and 'require_env' are the only two ways a recipe declines to
# run for want of its environment, and both die with a fixed phrase. Matching
# them is what separates "this machine lacks Samtools" from "this recipe is
# broken", so an absent tool is a counted skip and everything else fails.
function missing_prerequisite() {
    grep -qE "must be available|before generating fixtures" "${1}"
}


function check_regeneration() {
    local root="${1}"
    local before="${TEST_DIR_LOG}/fixture_paradigm/${root}.before"
    local after="${TEST_DIR_LOG}/fixture_paradigm/${root}.after"
    local log="${TEST_DIR_LOG}/fixture_paradigm/${root}.log"

    mkdir -p "$(dirname "${before}")"

    if ! run_recipe "${root}" "${log}"; then
        if missing_prerequisite "${log}"; then
            record_skip \
                "regeneration of ${root} needs a prerequisite this machine" \
                "does not provide: $(tail -1 "${log}")"
            return 0
        fi

        record_fail "fixture recipe failed to run: ${root}"
        tail -5 "${log}" >&2
        return 1
    fi

    inventory_fixture_root "${root}" > "${before}"

    if ! run_recipe "${root}" "${log}"; then
        record_fail "fixture recipe failed to rerun: ${root}"
        tail -5 "${log}" >&2
        return 1
    fi

    inventory_fixture_root "${root}" > "${after}"

    if ! cmp -s "${before}" "${after}"; then
        record_fail "fixture regeneration is not idempotent: ${root}"
        diff -u "${before}" "${after}" >&2 || true
        return 1
    fi

    if [[ ! -s "${before}" ]]; then
        record_fail "fixture recipe generated no output: ${root}"
        return 1
    fi

    return 0
}


# Print a digest of every generated file in one fixture root.
function inventory_fixture_root() {
    local root="${1}"

    find "${dir_fix}/${root}" -type f \
        ! -name README.md ! -name make.sh -print \
        | LC_ALL=C sort \
        | while IFS= read -r file; do
            printf '%s  %s\n' \
                "$(shasum -a 256 "${file}" | cut -d' ' -f1)" \
                "${file#"${dir_fix}/"}"
        done
}


# Require every generated 'cases.json' to be in canonical JSON source form.
#
# The ignore stanza takes these files out of the checker's discovery, so the
# repository inventory can no longer report them and their applicability entry
# in 'dev/config/rules.toml' would become a glob matching nothing. The checker
# accepts explicit paths and bypasses discovery, so the rule keeps reaching the
# generated output here rather than losing it to concealment.
function check_generated_json() {
    local root=""
    local file=""
    local -a arr_json=()
    local output=""

    for root in "$@"; do
        file="${dir_fix}/${root}/cases.json"
        [[ -f "${file}" ]] && arr_json+=( "${file}" )
    done

    if (( ${#arr_json[@]} == 0 )); then
        record_fail "no generated 'cases.json' fixtures were found to check"
        return 1
    fi

    if output="$(
        PYTHONDONTWRITEBYTECODE=1 \
            "${TEST_PYTHON}" -m dev.audit.json_source_form \
                --root "${ROOT_REPO}" "${arr_json[@]}" 2>&1
    )"; then
        record_pass \
            "${#arr_json[@]} generated 'cases.json' fixture(s) are in" \
            "canonical form"
        return 0
    fi

    record_fail "generated 'cases.json' fixtures are not in canonical form"
    printf '%s\n' "${output}" >&2
    return 1
}


print_section "${TEST_NAME}"

declare -a arr_roots=()
while IFS= read -r root; do
    arr_roots+=( "${root}" )
done < <(list_roots_filesystem)

check_population "${arr_roots[@]}" || finish

declare -a arr_converted=()
declare -a arr_unconverted=()
findings=0

for root in "${arr_roots[@]}"; do
    if [[ -r "${dir_fix}/${root}/make.sh" ]]; then
        arr_converted+=( "${root}" )
    else
        arr_unconverted+=( "${root}" )
    fi
done

check_declared_stanzas "${arr_converted[@]}" || findings=1

for root in "${arr_converted[@]}"; do
    check_ignore_stanza "${root}" || findings=1
    check_tracked_contents "${root}" || findings=1
    check_runner_wiring "${root}" || findings=1
    check_recipe_shape "${root}" || findings=1
    check_stale_sweep "${root}" || findings=1
    check_readme "${root}" || findings=1
    check_regeneration "${root}" || findings=1
done

if (( findings == 0 )); then
    record_pass \
        "${#arr_converted[@]} fixture root(s) satisfy every paradigm property"
fi

check_generated_json "${arr_converted[@]}" || findings=1

# Report the unconverted remainder without failing, so it stays visible instead
# of disappearing behind a green converted cohort.
if (( ${#arr_unconverted[@]} > 0 )); then
    record_warn \
        "${#arr_unconverted[@]} unconverted fixture root(s), tracked and" \
        "unread: ${arr_unconverted[*]}; their consuming tests build inputs" \
        "inline instead, so these roots are filled by the pass that moves" \
        "inline fixture construction out of tests"
fi

finish
