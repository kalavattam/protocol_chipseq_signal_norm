#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_help_output.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="help output"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


print_section "${TEST_NAME}"


#  Fail when one rendered help document repeats a top-level section
function assert_no_duplicate_rendered_sections() {
    local file="${1:?}"
    local lbl="${2:?}"
    local found=0
    local line

    while IFS= read -r line; do
        found=1
        record_fail \
            "duplicate rendered help section in ${lbl}: ${line}; see" \
            "$(print_relpath "${file}")"
    done < <(
        awk '
            BEGIN {
                heading = "Usage|Parameters|Returns|Notes|References|See Also|Examples"
            }
            /^-{3,}$/ && prev ~ "^(" heading ")$" {
                count[prev]++
            }
            { prev = $0 }
            END {
                for (sect in count) {
                    if (count[sect] > 1) {
                        print sect " appears " count[sect] " times"
                    }
                }
            }
        ' "${file}" | sort
    )

    if (( found == 0 )); then
        record_pass "${lbl} has no duplicate rendered sections"
    fi
}


#  Fail when maintainer TODO markers leak into rendered user-facing help
function assert_no_rendered_todo() {
    local file="${1:?}"
    local lbl="${2:?}"

    if grep -Eq -- '(^|[^A-Za-z])#?TODO([^A-Za-z]|$)' "${file}"; then
        record_fail \
            "rendered help contains TODO marker in ${lbl}; see" \
            "$(print_relpath "${file}")"
    else
        record_pass "${lbl} has no rendered TODO markers"
    fi
}


#  Fail when sourced helper/function inventories leak into rendered help
function assert_no_rendered_helper_inventory() {
    local file="${1:?}"
    local lbl="${2:?}"
    local found=0
    local line

    while IFS= read -r line; do
        found=1
        record_fail \
            "rendered help contains helper/function inventory in ${lbl}:" \
            "$(print_relpath "${file}"):${line}"
    done < <(
        perl -ne '
            if (/^\s*(?:-\s*)?(?:Functions|Function scripts|Sourced functions|Shell functions|Helpers|Helper functions|Sourced helpers)\s*:?\s*$/) {
                print "$.:$_";
                next;
            }

            if (/^\s*[+-]\s+[A-Za-z_][A-Za-z0-9_.\/-]*\s+##\s*[A-Za-z_][A-Za-z0-9_.\/-]*\s*##\s*$/) {
                print "$.:$_";
                next;
            }

            if (/^\s*[+-]\s+help\/[A-Za-z0-9_.\/-]+\s*$/) {
                print "$.:$_";
                next;
            }
        ' "${file}"
    )

    if (( found == 0 )); then
        record_pass "${lbl} has no rendered helper/function inventory"
    fi
}


#  Run the authoritative structural checker on one rendered help document
function assert_rendered_help_structure() {
    local file="${1:?}"
    local lbl="${2:?}"
    local output=""
    local line

    if output="$(
        PYTHONDONTWRITEBYTECODE=1 \
            python3 "${ROOT_REPO}/dev/audit/help_style.py" \
                --root "${ROOT_REPO}" \
                --rendered "${file}" 2>&1
    )"; then
        record_pass "${lbl} has valid rendered help structure"
        return 0
    fi

    while IFS= read -r line; do
        record_fail "rendered help structure in ${lbl}: ${line}"
    done <<< "${output}"
}


#  Protect the representative multiline Examples section byte-for-byte
function assert_combine_parts_examples_fixture() {
    local file="${1:?}"
    local expected="${ROOT_REPO}/tests/fixtures/help/combine_parts_scaling_factor.examples.txt"
    local observed="${TEST_DIR_LOG}/help/combine_parts_scaling_factor.examples.txt"
    local line

    if ! \
        PYTHONDONTWRITEBYTECODE=1 \
            python3 "${ROOT_REPO}/dev/audit/help_heredoc_reflow.py" \
                --extract-rendered-examples "${file}" > "${observed}"
    then
        record_fail "could not extract combine_parts_scaling_factor.sh Examples"
        return 0
    fi

    if cmp -s "${expected}" "${observed}"; then
        record_pass "combine_parts_scaling_factor.sh Examples match exact fixture"
        return 0
    fi

    record_fail "combine_parts_scaling_factor.sh Examples differ from exact fixture"
    while IFS= read -r line; do
        record_fail "combine_parts_scaling_factor.sh fixture diff: ${line}"
    done < <(diff -u "${expected}" "${observed}" || true)
}


#  Check short and detailed help rendering for top-level shell scripts
while IFS= read -r file; do
    rel="$(print_relpath "${file}")"
    log="${TEST_DIR_LOG}/help/${rel//\//__}.help.log"

    if \
        run_capture \
            "help ${rel}" \
            "${log}" \
            "${TEST_BASH}" "${file}" --help
    then
        record_pass "--help exits 0 for ${rel}"
        assert_no_duplicate_rendered_sections "${log}" "${rel} --help"
        assert_no_rendered_todo "${log}" "${rel} --help"
        assert_no_rendered_helper_inventory "${log}" "${rel} --help"
        assert_rendered_help_structure "${log}" "${rel} --help"

        if [[ "${rel}" == "bin/combine_parts_scaling_factor.sh" ]]; then
            assert_combine_parts_examples_fixture "${log}"
        fi

        warn_help_pattern_missing \
            "${log}" \
            '^Usage$' \
            "${rel} help has Usage"

        warn_help_pattern_missing \
            "${log}" \
            '^Parameters$' \
            "${rel} help has Parameters"

        warn_help_pattern_missing \
            "${log}" \
            '^Notes$' \
            "${rel} help has Notes"

        case "${rel}" in
            bin/execute_*.sh)
                assert_pattern_found \
                    "${log}" \
                    '--env_nam' \
                    "${rel} help documents --env_nam"
                ;;
        esac
    else
        record_fail "--help failed for ${rel}; see $(print_relpath "${log}")"
        continue
    fi

    if \
        grep -q -- '--details' "${log}"
    then
        dlog="${TEST_DIR_LOG}/help/${rel//\//__}.details.log"
        if \
            run_capture \
                "details ${rel}" \
                "${dlog}" \
                "${TEST_BASH}" "${file}" --details
        then
            record_pass "--details exits 0 for ${rel}"
            assert_no_duplicate_rendered_sections "${dlog}" "${rel} --details"
            assert_no_rendered_todo "${dlog}" "${rel} --details"
            assert_no_rendered_helper_inventory "${dlog}" "${rel} --details"
            assert_rendered_help_structure "${dlog}" "${rel} --details"

            warn_help_pattern_missing \
                "${dlog}" \
                '^Usage$' \
                "${rel} details has Usage"

            warn_help_pattern_missing \
                "${dlog}" \
                '^Parameters$' \
                "${rel} details has Parameters"

            warn_help_pattern_missing \
                "${dlog}" \
                '^Notes$' \
                "${rel} details has Notes"

            warn_help_pattern_missing \
                "${dlog}" \
                '^Examples$' \
                "${rel} details has Examples"
        else
            record_warn \
                "--details advertised but failed for ${rel}; see" \
                "$(print_relpath "${dlog}")"
        fi
    else
        record_skip "no detailed help flag advertised for ${rel}"
    fi
done < <(
    find "${ROOT_REPO}/bin" -maxdepth 1 -type f -name '*.sh' -print | sort
)

finish
