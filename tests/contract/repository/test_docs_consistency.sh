#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_docs_consistency.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="docs consistency"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


#  Fail when required style docs or focused style tests are missing
function scan_required_files() {
    local found=0
    local file

    for file in "${files_required[@]}"; do
        if [[ ! -f "${file}" ]]; then
            found=1
            record_fail \
                "required style file missing:" \
                "$(print_relpath "${file}")"
        fi
    done

    if (( found == 0 )); then
        record_pass "all required style docs and checks exist"
    fi
}


#  Fail when required owner IDs are missing from their maintained standards
function scan_owner_ids() {
    local found=0
    local file label pattern rule_id

    while IFS=$'\t' read -r label file rule_id; do
        [[ -n "${label}" ]] || continue
        printf -v pattern '`%s`' "${rule_id}"

        if ! grep -Fq -- "${pattern}" "${file}"; then
            found=1
            record_fail \
                "standard missing ${label} owner ${rule_id}:" \
                "$(print_relpath "${file}")"
        fi
    done << EOM
shared documentation schema	${doc_help}	HELP.SECTION.SCHEMA
public aliases	${doc_help}	HELP.ALIAS.PUBLIC
parameter vocabulary	${doc_help}	HELP.PARAMETER.VOCABULARY
canonical parameter descriptions	${doc_help}	PARAMETER.DESCRIPTIONS
documentation coverage	${doc_help}	DOC.COVERAGE
shell wrapper topology	${doc_shell}	SHELL.WRAPPER_TOPOLOGY
shell line length	${doc_shell}	SHELL.LINE_LENGTH
submit bootstrap	${doc_shell}	SHELL.SUBMIT.BOOTSTRAP
Python CLI parser	${doc_python}	PY.CLI.PARSER
Python NumPy docstrings	${doc_python}	PY.DOCSTRING.NUMPY
EOM

    if (( found == 0 )); then
        record_pass "standards contain the required stable owner IDs"
    fi
}


#  Fail if TODOs are centralized away from their associated code/docs
function scan_retired_todo_notebook() {
    local file="${ROOT_REPO}/docs/dev/maintainer_todos.md"

    if [[ -e "${file}" ]]; then
        record_fail \
            "retired centralized TODO notebook exists:" \
            "$(print_relpath "${file}")"
    else
        record_pass "centralized TODO notebook is absent"
    fi
}


print_section "${TEST_NAME}"

doc_help="${ROOT_REPO}/docs/standards/help.md"
doc_python="${ROOT_REPO}/docs/standards/python.md"
doc_shell="${ROOT_REPO}/docs/standards/shell.md"

chk_help="${ROOT_REPO}/tests/contract/repository/test_help_style.sh"
chk_python="${ROOT_REPO}/tests/contract/repository/test_python_style.sh"
chk_wrapper="${ROOT_REPO}/tests/contract/repository/test_shell_wrapper_style.sh"
chk_smoke="${ROOT_REPO}/tests/contract/repository/test_test_layout.sh"
chk_docs="${ROOT_REPO}/tests/contract/repository/test_docs_consistency.sh"
chk_doc_cov="${ROOT_REPO}/tests/contract/repository/test_doc_coverage.sh"
chk_param_docs="$(
    printf '%s' \
        "${ROOT_REPO}/tests/contract/repository/test_parameter_docs_consistency.sh"
)"

files_required=(
    "${doc_help}"
    "${doc_python}"
    "${doc_shell}"
    "${chk_help}"
    "${chk_python}"
    "${chk_wrapper}"
    "${chk_smoke}"
    "${chk_docs}"
    "${chk_doc_cov}"
    "${chk_param_docs}"
)

scan_required_files
scan_owner_ids
scan_retired_todo_notebook

finish
