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
shell diagnostic form	${doc_shell}	SHELL.DIAGNOSTIC.FORM
submit bootstrap	${doc_shell}	SHELL.SUBMIT.BOOTSTRAP
Python CLI parser	${doc_python}	PY.CLI.PARSER
Python CLI help layout	${doc_python}	PY.CLI.HELP.LAYOUT
Python callable annotations	${doc_python}	PY.TYPE.ANNOTATIONS
Python exact docstring layout	${doc_python}	PY.DOCSTRING.LAYOUT
Python NumPy docstrings	${doc_python}	PY.DOCSTRING.NUMPY
Python ordinary strings	${doc_python}	PY.STRING.QUOTES
Python comments	${doc_python}	PY.COMMENT.FORM
Python identifiers	${doc_python}	PY.NAMING.IDENTIFIERS
Python naming length	${doc_python}	PY.NAMING.LENGTH
Python source layout	${doc_python}	PY.SOURCE.LAYOUT
source layout candidates	${doc_source_layout}	SOURCE.LAYOUT.CANDIDATES
multiline delimited structures	${doc_source_layout}	SOURCE.DELIMITED.MULTILINE
EOM

    if (( found == 0 )); then
        record_pass "standards contain the required stable owner IDs"
    fi
}


#  Fail when the task router loses language routing, owners, or commands
function scan_python_task_router() {
    local found=0
    local required

    for required in \
        "Python implementation-body layout" \
        "Python implementation comments" \
        "Python implementation naming" \
        "Shell implementation-body layout, comments, or naming" \
        "R implementation-body layout, comments, or naming" \
        "Rust implementation-body layout, comments, or naming" \
        "New or revised obligation, vocabulary entry, checker facet, or rule ID" \
        "Python \`parse_args()\` help literal" \
        "Python docstring or API" \
        "Python CLI error handling" \
        "source_layout.md#semantic-paragraphs-and-density-sourcelayoutparagraphs" \
        "source_layout.md#comment-attachment-sourcecommentattachment" \
        "source_layout.md#grammatical-naming-and-migration-sourcenamingsemantics" \
        "shell.md#shell-source-form-shellsourceform" \
        "r.md#r-identifiers-rnamingidentifiers" \
        "rust.md#rust-identifiers-rustnamingidentifiers" \
        "governance.md#authoritative-standard-first-changes-govchangegoldenfirst" \
        "python.md#cli-help-literal-layout-pyclihelplayout" \
        "help.md#shared-python-docstring-schema-helpdocstringschema" \
        "python.md#error-and-status-ownership-pyerrorexit" \
        "--rule PY.SOURCE.LAYOUT <paths>" \
        "--rule PY.COMMENT.FORM <paths>" \
        "--rule PY.NAMING.IDENTIFIERS <paths>" \
        "--rule PY.CLI.HELP.LAYOUT <paths>" \
        "--rule PY.DOCSTRING.LAYOUT <paths>" \
        "if no language owner exists, route the proposed specialization" \
        "instead of borrowing Python rules" \
        "record the anti-accretion evidence"; do
        if ! grep -Fq -- "${required}" "${doc_router}"; then
            found=1
            record_fail "ordinary task router missing: ${required}"
        fi
    done

    if (( found == 0 )); then
        record_pass "ordinary Python task router retains owners and commands"
    fi
}


#  Fail when shared comment or delimiter semantics lack language realizations
function scan_language_layout_realizations() {
    local file found=0 label pattern

    while IFS=$'\t' read -r label file pattern; do
        [[ -n "${label}" ]] || continue

        if ! grep -Fq -- "${pattern}" "${file}"; then
            found=1
            record_fail "${label} missing language-owned text: ${pattern}"
        fi
    done << EOM
Python comment markers	${doc_python}	beginning each nonempty ordinary full-line or continuation comment
Python trailing comma	${doc_python}	trailing comma in expanded calls
Shell comment markers	${doc_shell}	Begin each nonempty ordinary comment line
Shell no trailing comma	${doc_shell}	do not acquire a trailing comma
R comment markers	${doc_r}	Ordinary R comments use
R documentation markers	${doc_r}	for roxygen2 documentation
R no trailing comma	${doc_r}	Multiline R calls do not use
Rust comment markers	${doc_rust}	Ordinary full-line, continuation, and inline comments
Rust documentation markers	${doc_rust}	outer documentation attached to the following item
Rust trailing comma	${doc_rust}	and a trailing comma
EOM

    if (( found == 0 )); then
        record_pass "language owners retain comment and delimiter realizations"
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


#  Fail when bounded public module docstrings retain broken private-doc paths
function scan_broken_python_see_also() {
    local found=0
    local file

    for file in "${python_see_also_modules[@]}"; do
        if grep -Fq -- "docs/dev/" "${file}"; then
            found=1
            record_fail \
                "public Python module retains a broken docs/dev reference:" \
                "$(print_relpath "${file}")"
        fi
    done

    if (( found == 0 )); then
        record_pass "bounded public Python docstrings omit broken docs/dev paths"
    fi
}


print_section "${TEST_NAME}"

doc_help="${ROOT_REPO}/docs/standards/help.md"
doc_python="${ROOT_REPO}/docs/standards/python.md"
doc_r="${ROOT_REPO}/docs/standards/r.md"
doc_router="${ROOT_REPO}/docs/standards/README.md"
doc_rust="${ROOT_REPO}/docs/standards/rust.md"
doc_shell="${ROOT_REPO}/docs/standards/shell.md"
doc_source_layout="${ROOT_REPO}/docs/standards/source_layout.md"

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
    "${doc_r}"
    "${doc_router}"
    "${doc_rust}"
    "${doc_shell}"
    "${doc_source_layout}"
    "${chk_help}"
    "${chk_python}"
    "${chk_wrapper}"
    "${chk_smoke}"
    "${chk_docs}"
    "${chk_doc_cov}"
    "${chk_param_docs}"
)

dir_cli="${ROOT_REPO}/src/protocol_chipseq_signal_norm/cli"
python_see_also_modules=(
    "${dir_cli}/calculate_scaling_factor_siqchip.py"
    "${dir_cli}/calculate_scaling_factor_spike.py"
    "${dir_cli}/compute_input_floor.py"
    "${dir_cli}/compute_pseudo.py"
    "${dir_cli}/compute_signal.py"
    "${dir_cli}/compute_signal_ratio.py"
)

scan_required_files
scan_owner_ids
scan_python_task_router
scan_language_layout_realizations
scan_retired_todo_notebook
scan_broken_python_see_also

finish
