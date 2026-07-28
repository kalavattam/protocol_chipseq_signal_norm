#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: make.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


#  Require Bash >= 4.4 before doing any work
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(make.sh):" \
        "this script must be run under Bash >= 4.4." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

#  Run in safe mode, exiting on errors, unset variables, and pipe failures
set -euo pipefail


#  Resolve repository, fixture, and checker paths
dir_fix="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_root="$(
    git -C "${dir_fix}" rev-parse --show-toplevel
)"
chk_policy="dev.audit.python_source_policy"

#  Source shared fixture-generation helpers
# shellcheck source=tests/support/fixture_helpers.sh
source "${dir_root}/tests/support/fixture_helpers.sh"


#  Define the tracked fixture cohorts
arr_positive=(
    "${dir_fix}/positive.py.fixture"
    "${dir_fix}/exceptions.py.fixture"
    "${dir_fix}/boundary.py.fixture"
)
negative="${dir_fix}/negative.py.fixture"
arr_expected_rules=(
    "HELP.PROSE.SENTENCES"
    "PY.CLI.HELP.LAYOUT"
    "PY.COMMENT.FORM"
    "PY.DOCSTRING.LAYOUT"
    "PY.DOCSTRING.NUMPY"
    "PY.NAMING.IDENTIFIERS"
    "PY.SOURCE.LAYOUT"
    "PY.STRING.QUOTES"
    "PY.TYPE.ANNOTATIONS"
    "SOURCE.DELIMITED.MULTILINE"
)


#  Require the managed project environment and checker dependencies
require_env "env_protocol" "for Python source-policy fixtures."
require_cmd git "to resolve the fixture repository."
require_cmd python "to validate Python source-policy fixtures."

#  Prove that canonical and exception fixtures have no deterministic findings
(
    cd "${dir_root}"
    python -m "${chk_policy}" \
        --root "${dir_root}" \
        "${arr_positive[@]}"
)

#  Prove that the negative fixture fails with every expected owner
if output="$(
    cd "${dir_root}"
    python -m "${chk_policy}" \
        --root "${dir_root}" \
        "${negative}" 2>&1
)"; then
    die "negative Python source-policy fixture passed unexpectedly."
fi

for rule_id in "${arr_expected_rules[@]}"; do
    if [[ "${output}" != *"${rule_id}:"* ]]; then
        die "negative fixture did not emit expected owner '${rule_id}'."
    fi
done

succeed "validated tracked Python source-policy fixtures under ${dir_fix}"
