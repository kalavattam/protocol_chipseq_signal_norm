#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: make.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


# Require Bash >= 4.4 before doing any work.
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
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

# Run in safe mode, exiting on errors, unset variables, and pipe failures.
set -euo pipefail


# Resolve paths relative to 'tests/fixtures'
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_fix="${dir_scr}"

# Source shared fixture-generation helpers.
# shellcheck source=tests/support/fixture_helpers.sh
source "${dir_scr}/../../support/fixture_helpers.sh"

# Declare every generated path up front. The directory names the verdict the
# checker must return for the document inside it: the canonical rendering is
# the only conforming input, and every other document reports at least one
# finding.
dir_acc="${dir_fix}/accepted"
dir_rej="${dir_fix}/rejected"
fil_canonical="${dir_acc}/canonical.json"
fil_overflow="${dir_rej}/inline_overflow.json"
fil_fits="${dir_rej}/expanded_fits.json"
fil_hybrid="${dir_rej}/hybrid_delimiter.json"
fil_indent="${dir_rej}/wrong_indent.json"
fil_tab="${dir_rej}/tab_indent.json"
fil_newline="${dir_rej}/no_trailing_newline.json"
fil_duplicate="${dir_rej}/duplicate_key.json"
fil_unreadable="${dir_rej}/unreadable.json"

# Remove stale outputs so regeneration is idempotent.
rm_files "${dir_fix}" \
    "${fil_canonical}" \
    "${fil_overflow}" \
    "${fil_fits}" \
    "${fil_hybrid}" \
    "${fil_indent}" \
    "${fil_tab}" \
    "${fil_newline}" \
    "${fil_duplicate}" \
    "${fil_unreadable}"

mkdirs "${dir_acc}" "${dir_rej}"


# Author every fixture literally. Each negative fixture departs from the
# canonical rendering along exactly one axis, so a finding about that axis
# cannot be confounded with a finding about another. The delimiter is quoted
# so no '$' in a fixture body reaches the shell.

# The canonical rendering itself, exercising both treatments the budget
# selects: a short structure stays inline, a long one expands, and a
# record-per-line array whose rows fit is preserved exactly as authored.
cat << 'EOM' > "${fil_canonical}"
{
  "schema_version": 1,
  "flags": {"strict": true, "verbose": false},
  "empty_object": {},
  "empty_array": [],
  "commands": [
    {"callable": "atria", "conceptual_names": ["Atria"]},
    {"callable": "bowtie2", "conceptual_names": ["Bowtie 2"]}
  ],
  "long_prose": "A scalar longer than the budget is not this rule's business, because JSON has no string continuation.",
  "expanded_because_it_does_not_fit": [
    "first entry that helps push this array past the budget",
    "second entry that helps push this array past the budget"
  ]
}
EOM

# One array packed onto its key's line, past the budget.
cat << 'EOM' > "${fil_overflow}"
{
  "supplied_artifacts": ["semantic_review/report.md", "findings.ndjson", "facts.ndjson", "limitations.ndjson"]
}
EOM

# One structure broken across lines that fits the budget inline.
cat << 'EOM' > "${fil_fits}"
{
  "schema_version": {
    "const": 2
  }
}
EOM

# One object opened inline and then continued vertically. The opening
# delimiter is not the last content on its line and the closing delimiter is
# not the first content on its line.
cat << 'EOM' > "${fil_hybrid}"
{
  "requirements": [
    {"name": "wrapper_roots", "edges": [
      {"from": "bin/execute.sh", "to": "lib/bash/core/source_helpers.sh"}
    ]}
  ]
}
EOM

# One expanded structure indented by three spaces instead of two.
cat << 'EOM' > "${fil_indent}"
{
   "paths": [
      "first entry that helps push this array past the column budget",
      "second entry that helps push this array past the column budget"
   ]
}
EOM

# One expanded structure indented with a tab. The tab is built with 'printf'
# rather than typed into a heredoc body, where it would be an invisible
# control character that any editor could silently convert to spaces.
tab="$(printf '\t')"

{
    echo '{'
    echo "${tab}\"paths\": ["
    echo "${tab}${tab}\"first entry pushing this array past the budget\","
    echo "${tab}${tab}\"second entry pushing this array past the budget\""
    echo "${tab}]"
    echo '}'
} > "${fil_tab}"

# One canonical document with its final newline removed. 'printf' is required
# because a heredoc always terminates its body with a newline.
printf '%s' \
    '{"schema_version": 1}' \
    > "${fil_newline}"

# One object declaring the same key twice. Rewriting this file would silently
# discard the first value, so the checker must refuse to read it.
cat << 'EOM' > "${fil_duplicate}"
{
  "schema_version": 1,
  "schema_version": 2
}
EOM

# One document that is not JSON at all.
cat << 'EOM' > "${fil_unreadable}"
{
  "schema_version": 1,
}
EOM


succeed "generated JSON source-form fixtures under ${dir_fix}"
