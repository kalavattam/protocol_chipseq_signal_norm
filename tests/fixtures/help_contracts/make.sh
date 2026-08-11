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

# Run in safe mode, exiting on errors, unset variables, and pipe failures.
set -euo pipefail


# Resolve paths relative to 'tests/fixtures'.
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_fix="${dir_scr}"

# Source shared fixture-generation helpers.
# shellcheck source=tests/support/fixture_helpers.sh
source "${dir_scr}/../../support/fixture_helpers.sh"

# Declare every generated path up front. 'cases.json' sits at the fixture root
# rather than under a verdict directory because it is the record of every
# verdict at once, not an input carrying one of them.
dir_acc="${dir_fix}/accepted"
fil_cases="${dir_fix}/cases.json"
fil_acc_shell="${dir_acc}/shell.sh"
fil_acc_python="${dir_acc}/python.py"
fil_acc_callable="${dir_acc}/python_callable.py"


# Remove stale outputs so regeneration is idempotent.
rm_files "${dir_fix}" \
    "${fil_cases}" \
    "${fil_acc_shell}" \
    "${fil_acc_python}" \
    "${fil_acc_callable}"

mkdirs "${dir_acc}"


# Author every fixture literally. The delimiter is quoted, so no '$' or
# backtick in a fixture body reaches the shell.

# One shell help surface. The fixture is itself a heredoc-bearing function, so
# its last two lines are its own 'EOM' terminator and closing brace. Those two
# lines are appended rather than written inside the heredoc below, because a
# body line reading 'EOM' would end this heredoc at the fixture's terminator.
# Switching this recipe to an 'EOF' delimiter would avoid that too, and is what
# 'SHELL.SOURCE.FORM' asks shell-text heredocs not to do.
cat << 'EOM' > "${fil_acc_shell}"
function demo() {
    cat << EOM
Usage
-----
  demo [--dry_run]

Parameters
----------
  -dr, --dry_run, --dry_run : flag
    Run in dry-run mode.
EOM

printf '%s\n' 'EOM' '}' >> "${fil_acc_shell}"

# Two Python help surfaces: one 'add_argument' help literal and one callable
# whose Examples section is source-language rather than a CLI invocation.
cat << 'EOM' > "${fil_acc_python}"
def build(parser):
    parser.add_argument(
        "--mode",
        help="Use `--mode dist` for distribution mode.",
    )
EOM

cat << 'EOM' > "${fil_acc_callable}"
def combine(left: float, right: float) -> tuple[float, float]:
    """
    Combine two values.

    Examples
    --------
    >>> combine(1.0, 3.0)
    (2.0, 2.0)

    >>> combine(2.0, 4.0)
    (3.0, 3.0)
    """

    return (left + right) / 2.0, (left + right) / 2.0
EOM

# The verdict record. It enumerates the accepted, rejected, boundary, and
# non-applicable readings the help-contract rules must produce, and it is
# written in the canonical form owned by 'JSON.SOURCE.FORM'.
cat << 'EOM' > "${fil_cases}"
{
  "schema_version": 3,
  "accepted": {
    "audience": "Installed help contains the runtime formula and conditions.",
    "token": "Use '--mode dist --method qntl_nz' for the default distribution method.",
    "parameter": "Maximum number of decimal places retained for finite emitted values.",
    "examples": 2,
    "usage_rows": "Registered current Usage rows preserve every approved group.",
    "rendered_cli": {
      "example_form": "rendered_cli_invocation",
      "representation_owner": "HELP.EXAMPLES",
      "lifecycle": "active_enforced"
    },
    "callable": {
      "example_form": "callable_source_language",
      "representation_owner": "PY.DOCSTRING.NUMPY",
      "lifecycle": "active_enforced",
      "authority": "authorized_fresh_pilot",
      "fingerprint_scope": "source_and_displayed_expected_result"
    },
    "reporter_binding": [{"source": "parsed_args"}, {"literal": null}]
  },
  "rejected": {
    "audience": "See the repository design document for required mode behavior.",
    "token": "Use `--mode dist --method qntl_nz` for the default distribution method.",
    "parameter": "Precision.",
    "examples": 0,
    "usage_rows": "Merged, split, or misordered registered Usage groups are rejected.",
    "callable": {
      "example_form": "callable_source_language",
      "representation_owner": "PY.DOCSTRING.NUMPY",
      "lifecycle": "deferred_migration",
      "deferred_record": "S3-MIG-001",
      "required_count": 2,
      "current_count": 0
    },
    "reporter_binding": {"source": "environment_variable"}
  },
  "boundary": {
    "parameter": "Reference FASTA file. Required for CRAM decoding and ignored otherwise.",
    "examples": 1,
    "token": "Ambiguous fragments require review.",
    "bash_block_meaning": "CLI invocation environment, not implementation language."
  },
  "non_applicable": {
    "token": "compute_input_floor --mode dist --method qntl_nz",
    "context": "structured_usage",
    "function": "is_even",
    "private_renderer": "not_public_api"
  }
}
EOM

succeed "generated help-contract fixtures under ${dir_fix}"
