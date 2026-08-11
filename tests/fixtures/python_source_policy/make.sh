#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: make.sh
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

# Declare every generated path up front. The directory names the verdict the
# checker must return for the module inside it; 'format/' holds the input and
# expected pairs for the help-literal producer, which claims a rewrite rather
# than a verdict.
dir_acc="${dir_fix}/accepted"
dir_bnd="${dir_fix}/boundary"
dir_rej="${dir_fix}/rejected"
dir_fmt="${dir_fix}/format"

fil_acc_canonical="${dir_acc}/canonical.py"
fil_acc_exceptions="${dir_acc}/exceptions.py"
fil_bnd_form="${dir_bnd}/source_form.py"
fil_rej_owners="${dir_rej}/deterministic_owners.py"
fil_fmt_input="${dir_fmt}/help_input.py"
fil_fmt_expected="${dir_fmt}/help_expected.py"
fil_fmt_uni_input="${dir_fmt}/help_unicode_input.py"
fil_fmt_uni_expected="${dir_fmt}/help_unicode_expected.py"


# Remove stale outputs so regeneration is idempotent.
rm_files "${dir_fix}" \
    "${fil_acc_canonical}" \
    "${fil_acc_exceptions}" \
    "${fil_bnd_form}" \
    "${fil_rej_owners}" \
    "${fil_fmt_input}" \
    "${fil_fmt_expected}" \
    "${fil_fmt_uni_input}" \
    "${fil_fmt_uni_expected}"

mkdirs "${dir_acc}" "${dir_bnd}" "${dir_rej}" "${dir_fmt}"


# Author every module literally. The delimiter is quoted, so no '$' in a
# fixture body reaches the shell.
#
# This recipe only writes. It does not run the checker over what it writes:
# 'make.sh' generates, and tests and checkers validate. The cohort claims this
# file used to assert now live in
# 'tests/unit/dev_audit/test_python_source_policy.py', where provoking an
# expected set of rule owners is a unit-level property rather than a
# precondition of having a fixture at all.

# Accepted: modules the checker must report nothing about
cat << 'EOM' > "${fil_acc_canonical}"
"""
Exercise canonical deterministic Python source forms.
"""

from __future__ import annotations

import argparse

PAIRS = (
    ("left", "right"),
    ("up", "down"),
)


def parse_args(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    """
    Parse fixture arguments.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        help=(
            "Write the selected result to one explicit output path while "
            "preserving the configured suffix."
        ),
    )

    return parser.parse_args(argv)


def result_pair() -> tuple[int, str]:
    """
    Return one canonical result pair.

    Returns
    -------
    result : tuple[
        int, str
    ]
        Canonical result values.
    """

    return 1, "value"


class Renderer:
    """
    Render one configured value.
    """

    def render(self, value: str) -> str:
        """
        Return the configured value.
        """

        # Preserve this value for the caller.
        result = [
            value,
        ]

        return result[0]
EOM

cat << 'EOM' > "${fil_acc_exceptions}"
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: exceptions.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Distributed under the MIT license.


r"""
Represent the literal \d regular-expression marker.
"""


def literal_values() -> tuple[str, str, str]:
    """
    Return bounded literal-content exceptions.
    """

    # noqa: RUF001
    quoted = 'The value contains "double quotes".'
    multiline = '''
The literal contains a """ delimiter.
'''
    pattern = r"\d+"

    return quoted, multiline, pattern
EOM

# Boundary: exact source-form edges — parameter markers, generated 79- and
# 80-column cases, and greedy docstring prose whose structural boundaries are
# not prose breaks. These must also report nothing.
cat << 'EOM' > "${fil_bnd_form}"
"""
Exercise exact deterministic source boundaries.
"""


def comments() -> None:
    """
    Exercise comment-width boundaries.
    """

    # This comment remains within the configured source-width boundary.
    value = 1  # This inline comment uses the canonical spacing form.
    _ = value


def parameters(
    first: int,
    /,
    second: int,
    *,
    third: int,
) -> None:
    """
    Exercise parameter-list delimiter boundaries.
    """

    return None


def prose_wrap(first: int, second: int) -> tuple[
    int,
    str,
]:
    """
    Exercise greedy docstring prose boundaries.

    A wrapped prose paragraph fills each physical line through the last
    complete word that fits, so the following word would reach column 80.

    Parameters
    ----------
    first : int
        Structural entry boundaries are not prose breaks, so this description
        never joins the entry header above it.
    second : int
        A bullet list wraps against its own continuation indentation:
        - the marker establishes a deeper continuation column, and the wrapped
          remainder aligns beneath the first word rather than the marker.

    Returns
    -------
    result : tuple[
        int,
        str,
    ]
        A multiline textual type is structure rather than prose.

    Examples
    --------
    >>> prose_wrap(1, 2)
    (1, 'value')
    """

    return first, "value"
EOM

# Rejected: bounded near misses for every deterministic owner. The module is
# deliberately malformed, and it is invisible to the checker because it is
# ignored, so it cannot be mistaken for maintained source violating the rules
# it exists to provoke.
cat << 'EOM' > "${fil_rej_owners}"
"""Reject noncanonical deterministic Python source forms."""

import argparse

def parse_args(argv):
    '''Parse bad fixture arguments'''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", "--outfile",
        help=(
            "Break this help prose "
            "far too early for greedy wrapping."
        )
    )
    parser.add_argument(
        "--static",
        help="This single static help literal is deliberately much too long for the governed source-width boundary.",
    )
    parser.add_argument(
        "--unseparated",
        help=(
            "Adjacent literals that lose the space the rendered prose"
            "requires, and a bracket pulled against the previous word"
            "(like this)."
        ),
    )
    parser.add_argument("--incomplete", help="incomplete help text")
    #  Bad marker.
    value = 'plain'
    return parser.parse_args(argv)


def BadName() -> None:
    """
    Exercise a noncanonical callable name.
    """

    return None


def bad_parameters(
    first: int, second: int
  ) -> None:
    """
        Exercise bad parameter and docstring indentation.
        """

    return None


def premature_prose() -> None:
    """
    Exercise insufficiently greedy docstring prose.

    Notes
    -----
    This sentence deliberately stops well before the governed width
    boundary.
    """

    return None


def bad_type_closer() -> tuple[int, str]:
    """
    Return one badly aligned result pair.

    Returns
    -------
    result : tuple[
        int, str
        ]
        Result values.
    """

    return 1, "value"
EOM

# Format: input and expected pairs for the help-literal producer, including the
# width-sensitive non-ASCII pair.
cat << 'EOM' > "${fil_fmt_input}"
def build(parser):
    parser.add_argument(
        "--value",
        help=(
            "This sentence uses complete words through "
            "the available column limit.\n"
        ),
    )
EOM

cat << 'EOM' > "${fil_fmt_expected}"
def build(parser):
    parser.add_argument(
        "--value",
        help=(
            "This sentence uses complete words through the available column "
            "limit.\n"
        ),
    )
EOM

cat << 'EOM' > "${fil_fmt_uni_input}"
import argparse


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--value",
        help=(
            "Unicode α β γ δ ε ζ η θ ι κ λ μ ν ξ ο π ρ σ τ υ "
            "must remain available for downstream summary output.\n"
        ),
    )
    return parser.parse_args(argv)
EOM

cat << 'EOM' > "${fil_fmt_uni_expected}"
import argparse


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--value",
        help=(
            "Unicode α β γ δ ε ζ η θ ι κ λ μ ν ξ ο π ρ σ τ υ must remain "
            "available for downstream summary output.\n"
        ),
    )
    return parser.parse_args(argv)
EOM


succeed "generated Python source-policy fixtures under ${dir_fix}"
