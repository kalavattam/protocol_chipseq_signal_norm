#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_shell_source_form.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5, Opus 5.5).
#
# Distributed under the MIT license.


"""
Tests for the deliberately narrow shell source-form checker.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from dev.audit.shell_source_form import (
    DIAGNOSTIC_RULE_ID,
    RULE_COMMENT,
    check_text,
)

ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tests" / "fixtures" / "shell_source_form"


def messages(text: str, path: str = "tool.sh") -> list[str]:
    """
    Return stable messages for one in-memory shell source.
    """

    return [finding.message for finding in check_text(text, path)]


def diagnostic_messages(text: str, path: str = "tool.sh") -> list[str]:
    """
    Return only direct-echo diagnostic-form messages.
    """

    return [
        finding.message
        for finding in check_text(text, path)
        if finding.rule_id == DIAGNOSTIC_RULE_ID
    ]


def test_accepts_recognized_naming_indentation_and_heredocs() -> None:
    """
    Canonical fixture-defined forms pass the bounded checker.
    """

    text = """
function parse_args() {
    local -a arr_values=()
    declare -ra arr_readonly=()
    cat << 'EOM'
literal text
EOM
    python << PY
print('ok')
PY
}
"""

    assert messages(text) == []


def test_rejects_recognized_noncanonical_forms() -> None:
    """
    The checker reports only its declared simple source forms.
    """

    text = """
function ParseArgs() {
  declare -a values=()
\tcat << EOF
literal text
EOF
}
"""

    output = messages(text)

    assert "recognized function name is not snake_case" in output
    assert "recognized array declaration should use the arr_ prefix" in output
    assert "leading indentation contains a tab" in output
    assert (
        "ordinary shell-text heredoc should use EOM instead of EOF" in output
    )


def test_function_names_allow_one_leading_underscore() -> None:
    """
    One leading underscore marks a private helper; two, or a capital after it,
    are rejected.
    """

    accepted = """
function _sort_qname_bam() {
    :
}
"""
    rejected = """
function __helper() {
    :
}

function _Helper() {
    :
}
"""

    assert messages(accepted) == []
    assert messages(rejected) == [
        "recognized function name is not snake_case",
        "recognized function name is not snake_case",
    ]


def test_requires_exactly_one_space_before_a_heredoc_delimiter() -> None:
    """
    Reject a missing or repeated heredoc operator separator.
    """

    joined = "\npython - <<PY\nprint('ok')\nPY\n"
    quoted = "\npython - <<'PY'\nprint('ok')\nPY\n"
    widened = "\npython - <<  PY\nprint('ok')\nPY\n"
    indented = "\npython - <<- PY\nprint('ok')\nPY\n"
    separator = (
        "heredoc operator and delimiter must be separated by exactly one space"
    )

    assert messages(joined) == [separator]
    assert messages(quoted) == [separator]
    assert messages(widened) == [separator]
    assert messages(indented) == []


def test_heredoc_separator_excludes_herestrings_and_shifts() -> None:
    """
    Keep herestrings and arithmetic shifts outside the separator facet.
    """

    herestring = '\nIFS="," read -r -a arr_values <<< "${csv_values}"\n'
    bare_herestring = "\nread -r value <<<word\n"
    spaced_shift = "\nvalue=$(( total << width ))\n"
    joined_shift = "\nvalue=$((total<<width))\n"
    guarded_shift = "\nif (( total<<width > limit )); then\n    :\nfi\n"

    assert messages(herestring) == []
    assert messages(bare_herestring) == []
    assert messages(spaced_shift) == []
    assert messages(joined_shift) == []
    assert messages(guarded_shift) == []


def test_rejects_local_array_without_arr_prefix() -> None:
    """
    Recognize the claimed simple local -a declaration boundary.
    """

    assert messages("local -a values=()\n") == [
        "recognized array declaration should use the arr_ prefix",
    ]


def test_excludes_posix_bootstrap_from_bash_only_forms() -> None:
    """
    The explicit POSIX bootstrap is not treated as maintained Bash.
    """

    observed = messages(
        "ParseArgs() {\n  :\n}\n",
        "install/scripts/install_envs_entrypoint.sh",
    )

    assert observed == []


def test_accepts_split_direct_echo_diagnostics() -> None:
    """
    Prefix-only first arguments and bounded prose continuations pass.
    """

    text = """
if [[ -z "${value}" ]]; then
    echo "error(test_helpers.sh):" \\
        "TEST_ARTIFACT_ROOT must be an absolute non-root path. This is" \\
        "additional text split across adjacent arguments." >&2
fi
"""

    assert diagnostic_messages(text) == []


def test_rejects_message_joined_to_diagnostic_prefix() -> None:
    """
    A diagnostic prefix and prose may not share the first argument.
    """

    text = (
        'echo "error(test_helpers.sh): TEST_ARTIFACT_ROOT must be an absolute '
        'non-root path." >&2\n'
    )

    assert diagnostic_messages(text) == [
        "recognized diagnostic source line exceeds 79 columns",
        "diagnostic prefix must be the first quoted argument by itself",
    ]


def test_rejects_bad_continuation_shape_and_width() -> None:
    """
    Continuation indentation, quoting, edge spaces, and width are bounded.
    """

    overlong = "x" * 80
    text = (
        'echo "warning(validate):" \\\n'
        '  " leading space" \\\n'
        f'    "{overlong}" trailing\n'
    )

    output = diagnostic_messages(text)

    continuation = (
        "diagnostic continuation must be indented one additional level"
    )
    nonempty = (
        "diagnostic prose argument must be nonempty without edge whitespace"
    )

    assert continuation in output
    assert nonempty in output
    assert "recognized diagnostic source line exceeds 79 columns" in output
    assert "diagnostic final line has an unrecognized source suffix" in output


def test_diagnostic_check_is_stable_and_ignores_literal_text() -> None:
    """
    Repeated scans are stable and heredoc-like text is not a direct command.
    """

    text = """
cat << 'EOM'
echo "error(literal): this is expected output."
EOM
value='echo "error(data): not executable source."'
printf 'error(%s): %s\n' "${name}" "${message}" >&2
echo_err "error(helper): helper-owned output."
"""
    first = check_text(text)

    assert first == check_text(text)
    assert [
        finding for finding in first if finding.rule_id == DIAGNOSTIC_RULE_ID
    ] == []


def test_rejects_echo_split_before_the_diagnostic_prefix() -> None:
    """
    The echo command and prefix belong together on the first source line.
    """

    text = """
echo \\
    "error(run_python.sh): command failed." >&2
"""

    assert diagnostic_messages(text) == [
        "diagnostic echo and prefix must share the first source line",
        "diagnostic prefix must be the first quoted argument by itself",
    ]


def test_diagnostic_width_boundary_is_79_columns() -> None:
    """
    A 79-column continuation passes and an 80-column continuation fails.
    """

    pass_line = f'    "{"x" * 69}" >&2'
    fail_line = f'    "{"x" * 70}" >&2'

    assert len(pass_line) == 79
    assert len(fail_line) == 80
    assert (
        diagnostic_messages(f'echo "debug(worker):" \\\n{pass_line}\n') == []
    )
    assert diagnostic_messages(f'echo "debug(worker):" \\\n{fail_line}\n') == [
        "recognized diagnostic source line exceeds 79 columns",
    ]


def test_rejects_empty_context_and_incomplete_continuation() -> None:
    """
    A context is required and a continued prefix needs diagnostic prose.
    """

    assert diagnostic_messages('echo "error(): message." >&2\n') == [
        "diagnostic script-or-function context must be nonempty",
    ]
    assert diagnostic_messages('echo "note(worker):" \\\n') == [
        "diagnostic continuation is incomplete",
    ]


def test_rejects_single_quoted_diagnostic_prefix() -> None:
    """
    Recognized direct diagnostics use the governed double-quoted form.
    """

    assert diagnostic_messages("echo 'error(worker): message.' >&2\n") == [
        "diagnostic prefix must use a double-quoted first argument",
    ]


def test_split_echo_preserves_rendered_spacing() -> None:
    """
    Adjacent quoted arguments preserve the prior visible diagnostic.
    """

    original = subprocess.run(
        ["bash", "-c", 'echo "error(worker): message text." >&2'],
        check=False,
        capture_output=True,
        text=True,
    )
    split = subprocess.run(
        [
            "bash",
            "-c",
            'echo "error(worker):" \\\n    "message text." >&2',
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert original.returncode == split.returncode == 0
    assert original.stdout == split.stdout == ""
    assert original.stderr == split.stderr


def test_rejects_superseded_ordinary_comment_markers() -> None:
    """
    Require '# ' for new ordinary comments.

    'SHELL.COMMENT.FORM' records '#  ' and '#+ ' as superseded rather than as
    compatibility forms, so source written now may not use them.
    """

    header = "\n" * 16
    accepted = f"{header}# One ordinary comment.\nvalue=1\n"
    wide = f"{header}#  One ordinary comment.\nvalue=1\n"
    continued = (
        f"{header}# One ordinary comment,\n#+ continued here.\nvalue=1\n"
    )
    expected = (
        "ordinary comment must begin with '# '; the '#  ' and '#+ ' forms "
        "are superseded"
    )

    assert expected not in messages(accepted)
    assert expected in messages(wide)
    assert expected in messages(continued)


def test_superseded_marker_check_spares_owned_syntax() -> None:
    """
    Leave the bounded header, directives, and heredoc bodies alone.

    Each carries separately owned syntax, so a retired ordinary marker inside
    them is not an ordinary comment at all.
    """

    header = "\n" * 16
    directive = f"{header}#  shellcheck disable=SC2034\nvalue=1\n"
    in_header = "#  Copyright row spacing.\n" + "\n" * 16 + "value=1\n"
    heredoc = f"{header}cat << 'EOM'\n#  Literal body row.\nEOM\n"
    expected = (
        "ordinary comment must begin with '# '; the '#  ' and '#+ ' forms "
        "are superseded"
    )

    assert expected not in messages(directive)
    assert expected not in messages(in_header)
    assert expected not in messages(heredoc)


def wrap_lines(text: str, path: str) -> list[int]:
    """
    Return the lines that the comment-wrap facet reports.
    """

    return [
        finding.line
        for finding in check_text(text, path)
        if finding.rule_id == RULE_COMMENT
        and finding.message.startswith("comment breaks before")
    ]


def test_comment_wrap_fixtures_accept_and_bound_greedy_prose() -> None:
    """
    Report nothing for greedy prose and for breaks on the fit-test edge.

    The accepted script holds every structural boundary the facet recognizes,
    and the boundary script holds breaks that are one column or one unit too
    long to move, so neither may produce any finding.
    """

    for name in ("accepted/comment_wrap.sh", "boundary/comment_wrap.sh"):
        text = (FIXTURES / name).read_text(encoding="utf-8")

        assert check_text(text, f"tests/{name}") == []


def test_comment_wrap_fixture_rejects_one_break_per_block() -> None:
    """
    Report exactly one premature break in each rejected block.

    Blocks are separated by blank lines. Each holds one break, so a facet that
    misses a case or reports a structural boundary changes the count of some
    block.
    """

    name = "rejected/comment_wrap.sh"
    text = (FIXTURES / name).read_text(encoding="utf-8")
    lines = text.splitlines()
    starts = [
        number
        for number, line in enumerate(lines, 1)
        if number > 16 and line.strip() and not lines[number - 2].strip()
    ]
    found = wrap_lines(text, f"tests/{name}")
    per_block = [
        sum(start <= line < end for line in found)
        for start, end in zip(
            starts,
            [*starts[1:], len(lines) + 1],
            strict=True,
        )
    ]

    assert per_block == [1, 1, 1, 1, 1]
