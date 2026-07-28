#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_shell_source_form.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""
Tests for the deliberately narrow shell source-form checker.
"""

from __future__ import annotations

import subprocess

from dev.audit.shell_source_form import DIAGNOSTIC_RULE_ID, check_text


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
