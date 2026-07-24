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


"""Tests for the deliberately narrow shell source-form checker."""

from __future__ import annotations

from dev.audit.shell_source_form import check_text


def messages(text: str, path: str = "tool.sh") -> list[str]:
    """Return stable messages for one in-memory shell source."""

    return [finding.message for finding in check_text(text, path)]


def test_accepts_recognized_naming_indentation_and_heredocs() -> None:
    """Canonical fixture-defined forms pass the bounded checker."""

    text = """function parse_args() {
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
    """The checker reports only its declared simple source forms."""

    text = """function ParseArgs() {
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
    assert "ordinary shell-text heredoc should use EOM instead of EOF" in output


def test_rejects_local_array_without_arr_prefix() -> None:
    """Recognize the claimed simple local -a declaration boundary."""

    assert messages("local -a values=()\n") == [
        "recognized array declaration should use the arr_ prefix"
    ]


def test_excludes_posix_bootstrap_from_bash_only_forms() -> None:
    """The explicit POSIX bootstrap is not treated as maintained Bash."""

    assert messages(
        "ParseArgs() {\n  :\n}\n",
        "install/scripts/install_envs_entrypoint.sh",
    ) == []
