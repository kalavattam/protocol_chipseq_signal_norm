#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_unknown_option_helpers.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Focused regressions for context-correct unknown-option diagnostics.
"""

from __future__ import annotations

import unittest
from pathlib import Path

from dev.audit.unknown_option_helpers import check_text, scan_repository

REPO_ROOT = Path(__file__).resolve().parents[3]

FUNCTION_HELP = """
function demo() {{
    local show_help
    show_help=$(cat << EOM
Usage
-----
  demo [--help]
Returns
-------
  Returns 0 or 1.
EOM
    )
    case "${{1:-}}" in
        *)
            {error}
            echo >&2
            echo "${{show_help}}" >&2
            return 1
            ;;
    esac
}}
"""


SCRIPT_MAIN = """
function main() {{
    case "${{1:-}}" in
        *)
            {error}
            help_demo >&2
            return 1
            ;;
    esac
}}
main "$@"
"""


class UnknownOptionHelpersTest(unittest.TestCase):
    """
    Distinguish reusable keyword interfaces from top-level CLIs.
    """

    def test_library_function_using_echo_err_func_passes(self) -> None:
        text = FUNCTION_HELP.format(
            error=(
                'echo_err_func "${FUNCNAME[0]}" '
                "\"unknown option/parameter passed: '${1}'.\""
            ),
        )

        self.assertEqual(check_text("lib/bash/core/demo.sh", text), [])

    def test_library_function_using_echo_err_fails(self) -> None:
        text = FUNCTION_HELP.format(
            error="echo_err \"unknown option/parameter passed: '${1}'.\"",
        )

        self.assertTrue(check_text("lib/bash/core/demo.sh", text))

    def test_top_level_main_uses_script_context(self) -> None:
        good = SCRIPT_MAIN.format(
            error="echo_err \"unknown option/parameter passed: '${1}'.\"",
        )
        bad = SCRIPT_MAIN.format(
            error=(
                'echo_err_func "${FUNCNAME[0]}" '
                "\"unknown option/parameter passed: '${1}'.\""
            ),
        )

        self.assertEqual(check_text("bin/demo.sh", good), [])
        self.assertTrue(check_text("bin/demo.sh", bad))

    def test_function_requires_funcname_and_exact_wording(self) -> None:
        missing_name = FUNCTION_HELP.format(
            error=(
                'echo_err_func "demo" '
                "\"unknown option/parameter passed: '${1}'.\""
            ),
        )
        wrong_words = FUNCTION_HELP.format(
            error=(
                'echo_err_func "${FUNCNAME[0]}" '
                "\"unknown parameter passed: '${1}'.\""
            ),
        )

        self.assertTrue(check_text("lib/bash/core/demo.sh", missing_name))
        self.assertTrue(check_text("lib/bash/core/demo.sh", wrong_words))

    def test_help_emission_and_return_are_required(self) -> None:
        missing_help = FUNCTION_HELP.replace(
            '            echo "${{show_help}}" >&2\n',
            "",
        )

        missing_status = FUNCTION_HELP.replace("            return 1\n", "")
        error = (
            'echo_err_func "${FUNCNAME[0]}" '
            "\"unknown option/parameter passed: '${1}'.\""
        )

        self.assertTrue(
            check_text(
                "lib/bash/core/demo.sh",
                missing_help.format(error=error),
            ),
        )
        self.assertTrue(
            check_text(
                "lib/bash/core/demo.sh",
                missing_status.format(error=error),
            ),
        )

    def test_current_repository_has_no_clear_helper_violations(self) -> None:
        findings, inventory = scan_repository(REPO_ROOT)

        self.assertTrue(inventory)
        self.assertEqual(findings, [])


if __name__ == "__main__":
    unittest.main()
