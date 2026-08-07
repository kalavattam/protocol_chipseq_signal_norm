#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_cli.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


import pytest

from protocol_chipseq_signal_norm.utilities import utils_cli
from protocol_chipseq_signal_norm.utilities.utils_cli import (
    CapArgumentParser,
    add_help_cap,
)


def test_cap_argument_parser_uses_capitalized_usage(
    capsys: pytest.CaptureFixture[str],
) -> None:
    parser = CapArgumentParser(prog="tool", description="Demo.")
    add_help_cap(parser)

    with pytest.raises(SystemExit) as excinfo:
        parser.parse_args(["--help"])

    assert excinfo.value.code == 0
    assert "Usage:" in capsys.readouterr().out


def test_cap_argument_parser_disables_abbreviation() -> None:
    parser = CapArgumentParser(prog="tool")
    parser.add_argument("--sample-name")

    with pytest.raises(SystemExit):
        parser.parse_args(["--sample", "x"])


def test_default_help_output_remains_exact() -> None:
    parser = CapArgumentParser(prog="tool", description="Demo.")
    add_help_cap(parser)
    parser.add_argument("--value", help="Set a value.")

    assert parser.format_help() == (
        "Usage:\n"
        "tool [-h] [--value VALUE]\n\n"
        "Demo.\n\n"
        "Options:\n"
        "  -h, --help     Show this help message and exit.\n"
        "                 \n"
        "  --value VALUE  Set a value.\n"
    )


def test_pilot_private_help_types_are_not_exported() -> None:
    assert "_HelpExample" not in utils_cli.__all__
    assert "_SectionedHelpConfig" not in utils_cli.__all__


def test_private_sectioned_help_preserves_semantic_lines() -> None:
    command_paragraph = (
        "Command prose remains on one intentionally supplied rendered line "
        "even when that physical output line extends beyond column 79."
    )
    second_command_paragraph = (
        "A second explicitly supplied paragraph remains separate."
    )
    option_paragraph = (
        "Ordinary option prose also remains on one intentionally supplied "
        "rendered line beyond the source-width target."
    )
    bullet = (
        "- nested: This deliberately long bullet remains one physical bullet "
        "line and is never automatically wrapped by the private renderer."
    )
    final_option_paragraph = (
        "A final option paragraph remains separated by one explicit blank."
    )
    parser = CapArgumentParser(
        prog="tool",
        description=f"{command_paragraph}\n\n{second_command_paragraph}",
        _sectioned_help=utils_cli._SectionedHelpConfig(
            usage_rows=(("value",),),
            examples=(
                utils_cli._HelpExample(
                    description="Run the private test command.",
                    command_lines=("tool --value example",),
                ),
            ),
        ),
    )
    parser.add_argument(
        "--value",
        help=(f"{option_paragraph}\n{bullet}\n\n{final_option_paragraph}"),
    )

    rendered = parser.format_help()

    assert len(f"  {command_paragraph}") > 79
    assert len(f"    {option_paragraph}") > 79
    assert len(f"      {bullet}") > 79
    assert rendered == (
        "Usage\n"
        "-----\n"
        "  tool\n"
        "    [--value VALUE]\n"
        "\n"
        f"  {command_paragraph}\n"
        "\n"
        f"  {second_command_paragraph}\n"
        "\n"
        "Options\n"
        "-------\n"
        "  --value VALUE\n"
        f"    {option_paragraph}\n"
        f"      {bullet}\n"
        "\n"
        f"    {final_option_paragraph}\n"
        "\n"
        "Examples\n"
        "--------\n"
        "  1. Run the private test command.\n"
        "    '''bash\n"
        "    tool --value example\n"
        "    '''\n"
    )
