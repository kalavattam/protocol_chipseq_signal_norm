#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_cli.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


import pytest

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
