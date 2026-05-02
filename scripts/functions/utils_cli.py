#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_cli.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5-series models) was used in development.
#
# Distributed under the MIT license.

"""
Script
------
utils_cli.py


Description
-----------
Small CLI-formatting helpers built around 'argparse'.

This module provides:
    - a help formatter that capitalizes section headings and prints 'Usage:'
      on its own line,
    - an 'ArgumentParser' subclass with preferred defaults for this codebase,
      and
    - a helper to add a standard '-h' / '--help' flag.


Variables
---------
__all__
    Public exports for the module: 'CapHelpFormatter', 'CapArgumentParser', and
    'add_help_cap'.


Classes
-------
CapHelpFormatter
    Subclass of 'argparse.RawTextHelpFormatter' that capitalizes section
    headings and prints 'Usage:' on its own line.

CapArgumentParser
    Subclass of 'argparse.ArgumentParser' configured for this codebase:
        - uses 'CapHelpFormatter',
        - suppresses implicit defaults in parsed namespaces,
        - disables automatic help injection, and
        - disables abbreviation of long options.


Functions
---------
add_help_cap()
    Add a standard '-h' / '--help' action with preferred wording and spacing.


Example
-------
'''python
from scripts.functions.utils_cli import CapArgumentParser, add_help_cap

parser = CapArgumentParser(description="Example parser.")
add_help_cap(parser)
parser.add_argument("--foo", type=int, required=True)
args = parser.parse_args()
'''


Notes
-----
- 'CapArgumentParser' sets 'add_help=False', so callers should usually invoke
  'add_help_cap(parser)' explicitly.
- 'CapHelpFormatter' preserves raw newlines in help text, as it subclasses
  'argparse.RawTextHelpFormatter'.
"""

from __future__ import annotations
import argparse

__all__ = ["CapHelpFormatter", "CapArgumentParser", "add_help_cap"]


class CapHelpFormatter(argparse.RawTextHelpFormatter):
    """
    RawText formatter with capitalized section headings and “Usage”.
    """
    def add_usage(self, usage, actions, groups, prefix=None):
        #  Have "Usage:" on its own line
        if prefix is None:
            prefix = "Usage:\n"
        return super().add_usage(usage, actions, groups, prefix)

    def start_section(self, heading):
        if heading:
            heading = heading[:1].upper() + heading[1:]
        return super().start_section(heading)


class CapArgumentParser(argparse.ArgumentParser):
    """
    Opinionated parser: capitalized help, RawText, suppress implicit defaults,
    and no implicit '-h' / '--help'.
    """
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("formatter_class", CapHelpFormatter)
        kwargs.setdefault("argument_default", argparse.SUPPRESS)
        kwargs.setdefault("add_help", False)
        kwargs.setdefault("allow_abbrev", False)
        super().__init__(*args, **kwargs)


def add_help_cap(parser: argparse.ArgumentParser) -> None:
    """
    Add a standard '-h' / '--help' with preferred wording and spacing.
    """
    parser.add_argument(
        "-h", "--help",
        action="help",
        help="Show this help message and exit.\n\n",
    )
