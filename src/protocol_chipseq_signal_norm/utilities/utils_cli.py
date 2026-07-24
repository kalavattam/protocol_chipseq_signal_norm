#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_cli.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6) were
# used in development and documentation.
#
# Distributed under the MIT license.


"""
Command-line parser helpers.

Notes
-----
'CapArgumentParser' and 'add_help_cap()' provide shared parser behavior for
user-facing Python CLIs.
"""

from __future__ import annotations

import argparse

__all__ = ["CapArgumentParser", "CapHelpFormatter", "add_help_cap"]


class CapHelpFormatter(argparse.RawTextHelpFormatter):
    """
    RawText formatter with capitalized section headings and “Usage”.
    """
    def add_usage(self, usage, actions, groups, prefix=None):
        """Add parser usage with the repository's preferred heading."""

        #  Have "Usage:" on its own line
        if prefix is None:
            prefix = "Usage:\n"
        return super().add_usage(usage, actions, groups, prefix)

    def start_section(self, heading):
        """Start a help section after capitalizing its heading."""

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
