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

    def add_usage(
        self,
        usage: str | None,
        actions: list[argparse.Action],
        groups: list[argparse._MutuallyExclusiveGroup],
        prefix: str | None = None,
    ) -> None:
        """
        Add parser usage with the repository's preferred heading.

        The formatter writes the usage block to its stream.

        Parameters
        ----------
        usage : str | None
            Usage text emitted by the customized formatter.
        actions : list[argparse.Action]
            Parser actions represented in the generated usage block.
        groups : list[argparse._MutuallyExclusiveGroup]
            Mutually exclusive parser groups represented in usage.
        prefix : str | None
            Usage-heading prefix, or None for the canonical heading.
        """

        if prefix is None:
            prefix = "Usage:\n"

        super().add_usage(usage, actions, groups, prefix)

    def start_section(self, heading: str | None) -> None:
        """
        Start a help section after capitalizing its heading.

        The formatter writes the section heading to its stream.

        Parameters
        ----------
        heading : str | None
            Help-section heading to emit.
        """

        if heading:
            heading = heading[:1].upper() + heading[1:]

        super().start_section(heading)


class CapArgumentParser(argparse.ArgumentParser):
    """
    Provide capitalized raw-text help without implicit defaults or help flags.
    """

    def __init__(self, *args: object, **kwargs: object) -> None:
        """
        Initialize a parser with the repository's stable help defaults.

        Parameters
        ----------
        *args : object
            Positional arguments forwarded to 'argparse.ArgumentParser'.
        **kwargs : object
            Parser options, with repository defaults supplied when absent.
        """

        kwargs.setdefault("formatter_class", CapHelpFormatter)
        kwargs.setdefault("argument_default", argparse.SUPPRESS)
        kwargs.setdefault("add_help", False)
        kwargs.setdefault("allow_abbrev", False)
        super().__init__(*args, **kwargs)


def add_help_cap(parser: argparse.ArgumentParser) -> None:
    """
    Add a standard '-h' / '--help' with preferred wording and spacing.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Argument parser whose help behavior is extended.
    """

    parser.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this help message and exit.\n\n",
    )
