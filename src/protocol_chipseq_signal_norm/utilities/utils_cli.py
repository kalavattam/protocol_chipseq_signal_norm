#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_cli.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6) were
# used in design, development, and documentation, with all output reviewed,
# edited, and approved by the author.
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
from dataclasses import dataclass

__all__ = ["CapArgumentParser", "CapHelpFormatter", "add_help_cap"]


@dataclass(frozen=True)
class _HelpExample:
    """
    Represent one pilot-private rendered help example.

    The command lines preserve the example owner's one-line or multiline
    intent. Rendering supplies delimiters, indentation, and continuations
    without changing command text.
    """

    description: str
    command_lines: tuple[str, ...]


@dataclass(frozen=True)
class _SectionedHelpConfig:
    """
    Configure the pilot-private sectioned help renderer.

    Each usage row names parser destinations explicitly. The renderer does not
    infer semantic categories from option names.
    """

    usage_rows: tuple[tuple[str, ...], ...]
    examples: tuple[_HelpExample, ...]


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

        sectioned_help = kwargs.pop("_sectioned_help", None)
        if sectioned_help is not None and not isinstance(
            sectioned_help,
            _SectionedHelpConfig,
        ):
            raise TypeError(
                "'_sectioned_help' must be a private sectioned-help config",
            )

        self._sectioned_help = sectioned_help
        kwargs.setdefault("formatter_class", CapHelpFormatter)
        kwargs.setdefault("argument_default", argparse.SUPPRESS)
        kwargs.setdefault("add_help", False)
        kwargs.setdefault("allow_abbrev", False)
        super().__init__(*args, **kwargs)

    def _action_by_destination(self, destination: str) -> argparse.Action:
        """
        Return the one visible action registered for a destination.
        """

        visible = [
            action
            for action in self._actions
            if action.dest == destination
            and action.help is not argparse.SUPPRESS
        ]
        if len(visible) != 1:
            raise ValueError(
                "sectioned help requires exactly one visible action for "
                f"destination '{destination}'",
            )

        return visible[0]

    def _format_usage_invocation(self, action: argparse.Action) -> str:
        """
        Return one bracketed invocation for a configured usage row.
        """

        formatter = self._get_formatter()
        long_options = [
            option
            for option in action.option_strings
            if option.startswith("--")
        ]
        option = long_options[-1] if long_options else action.option_strings[0]
        arguments = formatter._format_args(action, action.dest.upper())
        invocation = option if not arguments else f"{option} {arguments}"
        if action.required:
            return invocation

        return f"[{invocation}]"

    def _format_option_invocation(self, action: argparse.Action) -> str:
        """
        Return stable public aliases followed by one argument expression.
        """

        formatter = self._get_formatter()
        options = ", ".join(action.option_strings)
        arguments = formatter._format_args(action, action.dest.upper())

        return options if not arguments else f"{options} {arguments}"

    @staticmethod
    def _render_description(text: str) -> list[str]:
        """
        Render option prose and stable nested bullets below an invocation.
        """

        output: list[str] = []

        for source_line in text.strip().splitlines():
            stripped = source_line.strip()

            if not stripped:
                output.append("")

                continue

            if stripped.startswith("- "):
                output.append(f"      {stripped}")
            else:
                output.append(f"    {stripped}")

        return output

    def _render_sectioned_help(self) -> str:
        """
        Render the explicitly configured Shell-like pilot help document.
        """

        config = self._sectioned_help
        if config is None:
            raise RuntimeError("sectioned help was not configured")

        formatter = self._get_formatter()
        actions = {
            destination: self._action_by_destination(destination)
            for row in config.usage_rows
            for destination in row
        }
        lines = ["Usage", "-----", f"  {self.prog}"]

        for row in config.usage_rows:
            invocations = [
                self._format_usage_invocation(actions[destination])
                for destination in row
            ]
            lines.append("    " + " ".join(invocations))

        if self.description:
            lines.append("")

            for source_line in self.description.strip().splitlines():
                stripped = source_line.strip()
                lines.append("" if not stripped else f"  {stripped}")

        lines.extend(["", "Options", "-------"])

        for index, action in enumerate(
            actions[destination]
            for row in config.usage_rows
            for destination in row
        ):
            if index:
                lines.append("")

            invocation = self._format_option_invocation(action)
            lines.append(f"  {invocation}")
            help_text = formatter._expand_help(action)
            lines.extend(self._render_description(help_text))

        lines.extend(["", "Examples", "--------"])

        for index, example in enumerate(config.examples, 1):
            if index > 1:
                lines.append("")

            lines.append(f"  {index}. {example.description}")
            lines.append("    '''bash")

            if len(example.command_lines) == 1:
                lines.append(f"    {example.command_lines[0]}")
            else:
                for position, command_line in enumerate(example.command_lines):
                    continuation = (
                        " \\"
                        if position + 1 < len(example.command_lines)
                        else ""
                    )
                    indentation = "    " if position == 0 else "        "
                    lines.append(
                        f"{indentation}{command_line}{continuation}",
                    )

            lines.append("    '''")

        return "\n".join(lines) + "\n"

    def format_help(self) -> str:
        """
        Return default help or the explicitly configured private pilot form.
        """

        if self._sectioned_help is None:
            return super().format_help()

        return self._render_sectioned_help()


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
