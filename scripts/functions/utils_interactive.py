#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_interactive.py
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
utils_interactive.py


Description
-----------
Shared helper functions to support an “interactive mode” across Python scripts
in the codebase. These functions enable “CLI escape hatches” in the presence of
hardcoded 'interactive = True' or 'interactive = False' in the scripts.


Flags
-----
--interactive-help, --interactive_help
    Print a brief help message describing the interactive pre-parser flags,
    then exit(0) before the main parser runs.

--interactive
    Force interactive mode on (even if the script has a hardcoded
    'interactive = False').

--no-interactive, --no_interactive
    Force interactive mode off (even if the script has a hardcoded
    'interactive = True').

--interactive-echo, --interactive_echo
    When running in interactive mode, echo the derived paths and args (which
    is useful for debugging the 'set_interactive()').

--interactive-echo-only, --interactive_echo_only
    Echo the derived paths and args in interactive mode and exit without
    running the main workflow.

--interactive-echo-check, --interactive_echo_check
    Echo the derived paths and args in interactive mode, run any "ensure"/
    validation logic in 'set_interactive()', then exit without running the
    main workflow. Intended as a quick “checked dry run” of interactive mode.

--interactive-ensure, --interactive_ensure
    When running in interactive mode, perform stricter checks (e.g., assert
    files or directories exist, create output directories) inside
    'set_interactive()'.


Example
-------
Pre-parse these flags before the main parser to decide whether to synthesize
'args' via 'set_interactive()' or to require full CLI args:
'''python
from scripts.functions.utils_interactive import (
    echo_block,
    preparse_interactive,
    resolve_interactive
)

...

def main():
    argv = sys.argv[1:]
    pre_ns, remaining = preparse_interactive(argv)
    use_interactive = resolve_interactive(
        default=interactive,  # The module-level default
        pre_ns=pre_ns
    )

    if use_interactive:
        args = set_interactive(
            echo=getattr(pre_ns, "interactive_echo", False),
            ensure=getattr(pre_ns, "interactive_ensure", False),
        )
    else:
        args = parse_args(remaining)

    ...
'''


Notes
-----
- This module is meant to be pre-parsed prior to running a given script’s main
  parser. That way, the user can decide to assemble an 'argparse.Namespace'
  from 'set_interactive()' or to require full CLI arguments.
- These flags are intentionally defined in a minimal "pre-parser" with
  'add_help=False'. This keeps them invisible in the main '--help' output and
  avoids required-argument errors before deciding how to run.
"""

from __future__ import annotations

import argparse
import sys

from typing import Callable, Tuple

#  Require Python 3.10 across the codebase (PEP 585 list[str] is used)
assert sys.version_info >= (3, 10), "Python >= 3.10 required."

__all__ = [
    "echo_block",
    "help_interactive",
    "preparse_interactive",
    "resolve_interactive",
    "get_args_interactive"
]


def echo_block(title: str, obj, sort_keys: bool = False) -> None:
    """
    Pretty-print a titled block of 'key = value' pairs.

    Args:
        title : str
            Block title, printed as '## {title} ##'.
        obj :
            Either a mapping (e.g., 'dict') or an object with a '__dict__'
            (e.g., 'argparse.Namespace'). For the latter, 'vars(obj)' is used.
        sort_keys : bool
            - If True, keys are printed in alphabetical order.
            - If False (default), preserves the original insertion order of the
              mapping.

    Returns:
        None. Prints to stdout.

    Raises:
        None.
    """
    #  Normalize 'obj' to a mapping
    if hasattr(obj, "__dict__"):
        data = vars(obj)
    else:
        data = obj

    print(f"\n## {title} ##")

    if not data:
        print("(no entries)")
        return

    #  Use insertion order or sorted order depending on 'sort_keys'
    keys = sorted(data) if sort_keys else list(data.keys())

    w = max(len(k) for k in keys)
    for k in keys:
        print(f"{k:<{w}} = {data[k]}")


def help_interactive() -> str:
    """
    Return a compact help string for the interactive pre-parser flags.
    This is intentionally separate from the main parser’s '--help'.
    """
    return (
        "Interactive pre-parser flags (handled before the main parser):\n"
        "  --interactive-help       Show this help message and exit.\n"
        "  --interactive            Force interactive mode on (override 'interactive = False').\n"
        "  --no-interactive         Force interactive mode off (override 'interactive = True').\n"
        "  --interactive-echo       Force interactive mode on and “echo” interactive arguments before running.\n"
        "  --interactive-echo-only  Force interactive mode on and “echo” interactive arguments and exit without running main workflow.\n"
        "  --interactive-echo-check Force interactive mode on, “echo” interactive arguments, run any ensure/validation logic in\n"
        "                           'set_interactive()', and exit without running main workflow.\n"
        "  --interactive-ensure     Force interactive mode on and create any necessary output directories before running.\n"
        "\n"
        "Usage pattern:\n"
        "  pre_ns, remaining = preparse_interactive(sys.argv[1:])\n"
        "  use_interactive = resolve_interactive(default, pre_ns)\n"
        "  args = set_interactive(...) if use_interactive else parse_args(remaining)\n"
    )


def preparse_interactive(
    argv: list[str]
) -> tuple[argparse.Namespace, list[str]]:
    """
    Before running the main parser (which may enforce required args), pre-parse
    interactive-only flags, deciding between 'interactive' and 'CLI' (“command
    line interface”) modes.

    Args:
        argv : list[str]
            Typically 'sys.argv[1:]'

    Returns:
        (pre_ns, remaining_argv) : tuple[argparse.Namespace, list[str]]
            'pre_ns' contains booleans for the above flags; 'remaining_argv' is
            passed to the main parser

    Raises:
        None.
    """
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument(
        "--interactive-help", "--interactive_help",
        "--help-interactive", "--help_interactive",
        action="store_true"
    )
    parser.add_argument(
        "--interactive",
        action="store_true",
        dest="interactive",
        help="Run script in interactive mode."
    )
    parser.add_argument(
        "--no-interactive", "--no_interactive",
        "--interactive-no", "--interactive_no",
        action="store_true",
        dest="no_interactive"
    )
    parser.add_argument(
        "--interactive-echo", "--interactive_echo",
        "--echo-interactive", "--echo_interactive",
        action="store_true",
        dest="interactive_echo",
        help="Echo interactive arguments before running."
    )
    parser.add_argument(
        "--interactive-echo-only", "--interactive_echo_only",
        "--echo-only-interactive", "--echo_only_interactive",
        action="store_true",
        dest="interactive_echo_only",
        help=(
            "Echo interactive arguments and exit without executing main "
            "workflow."
        )
    )
    parser.add_argument(
        "--interactive-echo-check", "--interactive_echo_check",
        "--echo-check-interactive", "--echo_check_interactive",
        action="store_true",
        dest="interactive_echo_check",
        help=(
            "Echo interactive arguments, apply ensure/validation logic in "
            "interactive mode, and exit without executing main workflow."
        )
    )
    parser.add_argument(
        "--interactive-ensure", "--interactive_ensure",
        "--ensure-interactive", "--ensure_interactive",
        action="store_true",
        dest="interactive_ensure",
        help="Create any necessary output directories in interactive mode."
    )
    pre_ns, remaining = parser.parse_known_args(argv)

    #  If help is requested, print a concise message and exit before main parse
    if getattr(pre_ns, "interactive_help", False):
        print(help_interactive())
        sys.exit(0)

    return pre_ns, remaining


def resolve_interactive(
    default: bool,
    pre_ns: argparse.Namespace
) -> bool:
    """
    Decide whether to run interactive mode.

    Args:
        default : bool
            The script's module-level default (e.g., hardcoded
            'interactive = True' or 'interactive = False').
        pre_ns : argparse.Namespace
            Result from 'preparse_interactive()'.

    Returns:
        bool
            True if interactive mode should be used; otherwise False.

    Raises:
        ValueError
            If both '--interactive' and '--no_interactive' are provided.

    Precedence:
        1. Return True if any of '--interactive', '--interactive-echo',
           '--interactive-echo-only', '--interactive-echo-check', or
           '--interactive-ensure' is supplied.
        2. Else return False if '--no_interactive'.
        3. Otherwise, return 'default' (i.e., the script's hardcoded default).
    """
    force_on = any(
        getattr(pre_ns, name, False)
        for name in (
            "interactive",
            "interactive_echo",
            "interactive_echo_only",
            "interactive_echo_check",
            "interactive_ensure"
        )
    )
    force_off = getattr(pre_ns, "no_interactive", False)

    if force_on and force_off:
        raise ValueError(
            "Error: Cannot pass both '--interactive*' and '--no-interactive'."
        )

    if force_on:
        return True

    if force_off:
        return False

    return default


def get_args_interactive(
    argv: list[str] | None,
    interactive_default: bool,
    set_interactive: Callable[..., argparse.Namespace],
    parse_args: Callable[[list[str]], argparse.Namespace],
) -> Tuple[argparse.Namespace, bool]:
    """
    Shared helper to resolve interactive vs CLI args.

    Args:
        argv :
            Typically 'sys.argv[1:]' or 'None' (in which case 'sys.argv[1:]'
            is used inside this function).
        interactive_default : bool
            Module-level 'interactive' default for the script.
        set_interactive :
            Script-specific 'set_interactive(echo=..., ensure=...)' function.
        parse_args :
            Script-specific 'parse_args(remaining_argv)' function.

    Returns:
        (args, early_exit) :
            - args        argparse.Namespace containing the effective
                          arguments.
            - early_exit  True if interactive echo-only / echo-check requested
                          and the caller should return 0 immediately.

    Raises:
        ValueError
            If conflicting interactive flags are passed (via
            'resolve_interactive').
    """
    #  Normalize 'argv'
    argv_use = sys.argv[1:] if argv is None else argv

    #  Pre-parse interactive flags
    pre_ns, remaining = preparse_interactive(argv_use)
    use_interactive = resolve_interactive(interactive_default, pre_ns)

    if use_interactive:
        #  Echo if any of the 'echo' flags are set
        echo = any(
            getattr(pre_ns, name, False)
            for name in (
                "interactive_echo",
                "interactive_echo_only",
                "interactive_echo_check",
            )
        )

        #  Run 'ensure' logic if either of the following:
        #    - explicit '--interactive-ensure', or
        #    - '--interactive-echo-check' (checked dry run)
        ensure = (
            getattr(pre_ns, "interactive_ensure", False)
            or getattr(pre_ns, "interactive_echo_check", False)
        )

        args = set_interactive(echo, ensure)

        #  Echo-only / echo-check: caller should return 0 without running main
        if getattr(pre_ns, "interactive_echo_only", False) or getattr(
            pre_ns, "interactive_echo_check", False
        ):
            return args, True

        return args, False

    #  CLI mode: parse remaining 'argv'
    args = parse_args(remaining)
    return args, False
