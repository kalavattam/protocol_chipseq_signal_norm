#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: local_http.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Find a free local port, or report whether a local HTTP server answers.

'port' prints a TCP port on 127.0.0.1 that was free when asked. 'ready <url>'
exits 0 when the URL answers with status 200 within half a second, and 1
otherwise. The helpers in 'test_helpers.sh' that serve fixtures over HTTP run
it under the Python they were given, which is at least 3.11.
"""

from __future__ import annotations

import argparse
import socket
import urllib.request


def parse_args(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    """
    Parse the subcommand and, for 'ready', the URL.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Find a free local port, or report whether a local HTTP server "
            "answers."
        ),
    )
    sub = parser.add_subparsers(dest="cmd")
    sub.required = True
    sub.add_parser("port", help="Print a free TCP port on 127.0.0.1.")
    ready = sub.add_parser(
        "ready", help="Exit 0 if the URL answers with status 200."
    )
    ready.add_argument("url", help="URL to request.")

    return parser.parse_args(argv)


def print_free_port() -> int:
    """
    Print a TCP port on 127.0.0.1 that is free now.
    """

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(("127.0.0.1", 0))
        print(sock.getsockname()[1])

    return 0


def check_ready(url: str) -> int:
    """
    Return 0 if the URL answers with status 200 within half a second.
    """

    try:
        with urllib.request.urlopen(url, timeout=0.5) as response:
            return 0 if response.status == 200 else 1
    except Exception:
        return 1


def main(argv: list[str] | None = None) -> int:
    """
    Run the subcommand.
    """

    args = parse_args(argv)

    if args.cmd == "port":
        return print_free_port()

    return check_ready(args.url)


if __name__ == "__main__":
    raise SystemExit(main())
