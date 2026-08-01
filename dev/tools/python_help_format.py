#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: python_help_format.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.

"""
Format the proven adjacent-constant Python 'help=' literal subset.
"""

from __future__ import annotations

import argparse
import ast
import json
import re
from pathlib import Path

HELP_GROUP = re.compile(
    r"(?P<prefix>^(?P<base> +)help=\(\n)"
    r"(?P<body>(?:(?P=base)    \"(?:[^\"\\]|\\.)*\"\n){2,})"
    r"(?P<suffix>(?P=base)\),)",
    re.MULTILINE,
)
STRING_TOKEN = re.compile(r'"(?:[^"\\]|\\.)*"')


def _literal(value: str) -> str:
    """
    Serialize one Unicode-preserving one-line source literal.
    """

    return json.dumps(value, ensure_ascii=False)


def _chunks(value: str, indent: str) -> list[str]:
    """
    Split one rendered value greedily without changing any character.
    """

    chunks: list[str] = []

    for segment in value.splitlines(keepends=True) or [""]:
        units = re.findall(r"\S+\s*|\s+", segment)
        current = ""

        for unit in units:
            candidate = current + unit
            if current and len(indent + _literal(candidate)) > 79:
                chunks.append(current)
                current = unit
            else:
                current = candidate

        if current or not chunks:
            chunks.append(current)

    return chunks


def format_source(text: str) -> str:
    """
    Greedily format eligible adjacent literals and preserve AST value.
    """

    ast.parse(text)

    def replace(match: re.Match[str]) -> str:
        tokens = STRING_TOKEN.findall(match.group("body"))
        value = ast.literal_eval("(" + "\n".join(tokens) + ")")
        indent = match.group("base") + "    "
        body = "".join(
            f"{indent}{_literal(chunk)}\n" for chunk in _chunks(value, indent)
        )
        replacement = match.group("prefix") + body + match.group("suffix")

        before = ast.parse("x = (" + "\n".join(tokens) + ")")
        after_tokens = STRING_TOKEN.findall(body)
        after = ast.parse("x = (" + "\n".join(after_tokens) + ")")
        if ast.literal_eval(before.body[0].value) != ast.literal_eval(
            after.body[0].value,
        ):
            raise ValueError("help-literal formatting changed rendered value")

        return replacement

    formatted = HELP_GROUP.sub(replace, text)
    ast.parse(formatted)

    return formatted


def main(argv: list[str] | None = None) -> int:
    """
    Preview by default and write only with an explicit flag.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    changed = False

    for path in args.paths:
        absolute = path if path.is_absolute() else root / path
        original = absolute.read_text(encoding="utf-8")
        formatted = format_source(original)
        if formatted == original:
            continue

        changed = True
        if args.write:
            absolute.write_text(formatted, encoding="utf-8")
        else:
            print(absolute.relative_to(root))

    return 0 if args.write or not changed else 1


if __name__ == "__main__":
    raise SystemExit(main())
