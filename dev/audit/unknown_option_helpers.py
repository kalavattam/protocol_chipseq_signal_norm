#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: unknown_option_helpers.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""Check bounded, recognized unknown-option branches by interface owner."""

from __future__ import annotations

import argparse
import dataclasses
import json
import re
import sys
from pathlib import Path

from help_heredoc_reflow import FUNCTION_START, extract_help_heredocs

UNKNOWN_TEXT = "unknown option/parameter passed: '${1}'."
UNKNOWN_MARKER = re.compile(r"unknown (?:option/)?parameter passed: '\$\{1\}'\.?")
HELP_CALL = re.compile(r"(?:echo\s+\"\$\{show_help\}\"|\b(?:help|show_help|detail)_[A-Za-z0-9_]+\b)")
STATUS = re.compile(r"\b(?:return|exit)\s+1\b")
HELPER = re.compile(
    r"^\s*(?P<helper>echo_err_func|echo_err|error)(?:\s|$)", re.MULTILINE
)


@dataclasses.dataclass(frozen=True)
class Finding:
    path: str
    line: int
    owner: str
    rule_id: str
    message: str

    def format(self) -> str:
        return f"{self.rule_id}: {self.path}:{self.line}: owner={self.owner}; {self.message}"


def function_ranges(text: str) -> list[tuple[str, int, int]]:
    """Return bounded top-level function ranges using repository brace style."""

    lines = text.splitlines()
    ranges: list[tuple[str, int, int]] = []
    index = 0
    while index < len(lines):
        match = FUNCTION_START.match(lines[index])
        if match is None:
            index += 1
            continue
        depth = lines[index].count("{") - lines[index].count("}")
        cursor = index + 1
        while cursor < len(lines) and depth > 0:
            depth += lines[cursor].count("{") - lines[cursor].count("}")
            cursor += 1
        ranges.append((match.group("name"), index + 1, cursor))
        index = max(cursor, index + 1)
    return ranges


def lexical_owner(text: str, line: int) -> str:
    candidates = [row for row in function_ranges(text) if row[1] <= line <= row[2]]
    return candidates[-1][0] if candidates else "<file>"


def is_top_level_cli(path: str) -> bool:
    value = Path(path)
    return (
        (value.parent == Path("bin") and value.suffix == ".sh")
        or (value.parent == Path("install/scripts") and value.suffix == ".sh")
        or path == "tests/support/clean_artifacts.sh"
    )


def interface_class(path: str, owner: str, heredoc_owners: set[str]) -> str:
    """Classify a recognized parser by documented interface ownership."""

    if is_top_level_cli(path) and owner in {"<file>", "main", "parse_args", "check_args_light"}:
        return "script"
    if path.startswith("lib/bash/") or owner in heredoc_owners:
        return "reusable_function"
    return "script"


def helper_available(text: str, helper: str) -> bool:
    """Recognize the established format_outputs sourcing layer or local definition."""

    if re.search(rf"function\s+{re.escape(helper)}\s*\(\)", text):
        return True
    return "format_outputs" in text or "source_helpers" in text


def check_text(path: str, text: str) -> list[Finding]:
    """Check recognized unknown-option messages without parsing general shell."""

    lines = text.splitlines()
    heredoc_owners = {item.owner for item in extract_help_heredocs(text)}
    findings: list[Finding] = []
    for index, line in enumerate(lines, 1):
        if UNKNOWN_MARKER.search(line) is None:
            continue
        owner = lexical_owner(text, index)
        classification = interface_class(path, owner, heredoc_owners)
        start = max(0, index - 2)
        end = min(len(lines), index + 7)
        block = "\n".join(lines[start:end])
        helper_match = HELPER.search(block)
        current = helper_match.group("helper") if helper_match else "<missing>"
        required = "echo_err_func" if classification == "reusable_function" else "echo_err"
        identity = f"{path}::{owner}"
        if UNKNOWN_TEXT not in block:
            findings.append(Finding(path, index, identity, "SHELL.UNKNOWN_OPTION.MESSAGE", f"required exact message: {UNKNOWN_TEXT}"))
        if current != required and (
            helper_available(text, required)
            or current in {"echo_err", "echo_err_func"}
        ):
            findings.append(Finding(path, index, identity, "SHELL.UNKNOWN_OPTION.HELPER", f"{classification} interface requires {required}; found {current}"))
        if classification == "reusable_function" and "${FUNCNAME[0]}" not in block:
            findings.append(Finding(path, index, identity, "SHELL.UNKNOWN_OPTION.FUNCNAME", "function error requires ${FUNCNAME[0]}"))
        if HELP_CALL.search(block) is None:
            findings.append(Finding(path, index, identity, "SHELL.UNKNOWN_OPTION.HELP", "recognized branch must preserve help emission"))
        if STATUS.search(block) is None:
            findings.append(Finding(path, index, identity, "SHELL.UNKNOWN_OPTION.STATUS", "recognized branch must preserve return or exit status 1"))
    return findings


def shell_sources(root: Path) -> list[Path]:
    paths = [
        *(root / "bin").glob("*.sh"),
        *(root / "lib/bash").rglob("*.sh"),
        *(root / "install/scripts").glob("*.sh"),
    ]
    cleaner = root / "tests/support/clean_artifacts.sh"
    if cleaner.is_file():
        paths.append(cleaner)
    return sorted(paths)


def scan_repository(root: Path) -> tuple[list[Finding], list[dict[str, object]]]:
    findings: list[Finding] = []
    inventory: list[dict[str, object]] = []
    for target in shell_sources(root.resolve()):
        path = str(target.relative_to(root.resolve()))
        text = target.read_text(encoding="utf-8")
        path_findings = check_text(path, text)
        findings.extend(path_findings)
        for index, line in enumerate(text.splitlines(), 1):
            if UNKNOWN_MARKER.search(line) is None:
                continue
            owner = lexical_owner(text, index)
            classification = interface_class(
                path, owner, {item.owner for item in extract_help_heredocs(text)}
            )
            block = "\n".join(text.splitlines()[max(0, index - 2):index + 7])
            match = HELPER.search(block)
            required = "echo_err_func" if classification == "reusable_function" else "echo_err"
            inventory.append({
                "path": path,
                "line": index,
                "lexical_function": None if owner == "<file>" else owner,
                "documented_interface_owner": f"{path}::{owner}",
                "classification": classification,
                "current_helper": match.group("helper") if match else "<missing>",
                "required_helper": required,
                "message_text": UNKNOWN_TEXT if UNKNOWN_TEXT in block else "noncanonical",
                "help_emission": bool(HELP_CALL.search(block)),
                "status_behavior": "return_or_exit_1" if STATUS.search(block) else "missing",
                "helper_available": helper_available(text, required),
            })
    return findings, inventory


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--inventory-output", type=Path)
    args = parser.parse_args(argv)
    findings, inventory = scan_repository(args.root)
    for finding in findings:
        print(finding.format())
    if args.inventory_output:
        args.inventory_output.write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    if findings:
        print(f"SHELL.UNKNOWN_OPTION: {len(findings)} violation(s)")
        return 1
    print("SHELL.UNKNOWN_OPTION: pass")
    return 0


if __name__ == "__main__":
    sys.exit(main())
