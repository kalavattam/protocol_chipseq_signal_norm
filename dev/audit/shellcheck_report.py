#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: shellcheck_report.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Merge and summarize the JSON output of 'run_shellcheck.sh'.

- 'merge <output> <part>...': join the findings of ShellCheck's 'json1' batches
  into one '{"comments": [...]}' file.
- 'summary <summary> <executable> <version> <bash_raw> <bash_files>
  <bash_status> <sh_raw> <sh_files> <sh_status> <status>': write the run's
  summary JSON, and print '<finding count><tab><status>'. The status becomes 2,
  an infrastructure failure, when either raw file is missing or is not
  ShellCheck JSON.

Arguments are read by position, as 'run_shellcheck.sh' passes them.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

USAGE = (
    "usage: shellcheck_report.py merge <output> <part>...\n"
    "       shellcheck_report.py summary <summary> <executable> <version> "
    "<bash_raw> <bash_files> <bash_status> <sh_raw> <sh_files> <sh_status> "
    "<status>"
)


def merge_parts(output: str, parts: list[str]) -> None:
    """
    Join the findings of each batch into one file.
    """

    comments: list[object] = []

    for supplied in parts:
        value = json.loads(Path(supplied).read_text(encoding="utf-8"))

        if isinstance(value, dict) and isinstance(value.get("comments"), list):
            comments.extend(value["comments"])
        elif isinstance(value, list):
            comments.extend(value)
        else:
            raise ValueError(f"unexpected ShellCheck JSON in {supplied}")

    Path(output).write_text(
        json.dumps({"comments": comments}, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def load_findings(path: str) -> tuple[list[object], bool]:
    """
    Return a raw file's findings and whether it held ShellCheck JSON.
    """

    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return [], False

    if isinstance(value, dict) and isinstance(value.get("comments"), list):
        return value["comments"], True

    if isinstance(value, list):
        return value, True

    return [], False


def write_summary(
    summary_path: str,
    executable: str,
    version: str,
    bash_path: str,
    bash_files: str,
    bash_status: str,
    sh_path: str,
    sh_files: str,
    sh_status: str,
    requested_status: str,
) -> None:
    """
    Write the summary JSON and print the finding count and status.
    """

    bash_findings, bash_valid = load_findings(bash_path)
    sh_findings, sh_valid = load_findings(sh_path)
    status = int(requested_status)

    if not bash_valid or not sh_valid:
        status = 2

    disposition = (
        "clean"
        if status == 0
        else "completed_with_findings"
        if status == 1
        else "infrastructure_failure"
    )

    summary = {
        "schema_version": 1,
        "executable": executable,
        "version": version,
        "source_path": str(Path.cwd()),
        "flags": ["-x", "-P", str(Path.cwd())],
        "file_count": int(bash_files) + int(sh_files),
        "finding_count": len(bash_findings) + len(sh_findings),
        "status": status,
        "disposition": disposition,
        "languages": {
            "bash": {
                "shell": "bash",
                "file_count": int(bash_files),
                "finding_count": len(bash_findings),
                "status": int(bash_status),
                "raw_findings": Path(bash_path).name,
                "raw_valid": bash_valid,
            },
            "sh": {
                "shell": "sh",
                "file_count": int(sh_files),
                "finding_count": len(sh_findings),
                "status": int(sh_status),
                "raw_findings": Path(sh_path).name,
                "raw_valid": sh_valid,
            },
        },
    }

    Path(summary_path).write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"{summary['finding_count']}\t{status}")


def main(argv: list[str] | None = None) -> int:
    """
    Run the subcommand.
    """

    args = sys.argv[1:] if argv is None else argv

    if len(args) >= 2 and args[0] == "merge":
        merge_parts(args[1], args[2:])
    elif len(args) == 11 and args[0] == "summary":
        write_summary(*args[1:])
    else:
        raise SystemExit(USAGE)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
