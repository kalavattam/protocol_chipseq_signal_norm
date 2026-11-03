#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: check_submit_compute_signal.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Check how 'submit_compute_signal.sh' calls its dispatch helpers.

- 'arity <script>': for each 'run_comp_*' call site, compare the arguments it
  passes with the count its helper's '[[ $# -ne <n> ]]' guard expects. Prints
  one 'OK' or 'MISMATCH' line per call site, then 'CALLS <n>'.
- 'forwarding <script>': classify each 'cmd+=( --report_n_frg ... )' and
  '--report_n_ovlp' append. A bare append must sit under the 'DERIVE' guard,
  and a valued one under a non-empty guard. Prints one line per append, then
  'BARE <n>' and 'VALUED <n>'.

The calling test reads the printed report; both subcommands exit 0.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def parse_args(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    """
    Parse the subcommand and the script path.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Check how 'submit_compute_signal.sh' calls its dispatch helpers."
        ),
    )
    parser.add_argument("cmd", choices=("arity", "forwarding"))
    parser.add_argument("script", help="Path of 'submit_compute_signal.sh'.")

    return parser.parse_args(argv)


def report_arity(path: str) -> None:
    """
    Print each call site's argument count against its helper's guard.
    """

    source = Path(path).read_text(encoding="utf-8").splitlines()

    guards = {}
    current = None
    for line in source:
        named = re.match(r"function (run_comp_\w+)\(\)", line)
        if named:
            current = named.group(1)
        guard = re.search(r"\[\[ \$# -ne (\d+) \]\]", line)
        if guard and current:
            guards.setdefault(current, int(guard.group(1)))

    calls = []
    for index, line in enumerate(source):
        called = re.match(r"\s*(run_comp_\w+) \\\s*$", line)
        if not called:
            continue
        count, cursor = 0, index + 1
        while cursor < len(source):
            count += 1
            if not source[cursor].rstrip().endswith("\\"):
                break
            cursor += 1
        calls.append((called.group(1), index + 1, count))

    for name, line_no, count in calls:
        want = guards.get(name)
        status = "OK" if want == count else "MISMATCH"
        print(
            f"{status} {name} line {line_no}: passes {count}, "
            f"guard expects {want}"
        )

    print(f"CALLS {len(calls)}")


def report_forwarding(path: str) -> None:
    """
    Print how each report append is guarded.
    """

    source = Path(path).read_text(encoding="utf-8").splitlines()

    bare = 0
    valued = 0

    for index, line in enumerate(source):
        append = re.match(
            r"\s*cmd\+=\( --report_(n_frg|n_ovlp)(.*)\)\s*$",
            line,
        )

        if not append:
            continue

        label, tail = append.group(1), append.group(2).strip()
        var = f"report_{label}"
        line_no = index + 1
        window = [entry.strip() for entry in source[max(0, index - 2) : index]]

        derive = f'if [[ "${{{var}}}" == "DERIVE" ]]; then'
        valued_guard = f'elif [[ -n "${{{var}}}" ]]; then'

        if tail == "":
            bare += 1

            if derive in window:
                print(f"OK-DERIVE --report_{label} line {line_no}")
            else:
                print(f"BARE-UNGUARDED --report_{label} line {line_no}")
        elif tail == f'"${{{var}}}"':
            valued += 1

            if valued_guard in window:
                print(f"OK-VALUED --report_{label} line {line_no}")
            else:
                print(f"UNGUARDED --report_{label} line {line_no}")
        else:
            print(f"UNEXPECTED --report_{label} line {line_no}: {tail!r}")

    print(f"BARE {bare}")
    print(f"VALUED {valued}")


def main(argv: list[str] | None = None) -> int:
    """
    Print the subcommand's report.
    """

    args = parse_args(argv)

    if args.cmd == "arity":
        report_arity(args.script)
    else:
        report_forwarding(args.script)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
