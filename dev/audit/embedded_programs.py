#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: embedded_programs.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Count programs embedded in other source, and inputs that tests write.

One checker for 'SOURCE.EMBED.PROGRAMS' and its realizations, and for the
record-input and inline-stub parts of 'TEST.FIXTURES'. Each diagnostic is
counted against its own limit:
- 'SHELL.EMBED.PYTHON': Python fed to an interpreter through a heredoc, or a
  'python -c' program of more than one line.
- 'SHELL.EMBED.BASH_C': a 'bash -c' program of more than one line.
- 'SHELL.EMBED.AWK': an 'awk' program of more than three lines.
- 'SHELL.EMBED.SED': a 'sed' program of more than one line.
- 'SHELL.EMBED.SED_GNU': a 'sed' option or escape that BSD 'sed' lacks.
- 'SHELL.EMBED.PERL': any 'perl' program.
- 'PY.EMBED.STRING_INPUT': a string literal of four or more lines in a Python
  test, other than a docstring. Whether it is an input or an expected output is
  left to review.
- 'TEST.FIXTURES.RECORD_INPUT': a line in a test that writes a file without a
  'fixture-perturbation' marker.
- 'TEST.FIXTURES.PERTURBATION_REASON': a 'fixture-perturbation' marker with no
  reason.
- 'TEST.FIXTURES.INLINE_STUB': a 'chmod' in a fixture recipe on a path that no
  earlier 'cp' of a 'stub_*' file in that recipe wrote.

Lines inside a heredoc body, such as help text, or inside a multi-line quoted
string are not read as programs, and comment lines are skipped. Generated
fixtures are not read; their recipes are. The scan is lexical, so it covers the
common forms and leaves the rest to review.

A diagnostic whose count has reached zero blocks: a new finding of it fails the
run. The others are advisory and report their count, with no baseline file.
Prints one '<severity><tab><diagnostic><tab><path>:<line><tab><message>' row
per finding, then one 'COUNT<tab><diagnostic><tab><n><tab><mode>' row per
diagnostic, and exits 1 if a blocking diagnostic has a finding.
"""

from __future__ import annotations

import argparse
import ast
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path

from dev.audit.fixture_paths import is_fixture_path, is_fixture_recipe

RULE_SHELL_PYTHON = "SHELL.EMBED.PYTHON"
RULE_SHELL_BASH_C = "SHELL.EMBED.BASH_C"
RULE_SHELL_AWK = "SHELL.EMBED.AWK"
RULE_SHELL_SED = "SHELL.EMBED.SED"
RULE_SHELL_SED_GNU = "SHELL.EMBED.SED_GNU"
RULE_SHELL_PERL = "SHELL.EMBED.PERL"
RULE_PYTHON_STRING_INPUT = "PY.EMBED.STRING_INPUT"
RULE_RECORD_INPUT = "TEST.FIXTURES.RECORD_INPUT"
RULE_PERTURBATION_REASON = "TEST.FIXTURES.PERTURBATION_REASON"
RULE_INLINE_STUB = "TEST.FIXTURES.INLINE_STUB"

DIAGNOSTICS = (
    RULE_SHELL_PYTHON,
    RULE_SHELL_BASH_C,
    RULE_SHELL_AWK,
    RULE_SHELL_SED,
    RULE_SHELL_SED_GNU,
    RULE_SHELL_PERL,
    RULE_PYTHON_STRING_INPUT,
    RULE_RECORD_INPUT,
    RULE_PERTURBATION_REASON,
    RULE_INLINE_STUB,
)

# Diagnostics whose count has reached zero. Each blocks from then on; add one
# here when its count reaches zero.
BLOCKING = frozenset(
    {
        RULE_SHELL_PYTHON,
        RULE_SHELL_SED,
        RULE_SHELL_SED_GNU,
        RULE_PERTURBATION_REASON,
        RULE_INLINE_STUB,
    },
)

SHELL_ROOTS = ("bin", "lib", "install", "tests", "dev")
TEST_ROOTS = ("tests/unit", "tests/contract", "tests/integration")

MARK = "fixture-perturbation:"
MARK_RE = re.compile(r"#\s*fixture-perturbation:(?P<why>.*)$")
HEREDOC_RE = re.compile(r"<<-?\s*(['\"]?)(?P<tag>[A-Za-z_][A-Za-z0-9_]*)\1")
PYTHON_CALL_RE = re.compile(
    r"(?:\bpython3?\b|PYTHON[A-Z_]*\}?\"?|\$\{?py\}?\"?|py_cmd\[@\]\}\"?)"
    r"\s+(?:-I\s+)?(?P<flag>-c|-)(?:\s|$)",
)
BASH_C_RE = re.compile(r"(?:\bbash\b|\$\{?TEST_BASH\}?\"?)\s+-c\s+'")
AWK_RE = re.compile(r"\bg?awk\b(?P<rest>.*)$")
SED_RE = re.compile(r"\bsed\b(?P<rest>.*)$")
PERL_RE = re.compile(r"\bperl\s+-[a-z]*e")
SED_GNU_OPTION_RE = re.compile(r"(?:^|\s)(?:-[a-zA-Z]*[irz][a-zA-Z]*|--\w)")
SED_GNU_ESCAPE_RE = re.compile(r"\\[+|?]")
# A writer at command position: the start of a line, after any indentation, or
# after ';', '&&', '||', '{', '(', 'then', 'do', or 'else'.
CMD_START = r"(?:^\s*|[;{(]\s*|&&\s*|\|\|\s*|\b(?:then|do|else)\s+)"
WRITE_RE = re.compile(
    CMD_START
    + r"(?:(?:printf|echo|cat)\b[^#]*?(?<![0-9&])>>?\s*(?P<redir>[^&\s]\S*)"
    r"|:\s*>\s*(?P<trunc>\S+)"
    r"|(?:cp|touch|install)(?:\s|$))",
)
SINGLE_QUOTED_RE = re.compile(r"'[^']*'")
# A recipe installs a stub with 'cp' from a 'stub_*' file, then makes that same
# path executable; 'chmod' on any other path marks an inline stub.
STUB_COPY_RE = re.compile(
    CMD_START + r"cp\s+\S*/stub_\w+\.(?:sh|py)\"?\s+(?P<dest>\S+)\s*$",
)
CHMOD_RE = re.compile(CMD_START + r"chmod\s+(?P<args>.*)$")
HEREDOC_TO_FILE_RE = re.compile(r"<<-?.*>\s*\S|>\s*\S.*<<-?")
CAPTURE_RE = re.compile(r"/dev/null|>&2|log|err|stdout|stderr|_out\b")
PYTHON_WRITE_RE = re.compile(
    r"\.write_text\(|\.write_bytes\(|shutil\.copy"
    r"|\bopen\([^)]*['\"][wax]b?['\"]|gzip\.open\([^)]*['\"]w",
)

SIMPLE_STATEMENTS = (
    ast.Expr,
    ast.Assign,
    ast.AnnAssign,
    ast.AugAssign,
    ast.Return,
)


@dataclass(frozen=True)
class Finding:
    """
    One counted case.
    """

    rule_id: str
    path: str
    line: int
    message: str


def parse_args(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    """
    Parse the repository root and the output form.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Count programs embedded in other source, and inputs that tests "
            "write."
        ),
    )
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--json", action="store_true", help="Print findings as JSON."
    )

    return parser.parse_args(argv)


def iter_files(root: Path, bases: tuple[str, ...], suffix: str) -> list[Path]:
    """
    Return maintained files with a suffix under the given bases.
    """

    paths = set()

    for base in bases:
        top = root / base

        if top.is_dir():
            paths.update(top.rglob(f"*{suffix}"))

    return sorted(
        path
        for path in paths
        if "__pycache__" not in path.parts
        and not path.relative_to(root).as_posix().startswith("artifacts/")
        and path.is_file()
        and (
            not is_fixture_path(path, root)
            or is_fixture_recipe(path.relative_to(root).as_posix())
        )
    )


def is_test_source(rel: str) -> bool:
    """
    Return whether a path is a test, not a recipe, runner, or shared helper.
    """

    return rel.startswith(TEST_ROOTS) and Path(rel).name.startswith("test_")


def quoted_program(
    lines: list[str],
    index: int,
    start: int,
) -> tuple[int, int] | None:
    """
    Return a single-quoted program's line count and closing line.

    The quote opens at 'start' on line 'index'. Counts program lines only: the
    calling line counts when the program starts on it, and the closing line
    counts when program text comes before the quote.
    """

    head = lines[index][start + 1 :]
    if "'" in head:
        return 1, index

    count = 1 if head.strip() else 0

    for later in range(index + 1, len(lines)):
        text = lines[later]

        if "'" in text:
            if text[: text.index("'")].strip():
                count += 1

            return count, later

        count += 1

    return None


def program_quote(text: str) -> int | None:
    """
    Return where a tool's single-quoted program opens on a line.
    """

    for index, char in enumerate(text):
        if char == "'":
            return index

        if char == '"':
            return None

    return None


def next_quote_state(line: str, stack: list[str]) -> list[str]:
    """
    Return the quoting contexts still open at the end of a line.

    'stack' holds the contexts open at the line's start, innermost last: a
    single quote, a double quote, '$' for an ANSI-C '$'...'' string, or '(' for
    a command substitution or subshell, inside which the shell reads commands
    again. An unquoted '#' that starts a word ends the line as a comment.
    """

    stack = list(stack)
    index = 0

    while index < len(line):
        char = line[index]
        top = stack[-1] if stack else None

        if top in (None, "("):
            if char == "\\":
                index += 2

                continue

            if char == "#" and (index == 0 or line[index - 1].isspace()):
                return stack

            if char == "$" and line[index + 1 : index + 2] == "'":
                stack.append("$")
                index += 2

                continue

            if char == "$" and line[index + 1 : index + 2] == "(":
                stack.append("(")
                index += 2

                continue

            if char == "(":
                stack.append("(")
            elif char == ")" and top == "(":
                stack.pop()
            elif char in "'\"":
                stack.append(char)
        elif top == "'":
            if char == "'":
                stack.pop()
        elif top == "$":
            if char == "\\":
                index += 2

                continue

            if char == "'":
                stack.pop()
        elif top == '"':
            if char == "\\":
                index += 2

                continue

            if char == "$" and line[index + 1 : index + 2] == "(":
                stack.append("(")
                index += 2

                continue

            if char == '"':
                stack.pop()

        index += 1

    return stack


def scan_shell(root: Path, path: Path) -> list[Finding]:
    """
    Return the embedded programs and recorded inputs in one Shell file.
    """

    rel = path.relative_to(root).as_posix()
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    test = is_test_source(rel)
    recipe = is_fixture_recipe(rel)
    stub_paths: set[str] = set()
    found = []
    stack: list[str] = []
    tag = None

    def add(rule: str, line: int, message: str) -> None:
        found.append(Finding(rule, rel, line + 1, message))

    for index, line in enumerate(lines):
        if tag is not None:
            if line.strip() == tag:
                tag = None

            continue

        in_string = bool(stack) and stack[-1] != "("
        stack = next_quote_state(line, stack)
        stripped = line.strip()
        if in_string or not stripped or stripped[0] == "#":
            continue

        mark = MARK_RE.search(line)

        if mark and not mark.group("why").strip():
            add(RULE_PERTURBATION_REASON, index, "marker gives no reason")

        heredoc = HEREDOC_RE.search(line)
        python = PYTHON_CALL_RE.search(line)

        if heredoc and python and python.group("flag") == "-":
            add(RULE_SHELL_PYTHON, index, "Python fed through a heredoc")
        elif python and python.group("flag") == "-c":
            quote = line.find("'", python.end())

            if quote >= 0:
                result = quoted_program(lines, index, quote)

                if result and result[0] > 1:
                    add(
                        RULE_SHELL_PYTHON,
                        index,
                        f"'python -c' program of {result[0]} lines",
                    )

        bash_c = BASH_C_RE.search(line)

        if bash_c:
            result = quoted_program(lines, index, bash_c.end() - 1)

            if result and result[0] > 1:
                add(
                    RULE_SHELL_BASH_C,
                    index,
                    f"'bash -c' program of {result[0]} lines",
                )

        for regex, rule, limit, tool in (
            (AWK_RE, RULE_SHELL_AWK, 3, "awk"),
            (SED_RE, RULE_SHELL_SED, 1, "sed"),
        ):
            match = regex.search(line)
            if not match:
                continue

            quote = program_quote(match.group("rest"))
            if quote is None:
                continue

            start = match.start("rest") + quote
            result = quoted_program(lines, index, start)

            if result is None:
                continue

            if result[0] > limit:
                add(rule, index, f"'{tool}' program of {result[0]} lines")

            if tool == "sed":
                opts = match.group("rest")[:quote]
                end = line.find("'", start + 1)
                prog = line[start + 1 : end] if end > start else ""

                if SED_GNU_OPTION_RE.search(opts) or (
                    SED_GNU_ESCAPE_RE.search(prog)
                ):
                    add(RULE_SHELL_SED_GNU, index, "GNU-only 'sed' form")

        if PERL_RE.search(line):
            add(RULE_SHELL_PERL, index, "'perl' program")

        if recipe:
            copy = STUB_COPY_RE.search(line)

            if copy:
                stub_paths.add(copy.group("dest"))

            chmod = CHMOD_RE.search(line)

            if chmod and any(
                path not in stub_paths
                for path in chmod_paths(lines, index, chmod.group("args"))
            ):
                add(
                    RULE_INLINE_STUB,
                    index,
                    "'chmod' on a file not copied from a stub",
                )

        if test and not mark and not marked_above(lines, index):
            # A writer inside a quoted program, such as 'printf' in an 'awk'
            # program, is not the test writing a file.
            write = WRITE_RE.search(SINGLE_QUOTED_RE.sub("''", line))
            target = ""

            if write:
                target = write.group("redir") or write.group("trunc") or ""

            if HEREDOC_TO_FILE_RE.search(line) and not CAPTURE_RE.search(line):
                add(RULE_RECORD_INPUT, index, "heredoc written to a file")
            elif write and not CAPTURE_RE.search(target or line):
                add(RULE_RECORD_INPUT, index, "unmarked write in a test")

        if heredoc and (not stack or stack[-1] == "("):
            tag = heredoc.group("tag")

    return found


def chmod_paths(lines: list[str], index: int, args: str) -> list[str]:
    """
    Return the paths a 'chmod' names, across its backslash continuations.

    Drops the options and the mode, which come before the first path.
    """

    text = args
    later = index

    while text.rstrip().endswith("\\") and later + 1 < len(lines):
        later += 1
        text = text.rstrip()[:-1] + " " + lines[later]

    words = [word for word in text.split() if not word.startswith("-")]

    return words[1:]


def marked_above(lines: list[str], index: int) -> bool:
    """
    Return whether a continued command carries its marker on a later line.

    Follows the command's backslash continuations to its last line.
    """

    for later in range(index, len(lines)):
        if MARK in lines[later]:
            return True

        if not lines[later].rstrip().endswith("\\"):
            return False

    return False


def docstring_nodes(tree: ast.Module) -> set[int]:
    """
    Return the ids of string constants that are docstrings.
    """

    ids = set()

    for node in ast.walk(tree):
        if isinstance(
            node,
            (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef),
        ):
            body = node.body

            if (
                body
                and isinstance(body[0], ast.Expr)
                and isinstance(body[0].value, ast.Constant)
                and isinstance(body[0].value.value, str)
            ):
                ids.add(id(body[0].value))

    return ids


def scan_python_test(root: Path, path: Path) -> list[Finding]:
    """
    Return the long string literals and recorded inputs in one Python test.
    """

    rel = path.relative_to(root).as_posix()
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    found = []

    try:
        tree = ast.parse(text, filename=rel)
    except SyntaxError:
        return found

    docs = docstring_nodes(tree)

    # An f-string is one literal. Its own span is the same in every Python, but
    # the spans of the pieces inside it changed in 3.12, so they are skipped.
    pieces = {
        id(inner)
        for node in ast.walk(tree)
        if isinstance(node, ast.JoinedStr)
        for inner in ast.walk(node)
        if inner is not node
    }

    for node in ast.walk(tree):
        literal = isinstance(node, ast.JoinedStr) or (
            isinstance(node, ast.Constant) and isinstance(node.value, str)
        )

        if (
            literal
            and id(node) not in docs
            and id(node) not in pieces
            and node.end_lineno is not None
            and node.end_lineno - node.lineno + 1 >= 4
        ):
            found.append(
                Finding(
                    RULE_PYTHON_STRING_INPUT,
                    rel,
                    node.lineno,
                    "string literal of "
                    f"{node.end_lineno - node.lineno + 1} lines",
                ),
            )

    # A statement can span lines, and its marker may sit on any of them.
    marked_lines = set()

    for node in ast.walk(tree):
        if isinstance(node, SIMPLE_STATEMENTS) and node.end_lineno:
            span = range(node.lineno, node.end_lineno + 1)

            if any(MARK in lines[number - 1] for number in span):
                marked_lines.update(span)

    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("#"):
            continue

        mark = MARK_RE.search(line)

        if mark and not mark.group("why").strip():
            found.append(
                Finding(
                    RULE_PERTURBATION_REASON,
                    rel,
                    index + 1,
                    "marker gives no reason",
                ),
            )

        if PYTHON_WRITE_RE.search(line) and index + 1 not in marked_lines:
            found.append(
                Finding(
                    RULE_RECORD_INPUT,
                    rel,
                    index + 1,
                    "unmarked write in a test",
                ),
            )

    return found


def scan_repository(root: Path) -> list[Finding]:
    """
    Return every finding under the repository root.
    """

    found = []

    for path in iter_files(root, SHELL_ROOTS, ".sh"):
        found.extend(scan_shell(root, path))

    for path in iter_files(root, TEST_ROOTS, ".py"):
        if is_test_source(path.relative_to(root).as_posix()):
            found.extend(scan_python_test(root, path))

    return found


def main(argv: list[str] | None = None) -> int:
    """
    Print each finding and each diagnostic's count; exit 1 if any blocks.
    """

    args = parse_args(argv)
    root = args.root.resolve()
    found = scan_repository(root)
    counts = {rule: 0 for rule in DIAGNOSTICS}

    for item in found:
        counts[item.rule_id] += 1

    blocked = any(counts[rule] for rule in BLOCKING)

    if args.json:
        print(
            json.dumps(
                {
                    "findings": [asdict(item) for item in found],
                    "counts": counts,
                    "blocking": sorted(BLOCKING),
                },
                indent=2,
                sort_keys=True,
            ),
        )
    else:
        for item in found:
            severity = "FAIL" if item.rule_id in BLOCKING else "WARN"
            print(
                f"{severity}\t{item.rule_id}\t{item.path}:{item.line}\t"
                f"{item.message}",
            )

        for rule in DIAGNOSTICS:
            mode = "blocking" if rule in BLOCKING else "advisory"
            print(f"COUNT\t{rule}\t{counts[rule]}\t{mode}")

    return 1 if blocked else 0


if __name__ == "__main__":
    raise SystemExit(main())
