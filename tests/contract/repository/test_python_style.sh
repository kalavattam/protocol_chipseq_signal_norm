#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_python_style.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="python style"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


log_style="${TEST_DIR_LOG}/python_style/style.log"
mkdir -p "$(dirname "${log_style}")"

print_section "${TEST_NAME}"

if \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPATH="${ROOT_REPO}/src" \
    ROOT_REPO="${ROOT_REPO}" \
    python3 - << 'PY' > "${log_style}"
from __future__ import annotations

import ast
import os
import re
from pathlib import Path


ROOT = Path(os.environ["ROOT_REPO"])
SCRIPT_DIRS = [ROOT / "src/protocol_chipseq_signal_norm/cli"]
DOCSTRING_DIRS = [ROOT / "src/protocol_chipseq_signal_norm"]
SNAKE_RE = re.compile(r"^[a-z][a-z0-9_]*$")
GOOGLE_RE = re.compile(r"^\s*(Args|Returns|Raises):\s*$", re.M)
NOOP_RAISES_RE = re.compile(
    r"^\s*Raises\s*\n\s*-+\s*\n\s*(None|Nothing|No meaningful exceptions)\.?\s*$",
    re.I | re.M,
)
DOCSTRING_MINI_GRAMMAR_RE = re.compile(
    r"^\s*[A-Za-z_][A-Za-z0-9_]*\s+:\s+(?:enum:|csv:|spec\b|flt\b|num\b)",
    re.M,
)
DOCSTRING_CHOICE_TYPE_RE = re.compile(
    r"^\s*[A-Za-z_][A-Za-z0-9_]*\s+:\s+(\{[^}\n]+\})",
    re.M,
)
SECTION_SPACING_RE = re.compile(
    r"^(?P<indent>[ \t]*)"
    r"(?P<section>Usage|Parameters|Other Parameters|Returns|Yields|Receives|"
    r"Raises|Warns|Warnings|See Also|Notes|References|Examples|Attributes|"
    r"Methods|Script|Description|Arguments|Output|Testing|#TODO)\n"
    r"[ \t]*\n"
    r"(?P=indent)-+\n",
    re.M,
)
MODULE_SECTION_RE = re.compile(
    r"^(Usage|Parameters|Returns|Notes|References|See Also|Examples)\n-+\n",
    re.M,
)
MODULE_RETIRED_SECTION_RE = re.compile(
    r"^(Script|Description|Arguments|Output|Testing|#TODO)\n-+\n",
    re.M,
)
MODULE_ANY_SECTION_RE = re.compile(r"^([A-Z][A-Za-z ]+|#TODO)\n-+\n", re.M)
MODULE_DOC_MAX_LINES = 120
DOCSTRING_CLI_INVENTORY_RE = re.compile(
    r"Supported CLI options\s*:",
    re.I,
)


class Finding:
    def __init__(self, severity: str, path: Path, line: int, message: str) -> None:
        self.severity = severity
        self.path = path
        self.line = line
        self.message = message

    def emit(self) -> None:
        rel = self.path.relative_to(ROOT)
        print(f"{self.severity}:{rel}:{self.line}:{self.message}")


def ann_text(node: ast.AST | None) -> str:
    return "" if node is None else ast.unparse(node)


def literal_str(node: ast.AST) -> str | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    return None


def keyword(call: ast.Call, name: str) -> ast.keyword | None:
    for kw in call.keywords:
        if kw.arg == name:
            return kw
    return None


def is_nonempty_help(call: ast.Call) -> bool:
    kw = keyword(call, "help")
    if kw is None:
        return False
    text = literal_str(kw.value)
    if text is not None:
        return bool(text.strip())
    return True


def is_argparse_suppress(node: ast.AST) -> bool:
    return (
        isinstance(node, ast.Attribute)
        and node.attr == "SUPPRESS"
        and isinstance(node.value, ast.Name)
        and node.value.id == "argparse"
    )


def is_hidden_compat(call: ast.Call) -> bool:
    kw = keyword(call, "help")
    if kw is None or not is_argparse_suppress(kw.value):
        return False
    opts = option_strings(call)
    return bool(opts) and all(is_hyphen_duplicate(opt) for opt in opts)


def is_hyphen_duplicate(opt: str) -> bool:
    return opt.startswith("--") and "-" in opt[2:] and "_" not in opt[2:]


def add_argument_calls(func: ast.FunctionDef) -> list[ast.Call]:
    calls: list[ast.Call] = []
    for node in ast.walk(func):
        if not isinstance(node, ast.Call):
            continue
        callee = node.func
        if isinstance(callee, ast.Attribute) and callee.attr == "add_argument":
            calls.append(node)
    return calls


def option_strings(call: ast.Call) -> list[str]:
    opts: list[str] = []
    for arg in call.args:
        text = literal_str(arg)
        if text is None:
            break
        opts.append(text)
    return opts


def check_parse_args(path: Path, tree: ast.Module, findings: list[Finding]) -> None:
    funcs = {
        node.name: node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
    }
    parse_args = funcs.get("parse_args")
    if parse_args is None:
        findings.append(Finding("FAIL", path, 1, "missing parse_args()"))
        return

    if len(parse_args.args.args) != 1 or parse_args.args.args[0].arg != "argv":
        findings.append(
            Finding("FAIL", path, parse_args.lineno, "parse_args() must accept argv")
        )
    elif ann_text(parse_args.args.args[0].annotation) != "list[str] | None":
        findings.append(
            Finding(
                "FAIL",
                path,
                parse_args.lineno,
                "parse_args() argv must be typed list[str] | None",
            )
        )

    default_count = len(parse_args.args.defaults)
    if default_count != 1 or not isinstance(parse_args.args.defaults[0], ast.Constant):
        findings.append(
            Finding(
                "FAIL",
                path,
                parse_args.lineno,
                "parse_args() argv must default to None",
            )
        )
    elif parse_args.args.defaults[0].value is not None:
        findings.append(
            Finding(
                "FAIL",
                path,
                parse_args.lineno,
                "parse_args() argv must default to None",
            )
        )

    if ann_text(parse_args.returns) != "argparse.Namespace":
        findings.append(
            Finding(
                "FAIL",
                path,
                parse_args.lineno,
                "parse_args() must return argparse.Namespace",
            )
        )

    if not ast.get_docstring(parse_args):
        findings.append(
            Finding("FAIL", path, parse_args.lineno, "parse_args() needs a docstring")
        )

    has_cap_parser = any(
        isinstance(node, ast.Call)
        and (
            (isinstance(node.func, ast.Name) and node.func.id == "CapArgumentParser")
            or (
                isinstance(node.func, ast.Attribute)
                and node.func.attr == "CapArgumentParser"
            )
        )
        for node in ast.walk(parse_args)
    )
    if not has_cap_parser:
        findings.append(
            Finding("FAIL", path, parse_args.lineno, "parse_args() must use CapArgumentParser")
        )

    has_help_cap = any(
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "add_help_cap"
        for node in ast.walk(parse_args)
    )
    if not has_help_cap:
        findings.append(
            Finding("FAIL", path, parse_args.lineno, "parse_args() must call add_help_cap()")
        )

    for call in add_argument_calls(parse_args):
        opts = option_strings(call)
        if not opts:
            continue
        is_optional = all(opt.startswith("-") for opt in opts)
        hidden_compat = is_optional and is_hidden_compat(call)
        if is_optional and not hidden_compat:
            if "--ref" in opts:
                findings.append(
                    Finding("FAIL", path, call.lineno, "use --ref_fa, not --ref")
                )
            if "--ref_fa" in opts and "-rf" not in opts:
                findings.append(
                    Finding("FAIL", path, call.lineno, "use -rf, --ref_fa for reference FASTA")
                )
        if hidden_compat:
            if not opts or any(not is_hyphen_duplicate(opt) for opt in opts):
                findings.append(
                    Finding(
                        "FAIL",
                        path,
                        call.lineno,
                        "hidden compatibility aliases must be hyphenated long options",
                    )
                )
        elif not is_nonempty_help(call):
            findings.append(
                Finding("FAIL", path, call.lineno, "parser argument needs non-empty help=")
            )
        if is_optional:
            short_opts = [
                opt for opt in opts
                if re.match(r"^-[A-Za-z0-9][A-Za-z0-9_-]*$", opt)
            ]
            long_opts = [opt for opt in opts if opt.startswith("--")]
            if len(short_opts) > 2 or len(long_opts) > 2:
                findings.append(
                    Finding(
                        "FAIL",
                        path,
                        call.lineno,
                        "optional argument exceeds alias cap "
                        f"(short={len(short_opts)} long={len(long_opts)})",
                    )
                )
            if keyword(call, "dest") is None:
                findings.append(
                    Finding("FAIL", path, call.lineno, "optional argument needs explicit dest=")
                )
            else:
                dest = literal_str(keyword(call, "dest").value)  # type: ignore[union-attr]
                if dest == "ref":
                    findings.append(
                        Finding("FAIL", path, call.lineno, "use dest=\"ref_fa\", not dest=\"ref\"")
                    )
                if dest == "r" "nd":
                    findings.append(
                        Finding("FAIL", path, call.lineno, "use dest=\"dp\", not dest=\"r" "nd\"")
                    )
            if hidden_compat:
                continue
            if not any(re.match(r"^-[A-Za-z0-9][A-Za-z0-9_-]*$", opt) for opt in opts):
                findings.append(
                    Finding("FAIL", path, call.lineno, "optional argument needs a short form")
                )
            if not any(opt.startswith("--") for opt in opts):
                findings.append(
                    Finding("FAIL", path, call.lineno, "optional argument needs a long form")
                )
        else:
            name = opts[0]
            if not SNAKE_RE.match(name):
                findings.append(
                    Finding("FAIL", path, call.lineno, f"positional '{name}' is not snake_case")
                )


def check_main(path: Path, tree: ast.Module, findings: list[Finding]) -> None:
    funcs = {
        node.name: node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
    }
    main = funcs.get("main")
    if main is None:
        findings.append(Finding("FAIL", path, 1, "missing main()"))
        return

    if len(main.args.args) != 1 or main.args.args[0].arg != "argv":
        findings.append(Finding("FAIL", path, main.lineno, "main() must accept argv"))
    elif ann_text(main.args.args[0].annotation) != "list[str] | None":
        findings.append(
            Finding("FAIL", path, main.lineno, "main() argv must be typed list[str] | None")
        )

    if ann_text(main.returns) != "int":
        findings.append(Finding("FAIL", path, main.lineno, "main() must return int"))

    if not ast.get_docstring(main):
        findings.append(Finding("FAIL", path, main.lineno, "main() needs a docstring"))

    has_return_zero = any(
        isinstance(node, ast.Return)
        and isinstance(node.value, ast.Constant)
        and node.value.value == 0
        for node in ast.walk(main)
    )
    if not has_return_zero:
        findings.append(Finding("FAIL", path, main.lineno, "main() needs explicit return 0"))


def check_entrypoint(path: Path, tree: ast.Module, findings: list[Finding]) -> None:
    src = path.read_text()
    if 'argparse.ArgumentParser' in src:
        findings.append(Finding("FAIL", path, 1, "use CapArgumentParser, not argparse.ArgumentParser"))
    if "args.r" "nd" in src:
        findings.append(Finding("FAIL", path, 1, "use args.dp, not args.r" "nd"))
    if re.search(r"\brnd\s*=", src):
        findings.append(Finding("FAIL", path, 1, "use dp, not r" "nd, for decimal precision variables"))
    if "sys.exit(main())" in src:
        findings.append(Finding("FAIL", path, 1, "use raise SystemExit(main())"))

    has_guard = False
    has_system_exit = False
    for node in tree.body:
        if not isinstance(node, ast.If):
            continue
        test = ast.unparse(node.test)
        if test == "__name__ == '__main__'" or test == "__name__ == \"__main__\"":
            has_guard = True
            for stmt in ast.walk(node):
                if not isinstance(stmt, ast.Raise) or stmt.exc is None:
                    continue
                if ast.unparse(stmt.exc) == "SystemExit(main())":
                    has_system_exit = True
    if not has_guard:
        findings.append(Finding("FAIL", path, 1, "missing __main__ entrypoint guard"))
    elif not has_system_exit:
        findings.append(Finding("FAIL", path, 1, "entrypoint must raise SystemExit(main())"))


def check_docstrings(path: Path, tree: ast.Module, findings: list[Finding]) -> None:
    module_doc = ast.get_docstring(tree)
    if not module_doc:
        findings.append(Finding("WARN", path, 1, "module docstring is missing"))
    elif "\n" in module_doc and not MODULE_SECTION_RE.search(module_doc):
        findings.append(
            Finding("WARN", path, 1, "module docstring lacks NumPy/SciPy sections")
        )
    if module_doc:
        if DOCSTRING_CLI_INVENTORY_RE.search(module_doc):
            findings.append(
                Finding(
                    "WARN",
                    path,
                    1,
                    "module docstring contains a CLI option inventory; use parser help",
                )
            )
        for match in MODULE_RETIRED_SECTION_RE.finditer(module_doc):
            line = module_doc.count("\n", 0, match.start()) + 1
            findings.append(
                Finding(
                    "WARN",
                    path,
                    line,
                    f"module docstring uses retired section '{match.group(1)}'",
                )
            )
        sections: dict[str, int] = {}
        for match in MODULE_ANY_SECTION_RE.finditer(module_doc):
            section = match.group(1)
            sections[section] = sections.get(section, 0) + 1
        for section, count in sorted(sections.items()):
            if count > 1:
                findings.append(
                    Finding(
                        "WARN",
                        path,
                        1,
                        f"module docstring has duplicate section '{section}'",
                    )
                )
        line_count = module_doc.count("\n") + 1
        if line_count > MODULE_DOC_MAX_LINES:
            findings.append(
                Finding(
                    "WARN",
                    path,
                    1,
                    "module docstring is long; move method notes to docs/dev",
                )
            )

    for node in tree.body:
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        if isinstance(node, ast.FunctionDef) and node.name.startswith("_"):
            continue
        doc = ast.get_docstring(node)
        if not doc:
            findings.append(
                Finding("FAIL", path, node.lineno, f"{node.name} needs a docstring")
            )
            continue
        if DOCSTRING_CLI_INVENTORY_RE.search(doc):
            findings.append(
                Finding(
                    "WARN",
                    path,
                    node.lineno,
                    f"{node.name} docstring contains a CLI option inventory; use parser help",
                )
            )
        if GOOGLE_RE.search(doc):
            findings.append(
                Finding("WARN", path, node.lineno, f"{node.name} uses Google-style sections")
            )
        if NOOP_RAISES_RE.search(doc):
            findings.append(
                Finding("WARN", path, node.lineno, f"{node.name} has a no-op Raises section")
            )
        if DOCSTRING_MINI_GRAMMAR_RE.search(doc):
            findings.append(
                Finding(
                    "WARN",
                    path,
                    node.lineno,
                    f"{node.name} uses compact mini-grammar in a docstring type",
                )
            )
        for match in DOCSTRING_CHOICE_TYPE_RE.finditer(doc):
            choice = match.group(1)
            if '"' in choice or re.search(r",[^\s}]", choice):
                findings.append(
                    Finding(
                        "WARN",
                        path,
                        node.lineno,
                        f"{node.name} has a nonstandard choice-set type",
                    )
                )
                break


def check_docstring_section_spacing(path: Path, findings: list[Finding]) -> None:
    text = path.read_text()
    for match in SECTION_SPACING_RE.finditer(text):
        line = text.count("\n", 0, match.start()) + 1
        section = match.group("section")
        findings.append(
            Finding(
                "FAIL",
                path,
                line,
                f"{section} heading has a blank line before its underline",
            )
        )


def iter_user_scripts() -> list[Path]:
    files: list[Path] = []
    for directory in SCRIPT_DIRS:
        if not directory.exists():
            continue
        files.extend(
            path for path in directory.glob("*.py")
            if path.name != "__init__.py"
        )
    return sorted(files)


def iter_docstring_files() -> list[Path]:
    files: list[Path] = []
    for directory in DOCSTRING_DIRS:
        if not directory.exists():
            continue
        files.extend(
            path for path in directory.glob("*.py")
            if path.name != "__init__.py"
        )
    return sorted(set(files))


all_findings: list[Finding] = []
for path in iter_docstring_files():
    check_docstring_section_spacing(path, all_findings)

for path in iter_user_scripts():
    try:
        tree = ast.parse(path.read_text(), filename=str(path))
    except SyntaxError as exc:
        all_findings.append(
            Finding("FAIL", path, exc.lineno or 1, f"syntax error: {exc.msg}")
        )
        continue
    check_parse_args(path, tree, all_findings)
    check_main(path, tree, all_findings)
    check_entrypoint(path, tree, all_findings)
    check_docstrings(path, tree, all_findings)

for finding in all_findings:
    finding.emit()

if not any(f.severity == "FAIL" for f in all_findings):
    print("PASS:python CLI and docstring style")
PY
then
    :
else
    record_fail "python style scanner crashed; see $(print_relpath "${log_style}")"
fi

while IFS=: read -r sev rel line msg; do
    case "${sev}" in
        PASS)
            record_pass "${rel}"
            ;;
        WARN)
            record_warn "${rel}:${line}:${msg}"
            ;;
        FAIL)
            record_fail "${rel}:${line}:${msg}"
            ;;
    esac
done < "${log_style}"

finish
