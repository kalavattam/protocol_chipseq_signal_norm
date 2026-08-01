#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_parameter_docs_consistency.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="parameter description consistency"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


log_style="${TEST_DIR_LOG}/parameter_docs_consistency/style.log"
mkdir -p "$(dirname "${log_style}")"

print_section "${TEST_NAME}"

if \
    PYTHONDONTWRITEBYTECODE=1 \
    ROOT_REPO="${ROOT_REPO}" \
    python3 - << 'PY' > "${log_style}"
from __future__ import annotations

import ast
import copy
import json
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(os.environ["ROOT_REPO"])
DOC_HELP = ROOT / "docs" / "standards" / "help.md"
CONTRACTS = ROOT / os.environ.get(
    "HELP_CONTRACTS_CONFIG",
    "dev/config/help_contracts.json",
)
SCAN_DIRS = [
    ROOT / "bin",
    ROOT / "lib",
    ROOT / "install",
    ROOT / "src",
    ROOT / "scripts",
    ROOT / "tests" / "scripts",
    ROOT / "docs" / "standards",
]
OLD_IN = "in" + "file"
OLD_OUT = "out" + "file"
OLD_NAMES = [
    f"csv_{OLD_IN}",
    f"csv_{OLD_OUT}",
    f"arr_{OLD_IN}",
    f"arr_{OLD_OUT}",
    OLD_IN,
    OLD_OUT,
    f"{OLD_IN}s",
    f"{OLD_OUT}s",
]
SHARED_NAME_CANDIDATES = {
    "align_typ",
    "aligner",
    "aln_typ",
    "bt2_mode",
    "bwa_alg",
    "cfg_met",
    "chk_chr",
    "coef",
    "csv_min",
    "csv_mip",
    "csv_sin",
    "csv_sip",
    "dep_mip",
    "dep_sin",
    "dep_sip",
    "dir_fnc",
    "drp_nan",
    "eps",
    "eqn",
    "fil_A",
    "fil_B",
    "fil_aln",
    "floor",
    "fq_1",
    "fq_2",
    "index",
    "len_def",
    "len_min",
    "len_mip",
    "log_err",
    "log_out",
    "mapq",
    "mito",
    "mode_nz",
    "mtr",
    "pseudo",
    "pth_scr_py",
    "qname",
    "qntl_nz",
    "req_flg",
    "retain",
    "sfx_pe",
    "sfx_se",
    "siz_gen",
    "suffix_pe",
    "suffix_se",
    "tbl_met",
    "tg",
    "typ_out",
    "usr_frg",
}
RETIRED_RE = re.compile(
    r"\b(?:" + "|".join(re.escape(name) for name in OLD_NAMES) + r")\b"
    r"|--(?:csv_)?(?:" + OLD_IN + "|" + OLD_OUT + r")\b"
)
PARAM_ROW_RE = re.compile(
    r"^(?P<indent>[ \t]{2,})(?P<head>[^:\n]+?)\s+:\s+(?P<typ>[^\n]+)$"
)
SECTION_RE = re.compile(r"^[A-Z][A-Za-z ]+\n-+\n", re.M)


@dataclass(frozen=True)
class ParamDoc:
    path: Path
    line: int
    key: str
    desc: str
    symbol: str


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def norm_desc(text: str) -> str:
    text = re.sub(r"\s+", " ", text.strip())
    return text.rstrip()


def norm_compare(text: str) -> str:
    text = norm_desc(text)
    text = re.sub(
        r"\s*\(default:\s*(?:%\([^)]+\)s|[^)]+)\)",
        "",
        text,
    )
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def has_shared_core(desc: str, expected: str) -> bool:
    desc_cmp = norm_compare(desc)
    exp_cmp = norm_compare(expected)
    exp_base = exp_cmp.rstrip(".")
    if desc_cmp == exp_cmp:
        return True
    if desc_cmp.startswith(exp_cmp + " "):
        return True
    if desc_cmp.startswith(exp_cmp + "\n"):
        return True
    if desc_cmp.startswith(exp_base):
        nxt = desc_cmp[len(exp_base):len(exp_base) + 1]
        return nxt in {"", ".", ";", ":", ",", " "}
    return False


def opening_words(text: str, n: int = 3) -> str:
    words = re.findall(r"[A-Za-z0-9_]+", text.lower())
    return " ".join(words[:n])


def duplicated_registered_core(desc: str, expected: str) -> str | None:
    desc_cmp = norm_compare(desc)
    exp_cmp = norm_compare(expected)
    exp_base = re.escape(exp_cmp.rstrip("."))
    core_hits = re.findall(rf"\b{exp_base}\b", desc_cmp, flags=re.I)
    if len(core_hits) > 1:
        return f"repeats shared core '{expected}'"

    sentences = [
        sent.strip()
        for sent in re.split(r"(?<=[.!?])\s+", desc_cmp)
        if sent.strip()
    ]
    if len(sentences) < 2:
        return None
    first_words = opening_words(sentences[0])
    second_words = opening_words(sentences[1])
    if first_words and first_words == second_words:
        return (
            "repeats the opening phrase across adjacent sentences: "
            f"'{sentences[0]}' / '{sentences[1]}'"
        )
    return None


def canonical_from_row(head: str) -> str | None:
    head = head.strip()
    opts = re.findall(r"--([A-Za-z][A-Za-z0-9_]*)\b", head)
    if opts:
        return opts[-1]
    m_pos = re.match(r"\d+\s+([A-Za-z][A-Za-z0-9_]*)$", head)
    if m_pos:
        return m_pos.group(1)
    m_name = re.match(r"([A-Za-z][A-Za-z0-9_]*)$", head)
    if m_name:
        return m_name.group(1)
    return None


def read_canonical_table() -> dict[str, str]:
    registry: dict[str, str] = {}
    lines = DOC_HELP.read_text(encoding="utf-8").splitlines()
    start = next(
        (
            index
            for index, line in enumerate(lines)
            if line.startswith("## ")
            and "`PARAMETER.DESCRIPTIONS`" in line
        ),
        None,
    )
    if start is None:
        raise RuntimeError(
            "missing PARAMETER.DESCRIPTIONS owner section"
        )
    end = next(
        (
            index
            for index in range(start + 1, len(lines))
            if lines[index].startswith("## ")
        ),
        len(lines),
    )
    for line in lines[start + 1:end]:
        if line.startswith("Use `"):
            break
        if not line.startswith("| `"):
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) < 3:
            continue
        key = cells[0].strip("`")
        desc = norm_desc(cells[2])
        if key and desc and key != "Parameter":
            registry[key] = desc
    if not registry:
        raise RuntimeError(
            "PARAMETER.DESCRIPTIONS contains no canonical rows"
        )
    return registry


def family_registry_findings(
    table: dict[str, str],
    families: list[dict[str, object]],
) -> list[str]:
    """
    Return exact bidirectional canonical-table/family failures.
    """

    findings: list[str] = []
    shared = [
        family
        for family in families
        if family.get("status") == "registered_shared"
    ]
    local = [
        family
        for family in families
        if family.get("status") == "non_applicable_same_name"
    ]
    by_parameter: dict[str, list[dict[str, object]]] = defaultdict(list)
    for family in shared:
        by_parameter[str(family.get("parameter", ""))].append(family)

    for parameter, core in table.items():
        matches = by_parameter.get(parameter, [])
        if len(matches) != 1:
            findings.append(
                f"canonical row '{parameter}' maps to {len(matches)} "
                "registered_shared families"
            )
            continue
        if norm_desc(str(matches[0].get("canonical_core", ""))) != core:
            findings.append(
                f"canonical row '{parameter}' core differs from its family"
            )

    for parameter, matches in sorted(by_parameter.items()):
        if len(matches) != 1:
            findings.append(
                f"registered_shared family '{parameter}' occurs "
                f"{len(matches)} times"
            )
        if parameter not in table:
            findings.append(
                f"registered_shared family '{parameter}' has no canonical row"
            )

    for family in local:
        parameter = str(family.get("parameter", ""))
        if parameter in table:
            findings.append(
                f"non_applicable_same_name family '{parameter}' appears in "
                "the canonical table"
            )
        meanings = family.get("local_meanings")
        if not isinstance(meanings, list) or not meanings:
            findings.append(
                f"non_applicable_same_name family '{parameter}' lacks local "
                "meaning evidence"
            )

    return findings


def read_registered_families(
    table: dict[str, str],
) -> tuple[dict[str, str], list[dict[str, object]]]:
    """
    Load the shared family registry consumed by this semantic checker.
    """

    data = json.loads(CONTRACTS.read_text(encoding="utf-8"))
    families = data.get("parameter_families")
    if not isinstance(families, list):
        raise RuntimeError("help contracts lack parameter_families")
    for message in family_registry_findings(table, families):
        print(f"FAIL:{rel(CONTRACTS)}:1:{message}")
    registry = {
        str(family["parameter"]): norm_desc(str(family["canonical_core"]))
        for family in families
        if family.get("status") == "registered_shared"
    }
    return registry, families


def prove_family_fault_detection(
    table: dict[str, str],
    families: list[dict[str, object]],
) -> None:
    """
    Prove the bidirectional connection rejects representative divergence.
    """

    shared_index = next(
        index
        for index, family in enumerate(families)
        if family.get("status") == "registered_shared"
    )
    non_applicable_index = next(
        index
        for index, family in enumerate(families)
        if family.get("status") == "non_applicable_same_name"
    )
    mutations: list[tuple[str, dict[str, str], list[dict[str, object]]]] = []

    missing = copy.deepcopy(families)
    del missing[shared_index]
    mutations.append(("missing family", table, missing))

    duplicate = copy.deepcopy(families)
    duplicate.append(copy.deepcopy(duplicate[shared_index]))
    mutations.append(("duplicate family", table, duplicate))

    mismatched = copy.deepcopy(families)
    mismatched[shared_index]["canonical_core"] = "Injected mismatch."
    mutations.append(("core mismatch", table, mismatched))

    promoted_table = dict(table)
    local_parameter = str(families[non_applicable_index]["parameter"])
    promoted_table[local_parameter] = "Injected shared meaning."
    mutations.append(("non-applicable promotion", promoted_table, families))

    for label, candidate_table, candidate_families in mutations:
        if not family_registry_findings(
            candidate_table,
            candidate_families,
        ):
            print(
                f"FAIL:{rel(CONTRACTS)}:1:fault injection was not detected: "
                f"{label}"
            )


def iter_text_files() -> list[Path]:
    files: list[Path] = []
    for base in SCAN_DIRS:
        if not base.exists():
            continue
        for path in base.rglob("*"):
            if not path.is_file():
                continue
            if "artifacts/tests" in path.as_posix():
                continue
            if path.suffix not in {".sh", ".py", ".md"}:
                continue
            files.append(path)
    return sorted(set(files))


def scan_retired_names(paths: list[Path]) -> list[str]:
    findings: list[str] = []
    for path in paths:
        if path.name == "test_parameter_docs_consistency.sh":
            continue
        text = path.read_text(encoding="utf-8")
        for idx, line in enumerate(text.splitlines(), 1):
            if RETIRED_RE.search(line):
                findings.append(
                    f"FAIL:{rel(path)}:{idx}:"
                    "retired input/output parameter name: "
                    f"{line.strip()}"
                )
    return findings


def shell_param_docs(path: Path) -> list[ParamDoc]:
    docs: list[ParamDoc] = []
    lines = path.read_text(encoding="utf-8").splitlines()
    in_params = False
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.strip() == "Parameters":
            in_params = True
            i += 2
            continue
        if in_params and re.match(r"^[A-Z][A-Za-z ]+$", line.strip()):
            in_params = False
        if not in_params:
            i += 1
            continue
        m_row = PARAM_ROW_RE.match(line)
        if not m_row:
            i += 1
            continue
        key = canonical_from_row(m_row.group("head"))
        if key is None:
            i += 1
            continue
        desc_lines: list[str] = []
        j = i + 1
        while j < len(lines):
            nxt = lines[j]
            if PARAM_ROW_RE.match(nxt):
                break
            if re.match(r"^[A-Z][A-Za-z ]+$", nxt.strip()):
                break
            if nxt.startswith("    ") or nxt.startswith("\t"):
                stripped = nxt.strip()
                if stripped:
                    desc_lines.append(stripped)
                elif desc_lines:
                    break
            elif desc_lines:
                break
            j += 1
        if desc_lines:
            desc = norm_desc(" ".join(desc_lines))
            symbol = "<file>"
            for prior in reversed(lines[:i + 1]):
                match = re.match(
                    r"(?:function\s+)?([A-Za-z_][A-Za-z0-9_]*)"
                    r"\s*(?:\(\))?\s*\{",
                    prior.strip(),
                )
                if match:
                    symbol = match.group(1)
                    break
            docs.append(ParamDoc(path, i + 1, key, desc, symbol))
        i = max(j, i + 1)
    return docs


def literal_str(node: ast.AST) -> str | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.JoinedStr):
        return None
    try:
        value = ast.literal_eval(node)
    except Exception:
        return None
    return value if isinstance(value, str) else None


def keyword(call: ast.Call, name: str) -> ast.keyword | None:
    for kw in call.keywords:
        if kw.arg == name:
            return kw
    return None


def option_strings(call: ast.Call) -> list[str]:
    opts: list[str] = []
    for arg in call.args:
        text = literal_str(arg)
        if text is None:
            break
        opts.append(text)
    return opts


def python_param_docs(path: Path) -> list[ParamDoc]:
    docs: list[ParamDoc] = []
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except SyntaxError as exc:
        line = exc.lineno or 1
        print(f"FAIL:{rel(path)}:{line}:Python syntax error: {exc.msg}")
        return docs
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        callee = node.func
        if not (
            isinstance(callee, ast.Attribute)
            and callee.attr == "add_argument"
        ):
            continue
        opts = option_strings(node)
        if not opts or all(
            opt.startswith("--") and "-" in opt[2:]
            for opt in opts
        ):
            continue
        dest_kw = keyword(node, "dest")
        key = literal_str(dest_kw.value) if dest_kw is not None else None
        if key is None:
            longs = [opt[2:] for opt in opts if opt.startswith("--")]
            key = longs[-1].replace("-", "_") if longs else opts[0]
        help_kw = keyword(node, "help")
        if help_kw is None:
            continue
        desc = literal_str(help_kw.value)
        if desc is None or desc.strip() == "argparse.SUPPRESS":
            continue
        symbol = "<file>"
        for candidate in ast.walk(tree):
            if not isinstance(
                candidate,
                (ast.FunctionDef, ast.AsyncFunctionDef),
            ):
                continue
            end_line = getattr(candidate, "end_lineno", candidate.lineno)
            if candidate.lineno <= node.lineno <= end_line:
                symbol = candidate.name
                break
        docs.append(
            ParamDoc(path, node.lineno, key, norm_desc(desc), symbol)
        )
    return docs


def approved_realization(member: dict[str, object], fallback: str) -> str:
    """Return the exact registered natural realization, never a name guess."""

    explicit = member.get("approved_realization")
    if isinstance(explicit, str) and explicit:
        return norm_desc(explicit)
    evidence = str(member.get("evidence", ""))
    marker = ": "
    return norm_desc(evidence.rsplit(marker, 1)[-1]) if marker in evidence else fallback


canonical_table = read_canonical_table()
registry, families = read_registered_families(canonical_table)
prove_family_fault_detection(canonical_table, families)
paths = iter_text_files()
for finding in scan_retired_names(paths):
    print(finding)

docs: list[ParamDoc] = []
for path in paths:
    if path.suffix == ".py":
        docs.extend(python_param_docs(path))
    elif path.suffix == ".sh":
        docs.extend(shell_param_docs(path))

for family in families:
    if family.get("status") != "registered_shared":
        continue
    parameter = str(family["parameter"])
    expected = registry[parameter]
    for member in family.get("members", []):
        member_path = ROOT / str(member["path"])
        matching = [
            doc
            for doc in docs
            if doc.path == member_path
            and doc.key == parameter
            and doc.symbol == str(member["symbol"])
        ]
        if len(matching) != 1:
            print(
                f"FAIL:{rel(CONTRACTS)}:1:registered member "
                f"'{member['surface_id']}' maps to {len(matching)} "
                "semantic parameter rows"
            )
        else:
            approved = approved_realization(member, expected)
            duplicated = duplicated_registered_core(matching[0].desc, approved)
            if duplicated is not None:
                print(
                    f"FAIL:{rel(matching[0].path)}:{matching[0].line}:"
                    f"malformed description for {parameter}: {duplicated}"
                )
            elif norm_compare(matching[0].desc) != norm_compare(approved):
                print(
                    f"WARN:{rel(matching[0].path)}:{matching[0].line}:"
                    f"PARAMETER.DESCRIPTIONS registered-member divergence for "
                    f"'{member['surface_id']}'"
                )
PY
then
    while IFS=: read -r sev file line msg; do
        [[ -n "${sev}" ]] || continue
        case "${sev}" in
            FAIL)
                record_fail "${file}:${line}:${msg}"
                ;;
            WARN)
                record_warn "${file}:${line}:${msg}"
                ;;
            *)
                record_fail \
                    "unrecognized parameter-doc finding:" \
                    "${sev}:${file}:${line}:${msg}"
                ;;
        esac
    done < "${log_style}"
else
    record_fail "parameter description consistency scanner crashed"
fi

if [[ ! -s "${log_style}" ]]; then
    record_pass "no parameter description consistency findings"
fi

finish
