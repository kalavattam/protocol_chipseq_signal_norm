#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_numeric_emission_inventory.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.

"""
Test bounded static numeric emission inventory claims.
"""

from __future__ import annotations

import json
from pathlib import Path

from dev.audit.numeric_emission_inventory import inventory
from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[3]
SCHEMA = ROOT / "dev/schemas/numeric_emission_inventory.schema.json"


def fixture_inventory(tmp_path: Path) -> dict:
    (tmp_path / "src").mkdir()
    (tmp_path / "bin").mkdir()
    (tmp_path / "src/example.py").write_text(
        "\n".join(
            [
                "import csv",
                "import json",
                "import sys",
                "from pathlib import Path",
                'print(f"{value:.{dp}f}")',
                'sys.stdout.write(f"{value}\\n")',
                'with open("values.txt", "w") as stream:',
                '    stream.write(f"{value}\\n")',
                "    writer = csv.writer(stream)",
                "    writer.writerow([value])",
                "    json.dump(value, stream)",
                'Path("value.txt").write_text(str(value))',
                "getattr(target, method)(value)",
                "eval(expression)",
                "",
            ],
        ),
        encoding="utf-8",
    )
    (tmp_path / "bin/example.sh").write_text(
        "\n".join(
            [
                'printf \'%.*f\\n\' "${dp}" "${value}"',
                "builtin printf '%s\\n' \"${value}\"",
                'echo "${value}"',
                "printf '%s\\n' \"${value}\" > value.txt",
                "cat <<'EOF'",
                "1.25",
                "EOF",
                'eval "${command}"',
                'source "${dynamic_source}"',
                "awk '{ print $1 }' values.tsv",
                '"${command}" "${value}"',
                "",
            ],
        ),
        encoding="utf-8",
    )
    return inventory(tmp_path, {"python", "shell"})


def test_supported_sink_classes_are_enumerated_with_bounded_claim(
    tmp_path: Path,
) -> None:
    result = fixture_inventory(tmp_path)
    classes = {item["syntactic_class"] for item in result["sites"]}
    assert {
        "python.direct_print",
        "python.sys_stream_write",
        "python.bound_text_stream_write",
        "python.csv_writer_row",
        "python.json_dump",
        "python.path_write_text",
        "shell.direct_printf_or_builtin_printf",
        "shell.direct_echo",
        "shell.static_redirection",
        "shell.literal_heredoc_or_herestring",
    } <= classes
    assert result["inventory_scope"]["claim"] == (
        "complete_only_within_recognized_static_syntactic_classes"
    )
    assert {item["applicability"] for item in result["sites"]} == {
        "unresolved"
    }
    assert {item["exception"] for item in result["sites"]} == {"unresolved"}


def test_unsupported_constructs_are_explicit_review_candidates(
    tmp_path: Path,
) -> None:
    result = fixture_inventory(tmp_path)
    reasons = {
        item["unsupported_reason"] for item in result["review_candidates"]
    }
    assert {
        "dynamic_or_indirect_call_target",
        "eval_or_eval_like_construct",
        "sourced_code_outside_scanned_statement",
        "embedded_language",
        "function_or_command_variable_dispatch",
    } <= reasons


def test_numeric_formatting_is_evidence_only_at_the_same_sink(
    tmp_path: Path,
) -> None:
    result = fixture_inventory(tmp_path)
    formatted = [
        item for item in result["sites"] if item["formatting_evidence"]
    ]
    assert formatted
    assert all(item["review_status"] == "evidence_only" for item in formatted)
    assert all(item["dp_source"] == "unresolved" for item in formatted)


def test_inventory_validates_against_the_supplied_schema(
    tmp_path: Path,
) -> None:
    result = fixture_inventory(tmp_path)
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    assert list(Draft202012Validator(schema).iter_errors(result)) == []
