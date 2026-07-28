#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: source_policy.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""
Check ShellCheck source directives and submit bootstrap policy.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import re
import sys
from pathlib import Path, PurePosixPath

from dev.audit.help_heredoc_reflow import shell_paths

SOURCE_LINE = re.compile(
    r"^\s*(?:source|\.)\s+(?P<expression>.+?)\s*(?:\|\||$)",
)
DISABLE = re.compile(
    r"^\s*#\s*shellcheck\s+disable=(?P<codes>SC\d+(?:,SC\d+)*)\s*$",
)
DIRECTIVE = re.compile(r"^\s*#\s*shellcheck\s+source=(?P<path>\S+)\s*$")
HELP_BOOTSTRAP_PATHS = frozenset(
    {
        "bin/submit_align_fastqs.sh",
        "bin/submit_calculate_scaling_factor.sh",
        "bin/submit_compute_signal.sh",
        "bin/submit_convert_bam_bed.sh",
        "bin/submit_download_fastqs.sh",
        "bin/submit_filter_alignments.sh",
        "bin/submit_qsort_bam.sh",
        "bin/submit_trim_fastqs.sh",
    },
)
CANONICAL_BASH_SHEBANG = "#!/usr/bin/env bash"
SUBMIT_GLOB = "submit_*.sh"
SOURCE_CALL = re.compile(
    r"^\s*(?:source|\.)\s+[^#\n]*submit_[A-Za-z0-9_]+\.sh",
)
COMPAT_EXEC = re.compile(
    r'^exec "\$\{BASH\}" "\$\{dir_scr\}/'
    r'(?P<target>submit_[A-Za-z0-9_]+\.sh)" "\$@"$',
)
SUBMIT_MODE_MESSAGE = (
    "submit interfaces are interpreter/Slurm entry points and must remain "
    "non-executable"
)
COMPATIBILITY_BOOTSTRAP_MESSAGE = (
    "compatibility wrapper must delegate before helper or environment "
    "bootstrap"
)
GUARD_FRAGMENT_PREFIX = (
    "Bash guard missing canonical boundary/diagnostic fragment:"
)
DIRECT_CONTRACT_MESSAGE = (
    "worker-only interface must invoke main directly and is not sourceable"
)


@dataclasses.dataclass(frozen=True)
class SourceRecord:
    """
    One source statement's bounded resolution policy.
    """

    classification: str
    target: str | None
    required_suppression: str | None


@dataclasses.dataclass(frozen=True)
class Advisory:
    """
    One non-failing source-policy diagnostic.
    """

    path: str
    line: int
    message: str
    rule_id: str = "SHELL.SOURCE.DIRECTIVES"

    def format(self) -> str:
        """
        Render a stable source-policy diagnostic.
        """

        return f"{self.rule_id}: {self.path}:{self.line}: {self.message}"


@dataclasses.dataclass(frozen=True)
class SubmitInterface:
    """
    One source-derived submit or compatibility interface.
    """

    path: str
    stable_interface_identity: str
    current_shebang: str
    executable_bit: bool
    direct_execution_callers: tuple[str, ...]
    source_callers: tuple[str, ...]
    compatibility_delegators: tuple[str, ...]
    bash_guard_location: int | None
    guard_before_bash_4_4_syntax: bool
    helper_bootstrap_mechanism: str
    environment_bootstrap_mechanism: str
    slurm_local_parallel_behavior: str
    current_tests: tuple[str, ...]
    proposed_policy_disposition: str
    classification: str

    def as_dict(self) -> dict[str, object]:
        """
        Return a stable JSON-ready inventory row.
        """

        return dataclasses.asdict(self)


def _maintained_reference_paths(root: Path) -> list[Path]:
    """
    Return maintained text sources that may name submit interfaces.
    """

    paths = [
        *(root / "bin").glob("*.sh"),
        *(root / "lib/bash").rglob("*.sh"),
        *(root / "install/scripts").glob("*.sh"),
        *(root / "tests").rglob("*.sh"),
    ]

    return sorted(
        {path for path in paths if "artifacts/tests" not in path.parts},
    )


def _reference_inventory(
    root: Path,
    submit_path: Path,
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    """
    Return direct, source, compatibility, and test references.
    """

    name = submit_path.name
    direct: set[str] = set()
    sourced: set[str] = set()
    delegated: set[str] = set()
    tests: set[str] = set()

    for path in _maintained_reference_paths(root):
        if path == submit_path:
            continue

        text = path.read_text(encoding="utf-8")
        if name not in text:
            continue

        relative = path.relative_to(root).as_posix()
        matching = [line for line in text.splitlines() if name in line]

        if any(SOURCE_CALL.search(line) for line in matching):
            sourced.add(relative)
        elif any("exec " in line and name in line for line in matching):
            delegated.add(relative)
        elif relative.startswith("bin/execute_") and any(
            "scr_sub=" in line for line in matching
        ):
            direct.add(relative)

        if relative.startswith("tests/"):
            tests.add(relative)

    bootstrap_test = "tests/contract/interfaces/test_submit_bootstrap.sh"

    if (root / bootstrap_test).is_file():
        tests.add(bootstrap_test)

    return (
        tuple(sorted(direct)),
        tuple(sorted(sourced)),
        tuple(sorted(delegated)),
        tuple(sorted(tests)),
    )


def _first_code_line(lines: list[str]) -> int | None:
    """
    Return the first non-header executable source line index.
    """

    for index, line in enumerate(lines[1:], 1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        return index

    return None


def _guard_end(lines: list[str], start: int) -> int | None:
    """
    Return the canonical top-level guard's closing 'fi' index.
    """

    for index in range(start, len(lines)):
        if lines[index].strip() == "fi":
            return index

    return None


def discover_submit_interfaces(root: Path) -> list[SubmitInterface]:
    """
    Build the complete maintained submit-interface inventory.

    Parameters
    ----------
    root : Path
        Repository root containing submit wrappers and their callers.

    Returns
    -------
    interfaces : list[SubmitInterface]
        Source-derived startup, delegation, behavior, and test ownership rows.
    """

    root = root.resolve()
    rows: list[SubmitInterface] = []

    for path in sorted((root / "bin").glob(SUBMIT_GLOB)):
        text = path.read_text(encoding="utf-8")
        lines = text.splitlines()
        relative = path.relative_to(root).as_posix()
        direct, sourced, delegated, tests = _reference_inventory(root, path)
        exec_match = next(
            (
                COMPAT_EXEC.match(line.strip())
                for line in lines
                if COMPAT_EXEC.match(line.strip())
            ),
            None,
        )
        guard = next(
            (
                index
                for index, line in enumerate(lines)
                if line.strip().startswith('if [[ -z "${BASH_VERSION:-}"')
            ),
            None,
        )
        first_code = _first_code_line(lines)
        guard_before = exec_match is not None or (
            guard is not None
            and first_code is not None
            and first_code <= guard
        )
        execute_owner = (
            root / "bin" / path.name.replace("submit_", "execute_", 1)
        )

        if exec_match is not None:
            classification = "compatibility delegator"
            behavior = "exec delegation to the canonical guarded worker"
            helper = "none; canonical target owns helper bootstrap"
            environment = "none; canonical target owns environment bootstrap"
            disposition = (
                "retain as a non-executable compatibility delegator with the "
                "canonical Bash shebang and exact argument/status-preserving "
                "exec"
            )
        else:
            classification = "worker-only interface"
            modes: list[str] = []
            owner_text = (
                execute_owner.read_text(encoding="utf-8")
                if execute_owner.is_file()
                else text
            )

            if "sbatch" in owner_text or "SLURM_ARRAY_TASK_ID" in text:
                modes.append("Slurm")

            if "parallel" in owner_text:
                modes.append("GNU Parallel")

            if (
                '"${BASH}"' in owner_text
                or "run_local" in text
                or "run_job" in text
            ):
                modes.append("local")

            behavior = "/".join(modes) + " worker dispatch"
            helper = (
                "bootstrap help source, then "
                "source_helpers_submit/source_helpers"
            )
            environment = (
                "setup_env/handle_env after guard and helper bootstrap"
                if "setup_env" in text or "handle_env" in text
                else "caller-provided environment; no activation in worker"
            )
            disposition = (
                "retain as a non-executable worker entry point with the "
                "canonical Bash shebang, an immediate Bash >= "
                "4.4 guard, and explicit-interpreter/Slurm dispatch"
            )

        rows.append(
            SubmitInterface(
                path=relative,
                stable_interface_identity=path.stem,
                current_shebang=lines[0] if lines else "",
                executable_bit=os.access(path, os.X_OK),
                direct_execution_callers=direct,
                source_callers=sourced,
                compatibility_delegators=delegated,
                bash_guard_location=guard + 1 if guard is not None else None,
                guard_before_bash_4_4_syntax=guard_before,
                helper_bootstrap_mechanism=helper,
                environment_bootstrap_mechanism=environment,
                slurm_local_parallel_behavior=behavior,
                current_tests=tests,
                proposed_policy_disposition=disposition,
                classification=classification,
            ),
        )

    return rows


def _append_submit_interface_findings(
    row: SubmitInterface,
    findings: list[Advisory],
) -> None:
    """
    Record findings shared by worker and compatibility submit interfaces.

    Parameters
    ----------
    row : SubmitInterface
        Source-derived interface contract.
    findings : list[Advisory]
        Mutable deterministic finding collection.
    """

    if row.current_shebang != CANONICAL_BASH_SHEBANG:
        findings.append(
            Advisory(
                row.path,
                1,
                f"submit entry point requires '{CANONICAL_BASH_SHEBANG}'",
                "SHELL.SUBMIT.SHEBANG",
            ),
        )

    if row.executable_bit:
        findings.append(
            Advisory(
                row.path,
                1,
                SUBMIT_MODE_MESSAGE,
                "SHELL.SUBMIT.MODE",
            ),
        )

    if row.source_callers:
        findings.append(
            Advisory(
                row.path,
                1,
                ("worker/compatibility submit interface has source callers: ")
                + ", ".join(row.source_callers),
                "SHELL.SUBMIT.SOURCE_CONTRACT",
            ),
        )


def _append_compatibility_submit_findings(
    row: SubmitInterface,
    lines: list[str],
    exec_match: re.Match[str] | None,
    by_name: dict[str, SubmitInterface],
    findings: list[Advisory],
) -> None:
    """
    Validate one compatibility submit interface.

    Parameters
    ----------
    row : SubmitInterface
        Compatibility interface contract.
    lines : list[str]
        Physical source lines.
    exec_match : re.Match[str] | None
        Recognized compatibility delegation, if present.
    by_name : dict[str, SubmitInterface]
        Maintained submit interfaces keyed by basename.
    findings : list[Advisory]
        Mutable deterministic finding collection.
    """

    if exec_match is None or exec_match.group("target") not in by_name:
        findings.append(
            Advisory(
                row.path,
                len(lines),
                (
                    "compatibility wrapper must exec a maintained submit "
                    "target via '${BASH}' and exact '$@'"
                ),
                "SHELL.SUBMIT.COMPATIBILITY",
            ),
        )

    forbidden = re.compile(
        r"\b(?:source|conda|handle_env|setup_env)\b",
    )

    for index, line in enumerate(lines, 1):
        bootstraps = forbidden.search(line)
        is_comment = line.lstrip().startswith("#")

        if bootstraps and not is_comment:
            findings.append(
                Advisory(
                    row.path,
                    index,
                    COMPATIBILITY_BOOTSTRAP_MESSAGE,
                    "SHELL.SUBMIT.COMPATIBILITY",
                ),
            )


def _append_guard_contract_findings(
    row: SubmitInterface,
    lines: list[str],
    first_code: int,
    guard_end: int,
    findings: list[Advisory],
) -> None:
    """
    Validate the contents and placement of one worker's Bash guard.

    Parameters
    ----------
    row : SubmitInterface
        Worker interface contract.
    lines : list[str]
        Physical source lines.
    first_code : int
        Zero-based line containing the first executable statement.
    guard_end : int
        Zero-based line containing the guard's closing statement.
    findings : list[Advisory]
        Mutable deterministic finding collection.
    """

    guard_text = "\n".join(lines[first_code : guard_end + 1])
    required = (
        "BASH_VERSINFO[0] < 4",
        "BASH_VERSINFO[0] == 4",
        "BASH_VERSINFO[1] < 4",
        "requires Bash >= 4.4",
        "exit 1",
    )

    for fragment in required:
        if fragment.lower() not in guard_text.lower():
            findings.append(
                Advisory(
                    row.path,
                    first_code + 1,
                    f"{GUARD_FRAGMENT_PREFIX} {fragment}",
                    "SHELL.SUBMIT.GUARD_CONTRACT",
                ),
            )

    guard_index = row.bash_guard_location - 1
    pre_guard = "\n".join(lines[first_code:guard_index])
    bootstrap_before_guard = re.search(
        r"\b(?:source|conda|handle_env|setup_env|local -n|declare -gA)\b",
        pre_guard,
    )

    if bootstrap_before_guard:
        findings.append(
            Advisory(
                row.path,
                first_code + 1,
                (
                    "helper/environment bootstrap or Bash >= 4.4 syntax "
                    "precedes the version boundary"
                ),
                "SHELL.SUBMIT.GUARD_ORDER",
            ),
        )


def _append_worker_submit_findings(
    row: SubmitInterface,
    text: str,
    lines: list[str],
    first_code: int | None,
    findings: list[Advisory],
) -> None:
    """
    Validate one non-sourceable worker submit interface.

    Parameters
    ----------
    row : SubmitInterface
        Worker interface contract.
    text : str
        Complete worker source.
    lines : list[str]
        Physical source lines.
    first_code : int | None
        Zero-based first executable source line.
    findings : list[Advisory]
        Mutable deterministic finding collection.
    """

    guard_is_first = (
        row.bash_guard_location is not None
        and first_code is not None
        and lines[first_code].startswith("if [[")
    )

    if not guard_is_first:
        findings.append(
            Advisory(
                row.path,
                1,
                (
                    "worker must begin with the Bash >= 4.4 guard before "
                    "other source constructs"
                ),
                "SHELL.SUBMIT.GUARD_ORDER",
            ),
        )

        return

    guard_end = _guard_end(lines, first_code)

    if guard_end is None:
        findings.append(
            Advisory(
                row.path,
                first_code + 1,
                "Bash >= 4.4 guard has no top-level closing 'fi'",
                "SHELL.SUBMIT.GUARD_ORDER",
            ),
        )

        return

    _append_guard_contract_findings(
        row,
        lines,
        first_code,
        guard_end,
        findings,
    )

    invokes_main = re.search(r'^main "\$@"\s*$', text, re.MULTILINE)

    if not invokes_main:
        findings.append(
            Advisory(
                row.path,
                len(lines),
                DIRECT_CONTRACT_MESSAGE,
                "SHELL.SUBMIT.DIRECT_CONTRACT",
            ),
        )


def submit_bootstrap_findings(
    root: Path,
) -> tuple[list[Advisory], list[SubmitInterface]]:
    """
    Enforce the settled submit shebang and startup/bootstrap contract.

    Parameters
    ----------
    root : Path
        Repository root containing maintained submit interfaces.

    Returns
    -------
    findings, inventory : tuple[list[Advisory], list[SubmitInterface]]
        Deterministic startup findings and complete interface inventory.
    """

    root = root.resolve()
    findings: list[Advisory] = []

    inventory = discover_submit_interfaces(root)
    by_name = {Path(row.path).name: row for row in inventory}

    for row in inventory:
        path = root / row.path
        text = path.read_text(encoding="utf-8")
        lines = text.splitlines()
        first_code = _first_code_line(lines)
        exec_match = next(
            (
                COMPAT_EXEC.match(line.strip())
                for line in lines
                if COMPAT_EXEC.match(line.strip())
            ),
            None,
        )

        _append_submit_interface_findings(row, findings)

        if row.classification == "compatibility delegator":
            _append_compatibility_submit_findings(
                row,
                lines,
                exec_match,
                by_name,
                findings,
            )

            continue

        _append_worker_submit_findings(
            row,
            text,
            lines,
            first_code,
            findings,
        )

    return findings, inventory


def classify_source(path: str, line: str) -> SourceRecord:
    """
    Classify only source shapes used by the active repository tree.

    Parameters
    ----------
    path : str
        Repository-relative path containing the source statement.
    line : str
        Complete recognized Shell 'source' statement.

    Returns
    -------
    record : SourceRecord
        Static, dynamic, or runtime-selected source ownership.

    Raises
    ------
    ValueError
        If the line is not a recognized source statement.
    """

    expression = SOURCE_LINE.match(line)
    if expression is None:
        raise ValueError("line is not a recognized source statement")

    value = expression.group("expression")
    name = PurePosixPath(path).name

    if "conda info --base" in value or "source activate" in value:
        return SourceRecord("dynamic_activation", None, "SC1091")

    if re.fullmatch(r'"?\$\{?(?:fnc_src|path)\}?"?', value):
        return SourceRecord("runtime_selected", None, "SC1090")

    if "source_helpers.sh" in value:
        return SourceRecord("static", "lib/bash/core/source_helpers.sh", None)

    if (
        value.startswith('"$(')
        and path.startswith("tests/fixtures/")
        and name == "make.sh"
    ):
        return SourceRecord("static", "tests/support/fixture_helpers.sh", None)

    if value.startswith('"$(') and path.startswith("tests/"):
        return SourceRecord("static", "tests/support/test_helpers.sh", None)

    filter_workflow = path in {
        "lib/bash/workflows/filter_bam.sh",
        "lib/bash/workflows/filter_cram.sh",
    }

    if value.startswith('"$(') and filter_workflow:
        return SourceRecord(
            "static",
            "lib/bash/workflows/filter_alignment.sh",
            None,
        )

    if (
        value.startswith('"$(')
        and path == "lib/bash/help/help_execute_filter_bams.sh"
    ):
        return SourceRecord(
            "static",
            "lib/bash/help/help_execute_filter_alignments.sh",
            None,
        )

    if "fixture_helpers.sh" in value:
        return SourceRecord("static", "tests/support/fixture_helpers.sh", None)

    if "test_helpers.sh" in value:
        return SourceRecord("static", "tests/support/test_helpers.sh", None)

    if "filter_alignment.sh" in value:
        return SourceRecord(
            "static",
            "lib/bash/workflows/filter_alignment.sh",
            None,
        )

    if "help_submit_" in value:
        return SourceRecord(
            "static",
            f"lib/bash/help/help_{name}",
            None,
        )

    if name.startswith("submit_") and name.endswith(".sh"):
        return SourceRecord(
            "static",
            f"lib/bash/help/help_{name}",
            None,
        )

    return SourceRecord("unclassified", None, None)


def preceding_source_comments(
    lines: list[str],
    index: int,
) -> tuple[set[str], str | None]:
    """
    Return immediately attached disables and an exact source directive.
    """

    codes: set[str] = set()
    directive: str | None = None
    cursor = index - 1

    while cursor >= 0:
        line = lines[cursor]
        disabled = DISABLE.match(line)
        mapped = DIRECTIVE.match(line)

        if disabled:
            codes.update(disabled.group("codes").split(","))
        elif mapped:
            directive = mapped.group("path")
        else:
            break

        cursor -= 1

    return codes, directive


def suppression_mismatches(path: str, text: str) -> list[Advisory]:
    """
    Report stale source suppressions and incomplete static source maps.
    """

    findings: list[Advisory] = []
    lines = text.splitlines()

    for index, line in enumerate(lines):
        if SOURCE_LINE.match(line) is None:
            continue

        record = classify_source(path, line)
        disables, directive = preceding_source_comments(lines, index)

        if record.classification == "static":
            if directive != record.target:
                findings.append(
                    Advisory(
                        path,
                        index + 1,
                        "static source requires exact "
                        f"'shellcheck source={record.target}' before the "
                        f"statement",
                    ),
                )

            stale = sorted(
                code for code in disables if code in {"SC1090", "SC1091"}
            )

            if stale:
                findings.append(
                    Advisory(
                        path,
                        index + 1,
                        (
                            f"static source has stale {','.join(stale)} "
                            f"suppression; use source= mapping"
                        ),
                    ),
                )
        elif (
            record.required_suppression
            and record.required_suppression not in disables
        ):
            findings.append(
                Advisory(
                    path,
                    index + 1,
                    f"{record.classification} source requires exact "
                    f"{record.required_suppression} suppression",
                ),
            )

    return findings


def bootstrap_help_spacing_warnings(path: str, text: str) -> list[Advisory]:
    """
    Warn only for spacing before the eight submit-wrapper help sources.
    """

    if path not in HELP_BOOTSTRAP_PATHS:
        return []

    findings: list[Advisory] = []
    lines = text.splitlines()

    for index, line in enumerate(lines):
        if "source " not in line or "help_submit_" not in line:
            continue

        cursor = index - 1

        while cursor >= 0 and (
            DISABLE.match(lines[cursor]) or DIRECTIVE.match(lines[cursor])
        ):
            cursor -= 1

        blanks = 0

        while cursor >= 0 and not lines[cursor].strip():
            blanks += 1
            cursor -= 1

        if blanks != 1:
            findings.append(
                Advisory(
                    path,
                    index + 1,
                    "prefer exactly one blank line between the completed "
                    "bootstrap help-path assignment and its ShellCheck/source "
                    "statement",
                    "SHELL.HELP_SOURCE.SPACING",
                ),
            )

    return findings


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse the bounded source-policy checker command line.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    arguments : argparse.Namespace
        Parsed source-policy mode and output options.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--help-source-spacing", action="store_true")
    parser.add_argument("--submit-bootstrap", action="store_true")
    parser.add_argument("--submit-inventory-json", action="store_true")
    parser.add_argument("paths", nargs="*")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """
    Check current-diff source mapping or emit spacing advisories.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    status : int
        Zero when selected checks pass and one when findings remain.
    """

    args = parse_args(argv)
    root = args.root.resolve()

    if args.submit_bootstrap or args.submit_inventory_json:
        findings, inventory = submit_bootstrap_findings(root)

        if args.submit_inventory_json:
            print(
                json.dumps(
                    [row.as_dict() for row in inventory],
                    indent=2,
                    sort_keys=True,
                ),
            )
        else:
            for finding in findings:
                print(finding.format())

            counts: dict[str, int] = {}

            for row in inventory:
                counts[row.classification] = (
                    counts.get(row.classification, 0) + 1
                )

            print(f"SHELL.SUBMIT.BOOTSTRAP: interfaces={len(inventory)}")

            for classification in sorted(counts):
                print(
                    "SHELL.SUBMIT.BOOTSTRAP: "
                    f"{classification}={counts[classification]}",
                )

            print(f"SHELL.SUBMIT.BOOTSTRAP: findings={len(findings)}")

        return 1 if findings else 0

    paths = args.paths or shell_paths(root)
    warnings: list[Advisory] = []

    for path in paths:
        target = root / path
        if not target.is_file():
            continue

        text = target.read_text(encoding="utf-8")

        if args.help_source_spacing:
            warnings.extend(bootstrap_help_spacing_warnings(path, text))
        else:
            warnings.extend(suppression_mismatches(path, text))

    for warning in warnings:
        print(warning.format())

    if args.help_source_spacing:
        print(
            f"SHELL.HELP_SOURCE.SPACING: {len(warnings)} advisory warning(s)",
        )

        return 0

    if warnings:
        print(f"SHELL.SOURCE.DIRECTIVES: {len(warnings)} violation(s)")

        return 1

    print("SHELL.SOURCE.DIRECTIVES: pass (zero current-diff violations)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
