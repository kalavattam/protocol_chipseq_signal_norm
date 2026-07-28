#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_help_aliases.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""
Focused regressions for authoritative shell Parameter alias comparison.
"""

from __future__ import annotations

import unittest
from pathlib import Path

from dev.audit import help_aliases
from dev.audit.help_aliases import (
    AliasChunk,
    ParameterRow,
    alias_chunks,
    check_document,
    csv_short_alias_findings,
    expected_aliases,
    file_alias_chunks,
    normalize_document,
    parameter_rows,
    scan_repository,
    usage_aliases_by_owner,
    wrapper_documents,
)

REPO_ROOT = Path(__file__).resolve().parents[3]


class HelpAliasesTest(unittest.TestCase):
    """
    Prove bounded logical grouping, visibility, and normalization.
    """

    def test_combined_case_arm_splits_at_each_logical_short_alias(
        self,
    ) -> None:
        text = """
function parse_args() {
    case "${1}" in
        -dr|--dry|--dry[_-]run|-y|--yes)
            shift
            ;;
    esac
}
"""

        self.assertEqual(
            [chunk.aliases for chunk in alias_chunks(text)],
            [
                ("-dr", "--dry", "--dry_run", "--dry-run"),
                ("-y", "--yes"),
            ],
        )

    def test_literal_parser_acceptance_controls_public_visibility(
        self,
    ) -> None:
        row = ParameterRow("parse_args", 1, ("--dry_run",), "  ", "flag")
        chunk = AliasChunk(
            "parse_args",
            ("-dr", "--dry", "--dry_run", "--dry-run"),
        )
        public, hidden = expected_aliases(row, chunk)

        self.assertEqual(public, ("-dr", "--dry", "--dry_run"))
        self.assertEqual(hidden, ("--dry-run",))

    def test_chunk_abbreviations_are_hidden_compatibility_aliases(
        self,
    ) -> None:
        row = ParameterRow(
            "parse_args",
            1,
            ("-ck", "--chunk_size"),
            "  ",
            "int",
        )
        chunk = AliasChunk(
            "parse_args",
            (
                "-ck",
                "--chunk_size",
                "--chunk-size",
                "--chnk_size",
                "--chnk-size",
            ),
        )
        public, hidden = expected_aliases(row, chunk)

        self.assertEqual(public, ("-ck", "--chunk_size"))
        self.assertEqual(
            hidden,
            ("--chunk-size", "--chnk_size", "--chnk-size"),
        )

    def test_install_channel_inventory_separates_hidden_and_retired_aliases(
        self,
    ) -> None:
        """
        Installer inventory records canonical, hidden, and retired forms.
        """

        _, inventory = scan_repository(
            REPO_ROOT,
            [
                "install/scripts/install_envs_entrypoint.sh",
                "lib/bash/help/help_install_envs.sh",
            ],
        )
        selected_options = {"--channels", "--override_channels"}
        selected = {}

        for row in inventory:
            if row["logical_option"] not in selected_options:
                continue

            selected[(row["path"], row["logical_option"])] = row

        self.assertEqual(len(selected), 4)

        for path in (
            "install/scripts/install_envs_entrypoint.sh",
            "lib/bash/help/help_install_envs.sh",
        ):
            channels = selected[(path, "--channels")]

            self.assertEqual(
                channels["public_aliases"],
                ["-ch", "--channels"],
            )
            self.assertEqual(channels["hidden_aliases"], [])
            self.assertEqual(
                channels["rejected_or_retired_aliases"],
                ["--channel", "--channel-list", "--channel_list"],
            )

            override = selected[(path, "--override_channels")]

            self.assertEqual(
                override["public_aliases"],
                ["-oc", "--override_channels"],
            )
            self.assertEqual(
                override["hidden_aliases"],
                ["--override-channels"],
            )
            self.assertEqual(
                override["rejected_or_retired_aliases"],
                ["--override-channel", "--override_channel"],
            )

    def test_compression_inventory_exposes_only_dry_run(self) -> None:
        """
        Compression inventory records one public dry-run spelling.
        """

        _, inventory = scan_repository(
            REPO_ROOT,
            ["lib/bash/help/help_compress_remove_files.sh"],
        )
        row = next(
            row for row in inventory if row["logical_option"] == "--dry_run"
        )

        self.assertEqual(row["public_aliases"], ["-dr", "--dry_run"])
        self.assertEqual(row["hidden_aliases"], [])
        self.assertEqual(row["restored_aliases"], ["-dr"])
        self.assertEqual(
            row["rejected_or_retired_aliases"],
            [
                "--chk-exc",
                "--chk-exu",
                "--chk_exc",
                "--chk_exu",
                "--dry",
                "--dry-run",
                "-ce",
                "-cu",
            ],
        )

    def test_help_short_alias_is_public_without_historical_help_evidence(
        self,
    ) -> None:
        row = ParameterRow("parse_args", 1, ("--help",), "  ", "flag")
        public, hidden = expected_aliases(
            row,
            AliasChunk("parse_args", ("-h", "--hlp", "--help")),
        )

        self.assertEqual(public, ("-h", "--help"))
        self.assertEqual(hidden, ("--hlp",))

    def test_regex_help_parser_exposes_public_h_and_hides_hlp(self) -> None:
        text = """
function main() {
    if [[ "${1:-}" =~ ^(-h|--h[e]?lp)$ ]]; then
        return 0
    fi
}
"""

        self.assertEqual(
            [chunk.aliases for chunk in alias_chunks(text)],
            [("-h", "--hlp", "--help")],
        )

    def test_runtime_requirement_bullets_are_not_parser_aliases(self) -> None:
        text = """
function demo() {
    local show_help
    show_help=$(cat << EOM
Usage
-----
  demo

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - gzip or zcat (when processing a gzipped input)
EOM
    )
    case "${1:-}" in
        -i|--input)
            shift
            ;;
    esac
}
"""

        self.assertEqual(
            [chunk.aliases for chunk in alias_chunks(text)],
            [("-i", "--input")],
        )

    def test_accepted_help_aliases_require_a_parameters_row(self) -> None:
        text = """
function demo() {
    local show_help
    show_help=$(cat << EOM
Usage
-----
  demo
    [--help] value

Parameters
----------
  1  value : str
    Value to check.

Returns
-------
  Returns 0.
EOM
    )
    if [[ "${1:-}" =~ ^(-h|--h[e]?lp)$ ]]; then
        return 0
    fi
}
"""

        findings, _ = check_document(REPO_ROOT, "lib/bash/core/demo.sh", text)

        self.assertEqual(len(findings), 1)
        self.assertEqual(findings[0].expected, ("-h", "--help"))

    def test_regex_details_and_all_help_aliases_are_complete(self) -> None:
        text = """
function main() {
    if [[ "${1:-}" =~ ^(-d|--details)$ ]]; then
        return 0
    elif [[ "${1:-}" =~ ^(-ah|--all[_-]h[e]?lp)$ ]]; then
        return 0
    fi
}
"""
        chunks = alias_chunks(text)

        self.assertEqual(chunks[0].aliases, ("-d", "--details"))

        public, hidden = expected_aliases(
            ParameterRow("main", 1, ("--all_help",), "  ", "flag"),
            chunks[1],
        )

        self.assertEqual(public, ("-ah", "--all_help"))
        self.assertEqual(hidden, ("--all_hlp", "--all-hlp", "--all-help"))

    def test_documented_d_and_ah_short_aliases_are_public_when_accepted(
        self,
    ) -> None:
        row = ParameterRow("parse_args", 1, ("--debug",), "  ", "flag")
        public, _ = expected_aliases(
            row,
            AliasChunk("parse_args", ("-d", "--debug", "-ah", "--all_help")),
        )

        self.assertEqual(public, ("-d", "--debug", "-ah", "--all_help"))

    def test_function_body_keeps_parser_after_nested_brace_group(self) -> None:
        text = """
function parse_args() {
    require_optarg || {
        return 1
    }
    case "${1}" in
        -h|--help)
            return 0
            ;;
    esac
}
"""

        self.assertEqual(
            [chunk.aliases for chunk in alias_chunks(text)],
            [("-h", "--help")],
        )

    def test_top_level_parser_alias_fallback_uses_central_extractor(
        self,
    ) -> None:
        text = """
while [[ "$#" -gt 0 ]]; do
    case "${1}" in
        -dr|--dry|--dry[_-]run)
            shift
            ;;
        -o|--out|--outputs)
            shift
            ;;
    esac
done
"""

        self.assertEqual(
            [chunk.aliases for chunk in file_alias_chunks(text)],
            [
                ("-dr", "--dry", "--dry_run", "--dry-run"),
                ("-o", "--out", "--outputs"),
            ],
        )

    def test_wrapper_scope_excludes_function_owner_baseline(self) -> None:
        self.assertEqual(
            wrapper_documents(
                [
                    "lib/bash/help/help_execute_trim_fastqs.sh",
                    "lib/bash/workflows/align_fastqs.sh",
                    "lib/bash/help/help_find_files.sh",
                ],
            ),
            ["lib/bash/help/help_execute_trim_fastqs.sh"],
        )

    def test_positional_row_is_not_an_option_alias_row(self) -> None:
        text = """
function show_help() {
    cat << EOM
Parameters
----------
  1  fil_in : file
    Input file.
EOM
}
"""

        self.assertEqual(parameter_rows(text), [])

    def test_usage_inventory_records_only_observed_spellings(self) -> None:
        text = """
function demo() {
    cat << EOM
Usage
-----
  demo
    [--help] [--dry_run]
EOM
}
"""

        self.assertEqual(
            usage_aliases_by_owner(text)["demo"],
            ("--help", "--dry_run"),
        )

    def test_usage_may_omit_an_option_but_cannot_use_its_short_alias(
        self,
    ) -> None:
        omitted = """
function demo() {
    cat << EOM
Usage
-----
  demo value

Parameters
----------
  -h, --help : flag
    Display help.
EOM
    case "${1:-}" in
        -h|--help)
            return 0
            ;;
    esac
}
"""

        findings, _ = check_document(
            REPO_ROOT,
            "lib/bash/core/demo.sh",
            omitted,
        )

        self.assertEqual(findings, [])

        short = omitted.replace("demo value", "demo [-h] value")

        findings, _ = check_document(REPO_ROOT, "lib/bash/core/demo.sh", short)

        self.assertEqual(len(findings), 1)
        self.assertIn("Usage requires", findings[0].message)

    def test_shared_usage_variable_is_attributed_to_help_owners(self) -> None:
        text = """
usage=$(cat << 'USAGE'
Usage
-----
  demo [--help]
USAGE
)
function demo() {
    cat << EOM
${usage}

Parameters
----------
  -h, --help : flag
    Display help.
EOM
}
"""

        self.assertEqual(usage_aliases_by_owner(text)["demo"], ("--help",))

    def test_normalization_adds_aliases_on_one_row(self) -> None:
        text = """
Parameters
----------
  --dry_run : flag
""".removeprefix("\n")
        normalized = normalize_document(
            text,
            {3: ("-dr", "--dry", "--dry_run")},
        )

        self.assertIn("  -dr, --dry, --dry_run : flag", normalized)

    def test_normalization_is_idempotent(self) -> None:
        text = """
Parameters
----------
  -dr, --dry, --dry_run : flag
"""
        aliases = {3: ("-dr", "--dry", "--dry_run")}

        self.assertEqual(
            normalize_document(normalize_document(text, aliases), aliases),
            normalize_document(text, aliases),
        )

    def delegated_filter_inventory(self) -> list[dict[str, object]]:
        """
        Return the source-derived filter delegations under test.
        """

        inventory = getattr(help_aliases, "delegated_parser_inventory", None)

        self.assertIsNotNone(
            inventory,
            "authoritative delegated-parser inventory API is missing",
        )
        assert inventory is not None

        filter_owners = {"filter_alignment_sc", "filter_alignment_sp"}
        rows = inventory(
            REPO_ROOT,
            ["lib/bash/workflows/filter_alignment.sh"],
        )

        return [
            row for row in rows if row["documented_owner"] in filter_owners
        ]

    def delegated_filter_row(self, owner: str) -> dict[str, object]:
        """
        Return one exact public filter-owner delegation.
        """

        rows = [
            row
            for row in self.delegated_filter_inventory()
            if row["documented_owner"] == owner
            and row["parser_owner"] == "_parse_args_filter_alignment"
        ]

        self.assertEqual(len(rows), 1)

        return rows[0]

    def test_public_owner_inherits_directly_delegated_parser_contract(
        self,
    ) -> None:
        source = """
function shared_parser() {
    case "${1}" in
        -t|--threads)
            shift
            ;;
    esac
}
function public_owner() {
    local show_help
    show_help=$(cat << EOM
Usage
-----
  public_owner --threads <int>

Parameters
----------
  -t, --threads : int
    Number of threads.
EOM
    )
    shared_parser "$@"
}
"""

        bindings = getattr(help_aliases, "delegated_parser_bindings", None)

        self.assertIsNotNone(
            bindings,
            "delegated binding discovery API is missing",
        )
        assert bindings is not None

        observed = bindings(
            REPO_ROOT,
            "lib/bash/core/demo.sh",
            source,
        )

        self.assertEqual(len(observed), 1)
        self.assertEqual(observed[0].documented_owner, "public_owner")
        self.assertEqual(observed[0].parser_owner, "shared_parser")
        self.assertEqual(
            [chunk.aliases for chunk in observed[0].applicable_chunks],
            [("-t", "--threads")],
        )

    def test_omitted_delegated_public_short_alias_fails(self) -> None:
        text = """
function shared_parser() {
    case "${1}" in
        -t|--threads)
            ;;
    esac
}
function public_owner() {
    local show_help
    show_help=$(cat << EOM
Usage
-----
  public_owner --threads <int>

Parameters
----------
  --threads : int
    Number of threads.
EOM
    )
    shared_parser "$@"
}
"""

        findings, _ = check_document(
            REPO_ROOT,
            "lib/bash/core/demo.sh",
            text,
        )

        self.assertEqual(len(findings), 1)
        self.assertEqual(findings[0].owner, "public_owner")
        self.assertEqual(findings[0].expected, ("-t", "--threads"))
        self.assertIn("delegated", findings[0].message)

    def test_filter_alignment_sc_receives_only_shared_public_contract(
        self,
    ) -> None:
        row = self.delegated_filter_row("filter_alignment_sc")

        self.assertEqual(row["fixed_mode"], "sc")
        self.assertEqual(
            row["public_aliases"],
            [
                "-t",
                "--thr",
                "--threads",
                "-fi",
                "--fil_in",
                "-fo",
                "--fil_out",
                "-m",
                "--mit",
                "--mito",
                "-cc",
                "--chk_chr",
                "-rf",
                "--ref_fa",
            ],
        )
        self.assertEqual(
            row["rejected_aliases"],
            ["-tg", "--tg", "-mr", "--mr", "--mtr"],
        )

    def test_filter_alignment_sp_receives_shared_and_sp_public_contract(
        self,
    ) -> None:
        row = self.delegated_filter_row("filter_alignment_sp")

        self.assertEqual(row["fixed_mode"], "sp")
        self.assertEqual(
            row["public_aliases"],
            [
                "-t",
                "--thr",
                "--threads",
                "-fi",
                "--fil_in",
                "-fo",
                "--fil_out",
                "-m",
                "--mit",
                "--mito",
                "-cc",
                "--chk_chr",
                "-rf",
                "--ref_fa",
                "-tg",
                "--tg",
                "-mr",
                "--mr",
                "--mtr",
            ],
        )
        self.assertEqual(row["rejected_aliases"], [])

    def test_filter_alignment_sc_rejected_mode_alias_is_an_overclaim(
        self,
    ) -> None:
        text = (
            REPO_ROOT / "lib/bash/workflows/filter_alignment.sh"
        ).read_text()
        marker = """
  -cc, --chk_chr : flag
    Check chromosomes in output alignment files.
"""
        overclaim = """
  -tg, --tg : flag
    Retain SP_II_TG chromosome.

"""

        text = text.replace(marker, overclaim + marker, 1)

        findings, _ = check_document(
            REPO_ROOT,
            "lib/bash/workflows/filter_alignment.sh",
            text,
        )

        rejected = [
            finding
            for finding in findings
            if finding.owner == "filter_alignment_sc"
            and "rejected by delegated fixed-mode" in finding.message
        ]

        self.assertEqual(len(rejected), 1)
        self.assertEqual(rejected[0].documented, ("-tg", "--tg"))

    def test_delegated_hidden_pattern_aliases_are_not_public_rows(
        self,
    ) -> None:
        for owner in ("filter_alignment_sc", "filter_alignment_sp"):
            row = self.delegated_filter_row(owner)

            self.assertEqual(row["hidden_aliases"], ["--chk-chr", "--ref-fa"])
            self.assertNotIn("--chk-chr", row["documented_aliases"])
            self.assertNotIn("--ref-fa", row["documented_aliases"])

    def test_unsupported_and_transposed_filter_spellings_stay_rejected(
        self,
    ) -> None:
        for owner in ("filter_alignment_sc", "filter_alignment_sp"):
            row = self.delegated_filter_row(owner)
            accepted = set(row["accepted_aliases"])

            for rejected in ("--thread", "--fil-in", "--ref_af", "--chr_chk"):
                self.assertNotIn(rejected, accepted)
                self.assertNotIn(rejected, row["documented_aliases"])

    def test_delegated_usage_remains_canonical_long_only(self) -> None:
        text = (
            REPO_ROOT / "lib/bash/workflows/filter_alignment.sh"
        ).read_text()
        usage = usage_aliases_by_owner(text)
        expected = {
            "filter_alignment_sc": (
                "--help",
                "--threads",
                "--fil_in",
                "--fil_out",
                "--ref_fa",
                "--mito",
                "--chk_chr",
            ),
            "filter_alignment_sp": (
                "--help",
                "--threads",
                "--fil_in",
                "--fil_out",
                "--ref_fa",
                "--mito",
                "--tg",
                "--mtr",
                "--chk_chr",
            ),
        }

        for owner, aliases in expected.items():
            self.assertEqual(usage[owner], aliases)
            self.assertTrue(
                all(alias.startswith("--") for alias in usage[owner]),
            )

    def test_internal_filter_parser_documents_options_after_delimiter(
        self,
    ) -> None:
        text = (
            REPO_ROOT / "lib/bash/workflows/filter_alignment.sh"
        ).read_text()
        rows = [
            row.aliases
            for row in parameter_rows(text)
            if row.owner == "_parse_args_filter_alignment"
        ]

        self.assertEqual(
            rows,
            [
                ("-h", "--help"),
                ("-t", "--thr", "--threads"),
                ("-fi", "--fil_in"),
                ("-fo", "--fil_out"),
                ("-m", "--mit", "--mito"),
                ("-cc", "--chk_chr"),
                ("-rf", "--ref_fa"),
                ("-tg", "--tg"),
                ("-mr", "--mr", "--mtr"),
            ],
        )

        body = help_aliases.function_bodies(text)[
            "_parse_args_filter_alignment"
        ]

        self.assertIn('if [[ "${1:-}" == "--" ]]', body)

        for owner in ("filter_alignment_sc", "filter_alignment_sp"):
            calls = help_aliases.parser_calls(
                text,
                owner,
                "_parse_args_filter_alignment",
            )

            self.assertEqual(len(calls), 1)
            self.assertIn("--", calls[0][1])

    def test_help_heredoc_masking_survives_delegated_discovery(self) -> None:
        source = """
function shared_parser() {
    local show_help
    show_help=$(cat << EOM
Usage
-----
  shared_parser [--fake]

Parameters
----------
  --fake : flag
    Text that must not become parser evidence.
EOM
    )
    case "${1}" in
        -t|--threads)
            ;;
    esac
}
function public_owner() {
    local show_help
    show_help=$(cat << EOM
Usage
-----
  public_owner --threads <int>

Parameters
----------
  -t, --threads : int
    Number of threads.
EOM
    )
    shared_parser "$@"
}
"""

        self.assertEqual(
            [chunk.aliases for chunk in alias_chunks(source)],
            [("-t", "--threads")],
        )

        bindings = getattr(help_aliases, "delegated_parser_bindings", None)

        self.assertIsNotNone(
            bindings,
            "delegated binding discovery API is missing",
        )
        assert bindings is not None

        observed = bindings(
            REPO_ROOT,
            "lib/bash/core/demo.sh",
            source,
        )

        self.assertEqual(
            [
                chunk.aliases
                for binding in observed
                for chunk in binding.applicable_chunks
            ],
            [("-t", "--threads")],
        )

    def test_delegated_inventory_is_deterministic(self) -> None:
        inventory = getattr(help_aliases, "delegated_parser_inventory", None)

        self.assertIsNotNone(inventory, "delegated inventory API is missing")
        assert inventory is not None

        first = inventory(REPO_ROOT)
        second = inventory(REPO_ROOT)

        self.assertEqual(first, second)
        self.assertEqual(
            first,
            sorted(
                first,
                key=lambda row: (
                    str(row["documentation_path"]),
                    str(row["documented_owner"]),
                    str(row["parser_path"]),
                    str(row["parser_owner"]),
                ),
            ),
        )

    def test_directly_owned_parser_help_relationship_is_unchanged(
        self,
    ) -> None:
        text = """
function demo() {
    cat << EOM
Usage
-----
  demo --threads <int>

Parameters
----------
  -t, --threads : int
    Number of threads.
EOM
    case "${1}" in
        -t|--threads) ;;
    esac
}
"""

        bindings = getattr(help_aliases, "delegated_parser_bindings", None)

        self.assertIsNotNone(
            bindings,
            "delegated binding discovery API is missing",
        )
        assert bindings is not None
        self.assertEqual(
            bindings(REPO_ROOT, "lib/bash/core/demo.sh", text),
            [],
        )

        findings, _ = check_document(
            REPO_ROOT,
            "lib/bash/core/demo.sh",
            text,
        )

        self.assertEqual(findings, [])

    def test_all_current_help_owners_have_zero_unexplained_alias_gaps(
        self,
    ) -> None:
        paths = sorted(
            str(path.relative_to(REPO_ROOT))
            for base in (
                REPO_ROOT / "scripts",
                REPO_ROOT / "install" / "scripts",
                REPO_ROOT / "tests",
            )
            for path in base.rglob("*.sh")
            if "outputs" not in path.relative_to(REPO_ROOT).parts
        )

        findings, _ = scan_repository(REPO_ROOT, paths)

        self.assertEqual(findings, [])

    def test_public_csv_short_aliases_must_begin_with_c(self) -> None:
        rows = [
            {
                "path": "demo.sh",
                "owner": "parse_args",
                "logical_option": "--csv_value",
                "public_aliases": ["-v", "--csv_value"],
            },
        ]

        findings = csv_short_alias_findings(rows)

        self.assertEqual(len(findings), 1)
        self.assertIn("-v", findings[0].message)

    def test_csv_option_without_short_alias_does_not_acquire_one(self) -> None:
        rows = [
            {
                "path": "demo.sh",
                "owner": "parse_args",
                "logical_option": "--csv_value",
                "public_aliases": ["--csv_value"],
            },
        ]

        self.assertEqual(csv_short_alias_findings(rows), [])

    def test_submit_compute_signal_csv_fragment_aliases_are_settled(
        self,
    ) -> None:
        source = (REPO_ROOT / "bin/submit_compute_signal.sh").read_text()
        help_text = (
            REPO_ROOT / "lib/bash/help/help_submit_compute_signal.sh"
        ).read_text()
        chunks = [
            chunk.aliases
            for chunk in alias_chunks(source, {"parse_args"})
            if "--csv_usr_frg" in chunk.aliases
        ]

        self.assertEqual(
            chunks,
            [("-cuf", "--csv_usr_frg", "--csv-usr-frg")],
        )

        rows = [
            row.aliases
            for row in parameter_rows(help_text)
            if "--csv_usr_frg" in row.aliases
        ]

        self.assertEqual(rows, [("-cuf", "--csv_usr_frg")])

        usage = usage_aliases_by_owner(help_text)["help_submit_compute_signal"]

        self.assertIn("--csv_usr_frg", usage)
        self.assertNotIn("-cuf", usage)

        parser_body = help_aliases.function_bodies(source)["parse_args"]

        for retired in ("-uf", "--usr_frg", "--usr-frg"):
            self.assertNotIn(retired, parser_body)
            self.assertNotIn(retired, help_text)


if __name__ == "__main__":
    unittest.main()
