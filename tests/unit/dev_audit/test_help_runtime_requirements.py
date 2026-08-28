#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_help_runtime_requirements.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Focused regressions for Runtime-requirements extraction and policy.
"""

from __future__ import annotations

import json
import re
import tempfile
import unittest
from pathlib import Path

from dev.audit.help_heredoc_reflow import Heredoc
from dev.audit.help_runtime_requirements import (
    CHECKSUM_REQUIREMENT,
    CONDITION,
    DECOMPRESSION_REQUIREMENT,
    SETTLED_OWNER_REQUIREMENTS,
    Evidence,
    analyze_owner,
    candidate,
    evidence_for_owner,
    owner_units,
    rewrite_runtime_blocks,
    runtime_entries,
    runtime_policy_messages,
    settled_requirements,
)
from dev.audit.help_style import exact_callable_message
from dev.audit.shell_help_pilot import validate_command_registry
from dev.audit.shell_validation import MINIMUM_BASH


def document(body: str) -> Heredoc:
    return Heredoc(
        "demo",
        "EOM",
        1,
        len(body.splitlines()) + 1,
        tuple(enumerate(body.splitlines(), 1)),
    )


class RuntimeRequirementsTest(unittest.TestCase):
    def policy_rule_ids(self, *entries: str) -> set[str]:
        """
        Return policy rule IDs for one synthetic Runtime owner.
        """

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = Path(temp_dir) / "runtime_policy.sh"
            source.write_text(
                "function demo() {\n    :\n}\n",
                encoding="utf-8",
            )

            body = "Notes\n-----\n  Runtime requirements:\n" + "\n".join(
                f"    - {entry}" for entry in entries
            )
            identity = f"{source.relative_to(root)}::demo"

            _, findings, _ = analyze_owner(
                identity,
                ("function", source, source, document(body)),
                root,
                {
                    "bash",
                    "bowtie2",
                    "conda",
                    "parallel",
                    "python",
                    "python3",
                    "Rscript",
                    "samtools",
                },
                {
                    "Bowtie2": {"bowtie2"},
                    "Conda": {"conda"},
                    "Python": {"python"},
                    "SAMtools": {"samtools"},
                },
            )

        return {finding.rule_id for finding in findings}

    def test_singleton_is_unbulleted(self) -> None:
        label, entries, grammar = runtime_entries(
            document("Notes\n-----\n  Runtime requirements:\n    bash >= 4.4"),
        )

        self.assertEqual(label, 3)
        self.assertEqual(entries, [(4, "bash >= 4.4")])
        self.assertEqual(grammar, [])

    def test_categories_are_rejected(self) -> None:
        _, _, grammar = runtime_entries(
            document(
                (
                    "Notes\n-----\n  Runtime requirements:\n    - "
                    "Programs\n      + bash >= 4.4"
                ),
            ),
        )

        self.assertTrue(grammar)

    def test_runtime_requirements_outside_notes_are_rejected(self) -> None:
        _, _, grammar = runtime_entries(
            document(
                "Returns\n-------\n  Returns 0.\n\n"
                "  Runtime requirements:\n    bash >= 4.4",
            ),
        )

        self.assertIn(
            "Runtime requirements must occur inside the active top-level "
            "Notes section",
            grammar,
        )

    def test_runtime_requirements_duplicate_is_rejected(self) -> None:
        _, _, grammar = runtime_entries(
            document(
                "Notes\n-----\n  Runtime requirements:\n    bash >= 4.4\n\n"
                "  Runtime requirements:\n    samtools",
            ),
        )

        self.assertIn("Runtime requirements may appear only once", grammar)

    def test_runtime_rewrite_preserves_invalid_parent_and_is_idempotent(
        self,
    ) -> None:
        text = (
            "function demo() {\n"
            "    cat << 'EOM'\n"
            "Returns\n"
            "-------\n"
            "  Returns 0.\n"
            "\n"
            "  Runtime requirements:\n"
            "    samtools\n"
            "    bash >= 4.4\n"
            "EOM\n"
            "}\n"
        )
        once, changed = rewrite_runtime_blocks(text)
        twice, changed_twice = rewrite_runtime_blocks(once)

        self.assertEqual(changed, 1)
        self.assertEqual(changed_twice, 0)
        self.assertEqual(once, twice)
        self.assertIn("Returns\n-------", once)
        self.assertNotIn("Notes\n-----", once)

    def test_bash_boundary_policy_is_4_4(self) -> None:
        self.assertEqual(MINIMUM_BASH, (4, 4))
        self.assertFalse(MINIMUM_BASH <= (4, 3))
        self.assertTrue(MINIMUM_BASH <= (4, 4))
        self.assertTrue(MINIMUM_BASH <= (4, 99))
        self.assertTrue(MINIMUM_BASH <= (5, 0))

    def test_semantic_candidate_signature_is_stable(self) -> None:
        evidence = Evidence(
            "demo",
            "demo.sh",
            4,
            "unknown_callable",
            "tool",
            "command -v tool",
            "direct",
            "conditional",
            "available in PATH",
            True,
            "low",
        )

        self.assertEqual(
            candidate("demo", evidence, ["bash >= 4.4"]),
            candidate("demo", evidence, ["bash >= 4.4"]),
        )

    def test_executable_requirements_use_exact_cli_casing(self) -> None:
        callables = {
            "bowtie2",
            "bwa",
            "bwa-mem2",
            "python",
            "Rscript",
            "samtools",
        }
        concepts = {
            "Bowtie2": {"bowtie2"},
            "BWA": {"bwa"},
            "BWA-MEM2": {"bwa-mem2"},
            "Python": {"python"},
            "SAMtools": {"samtools"},
        }

        self.assertIsNone(
            exact_callable_message("bowtie2", callables, concepts),
        )
        self.assertIsNone(exact_callable_message("bwa", callables, concepts))
        self.assertIsNone(
            exact_callable_message("bwa-mem2", callables, concepts),
        )
        self.assertIsNone(
            exact_callable_message("python", callables, concepts),
        )
        self.assertIsNone(
            exact_callable_message("Rscript", callables, concepts),
        )
        self.assertIsNone(
            exact_callable_message("samtools", callables, concepts),
        )

    def test_unsupported_display_capitalization_is_rejected(self) -> None:
        callables = {
            "bowtie2",
            "bwa",
            "bwa-mem2",
            "python",
            "Rscript",
            "samtools",
        }
        concepts = {
            "Bowtie2": {"bowtie2"},
            "BWA": {"bwa"},
            "BWA-MEM2": {"bwa-mem2"},
            "Python": {"python"},
            "SAMtools": {"samtools"},
        }

        for observed, expected in (
            ("Bowtie2", "bowtie2"),
            ("BWA", "bwa"),
            ("BWA-MEM2", "bwa-mem2"),
            ("Python", "python"),
            ("rscript", "Rscript"),
            ("Samtools", "samtools"),
        ):
            message = exact_callable_message(observed, callables, concepts)

            self.assertIsNotNone(message)
            self.assertIn(f"'{expected}'", str(message))

    def test_python_runtime_floor_accepts_exact_executable_spellings(
        self,
    ) -> None:
        for requirement in (
            "python >= 3.11",
            "python3 >= 3.11",
            "Caller-supplied Python interpreter >= 3.11",
        ):
            rules = self.policy_rule_ids("bash >= 4.4", requirement)

            self.assertNotIn("RUNTIME.PYTHON.MINIMUM", rules)

    def test_python_runtime_floor_rejects_missing_or_old_minimum(self) -> None:
        for requirement in ("python", "python >= 3.10", "Python >= 3.11"):
            rules = self.policy_rule_ids("bash >= 4.4", requirement)

            self.assertIn("RUNTIME.PYTHON.MINIMUM", rules)

    def test_runtime_conditions_accept_straight_single_quoted_api_syntax(
        self,
    ) -> None:
        for requirement in (
            "parallel (when '--slurm' is not specified)",
            "bowtie2 (when '--aligner bowtie2' is specified)",
        ):
            rules = self.policy_rule_ids("bash >= 4.4", requirement)

            self.assertNotIn("RUNTIME.INLINE_API", rules)

    def test_runtime_conditions_reject_unquoted_or_malformed_api_syntax(
        self,
    ) -> None:
        for requirement in (
            "parallel (when --slurm is not specified)",
            "parallel (when `--slurm` is not specified)",
            "bowtie2 (when '--aligner' bowtie2 is specified)",
            "bwa-mem2 (when -'-aligner bwa-mem2' is specified)",
        ):
            rules = self.policy_rule_ids("bash >= 4.4", requirement)

            self.assertIn("RUNTIME.INLINE_API", rules)

    def test_sentence_style_resource_capitalization_accepts_prose(
        self,
    ) -> None:
        for requirement in (
            "A compatible Conda environment providing the listed tools",
            "Caller-supplied command executable",
            "Reference FASTA and required index",
        ):
            rules = self.policy_rule_ids("bash >= 4.4", requirement)

            self.assertNotIn("RUNTIME.CAPITALIZATION", rules)

    def test_sentence_style_resource_capitalization_rejects_lowercase_prose(
        self,
    ) -> None:
        for requirement in (
            "a compatible Conda environment providing the listed tools",
            "caller-supplied command executable",
            "reference FASTA and required index",
        ):
            rules = self.policy_rule_ids("bash >= 4.4", requirement)

            self.assertIn("RUNTIME.CAPITALIZATION", rules)

    def test_lowercase_executable_names_are_not_sentence_capitalized(
        self,
    ) -> None:
        for requirement in ("samtools", "bowtie2", "Rscript"):
            rules = self.policy_rule_ids("bash >= 4.4", requirement)

            self.assertNotIn("RUNTIME.CAPITALIZATION", rules)

    def test_command_registry_recognizes_rmdir_as_an_executable(self) -> None:
        root = Path(__file__).resolve().parents[3]
        registry = json.loads(
            (root / "dev/config/command_names.json").read_text(
                encoding="utf-8",
            ),
        )
        callables, _ = validate_command_registry(registry)
        rule_ids = {
            rule_id
            for rule_id, _ in runtime_policy_messages("rmdir", callables)
        }

        self.assertNotIn("RUNTIME.CAPITALIZATION", rule_ids)

    def test_conda_environment_terminology_accepts_resource_and_activation(
        self,
    ) -> None:
        for requirement in (
            "A compatible Conda environment",
            "conda (when the requested environment is not active)",
        ):
            rules = self.policy_rule_ids("bash >= 4.4", requirement)

            self.assertNotIn("RUNTIME.CONDA_ENVIRONMENT", rules)

    def test_conda_environment_terminology_rejects_mamba_or_lowercase_resource(
        self,
    ) -> None:
        for requirement in (
            "A compatible Mamba environment",
            "A compatible Mamba or Conda environment",
            "a compatible Conda environment",
        ):
            rules = self.policy_rule_ids("bash >= 4.4", requirement)

            self.assertIn("RUNTIME.CONDA_ENVIRONMENT", rules)

    def test_redundant_environment_tool_composite_is_rejected(self) -> None:
        rules = self.policy_rule_ids(
            "bash >= 4.4",
            (
                "samtools or conda with a compatible environment providing "
                "samtools"
            ),
        )

        self.assertIn("RUNTIME.REDUNDANT_COMPOSITE", rules)

    def test_genuine_alternative_provider_remains_one_entry(self) -> None:
        requirement = "sha256sum or shasum (when '--dry_run' is not specified)"
        _, entries, grammar = runtime_entries(
            document(
                f"Notes\n-----\n  Runtime requirements:\n    {requirement}",
            ),
        )

        self.assertEqual(entries, [(4, requirement)])
        self.assertEqual(grammar, [])

    def test_runtime_entries_reject_terminal_periods(self) -> None:
        rules = self.policy_rule_ids("bash >= 4.4", "samtools.")

        self.assertIn("RUNTIME.TERMINAL_PUNCTUATION", rules)

    def test_runtime_list_accepts_casefolded_alphabetical_order(self) -> None:
        entries = (
            "A compatible Conda environment providing the listed tools",
            "bash >= 4.4",
            "conda",
            "python >= 3.11",
            "Reference FASTA and required index",
            "samtools",
        )

        self.assertNotIn("RUNTIME.ORDER", self.policy_rule_ids(*entries))

    def test_runtime_list_rejects_out_of_order_entries(self) -> None:
        entries = (
            "bash >= 4.4",
            "A compatible Conda environment providing the listed tools",
            "samtools",
            "python >= 3.11",
        )

        self.assertIn("RUNTIME.ORDER", self.policy_rule_ids(*entries))

    def test_runtime_sort_key_uses_exact_text_as_casefold_tie_breaker(
        self,
    ) -> None:
        from dev.audit.help_runtime_requirements import runtime_sort_key

        values = ["alpha resource", "Alpha resource"]

        self.assertEqual(
            sorted(values, key=runtime_sort_key),
            ["Alpha resource", "alpha resource"],
        )

    def test_runtime_normalizer_preserves_case_conditions_and_singletons(
        self,
    ) -> None:
        from dev.audit.help_runtime_requirements import (
            normalize_runtime_entries,
        )

        entries = [
            "samtools",
            "Reference FASTA (when processing CRAM)",
            "bash >= 4.4",
        ]

        self.assertEqual(
            normalize_runtime_entries(entries),
            [
                "bash >= 4.4",
                "Reference FASTA (when processing CRAM)",
                "samtools",
            ],
        )
        self.assertEqual(normalize_runtime_entries(["Rscript"]), ["Rscript"])

    def test_resource_prose_is_not_treated_as_an_executable(self) -> None:
        registry = {
            "schema_version": 1,
            "commands": [
                {"callable": "python", "conceptual_names": ["Python"]},
                {"callable": "samtools", "conceptual_names": ["SAMtools"]},
            ],
        }
        callables, concepts = validate_command_registry(registry)
        resource = "Reference FASTA and required index (when processing CRAM)"

        self.assertIsNone(
            exact_callable_message(resource, callables, concepts),
        )

    def test_checksum_requirement_is_one_alternative_entry(self) -> None:
        self.assertEqual(
            CHECKSUM_REQUIREMENT,
            "sha256sum or shasum (when '--dry_run' is not specified)",
        )

    def test_settled_owner_requirements_distinguish_execution_from_validation(
        self,
    ) -> None:
        self.assertEqual(
            settled_requirements(
                "lib/bash/core/check_env.sh::check_pgrm_path",
            ),
            ("bash >= 4.4",),
        )
        self.assertEqual(
            settled_requirements(
                "lib/bash/core/construct_find.sh::construct_find",
            ),
            ("bash >= 4.4",),
        )
        self.assertEqual(
            settled_requirements("lib/bash/core/check_unity.sh::check_unity"),
            (
                "bash >= 4.4",
                "awk",
                "gunzip (when processing a gzipped bedGraph)",
            ),
        )

    def test_settled_requirements_preserve_alternatives_and_conditions(
        self,
    ) -> None:
        self.assertEqual(
            settled_requirements(
                "lib/bash/dispatch/manage_parallel.sh::determine_cores",
            ),
            ("bash >= 4.4", "nproc, sysctl, or getconf"),
        )

    def test_c4_alignment_keeps_conditional_providers_without_versions(
        self,
    ) -> None:
        self.assertEqual(
            settled_requirements(
                "lib/bash/workflows/align_fastqs.sh::align_fastqs",
            ),
            (
                "bash >= 4.4",
                "bowtie2 (when '--aligner bowtie2' is specified)",
                "bwa (when '--aligner bwa' is specified)",
                "bwa-mem2 (when '--aligner bwa-mem2' is specified)",
                "grep",
                "Reference FASTA and required index (when writing CRAM)",
                "samtools",
            ),
        )

    def test_c4_filter_contracts_couple_cram_resources_conditionally(
        self,
    ) -> None:
        self.assertEqual(
            settled_requirements(
                (
                    "lib/bash/workflows/filter_alignment.sh::_validate_args_fi"
                    "lter_alignment"
                ),
            ),
            ("bash >= 4.4", "dirname"),
        )
        self.assertEqual(
            settled_requirements(
                "lib/bash/workflows/filter_alignment.sh::_check_chr_alignment",
            ),
            (
                "bash >= 4.4",
                "awk",
                "Reference FASTA and required index (when processing CRAM)",
                "samtools",
                "sort",
                "uniq",
            ),
        )
        self.assertEqual(
            settled_requirements(
                (
                    "lib/bash/workflows/filter_alignment.sh::_finalize_alignme"
                    "nt_filter"
                ),
            ),
            (
                "bash >= 4.4",
                "awk (when 'chk_chr' is true)",
                (
                    "Reference FASTA and required index (when 'chk_chr' is "
                    "true and processing CRAM)"
                ),
                "samtools",
                "sort (when 'chk_chr' is true)",
                "uniq (when 'chk_chr' is true)",
            ),
        )
        self.assertIn(
            "Input BAM or CRAM index",
            settled_requirements(
                "lib/bash/workflows/filter_alignment.sh::filter_alignment_sc",
            ),
        )
        self.assertNotIn(
            "Input BAM or CRAM index",
            settled_requirements(
                "lib/bash/workflows/filter_alignment.sh::filter_alignment_sp",
            ),
        )

    def test_c4_bam_cram_distinguishes_validation_and_resources(
        self,
    ) -> None:
        self.assertEqual(
            settled_requirements(
                (
                    "lib/bash/workflows/process_alignments.sh::_validate_proce"
                    "ss_alignments"
                ),
            ),
            ("bash >= 4.4", "dirname"),
        )
        self.assertEqual(
            settled_requirements(
                (
                    "lib/bash/workflows/process_alignments.sh::qsort_file_alig"
                    "nments"
                ),
            ),
            (
                "bash >= 4.4",
                "Reference FASTA and required index (when processing CRAM)",
                "rm (when processing CRAM)",
                "samtools",
            ),
        )
        self.assertEqual(
            settled_requirements(
                (
                    "lib/bash/workflows/process_alignments.sh::convert_alignme"
                    "nts_bed_python"
                ),
            ),
            (
                "bash >= 4.4",
                "python >= 3.11",
                "Reference FASTA and required index (when processing CRAM)",
            ),
        )

    def test_c4_zcat_is_one_conditional_alternative_provider_contract(
        self,
    ) -> None:
        self.assertEqual(
            DECOMPRESSION_REQUIREMENT,
            "gzip or zcat (when processing a gzipped bedGraph)",
        )
        self.assertEqual(
            settled_requirements(
                "lib/bash/workflows/process_region.sh::check_region_bdg",
            ),
            (
                "bash >= 4.4",
                "cat (when processing a plain-text bedGraph)",
                "gzip or zcat (when processing a gzipped bedGraph)",
            ),
        )

        root = Path(__file__).resolve().parents[3]
        identity = "lib/bash/workflows/process_region.sh::check_region_bdg"
        kind, source, _, _ = owner_units(root)[identity]

        evidence = evidence_for_owner(identity, source, kind, root)
        alternatives = [
            item for item in evidence if item.key == "gzip_or_zcat"
        ]

        self.assertEqual(
            [item.invocation for item in alternatives],
            ["command -v gzip", "command -v zcat"],
        )
        self.assertTrue(
            all(
                item.category == "callable_alternative"
                for item in alternatives
            ),
        )

    def test_c4_atria_runs_generated_command_and_excludes_python(
        self,
    ) -> None:
        self.assertEqual(
            settled_requirements(
                "lib/bash/workflows/process_sequences.sh::trim_fastqs_atria",
            ),
            ("bash >= 4.4", "atria", "dirname"),
        )

        converter = settled_requirements(
            (
                "lib/bash/workflows/process_alignments.sh::convert_alignments_"
                "bed_python"
            ),
        )

        self.assertNotIn("pth_scr_py", converter)

    def test_c4_table_contracts_preserve_direct_and_delegated_parser_tools(
        self,
    ) -> None:
        self.assertEqual(
            settled_requirements(
                "lib/bash/workflows/process_tables.sh::_parse_table_core",
            ),
            (
                "bash >= 4.4",
                "awk",
                "basename (when 'mode' is 'complex')",
                "gawk",
            ),
        )
        self.assertEqual(
            settled_requirements(
                "lib/bash/workflows/process_tables.sh::parse_table_simple",
            ),
            ("bash >= 4.4", "awk", "gawk"),
        )

    def test_c4_settled_conditions_use_exact_grammar_and_no_tool_versions(
        self,
    ) -> None:
        fourth_cohort_prefixes = (
            "lib/bash/workflows/align_fastqs.sh::",
            "lib/bash/workflows/filter_alignment.sh::",
            "lib/bash/help/help_execute_compute_signal.sh::",
            "lib/bash/workflows/process_alignments.sh::",
            "lib/bash/workflows/process_region.sh::",
            "lib/bash/workflows/process_sequences.sh::",
            "lib/bash/workflows/process_tables.sh::",
        )
        contracts = [
            requirement
            for owner, requirements in SETTLED_OWNER_REQUIREMENTS.items()
            if owner.startswith(fourth_cohort_prefixes)
            for requirement in requirements
        ]

        self.assertTrue(contracts)

        for requirement in contracts:
            if " (when " in requirement:
                self.assertIsNotNone(CONDITION.fullmatch(requirement))

            if not re.match(r"^(?:bash|python3?) >= ", requirement):
                self.assertNotRegex(requirement, r"(?:>=|<=|==|>|<)\s*\d")

        self.assertEqual(
            settled_requirements(
                "lib/bash/core/format_outputs.sh::format_print_cmd",
            ),
            (
                "bash >= 4.4",
                "sed (when formatting a non-Slurm command)",
                "awk (when formatting an 'sbatch' command)",
            ),
        )

    def test_c5_scaling_contracts_distinguish_direct_and_delegated_tools(
        self,
    ) -> None:
        self.assertEqual(
            settled_requirements(
                (
                    "lib/bash/workflows/calculate_scaling_factor.sh::_calculat"
                    "e_frag_avg"
                ),
            ),
            (
                "bash >= 4.4",
                "awk (when deriving paired-end fragment length)",
                (
                    "Reference FASTA and required index (when deriving "
                    "fragment length from CRAM)"
                ),
                "samtools",
            ),
        )
        self.assertEqual(
            settled_requirements(
                (
                    "lib/bash/workflows/calculate_scaling_factor.sh::_compute_"
                    "dep_all"
                ),
            ),
            ("bash >= 4.4", "bc"),
        )

    def test_c5_generated_options_do_not_inherit_executed_command_requirements(
        self,
    ) -> None:
        self.assertEqual(
            settled_requirements("bin/submit_compute_signal.sh::set_args_opt"),
            ("bash >= 4.4",),
        )
        self.assertEqual(
            settled_requirements(
                "bin/submit_compute_signal.sh::run_dry_or_wet",
            ),
            (
                "bash >= 4.4",
                "Command-array executable (when 'dry_run' is false)",
                "dirname",
            ),
        )

    def test_c5_signal_tracks_python_environment_and_resources(
        self,
    ) -> None:
        self.assertEqual(
            settled_requirements("bin/submit_compute_signal.sh::run_comp_sig"),
            (
                "bash >= 4.4",
                "A compatible Conda environment providing the listed tools",
                "dirname",
                (
                    "Input BAM or CRAM index (when the selected signal engine "
                    "uses indexed access)"
                ),
                "python >= 3.11",
                "Reference FASTA and required index (when processing CRAM)",
            ),
        )
        self.assertEqual(
            settled_requirements("bin/submit_compute_signal.sh::run_comp_rat"),
            (
                "bash >= 4.4",
                "A compatible Conda environment providing the listed tools",
                "dirname",
                "python >= 3.11",
            ),
        )

    def test_c5_submit_worker_contract_does_not_claim_sbatch(self) -> None:
        contract = settled_requirements(
            "bin/submit_compute_signal.sh::run_task_sig",
        )

        self.assertIn("ln (when run as a Slurm array task)", contract)
        self.assertIn(
            "rm (when run as a wet Slurm array task)",
            contract,
        )
        self.assertNotIn("sbatch", contract)

    def test_c5_dynamic_filter_dispatch_retains_complete_reachable_contract(
        self,
    ) -> None:
        self.assertEqual(
            settled_requirements(
                "bin/submit_filter_alignments.sh::run_filtering",
            ),
            (
                "bash >= 4.4",
                "awk",
                "basename",
                "dirname",
                (
                    "Input BAM or CRAM index (when filtering S. cerevisiae "
                    "alignments)"
                ),
                "mv (when writing S. cerevisiae BAM output)",
                "Reference FASTA and required index (when processing CRAM)",
                "rm",
                "samtools",
                "sort (when 'chk_chr' is true)",
                "uniq (when 'chk_chr' is true)",
            ),
        )

    def test_c6_direct_python3_execution_is_literal_callable_evidence(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = Path(temp_dir) / "direct_python3.sh"
            source.write_text(
                "function direct_python3() {\n    python3 -c 'print(1)'\n}\n",
                encoding="utf-8",
            )

            evidence = evidence_for_owner(
                f"{source.relative_to(root)}::direct_python3",
                source,
                "function",
                root,
            )

        direct = [item for item in evidence if item.key == "python3"]

        self.assertEqual(len(direct), 1)
        self.assertEqual(direct[0].category, "callable")
        self.assertEqual(direct[0].invocation, "python3")

    def test_c6_python3_in_help_fixture_text_is_not_runtime_evidence(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = Path(temp_dir) / "fixture_text.sh"
            source.write_text(
                """
function fixture_text() {
    local show_help
    show_help=$(cat << 'EOM'
Usage
-----
  fixture_text

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Keep a command only as fixture text.
    '''bash
    printf '%s\\n' python3
    '''
EOM
    )
    printf '%s\\n' fixture
}
""",
                encoding="utf-8",
            )

            evidence = evidence_for_owner(
                f"{source.relative_to(root)}::fixture_text",
                source,
                "function",
                root,
            )

        self.assertFalse(any(item.key == "python3" for item in evidence))

    def test_c6_python_discovery_and_loopback_are_separate(
        self,
    ) -> None:
        self.assertEqual(
            settled_requirements("tests/support/test_helpers.sh::find_python"),
            ("bash >= 4.4", "python >= 3.11 or python3 >= 3.11"),
        )
        self.assertEqual(
            settled_requirements(
                "tests/support/test_helpers.sh::find_python_loopback",
            ),
            (
                "bash >= 4.4",
                (
                    "python >= 3.11 or python3 >= 3.11, each with loopback "
                    "socket support"
                ),
            ),
        )

    def test_c6_dynamic_interpreter_contract_does_not_invent_one_spelling(
        self,
    ) -> None:
        self.assertEqual(
            settled_requirements(
                "tests/support/test_helpers.sh::find_port_free",
            ),
            ("bash >= 4.4", "Caller-supplied Python interpreter >= 3.11"),
        )
        self.assertEqual(
            settled_requirements(
                "tests/support/test_helpers.sh::wait_http_local",
            ),
            (
                "bash >= 4.4",
                (
                    "Caller-supplied Python interpreter >= 3.11 with 'urllib' "
                    "support"
                ),
                "sleep",
            ),
        )

    def test_c6_tested_child_requirements_are_not_helper_requirements(
        self,
    ) -> None:
        expected = {
            "tests/support/test_helpers.sh::run_case_filter": (
                "bash >= 4.4",
                "dirname",
                "mkdir",
            ),
            "tests/support/test_helpers.sh::run_case_compute_signal": (
                "bash >= 4.4",
                "dirname",
                "mkdir",
            ),
            "tests/support/test_helpers.sh::run_case_compute_signal_ratio": (
                "bash >= 4.4",
                "dirname",
                "mkdir",
            ),
            "tests/support/test_helpers.sh::run_case_scaling_factor_execute": (
                "bash >= 4.4",
                "dirname",
                "grep",
                "mkdir",
            ),
            (
                "tests/support/test_helpers.sh::run_case_scaling_factor_submit"
                "_part"
            ): (
                "bash >= 4.4",
                "cat",
                "dirname",
                "grep",
                "mkdir",
                "wc",
            ),
        }

        for owner, contract in expected.items():
            contract = settled_requirements(owner)

            self.assertEqual(contract, expected[owner])
            self.assertNotIn("python", contract)
            self.assertNotIn("samtools", contract)

    def test_c6_optional_integration_dependencies_keep_source_conditions(
        self,
    ) -> None:
        self.assertEqual(
            settled_requirements(
                "tests/support/test_helpers.sh::require_env_parallel",
            ),
            (
                "A compatible Conda environment providing parallel",
                "bash >= 4.4",
                "conda (when the requested environment is not active)",
                "dirname",
                "mkdir",
                "parallel",
            ),
        )

    def test_c6_python_candidates_are_fully_dispositioned(self) -> None:
        root = Path(__file__).resolve().parents[3]

        for identity in (
            "tests/support/test_helpers.sh::find_python",
            "tests/support/test_helpers.sh::find_python_loopback",
        ):
            kind, source, _, _ = owner_units(root)[identity]

            evidence = evidence_for_owner(identity, source, kind, root)

            self.assertFalse(
                [
                    item
                    for item in evidence
                    if item.category == "unknown_callable"
                ],
            )


if __name__ == "__main__":
    unittest.main()
