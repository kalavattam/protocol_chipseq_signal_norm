#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: help_runtime_requirements.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""Audit source-grounded Runtime requirements in maintained shell help.

The owner universe deliberately comes from 'help_examples'. This module
adds requirement evidence and grammar checks; it must not grow a competing
owner or compatibility-wrapper parser.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import re
import sys
from pathlib import Path

from help_aliases import alias_chunks
from help_examples import function_help_units, script_help_units
from help_heredoc_reflow import Heredoc, extract_help_heredocs
from help_style import exact_callable_message, load_registry
from shell_validation import MINIMUM_BASH
from source_policy import classify_source

RULE_OWNER_SECTION = "RUNTIME.OWNER.SECTION"
RULE_GRAMMAR = "RUNTIME.GRAMMAR"
RULE_CARDINALITY = "RUNTIME.CARDINALITY"
RULE_CALLABLE = "RUNTIME.CALLABLE"
RULE_VERSION = "RUNTIME.VERSION"
RULE_CONDITION = "RUNTIME.CONDITION"
RULE_MISSING = "RUNTIME.MISSING"
RULE_OVERCLAIM = "RUNTIME.OVERCLAIM"
RULE_INTERNAL = "RUNTIME.INTERNAL_EXPOSED"
RULE_ENV_DEFAULT = "RUNTIME.ENV.DEFAULT_EXPOSED"
RULE_PYTHON_MINIMUM = "RUNTIME.PYTHON.MINIMUM"
RULE_INLINE_API = "RUNTIME.INLINE_API"
RULE_CAPITALIZATION = "RUNTIME.CAPITALIZATION"
RULE_PUNCTUATION = "RUNTIME.TERMINAL_PUNCTUATION"
RULE_CONDA_ENVIRONMENT = "RUNTIME.CONDA_ENVIRONMENT"
RULE_REDUNDANT_COMPOSITE = "RUNTIME.REDUNDANT_COMPOSITE"
RULE_ORDER = "RUNTIME.ORDER"

RUNTIME_LABEL = re.compile(r"^(?P<indent> *)Runtime requirements:\s*$")
ITEM = re.compile(r"^(?P<indent> *)(?P<marker>- )?(?P<text>\S.*)$")
CONDITION = re.compile(r"^(?P<requirement>.+?) \(when (?P<condition>[^()]+)\)$")
VERSION = re.compile(r"^(?P<name>[A-Za-z][A-Za-z0-9_.-]*)\s*(?P<op>>=|<=|==|>|<)\s*(?P<value>\S+)$")
COMMAND = re.compile(r"^\s*(?P<name>[A-Za-z][A-Za-z0-9_.-]*)\b")
COMMAND_V = re.compile(r"\bcommand\s+-v\s+(?P<name>[A-Za-z][A-Za-z0-9_.-]*)")
INTERNAL = re.compile(r"(?:^|\s)(?:bin/|lib/bash/|src/|install/scripts/|tests/|[A-Za-z_][A-Za-z0-9_]*\.py\b)")
FUNCTION_START = re.compile(r"^\s*(?:function\s+)?(?P<name>[A-Za-z_][A-Za-z0-9_]*)\s*\(\s*\)\s*\{")

# This deliberately small set prevents a regex scan from turning shell syntax,
# local functions, and arbitrary variables into asserted requirements.  Unknown
# literal command heads are emitted as semantic candidates instead.
KNOWN_CALLABLES = frozenset({
    "awk", "basename", "bash", "bc", "cat", "conda", "cp", "curl", "cut", "getconf", "gunzip",
    "dirname", "find", "gawk", "git", "grep", "gzip", "head", "ln", "mamba",
    "mkdir", "mktemp", "mv", "nproc", "parallel", "pbzip2", "pigz", "python", "python3", "realpath",
    "rm", "Rscript", "samtools", "sbatch", "sed", "sort", "sysctl", "tail", "tar", "tr",
    "uniq", "wget", "zcat",
})
CHECKSUM_PROVIDERS = frozenset({"sha256sum", "shasum"})
CHECKSUM_REQUIREMENT = "sha256sum or shasum (when '--dry_run' is not specified)"
DECOMPRESSION_PROVIDERS = frozenset({"gzip", "zcat"})
DECOMPRESSION_REQUIREMENT = "gzip or zcat (when processing a gzipped bedGraph)"
SHELL_WORDS = frozenset({
    "if", "then", "else", "elif", "fi", "for", "while", "do", "done", "case",
    "esac", "function", "local", "readonly", "return", "exit", "export", "source",
    "set", "echo", "printf", "test", "true", "false", "eval", "exec", "time",
})
PYTHON_MINIMUM = "3.11"
REDUNDANT_COMPOSITE = re.compile(
    r"\bor conda with a compatible environment providing\b",
    re.I,
)
CONDA_ENVIRONMENT = re.compile(
    r"\b(?:Mamba(?: or Conda)?|Conda or Mamba) environment\b",
    re.I,
)
OPTION_TOKEN = re.compile(r"--[A-Za-z0-9][A-Za-z0-9_-]*")
SNAKE_API_TOKEN = re.compile(r"\b[a-z][a-z0-9]*_[a-z0-9_]+\b")
SINGLE_QUOTED = re.compile(r"'[^'\n]+'")
PYTHON_EXECUTABLE = re.compile(r"^(?P<name>python3?|Python)(?:\s+(?P<rest>.*))?$")
CALLER_PYTHON = re.compile(r"^Caller-supplied Python interpreter\b")
OLD_ENVIRONMENT = "A compatible Mamba or Conda environment"
CONDA_ACTIVATION = "conda (when the requested environment is not active)"
SAMTOOLS_COMPOSITE = (
    "samtools or conda with a compatible environment providing samtools"
)
SAMTOOLS_FALLBACK = (
    "A compatible Conda environment providing samtools "
    "(when 'samtools' is unavailable)",
    "conda (when 'samtools' is unavailable)",
    "samtools",
)

# Generic command extraction cannot safely distinguish a validator's dynamic
# input from an executable it runs, nor can it prove generated arrays execute.
# These settled entries therefore retain the exact, source-reviewed contracts
# for the completed validation and orchestration cohorts.  The table is keyed
# by stable owner identity rather than filename branches.
SETTLED_OWNER_REQUIREMENTS: dict[str, tuple[str, ...]] = {
    "lib/bash/core/check_source.sh::err_source_only": (
        "bash >= 4.4", "basename",
    ),
    "lib/bash/core/check_unity.sh::check_unity": (
        "bash >= 4.4", "awk", "gunzip (when processing a gzipped bedGraph)",
    ),
    "lib/bash/core/check_env.sh::check_env_installed": (
        "bash >= 4.4", "awk", "conda",
    ),
    "lib/bash/core/format_outputs.sh::format_print_cmd": (
        "bash >= 4.4",
        "sed (when formatting a non-Slurm command)",
        "awk (when formatting an 'sbatch' command)",
    ),
    "lib/bash/core/handle_env.sh::_handle_env_deactivate": (
        "bash >= 4.4", "basename", "conda (when deactivating an active environment)",
    ),
    "lib/bash/core/handle_env.sh::_handle_env_activate": (
        "bash >= 4.4", "basename", "conda", "grep",
    ),
    "lib/bash/core/handle_env.sh::handle_env": (
        "bash >= 4.4",
        "basename (when activating or deactivating an environment)",
        "conda (when activating or deactivating an environment)",
        "grep (when activating an environment)",
    ),
    "lib/bash/dispatch/manage_parallel.sh::determine_cores": (
        "bash >= 4.4", "nproc, sysctl, or getconf",
    ),
    "lib/bash/dispatch/manage_parallel.sh::set_params_parallel": (
        "bash >= 4.4", "nproc, sysctl, or getconf",
    ),
    "lib/bash/dispatch/manage_slurm.sh::set_logs_slurm": (
        "bash >= 4.4", "ln",
    ),
    "lib/bash/dispatch/run_python.sh::_resolve_dir_rep_run_py": (
        "bash >= 4.4", "dirname",
    ),
    "lib/bash/dispatch/run_python.sh::to_module": (
        "bash >= 4.4", "dirname",
    ),
    "lib/bash/dispatch/run_python.sh::run_py": (
        "bash >= 4.4", "dirname", "python >= 3.11",
    ),
    "lib/bash/core/source_helpers.sh::source_once": (
        "bash >= 4.4", "basename", "dirname",
    ),
    "lib/bash/core/source_helpers.sh::source_helpers": (
        "bash >= 4.4", "basename", "dirname",
    ),
    "lib/bash/core/wrap_cmd.sh::get_submit_logs": (
        "bash >= 4.4", "basename",
    ),
    "lib/bash/workflows/align_fastqs.sh::align_fastqs": (
        "bash >= 4.4",
        "bowtie2 (when '--aligner bowtie2' is specified)",
        "bwa (when '--aligner bwa' is specified)",
        "bwa-mem2 (when '--aligner bwa-mem2' is specified)",
        "grep",
        "Reference FASTA and required index (when writing CRAM)",
        "samtools",
    ),
    "lib/bash/workflows/filter_alignment.sh::_validate_args_filter_alignment": (
        "bash >= 4.4", "dirname",
    ),
    "lib/bash/workflows/filter_alignment.sh::_check_chr_alignment": (
        "bash >= 4.4",
        "awk",
        "Reference FASTA and required index (when processing CRAM)",
        "samtools",
        "sort",
        "uniq",
    ),
    "lib/bash/workflows/filter_alignment.sh::_finalize_alignment_filter": (
        "bash >= 4.4",
        "awk (when 'chk_chr' is true)",
        "Reference FASTA and required index (when 'chk_chr' is true and processing CRAM)",
        "samtools",
        "sort (when 'chk_chr' is true)",
        "uniq (when 'chk_chr' is true)",
    ),
    "lib/bash/workflows/filter_alignment.sh::_filter_sam_chr": (
        "bash >= 4.4", "awk",
    ),
    "lib/bash/workflows/filter_alignment.sh::filter_alignment_sc": (
        "bash >= 4.4",
        "awk",
        "basename",
        "dirname",
        "Input BAM or CRAM index",
        "mv (when writing BAM)",
        "Reference FASTA and required index (when processing CRAM)",
        "rm",
        "samtools",
        "sort (when '--chk_chr' is specified)",
        "uniq (when '--chk_chr' is specified)",
    ),
    "lib/bash/workflows/filter_alignment.sh::filter_alignment_sp": (
        "bash >= 4.4",
        "awk",
        "basename",
        "dirname",
        "Reference FASTA and required index (when processing CRAM)",
        "rm",
        "samtools",
        "sort (when '--chk_chr' is specified)",
        "uniq (when '--chk_chr' is specified)",
    ),
    "lib/bash/help/help_execute_compute_signal.sh::help_execute_compute_signal": (
        "bash >= 4.4",
        "A compatible Conda environment providing the listed tools",
        "basename",
        "dirname",
        "parallel (when '--slurm' is not specified and multiple jobs are run)",
        "python >= 3.11",
        "rm (when '--slurm' is not specified and multiple jobs are run)",
        "sbatch (when '--slurm' is specified)",
        "tr",
    ),
    "lib/bash/workflows/process_alignments.sh::_validate_process_alignments": (
        "bash >= 4.4", "dirname",
    ),
    "lib/bash/workflows/process_alignments.sh::qsort_file_alignments": (
        "bash >= 4.4",
        "Reference FASTA and required index (when processing CRAM)",
        "rm (when processing CRAM)",
        "samtools",
    ),
    "lib/bash/workflows/process_alignments.sh::convert_alignments_bed_awk": (
        "bash >= 4.4",
        "awk",
        "gzip",
        "Reference FASTA and required index (when processing CRAM)",
        "samtools",
        "sort",
    ),
    "lib/bash/workflows/process_alignments.sh::convert_alignments_bed_python": (
        "bash >= 4.4",
        "python >= 3.11",
        "Reference FASTA and required index (when processing CRAM)",
    ),
    "lib/bash/workflows/process_region.sh::check_region_bam": (
        "bash >= 4.4", "samtools",
    ),
    "lib/bash/workflows/process_region.sh::check_region_bdg": (
        "bash >= 4.4",
        "cat (when processing a plain-text bedGraph)",
        "gzip or zcat (when processing a gzipped bedGraph)",
    ),
    "lib/bash/workflows/process_sequences.sh::check_seq_type": (
        "bash >= 4.4", "grep", "samtools",
    ),
    "lib/bash/workflows/process_sequences.sh::check_string_fastqs": (
        "bash >= 4.4", "awk",
    ),
    "lib/bash/workflows/process_sequences.sh::parse_fastq_entry": (
        "bash >= 4.4", "awk", "basename",
    ),
    "lib/bash/workflows/process_sequences.sh::trim_fastqs_atria": (
        "bash >= 4.4", "atria", "dirname",
    ),
    "lib/bash/workflows/process_sequences.sh::pair_fastqs": (
        "bash >= 4.4", "awk",
    ),
    "lib/bash/workflows/process_tables.sh::check_table": (
        "bash >= 4.4", "awk",
    ),
    "lib/bash/workflows/process_tables.sh::check_table_column": (
        "bash >= 4.4", "awk", "tr",
    ),
    "lib/bash/workflows/process_tables.sh::extract_field_str": (
        "bash >= 4.4", "awk",
    ),
    "lib/bash/workflows/process_tables.sh::_validate_args_table": (
        "bash >= 4.4", "awk",
    ),
    "lib/bash/workflows/process_tables.sh::_parse_table_core": (
        "bash >= 4.4",
        "awk",
        "basename (when 'mode' is 'complex')",
        "gawk",
    ),
    "lib/bash/workflows/process_tables.sh::parse_table": (
        "bash >= 4.4", "awk", "basename", "gawk",
    ),
    "lib/bash/workflows/process_tables.sh::parse_table_simple": (
        "bash >= 4.4", "awk", "gawk",
    ),
    "lib/bash/workflows/calculate_scaling_factor.sh::_set_ref_arg_cram": (
        "bash >= 4.4",
        "Reference FASTA and required index (when processing CRAM)",
    ),
    "lib/bash/workflows/calculate_scaling_factor.sh::_detect_typ_aln": (
        "bash >= 4.4",
        "cut",
        "head",
        "Reference FASTA and required index (when processing CRAM)",
        "samtools",
    ),
    "lib/bash/workflows/calculate_scaling_factor.sh::_resolve_typ_fil": (
        "bash >= 4.4",
        "cut (when 'pref' is 'auto' or empty)",
        "head (when 'pref' is 'auto' or empty)",
        "Reference FASTA and required index (when 'pref' is 'auto' or empty and processing CRAM)",
        "samtools (when 'pref' is 'auto' or empty)",
    ),
    "lib/bash/workflows/calculate_scaling_factor.sh::_count_alignments": (
        "bash >= 4.4",
        "Reference FASTA and required index (when processing CRAM)",
        "samtools",
    ),
    "lib/bash/workflows/calculate_scaling_factor.sh::_calculate_frag_avg": (
        "bash >= 4.4",
        "awk (when deriving paired-end fragment length)",
        "Reference FASTA and required index (when deriving fragment length from CRAM)",
        "samtools",
    ),
    "lib/bash/workflows/calculate_scaling_factor.sh::_compute_scl_fct": (
        "bash >= 4.4", "python >= 3.11",
    ),
    "lib/bash/workflows/calculate_scaling_factor.sh::_parse_metadata": (
        "bash >= 4.4", "python >= 3.11",
    ),
    "lib/bash/workflows/calculate_scaling_factor.sh::_calculate_dep_fct": (
        "bash >= 4.4", "bc",
    ),
    "lib/bash/workflows/calculate_scaling_factor.sh::_calculate_dep_arr": (
        "bash >= 4.4", "bc",
    ),
    "lib/bash/workflows/calculate_scaling_factor.sh::_compute_dep_all": (
        "bash >= 4.4", "bc",
    ),
    "lib/bash/workflows/calculate_scaling_factor.sh::process_samp_siq": (
        "bash >= 4.4",
        "A compatible Conda environment providing the listed tools",
        "awk (when deriving paired-end fragment length)",
        "cut (when auto-detecting alignment layout)",
        "dirname",
        "head (when auto-detecting alignment layout)",
        "python >= 3.11",
        "Reference FASTA and required index (when processing CRAM)",
        "samtools (when alignment depth, fragment length, or layout must be derived)",
    ),
    "lib/bash/workflows/calculate_scaling_factor.sh::process_samp_spike": (
        "bash >= 4.4",
        "A compatible Conda environment providing the listed tools",
        "cut (when auto-detecting alignment layout)",
        "dirname",
        "head (when auto-detecting alignment layout)",
        "python >= 3.11",
        "Reference FASTA and required index (when processing CRAM)",
        "samtools (when alignment depth or layout must be derived)",
    ),
    "bin/submit_align_fastqs.sh::parse_entry_align_fastq": (
        "bash >= 4.4", "awk", "basename",
    ),
    "bin/submit_align_fastqs.sh::run_alignment": (
        "bash >= 4.4",
        "bowtie2 (when 'aligner' is 'bowtie2')",
        "bwa (when 'aligner' is 'bwa')",
        "bwa-mem2 (when 'aligner' is 'bwa-mem2')",
        "grep",
        "Reference FASTA and required index (when writing CRAM)",
        "samtools",
    ),
    "bin/submit_calculate_scaling_factor.sh::derive_samp_sf": (
        "bash >= 4.4", "basename",
    ),
    "bin/submit_compute_signal.sh::process_io": (
        "bash >= 4.4", "basename",
    ),
    "bin/submit_compute_signal.sh::run_dry_or_wet": (
        "bash >= 4.4",
        "Command-array executable (when 'dry_run' is false)",
        "dirname",
    ),
    "bin/submit_compute_signal.sh::run_comp_sig": (
        "bash >= 4.4",
        "A compatible Conda environment providing the listed tools",
        "dirname",
        "Input BAM or CRAM index (when the selected signal engine uses indexed access)",
        "python >= 3.11",
        "Reference FASTA and required index (when processing CRAM)",
    ),
    "bin/submit_compute_signal.sh::run_comp_rat": (
        "bash >= 4.4",
        "A compatible Conda environment providing the listed tools",
        "dirname",
        "python >= 3.11",
    ),
    "bin/submit_compute_signal.sh::task_pro": (
        "bash >= 4.4",
        "basename",
        "ln (when run as a Slurm array task)",
    ),
    "bin/submit_compute_signal.sh::task_epi": (
        "bash >= 4.4", "rm (when run as a wet Slurm array task)",
    ),
    "bin/submit_compute_signal.sh::run_task_sig": (
        "bash >= 4.4",
        "A compatible Conda environment providing the listed tools",
        "basename",
        "dirname",
        "Input BAM or CRAM index (when the selected signal engine uses indexed access)",
        "ln (when run as a Slurm array task)",
        "python >= 3.11",
        "Reference FASTA and required index (when processing CRAM)",
        "rm (when run as a wet Slurm array task)",
    ),
    "bin/submit_compute_signal.sh::run_task_rat": (
        "bash >= 4.4",
        "A compatible Conda environment providing the listed tools",
        "basename",
        "dirname",
        "ln (when run as a Slurm array task)",
        "python >= 3.11",
        "rm (when run as a wet Slurm array task)",
    ),
    "bin/submit_compute_signal.sh::run_task_coord": (
        "bash >= 4.4",
        "A compatible Conda environment providing the listed tools",
        "basename",
        "dirname",
        "ln (when run as a Slurm array task)",
        "python >= 3.11",
        "Reference FASTA and required index (when processing CRAM)",
        "rm (when run as a wet Slurm array task)",
    ),
    "bin/submit_filter_alignments.sh::parse_filter_alignment_entry": (
        "bash >= 4.4", "basename",
    ),
    "bin/submit_filter_alignments.sh::run_filtering": (
        "bash >= 4.4",
        "awk",
        "basename",
        "dirname",
        "Input BAM or CRAM index (when filtering S. cerevisiae alignments)",
        "mv (when writing S. cerevisiae BAM output)",
        "Reference FASTA and required index (when processing CRAM)",
        "rm",
        "samtools",
        "sort (when 'chk_chr' is true)",
        "uniq (when 'chk_chr' is true)",
    ),
    "bin/submit_trim_fastqs.sh::parse_entry_trim_fastq": (
        "bash >= 4.4", "awk", "basename",
    ),
    "tests/support/fixture_helpers.sh::gzip_n": (
        "bash >= 4.4", "gzip",
    ),
    "tests/support/fixture_helpers.sh::rm_file": (
        "bash >= 4.4", "rm",
    ),
    "tests/support/fixture_helpers.sh::rm_files": (
        "bash >= 4.4", "rm",
    ),
    "tests/support/test_helpers.sh::assert_cram_count": (
        "A compatible Conda environment providing samtools (when 'samtools' is unavailable)",
        "bash >= 4.4",
        "conda (when 'samtools' is unavailable)",
        "dirname",
        "grep",
        "mkdir",
        "Reference FASTA for CRAM decoding",
        "samtools",
    ),
    "tests/support/test_helpers.sh::assert_fastq_gzip": (
        "bash >= 4.4",
        "awk",
        "cmp (when a source fixture is supplied)",
        "grep",
        "gzip",
    ),
    "tests/support/test_helpers.sh::assert_file_exact_line": (
        "bash >= 4.4", "cat", "wc",
    ),
    "tests/support/test_helpers.sh::assert_files_equal": (
        "bash >= 4.4", "cmp",
    ),
    "tests/support/test_helpers.sh::assert_filter_alignments_pg_header": (
        "A compatible Conda environment providing samtools (when 'samtools' is unavailable)",
        "bash >= 4.4",
        "conda (when 'samtools' is unavailable)",
        "dirname",
        "grep",
        "mkdir",
        "Reference FASTA (when inspecting CRAM)",
        "samtools",
    ),
    "tests/support/test_helpers.sh::assert_pattern_absent": (
        "bash >= 4.4", "grep",
    ),
    "tests/support/test_helpers.sh::assert_pattern_found": (
        "bash >= 4.4", "grep",
    ),
    "tests/support/test_helpers.sh::assert_scaling_factor_header": (
        "bash >= 4.4", "grep",
    ),
    "tests/support/test_helpers.sh::build_filter_alignments_fixture_bam": (
        "A compatible Conda environment providing samtools (when 'samtools' is unavailable)",
        "bash >= 4.4",
        "conda (when 'samtools' is unavailable)",
        "dirname",
        "mkdir",
        "samtools",
    ),
    "tests/support/test_helpers.sh::build_filter_alignments_fixture_cram": (
        "A compatible Conda environment providing samtools (when 'samtools' is unavailable)",
        "bash >= 4.4",
        "conda (when 'samtools' is unavailable)",
        "dirname",
        "mkdir",
        "Reference FASTA for CRAM writing",
        "samtools",
    ),
    "tests/support/test_helpers.sh::check_env_cmds": (
        "A compatible Conda environment",
        "bash >= 4.4",
        "conda (when the requested environment is not active)",
    ),
    "tests/support/test_helpers.sh::find_port_free": (
        "bash >= 4.4", "Caller-supplied Python interpreter >= 3.11",
    ),
    "tests/support/test_helpers.sh::find_python": (
        "bash >= 4.4", "python >= 3.11 or python3 >= 3.11",
    ),
    "tests/support/test_helpers.sh::find_python_loopback": (
        "bash >= 4.4",
        "python >= 3.11 or python3 >= 3.11, each with loopback socket support",
    ),
    "tests/support/test_helpers.sh::require_env_atria": (
        "A compatible Conda environment providing the listed tools",
        "bash >= 4.4",
        "atria",
        "conda (when the requested environment is not active)",
        "dirname",
        "gzip",
        "mkdir",
        "pbzip2",
        "pigz",
    ),
    "tests/support/test_helpers.sh::require_env_download": (
        "A compatible Conda environment providing the listed tools",
        "bash >= 4.4",
        "conda (when the requested environment is not active)",
        "dirname",
        "gzip",
        "mkdir",
        "wget",
    ),
    "tests/support/test_helpers.sh::require_env_parallel": (
        "A compatible Conda environment providing parallel",
        "bash >= 4.4",
        "conda (when the requested environment is not active)",
        "dirname",
        "mkdir",
        "parallel",
    ),
    "tests/support/test_helpers.sh::require_env_project": (
        "A compatible Conda environment",
        "bash >= 4.4",
        "awk (when no non-base project environment is active)",
        "conda (when no non-base project environment is active)",
    ),
    "tests/support/test_helpers.sh::run_capture": (
        "bash >= 4.4", "Caller-supplied command executable", "dirname", "mkdir",
    ),
    "tests/support/test_helpers.sh::run_case_compute_signal": (
        "bash >= 4.4", "dirname", "mkdir",
    ),
    "tests/support/test_helpers.sh::run_case_compute_signal_ratio": (
        "bash >= 4.4", "dirname", "mkdir",
    ),
    "tests/support/test_helpers.sh::run_case_filter": (
        "bash >= 4.4", "dirname", "mkdir",
    ),
    "tests/support/test_helpers.sh::run_case_scaling_factor_execute": (
        "bash >= 4.4", "dirname", "grep", "mkdir",
    ),
    "tests/support/test_helpers.sh::run_case_scaling_factor_submit_part": (
        "bash >= 4.4", "cat", "dirname", "grep", "mkdir", "wc",
    ),
    "tests/support/test_helpers.sh::run_samtools": (
        "A compatible Conda environment providing samtools (when 'samtools' is unavailable)",
        "bash >= 4.4",
        "conda (when 'samtools' is unavailable)",
        "samtools",
    ),
    "tests/support/test_helpers.sh::wait_http_local": (
        "bash >= 4.4",
        "Caller-supplied Python interpreter >= 3.11 with 'urllib' support",
        "sleep",
    ),
    "tests/support/test_helpers.sh::warn_help_pattern_missing": (
        "bash >= 4.4", "grep",
    ),
}


def settled_requirements(identity: str) -> tuple[str, ...]:
    """Return source-reviewed requirements where generic inference is unsafe."""

    return SETTLED_OWNER_REQUIREMENTS.get(identity, ("bash >= 4.4",))


@dataclasses.dataclass(frozen=True)
class Finding:
    """One deterministic Runtime-requirements mismatch."""

    rule_id: str
    owner: str
    path: str
    line: int
    message: str
    expected: tuple[str, ...] = ()
    observed: tuple[str, ...] = ()

    def as_dict(self) -> dict[str, object]:
        return dataclasses.asdict(self)

    def format(self) -> str:
        return f"{self.rule_id}: {self.path}:{self.line}: owner={self.owner}; {self.message}"


@dataclasses.dataclass(frozen=True)
class Evidence:
    """One bounded source observation for an owner."""

    owner: str
    path: str
    line: int
    category: str
    key: str
    invocation: str
    relation: str
    condition_kind: str
    condition: str
    public_interface: bool
    confidence: str

    def as_dict(self) -> dict[str, object]:
        return dataclasses.asdict(self)


def owner_units(root: Path) -> dict[str, tuple[str, Path, Path, Heredoc]]:
    """Return every Examples-owned help document without rediscovering owners."""

    result: dict[str, tuple[str, Path, Path, Heredoc]] = {}
    for identity, (source, documentation, heredoc) in script_help_units(root).items():
        result[identity] = ("script", source, documentation, heredoc)
    for identity, (source, heredoc) in function_help_units(root).items():
        result[identity] = ("function", source, source, heredoc)
    return dict(sorted(result.items()))


def runtime_entries(heredoc: Heredoc) -> tuple[int | None, list[tuple[int, str]], list[str]]:
    """Extract a flat Runtime requirements list from one rendered help unit."""

    rows = list(heredoc.lines)
    label: int | None = None
    entries: list[tuple[int, str]] = []
    grammar: list[str] = []
    for index, (number, line) in enumerate(rows):
        label_match = RUNTIME_LABEL.fullmatch(line)
        if label_match is None:
            continue
        if label is not None:
            grammar.append("Runtime requirements may appear only once")
            continue
        label = number
        base = len(label_match.group("indent"))
        for content_number, content in rows[index + 1 :]:
            if not content.strip():
                break
            indent = len(content) - len(content.lstrip(" "))
            if indent <= base:
                break
            match = ITEM.fullmatch(content)
            if match is None or indent != base + 2:
                grammar.append(f"line {content_number} is not a flat entry")
                continue
            text = match.group("text")
            if text.endswith(":") or text.startswith(("+ ", "- Programs", "- Environment", "- Files", "- Shell scripts", "- Python scripts")):
                grammar.append(f"line {content_number} uses a category or continuation")
                continue
            entries.append((content_number, (match.group("marker") or "") + text))
        break
    return label, entries, grammar


def runtime_sort_key(text: str) -> tuple[str, str]:
    """Return the case-insensitive Runtime ordering key with an exact tie-breaker."""

    return text.casefold(), text


def normalize_runtime_entries(entries: list[str]) -> list[str]:
    """Sort Runtime entries without changing their displayed text."""

    return sorted(entries, key=runtime_sort_key)


def normalize_condition_wording(condition: str) -> str:
    """Normalize API delimiters in one already source-reviewed condition."""

    replacements = {
        "aligner is bwa-mem2": "'aligner' is 'bwa-mem2'",
        "aligner is bowtie2": "'aligner' is 'bowtie2'",
        "aligner is bwa": "'aligner' is 'bwa'",
        "mode is complex": "'mode' is 'complex'",
        "pref is auto or empty": "'pref' is 'auto' or empty",
        "formatting an sbatch command": "formatting an 'sbatch' command",
    }
    for old, new in replacements.items():
        condition = condition.replace(old, new)
    condition = re.sub(
        r"(?<!')(?P<api>--aligner\s+(?:bwa-mem2|bowtie2|bwa)|"
        r"--[A-Za-z0-9][A-Za-z0-9_-]*)(?!')",
        lambda match: f"'{match.group('api')}'",
        condition,
    )
    ranges = quoted_ranges(condition)
    matches = list(SNAKE_API_TOKEN.finditer(condition))
    for match in reversed(matches):
        if position_is_quoted(match.start(), ranges):
            continue
        condition = (
            condition[: match.start()]
            + f"'{match.group(0)}'"
            + condition[match.end() :]
        )
    return condition


def normalize_runtime_wording(value: str) -> str:
    """Normalize one high-confidence Runtime prose fragment."""

    head, condition = requirement_parts(value)
    if head.startswith(OLD_ENVIRONMENT):
        head = head.replace(OLD_ENVIRONMENT, "A compatible Conda environment", 1)
    if head.startswith("a compatible Conda environment"):
        head = "A" + head[1:]
    if head in {
        "A compatible Conda environment",
        "A compatible Conda environment providing parallel",
        "A compatible Conda environment providing the listed tools",
    } and condition in {
        "the requested environment is not active",
        "no non-base project environment is active",
    }:
        condition = None
    if head == "python":
        head = "python >= 3.11"
    elif head == "python3":
        head = "python3 >= 3.11"
    elif head == "python or python3":
        head = "python >= 3.11 or python3 >= 3.11"
    elif head == "python or python3 with loopback socket support":
        head = (
            "python >= 3.11 or python3 >= 3.11, "
            "each with loopback socket support"
        )
    if head == "caller-supplied Python interpreter":
        head = "Caller-supplied Python interpreter >= 3.11"
    elif head == "caller-supplied Python interpreter with urllib support":
        head = "Caller-supplied Python interpreter >= 3.11 with 'urllib' support"
    elif head == "caller-supplied command executable":
        head = "Caller-supplied command executable"
    elif head == "command-array executable":
        head = "Command-array executable"
    elif head.startswith("input BAM or CRAM index"):
        head = "Input" + head[len("input") :]
    if condition is not None:
        condition = normalize_condition_wording(condition)
        return f"{head} (when {condition})"
    return head


def normalize_runtime_contract(entries: list[str]) -> list[str]:
    """Normalize and sort one source-reviewed flat Runtime contract."""

    normalized: list[str] = []
    add_conda_activation = any(entry.startswith(OLD_ENVIRONMENT) for entry in entries)
    manager_alternative = (
        "mamba" in entries
        and "conda (when mamba is unavailable)" in entries
    )
    for entry in entries:
        if manager_alternative and entry in {
            "mamba",
            "conda (when mamba is unavailable)",
        }:
            continue
        if entry == SAMTOOLS_COMPOSITE:
            normalized.extend(SAMTOOLS_FALLBACK)
            continue
        normalized.append(normalize_runtime_wording(entry))
    if manager_alternative:
        normalized.append("mamba or conda")
    if add_conda_activation:
        normalized.append(CONDA_ACTIVATION)
    return normalize_runtime_entries(list(dict.fromkeys(normalized)))


def rewrite_runtime_blocks(text: str) -> tuple[str, int]:
    """Normalize recognized Runtime blocks without touching other shell text."""

    lines = text.splitlines()
    replacements: list[tuple[int, int, list[str]]] = []
    for heredoc in extract_help_heredocs(text):
        label, entries_with_lines, _ = runtime_entries(heredoc)
        if label is None or not entries_with_lines:
            continue
        values = [
            entry[2:] if entry.startswith("- ") else entry
            for _, entry in entries_with_lines
        ]
        normalized = normalize_runtime_contract(values)
        indent = len(lines[label - 1]) - len(lines[label - 1].lstrip(" ")) + 2
        marker = "- " if len(normalized) > 1 else ""
        rendered = [" " * indent + marker + entry for entry in normalized]
        start = entries_with_lines[0][0] - 1
        end = entries_with_lines[-1][0]
        if lines[start:end] != rendered:
            replacements.append((start, end, rendered))
    for start, end, rendered in reversed(replacements):
        lines[start:end] = rendered
    suffix = "\n" if text.endswith("\n") else ""
    return "\n".join(lines) + suffix, len(replacements)


def rewrite_repository(root: Path) -> tuple[int, int]:
    """Normalize every Runtime documentation owner in place."""

    paths = sorted({unit[2] for unit in owner_units(root).values()})
    changed_files = 0
    changed_blocks = 0
    for path in paths:
        before = path.read_text(encoding="utf-8")
        after, count = rewrite_runtime_blocks(before)
        if after == before:
            continue
        path.write_text(after, encoding="utf-8")
        changed_files += 1
        changed_blocks += count
    return changed_files, changed_blocks


def requirement_parts(value: str) -> tuple[str, str | None]:
    """Return one requirement head and its optional source-grounded condition."""

    match = CONDITION.fullmatch(value)
    if match is None:
        return value, None
    return match.group("requirement"), match.group("condition")


def quoted_ranges(text: str) -> list[tuple[int, int]]:
    """Return straight-single-quoted ranges in one prose fragment."""

    return [(match.start(), match.end()) for match in SINGLE_QUOTED.finditer(text)]


def position_is_quoted(position: int, ranges: list[tuple[int, int]]) -> bool:
    """Return whether one character position is within a quoted range."""

    return any(start < position < end for start, end in ranges)


def is_executable_requirement(value: str, callables: set[str]) -> bool:
    """Classify an entry as exact executable/API syntax rather than prose."""

    head, _ = requirement_parts(value)
    known = set(KNOWN_CALLABLES) | set(CHECKSUM_PROVIDERS) | set(callables)
    tokens = re.findall(r"[A-Za-z][A-Za-z0-9_.-]*", head)
    if not tokens:
        return False
    first = tokens[0]
    if first in known:
        return True
    provider_tokens = [token for token in tokens if token in known]
    return bool(provider_tokens) and re.search(r"\b(?:or|and)\b", head) is not None


def python_minimum_message(value: str) -> str | None:
    """Return a diagnostic for a Python interpreter requirement below policy."""

    head, _ = requirement_parts(value)
    if CALLER_PYTHON.match(head):
        if re.search(r"\bPython interpreter >= 3\.11\b", head) is None:
            return "caller-supplied interpreters must require 'Python interpreter >= 3.11'"
        return None
    if head.startswith("caller-supplied Python interpreter"):
        return "caller-supplied interpreter prose must begin with 'Caller-supplied'"
    if re.match(r"^(?:python|python3|Python)\b", head) is None:
        return None
    if head.startswith("Python"):
        return "an exact python executable requirement must use lowercase callable spelling"
    providers = re.findall(r"\b(?:python3|python)(?:\s*>=\s*[^\s,]+)?", head)
    if not providers:
        return None
    expected = {
        provider.split()[0] + f" >= {PYTHON_MINIMUM}"
        for provider in providers
    }
    observed = {provider.strip() for provider in providers}
    if observed != expected:
        return "every python/python3 provider must use the minimum '>= 3.11'"
    return None


def inline_api_message(
    value: str,
    callables: set[str],
) -> str | None:
    """Return a diagnostic for malformed or unquoted API syntax in a condition."""

    _, condition = requirement_parts(value)
    if condition is None:
        return None
    if "`" in condition:
        return "Runtime conditions use straight single quotes, not Markdown backticks"
    if "-'" in condition or re.search(r"'--[^']+'\s+[A-Za-z0-9_.-]+\s+is\b", condition):
        return "quote one complete API expression with straight single quotes"
    ranges = quoted_ranges(condition)
    for pattern in (OPTION_TOKEN, SNAKE_API_TOKEN):
        for match in pattern.finditer(condition):
            if not position_is_quoted(match.start(), ranges):
                return f"quote API expression '{match.group(0)}' with straight single quotes"
    known = set(KNOWN_CALLABLES) | set(CHECKSUM_PROVIDERS) | set(callables)
    for match in re.finditer(r"\b[A-Za-z][A-Za-z0-9_.-]*\b", condition):
        if match.group(0) in known and not position_is_quoted(match.start(), ranges):
            return f"quote API expression '{match.group(0)}' with straight single quotes"
    return None


def runtime_policy_messages(
    value: str,
    callables: set[str],
) -> list[tuple[str, str]]:
    """Return generic structured-prose diagnostics for one Runtime entry."""

    messages: list[tuple[str, str]] = []
    python_message = python_minimum_message(value)
    if python_message:
        messages.append((RULE_PYTHON_MINIMUM, python_message))
    api_message = inline_api_message(value, callables)
    if api_message:
        messages.append((RULE_INLINE_API, api_message))
    if value.endswith("."):
        messages.append((RULE_PUNCTUATION, "Runtime entries do not end with a period"))
    if CONDA_ENVIRONMENT.search(value):
        messages.append(
            (
                RULE_CONDA_ENVIRONMENT,
                "environment resources use 'A compatible Conda environment' wording",
            )
        )
    elif value.startswith("a compatible Conda environment"):
        messages.append(
            (
                RULE_CONDA_ENVIRONMENT,
                "Conda environment resources begin with capitalized 'A compatible'",
            )
        )
    if REDUNDANT_COMPOSITE.search(value):
        messages.append(
            (
                RULE_REDUNDANT_COMPOSITE,
                "separate the Conda environment, activation executable, and tool entries",
            )
        )
    if (
        value
        and value[0].islower()
        and not is_executable_requirement(value, callables)
    ):
        messages.append(
            (
                RULE_CAPITALIZATION,
                "prose/resource requirements begin with sentence-style capitalization",
            )
        )
    return messages


def source_without_help(path: Path, owner: str | None = None) -> list[tuple[int, str]]:
    """Return source rows after masking only recognized rendered help heredocs."""

    text = path.read_text(encoding="utf-8")
    masked = set()
    for heredoc in extract_help_heredocs(text):
        masked.update(range(heredoc.start_line, heredoc.end_line + 1))
    rows = [(number, line) for number, line in enumerate(text.splitlines(), 1) if number not in masked]
    if owner is None:
        return rows
    for index, (_, line) in enumerate(rows):
        match = FUNCTION_START.match(line)
        if match is None or match.group("name") != owner:
            continue
        depth = line.count("{") - line.count("}")
        end = index + 1
        while end < len(rows) and depth > 0:
            depth += rows[end][1].count("{") - rows[end][1].count("}")
            end += 1
        return rows[index:end]
    return []


def evidence_for_owner(identity: str, source: Path, kind: str, root: Path) -> list[Evidence]:
    """Collect literal callable evidence, retaining ambiguity as review data."""

    evidence = [Evidence(identity, str(source.relative_to(root)), 1, "shell", "bash", "bash", "inherited", "always", "", True, "settled")]
    owner = identity.rsplit("::", 1)[-1] if kind == "function" else None
    rows = source_without_help(source, owner)
    functions = {
        match.group("name")
        for _, line in source_without_help(source)
        if (match := FUNCTION_START.match(line)) is not None
    }
    for number, line in rows:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        for match in COMMAND_V.finditer(stripped):
            name = match.group("name")
            if name in CHECKSUM_PROVIDERS and source.name == "install_atria.sh":
                evidence.append(Evidence(identity, str(source.relative_to(root)), number, "callable_alternative", "sha256sum_or_shasum", f"command -v {name}", "direct", "conditional", "when --dry_run is not specified", True, "settled"))
            elif identity == "lib/bash/workflows/process_region.sh::check_region_bdg" and name in DECOMPRESSION_PROVIDERS:
                evidence.append(Evidence(identity, str(source.relative_to(root)), number, "callable_alternative", "gzip_or_zcat", f"command -v {name}", "direct", "conditional", "when processing a gzipped bedGraph", True, "settled"))
            elif name in {"julia", "atria"} and source.name == "install_atria.sh":
                evidence.append(Evidence(identity, str(source.relative_to(root)), number, "optional_reuse_target", name, f"command -v {name}", "direct", "conditional", "when --if_exists reuse is specified", False, "settled"))
            else:
                category = "callable" if name in KNOWN_CALLABLES else "unknown_callable"
                evidence.append(Evidence(identity, str(source.relative_to(root)), number, category, name, f"command -v {name}", "direct", "conditional", "available in PATH", True, "high" if category == "callable" else "low"))
        command = COMMAND.search(stripped)
        if command is None:
            continue
        name = command.group("name")
        if name in SHELL_WORDS or name in {"command"} or name in functions:
            continue
        if name in KNOWN_CALLABLES:
            evidence.append(Evidence(identity, str(source.relative_to(root)), number, "callable", name, name, "direct", "always", "", True, "medium"))
        if stripped.startswith(("source ", ". ")):
            try:
                record = classify_source(str(source.relative_to(root)), stripped)
            except ValueError:
                continue
            if record.classification != "static":
                evidence.append(Evidence(identity, str(source.relative_to(root)), number, "source", record.classification, stripped, "direct", "unclassified", "", False, "low"))
    return evidence


def candidate(owner: str, evidence: Evidence, documented: list[str]) -> dict[str, object]:
    """Render a stable unresolved semantic candidate."""

    payload = {"owner": owner, "evidence": evidence.as_dict(), "documented": documented}
    signature = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return {
        "signature": f"sha256:{signature}", "rule_id": "RUNTIME.SEMANTIC.UNKNOWN_CALLABLE",
        "owner": owner, "uncertainty_class": "literal-command-classification",
        "evidence": evidence.as_dict(), "current_text": documented,
        "disposition_options": ["document as callable", "classify as shell/local syntax", "exclude as internal"],
        "preferred_disposition": "classify before asserting a public requirement",
    }


def analyze_owner(
    identity: str,
    unit: tuple[str, Path, Path, Heredoc],
    root: Path,
    callables: set[str] | None = None,
    concepts: dict[str, set[str]] | None = None,
) -> tuple[dict[str, object], list[Finding], list[dict[str, object]]]:
    """Analyze one discovered owner and return JSON-ready inventory data."""

    kind, source, documentation, heredoc = unit
    label, entries_with_lines, grammar = runtime_entries(heredoc)
    entries = [text for _, text in entries_with_lines]
    evidence = evidence_for_owner(identity, source, kind, root)
    expected_requirements = settled_requirements(identity)
    findings: list[Finding] = []
    path = str(documentation.relative_to(root))
    line = label or heredoc.start_line
    if grammar:
        findings.extend(Finding(RULE_GRAMMAR, identity, path, line, message) for message in grammar)
    if label is None:
        findings.append(Finding(RULE_OWNER_SECTION, identity, path, heredoc.start_line, "current owners inherit bash >= 4.4 and require a Runtime requirements section"))
    else:
        markers = [text.startswith("- ") for text in entries]
        if len(entries) == 1 and markers and entries[0].startswith("- "):
            findings.append(Finding(RULE_CARDINALITY, identity, path, line, "one requirement must be unbulleted"))
        if len(entries) > 1 and not all(markers):
            findings.append(Finding(RULE_CARDINALITY, identity, path, line, "two or more requirements must use '-' bullets"))
        normalized = [text[2:] if text.startswith("- ") else text for text in entries]
        expected_order = normalize_runtime_entries(normalized)
        if len(normalized) > 1 and normalized != expected_order:
            findings.append(
                Finding(
                    RULE_ORDER,
                    identity,
                    path,
                    line,
                    "multi-item Runtime requirements must be sorted case-insensitively by complete displayed text",
                    tuple(expected_order),
                    tuple(normalized),
                )
            )
        if f"bash >= {MINIMUM_BASH[0]}.{MINIMUM_BASH[1]}" not in normalized:
            findings.append(Finding(RULE_MISSING, identity, path, line, "missing inherited requirement 'bash >= 4.4'", ("bash >= 4.4",), tuple(normalized)))
        for requirement in expected_requirements:
            if requirement not in normalized:
                findings.append(Finding(RULE_MISSING, identity, path, line, f"missing source-reviewed requirement '{requirement}'", (requirement,), tuple(normalized)))
        if any(item.key == "sha256sum_or_shasum" for item in evidence) and CHECKSUM_REQUIREMENT not in normalized:
            findings.append(Finding(RULE_MISSING, identity, path, line, "missing alternative SHA-256 provider requirement", (CHECKSUM_REQUIREMENT,), tuple(normalized)))
        if any(item.key == "gzip_or_zcat" for item in evidence) and DECOMPRESSION_REQUIREMENT not in normalized:
            findings.append(Finding(RULE_MISSING, identity, path, line, "missing alternative bedGraph decompression provider requirement", (DECOMPRESSION_REQUIREMENT,), tuple(normalized)))
        for entry_line, entry in entries_with_lines:
            value = entry[2:] if entry.startswith("- ") else entry
            callable_message = exact_callable_message(
                value,
                callables or set(),
                concepts or {},
            )
            if callable_message is not None:
                findings.append(
                    Finding(
                        RULE_CALLABLE,
                        identity,
                        path,
                        entry_line,
                        callable_message,
                    )
                )
            for rule_id, message in runtime_policy_messages(
                value,
                callables or set(),
            ):
                findings.append(
                    Finding(rule_id, identity, path, entry_line, message)
                )
            requirement_head, _ = requirement_parts(value)
            parsed = VERSION.fullmatch(requirement_head)
            if parsed and parsed.group("name") not in {"bash", "python", "python3"}:
                findings.append(Finding(RULE_VERSION, identity, path, entry_line, "non-Bash/Python version claims need source-grounded validation"))
            if parsed and parsed.group("name") == "bash" and requirement_head != "bash >= 4.4":
                findings.append(Finding(RULE_VERSION, identity, path, entry_line, "Bash must use the canonical minimum 'bash >= 4.4'"))
            if "env_protocol" in value:
                findings.append(Finding(RULE_ENV_DEFAULT, identity, path, entry_line, "env_protocol is an --env_nam default, not a Runtime requirement"))
            if INTERNAL.search(value):
                findings.append(Finding(RULE_INTERNAL, identity, path, entry_line, "repository-local helpers, scripts, and modules are not Runtime requirements"))
            condition = CONDITION.fullmatch(value)
            if " when " in value and condition is None:
                findings.append(Finding(RULE_CONDITION, identity, path, entry_line, "conditions must use '<requirement> (when <source-grounded trigger>)'"))
        if identity in SETTLED_OWNER_REQUIREMENTS:
            for entry_line, entry in entries_with_lines:
                value = entry[2:] if entry.startswith("- ") else entry
                if value not in expected_requirements:
                    findings.append(Finding(RULE_OVERCLAIM, identity, path, entry_line, "entry is not part of the source-reviewed owner contract", expected_requirements, tuple(normalized)))
    candidates = [candidate(identity, item, entries) for item in evidence if item.category == "unknown_callable"]
    public_aliases = sorted(
        {
            alias
            for chunk in alias_chunks(source.read_text(encoding="utf-8"))
            for alias in chunk.aliases
            if alias.startswith("--") and alias != "--hlp"
        }
    )
    inventory = {
        "identity": identity, "kind": kind, "source_path": str(source.relative_to(root)),
        "documentation_path": path, "help_start_line": heredoc.start_line,
        "runtime_label_line": label, "documented_entries": entries,
        "public_aliases": public_aliases,
        "evidence": [item.as_dict() for item in evidence],
    }
    return inventory, findings, candidates


def audit(root: Path) -> tuple[list[dict[str, object]], list[Finding], list[dict[str, object]]]:
    """Audit all current Examples-owned Runtime-requirements documents."""

    inventories: list[dict[str, object]] = []
    findings: list[Finding] = []
    candidates: list[dict[str, object]] = []
    callables, concepts = load_registry(root / "dev/config/command_names.json")
    for identity, unit in owner_units(root).items():
        inventory, owner_findings, owner_candidates = analyze_owner(
            identity,
            unit,
            root,
            callables,
            concepts,
        )
        inventories.append(inventory)
        findings.extend(owner_findings)
        candidates.extend(owner_candidates)
    return inventories, findings, candidates


def write_json(path: Path, value: object) -> None:
    """Write deterministic, newline-terminated JSON."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--inventory-output", type=Path)
    parser.add_argument("--evidence-output", type=Path)
    parser.add_argument("--semantic-output", type=Path)
    parser.add_argument(
        "--fix",
        action="store_true",
        help="normalize source-reviewed Runtime blocks before auditing",
    )
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.root.resolve()
    if args.fix:
        changed_files, changed_blocks = rewrite_repository(root)
        print(
            "Runtime normalization: "
            f"{changed_files} file(s), {changed_blocks} block(s) changed"
        )
    inventory, findings, candidates = audit(root)
    if args.inventory_output:
        write_json(args.inventory_output, {"owners": inventory, "findings": [item.as_dict() for item in findings]})
    if args.evidence_output:
        write_json(args.evidence_output, {"owners": [{"identity": item["identity"], "evidence": item["evidence"]} for item in inventory]})
    if args.semantic_output:
        write_json(args.semantic_output, {"candidates": candidates})
    for finding in findings:
        print(finding.format())
    print(f"Runtime requirements: {len(inventory)} owner(s), {len(findings)} finding(s), {len(candidates)} semantic candidate(s)")
    return 1 if args.strict and (findings or candidates) else 0


if __name__ == "__main__":
    sys.exit(main())
