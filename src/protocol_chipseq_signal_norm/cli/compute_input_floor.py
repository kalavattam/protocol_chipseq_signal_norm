#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: compute_input_floor.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
#   GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


"""
Compute denominator floors for ratio-based ChIP-seq signal workflows.

The CLI accepts input, mode, genome-size, bin-size, distribution-filtering, and
formatting options. It prints one rounded floating-point floor recommendation
to stdout.

Examples
--------
python -m protocol_chipseq_signal_norm.cli.compute_input_floor \\
    [--mode <mode>] [options]
"""

from __future__ import annotations

import argparse
import math

# Import pysam lazily in `_count_alignment_records()`.
import signal
import sys
from contextlib import redirect_stdout, suppress

from protocol_chipseq_signal_norm.utilities.utils_check import (
    check_exists,
    validate_comparison,
)
from protocol_chipseq_signal_norm.utilities.utils_cli import (
    CapArgumentParser,
    _HelpExample,
    _SectionedHelpConfig,
    add_help_cap,
)
from protocol_chipseq_signal_norm.utilities.utils_format import format_value
from protocol_chipseq_signal_norm.utilities.utils_io import (
    DEF_SKP_PFX,
    is_header,
    open_in,
    parse_skp_pfx,
)
from protocol_chipseq_signal_norm.utilities.utils_stabilizer import (
    determine_coef_eff,
    iter_vals_bdg,
    pick_stabilizer,
)

with suppress(AttributeError, ValueError):
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

assert sys.version_info >= (3, 11), "Python >= 3.11 required."

# These default SAM allowlists select paired-end and single-end alignments.
PAIRED_FLAGS = {99, 1123, 163, 1187}
SINGLE_FLAGS = {0, 16, 1024, 1040}

_CANONICAL_INPUT_FORMATS = ("bam", "cram", "bed", "bedgraph")
_FORMAT_HINT_ALIASES = {
    "bam": "bam",
    "cram": "cram",
    "bed": "bed",
    "bedgraph": "bedgraph",
    "bdg": "bedgraph",
    "bg": "bedgraph",
}


class InputFloorError(Exception):
    """
    Represent an anticipated input-floor failure.
    """


class FlagParseError(InputFloorError):
    """
    Represent an invalid SAM FLAG list.
    """


class AlignmentReadError(InputFloorError):
    """
    Represent an unreadable alignment input.
    """


class InputFloorValidationError(InputFloorError):
    """
    Represent an invalid input-floor computation request.
    """


def _note_ignored(option: str, reason: str, value: object) -> None:
    """
    Print a standard “ignored option” note when a mode does not use an arg.

    Parameters
    ----------
    option : str
        Option name as shown to the user (e.g., '--fil_in').
    reason : str
        Short reason string (e.g., \"in '--mode norm'\").
    value : object
        Value that triggered the note (printed only when useful).
    """

    if value is None:
        return

    print(
        f"Note: '{option}' is ignored {reason}: got '{value}'.",
        file=sys.stderr,
    )


def parse_flag_csv(csv_text: str, label: str) -> set[int]:
    """
    Parse a comma-separated string of SAM FLAGs into a set of ints.

    Decimal values and hexadecimal values with a '0x' prefix are accepted.

    Parameters
    ----------
    csv_text : str
        Comma-separated FLAG values. Whitespace around commas is ignored, and
        empty tokens are skipped.
    label : str
        Short option name used in error messages.

    Returns
    -------
    flags : set[int]
        Non-empty set of integer FLAG values.

    Raises
    ------
    FlagParseError
        If a token is invalid, outside the unsigned 16-bit range, or no FLAG
        values are present.
    """

    flags: set[int] = set()

    for token in csv_text.split(","):
        stripped = token.strip()

        if not stripped:
            continue

        try:
            integer_base = 16 if stripped.lower().startswith("0x") else 10
            flag = int(stripped, integer_base)

            if not (0 <= flag <= 0xFFFF):
                raise FlagParseError(
                    f"Error: FLAG out of range (0..65535): {stripped}",
                )

            flags.add(flag)
        except ValueError as error:
            raise FlagParseError(
                f"Error: Invalid FLAG '{stripped}' in '{label}'. Use decimal "
                "(e.g., 99) or hex (e.g., 0x63).",
            ) from error

    if not flags:
        raise FlagParseError(
            f"Error: No valid FLAGs parsed for '{label}'.",
        )

    return flags


def _count_alignment_records(
    alignment_path: str,
    paired_flags: set[int] | None = None,
    single_flags: set[int] | None = None,
    *,
    alignment_format: str,
    ref_fa: str | None = None,
) -> int:
    """
    Count BAM or CRAM alignment records selected by SAM FLAG.

    Duplicate-flagged reads remain eligible because rDNA copy number and
    high-coverage repeats can represent biological signal rather than only PCR
    or optical artifacts. Paired-end records serve as fragment proxies.

    Parameters
    ----------
    alignment_path : str
        Alignment path, or '-' for standard input. Streaming does not require
        an index.
    paired_flags : set[int] | None
        Paired-end FLAG allowlist. The default is 'PAIRED_FLAGS'.
    single_flags : set[int] | None
        Single-end FLAG allowlist. The default is 'SINGLE_FLAGS'.
    alignment_format : str
        Resolved alignment format, either 'bam' or 'cram'.
    ref_fa : str | None
        Reference FASTA used for CRAM decoding.

    Returns
    -------
    count : int
        Number of records whose FLAG occurs in either allowlist.

    Raises
    ------
    AlignmentReadError
        If 'pysam' is unavailable or the input cannot be opened or read.
    """

    try:
        import pysam
    except ImportError as error:
        raise AlignmentReadError(
            "Error: 'pysam' is required for alignment input.",
        ) from error

    if alignment_format not in {"bam", "cram"}:
        raise AlignmentReadError(
            "Error: Alignment format must be 'bam' or 'cram' "
            f"(got '{alignment_format}').",
        )

    if alignment_format == "cram" and not ref_fa:
        raise AlignmentReadError(
            "Error: CRAM input requires a reference FASTA.",
        )

    if paired_flags is None:
        paired_flags = PAIRED_FLAGS

    if single_flags is None:
        single_flags = SINGLE_FLAGS

    n_in = 0
    read_mode = "rc" if alignment_format == "cram" else "rb"
    alignment_source = (
        sys.stdin.buffer if alignment_path == "-" else alignment_path
    )
    alignment_options = {}

    if alignment_format == "cram":
        alignment_options["reference_filename"] = ref_fa

    try:
        alignment_handle = pysam.AlignmentFile(
            alignment_source,
            read_mode,
            **alignment_options,
        )

        with alignment_handle as alignment_stream:
            for read in alignment_stream.fetch(until_eof=True):
                if read.flag in paired_flags or read.flag in single_flags:
                    n_in += 1
    except (FileNotFoundError, OSError, ValueError) as error:
        raise AlignmentReadError(
            f"Error: Cannot process {alignment_format} alignment "
            f"'{alignment_path}': {error}",
        ) from error

    return n_in


def _count_bed_records(bed_path: str, skp_pfx: tuple[str, ...]) -> int:
    """
    Count data records in a BED input.

    Blank rows and configured header or metadata rows are skipped.

    Parameters
    ----------
    bed_path : str
        BED path, which may be gzip-compressed.
    skp_pfx : tuple[str, ...]
        Header and metadata prefixes to skip.

    Returns
    -------
    count : int
        Number of data records.

    Raises
    ------
    AlignmentReadError
        If the input cannot be opened or read.
    """

    try:
        n = 0

        with open_in(bed_path) as f:
            for line in f:
                if is_header(line, skp_pfx) or not line.strip():
                    continue

                n += 1

        return n
    except (FileNotFoundError, PermissionError) as error:
        raise AlignmentReadError(
            f"Error: Cannot open bed file '{bed_path}': {error}",
        ) from error
    except Exception as error:
        raise AlignmentReadError(
            f"Error: Unexpected error with '{bed_path}': {error}",
        ) from error


def _canonicalize_input_format_hint(hint: str | None) -> str | None:
    """
    Canonicalize an explicit project-owned input-format hint.

    Parameters
    ----------
    hint : str | None
        Explicit format hint in any letter case.

    Returns
    -------
    canonical_hint : str | None
        One of 'bam', 'cram', 'bed', or 'bedgraph'. An unknown hint is
        case-normalized for the parser or caller to reject.
    """

    if hint is None:
        return None

    lowercase_hint = hint.casefold()

    return _FORMAT_HINT_ALIASES.get(lowercase_hint, lowercase_hint)


def infer_input_format(path: str, hint: str | None = None) -> str:
    """
    Infer file format label from the path suffix.

    Parameters
    ----------
    path : str
        Input file path, or '-' for standard input.
    hint : str | None
        Optional case-insensitive format hint used only when 'path == "-"' is
        true. Accepts 'bam', 'cram', 'bed', 'bedGraph', 'bdg', or 'bg'.

    Returns
    -------
    format_name : str
        Canonical 'bam', 'cram', 'bed', or 'bedgraph' for recognized input;
        otherwise, 'other'.
    """

    if path == "-":
        canonical_hint = (
            hint
            if hint in _CANONICAL_INPUT_FORMATS
            else _canonicalize_input_format_hint(hint)
        )

        if canonical_hint in _CANONICAL_INPUT_FORMATS:
            return canonical_hint

    lowercase_path = path.lower()
    if lowercase_path.endswith(".bam"):
        return "bam"

    if lowercase_path.endswith(".cram"):
        return "cram"

    if lowercase_path.endswith(".bed") or lowercase_path.endswith(".bed.gz"):
        return "bed"

    if (
        lowercase_path.endswith(".bedgraph")
        or lowercase_path.endswith(".bedgraph.gz")
        or lowercase_path.endswith(".bdg")
        or lowercase_path.endswith(".bdg.gz")
        or lowercase_path.endswith(".bg")
        or lowercase_path.endswith(".bg.gz")
    ):
        return "bedgraph"

    return "other"


def _validate_positive_mode_dimensions(
    mode: str,
    siz_bin: int,
    siz_gen: int,
) -> None:
    """
    Validate positive dimensions for fragment and normalized calculations.

    Parameters
    ----------
    mode : str
        Floor-computation mode.
    siz_bin : int
        Target signal-bin width in base pairs.
    siz_gen : int
        Effective genome size in base pairs.

    Raises
    ------
    InputFloorValidationError
        If a dimension used by 'frag' or 'norm' is nonpositive.
    """

    if mode not in {"frag", "norm"}:
        return

    if siz_bin <= 0:
        raise InputFloorValidationError(
            f"Error: 'siz_bin' must be positive (got siz_bin={siz_bin}).",
        )

    if siz_gen <= 0:
        raise InputFloorValidationError(
            f"Error: 'siz_gen' must be positive (got siz_gen={siz_gen}).",
        )


def _validate_mode_dimensions(
    mode: str,
    siz_bin: int,
    siz_gen: int,
) -> None:
    """
    Validate dimensions used by fragment and normalized floor calculations.

    Parameters
    ----------
    mode : str
        Floor-computation mode.
    siz_bin : int
        Target signal-bin width in base pairs.
    siz_gen : int
        Effective genome size in base pairs.

    Raises
    ------
    InputFloorValidationError
        If a dimension used by 'frag' or 'norm' is nonpositive, equal, or
        reversed.
    """

    if mode not in {"frag", "norm"}:
        return

    _validate_positive_mode_dimensions(mode, siz_bin, siz_gen)

    if siz_bin >= siz_gen:
        raise InputFloorValidationError(
            "Error: 'siz_bin' must be smaller than 'siz_gen' (got "
            f"siz_bin={siz_bin}, siz_gen={siz_gen}).",
        )


def _validate_cram_reference(
    format_name: str,
    ref_fa: str | None,
) -> None:
    """
    Require one existing reference FASTA for resolved CRAM input.

    Parameters
    ----------
    format_name : str
        Canonical input format.
    ref_fa : str | None
        Reference FASTA path supplied for CRAM decoding.

    Raises
    ------
    InputFloorValidationError
        If CRAM input has no readable reference FASTA.
    """

    if format_name != "cram":
        return

    if not ref_fa:
        raise InputFloorValidationError(
            "Error: CRAM input requires '--ref_fa' for reference-backed "
            "decoding.",
        )

    try:
        check_exists(ref_fa, kind="file", label="Reference FASTA")
    except FileNotFoundError as error:
        raise InputFloorValidationError(f"Error: {error}") from None


def compute_input_floor(
    fil_in: str,
    siz_bin: int,
    siz_gen: int,
    mode: str = "dist",
    method: str = "qntl_nz",
    qntl_nz: float = 1.0,
    coef: float | None = None,
    floor: float = 0.0,
    eps: float = 0.0,
    mode_nz: str = "closed",
    paired_flags: set[int] | None = None,
    single_flags: set[int] | None = None,
    skp_pfx: tuple[str, ...] = DEF_SKP_PFX,
    infmt: str | None = None,
    ref_fa: str | None = None,
) -> float:
    """
    Compute a denominator floor for input normalization.

    The 'mode' parameter selects 'dist', 'frag', or 'norm'. Use 'dist' for new
    analyses. Use 'frag' or 'norm' when reproducing the fragment-normalization
    or normalized-coverage floor calculations used in the Dickson/siQ-ChIP and
    *Bio-protocol* workflows.

    Parameters
    ----------
    fil_in : str
        Input path. 'dist' accepts 'bedGraph', 'bdg', or 'bg', optionally with
        '.gz'; 'frag' accepts 'bam', 'cram', 'bed', or 'bed.gz'; 'norm' ignores
        'fil_in'. For 'dist' and 'frag', '-' reads standard input and requires
        'infmt'.
    siz_bin : int
        Target signal-bin width in base pairs. Used only by 'frag' and 'norm';
        callers must supply a positive value smaller than 'siz_gen'.
    siz_gen : int
        Effective genome size in base pairs. Used only by 'frag' and 'norm';
        callers must supply a positive value.
    mode : str
        One of 'dist', 'frag', or 'norm'. 'dist' reads bedGraph column four,
        filters finite values, applies 'eps' and 'mode_nz', retains 'v_i > 0',
        applies 'method', then applies 'dep_min := max(dep_min, floor)' without
        using or validating 'siz_bin' or 'siz_gen'. 'frag' returns
        '((n * b) / g) / (1 - (b / g))', where 'n' is the counted
        alignment-record total. 'norm' returns '(b / g) / (1 - (b / g))'. Here,
        'b = siz_bin' and 'g = siz_gen'.
    method : str
        Distribution rule used only by 'dist'. 'qntl_nz' selects the quantile;
        'frc_mdn_nz' returns 'coef * median(v_i)'; 'frc_avg_nz' returns
        'coef * mean(v_i)'; and 'min_nz' returns 'coef * min(v_i)' over the
        filtered positive values.
    qntl_nz : float
        Quantile percentage in '[0, 100]', used only by 'dist/qntl_nz'. For
        sorted filtered values and 'q = qntl_nz / 100', select
        'i = floor(q * (N - 1))', clamped to '[0, N - 1]', then
        'dep_min = sorted_vals[i]'.
    coef : float | None
        Nonnegative coefficient for the three coefficient-based 'dist' methods.
        'None' matches 'compute_pseudo.py' defaults: '0.01' for 'frc_*' and
        '1.0' for 'min_nz'. Ignored by 'qntl_nz'.
    floor : float
        Nonnegative lower bound applied after the raw 'dist' statistic.
    eps : float
        Nonnegative zero tolerance used only by 'dist'.
    mode_nz : str
        Epsilon rule used only by 'dist': 'closed' drops 'abs(v_i) <= eps';
        'open' drops 'abs(v_i) < eps'; and 'off' disables epsilon filtering.
        Positive-only filtering follows.
    paired_flags : set[int] | None
        Optional paired-end main-alignment FLAG allowlist for BAM or CRAM
        counting. 'None' uses the defaults.
    single_flags : set[int] | None
        Optional single-end main-alignment FLAG allowlist for BAM or CRAM
        counting. 'None' uses the defaults.
    skp_pfx : tuple[str, ...]
        Prefixes skipped as header or metadata rows in BED and bedGraph inputs.
    infmt : str | None
        Required case-insensitive format hint when 'fil_in' is '-'. Accepts
        'bam', 'cram', 'bed', 'bedGraph', 'bdg', or 'bg' and resolves to a
        canonical format; ignored for named paths.
    ref_fa : str | None
        Reference FASTA required for CRAM decoding and ignored otherwise.

    Returns
    -------
    dep_min : float
        Unrounded denominator floor. Callers apply it as
        'denominator := max(denominator, dep_min)'.

    Raises
    ------
    InputFloorValidationError
        If the mode, 'frag'/'norm' dimensions, input format, filtered values,
        or computed floor are invalid.
    AlignmentReadError
        If an alignment input cannot be read.
    """

    if mode not in {"dist", "frag", "norm"}:
        raise InputFloorValidationError(
            f"Error: Invalid mode '{mode}'. Use 'dist', 'frag', or 'norm'.",
        )

    _validate_mode_dimensions(mode, siz_bin, siz_gen)

    if mode == "dist":
        format_name = infer_input_format(fil_in, infmt)
        if format_name != "bedgraph":
            raise InputFloorValidationError(
                "Error: '--mode dist' requires a bedGraph-like input file "
                "(bedGraph, bdg, or bg, optionally with .gz).",
            )

        # Ratio denominators occupy a positive domain, so the iterator policy
        # and explicit guard retain only positive values for the floor.
        vals = list(
            iter_vals_bdg(
                fil_in,
                eps=eps,
                mode_nz=mode_nz,
                skp_pfx=skp_pfx,
                nz_policy="pos",
            ),
        )

        vals = [v for v in vals if v > 0.0]

        if not vals:
            raise InputFloorValidationError(
                "Error: No positive finite values found after filtering "
                "(cannot compute distribution-based 'dep_min').",
            )

        coef_eff = determine_coef_eff(method, coef)

        dep_min = pick_stabilizer(
            vals,
            method=method,
            coef=coef_eff,
            qntl_pct=qntl_nz,
            floor=floor,
            qntl_rule="floor",
        )

        if not math.isfinite(dep_min):
            raise InputFloorValidationError(
                "Error: 'dep_min' is non-finite; check input values and "
                "filters.",
            )

        return dep_min

    b_over_g = siz_bin / siz_gen

    # Normalized mode does not depend on the alignment count.
    if mode == "norm":
        dep_min = b_over_g / (1.0 - b_over_g)
        return dep_min

    format_name = infer_input_format(fil_in, infmt)

    if format_name in {"bam", "cram"}:
        _validate_cram_reference(format_name, ref_fa)
        n_in = _count_alignment_records(
            fil_in,
            paired_flags=paired_flags,
            single_flags=single_flags,
            alignment_format=format_name,
            ref_fa=ref_fa,
        )
    elif format_name == "bed":
        n_in = _count_bed_records(fil_in, skp_pfx)
    elif format_name == "bedgraph":
        raise InputFloorValidationError(
            "Error: '--mode frag' expects bam, cram, or bed/bed.gz "
            "(alignment records), not bedGraph.",
        )
    else:
        raise InputFloorValidationError(
            f"Error: Unsupported file type: {fil_in}. Provide a bam, cram, "
            "or bed/bed.gz file.",
        )

    dep_min = ((n_in * siz_bin) / siz_gen) / (1.0 - b_over_g)
    return dep_min


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    arguments : argparse.Namespace
        Parsed alignment, binning, and input-floor options.

    Raises
    ------
    SystemExit
        If argument parsing fails or help is requested.
    """

    parser = CapArgumentParser(
        description=(
            "Compute 'dep_min', a floor for positive input denominators that "
            "helps prevent extreme or erroneous IP/input ratios. The current "
            "repository consumer, 'compute_signal_ratio', applies it as "
            "'denominator := max(denominator, dep_min)'. Use 'dist' for new "
            "analyses. Use 'frag' or 'norm' when reproducing the "
            "fragment-normalization or normalized-coverage floor calculations "
            "used in the Dickson/siQ-ChIP and *Bio-protocol* workflows."
        ),
        prog="compute_input_floor",
        _sectioned_help=_SectionedHelpConfig(
            usage_rows=(
                ("help", "verbose"),
                ("mode",),
                ("fil_in", "infmt", "ref_fa", "skp_pfx"),
                ("method", "qntl_nz", "coef", "eps", "mode_nz", "floor"),
                ("siz_bin", "siz_gen", "flags_pe", "flags_se"),
                ("dp",),
            ),
            examples=(
                _HelpExample(
                    description=(
                        "Compute a first-percentile floor from bedGraph "
                        "values."
                    ),
                    command_lines=(
                        "compute_input_floor",
                        "--mode dist",
                        "--fil_in signal.bdg",
                        "--method qntl_nz",
                        "--qntl_nz 1",
                    ),
                ),
                _HelpExample(
                    description=(
                        "Compute a normalized floor from explicit dimensions."
                    ),
                    command_lines=(
                        "compute_input_floor --mode norm --siz_bin 30 "
                        "--siz_gen 12157105",
                    ),
                ),
            ),
        ),
    )
    add_help_cap(parser)
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Run script in verbose mode.",
    )

    parser.add_argument(
        "-md",
        "--mode",
        dest="mode",
        choices=["dist", "frag", "norm"],
        default="dist",
        help=(
            "Floor-computation mode (default: '%(default)s').\n"
            "- dist: Recommended for new analyses. Read bedGraph column four, "
            "skip non-data rows, filter to finite values, apply "
            "'--eps'/'--mode_nz', retain 'v_i > 0' because ratio denominators "
            "occupy a positive domain, and summarize with '--method'. This "
            "mode does not use, validate, compare, infer, or warn about "
            "'--siz_bin' or '--siz_gen'. See the method-specific options for "
            "calculations and bounds.\n"
            "- frag: Reproduce the fragment-normalization floor calculation "
            "used in the Dickson/siQ-ChIP and *Bio-protocol* workflows. Count "
            "BAM, CRAM, or BED alignment records and compute 'dep_min = ((n * "
            "b) / g) / [1 - (b / g)]'. Here, 'n' is the counted-record total, "
            "'b = siz_bin', and 'g = siz_gen'.\n"
            "- norm: Reproduce the normalized-coverage floor calculation used "
            "in the Dickson/siQ-ChIP and *Bio-protocol* workflows. Compute "
            "'dep_min = (b / g) / [1 - (b / g)]' from 'b = siz_bin' and 'g = "
            "siz_gen'; '--fil_in' is ignored. The command returns one scalar "
            "'dep_min'."
        ),
    )

    parser.add_argument(
        "-fi",
        "--fil_in",
        dest="fil_in",
        type=str,
        required=False,
        default=None,
        help=(
            "Input path. For 'dist', provide bedGraph, bdg, or bg, optionally "
            "with '.gz'. For 'frag', provide BAM, CRAM, or BED/BED.GZ. For "
            "'dist' and 'frag', '-' reads stdin and requires '--infmt'. "
            "'norm' ignores '--fil_in'."
        ),
    )
    parser.add_argument(
        "-if",
        "--infmt",
        dest="infmt",
        type=_canonicalize_input_format_hint,
        choices=_CANONICAL_INPUT_FORMATS,
        metavar="{bam,cram,bed,bedGraph,bdg,bg}",
        default=None,
        help=(
            "Case-insensitive input-format hint for 'dist' or 'frag'. "
            "Required when '--fil_in -' reads stdin and ignored for named "
            "paths. Choose 'bam', 'cram', 'bed', 'bedGraph', 'bdg', or 'bg'; "
            "accepted values resolve to 'bam', 'cram', 'bed', or 'bedgraph'."
        ),
    )
    parser.add_argument(
        "-rf",
        "--ref_fa",
        dest="ref_fa",
        default=None,
        help=(
            "Reference FASTA file required for CRAM decoding. Ignored for "
            "other input formats."
        ),
    )

    parser.add_argument(
        "-sp",
        "--skp_pfx",
        dest="skp_pfx",
        type=str,
        default=",".join(DEF_SKP_PFX),
        help=(
            "Comma-separated header prefixes skipped in BED/BED.GZ and "
            "bedGraph-like input; an empty string disables skipping (default: "
            "'%(default)s')."
        ),
    )

    parser.add_argument(
        "-m",
        "--method",
        dest="method",
        choices=("qntl_nz", "frc_mdn_nz", "frc_avg_nz", "min_nz"),
        default="qntl_nz",
        help=(
            "Distribution method used only in '--mode dist' (default: "
            "'%(default)s').\n"
            "- qntl_nz: 'dep_min = Q_q({v_i : v_i > 0})'. See '--qntl_nz' for "
            "the selection rule.\n"
            "- frc_mdn_nz: 'dep_min = coef * median({v_i : v_i > 0})'.\n"
            "- frc_avg_nz: 'dep_min = coef * mean({v_i : v_i > 0})'.\n"
            "- min_nz: 'dep_min = coef * min({v_i : v_i > 0})'.\n"
            "After this statistic, '--floor' applies the lower bound."
        ),
    )
    parser.add_argument(
        "-qn",
        "--qntl_nz",
        dest="qntl_nz",
        type=float,
        default=1.0,
        help=(
            "Quantile percentage in '[0, 100]' used only by '--mode dist "
            "--method qntl_nz'. For sorted filtered values and 'q = qntl_nz / "
            "100', select 'i = floor(q * (N - 1))', clamped to '[0, N - 1]' "
            "and then 'dep_min = sorted_vals[i]' (default: %(default)s)."
        ),
    )
    parser.add_argument(
        "-c",
        "--coef",
        dest="coef",
        type=float,
        default=None,
        help=(
            "Nonnegative coefficient used only by the 'frc_mdn_nz', "
            "'frc_avg_nz', and 'min_nz' distribution methods. If omitted, "
            "defaults match 'compute_pseudo.py': 0.01 for 'frc_*' and 1.0 for "
            "'min_nz'. The 'qntl_nz' method ignores it."
        ),
    )
    parser.add_argument(
        "-e",
        "--eps",
        dest="eps",
        type=float,
        default=0.0,
        help=(
            "Nonnegative zero tolerance used only in '--mode dist' (default: "
            "%(default)s)."
        ),
    )
    parser.add_argument(
        "-mn",
        "--mode_nz",
        dest="mode_nz",
        choices=("closed", "open", "off"),
        default="closed",
        help=(
            "Epsilon rule used only in '--mode dist': 'closed' drops "
            "'abs(v_i) <= eps', 'open' drops 'abs(v_i) < eps', and 'off' "
            "disables epsilon filtering. Positive-only filtering follows "
            "(default: '%(default)s')."
        ),
    )
    parser.add_argument(
        "-f",
        "--floor",
        dest="floor",
        type=float,
        default=0.0,
        help=(
            "Nonnegative lower bound applied after the raw '--mode dist' "
            "statistic as 'dep_min := max(dep_min, floor)' (default: "
            "%(default)s)."
        ),
    )

    parser.add_argument(
        "-sb",
        "--siz_bin",
        dest="siz_bin",
        type=int,
        default=10,
        help=(
            "Target signal-bin width in base pairs for 'frag' and 'norm'; it "
            "must be positive and smaller than '--siz_gen'. 'dist' ignores "
            "this value (default: %(default)s)."
        ),
    )
    parser.add_argument(
        "-sg",
        "--siz_gen",
        dest="siz_gen",
        type=int,
        default=12157105,
        help=(
            "Positive effective genome size in base pairs for 'frag' and "
            "'norm'. 'dist' ignores this value (default: %(default)s, "
            "appropriate for S. cerevisiae when retaining multi-mapping "
            "alignments)."
        ),
    )

    parser.add_argument(
        "-fp",
        "--flags_pe",
        dest="flags_pe",
        type=str,
        default=None,
        help=(
            "Comma-separated SAM FLAGs for paired-end main alignments in "
            "'frag' BAM/CRAM input. BED input ignores them. Accepts decimal "
            "or hexadecimal values (default: 99,1123,163,1187; e.g. 0x63)."
        ),
    )
    parser.add_argument(
        "--flags-pe",
        dest="flags_pe",
        type=str,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-fs",
        "--flags_se",
        dest="flags_se",
        type=str,
        default=None,
        help=(
            "Comma-separated SAM FLAGs for single-end main alignments in "
            "'frag' BAM/CRAM input. BED input ignores them. Accepts decimal "
            "or hexadecimal values (default: 0,16,1024,1040; e.g. 0x400)."
        ),
    )
    parser.add_argument(
        "--flags-se",
        dest="flags_se",
        type=str,
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "-dp",
        "--dp",
        dest="dp",
        type=int,
        default=24,
        help=(
            "Maximum number of decimal places retained for finite emitted "
            "values. Non-informative trailing zeros and a trailing decimal "
            "point are removed, and negative zero is emitted as '0' (default: "
            "%(default)s)."
        ),
    )

    argv_parse = sys.argv[1:] if argv is None else argv

    if not argv_parse:
        parser.print_help(sys.stderr)
        raise SystemExit(0)

    return parser.parse_args(argv_parse)


def _validate_norm_arguments(
    args: argparse.Namespace,
    flags_pe: str | None,
    flags_se: str | None,
) -> None:
    """
    Report arguments that normalization mode intentionally ignores.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed input-floor arguments.
    flags_pe : str | None
        Raw paired-end flag selection.
    flags_se : str | None
        Raw single-end flag selection.
    """

    ignored = (
        ("--fil_in", args.fil_in, args.fil_in is not None),
        ("--infmt", args.infmt, args.infmt is not None),
        ("--ref_fa", args.ref_fa, args.ref_fa is not None),
        ("--method", args.method, args.method != "qntl_nz"),
        ("--qntl_nz", args.qntl_nz, args.qntl_nz != 1.0),
        ("--coef", args.coef, args.coef is not None),
        ("--floor", args.floor, args.floor != 0.0),
        ("--eps", args.eps, args.eps != 0.0),
        ("--mode_nz", args.mode_nz, args.mode_nz != "closed"),
        ("--flags_pe", flags_pe, bool(flags_pe)),
        ("--flags_se", flags_se, bool(flags_se)),
    )

    for option, value, supplied in ignored:
        if supplied:
            _note_ignored(option, "in '--mode norm'", value)


def _validate_data_arguments(args: argparse.Namespace) -> str:
    """
    Validate input format and mode-specific numeric arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed input-floor arguments for distance or fragment mode.

    Returns
    -------
    format_name : str
        Canonical input format.

    Raises
    ------
    FileNotFoundError
        If a named input path does not exist.
    InputFloorValidationError
        If the required input or its format is invalid.
    ValueError
        If a numeric distance-mode argument is invalid.
    """

    if not args.fil_in:
        raise InputFloorValidationError(
            "Error: '--fil_in' is required when '--mode dist' or "
            "'--mode frag'.",
        )

    format_name = infer_input_format(args.fil_in, args.infmt)

    if args.mode == "dist":
        if format_name != "bedgraph":
            if args.fil_in == "-":
                message = (
                    "Error: When '--fil_in -' is used with '--mode dist', "
                    "provide '--infmt {bedGraph,bdg,bg}'."
                )
            else:
                message = (
                    f"Error: Unsupported file type for '--mode dist': "
                    f"'{args.fil_in}'. Provide bedGraph-like input "
                    "(bedGraph, bdg, or bg, optionally with .gz)."
                )

            raise InputFloorValidationError(message)

        validate_comparison(args.floor, "ge", 0.0, "floor", allow_none=False)
        validate_comparison(args.eps, "ge", 0.0, "eps", allow_none=False)
        validate_comparison(args.coef, "ge", 0.0, "coef", allow_none=True)

        if args.method == "qntl_nz" and (
            not math.isfinite(args.qntl_nz) or not 0.0 <= args.qntl_nz <= 100.0
        ):
            raise InputFloorValidationError(
                "Error: '--qntl_nz' must be finite and in [0, 100].",
            )

        if args.method != "qntl_nz" and args.qntl_nz != 1.0:
            print(
                "Note: '--qntl_nz' is ignored unless '--method qntl_nz'.",
                file=sys.stderr,
            )

        if args.method == "qntl_nz" and args.coef is not None:
            print(
                "Note: '--coef' is ignored for '--method qntl_nz'.",
                file=sys.stderr,
            )

    elif format_name not in {"bam", "cram", "bed"}:
        if args.fil_in == "-":
            message = (
                "Error: When '--fil_in -' is used with '--mode frag', provide "
                "'--infmt {bam,cram,bed}'."
            )
        else:
            input_path = args.fil_in
            mode = "--mode frag"
            detail = f"Unsupported file type for '{mode}': '{input_path}'."
            message = f"Error: {detail} Provide bam, cram, or bed/bed.gz."

        raise InputFloorValidationError(message)

    if args.fil_in != "-":
        check_exists(
            args.fil_in,
            kind="file",
            label=format_name.upper(),
        )

    _validate_cram_reference(format_name, args.ref_fa)

    return format_name


def _validate_input_floor_arguments(
    args: argparse.Namespace,
) -> tuple[str | None, str | None, str | None]:
    """
    Validate shared and mode-specific input-floor arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed input-floor arguments.

    Returns
    -------
    format_name, flags_pe, flags_se : tuple[
        str | None, str | None, str | None
    ]
        Canonical input format and raw paired- and single-end flag values.
    """

    _validate_mode_dimensions(
        args.mode,
        args.siz_bin,
        args.siz_gen,
    )

    validate_comparison(args.dp, "ge", 0, "dp", allow_none=False)

    flags_pe = getattr(args, "flags_pe", None)
    flags_se = getattr(args, "flags_se", None)

    if args.mode == "norm":
        _validate_norm_arguments(args, flags_pe, flags_se)
        format_name = None
    else:
        format_name = _validate_data_arguments(args)

    return format_name, flags_pe, flags_se


def _parse_fragment_flags(
    args: argparse.Namespace,
    flags_pe: str | None,
    flags_se: str | None,
) -> tuple[set[int] | None, set[int] | None]:
    """
    Parse fragment-mode alignment flags and report ignored selections.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed input-floor arguments.
    flags_pe : str | None
        Raw paired-end flag selection.
    flags_se : str | None
        Raw single-end flag selection.

    Returns
    -------
    paired_flags, single_flags : tuple[set[int] | None, set[int] | None]
        Parsed paired- and single-end flag allowlists.
    """

    paired_flags = None
    single_flags = None

    if args.mode == "frag":
        if flags_pe:
            paired_flags = parse_flag_csv(flags_pe, "flags_pe")

        if flags_se:
            single_flags = parse_flag_csv(flags_se, "flags_se")
    else:
        if flags_pe:
            _note_ignored("--flags_pe", "unless '--mode frag'", flags_pe)

        if flags_se:
            _note_ignored("--flags_se", "unless '--mode frag'", flags_se)

    return paired_flags, single_flags


def _print_input_floor_arguments(
    args: argparse.Namespace,
    paired_flags: set[int] | None,
    single_flags: set[int] | None,
    skp_pfx: tuple[str, ...],
) -> None:
    """
    Print the verbose input-floor argument banner.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed input-floor arguments.
    paired_flags : set[int] | None
        Parsed paired-end flag allowlist.
    single_flags : set[int] | None
        Parsed single-end flag allowlist.
    skp_pfx : tuple[str, ...]
        Parsed input-header prefixes.
    """

    if not args.verbose:
        return

    paired_flags_text = (
        f"None (default: {PAIRED_FLAGS})"
        if paired_flags is None
        else str(paired_flags)
    )
    single_flags_text = (
        f"None (default: {SINGLE_FLAGS})"
        if single_flags is None
        else str(single_flags)
    )

    with redirect_stdout(sys.stderr):
        print("####################################")
        print("## Arguments: compute_input_floor ##")
        print("####################################")
        print("")
        print("--verbose")
        print(f"--mode {args.mode}")
        print(f"--fil_in {args.fil_in}")

        if args.infmt is not None:
            print(f"--infmt {args.infmt}")

        if args.ref_fa is not None:
            print(f"--ref_fa {args.ref_fa}")

        print(f"--skp_pfx {skp_pfx}")

        if args.mode == "dist":
            print(f"--method {args.method}")

            if args.method == "qntl_nz":
                print(f"--qntl_nz {args.qntl_nz}")

            print(f"--coef {args.coef}")
            print(f"--eps {args.eps}")
            print(f"--mode_nz {args.mode_nz}")
            print(f"--floor {args.floor}")

        print(f"--siz_bin {args.siz_bin}")
        print(f"--siz_gen {args.siz_gen}")

        if args.mode == "frag":
            print(f"--flags_pe {paired_flags_text}")
            print(f"--flags_se {single_flags_text}")

        print(f"--dp {args.dp}")
        print("")
        print("")


def _compute_input_floor_from_args(
    args: argparse.Namespace,
    paired_flags: set[int] | None,
    single_flags: set[int] | None,
    skp_pfx: tuple[str, ...],
) -> float:
    """
    Compute one input floor from validated command-line arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Validated input-floor arguments.
    paired_flags : set[int] | None
        Parsed paired-end flag allowlist.
    single_flags : set[int] | None
        Parsed single-end flag allowlist.
    skp_pfx : tuple[str, ...]
        Parsed input-header prefixes.

    Returns
    -------
    floor : float
        Computed minimum input depth.
    """

    data_mode = args.mode in {"dist", "frag"}

    return compute_input_floor(
        fil_in=args.fil_in if data_mode else "-",
        siz_bin=args.siz_bin,
        siz_gen=args.siz_gen,
        mode=args.mode,
        method=args.method,
        qntl_nz=args.qntl_nz,
        coef=args.coef,
        floor=args.floor,
        eps=args.eps,
        mode_nz=args.mode_nz,
        paired_flags=paired_flags,
        single_flags=single_flags,
        skp_pfx=skp_pfx,
        infmt=args.infmt if data_mode else None,
        ref_fa=args.ref_fa if data_mode else None,
    )


def main(argv: list[str] | None = None) -> int:
    """
    Execute the primary control flow for the script.

    Parameters
    ----------
    argv : list[str] | None
        Arguments to parse. The process arguments are used by default.

    Returns
    -------
    status : int
        Zero on success and one for an anticipated computation failure.

    Raises
    ------
    SystemExit
        For parser-controlled help and argument-validation failures.

    Notes
    -----
    - Emits a note to stderr if '--flags_pe' / '--flags_se' are supplied for
      bed/bed.gz input (as flags are ignored for bed inputs).
    - Prints human-readable error messages to stderr on failure.
    - BrokenPipeError is handled in the '__main__' wrapper.
    """

    args = parse_args(argv)

    try:
        (
            format_name,
            flags_pe,
            flags_se,
        ) = _validate_input_floor_arguments(args)
    except (FileNotFoundError, InputFloorError, ValueError) as error:
        print(str(error), file=sys.stderr)

        return 1

    skp_pfx = parse_skp_pfx(args.skp_pfx, default=DEF_SKP_PFX)

    if args.mode in {"frag", "norm"} and args.siz_bin > (0.5 * args.siz_gen):
        print(
            "Warning: 'siz_bin' is a large fraction of 'siz_gen'; 'dep_min' "
            "may be very large.",
            file=sys.stderr,
        )

    try:
        paired_flags, single_flags = _parse_fragment_flags(
            args,
            flags_pe,
            flags_se,
        )
    except FlagParseError as error:
        print(str(error), file=sys.stderr)

        return 1

    has_alignment_flags = bool(paired_flags or single_flags)

    if args.mode == "frag" and format_name == "bed" and has_alignment_flags:
        print(
            "Note: '--flags_pe' / '--flags_se' are ignored for bed inputs.",
            file=sys.stderr,
        )

    _print_input_floor_arguments(args, paired_flags, single_flags, skp_pfx)

    try:
        dep_min = _compute_input_floor_from_args(
            args,
            paired_flags,
            single_flags,
            skp_pfx,
        )
    except InputFloorError as error:
        print(str(error), file=sys.stderr)

        return 1
    except Exception as error:
        print(f"Error: Unexpected error: {error}", file=sys.stderr)

        return 1

    if args.mode == "frag" and dep_min == 0:
        print(
            "Warning: No alignments counted in input; dep_min will be 0 in "
            "'frag' mode.",
            file=sys.stderr,
        )

    print(format_value(dep_min, args.dp))

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        with suppress(Exception):
            sys.stdout.close()

        with suppress(Exception):
            sys.stderr.close()

        raise SystemExit(0) from None
