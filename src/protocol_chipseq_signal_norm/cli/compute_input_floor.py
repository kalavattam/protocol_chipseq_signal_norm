#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: compute_input_floor.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in development and documentation.
#
# Distributed under the MIT license.


"""
Compute denominator floors for ratio-based ChIP-seq signal workflows.

The CLI accepts input, mode, genome-size, bin-size, distribution-filtering, and
formatting options. It prints one rounded floating-point floor recommendation
to stdout.

Examples
--------
python -m protocol_chipseq_signal_norm.cli.compute_input_floor [options] \\
    --mode <mode>
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
    add_help_cap,
)
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


# TODO: Add unit tests for track/blank-line parsing, 'siz_bin >= siz_gen',
# hexadecimal flag parsing, and deepTools-style normalization extensions.


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
    high-coverage repeats can represent biological signal rather than only
    PCR or optical artifacts. Paired-end records serve as fragment proxies.

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


def infer_input_format(path: str, hint: str | None = None) -> str:
    """
    Infer file format label from the path suffix.

    Parameters
    ----------
    path : str
        Input file path. Path to the input file, or "-" for stdin.

    hint : str | None = None
        Optional format hint used only when 'path == "-"'. Must be one of
        {'bam', 'cram', 'bed', 'bedGraph', 'bedgraph', 'bdg', 'bg'} if
        provided.

    Returns
    -------
    format_name : str
        - "bam" for bam.
        - "cram" for cram.
        - "bed" for bed (optionally with .gz).
        - "bedgraph" for bedGraph, bedgraph, bdg, or bg (optionally with
          .gz).
        - "other" for otherwise.
        - When 'path == "-"' and a valid 'hint' is given, the hint value is
          returned.
    """

    if path == "-":
        # Standard input requires a format hint.
        if hint in {"bam", "cram", "bed"}:
            return hint

        if hint in {"bedgraph", "bedGraph", "bdg", "bg"}:
            return "bedgraph"

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
    Compute depth factor for input normalization.

    Parameters
    ----------
    fil_in : str
        Input file path. Path to input file:
            - mode=dist  bedGraph-like (bedGraph, bedgraph, bdg, or bg,
                         optionally with .gz)
            - mode=frag  bam, cram, or bed/bed.gz
            - mode=norm  fil_in is ignored

    siz_bin : int
        Output bin size in base pairs (i.e., the bedGraph bin width).

    siz_gen : int
        Effective genome size of the organism in base pairs.

    mode : str
        One of {'dist', 'frag', 'norm'}:

            - dist:
                Compute a distribution-based denominator floor from a
                bedGraph-like input (column 4). Values are
                    (1) filtered to finite values, then
                    (2) filtered by 'eps'/'mode_nz', then
                    (3) restricted to positive values only ('v_i > 0').

                The final floor is selected by 'method' and optional knobs:
                    + method=qntl_nz:

                        dep_min = Q_q({v_i : v_i > 0}),

                        where 'q' is '--qntl_nz' expressed as a fraction
                        (qntl_nz/100). Quantile selection uses a simple,
                        floor-based nearest-rank rule on sorted values:

                            i = floor(q * (N - 1)) (clamped to [0, N - 1])
                            dep_min = sorted_vals[i]

                    + method=frc_mdn_nz:

                        dep_min = coef × median({v_i : v_i > 0})

                    + method=frc_avg_nz:

                        dep_min = coef × mean({v_i : v_i > 0})

                    + method=min_nz:

                        dep_min = coef × min({v_i : v_i > 0})

                After the method-specific value is computed ('mode=dist'):

                    dep_min := max(dep_min, floor)

            - frag:

                dep_min = [(n * b) / g] / [1 - (b / g)],

                where 'n' is the number of counted alignment records (bam,
                cram, or bed/bed.gz). This matches the fragment-normalized
                signal derivation used in the siQ-ChIP code paths.

            - norm:

                dep_min = (b / g) / [1 - (b / g)].

                This matches “normalized coverage” derivations
                (siQ-ChIP-style), and depends only on 'b' and 'g'; fil_in
                is ignored.

    method : str
        Distribution-based rule used only when mode='dist'. One of:
        {'qntl_nz', 'frc_mdn_nz', 'frc_avg_nz', 'min_nz'}.

    qntl_nz : float
        Quantile in percent (0..100) used only when mode='dist' and
        method='qntl_nz'. The quantile is computed on positive values only
        (v_i > 0) after eps/mode_nz filtering.

    coef : float | None
        Coefficient used only when mode='dist' and method is one of
        {'frc_mdn_nz', 'frc_avg_nz', 'min_nz'}. If None, defaults match
        compute_pseudo.py (0.01 for frc_*; 1.0 for min_nz). Ignored for
        method='qntl_nz'.

    floor : float
        Nonnegative lower bound applied after computing the method’s raw
        value in mode='dist':
            dep_min := max(dep_min, floor)

    eps : float
        Zero tolerance epsilon used only in mode='dist' (see mode_nz).

    mode_nz : str
        Epsilon/zero-handling mode used only in mode='dist':
            - "closed"  drop values with |v_i| <= eps
            - "open"    drop values with |v_i| <  eps
            - "off"     disable epsilon-based filtering
        After epsilon filtering, distribution-based methods are restricted
        to positive values only (v_i > 0).


    paired_flags : set[int] | None
        Optional allow-list of FLAG integers for paired-end “main”
        alignments when counting a bam or cram. If None, defaults are used.

    single_flags : set[int] | None
        Optional allow-list of FLAG integers for single-end “main”
        alignments when counting a bam or cram. If None, defaults are used.

    skp_pfx : tuple[str, ...]
        Prefixes to skip as header/meta lines bed/bed.gz and bedGraph-like
        inputs.

    infmt : str | None
        Optional format hint forwarded to 'infer_input_format' when
        'fil_in == "-".
        Accepts "bam", "cram", "bed", or "bedgraph" (and bedGraph-like
        aliases such as "bedGraph", "bdg", "bg"). Ignored otherwise.
    ref_fa : str | None
        Reference FASTA required for CRAM decoding. Ignored for other input
        formats.

    Returns
    -------
    floor : float
        Minimum denominator floor. Distribution mode applies the selected
        statistic and lower bound; fragment and normalized modes apply their
        documented bin-to-genome-size formulas.

    Raises
    ------
    InputFloorValidationError
        If the mode, dimensions, input format, filtered values, or computed
        floor are invalid.
    AlignmentReadError
        If an alignment input cannot be read.
    """

    if mode not in {"dist", "frag", "norm"}:
        raise InputFloorValidationError(
            f"Error: Invalid mode '{mode}'. Use 'dist', 'frag', or 'norm'.",
        )

    if siz_bin >= siz_gen:
        raise InputFloorValidationError(
            "Error: 'siz_bin' must be smaller than 'siz_gen' (got "
            f"siz_bin={siz_bin}, siz_gen={siz_gen}).",
        )

    if mode == "dist":
        format_name = infer_input_format(fil_in, infmt)
        if format_name != "bedgraph":
            raise InputFloorValidationError(
                "Error: '--mode dist' requires a bedGraph-like input file "
                "(bedGraph, bedgraph, bdg, or bg, optionally with .gz).",
            )

        # Read values with “positive-only” nonzero policy.
        # This keeps the floor in the same conceptual space as denominators.
        vals = list(
            iter_vals_bdg(
                fil_in,
                eps=eps,
                mode_nz=mode_nz,
                skp_pfx=skp_pfx,
                nz_policy="pos",
            ),
        )

        # Positive-only bins keep `dep_min` in denominator units.
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
            "Compute the minimum input depth ('dep_min'), a “denominator "
            "floor” used to avoid extreme or erroneous divisions when "
            "normalizing by input (IP ÷ input).\n"
            "\n"
            "Downstream ratio code (e.g., 'compute_signal_ratio.py') can use "
            "the computed value to “clamp” denominators below 'dep_min' up to "
            "'dep_min'.\n"
            "\n"
            "Note the following variables in below equations and "
            "expressions:\n"
            "   - n: ...\n"
            "   - g: ...\n"
            "   - etc.\n"
            "\n"
        ),
    )
    add_help_cap(parser)
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Run script in verbose mode.\n\n",
    )

    parser.add_argument(
        "-fi",
        "--fil_in",
        dest="fil_in",
        type=str,
        required=False,
        default=None,
        help=(
            "Input file path. Path to input file, or use '-' for stdin.\n"
            "    - '--mode dist' requires a bedGraph-like file: bedGraph, "
            "bedgraph, bdg, or bg (optionally with .gz).\n"
            "    - '--mode frag' requires bam, cram, or bed/bed.gz.\n"
            "    - '--mode norm' fil_in is ignored.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-if",
        "--infmt",
        dest="infmt",
        choices=["bam", "cram", "bed", "bedGraph", "bedgraph", "bdg", "bg"],
        default=None,
        help=(
            "Input format hint. Required only when '--fil_in -' (stdin). "
            "Choose 'bam', 'cram', 'bed', 'bedGraph', 'bedgraph', 'bdg', or "
            "'bg'.\n\n"
        ),
    )
    parser.add_argument(
        "-rf",
        "--ref_fa",
        dest="ref_fa",
        default=None,
        help=(
            "Reference FASTA file required for CRAM decoding. Ignored for "
            "other input formats.\n\n"
        ),
    )

    parser.add_argument(
        "-sb",
        "--siz_bin",
        dest="siz_bin",
        type=int,
        default=10,
        help=(
            "Bin size in base pairs for bedGraph input (default: %(default)s)."
            "\n\n"
        ),
    )
    parser.add_argument(
        "-sg",
        "--siz_gen",
        dest="siz_gen",
        type=int,
        default=12157105,
        help=(
            "Effective genome size of the model organism (default: "
            "%(default)s [which is appropriate for S. cerevisiae when "
            "retaining multi-mapping alignments]).\n"
            "\n"
        ),
    )

    parser.add_argument(
        "-fp",
        "--flags-pe",
        dest="flags_pe",
        type=str,
        help=(
            "Comma-separated string of SAM FLAGs to count as paired-end main "
            "alignments (default: 99,1123,163,1187). Accepts decimal or hex "
            "(e.g., 0x63).\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-fs",
        "--flags-se",
        dest="flags_se",
        type=str,
        help=(
            "Comma-separated string of SAM FLAGs to count as single-end main "
            "alignments (default: 0,16,1024,1040). Accepts decimal or hex "
            "(e.g., 0x400).\n"
            "\n"
        ),
    )

    parser.add_argument(
        "-md",
        "--mode",
        dest="mode",
        choices=["dist", "frag", "norm"],
        default="dist",
        help=(
            "Workflow mode. Normalization mode (default: '%(default)s').\n"
            "    - dist: Compute a distribution-based denominator floor from "
            "a bedGraph-like input track (column 4). In this mode:\n"
            "        + Values are restricted to finite values ('NaN' / 'inf' "
            "ignored).\n"
            "        + Values are filtered by '--eps' / '--mode_nz' to define "
            "which bins are treated as “nonzero”.\n"
            "        + Values are then restricted to positive-only bins ('v_i "
            "> 0'), so the resulting floor is in the same conceptual space as "
            "denominators in IP ÷ input.\n"
            "        + '--method' selects how the filtered values are "
            "summarized into a single 'dep_min'.\n"
            "        + '--qntl_nz' is used only with '--method qntl_nz'.\n"
            "        + '--coef' is used only with '--method frc_mdn_nz', "
            "'--method frc_avg_nz', and '--method min_nz'.\n"
            "        + '--floor' is a nonnegative lower bound applied after "
            "the method-specific value is computed:\n"
            "\n"
            "              dep_min := max(dep_min, floor)\n"
            "\n"
            "    - frag: Compute a siQ-ChIP-style fragment-normalized factor "
            "from an alignment-record file (bam, cram, or bed/bed.gz):\n"
            "\n"
            "        dep_min = ((n * b) / g) / [1 - (b / g)]\n"
            "\n"
            "    - norm: Compute a siQ-ChIP-style “normalized-coverage” "
            "factor depending only on bin size ('b') and genome size ('g'):\n"
            "\n"
            "        dep_min = (b / g) / [1 - (b / g)],\n"
            "\n"
            "    The formula follows from this algebraic simplification.\n"
            "\n"
            "        ([(n * b) / g] / [1 - (b / g)]) / n = (b / g) / [1 - (b "
            "/ g)].\n"
            "\n"
            "    Here, normalized coverage means the corresponding genome-"
            "wide signal has been scaled so its integral (i.e., sum over all "
            "bins) is 1. Under this convention, the expected per-bin depth "
            "depends only on 'b' and 'g', so no input file is needed and "
            "'--fil_in' is ignored.\n"
            "\n"
            "Summary of calculations:\n"
            "    - dist (bedGraph column 4 values 'v_i'; computed on positive-"
            "only bins):\n"
            "        + method=qntl_nz\n"
            "\n"
            "            dep_min = q-th percentile of {v_i : v_i > 0}\n"
            "\n"
            "        + method=frc_mdn_nz\n"
            "\n"
            "            dep_min = coef × median({v_i : v_i > 0})\n"
            "\n"
            "        + method=frc_avg_nz\n"
            "\n"
            "            dep_min = coef × mean({v_i : v_i > 0})\n"
            "\n"
            "        + method=min_nz\n"
            "\n"
            "            dep_min = coef × min({v_i : v_i > 0})\n"
            "\n"
            "        + then\n"
            "\n"
            "            dep_min := max(dep_min, floor)\n"
            "\n"
            "    - frag:\n"
            "\n"
            "        dep_min = ((n * b) / g) / [1 - (b / g)]\n"
            "\n"
            "    - norm:\n"
            "\n"
            "        dep_min = (b / g) / [1 - (b / g)]\n"
            "\n"
            "Notes:\n"
            "    - 'frag' and 'norm' are included here for consistency with "
            "the siQ-ChIP derivations/code paths by Brad Dickson.\n"
            "    - Intended downstream behavior:\n"
            "\n"
            "        denominator := max(denominator, dep_min)\n"
            "\n"
        ),
    )

    parser.add_argument(
        "-m",
        "--method",
        dest="method",
        choices=("qntl_nz", "frc_mdn_nz", "frc_avg_nz", "min_nz"),
        default="qntl_nz",
        help=(
            "Workflow method. Distribution-based method used only in '--mode "
            "dist' (default: %(default)s):\n"
            "    - qntl_nz     dep_min = q-th percentile of positive bins\n"
            "    - frc_mdn_nz  dep_min = coef × median of positive bins\n"
            "    - frc_avg_nz  dep_min = coef × mean of positive bins\n"
            "    - min_nz      dep_min = coef × minimum positive bin\n"
            "\n"
            "Notes:\n"
            "    - The input distribution is taken from bedGraph column 4 "
            "after skipping header/meta lines and non-data lines.\n"
            "    - Values are filtered to finite values, filtered by '--eps' "
            "/ '--mode_nz', then restricted to positive-only values ('v_i > "
            "0').\n"
            "    - '--qntl_nz' is used only with '--method qntl_nz'.\n"
            "    - '--coef' is used only with the 'frc_*' and 'min_nz' "
            "methods.\n"
            "    - After the method value is computed in 'mode=dist', "
            "'--floor' is applied as follows:\n"
            "\n"
            "        dep_min := max(dep_min, floor).\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-qn",
        "--qntl_nz",
        dest="qntl_nz",
        type=float,
        default=1.0,
        help=(
            "Quantile in percent used only when '--mode dist' AND '--method "
            "qntl_nz' (default: %(default)s).\n\n"
        ),
    )
    parser.add_argument(
        "-c",
        "--coef",
        dest="coef",
        type=float,
        default=None,
        help=(
            "Coefficient for '--method frc_mdn_nz', '--method frc_avg_nz', "
            "and '--method min_nz'. If omitted, defaults match "
            "'compute_pseudo.py' (0.01 for 'frc_*'; 1.0 for 'min_nz').\n\n"
        ),
    )
    parser.add_argument(
        "-f",
        "--floor",
        dest="floor",
        type=float,
        default=0.0,
        help=(
            "Lower bound applied to the computed floor in '--mode dist' "
            "(default: %(default)s).\n\n"
        ),
    )
    parser.add_argument(
        "-e",
        "--eps",
        dest="eps",
        type=float,
        default=0.0,
        help=(
            "Zero tolerance epsilon used only in '--mode dist' (default: "
            "%(default)s).\n\n"
        ),
    )
    parser.add_argument(
        "-mn",
        "--mode_nz",
        dest="mode_nz",
        choices=("closed", "open", "off"),
        default="closed",
        help=(
            "Epsilon/zero-handling mode used only in '--mode dist' (default: "
            "%(default)s).\n\n"
        ),
    )

    parser.add_argument(
        "-dp",
        "--dp",
        dest="dp",
        type=int,
        default=24,
        help=(
            "Maximum number of decimal places retained for finite emitted "
            "values. Number of decimal places for rounding result (default: "
            "%(default)s).\n\n"
        ),
    )
    parser.add_argument(
        "-sp",
        "--skp_pfx",
        dest="skp_pfx",
        type=str,
        default=",".join(DEF_SKP_PFX),
        help=(
            "Comma-separated list of header prefixes to skip in bed/bed.gz or "
            "bedGraph-like inputs; to disable skipping, pass an empty string "
            "(default: %(default)s).\n\n"
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
        ("--flags-pe", flags_pe, bool(flags_pe)),
        ("--flags-se", flags_se, bool(flags_se)),
        ("--method", args.method, args.method != "qntl_nz"),
        ("--qntl_nz", args.qntl_nz, args.qntl_nz != 1.0),
        ("--coef", args.coef, args.coef is not None),
        ("--floor", args.floor, args.floor != 0.0),
        ("--eps", args.eps, args.eps != 0.0),
        ("--mode_nz", args.mode_nz, args.mode_nz != "closed"),
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
                    "provide '--infmt {bedgraph,bedGraph,bdg,bg}'."
                )
            else:
                message = (
                    f"Error: Unsupported file type for '--mode dist': "
                    f"'{args.fil_in}'. Provide bedGraph-like input "
                    "(bedgraph, bedGraph, bdg, or bg, optionally with .gz)."
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

    validate_comparison(args.siz_bin, "gt", 0, "siz_bin", allow_none=False)
    validate_comparison(args.siz_gen, "gt", 0, "siz_gen", allow_none=False)
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
            paired_flags = parse_flag_csv(flags_pe, "flags-pe")

        if flags_se:
            single_flags = parse_flag_csv(flags_se, "flags-se")
    else:
        if flags_pe:
            _note_ignored("--flags-pe", "unless '--mode frag'", flags_pe)

        if flags_se:
            _note_ignored("--flags-se", "unless '--mode frag'", flags_se)

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
        print(f"--fil_in {args.fil_in}")

        if args.infmt is not None:
            print(f"--infmt {args.infmt}")

        if args.ref_fa is not None:
            print(f"--ref_fa {args.ref_fa}")

        print(f"--siz_bin {args.siz_bin}")
        print(f"--siz_gen {args.siz_gen}")

        if args.mode == "frag":
            print(f"--flags-pe {paired_flags_text}")
            print(f"--flags-se {single_flags_text}")

        print(f"--mode {args.mode}")

        if args.mode == "dist":
            print(f"--method {args.method}")

            if args.method == "qntl_nz":
                print(f"--qntl_nz {args.qntl_nz}")

            print(f"--coef {args.coef}")
            print(f"--floor {args.floor}")
            print(f"--eps {args.eps}")
            print(f"--mode_nz {args.mode_nz}")

        print(f"--dp {args.dp}")
        print(f"--skp_pfx {skp_pfx}")
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
    - Emits a note to stderr if '--flags-pe' / '--flags-se' are supplied
      for bed/bed.gz input (as flags are ignored for bed inputs).
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

    if args.siz_bin > (0.5 * args.siz_gen):
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
            "Note: '--flags-pe' / '--flags-se' are ignored for bed inputs.",
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

    print(f"{dep_min:.{args.dp}f}")

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
