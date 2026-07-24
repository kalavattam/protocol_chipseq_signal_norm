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

Usage
-----
python -m protocol_chipseq_signal_norm.cli.compute_input_floor [options] --mode <mode>

Parameters
----------
Input, mode, genome-size, bin-size, distribution-filtering, and formatting
options are parsed by 'parse_args()'.

Returns
-------
Prints one rounded floating-point floor recommendation to stdout.

See Also
--------
docs/dev/signal_stabilization.md
    Maintainer notes on distribution-, fragment-, and normalization-based floor
    computations.
"""

from __future__ import annotations

import argparse
import math

# import pysam  # Lazy-imported in 'count_aln_bam()'
import signal
import sys
from contextlib import redirect_stdout, suppress

from protocol_chipseq_signal_norm.utilities.utils_check import (
    check_cmp,
    check_exists,
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

#  Define global default SAM FLAG values
FLG_PE = {99, 1123, 163, 1187}  # Proper paired-end alignments
FLG_SE = {0, 16, 1024, 1040}    # Single-end alignments

# TODO: Add unit tests for track/blank-line parsing, 'siz_bin >= siz_gen',
# hexadecimal flag parsing, and deepTools-style normalization extensions.


def note_ignored(opt: str, reason: str, val) -> None:
    """
    Print a standard “ignored option” note when a mode does not use an arg.

    Parameters
    ----------
        opt : str
            Option name as shown to the user (e.g., '--fil_in').
        reason : str
            Short reason string (e.g., \"in '--mode norm'\").
        val :
            Value that triggered the note (printed only when useful).
    """
    if val is None:
        return

    print(f"Note: '{opt}' is ignored {reason}: got '{val}'.", file=sys.stderr)


def parse_flag_csv(csv_str: str, label: str) -> set[int]:
    """
    Parse a comma-separated string of SAM FLAGs into a set of ints. Accepts
    decimal (e.g., 99) or hex with 0x prefix (e.g., 0x63).

    Parameters
    ----------
        csv_str : str
            Comma-separated list of FLAG values to allow. Tokens may be decimal
            (e.g., '99') or hexadecimal with a '0x' prefix (e.g., '0x63').
            Whitespace around commas is ignored, and empty tokens are skipped.
        label : str
            Short name used in error messages to identify which option the
            flags came from (e.g., 'flags-pe' or 'flags-se').

    Returns
    -------
        flags : set[int]
            A non-empty set of integer FLAG values.

    Exits:
        sys.exit(1) on invalid tokens (ValueError) or if no flags parsed.
    """
    flags: set[int] = set()
    for tok in csv_str.split(','):
        s = tok.strip()
        if not s:
            continue
        try:
            base = 16 if s.lower().startswith('0x') else 10
            val = int(s, base)

            if not (0 <= val <= 0xFFFF):
                print(
                    f"Error: FLAG out of range (0..65535): {s}",
                    file=sys.stderr
                )
                sys.exit(1)
            flags.add(val)
        except ValueError:
            print(
                f"Error: Invalid FLAG '{s}' in '{label}'. Use decimal (e.g., "
                "99) or hex (e.g., 0x63).",
                file=sys.stderr
            )
            sys.exit(1)
    if not flags:
        print(
            f"Error: No valid FLAGs parsed for '{label}'.",
            file=sys.stderr
        )
        sys.exit(1)
    return flags


# TODO: rename b/c works with CRAM too
def count_aln_bam(
    pth_bam: str,
    flg_pe: set[int] | None = None,
    flg_se: set[int] | None = None
) -> int:
    """
    Count alignments in a bam using explicit FLAG allow-lists. By default,
    those are the following:
        - Paired-end alignments: {99, 1123, 163, 1187}
        - Single-end alignments: {0, 16, 1024, 1040}

    Note: By default, we intentionally count reads flagged as duplicates (e.g.,
    0x400/1024 and 0x410/1040 in the single-end set, 0x463/1123 and 0x4A3/1187
    in the paired-end set). This is necessary for rDNA-focused analyses where
    duplicate marking reflects genuine biological copy number and high-coverage
    repeats, not (just) PCR/optical artifacts.

    Assumption for paired-end alignments: counts represent read alignment-
    inferred fragments. In practice, using the {99,163,1123,1187} allow-list
    approximates counting primary, properly paired alignments as fragment
    proxies.

    Parameters
    ----------
        pth_bam : str
            Path to bam file (coordinate- or name-sorted). A bam index is not
            required for streaming counts, as we stream the file with
            'until_eof=True'.
        flg_pe : set[int] | None
            Allow-list of FLAGs to count as paired-end sequenced “main”
            alignments. If None, defaults to {99, 1123, 163, 1187}.
        flg_se : set[int] | None
            Allow-list of FLAG integers to count as single-end “main”
            alignments. If None, defaults to {0, 16, 1024, 1040}.

    Returns
        n : int
            Total count of alignments whose FLAG is in either allow-list.

    Exits:
        sys.exit(1) if the bam cannot be opened/read (FileNotFoundError,
        OSError, ValueError).
    """
    try:
        import pysam
    except ImportError:
        print("Error: 'pysam' is required for bam input.", file=sys.stderr)
        sys.exit(1)

    if flg_pe is None:
        flg_pe = FLG_PE
    if flg_se is None:
        flg_se = FLG_SE

    n_in = 0
    try:
        if pth_bam == "-":
            fh = pysam.AlignmentFile(sys.stdin.buffer, "rb")
        else:
            fh = pysam.AlignmentFile(pth_bam, "rb")
        with fh as fil_bam:
            for read in fil_bam.fetch(until_eof=True):
                if read.flag in flg_pe or read.flag in flg_se:
                    n_in += 1
    except (FileNotFoundError, OSError, ValueError) as e:
        print(
            f"Error: Cannot process bam file '{pth_bam}': {e}",
            file=sys.stderr
        )
        sys.exit(1)

    return n_in


def count_aln_bed(pth_bed: str, skp_pfx: tuple[str, ...]) -> int:
    """
    Count the number of records (lines) in a bed/bed.gz file, skipping blank
    lines and header/meta lines. Header prefixes are configurable via
    '--skp_pfx' (default: '#,track,browser').

    Parameters
    ----------
        pth_bed : str
            Path to bed file (can be gzipped).

    Returns
    -------
        n : int
            Total count of alignments (lines).

    Exits:
        System exit (status 1) on open/read errors (FileNotFoundError,
        PermissionError, Exception).
    """
    try:
        n = 0
        with open_in(pth_bed) as f:
            for line in f:
                if is_header(line, skp_pfx) or not line.strip():
                    continue
                n += 1
        return n
    except (FileNotFoundError, PermissionError) as e:
        print(
            f"Error: Cannot open bed file '{pth_bed}': {e}",
            file=sys.stderr
        )
        sys.exit(1)
    except Exception as e:
        print(
            f"Error: Unexpected error with '{pth_bed}': {e}",
            file=sys.stderr
        )
        sys.exit(1)


def infer_fmt(p: str, hint: str | None = None) -> str:
    """
    Infer file format label from the path suffix.

    Parameters
    ----------
        p : str
            Input file path. Path to the input file, or "-" for stdin.

        hint : str | None = None
            Optional format hint used only when 'p == "-"'. Must be one of
            {'bam', 'bed', 'bedGraph', 'bedgraph', 'bdg', 'bg'} if provided.

    Returns
    -------
        str
            - "bam" for bam.
            - "bed" for bed (optionally with .gz).
            - "bedgraph" for bedGraph, bedgraph, bdg, or bg (optionally with
              .gz).
            - "other" for otherwise.
            - When 'p == "-"' and a valid 'hint' is given, the hint value is
              returned.

    """
    if p == "-":
        #  For stdin, accept a format hint and normalize bedGraph aliases
        if hint in {"bam", "bed"}:
            return hint
        if hint in {"bedgraph", "bedGraph", "bdg", "bg"}:
            return "bedgraph"

    q = p.lower()
    if q.endswith(".bam"):
        return "bam"

    if q.endswith(".bed") or q.endswith(".bed.gz"):
        return "bed"

    if (
        q.endswith(".bedgraph") or q.endswith(".bedgraph.gz") or
        q.endswith(".bdg") or q.endswith(".bdg.gz") or
        q.endswith(".bg") or q.endswith(".bg.gz")
    ):
        return "bedgraph"

    return "other"


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
    flg_pe: set[int] | None = None,
    flg_se: set[int] | None = None,
    skp_pfx: tuple[str, ...] = DEF_SKP_PFX,
    infmt: str | None = None
) -> float:
    """
    Compute depth factor for input normalization.

    Parameters
    ----------
        fil_in : str
            Input file path. Path to input file:
                - mode=dist  bedGraph-like (bedGraph, bedgraph, bdg, or bg,
                             optionally with .gz)
                - mode=frag  bam or bed/bed.gz
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

                    where 'n' is the number of counted alignment records (bam
                    or bed/bed.gz). This matches the fragment-normalized signal
                    derivation used in the siQ-ChIP code paths.

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


        flg_pe : set[int] | None
            Optional allow-list of FLAG integers for paired-end “main”
            alignments when counting a bam. If None, defaults are used.

        flg_se : set[int] | None
            Optional allow-list of FLAG integers for single-end “main”
            alignments when counting a bam. If None, defaults are used.

        skp_pfx : tuple[str, ...]
            Prefixes to skip as header/meta lines bed/bed.gz and bedGraph-like
            inputs.

        infmt : str | None
            Optional format hint forwarded to 'infer_fmt' when 'fil_in == "-"'.
            Accepts "bam", "bed", or "bedgraph" (and bedGraph-like aliases such
            as "bedGraph", "bdg", "bg"). Ignored otherwise.

    Returns
    -------
        dep_min : float
            The minimum input depth to use as a denominator floor.

            In 'dist' mode, 'dep_min' is computed from a bedGraph-like track
            (column 4) using one of the following rules on filtered values.
            Values are restricted to finite values, filtered by 'eps' /
            'mode_nz', then restricted to positive values only ('v_i > 0'):

                - method=qntl_nz:

                    dep_min = Q_q({v_i : v_i > 0}),

                    where 'q = qntl_nz / 100' and 'qntl_nz' is in [0, 100].
                    Quantile selection uses a floor-based nearest-rank rule:

                        i = floor(q * (N - 1)) (clamped to [0, N - 1])
                        dep_min = sorted_vals[i]

                - method=frc_mdn_nz:

                    dep_min = coef × median({v_i : v_i > 0})

                - method=frc_avg_nz:

                    dep_min = coef × mean({v_i : v_i > 0})

                - method=min_nz:

                    dep_min = coef × min({v_i : v_i > 0})

            After the method-specific value is computed:

                dep_min := max(dep_min, floor)

            In 'frag' mode,

                ((n * b) / g) / (1 - (b / g)).

            In 'norm' mode,

                (b / g) / [1 - (b / g)].

            'norm' mode via algebraic simplification:

                ([(n * b) / g] / [1 - (b / g)]) / n = (b / g) / [1 - (b / g)].

    Exits:
        System exit (status 1) on the following:
            - invalid mode (not 'dist', 'frag', or 'norm');
            - 'siz_bin >= siz_gen';
            - unsupported fil_in extension (dist: bedGraph-like; frag: bam or
              bed/bed.gz);
            - I/O/open errors while counting (propagated as error messages);
            - invalid FLAG tokens if provided (handled upstream in
              'parse_flag_csv').
    """
    if mode not in {"dist", "frag", "norm"}:
        print(
            f"Error: Invalid mode '{mode}'. Use 'dist', 'frag', or 'norm'.",
            file=sys.stderr
        )
        sys.exit(1)

    #  Check that 'siz_bin' is not greater than or equal to 'siz_gen'
    if siz_bin >= siz_gen:
        print(
            "Error: 'siz_bin' must be smaller than 'siz_gen' (got "
            f"siz_bin={siz_bin}, siz_gen={siz_gen}).",
            file=sys.stderr
        )
        sys.exit(1)

    if mode == "dist":
        fmt = infer_fmt(fil_in, infmt)
        if fmt != "bedgraph":
            print(
                "Error: '--mode dist' requires a bedGraph-like input file "
                "(bedGraph, bedgraph, bdg, or bg, optionally with .gz).",
                file=sys.stderr
            )
            sys.exit(1)

        #  Read values with “positive-only” nonzero policy.
        #  This keeps the floor in the same conceptual space as denominators.
        vals = list(
            iter_vals_bdg(
                fil_in,
                eps=eps,
                mode_nz=mode_nz,
                skp_pfx=skp_pfx,
                nz_policy="pos",
            )
        )

        #  Enforce positive-only bins for denominator floors, regardless of
        #  epsilon mode, thereby keeping 'dep_min' in denominator units
        vals = [v for v in vals if v > 0.0]

        if not vals:
            print(
                "Error: No positive finite values found after filtering "
                "(cannot compute distribution-based 'dep_min').",
                file=sys.stderr
            )
            sys.exit(1)

        coef_eff = determine_coef_eff(method, coef)

        dep_min = pick_stabilizer(
            vals,
            method=method,
            coef=coef_eff,
            qntl_pct=qntl_nz,
            floor=floor,
            qntl_rule="floor",
        )

        #  If dep_min is nan (should only happen if vals empty), hard error
        if not math.isfinite(dep_min):
            print(
                "Error: 'dep_min' is non-finite; check input values and "
                "filters.",
                file=sys.stderr
            )
            sys.exit(1)

        return dep_min

    #  Compute ratio of bin size to effective genome size
    b_over_g = siz_bin / siz_gen

    #  Follow fast processing path for 'norm', as 'dep_min' is independent of
    #  'n' (see "General notes")
    if mode == "norm":
        dep_min = b_over_g / (1.0 - b_over_g)
        return dep_min

    #  Otherwise, follow processing path in which we tally 'n' ('frag')
    fmt = infer_fmt(fil_in, infmt)
    if fmt == "bam":
        n_in = count_aln_bam(fil_in, flg_pe=flg_pe, flg_se=flg_se)
    elif fmt == "bed":
        if (flg_pe or flg_se):
            print(
                "Note: '--flags-pe' / '--flags-se' are ignored for bed "
                "inputs.",
                file=sys.stderr
            )
        n_in = count_aln_bed(fil_in, skp_pfx)
    elif fmt == "bedgraph":
        print(
            "Error: '--mode frag' expects bam or bed/bed.gz (alignment "
            "records), not bedGraph.",
            file=sys.stderr
        )
        sys.exit(1)
    else:
        print(
            f"Error: Unsupported file type: {fil_in}. Provide a bam or "
            "bed/bed.gz file.",
            file=sys.stderr
        )
        sys.exit(1)

    #  Compute depth factor for 'frag'
    if n_in == 0:
        print(
            "Warning: No alignments counted in input; dep_min will be 0 in "
            "'frag' mode.",
            file=sys.stderr
        )
    dep_min = ((n_in * siz_bin) / siz_gen) / (1.0 - b_over_g)
    return dep_min


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments.
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
        )
    )
    add_help_cap(parser)
    parser.add_argument(
        "-v", "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Run script in verbose mode.\n\n"
    )
    parser.add_argument(
        "-fi", "--fil_in",
        dest="fil_in",
        type=str,
        required=False,
        default=None,
        help=(
            "Input file path. Path to input file, or use '-' for stdin.\n"
            "    - '--mode dist' requires a bedGraph-like file: bedGraph, "
            "bedgraph, bdg, or bg (optionally with .gz).\n"
            "    - '--mode frag' requires bam or bed/bed.gz.\n"
            "    - '--mode norm' fil_in is ignored.\n"
            "\n"
        )
    )
    parser.add_argument(
        "-if", "--infmt",
        dest="infmt",
        choices=["bam", "bed", "bedGraph", "bedgraph", "bdg", "bg"],
        default=None,
        help=(
            "Input format hint. Required only when '--fil_in -' (stdin). "
            "Choose 'bam', 'bed', 'bedGraph', 'bedgraph', 'bdg', or 'bg'.\n"
            "\n"
        )
    )
    parser.add_argument(
        "-sb", "--siz_bin",
        dest="siz_bin",
        type=int,
        default=10,
        help="Bin size in base pairs for bedGraph input (default: %(default)s).\n\n"
    )
    parser.add_argument(
        "-sg", "--siz_gen",
        dest="siz_gen",
        type=int,
        default=12157105,
        help=(
            "Effective genome size of the model organism (default: "
            "%(default)s [which is appropriate for S. cerevisiae when "
            "retaining multi-mapping alignments]).\n"
            "\n"
        )
    )
    parser.add_argument(
        "-fp", "--flags-pe",
        dest="flags_pe",
        type=str,
        help=(
            "Comma-separated string of SAM FLAGs to count as paired-end main "
            "alignments (default: 99,1123,163,1187). Accepts decimal or hex "
            "(e.g., 0x63).\n"
            "\n"
        )
    )
    parser.add_argument(
        "-fs", "--flags-se",
        dest="flags_se",
        type=str,
        help=(
            "Comma-separated string of SAM FLAGs to count as single-end main "
            "alignments (default: 0,16,1024,1040). Accepts decimal or hex "
            "(e.g., 0x400).\n"
            "\n"
        )
    )
    parser.add_argument(
        "-md", "--mode",
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
            "        + Values are then restricted to positive-only bins "
            "('v_i > 0'), so the resulting floor is in the same conceptual "
            "space as denominators in IP ÷ input.\n"
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
            "from an alignment-record file (bam or bed/bed.gz):\n"
            "\n"
            "        dep_min = ((n * b) / g) / [1 - (b / g)]\n"
            "\n"
            "    - norm: Compute a siQ-ChIP-style “normalized-coverage” "
            "factor depending only on bin size ('b') and genome size ('g'):\n"
            "\n"
            "        dep_min = (b / g) / [1 - (b / g)],\n"
            "\n"
            "    derived from the following algebraic simplification:\n"
            "\n"
            "        ([(n * b) / g] / [1 - (b / g)]) / n = (b / g) / [1 - "
            "(b / g)].\n"
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
        )
    )
    parser.add_argument(
        "-m", "--method",
        dest="method",
        choices=("qntl_nz", "frc_mdn_nz", "frc_avg_nz", "min_nz"),
        default="qntl_nz",
        help=(
            "Workflow method. Distribution-based method used only in '--mode dist' "
            "(default: %(default)s):\n"
            "    - qntl_nz     dep_min = q-th percentile of positive bins\n"
            "    - frc_mdn_nz  dep_min = coef × median of positive bins\n"
            "    - frc_avg_nz  dep_min = coef × mean of positive bins\n"
            "    - min_nz      dep_min = coef × minimum positive bin\n"
            "\n"
            "Notes:\n"
            "    - The input distribution is taken from bedGraph column 4 "
            "after skipping header/meta lines and non-data lines.\n"
            "    - Values are filtered to finite values, filtered by '--eps' "
            "/ '--mode_nz', then restricted to positive-only values "
            "('v_i > 0').\n"
            "    - '--qntl_nz' is used only with '--method qntl_nz'.\n"
            "    - '--coef' is used only with the 'frc_*' and 'min_nz' "
            "methods.\n"
            "    - After the method value is computed in 'mode=dist', "
            "'--floor' is applied as follows:\n"
            "\n"
            "        dep_min := max(dep_min, floor).\n"
            "\n"
        )
    )
    parser.add_argument(
        "-qn", "--qntl_nz",
        dest="qntl_nz",
        type=float,
        default=1.0,
        help=(
            "Quantile in percent used only when '--mode dist' AND "
            "'--method qntl_nz' (default: %(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-c", "--coef",
        dest="coef",
        type=float,
        default=None,
        help=(
            "Coefficient for '--method frc_mdn_nz', '--method frc_avg_nz', "
            "and '--method min_nz'. If omitted, defaults match "
            "'compute_pseudo.py' (0.01 for 'frc_*'; 1.0 for 'min_nz').\n\n"
        )
    )
    parser.add_argument(
        "-f", "--floor",
        dest="floor",
        type=float,
        default=0.0,
        help=(
            "Lower bound applied to the computed floor in '--mode dist' "
            "(default: %(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-e", "--eps",
        dest="eps",
        type=float,
        default=0.0,
        help=(
            "Zero tolerance epsilon used only in '--mode dist' (default: "
            "%(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-mn", "--mode_nz",
        dest="mode_nz",
        choices=("closed", "open", "off"),
        default="closed",
        help=(
            "Epsilon/zero-handling mode used only in '--mode dist' "
            "(default: %(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-dp", "--dp",
        dest="dp",
        type=int,
        default=24,
        help=(
            "Maximum number of decimal places retained for finite emitted values. Number of decimal places for rounding result (default: "
            "%(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-sp", "--skp_pfx",
        dest="skp_pfx",
        type=str,
        default=",".join(DEF_SKP_PFX),
        help=(
            "Comma-separated list of header prefixes to skip in "
            "bed/bed.gz or bedGraph-like inputs; to disable skipping, pass an "
            "empty string (default: %(default)s).\n\n"
        )
    )

    #  If no arguments are provided, use 'argv' to display help and exit
    argv_parse = sys.argv[1:] if argv is None else argv
    if not argv_parse:
        parser.print_help(sys.stderr)
        raise SystemExit(0)

    return parser.parse_args(argv_parse)


def main(argv: list[str] | None = None) -> int:
    """
    Execute the primary control flow for the script.

    Parameters
    ----------
        None.

    Returns
    -------
        None. On success, prints the “floor” (“input depth factor”) to stdout.

    Raises
    ------
        SystemExit
            On validation failures or when required dependencies are missing.

    Notes
    -----
        - Emits a note to stderr if '--flags-pe' / '--flags-se' are supplied
          for bed/bed.gz input (as flags are ignored for bed inputs).
        - Prints human-readable error messages to stderr on failure.
        - BrokenPipeError is handled in the '__main__' wrapper.
    """
    #  Parse CLI arguments
    args = parse_args(argv)

    #  Check arguments
    try:
        #  Always validate numeric args
        check_cmp(args.siz_bin, "gt", 0, "siz_bin", allow_none=False)
        check_cmp(args.siz_gen, "gt", 0, "siz_gen", allow_none=False)
        check_cmp(args.dp, "ge", 0, "dp", allow_none=False)

        #  Avoid accessing 'args.flags_(pe|se)' directly
        flags_pe = getattr(args, "flags_pe", None)
        flags_se = getattr(args, "flags_se", None)

        if args.mode == "norm":
            #  UX note: in 'norm' mode, fil_in/infmt/flags are ignored, as are
            #           'dist'-mode knobs
            if args.fil_in is not None:
                note_ignored("--fil_in", "in '--mode norm'", args.fil_in)
            if args.infmt is not None:
                note_ignored("--infmt", "in '--mode norm'", args.infmt)

            if flags_pe:
                note_ignored("--flags-pe", "in '--mode norm'", args.flags_pe)
            if flags_se:
                note_ignored("--flags-se", "in '--mode norm'", args.flags_se)

            if args.method != "qntl_nz":
                note_ignored("--method", "in '--mode norm'", args.method)
            if args.qntl_nz != 1.0:
                note_ignored("--qntl_nz", "in '--mode norm'", args.qntl_nz)
            if args.coef is not None:
                note_ignored("--coef", "in '--mode norm'", args.coef)
            if args.floor != 0.0:
                note_ignored("--floor", "in '--mode norm'", args.floor)
            if args.eps != 0.0:
                note_ignored("--eps", "in '--mode norm'", args.eps)
            if args.mode_nz != "closed":
                note_ignored("--mode_nz", "in '--mode norm'", args.mode_nz)

        else:
            #  In 'dist' and 'frag' modes, '--fil_in' is required
            if not args.fil_in:
                print(
                    "Error: '--fil_in' is required when '--mode dist' or "
                    "'--mode frag'.",
                    file=sys.stderr
                )
                raise SystemExit(1)

            fmt = infer_fmt(args.fil_in, args.infmt)

            if args.mode == "dist":
                #  Validate input type for 'mode=dist'
                if fmt != "bedgraph":
                    if args.fil_in == "-":
                        print(
                            "Error: When '--fil_in -' is used with '--mode "
                            "dist', provide '--infmt {bedgraph,bedGraph,bdg,"
                            "bg}'.",
                            file=sys.stderr
                        )
                    else:
                        print(
                            f"Error: Unsupported file type for '--mode dist': "
                            f"'{args.fil_in}'. Provide bedGraph-like input "
                            "(bedgraph, bedGraph, bdg, or bg, optionally with "
                            ".gz).",
                            file=sys.stderr
                        )
                    raise SystemExit(1)

                check_cmp(args.floor, "ge", 0.0, "floor", allow_none=False)
                check_cmp(args.eps,   "ge", 0.0, "eps",   allow_none=False)
                check_cmp(args.coef,  "ge", 0.0, "coef",  allow_none=True)

                if args.method == "qntl_nz" and (
                    not math.isfinite(args.qntl_nz)
                    or not 0.0 <= args.qntl_nz <= 100.0
                ):
                    raise SystemExit(
                        "Error: '--qntl_nz' must be finite and in [0, 100]."
                    )

                if args.method != "qntl_nz" and args.qntl_nz != 1.0:
                    print(
                        "Note: '--qntl_nz' is ignored unless '--method "
                        "qntl_nz'.",
                        file=sys.stderr
                    )

                if args.method == "qntl_nz" and args.coef is not None:
                    print(
                        "Note: '--coef' is ignored for '--method qntl_nz'.",
                        file=sys.stderr
                    )

            elif args.mode == "frag":
                if fmt not in {"bam", "bed"}:
                    if args.fil_in == "-":
                        print(
                            "Error: When '--fil_in -' is used with "
                            "'--mode frag', provide '--infmt {bam,bed}'.",
                            file=sys.stderr
                        )
                    else:
                        print(
                            f"Error: Unsupported file type for '--mode frag': "
                            f"'{args.fil_in}'. Provide bam or bed/bed.gz.",
                            file=sys.stderr
                        )
                    raise SystemExit(1)

            if args.fil_in != "-":
                check_exists(args.fil_in, kind="file", label=fmt.upper())

    except ValueError as e:
        raise SystemExit(str(e)) from None

    #  Parse header-skip prefixes from CLI (empty string: no skipping)
    skp_pfx = parse_skp_pfx(args.skp_pfx, default=DEF_SKP_PFX)

    if args.siz_bin > (0.5 * args.siz_gen):
        print(
            "Warning: 'siz_bin' is a large fraction of 'siz_gen'; 'dep_min' "
            "may be very large.",
            file=sys.stderr
        )

    #  Parse user-supplied flags if present
    flg_pe = None
    flg_se = None

    if args.mode == "frag":
        if flags_pe:
            flg_pe = parse_flag_csv(flags_pe, "flags-pe")
        if flags_se:
            flg_se = parse_flag_csv(flags_se, "flags-se")
    else:
        #  UX note: flags are ignored outside 'frag' mode
        if flags_pe:
            note_ignored("--flags-pe", "unless '--mode frag'", flags_pe)
        if flags_se:
            note_ignored("--flags-se", "unless '--mode frag'", flags_se)

    if flg_pe is None:
        flg_pe_str = f"None (default: {FLG_PE})"
    else:
        flg_pe_str = str(flg_pe)

    if flg_se is None:
        flg_se_str = f"None (default: {FLG_SE})"
    else:
        flg_se_str = str(flg_se)

    #  If '--verbose', print banner
    if args.verbose:
        with redirect_stdout(sys.stderr):
            print("####################################")
            print("## Arguments: compute_input_floor ##")
            print("####################################")
            print("")
            print("--verbose")
            print(f"--fil_in {args.fil_in}")
            if args.infmt is not None:
                print(f"--infmt {args.infmt}")
            print(f"--siz_bin {args.siz_bin}")
            print(f"--siz_gen {args.siz_gen}")
            if args.mode == "frag":
                print(f"--flags-pe {flg_pe_str}")
                print(f"--flags-se {flg_se_str}")
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

    #  Compute depth factor
    try:
        dep_min = compute_input_floor(
            fil_in=args.fil_in if args.mode in {"dist", "frag"} else "-",
            siz_bin=args.siz_bin,
            siz_gen=args.siz_gen,
            mode=args.mode,
            method=args.method,
            qntl_nz=args.qntl_nz,
            coef=args.coef,
            floor=args.floor,
            eps=args.eps,
            mode_nz=args.mode_nz,
            flg_pe=flg_pe,
            flg_se=flg_se,
            skp_pfx=skp_pfx,
            infmt=args.infmt if args.mode in {"dist", "frag"} else None
        )
    except Exception as e:
        #  Let explicit SystemExit from helpers “bubble” through; only catch
        #  the truly unexpected
        print(f"Error: Unexpected error: {e}", file=sys.stderr)
        raise SystemExit(1) from e

    #  Print result
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
