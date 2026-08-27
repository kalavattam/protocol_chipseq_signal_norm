#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: compute_pseudo.py
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
Compute pseudocount recommendations for bedGraph signal tracks.

The CLI accepts input tracks and a method. '--method edger' derives the
pseudocount from library sizes, taking a target normalization, a prior count,
and optionally the sizes themselves; the four distribution-based methods derive
it from the value distribution, taking a coefficient, nonzero filtering, and a
symmetry mode. It prints one pseudocount or an A:B pair to stdout, or a
deepTools argument string under '--prt_arg', and with '--prt_jsn' also prints a
JSON summary.

Examples
--------
python -m protocol_chipseq_signal_norm.cli.compute_pseudo \\
    --fil_A <file> [--fil_B <file>] [options]
"""

from __future__ import annotations

import argparse
import json
import math
import signal
import sys
from contextlib import redirect_stdout, suppress

from protocol_chipseq_signal_norm.utilities.utils_bdg import sum_counts_bdg
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
    ensure_single_stdin,
    parse_skp_pfx,
)
from protocol_chipseq_signal_norm.utilities.utils_stabilizer import (
    NORM_CHOICES,
    canonicalize_norm,
    compute_pseudo_edger,
    compute_stats_robust,
    determine_coef_eff,
    iter_vals_bdg,
    pick_stabilizer,
)

with suppress(AttributeError, ValueError):
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

assert sys.version_info >= (3, 11), "Python >= 3.11 required."

# The distribution-shaping options, each with every spelling argparse accepts
# for it. 'edger' consults none of them, and silently ignoring one lets a
# wrong mental model survive.
OPT_IGNORED_EDGER = {
    "coef": ("-c", "--coef"),
    "qntl_nz": ("-q", "--qntl_nz"),
    "floor": ("-fl", "--floor"),
    "eps": ("-e", "--eps"),
    "mode_nz": ("-mz", "--mode_nz"),
    "sym": ("-s", "--sym"),
}

# Every short option 'parse_args' registers. Resolving an attached value needs
# the whole set rather than the six above: '-sfA0.5' must resolve to '-sfA',
# not to '-s' carrying 'fA0.5', and only the longest registered prefix tells
# the two apart. 'test_pseudo.py' asserts both constants against the parser,
# so neither can drift away from it.
OPT_SHORT_ALL = (
    "-h",
    "-v",
    "-fA",
    "-fB",
    "-sp",
    "-m",
    "-q",
    "-c",
    "-fl",
    "-e",
    "-mz",
    "-s",
    "-nm",
    "-pc",
    "-sb",
    "-lA",
    "-lB",
    "-sfA",
    "-sfB",
    "-gA",
    "-gB",
    "-dp",
    "-pj",
    "-pa",
)


# TODO: Extend compressed-input handling to '.bgz' and '.bgzf' here and in
# related bedGraph parsers. Revisit JSON summarization and whether a shared
# warning/error-text helper should be exported to other Python CLIs.
def combine_pseudo_sym(
    pseudo_a: float,
    pseudo_b: float,
    mode: str = "none",
) -> tuple[float, float]:
    """
    Combine per-track pseudocounts into a symmetric pair if requested.

    Parameters
    ----------
    pseudo_a : float
        Pseudocount proposed for track A (may be finite or nonfinite).
    pseudo_b : float
        Pseudocount proposed for track B (may be finite or nonfinite).
    mode : str
        Symmetrization rule; see notes below.

    Returns
    -------
    pseudo_a, pseudo_b : tuple[float, float]
        A pair (pseudo_A, pseudo_B) after applying the symmetrization policy.

    Raises
    ------
    ValueError
        If both pseudocounts are finite and 'mode' is not a recognized
        symmetrization rule.

    Notes
    -----
    - 'none' returns both inputs as-is, including nonfinite values. With two
      finite inputs, 'max' and 'min' select the respective input; 'arith',
      'geom', and 'harm' use the arithmetic, geometric, and harmonic means; and
      'use_A' and 'use_B' copy the selected input to both outputs.
    - For finite inputs, 'geom' falls back to the minimum with a stderr warning
      if either input is negative; 'harm' does the same if either input is
      nonpositive.
    - For every mode except 'none':
        + If exactly one is finite, then mirror the finite value while issuing
          a warning.
        + If both are nonfinite, then return as-is (i.e., let the user decide
          how to proceed) while issuing a warning. These nonfinite paths do not
          validate 'mode'.

    Examples
    --------
    >>> combine_pseudo_sym(1.0, 3.0, mode="arith")
    (2.0, 2.0)

    >>> combine_pseudo_sym(1.0, 3.0, mode="use_A")
    (1.0, 1.0)
    """

    if mode == "none":
        return pseudo_a, pseudo_b

    finite_a = math.isfinite(pseudo_a)
    finite_b = math.isfinite(pseudo_b)

    if finite_a and finite_b:
        if mode == "max":
            pseudo = max(pseudo_a, pseudo_b)
        elif mode == "min":
            pseudo = min(pseudo_a, pseudo_b)
        elif mode == "arith":
            pseudo = 0.5 * (pseudo_a + pseudo_b)
        elif mode == "geom":
            if pseudo_a < 0.0 or pseudo_b < 0.0:
                print(
                    "Geometric mean undefined for negative values; falling "
                    "back to min(pseudo_A, pseudo_B).",
                    file=sys.stderr,
                )
                pseudo = min(pseudo_a, pseudo_b)
            else:
                pseudo = math.sqrt(pseudo_a * pseudo_b)
        elif mode == "harm":
            if pseudo_a <= 0.0 or pseudo_b <= 0.0:
                print(
                    "Harmonic mean undefined for nonpositive values; falling "
                    "back to min(pseudo_A, pseudo_B).",
                    file=sys.stderr,
                )
                pseudo = min(pseudo_a, pseudo_b)
            else:
                pseudo = 2.0 / (1.0 / pseudo_a + 1.0 / pseudo_b)
        elif mode == "use_A":
            pseudo = pseudo_a
        elif mode == "use_B":
            pseudo = pseudo_b
        else:
            raise ValueError(f"Error: Unknown --sym: {mode!r}")

        return pseudo, pseudo

    if finite_a and not finite_b:
        print(
            "pseudo_B is nonfinite; mirroring pseudo_A in symmetric mode "
            f"{mode!r}.",
            file=sys.stderr,
        )

        return pseudo_a, pseudo_a

    if finite_b and not finite_a:
        print(
            "pseudo_A is nonfinite; mirroring pseudo_B in symmetric mode "
            f"{mode!r}.",
            file=sys.stderr,
        )

        return pseudo_b, pseudo_b

    print("Both pseudocounts are nonfinite; returning as-is.", file=sys.stderr)

    return pseudo_a, pseudo_b


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
        Parsed bedGraph inputs and pseudocount policy.

    Raises
    ------
    SystemExit
        If argument parsing fails or help is requested.
    """

    parser = CapArgumentParser(
        description=(
            "Compute data-driven pseudocount(s) for one (or two) bedGraph "
            "track(s)."
        ),
        prog="compute_pseudo",
        _sectioned_help=_SectionedHelpConfig(
            usage_rows=(
                ("help", "verbose"),
                ("fil_A", "fil_B", "skp_pfx"),
                (
                    "method",
                    "qntl_nz",
                    "coef",
                    "floor",
                    "eps",
                    "mode_nz",
                    "sym",
                ),
                (
                    "normalization",
                    "prior_count",
                    "siz_bin",
                    "lib_A",
                    "lib_B",
                    "sf_A",
                    "sf_B",
                    "frg_A",
                    "frg_B",
                ),
                ("dp", "prt_jsn", "prt_arg"),
            ),
            examples=(
                _HelpExample(
                    description=(
                        "Compute a default pseudocount from one bedGraph "
                        "track."
                    ),
                    command_lines=("compute_pseudo --fil_A signal_A.bdg",),
                ),
                _HelpExample(
                    description=(
                        "Compute first-percentile pseudocounts and symmetrize "
                        "them by maximum."
                    ),
                    command_lines=(
                        "compute_pseudo",
                        "--fil_A signal_A.bdg",
                        "--fil_B signal_B.bdg",
                        "--method qntl_nz",
                        "--qntl_nz 1",
                        "--sym max",
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
        help="Run script in verbose mode.\n\n",
    )

    parser.add_argument(
        "-fA",
        "--fil_A",
        dest="fil_A",
        required=True,
        help=(
            "First bedGraph input file, file A. Use '-' for stdin; '.gz' is "
            "handled.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-fB",
        "--fil_B",
        dest="fil_B",
        help=(
            "Second bedGraph input file, file B. This file is optional for "
            "single-file pseudocount modes; use '-' for stdin; '.gz' is "
            "handled.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-sp",
        "--skp_pfx",
        dest="skp_pfx",
        type=str,
        default=",".join(DEF_SKP_PFX),
        help=(
            "Comma-separated list of header prefixes to skip in bedGraph "
            "file(s) (default: %(default)s).\n"
            "\n"
        ),
    )

    parser.add_argument(
        "-m",
        "--method",
        dest="method",
        choices=("edger", "frc_mdn_nz", "qntl_nz", "frc_avg_nz", "min_nz"),
        default="edger",
        help=(
            "Workflow method to compute per-track pseudocount (default: "
            "%(default)s):\n"
            "    - edger       edgeR's prior.count rule; see "
            "'--normalization'\n"
            "    - frc_mdn_nz  value = coef × median of nonzero bins\n"
            "    - qntl_nz     value = q-th percentile of nonzero bins\n"
            "    - frc_avg_nz  value = coef × mean of nonzero bins\n"
            "    - min_nz      value = coef × minimum nonzero bin\n"
            "\n"
            "Notes:\n"
            "    - '--method edger' takes one track or two, as edgeR does.\n"
            "    - With one track, edgeR's per-sample prior scaling "
            "degenerates to a no-op and the mean library size is that track's "
            "own.\n"
            "    - A given track's one-track pseudocount is therefore not its "
            "two-track pseudocount; both are correct for their own frame.\n"
            "    - If '--method qntl_nz', set percentile with '--qntl_nz' "
            "[decimals OK (e.g., 0.1 = 0.1th percentile); nearest-rank "
            "determined via 'round'].\n"
            "    - 'nonzero' means '|x| > eps' with '--mode_nz closed', '|x| "
            ">= eps' with '--mode_nz open', and every finite value with "
            "'--mode_nz off'.\n"
            "    - If '--coef' is omitted, then defaults to '--coef 0.01' for "
            "'--method frc_*' and '--coef 1.0' for '--method min_nz'.\n"
            "    - '--method min_nz' typically needs a larger coef (e.g., "
            "0.1–1.0) in comparison to '--method frc_*' (e.g., 0.01).\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-q",
        "--qntl_nz",
        dest="qntl_nz",
        type=float,
        default=1.0,
        help=(
            "Quantile in percent for '--method qntl_nz' (0..100). Decimals "
            "are allowed (e.g., 0.5 = 0.5th percentile). Ignored if "
            "'--method' is not 'qntl_nz' (default: %(default)s).\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-c",
        "--coef",
        dest="coef",
        type=float,
        default=None,
        help=(
            "Coefficient for median, mean, and min methods. If not specified, "
            "then defaults to 0.01 for '--method frc_mdn_nz' and '--method "
            "frc_avg_nz', or 1.0 for '--method min_nz'.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-fl",
        "--floor",
        dest="floor",
        type=float,
        default=0.0,
        help=(
            "Lower bound for computation of pseudocount(s) (default: "
            "%(default)s).\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-e",
        "--eps",
        dest="eps",
        type=float,
        default=0.0,
        help=(
            "Zero tolerance epsilon in computation of pseudocount(s) "
            "(default: %(default)s).\n"
            "\n"
            "When '--mode_nz closed' (default), values with '|x| <= ε' are "
            "treated as zero and excluded from statistics; with '--mode_nz "
            "open', values with '|x| < ε' are excluded; with '--mode_nz off', "
            "ε-based filtering is disabled.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-mz",
        "--mode_nz",
        dest="mode_nz",
        choices=("closed", "open", "off"),
        default="closed",
        help=(
            "Epsilon/zero-handling mode for selecting 'nonzero' bins:\n"
            "    - closed (default)  drop values with '|x| <= eps'.\n"
            "    - open              drop values with '|x| < eps'.\n"
            "    - off               disable 'eps'-based zero filtering.\n"
            "\n"
        ),
    )

    parser.add_argument(
        "-s",
        "--sym",
        dest="sym",
        choices=(
            "none",
            "max",
            "min",
            "arith",
            "geom",
            "harm",
            "use_A",
            "use_B",
        ),
        default="none",
        help=(
            "Symmetrize pseudocounts across A and B.\n"
            "    - none returns A and B as computed; max/min select the "
            "larger/smaller value; arith/geom/harm use "
            "arithmetic/geometric/harmonic means; use_A/use_B copy A/B to "
            "both outputs.\n"
            "    - With finite values, geom falls back to min for a negative "
            "input and harm falls back to min for a nonpositive input; each "
            "fallback writes a warning to stderr.\n"
            "    - In any non-none mode, one finite value is mirrored with a "
            "warning; two nonfinite values are retained with a warning.\n"
            "    - If both A AND B are given, apply the chosen rule to "
            "'(pseudo_A, pseudo_B)', then print 'pseudo_A:pseudo_B'.\n"
            "    - If only A is given AND '--sym' is provided AND NOT 'none', "
            "mirror A to B and print: 'pseudo_A:pseudo_B'.\n"
            "    - If only A is given AND '--sym' is omitted OR '--sym none', "
            "print a single value: 'pseudo_A'.\n"
            "\n"
        ),
    )

    parser.add_argument(
        "-nm",
        "--normalization",
        dest="normalization",
        choices=NORM_CHOICES,
        default="CPM",
        help=(
            "Target deepTools normalization (default: %(default)s).\n"
            "    - CPM, BPM  exact; BPM is CPM in deepTools\n"
            "    - RPKM      exact; needs '--siz_bin'\n"
            "    - None      correct up to a constant log offset; not edgeR\n"
            "    - RPGC      needs '--sf_A' and '--sf_B'; not edgeR\n"
            "    - norm      normalized coverage; not edgeR. Needs '--frg_A' "
            "and '--frg_B'; pass the raw-count tracks as '--fil_A' and "
            "'--fil_B', not the normalized coverage tracks\n"
            "\n"
            "Applies to '--method edger' only; ignored otherwise.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-pc",
        "--prior_count",
        dest="prior_count",
        type=float,
        default=2.0,
        help=(
            "edgeR 'prior.count' before library-size scaling (default: "
            "%(default)s).\n"
            "\n"
            "Applies to '--method edger' only; ignored otherwise.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-sb",
        "--siz_bin",
        dest="siz_bin",
        type=int,
        default=None,
        help=(
            "Bin width in bp used to write the track(s).\n"
            "\n"
            "Applies to '--method edger' only; ignored otherwise.\n"
            "\n"
            "Inferred from the track when omitted, and cross-checked against "
            "it when given, so a value the track contradicts is refused "
            "rather than silently rescaling the library size.\n"
            "\n"
            "Required only when no track is read, i.e., when both '--lib_A' "
            "and '--lib_B' are supplied and '--normalization RPKM' needs a "
            "width for its scale factor.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-lA",
        "--lib_A",
        dest="lib_A",
        type=float,
        default=None,
        help=(
            "Library size for track A: the column sum of the bin matrix, "
            "which is edgeR's 'lib.size'. Computed from '--fil_A' when "
            "omitted, which requires a non-normalized track. Supplying it "
            "skips that read and changes nothing else.\n"
            "\n"
            "Applies to '--method edger' only; ignored otherwise.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-lB",
        "--lib_B",
        dest="lib_B",
        type=float,
        default=None,
        help=(
            "Library size for track B; see '--lib_A'. Omit this and '--fil_B' "
            "to compute a single-track pseudocount.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-sfA",
        "--sf_A",
        dest="sf_A",
        type=float,
        default=None,
        help=(
            "Scaling factor read from 'bamCoverage --verbose'. Generate that "
            "run with '--exactScaling'; otherwise, the factor is a sampled "
            "estimate and so is every pseudocount derived from it.\n"
            "\n"
            "Applies to '--method edger --normalization RPGC' only, where it "
            "is required.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-sfB",
        "--sf_B",
        dest="sf_B",
        type=float,
        default=None,
        help=("deepTools scale factor for track B; see '--sf_A'.\n\n"),
    )
    parser.add_argument(
        "-gA",
        "--frg_A",
        dest="frg_A",
        type=float,
        default=None,
        help=(
            "Fragment count for track A: the number of fragments "
            "'compute_signal' divided by, not the number of alignment "
            "records. The two coincide only when exactly one alignment per "
            "fragment survives filtering. A normalized-coverage track sums to "
            "1 by construction, so its own total cannot supply this. Get it "
            "from, e.g., the fragment count 'compute_signal' reports, or as "
            "1e6 divided by the CPM scale factor that 'bamCoverage --verbose' "
            "reports.\n"
            "\n"
            "Applies to '--method edger' with normalized coverage "
            "('--normalization norm', aliases 'nc', 'n', 'nrm', 'normalized') "
            "only, where it is required.\n"
        ),
    )
    parser.add_argument(
        "-gB",
        "--frg_B",
        dest="frg_B",
        type=float,
        default=None,
        help=("Fragment count for track B; see '--frg_A'.\n\n"),
    )
    parser.add_argument(
        "-dp",
        "--dp",
        dest="dp",
        type=int,
        default=24,
        help=(
            "Maximum number of decimal places retained for finite emitted "
            "values (default: %(default)s).\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-pj",
        "--prt_jsn",
        dest="prt_jsn",
        action="store_true",
        default=False,
        help="Print a JSON summary to stdout.\n\n",
    )
    parser.add_argument(
        "-pa",
        "--prt_arg",
        dest="prt_arg",
        action="store_true",
        default=False,
        help=(
            "Print a ready-to-paste deepTools argument string instead of the "
            "bare pseudocount pair.\n"
            "\n"
            "Applies to '--method edger' only; ignored otherwise.\n"
            "\n"
        ),
    )

    argv_parse = sys.argv[1:] if argv is None else argv

    if not argv_parse:
        parser.print_help(sys.stderr)
        raise SystemExit(0)

    return parser.parse_args(argv_parse)


def _print_pseudo_arguments(
    args: argparse.Namespace,
    coef_eff: float | None,
    skp_pfx: tuple[str, ...],
    siz_bin: int | None = None,
) -> None:
    """
    Print the verbose argument banner in the reviewed semantic order.
    """

    with redirect_stdout(sys.stderr):
        print("####################################")
        print("## Arguments for 'compute_pseudo' ##")
        print("####################################")
        print("")
        print("--verbose")
        print(f"--fil_A   {args.fil_A}")

        if getattr(args, "fil_B", None):
            print(f"--fil_B   {args.fil_B}")

        print(f"--skp_pfx {skp_pfx}")
        print(f"--method  {args.method}")

        if args.method == "edger":
            print(f"--normalization {args.normalization}")
            print(f"--prior_count {args.prior_count}")

            if siz_bin is None and args.siz_bin is None:
                print("--siz_bin (unset)")
            elif siz_bin is None:
                print(f"--siz_bin {args.siz_bin}")
            elif args.siz_bin is None:
                print(f"--siz_bin {siz_bin}  ## inferred from track ##")
            else:
                print(f"--siz_bin {siz_bin}")

            if args.lib_A is not None:
                print(f"--lib_A   {args.lib_A}")

            if args.lib_B is not None:
                print(f"--lib_B   {args.lib_B}")

            if args.sf_A is not None:
                print(f"--sf_A    {args.sf_A}")

            if args.sf_B is not None:
                print(f"--sf_B    {args.sf_B}")

            if args.frg_A is not None:
                print(f"--frg_A   {args.frg_A}")

            if args.frg_B is not None:
                print(f"--frg_B   {args.frg_B}")

            if args.prt_arg:
                print("--prt_arg")
        else:
            if args.method == "qntl_nz":
                print(f"--qntl_nz {args.qntl_nz}")

            if coef_eff is not None and coef_eff != args.coef:
                print(f"--coef    {args.coef}  ## coef_eff = {coef_eff} ##")
            else:
                print(f"--coef    {args.coef}")

            print(f"--floor   {args.floor}")
            print(f"--eps     {args.eps}")
            print(f"--mode_nz {args.mode_nz}")
            print(f"--sym     {args.sym}")

        print(f"--dp      {args.dp}")

        if args.prt_jsn:
            print("--prt_jsn")

        print("")
        print("")


def _is_one_track(args: argparse.Namespace) -> bool:
    """
    Report whether the edgeR request describes a single track.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments.

    Returns
    -------
    one_track : bool
        True when neither '--fil_B' nor '--lib_B' was supplied.

    Notes
    -----
    edgeR itself has no two-sample requirement, so neither does this. In
    'add_prior_count.c:88' 'compute_offsets' averages the library sizes over
    all columns and scales each sample's prior by 'offset[lib] / ave_lib'.
    With one column that ratio is exactly 1: the per-sample scaling degenerates
    to a no-op, the prior stays at its nominal 'prior.count', and the
    denominator becomes 'L + 2 * prior.count'. Nothing guards the path --
    'cpm.default' checks only for a zero-length dimension, and 'rpkm.default'
    is 'cpm.default' followed by a length division, so it inherits the same
    behavior.

    Reproducing that here therefore needs no separate estimator: passing the
    one library size as both 'lib_a' and 'lib_b' makes 'L_bar' equal 'L_A',
    which is exactly what 'ave_lib' becomes when 'nlib' is 1. Confirmed against
    edgeR 4.4.0 and pinned by
    'test_compute_pseudo_edger_reproduces_edger_for_one_library'.

    The consequence a caller must know: a track's one-track pseudocount is not
    its two-track pseudocount, because 'L_bar' differs. That is edgeR's own
    behavior -- 'cpm(log = TRUE)' on a one-column 'DGEList' differs from the
    same column inside a two-column one -- and not an artifact here.
    """

    return not getattr(args, "fil_B", None) and args.lib_B is None


def _resolve_ignored(token: str) -> str | None:
    """
    Report which ignored option a supplied token names, if any.

    Parameters
    ----------
    token : str
        One command-line token as supplied.

    Returns
    -------
    flag : str | None
        The option's long form, or None when the token names something else.

    Notes
    -----
    Three spellings reach the same option: the bare flag, an inline value after
    '=', and, for a single-character short option, a value attached directly as
    in '-c0.05'. Only the third needs care, because a short option is a prefix
    of longer ones -- '-s' of '-sp', '-sb', '-sfA', and '-sfB' -- so the token
    resolves to the longest registered option that prefixes it, which is how
    argparse itself decides. Abbreviation is off ('allow_abbrev=False'), so no
    partial long form has to be recognized.
    """

    name = token.split("=", 1)[0]

    for spellings in OPT_IGNORED_EDGER.values():
        if name in spellings:
            return spellings[-1]

    if not token.startswith("-") or token.startswith("--"):
        return None

    prefixes = [opt for opt in OPT_SHORT_ALL if token.startswith(opt)]

    if not prefixes:
        return None

    longest = max(prefixes, key=len)

    for spellings in OPT_IGNORED_EDGER.values():
        if longest in spellings:
            return spellings[-1]

    return None


def _warn_inapplicable(
    args: argparse.Namespace,
    argv: list[str] | None,
) -> None:
    """
    Warn when arguments that do not apply to the chosen method are passed.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments.
    argv : list[str] | None
        Arguments as supplied, or None to read the process arguments.

    Notes
    -----
    'edger' derives its pseudocount from library sizes, so it consults none of
    the distribution-shaping options. Silently ignoring them lets a wrong
    mental model survive: '--coef 0.05 --method edger' currently runs, ignores
    the coefficient, and reports nothing.

    Detection reads the supplied tokens rather than comparing against defaults,
    so an explicitly passed default is still reported. It resolves every
    spelling argparse accepts, short and long alike; reading long forms only
    would miss '-c 0.05', which is the very case the note exists for.
    """

    supplied = sys.argv[1:] if argv is None else argv
    seen: list[str] = []

    for token in supplied:
        flag = _resolve_ignored(token)

        if flag is not None and flag not in seen:
            seen.append(flag)

    if seen:
        print(
            f"Note: {', '.join(seen)} do not apply to '--method edger' and "
            "were ignored; it derives the pseudocount from library sizes, not "
            "from the value distribution.",
            file=sys.stderr,
        )


def _run_edger(
    args: argparse.Namespace,
    skp_pfx: tuple[str, ...],
) -> int:
    """
    Emit scale factors and pseudocounts from edgeR's prior rule.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments carrying the edgeR policy.
    skp_pfx : tuple[str, ...]
        Header prefixes to skip when reading library sizes from tracks.

    Returns
    -------
    status : int
        Zero. Every failure leaves through 'SystemExit' instead.

    Raises
    ------
    SystemExit
        For a nonpositive library size, an unusable normalization request, an
        unreadable track, a '--siz_bin' the track contradicts, or a
        '--normalization RPKM' run with no width and no track to infer one
        from.

    Notes
    -----
    Library sizes come from '--lib_A' and '--lib_B' when supplied, and are
    otherwise summed from the tracks. Summing requires a non-normalized track,
    as a CPM or RPKM bedGraph sums to a normalized total, not to a library
    size, and would silently rescale every value this function returns.

    Reading a track also resolves the bin width, so '--siz_bin' is needed only
    when both library sizes are supplied, which reads no track, and
    '--normalization RPKM' still wants a width for its scale factor.

    The result line reaches stdout:
      - the pseudocount pair as 'A:B',
      - a single value in single-track mode, or
      - a deepTools argument string under '--prt_arg'.

    The pair is the spec 'compute_signal_ratio' accepts for '--pseudo' and
    '--scl_fct'. A JSON summary follows on stdout under '--prt_jsn'.
    """

    one_track = _is_one_track(args)

    lib_a = args.lib_A
    lib_b = args.lib_B
    frg_a, frg_b = args.frg_A, args.frg_B
    sf_a, sf_b = args.sf_A, args.sf_B

    siz_bin = args.siz_bin

    # Reading a track resolves the bin width, so a run that recovers either
    # library size needs no '--siz_bin' at all. Only a run given both sizes
    # reads nothing, and only 'RPKM' then still needs a width.
    try:
        if lib_a is None:
            counts_a = sum_counts_bdg(args.fil_A, siz_bin, skp_pfx)
            lib_a, siz_bin = counts_a.total, counts_a.siz_bin

        if one_track:
            # Mirror A onto B so 'L_bar' collapses to 'L_A'. That reproduces
            # edgeR's 'nlib == 1' behavior exactly rather than approximating
            # it; see '_is_one_track' for the source reading it comes from.
            lib_b = lib_a
            frg_b = frg_a
            sf_b = sf_a
        elif lib_b is None:
            counts_b = sum_counts_bdg(args.fil_B, siz_bin, skp_pfx)
            lib_b, siz_bin = counts_b.total, counts_b.siz_bin
    except (OSError, ValueError) as e:
        if args.verbose:
            _print_pseudo_arguments(args, None, skp_pfx, siz_bin)

        raise SystemExit(str(e)) from None

    if args.verbose:
        _print_pseudo_arguments(args, None, skp_pfx, siz_bin)

    if siz_bin is None and canonicalize_norm(args.normalization) == "RPKM":
        raise SystemExit(
            "'--siz_bin' is required for '--normalization RPKM' when both "
            "'--lib_A' and '--lib_B' are supplied, because no track is read "
            "to infer the bin width from.",
        )

    try:
        result = compute_pseudo_edger(
            lib_a=lib_a,
            lib_b=lib_b,
            prior_count=args.prior_count,
            norm=args.normalization,
            siz_bin=siz_bin,
            scale_a=sf_a,
            scale_b=sf_b,
            frg_a=frg_a,
            frg_b=frg_b,
        )
    except ValueError as e:
        raise SystemExit(str(e)) from None

    if not result["is_edger"]:
        print(
            f"Note: '--normalization {args.normalization}' does not "
            f"reproduce edgeR's estimator: {result['note']}.",
            file=sys.stderr,
        )

    def render(value: float) -> str:
        return format_value(value, args.dp)

    scale_a = result["scale_A"]
    scale_b = result["scale_B"]
    pseudo_a = result["pseudo_A"]
    pseudo_b = result["pseudo_B"]

    if args.verbose:
        with redirect_stdout(sys.stderr):
            print(f"lib_A          {format_value(lib_a, args.dp)}")
            print(f"lib_B          {format_value(lib_b, args.dp)}")
            print(
                "prior_scaled_A "
                f"{format_value(result['prior_scaled_A'], args.dp)}",
            )
            print(
                "prior_scaled_B "
                f"{format_value(result['prior_scaled_B'], args.dp)}",
            )

            if "k_A" in result:
                print(f"k_A            {format_value(result['k_A'], args.dp)}")
                print(f"k_B            {format_value(result['k_B'], args.dp)}")
            print(f"is_edger       {result['is_edger']}")
            print("")

    if args.prt_arg:
        print(
            f"--scaleFactors {scale_a:.{args.dp}f}:{scale_b:.{args.dp}f} "
            f"--pseudocount {render(pseudo_a)} {render(pseudo_b)}",
        )
    elif one_track:
        print(render(pseudo_a))
    else:
        print(f"{render(pseudo_a)}:{render(pseudo_b)}")

    if args.prt_jsn:
        out = {
            "fil_A": args.fil_A,
            "fil_B": getattr(args, "fil_B", None),
            "method": "edger",
            "params": {
                "normalization": args.normalization,
                "prior_count": args.prior_count,
                "siz_bin": siz_bin,
                "dp": args.dp,
                "skp_pfx": list(skp_pfx),
            },
            "lib_sizes": {"A": lib_a, "B": lib_b},
            "k": (
                {"A": result["k_A"], "B": result["k_B"]}
                if "k_A" in result
                else None
            ),
            # Not derivable under 'norm': the pseudocount is symmetric there, so
            # 'pseudo_i / scale_i' returns it rather than the prior.
            "prior_scaled": {
                "A": result["prior_scaled_A"],
                "B": result["prior_scaled_B"],
            },
            "scale_factors": {"A": scale_a, "B": scale_b},
            "pseudocounts": {
                "pseudo_A": pseudo_a,
                "pseudo_B": pseudo_b,
                "pseudo_A_str": render(pseudo_a),
                "pseudo_B_str": render(pseudo_b),
            },
            "is_edger": result["is_edger"],
            "note": result["note"],
            # The B fields mirror A here rather than being dropped, so the
            # schema does not change shape between one- and two-track runs.
            # 'one_track' is what tells a consumer why they are equal.
            "one_track": one_track,
        }

        try:
            print(json.dumps(out, separators=(",", ":"), allow_nan=False))
        except ValueError:
            print(
                "Strict JSON disallows nan and inf; check '--lib_A' and "
                "'--lib_B', or just skip '--prt_jsn'.",
                file=sys.stderr,
            )

    return 0


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
        Zero on success after printing one pseudocount, an 'A:B' pair, and
        optionally a single-line JSON summary.

    Raises
    ------
    SystemExit
        For help, parser rejection, invalid numeric arguments, or missing and
        conflicting input paths.

    Notes
    -----
    The four distribution-based methods filter inputs through '--eps' and
    '--mode_nz' before computing, and can symmetrize the pair with '--sym';
    '--method edger' consults none of those, deriving the pseudocount from
    library sizes instead. Malformed or nonnumeric bedGraph rows are skipped by
    the row iterator. With '--prt_jsn', the command prints a strict one-line
    JSON summary only when all serialized values are finite; otherwise it warns
    on stderr, omits the JSON line, and returns zero. Warnings about empty
    inputs, zero pseudocounts, or symmetrization are written to stderr.
    """

    args = parse_args(argv)

    paths = [
        p for p in (args.fil_A, getattr(args, "fil_B", None)) if p is not None
    ]

    try:
        ensure_single_stdin(paths)
    except ValueError as e:
        raise SystemExit(str(e)) from None

    try:
        for label, p in (
            ("A", args.fil_A),
            ("B", getattr(args, "fil_B", None)),
        ):
            if p is None or p == "-":
                continue

            check_exists(p, kind="file", label=f"bedGraph {label}")
    except FileNotFoundError as e:
        raise SystemExit(str(e)) from None

    try:
        if args.method == "qntl_nz" and (
            not math.isfinite(args.qntl_nz)
            or not (0.0 <= args.qntl_nz <= 100.0)
        ):
            raise ValueError("'--qntl_nz' must be finite and in [0, 100].")

        if args.method == "edger":
            one_track = _is_one_track(args)

            validate_comparison(
                args.prior_count,
                "ge",
                0.0,
                "prior_count",
                allow_none=False,
            )
            validate_comparison(
                args.siz_bin,
                "gt",
                0,
                "siz_bin",
                allow_none=True,
            )

            if one_track:
                need_sf = args.sf_A is None
                need_aln = args.frg_A is None
                both = "'--sf_A'"
                both_frg = "'--frg_A'"
            else:
                need_sf = args.sf_A is None or args.sf_B is None
                need_aln = args.frg_A is None or args.frg_B is None
                both = "both '--sf_A' and '--sf_B'"
                both_frg = "both '--frg_A' and '--frg_B'"

            if canonicalize_norm(args.normalization) == "norm" and need_aln:
                raise ValueError(
                    f"'--normalization {args.normalization}' requires "
                    f"{both_frg}; a normalized-coverage track sums to 1, so "
                    "the fragment count cannot be recovered from it."
                )

            if args.normalization == "RPGC" and need_sf:
                raise ValueError(
                    f"'--normalization RPGC' requires {both}; read them from "
                    "'bamCoverage --verbose'.",
                )

            if one_track and args.prt_arg:
                # The deepTools 'bamCompare' options '--scaleFactors' and
                # '--pseudocount' take a pair each, so there is no one-track
                # spelling of them. 'bamCoverage --scaleFactor' exists but
                # accepts no pseudocount, so emitting it would silently drop
                # the value that was asked for.
                raise ValueError(
                    "'--prt_arg' writes the two-track 'bamCompare' argument "
                    "string, which has no single-track form. Drop "
                    "'--prt_arg' to print the pseudocount, or supply track B."
                )

        validate_comparison(args.coef, "ge", 0.0, "coef", allow_none=True)
        validate_comparison(args.floor, "ge", 0.0, "floor", allow_none=False)
        validate_comparison(args.eps, "ge", 0.0, "eps", allow_none=False)
        validate_comparison(args.dp, "ge", 0, "dp", allow_none=False)
    except ValueError as e:
        raise SystemExit(str(e)) from None

    skp_pfx = parse_skp_pfx(args.skp_pfx, default=DEF_SKP_PFX)

    if args.method == "edger":
        _warn_inapplicable(args, argv)

        return _run_edger(args, skp_pfx)

    coef_eff = determine_coef_eff(args.method, args.coef)

    mode_nz = args.mode_nz

    if args.verbose:
        _print_pseudo_arguments(args, coef_eff, skp_pfx)

    vals_a = list(
        iter_vals_bdg(
            args.fil_A,
            eps=args.eps,
            mode_nz=mode_nz,
            skp_pfx=skp_pfx,
        ),
    )

    if not vals_a:
        print(
            "No finite values in A after filtering; pseudocount will be "
            "'nan'. Check '--eps', '--mode_nz', and/or '--skp_pfx'.",
            file=sys.stderr,
        )

    pseudo_a = pick_stabilizer(
        vals_a,
        method=args.method,
        coef=coef_eff,
        qntl_pct=args.qntl_nz,
        floor=args.floor,
        qntl_rule="round",
    )

    pseudo_b = float("nan")
    vals_b = []

    if getattr(args, "fil_B", None):
        vals_b = list(
            iter_vals_bdg(
                args.fil_B,
                eps=args.eps,
                mode_nz=mode_nz,
                skp_pfx=skp_pfx,
            ),
        )

        if not vals_b:
            print(
                "No finite values in B after filtering; pseudocount will be "
                "'nan'. Check '--eps', '--mode_nz', and/or '--skp_pfx'.",
                file=sys.stderr,
            )

        pseudo_b = pick_stabilizer(
            vals_b,
            method=args.method,
            coef=coef_eff,
            qntl_pct=args.qntl_nz,
            floor=args.floor,
            qntl_rule="round",
        )
    else:
        # If only one file is supplied, then mirror A to B if and only if the
        # user explicitly requested a non-'none' '--sym' rule.
        if (args.sym != "none") and math.isfinite(pseudo_a):
            pseudo_b = pseudo_a

        # Most symmetry modes with only fil_A are mathematically degenerate.
        file_b_missing = getattr(args, "fil_B", None) is None
        symmetry_requires_b = args.sym not in {"none", "use_A"}

        if file_b_missing and symmetry_requires_b:
            print(
                f"Note: '--sym {args.sym}' was used with only '--fil_A'; "
                "symmetrization degenerates to 'pseudo_A:pseudo_A'.",
                file=sys.stderr,
            )

    # Warn if a computed pseudo is exactly 0 (can cause -inf in log2 ratios).
    for tag, pseudo in (("A", pseudo_a), ("B", pseudo_b)):
        if math.isfinite(pseudo) and pseudo == 0.0:
            print(
                f"Pseudocount for {tag} is 0.0; log2 ratios may produce -inf "
                "at zeros. Consider a positive '--floor' or a larger "
                "'--coef'.",
                file=sys.stderr,
            )

    if args.sym != "none":
        pseudo_a, pseudo_b = combine_pseudo_sym(
            pseudo_a,
            pseudo_b,
            mode=args.sym,
        )

    # Symmetrization can introduce new zero values.
    for tag, pseudo in (("A", pseudo_a), ("B", pseudo_b)):
        if math.isfinite(pseudo) and pseudo == 0.0:
            print(
                f"Pseudocount for {tag} is 0.0 after symmetrization; log2 "
                "ratios may produce -inf at zeros. Consider a positive "
                "'--floor' or larger '--coef'.",
                file=sys.stderr,
            )

    want_pair = (getattr(args, "fil_B", None) is not None) or (
        args.sym != "none"
    )

    def format_pseudocount(value: float) -> str:
        return format_value(value, args.dp)

    if want_pair:
        rendered_a = format_pseudocount(pseudo_a)
        rendered_b = format_pseudocount(
            pseudo_b if math.isfinite(pseudo_b) else pseudo_a,
        )
        print(f"{rendered_a}:{rendered_b}")
    else:
        print(f"{format_pseudocount(pseudo_a)}")

    # Keep the JSON summary local until other CLIs share a reviewed contract.
    if args.prt_jsn:
        pseudocounts = {
            "pseudo_A": pseudo_a,
            "pseudo_A_str": format_pseudocount(pseudo_a),
        }

        if want_pair:
            value_b = pseudo_b if math.isfinite(pseudo_b) else pseudo_a
            pseudocounts.update(
                {
                    "pseudo_B": value_b,
                    "pseudo_B_str": format_pseudocount(value_b),
                },
            )

        params = {"coef": args.coef}

        if coef_eff is not None:
            params["coef_eff"] = coef_eff

        params.update(
            {
                "qntl_nz": args.qntl_nz,
                "floor": args.floor,
                "eps": args.eps,
                "mode_nz": mode_nz,
                "sym": args.sym,
                "dp": args.dp,
                "skp_pfx": list(skp_pfx),
            },
        )

        out = {
            "fil_A": args.fil_A,
            "fil_B": getattr(args, "fil_B", None),
            "method": args.method,
            "params": params,
            "stats": {
                "A": compute_stats_robust(vals_a),
                "B": compute_stats_robust(vals_b) if vals_b else None,
            },
            "pseudocounts": pseudocounts,
        }

        # Strict JSON rejects nonfinite values, so use a stable fallback.
        try:
            print(json.dumps(out, separators=(",", ":"), allow_nan=False))
        except ValueError:
            print(
                "Strict JSON disallows nan and inf; adjust '--floor' and "
                "'--coef', or just skip '--prt_jsn'.",
                file=sys.stderr,
            )

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
