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
Compute pseudocounts for bedGraph signal tracks.

The CLI accepts input tracks and a method. '--method edger' derives the
pseudocount from fragment-bin overlap counts, taking a target signal type, a
prior count, and optionally the counts themselves; the four distribution-based
methods derive it from the value distribution, taking a coefficient, nonzero
filtering, and a symmetry mode. It prints one pseudocount or an A:B pair to
stdout, and with '--prt_jsn' also prints a JSON summary.

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
    canonicalize_typ_sig,
    compute_pseudo_edger,
    compute_stats_robust,
    determine_coef_eff,
    iter_vals_bdg,
    pick_stabilizer,
)

with suppress(AttributeError, ValueError):
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

assert sys.version_info >= (3, 11), "Python >= 3.11 required."

# The signal types this tool serves: the project's own, as 'compute_signal'
# emits them. 'TYP_SIG_CANON' spans both tools, so each CLI restricts its own
# choices rather than carrying its own vocabulary. The deepTools signal types
# belong to 'compute_pseudo_deeptools'.
TYP_SIG_SERVED = ("unadj", "frag", "norm", "nc", "count", "cpm")

# Signal types whose prior is denominated against the fragment count, so they
# need '--n_frg_A' and '--n_frg_B'. 'unadj' additionally needs its own column
# sum, which the track supplies and no report flag writes.
SUB_FRACTIONAL = ("unadj", "frag", "norm")

# Signal types this tool no longer accepts and where they went. A bare argparse
# rejection would leave a user guessing which of the two tools to reach for.
SUB_MOVED = ("CPM", "BPM", "RPKM", "RPGC", "None")

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
# the whole set rather than the six above: '-sb10' must resolve to '-sb', not
# to '-s' carrying 'b10', and only the longest registered prefix tells the two
# apart. 'test_pseudo.py' asserts both constants against the parser, so neither
# can drift away from it.
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
    "-ts",
    "-pc",
    "-sb",
    "-noA",
    "-noB",
    "-nfA",
    "-nfB",
    "-dp",
    "-pj",
)


# TODO: Extend compressed-input handling to '.bgz' and '.bgzf' here and in
# related bedGraph parsers. Revisit JSON summarization and whether a shared
# warning/error-text helper should be exported to other Python CLIs.
def _check_typ_sig(value: str) -> str:
    """
    Accept a served signal type, naming the other tool for a moved one.

    Parameters
    ----------
    value : str
        One '--typ_sig' value as supplied.

    Returns
    -------
    typ_sig : str
        The value unchanged, for the 'choices' check to validate.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value names a deepTools signal type, which moved to
        'compute_pseudo_deeptools' when the tools split.

    Notes
    -----
    argparse applies 'type' before 'choices', so a moved name is intercepted
    here and never reaches the bare rejection. Leaving it to 'choices' instead
    would either list the moved names in the advertised set or drop the
    migration hint, and the split is recent enough that a user will hit it.
    """

    if value in SUB_MOVED:
        raise argparse.ArgumentTypeError(
            f"'{value}' is a deepTools signal type; use "
            f"'compute_pseudo_deeptools', which serves {', '.join(SUB_MOVED)} "
            "through the supplied-scale-factor workflow.",
        )

    return value


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
                    "typ_sig",
                    "prior_count",
                    "siz_bin",
                    "n_ovlp_A",
                    "n_ovlp_B",
                    "n_frg_A",
                    "n_frg_B",
                ),
                ("dp", "prt_jsn"),
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
        help=(
            "First bedGraph input file, file A. Use '-' for stdin; '.gz' is "
            "handled.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "--fil-A",
        dest="fil_A",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-fB",
        "--fil_B",
        dest="fil_B",
        default=None,
        help=(
            "Second bedGraph input file, file B. This file is optional for "
            "single-file pseudocount modes; use '-' for stdin; '.gz' is "
            "handled.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "--fil-B",
        dest="fil_B",
        help=argparse.SUPPRESS,
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
        "--skp-pfx",
        dest="skp_pfx",
        type=str,
        help=argparse.SUPPRESS,
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
            "\n"
            "| method     | description                             |\n"
            "| :---       | :---                                    |\n"
            "| edger      | edgeR's prior.count rule                |\n"
            "| frc_mdn_nz | value = coef × median of nonzero bins   |\n"
            "| qntl_nz    | value = q-th percentile of nonzero bins |\n"
            "| frc_avg_nz | value = coef × mean of nonzero bins     |\n"
            "| min_nz     | value = coef × minimum nonzero bin      |\n"
            "\n"
            "Notes:\n"
            "  - '--method edger' takes one track or two, as edgeR does.\n"
            "  - With one track, edgeR's per-sample prior scaling degenerates "
            "to a no-op and the mean fragment-bin overlap count is that "
            "track's own.\n"
            "  - For a given track, its single-track pseudocount is not its "
            "two-track pseudocount; both are correct for their own frame.\n"
            "  - If '--method qntl_nz', set percentile with '--qntl_nz' "
            "[decimals OK (e.g., 0.1 = 0.1th percentile); nearest-rank "
            "determined via 'round'].\n"
            "  - 'nonzero' means '|x| > eps' with '--mode_nz closed', "
            "'|x| >= eps' with '--mode_nz open', and every finite value with "
            "'--mode_nz off'.\n"
            "  - If '--coef' is omitted, then defaults to '--coef 0.01' for "
            "'--method frc_*' and '--coef 1.0' for '--method min_nz'.\n"
            "  - '--method min_nz' typically needs a larger coef (e.g., "
            "0.1–1.0) in comparison to '--method frc_*' (e.g., 0.01).\n"
            "\n"
            "For '--method edger', see '--typ_sig' for more details.\n"
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
        "--qntl-nz",
        dest="qntl_nz",
        type=float,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-c",
        "--coef",
        dest="coef",
        type=float,
        default=None,
        help=(
            "Coefficient for median, mean, and min methods. If not specified, "
            "then defaults to 0.01 for '--method frc_mdn_nz' and "
            "'--method frc_avg_nz', or 1.0 for '--method min_nz'.\n"
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
            "treated as zero and excluded from statistics; with "
            "'--mode_nz open', values with '|x| < ε' are excluded; with "
            "'--mode_nz off', ε-based filtering is disabled.\n"
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
            "\n"
            "|                  | description                         |\n"
            "| :---             | :---                                |\n"
            "| closed (default) | drop values with '|x| <= eps'.      |\n"
            "| open             | drop values with '|x| < eps'.       |\n"
            "| off              | disable 'eps'-based zero filtering. |\n"
            "\n"
        ),
    )
    parser.add_argument(
        "--mode-nz",
        dest="mode_nz",
        choices=("closed", "open", "off"),
        help=argparse.SUPPRESS,
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
            "  - none returns A and B as computed; max/min select the "
            "larger/smaller value; arith/geom/harm use "
            "arithmetic/geometric/harmonic means; use_A/use_B copy A/B to "
            "both outputs.\n"
            "  - With finite values, geom falls back to min for a negative "
            "input and harm falls back to min for a nonpositive input; each "
            "fallback writes a warning to stderr.\n"
            "  - In any non-none mode, one finite value is mirrored with a "
            "warning; two nonfinite values are retained with a warning.\n"
            "  - If both A AND B are given, apply the chosen rule to "
            "'(pseudo_A, pseudo_B)', then print 'pseudo_A:pseudo_B'.\n"
            "  - If only A is given AND '--sym' is provided AND NOT 'none', "
            "mirror A to B and print: 'pseudo_A:pseudo_B'.\n"
            "  - If only A is given AND '--sym' is omitted OR '--sym none', "
            "print a single value: 'pseudo_A'.\n"
            "\n"
        ),
    )

    parser.add_argument(
        "-ts",
        "--typ_sig",
        dest="typ_sig",
        type=_check_typ_sig,
        choices=TYP_SIG_SERVED,
        default="norm",
        help=(
            "Target signal type: which '--method' wrote the track (default: "
            "%(default)s). Applies to '--method edger' only; ignored "
            "otherwise.\n"
            "\n"
            "|       | deposition          | column sum  | n_frg | n_ovlp |\n"
            "| :---  | :---                | :---        | :---  | :---   |\n"
            "| unadj | overlap in bp       | frag bp     | yes   | yes    |\n"
            "| frag  | overlap / length    | frag count  | yes   | infer  |\n"
            "| norm  | overlap / len / 'N' | one         | yes   | infer  |\n"
            "| count | one per touched bin | overlaps    | no    | infer  |\n"
            "| cpm   | touches x 1e6 / 'L' | one million | no    | infer  |\n"
            "\n"
            "For fractional signal types 'unadj', 'frag', and 'norm', a "
            "partly covered bin takes a share of the overlap; for whole-count "
            "signal types 'count' and 'cpm', a partly covered bin takes one "
            "count.\n"
            "\n"
            "Take the counts the last two columns ask for from "
            "'compute_signal --report_n_frg --report_n_ovlp', which writes "
            "them beside the track it is already producing.\n"
            "\n"
            "Where that column reads 'infer', '--n_ovlp_A' and '--n_ovlp_B' "
            "may be omitted and the tool sums '--fil_A' and '--fil_B' in "
            "their place. Only a 'count' track sums to 'L', so 'count' tracks "
            "must be used for inference regardless of chosen signal type.\n"
            "\n"
            "Here, 'n_ovlp' ('L') adds up, across all fragments, how many "
            "bins each fragment spans (it is edgeR's 'lib.size').\n"
            "\n"
            "Additional notes:\n"
            "  - 'nc' is an alias for 'norm'.\n"
            "  - For deepTools signal types, use 'compute_pseudo_deeptools'.\n"
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
            "edgeR 'prior.count', before the per-sample scaling each signal "
            "type applies (default: %(default)s).\n"
            "\n"
            "Applies to '--method edger' only; ignored otherwise.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "--prior-count",
        dest="prior_count",
        type=float,
        help=argparse.SUPPRESS,
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
            "rather than silently rescaling the fragment-bin overlap count.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "--siz-bin",
        dest="siz_bin",
        type=int,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-noA",
        "--n_ovlp_A",
        dest="n_ovlp_A",
        type=float,
        default=None,
        help=(
            "Fragment-bin overlap count 'L' for track A: the column sum of "
            "the bin matrix, which is edgeR's 'lib.size'. It adds up, across "
            "all fragments, how many bins each fragment spans (note this is "
            "not the fragment count, which is '--n_frg_A').\n"
            "\n"
            "Written by 'compute_signal --report_n_ovlp', as "
            "'<track>.n_ovlp.txt' when that flag is given without a path. "
            "Computed from '--fil_A' when omitted, which reads the track's "
            "column sum as 'L' and so holds only for a 'count' track; a sum "
            "below '--n_frg_A' is impossible for one and is refused.\n"
            "\n"
            "Applies to '--method edger' only; ignored otherwise.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "--n-ovlp-A",
        "--n_ovlp-A",
        "--n-ovlp_A",
        dest="n_ovlp_A",
        type=float,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-noB",
        "--n_ovlp_B",
        dest="n_ovlp_B",
        type=float,
        default=None,
        help=(
            "Fragment-bin overlap count 'L' for track B; see '--n_ovlp_A'. "
            "Omit this and '--fil_B' to compute a single-track pseudocount.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "--n-ovlp-B",
        "--n_ovlp-B",
        "--n-ovlp_B",
        dest="n_ovlp_B",
        type=float,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-nfA",
        "--n_frg_A",
        dest="n_frg_A",
        type=float,
        default=None,
        help=(
            "Fragment count for track A: the number of fragments "
            "'compute_signal' divided by, not the number of alignment "
            "records. The two coincide only when exactly one alignment per "
            "fragment survives filtering. A fractional track's own column "
            "total is a deposition total rather than a count, so it cannot "
            "supply this.\n"
            "\n"
            "Written by 'compute_signal --report_n_frg', as "
            "'<track>.n_frg.txt' when that flag is given without a path.\n"
            "\n"
            "Applies to '--method edger' with a fractional signal type "
            "('unadj', 'frag', 'norm', and its alias 'nc') only, where it is "
            "required.\n"
        ),
    )
    parser.add_argument(
        "--n-frg-A",
        "--n_frg-A",
        "--n-frg_A",
        dest="n_frg_A",
        type=float,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-nfB",
        "--n_frg_B",
        dest="n_frg_B",
        type=float,
        default=None,
        help=("Fragment count for track B; see '--n_frg_A'.\n\n"),
    )
    parser.add_argument(
        "--n-frg-B",
        "--n_frg-B",
        "--n-frg_B",
        dest="n_frg_B",
        type=float,
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
        "--prt-jsn",
        dest="prt_jsn",
        action="store_true",
        help=argparse.SUPPRESS,
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
            print(f"--typ_sig {args.typ_sig}")
            print(f"--prior_count {args.prior_count}")

            if siz_bin is None and args.siz_bin is None:
                print("--siz_bin (unset)")
            elif siz_bin is None:
                print(f"--siz_bin {args.siz_bin}")
            elif args.siz_bin is None:
                print(f"--siz_bin {siz_bin}  ## inferred from track ##")
            else:
                print(f"--siz_bin {siz_bin}")

            if args.n_ovlp_A is not None:
                print(f"--n_ovlp_A {args.n_ovlp_A}")

            if args.n_ovlp_B is not None:
                print(f"--n_ovlp_B {args.n_ovlp_B}")

            if args.n_frg_A is not None:
                print(f"--n_frg_A {args.n_frg_A}")

            if args.n_frg_B is not None:
                print(f"--n_frg_B {args.n_frg_B}")

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
        True when neither '--fil_B' nor '--n_ovlp_B' was supplied.

    Notes
    -----
    edgeR itself has no two-sample requirement, so neither does this. In
    'add_prior_count.c:88' 'compute_offsets' averages the fragment-bin overlap
    counts over all columns and scales each sample's prior by
    'offset[lib] / ave_lib'. With one column that ratio is exactly 1: the
    per-sample scaling degenerates to a no-op, the prior stays at its nominal
    'prior.count', and the denominator becomes 'L + 2 * prior.count'. Nothing
    guards the path: 'cpm.default' checks only for a zero-length dimension, and
    'rpkm.default' is 'cpm.default' followed by a length division, so it
    inherits the same behavior.

    Thus, reproducing that here needs no separate estimator: passing the one
    fragment-bin overlap count as both 'n_ovlp_a' and 'n_ovlp_b' makes 'L_bar'
    equal 'L_A', which is exactly what 'ave_lib' becomes when 'nlib' is 1.
    Confirmed against edgeR 4.4.0, re-confirmed against 4.8.2, and pinned by
    'test_compute_pseudo_edger_reproduces_edger_for_one_track'.

    The consequence a caller must know: a track's single-track pseudocount is
    not its two-track pseudocount, because 'L_bar' is that track's own count in
    one mode and the mean of both in the other. That is edgeR's own behavior,
    not an artifact here: 'cpm(log = TRUE)' on a one-column 'DGEList' differs
    from the same column inside a two-column one.
    """

    return not getattr(args, "fil_B", None) and args.n_ovlp_B is None


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
    of longer ones ('-s' of '-sp', '-su', and '-sb'), so the token resolves to
    the longest registered option that prefixes it, which is how argparse
    itself decides. Abbreviation is off ('allow_abbrev=False'), so no partial
    long form has to be recognized.
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
    edgeR derives its pseudocount from fragment-bin overlap counts, so it
    consults none of the distribution-shaping options. Silently ignoring them
    lets a wrong mental model survive: '--coef 0.05 --method edger' currently
    runs, ignores the coefficient, and reports nothing.

    Detection reads the supplied tokens rather than comparing against defaults,
    so an explicitly passed default is still reported. It resolves every
    spelling argparse accepts, short and long alike.
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
            "were ignored; it derives the pseudocount from fragment-bin "
            "overlap counts, not from the value distribution.",
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
        Header prefixes to skip when reading fragment-bin overlap counts from
        tracks.

    Returns
    -------
    status : int
        Zero. Every failure leaves through 'SystemExit' instead.

    Raises
    ------
    SystemExit
        For a nonpositive fragment-bin overlap count, an unusable signal type
        request, an unreadable track, a '--siz_bin' the track contradicts, a
        fractional signal type without both fragment counts, an 'unadj' run
        whose column total is absent or below the fragment count, or an
        inferred overlap count that falls below the fragment count.

    Notes
    -----
    Fragment-bin overlap counts come from '--n_ovlp_A' and '--n_ovlp_B' when
    supplied, and are otherwise summed from the tracks. Summing requires a
    whole-count track, since a fractional one sums to its own total rather than
    to 'L' and would silently rescale every value this function returns.

    Reading a track also resolves the bin width, so '--siz_bin' is never
    required: no signal type this tool serves needs a width of its own, and the
    option exists to cross-check the track rather than to supply a missing
    number.

    One result line reaches stdout: a single value in single-track mode, or the
    pseudocount pair as 'A:B'. The pair is the spec 'compute_signal_ratio'
    accepts for '--pseudo' and '--scl_fct'. A JSON summary follows on stdout
    under '--prt_jsn'.
    """

    one_track = _is_one_track(args)

    n_ovlp_a = args.n_ovlp_A
    n_ovlp_b = args.n_ovlp_B
    n_frg_a, n_frg_b = args.n_frg_A, args.n_frg_B

    siz_bin = args.siz_bin

    # Reading a track resolves the bin width, so a run that recovers either
    # overlap count needs no '--siz_bin' at all, and a run given both of them
    # reads nothing and needs none either.
    try:
        if n_ovlp_a is None:
            counts_a = sum_counts_bdg(args.fil_A, siz_bin, skp_pfx)
            n_ovlp_a, siz_bin = counts_a.total, counts_a.siz_bin

        if one_track:
            # Mirror A onto B so 'L_bar' collapses to 'L_A'. That reproduces
            # edgeR's 'nlib == 1' behavior exactly rather than approximating
            # it; see '_is_one_track' for the source reading it comes from.
            n_ovlp_b = n_ovlp_a
            n_frg_b = n_frg_a
        elif n_ovlp_b is None:
            counts_b = sum_counts_bdg(args.fil_B, siz_bin, skp_pfx)
            n_ovlp_b, siz_bin = counts_b.total, counts_b.siz_bin
    except (OSError, ValueError) as e:
        if args.verbose:
            _print_pseudo_arguments(args, None, skp_pfx, siz_bin)

        raise SystemExit(str(e)) from None

    # Inferring 'L' reads the track's column sum, which is 'L' only for a
    # whole-count track. Every fragment touches at least one bin, so a sum
    # below 'N' proves the track is something else, and inferring from it would
    # rescale the prior in silence. This tests impossibility, not plausibility:
    # a real whole-count track cannot trip this, but one summing too high still
    # slips through.
    for label, inferred, total, frg in (
        ("A", args.n_ovlp_A is None, n_ovlp_a, n_frg_a),
        ("B", args.n_ovlp_B is None and not one_track, n_ovlp_b, n_frg_b),
    ):
        if inferred and frg is not None and total < frg:
            raise SystemExit(
                f"Inferred '--n_ovlp_{label}' is {total}, below "
                f"'--n_frg_{label}' at {frg}; the summed track cannot be a "
                f"whole-count track, so '--n_ovlp_{label}' is required for "
                f"'--typ_sig {args.typ_sig}'.",
            )

    if args.verbose:
        _print_pseudo_arguments(args, None, skp_pfx, siz_bin)

    total_a, total_b = None, None

    if canonicalize_typ_sig(args.typ_sig) == "unadj":
        # The track's own column sum is the total fragment base pairs, which is
        # what the prior is denominated in. It is read here rather than taken
        # from '--n_ovlp_A', which carries the overlap count 'k' needs and is a
        # different quantity entirely.
        try:
            total_a = sum_counts_bdg(args.fil_A, siz_bin, skp_pfx).total
            total_b = (
                total_a
                if one_track
                else sum_counts_bdg(args.fil_B, siz_bin, skp_pfx).total
            )
        except (OSError, ValueError) as e:
            raise SystemExit(str(e)) from None

    try:
        result = compute_pseudo_edger(
            n_ovlp_a=n_ovlp_a,
            n_ovlp_b=n_ovlp_b,
            total_a=total_a,
            total_b=total_b,
            prior_count=args.prior_count,
            typ_sig=args.typ_sig,
            siz_bin=siz_bin,
            n_frg_a=n_frg_a,
            n_frg_b=n_frg_b,
        )
    except ValueError as e:
        raise SystemExit(str(e)) from None

    if not result["is_edger"]:
        print(
            f"Note: '--typ_sig {args.typ_sig}' does not reproduce edgeR's "
            f"estimator: {result['note']}.",
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
            print(f"n_ovlp_A        {format_value(n_ovlp_a, args.dp)}")
            print(f"n_ovlp_B        {format_value(n_ovlp_b, args.dp)}")
            prior_a = format_value(result["prior_scaled_A"], args.dp)
            prior_b = format_value(result["prior_scaled_B"], args.dp)

            print(f"prior_scaled_A {prior_a}")
            print(f"prior_scaled_B {prior_b}")

            if "k_A" in result:
                print(f"k_A            {format_value(result['k_A'], args.dp)}")
                print(f"k_B            {format_value(result['k_B'], args.dp)}")

            print(f"is_edger       {result['is_edger']}")
            print("")

    if one_track:
        print(render(pseudo_a))
    else:
        print(f"{render(pseudo_a)}:{render(pseudo_b)}")

    if args.prt_jsn:
        out = {
            "fil_A": args.fil_A,
            "fil_B": getattr(args, "fil_B", None),
            "method": "edger",
            "params": {
                "typ_sig": args.typ_sig,
                "prior_count": args.prior_count,
                "siz_bin": siz_bin,
                "dp": args.dp,
                "skp_pfx": list(skp_pfx),
            },
            "n_ovlp": {"A": n_ovlp_a, "B": n_ovlp_b},
            "k": (
                {"A": result["k_A"], "B": result["k_B"]}
                if "k_A" in result
                else None
            ),

            # Not derivable under 'unadj', 'frag', 'norm' or 'cpm': those
            # pseudocounts come from a closed form rather than from
            # 'scale_i * prior_scaled_i', so dividing does not recover the
            # prior. For that reason, emitted rather than left to the consumer.
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
                "Strict JSON disallows nan and inf; check '--n_ovlp_A' and "
                "'--n_ovlp_B', or just skip '--prt_jsn'.",
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
    fragment-bin overlap counts instead. Malformed or nonnumeric bedGraph rows
    are skipped by the row iterator. With '--prt_jsn', the command prints a
    strict one-line JSON summary only when all serialized values are finite;
    otherwise it warns on stderr, omits the JSON line, and returns zero.
    Warnings about empty inputs, zero pseudocounts, or symmetrization are
    written to stderr.
    """

    args = parse_args(argv)

    # A hidden hyphen spelling is a separate action that argparse's own
    # 'required' check cannot see, so '--fil-A' alone would be rejected for
    # want of '--fil_A'. Check the parsed value here instead.
    if getattr(args, "fil_A", None) is None:
        raise SystemExit(
            "'--fil_A' is required. Supply the first bedGraph input path.",
        )

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
                need_aln = args.n_frg_A is None
                both_frg = "'--n_frg_A'"
            else:
                need_aln = args.n_frg_A is None or args.n_frg_B is None
                both_frg = "both '--n_frg_A' and '--n_frg_B'"

            if canonicalize_typ_sig(args.typ_sig) == "unadj" and (
                args.n_ovlp_A is None
                or (not one_track and args.n_ovlp_B is None)
            ):
                # The track this run reads supplies the signal type's own
                # column total, the fragment base pairs. Summing it for
                # 'n_ovlp' too would put that total where 'L' belongs and scale
                # every prior by 'total / L' without saying so.
                raise ValueError(
                    "'--typ_sig unadj' requires '--n_ovlp_A' and "
                    "'--n_ovlp_B'; the track supplies its own column total, "
                    "not the fragment-bin overlap count that 'k' needs.",
                )

            if (
                canonicalize_typ_sig(args.typ_sig) in SUB_FRACTIONAL
                and need_aln
            ):
                raise ValueError(
                    f"'--typ_sig {args.typ_sig}' requires {both_frg}; the "
                    "fractional signal types are denominated against the "
                    "fragment count, which their own totals cannot supply.",
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
