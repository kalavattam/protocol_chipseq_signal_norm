#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: compute_pseudo_deeptools.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Compute pseudocounts for deepTools bedGraph tracks.

The CLI applies edgeR's 'prior.count' rule for a track written by deepTools,
taking the track's substrate, a prior count, and optionally the fragment-bin
overlap counts themselves. It prints one pseudocount or an A:B pair to stdout
or, under '--prt_arg', an argument string for use with 'bamCompare'; with
'--prt_jsn', it prints a JSON summary.

This is a deepTools-interop version of 'compute_pseudo' in which the arithmetic
is shared and the substrate vocabulary differs. For this project's own
substrates, use 'compute_pseudo'.

Examples
--------
python -m protocol_chipseq_signal_norm.cli.compute_pseudo_deeptools \\
    --fil_A <file> [--fil_B <file>] [options]
"""

from __future__ import annotations

import argparse
import json
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
    canonicalize_substrate,
    compute_pseudo_edger,
)

with suppress(AttributeError, ValueError):
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

assert sys.version_info >= (3, 11), "Python >= 3.11 required."

# These are the deepTools substrates this tool serves. 'SUBSTRATE_CANON' spans
# both tools, so each CLI restricts its own choices rather than carrying its
# own vocabulary. This tool's own substrates are the ones deepTools
# '--normalizeUsing' writes.
# TODO: Potential change of SUBSTRATE, SUB, etc. to something clearer, here and
# elsewhere.
SUB_CHOICES = ("CPM", "BPM", "RPKM", "RPGC", "None")


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
        Parsed bedGraph inputs and deepTools pseudocount policy.

    Raises
    ------
    SystemExit
        If argument parsing fails or help is requested.
    """

    parser = CapArgumentParser(
        description=(
            "Compute edgeR-equivalent pseudocount(s) for one (or two) "
            "deepTools bedGraph track(s)."
        ),
        prog="compute_pseudo_deeptools",
        _sectioned_help=_SectionedHelpConfig(
            usage_rows=(
                ("help", "verbose"),
                ("fil_A", "fil_B", "skp_pfx"),
                (
                    "substrate",
                    "prior_count",
                    "siz_bin",
                    "n_ovlp_A",
                    "n_ovlp_B",
                    "sf_A",
                    "sf_B",
                ),
                ("dp", "prt_jsn", "prt_arg"),
            ),
            examples=(
                _HelpExample(
                    description=(
                        "Compute a CPM pseudocount pair from two deepTools "
                        "tracks."
                    ),
                    command_lines=(
                        "compute_pseudo_deeptools",
                        "--fil_A signal_A.bedGraph",
                        "--fil_B signal_B.bedGraph",
                    ),
                ),
                _HelpExample(
                    description=(
                        "Write a ready-to-paste 'bamCompare' argument string "
                        "from supplied fragment-bin overlap counts."
                    ),
                    command_lines=(
                        "compute_pseudo_deeptools",
                        "--fil_A signal_A.bedGraph",
                        "--fil_B signal_B.bedGraph",
                        "--substrate RPKM",
                        "--n_ovlp_A 41318705",
                        "--n_ovlp_B 39204118",
                        "--siz_bin 10",
                        "--prt_arg",
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
            "single-track mode; use '-' for stdin; '.gz' is handled.\n"
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
        "-su",
        "--substrate",
        dest="substrate",
        choices=SUB_CHOICES,
        default="CPM",
        help=(
            "Target substrate: the deepTools '--normalizeUsing' the track was "
            "written with (default: %(default)s).\n"
            "\n"
            "|          | description                                    |\n"
            "| :---     | :---                                           |\n"
            "| CPM, BPM | exact; BPM reduces to CPM with fixed bin width |\n"
            "| RPKM     | exact; needs '--siz_bin'                       |\n"
            "| None     | correct up to a constant log offset; not edgeR |\n"
            "| RPGC     | needs '--sf_A' and '--sf_B'; not edgeR         |\n"
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
            "edgeR 'prior.count' before scaling by the fragment-bin overlap "
            "count (default: %(default)s).\n"
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
            "Inferred from the track when omitted, and cross-checked against "
            "it when given, so a value the track contradicts is refused "
            "rather than silently rescaling the fragment-bin overlap count.\n"
            "\n"
            "Required only when no track is read, i.e., when both "
            "'--n_ovlp_A' and '--n_ovlp_B' are supplied and "
            "'--substrate RPKM' needs a width for its scale factor.\n"
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
            "all fragments, how many bins each fragment spans.\n"
            "\n"
            "Written by 'compute_signal --report_n_ovlp', as "
            "'<track>.n_ovlp.txt' when that flag is given without a path. "
            "Computed from '--fil_A' when omitted, which requires a track "
            "written with '--normalizeUsing None'; supplying it skips that "
            "read and changes nothing else.\n"
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
        "-sfA",
        "--sf_A",
        dest="sf_A",
        type=float,
        default=None,
        help=(
            "Scaling factor read from 'bamCoverage --verbose'. Applies to "
            "'--substrate RPGC' only, for which it is required.\n"
            "\n"
            "Generate with '--exactScaling'; otherwise, the factor is a "
            "sampled estimate and so is every pseudocount derived from it.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "--sf-A",
        dest="sf_A",
        type=float,
        help=argparse.SUPPRESS,
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
        "--sf-B",
        dest="sf_B",
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
        ),
    )
    parser.add_argument(
        "--prt-arg",
        dest="prt_arg",
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
    skp_pfx: tuple[str, ...],
    siz_bin: int | None = None,
) -> None:
    """
    Print the verbose argument banner in the reviewed semantic order.
    """

    with redirect_stdout(sys.stderr):
        print("#############################################")
        print("## Arguments for 'compute_pseudo_deeptools' ##")
        print("#############################################")
        print("")
        print("--verbose")
        print(f"--fil_A   {args.fil_A}")

        if getattr(args, "fil_B", None):
            print(f"--fil_B   {args.fil_B}")

        print(f"--skp_pfx {skp_pfx}")
        print(f"--substrate {args.substrate}")
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

        if args.sf_A is not None:
            print(f"--sf_A    {args.sf_A}")

        if args.sf_B is not None:
            print(f"--sf_B    {args.sf_B}")

        if args.prt_arg:
            print("--prt_arg")

        print(f"--dp      {args.dp}")

        if args.prt_jsn:
            print("--prt_jsn")

        print("")
        print("")


def _is_one_track(args: argparse.Namespace) -> bool:
    """
    Report whether the request describes a single track.

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
    'prior.count', and the denominator becomes 'L + 2 * prior.count'.

    Reproducing that here therefore needs no separate estimator: passing the
    one fragment-bin overlap count as both 'n_ovlp_a' and 'n_ovlp_b' makes
    'L_bar' equal 'L_A', which is what 'ave_lib' becomes when 'nlib' is 1.

    The consequence a caller must know: a track's single-track pseudocount is
    not its two-track pseudocount, because 'L_bar' is that track's own count in
    one mode and the mean of both in the other. That is edgeR's own behavior,
    not an artifact here.
    """

    return not getattr(args, "fil_B", None) and args.n_ovlp_B is None


def _run_edger(
    args: argparse.Namespace,
    skp_pfx: tuple[str, ...],
) -> int:
    """
    Emit scale factors and pseudocounts from edgeR's prior rule.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments carrying the deepTools substrate policy.
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
        For a nonpositive fragment-bin overlap count, an unusable substrate
        request, an unreadable track, a '--siz_bin' the track contradicts, a
        '--substrate RPKM' run with no width and no track to infer one from, or
        a '--substrate RPGC' run without '--sf_A'.

    Notes
    -----
    Fragment-bin overlap counts come from '--n_ovlp_A' and '--n_ovlp_B' when
    supplied, and are otherwise summed from the tracks. Summing requires a
    whole-count track, since a CPM or RPKM bedGraph sums to its own normalized
    total rather than to 'L', and would silently rescale every value this
    function returns.

    Reading a track also resolves the bin width, so '--siz_bin' is needed only
    when both fragment-bin overlap counts are supplied, which reads no track,
    and '--substrate RPKM' still wants a width for its scale factor.

    One result line reaches stdout, the first of these that applies:
      - a deepTools argument string under '--prt_arg',
      - a single value in single-track mode, or
      - the pseudocount pair as 'A:B'.

    A JSON summary follows on stdout under '--prt_jsn'.
    """

    one_track = _is_one_track(args)

    n_ovlp_a = args.n_ovlp_A
    n_ovlp_b = args.n_ovlp_B
    sf_a, sf_b = args.sf_A, args.sf_B

    siz_bin = args.siz_bin

    # Reading a track resolves the bin width, so a run that recovers either
    # overlap count needs no '--siz_bin' at all. Only a run given both of them
    # reads nothing, and only 'RPKM' then still needs a width.
    try:
        if n_ovlp_a is None:
            counts_a = sum_counts_bdg(args.fil_A, siz_bin, skp_pfx)
            n_ovlp_a, siz_bin = counts_a.total, counts_a.siz_bin

        if one_track:
            # Mirror A onto B so 'L_bar' collapses to 'L_A'. That reproduces
            # edgeR's 'nlib == 1' behavior exactly rather than approximating
            # it; see '_is_one_track' for the source reading it comes from.
            n_ovlp_b = n_ovlp_a
            sf_b = sf_a
        elif n_ovlp_b is None:
            counts_b = sum_counts_bdg(args.fil_B, siz_bin, skp_pfx)
            n_ovlp_b, siz_bin = counts_b.total, counts_b.siz_bin
    except (OSError, ValueError) as e:
        if args.verbose:
            _print_pseudo_arguments(args, skp_pfx, siz_bin)

        raise SystemExit(str(e)) from None

    if args.verbose:
        _print_pseudo_arguments(args, skp_pfx, siz_bin)

    if siz_bin is None and canonicalize_substrate(args.substrate) == "RPKM":
        raise SystemExit(
            "'--siz_bin' is required for '--substrate RPKM' when both "
            "'--n_ovlp_A' and '--n_ovlp_B' are supplied, because no track is "
            "read to infer the bin width from.",
        )

    try:
        result = compute_pseudo_edger(
            n_ovlp_a=n_ovlp_a,
            n_ovlp_b=n_ovlp_b,
            prior_count=args.prior_count,
            substrate=args.substrate,
            siz_bin=siz_bin,
            scale_a=sf_a,
            scale_b=sf_b,
        )
    except ValueError as e:
        raise SystemExit(str(e)) from None

    if not result["is_edger"]:
        print(
            f"Note: '--substrate {args.substrate}' does not reproduce edgeR's "
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
                "substrate": args.substrate,
                "prior_count": args.prior_count,
                "siz_bin": siz_bin,
                "dp": args.dp,
                "skp_pfx": list(skp_pfx),
            },
            "n_ovlp": {"A": n_ovlp_a, "B": n_ovlp_b},

            # The fractional substrates own 'k_A' and 'k_B', and this tool
            # serves none of them, so the key is absent rather than null.
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
    Every run takes edgeR's prior rule; there is no '--method', because the
    distribution-based estimators read the value distribution alone and are
    therefore substrate-agnostic. Run those through 'compute_pseudo', which
    accepts a deepTools track for them without needing to know its substrate.

    With '--prt_jsn', the command prints a strict one-line JSON summary only
    when all serialized values are finite; otherwise it warns on stderr, omits
    the JSON line, and returns zero.
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
        one_track = _is_one_track(args)

        validate_comparison(
            args.prior_count,
            "ge",
            0.0,
            "prior_count",
            allow_none=False,
        )
        validate_comparison(args.siz_bin, "gt", 0, "siz_bin", allow_none=True)
        validate_comparison(args.dp, "ge", 0, "dp", allow_none=False)

        if one_track:
            need_sf = args.sf_A is None
            both = "'--sf_A'"
        else:
            need_sf = args.sf_A is None or args.sf_B is None
            both = "both '--sf_A' and '--sf_B'"

        if args.substrate == "RPGC" and need_sf:
            raise ValueError(
                f"'--substrate RPGC' requires {both}; read them from "
                "'bamCoverage --verbose'.",
            )

        if one_track and args.prt_arg:
            # The deepTools 'bamCompare' options '--scaleFactors' and
            # '--pseudocount' take a pair each, so there is no single-track
            # spelling of them. 'bamCoverage --scaleFactor' exists but accepts
            # no pseudocount, so emitting it would silently drop the value that
            # was asked for.
            raise ValueError(
                "'--prt_arg' writes the two-track 'bamCompare' argument "
                "string, which has no single-track form. Drop '--prt_arg' to "
                "print the pseudocount, or supply track B.",
            )
    except ValueError as e:
        raise SystemExit(str(e)) from None

    skp_pfx = parse_skp_pfx(args.skp_pfx, default=DEF_SKP_PFX)

    _run_edger(args, skp_pfx)

    # Reaching this line is the one success path, because '_run_edger' returns
    # zero or leaves through 'SystemExit'.
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
