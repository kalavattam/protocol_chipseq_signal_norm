#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: compute_pseudo.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in development and documentation.
#
# Distributed under the MIT license.


"""
Compute pseudocount recommendations for bedGraph signal tracks.

The CLI accepts input tracks, statistic method, coefficient, nonzero filtering,
symmetry mode, and formatting options. It prints one pseudocount or an A:B
pseudocount pair to stdout and, with '--prt_jsn', also prints a JSON summary.

Examples
--------
python -m protocol_chipseq_signal_norm.cli.compute_pseudo [options] \\
    --fil_A <file> [--fil_B <file>]
"""

from __future__ import annotations

import argparse
import json
import math
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
from protocol_chipseq_signal_norm.utilities.utils_io import (
    DEF_SKP_PFX,
    ensure_single_stdin,
    parse_skp_pfx,
)
from protocol_chipseq_signal_norm.utilities.utils_stabilizer import (
    compute_stats_robust,
    determine_coef_eff,
    iter_vals_bdg,
    pick_stabilizer,
)

with suppress(AttributeError, ValueError):
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

assert sys.version_info >= (3, 11), "Python >= 3.11 required."


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
        A pair (pseudo_A, pseudo_B) after applying the symmetrization
        policy.

    Raises
    ------
    ValueError
        If both pseudocounts are finite and 'mode' is not a recognized
        symmetrization rule.

    Notes
    -----
    - 'none' returns both inputs as-is, including nonfinite values. With two
      finite inputs, 'max' and 'min' select the respective input; 'arith',
      'geom', and 'harm' use the arithmetic, geometric, and harmonic means;
      and 'use_A' and 'use_B' copy the selected input to both outputs.
    - For finite inputs, 'geom' falls back to the minimum with a stderr
      warning if either input is negative; 'harm' does the same if either
      input is nonpositive.
    - For every mode except 'none':
        + If exactly one is finite, then mirror the finite value while
          issuing a warning.
        + If both are nonfinite, then return as-is (i.e., let the user
          decide how to proceed) while issuing a warning. These nonfinite
          paths do not validate 'mode'.

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
                        "Compute first-percentile pseudocounts and "
                        "symmetrize them by maximum."
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
        choices=("frc_mdn_nz", "qntl_nz", "frc_avg_nz", "min_nz"),
        default="frc_mdn_nz",
        help=(
            "Workflow method. Method to compute per-track pseudocount "
            "(default: %(default)s):\n"
            "    - frc_mdn_nz  value = coef × median of nonzero bins\n"
            "    - qntl_nz     value = q-th percentile of nonzero bins\n"
            "    - frc_avg_nz  value = coef × mean of nonzero bins\n"
            "    - min_nz      value = coef × minimum nonzero bin\n"
            "\n"
            "Notes:\n"
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

    argv_parse = sys.argv[1:] if argv is None else argv

    if not argv_parse:
        parser.print_help(sys.stderr)
        raise SystemExit(0)

    return parser.parse_args(argv_parse)


def _print_pseudo_arguments(
    args: argparse.Namespace,
    coef_eff: float | None,
    skp_pfx: tuple[str, ...],
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
        print(f"--dp     {args.dp}")

        if args.prt_jsn:
            print("--prt_jsn")

        print("")
        print("")


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
    Inputs are filtered through '--eps' and '--mode_nz' before pseudocount
    computation; malformed or nonnumeric bedGraph rows are skipped by the row
    iterator. Pairwise values can be symmetrized with '--sym'. With
    '--prt_jsn', the command prints a strict one-line JSON summary only when
    all serialized values are finite; otherwise it warns on stderr, omits the
    JSON line, and returns zero. Warnings about empty inputs, zero
    pseudocounts, or symmetrization are written to stderr.
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

        validate_comparison(args.coef, "ge", 0.0, "coef", allow_none=True)
        validate_comparison(args.floor, "ge", 0.0, "floor", allow_none=False)
        validate_comparison(args.eps, "ge", 0.0, "eps", allow_none=False)
        validate_comparison(args.dp, "ge", 0, "dp", allow_none=False)
    except ValueError as e:
        raise SystemExit(str(e)) from None

    coef_eff = determine_coef_eff(args.method, args.coef)

    skp_pfx = parse_skp_pfx(args.skp_pfx, default=DEF_SKP_PFX)

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
        return "nan" if not math.isfinite(value) else f"{value:.{args.dp}f}"

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
