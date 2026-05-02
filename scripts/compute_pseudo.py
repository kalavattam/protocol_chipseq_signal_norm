#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.

"""
Script
------
compute_pseudo.py


Description
-----------
Computes pseudocounts for bedGraph signal tracks using a few simple, defensible
rules (see details below). Intended to be run on tracks that are already
scaled/normalized for ratio computation, e.g., spike-in-scaled signal and siQ-
ChIP-scaled signal (PMID 40364978), as well as scaling methods (SES, readCount)
and normalization methods ({F,R}PKM, CPM, {T,B}PM, {F,R}PGC) implemented in
deepTools.


Usage
-----
python -m scripts.compute_pseudo \
    [--help] [--verbose] \
    --fil_A <str_file.(bedgraph|bdg|bg)(.gz)> \
    [--fil_B <str_file.(bedgraph|bdg|bg)(.gz)>] \
    [--method {frc_mdn_nz,qntl_nz,frc_avg_nz,min_nz}] [--coef <float>] \
    [--qntl_nz <0..100>] [--floor <float>] [--eps <float>] \
    [--mode_nz {closed,open,off}] \
    [--sym {none,max,min,arith,geom,harm,use_A,use_B}] [--rnd <int_decimals>] \
    [--skp_pfx <csv>] [--prt_jsn]


Output
------
- Prints a one-line recommendation to stdout: '<pseudo_A>' or
  '<pseudo_A>:<pseudo_B>'
- With "--prt_jsn", also prints a JSON summary to stdout.
- Use "--rnd" to control rounding displayed in the recommendation.
- Supports plain text, ".gz", and "-" (stdin) via 'utils_io.open_in()'.


Examples
--------
1. Default (median-based, coef=0.01), two tracks
'''bash
python -m scripts.compute_pseudo \
    --fil_A ip.cpm.bdg.gz \
    --fil_B in.cpm.bdg.gz
'''

2. Quantile method: 1st percentile; make symmetric via geometric mean
'''bash
python compute_pseudo.py \
    --fil_A ip.rpkm.bdg \
    --fil_B in.rpkm.bdg \
    --method qntl_nz \
    --qntl_nz 1 \
    --sym geom
'''

3. Minimum-nonzero method with larger coef and epsilon for float noise
'''bash
python -m scripts.compute_pseudo \
    --fil_A ip.cpm.bdg.gz \
    --method min_nz \
    --coef 0.5 \
    --eps 1e-12
'''

4. One track only; mirror via '--sym max' (prints 'pseudo_A:pseudo_B')
'''bash
python -m scripts.compute_pseudo \
    --fil_A ip.cpm.bdg.gz \
    --sym max
'''

5. Also emit a JSON summary (last line is JSON)
'''bash
python -m scripts.compute_pseudo \
    --fil_A ip.cpm.bdg.gz \
    --prt_jsn
'''


Methods (per track)
-------------------
- "frc_mdn_nz"  c = coef * median of nonzero bins     [default; typ. coef=0.01]
- "qntl_nz"     c = q-th percentile of nonzero bins   [q in 0..100]‡
- "frc_avg_nz"  c = coef * mean of nonzero bins       [typ. coef=0.01]
- "min_nz"      c = coef * minimum nonzero bin value  [typ. coef=0.1–1.0]

‡[decimals OK; nearest rank via 'round']


Combination ('--sym') across A and B
------------------------------------
- "none"   keep asymmetric pseudo_A:pseudo_B                  [default]
- "max"    use max(pseudo_A, pseudo_B) for both tracks
- "min"    use min(pseudo_A, pseudo_B) for both tracks
- "arith"  use (pseudo_A + pseudo_B) / 2 for both tracks
- "geom"   use sqrt(pseudo_A * pseudo_B) for both tracks
- "harm"   use 2 / (1/pseudo_A + 1/pseudo_B) for both tracks
- "use_A"  copy pseudo_A to both tracks
- "use_B"  copy pseudo_B to both tracks

'--sym geom' assumes nonnegative pseudocounts; if either is 0, the result will
be 0. '--sym harm' requires both pseudocounts are greater than 0; if not, the
processing code falls back to the '--sym min' computation with a warning.


Zero handling
-------------
- Parsing uses 'utils_bdg.iter_bdg_rows', which skips headers, blank lines,
  and malformed rows, and yields 'v_num=None' for non-finite tokens; those are
  ignored when computing statistics.
- By default, statistics are computed on “nonzero” bins only, where nonzero
  means '|x| > eps' ('eps=0.0' by default). Set '--eps' to, e.g., treat tiny
  float noise as zero.
- Use '--mode_nz' to control how the epsilon filter is applied:
    + '--mode_nz closed' (default): drop values with '|x| <= eps'.
    + '--mode_nz open': drop values with '|x| < eps' (so, with 'eps=0', exact
      zeros are included).
    + '--mode_nz off': disable epsilon-based filtering entirely; all finite
      values are used.
- Practical notes on zero handling:
    + With very small 'eps' and quantile/min methods, low pseudocounts are
      common; consider a positive '--floor' or larger '--coef' if needed.


General notes
-------------
- If unsure, start with '--method frc_mdn_nz' (default), leave '--coef'
  unset (auto 0.01), and retain the default '--mode_nz closed'. Also, consider
  using a small '--eps' to treat float noise as zero; the value to use depends
  on the input data value range.
- A note on '--method min_nz': The minimum nonzero bin is often at the noise
  floor.
    + To keep the pseudocount useful (but not distortive), the coefficient
      usually needs to be an order of magnitude or two larger than for the
      median and mean methods.
    + As a rule of thumb, start with 'coef = 0.1–1.0' for '--method min_nz' vs
      'coef = 0.01' for '--method frc_mdn_nz' and '--method frc_avg_nz'.


Performance notes
-----------------
- Materializes all signal values per track to compute medians, means, or
  quantiles; for very large bedGraph files, this can be memory- and time-
  intensive.
- Ideas for improving performance:
    + Add a '--sample N' option to subsample values uniformly for faster,
      approximate stats (i.e., similar to the default behavior of deepTools).
    + Implement streaming quantile estimators (e.g., P^2 or t-digest) to avoid
      full sorts; median/quantiles would then be approximate but scalable.
    + Keep current behavior as the default, but enable “approximate modes”
      behind explicit flags.


#TODO
-----
- Regarding the previous 'module_tag' / 'warn_pseudo()' approach to generating
  error text (which has since been refactored out): Flesh this out, then export
  it to other Python scripts?
- Extend '.gz' handling to '.bgz' and '.bgzf'. (Here and in other scripts.)
- JSON summarization? (I have roughly implemented it here.)
"""

from __future__ import annotations

import argparse
import json
import math
import signal
import sys

from contextlib import redirect_stdout

from scripts.functions.utils_check import (
    check_cmp,
    check_exists,
)
from scripts.functions.utils_cli import (
    add_help_cap,
    CapArgumentParser,
)
from scripts.functions.utils_interactive import (
    echo_block,
    get_args_interactive,
)
from scripts.functions.utils_io import (
    DEF_SKP_PFX,
    ensure_single_stdin,
    parse_skp_pfx,
)
from scripts.functions.utils_stabilizer import (
    compute_stats_robust,
    determine_coef_eff,
    iter_vals_bdg,
    pick_stabilizer,
)

try:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except (AttributeError, ValueError):
    pass

assert sys.version_info >= (3, 10), "Python >= 3.10 required."

#  Run script in “interactive mode” (true) or “CLI mode” (false)
interactive = False


#  Define functions
def set_interactive(
    echo: bool = False, ensure: bool = False
) -> argparse.Namespace:
    """
    Set paths and parameters for interactive mode.
    """
    #  Set general paths
    dir_bas = "/home/kalavatt/tsukiyamalab/Kris"
    dir_rep = f"{dir_bas}/protocol_chipseq_signal_norm"
    dir_dat = f"{dir_rep}/data"
    dir_pro = f"{dir_dat}/processed"
    dir_exp = f"{dir_pro}/2025-09-01"
    dir_sig = f"{dir_exp}/compute_signal"

    aln = "bowtie2"
    atp = "global"
    flg = 2
    mpq = 1
    det = f"{aln}_{atp}_flag-{flg}_mapq-{mpq}"
    sig = "norm"  # "alpha", "norm", "raw", or "spike"
    spc = "sc"
    ext = "bdg.gz"

    dir_nrm = f"{dir_sig}/{det}/{sig}"
    fil_A = f"IP_SIR2-Mcm4-flag_H3_G1-1st_7562.{spc}.{ext}"
    fil_B = f"in_SIR2-Mcm4-flag_H3_G1-1st_7562.{spc}.{ext}"
    pth_A = f"{dir_nrm}/{fil_A}"
    pth_B = f"{dir_nrm}/{fil_B}"

    #  Check that paths exist (optional for debugging/development purposes)
    if ensure:
        check_exists(dir_nrm, kind="dir")
        check_exists(pth_A, kind="file")
        if pth_B is not None:
            check_exists(pth_B, kind="file")

    #  Set argument values wrapped in argparse.Namespace
    ns = argparse.Namespace(
        verbose=True,
        fil_A=pth_A,
        fil_B=pth_B,          # Optional second track, e.g., for per-track
                              # stats and/or true cross-track symmetrization
                              # via '--sym'
        method="frc_mdn_nz",  # Default method: 'frc_mdn_nz'
        qntl_nz=1.0,          # Only used if '--method qntl_nz'
        coef=None,            # Let script pick 0.01 or 1.0 (method-dep.)
        floor=0.0,            # Lower bound on pseudo; set >0 to avoid 0
        eps=0.0,              # Consider a tiny ε to tolerate float noise
        mode_nz="closed",     # Epsilon/zero-handling mode
        sym="none",           # Symmetrization; 'none' keeps A/B separate
        rnd=24,               # Decimal precision for computed pseudo
        skp_pfx=",".join(DEF_SKP_PFX),  # bedGraph header/comments to skip
        prt_jsn=False         # Set True to also emit a JSON summary line
    )

    #  “Echo” the interactive argument assignments (etc.) if specified
    if echo:
        echo_block("paths, etc.", {
            "dir_bas": dir_bas, "dir_rep": dir_rep, "dir_dat": dir_dat,
            "dir_pro": dir_pro, "dir_exp": dir_exp, "dir_sig": dir_sig,
            "dir_nrm": dir_nrm,
            "fil_A":   fil_A,   "fil_B":   fil_B,
            "pth_A":   pth_A,   "pth_B":   pth_B,
        })
        echo_block("args", vars(ns))

    #  Return the argparse.Namespace-wrapped arguments
    return ns


def combine_pseudo_sym(
    pseudo_A: float,
    pseudo_B: float,
    mode: str = "none"
) -> tuple[float, float]:
    """
    Combine per-track pseudocounts into a symmetric pair if requested.

    Args:
        pseudo_A : float
            Pseudocount proposed for track A (may be finite or nonfinite).
        pseudo_B : float
            Pseudocount proposed for track B (may be finite or nonfinite).
        mode : {"none","max","min","arith","geom","harm","use_A","use_B"}
            Symmetrization rule; see notes below.

    Returns:
        tuple[float, float]
            A pair (pseudo_A, pseudo_B) after applying the symmetrization
            policy.

    Raises:
        ValueError
            If 'mode' is not one of the recognized options.

    Notes:
        - If mode == "none", then return as-is, even if a nonfinite value is
          present.
        - Symmetric modes:
            + If both are finite, then combine as requested.
            + If exactly one is finite, then mirror the finite value while
              issuing a warning.
            + If both are nonfinite, then return as-is (i.e., let the user
              decide how to proceed) while issuing a warning.
    """
    if mode == "none":
        return pseudo_A, pseudo_B

    fnt_A = math.isfinite(pseudo_A)
    fnt_B = math.isfinite(pseudo_B)

    if fnt_A and fnt_B:
        if mode == "max":
            pseudo = max(pseudo_A, pseudo_B)
        elif mode == "min":
            pseudo = min(pseudo_A, pseudo_B)
        elif mode == "arith":
            pseudo = 0.5 * (pseudo_A + pseudo_B)
        elif mode == "geom":
            if pseudo_A < 0.0 or pseudo_B < 0.0:
                print(
                    "Geometric mean undefined for negative values; falling "
                    "back to min(pseudo_A, pseudo_B).",
                    file=sys.stderr
                )
                pseudo = min(pseudo_A, pseudo_B)
            else:
                pseudo = math.sqrt(pseudo_A * pseudo_B)
        elif mode == "harm":
            if pseudo_A <= 0.0 or pseudo_B <= 0.0:
                print(
                    "Harmonic mean undefined for nonpositive values; falling "
                    "back to min(pseudo_A, pseudo_B).",
                    file=sys.stderr
                )
                pseudo = min(pseudo_A, pseudo_B)
            else:
                pseudo = 2.0 / (1.0 / pseudo_A + 1.0 / pseudo_B)
        elif mode == "use_A":
            pseudo = pseudo_A
        elif mode == "use_B":
            pseudo = pseudo_B
        else:
            raise ValueError(f"Error: Unknown --sym: {mode!r}")
        return pseudo, pseudo

    if fnt_A and not fnt_B:
        print(
            "pseudo_B is nonfinite; mirroring pseudo_A in symmetric mode "
            f"{mode!r}.",
            file=sys.stderr
        )
        return pseudo_A, pseudo_A

    if fnt_B and not fnt_A:
        print(
            "pseudo_A is nonfinite; mirroring pseudo_B in symmetric mode "
            f"{mode!r}.",
            file=sys.stderr
        )
        return pseudo_B, pseudo_B

    print("Both pseudocounts are nonfinite; returning as-is.", file=sys.stderr)
    return pseudo_A, pseudo_B


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = CapArgumentParser(
        description=(
            "Compute data-driven pseudocount(s) for one (or two) bedGraph "
            "track(s)."
        )
    )
    add_help_cap(parser)
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        default=False,
        help="Increase output verbosity.\n\n"
    )
    parser.add_argument(
        "-fA", "--fil_A",
        required=True,
        help=(
            "Path to the first bedGraph file (file 'A'), or '-' for stdin; "
            ".gz is handled.\n\n"
        )
    )
    parser.add_argument(
        "-fB", "--fil_B",
        help=(
            "Path to an optional second bedGraph file (file 'B'), or '-' for "
            "stdin; .gz is handled.\n\n"
        )
    )
    parser.add_argument(
        "-m", "--method",
        choices=("frc_mdn_nz", "qntl_nz", "frc_avg_nz", "min_nz"),
        default="frc_mdn_nz",
        help=(
            "Method to compute per-track pseudocount (default: %(default)s):\n"
            "    - frc_mdn_nz  value = coef × median of nonzero bins\n"
            "    - qntl_nz     value = q-th percentile of nonzero bins\n"
            "    - frc_avg_nz  value = coef × mean of nonzero bins\n"
            "    - min_nz      value = coef × minimum nonzero bin\n"
            "\n"
            "Notes:\n"
            "    - If '--method qntl_nz', set percentile with '--qntl_nz' "
            "[decimals OK (e.g., 0.1 = 0.1th percentile); nearest-rank "
            "determined via 'round'].\n"
            "    - 'nonzero' means '|x| > eps' (via '--eps' as defined by "
            "'--mode_nz').\n"
            "    - If '--coef' is omitted, then defaults to '--coef 0.01' for "
            "'--method frc_*' and '--coef 1.0' for '--method min_nz'.\n"
            "    - '--method min_nz' typically needs a larger coef "
            "(e.g., 0.1–1.0) in comparison to '--method frc_*' (e.g., "
            "0.01).\n\n"
        )
    )
    parser.add_argument(
        "-q", "--qntl_nz",
        type=float,
        default=1.0,
        help=(
            "Quantile in percent for '--method qntl_nz' (0..100). Decimals "
            "are allowed (e.g., 0.5 = 0.5th percentile). Ignored if "
            "'--method' is not 'qntl_nz' (default: %(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-c", "--coef",
        type=float,
        default=None,
        help=(
            "Coefficient for median, mean, and min methods. If not specified, "
            "then defaults to 0.01 for '--method frc_mdn_nz' and '--method "
            "frc_avg_nz', or 1.0 for '--method min_nz'.\n\n"
        )
    )
    parser.add_argument(
        "-fl", "--floor",
        type=float,
        default=0.0,
        help=(
            "Lower bound for computation of pseudocount(s) (default: "
            "%(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-e", "--eps",
        type=float,
        default=0.0,
        help=(
            "Zero tolerance ε (epsilon) in computation of pseudocount(s) "
            "(default: %(default)s).\n"
            "\n"
            "When '--mode_nz closed' (default), values with '|x| <= ε' are "
            "treated as zero and excluded from statistics; with "
            "'--mode_nz open', values with '|x| < ε' are excluded; with "
            "'--mode_nz off', ε-based filtering is disabled.\n\n"
        )
    )
    parser.add_argument(
        "-mz", "--mode_nz",
        choices=("closed", "open", "off"),
        default="closed",
        help=(
            "Epsilon/zero-handling mode for selecting 'nonzero' bins:\n"
            "    - closed (default): drop values with '|x| <= eps'.\n"
            "    - open            : drop values with '|x| < eps'.\n"
            "    - off             : disable 'eps'-based zero filtering.\n\n"
        )
    )
    parser.add_argument(
        "-s", "--sym",
        choices=(
            "none", "max", "min", "arith", "geom", "harm", "use_A", "use_B"
        ),
        default="none",
        help=(
            "Symmetrize pseudocounts across A and B.\n"
            "    - If both A AND B are given, apply the chosen rule to "
            "'(pseudo_A, pseudo_B)', then print 'pseudo_A:pseudo_B'.\n"
            "    - If only A is given AND '--sym' is provided AND NOT 'none', "
            "mirror A to B and print: 'pseudo_A:pseudo_B'.\n"
            "    - If only A is given AND '--sym' is omitted OR '--sym none', "
            "print a single value: 'pseudo_A'.\n\n"
        )
    )
    parser.add_argument(
        "-r", "--rnd",
        type=int,
        default=24,
        help=(
            "Number of decimal places for rounding output pseudocount(s) "
            "(default: %(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-sp", "--skp_pfx",
        type=str,
        default=",".join(DEF_SKP_PFX),
        help=(
            "Comma-separated prefixes to skip in bedGraph file(s) (default: "
            "%(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-pj", "--prt_jsn",
        action="store_true",
        default=False,
        help="Print a JSON summary to stdout.\n\n"
    )

    #  If no arguments are provided, use 'argv' to display help and exit
    argv_parse = sys.argv[1:] if argv is None else argv
    if not argv_parse:
        parser.print_help(sys.stderr)
        sys.exit(0)

    return parser.parse_args(argv_parse)


def main(argv: list[str] | None = None) -> None:
    """
    Execute the primary control flow for the script.

    Args:
        None.

    Returns:
        None. On success, prints one pseudocount (or 'A:B') to stdout and, if
        requested, a single-line JSON summary.

    Exits:
        - 0 on success or when showing help with no arguments.
        - 1 on validation or computation errors, including (non-exhaustive):
            + Non-finite or negative '--coef', '--floor', or '--eps'
            + Invalid '--method' or out-of-range/NaN '--qntl_nz'
            + File I/O errors or malformed bedGraph lines
            + Strict JSON serialization failure when '--prt_jsn' is set

    Behavior:
        - Parses arguments (or uses interactive defaults if enabled).
        - Validates numeric options ('--coef', '--floor', '--eps';
          non-negative).
        - Reads bedGraph(s), extracting column 4 values with optional header
          skipping and zero filtering via '--eps' / '--mode_nz'.
        - Computes per-track pseudocount(s) according to '--method' and
          options.
        - Optionally symmetrizes A/B via '--sym'.
        - Prints either a single value or 'A:B' to stdout with '--rnd'
          rounding.
        - With '--prt_jsn', also prints a one-line JSON summary (strict JSON;
          NaN/inf are rejected).

    Notes:
        - Prints human-readable error messages to stderr on failure.
        - Warnings to stderr for empty/non-finite inputs, zero pseudocounts,
          and incompatible symmetrization cases.
    """
    #  Handle interactive or CLI arguments
    args, early_exit = get_args_interactive(
        argv, interactive, set_interactive, parse_args
    )
    if early_exit:
        return

    #  Check input file(s) and stdin usage
    paths = [
        p for p in (args.fil_A, getattr(args, "fil_B", None)) if p is not None
    ]

    #  Allow at most one "-" across 'fil_A' or 'fil_B'
    try:
        ensure_single_stdin(paths)
    except ValueError as e:
        #  Convert to a clean CLI error without a traceback
        raise SystemExit(str(e))

    #  For non-stdin paths, ensure the files exist
    try:
        for label, p in (
            ("A", args.fil_A),
            ("B", getattr(args, "fil_B", None))
        ):
            if p is None or p == "-":
                continue
            check_exists(p, kind="file", label=f"bedGraph {label}")
    except FileNotFoundError as e:
        #  Convert to a clean CLI error without a traceback
        raise SystemExit(str(e))

    #  Check numeric arguments
    try:
        if args.method == "qntl_nz" and (
            not math.isfinite(args.qntl_nz) or
            not (0.0 <= args.qntl_nz <= 100.0)
        ):
            raise ValueError("'--qntl_nz' must be finite and in [0, 100].")

        check_cmp(args.coef,  "ge", 0.0, "coef",  allow_none=True)
        check_cmp(args.floor, "ge", 0.0, "floor", allow_none=False)
        check_cmp(args.eps,   "ge", 0.0, "eps",   allow_none=False)
        check_cmp(args.rnd,   "ge", 0,   "rnd",   allow_none=False)
    except ValueError as e:
        raise SystemExit(str(e))

    #  Resolve effective coefficient once for the chosen method
    coef_eff = determine_coef_eff(args.method, args.coef)

    #  Parse bedGraph line prefixes to skip
    skp_pfx = parse_skp_pfx(args.skp_pfx, default=DEF_SKP_PFX)

    #  Handle inclusion of zero/epsilon in statistics
    mode_nz = args.mode_nz

    #  Optionally emit a verbose argument banner
    if args.verbose:
        with redirect_stdout(sys.stderr):
            print("####################################")
            print("## Arguments for 'compute_pseudo' ##")
            print("####################################")
            print("")
            print("--verbose")
            print(f"--fil_A   {args.fil_A}")
            if getattr(args, "fil_B", None):
                print(f"--fil_B   {args.fil_B}")
            print(f"--method  {args.method}")
            if args.method == "qntl_nz":
                print(f"--qntl_nz {args.qntl_nz}")
            if coef_eff is not None and coef_eff != args.coef:
                print(f"--coef    {args.coef}  ## coef_eff = {coef_eff} ##")
            else:
                print(f"--coef    {args.coef}")
            print(f"--floor   {args.floor}")
            print(f"--eps     {args.eps}")
            print(f"--mode_nz {mode_nz}")
            print(f"--sym     {args.sym}")
            print(f"--rnd     {args.rnd}")
            print(f"--skp_pfx {skp_pfx}")
            if args.prt_jsn:
                print("--prt_jsn")
            print("")
            print("")

    vals_A = list(
        iter_vals_bdg(
            args.fil_A,
            eps=args.eps,
            mode_nz=mode_nz,
            skp_pfx=skp_pfx
        )
    )
    if not vals_A:
        print(
            "No finite values in A after filtering; pseudocount will be "
            "'nan'. Check '--eps', '--mode_nz', and/or '--skp_pfx'.",
            file=sys.stderr
        )
    pseudo_A = pick_stabilizer(
        vals_A,
        method=args.method,
        coef=coef_eff,
        qntl_pct=args.qntl_nz,
        floor=args.floor,
        qntl_rule="round",
    )

    pseudo_B = float("nan")
    vals_B = []
    if getattr(args, "fil_B", None):
        vals_B = list(
            iter_vals_bdg(
                args.fil_B,
                eps=args.eps,
                mode_nz=mode_nz,
                skp_pfx=skp_pfx
            )
        )
        if not vals_B:
            print(
                "No finite values in B after filtering; pseudocount will be "
                "'nan'. Check '--eps', '--mode_nz', and/or '--skp_pfx'.",
                file=sys.stderr
            )
        pseudo_B = pick_stabilizer(
            vals_B,
            method=args.method,
            coef=coef_eff,
            qntl_pct=args.qntl_nz,
            floor=args.floor,
            qntl_rule="round",
        )
    else:
        #  If only one file is supplied, then mirror A to B if and only if the
        #  user explicitly requested a non-'none' '--sym' rule
        if (args.sym != "none") and math.isfinite(pseudo_A):
            pseudo_B = pseudo_A

        #  NOTE: using most '--sym' modes with only 'fil_A' is mathematically
        #        degenerate ('pseudo_A:pseudo_A'); allow it, but inform the
        #        user
        fil_B = getattr(args, "fil_B", None)  # MAYBE: OK to to delete this?
        if fil_B is None and args.sym not in {"none", "use_A"}:
            print(
                f"Note: '--sym {args.sym}' was used with only '--fil_A'; "
                "symmetrization degenerates to 'pseudo_A:pseudo_A'.",
                file=sys.stderr,
            )

    #  Warn if a computed pseudo is exactly 0 (can cause -inf in log2 ratios)
    for tag, pseudo in (("A", pseudo_A), ("B", pseudo_B)):
        if math.isfinite(pseudo) and pseudo == 0.0:
            print(
                f"Pseudocount for {tag} is 0.0; log2 ratios may produce -inf "
                "at zeros. Consider a positive '--floor' or a larger "
                "'--coef'.",
                file=sys.stderr
            )

    #  Apply symmetrization only if the user requested a non-'none' mode for
    #  '--sym'
    if args.sym != "none":
        pseudo_A, pseudo_B = combine_pseudo_sym(
            pseudo_A, pseudo_B, mode=args.sym
        )

    #  Re-check after symmetrization (as '--sym' modes can introduce zeros)
    for tag, pseudo in (("A", pseudo_A), ("B", pseudo_B)):
        if math.isfinite(pseudo) and pseudo == 0.0:
            print(
                f"Pseudocount for {tag} is 0.0 after symmetrization; log2 "
                "ratios may produce -inf at zeros. Consider a positive "
                "'--floor' or larger '--coef'.",
                file=sys.stderr
            )

    #  Prepare printed recommendation
    want_pair = (
        getattr(args, "fil_B", None) is not None
    ) or (
        args.sym != "none"
    )

    #  Set up formatter (helper function) for CLI output: print 'nan' for non-
    #  finite, otherwise print values with a specified decimal precision
    def fmt(x: float) -> str:
        return "nan" if not math.isfinite(x) else f"{x:.{args.rnd}f}"

    #  Emit the primary recommendation line; if a pair was requested (two files
    #  or explicit '--sym'), print A:B; otherwise (single file, '--sym none'),
    #  print a single value
    if want_pair:
        out_A = fmt(pseudo_A)
        out_B = fmt(pseudo_B if math.isfinite(pseudo_B) else pseudo_A)
        print(f"{out_A}:{out_B}")
    else:
        print(f"{fmt(pseudo_A)}")

    #  Optionally emit a single-line JSON summary for potential logging/parsing
    #
    #  NOTE: workshopped JSON summarization in this script; may add to other
    #        other repo scripts in the future
    if args.prt_jsn:
        #  Include the computed pseudocounts; for a requested pair, include B
        pseudocounts = {"pseudo_A": pseudo_A, "pseudo_A_str": fmt(pseudo_A)}
        if want_pair:
            p_B = pseudo_B if math.isfinite(pseudo_B) else pseudo_A
            pseudocounts.update({"pseudo_B": p_B, "pseudo_B_str": fmt(p_B)})

        #  Include the provided parameters; for a requested pair, include B
        params = {"coef": args.coef}       # As provided by user (may be None)

        if coef_eff is not None:           # Ensure eff is directly after coef
            params["coef_eff"] = coef_eff  # Actual coef used by method

        params.update({
            "qntl_nz": args.qntl_nz,
            "floor": args.floor,
            "eps": args.eps,
            "mode_nz": mode_nz,
            "sym": args.sym,
            "rnd": args.rnd,
            "skp_pfx": list(skp_pfx)
        })

        #  Include output details
        out = {
            "fil_A": args.fil_A,
            "fil_B": getattr(args, "fil_B", None),
            "method": args.method,
            "params": params,
            "stats": {
                "A": compute_stats_robust(vals_A),
                "B": compute_stats_robust(vals_B) if vals_B else None
            },
            "pseudocounts": pseudocounts
        }

        #  NOTE: strict parsers reject nan, inf, etc.; will need to update
        #        and/or extend the below code if I ever plan to parse JSON
        #        output
        try:
            #  Print as a single line to keep logs tidy
            print(json.dumps(out, separators=(",", ":"), allow_nan=False))
        except ValueError:
            print(
                "Strict JSON disallows nan and inf; adjust '--floor' and "
                "'--coef', or just skip '--prt_jsn'.",
                file=sys.stderr
            )


if __name__ == "__main__":
    try:
        # FIXME: 'main()' sometimes effectively returns an int value via
        #        'sys.exit(main() or 0)'; change (above)
        #        'def main(argv: list[str] | None = None) -> None:' to
        #        'def main(argv: list[str] | None = None) -> int | None:', or
        #        change (here)
        #        '''
        #        sys.exit(main() or 0)
        #        ''' to
        #        '''
        #        main()
        #        sys.exit(0)
        #        '''
        #
        # NOTE: doing the above makes “run directly” differ in style from other
        #       Python scripts in repo; if make changes here, then need to
        #       consider making changes elsewhere
        sys.exit(main() or 0)
    except BrokenPipeError:
        try:
            sys.stdout.close()
        except Exception:
            pass
        try:
            sys.stderr.close()
        except Exception:
            pass
        sys.exit(0)
