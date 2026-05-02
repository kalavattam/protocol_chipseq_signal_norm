#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2025 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Distributed under the MIT license.

"""
Script
------
merge_bins_bdg.py


Description
------------
Merge adjacent same-valued bedGraph bins on each chromosome when intervals are
contiguous (i.e., when the previous end equals the current start) and values
“match.”

Matching modes (priority):
    1. '--rnd N' and '--eps E' together
       Round each finite value to N decimals, then merge when |Δ| <= E on
       those rounded numbers (printed with N decimals).

    2. '--rnd N' only
       Compare on N-decimal rounded strings (printed with N decimals).

    3. '--eps E' only
       Compare raw numeric values within |Δ| <= E (print preserves the first
       token).

    4. <default>
       Exact text match of the 4th field (after “canonicalizing” NaN/±inf).

Non-finite tokens ("nan", "inf", "-inf") never merge with anything else.
Header/metadata lines (e.g., '#', 'track', 'browser') are forwarded verbatim
and always terminate the current run (they don’t participate in merging). Data
are streamed, and .gz in/out and '-' for stdin/stdout are supported.


Usage
-----
merge_bins_bdg.py \
    -i <in.bdg[.gz]|-> \
    -o <out.bdg[.gz]|-> \
    [--rnd N] \
    [--eps E] \
    [--skp_pfx "#,track,browser"]

Typical patterns:
- Exact text match (default): merges only when column 4 tokens are identical.
- Rounded match: '--rnd N' merges values that round to the same N-decimal
  string.
- Numeric-tolerance match: '--eps E' merges values within |Δ| <= E.
  If '--rnd N' is also given, values are first rounded to N decimals and then
  compared with |Δ| <= E on the rounded numbers.

Notes:
- Headers recognized by 'utils_io.is_header' (using prefixes from '--skp_pfx')
  and blank/whitespace-only lines are forwarded as-is and always terminate the
  current run.
- Input is expected to be sorted by (chrom, start). Unsorted input may result
  in unexpected merges and the prevention of expected merges.


Output
------
Writes a merged bedGraph to '-o' with the following lines:
    chrom <tab> start <tab> end <tab> value

The value field is printed as follows:
- With '--rnd N', the first finite numeric token in a run is printed to N
  decimals.
- Without '--rnd', the first token of the run is preserved verbatim.
- Non-finite tokens never merge; they are emitted one per input interval.

Exit behavior:
- Returns 0 on success.
- Raises ValueError on malformed data lines (e.g., < 4 fields, non-integer
  coords, non-positive width).
- Prints an error and returns 1 when '--eps' is negative; rounding is validated
  to be '--rnd >= 0'.


Examples
--------
1. Merge by exact text (default behavior)
```bash
python merge_bins_bdg.py -i track.raw.bdg.gz -o track.merged.bdg.gz
```

2. Merge near-equal numeric values via rounding to 6 decimals
```bash
python merge_bins_bdg.py -i track.raw.bdg -o track.rnd6.bdg --rnd 6
```

3. Merge with numeric tolerance 1e-6 (with no rounding of the printed token)
```bash
python merge_bins_bdg.py -i track.raw.bdg -o track.eps1e-6.bdg --eps 1e-6
```

4. Stream from stdin to stdout, keeping only default header prefixes
```bash
zcat track.raw.bdg.gz \
    | python merge_bins_bdg.py -i - -o - --rnd 4 --skp_pfx "#,track,browser" \
    > track.rnd4.merged.bdg
```

5. Disable header skipping (treat all lines as data by passing an empty string
   to '--skp_pfx')
```bash
python merge_bins_bdg.py -i track.raw.bdg -o track.merged.bdg --skp_pfx ""
```

6. Stream from Samtools to signal computation to bin merger to gzip-compression
```bash
samtools view -b input.bam \
    | python compute_signal.py -i - -o - --outfmt bdg -sb 10 -me norm \
    | python merge_bins_bdg.py -i - -o - --rnd 6 \
    | gzip > sample.norm.rnd6.merged.bdg.gz
```

7. Use process substitution with pigz (bash) for fast de/compression
```bash
bash -c '
    python merge_bins_bdg.py \
        -i <(pigz -dc track.raw.bdg.gz) \
        -o >(pigz -c > track.merged.bdg.gz) \
        --eps 1e-6
'
```


General notes
-------------
- Merges adjacent bins on the same chrom when the end of one equals the start
  of the next and values “match.”
- Default “match” is exact text equality of the signal field (column 4; so,
  e.g., 0.5 and 0.500000 are different).
- Optional '--rnd <int>' round/prints to '<int>' decimals before comparing
  (robust to tiny jitter).
- Optional '--eps <float>' compares numerically within tolerance (prints using
  '--rnd <int>' if given, otherwise preserves the initial file's original
  text).
- If both '--rnd N' and '--eps E' are provided, values are first rounded to N
  decimals and then merged when |Δ| <= E on those rounded numbers. Printing
  uses N decimals.
- Tip: with rnd=N, eps < 0.5·10^-N behaves like exact rounded-string equality;
  eps ≥ 10^-N merges bins that differ by one rounding unit.
- By default, passes through 'track', 'browser', '#' lines unchanged.
- Canonicalization of non-finite tokens (NaN/±inf) is via
  'utils_bdg.canon_nonfinite'.
- Works with '.gz' in/out. Full streaming; no buffering of whole files.
- Additional general notes are described under 'Description', 'Usage', and
  'Output' above.


Performance notes
-----------------
- Streaming, single-pass algorithm. Memory is O(1) with respect to file size.
- CPU cost is O(N) over intervals with a constant amount of work per record.
- If there’s any bottleneck, it’s typically I/O, and gzip input/output adds
  decompression/compression overhead.
- For increased throughput, colocate data
  on fast storage and, when possible, avoid double-compressing pipes.
- Rounding and numeric parsing are only applied to column 4; header detection
  is prefix-based and O(1) per line.


#TODO
-----
- Test with, e.g., a 5-line bedGraph with exact ties, near-ties, and a
  non-finite token to validate all three modes and the combined mode.
- Test with one header line mid-stream to confirm it flushes and re-starts a
  run.
"""

from __future__ import annotations

from contextlib import redirect_stdout

from scripts.functions.utils_bdg import canon_nonfinite, try_float
from scripts.functions.utils_check import (
    check_cmp, check_exists, check_parse_outfile, check_writable
)
from scripts.functions.utils_cli import add_help_cap, CapArgumentParser
from scripts.functions.utils_io import (
    DEF_SKP_PFX, is_header, open_in, open_out, parse_skp_pfx
)

from typing import TextIO

import argparse
import math
import signal
import sys

try:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except (AttributeError, ValueError):
    pass

assert sys.version_info >= (3, 10), "Python >= 3.10 required."


def format_value(
    tok: str, rnd: int | None, cache: float | None
) -> str:
    """
    Decide what to print for a run’s value.

    Args:
        tok : str
            Original token from column 4 for the first bin in the run.
        rnd : int | None
            If not None, number of decimals to print when a finite numeric
             cache is available.
        cache : float | None
            Parsed float for the first token in the run (if parseable).

    Returns:
        tok | nf : str
            String to print in column 4 for the merged interval.

    Raises:
        None.

    Notes:
        - If 'rnd' is set and 'cache' is finite, print formatted to 'rnd'.
        - Else if the original token is a recognized non-finite, return the
          canonical token.
        - Else return the original token unchanged.
    """
    nf = canon_nonfinite(tok)
    if (rnd is not None) and (cache is not None) and math.isfinite(cache):
        v = round(cache, rnd)
        if v == 0.0:  # Avoid printing “-0.0”
            v = 0.0
        return f"{v:.{rnd}f}"

    if nf is not None:
        return nf

    return tok


def merge_bins(
    fil_in: str,
    fil_out: str,
    *,
    rnd: int | None = None,
    eps: float | None = None,
    skp_pfx: tuple[str, ...]
):
    """
    Stream through a (sorted) bedGraph, merging adjacent bins when values
    match.

    Args:
        fil_in : str
            Input bedGraph path ('.gz' ok) or '-' for stdin.
        fil_out : str
            Output bedGraph path ('.gz' ok) or '-' for stdout.
        rnd : int | None = None
            If not None, compare by rounded-to-N-decimal strings and print to N
            decimals.
        eps : float | None = None
            If not None and 'rnd' is None, compare numerically within
            |Δ| <= eps.
        skp_pfx : tuple[str, ...]
            Header/metadata prefixes to skip (after left-stripping).

    Returns:
        None (writes merged bedGraph to 'fil_out').

    Raises:
        ValueError
            If '--eps' < 0, or on malformed data lines (non-integer coords,
            non-positive width, or < 4 fields). Headers pass through.

    Notes:
        - Matching modes (priority):
            + if rnd is not None: compare on rounded string to 'rnd' places.
            + elif eps is not None: compare numeric values within |Δ| <= eps.
            + else: compare exact text of the 4th field (after canonicalizing
              nan/inf).
        - Non-finite values ('nan', 'inf', '-inf') never merge with anything
          else.
    """
    #  MAYBE: these checks are duplicated in 'main()'; cut them?
    #  (Silence but retain for potential programmatic calls)
    # if rnd is not None and rnd < 0:
    #     raise ValueError("Error: 'rnd' must be >= 0.")
    # if eps is not None and eps < 0:
    #     raise ValueError("Error: 'eps' must be >= 0.")

    #  Set the “rolling state” for the current merge “run;” a run extends while
    #  (a) chrom stays the same, (b) bins are contiguous, and (c) the value
    #  matches under the chosen rule (rounded/eps/text)
    last_chrom: str | None = None
    last_start: int | None = None
    last_end: int | None = None
    last_tok: str | None = None    # Original token for printing if needed
    last_num: float | None = None  # Numeric cache (float) if parseable
    last_key: tuple | None = None  # ("rounded", str) | ("eps", float) | ("text", str) | ("nonfin", "nan"/"inf"/"-inf")
    last_nonfin: bool = False      # True if last value is non-finite

    def flush(fo: TextIO) -> None:
        """
        Emit the current run (if any) to 'fo' and reset the rolling state.
        """
        nonlocal last_chrom, last_start, last_end
        nonlocal last_tok, last_num, last_key, last_nonfin
        if last_chrom is None:
            return
        out_val = format_value(last_tok, rnd, last_num)
        fo.write(f"{last_chrom}\t{last_start}\t{last_end}\t{out_val}\n")
        last_chrom = last_start = last_end = last_tok = None
        last_num = None
        last_key = None
        last_nonfin = False

    #  Stream input to output; headers are forwarded verbatim and terminate
    #  runs
    with open_in(fil_in) as fi, open_out(fil_out) as fo:
        for raw in fi:
            #  Headers pass through and break runs
            if is_header(raw, skp_pfx):
                flush(fo)
                fo.write(raw)
                continue

            #  Blank/whitespace-only lines also pass through and break runs
            if not raw.strip():
                flush(fo)
                fo.write(raw)
                continue

            #  Validate data lines (need at least 4 fields)
            parts = raw.rstrip("\n").split()
            if len(parts) < 4:
                flush(fo)
                raise ValueError(
                    f"Malformed bedGraph line; needs >=4 fields: {raw!r}"
                )

            #  Parse and check coordinates: reject non-integers and
            #  non-positive spans, as bedGraph requires start < end
            chrom, s_str, e_str, v_str = parts[0], parts[1], parts[2], parts[3]
            try:
                s = int(s_str)
                e = int(e_str)
            except ValueError:
                flush(fo)
                raise ValueError(f"Non-integer start/end: {raw!r}")
            if e <= s:
                flush(fo)
                raise ValueError(f"Non-positive interval width: {raw!r}")

            #  Build comparison key
            nf = canon_nonfinite(v_str)
            v_num = None if nf is not None else try_float(v_str)

            # MAYBE: extract the “build comparison key” block into a small
            #        helper function
            if nf is not None:
                #  Never merges across lines
                cmp_key = ("nonfin",)
            elif (
                (rnd is not None) and
                (v_num is not None) and
                math.isfinite(v_num)
            ):
                #  Round once, using the rounded value for both comparisons and
                #  printing
                v_rnd = round(v_num, rnd)
                if eps is not None:
                    #  Both '--rnd' and '--eps': compare on rounded floats
                    #  within eps
                    cmp_key = ("rnd_eps", v_rnd)
                else:
                    # '--rnd' only: compare on rounded-string equality
                    cmp_key = ("rounded", f"{v_rnd:.{rnd}f}")
            else:
                if (
                    (eps is not None) and
                    (v_num is not None) and
                    math.isfinite(v_num)
                ):
                    #  '--eps' only: Compare on raw numeric values
                    cmp_key = ("eps", v_num)
                else:
                    #  Default: Compare on exact token text (non-finites
                    #  already handled)
                    cmp_key = ("text", v_str)

            #  Start a new run if (a) chrom changes, (b) bins are not
            #  contiguous, or (c) the mode-dependent comparison key no longer
            #  matches
            start_new = (
                last_chrom is None or
                chrom != last_chrom or
                s != last_end or
                nf is not None or last_nonfin or
                last_key is None or
                cmp_key[0] != last_key[0] or
                (
                    cmp_key[0] == "rounded" and
                    last_key[1] != cmp_key[1]
                ) or
                (
                    cmp_key[0] == "text" and
                    last_key[1] != cmp_key[1]
                ) or
                (
                    #  MAYBE: 'last_key[0] == "eps"' redundant? Cut it?
                    cmp_key[0] == "eps" and
                    last_key[0] == "eps" and
                    abs(last_key[1] - cmp_key[1]) > (eps or 0.0)
                ) or
                (
                    #  MAYBE: 'last_key[0] == "rnd_eps"' redundant? Cut it?
                    cmp_key[0] == "rnd_eps" and
                    last_key[0] == "rnd_eps" and
                    abs(last_key[1] - cmp_key[1]) > (eps or 0.0)
                )
            )

            if start_new:
                flush(fo)
                last_chrom = chrom
                last_start = s
                last_end = e
                last_tok = v_str
                last_num = v_num
                last_nonfin = (nf is not None)
                last_key = cmp_key
            else:
                #  Extend current run
                last_end = e

        #  “Flush” the tail
        flush(fo)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = CapArgumentParser(
        description="Merge adjacent same-valued bedGraph bins (streaming)."
    )
    add_help_cap(parser)
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        default=False,
        help="Increase output verbosity.\n\n"
    )
    parser.add_argument(
        "-i", "--infile",
        required=True,
        help=(
            "Path to the input bedGraph file (.gz is handled), or '-' for "
            "stdin.\n\n"
        )
    )
    parser.add_argument(
        "-o", "--outfile",
        required=True,
        help=(
            "Path to the output bedGraph file (.gz is handled), or '-' for "
            "stdout.\n\n"
        )
    )
    parser.add_argument(
        "-r", "--rnd",
        type=int,
        default=None,
        help=(
            "Number of decimals in the printed total.\n"
            "\n"
            "Controls the value comparisons by N-decimal rounded strings. "
            "With '--eps', rounding is applied first, then |Δ| (absolute "
            "difference) <= eps is tested on the rounded numbers.\n\n"
        )
    )
    parser.add_argument(
        "-e", "--eps",
        type=float,
        default=None,
        help=(
            "Numeric tolerance for matching (merge if |Δ| <= eps).\n"
            "  - Alone: compare raw numeric values with |Δ| <= eps.\n"
            "  - With '--rnd N': compare on rounded values (to N decimals) "
            "with |Δ| <= eps.\n\n"
        )
    )
    parser.add_argument(
        "-sp", "--skp_pfx",
        type=str,
        default=",".join(DEF_SKP_PFX),
        help=(
            "Comma-separated list of bedGraph prefixes to skip as headers/"
            "metadata; to disable skipping, pass an empty string (default: "
            "%(default)s).\n\n"
        )
    )

    argv_parse = sys.argv[1:] if argv is None else argv
    if not argv_parse:
        parser.print_help(sys.stderr)
        sys.exit(0)

    return parser.parse_args(argv_parse)


def main(argv: list[str] | None = None) -> int:
    """
    Execute the primary control flow for the script.

    Args:
        argv : list[str] | None = None
            Optional argument vector for testing; defaults to sys.argv[1:].

    Returns:
        Process exit code (int).

    Raises:
        ValueError
            Propagated from 'merge_bins' when input data lines are malformed.
    """
    args = parse_args(argv)

    #  Perform argument checks
    try:
        if args.infile != "-":
            check_exists(args.infile, "file", "bedGraph")

        if args.outfile != "-":
            outfile, _, _ = check_parse_outfile(
                args.outfile, ("bedgraph", "bdg", "bg")
            )
            check_writable(outfile, "file")
        else:
            outfile = "-"

        #  MAYBE: 'rnd' check in 'main()' and 'merge_bins()'; pick one?
        check_cmp(args.rnd, "ge", 0, "rnd", allow_none=True)

        #  MAYBE: 'eps' check in 'main()' and 'merge_bins()'; pick one?
        check_cmp(args.eps, "ge", 0, "eps", allow_none=True)

    except ValueError as e:
        raise SystemExit(str(e))

    #  Standardize header-skip prefixes
    skp_pfx = parse_skp_pfx(args.skp_pfx, default=DEF_SKP_PFX)

    #  Print verbose output
    if args.verbose:
        with redirect_stdout(sys.stderr):
            print("####################################")
            print("## Arguments for 'merge_bins_bdg' ##")
            print("####################################")
            print("")
            print(f"--verbose {args.verbose}")
            print(f"--infile  {args.infile}")
            print(f"--outfile {outfile}")
            if args.rnd is not None:
                print(f"--rnd     {args.rnd}")
            if args.eps is not None:
                print(f"--eps     {args.eps}")
            if args.skp_pfx:
                print(f"--skp_pfx {skp_pfx}")
            print("")
            print("")

    try:
        merge_bins(
            fil_in=args.infile,
            fil_out=args.outfile,
            rnd=args.rnd,
            eps=args.eps,
            skp_pfx=skp_pfx
        )
    except (ValueError, FileNotFoundError, PermissionError, OSError) as e:
        print(str(e), file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    try:
        sys.exit(main())
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
