#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: compute_signal_ratio.py
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
compute_signal_ratio.py


Description
-----------
Compute per-bin ratios between two bedGraph tracks [e.g., (IP ÷ input),
log2(IP ÷ input), and reciprocals thereof] with a deepTools bamCompare-like
order of operations.


Usage
-----
python -m compute_signal_ratio \
    [--verbose] \
    --fil_A <str_file_A.(bedGraph|bedgraph|bdg|bg)(.gz)|-> \
    --fil_B <str_file_B.(bedGraph|bedgraph|bdg|bg)(.gz)|-> \
    --fil_out <str_file.(bedGraph|bedgraph|bdg|bg)(.gz)|-> \
    [--method <str_method>] \
    [--scl_fct <flt_A>[:<flt_B>]] \
    [--pseudo <flt_A>[:<flt_B>]] \
    [--dep_min <flt_floor>] \
    [--skip_00 {pre_scale,post_scale}] \
    [--eps <flt_zero_tol>] \
    [--drp_nan] \
    [--track] \
    [--rnd <int_dec>] \
    [--skp_pfx <str,...>]


Arguments
---------
 -v, --verbose
    Run script in verbose mode (stderr banner of parsed args).

-fA, --fil_A <str>
    Path to first bedGraph input (file A; e.g., IP), or '-' to read from
    stdin. Gzip-compressed input ('.gz') is supported.

-fB, --fil_B <str>
    Path to second bedGraph input (file B; e.g., input), or '-' to read from
    stdin. Gzip-compressed input ('.gz') is supported.

-fo, --fil_out <str>
    Path to output bedGraph, or '-' to write to stdout. If a real output path
    ends with '.gz', output is gzip-compressed. (stdout output is intended for
    direct use of 'compute_signal_ratio.py'; the higher-level Shell wrappers do
    not support this mode.)

-me, --method <str>
    Ratio-computation subtype. Accepted aliases are:
        - 'r', 'raw', 'u', 'unadj', 'unadjusted', 's', 'smp', 'simple': compute
          the simple unadjusted ratio 'A / B'. Internally standardized to
          'method=unadj'.
        - '2', 'l2', 'lg2', 'log2': compute 'log2(A / B)'. Internally
          standardized to 'method=log2'.
        - 'rr', 'raw_r', 'ur', 'unadj_r', 'unadjusted_r', 'sr', 'smp_r',
          'simple_r': compute the reciprocal of the simple ratio, i.e.,
          'B / A'. Internally standardized to 'method=unadj_r'.
        - '2r', 'l2r', 'l2_r', 'lg2_r', 'log2_r': compute the reciprocal of the
          log2 ratio, i.e., 'log2(B / A)' = '-log2(A / B)'. Internally
          standardized to 'method=log2_r'.

-sf, --scl_fct <flt_A>[:<flt_B>]
    Optional per-file multiplicative scaling factors (default: None [treated as
    '1.0:1.0']). If only A is given, then 'B=1.0'. Must be > 0. Examples:
    '1.0:1.0', '1.3:0.95', etc.

-ps, --pseudo <flt_A>[:<flt_B>]
    Optional per-file pseudocount summands added after any scaling (default:
    0:0). Primarily useful for log2-ratio methods, where it helps avoid
    undefined values such as 'log2(A / 0)'. For linear ratios, '--dep_min' may
    be preferable when the main concern is bounding low-depth extremes. Also,
    using '--pseudo' together with '--dep_min' is allowed but usually makes
    low-depth stabilization harder to interpret.

-dm, --dep_min <flt>
    Optional minimum denominator clamp ("minimum input depth"; default: None
    [no clamp]). Applied after any scaling and/or pseudocounts:
    'B := max(B, dep_min)'. After clamping, bins with '|B| <= ε' yield 'nan'.
    In practice choose 'dep_min > ε'. Also, using '--dep_min' together with
    '--pseudo' is allowed, but can make low-depth stabilization harder to
    interpret.

 -e, --eps <flt>
    Epsilon (ε) tolerance for zero checks (default: 0.0 [deepTools-like]).
    Values satisfying '|value| <= ε' are treated as zero. Used in two places:
    (i) for '--skip_00' zero/zero tests, and (ii) as a final zero guard on the
    denominator after optional clamping.

-s0, --skip_00 {pre_scale,post_scale}
    Optionally drop bins where both values are treated as zero ('|A| <= ε' and
    '|B| <= ε').

    Modes:
        - pre_scale   Test raw values before scaling and pseudocount addition.
                      This is deepTools-like behavior.
        - post_scale  Test after scaling but before pseudocount addition; ε is
                      interpreted in scaled units when scaling ≠ 1.

    Notes:
        - If ε = 0 or scaling is neutral (1:1), modes are equivalent.
        - Pseudocounts do not affect this check.

-dn, --drp_nan
    Drop non-finite values ('nan', 'inf', '-inf') from the main output file.
    The '.track' sidecar (if requested via '--track') always drops non-finite
    values. Default: False.

-tr, --track
    Also write a '.track' sidecar bedGraph that omits non-finite values ('nan',
    'inf', '-inf'), which is convenient for genome browsers (e.g., IGV).

-dp, --dp, --rnd, --round, --decimals, --digits <int>
    Maximum number of decimal places retained for emitted numeric values in
    finite rows. After rounding, non-informative trailing zeros are stripped.
    Default: 24.

-sp, --skp_pfx <csv>
    Comma-separated prefixes to treat as header/metadata lines in bedGraph
    tracks. Lines that start with any listed prefix are skipped. Default:
    "#,track,browser".


Output
------
- bedGraph ('.bedGraph', '.bedGraph.gz', '.bedgraph', '.bedgraph.gz', '.bdg',
  '.bdg.gz', '.bg', '.bg.gz'):
    + One line per emitted bin: chrom  start  end  value.
    + Values are ratios (A / B), log2(A / B), reciprocals, or reciprocal log2
      depending on '--method'.
    + Finite values are rounded to at most '--rnd' decimal places; after
      rounding, non-informative trailing zeros and any trailing decimal point
      are stripped.
    + Non-finite values are written as 'nan', 'inf', or '-inf' unless
      '--drp_nan' is set.
- Sidecar track (when '--track'):
    + Same schema as the main bedGraph.
    + All non-finite values are omitted unconditionally for genome browser
      compatibility.
    + Filename is derived by inserting '.track' before the main extension and
      preserving '.gz' if present.
- Compression:
    + If the output path ends with '.gz', the corresponding file(s) are
      gzip-compressed in text mode.


Examples
--------
1. Compute a log2 ratio track with zero/zero skipping after scaling (contrast
   with deepTools-like pre-scale behavior), and gzip the output
'''bash
python -m scripts.compute_signal_ratio \
    -fA IP_sample.sc.bdg.gz \
    -fB in_sample.sc.bdg.gz \
    -fo log2_sample.sc.bdg.gz \
    -me log2 \
    -s0 post_scale
'''

2. Compute a linear ratio track with a denominator clamp and write a '.track'
   sidecar
'''bash
python -m scripts.compute_signal_ratio \
    -fA IP_sample.sc.bdg \
    -fB in_sample.sc.bdg \
    -fo ratio_sample.sc.bedGraph \
    -dm 5.0 \
    -tr
'''

3. Compute a log2 ratio track with distinct scaling factors and small
   pseudocounts
'''bash
python -m scripts.compute_signal_ratio \
    -fA IP_sample.sc.bdg.gz \
    -fB in_sample.sc.bdg.gz \
    -fo log2_scaled.sc.bdg.gz \
    -sf 1.25:0.93 \
    -ps 1e-6:1e-6 \
    -me log2
'''

4. Return a track that is the reciprocal of the linear ratio, dropping
   non-finite values from the main output file
'''bash
python -m scripts.compute_signal_ratio \
    -fA IP_sample.sc.bdg \
    -fB in_sample.sc.bdg \
    -fo recip_ratio.sc.bedgraph \
    -me unadj_r \
    -dn
'''

5. Compute the reciprocal of a log2 ratio track [i.e., 'log2(B / A)' =
   '-log2(A / B)']
'''bash
python -m scripts.compute_signal_ratio \
    -fA IP_sample.sc.bdg.gz \
    -fB in_sample.sc.bdg.gz \
    -fo log2_recip_sample.sc.bdg.gz \
    -me log2_r
'''

6. Skip custom header prefixes (e.g., unusual bedGraph headers)
'''bash
python -m scripts.compute_signal_ratio \
    -fA A.bdg -fB B.bdg -fo out.bdg \
    -sp "#,track,browser,customHeader"
'''

7. Example using stdin and writing to stdout (not runnable as written; not
   supported by Shell wrappers)
'''bash
#  Example: upstream 'compute_signal' call abbreviated
samtools view -b input.bam \
    | python -m scripts.compute_signal ... \
    | python -m scripts.compute_signal_ratio \
        -fA - \
        -fB input.sc.bdg.gz \
        -fo -
'''


Order of operations (per bin)
-----------------------------
1. Optional skip of zero-zero bins ('0 / 0' or 'ε / ε') before (optional)
   scaling ‡
    - This step precedes optional multiplicative scaling.
    - If '--skip_00 pre_scale', drop bins where both raw values are treated
      as zero: '|A| <= ε' and '|B| <= ε' (ε set by '--eps'; default 0.0 for
      exact zero).
    - This mirrors the intent of deepTools 'bamCompare --skipZeroOverZero',
      which is to ignore unscaled bins in both the dividend and divisor with
      signal values of zero.
    - '--eps' ("epsilon" or ε) controls zero tolerance [default: 0.0 (exact
      zero)]; e.g., the user can set a small ε (e.g., 1e-12) to treat float
      noise as zero.

2. Optional per-file scaling
    - Multiply the first file (file A; e.g., "IP") by scaling factor A, and
      multiply the second file (file B; e.g., "input") by scaling factor B (via
      '--scl_fct A[:B]').
    - If only A is provided, B defaults to 1.0.

3. Optional skip of zero-zero bins ('0 / 0' or 'ε / ε') after (optional)
   scaling ‡
    - This step takes place after optional multiplicative scaling but prior to
      optional pseudocount addition.
    - If '--skip_00 post_scale', then drop bins where both scaled values are
      treated as zero: '|scl_A × A| <= ε' and '|scl_B × B| <= ε' (ε set by
      '--eps'; default 0.0 for exact zero).
    - With 'ε > 0' and 'scaling ≠ 1', this interprets ε in scaled units (i.e.,
      the units that are actually divided).

4. Optional per-file pseudocount addition
    - Add pseudocount A to file A (numerator; e.g., IP; may or may not be
      scaled) and add pseudocount B to file B (denominator; e.g., input; may or
      may not be scaled) via '--pseudo A[:B]'.
    - Pseudocounts default to '0:0'.
    - Pseudocounts are applied after optional scaling and optional skipping of
      true-zero bins.
    - Pseudocounts are used mainly for log2 ratios, where they prevent
      undefined values such as 'log2(A / 0)', 'log2(0 / B)', and 'log2(0 / 0)'.
    - In log2-ratio mode,
        + pseudocounts address the zero-value problem directly and are a
          standard way to make log ratios finite.
        + For log2 ratios, '--pseudo A[:B]' is often preferable to
          '--dep_min <flt>' (for linear ratios, '--dep_min <flt>' may be
          preferable to '--pseudo A[:B]').


5. Optional denominator clamp (i.e., threshold or floor)
    - If '--dep_min <flt>', 'B := max(B, dep_min)' after any optional
      multiplicative scaling and pseudocount addition.
    - This is the "minimum input depth" behavior described in PMID 40364978. It
      thresholds very small denominators, preventing extreme divisions such as
      'A / ("very small B")' and 'A / 0', thereby stabilizing ratios.
    - In linear-ratio mode,
        + pseudocounts can alter low-depth bins, sometimes substantially, by
          shifting both numerator and denominator away from their observed
          values.
        + When the main concern is preventing extreme linear ratios caused by
          very small denominator values, '--dep_min <flt>' may be preferable to
          '--pseudo A[:B]' (by contrast, for log2 ratios, '--pseudo A[:B]' is
          often preferable to '--dep_min <flt>').
    - Unlike the other operations here, the use of an optional denominator
      clamp is not available in deepTools bamCompare.

6. A zero guard is applied: If '|B| <= ε', the bin yields 'nan'.
    - This takes place after any optional denominator clamping.
    - In practice, choose 'dep_min > ε'.

7. Division
    - Compute 'ratio = A / B' (after optional steps 1–6).

8. Optional log2 transformation
    - If '--method' selects a log2-based ratio, transform to 'log2(A / B)'.
      ('log2(0)' becomes '-inf'; negative ratio becomes 'nan')

9. Optional reciprocal computation
    - If '--method' selects a reciprocal linear ratio, return '1 / (A / B)'.
    - If '--method' selects a reciprocal log2 ratio, return '-log2(A / B)'.

‡ If specified, only optional #1 OR optional #3 can be applied, not both.


I/O handling
------------
- In the streaming merge, missing bins are treated as 0.0 for missing data.
- Main output may optionally drop non-finite values ('nan', 'inf', '-inf') with
  '--drp_nan'.
- Input files may be regular paths or '-' for stdin. At most one of '--fil_A'
  or '--fil_B' may be '-', as stdin can only be consumed once.
- When either input is '-', the light bin-size consistency check (see 'General
  notes') is skipped, as stdin cannot be rewound without buffering. In that
  case, the script assumes that A and B are already on a compatible binning
  grid (same chromosomes, bin starts, and widths).
- The output can be written to a path or to stdout using '--fil_out -'; in all
  cases the format is bedGraph. (This stdout workflow is intended for direct
  use of 'compute_signal_ratio.py'; the higher-level Shell wrappers do not
  support stdout output.)
- If '--track' is given, a '.track' sidecar is written that excludes 'nan',
  'inf', and '-inf', supporting usability with genome browsers such as IGV.
  This requires a real outfile path; it is not available when '--fil_out -' is
  used.


General notes
-------------
- Coordinates are 0-based, half-open '[start, end)' in bedGraph.
- Chromosome labels are preserved; e.g., no automatic 'chr' standardization.
- Sorting for emission uses Roman-numeral ordering for I..XVI first (case-
  insensitive, with or without 'chr' prefix), then lexical for all others.
- Streaming merge treats missing bins in either file as '0.0' for the given
  file.
- Header/meta lines are skipped by prefix matching (configurable via
  '--skp_pfx').
- A light bin-size consistency check compares the first few data lines of both
  files and fails fast on mismatch or malformed rows. This check is skipped
  when either input is read from stdin ('-'), as stdin cannot be safely read
  twice without buffering.
- On the handling of non-finite values:
    + Main output: written as 'nan', 'inf', '-inf' unless '--drp_nan'.
    + Track output (with '--track'): non-finite are always omitted.
- Emission order reflects the merge order defined by 'key_bin' (Roman I..XVI
  first, then numeric, then mitochondria, then lexical) assuming both inputs
  are individually sorted by the same order.
- The streaming merge matches bins by chromosome and start coordinate (via
  'key_bin'). Thus, it assumes that corresponding bins in A and B also share
  the same end coordinate. A light upfront consistency check is performed, but
  strict whole-file validation of bin boundaries is not currently implemented.


Performance notes
-----------------
- The merge is streamed line-by-line with O(1) memory relative to file size.
- I/O is the main bottleneck.
    + Using local SSDs or some other fast storage improves throughput.
    + 'gzip' adds CPU cost (which negatively affects throughput) in exchange
      for smaller files. (Worth it, I think.)


#TODO
-----
- In 'check_bin_size()', I currently compare up to 5 paired data lines. This is
  a fast heuristic, but I should implement an option for stricter validation,
  e.g., by adding a '--strict_bins' (or something like that) flag to scan more
  or all lines, leaving the quick check as the default behavior. (However,
  currently not a priority.)
- Strengthen bin-grid validation beyond the current light upfront check. The
  streaming merge matches bins by chromosome and start coordinate (via
  'key_bin') and currently assumes that corresponding bins in A and B also
  share the same end coordinate. A stricter validation mode should verify this
  assumption across more or all rows before merging.
"""

from __future__ import annotations

from contextlib import (
    nullcontext,
    redirect_stdout
)

from scripts.functions.utils_bdg import (
    check_bin_size,
    generate_track_name,
    key_bin
)
from scripts.functions.utils_check import (
    check_cmp,
    check_exists,
    check_parse_outfile,
    check_writable
)
from scripts.functions.utils_cli import (
    add_help_cap,
    CapArgumentParser
)
from scripts.functions.utils_io import (
    DEF_SKP_PFX,
    open_in,
    open_out,
    parse_skp_pfx,
    read_data_line
)

import argparse
import math
import signal
import sys

try:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except (AttributeError, ValueError):
    pass

assert sys.version_info >= (3, 10), "Python >= 3.10 required."

#  Set allowed outfile extensions
ONLY_BDG = ("bedGraph", "bedgraph", "bdg", "bg")

#  Accepted '--method' values and their canonical internal names
METHOD_CANON = {
    #  Simple unadjusted ratio: A / B
    "r": "unadj",
    "raw": "unadj",
    "u": "unadj",
    "unadj": "unadj",
    "unadjusted": "unadj",
    "s": "unadj",
    "smp": "unadj",
    "simple": "unadj",

    #  Log2 ratio: log2(A / B)
    "2": "log2",
    "l2": "log2",
    "lg2": "log2",
    "log2": "log2",

    #  Reciprocal of simple ratio: B / A
    "rr": "unadj_r",
    "raw_r": "unadj_r",
    "ur": "unadj_r",
    "unadj_r": "unadj_r",
    "unadjusted_r": "unadj_r",
    "sr": "unadj_r",
    "smp_r": "unadj_r",
    "simple_r": "unadj_r",

    #  Reciprocal of log2 ratio: log2(B / A) = -log2(A / B)
    "2r": "log2_r",
    "l2r": "log2_r",
    "l2_r": "log2_r",
    "lg2_r": "log2_r",
    "log2_r": "log2_r",
}
METHOD_CHOICES = tuple(METHOD_CANON.keys())


def parse_pair(val: str, def_sec: float) -> tuple[float, float]:
    """
    Parse 'A' or 'A:B' into two floats. If only 'A' is provided, the second
    value defaults to a user-assigned default second value, 'def_sec'.

    Args:
        val : str
            String, e.g., '2.0', '2.0:1.5', etc.
        def_sec : float
            Default for the second value when only 'A' is provided.

    Returns:
        (float, float)
            Tuple (A, B).

    Raises:
        argparse.ArgumentTypeError
            If the string is not of the form 'A' or 'A:B' with numeric parts.
    """
    parts = val.split(":")
    if len(parts) == 1:
        try:
            return float(parts[0]), def_sec
        except ValueError as e:
            raise argparse.ArgumentTypeError("Expected number for 'A'.") from e
    if len(parts) == 2:
        try:
            return float(parts[0]), float(parts[1])
        except ValueError as e:
            raise argparse.ArgumentTypeError("Expected numeric 'A:B'.") from e
    raise argparse.ArgumentTypeError("Expected 'A' or 'A:B'.")


def calc_rat_bin(
    sig_A: float,
    sig_B: float,
    scl_A: float | None,
    scl_B: float | None,
    psc_A: float | None,
    psc_B: float | None,
    dep_min: float | None,
    log2: bool,
    recip: bool,
    skip_00: str | None,
    eps: float = 0.0
) -> float | None:
    """
    Compute 'A / B' (optionally log2-transformed and/or the reciprocal) with
    optional per-file multiplicative scaling and optional pseudocount addition,
    as well as an optional "denominator clamp" (e.g., a "minimum input depth
    value" as described in PMID 40364978).

    See module docstring for more details.

    Args:
        sig_A : float
            First file (e.g., IP) signal for a bin (A).
        sig_B : float
            Second file (e.g., input) signal for a bin (B).
        scl_A : float | None
            Per-file multiplicative scale factor for A. If None or 1.0, treated
            as neutral.
        scl_B : float | None
            Per-file multiplicative scale factor for B. If None or 1.0, treated
            as neutral.
        psc_A : float | None
            Pseudocount added to A (post-scaling). If None or 0.0, skipped.
        psc_B : float | None
            Pseudocount added to B (post-scaling). If None or 0.0, skipped.
        dep_min : float | None
            Minimum allowed denominator after any optional scaling and/or
            pseudocount addition. If provided, B := max(B, dep_min) to avoid
            extreme and undefined (e.g., 'n / 0') divisions.
        log2 : bool
            If True, return 'log2(A / B)'. If False, return linear 'A / B'.
        recip : bool
            If True, return the reciprocal of the computed ratio:
                - linear: '1 / (A / B)'
                - log2: '-log2(A / B)' (since 'log2(1 / x) = -log2(x)')
        skip_00 : str | None
            Optional zero-bin ('0 / 0' or 'ε / ε') drop stage. One of
            "pre_scale", "post_scale", or None.
                - "pre_scale": Test on raw values (before optional scaling
                               and/or pseudocounts addition).
                - "post_scale": Test after optional scaling, before optional
                                pseudocount addition.
                - None: Do not drop '0 / 0' bins.
        eps : float, default 0.0
            Tolerance (epsilon value, ε) for treating values as zero in the
            pre-pseudocount zero-zero check. Use 0.0 for exact-zero behavior,
            similar to what is done in deepTools. A tiny value (e.g., say,
            '1e-12') can be used to ignore "float noise."

    Returns:
        ratio | xfrm | -xfrm : float | None
            The computed value (ratio or log2 ratio, possibly reciprocated).

            May also return the following:
                - None  When 'skip_00' is specified and the bin is '0 / 0' or
                        'ε / ε'.
                - -inf  When 'log2' is requested and 'ratio == 0'.
                -  inf  When 'reciprocal' is requested on a zero linear ratio.
                -  nan  When the computation is undefined (e.g., negative ratio
                        for 'log2', or 'B == 0' even after clamping).

            Caller skips writing the bin when a return value is None.

    Raises:
        None.
            This function is defensive: domain issues are mapped to sentinel
            return values instead of exceptions. Specifically,
                - '0 / 0' (or 'ε / ε') bins may return None (i.e., when
                  '--skip_00' applies).
                - denominator underflow after optional clamping yields
                  'float('nan')'.
                - with 'log2=True', 'ratio == 0' yields '-inf', and 'ratio < 0'
                  yields 'nan'.
                - with 'log2=False' and 'recip=True', 'ratio == 0' yields
                  'inf'.

            Callers should validate option domains (e.g., 'scl_fct > 0',
            'eps >= 0') prior to calling. A TypeError may still propagate if
            non-numeric inputs are passed.

    Notes:
        Order of operations:
            1. Optionally skip zero-zero bins: '0 / 0' or 'ε / ε' (deepTools-
               like). ‡
            2. Optionally scale each file.
            3. Optionally skip scaled zero-zero bins:
               '[(sf_A × 0) / (sf_B × 0)]' or '[(sf_A × ε) / (sf_B × ε)]'. ‡
            4. Optionally add pseudocounts.
            5. Optionally clamp denominator by 'dep_min'.
            6. If the denominator is <= ε, treat the bin as undefined.
            7. Divide 'A / B'.
            8. Optionally perform log2 transformation.
            9. Optionally compute reciprocal of #7 (linear) or #8 (log2).

            ‡ Either of optional #1 or optional #3, not both.
    """
    #  1. Optionally skip on '0 / 0' or 'ε / ε' prior to optional scaling
    if skip_00 == "pre_scale":
        if abs(sig_A) <= eps and abs(sig_B) <= eps:
            return None

    #  2. Optionally scale each file, skipping neutral factors
    num = sig_A if (scl_A is None or scl_A == 1.0) else (scl_A * sig_A)
    den = sig_B if (scl_B is None or scl_B == 1.0) else (scl_B * sig_B)

    #  3. Optionally skip on '0 / 0' or 'ε / ε' after optional scaling, prior
    #     to optional pseudocount addition
    if skip_00 == "post_scale":
        if abs(num) <= eps and abs(den) <= eps:
            return None

    #  4. Optionally add pseudocounts, skipping neutral summands
    if psc_A not in (None, 0.0):
        num += psc_A
    if psc_B not in (None, 0.0):
        den += psc_B

    #  5. Optionally clamp the denominator (if a minimum denominator threshold
    #     is specified)
    if dep_min is not None:
        if den < dep_min:
            den = dep_min

    #  6. Perform a zero guard: If the denominator (optionally clamped or not)
    #     remains effectively zero ('|den| <= ε'), treat the bin as undefined
    if abs(den) <= eps:  # Advice: pick 'dep_min > ε' so this never triggers
        return float("nan")

    #  7. Compute ratio
    ratio = num / den

    #  8. Optionally compute a log2 transformation of the ratio
    if log2:
        if ratio > 0.0:
            xfrm = math.log2(ratio)
        elif ratio == 0.0:
            xfrm = float("-inf")
        else:
            xfrm = float("nan")

        #  9. Optionally compute a reciprocal on the log2 value (log space)
        return -xfrm if recip else xfrm

    #  9. Optionally compute a reciprocal on the ratio (linear space)
    if recip:
        return (1.0 / ratio) if ratio != 0.0 else float("inf")

    return ratio


def comp_sig_rat(
    fil_A: str,
    fil_B: str,
    fil_out: str,
    scl_A: float,
    scl_B: float,
    psc_A: float,
    psc_B: float,
    dep_min: float | None,
    rnd: int,
    log2: bool,
    recip: bool,
    skip_00: str | None,
    eps: float,
    track: bool,
    drp_nan: bool,
    skp_pfx: tuple[str, ...] = DEF_SKP_PFX
) -> None:
    """
    Compute the ratio, log2 ratio, reciprocal ratio, or reciprocal log2 ratio
    of the first ('A'; e.g., IP) and second ('B'; e.g., input) bedGraph tracks,
    with optional per-file multiplicative scaling, addition of pseudocount(s),
    and denominator clamping (e.g., "minimum input depth"). Streams the merge.
    Handles cases where bins are missing in one of the files. The merge matches
    bins by chromosome and start coordinate and therefore assumes that matched
    bins also share the same end coordinate.

    See module docstring for more details.

    Args:
        fil_A : str
            Path to first (e.g., IP) bedGraph file (can be gzipped; A).
        fil_B : str
            Path to second (e.g., input) bedGraph file (can be gzipped; B).
        fil_out : str
            Output file path (if '.gz' extension, gzip compression).
        scl_A : float
            Scale factor for A.
        scl_B : float
            Scale factor for B.
        psc_A : float
            Pseudocount added to A after scaling (0.0 means no addition).
        psc_B : float
            Pseudocount added to B after scaling (0.0 means no addition).
        dep_min : float | None
            Denominator clamp ("minimum input depth") to avoid extreme and/or
            erroneous division.
        rnd : int
            Maximum number of decimal places retained for finite emitted
            values.
        log2 : bool
            If True, return 'log2(A / B)'. If False, return linear 'A / B'.
        recip : bool
            If True, return reciprocal of the final value:
                - linear path: '1 / (A / B)'
                - log2 path: '-log2(A / B)'
        skip_00 : str | None
            Optional zero-bin ('0 / 0' or 'ε / ε') drop stage. One of
            "pre_scale", "post_scale", or None.
                - "pre_scale": Test on raw values (before optional scaling
                               and/or pseudocounts addition).
                - "post_scale": Test after optional scaling, before optional
                                pseudocount addition.
                - None: Do not drop '0 / 0' bins.
        eps : float
            Epsilon for zero tests used in the pre-optional-pseudocount '0 / 0'
            drop and for post-clamp denominator zero checks.
        track : bool
            If True, write a '.track' sidecar omitting non-finite values
            ('inf', '-inf', and 'nan').
        drp_nan : bool
            If True, omit non-finite values from the main output as well.
        skp_pfx : tuple[str, ...]
            Prefixes to skip as bedGraph header/meta lines.

    Returns:
        None. Writes one or two bedGraphs to disk.
            - Standard bedGraph file: Contains all computed ratios unless
              '--drp_nan' is invoked, in which case rows with 'inf', '-inf',
              and 'nan' values are omitted. Finite values are rounded to at
              most 'rnd' decimal places and then have non-informative trailing
              zeros stripped.
            - If track is 'True', a second file with '.track' before the
              extension is created, excluding rows with 'inf', '-inf', and
              'nan' values.

    Raises:
        FileNotFoundError
            If an input file or the output directory is missing.
        PermissionError
            If the output directory is not writable.
        ValueError
            On malformed bedGraph lines, inconsistent bin sizes, invalid
            numeric arguments, or non-finite formatting issues.
        OSError
            On I/O errors opening, reading, and/or writing files.
    """
    #  Check bin size consistency before proceeding
    if fil_A != "-" and fil_B != "-":
        check_bin_size(fil_A=fil_A, fil_B=fil_B, skp_pfx=skp_pfx)

    #  Generate file name for optional track (insert '.track' before the ext)
    fil_trk = generate_track_name(fil_out) if track else None

    #  Open input/output safely; ensure sidecar (if any) is closed on error
    with (
        open_in(fil_A) as opn_A,
        open_in(fil_B) as opn_B,
        open_out(fil_out) as f_out,
        (open_out(fil_trk) if track else nullcontext()) as f_trk
    ):
        def write_line(chrom, start, end, val):
            """
            Helper function: Write one bedGraph row. Optionally drop non-finite
            values ('inf', '-inf', 'nan') in main output; if enabled, always
            drop non-finite values in '.track' sidecar output.
            """
            if val is None:
                return  # Skip '0 / 0' bins (assessed prior to pseudo addition)

            #  Optionally drop non-finite values in main output:
            if drp_nan and not math.isfinite(val):
                return

            #  Safely format value for main output
            if math.isfinite(val):
                v = round(val, rnd)
                if v == 0.0:
                    v = 0.0  # Collapse "-0"

                #  Format with at most 'rnd' decimal places, then strip only
                #  non-informative trailing zeros and any trailing decimal
                #  point
                s_val = f"{v:.{rnd}f}"
                if "." in s_val:
                    s_val = s_val.rstrip("0").rstrip(".")
                if s_val == "-0":
                    s_val = "0"
            else:
                if math.isnan(val):
                    s_val = "nan"
                else:
                    s_val = "-inf" if val < 0 else "inf"

            f_out.write(f"{chrom}\t{start}\t{end}\t{s_val}\n")

            #  Write track 'sidecar': drop all non-finite values
            if (f_trk is not None) and math.isfinite(val):
                val_trk = round(val, rnd)
                if val_trk == 0.0:
                    val_trk = 0.0

                s_trk = f"{val_trk:.{rnd}f}"
                if "." in s_trk:
                    s_trk = s_trk.rstrip("0").rstrip(".")
                if s_trk == "-0":
                    s_trk = "0"

                f_trk.write(f"{chrom}\t{start}\t{end}\t{s_trk}\n")

        #  Prime the data streams
        lin_A = read_data_line(fh=opn_A, skp_pfx=skp_pfx)
        lin_B = read_data_line(fh=opn_B, skp_pfx=skp_pfx)

        #  Merge the two bedGraph streams in order; if a bin exists in only one
        #  file, we treat the missing partner as 0.0 for that bin, making it so
        #  that we can still compute 'A / B'
        while lin_A or lin_B:  # Continue until both files are exhausted
            fld_A = lin_A.split() if lin_A else None
            fld_B = lin_B.split() if lin_B else None

            #  Fail fast on malformed bedGraph rows
            if fld_A and len(fld_A) < 4:
                raise ValueError(
                    "Malformed bedGraph line in first file (e.g., IP): "
                    f"{lin_A!r}"
                )
            if fld_B and len(fld_B) < 4:
                raise ValueError(
                    "Malformed bedGraph line in second file (e.g., input): "
                    f"{lin_B!r}"
                )

            #  Parse bedGraph format, handling missing bins
            if fld_A:
                try:
                    chr_A = fld_A[0]
                    start_A = int(fld_A[1])
                    end_A = int(fld_A[2])
                    sig_A = float(fld_A[3])
                except ValueError as e:
                    raise ValueError(
                        "Non-numeric start/end/signal in first-file (e.g., "
                        f"IP) line: {lin_A!r}"
                    ) from e
            else:
                chr_A, start_A, end_A, sig_A = None, None, None, 0.0

            if fld_B:
                try:
                    chr_B = fld_B[0]
                    start_B = int(fld_B[1])
                    end_B = int(fld_B[2])
                    sig_B = float(fld_B[3])
                except ValueError as e:
                    raise ValueError(
                        "Non-numeric start/end/signal in second-file (e.g., "
                        f"input) line: {lin_B!r}"
                    ) from e
            else:
                chr_B, start_B, end_B, sig_B = None, None, None, 0.0

            #  Post-parsing sanity check for widths
            if fld_A and start_A >= end_A:
                raise ValueError(
                    "Non-positive width in first-file (e.g., IP) line: "
                    f"{lin_A!r}"
                )
            if fld_B and start_B >= end_B:
                raise ValueError(
                    "Non-positive width in second-file (e.g., input) line: "
                    f"{lin_B!r}"
                )

            #  Choose which bin to emit next
            if (fld_A is not None) and (fld_B is not None):
                # TODO: The merge key currently uses chromosome and start
                #       coordinate only (via 'key_bin'). For bins treated as
                #       matches ('key_A == key_B'), the end coordinate should
                #       also be checked. When stricter bin validation is added,
                #       enforce matching '(chrom, start, end)' across more (or
                #       all) rows rather than assuming that equal
                #       '(chrom, start)' implies a true bin match.
                key_A = key_bin(chr_A, start_A)
                key_B = key_bin(chr_B, start_B)

                if key_A == key_B:
                    #  Same bin in both files: compute ratio with A and B as-is
                    val = calc_rat_bin(
                        sig_A=sig_A,
                        sig_B=sig_B,
                        scl_A=scl_A,
                        scl_B=scl_B,
                        psc_A=psc_A,
                        psc_B=psc_B,
                        dep_min=dep_min,
                        log2=log2,
                        recip=recip,
                        skip_00=skip_00,
                        eps=eps
                    )
                    write_line(chr_A, start_A, end_A, val)
                    lin_A = read_data_line(fh=opn_A, skp_pfx=skp_pfx)
                    lin_B = read_data_line(fh=opn_B, skp_pfx=skp_pfx)

                elif key_A < key_B:
                    #  File A is "ahead", meaning this bin exists only in A at
                    #  the moment; thus, treat B as 0.0 for this bin and
                    #  compute 'A / 0' subject to guards
                    val = calc_rat_bin(
                        sig_A=sig_A,
                        sig_B=0.0,
                        scl_A=scl_A,
                        scl_B=scl_B,
                        psc_A=psc_A,
                        psc_B=psc_B,
                        dep_min=dep_min,
                        log2=log2,
                        recip=recip,
                        skip_00=skip_00,
                        eps=eps
                    )
                    write_line(chr_A, start_A, end_A, val)
                    lin_A = read_data_line(fh=opn_A, skp_pfx=skp_pfx)

                else:
                    #  File B is "ahead", meaning this bin exists only in B
                    #  right now; thus, treat A as 0.0 for this bin and compute
                    #  '0 / B' subject to guards
                    val = calc_rat_bin(
                        sig_A=0.0,
                        sig_B=sig_B,
                        scl_A=scl_A,
                        scl_B=scl_B,
                        psc_A=psc_A,
                        psc_B=psc_B,
                        dep_min=dep_min,
                        log2=log2,
                        recip=recip,
                        skip_00=skip_00,
                        eps=eps
                    )
                    write_line(chr_B, start_B, end_B, val)
                    lin_B = read_data_line(fh=opn_B, skp_pfx=skp_pfx)

            elif fld_A is not None:
                #  Only A has bins left; B is exhausted, so treat B as 0.0 for
                #  remaining bins
                val = calc_rat_bin(
                    sig_A=sig_A,
                    sig_B=0.0,
                    scl_A=scl_A,
                    scl_B=scl_B,
                    psc_A=psc_A,
                    psc_B=psc_B,
                    dep_min=dep_min,
                    log2=log2,
                    recip=recip,
                    skip_00=skip_00,
                    eps=eps
                )
                write_line(chr_A, start_A, end_A, val)
                lin_A = read_data_line(fh=opn_A, skp_pfx=skp_pfx)

            else:
                #  Only B has bins left; A is exhausted, so treat A as 0.0 for
                #  remaining bins
                val = calc_rat_bin(
                    sig_A=0.0,
                    sig_B=sig_B,
                    scl_A=scl_A,
                    scl_B=scl_B,
                    psc_A=psc_A,
                    psc_B=psc_B,
                    dep_min=dep_min,
                    log2=log2,
                    recip=recip,
                    skip_00=skip_00,
                    eps=eps
                )
                write_line(chr_B, start_B, end_B, val)
                lin_B = read_data_line(fh=opn_B, skp_pfx=skp_pfx)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = CapArgumentParser(
        description=(
            "Compute per-bin ratios between two bedGraphs [A (numerator) and "
            "B (denominator)] with optional per-file multiplicative scaling, "
            "pseudocount addition, 'dep_min' “clamping” (denominator "
            "thresholding, i.e., a divisor “floor”), log2 transformation, "
            "and reciprocal computation.\n"
            "\n"
            "Can assign an error tolerance (ε [an “epsilon value”]) for "
            "treating values (e.g., “float noise”) as zero. deepTools "
            "bamCompare-like behavior is 'ε = 0.0' (the default).\n"
            "\n"
            "Order of operations is generally deepTools bamCompare-like:\n"
            "    1. Optionally skip zero-zero bins: '0 / 0' or 'ε / ε' "
            "(deepTools-like). ‡\n"
            "    2. Optionally scale each file.\n"
            "    3. Optionally skip scaled zero-zero bins: '[(sf_A × 0) / "
            "(sf_B × 0)]' or '[(sf_A × ε) / (sf_B × ε)]'. ‡\n"
            "    4. Optionally add pseudocounts.\n"
            "    5. Optionally clamp divisor (denominator) by 'dep_min' (see "
            "below).\n"
            "    6. If the divisor (denominator) is <= ε, treat the bin as "
            "undefined.\n"
            "    7. Divide 'A / B'.\n"
            "    8. Optionally perform log2 transformation.\n"
            "    9. Optionally compute reciprocal of #7 (linear) or #8 "
            "(log2).\n"
            "\n"
            "‡ Either of optional #1 or optional #3 can be applied, not "
            "both.\n"
            "\n"
            "(See module docstring for more details.)"
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
        "-ia", "-iA", "-fa", "-fA",
        "--infile_a", "--infile_A", "--fil_a", "--fil_A",
        dest="fil_A",
        required=True,
        type=str,
        help=(
            "Path to the first bedGraph input file (file 'A'; e.g., IP), or "
            "'-' to read from stdin. Supports plain text and '.gz'.\n"
            "\n"
            "Notes:\n"
            "  - At most one of '--fil_A' or '--fil_B' may be '-', since "
            "stdin can be consumed only once.\n"
            "  - When '-' is used for either input, the pre-merge bin-size "
            "consistency check is skipped (see 'I/O handling' in the module "
            "docstring).\n"
            "  - stdin/stdout workflows using '-' are supported when calling "
            "'compute_signal_ratio.py' directly; the higher-level Shell "
            "wrappers do not support '-' input/output.\n\n"
        )
    )
    parser.add_argument(
        "-ib", "-iB", "-fb", "-fB",
        "--infile_b", "--infile_B", "--fil_b", "--fil_B",
        dest="fil_B",
        required=True,
        type=str,
        help=(
            "Path to the second bedGraph input file (file 'B'; e.g., input), "
            "or '-' to read from stdin. Supports plain text and '.gz'.\n"
            "\n"
            "See '--fil_A' for stdin notes.\n\n"
        )
    )
    parser.add_argument(
        "-o", "-fo", "--outfile", "--fil_out",
        dest="fil_out",
        required=True,
        type=str,
        help=(
            "Path to the bedGraph outfile, or '-' to write to stdout. "
            "When a real path is provided, accepted extensions are "
            "'.bedGraph', '.bedgraph', '.bdg', and '.bg', each optionally "
            "followed by '.gz'.\n\n"
        )
    )
    parser.add_argument(
        "-me", "--method",
        choices=METHOD_CHOICES,
        default="unadj",
        help=(
            "Ratio-computation subtype (default: '%(default)s').\n"
            "  - Unadjusted aliases: 'r', 'raw', 'u', 'unadj', 'unadjusted', "
            "'s', 'smp', 'simple'. Internally standardized to 'unadj'.\n"
            "  - Log2 aliases: '2', 'l2', 'lg2', 'log2'. Internally "
            "standardized to 'log2'.\n"
            "  - Reciprocal-unadjusted aliases: 'rr', 'raw_r', 'ur', "
            "'unadj_r', 'unadjusted_r', 'sr', 'smp_r', 'simple_r'. Internally "
            "standardized to 'unadj_r'.\n"
            "  - Reciprocal-log2 aliases: '2r', 'l2r', 'l2_r', 'lg2_r', "
            "'log2_r'. Internally standardized to 'log2_r'.\n\n"
        )
    )
    parser.add_argument(
        "-sf", "--scl_fct",
        type=str,
        default=None,
        help=(
            "Per-file scale factor(s) 'A[:B]'. If only 'A' is given, 'B' "
            "defaults to 1.0.\n\n"
        )
    )
    parser.add_argument(
        "-ps", "--pseudo",
        type=str,
        default="0:0",
        help=(
            "Per-file pseudocount(s) 'A[:B]' added after scaling. (default: "
            "%(default)s)\n"
            "\n"
            "Note: This is primarily useful for log2-ratio methods, where it "
            "helps avoid undefined values such as 'log2(A / 0)'; for linear "
            "ratios, '--dep_min' may be preferable when the main concern is "
            "bounding low-depth extremes.\n"
            "\n"
            "Using '--pseudo' together with '--dep_min' is allowed, but "
            "usually makes low-depth stabilization harder to interpret.\n\n"
        )
    )
    parser.add_argument(
        "-dm", "--dep_min",
        type=float,
        default=None,
        help=(
            "Minimum allowed denominator threshold or “clamp” (i.e., “minimum "
            "input depth” as described in PMID 40364978). This is a “floor” "
            "(threshold) ensuring denominators do not fall below this value.\n"
            "\n"
            "'dep_min' is applied after any optional scaling and/or "
            "pseudocount addition ('B := max(B, dep_min)').\n"
            "\n"
            "With or without clamping, a zero-guard with ε ('--eps') is "
            "applied after this optional step in the order of operations; if "
            "'--dep_min <flt>' is specified, choose 'dep_min > ε'.\n"
            "\n"
            "Using '--dep_min' together with '--pseudo' is allowed, but "
            "usually makes low-depth stabilization harder to interpret.\n\n"
        )
    )
    parser.add_argument(
        "-e", "--eps",
        type=float,
        default=0.0,
        help=(
            "Zero tolerance ε (epsilon) for zero checks: Values satisfying "
            "'|value| <= ε' are treated as zero (default: %(default)s).\n"
            "\n"
            "Applied in two places:\n"
            "  - '--skip_00':\n"
            "    + If '--skip_00 pre_scale', check raw values (before any "
            "optional scaling and/or pseudocount addition).\n"
            "    + If '--skip_00 post_scale', check scaled values (before any "
            "optional pseudocount addition).\n"
            "  - Denominator guard:\n"
            "    + If '|B| <= ε', then the bin is emitted as 'nan'.\n"
            "    + If optional denominator clamping is specified, the "
            "denominator guard takes place after that.\n"
            "\n"
            "Notes:\n"
            "  - With 'ε = 0.0', '--skip_00 pre_scale' and '--skip_00 "
            "post_scale' are equivalent.\n"
            "  - With 'ε > 0' and 'scaling != 1', '--skip_00 post_scale' "
            "interprets ε in scaled units (i.e., the same units that are "
            "divided).\n\n"
        )
    )
    parser.add_argument(
        "-s0", "--skip_00",
        choices=["pre_scale", "post_scale"],
        default=None,
        help=(
            "Optionally drop bins where both values are treated as zero with "
            "ε ('--eps').\n"
            "  - '--skip_00 pre_scale':  Test on raw values (deepTools "
            "bamCompare-like behavior).\n"
            "  - '--skip_00 post_scale': Test after scaling.\n"
            "  - Omit entirely to disable the '0 / 0' (or 'ε / ε') drop.\n"
            "\n"
            "Notes:\n"
            "  - If not omitted, then 'pre_scale' and 'post_scale' modes are "
            "equivalent if 'ε = 0' or if no scaling is applied (i.e., 'scl_A "
            "= scl_B = 1').\n"
            "  - Pseudocounts do not affect this check in either mode.\n\n"
        )
    )
    parser.add_argument(
        "-dn", "--drp_nan",
        action="store_true",
        default=False,
        help=(
            "Skip non-finite 'nan', 'inf', and '-inf' (non-finite) values in "
            "the main output.\n"
            "\n"
            "Note: '--track' output will still drop 'inf', '-inf', and "
            "'nan'.\n\n"
        )
    )
    parser.add_argument(
        "-tr", "--track",
        action="store_true",
        default=False,
        help=(
            "Generate an additional bedGraph file where rows with 'inf', "
            "'-inf', and 'nan' values are excluded. The new file will have "
            "'.track' before the extension.\n"
            "\n"
            "Notes:\n"
            "  - This file is safe for use in genome browsers such as IGV.\n"
            "  - '--track' is not available when '--fil_out -' (stdout) is "
            "used.\n\n"
        )
    )
    parser.add_argument(
        "-dp", "--dp", "--rnd", "--round", "--decimals", "--digits",
        dest="rnd",
        type=int,
        default=24,
        help=(
            "Maximum number of decimal places retained for finite emitted "
            "ratio values (default: %(default)s). After rounding, "
            "non-informative trailing zeros are stripped.\n\n"
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

    #  If no arguments are provided, display help and exit
    argv_parse = sys.argv[1:] if argv is None else argv
    if not argv_parse:
        parser.print_help(sys.stderr)
        sys.exit(0)

    return parser.parse_args(argv_parse)


def main(argv: list[str] | None = None) -> int:
    """
    Execute the primary control flow for the script.

    Args:
        argv : list[str] | None
            Optional list of command-line arguments. When None (the default),
            'sys.argv[1:]' is used (the typical CLI entry point).

    Returns:
        int.
            On success, returns 0 and writes the main bedGraph, along with an
            optional '.track' sidecar when requested.

    Side effects:
        - When '--verbose' is set, prints a human-readable banner of the
          resolved arguments to stderr.
        - Prints error messages to stderr on validation/I/O failures.
        - Exits the interpreter (sys.exit) with status codes below.

    Exits:
        - 0 on success or when showing help with no arguments.
        - 1 on validation or computation errors, including (non-exhaustive):
            + missing inputs or unwritable output directory,
            + malformed bedGraph rows or inconsistent bin sizes,
            + invalid numeric arguments (e.g., 'scl_fct <= 0', 'dep_min <= 0',
              'eps < 0', 'rnd < 0'), and
            + read/write I/O errors.
    """
    #  Parse CLI arguments
    args = parse_args(argv)

    #  At most one input may be read from stdin ('-'; stdin can’t be consumed
    #  twice)
    n_stdin = sum(1 for f in (args.fil_A, args.fil_B) if f == "-")
    if n_stdin > 1:
        raise SystemExit(
            "At most one of '--fil_A' or '--fil_B' may be '-', since stdin "
            "can only be consumed once."
        )

    #  Warn when stdin is used: bin-size consistency check is skipped (as stdin
    #  cannot be rewound without buffering)
    if n_stdin == 1 and args.verbose:
        print(
            "Warning: Bin-size consistency check is skipped when either "
            "'--fil_A' or '--fil_B' is '-', since stdin cannot be safely "
            "read twice. Ensure that the A and B bedGraphs share a compatible "
            "binning grid (same chromosomes, bin starts, and widths).",
            file=sys.stderr
        )

    #  Validate infile existence
    if args.fil_A != "-":
        try:
            check_exists(args.fil_A, "file", "First file (A)")
        except FileNotFoundError as e:
            #  Print a one-line message and exit cleanly
            raise SystemExit(str(e))

    if args.fil_B != "-":
        try:
            check_exists(args.fil_B, "file", "Second file (B)")
        except FileNotFoundError as e:
            #  Print a one-line message and exit cleanly
            raise SystemExit(str(e))

    #  Check and standardize outfile (bedGraph only in this script); validate
    #  writability too; when '--fil_out -', write plain bedGraph to stdout
    if args.fil_out == "-":
        outfile = "-"
    else:
        try:
            outfile, _, _ = check_parse_outfile(args.fil_out, ONLY_BDG)
            check_writable(outfile, kind="file")
        except (
            ValueError, FileNotFoundError, PermissionError, IsADirectoryError
        ) as e:
            raise SystemExit(str(e))

    #  Disallow '--track' when writing main output to stdout
    if outfile == "-" and args.track:
        raise SystemExit(
            "Cannot use '--track' when '--fil_out -', as the '.track' sidecar "
            "requires a real output path."
        )

    #  Parse header-skip prefixes from CLI (empty string: no skipping)
    skp_pfx = parse_skp_pfx(args.skp_pfx, default=DEF_SKP_PFX)

    #  Parse per-file scale factors and pseudocounts
    try:
        if args.scl_fct is not None:
            scl_A, scl_B = parse_pair(args.scl_fct, 1.0)
        else:
            #  If not specified, default to '1.0:1.0'
            scl_A, scl_B = 1.0, 1.0

        psc_A, psc_B = parse_pair(args.pseudo, 0.0)
    except argparse.ArgumentTypeError as e:
        #  Surface malformed "A" or "A:B" string as a one-line error
        raise SystemExit(str(e))

    #  Validate 'skip_00'
    if args.skip_00 not in {None, "pre_scale", "post_scale"}:
        raise SystemExit(
            "Invalid '--skip_00'; expected None, 'pre_scale', or 'post_scale'."
        )

    #  Validate numeric arguments
    try:
        for lbl, v in (("scl_fct:A", scl_A), ("scl_fct:B", scl_B)):
            check_cmp(v, "gt", 0.0, lbl, allow_none=False)              # > 0
        for lbl, v in (("pseudo:A", psc_A), ("pseudo:B", psc_B)):
            check_cmp(v, "ge", 0.0, lbl, allow_none=False)              # >= 0
        check_cmp(args.dep_min, "gt", 0.0, "dep_min", allow_none=True)  # > 0
        check_cmp(args.eps, "ge", 0.0, "eps", allow_none=False)         # >= 0
        check_cmp(args.rnd, "ge", 0, "rnd", allow_none=False)           # >= 0

    except ValueError as e:
        raise SystemExit(str(e))

    #  Standardize '--method' to canonical internal name, then derive the
    #  internal computation flags from it
    mthd_in = args.method
    args.method = METHOD_CANON[args.method]

    log2 = args.method in {"log2", "log2_r"}
    recip = args.method in {"unadj_r", "log2_r"}

    #  Warn when the clamp threshold is at/below the zero tolerance
    if (
        args.verbose and args.dep_min is not None and
        args.dep_min <= args.eps
    ):
        #  If 'dep_min <= ε', the post-clamp denominator guard ('|B| <= ε')
        #  will still emit NaNs, effectively negating the clamp
        print(
            f"Warning: 'dep_min' ({args.dep_min}) <= 'ε' ({args.eps}), so "
            "denominator guard will emit 'nan' whenever '|B| <= ε'. Consider "
            "choosing 'dep_min > ε'.",
            file=sys.stderr
        )

    #  Warn when both stabilization strategies are supplied together
    if (
        args.verbose and args.dep_min is not None
        and (psc_A > 0.0 or psc_B > 0.0)
    ):
        print(
            "Warning: Both '--dep_min' and '--pseudo' were supplied. This is "
            "allowed, but interpretability may be reduced because both "
            "arguments stabilize low-depth ratio behavior in different ways.",
            file=sys.stderr
        )

    #  Warn when zero-zero skipping is combined with non-zero pseudocounts
    #  and a positive epsilon value
    if (
        args.verbose
        and args.skip_00 is not None
        and (psc_A > 0.0 or psc_B > 0.0)
        and args.eps > 0.0
    ):
        print(
            "Warning: '--skip_00', '--eps', and non-zero '--pseudo' were all "
            "supplied. This is allowed, but note that the zero-zero skip is "
            "applied before pseudocount addition, so pseudocounts do not "
            "rescue bins dropped by '--skip_00'.",
            file=sys.stderr
        )

    #  Print verbose output
    if args.verbose:
        with redirect_stdout(sys.stderr):
            print("##########################################")
            print("## Arguments for 'compute_signal_ratio' ##")
            print("##########################################")
            print("")
            print(f"--fil_A    {args.fil_A}")
            print(f"--fil_B    {args.fil_B}")
            print(f"--fil_out  {outfile}")
            if mthd_in != args.method:
                print(
                    f"--method   {mthd_in}  (standardized internally to "
                    f"{args.method})"
                )
            else:
                print(f"--method   {args.method}")
            if args.scl_fct is not None:
                print(f"--scl_fct {scl_A}:{scl_B}")
            print(f"--pseudo   {psc_A}:{psc_B}")
            if args.dep_min is not None:
                print(f"--dep_min  {args.dep_min}")
            if args.skip_00:
                print(f"--skip_00  {args.skip_00}")
            print(f"--eps      {args.eps}")
            if args.drp_nan:
                print("--drp_nan")
            print(f"--track    {args.track}")
            print(f"--rnd      {args.rnd}")
            print(f"--skp_pfx  {skp_pfx}")
            print("")
            print("")

    #  Call function for bin-wise ratio computations
    try:
        comp_sig_rat(
            fil_A=args.fil_A,
            fil_B=args.fil_B,
            fil_out=outfile,
            scl_A=scl_A,
            scl_B=scl_B,
            psc_A=psc_A,
            psc_B=psc_B,
            dep_min=args.dep_min,
            rnd=args.rnd,
            log2=log2,
            recip=recip,
            skip_00=args.skip_00,
            eps=args.eps,
            track=args.track,
            drp_nan=args.drp_nan,
            skp_pfx=skp_pfx
        )

    except (ValueError, FileNotFoundError, PermissionError, OSError) as e:
        raise SystemExit(str(e))

    return 0


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
