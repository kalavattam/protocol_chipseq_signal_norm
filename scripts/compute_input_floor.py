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
compute_input_floor.py


Description
-----------
Computes the minimum input depth ('dep_min') used as a “denominator floor” (or
“clamp”) when dividing IP by input in ratio-based ChIP-seq signal calculations.
The floor is intended to stabilize downstream divisions ('IP ÷ input') by
preventing extreme or erroneous ratios in low-depth bins.

The 'dep_min' computation depends on one or more of the following values:
    - Distribution-based ('mode=dist') quantities (bedGraph-like input):
        + v_i
            Per-bin signal values from a bedGraph-like file (column 4). Only
            finite values are used. For distribution-based methods, statistics
            are computed on a filtered subset of the values (see below).

        + eps, mode_nz
            Epsilon-based “zero handling” used to define which bins are treated
            as “nonzero” before applying the distribution rule.
                - 'mode_nz=closed'  drop values with '|v_i| <= eps'
                - 'mode_nz=open'    drop values with '|v_i| <  eps'
                - 'mode_nz=off'     disable epsilon-based filtering

            In this script, distribution-based floors are further restricted to
            positive values only ('v_i > 0'), keeping the floor in the same
            conceptual space as denominators in IP ÷ input ratio tracks.

        + method
            Distribution rule used to convert the filtered values into a
            single recommended denominator floor (see “Three computation
            modes” below).

        + qntl_nz
            Quantile in percent (0..100) used only when 'method=qntl_nz'.

        + coef
            Multiplicative coefficient used only for the
            “fraction-of-statistic” methods ('frc_mdn_nz', 'frc_avg_nz') and
            for 'min_nz'; ignored for 'qntl_nz'. If omitted, defaults match
            'compute_pseudo.py' (0.01 for 'frc_*'; 1.0 for 'min_nz').

        + floor
            A nonnegative lower bound applied after computing the method’s raw
            value:

                dep_min := max(dep_min, floor)

    - Fragment ('mode=frag') and normalization ('mode=norm')-based quantities:
        + n
            The number of counted alignments (reads or, preferably, fragments
            inferred from reads), as counted from an alignment-record file
            (bam or bed/bed.gz). In this script, n is not inferred from
            bedGraph-like inputs.

        + b
            The output bin size in bp ('--siz_bin').

        + g
            The effective genome size in bp ('--siz_gen').

Three computation modes are supported:
    - mode=dist
        + Compute a distribution-based denominator floor from a bedGraph-like
          input (column 4).

          Parsing and filtering pipeline for distribution-based methods:
            (1) Parse bedGraph column 4 values and keep only finite values.
            (2) Apply epsilon/zero handling ('--eps', '--mode_nz') to define
                which bins are treated as “nonzero”.
            (3) Restrict to positive values only ('v_i > 0'). (In practice,
                '--eps' influences which small positive values survive step
                (2): with '--mode_nz closed', values with 'v_i <= eps' are
                removed; with '--mode_nz open', values with 'v_i < eps' are
                removed.)

            - 'method=qntl_nz':

                dep_min = Q_q({v_i : v_i > 0}),

                where q = 'qntl_nz / 100' and 'qntl_nz' is provided by
                '--qntl_nz' in [0, 100]. Quantile selection uses a floor-based
                nearest-rank rule on sorted values:

                    i = floor(q * (N - 1)) (clamped to [0, N - 1])
                    dep_min = sorted_vals[i]

            - 'method=frc_mdn_nz':

                dep_min = coef × median({v_i : v_i > 0})

            - 'method=frc_avg_nz':

                dep_min = coef × mean({v_i : v_i > 0})

            - 'method=min_nz':

                dep_min = coef × min({v_i : v_i > 0})


        + After the method-specific value is computed ('mode=dist' only):

            dep_min := max(dep_min, floor)

    - mode=frag

        dep_min = [(n * b) / g] / [1 - (b / g)]

    - mode=norm

        dep_min = (b / g) / [1 - (b / g)]

Mode descriptions:
    - 'mode=dist'
        Computes a distribution-based denominator floor from a bedGraph-like
        input track (column 4). Values are filtered to finite values, then
        filtered by '--eps' / '--mode_nz', then restricted to positive-only
        values ('v_i > 0'). The resulting set is summarized into a single
        recommendation using '--method' ('qntl_nz', 'frc_mdn_nz', 'frc_avg_nz',
        'min_nz').

        If '--floor' is provided (>= 0), it is applied after the
        method-specific value is computed:

            dep_min := max(dep_min, floor)

        Note: '--method', '--qntl_nz', '--coef', '--floor', '--eps', and
        '--mode_nz' are only used in 'mode=dist'.

    - 'mode=frag'
        Computes 'dep_min' from the alignment-record count 'n' and sizes 'b'
        and 'g' according to the siQ-ChIP-style fragment-normalization
        derivation:

            dep_min = [(n * b) / g] / [1 - (b / g)]

        Requires an alignment-record input (bam or bed/bed.gz).

    - 'mode=norm'
        Computes 'dep_min' for “normalized coverage” (siQ-ChIP-style), which
        depends only on 'b' and 'g':

            dep_min = (b / g) / [1 - (b / g)]

        Here, normalized coverage means the corresponding genome-wide signal
        is scaled so its integral (i.e., sum over all bins) is 1. No input file
        is needed for this calculation (see algebraic simplification below);
        '--infile' is ignored.

Input types:
    - 'mode=dist'
        bedGraph-like (bedGraph, bedgraph, bdg, or bg; optionally, with .gz).

    - 'mode=frag'
        bam or bed/bed.gz (records counted as alignments/fragments).

    - 'mode=norm'
        No infile needed (infile ignored).

Note that bam counting uses the following default FLAG sets:
    - paired-end  '{99, 1123, 163, 1187}'
    - single-end  '{0, 16, 1024, 1040}'
These reflect common “main” alignments while tolerating duplicate flags; they
are configurable via optional arguments '--flags-pe' and '--flags-se'.


Usage
-----
python -m scripts.compute_input_floor \\
    [--help] [--verbose] \\
    --infile <str_fil.(bam|bed(.gz)?|(bedGraph|bedgraph|bdg|bg)(.gz)?)> \\
    --siz_bin <int_bp> --siz_gen <int_bp> \\
    --mode {dist,frag,norm} [--dp <int_decimals>] [--skp_pfx <str,...>]


Output
------
Writes a single floating-point number (rounded to '--dp' decimals) to stdout.
Non-zero exit on error (missing file, unsupported extension, bad parameters).


Examples
--------
1. bam input, 10-bp bins, S. cerevisiae effective genome size, “norm” factor
'''bash
python -m scripts.compute_input_floor \
    -i input.bam \
    -sb 10 \
    -sg 12157105 \
    -md norm
'''

2. bed.gz input, 50-bp bins, custom genome size, “frag” factor to eight decimal
   places
'''bash
python -m scripts.compute_input_floor \
    -i input.bed.gz \
    -sb 50 \
    -sg 2864785220 \
    -md frag \
    -dp 8
'''

3. bam with custom FLAG allow-lists (hex OK)
'''bash
python -m scripts.compute_input_floor \
    -i in.bam \
    -sb 10 \
    -sg 12157105 \
    -md frag \
    -fp 99,0x463,163,0x4A3 \
    -fs 0,16,0x400,0x410
'''

4. Compute 'dep_min' once from an input library and pass it to a signal ratio
   computation.
'''bash
dep_min=$(
    python -m scripts.compute_input_floor \
        -i input.bam \
        -sb 10 \
        -sg 12157105 \
        -md norm \
        -dp 12
)

#  Generate a log2(IP/input) track while clamping 'fB' with 'dep_min'
python -m scripts.compute_input_floor \
    -fA IP.sc.bdg.gz \
    -fB in.sc.bdg.gz \
    -fo log2.sc.bdg.gz \
    -l2 \
    -dm "${dep_min}"
'''

5. Pipe a bam file (S288C) to 'compute_input_floor.py'
'''bash
samtools view -b in.bam \
    | python -m scripts.compute_input_floor \
        -i - --infmt bam -sb 10 -sg 12157105 -md frag
'''

6. Pipe a bed.gz (GRCh38) file to 'compute_input_floor.py'
'''bash
zcat in.bed.gz \
    | python -m scripts.compute_input_floor \
        -i - --infmt bed -sb 50 -sg 2913022398 -md frag
'''

7. Compute a distribution-based floor: quantile; 1st percentile of positive
   bins.
'''bash
python -m scripts.compute_input_floor \
    --infile in.sc.bdg.gz \
    --mode dist \
    --method qntl_nz \
    --qntl_nz 1
'''

8. Compute a distribution-based floor: fraction of median; default 'coef=0.01'.
'''bash
python -m scripts.compute_input_floor \
    --infile in.sc.bdg.gz \
    --mode dist \
    --method frc_mdn_nz
'''

9. Compute a distribution-based floor: fraction of mean; explicit coef; treat
   tiny float noise as zero via 'eps'.
'''bash
python -m scripts.compute_input_floor \
    --infile in.sc.bdg.gz \
    --mode dist \
    --method frc_avg_nz \
    --coef 0.01 \
    --eps 1e-12
'''

10. Compute a distribution-based floor (fraction of minimum positive bin; clamp
    with '--floor' to avoid a 0-like floor.
'''bash
python -m scripts.compute_input_floor \
    --infile in.sc.bdg.gz \
    --mode dist \
    --method min_nz \
    --coef 0.5 \
    --floor 1e-12
'''


General notes
-------------
- Effective genome size ('g') should match the mappable/usable size for the
  library and reference build used downstream.
- Choose bin size ('b') to match the binning used to compute signal. [TODO: Is
  this checked for / warned against / errored?]
- The 'mode=norm' mode mirrors unity-normalized coverage workflows where total
  genome-wide signal integrates to 1; 'mode=frag' matches fragment-length
  normalization sans the unity constraint.
- In 'mode=norm', the algebra simplifies:
      dep_min = ([(n * b) / g] / [1 - (b / g)]) / n = (b / g) / [1 - (b / g)].
    + As '(b / g)' is expected to be a very small number for typical bin sizes
      in ChIP-seq analyses, '(b / g) / [1 - (b / g)] ≈ (b / g)'.
    + Since 'mode=norm' does not depend on the read/fragment count 'n', we do
      not tally alignments in that mode.


Performance notes
-----------------
- bam counting streams via 'pysam' and scales linearly with the number of
  alignments; bed/bed.gz counting streams line-by-line.
- I/O tends to dominate runtime in 'mode=frag' mode; colocating input on fast
  storage can help with this.
- In 'mode=norm', no alignment counting or file I/O is performed; the result
  depends only on 'b' and 'g', so runtime is much quicker.


#TODO
-----
- Unit tests:
    + a small test that feeds a short bed with a track line and blanks;
    + a test that asserts the 'siz_bin >= siz_gen' check; and
    + a flag-parsing test with hex tokens.
    + Maybe/probably others.
- Extend 'compute_input_floor()' for deepTools normalizations, including
  (non-exhaustive) CPM, {R,F}PKM, {T,B}PM, and {R,F}PGC (and perhaps other non-
  deepTools normalizations).
- Add optional stdin and auto-compressed input support for bed/bed.gz files via
  a small helper [e.g., 'open_text_auto(path: str)' that accepts "-" and
  chooses 'gzip'/plain automatically], and wire it into 'count_aln_bed()'.
- Perhaps something similar to the above for bam files?
"""

from __future__ import annotations

from contextlib import redirect_stdout

import argparse
import math
# import pysam  # Lazy-imported in 'count_aln_bam()'
import signal
import sys

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
    is_header,
    open_in,
    parse_skp_pfx,
)
from scripts.functions.utils_stabilizer import (
    determine_coef_eff,
    iter_vals_bdg,
    pick_stabilizer,
)

try:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except (AttributeError, ValueError):
    pass

assert sys.version_info >= (3, 10), "Python >= 3.10 required."

#  Define global default SAM FLAG values
FLG_PE = {99, 1123, 163, 1187}  # Proper paired-end alignments
FLG_SE = {0, 16, 1024, 1040}    # Single-end alignments

#  Run script in “interactive mode” (true) or “CLI mode” (false)
interactive = False


def set_interactive(
    echo: bool = False, ensure: bool = False
) -> argparse.Namespace:
    """
    Return a pre-filled argument namespace for interactive runs.

    Args:
        echo : bool
            If True, echo chosen parameters to stderr.
        ensure : bool
            If True, perform light validation of the interactive defaults
            before returning.

    Returns:
        argparse.Namespace
            Namespace mirroring the CLI options.
    """

    #  Set general paths
    dir_bas = "/home/kalavatt/tsukiyamalab/Kris"
    dir_rep = f"{dir_bas}/protocol_chipseq_signal_norm"
    dir_dat = f"{dir_rep}/data"
    dir_pro = f"{dir_dat}/processed"
    dir_exp = f"{dir_pro}/2025-08-19"
    dir_aln = f"{dir_exp}/align_fastqs"

    aln = "bowtie2"
    atp = "global"
    flg = 2
    mpq = 1
    det = f"{aln}_{atp}_flag-{flg}_mapq-{mpq}"
    spc = "sc"

    dir_bam = f"{dir_aln}/{det}/{spc}"
    fil_bam = "in_WT_H3_Q_5781.sc.bam"
    pth_bam = f"{dir_bam}/{fil_bam}"

    #  Set argument values
    verbose = True
    infile = pth_bam
    infmt = None
    siz_bin = 10
    siz_gen = 12157105  # Match 'parse_args' default
    flags_pe = None
    flags_se = None
    mode = "frag"       # "norm" or "dist"

    #  Distribution-mode knobs (used only when 'mode=dist')
    method = "qntl_nz"
    qntl_nz = 1.0
    coef = None
    floor = 0.0
    eps = 0.0
    mode_nz = "closed"

    dp = 24            # Match 'parse_args' default
    skp_pfx = None

    #  Wrap the arguments in argparse.Namespace
    ns = argparse.Namespace(
        verbose=verbose,
        infile=infile,
        infmt=infmt,
        siz_bin=siz_bin,
        siz_gen=siz_gen,
        flags_pe=flags_pe,
        flags_se=flags_se,
        mode=mode,

        #  'mode=dist' args
        method=method,
        qntl_nz=qntl_nz,
        coef=coef,
        floor=floor,
        eps=eps,
        mode_nz=mode_nz,

        dp=dp,
        skp_pfx=skp_pfx
    )

    #  “Echo” the interactive argument assignments (etc.) if specified
    if echo:
        echo_block("paths, etc.", {
            "dir_bas": dir_bas, "dir_rep": dir_rep, "dir_dat": dir_dat,
            "dir_pro": dir_pro, "dir_exp": dir_exp, "dir_aln": dir_aln,
            "dir_bam": dir_bam, "fil_bam": fil_bam, "pth_bam": pth_bam
        })
        echo_block("args", vars(ns))

    #  Return the argparse.Namespace-wrapped arguments
    return ns


def note_ignored(opt: str, reason: str, val) -> None:
    """
    Print a standard “ignored option” note when a mode does not use an arg.

    Args:
        opt : str
            Option name as shown to the user (e.g., '--infile').
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

    Args:
        csv_str : str
            Comma-separated list of FLAG values to allow. Tokens may be decimal
            (e.g., '99') or hexadecimal with a '0x' prefix (e.g., '0x63').
            Whitespace around commas is ignored, and empty tokens are skipped.
        label : str
            Short name used in error messages to identify which option the
            flags came from (e.g., 'flags-pe' or 'flags-se').

    Returns:
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

    Args:
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

    Args:
        pth_bed : str
            Path to bed file (can be gzipped).

    Returns:
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

    Args:
        p : str
            Path to the input file, or "-" for stdin.

        hint : str | None = None
            Optional format hint used only when 'p == "-"'. Must be one of
            {"bam", "bed", "bedGraph", "bedgraph", "bdg", "bg"} if provided.

    Returns:
        str
            - "bam" for bam.
            - "bed" for bed (optionally with .gz).
            - "bedgraph" for bedGraph, bedgraph, bdg, or bg (optionally with
              .gz).
            - "other" for otherwise.
            - When 'p == "-"' and a valid 'hint' is given, the hint value is
              returned.

    Raises:
        None.
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
    infile: str,
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

    Args:
        infile : str
            Path to input file:
                - mode=dist  bedGraph-like (bedGraph, bedgraph, bdg, or bg,
                             optionally with .gz)
                - mode=frag  bam or bed/bed.gz
                - mode=norm  infile is ignored

        siz_bin : int
            Output bin size in base pairs (i.e., the bedGraph bin width).

        siz_gen : int
            Effective genome size of the organism in base pairs.

        mode : str
            One of {'dist','frag','norm'}:

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
                    (siQ-ChIP-style), and depends only on 'b' and 'g'; infile
                    is ignored.

        method : str
            Distribution-based rule used only when mode='dist'. One of:
            {'qntl_nz','frc_mdn_nz','frc_avg_nz','min_nz'}.

        qntl_nz : float
            Quantile in percent (0..100) used only when mode='dist' and
            method='qntl_nz'. The quantile is computed on positive values only
            (v_i > 0) after eps/mode_nz filtering.

        coef : float | None
            Coefficient used only when mode='dist' and method is one of
            {'frc_mdn_nz','frc_avg_nz','min_nz'}. If None, defaults match
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
            Optional format hint forwarded to 'infer_fmt' when 'infile == "-"'.
            Accepts "bam", "bed", or "bedgraph" (and bedGraph-like aliases such
            as "bedGraph", "bdg", "bg"). Ignored otherwise.

    Returns:
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
            - unsupported infile extension (dist: bedGraph-like; frag: bam or
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
        fmt = infer_fmt(infile, infmt)
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
                infile,
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
    fmt = infer_fmt(infile, infmt)
    if fmt == "bam":
        n_in = count_aln_bam(infile, flg_pe=flg_pe, flg_se=flg_se)
    elif fmt == "bed":
        if (flg_pe or flg_se):
            print(
                "Note: '--flags-pe' / '--flags-se' are ignored for bed "
                "inputs.",
                file=sys.stderr
            )
        n_in = count_aln_bed(infile, skp_pfx)
    elif fmt == "bedgraph":
        print(
            "Error: '--mode frag' expects bam or bed/bed.gz (alignment "
            "records), not bedGraph.",
            file=sys.stderr
        )
        sys.exit(1)
    else:
        print(
            f"Error: Unsupported file type: {infile}. Provide a bam or "
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
        action="store_true",
        default=False,
        help="Increase output verbosity.\n\n"
    )
    parser.add_argument(
        "-i", "--infile",
        type=str,
        required=False,
        default=None,
        help=(
            "Path to input file, or use '-' for stdin.\n"
            "    - '--mode dist' requires a bedGraph-like file: bedGraph, "
            "bedgraph, bdg, or bg (optionally with .gz).\n"
            "    - '--mode frag' requires bam or bed/bed.gz.\n"
            "    - '--mode norm' infile is ignored.\n"
            "\n"
        )
    )
    parser.add_argument(
        "-if", "--infmt",
        choices=["bam", "bed", "bedGraph", "bedgraph", "bdg", "bg"],
        default=None,
        help=(
            "Input format hint. Required only when '--infile -' (stdin). "
            "Choose 'bam', 'bed', 'bedGraph', 'bedgraph', 'bdg', or 'bg'.\n"
            "\n"
        )
    )
    parser.add_argument(
        "-sb", "--siz_bin",
        type=int,
        default=10,
        help="bedGraph bin size in base pairs (default: %(default)s).\n\n"
    )
    parser.add_argument(
        "-sg", "--siz_gen",
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
        choices=["dist", "frag", "norm"],
        default="dist",
        help=(
            "Normalization mode (default: '%(default)s').\n"
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
            "'--infile' is ignored.\n"
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
        "--method",
        choices=("qntl_nz", "frc_mdn_nz", "frc_avg_nz", "min_nz"),
        default="qntl_nz",
        help=(
            "Distribution-based method used only in '--mode dist' "
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
        "--qntl_nz",
        type=float,
        default=1.0,
        help=(
            "Quantile in percent (0..100) used only when '--mode dist' AND "
            "'--method qntl_nz' (default: %(default)s).\n\n"
        )
    )
    parser.add_argument(
        "--coef",
        type=float,
        default=None,
        help=(
            "Coefficient for '--method frc_mdn_nz', '--method frc_avg_nz', "
            "and '--method min_nz'. If omitted, defaults match "
            "'compute_pseudo.py' (0.01 for 'frc_*'; 1.0 for 'min_nz').\n\n"
        )
    )
    parser.add_argument(
        "--floor",
        type=float,
        default=0.0,
        help=(
            "Lower bound applied to the computed floor in '--mode dist' "
            "(default: %(default)s).\n\n"
        )
    )
    parser.add_argument(
        "--eps",
        type=float,
        default=0.0,
        help=(
            "Zero tolerance epsilon used only in '--mode dist' (default: "
            "%(default)s).\n\n"
        )
    )
    parser.add_argument(
        "--mode_nz",
        choices=("closed", "open", "off"),
        default="closed",
        help=(
            "Epsilon/zero-handling mode used only in '--mode dist' "
            "(default: %(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-dp", "--dp", "--rnd", "--round", "--decimals", "--digits",
        type=int,
        default=24,
        help=(
            "Number of decimal places for rounding result (default: "
            "%(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-sp", "--skp_pfx",
        type=str,
        default=",".join(DEF_SKP_PFX),
        help=(
            "Comma-separated list of prefixes to skip as header/meta lines in "
            "bed/bed.gz or bedGraph-like inputs; to disable skipping, pass an "
            "empty string (default: %(default)s).\n\n"
        )
    )

    #  If no arguments are provided, use 'argv' to display help and exit
    argv_parse = sys.argv[1:] if argv is None else argv
    if not argv_parse:
        parser.print_help(sys.stderr)
        sys.exit(0)

    return parser.parse_args(argv_parse)


def main(argv: list[str] | None = None) -> int:
    """
    Execute the primary control flow for the script.

    Args:
        None.

    Returns:
        None. On success, prints the “floor” (“input depth factor”) to stdout.

    Raises:
        SystemExit
            On validation failures or when required dependencies are missing.

    Notes:
        - Emits a note to stderr if '--flags-pe' / '--flags-se' are supplied
          for bed/bed.gz input (as flags are ignored for bed inputs).
        - Prints human-readable error messages to stderr on failure.
        - BrokenPipeError is handled in the '__main__' wrapper.
    """
    #  Handle interactive or CLI arguments
    args, early_exit = get_args_interactive(
        argv, interactive, set_interactive, parse_args
    )
    if early_exit:
        return 0

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
            #  UX note: in 'norm' mode, infile/infmt/flags are ignored, as are
            #           'dist'-mode knobs
            if args.infile is not None:
                note_ignored("--infile", "in '--mode norm'", args.infile)
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
            #  In 'dist' and 'frag' modes, '--infile' is required
            if not args.infile:
                print(
                    "Error: '--infile' is required when '--mode dist' or "
                    "'--mode frag'.",
                    file=sys.stderr
                )
                raise SystemExit(1)

            fmt = infer_fmt(args.infile, args.infmt)

            if args.mode == "dist":
                #  Validate input type for 'mode=dist'
                if fmt != "bedgraph":
                    if args.infile == "-":
                        print(
                            "Error: When '--infile -' is used with '--mode "
                            "dist', provide '--infmt {bedgraph,bedGraph,bdg,"
                            "bg}'.",
                            file=sys.stderr
                        )
                    else:
                        print(
                            f"Error: Unsupported file type for '--mode dist': "
                            f"'{args.infile}'. Provide bedGraph-like input "
                            "(bedgraph, bedGraph, bdg, or bg, optionally with "
                            ".gz).",
                            file=sys.stderr
                        )
                    raise SystemExit(1)

                check_cmp(args.floor, "ge", 0.0, "floor", allow_none=False)
                check_cmp(args.eps,   "ge", 0.0, "eps",   allow_none=False)
                check_cmp(args.coef,  "ge", 0.0, "coef",  allow_none=True)

                if args.method == "qntl_nz":
                    if (
                        not math.isfinite(args.qntl_nz)
                    ) or not (
                        0.0 <= args.qntl_nz <= 100.0
                    ):
                        raise SystemExit(
                            "Error: '--qntl_nz' must be finite and in [0, "
                            "100]."
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
                    if args.infile == "-":
                        print(
                            "Error: When '--infile -' is used with "
                            "'--mode frag', provide '--infmt {bam,bed}'.",
                            file=sys.stderr
                        )
                    else:
                        print(
                            f"Error: Unsupported file type for '--mode frag': "
                            f"'{args.infile}'. Provide bam or bed/bed.gz.",
                            file=sys.stderr
                        )
                    raise SystemExit(1)

            if args.infile != "-":
                check_exists(args.infile, kind="file", label=fmt.upper())

    except ValueError as e:
        raise SystemExit(str(e))

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
            print(f"--infile {args.infile}")
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
            infile=args.infile if args.mode in {"dist", "frag"} else "-",
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
        raise SystemExit(1)

    #  Print result
    print(f"{dep_min:.{args.dp}f}")


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
