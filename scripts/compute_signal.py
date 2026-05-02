#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: compute_signal.py
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


"""
Script
------
compute_signal.py


Description
-----------
Calculates binned signal from a BAM file, typically outputting bedGraph
signal tracks. Users can compute raw (unadjusted) signal, signal normalized by
fragment length (as in PMID 37160995), or "normalized coverage" (i.e., signal
adjusted for both fragment length and total fragment count, as in PMID
37160995). An optional scaling factor can also be applied.

The output format is determined by the '--outfile' extension when '--outfile'
is a real path:
    - extensions '.bedGraph', '.bedGraph.gz', '.bedgraph', '.bedgraph.gz',
      '.bdg', '.bdg.gz', '.bg', and '.bg.gz' result in bedGraph-formatted
      output.
    - extensions '.bed' and '.bed.gz' output processed fragment coordinates,
      not signal, in BED-like format.
    - extension '.gz' results in additional gzip compression.

When '--outfile -' (stdout) is used, '--outfmt' is required to tell the script
which format to emit: 'bedGraph', 'bedgraph', 'bdg', 'bg', or 'bed'. (This
stdout workflow is intended for direct use of 'compute_signal.py'; the higher-
level Shell wrappers do not support stdout output.)

BED output is useful when generating fragment-coordinate input for the original
siQ-ChIP implementation: https://github.com/BradleyDickson/siQ-ChIP.

The script supports parallel processing, with per-chromosome signal
calculations distributed across multiple CPU threads.


Usage
-----
python -m scripts.compute_signal \
    [--help] [--verbose] \
    [--threads <int_no_cores>] \
    --infile <str_file.bam|-> \
    --outfile <str_file.(bedGraph|bedgraph|bdg|bg|bed)(.gz)|-> \
    [--outfmt {bedGraph,bedgraph,bdg,bg,bed}] \
    [--method <str_method>] \
    [--siz_bin <int_bin_size_bp>] \
    [--scl_fct <flt_scale_factor>] \
    [--usr_frg <int_frag_len_bp>] \
    [--rnd <int_decimals>]


Arguments
---------
 -v, --verbose
    Run script in verbose mode (stderr banner of parsed args).

 -t, --threads <int>
    Number of threads for parallel processing (default: 1).

 -i, --infile <str>
    Path to the input BAM file, or '-' to read from stdin.

 -o, --outfile <str>
    Path to the output file, or '-' to write to stdout. Supported formats are
    '.bedGraph', '.bedgraph', '.bdg', '.bg', and '.bed', each optionally with a
    '.gz' suffix for gzip compression. (stdout output is intended for direct
    use of 'compute_signal.py'; the higher-level Shell wrappers do not support
    this mode.)

-of, --outfmt {(bedGraph|bedgraph|bdg|bg),bed}
    Output format hint used only when '--outfile -' (stdout). Must be one of
    'bedGraph', 'bedgraph', 'bdg', 'bg', or 'bed'. Ignored when '--outfile' is
    a real path, in which case the format is inferred from the outfile
    extension. (This option is mainly for direct use of 'compute_signal.py',
    since the higher-level Shell wrappers do not support stdout output.)

-me, --method {...}
    Specify signal calculation type (default: 'norm'). Options:
        - 'r', 'raw', 'u', 'unadj', 'unadjusted', 's', 'smp', 'simple': Compute
          unadjusted signal. Internally, all of these values are standardized
          to 'method=unadj'.
        - 'f', 'frg', 'frag', 'frg_len', 'frag_len', 'l', 'len', 'len_frg',
          'len_frag': Normalize signal by fragment length (as described in PMID
          37160995). Internally, all of these values are standardized to
          'method=frag'.
        - 'n', 'nrm', 'norm', 'normalized': Normalize signal by fragment length
          and total fragment count so that the genome-wide summed signal is
          approximately 1 (up to small floating-point error), i.e., "normalized
          coverage" (as described in PMID 37160995). Internally, all of these
          values are standardized to 'method=norm'.
        - '--method' is ignored for BED fragment coordinate output.

-sb, --siz_bin <int>
    Bin size for signal calculation in base pairs (default: 10). Ignored for
    BED output.

-sf, --scl_fct <flt>
    Optional scaling factor to apply to the signal. Can be used with any
    '--method'. Ignored for BED output.

-uf, --usr_frg <int>
    Optional fixed fragment length to use instead of (a) BAM query alignment
    lengths for single-end alignments or (b) template lengths for paired-end
    alignments. Affects both bedGraph and BED outputs. For single-end data,
    '--usr_frg' is generally the preferred way to obtain fragment-length-aware
    signal.

-dp, --dp, --rnd, --round, --decimals, --digits <int>
    Maximum number of decimal places retained for bedGraph values (default:
    24). After rounding, non-informative trailing zeros are stripped. Ignored
    for BED output.


Output
------
- bedGraph ('.bedGraph', '.bedGraph.gz', '.bedgraph', '.bedgraph.gz', '.bdg',
  '.bdg.gz', '.bg', '.bg.gz'):
    + Contains binned signal values.
    + Uses '--rnd' as the maximum decimal precision for finite values.
    + After rounding, non-informative trailing zeros are stripped.
- BED ('.bed', '.bed.gz'):
    + Contains processed fragment coordinates instead of binned signal.
    + Ignores '--siz_bin', '--method', '--scl_fct', and '--rnd'.
    + Honors '--usr_frg' (affects interval lengths/positions).
    + Column 4 contains the fragment length used for that interval.
- Compression: If '.gz' is in the filename, output is automatically
  gzip-compressed.


Examples
--------
1. Compute normalized coverage at 10 bp bins, gzipping bedGraph output
'''bash
python -m scripts.compute_signal \
    -i sample.bam \
    -o sample.bdg.gz \
    -sb 10 \
    -me norm
'''

2. Compute fragment-length normalized signal, 25-bp bins, scaled by 1.732
'''bash
python -m scripts.compute_signal \
    -i sample.bam \
    -o sample.bedgraph \
    -sb 25 \
    -me frag \
    -sf 1.732
'''

3. Compute raw (unadjusted) signal, 50-bp bins, with fixed fragment length 150
'''bash
python -m scripts.compute_signal \
    -i sample.bam \
    -o sample.bg \
    -sb 50 \
    -me unadj \
    -uf 150
'''

4. Export processed fragment coordinates in BED format for use with siQ-ChIP
'''bash
python -m scripts.compute_signal -i sample.bam -o sample.bed.gz
'''

5. Stream BAM from samtools and write bedGraph of normalized coverage to
   stdout (then gzip; not supported by Shell wrappers)
'''bash
samtools view -b input.bam \
    | python -m scripts.compute_signal -i - -o - --outfmt bdg -sb 10 -me norm \
    | gzip > sample.norm.bdg.gz
'''

6. Write BED to stdout (not supported by Shell wrappers)
'''bash
python -m scripts.compute_signal -i sample.bam -o - --outfmt bed -uf 150 \
    > sample.bed
'''


General notes
-------------
- The input BAM file does not need to be coordinate-sorted or queryname-sorted,
  as read order is not relevant for the parsing implemented here.
    + No BAM index is required (streaming is handled with
      'pysam.fetch(until_eof=True)'); random access is not used.
    + '--infile -' is also supported and behaves the same way.
    + For paired-end data, fragments are inferred from leftmost mates with
      TLEN > 0; for single-end data, intervals are extended from the 5' ends
      of alignments.
    + Binning accumulates fragments into per-chromosome dict buckets, and
      sorting occurs at emit.
- Coordinates are 0-based, half-open [start, end) in both BED and bedGraph.
- Chromosome naming is preserved from the BAM.
- Roman numeral sorting is applied when emitting both BED and bedGraph (I..XVI
  first), reflecting that this code was written with S. cerevisiae and S. pombe
  ChIP-seq data in mind.
    + After Roman numerals, purely numeric names, then mitochondrial aliases,
      then lexical scaffolds are used; within each chromosome, records are
      sorted by start.
- For single-end read alignments, using '--usr_frg' yields aligner/read-length-
  agnostic coverage analogous to '--extendReads' in deepTools.
- bedGraph output rounds finite values to at most '--rnd' decimal places, then
  strips non-informative trailing zeros and any trailing decimal point.


Performance notes
-----------------
- When '--threads > 1', work is per-chromosome via 'ProcessPoolExecutor'.
- I/O dominates for large BAMs.
    + Placing BAM and output on fast local storage can help with this.
    + Too many threads can increase contention on slow disks.
- Memory use scales with the number of fragments retained per chromosome; for
  very dense chromosome data, chunked processing may be worth implementing.


#TODO
-----
- May implement chunked processing; it may be faster than the current
  per-chromosome approach and may also help reduce memory use for large,
  data-dense chromosomes.
- May reimplement the core algorithm in a compiled language.
- Currently, bins can extend beyond chromosome ends. This can cause problems in
  downstream analyses (e.g., conversion of bedGraph to bigWig with Kent tools)
  and should be addressed by trimming bins to chromosome boundaries.
- Test with model organisms with non-Roman-numeral chromosome names, both with
  and without 'chr' prefixes.
"""

from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout

from scripts.functions.utils_bdg import write_bdg
from scripts.functions.utils_check import (
    ALLOWED,
    check_cmp,
    check_exists,
    check_parse_outfile,
    check_writable
)
from scripts.functions.utils_chrom import sort_chrom
from scripts.functions.utils_cli import (
    add_help_cap,
    CapArgumentParser
)
from scripts.functions.utils_io import open_out

import argparse
import os
import pysam
import signal
import sys

try:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except (AttributeError, ValueError):
    pass

assert sys.version_info >= (3, 10), "Python >= 3.10 required."

#  Accepted '--method' values and their canonical internal names
METHOD_CANON = {
    #  Unadjusted signal
    "r": "unadj",
    "raw": "unadj",
    "u": "unadj",
    "unadj": "unadj",
    "unadjusted": "unadj",
    "s": "unadj",
    "smp": "unadj",
    "simple": "unadj",

    #  Fragment-length-normalized signal
    "f": "frag",
    "frg": "frag",
    "frag": "frag",
    "frg_len": "frag",
    "frag_len": "frag",
    "l": "frag",
    "len": "frag",
    "len_frg": "frag",
    "len_frag": "frag",

    #  Fragment-length-and-depth-normalized signal ("normalized coverage")
    "n": "norm",
    "nrm": "norm",
    "norm": "norm",
    "normalized": "norm",
}
METHOD_CHOICES = tuple(METHOD_CANON.keys())


def parse_bam(
    path_bam: str,
    usr_frg: int | None = None,
    allow_sec: bool = True,
    allow_sup: bool = False,
    allow_dup: bool = True
) -> dict[str, list[tuple[int, int, int]]]:
    """
    Parse a BAM into chromosome-grouped fragment intervals.

    Args:
        path_bam : str
            Input BAM path.
        usr_frg : int | None
            Optional fixed fragment length override. If provided:
                - paired-end read alignments: overrides TLEN for leftmost
                                              anchors.
                - single-end read alignments: used for both strands from the 5’
                                              end.
        allow_sec : bool
            Include secondary (multi-mapping) alignments. Default True.
        allow_sup : bool
            Include supplementary alignments. Default False.
        allow_dup : bool
            Include duplicate-marked alignments. Default True.

    Returns:
        frg_tup : dict[str, list[tuple[int, int, int]]]
            A dictionary mapping chromosomes to lists of (start, end, frg_len),
            where [start, end) is a half-open interval in 0-based coordinates
            and frg_len is the fragment length used for that interval.

    Raises:
        FileNotFoundError
            If 'path_bam' does not exist.
        ValueError
            If 'usr_frg' is provided but is <= 0 for paired-end or single-end
            extension.

    Notes:
        The output uses half-open intervals; downstream binning iterates
        positions in [start, end) (i.e., range(start, end)).

    Policies, semantics:
        - One fragment per proper paired-end alignment pair, anchored on the
          leftmost mate:
            + A read alignment is treated as a leftmost paired-end anchor iff
              all of the following are true:
              '''
                  read.is_paired
              and read.is_proper_pair
              and read.reference_id == read.next_reference_id
              and read.template_length > 0
              '''
            + For such reads, the fragment starts at 'read.reference_start',
              and the fragment length is TLEN unless overridden by 'usr_frg'.
            + In the current implementation, only leftmost proper-pair anchors
              with 'TLEN > 0' are treated as paired-end fragment anchors.

            NOTE: This intentionally differs from deepTools, where proper pairs
                  always use the observed TLEN (also, deepTools applies a
                  distance guard). Here, the user may override TLEN with
                  'usr_frg' by design (also, no distance guard is in place).

        - Single-end read alignments ('read.is_paired == False') are extended
          strand-aware from the 5’ end to length 'frg_len', where
          'frg_len = usr_frg' if provided, else 'read.query_alignment_length'.
            + Forward strand (not 'read.is_reverse'):
                  [start, end) = [reference_start, reference_start + frg_len)
            + Reverse strand ('read.is_reverse'):
                  [start, end) = [reference_end - frg_len, reference_end)

            NOTE: A fixed fragment length ('usr_frg') is preferred, mirroring
                  deepTools’ '--extendReads', as it produces aligner- and read-
                  length–independent coverage. However, we do not enforce this;
                  when 'usr_frg' is not provided, we fall back to
                  'read.query_alignment_length' (aligned span excluding soft
                  clips) and extend from the 5’ end.

        - Filtering toggles:
            + allow_sec: Include/exclude secondary alignments (multi-mappers).
            + allow_sup: Include/exclude supplementary alignments.
            + allow_dup: Include/exclude duplicate-marked alignments.

            These toggles apply uniformly to both paired-end and single-end
            alignments. Keeping 'allow_dup=True' retains duplicate-marked
            proper-pair alignments (e.g., FLAGs 1123 and 1187, which correspond
            to duplicate-marked versions of 99 and 163). If the current data
            lack secondary alignments (which is expected in alignment output
            from the Tsukiyama Lab Bio-protocol workflow), setting
            'allow_sec=True' does nothing.

        - Coordinate handling:
            Intervals are clamped to chromosome bounds '[0, chrom_len]'; zero-
            or negative-length intervals after clamping are skipped.

    Intentional differences from deepTools:
        1. Proper pairs: the user may override TLEN with 'usr_frg'; deepTools
           does not allow this.
        2. No distance guard is applied here, whereas deepTools can fall back
           to single-end-style extension for far or discordant pairs.
        3. Default length for single-end alignments: use
           'read.query_alignment_length' unless 'usr_frg' is provided;
           deepTools’ “extend mode” requires a fixed length.
    """
    frg_tup = defaultdict(list)

    try:
        with pysam.AlignmentFile(path_bam, "rb") as bam_file:
            chrom_sizes = dict(zip(bam_file.references, bam_file.lengths))

            #  Stream the entire file; works for stdin (no index available)
            for read in bam_file.fetch(until_eof=True):
                #  Handle basic exclusions
                if read.is_unmapped or read.reference_id < 0:
                    continue

                #  Handle policy toggles
                if not allow_sec and read.is_secondary:
                    continue
                if not allow_sup and read.is_supplementary:
                    continue
                if not allow_dup and read.is_duplicate:
                    continue

                chrom = bam_file.get_reference_name(read.reference_id)
                chrom_len = chrom_sizes.get(chrom, None)

                #  Handle paired-end alignments: one fragment is emitted per
                #  leftmost anchor in a proper pair
                is_leftmost_pe = (
                    read.is_paired
                    and read.is_proper_pair
                    and read.reference_id == read.next_reference_id
                    and read.template_length > 0
                )

                if is_leftmost_pe:
                    start = read.reference_start
                    tlen = read.template_length

                    frg_len = usr_frg if usr_frg is not None else tlen

                    #  Enforce that user-specified length must be > 0
                    if usr_frg is not None and frg_len <= 0:
                        raise ValueError(
                            "usr_frg must be > 0 for paired-end extension."
                        )

                    frg_start = start
                    frg_end = start + frg_len

                #  Handle single-end alignments: strand-aware extension from
                #  the 5’ end
                elif not read.is_paired:
                    frg_len = (
                        usr_frg
                        if usr_frg is not None
                        else read.query_alignment_length
                    )  # 'query_alignment_length': strictly aligned bases

                    #  Enforce that user-specified fragment length must be > 0,
                    #  and derived zero-length reads get skipped
                    if usr_frg is not None and frg_len <= 0:
                        raise ValueError(
                            "usr_frg must be > 0 for single-end extension."
                        )
                    if usr_frg is None and frg_len <= 0:
                        continue  # Rare: align w/no reference-consuming ops

                    if read.is_reverse:
                        #  5’ end is at 'reference_end' on reverse strand
                        ref_end = (
                            read.reference_end
                            if read.reference_end is not None
                            else read.reference_start
                        )
                        frg_end = ref_end
                        frg_start = frg_end - frg_len
                    else:
                        frg_start = read.reference_start
                        frg_end = frg_start + frg_len
                else:
                    #  Non-leftmost mate of a paired-end pair is ignored by
                    #  design
                    continue

                #  Clamp to contig bounds (if known)
                if chrom_len is not None:
                    if frg_start < 0:
                        frg_start = 0
                    if frg_end > chrom_len:
                        frg_end = chrom_len

                #  Skip degenerate intervals after clamping
                if frg_end <= frg_start:
                    continue

                frg_tup[chrom].append((frg_start, frg_end, frg_len))

    except FileNotFoundError:
        print(f"Error: BAM file '{path_bam}' not found.", file=sys.stderr)
        raise
    except ValueError:
        #  Let caller or outer wrapper handle ValueError details
        raise
    except Exception as e:
        print(
            f"Unexpected error with BAM file '{path_bam}': {e}",
            file=sys.stderr
        )
        raise

    return frg_tup


def calc_sig_chrom(
    chrom: str,
    frg_tup: list[tuple[int, int, int]],
    frg_tot: int,
    siz_bin: int,
    is_len: bool,
    is_norm: bool,
    scl_fct: float | None = None
) -> dict[tuple[str, int], float]:
    """
    Compute binned signal for one chromosome using exact fragment–bin overlap.

    Function respects half-open '[start, end)' fragments and avoids per-base
    loops. Overlap is measured in bases, and output is per-bin sums of per-base
    contributions.

    If provided, 'scl_fct' is applied after any optional normalization.
    ('scl_fct > 0' is required to avoid silent zeroing.)

    Args:
        chrom : str
            Chromosome name.
        frg_tup : list[tuple[int, int, int]]
            List of fragment tuples '(start, end, frg_len)'.
        frg_tot : int
            Total number of fragments (used when 'is_norm=True').
        siz_bin : int
            Bin size in base pairs.
        is_len : bool
            If 'True', normalize by fragment length.
        is_norm : bool
            If 'True', normalize by both fragment length and total fragment
            count so that the genome-wide summed signal is approximately 1.
        scl_fct : float | None
            Optional scaling factor applied to signal.

    Returns:
        sig_bin : dict[tuple[str, int], float]
            A dictionary of binned signal data, where keys are
            '(chrom, bin_start)' and values are per-bin signal scores.

    Raises:
        ValueError
            - If 'siz_bin' <= 0.
            - If 'is_len' or 'is_norm' is True and any 'frg_len' <= 0.
            - If 'is_norm' is True and 'frg_tot' <= 0.
            - If a provided 'scl_fct' <= 0.

    Notes:
        - If 'is_len=True' or 'is_norm=True', each fragment contributes
          '1 / frg_len' per covered base.
        - If 'is_norm=True', signal is additionally divided by 'frg_tot' so
          that genome-wide summed signal is approximately 1.
        - If 'scl_fct' is provided, signal values are scaled accordingly.
        - Signal is accumulated per bin according to the number of overlapping
          bases contributed by each fragment.
    """
    #  Check for invalid bin size
    check_cmp(siz_bin, "gt", 0, "siz_bin", allow_none=False)

    sig_bin = defaultdict(float)

    #  Accumulate by overlap with each bin
    for frg_start, frg_end, frg_len in frg_tup:
        if frg_end <= frg_start:
            continue  # Already filtered upstream, but be safe

        if (is_len or is_norm) and frg_len <= 0:
            #  This shouldn't happen given earlier guards, but fail-fast if it
            #  does
            raise ValueError("'frg_len' must be > 0 when using normalization.")

        #  Compute per-base contribution
        per_bas = (1.0 / frg_len) if (is_len or is_norm) else 1.0

        #  Determine first and last bins “touched” by a half-open fragment:
        #  [frag_start, frag_end)
        fst_bin_start = (frg_start // siz_bin) * siz_bin
        lst_bin_start = ((frg_end - 1) // siz_bin) * siz_bin

        #  Iterate over every bin intersecting [frg_start, frg_end); fst/lst
        #  are bin anchors (multiples of siz_bin, lst inclusive)
        for bin_start in range(fst_bin_start, lst_bin_start + 1, siz_bin):
            bin_end = bin_start + siz_bin
            overlap = min(frg_end, bin_end) - max(frg_start, bin_start)
            if overlap > 0:
                sig_bin[(chrom, bin_start)] += per_bas * overlap

    # #  Optional signal-conservation sanity check
    # if False:  # MAYBE: change to a flag like 'if off_by_one:'
    #     total = sum(sig_bin.values())
    #     expected = (
    #         #  Each fragment should contribute exactly 1 before '÷ frg_tot'
    #         len(frg_tup) if (is_len or is_norm)
    #         #  Otherwise, compute sum of fragment lengths
    #         else sum(e - s for s, e, _ in frg_tup)
    #     )
    #
    #     #  Allow tiny amount of float jitter
    #     assert \
    #         abs(total - expected) < 1e-6 * max(1.0, expected), \
    #         (total, expected)

    #  If requested, divide by total fragment count so that genome-wide summed
    #  signal is approximately 1
    if is_norm:
        if frg_tot <= 0:
            raise ValueError(
                "Normalization requires non-zero total fragments."
            )
        for k in sig_bin:
            sig_bin[k] /= frg_tot

    #  Apply scaling factor if specified
    if scl_fct is not None:
        check_cmp(scl_fct, "gt", 0, "scl_fct")
        for k in sig_bin:
            sig_bin[k] *= scl_fct

    return sig_bin


def calc_sig_task(data):
    """
    Unpack one worker-task tuple and dispatch to 'calc_sig_chrom()'.
    """
    return calc_sig_chrom(*data)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments for computing binned signal from a BAM file.

    Args:
        argv : list[str] | None
            Optional argument vector to parse. If None, use 'sys.argv[1:]'.

    Returns:
        argparse.Namespace
            Parsed command-line arguments.

    Raises:
        SystemExit
            Raised by argparse when '--help' is shown or when required
            arguments or valid choices are missing.

    Notes:
        - Supported CLI options:
            -v, --verbose
                Enable verbose output.
            -t, --threads
                Number of threads for parallel processing (default: 1).
            -i, --infile
                Path to input BAM file or '-'.
            -o, --outfile
                Path to output file or '-'. The extension determines format:
                    - '.bdg', '.bg', '.bedgraph': bedGraph of binned signal.
                    - '.bed': BED-like output of processed fragment
                      coordinates, not signal.
                    - Adding '.gz' to either output type requests gzip-
                      compressed output (e.g., 'output.bdg.gz',
                      'output.bed.gz').
            -of, --outfmt
                Output format hint used only when '--outfile -' (stdout). Must
                be one of 'bedGraph', 'bedgraph', 'bdg', 'bg', or 'bed'.
                Ignored when '--outfile' is a real path, in which case the
                format is inferred from the outfile extension. (This option is
                for direct use of 'compute_signal.py', since the Shell wrappers
                do not support stdout output.)
            -me, --method
                Signal calculation method (default: 'norm'):
                    - unadjusted aliases: 'r', 'raw', 'unadj', 'u',
                      'unadjusted', 's', 'smp', 'simple'
                    - fragment-length-normalized aliases: 'f', 'frg', 'frag',
                      'frg_len', 'frag_len', 'l', 'len', 'len_frg', 'len_frag'
                    - normalized-coverage aliases: 'n', 'nrm', 'norm',
                      'normalized'
                Internally, these are standardized to 'unadj', 'frag', or
                'norm'. Ignored for '.bed'.
            -sb, --siz_bin
                Bin size in base pairs (default: 10). Ignored for '.bed'.
            -sf, --scl_fct
                Optional scaling factor. Ignored for '.bed'.
            -uf, --usr_frg
                Fixed fragment length instead of read/template lengths.
            -dp, --dp, --rnd, --round, --decimals, --digits
                Maximum decimal places retained for finite bedGraph values
                (default: 24). After rounding, non-informative trailing zeros
                are stripped. Ignored for '.bed'.
        - Output format behavior:
            - '.bedGraph', '.bedgraph', '.bdg', and '.bg', optionally followed
              by '.gz', produce bedGraph output.
            - '.bed', optionally followed by '.gz', produces BED-like processed
              fragment-coordinate output rather than bedGraph signal.
            - When '--outfile -' is used, '--outfmt' must be provided because
              the output format cannot be inferred from a filename extension.
    """
    parser = CapArgumentParser(
        description=(
            "Compute binned signal from a BAM file in bedGraph format, "
            "optionally applying normalization.\n"
            "\n"
            "Alternatively, extract and output processed fragment coordinates "
            "in a BED-like format, which can be used as input to the original "
            "implementation of siQ-ChIP "
            "(https://github.com/BradleyDickson/siQ-ChIP) or one updated for "
            "use with S. cerevisiae data "
            "(https://github.com/kalavattam/siQ-ChIP/tree/protocol)."
        )
    )
    add_help_cap(parser)
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        default=False,
        help="Increase output verbosity (stderr banner of parsed args).\n\n"
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=1,
        help=(
            "Number of threads for parallel processing (>= 1; default: "
            "%(default)s).\n"
            "\n"
            "When '--threads > 1', different chromosomes are processed in "
            "parallel.\n\n"
        )
    )
    parser.add_argument(
        "-i", "-f", "--infile", "--file",
        dest="infile",
        required=True,
        help="Path to the input BAM file, or '-' for stdin.\n\n"
    )
    parser.add_argument(
        "-o", "-fo", "--outfile", "--fil_out",
        dest="outfile",
        required=True,
        help=(
            "Path to the output file with extension, or '-' for stdout. "
            "Supported output types are bedGraph ('bedGraph', 'bedgraph', "
            "'bdg', 'bg') and BED ('bed').\n"
            "\n"
            "Append '.gz' for gzip compression, e.g., 'output.bdg.gz'.\n"
            "\n"
            "Note: requesting BED output causes the script to write processed "
            "fragment coordinates in a BED-like format, and '--siz_bin', "
            "'--method', '--scl_fct', and '--rnd' are ignored.\n\n"
        )
    )
    parser.add_argument(
        "-of", "--outfmt",
        choices=["bedGraph", "bedgraph", "bdg", "bg", "bed"],
        default=None,
        help=(
            "Output format hint. Required only when '--outfile -' (stdout). "
            "Choose 'bedGraph', 'bedgraph', 'bdg', 'bg', or 'bed'.\n"
            "\n"
            "Note: stdout output ('--outfile -') is supported only when "
            "calling 'compute_signal.py' directly; the Shell wrappers do not "
            "support this mode.\n\n"
        ),
    )
    parser.add_argument(
        "-me", "--method",
        choices=METHOD_CHOICES,
        default='norm',
        help=(
            "Specify signal calculation type (default: '%(default)s').\n"
            "  - Unadjusted aliases: 'r', 'raw', 'u', 'unadj', 'unadjusted', "
            "'s', 'smp', 'simple'. Internally standardized to 'unadj'.\n"
            "  - Fragment-length-normalized aliases: 'f', 'frg', 'frag', "
            "'frg_len', 'frag_len', 'l', 'len', 'len_frg', 'len_frag'. "
            "Internally standardized to 'frag'.\n"
            "  - Normalized-coverage aliases: 'n', 'nrm', 'norm', "
            "'normalized'. Internally standardized to 'norm'.\n"
            "\n"
            "Note: 'norm' normalizes for both fragment length and total "
            "fragment count so that the genome-wide summed signal is "
            "approximately 1.\n\n"
        )
    )
    parser.add_argument(
        "-sb", "--siz_bin",
        type=int,
        default=10,
        help=(
            "Bin size for signal calculation in base pairs (default: "
            "%(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-sf", "--scl_fct",
        type=float,
        default=None,
        help=(
            "Optional scaling factor to apply to the signal (default: "
            "%(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-uf", "--usr_frg",
        type=int,
        default=None,
        help=(
            "Optional fixed fragment length to use instead of read lengths "
            "(single-end alignments) or template lengths (paired-end "
            "alignments; default: %(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-dp", "--dp", "--rnd", "--round", "--decimals", "--digits",
        dest="rnd",
        type=int,
        default=24,
        help=(
            "Maximum number of decimal places retained for binned signal "
            "values (default: %(default)s). After rounding, non-informative "
            "trailing zeros are stripped.\n\n"
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
            'sys.argv[1:]' is used.

    Returns:
        int.
            On success, returns 0 and writes either a bedGraph of binned signal
            or a BED-like file of processed fragment coordinates (optionally
            gzipped) to '--outfile'.

    Side effects:
        - If '--verbose' is set, prints a human-readable banner of the parsed
          arguments to stderr.
        - Prints human-readable error messages to stderr on validation or I/O
          failures.

    Exits:
        - 0 on success or when showing help with no arguments.
        - 1 on validation or computation errors, including (non-exhaustive):
            + '--siz_bin <= 0' or '--usr_frg <= 0',
            + unsupported/invalid '--outfile' extension,
            + normalization requested with zero fragments,
            + provided '--scl_fct <= 0', and
            + BAM parsing errors (file not found/unreadable).
    """
    #  Parse CLI arguments
    args = parse_args(argv)

    #  Check that input BAM exists
    try:
        if args.infile != "-":
            check_exists(args.infile, kind="file", label="BAM")
    except FileNotFoundError as e:
        #  Print a one-line message and exit cleanly
        raise SystemExit(str(e))

    #  Check and parse output file
    if args.outfile == "-":
        if not args.outfmt:
            raise SystemExit(
                "When '--outfile -' is used, '--outfmt' must also be provided "
                "('bedGraph', 'bedgraph', 'bdg', 'bg', or 'bed')."
            )
        if args.outfmt not in ALLOWED:
            raise SystemExit(
                "Invalid '--outfmt'; choose one of 'bedGraph', 'bedgraph', "
                "'bdg', 'bg', or 'bed'."
            )
        out_fmt = args.outfmt
        outfile = "-"  # Pass through to writers (uncompressed)
    else:
        try:
            #  Parse output path and validate extension
            outfile, out_fmt, _ = check_parse_outfile(args.outfile, ALLOWED)

            #  Check file creatability/writability
            check_writable(outfile, "file")
        except (
            ValueError, FileNotFoundError, PermissionError, IsADirectoryError
        ) as e:
            raise SystemExit(str(e))

    #  Check numeric arguments
    try:
        #  Check threads >= 1
        check_cmp(args.threads, "ge", 1, "threads", allow_none=False)

        #  Only require bin size and rounding for bedGraph paths
        if out_fmt != "bed":
            check_cmp(args.siz_bin, "gt", 0, "siz_bin", allow_none=False)
            check_cmp(args.rnd, "ge", 0, "rnd", allow_none=False)

        #  Handle optional scalars: allow None but must be > 0 when provided
        check_cmp(args.scl_fct, "gt", 0, "scl_fct", allow_none=True)
        check_cmp(args.usr_frg, "gt", 0, "usr_frg", allow_none=True)

    except ValueError as e:
        raise SystemExit(str(e))

    #  Preserve the user-supplied '--method' token, then standardize it to the
    #  canonical internal name
    mthd_in = args.method
    args.method = METHOD_CANON[args.method]

    #  Determine whether to normalize based on canonicalized method
    is_len = (args.method == "frag")
    is_norm = (args.method == "norm")

    #  Print verbose output
    if args.verbose:
        with redirect_stdout(sys.stderr):
            print("")
            print("####################################")
            print("## Arguments for 'compute_signal' ##")
            print("####################################")
            print("")
            print("--verbose")
            print(f"--threads {args.threads}")
            print(f"--infile  {args.infile}")
            print(f"--outfile {outfile}")
            if outfile == "-":
                print(f"--outfmt  {out_fmt}")
            if out_fmt == "bed":
                print(f"--usr_frg {args.usr_frg}")
                print(
                    "\n\n(BED output mode: signal computation arguments "
                    "ignored)\n"
                )
            else:
                if mthd_in != args.method:
                    print(
                        f"--method  {mthd_in}  (standardized internally to "
                        f"{args.method})"
                    )
                else:
                    print(f"--method  {args.method}")
                print(f"--siz_bin {args.siz_bin}")
                print(f"--scl_fct {args.scl_fct}")
                print(f"--usr_frg {args.usr_frg}")
                print(f"--rnd     {args.rnd}")
            print("")
            print("")

    #  Parse and process BAM file
    try:
        frg_tup = parse_bam(args.infile, args.usr_frg)
    except FileNotFoundError:
        raise SystemExit(f"BAM not found: {args.infile}")
    except ValueError as e:
        raise SystemExit(str(e))
    except OSError as e:
        raise SystemExit(f"I/O error while reading BAM: {e}")

    try:
        if out_fmt == 'bed':
            #  Write fragments processed by 'parse_bam()' in a BED-like format,
            #  omitting signal calculations, etc.
            with open_out(outfile) as bed_file:
                for chrom in sorted(frg_tup.keys(), key=sort_chrom):
                    for start, end, length in sorted(
                        frg_tup[chrom], key=lambda t: t[0]
                    ):
                        bed_file.write(f"{chrom}\t{start}\t{end}\t{length}\n")
            return 0

        #  Otherwise, compute and write bedGraph signal from fragments
        #  processed by 'parse_bam()'
        frg_tot = sum(len(entries) for entries in frg_tup.values())
        if is_norm and frg_tot == 0:
            raise ValueError(
                "Normalization requires non-zero total fragments. Check the "
                "BAM file."
            )

        #  Prepare and execute parallel tasks (when '--threads > 1')
        tsk_dat = [
            (
                chrom, entries, frg_tot, args.siz_bin, is_len, is_norm,
                args.scl_fct
            )
            for chrom, entries in frg_tup.items()
        ]

        #  Compute per-chromosome signal, in parallel if '--threads > 1'
        if args.threads == 1:
            results = (calc_sig_task(d) for d in tsk_dat)
        else:
            with ProcessPoolExecutor(
                max_workers=min(args.threads, os.cpu_count() or 1)
            ) as executor:
                results = executor.map(calc_sig_task, tsk_dat)

        #  Combine results from all processes
        sig_cmb = defaultdict(float)
        for result in results:
            for key, value in result.items():
                sig_cmb[key] += value

        #  Write bedGraph output
        write_bdg(sig_cmb, outfile, args.siz_bin, args.rnd)

        return 0

    except ValueError as e:
        raise SystemExit(str(e))

    except OSError as e:
        raise SystemExit(f"I/O error: {e}")


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
