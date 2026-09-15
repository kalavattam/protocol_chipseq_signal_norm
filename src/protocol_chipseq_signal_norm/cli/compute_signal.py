#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: compute_signal.py
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
#   GPT-5.6);
# - Anthropic Claude Code (Opus 5, Fable 5).
#
# Distributed under the MIT license.


"""
Calculate binned signal or fragment-coordinate output from BAM/CRAM input.

The CLI accepts input, reference, and output paths, threading, the signal
method, bin and window sizes, the processing engine, fragment length, scaling,
reporting, and value formatting.

It writes bedGraph-like signal tracks, BED-like fragment-coordinate records, or
fragment ('--report_n_frg') and bin ('--report_n_bin') counts with or without a
track.

When writing bedGraph output, finite values are rounded to at most '--dp'
decimal places and trailing zeros are stripped.

Examples
--------
python -m protocol_chipseq_signal_norm.cli.compute_signal \\
    --fil_in <file> --fil_out <file> [options]
python -m protocol_chipseq_signal_norm.cli.compute_signal \\
    --fil_in <file> --fil_out <file> --report_n_frg --report_n_bin [options]
python -m protocol_chipseq_signal_norm.cli.compute_signal \\
    --fil_in <file> --report_n_frg <file> --report_n_bin <file> [options]
"""

from __future__ import annotations

import argparse
import json
import math
import os
import signal
import sys
import time
from collections import defaultdict
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import redirect_stdout, suppress
from typing import Any

import numpy as np
import pysam

from protocol_chipseq_signal_norm.utilities.utils_check import (
    ALLOWED_OUTPUT_FORMATS,
    check_exists,
    check_writable,
    validate_comparison,
    validate_output_path,
)
from protocol_chipseq_signal_norm.utilities.utils_chrom import sort_chrom
from protocol_chipseq_signal_norm.utilities.utils_cli import (
    CapArgumentParser,
    add_help_cap,
)
from protocol_chipseq_signal_norm.utilities.utils_io import open_out

with suppress(AttributeError, ValueError):
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

assert sys.version_info >= (3, 11), "Python >= 3.11 required."

# TODO: In the performance pass, reassess a compiled core. Separately, assert
# emitted terminal-bin end coordinates never exceed chromosome size, which is a
# key property required by downstream bigWig converters.

# Map accepted '--method' values to canonical internal names.
# fmt: off
METHOD_CANON = {
    # Compute unadjusted signal.
    "r": "unadj",
    "raw": "unadj",
    "u": "unadj",
    "unadj": "unadj",
    "unadjusted": "unadj",
    "s": "unadj",
    "smp": "unadj",
    "simple": "unadj",

    # Normalize signal by fragment length.
    "f": "frag",
    "frg": "frag",
    "frag": "frag",
    "frg_len": "frag",
    "frag_len": "frag",
    "l": "frag",
    "len": "frag",
    "len_frg": "frag",
    "len_frag": "frag",

    # Normalize signal by fragment length and depth.
    "n": "norm",
    "nc": "norm",
    "nrm": "norm",
    "norm": "norm",
    "normalized": "norm",
}
# fmt: on
METHOD_CHOICES = tuple(METHOD_CANON.keys())
ENGINE_CHOICES = ("chrom", "window")
ENGINE_STRAT = {
    "chrom": "idx_chrom",
    "window": "idx_win",
}
MODE_EXEC_CHOICES = ("map", "as_completed")
STRAT_BED_CHOICES = (
    "auto",
    "serial",
    "idx_chrom",
    "idx_win",
)
STRAT_WRITER_CHOICES = ("serial", "parallel_ordered")

# Stand in for a report path the user did not spell out, which is what a bare
# '--report_n_frg' or '--report_n_bin' supplies. A NUL byte cannot occur in a
# path, so the sentinel can never collide with one a user typed.
REPORT_DERIVE = "\x00derive"


def get_siz_chr(
    fil_aln: str,
    ref_fa: str | None = None,
) -> dict[str, int]:
    """
    Read chromosome sizes from a BAM or CRAM header.
    """

    kwargs = {}

    if ref_fa is not None:
        kwargs["reference_filename"] = ref_fa

    with pysam.AlignmentFile(fil_aln, "rb", **kwargs) as alignment_file:
        return {
            chrom: size
            for chrom, size in zip(
                alignment_file.references,
                alignment_file.lengths,
                strict=True,
            )
            if size is not None and size > 0
        }


def resolve_siz_chr(siz_chr_hdr: dict[str, int]) -> dict[str, int]:
    """
    Return chromosome sizes from a BAM/CRAM header, rejecting an empty set.
    """

    if not siz_chr_hdr:
        raise ValueError(
            "Chromosome sizes are required to trim bedGraph bins. Provide a "
            "BAM/CRAM whose header carries usable sequence lengths ('LN' "
            "tags).",
        )

    return dict(siz_chr_hdr)


def resolve_report_path(
    value: str | None,
    fil_out: str | None,
    label: str,
    flag: str,
) -> str | None:
    """
    Resolve a report path, deriving it from the output path when unspelled.

    Parameters
    ----------
    value : str | None
        Report path as parsed, 'REPORT_DERIVE' for a bare flag, or None.
    fil_out : str | None
        Validated output path or None in report-only mode.
    label : str
        Count label placed before '.txt', either 'n_frg' or 'n_bin'.
    flag : str
        Option spelling named in the error, e.g., '--report_n_frg'.

    Returns
    -------
    path_report : str | None
        Path to write, unchanged unless it was derived.

    Raises
    ------
    ValueError
        If a bare flag is given without an output path to derive from.
    """

    if value != REPORT_DERIVE:
        return value

    if fil_out is None:
        raise ValueError(
            f"'{flag}' was given without a path, so the report path is "
            f"derived from '--fil_out', but no '--fil_out' was given. Supply "
            f"a path to '{flag}' or an output path to '--fil_out'.",
        )

    # Strip one '.gz' and then one extension, matching how the wrappers derive
    # the same name, so a sample's counts sit beside its track either way.
    base = fil_out[:-3] if fil_out.endswith(".gz") else fil_out

    return f"{os.path.splitext(base)[0]}.{label}.txt"


def start_profile(
    args: argparse.Namespace,
    fil_out: str,
    fmt_out: str,
) -> dict | None:
    """
    Initialize optional timing metadata for internal profiling.
    """

    if args.profile_json is None:
        return None

    return {
        "version": 1,
        "script": "compute_signal.py",
        "engine": args.engine,
        "mode_exec": args.mode_exec,
        "strat_bed": args.strat_bed,
        "strat_writer": args.strat_writer,
        "wrk_writer": args.wrk_writer,
        "threads": args.threads,
        "fil_in": args.fil_in,
        "ref_fa": args.ref_fa,
        "fil_out": fil_out,
        "fmt_out": fmt_out,
        "method": args.method,
        "siz_bin": args.siz_bin,
        "siz_win": args.siz_win,
        "phases_s": {},
        "tasks": [],
    }


def record_phase(profile: dict | None, name: str, start: float) -> None:
    """
    Record elapsed seconds for one named phase when profiling is active.
    """

    if profile is not None:
        profile["phases_s"][name] = time.perf_counter() - start


def write_profile(path: str | None, profile: dict | None) -> None:
    """
    Write optional profiling metadata as JSON.
    """

    if path is None or profile is None:
        return

    with open(path, "w", encoding="utf-8") as handle:
        json.dump(profile, handle, indent=2, sort_keys=True)
        handle.write("\n")


def count_result_bins(result: Any) -> int:
    """
    Count signal bins in worker results and merged sparse containers.

    Parameters
    ----------
    result : Any
        Tagged sparse signal result, optionally paired with an integer fragment
        count.

    Returns
    -------
    count : int
        Number of represented signal bins.
    """

    if isinstance(result, tuple) and result:
        if len(result) == 2 and isinstance(result[1], int):
            return count_result_bins(result[0])

        tag = result[0]
        if tag in ("direct_sparse_np", "sparse_merged"):
            return sum(int(starts.size) for _, starts, _ in result[1])

    return len(result)


def est_result_bytes(result: Any) -> int:
    """
    Estimate NumPy-array payload bytes for a tagged sparse signal result.

    Parameters
    ----------
    result : Any
        Tagged sparse signal result.

    Returns
    -------
    byte_count : int
        Sum of owned NumPy array payload bytes.
    """

    if isinstance(result, tuple) and result:
        tag = result[0]
        if tag in ("direct_sparse_np", "sparse_merged"):
            return sum(
                int(starts.nbytes + values.nbytes)
                for _, starts, values in result[1]
            )

    return 0


def check_idx_aln(
    fil_aln: str,
    ref_fa: str | None,
) -> None:
    """
    Validate prerequisites for indexed signal engines.

    Parameters
    ----------
    fil_aln : str
        BAM or CRAM path that must have a corresponding index.
    ref_fa : str | None
        Reference FASTA required for indexed CRAM access.

    Raises
    ------
    ValueError
        If the path, index, or required CRAM reference is absent.
    """

    if fil_aln == "-":
        raise ValueError(
            "Indexed signal engines require a named BAM or CRAM path.",
        )

    if fil_aln.lower().endswith(".cram") and ref_fa is None:
        raise ValueError(
            "Indexed CRAM signal engines require '--ref_fa' so every worker "
            "can open the CRAM deterministically.",
        )

    kwargs = {}

    if ref_fa is not None:
        kwargs["reference_filename"] = ref_fa

    with pysam.AlignmentFile(fil_aln, "rb", **kwargs) as alignment_file:
        if not alignment_file.has_index():
            raise ValueError(
                "Indexed signal engines require an alignment index "
                "(.bai for BAM or .crai for CRAM).",
            )


def has_idx_aln(
    fil_aln: str,
    ref_fa: str | None = None,
) -> bool:
    """
    Return True when an alignment can be accessed through coordinate indexes.
    """

    kwargs = {}

    if ref_fa is not None:
        kwargs["reference_filename"] = ref_fa

    try:
        with pysam.AlignmentFile(
            fil_aln,
            "rb",
            **kwargs,
        ) as alignment_file:
            return bool(alignment_file.has_index())
    except (OSError, ValueError):
        return False


def keep_read(
    read: pysam.AlignedSegment,
    allow_sec: bool = True,
    allow_supp: bool = False,
    allow_dup: bool = True,
) -> bool:
    """
    Return whether a read passes the shared alignment-filter policy.
    """

    if read.is_unmapped or read.reference_id < 0:
        return False

    if not allow_sec and read.is_secondary:
        return False

    if not allow_supp and read.is_supplementary:
        return False

    return not (not allow_dup and read.is_duplicate)


def read_to_frg(
    read: pysam.AlignedSegment,
    get_reference_name: Callable[[int], str],
    siz_chr: dict[str, int],
    usr_frg: int | None = None,
) -> tuple[str, int, int, int] | None:
    """
    Convert one accepted alignment record to a processed fragment interval.

    Parameters
    ----------
    read : pysam.AlignedSegment
        Accepted alignment record.
    get_reference_name : Callable[[int], str]
        Resolver for numeric reference identifiers.
    siz_chr : dict[str, int]
        Maximum coordinate for each accepted chromosome.
    usr_frg : int | None
        Optional fragment length used to extend alignments.

    Returns
    -------
    fragment : tuple[str, int, int, int] | None
        Chromosome, clipped start, clipped end, and fragment length, or 'None'
        when the alignment does not produce an interval.

    Raises
    ------
    ValueError
        If extension length or chromosome-size data is invalid.
    """

    chrom = get_reference_name(read.reference_id)
    chrom_len = siz_chr.get(chrom)
    if chrom_len is None:
        raise ValueError(
            f"Chromosome {chrom!r} is missing from the alignment header's "
            "usable sequence lengths; fix the header before rerunning.",
        )

    # Handle paired-end alignments: one fragment is emitted per leftmost anchor
    # in a proper pair.
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

        if usr_frg is not None and frg_len <= 0:
            raise ValueError("usr_frg must be > 0 for paired-end extension.")

        frg_start = start
        frg_end = start + frg_len

    elif not read.is_paired:
        frg_len = (
            usr_frg if usr_frg is not None else read.query_alignment_length
        )

        if usr_frg is not None and frg_len <= 0:
            raise ValueError("usr_frg must be > 0 for single-end extension.")

        if usr_frg is None and frg_len <= 0:
            return None

        if read.is_reverse:
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
        return None

    if frg_start < 0:
        frg_start = 0

    if frg_end > chrom_len:
        frg_end = chrom_len

    if frg_end <= frg_start:
        return None

    return chrom, frg_start, frg_end, frg_len


def iter_aln_frg(
    fil_aln: str,
    siz_chr: dict[str, int],
    usr_frg: int | None = None,
    ref_fa: str | None = None,
    allow_sec: bool = True,
    allow_supp: bool = False,
    allow_dup: bool = True,
) -> Iterator[tuple[str, int, int, int]]:
    """
    Stream processed fragment intervals from a BAM or CRAM file.

    Parameters
    ----------
    fil_aln : str
        BAM or CRAM input path.
    siz_chr : dict[str, int]
        Maximum coordinate for each accepted chromosome.
    usr_frg : int | None
        Optional fragment length used to extend alignments.
    ref_fa : str | None
        Reference FASTA used for CRAM decoding.
    allow_sec, allow_supp, allow_dup : bool
        Whether to retain secondary, supplementary, and duplicate records.

    Yields
    ------
    fragment : tuple[str, int, int, int]
        Chromosome, clipped start, clipped end, and fragment length.
    """

    kwargs = {}

    if ref_fa is not None:
        kwargs["reference_filename"] = ref_fa

    with pysam.AlignmentFile(fil_aln, "rb", **kwargs) as alignment_file:
        for read in alignment_file.fetch(until_eof=True):
            if not keep_read(
                read,
                allow_sec,
                allow_supp,
                allow_dup,
            ):
                continue

            fragment = read_to_frg(
                read,
                alignment_file.get_reference_name,
                siz_chr,
                usr_frg,
            )
            if fragment is not None:
                yield fragment


def iter_idx_frg(
    fil_aln: str,
    chrom: str,
    siz_chr: dict[str, int],
    usr_frg: int | None = None,
    ref_fa: str | None = None,
    start: int | None = None,
    end: int | None = None,
    allow_sec: bool = True,
    allow_supp: bool = False,
    allow_dup: bool = True,
) -> Iterator[tuple[str, int, int, int]]:
    """
    Stream processed fragments from one indexed chromosome or window fetch.

    Parameters
    ----------
    fil_aln : str
        Indexed BAM or CRAM input path.
    chrom : str
        Chromosome to fetch.
    siz_chr : dict[str, int]
        Maximum coordinate for each accepted chromosome.
    usr_frg : int | None
        Optional fragment length used to extend alignments.
    ref_fa : str | None
        Reference FASTA used for CRAM decoding.
    start, end : int | None
        Optional half-open fetch bounds.
    allow_sec, allow_supp, allow_dup : bool
        Whether to retain secondary, supplementary, and duplicate records.

    Yields
    ------
    fragment : tuple[str, int, int, int]
        Chromosome, clipped start, clipped end, and fragment length.
    """

    kwargs = {}

    if ref_fa is not None:
        kwargs["reference_filename"] = ref_fa

    with pysam.AlignmentFile(fil_aln, "rb", **kwargs) as alignment_file:
        if start is None and end is None:
            reads = alignment_file.fetch(chrom)
        else:
            reads = alignment_file.fetch(chrom, start, end)

        for read in reads:
            if not keep_read(
                read,
                allow_sec,
                allow_supp,
                allow_dup,
            ):
                continue

            if (
                start is not None
                and end is not None
                and (
                    read.reference_start < start or read.reference_start >= end
                )
            ):
                # Partition by read anchor, not computed fragment start. This
                # prevents double-counting while preserving reverse-strand
                # single-end fragments that extend left of a window.
                continue

            fragment = read_to_frg(
                read,
                alignment_file.get_reference_name,
                siz_chr,
                usr_frg,
            )
            if fragment is not None:
                yield fragment


def count_frgs_and_bins(
    fragments: Iterator[tuple[str, int, int, int]],
    siz_bin: int,
) -> tuple[int, int]:
    """
    Count fragments ('N') and the bins they span in total ('L').

    Parameters
    ----------
    fragments : Iterator[tuple[str, int, int, int]]
        Processed fragments as '(chrom, start, end, length)', from
        'iter_aln_frg' or 'iter_idx_frg'.
    siz_bin : int
        Bin size in base pairs.

    Returns
    -------
    n_frg, n_bin : tuple[int, int]
        The number of fragments ('N') and the summed count of bins the
        fragments touch ('L').

    Notes
    -----
    - 'L' counts a fragment once per bin it touches.
      + That is deliberately not what the signal accumulation does: the
        accumulators add base-pair overlap, so an unadjusted track sums to
        total base pairs rather than to total spanned bins, while a
        fragment-length-normalized track sums to 'N'.
      + Both are measures of span, but only the bin count makes 'k = L / N'
        come out in bins, which is the unit 'compute_pseudo' needs for a
        pseudocount added to per-bin coverage.
    - Because the fragments arrive from the same iterators the signal path
      uses, 'N' is the denominator a '--method norm' run divides by, and both
      numbers describe the same population of read alignments.
    - Fragments with a nonpositive span are skipped, so they raise neither
      count.
    - The carried length field is ignored; both counts come from the start and
      end coordinates alone.
    - A half-open fragment '[start, end)' touches bins
      'start // siz_bin' through '(end - 1) // siz_bin' inclusive.
    """

    n_frg = 0
    n_bin = 0

    for _chrom, frg_start, frg_end, _length in fragments:
        if frg_end <= frg_start:
            continue

        n_frg += 1
        n_bin += ((frg_end - 1) // siz_bin) - (frg_start // siz_bin) + 1

    return n_frg, n_bin


def collect_frg_arr(
    fragments: Iterator[tuple[str, int, int, int]],
) -> tuple[dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]], int]:
    """
    Collect processed fragments into per-chromosome NumPy arrays.
    """

    starts: dict[str, list[int]] = defaultdict(list)
    ends: dict[str, list[int]] = defaultdict(list)
    lengths: dict[str, list[int]] = defaultdict(list)
    n_frg = 0

    for chrom, frg_start, frg_end, frg_len in fragments:
        starts[chrom].append(frg_start)
        ends[chrom].append(frg_end)
        lengths[chrom].append(frg_len)
        n_frg += 1

    arrays: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    for chrom in starts:
        arrays[chrom] = (
            np.asarray(starts[chrom], dtype=np.int64),
            np.asarray(ends[chrom], dtype=np.int64),
            np.asarray(lengths[chrom], dtype=np.float64),
        )

    return arrays, n_frg


def collect_bed_arr(
    fragments: Iterator[tuple[str, int, int, int]],
) -> tuple[str, list[tuple[str, np.ndarray, np.ndarray, np.ndarray]], int]:
    """
    Collect processed fragments into per-chromosome BED output arrays.
    """

    starts: dict[str, list[int]] = defaultdict(list)
    ends: dict[str, list[int]] = defaultdict(list)
    lengths: dict[str, list[int]] = defaultdict(list)
    n_frg = 0

    for chrom, frg_start, frg_end, frg_len in fragments:
        starts[chrom].append(frg_start)
        ends[chrom].append(frg_end)
        lengths[chrom].append(frg_len)
        n_frg += 1

    arrays = [
        (
            chrom,
            np.asarray(starts[chrom], dtype=np.int64),
            np.asarray(ends[chrom], dtype=np.int64),
            np.asarray(lengths[chrom], dtype=np.int64),
        )
        for chrom in starts
    ]

    return "bed_np", arrays, n_frg


def parse_bam(
    fil_aln: str,
    usr_frg: int | None = None,
    ref_fa: str | None = None,
    siz_chr: dict[str, int] | None = None,
    allow_sec: bool = True,
    allow_supp: bool = False,
    allow_dup: bool = True,
) -> dict[str, list[tuple[int, int, int]]]:
    """
    Parse a BAM or CRAM into chromosome-grouped fragment intervals.

    Parameters
    ----------
    fil_aln : str
        Input BAM or CRAM path.
    usr_frg : int | None
        Optional fixed fragment length override. If provided:
            - paired-end read alignments: overrides TLEN for leftmost anchors.
            - single-end read alignments: used for both strands from 5’ end.
    ref_fa : str | None
        Reference FASTA file for CRAM decoding.
    siz_chr : dict[str, int] | None
        Optional chromosome sizes used to filter and clamp fragments.
    allow_sec : bool
        Include secondary (multi-mapping) alignments. Default True.
    allow_supp : bool
        Include supplementary alignments. Default False.
    allow_dup : bool
        Include duplicate-marked alignments. Default True.

    Returns
    -------
    frg_tup : dict[str, list[tuple[int, int, int]]]
        Chromosome-keyed lists of start, end, and fragment-length tuples. Every
        interval is half-open and uses zero-based coordinates.

    Raises
    ------
    FileNotFoundError
        If 'fil_aln' does not exist.
    ValueError
        If '--usr_frg' is nonpositive for paired-end or single-end extension.

    Notes
    -----
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
              and the fragment length is TLEN unless overridden by '--usr_frg'.
            + In the current implementation, only leftmost proper-pair anchors
              with 'TLEN > 0' are treated as paired-end fragment anchors.

            NOTE: This intentionally differs from deepTools, where proper pairs
                  always use the observed TLEN (also, deepTools applies a
                  distance guard). Here, the user may override TLEN with
                  '--usr_frg' by design (also, no distance guard is in place).

        - Single-end read alignments ('read.is_paired == False') are extended
          strand-aware from the 5’ end to the selected fragment length. This is
          '--usr_frg' when provided and 'read.query_alignment_length'
          otherwise.
            + Forward strand (not 'read.is_reverse'):
                  [start, end) = [reference_start, reference_start + length)
            + Reverse strand ('read.is_reverse'):
                  [start, end) = [reference_end - length, reference_end)

            NOTE: A fixed fragment length ('--usr_frg') is preferred, mirroring
                  deepTools’ '--extendReads', as it produces aligner- and
                  read-length–independent coverage. We do not enforce this;
                  when '--usr_frg' is not provided, we fall back to
                  'read.query_alignment_length' (aligned span excluding soft
                  clips) and extend from the 5’ end.

        - Filtering toggles:
            + 'allow_sec': Include or exclude secondary alignments.
            + 'allow_supp': Include or exclude supplementary alignments.
            + 'allow_dup': Include or exclude duplicate-marked alignments.

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
        1. Proper pairs: the user may override TLEN with '--usr_frg'; deepTools
           does not allow this.
        2. No distance guard is applied here, whereas deepTools can fall back
           to single-end-style extension for far or discordant pairs.
        3. Default length for single-end alignments: use
           'read.query_alignment_length' unless '--usr_frg' is provided;
           deepTools’ “extend mode” requires a fixed length.
    """

    frg_tup = defaultdict(list)

    try:
        if siz_chr is None:
            siz_chr_hdr = get_siz_chr(fil_aln, ref_fa)
            siz_chr = resolve_siz_chr(siz_chr_hdr)

        for (
            chrom,
            frg_start,
            frg_end,
            frg_len,
        ) in iter_aln_frg(
            fil_aln=fil_aln,
            siz_chr=siz_chr,
            usr_frg=usr_frg,
            ref_fa=ref_fa,
            allow_sec=allow_sec,
            allow_supp=allow_supp,
            allow_dup=allow_dup,
        ):
            frg_tup[chrom].append(
                (frg_start, frg_end, frg_len),
            )

    except FileNotFoundError:
        print(
            f"Error: BAM or CRAM file '{fil_aln}' not found.",
            file=sys.stderr,
        )

        raise
    except ValueError:
        raise
    except Exception as error:
        print(
            f"Unexpected error with BAM or CRAM file '{fil_aln}': {error}",
            file=sys.stderr,
        )

        raise

    return frg_tup


def calc_sig_chrom_direct_sparse_np(
    chrom: str,
    starts: np.ndarray,
    ends: np.ndarray,
    lengths: np.ndarray,
    chrom_size: int,
    siz_bin: int,
    is_len: bool,
) -> object:
    """
    Compute binned signal and return sparse NumPy arrays without a dict.

    Parameters
    ----------
    chrom : str
        Chromosome identifier for the selected signal calculation.
    starts : np.ndarray
        Zero-based interval start coordinates.
    ends : np.ndarray
        Zero-based interval end coordinates.
    lengths : np.ndarray
        Fragment lengths aligned with the interval coordinates.
    chrom_size : int
        Chromosome length in base pairs.
    siz_bin : int
        Positive signal-bin width in base pairs.
    is_len : bool
        Whether signal weights are normalized by fragment length.

    Returns
    -------
    result : object
        Strategy tag and sparse chromosome result payload.

    Raises
    ------
    ValueError
        If a supplied value violates the validated contract.
    """

    validate_comparison(siz_bin, "gt", 0, "siz_bin", allow_none=False)

    if chrom_size <= 0:
        raise ValueError(f"Chromosome size must be > 0 for {chrom!r}.")

    n_bins = math.ceil(chrom_size / siz_bin)
    if starts.size == 0:
        return "direct_sparse_np", []

    valid = ends > starts

    if not np.all(valid):
        starts = starts[valid]
        ends = ends[valid]
        lengths = lengths[valid]

    if starts.size == 0:
        return "direct_sparse_np", []

    if is_len:
        if np.any(lengths <= 0):
            raise ValueError("'frg_len' must be > 0 when using normalization.")

        weights = 1.0 / lengths
    else:
        weights = np.ones(starts.shape, dtype=np.float64)

    sig = np.zeros(n_bins, dtype=np.float64)
    start_bins = starts // siz_bin
    end_bins = (ends - 1) // siz_bin

    same = start_bins == end_bins

    if np.any(same):
        same_values = (ends[same] - starts[same]) * weights[same]

        np.add.at(sig, start_bins[same], same_values)

    multi = ~same

    if np.any(multi):
        starts_subset = starts[multi]
        ends_subset = ends[multi]
        start_bins_subset = start_bins[multi]
        end_bins_subset = end_bins[multi]
        weights_subset = weights[multi]

        left_overlap = (
            (start_bins_subset + 1) * siz_bin - starts_subset
        ) * weights_subset
        right_overlap = (
            ends_subset - end_bins_subset * siz_bin
        ) * weights_subset

        np.add.at(sig, start_bins_subset, left_overlap)
        np.add.at(sig, end_bins_subset, right_overlap)

        has_interior = end_bins_subset > start_bins_subset + 1

        if np.any(has_interior):
            first = start_bins_subset[has_interior] + 1
            last_exclusive = end_bins_subset[has_interior]
            value = siz_bin * weights_subset[has_interior]

            diff = np.zeros(n_bins + 1, dtype=np.float64)
            np.add.at(diff, first, value)
            np.add.at(diff, last_exclusive, -value)

            # The 'int64' interior sum is exact; 'float64' rounding can leave
            # ~1e-20 residue in true-zero bins that 'sig != 0.0' would emit.
            # Masking by the integer sum restores exact zero, tolerance-free.
            interior_diff = np.zeros(n_bins + 1, dtype=np.int64)
            np.add.at(interior_diff, first, 1)
            np.add.at(interior_diff, last_exclusive, -1)
            interior = np.cumsum(interior_diff[:-1]) > 0

            sig += np.where(interior, np.cumsum(diff[:-1]), 0.0)

    idx = np.flatnonzero(sig != 0.0)

    if idx.size == 0:
        return "direct_sparse_np", []

    return (
        "direct_sparse_np",
        [
            (
                chrom,
                (idx * siz_bin).astype(np.int64, copy=False),
                sig[idx].astype(np.float64, copy=False),
            ),
        ],
    )


def calc_sig_idx_fetch_task(
    data: tuple[object, ...],
) -> tuple[object, int]:
    """
    Fetch one indexed region, parse fragments, and compute unscaled signal.

    Parameters
    ----------
    data : tuple[object, ...]
        Serialized indexed-fetch task arguments.

    Returns
    -------
    sig, n_frg : tuple[object, int]
        Tagged sparse chromosome result and the accepted fragment count.
    """

    (
        fil_aln,
        chrom,
        siz_chr,
        usr_frg,
        ref_fa,
        start,
        end,
        siz_bin,
        is_len,
    ) = data

    frg_arr, n_frg = collect_frg_arr(
        iter_idx_frg(
            fil_aln=fil_aln,
            chrom=chrom,
            siz_chr=siz_chr,
            usr_frg=usr_frg,
            ref_fa=ref_fa,
            start=start,
            end=end,
        ),
    )

    arrays = frg_arr.get(chrom)

    if arrays is None:
        sig = ("direct_sparse_np", [])
    else:
        sig = calc_sig_chrom_direct_sparse_np(
            chrom=chrom,
            starts=arrays[0],
            ends=arrays[1],
            lengths=arrays[2],
            chrom_size=siz_chr[chrom],
            siz_bin=siz_bin,
            is_len=is_len,
        )

    return sig, n_frg


def collect_bed_idx_task(data: tuple[object, ...]) -> object:
    """
    Fetch one indexed region and return compact BED row arrays.
    """

    (
        fil_aln,
        chrom,
        siz_chr,
        usr_frg,
        ref_fa,
        start,
        end,
    ) = data

    return collect_bed_arr(
        iter_idx_frg(
            fil_aln=fil_aln,
            chrom=chrom,
            siz_chr=siz_chr,
            usr_frg=usr_frg,
            ref_fa=ref_fa,
            start=start,
            end=end,
        ),
    )


def count_bed_array_rows(result: Any) -> int:
    """
    Count rows in one compact BED array result.
    """

    if not isinstance(result, tuple) or len(result) < 2:
        return 0

    return sum(
        int(starts.size) for _chrom, starts, _ends, _lengths in result[1]
    )


def est_bed_bytes(result: Any) -> int:
    """
    Estimate NumPy payload bytes in one compact BED array result.
    """

    if not isinstance(result, tuple) or len(result) < 2:
        return 0

    return sum(
        int(starts.nbytes + ends.nbytes + lengths.nbytes)
        for _chrom, starts, ends, lengths in result[1]
    )


# TODO: Revisit in the post-M8-c performance pass. This wrapper, the
# '--profile_json' layer, and the four hidden operator knobs ('--mode_exec',
# '--strat_bed', '--strat_writer', '--wrk_writer') exist for that pass; its
# exit decision is to adopt them as supported surface or to move them to
# 'dev/algorithm_testing/' (decided 2026-09-12).
def calc_sig_profile_task(data: tuple[object, ...]) -> object:
    """
    Compute one profiled worker task and return timing metadata.

    Parameters
    ----------
    data : tuple[object, ...]
        Task kind, payload, and profiling metadata.

    Returns
    -------
    result : object
        Worker result paired with timing and process metadata.

    Raises
    ------
    ValueError
        If the task kind is unknown.
    """

    kind, task_id, payload = data
    start = time.perf_counter()

    if kind in ("idx_chrom", "idx_win"):
        result = calc_sig_idx_fetch_task(payload)
    elif kind == "indexed_bed":
        result = collect_bed_idx_task(payload)
    else:
        raise ValueError(f"Unknown profiled task kind: {kind!r}.")

    end = time.perf_counter()
    result_count = (
        count_bed_array_rows(result)
        if kind == "indexed_bed"
        else count_result_bins(result)
    )

    return (
        task_id,
        {
            "worker_start_s": start,
            "worker_end_s": end,
            "worker_s": end - start,
            "result_bins": result_count,
        },
        result,
    )


def iter_task_results(
    kind: str,
    task_func: Callable[[object], object],
    task_data: Iterable[object],
    threads: int,
    mode_exec: str,
    task_profiles: list[dict] | None = None,
) -> Iterator[object]:
    """
    Yield task results, optionally recording per-task worker timings.

    Parameters
    ----------
    kind : str
        Worker task kind.
    task_func : Callable[[object], object]
        Callable executed for each task payload.
    task_data : Iterable[object]
        Ordered task payloads.
    threads : int
        Maximum worker count.
    mode_exec : str
        Ordered-map or completion-order execution mode.
    task_profiles : list[dict] | None
        Optional collection receiving worker timing records.

    Yields
    ------
    result : object
        Task results in the order defined by 'mode_exec'.
    """

    if threads == 1:
        for task_id, payload in enumerate(task_data):
            if task_profiles is None:
                yield task_func(payload)
            else:
                start = time.perf_counter()
                result = task_func(payload)
                end = time.perf_counter()
                received = time.perf_counter()
                task_profiles[task_id].update(
                    {
                        "worker_start_s": start,
                        "worker_end_s": end,
                        "worker_s": end - start,
                        "recv_s": received,
                        "recv_lag_s": received - end,
                        "result_bins": (
                            count_bed_array_rows(result)
                            if kind == "indexed_bed"
                            else count_result_bins(result)
                        ),
                    },
                )

                yield result

        return

    max_workers = min(threads, os.cpu_count() or 1)

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        if task_profiles is None:
            if mode_exec == "map":
                yield from executor.map(task_func, task_data)
            else:
                future_to_id = {
                    executor.submit(task_func, payload): task_id
                    for task_id, payload in enumerate(task_data)
                }

                for future in as_completed(future_to_id):
                    yield future.result()

            return

        profiled_data = [
            (kind, task_id, payload)
            for task_id, payload in enumerate(task_data)
        ]

        if mode_exec == "map":
            for task_id, timing, result in executor.map(
                calc_sig_profile_task,
                profiled_data,
            ):
                received = time.perf_counter()
                timing["recv_s"] = received
                timing["recv_lag_s"] = received - timing["worker_end_s"]
                task_profiles[task_id].update(timing)

                yield result
        else:
            future_to_id = {
                executor.submit(calc_sig_profile_task, payload): payload[1]
                for payload in profiled_data
            }

            for order, future in enumerate(as_completed(future_to_id)):
                task_id, timing, result = future.result()
                received = time.perf_counter()
                timing["completion_order"] = order
                timing["recv_s"] = received
                timing["recv_lag_s"] = received - timing["worker_end_s"]
                task_profiles[task_id].update(timing)

                yield result


def summarize_profiles(profile: dict | None) -> None:
    """
    Add derived worker timing summaries to an active profile.
    """

    if profile is None:
        return

    worker_s = [
        task["worker_s"]
        for task in profile.get("tasks", [])
        if "worker_s" in task
    ]
    recv_lag_s = [
        task["recv_lag_s"]
        for task in profile.get("tasks", [])
        if "recv_lag_s" in task
    ]

    if worker_s:
        profile["worker_sum_s"] = sum(worker_s)
        profile["worker_max_s"] = max(worker_s)

    if recv_lag_s:
        profile["recv_lag_sum_s"] = sum(recv_lag_s)
        profile["recv_lag_max_s"] = max(recv_lag_s)


def apply_sig_adj(
    sig_bin: Any,
    frg_tot: int,
    is_norm: bool,
    scl_fct: float | None = None,
) -> None:
    """
    Apply global depth normalization and scaling to already merged signal.

    Parameters
    ----------
    sig_bin : Any
        Mutable merged sparse signal container.
    frg_tot : int
        Total accepted fragments used for normalization.
    is_norm : bool
        Whether to apply depth normalization.
    scl_fct : float | None
        Optional multiplicative scale factor.

    Raises
    ------
    ValueError
        If normalization is requested with no accepted fragments.
    """

    if is_norm:
        if frg_tot <= 0:
            raise ValueError(
                "Normalization requires non-zero total fragments.",
            )

        for values in iter_sig_arr(sig_bin):
            values /= frg_tot

    if scl_fct is not None:
        validate_comparison(scl_fct, "gt", 0, "scl_fct")

        for values in iter_sig_arr(sig_bin):
            values *= scl_fct


def is_sparse_sig(sig_bin: object) -> bool:
    """
    Return True for the merged sparse-array signal container.
    """

    return (
        isinstance(sig_bin, tuple)
        and len(sig_bin) == 2
        and sig_bin[0] == "sparse_merged"
    )


def iter_sig_arr(sig_bin: Any) -> Iterator[np.ndarray]:
    """
    Yield mutable value arrays from the merged sparse signal container.
    """

    if is_sparse_sig(sig_bin):
        for _, _, values in sig_bin[1]:
            yield values


def n_sig_rows(sig_bin: Any) -> int:
    """
    Count output bedGraph rows in the merged sparse signal container.
    """

    return sum(int(starts.size) for _, starts, _ in sig_bin[1])


def is_parts_ordered(starts_parts: list[np.ndarray]) -> bool:
    """
    Return True when sparse parts are individually and collectively ordered.
    """

    last_seen = None

    for starts in starts_parts:
        if starts.size == 0:
            continue

        if starts.size > 1 and np.any(starts[1:] < starts[:-1]):
            return False

        first = int(starts[0])
        last = int(starts[-1])
        if last_seen is not None and first < last_seen:
            return False

        last_seen = last

    return True


def coalesce_sparse(
    starts: np.ndarray,
    values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Sum adjacent duplicate starts in already ordered sparse arrays.
    """

    if starts.size == 0:
        return starts, values

    boundaries = np.empty(starts.size, dtype=bool)
    boundaries[0] = True
    boundaries[1:] = starts[1:] != starts[:-1]
    first_idx = np.flatnonzero(boundaries)
    summed = np.add.reduceat(values, first_idx)
    unique_starts = starts[first_idx]
    keep = summed != 0.0

    return (
        unique_starts[keep].astype(np.int64, copy=False),
        summed[keep].astype(np.float64, copy=False),
    )


def coalesce_parts(
    parts_chr: dict[str, list[tuple[np.ndarray, np.ndarray]]],
    optimize: bool = True,
    stats: dict[str, int] | None = None,
) -> object:
    """
    Sort and sum sparse signal arrays by chromosome without building a dict.

    Parameters
    ----------
    parts_chr : dict[str, list[tuple[np.ndarray, np.ndarray]]]
        Sparse index/value array parts grouped by chromosome.
    optimize : bool
        Whether to compact the final array representation.
    stats : dict[str, int] | None
        Optional collection receiving coalescing metrics.

    Returns
    -------
    signal : object
        Ordered compact signal arrays grouped by chromosome.
    """

    if stats is not None:
        stats.setdefault("single_chroms", 0)
        stats.setdefault("ordered_chroms", 0)
        stats.setdefault("sorted_chroms", 0)
        stats.setdefault("parts", 0)
        stats.setdefault("rows_before", 0)
        stats.setdefault("rows_after", 0)

    merged = []

    for chrom in sorted(parts_chr, key=sort_chrom):
        starts_parts = []
        values_parts = []

        for starts, values in parts_chr[chrom]:
            if starts.size == 0:
                continue

            starts_parts.append(starts.astype(np.int64, copy=False))
            values_parts.append(values.astype(np.float64, copy=False))

        if not starts_parts:
            continue

        if stats is not None:
            stats["parts"] += len(starts_parts)
            stats["rows_before"] += sum(
                int(starts.size) for starts in starts_parts
            )

        if optimize and len(starts_parts) == 1:
            starts = starts_parts[0]
            values = values_parts[0]
            keep = values != 0.0

            if np.any(keep):
                out_starts = starts[keep].astype(np.int64, copy=False)
                out_values = values[keep].astype(np.float64, copy=False)
                merged.append((chrom, out_starts, out_values))

                if stats is not None:
                    stats["single_chroms"] += 1
                    stats["rows_after"] += int(out_starts.size)

            continue

        is_ordered = optimize and is_parts_ordered(starts_parts)
        starts = np.concatenate(starts_parts)
        values = np.concatenate(values_parts)
        if starts.size == 0:
            continue

        if not is_ordered:
            order = np.argsort(starts, kind="mergesort")
            starts = starts[order]
            values = values[order]

            if stats is not None:
                stats["sorted_chroms"] += 1
        elif stats is not None:
            stats["ordered_chroms"] += 1

        out_starts, out_values = coalesce_sparse(starts, values)

        if out_starts.size > 0:
            merged.append((chrom, out_starts, out_values))

            if stats is not None:
                stats["rows_after"] += int(out_starts.size)

    return "sparse_merged", merged


def merge_sig(
    results: Iterable[Any],
    siz_bin: int,
    profile: dict[str, Any] | None = None,
) -> tuple[Any, int, int, int]:
    """
    Merge worker signal results and return signal, fragments, bins, and bytes.

    Parameters
    ----------
    results : Iterable[Any]
        Per-task pairs of a tagged sparse signal result and the fragment count
        that task accepted.
    siz_bin : int
        Positive signal-bin width in base pairs.
    profile : dict[str, Any] | None
        Optional mutable timing and merge-statistics payload.

    Returns
    -------
    combined_signal, frg_tot, n_bin, payload_bytes : tuple[
        Any, int, int, int
    ]
        Merged sparse signal container and execution statistics.
    """

    sparse_parts: dict[str, list[tuple[np.ndarray, np.ndarray]]] = defaultdict(
        list,
    )

    frg_tot = 0
    n_bin = 0
    payload_bytes = 0

    meta_s = 0.0
    accum_s = 0.0

    stats_coal: dict[str, int] = {}

    for result in results:
        time_meta = time.perf_counter()
        signal_result, n_frg = result
        frg_tot += n_frg
        n_bin += count_result_bins(signal_result)
        payload_bytes += est_result_bytes(signal_result)

        if profile is not None:
            meta_s += time.perf_counter() - time_meta

        time_accum = time.perf_counter()

        for chrom, starts, values in signal_result[1]:
            sparse_parts[chrom].append(
                (starts.astype(np.int64, copy=False), values),
            )

        if profile is not None:
            accum_s += time.perf_counter() - time_accum

    time_coal = time.perf_counter()
    combined_signal = coalesce_parts(
        sparse_parts,
        optimize=True,
        stats=stats_coal,
    )

    if profile is not None:
        phases = profile["phases_s"]
        phases["merge_result_metadata"] = meta_s
        phases["merge_accumulate_results"] = accum_s
        phases["merge_array_coalesce"] = time.perf_counter() - time_coal
        profile["stats_coal"] = stats_coal

    return (
        combined_signal,
        frg_tot,
        n_bin,
        payload_bytes,
    )


def format_bdg_value(value: float, dp: int) -> str:
    """
    Format one bedGraph value with the same rounding used by write_bdg().
    """

    rounded_value = round(float(value), dp)

    if rounded_value == 0.0:
        rounded_value = 0.0

    rendered_value = f"{rounded_value:.{dp}f}"

    if "." in rendered_value:
        rendered_value = rendered_value.rstrip("0").rstrip(".")

    if rendered_value == "-0":
        rendered_value = "0"

    return rendered_value


def iter_bed_rows(
    results: Iterable[object],
) -> Iterator[tuple[str, int, int, int]]:
    """
    Yield deterministic BED rows from compact worker array results.

    Parameters
    ----------
    results : Iterable[object]
        Compact array-backed worker results.

    Yields
    ------
    row : tuple[str, int, int, int]
        Chromosome, start, end, and fragment length in coordinate order.
    """

    rows_chr: dict[str, list[tuple[np.ndarray, np.ndarray, np.ndarray]]] = (
        defaultdict(list)
    )

    for result in results:
        if not isinstance(result, tuple) or len(result) < 2:
            continue

        for chrom, starts, ends, lengths in result[1]:
            if starts.size:
                rows_chr[chrom].append((starts, ends, lengths))

    for chrom in sorted(rows_chr, key=sort_chrom):
        starts = np.concatenate([part[0] for part in rows_chr[chrom]])
        ends = np.concatenate([part[1] for part in rows_chr[chrom]])
        lengths = np.concatenate([part[2] for part in rows_chr[chrom]])
        order = np.lexsort((lengths, ends, starts))

        for idx in order:
            yield (chrom, int(starts[idx]), int(ends[idx]), int(lengths[idx]))


def write_bed_results(
    results: Iterable[object],
    fil_out: str,
) -> int:
    """
    Write compact BED worker results and return the emitted row count.
    """

    row_count = 0

    with open_out(fil_out) as bed_file:
        for chrom, start, end, length in iter_bed_rows(results):
            bed_file.write(f"{chrom}\t{start}\t{end}\t{length}\n")
            row_count += 1

    return row_count


def write_bed_frg(
    frg_tup: dict[str, list[tuple[int, int, int]]],
    fil_out: str,
) -> int:
    """
    Write chromosome-grouped fragment tuples in deterministic BED order.
    """

    row_count = 0

    with open_out(fil_out) as bed_file:
        for chrom in sorted(frg_tup.keys(), key=sort_chrom):
            for start, end, length in sorted(
                frg_tup[chrom],
                key=lambda t: t[0],
            ):
                bed_file.write(f"{chrom}\t{start}\t{end}\t{length}\n")
                row_count += 1

    return row_count


def format_bdg_rows(
    chrom: str,
    starts: np.ndarray,
    values: np.ndarray,
    siz_bin: int,
    dp: int,
    chrom_size: int | None,
) -> str:
    """
    Format one sparse-array bedGraph shard.
    """

    lines = []

    for start, value in zip(starts, values, strict=True):
        bin_start = int(start)
        bin_end = bin_start + siz_bin

        if chrom_size is not None:
            if bin_start < 0 or bin_start >= chrom_size:
                raise ValueError(
                    "bedGraph bin start is outside chromosome bounds: "
                    f"{chrom}:{bin_start} (chromosome size {chrom_size}).",
                )

            bin_end = min(bin_end, chrom_size)

        lines.append(
            f"{chrom}\t{bin_start}\t{bin_end}\t"
            f"{format_bdg_value(float(value), dp)}\n",
        )

    return "".join(lines)


def format_bdg_rows_task(data: tuple[object, ...]) -> str:
    """
    ProcessPool-friendly wrapper for sparse-array bedGraph formatting.
    """

    return format_bdg_rows(*data)


def iter_sig_chunks(
    sig_bin: Any,
    n_chunks: int,
) -> Iterator[tuple[str, np.ndarray, np.ndarray]]:
    """
    Yield ordered sparse signal chunks for parallel formatting.

    Parameters
    ----------
    sig_bin : Any
        Merged sparse signal container.
    n_chunks : int
        Approximate number of ordered output chunks.

    Yields
    ------
    chunk : tuple[str, np.ndarray, np.ndarray]
        Chromosome with matching bin-start and value arrays.
    """

    total_rows = max(n_sig_rows(sig_bin), 1)
    chunk_size = max(1, math.ceil(total_rows / max(n_chunks, 1)))

    for chrom, starts, values in sig_bin[1]:
        for i in range(0, starts.size, chunk_size):
            yield (
                chrom,
                starts[i : i + chunk_size],
                values[i : i + chunk_size],
            )


def write_bdg_sparse(
    sig_bin: Any,
    fil_out: str,
    siz_bin: int,
    dp: int,
    siz_chr: dict[str, int] | None,
) -> None:
    """
    Write a sparse-array signal directly to bedGraph.
    """

    with open_out(fil_out) as handle:
        for chrom, starts, values in iter_sig_chunks(sig_bin, 1):
            handle.write(
                format_bdg_rows(
                    chrom,
                    starts,
                    values,
                    siz_bin,
                    dp,
                    siz_chr.get(chrom) if siz_chr is not None else None,
                ),
            )


def write_bdg_par(
    sig_bin: Any,
    fil_out: str,
    siz_bin: int,
    dp: int,
    siz_chr: dict[str, int] | None,
    workers: int,
) -> None:
    """
    Format sorted bedGraph shards in worker processes and write in order.

    Parameters
    ----------
    sig_bin : Any
        Merged sparse signal container.
    fil_out : str
        Output bedGraph path.
    siz_bin : int
        Signal-bin width.
    dp : int
        Decimal precision for rendered values.
    siz_chr : dict[str, int] | None
        Optional chromosome sizes used to clip final intervals.
    workers : int
        Number of formatting worker processes.
    """

    validate_comparison(
        workers,
        "ge",
        1,
        "wrk_writer",
        allow_none=False,
    )

    if workers == 1:
        write_bdg_sparse(
            sig_bin,
            fil_out,
            siz_bin,
            dp,
            siz_chr,
        )

        return

    task_data = [
        (
            chrom,
            starts,
            values,
            siz_bin,
            dp,
            siz_chr.get(chrom) if siz_chr is not None else None,
        )
        for chrom, starts, values in iter_sig_chunks(
            sig_bin,
            workers,
        )
    ]

    try:
        with ProcessPoolExecutor(
            max_workers=min(workers, os.cpu_count() or 1),
        ) as ex:
            formatted = list(ex.map(format_bdg_rows_task, task_data))
    except OSError:
        write_bdg_sparse(
            sig_bin,
            fil_out,
            siz_bin,
            dp,
            siz_chr,
        )

        return

    with open_out(fil_out) as handle:
        for shard in formatted:
            handle.write(shard)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse arguments for BAM or CRAM signal computation.

    Parameters
    ----------
    argv : list[str] | None
        Optional argument vector to parse. If None, use 'sys.argv[1:]'.

    Returns
    -------
    arguments : argparse.Namespace
        Parsed command-line arguments.

    Raises
    ------
    SystemExit
        Raised by argparse when '--help' is shown or when required arguments or
        valid choices are missing.

    Notes
    -----
    Parser option names, defaults, and accepted choices are documented in the
    'add_argument()' definitions below and in rendered '--help' output.
    """

    parser = CapArgumentParser(
        description=(
            "Compute binned signal from a BAM or CRAM file in bedGraph "
            "format, optionally applying normalization.\n"
            "\n"
            "Alternatively, extract and output processed fragment coordinates "
            "in a BED-like format, which can be used as input to the original "
            "implementation of siQ-ChIP "
            "(https://github.com/BradleyDickson/siQ-ChIP) or one updated for "
            "use with S. cerevisiae data "
            "(https://github.com/kalavattam/siQ-ChIP/tree/protocol).\n"
            "\n"
            "The fragment count 'N' and the spanned-bin count 'L' can be "
            "written alongside either output or on their own: when "
            "'--fil_out' is omitted, the run counts and writes only the "
            "requested reports. Both are inputs to 'compute_pseudo', where "
            "'k = L / N' puts the pseudocount in bins, the unit per-bin "
            "coverage is added in."
        ),
    )
    add_help_cap(parser)

    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Run script in verbose mode (stderr banner of parsed args).\n\n",
    )
    parser.add_argument(
        "-t",
        "--threads",
        dest="threads",
        type=int,
        default=1,
        help=(
            "Number of threads to use for parallel processing (>= 1; default: "
            "%(default)s).\n"
            "\n"
            "When '--threads > 1', different chromosomes are processed in "
            "parallel.\n"
            "\n"
        ),
    )

    parser.add_argument(
        "-fi",
        "--fil_in",
        dest="fil_in",
        default=None,
        help="Input file path for the BAM or CRAM file.\n\n",
    )
    parser.add_argument(
        "--fil-in",
        dest="fil_in",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-rf",
        "--ref_fa",
        dest="ref_fa",
        default=None,
        help=(
            "Reference FASTA file for CRAM decoding. Required for CRAM inputs "
            "unless the reference is otherwise available to htslib/pysam.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "--ref-fa",
        dest="ref_fa",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-fo",
        "--fil_out",
        dest="fil_out",
        required=False,
        default=None,
        help=(
            "Output file path. Required unless '--report_n_frg' or "
            "'--report_n_bin' is given, in which case it may be omitted to "
            "count without writing a track.\n"
            "\n"
            "Supported output types are bedGraph ('bedGraph', 'bedgraph', "
            "'bdg', 'bg') and BED ('bed').\n"
            "\n"
            "Append '.gz' for gzip compression, e.g., 'output.bdg.gz'.\n"
            "\n"
            "Note: requesting BED output causes the script to write processed "
            "fragment coordinates in a BED-like format, and '--method', "
            "'--scl_fct', and '--dp' are ignored. '--siz_bin' is ignored too, "
            "unless '--report_n_bin' is given, since 'L' counts bins.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "--fil-out",
        dest="fil_out",
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "-me",
        "--method",
        dest="method",
        choices=METHOD_CHOICES,
        default="norm",
        help=(
            "Workflow method. Specify signal calculation type (default: "
            "'%(default)s').\n"
            "  - Unadjusted aliases: 'r', 'raw', 'u', 'unadj', 'unadjusted', "
            "'s', 'smp', 'simple'. Internally standardized to 'unadj'.\n"
            "  - Fragment-length-normalized aliases: 'f', 'frg', 'frag', "
            "'frg_len', 'frag_len', 'l', 'len', 'len_frg', 'len_frag'. "
            "Internally standardized to 'frag'.\n"
            "  - Normalized-coverage aliases: 'n', 'nc', 'nrm', 'norm', "
            "'normalized'. Internally standardized to 'norm'.\n"
            "\n"
            "Note: 'norm' normalizes for both fragment length and total "
            "fragment count so that the genome-wide summed signal is "
            "approximately 1.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-sb",
        "--siz_bin",
        dest="siz_bin",
        type=int,
        default=10,
        help=(
            "Bin size in base pairs for signal calculation (default: "
            "%(default)s).\n"
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
        "-eg",
        "--engine",
        dest="engine",
        choices=ENGINE_CHOICES,
        default="chrom",
        help=(
            "Processing engine for bedGraph output (default: %(default)s).\n"
            "  - 'chrom' parallelizes indexed chromosome fetching and "
            "array-backed signal calculation.\n"
            "  - 'window' parallelizes indexed coordinate-window fetching for "
            "finer load balance on, e.g., large BAM inputs.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "-sw",
        "--siz_win",
        dest="siz_win",
        type=int,
        default=100000,
        help=(
            "Window size in base pairs for '--engine window' indexed fetch "
            "tasks (default: %(default)s). Ignored by '--engine chrom'.\n"
            "\n"
            "Each chromosome is split into windows of this size, and one "
            "worker task fetches and bins each window. Smaller windows give "
            "finer load balance across threads at the cost of more fetch "
            "overhead.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "--siz-win",
        dest="siz_win",
        type=int,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-sf",
        "--scl_fct",
        dest="scl_fct",
        type=float,
        default=None,
        help=(
            "Scaling factor to apply to the signal (default: %(default)s).\n\n"
        ),
    )
    parser.add_argument(
        "--scl-fct",
        dest="scl_fct",
        type=float,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-uf",
        "--usr_frg",
        dest="usr_frg",
        type=int,
        default=None,
        help=(
            "Fixed fragment length to use instead of read lengths (single-end "
            "alignments) or template lengths (paired-end alignments; default: "
            "%(default)s).\n"
            "\n"
        ),
    )
    parser.add_argument(
        "--usr-frg",
        dest="usr_frg",
        type=int,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-rnf",
        "--report_n_frg",
        dest="report_n_frg",
        nargs="?",
        default=None,
        const=REPORT_DERIVE,
        help=(
            "Write the fragment count 'N' to this file as a single integer on "
            "one line (default: %(default)s).\n"
            "\n"
            "Given without a path, the file is written beside '--fil_out', "
            "with its extension (and any '.gz') replaced by '.n_frg.txt'. If "
            "'--fil_out' is not used, then an output path must be given.\n"
            "\n"
            "The count covers the fragments the signal path uses, taken from "
            "the same iterator under the same '--usr_frg' and filter "
            "settings, so 'N' is identical for every '--method' as well as in "
            '"report-only mode".\n'
            "\n"
            "For '--method norm', 'N' is the denominator the signal is "
            "divided by, and a '--method frag' track's value column sums to "
            "'N'. Other tracks do not, and 'N' is counted from the fragments "
            "rather than summed off any track.\n"
            "\n"
        ),
    )
    parser.add_argument(
        "--report-n-frg",
        dest="report_n_frg",
        nargs="?",
        const=REPORT_DERIVE,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-rnb",
        "--report_n_bin",
        dest="report_n_bin",
        nargs="?",
        default=None,
        const=REPORT_DERIVE,
        help=(
            "Write the spanned-bin count 'L' to this file as a single integer "
            "on one line (default: %(default)s).\n"
            "\n"
            "Given without a path, the file is written beside '--fil_out', "
            "with its extension (and any '.gz') replaced by '.n_bin.txt'. If "
            "'--fil_out' is not used, then an output path must be given.\n"
            "\n"
            "Depends on '--siz_bin' and '--usr_frg', which set how each "
            "fragment is divided and how far it reaches, respectively.\n"
            "\n"
            "Bins are counted once per fragment that touches them, so 'L' "
            "totals the bins spanned by the same fragments counted by 'N' "
            "(see '--report_n_frg' above).\n"
            "\n"
        ),
    )
    parser.add_argument(
        "--report-n-bin",
        dest="report_n_bin",
        nargs="?",
        const=REPORT_DERIVE,
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
            "values (default: %(default)s). After rounding, non-informative "
            "trailing zeros are stripped.\n"
            "\n"
        ),
    )

    parser.add_argument(
        "-pj",
        "--profile_json",
        "--profile-json",
        dest="profile_json",
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-mx",
        "--mode_exec",
        "--mode-exec",
        dest="mode_exec",
        choices=MODE_EXEC_CHOICES,
        default="map",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-bs",
        "--strat_bed",
        "--strat-bed",
        dest="strat_bed",
        choices=STRAT_BED_CHOICES,
        default="auto",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-ws",
        "--strat_writer",
        "--strat-writer",
        dest="strat_writer",
        choices=STRAT_WRITER_CHOICES,
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-ww",
        "--wrk_writer",
        "--wrk-writer",
        dest="wrk_writer",
        type=int,
        default=None,
        help=argparse.SUPPRESS,
    )

    argv_parse = sys.argv[1:] if argv is None else argv

    if not argv_parse:
        parser.print_help(sys.stderr)
        raise SystemExit(0)

    return parser.parse_args(argv_parse)


def _write_bed(
    args: argparse.Namespace,
    fil_out: str,
    siz_chr_hdr: dict[str, int],
    siz_chr: dict[str, int],
    profile: dict | None,
    time_tot: float,
) -> int:
    """
    Compute and write fragment-coordinate BED output.

    Parameters
    ----------
    args : argparse.Namespace
        Validated command-line arguments.
    fil_out : str
        Writable BED output path.
    siz_chr_hdr : dict[str, int]
        Alignment-header chromosome lengths in source order.
    siz_chr : dict[str, int]
        Resolved chromosome lengths selected for output.
    profile : dict | None
        Optional mutable execution profile.
    time_tot : float
        Monotonic start time for total-duration reporting.

    Returns
    -------
    status : int
        Zero after serial or indexed BED output succeeds.

    Raises
    ------
    OSError
        If an explicitly requested indexed strategy cannot read the alignment.
    """

    strat_bed = args.strat_bed

    if strat_bed == "auto":
        if args.threads > 1 and has_idx_aln(
            args.fil_in,
            args.ref_fa,
        ):
            strat_bed = "idx_chrom"
        else:
            strat_bed = "serial"

    if profile is not None:
        profile["strat_bed_resolved"] = strat_bed

    if strat_bed == "serial":
        time_phase = time.perf_counter()
        fragments = parse_bam(
            args.fil_in,
            args.usr_frg,
            args.ref_fa,
            siz_chr=siz_chr,
        )
        record_phase(profile, "parse_bed_fragments", time_phase)

        time_phase = time.perf_counter()
        row_count = write_bed_frg(fragments, fil_out)
        record_phase(profile, "write_bed", time_phase)

        if profile is not None:
            profile["n_rows_out"] = row_count
    else:
        time_phase = time.perf_counter()
        check_idx_aln(args.fil_in, args.ref_fa)
        record_phase(
            profile,
            "validate_indexed_bed_alignment",
            time_phase,
        )

        time_phase = time.perf_counter()
        task_chr = [chrom for chrom in siz_chr_hdr if chrom in siz_chr]

        if strat_bed == "idx_chrom":
            tasks = [
                (
                    args.fil_in,
                    chrom,
                    siz_chr,
                    args.usr_frg,
                    args.ref_fa,
                    None,
                    None,
                )
                for chrom in task_chr
            ]
        else:
            tasks = [
                (
                    args.fil_in,
                    chrom,
                    siz_chr,
                    args.usr_frg,
                    args.ref_fa,
                    start,
                    min(
                        start + args.siz_win,
                        siz_chr[chrom],
                    ),
                )
                for chrom in task_chr
                for start in range(
                    0,
                    siz_chr[chrom],
                    args.siz_win,
                )
            ]

        if profile is not None:
            profile["tasks"] = [
                {
                    "task_id": task_id,
                    "kind": "indexed_bed",
                    "strat_bed": strat_bed,
                    "chrom": task[1],
                    "start": task[5],
                    "end": task[6],
                    "chrom_size": int(siz_chr[task[1]]),
                }
                for task_id, task in enumerate(tasks)
            ]
            profile["n_tasks"] = len(tasks)

        record_phase(profile, "build_indexed_bed_tasks", time_phase)
        time_phase = time.perf_counter()

        try:
            results = list(
                iter_task_results(
                    "indexed_bed",
                    collect_bed_idx_task,
                    tasks,
                    args.threads,
                    args.mode_exec,
                    profile["tasks"] if profile is not None else None,
                ),
            )
        except OSError:
            if args.strat_bed != "auto":
                raise

            if profile is not None:
                profile["strat_bed_resolved"] = "serial_fallback"

            fragments = parse_bam(
                args.fil_in,
                args.usr_frg,
                args.ref_fa,
                siz_chr=siz_chr,
            )
            record_phase(profile, "parse_bed_fragments", time_phase)

            time_phase = time.perf_counter()
            row_count = write_bed_frg(fragments, fil_out)
            record_phase(profile, "write_bed", time_phase)

            if profile is not None:
                profile["n_rows_out"] = row_count
                profile["phases_s"]["total"] = time.perf_counter() - time_tot
                write_profile(args.profile_json, profile)

            return 0

        record_phase(
            profile,
            "receive_indexed_bed_results",
            time_phase,
        )
        summarize_profiles(profile)

        if profile is not None:
            profile["result_payload_bytes"] = sum(
                est_bed_bytes(result) for result in results
            )

        time_phase = time.perf_counter()
        row_count = write_bed_results(results, fil_out)
        record_phase(
            profile,
            "merge_sort_write_indexed_bed",
            time_phase,
        )

        if profile is not None:
            profile["n_rows_out"] = row_count

    if profile is not None:
        profile["phases_s"]["total"] = time.perf_counter() - time_tot
        write_profile(args.profile_json, profile)

    return 0


def _write_bdg(
    args: argparse.Namespace,
    sig_bin: Any,
    fil_out: str,
    siz_chr: dict[str, int],
    profile: dict | None,
    time_tot: float,
) -> int:
    """
    Write the merged bedGraph signal.

    Parameters
    ----------
    args : argparse.Namespace
        Validated command-line arguments and output strategy.
    sig_bin : Any
        Merged sparse signal container.
    fil_out : str
        Writable bedGraph output path.
    siz_chr : dict[str, int]
        Resolved chromosome lengths keyed by chromosome.
    profile : dict | None
        Optional mutable execution profile.
    time_tot : float
        Monotonic start time for total-duration reporting.

    Returns
    -------
    status : int
        Zero after bedGraph output succeeds.
    """

    time_phase = time.perf_counter()

    if args.strat_writer == "parallel_ordered":
        write_bdg_par(
            sig_bin,
            fil_out,
            args.siz_bin,
            args.dp,
            siz_chr,
            args.wrk_writer,
        )
    else:
        write_bdg_sparse(
            sig_bin,
            fil_out,
            args.siz_bin,
            args.dp,
            siz_chr,
        )

    if profile is not None:
        profile["output_written"] = True
        profile["n_rows_out"] = n_sig_rows(sig_bin)

    record_phase(profile, "write_bedgraph", time_phase)

    if profile is not None:
        profile["phases_s"]["total"] = time.perf_counter() - time_tot
        write_profile(args.profile_json, profile)

    return 0


def _sig_results(
    args: argparse.Namespace,
    siz_chr_hdr: dict[str, int],
    siz_chr: dict[str, int],
    is_len: bool,
    is_norm: bool,
    profile: dict | None,
) -> Iterable[Any]:
    """
    Build indexed worker tasks for a public signal engine.

    Parameters
    ----------
    args : argparse.Namespace
        Validated command-line arguments and public engine selection.
    siz_chr_hdr : dict[str, int]
        Alignment-header chromosome lengths in source order.
    siz_chr : dict[str, int]
        Resolved chromosome lengths selected for signal computation.
    is_len : bool
        Whether signal weights use fragment length.
    is_norm : bool
        Whether signal is normalized by total fragments.
    profile : dict | None
        Optional mutable task and timing profile.

    Returns
    -------
    results : Iterable[Any]
        Ordered worker pairs of a tagged sparse signal result and the fragment
        count that task accepted.

    Raises
    ------
    ValueError
        If the alignment lacks the index the engines require.
    """

    strategy = ENGINE_STRAT[args.engine]

    time_phase = time.perf_counter()
    check_idx_aln(args.fil_in, args.ref_fa)
    record_phase(profile, "validate_indexed_alignment", time_phase)

    time_phase = time.perf_counter()
    task_chr = [chrom for chrom in siz_chr_hdr if chrom in siz_chr]

    if strategy == "idx_chrom":
        tasks = [
            (
                args.fil_in,
                chrom,
                siz_chr,
                args.usr_frg,
                args.ref_fa,
                None,
                None,
                args.siz_bin,
                is_len or is_norm,
            )
            for chrom in task_chr
        ]
    else:
        tasks = [
            (
                args.fil_in,
                chrom,
                siz_chr,
                args.usr_frg,
                args.ref_fa,
                start,
                min(
                    start + args.siz_win,
                    siz_chr[chrom],
                ),
                args.siz_bin,
                is_len or is_norm,
            )
            for chrom in task_chr
            for start in range(
                0,
                siz_chr[chrom],
                args.siz_win,
            )
        ]

    if profile is not None:
        profile["fragments_total"] = None
        profile["tasks"] = [
            {
                "task_id": task_id,
                "kind": strategy,
                "chrom": task[1],
                "start": task[5],
                "end": task[6],
                "chrom_size": int(siz_chr[task[1]]),
                "n_bins": math.ceil(siz_chr[task[1]] / args.siz_bin),
            }
            for task_id, task in enumerate(tasks)
        ]
        profile["n_tasks"] = len(tasks)

    record_phase(profile, "build_tasks", time_phase)

    return iter_task_results(
        strategy,
        calc_sig_idx_fetch_task,
        tasks,
        args.threads,
        args.mode_exec,
        profile["tasks"] if profile is not None else None,
    )


def main(argv: list[str] | None = None) -> int:
    """
    Execute the primary control flow for the script.

    Parameters
    ----------
    argv : list[str] | None
        Optional list of command-line arguments. When None (the default),
        'sys.argv[1:]' is used.

    Returns
    -------
    status : int
        Zero on success after writing binned bedGraph signal or processed
        fragment coordinates to '--fil_out'.

    Raises
    ------
    SystemExit
        For help, invalid dimensions or scaling, unsupported output paths,
        zero-fragment normalization, or BAM/CRAM read failures.

    Notes
    -----
    Verbose argument banners and human-readable failure diagnostics are written
    to stderr.
    """

    time_tot = time.perf_counter()

    args = parse_args(argv)

    if args.fil_in == "-":
        raise SystemExit(
            "'--fil_in -' is not supported; provide a BAM or CRAM path.",
        )

    if args.fil_out == "-":
        raise SystemExit(
            "'--fil_out -' is not supported; provide an output path.",
        )

    # Neither rule fits argparse: '--fil_out' is optional only in report-only
    # mode, and a hidden hyphen spelling is a separate action that argparse's
    # own 'required' check cannot see. Check the parsed values here instead.
    if args.fil_in is None:
        raise SystemExit(
            "'--fil_in' is required. Supply the BAM or CRAM input path.",
        )

    report_only = args.fil_out is None

    if report_only and args.report_n_frg is None and args.report_n_bin is None:
        raise SystemExit(
            "'--fil_out' is required unless '--report_n_frg' or "
            "'--report_n_bin' is given. Supply an output path to write a "
            "track, or a report path to count without writing one.",
        )

    try:
        check_exists(args.fil_in, kind="file", label="Alignment file")

        if args.ref_fa is not None:
            check_exists(args.ref_fa, kind="file", label="Reference FASTA")

    except FileNotFoundError as error:
        raise SystemExit(str(error)) from None

    try:
        if report_only:
            fil_out, fmt_out = None, None
        else:
            fil_out, fmt_out, _ = validate_output_path(
                args.fil_out,
                ALLOWED_OUTPUT_FORMATS,
            )

            check_writable(fil_out, "file")

        if args.profile_json is not None:
            check_writable(args.profile_json, "file")

        # Resolve before the writability check below, so a path the user left
        # to derivation is validated exactly like one they spelled out.
        args.report_n_frg = resolve_report_path(
            args.report_n_frg,
            fil_out,
            "n_frg",
            "--report_n_frg",
        )
        args.report_n_bin = resolve_report_path(
            args.report_n_bin,
            fil_out,
            "n_bin",
            "--report_n_bin",
        )

        for path_report in (args.report_n_frg, args.report_n_bin):
            if path_report is not None:
                check_writable(path_report, "file")
    except (
        ValueError,
        FileNotFoundError,
        PermissionError,
        IsADirectoryError,
    ) as error:
        raise SystemExit(str(error)) from None

    if fmt_out != "bed":
        if args.strat_writer is None:
            args.strat_writer = "parallel_ordered"

        if args.wrk_writer is None:
            args.wrk_writer = (
                1
                if args.strat_writer == "serial"
                else min(4, os.cpu_count() or 1)
            )
    else:
        if args.strat_writer is None:
            args.strat_writer = "serial"

        if args.wrk_writer is None:
            args.wrk_writer = 1

    try:
        validate_comparison(args.threads, "ge", 1, "threads", allow_none=False)

        if fmt_out != "bed":
            validate_comparison(
                args.siz_bin,
                "gt",
                0,
                "siz_bin",
                allow_none=False,
            )
            validate_comparison(
                args.siz_win,
                "gt",
                0,
                "siz_win",
                allow_none=False,
            )
            validate_comparison(
                args.wrk_writer,
                "ge",
                1,
                "wrk_writer",
                allow_none=False,
            )
            validate_comparison(
                args.wrk_writer,
                "le",
                4,
                "wrk_writer",
                allow_none=False,
            )
            validate_comparison(args.dp, "ge", 0, "dp", allow_none=False)

            if args.strat_writer == "serial" and args.wrk_writer != 1:
                raise ValueError(
                    "'--wrk_writer' must be 1 when '--strat_writer serial'.",
                )
        else:
            validate_comparison(
                args.siz_win,
                "gt",
                0,
                "siz_win",
                allow_none=False,
            )

        validate_comparison(args.scl_fct, "gt", 0, "scl_fct", allow_none=True)
        validate_comparison(
            args.usr_frg,
            "gt",
            0,
            "usr_frg",
            allow_none=True,
        )

    except ValueError as error:
        raise SystemExit(str(error)) from None

    mthd_in = args.method
    args.method = METHOD_CANON[args.method]

    is_len = args.method == "frag"
    is_norm = args.method == "norm"
    profile = start_profile(args, fil_out, fmt_out)

    if args.verbose:
        with redirect_stdout(sys.stderr):
            print("")
            print("####################################")
            print("## Arguments for 'compute_signal' ##")
            print("####################################")
            print("")
            print("--verbose")
            print(f"--threads  {args.threads}")
            print(f"--fil_in   {args.fil_in}")
            print(f"--ref_fa   {args.ref_fa}")
            print(f"--fil_out  {fil_out}")

            if fmt_out == "bed":
                print(f"--usr_frg  {args.usr_frg}")
                print(f"--report_n_frg {args.report_n_frg}")
                print(f"--report_n_bin {args.report_n_bin}")
                print(
                    "\n\n(BED output mode: signal computation arguments "
                    "ignored)\n",
                )
            else:
                if mthd_in != args.method:
                    mthd_msg = (
                        f"{mthd_in}  (standardized internally to "
                        f"{args.method})"
                    )

                    print(f"--method   {mthd_msg}")
                else:
                    print(f"--method   {args.method}")

                print(f"--siz_bin  {args.siz_bin}")
                print(f"--engine   {args.engine}")
                print(f"--siz_win  {args.siz_win}")
                print(f"--scl_fct  {args.scl_fct}")
                print(f"--usr_frg  {args.usr_frg}")
                print(f"--report_n_frg {args.report_n_frg}")
                print(f"--report_n_bin {args.report_n_bin}")
                print(f"--dp       {args.dp}")

            print("")
            print("")

    try:
        time_phase = time.perf_counter()
        siz_chr_hdr = get_siz_chr(args.fil_in, args.ref_fa)
        siz_chr = resolve_siz_chr(siz_chr_hdr)
        record_phase(profile, "resolve_siz_chr", time_phase)

        if profile is not None:
            profile["n_chrom_sizes"] = len(siz_chr)

        if args.report_n_frg is not None or args.report_n_bin is not None:
            time_phase = time.perf_counter()

            # Counted from the same iterator the signal path consumes, with the
            # same '--usr_frg', so 'N' here is the very number a
            # '--method norm' run divides by rather than an estimate of it.
            n_frg, n_bin = count_frgs_and_bins(
                iter_aln_frg(
                    fil_aln=args.fil_in,
                    siz_chr=siz_chr,
                    usr_frg=args.usr_frg,
                    ref_fa=args.ref_fa,
                ),
                args.siz_bin,
            )

            record_phase(profile, "count_frgs_and_bins", time_phase)

            if profile is not None:
                profile["report_n_frg"] = n_frg
                profile["report_n_bin"] = n_bin

            for path_report, value, label in (
                (args.report_n_frg, n_frg, "N"),
                (args.report_n_bin, n_bin, "L"),
            ):
                if path_report is None:
                    continue

                with open(path_report, "w") as handle:
                    handle.write(f"{value}\n")

                if args.verbose:
                    print(
                        f"{label} {value} -> {path_report}",
                        file=sys.stderr,
                    )

            if report_only:
                if args.profile_json is not None:
                    write_profile(args.profile_json, profile)

                return 0
    except FileNotFoundError:
        raise SystemExit(f"Alignment file not found: {args.fil_in}") from None
    except ValueError as error:
        raise SystemExit(str(error)) from None
    except OSError as error:
        raise SystemExit(
            f"I/O error while reading BAM or CRAM: {error}",
        ) from None

    try:
        if fmt_out == "bed":
            _write_bed(
                args,
                fil_out,
                siz_chr_hdr,
                siz_chr,
                profile,
                time_tot,
            )

            return 0

        # Otherwise, compute and write bedGraph signal. The two public engines
        # are 'chrom' and 'window', and both map to indexed fetch tasks.
        results = _sig_results(
            args,
            siz_chr_hdr,
            siz_chr,
            is_len,
            is_norm,
            profile,
        )

        # Receive worker results, then merge in the parent process.
        time_collect_merge = time.perf_counter()
        time_phase = time.perf_counter()
        result_list = list(results)
        record_phase(profile, "receive_worker_results", time_phase)
        summarize_profiles(profile)

        time_phase = time.perf_counter()
        (
            combined_signal,
            frg_tot,
            n_bin_pre,
            payload_bytes,
        ) = merge_sig(
            results=result_list,
            siz_bin=args.siz_bin,
            profile=profile,
        )

        record_phase(profile, "parent_merge_results", time_phase)

        time_phase = time.perf_counter()
        apply_sig_adj(
            combined_signal,
            frg_tot,
            is_norm,
            args.scl_fct,
        )
        record_phase(profile, "apply_sig_adj", time_phase)

        record_phase(profile, "collect_and_merge_results", time_collect_merge)

        if profile is not None:
            profile["fragments_total"] = frg_tot
            profile["result_bins_before_merge"] = n_bin_pre
            profile["result_bins_after_merge"] = n_sig_rows(
                combined_signal,
            )
            profile["result_payload_bytes"] = payload_bytes

        _write_bdg(
            args,
            combined_signal,
            fil_out,
            siz_chr,
            profile,
            time_tot,
        )

        return 0

    except ValueError as error:
        raise SystemExit(str(error)) from None

    except OSError as error:
        raise SystemExit(f"I/O error: {error}") from None


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        with suppress(Exception):
            sys.stdout.close()

        with suppress(Exception):
            sys.stderr.close()

        raise SystemExit(0) from None
