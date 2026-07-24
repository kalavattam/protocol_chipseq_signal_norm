#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: compute_signal.py
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in development and documentation.
#
# Distributed under the MIT license.


"""
Calculate binned signal or fragment-coordinate output from BAM/CRAM input.

Usage
-----
python -m protocol_chipseq_signal_norm.cli.compute_signal [options] --fil_in <file> --fil_out <file>

Parameters
----------
Input, output, signal-mode, engine, scaling, and formatting options are parsed
by 'parse_args()'.

Returns
-------
Writes bedGraph-like signal tracks or BED-like fragment-coordinate records.
When writing bedGraph, finite values are rounded to at most '--dp' decimal
places and trailing zeros are stripped.

See Also
--------
docs/dev/compute_signal.md
    Maintainer notes on signal modes, engines, coordinate assumptions, and
    performance.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import redirect_stdout, suppress
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import argparse
import hashlib
import json
import math
import os
import signal
import time
from collections.abc import Iterator

import numpy as np
import pysam

from protocol_chipseq_signal_norm.utilities.utils_bdg import (
    key_bin,
    load_chr_sizes,
    write_bdg,
)
from protocol_chipseq_signal_norm.utilities.utils_check import (
    ALLOWED,
    check_cmp,
    check_exists,
    check_parse_fil_out,
    check_writable,
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

# TODO: Reassess chunked processing, compiled-core options, and downstream
# bigWig compatibility after the current indexed engines are benchmarked.

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
ENGINE_CHOICES = ("chrom", "window")
PUBLIC_ENGINE_STRATEGY = {
    "chrom": "indexed_chrom",
    "window": "indexed_window",
}
EXECUTOR_MODE_CHOICES = ("map", "as_completed")
PROTOTYPE_PARSE_STRATEGY_CHOICES = (
    "serial",
    "indexed_chrom",
    "indexed_window"
)
PROTOTYPE_BED_STRATEGY_CHOICES = (
    "auto",
    "serial",
    "indexed_chrom",
    "indexed_window"
)
PROTOTYPE_RESULT_FORMAT_CHOICES = (
    "dict",
    "sparse_np",
    "dense_np",
    "direct_sparse_np",
    "direct_sparse_idx_np",
    "direct_sparse_local_np",
    "direct_sparse_touched_np",
    "direct_sparse_bincount_np",
    "direct_sparse_touched_bincount_np",
    "direct_sparse_local_bincount_np",
    "direct_dense_np",
    "event_np"
)
PROTOTYPE_MERGE_STRATEGY_CHOICES = (
    "dict_merge",
    "chrom_array_merge",
    "vectorized_merge",
    "array_sparse_merge_legacy",
    "array_sparse_merge",
    "array_dense_merge",
    "event_diff_merge"
)
PROTOTYPE_WRITER_STRATEGY_CHOICES = ("serial", "parallel_ordered")
PROTOTYPE_WRITE_MODE_CHOICES = ("full", "digest", "profile_only")
DIRECT_SPARSE_RESULT_FORMATS = (
    "direct_sparse_np",
    "direct_sparse_idx_np",
    "direct_sparse_local_np",
    "direct_sparse_touched_np",
    "direct_sparse_bincount_np",
    "direct_sparse_touched_bincount_np",
    "direct_sparse_local_bincount_np"
)


def get_alignment_chrom_sizes(
    path_bam: str,
    ref_fa: str | None = None
) -> dict[str, int]:
    """
    Read chromosome sizes from a BAM or CRAM header.
    """
    kwargs = {}
    if ref_fa is not None:
        kwargs["reference_filename"] = ref_fa

    with pysam.AlignmentFile(path_bam, "rb", **kwargs) as bam_file:
        return {
            chrom: size
            for chrom, size in zip(
                bam_file.references, bam_file.lengths, strict=True
            )
            if size is not None and size > 0
        }


def resolve_chrom_sizes(
    header_sizes: dict[str, int],
    path_chr_sizes: str | None = None
) -> dict[str, int]:
    """
    Resolve chromosome sizes from BAM/CRAM headers plus an optional TSV.
    """
    chrom_sizes = dict(header_sizes)

    if path_chr_sizes is not None:
        file_sizes = load_chr_sizes(path_chr_sizes)
        conflicts = [
            (chrom, chrom_sizes[chrom], size)
            for chrom, size in file_sizes.items()
            if chrom in chrom_sizes and chrom_sizes[chrom] != size
        ]
        if conflicts:
            chrom, header_size, file_size = conflicts[0]
            raise ValueError(
                "Conflicting chromosome sizes for "
                f"{chrom!r}: BAM/CRAM header has {header_size}, "
                f"but chr.sizes file has {file_size}."
            )
        chrom_sizes.update(file_sizes)

    if not chrom_sizes:
        raise ValueError(
            "Chromosome sizes are required to trim bedGraph bins. Provide a "
            "BAM/CRAM with sequence lengths in the header or pass "
            "'--chr_sizes' with a UCSC-style two-column TSV."
        )

    return chrom_sizes


def start_profile(
    args: argparse.Namespace, fil_out: str, out_fmt: str
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
        "executor_mode": args.executor_mode,
        "prototype_parse_strategy": args.prototype_parse_strategy,
        "prototype_bed_strategy": args.prototype_bed_strategy,
        "prototype_result_format": args.prototype_result_format,
        "prototype_merge_strategy": args.prototype_merge_strategy,
        "prototype_writer_strategy": args.prototype_writer_strategy,
        "prototype_writer_workers": args.prototype_writer_workers,
        "prototype_write_mode": args.prototype_write_mode,
        "threads": args.threads,
        "fil_in": args.fil_in,
        "ref_fa": args.ref_fa,
        "chr_sizes": args.chr_sizes,
        "fil_out": fil_out,
        "out_fmt": out_fmt,
        "method": args.method,
        "siz_bin": args.siz_bin,
        "chunk_size": args.chunk_size,
        "indexed_window_size": args.indexed_window_size,
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


def signal_dict_to_result(
    sig: dict[tuple[str, int], float],
    result_format: str,
    siz_bin: int
):
    """
    Convert a signal dict to one private benchmark result representation.
    """
    if result_format == "dict":
        return sig

    by_chrom: dict[str, list[tuple[int, float]]] = defaultdict(list)
    for (chrom, start), value in sig.items():
        by_chrom[chrom].append((start, value))

    if result_format in ("sparse_np", "direct_sparse_np"):
        return (
            result_format,
            [
                (
                    chrom,
                    np.asarray([start for start, _ in rows], dtype=np.int64),
                    np.asarray([value for _, value in rows], dtype=np.float64)
                )
                for chrom, rows in by_chrom.items()
            ]
        )

    if result_format == "dense_np":
        dense_parts = []
        for chrom, rows in by_chrom.items():
            starts = [start for start, _ in rows]
            first = min(starts)
            last = max(starts)
            values = np.zeros(
                ((last - first) // siz_bin) + 1, dtype=np.float64
            )
            for start, value in rows:
                values[(start - first) // siz_bin] += value
            dense_parts.append((chrom, first, values))
        return "dense_np", dense_parts

    raise ValueError(f"Unknown prototype result format: {result_format!r}.")


def count_result_bins(result) -> int:
    """
    Count signal bins for standard dict results and prototype tuple results.
    """
    if isinstance(result, tuple) and result:
        if len(result) == 2 and isinstance(result[1], int):
            return count_result_bins(result[0])
        tag = result[0]
        if tag in ("sparse_np", "direct_sparse_np", "direct_sparse_idx_np"):
            return sum(int(starts.size) for _, starts, _ in result[1])
        if tag == "array_sparse_merged":
            return sum(int(starts.size) for _, starts, _ in result[1])
        if tag == "direct_dense_np":
            return sum(
                int(np.count_nonzero(touched & (values != 0.0)))
                for _, values, touched in result[1]
            )
        if tag == "array_dense_merged":
            return sum(
                int(np.count_nonzero(touched & (values != 0.0)))
                for _, values, touched in result[2]
            )
        if tag == "event_np":
            return sum(
                int(edge_bins.size + diff_bins.size + touch_bins.size)
                for (
                    _,
                    _,
                    edge_bins,
                    _,
                    diff_bins,
                    _,
                    touch_bins,
                    _
                ) in result[1]
            )
        if tag == "dense_np":
            return sum(
                int(np.count_nonzero(values)) for _, _, values in result[1]
            )
        if isinstance(tag, dict):
            return len(tag)
    return len(result)


def estimate_result_payload_bytes(result) -> int:
    """
    Estimate NumPy-array payload bytes for private benchmark result formats.
    """
    if isinstance(result, tuple) and result:
        tag = result[0]
        if tag in ("sparse_np", "direct_sparse_np", "direct_sparse_idx_np"):
            return sum(
                int(starts.nbytes + values.nbytes)
                for _, starts, values in result[1]
            )
        if tag == "array_sparse_merged":
            return sum(
                int(starts.nbytes + values.nbytes)
                for _, starts, values in result[1]
            )
        if tag == "direct_dense_np":
            return sum(
                int(values.nbytes + touched.nbytes)
                for _, values, touched in result[1]
            )
        if tag == "array_dense_merged":
            return sum(
                int(values.nbytes + touched.nbytes)
                for _, values, touched in result[2]
            )
        if tag == "event_np":
            return sum(
                int(
                    edge_bins.nbytes
                    + edge_values.nbytes
                    + diff_bins.nbytes
                    + diff_values.nbytes
                    + touch_bins.nbytes
                    + touch_values.nbytes
                )
                for (
                    _,
                    _,
                    edge_bins,
                    edge_values,
                    diff_bins,
                    diff_values,
                    touch_bins,
                    touch_values
                ) in result[1]
            )
        if tag == "dense_np":
            return sum(int(values.nbytes) for _, _, values in result[1])
        if isinstance(tag, dict):
            return 0
    return 0


def check_indexed_alignment(
    path_bam: str,
    ref_fa: str | None,
    strategy: str
) -> None:
    """
    Validate prerequisites for indexed signal engines.
    """
    if strategy == "serial":
        return

    if path_bam == "-":
        raise ValueError(
            "Indexed signal engines require a named BAM or CRAM path."
        )

    if path_bam.lower().endswith(".cram") and ref_fa is None:
        raise ValueError(
            "Indexed CRAM signal engines require '--ref_fa' so every "
            "worker can open the CRAM deterministically."
        )

    kwargs = {}
    if ref_fa is not None:
        kwargs["reference_filename"] = ref_fa

    with pysam.AlignmentFile(path_bam, "rb", **kwargs) as bam_file:
        if not bam_file.has_index():
            raise ValueError(
                "Indexed signal engines require an alignment index "
                "(.bai for BAM or .crai for CRAM)."
            )


def has_indexed_alignment(path_bam: str, ref_fa: str | None = None) -> bool:
    """
    Return True when an alignment can be accessed through coordinate indexes.
    """
    kwargs = {}
    if ref_fa is not None:
        kwargs["reference_filename"] = ref_fa

    try:
        with pysam.AlignmentFile(path_bam, "rb", **kwargs) as bam_file:
            return bool(bam_file.has_index())
    except (OSError, ValueError):
        return False


def keep_read(
    read,
    allow_sec: bool = True,
    allow_sup: bool = False,
    allow_dup: bool = True
) -> bool:
    """
    Return whether a read passes the shared alignment-filter policy.
    """
    if read.is_unmapped or read.reference_id < 0:
        return False
    if not allow_sec and read.is_secondary:
        return False
    if not allow_sup and read.is_supplementary:
        return False
    return not (not allow_dup and read.is_duplicate)


def read_to_fragment(
    read,
    get_reference_name,
    chrom_sizes: dict[str, int],
    usr_frg: int | None = None
) -> tuple[str, int, int, int] | None:
    """
    Convert one accepted alignment record to a processed fragment interval.
    """
    chrom = get_reference_name(read.reference_id)
    chrom_len = chrom_sizes.get(chrom)
    if chrom_len is None:
        raise ValueError(
            f"Chromosome {chrom!r} is missing from chromosome sizes. "
            "Provide '--chr_sizes' with a UCSC-style two-column TSV."
        )

    #  Handle paired-end alignments: one fragment is emitted per leftmost
    #  anchor in a proper pair.
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
            usr_frg
            if usr_frg is not None
            else read.query_alignment_length
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


def iter_bam_fragments(
    path_bam: str,
    chrom_sizes: dict[str, int],
    usr_frg: int | None = None,
    ref_fa: str | None = None,
    allow_sec: bool = True,
    allow_sup: bool = False,
    allow_dup: bool = True
) -> Iterator[tuple[str, int, int, int]]:
    """
    Stream processed fragment intervals from a BAM or CRAM file.
    """
    kwargs = {}
    if ref_fa is not None:
        kwargs["reference_filename"] = ref_fa

    with pysam.AlignmentFile(path_bam, "rb", **kwargs) as bam_file:
        for read in bam_file.fetch(until_eof=True):
            if not keep_read(read, allow_sec, allow_sup, allow_dup):
                continue

            fragment = read_to_fragment(
                read,
                bam_file.get_reference_name,
                chrom_sizes,
                usr_frg
            )
            if fragment is not None:
                yield fragment


def iter_indexed_fragments(
    path_bam: str,
    chrom: str,
    chrom_sizes: dict[str, int],
    usr_frg: int | None = None,
    ref_fa: str | None = None,
    start: int | None = None,
    end: int | None = None,
    allow_sec: bool = True,
    allow_sup: bool = False,
    allow_dup: bool = True
) -> Iterator[tuple[str, int, int, int]]:
    """
    Stream processed fragments from one indexed chromosome or window fetch.
    """
    kwargs = {}
    if ref_fa is not None:
        kwargs["reference_filename"] = ref_fa

    with pysam.AlignmentFile(path_bam, "rb", **kwargs) as bam_file:
        if start is None and end is None:
            reads = bam_file.fetch(chrom)
        else:
            reads = bam_file.fetch(chrom, start, end)

        for read in reads:
            if not keep_read(read, allow_sec, allow_sup, allow_dup):
                continue

            if (
                start is not None
                and end is not None
                and (read.reference_start < start or read.reference_start >= end)
            ):
                #  Partition by read anchor, not computed fragment start.
                #  This prevents double-counting while preserving reverse-
                #  strand single-end fragments that extend left of a window.
                continue

            fragment = read_to_fragment(
                read,
                bam_file.get_reference_name,
                chrom_sizes,
                usr_frg
            )
            if fragment is not None:
                yield fragment


def collect_fragment_arrays_from_iter(
    fragments: Iterator[tuple[str, int, int, int]]
) -> tuple[dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]], int]:
    """
    Collect processed fragments into per-chromosome NumPy arrays.
    """
    starts: dict[str, list[int]] = defaultdict(list)
    ends: dict[str, list[int]] = defaultdict(list)
    lengths: dict[str, list[int]] = defaultdict(list)
    frg_tot = 0

    for chrom, frg_start, frg_end, frg_len in fragments:
        starts[chrom].append(frg_start)
        ends[chrom].append(frg_end)
        lengths[chrom].append(frg_len)
        frg_tot += 1

    arrays: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for chrom in starts:
        arrays[chrom] = (
            np.asarray(starts[chrom], dtype=np.int64),
            np.asarray(ends[chrom], dtype=np.int64),
            np.asarray(lengths[chrom], dtype=np.float64)
        )

    return arrays, frg_tot


def collect_bed_arrays_from_iter(
    fragments: Iterator[tuple[str, int, int, int]]
) -> tuple[str, list[tuple[str, np.ndarray, np.ndarray, np.ndarray]], int]:
    """
    Collect processed fragments into per-chromosome BED output arrays.
    """
    starts: dict[str, list[int]] = defaultdict(list)
    ends: dict[str, list[int]] = defaultdict(list)
    lengths: dict[str, list[int]] = defaultdict(list)
    n_fragments = 0

    for chrom, frg_start, frg_end, frg_len in fragments:
        starts[chrom].append(frg_start)
        ends[chrom].append(frg_end)
        lengths[chrom].append(frg_len)
        n_fragments += 1

    arrays = [
        (
            chrom,
            np.asarray(starts[chrom], dtype=np.int64),
            np.asarray(ends[chrom], dtype=np.int64),
            np.asarray(lengths[chrom], dtype=np.int64)
        )
        for chrom in starts
    ]
    return "bed_np", arrays, n_fragments


def iter_fragment_chunks(
    path_bam: str,
    chrom_sizes: dict[str, int],
    chunk_size: int,
    usr_frg: int | None = None,
    ref_fa: str | None = None
) -> Iterator[list[tuple[str, int, int, int]]]:
    """
    Stream processed fragments in fixed-size chunks.
    """
    check_cmp(chunk_size, "gt", 0, "chunk_size", allow_none=False)

    chunk: list[tuple[str, int, int, int]] = []
    for fragment in iter_bam_fragments(
        path_bam=path_bam,
        chrom_sizes=chrom_sizes,
        usr_frg=usr_frg,
        ref_fa=ref_fa
    ):
        chunk.append(fragment)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []

    if chunk:
        yield chunk


def count_fragments(
    path_bam: str,
    chrom_sizes: dict[str, int],
    usr_frg: int | None = None,
    ref_fa: str | None = None
) -> int:
    """
    Count processed fragments after filtering and coordinate clamping.
    """
    return sum(
        1
        for _ in iter_bam_fragments(
            path_bam=path_bam,
            chrom_sizes=chrom_sizes,
            usr_frg=usr_frg,
            ref_fa=ref_fa
        )
    )


def parse_bam(
    path_bam: str,
    usr_frg: int | None = None,
    ref_fa: str | None = None,
    chrom_sizes: dict[str, int] | None = None,
    allow_sec: bool = True,
    allow_sup: bool = False,
    allow_dup: bool = True
) -> dict[str, list[tuple[int, int, int]]]:
    """
    Parse a BAM or CRAM into chromosome-grouped fragment intervals.

    Parameters
    ----------
        path_bam : str
            Input BAM or CRAM path.
        usr_frg : int | None
            Optional fixed fragment length override. If provided:
                - paired-end read alignments: overrides TLEN for leftmost
                                              anchors.
                - single-end read alignments: used for both strands from the 5’
                                              end.
        ref_fa : str | None
            Reference FASTA file for CRAM decoding.
        allow_sec : bool
            Include secondary (multi-mapping) alignments. Default True.
        allow_sup : bool
            Include supplementary alignments. Default False.
        allow_dup : bool
            Include duplicate-marked alignments. Default True.

    Returns
    -------
        frg_tup : dict[str, list[tuple[int, int, int]]]
            A dictionary mapping chromosomes to lists of (start, end, frg_len),
            where [start, end) is a half-open interval in 0-based coordinates
            and frg_len is the fragment length used for that interval.

    Raises
    ------
        FileNotFoundError
            If 'path_bam' does not exist.
        ValueError
            If 'usr_frg' is provided but is <= 0 for paired-end or single-end
            extension.

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
        if chrom_sizes is None:
            header_sizes = get_alignment_chrom_sizes(path_bam, ref_fa)
            chrom_sizes = resolve_chrom_sizes(header_sizes)

        for chrom, frg_start, frg_end, frg_len in iter_bam_fragments(
            path_bam=path_bam,
            chrom_sizes=chrom_sizes,
            usr_frg=usr_frg,
            ref_fa=ref_fa,
            allow_sec=allow_sec,
            allow_sup=allow_sup,
            allow_dup=allow_dup
        ):
            frg_tup[chrom].append((frg_start, frg_end, frg_len))

    except FileNotFoundError:
        print(
            f"Error: BAM or CRAM file '{path_bam}' not found.",
            file=sys.stderr
        )
        raise
    except ValueError:
        #  Let caller or outer wrapper handle ValueError details
        raise
    except Exception as e:
        print(
            f"Unexpected error with BAM or CRAM file '{path_bam}': {e}",
            file=sys.stderr
        )
        raise

    return frg_tup


def collect_fragment_arrays(
    path_bam: str,
    chrom_sizes: dict[str, int],
    usr_frg: int | None = None,
    ref_fa: str | None = None,
    allow_sec: bool = True,
    allow_sup: bool = False,
    allow_dup: bool = True
) -> tuple[dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]], int]:
    """
    Parse a BAM or CRAM once into per-chromosome NumPy fragment arrays.
    """
    return collect_fragment_arrays_from_iter(
        iter_bam_fragments(
            path_bam=path_bam,
            chrom_sizes=chrom_sizes,
            usr_frg=usr_frg,
            ref_fa=ref_fa,
            allow_sec=allow_sec,
            allow_sup=allow_sup,
            allow_dup=allow_dup
        )
    )


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

    Parameters
    ----------
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
            Scaling factor applied to signal.

    Returns
    -------
        sig_bin : dict[tuple[str, int], float]
            A dictionary of binned signal data, where keys are
            '(chrom, bin_start)' and values are per-bin signal scores.

    Raises
    ------
        ValueError
            - If 'siz_bin' <= 0.
            - If 'is_len' or 'is_norm' is True and any 'frg_len' <= 0.
            - If 'is_norm' is True and 'frg_tot' <= 0.
            - If a provided 'scl_fct' <= 0.

    Notes
    -----
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


def calc_sig_chrom_array(
    chrom: str,
    starts: np.ndarray,
    ends: np.ndarray,
    lengths: np.ndarray,
    chrom_size: int,
    frg_tot: int,
    siz_bin: int,
    is_len: bool,
    is_norm: bool,
    scl_fct: float | None = None
) -> dict[tuple[str, int], float]:
    """
    Compute binned signal for one chromosome with vectorized range additions.
    """
    check_cmp(siz_bin, "gt", 0, "siz_bin", allow_none=False)

    if chrom_size <= 0:
        raise ValueError(f"Chromosome size must be > 0 for {chrom!r}.")

    n_bins = math.ceil(chrom_size / siz_bin)
    if starts.size == 0:
        return {}

    valid = ends > starts
    if not np.all(valid):
        starts = starts[valid]
        ends = ends[valid]
        lengths = lengths[valid]

    if starts.size == 0:
        return {}

    if is_len or is_norm:
        if np.any(lengths <= 0):
            raise ValueError(
                "'frg_len' must be > 0 when using normalization."
            )
        weights = 1.0 / lengths
    else:
        weights = np.ones(starts.shape, dtype=np.float64)

    sig = np.zeros(n_bins, dtype=np.float64)
    touched = np.zeros(n_bins, dtype=bool)
    start_bins = starts // siz_bin
    end_bins = (ends - 1) // siz_bin

    same = start_bins == end_bins
    if np.any(same):
        np.add.at(
            sig,
            start_bins[same],
            (ends[same] - starts[same]) * weights[same]
        )
        touched[start_bins[same]] = True

    multi = ~same
    if np.any(multi):
        s = starts[multi]
        e = ends[multi]
        sb = start_bins[multi]
        eb = end_bins[multi]
        w = weights[multi]

        left_overlap = ((sb + 1) * siz_bin - s) * w
        right_overlap = (e - eb * siz_bin) * w
        np.add.at(sig, sb, left_overlap)
        np.add.at(sig, eb, right_overlap)
        touched[sb] = True
        touched[eb] = True

        has_interior = eb > sb + 1
        if np.any(has_interior):
            diff = np.zeros(n_bins + 1, dtype=np.float64)
            touched_diff = np.zeros(n_bins + 1, dtype=np.int64)
            first = sb[has_interior] + 1
            last_exclusive = eb[has_interior]
            value = siz_bin * w[has_interior]
            np.add.at(diff, first, value)
            np.add.at(diff, last_exclusive, -value)
            np.add.at(touched_diff, first, 1)
            np.add.at(touched_diff, last_exclusive, -1)
            sig += np.cumsum(diff[:-1])
            touched |= np.cumsum(touched_diff[:-1]) > 0

    if is_norm:
        if frg_tot <= 0:
            raise ValueError(
                "Normalization requires non-zero total fragments."
            )
        sig /= frg_tot

    if scl_fct is not None:
        check_cmp(scl_fct, "gt", 0, "scl_fct")
        sig *= scl_fct

    idx = np.flatnonzero(touched & (sig != 0.0))
    return {
        (chrom, int(i) * siz_bin): float(sig[i])
        for i in idx
    }


def calc_sig_chrom_direct_sparse_np(
    chrom: str,
    starts: np.ndarray,
    ends: np.ndarray,
    lengths: np.ndarray,
    chrom_size: int,
    siz_bin: int,
    is_len: bool,
    return_bin_indices: bool = False,
    local_span: bool = False,
    use_bincount: bool = False,
    emit_touched_only: bool = False
):
    """
    Compute binned signal and return sparse NumPy arrays without a dict.
    """
    check_cmp(siz_bin, "gt", 0, "siz_bin", allow_none=False)

    if chrom_size <= 0:
        raise ValueError(f"Chromosome size must be > 0 for {chrom!r}.")

    genome_n_bins = math.ceil(chrom_size / siz_bin)
    if starts.size == 0:
        return (
            "direct_sparse_idx_np" if return_bin_indices else "direct_sparse_np",
            []
        )

    valid = ends > starts
    if not np.all(valid):
        starts = starts[valid]
        ends = ends[valid]
        lengths = lengths[valid]

    if starts.size == 0:
        return (
            "direct_sparse_idx_np" if return_bin_indices else "direct_sparse_np",
            []
        )

    if is_len:
        if np.any(lengths <= 0):
            raise ValueError(
                "'frg_len' must be > 0 when using normalization."
            )
        weights = 1.0 / lengths
    else:
        weights = np.ones(starts.shape, dtype=np.float64)

    start_bins_global = starts // siz_bin
    end_bins_global = (ends - 1) // siz_bin
    if local_span:
        first_bin = int(np.min(start_bins_global))
        last_bin = int(np.max(end_bins_global))
    else:
        first_bin = 0
        last_bin = genome_n_bins - 1
    first_bin = max(0, first_bin)
    last_bin = min(genome_n_bins - 1, last_bin)
    n_bins = last_bin - first_bin + 1
    if n_bins <= 0:
        return (
            "direct_sparse_idx_np" if return_bin_indices else "direct_sparse_np",
            []
        )

    sig = np.zeros(n_bins, dtype=np.float64)
    touched = np.zeros(n_bins, dtype=bool) if emit_touched_only else None
    start_bins = start_bins_global - first_bin
    end_bins = end_bins_global - first_bin

    same = start_bins == end_bins
    if np.any(same):
        same_values = (ends[same] - starts[same]) * weights[same]
        if use_bincount:
            sig += np.bincount(
                start_bins[same],
                weights=same_values,
                minlength=n_bins
            )[:n_bins]
        else:
            np.add.at(sig, start_bins[same], same_values)
        if touched is not None:
            touched[start_bins[same]] = True

    multi = ~same
    if np.any(multi):
        s = starts[multi]
        e = ends[multi]
        sb = start_bins[multi]
        eb = end_bins[multi]
        sb_global = start_bins_global[multi]
        eb_global = end_bins_global[multi]
        w = weights[multi]

        left_overlap = ((sb_global + 1) * siz_bin - s) * w
        right_overlap = (e - eb_global * siz_bin) * w
        if use_bincount:
            sig += np.bincount(sb, weights=left_overlap, minlength=n_bins)[:n_bins]
            sig += np.bincount(eb, weights=right_overlap, minlength=n_bins)[:n_bins]
        else:
            np.add.at(sig, sb, left_overlap)
            np.add.at(sig, eb, right_overlap)
        if touched is not None:
            touched[sb] = True
            touched[eb] = True

        has_interior = eb > sb + 1
        if np.any(has_interior):
            first = sb[has_interior] + 1
            last_exclusive = eb[has_interior]
            value = siz_bin * w[has_interior]
            touched_diff = None
            if use_bincount:
                diff = np.bincount(
                    first,
                    weights=value,
                    minlength=n_bins + 1
                ).astype(np.float64, copy=False)
                diff -= np.bincount(
                    last_exclusive,
                    weights=value,
                    minlength=n_bins + 1
                ).astype(np.float64, copy=False)
            else:
                diff = np.zeros(n_bins + 1, dtype=np.float64)
                np.add.at(diff, first, value)
                np.add.at(diff, last_exclusive, -value)
            if touched is not None:
                touched_diff = np.zeros(n_bins + 1, dtype=np.int64)
                np.add.at(touched_diff, first, 1)
                np.add.at(touched_diff, last_exclusive, -1)
            sig += np.cumsum(diff[:-1])
            if touched_diff is not None:
                touched |= np.cumsum(touched_diff[:-1]) > 0

    if touched is None:
        idx = np.flatnonzero(sig != 0.0)
    else:
        idx = np.flatnonzero(touched & (sig != 0.0))
    if idx.size == 0:
        return (
            "direct_sparse_idx_np" if return_bin_indices else "direct_sparse_np",
            []
        )

    if return_bin_indices:
        index_dtype = (
            np.int32
            if n_bins <= np.iinfo(np.int32).max
            else np.int64
        )
        return (
            "direct_sparse_idx_np",
            [
                (
                    chrom,
                    (idx + first_bin).astype(index_dtype, copy=False),
                    sig[idx].astype(np.float64, copy=False)
                )
            ]
        )

    return (
        "direct_sparse_np",
        [
            (
                chrom,
                ((idx + first_bin) * siz_bin).astype(np.int64, copy=False),
                sig[idx].astype(np.float64, copy=False)
            )
        ]
    )


def calc_sig_chrom_direct_dense_np(
    chrom: str,
    starts: np.ndarray,
    ends: np.ndarray,
    lengths: np.ndarray,
    chrom_size: int,
    siz_bin: int,
    is_len: bool
):
    """
    Compute binned signal and return one dense chromosome array.
    """
    check_cmp(siz_bin, "gt", 0, "siz_bin", allow_none=False)

    if chrom_size <= 0:
        raise ValueError(f"Chromosome size must be > 0 for {chrom!r}.")

    n_bins = math.ceil(chrom_size / siz_bin)
    if starts.size == 0:
        return "direct_dense_np", []

    valid = ends > starts
    if not np.all(valid):
        starts = starts[valid]
        ends = ends[valid]
        lengths = lengths[valid]

    if starts.size == 0:
        return "direct_dense_np", []

    if is_len:
        if np.any(lengths <= 0):
            raise ValueError(
                "'frg_len' must be > 0 when using normalization."
            )
        weights = 1.0 / lengths
    else:
        weights = np.ones(starts.shape, dtype=np.float64)

    values = np.zeros(n_bins, dtype=np.float64)
    touched = np.zeros(n_bins, dtype=bool)
    start_bins = starts // siz_bin
    end_bins = (ends - 1) // siz_bin

    same = start_bins == end_bins
    if np.any(same):
        np.add.at(
            values,
            start_bins[same],
            (ends[same] - starts[same]) * weights[same]
        )
        touched[start_bins[same]] = True

    multi = ~same
    if np.any(multi):
        s = starts[multi]
        e = ends[multi]
        sb = start_bins[multi]
        eb = end_bins[multi]
        w = weights[multi]

        np.add.at(values, sb, ((sb + 1) * siz_bin - s) * w)
        np.add.at(values, eb, (e - eb * siz_bin) * w)
        touched[sb] = True
        touched[eb] = True

        has_interior = eb > sb + 1
        if np.any(has_interior):
            diff = np.zeros(n_bins + 1, dtype=np.float64)
            touched_diff = np.zeros(n_bins + 1, dtype=np.int64)
            first = sb[has_interior] + 1
            last_exclusive = eb[has_interior]
            value = siz_bin * w[has_interior]
            np.add.at(diff, first, value)
            np.add.at(diff, last_exclusive, -value)
            np.add.at(touched_diff, first, 1)
            np.add.at(touched_diff, last_exclusive, -1)
            values += np.cumsum(diff[:-1])
            touched |= np.cumsum(touched_diff[:-1]) > 0

    if not np.any(touched & (values != 0.0)):
        return "direct_dense_np", []

    return "direct_dense_np", [(chrom, values, touched)]


def calc_sig_chrom_event_np(
    chrom: str,
    starts: np.ndarray,
    ends: np.ndarray,
    lengths: np.ndarray,
    chrom_size: int,
    siz_bin: int,
    is_len: bool
):
    """
    Return bin edge and range-add events for parent-side materialization.
    """
    check_cmp(siz_bin, "gt", 0, "siz_bin", allow_none=False)

    if chrom_size <= 0:
        raise ValueError(f"Chromosome size must be > 0 for {chrom!r}.")

    n_bins = math.ceil(chrom_size / siz_bin)
    if starts.size == 0:
        return "event_np", []

    valid = ends > starts
    if not np.all(valid):
        starts = starts[valid]
        ends = ends[valid]
        lengths = lengths[valid]

    if starts.size == 0:
        return "event_np", []

    if is_len:
        if np.any(lengths <= 0):
            raise ValueError(
                "'frg_len' must be > 0 when using normalization."
            )
        weights = 1.0 / lengths
    else:
        weights = np.ones(starts.shape, dtype=np.float64)

    start_bins = starts // siz_bin
    end_bins = (ends - 1) // siz_bin

    edge_bin_parts = []
    edge_value_parts = []
    diff_bin_parts = []
    diff_value_parts = []
    touch_bin_parts = []
    touch_value_parts = []

    same = start_bins == end_bins
    if np.any(same):
        same_bins = start_bins[same].astype(np.int64, copy=False)
        edge_bin_parts.append(same_bins)
        edge_value_parts.append(
            ((ends[same] - starts[same]) * weights[same]).astype(
                np.float64,
                copy=False
            )
        )

    multi = ~same
    if np.any(multi):
        s = starts[multi]
        e = ends[multi]
        sb = start_bins[multi]
        eb = end_bins[multi]
        w = weights[multi]

        edge_bin_parts.extend([
            sb.astype(np.int64, copy=False),
            eb.astype(np.int64, copy=False)
        ])
        edge_value_parts.extend([
            (((sb + 1) * siz_bin - s) * w).astype(np.float64, copy=False),
            ((e - eb * siz_bin) * w).astype(np.float64, copy=False)
        ])
        has_interior = eb > sb + 1
        if np.any(has_interior):
            first = (sb[has_interior] + 1).astype(np.int64, copy=False)
            last_exclusive = eb[has_interior].astype(np.int64, copy=False)
            value = (siz_bin * w[has_interior]).astype(
                np.float64,
                copy=False
            )
            diff_bin_parts.extend([first, last_exclusive])
            diff_value_parts.extend([value, -value])
            touch_bin_parts.extend([first, last_exclusive])
            touch_value_parts.extend([
                np.ones(first.shape, dtype=np.int64),
                -np.ones(last_exclusive.shape, dtype=np.int64)
            ])

    edge_bins = (
        np.concatenate(edge_bin_parts)
        if edge_bin_parts else np.empty(0, dtype=np.int64)
    )
    edge_values = (
        np.concatenate(edge_value_parts)
        if edge_value_parts else np.empty(0, dtype=np.float64)
    )
    diff_bins = (
        np.concatenate(diff_bin_parts)
        if diff_bin_parts else np.empty(0, dtype=np.int64)
    )
    diff_values = (
        np.concatenate(diff_value_parts)
        if diff_value_parts else np.empty(0, dtype=np.float64)
    )
    touch_bins = (
        np.concatenate(touch_bin_parts)
        if touch_bin_parts else np.empty(0, dtype=np.int64)
    )
    touch_values = (
        np.concatenate(touch_value_parts)
        if touch_value_parts else np.empty(0, dtype=np.int64)
    )

    if edge_bins.size == 0 and diff_bins.size == 0:
        return "event_np", []

    return (
        "event_np",
        [
            (
                chrom,
                n_bins,
                edge_bins,
                edge_values,
                diff_bins,
                diff_values,
                touch_bins,
                touch_values
            )
        ]
    )


def calc_sig_task(data):
    """
    Unpack one worker-task tuple and dispatch to 'calc_sig_chrom()'.
    """
    return calc_sig_chrom(*data)


def calc_sig_array_task(data):
    """
    Unpack one worker-task tuple and dispatch to 'calc_sig_chrom_array()'.
    """
    return calc_sig_chrom_array(*data)


def calc_sig_indexed_fetch_task(data):
    """
    Fetch one indexed region, parse fragments, and compute unscaled signal.
    """
    (
        path_bam,
        chrom,
        chrom_sizes,
        usr_frg,
        ref_fa,
        start,
        end,
        siz_bin,
        is_len,
        result_format
    ) = data

    frg_arr, frg_tot = collect_fragment_arrays_from_iter(
        iter_indexed_fragments(
            path_bam=path_bam,
            chrom=chrom,
            chrom_sizes=chrom_sizes,
            usr_frg=usr_frg,
            ref_fa=ref_fa,
            start=start,
            end=end
        )
    )

    arrays = frg_arr.get(chrom)
    if arrays is None:
        if result_format in DIRECT_SPARSE_RESULT_FORMATS:
            tag = (
                "direct_sparse_idx_np"
                if result_format == "direct_sparse_idx_np"
                else "direct_sparse_np"
            )
            sig = (tag, [])
        elif result_format in ("direct_dense_np", "event_np"):
            sig = (result_format, [])
        else:
            sig = {}
    elif result_format in DIRECT_SPARSE_RESULT_FORMATS:
        sig = calc_sig_chrom_direct_sparse_np(
            chrom=chrom,
            starts=arrays[0],
            ends=arrays[1],
            lengths=arrays[2],
            chrom_size=chrom_sizes[chrom],
            siz_bin=siz_bin,
            is_len=is_len,
            return_bin_indices=result_format == "direct_sparse_idx_np",
            local_span=result_format in (
                "direct_sparse_local_np",
                "direct_sparse_local_bincount_np"
            ),
            use_bincount=result_format in (
                "direct_sparse_bincount_np",
                "direct_sparse_touched_bincount_np",
                "direct_sparse_local_bincount_np"
            ),
            emit_touched_only=result_format in (
                "direct_sparse_touched_np",
                "direct_sparse_touched_bincount_np"
            )
        )
    elif result_format == "direct_dense_np":
        sig = calc_sig_chrom_direct_dense_np(
            chrom=chrom,
            starts=arrays[0],
            ends=arrays[1],
            lengths=arrays[2],
            chrom_size=chrom_sizes[chrom],
            siz_bin=siz_bin,
            is_len=is_len
        )
    elif result_format == "event_np":
        sig = calc_sig_chrom_event_np(
            chrom=chrom,
            starts=arrays[0],
            ends=arrays[1],
            lengths=arrays[2],
            chrom_size=chrom_sizes[chrom],
            siz_bin=siz_bin,
            is_len=is_len
        )
    else:
        sig = calc_sig_chrom_array(
            chrom=chrom,
            starts=arrays[0],
            ends=arrays[1],
            lengths=arrays[2],
            chrom_size=chrom_sizes[chrom],
            frg_tot=1,
            siz_bin=siz_bin,
            is_len=is_len,
            is_norm=False,
            scl_fct=None
        )

    if result_format in (*DIRECT_SPARSE_RESULT_FORMATS,
        "direct_dense_np",
        "event_np"
    ):
        return sig, frg_tot

    return signal_dict_to_result(sig, result_format, siz_bin), frg_tot


def collect_bed_indexed_task(data):
    """
    Fetch one indexed region and return compact BED row arrays.
    """
    (
        path_bam,
        chrom,
        chrom_sizes,
        usr_frg,
        ref_fa,
        start,
        end
    ) = data

    return collect_bed_arrays_from_iter(
        iter_indexed_fragments(
            path_bam=path_bam,
            chrom=chrom,
            chrom_sizes=chrom_sizes,
            usr_frg=usr_frg,
            ref_fa=ref_fa,
            start=start,
            end=end
        )
    )


def count_bed_array_rows(result) -> int:
    """
    Count rows in one compact BED array result.
    """
    if not isinstance(result, tuple) or len(result) < 2:
        return 0
    return sum(int(starts.size) for _chrom, starts, _ends, _lengths in result[1])


def estimate_bed_payload_bytes(result) -> int:
    """
    Estimate NumPy payload bytes in one compact BED array result.
    """
    if not isinstance(result, tuple) or len(result) < 2:
        return 0
    return sum(
        int(starts.nbytes + ends.nbytes + lengths.nbytes)
        for _chrom, starts, ends, lengths in result[1]
    )


def calc_sig_chunk_task(data):
    """
    Unpack one worker-task tuple and compute signal for a fragment chunk.
    """
    fragments, frg_tot, siz_bin, is_len, is_norm, scl_fct = data

    by_chrom = defaultdict(list)
    for chrom, start, end, frg_len in fragments:
        by_chrom[chrom].append((start, end, frg_len))

    sig_chunk = defaultdict(float)
    for chrom, entries in by_chrom.items():
        sig_chrom = calc_sig_chrom(
            chrom=chrom,
            frg_tup=entries,
            frg_tot=frg_tot,
            siz_bin=siz_bin,
            is_len=is_len,
            is_norm=is_norm,
            scl_fct=scl_fct
        )
        for key, value in sig_chrom.items():
            sig_chunk[key] += value

    return sig_chunk


def calc_sig_profile_task(data):
    """
    Compute one profiled worker task and return timing metadata.
    """
    kind, task_id, payload = data
    start = time.perf_counter()
    if kind == "chrom":
        result = calc_sig_array_task(payload)
    elif kind in ("indexed_chrom", "indexed_window"):
        result = calc_sig_indexed_fetch_task(payload)
    elif kind == "indexed_bed":
        result = collect_bed_indexed_task(payload)
    elif kind == "chunk":
        result = calc_sig_chunk_task(payload)
    else:
        raise ValueError(f"Unknown profiled task kind: {kind!r}.")
    end = time.perf_counter()
    result_count = (
        count_bed_array_rows(result)
        if kind == "indexed_bed"
        else count_result_bins(result)
    )

    return task_id, {
        "worker_start_s": start,
        "worker_end_s": end,
        "worker_elapsed_s": end - start,
        "result_bins": result_count,
    }, result


def iter_task_results(
    kind: str,
    task_func,
    task_data,
    threads: int,
    executor_mode: str,
    task_profiles: list[dict] | None = None
):
    """
    Yield task results, optionally recording per-task worker timings.
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
                task_profiles[task_id].update({
                    "worker_start_s": start,
                    "worker_end_s": end,
                    "worker_elapsed_s": end - start,
                    "parent_received_s": received,
                    "parent_receive_lag_s": received - end,
                    "result_bins": (
                        count_bed_array_rows(result)
                        if kind == "indexed_bed"
                        else count_result_bins(result)
                    ),
                })
                yield result
        return

    max_workers = min(threads, os.cpu_count() or 1)
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        if task_profiles is None:
            if executor_mode == "map":
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
        if executor_mode == "map":
            for task_id, timing, result in executor.map(
                calc_sig_profile_task,
                profiled_data
            ):
                received = time.perf_counter()
                timing["parent_received_s"] = received
                timing["parent_receive_lag_s"] = (
                    received - timing["worker_end_s"]
                )
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
                timing["parent_received_s"] = received
                timing["parent_receive_lag_s"] = (
                    received - timing["worker_end_s"]
                )
                task_profiles[task_id].update(timing)
                yield result


def summarize_task_profiles(profile: dict | None) -> None:
    """
    Add derived worker timing summaries to an active profile.
    """
    if profile is None:
        return

    worker_elapsed = [
        task["worker_elapsed_s"]
        for task in profile.get("tasks", [])
        if "worker_elapsed_s" in task
    ]
    receive_lag = [
        task["parent_receive_lag_s"]
        for task in profile.get("tasks", [])
        if "parent_receive_lag_s" in task
    ]
    if worker_elapsed:
        profile["worker_elapsed_sum_s"] = sum(worker_elapsed)
        profile["worker_elapsed_max_s"] = max(worker_elapsed)
    if receive_lag:
        profile["parent_receive_lag_sum_s"] = sum(receive_lag)
        profile["parent_receive_lag_max_s"] = max(receive_lag)


def apply_signal_adjustments(
    signal_bins: dict[tuple[str, int], float],
    frg_tot: int,
    is_norm: bool,
    scl_fct: float | None = None
) -> None:
    """
    Apply global depth normalization and scaling to already merged signal.
    """
    if is_array_signal(signal_bins):
        if is_norm:
            if frg_tot <= 0:
                raise ValueError(
                    "Normalization requires non-zero total fragments."
                )
            for values in iter_signal_value_arrays(signal_bins):
                values /= frg_tot

        if scl_fct is not None:
            check_cmp(scl_fct, "gt", 0, "scl_fct")
            for values in iter_signal_value_arrays(signal_bins):
                values *= scl_fct
        return

    if is_norm:
        if frg_tot <= 0:
            raise ValueError(
                "Normalization requires non-zero total fragments."
            )
        for key in signal_bins:
            signal_bins[key] /= frg_tot

    if scl_fct is not None:
        check_cmp(scl_fct, "gt", 0, "scl_fct")
        for key in signal_bins:
            signal_bins[key] *= scl_fct


def is_array_sparse_signal(signal_bins) -> bool:
    """
    Return True for the merged sparse-array signal container.
    """
    return (
        isinstance(signal_bins, tuple)
        and len(signal_bins) == 2
        and signal_bins[0] == "array_sparse_merged"
    )


def is_array_dense_signal(signal_bins) -> bool:
    """
    Return True for the merged dense-array signal container.
    """
    return (
        isinstance(signal_bins, tuple)
        and len(signal_bins) == 3
        and signal_bins[0] == "array_dense_merged"
    )


def is_array_signal(signal_bins) -> bool:
    """
    Return True for merged signal containers that avoid the final dict.
    """
    return is_array_sparse_signal(signal_bins) or is_array_dense_signal(signal_bins)


def iter_signal_value_arrays(signal_bins):
    """
    Yield mutable value arrays from merged array-backed signal containers.
    """
    if is_array_sparse_signal(signal_bins):
        for _, _, values in signal_bins[1]:
            yield values
    elif is_array_dense_signal(signal_bins):
        for _, values, _ in signal_bins[2]:
            yield values


def signal_row_count(signal_bins) -> int:
    """
    Count output bedGraph rows without assuming a dict-backed signal.
    """
    if is_array_sparse_signal(signal_bins):
        return sum(int(starts.size) for _, starts, _ in signal_bins[1])
    if is_array_dense_signal(signal_bins):
        return sum(
            int(np.count_nonzero(touched & (values != 0.0)))
            for _, values, touched in signal_bins[2]
        )
    return len(signal_bins)


def sparse_parts_are_ordered(
    starts_parts: list[np.ndarray]
) -> bool:
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


def coalesce_ordered_sparse_arrays(
    starts: np.ndarray,
    values: np.ndarray
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
        summed[keep].astype(np.float64, copy=False)
    )


def coalesce_sparse_signal_parts(
    parts_by_chrom: dict[str, list[tuple[np.ndarray, np.ndarray]]],
    optimize: bool = True,
    stats: dict[str, int] | None = None
):
    """
    Sort and sum sparse signal arrays by chromosome without building a dict.
    """
    if stats is not None:
        stats.setdefault("single_chroms", 0)
        stats.setdefault("ordered_chroms", 0)
        stats.setdefault("sorted_chroms", 0)
        stats.setdefault("parts", 0)
        stats.setdefault("rows_before", 0)
        stats.setdefault("rows_after", 0)

    merged = []
    for chrom in sorted(parts_by_chrom, key=sort_chrom):
        starts_parts = []
        values_parts = []
        for starts, values in parts_by_chrom[chrom]:
            if starts.size == 0:
                continue
            starts_parts.append(starts.astype(np.int64, copy=False))
            values_parts.append(values.astype(np.float64, copy=False))

        if not starts_parts:
            continue

        if stats is not None:
            stats["parts"] += len(starts_parts)
            stats["rows_before"] += sum(int(starts.size) for starts in starts_parts)

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

        is_ordered = optimize and sparse_parts_are_ordered(starts_parts)
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

        out_starts, out_values = coalesce_ordered_sparse_arrays(starts, values)
        if out_starts.size > 0:
            merged.append((
                chrom,
                out_starts,
                out_values
            ))
            if stats is not None:
                stats["rows_after"] += int(out_starts.size)

    return "array_sparse_merged", merged


def sparse_payload_to_starts(
    tag: str,
    coords: np.ndarray,
    siz_bin: int
) -> np.ndarray:
    """
    Convert sparse payload coordinates to bedGraph bin starts.
    """
    if tag == "direct_sparse_idx_np":
        return (coords.astype(np.int64, copy=False) * siz_bin)
    return coords.astype(np.int64, copy=False)


def materialize_event_signal_parts(
    parts_by_chrom: dict[str, list[tuple]],
    siz_bin: int
):
    """
    Combine event arrays and materialize final sparse signal arrays.
    """
    merged = []
    for chrom in sorted(parts_by_chrom, key=sort_chrom):
        n_bins = max(part[0] for part in parts_by_chrom[chrom])
        values = np.zeros(n_bins, dtype=np.float64)
        touched = np.zeros(n_bins, dtype=bool)

        edge_bins_parts = [part[1] for part in parts_by_chrom[chrom]]
        edge_values_parts = [part[2] for part in parts_by_chrom[chrom]]
        edge_bins = np.concatenate(edge_bins_parts) if edge_bins_parts else (
            np.empty(0, dtype=np.int64)
        )
        edge_values = (
            np.concatenate(edge_values_parts)
            if edge_values_parts else np.empty(0, dtype=np.float64)
        )
        if edge_bins.size > 0:
            np.add.at(values, edge_bins, edge_values)
            touched[edge_bins] = True

        diff_bins_parts = [part[3] for part in parts_by_chrom[chrom]]
        diff_values_parts = [part[4] for part in parts_by_chrom[chrom]]
        diff_bins = np.concatenate(diff_bins_parts) if diff_bins_parts else (
            np.empty(0, dtype=np.int64)
        )
        diff_values = (
            np.concatenate(diff_values_parts)
            if diff_values_parts else np.empty(0, dtype=np.float64)
        )
        if diff_bins.size > 0:
            diff = np.zeros(n_bins + 1, dtype=np.float64)
            np.add.at(diff, diff_bins, diff_values)
            values += np.cumsum(diff[:-1])

        touch_bins_parts = [part[5] for part in parts_by_chrom[chrom]]
        touch_values_parts = [part[6] for part in parts_by_chrom[chrom]]
        touch_bins = np.concatenate(touch_bins_parts) if touch_bins_parts else (
            np.empty(0, dtype=np.int64)
        )
        touch_values = (
            np.concatenate(touch_values_parts)
            if touch_values_parts else np.empty(0, dtype=np.int64)
        )
        if touch_bins.size > 0:
            touch_diff = np.zeros(n_bins + 1, dtype=np.int64)
            np.add.at(touch_diff, touch_bins, touch_values)
            touched |= np.cumsum(touch_diff[:-1]) > 0

        idx = np.flatnonzero(touched & (values != 0.0))
        if idx.size > 0:
            merged.append((
                chrom,
                (idx * siz_bin).astype(np.int64, copy=False),
                values[idx].astype(np.float64, copy=False)
            ))

    return "array_sparse_merged", merged


def merge_signal_results(
    results,
    merge_strategy: str,
    chrom_sizes: dict[str, int],
    siz_bin: int,
    is_prototype_strategy: bool,
    profile: dict | None = None
) -> tuple[defaultdict[tuple[str, int], float], int, int, int]:
    """
    Merge worker signal results and return signal, fragments, bins, and bytes.
    """
    sig_cmb: defaultdict[tuple[str, int], float] = defaultdict(float)
    array_cmb: dict[str, np.ndarray] = {}
    array_touched: dict[str, np.ndarray] = {}
    sparse_parts: dict[str, list[tuple[np.ndarray, np.ndarray]]] = defaultdict(list)
    event_parts: dict[str, list[tuple]] = defaultdict(list)
    frg_tot = 0
    n_result_bins = 0
    result_payload_bytes = 0
    metadata_s = 0.0
    accumulate_s = 0.0
    array_to_dict_s = 0.0
    array_coalesce_s = 0.0
    array_dense_finalize_s = 0.0
    event_materialize_s = 0.0
    sparse_coalesce_stats: dict[str, int] = {}

    def ensure_array(chrom: str) -> np.ndarray:
        """
        Ensure a dense chromosome array exists for one chromosomal key.

        Parameters
        ----------
        chrom
            Chromosome name used as the lookup key.

        Returns
        -------
        numpy.ndarray
            Mutable float64 array sized to the chromosome span under the
            current bin width.
        """
        if chrom not in array_cmb:
            array_cmb[chrom] = np.zeros(
                math.ceil(chrom_sizes[chrom] / siz_bin),
                dtype=np.float64
            )
        return array_cmb[chrom]

    def ensure_touched(chrom: str) -> np.ndarray:
        """
        Ensure a boolean touch mask exists for one chromosomal key.

        Parameters
        ----------
        chrom
            Chromosome name used as the lookup key.

        Returns
        -------
        numpy.ndarray
            Mutable boolean mask sized to the chromosome span under the
            current bin width.
        """
        if chrom not in array_touched:
            array_touched[chrom] = np.zeros(
                math.ceil(chrom_sizes[chrom] / siz_bin),
                dtype=bool
            )
        return array_touched[chrom]

    for result in results:
        time_meta = time.perf_counter()
        if is_prototype_strategy:
            sig_result, n_fragments = result
            frg_tot += n_fragments
        else:
            sig_result = result

        n_result_bins += count_result_bins(sig_result)
        result_payload_bytes += estimate_result_payload_bytes(sig_result)
        if profile is not None:
            metadata_s += time.perf_counter() - time_meta

        time_accumulate = time.perf_counter()
        if merge_strategy in ("array_sparse_merge", "array_sparse_merge_legacy"):
            if isinstance(sig_result, tuple) and sig_result:
                tag = sig_result[0]
                if tag in (
                    "sparse_np",
                    "direct_sparse_np",
                    "direct_sparse_idx_np"
                ):
                    for chrom, coords, values in sig_result[1]:
                        sparse_parts[chrom].append((
                            sparse_payload_to_starts(tag, coords, siz_bin),
                            values
                        ))
                elif tag == "dense_np":
                    for chrom, first, values in sig_result[1]:
                        idx = np.flatnonzero(values != 0.0)
                        sparse_parts[chrom].append((
                            first + idx.astype(np.int64) * siz_bin,
                            values[idx]
                        ))
                else:
                    by_chrom: dict[str, list[tuple[int, float]]] = defaultdict(list)
                    for (chrom, start), value in sig_result.items():
                        by_chrom[chrom].append((start, value))
                    for chrom, rows in by_chrom.items():
                        sparse_parts[chrom].append((
                            np.asarray([start for start, _ in rows], dtype=np.int64),
                            np.asarray([value for _, value in rows], dtype=np.float64)
                        ))
            else:
                by_chrom: dict[str, list[tuple[int, float]]] = defaultdict(list)
                for (chrom, start), value in sig_result.items():
                    by_chrom[chrom].append((start, value))
                for chrom, rows in by_chrom.items():
                    sparse_parts[chrom].append((
                        np.asarray([start for start, _ in rows], dtype=np.int64),
                        np.asarray([value for _, value in rows], dtype=np.float64)
                    ))
            if profile is not None:
                accumulate_s += time.perf_counter() - time_accumulate
            continue

        if merge_strategy == "array_dense_merge":
            if isinstance(sig_result, tuple) and sig_result:
                tag = sig_result[0]
                if tag == "direct_dense_np":
                    for chrom, values, touched in sig_result[1]:
                        arr = ensure_array(chrom)
                        arr[:values.size] += values
                        ensure_touched(chrom)[:touched.size] |= touched
                elif tag == "dense_np":
                    for chrom, first, values in sig_result[1]:
                        arr = ensure_array(chrom)
                        first_idx = first // siz_bin
                        arr[first_idx:first_idx + values.size] += values
                        ensure_touched(chrom)[
                            first_idx:first_idx + values.size
                        ] |= values != 0.0
                elif tag in (
                    "sparse_np",
                    "direct_sparse_np",
                    "direct_sparse_idx_np"
                ):
                    for chrom, coords, values in sig_result[1]:
                        if tag == "direct_sparse_idx_np":
                            idx = coords.astype(np.int64, copy=False)
                        else:
                            idx = coords // siz_bin
                        np.add.at(ensure_array(chrom), idx, values)
                        ensure_touched(chrom)[idx] = True
                else:
                    for (chrom, start), value in sig_result.items():
                        ensure_array(chrom)[start // siz_bin] += value
                        ensure_touched(chrom)[start // siz_bin] = True
            else:
                for (chrom, start), value in sig_result.items():
                    ensure_array(chrom)[start // siz_bin] += value
                    ensure_touched(chrom)[start // siz_bin] = True
            if profile is not None:
                accumulate_s += time.perf_counter() - time_accumulate
            continue

        if merge_strategy == "event_diff_merge":
            if isinstance(sig_result, tuple) and sig_result:
                tag = sig_result[0]
                if tag == "event_np":
                    for (
                        chrom,
                        n_bins,
                        edge_bins,
                        edge_values,
                        diff_bins,
                        diff_values,
                        touch_bins,
                        touch_values
                    ) in sig_result[1]:
                        event_parts[chrom].append((
                            n_bins,
                            edge_bins,
                            edge_values,
                            diff_bins,
                            diff_values,
                            touch_bins,
                            touch_values
                        ))
                elif tag in (
                    "sparse_np",
                    "direct_sparse_np",
                    "direct_sparse_idx_np"
                ):
                    for chrom, coords, values in sig_result[1]:
                        sparse_parts[chrom].append((
                            sparse_payload_to_starts(tag, coords, siz_bin),
                            values
                        ))
                else:
                    for (chrom, start), value in sig_result.items():
                        sparse_parts[chrom].append((
                            np.asarray([start], dtype=np.int64),
                            np.asarray([value], dtype=np.float64)
                        ))
            else:
                for (chrom, start), value in sig_result.items():
                    sparse_parts[chrom].append((
                        np.asarray([start], dtype=np.int64),
                        np.asarray([value], dtype=np.float64)
                    ))
            if profile is not None:
                accumulate_s += time.perf_counter() - time_accumulate
            continue

        if merge_strategy == "dict_merge":
            if isinstance(sig_result, tuple) and sig_result:
                tag = sig_result[0]
                if tag in (
                    "sparse_np",
                    "direct_sparse_np",
                    "direct_sparse_idx_np"
                ):
                    for chrom, coords, values in sig_result[1]:
                        starts = sparse_payload_to_starts(tag, coords, siz_bin)
                        for start, value in zip(starts, values, strict=True):
                            sig_cmb[(chrom, int(start))] += float(value)
                elif tag == "direct_dense_np":
                    for chrom, values, touched in sig_result[1]:
                        idx = np.flatnonzero(touched & (values != 0.0))
                        for i in idx:
                            sig_cmb[(chrom, int(i) * siz_bin)] += float(
                                values[i]
                            )
                elif tag == "dense_np":
                    for chrom, first, values in sig_result[1]:
                        for offset, value in enumerate(values):
                            if value != 0.0:
                                sig_cmb[(chrom, first + offset * siz_bin)] += (
                                    float(value)
                                )
                else:
                    for key, value in sig_result.items():
                        sig_cmb[key] += value
            else:
                for key, value in sig_result.items():
                    sig_cmb[key] += value
            if profile is not None:
                accumulate_s += time.perf_counter() - time_accumulate
            continue

        if isinstance(sig_result, tuple) and sig_result:
            tag = sig_result[0]
            if tag in (
                "sparse_np",
                "direct_sparse_np",
                "direct_sparse_idx_np"
            ):
                for chrom, coords, values in sig_result[1]:
                    arr = ensure_array(chrom)
                    if tag == "direct_sparse_idx_np":
                        idx = coords.astype(np.int64, copy=False)
                    else:
                        idx = coords // siz_bin
                    np.add.at(arr, idx, values)
            elif tag == "dense_np":
                for chrom, first, values in sig_result[1]:
                    arr = ensure_array(chrom)
                    first_idx = first // siz_bin
                    arr[first_idx:first_idx + values.size] += values
            else:
                for (chrom, start), value in sig_result.items():
                    ensure_array(chrom)[start // siz_bin] += value
        else:
            for (chrom, start), value in sig_result.items():
                ensure_array(chrom)[start // siz_bin] += value
        if profile is not None:
            accumulate_s += time.perf_counter() - time_accumulate

    dict_converting_merges = (
        "chrom_array_merge",
        "vectorized_merge"
    )
    if merge_strategy in dict_converting_merges:
        time_convert = time.perf_counter()
        for chrom, arr in array_cmb.items():
            idx = np.flatnonzero(arr != 0.0)
            for i in idx:
                sig_cmb[(chrom, int(i) * siz_bin)] = float(arr[i])
        if profile is not None:
            array_to_dict_s = time.perf_counter() - time_convert

    if merge_strategy in ("array_sparse_merge", "array_sparse_merge_legacy"):
        time_coalesce = time.perf_counter()
        sig_cmb = coalesce_sparse_signal_parts(
            sparse_parts,
            optimize=merge_strategy == "array_sparse_merge",
            stats=sparse_coalesce_stats
        )
        if profile is not None:
            array_coalesce_s = time.perf_counter() - time_coalesce

    if merge_strategy == "array_dense_merge":
        time_dense = time.perf_counter()
        sig_cmb = (
            "array_dense_merged",
            siz_bin,
            [
                (chrom, array_cmb[chrom], array_touched[chrom])
                for chrom in sorted(array_cmb, key=sort_chrom)
                if np.any(array_touched[chrom] & (array_cmb[chrom] != 0.0))
            ]
        )
        if profile is not None:
            array_dense_finalize_s = time.perf_counter() - time_dense

    if merge_strategy == "event_diff_merge":
        time_event = time.perf_counter()
        sig_cmb = materialize_event_signal_parts(event_parts, siz_bin)
        if sparse_parts:
            _, sparse_merged = coalesce_sparse_signal_parts(sparse_parts)
            combined_parts: dict[
                str,
                list[tuple[np.ndarray, np.ndarray]]
            ] = defaultdict(list)
            for chrom, starts, values in sig_cmb[1] + sparse_merged:
                combined_parts[chrom].append((starts, values))
            sig_cmb = coalesce_sparse_signal_parts(combined_parts)
        if profile is not None:
            event_materialize_s = time.perf_counter() - time_event

    if profile is not None:
        phases = profile["phases_s"]
        phases["merge_result_metadata"] = metadata_s
        phases["merge_accumulate_results"] = accumulate_s
        phases["merge_array_to_dict"] = array_to_dict_s
        phases["merge_array_coalesce"] = array_coalesce_s
        phases["merge_array_dense_finalize"] = array_dense_finalize_s
        phases["merge_event_materialize"] = event_materialize_s
        profile["sparse_coalesce_stats"] = sparse_coalesce_stats

    return sig_cmb, frg_tot, n_result_bins, result_payload_bytes


def format_bdg_rows(
    rows: list[tuple[tuple[str, int], float]],
    siz_bin: int,
    dp: int,
    chrom_sizes: dict[str, int] | None
) -> str:
    """
    Format a bedGraph row shard with the same row semantics as write_bdg().
    """
    lines = []
    for (chrom, bin_start), value in rows:
        bin_end = bin_start + siz_bin
        if chrom_sizes is not None:
            chrom_size = chrom_sizes.get(chrom)
            if chrom_size is None:
                raise ValueError(
                    f"Missing chromosome size for bedGraph row: {chrom!r}."
                )
            if bin_start < 0 or bin_start >= chrom_size:
                raise ValueError(
                    "bedGraph bin start is outside chromosome bounds: "
                    f"{chrom}:{bin_start} (chromosome size {chrom_size})."
                )
            bin_end = min(bin_end, chrom_size)

        lines.append(
            f"{chrom}\t{bin_start}\t{bin_end}\t"
            f"{format_bdg_value(value, dp)}\n"
        )

    return "".join(lines)


def format_bdg_value(value: float, dp: int) -> str:
    """
    Format one bedGraph value with the same rounding used by write_bdg().
    """
    v = round(float(value), dp)
    if v == 0.0:
        v = 0.0

    s_val = f"{v:.{dp}f}"
    if "." in s_val:
        s_val = s_val.rstrip("0").rstrip(".")
    if s_val == "-0":
        s_val = "0"
    return s_val


def iter_sorted_bed_array_rows(results):
    """
    Yield deterministic BED rows from compact worker array results.
    """
    by_chrom: dict[str, list[tuple[np.ndarray, np.ndarray, np.ndarray]]] = (
        defaultdict(list)
    )
    for result in results:
        if not isinstance(result, tuple) or len(result) < 2:
            continue
        for chrom, starts, ends, lengths in result[1]:
            if starts.size:
                by_chrom[chrom].append((starts, ends, lengths))

    for chrom in sorted(by_chrom, key=sort_chrom):
        starts = np.concatenate([part[0] for part in by_chrom[chrom]])
        ends = np.concatenate([part[1] for part in by_chrom[chrom]])
        lengths = np.concatenate([part[2] for part in by_chrom[chrom]])
        order = np.lexsort((lengths, ends, starts))
        for idx in order:
            yield (
                chrom,
                int(starts[idx]),
                int(ends[idx]),
                int(lengths[idx])
            )


def write_bed_array_results(results, fil_out: str) -> int:
    """
    Write compact BED worker results and return the emitted row count.
    """
    row_count = 0
    with open_out(fil_out) as bed_file:
        for chrom, start, end, length in iter_sorted_bed_array_rows(results):
            bed_file.write(f"{chrom}\t{start}\t{end}\t{length}\n")
            row_count += 1
    return row_count


def write_bed_fragment_dict(
    frg_tup: dict[str, list[tuple[int, int, int]]],
    fil_out: str
) -> int:
    """
    Write chromosome-grouped fragment tuples in deterministic BED order.
    """
    row_count = 0
    with open_out(fil_out) as bed_file:
        for chrom in sorted(frg_tup.keys(), key=sort_chrom):
            for start, end, length in sorted(frg_tup[chrom], key=lambda t: t[0]):
                bed_file.write(f"{chrom}\t{start}\t{end}\t{length}\n")
                row_count += 1
    return row_count


def format_bdg_array_rows(
    chrom: str,
    starts: np.ndarray,
    values: np.ndarray,
    siz_bin: int,
    dp: int,
    chrom_size: int | None
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
                    f"{chrom}:{bin_start} (chromosome size {chrom_size})."
                )
            bin_end = min(bin_end, chrom_size)

        lines.append(
            f"{chrom}\t{bin_start}\t{bin_end}\t"
            f"{format_bdg_value(float(value), dp)}\n"
        )
    return "".join(lines)


def format_bdg_rows_task(data) -> str:
    """
    ProcessPool-friendly wrapper for formatting one bedGraph shard.
    """
    return format_bdg_rows(*data)


def format_bdg_array_rows_task(data) -> str:
    """
    ProcessPool-friendly wrapper for sparse-array bedGraph formatting.
    """
    return format_bdg_array_rows(*data)


def iter_array_signal_chunks(signal_bins, target_chunks: int):
    """
    Yield ordered array-backed signal chunks for parallel formatting.
    """
    total_rows = max(signal_row_count(signal_bins), 1)
    chunk_size = max(1, math.ceil(total_rows / max(target_chunks, 1)))
    if is_array_sparse_signal(signal_bins):
        for chrom, starts, values in signal_bins[1]:
            for i in range(0, starts.size, chunk_size):
                yield chrom, starts[i:i + chunk_size], values[i:i + chunk_size]
    elif is_array_dense_signal(signal_bins):
        dense_siz_bin = signal_bins[1]
        for chrom, values, touched in signal_bins[2]:
            idx = np.flatnonzero(touched & (values != 0.0))
            starts = (idx * dense_siz_bin).astype(np.int64, copy=False)
            for i in range(0, idx.size, chunk_size):
                chunk_starts = starts[i:i + chunk_size]
                chunk_values = values[idx[i:i + chunk_size]]
                yield (
                    chrom,
                    chunk_starts,
                    chunk_values
                )


def digest_bdg_rows(
    cvg,
    siz_bin: int,
    dp: int,
    chrom_sizes: dict[str, int] | None
) -> tuple[int, str]:
    """
    Return row count and SHA-256 for deterministic bedGraph-formatted rows.
    """
    if is_array_signal(cvg):
        digest = hashlib.sha256()
        row_count = 0
        for chrom, starts, values in iter_array_signal_chunks(cvg, 1):
            text = format_bdg_array_rows(
                chrom,
                starts,
                values,
                siz_bin,
                dp,
                chrom_sizes.get(chrom) if chrom_sizes is not None else None
            )
            row_count += text.count("\n")
            digest.update(text.encode("utf-8"))
        return row_count, digest.hexdigest()

    digest = hashlib.sha256()
    row_count = 0
    items = sorted(cvg.items(), key=lambda kv: key_bin(kv[0][0], kv[0][1]))
    chunk_size = 100000
    for i in range(0, len(items), chunk_size):
        text = format_bdg_rows(
            items[i:i + chunk_size],
            siz_bin,
            dp,
            chrom_sizes
        )
        row_count += text.count("\n")
        digest.update(text.encode("utf-8"))

    return row_count, digest.hexdigest()


def write_bdg_array_sparse(
    signal_bins,
    fil_out: str,
    siz_bin: int,
    dp: int,
    chrom_sizes: dict[str, int] | None
) -> None:
    """
    Write a sparse-array signal directly to bedGraph.
    """
    with open_out(fil_out) as handle:
        for chrom, starts, values in iter_array_signal_chunks(signal_bins, 1):
            handle.write(format_bdg_array_rows(
                chrom,
                starts,
                values,
                siz_bin,
                dp,
                chrom_sizes.get(chrom) if chrom_sizes is not None else None
            ))


def write_bdg_parallel_ordered(
    cvg,
    fil_out: str,
    siz_bin: int,
    dp: int,
    chrom_sizes: dict[str, int] | None,
    workers: int
) -> None:
    """
    Format sorted bedGraph shards in worker processes and write in order.
    """
    check_cmp(workers, "ge", 1, "prototype_writer_workers", allow_none=False)
    if is_array_signal(cvg):
        if workers == 1:
            write_bdg_array_sparse(cvg, fil_out, siz_bin, dp, chrom_sizes)
            return

        task_data = [
            (
                chrom,
                starts,
                values,
                siz_bin,
                dp,
                chrom_sizes.get(chrom) if chrom_sizes is not None else None
            )
            for chrom, starts, values in iter_array_signal_chunks(
                cvg,
                workers
            )
        ]
        try:
            with ProcessPoolExecutor(
                max_workers=min(workers, os.cpu_count() or 1)
            ) as ex:
                formatted = list(ex.map(format_bdg_array_rows_task, task_data))
        except OSError:
            write_bdg_array_sparse(cvg, fil_out, siz_bin, dp, chrom_sizes)
            return

        with open_out(fil_out) as handle:
            for shard in formatted:
                handle.write(shard)
        return

    if workers == 1:
        write_bdg(cvg, fil_out, siz_bin, dp, chrom_sizes=chrom_sizes)
        return

    items = sorted(cvg.items(), key=lambda kv: key_bin(kv[0][0], kv[0][1]))
    if not items:
        with open_out(fil_out):
            return

    chunk_size = math.ceil(len(items) / workers)
    chunks = [
        items[i:i + chunk_size]
        for i in range(0, len(items), chunk_size)
    ]
    task_data = [
        (chunk, siz_bin, dp, chrom_sizes)
        for chunk in chunks
    ]

    try:
        with ProcessPoolExecutor(
            max_workers=min(workers, os.cpu_count() or 1)
        ) as ex:
            formatted = list(ex.map(format_bdg_rows_task, task_data))
    except OSError:
        write_bdg(cvg, fil_out, siz_bin, dp, chrom_sizes=chrom_sizes)
        return

    with open_out(fil_out) as handle:
        for shard in formatted:
            handle.write(shard)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments for computing binned signal from a BAM or CRAM
    file.

    Parameters
    ----------
        argv : list[str] | None
            Optional argument vector to parse. If None, use 'sys.argv[1:]'.

    Returns
    -------
        argparse.Namespace
            Parsed command-line arguments.

    Raises
    ------
        SystemExit
            Raised by argparse when '--help' is shown or when required
            arguments or valid choices are missing.

    Notes
    -----
        Parser option names, defaults, and accepted choices are documented in
        the 'add_argument()' definitions below and in rendered '--help' output.
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
            "(https://github.com/kalavattam/siQ-ChIP/tree/protocol)."
        )
    )
    add_help_cap(parser)
    parser.add_argument(
        "-v", "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Run script in verbose mode (stderr banner of parsed args).\n\n"
    )
    parser.add_argument(
        "-t", "--threads",
        dest="threads",
        type=int,
        default=1,
        help=(
            "Number of threads to use for parallel processing (>= 1; default: "
            "%(default)s).\n"
            "\n"
            "When '--threads > 1', different chromosomes are processed in "
            "parallel.\n\n"
        )
    )
    parser.add_argument(
        "-fi", "--fil_in",
        dest="fil_in",
        required=True,
        help="Input file path for the BAM or CRAM file.\n\n"
    )
    parser.add_argument(
        "-rf", "--ref_fa",
        dest="ref_fa",
        default=None,
        help=(
            "Reference FASTA file for CRAM decoding. Required for CRAM inputs "
            "unless the reference is otherwise available to htslib/pysam.\n\n"
        )
    )
    parser.add_argument(
        "-cs", "--chr_sizes", "--chrom_sizes",
        dest="chr_sizes",
        default=None,
        help=(
            "Chromosome sizes file in UCSC-style TSV format with chromosome "
            "name and positive integer size columns. Header sizes from "
            "BAM/CRAM are used when available; this file can supplement "
            "missing header sizes, but conflicting sizes are rejected.\n\n"
        )
    )
    parser.add_argument(
        "-fo", "--fil_out",
        dest="fil_out",
        required=True,
        help=(
            "Output file path. Supported output types are bedGraph "
            "('bedGraph', 'bedgraph', 'bdg', 'bg') and BED ('bed').\n"
            "\n"
            "Append '.gz' for gzip compression, e.g., 'output.bdg.gz'.\n"
            "\n"
            "Note: requesting BED output causes the script to write processed "
            "fragment coordinates in a BED-like format, and '--siz_bin', "
            "'--method', '--scl_fct', and '--dp' are ignored.\n\n"
        )
    )
    parser.add_argument(
        "-me", "--method",
        dest="method",
        choices=METHOD_CHOICES,
        default='norm',
        help=(
            "Workflow method. Specify signal calculation type (default: "
            "'%(default)s').\n"
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
        dest="siz_bin",
        type=int,
        default=10,
        help=(
            "Bin size in base pairs for signal calculation (default: "
            "%(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-eg", "--engine",
        dest="engine",
        choices=ENGINE_CHOICES,
        default="chrom",
        help=(
            "Processing engine for bedGraph output (default: %(default)s). "
            "'chrom' parallelizes indexed chromosome fetching and "
            "array-backed signal calculation. 'window' parallelizes indexed "
            "coordinate-window fetching for finer load balance on large BAM "
            "inputs.\n\n"
        )
    )
    parser.add_argument(
        "-ck", "--chunk_size", "--chunk-size",
        dest="chunk_size",
        type=int,
        default=100000,
        help=(
            "Number of records to process per chunk, retained for "
            "compatibility with older workflows (default: %(default)s). "
            "Ignored by public engines.\n\n"
        )
    )
    parser.add_argument(
        "-sf", "--scl_fct",
        dest="scl_fct",
        type=float,
        default=None,
        help=(
            "Scaling factor to apply to the signal (default: %(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-uf", "--usr_frg",
        dest="usr_frg",
        type=int,
        default=None,
        help=(
            "Fixed fragment length to use instead of read lengths (single-end "
            "alignments) or template lengths (paired-end alignments; default: "
            "%(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-dp", "--dp",
        dest="dp",
        type=int,
        default=24,
        help=(
            "Maximum number of decimal places retained for finite emitted "
            "values (default: %(default)s). After rounding, non-informative "
            "trailing zeros are stripped.\n\n"
        )
    )
    parser.add_argument(
        "-pj", "--profile_json",
        dest="profile_json",
        default=None,
        help=argparse.SUPPRESS
    )
    parser.add_argument(
        "-em", "--executor_mode",
        dest="executor_mode",
        choices=EXECUTOR_MODE_CHOICES,
        default="map",
        help=argparse.SUPPRESS
    )
    parser.add_argument(
        "-pps", "--prototype_parse_strategy",
        dest="prototype_parse_strategy",
        choices=PROTOTYPE_PARSE_STRATEGY_CHOICES,
        default=None,
        help=argparse.SUPPRESS
    )
    parser.add_argument(
        "-pbs", "--prototype_bed_strategy",
        dest="prototype_bed_strategy",
        choices=PROTOTYPE_BED_STRATEGY_CHOICES,
        default="auto",
        help=argparse.SUPPRESS
    )
    parser.add_argument(
        "-ws", "--indexed_window_size",
        dest="indexed_window_size",
        type=int,
        default=100000,
        help=argparse.SUPPRESS
    )
    parser.add_argument(
        "-prf", "--prototype_result_format",
        dest="prototype_result_format",
        choices=PROTOTYPE_RESULT_FORMAT_CHOICES,
        default=None,
        help=argparse.SUPPRESS
    )
    parser.add_argument(
        "-pms", "--prototype_merge_strategy",
        dest="prototype_merge_strategy",
        choices=PROTOTYPE_MERGE_STRATEGY_CHOICES,
        default=None,
        help=argparse.SUPPRESS
    )
    parser.add_argument(
        "-pws", "--prototype_writer_strategy",
        dest="prototype_writer_strategy",
        choices=PROTOTYPE_WRITER_STRATEGY_CHOICES,
        default=None,
        help=argparse.SUPPRESS
    )
    parser.add_argument(
        "-pww", "--prototype_writer_workers",
        dest="prototype_writer_workers",
        type=int,
        default=None,
        help=argparse.SUPPRESS
    )
    parser.add_argument(
        "-pwm", "--prototype_write_mode",
        dest="prototype_write_mode",
        choices=PROTOTYPE_WRITE_MODE_CHOICES,
        default="full",
        help=argparse.SUPPRESS
    )

    #  If no arguments are provided, display help and exit
    argv_parse = sys.argv[1:] if argv is None else argv
    if not argv_parse:
        parser.print_help(sys.stderr)
        raise SystemExit(0)

    return parser.parse_args(argv_parse)


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
        int.
            On success, returns 0 and writes either a bedGraph of binned signal
            or a BED-like file of processed fragment coordinates (optionally
            gzipped) to '--fil_out'.

    Side effects:
        - If '--verbose' is set, prints a human-readable banner of the parsed
          arguments to stderr.
        - Prints human-readable error messages to stderr on validation or I/O
          failures.

    Exits:
        - 0 on success or when showing help with no arguments.
        - 1 on validation or computation errors, including (non-exhaustive):
            + '--siz_bin <= 0' or '--usr_frg <= 0',
            + unsupported/invalid '--fil_out' extension,
            + normalization requested with zero fragments,
            + provided '--scl_fct <= 0', and
            + BAM or CRAM parsing errors (file not found/unreadable).
    """
    time_total_start = time.perf_counter()

    #  Parse CLI arguments
    args = parse_args(argv)

    if args.fil_in == "-":
        raise SystemExit(
            "'--fil_in -' is no longer supported; provide a BAM or CRAM path."
        )
    if args.fil_out == "-":
        raise SystemExit(
            "'--fil_out -' is no longer supported; provide an output path."
        )

    #  Check that input alignment/reference files exist
    try:
        check_exists(args.fil_in, kind="file", label="Alignment file")
        if args.ref_fa is not None:
            check_exists(args.ref_fa, kind="file", label="Reference FASTA")
        if args.chr_sizes is not None:
            check_exists(args.chr_sizes, kind="file", label="chr.sizes file")
    except FileNotFoundError as e:
        #  Print a one-line message and exit cleanly
        raise SystemExit(str(e)) from None

    #  Check and parse output file
    try:
        #  Parse output path and validate extension
        fil_out, out_fmt, _ = check_parse_fil_out(args.fil_out, ALLOWED)

        #  Check file creatability/writability
        check_writable(fil_out, "file")
        if args.profile_json is not None:
            check_writable(args.profile_json, "file")
    except (
        ValueError, FileNotFoundError, PermissionError, IsADirectoryError
    ) as e:
        raise SystemExit(str(e)) from None

    if out_fmt != "bed":
        args.prototype_parse_strategy = PUBLIC_ENGINE_STRATEGY[args.engine]
        if args.prototype_result_format is None:
            args.prototype_result_format = "direct_sparse_np"
        if args.prototype_merge_strategy is None:
            args.prototype_merge_strategy = "array_sparse_merge"
        if args.prototype_writer_strategy is None:
            args.prototype_writer_strategy = "parallel_ordered"
        if args.prototype_writer_workers is None:
            args.prototype_writer_workers = (
                1
                if args.prototype_writer_strategy == "serial"
                else min(4, os.cpu_count() or 1)
            )
    else:
        if args.prototype_parse_strategy is None:
            args.prototype_parse_strategy = "serial"
        if args.prototype_result_format is None:
            args.prototype_result_format = "dict"
        if args.prototype_merge_strategy is None:
            args.prototype_merge_strategy = "dict_merge"
        if args.prototype_writer_strategy is None:
            args.prototype_writer_strategy = "serial"
        if args.prototype_writer_workers is None:
            args.prototype_writer_workers = 1

    #  Check numeric arguments
    try:
        #  Check threads >= 1
        check_cmp(args.threads, "ge", 1, "threads", allow_none=False)

        #  Only require bin size and rounding for bedGraph paths
        if out_fmt != "bed":
            check_cmp(args.siz_bin, "gt", 0, "siz_bin", allow_none=False)
            check_cmp(
                args.chunk_size, "gt", 0, "chunk_size", allow_none=False
            )
            check_cmp(
                args.indexed_window_size,
                "gt",
                0,
                "indexed_window_size",
                allow_none=False
            )
            check_cmp(
                args.prototype_writer_workers,
                "ge",
                1,
                "prototype_writer_workers",
                allow_none=False
            )
            check_cmp(
                args.prototype_writer_workers,
                "le",
                4,
                "prototype_writer_workers",
                allow_none=False
            )
            check_cmp(args.dp, "ge", 0, "dp", allow_none=False)

            if (
                args.prototype_writer_strategy == "serial"
                and args.prototype_writer_workers != 1
            ):
                raise ValueError(
                    "'--prototype_writer_workers' must be 1 when "
                    "'--prototype_writer_strategy serial'."
                )
        else:
            check_cmp(
                args.indexed_window_size,
                "gt",
                0,
                "indexed_window_size",
                allow_none=False
            )

        #  Handle optional scalars: allow None but must be > 0 when provided
        check_cmp(args.scl_fct, "gt", 0, "scl_fct", allow_none=True)
        check_cmp(args.usr_frg, "gt", 0, "usr_frg", allow_none=True)

    except ValueError as e:
        raise SystemExit(str(e)) from None

    #  Preserve the user-supplied '--method' token, then standardize it to the
    #  canonical internal name
    mthd_in = args.method
    args.method = METHOD_CANON[args.method]

    #  Determine whether to normalize based on canonicalized method
    is_len = (args.method == "frag")
    is_norm = (args.method == "norm")
    profile = start_profile(args, fil_out, out_fmt)

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
            print(f"--fil_in  {args.fil_in}")
            print(f"--ref_fa  {args.ref_fa}")
            print(f"--chr_sizes {args.chr_sizes}")
            print(f"--fil_out {fil_out}")
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
                print(f"--engine  {args.engine}")
                print(f"--chunk_size {args.chunk_size}")
                print(f"--scl_fct {args.scl_fct}")
                print(f"--usr_frg {args.usr_frg}")
                print(f"--dp     {args.dp}")
            print("")
            print("")

    try:
        time_phase = time.perf_counter()
        header_sizes = get_alignment_chrom_sizes(args.fil_in, args.ref_fa)
        chrom_sizes = resolve_chrom_sizes(header_sizes, args.chr_sizes)
        record_phase(profile, "resolve_chrom_sizes", time_phase)
        if profile is not None:
            profile["n_chrom_sizes"] = len(chrom_sizes)
    except FileNotFoundError:
        raise SystemExit(
            f"Alignment file not found: {args.fil_in}"
        ) from None
    except ValueError as e:
        raise SystemExit(str(e)) from None
    except OSError as e:
        raise SystemExit(
            f"I/O error while reading BAM or CRAM: {e}"
        ) from None

    try:
        if out_fmt == 'bed':
            bed_strategy = args.prototype_bed_strategy
            if bed_strategy == "auto":
                if (
                    args.threads > 1
                    and has_indexed_alignment(args.fil_in, args.ref_fa)
                ):
                    bed_strategy = "indexed_chrom"
                else:
                    bed_strategy = "serial"
            if profile is not None:
                profile["prototype_bed_strategy_resolved"] = bed_strategy

            if bed_strategy == "serial":
                time_phase = time.perf_counter()
                frg_tup = parse_bam(
                    args.fil_in,
                    args.usr_frg,
                    args.ref_fa,
                    chrom_sizes=chrom_sizes
                )
                record_phase(profile, "parse_bed_fragments", time_phase)

                time_phase = time.perf_counter()
                row_count = write_bed_fragment_dict(frg_tup, fil_out)
                record_phase(profile, "write_bed", time_phase)
                if profile is not None:
                    profile["output_row_count"] = row_count
            else:
                time_phase = time.perf_counter()
                check_indexed_alignment(args.fil_in, args.ref_fa, bed_strategy)
                record_phase(
                    profile, "validate_indexed_bed_alignment", time_phase
                )

                time_phase = time.perf_counter()
                chrom_task_names = [
                    chrom
                    for chrom in header_sizes
                    if chrom in chrom_sizes
                ]
                if bed_strategy == "indexed_chrom":
                    tsk_dat = [
                        (
                            args.fil_in,
                            chrom,
                            chrom_sizes,
                            args.usr_frg,
                            args.ref_fa,
                            None,
                            None
                        )
                        for chrom in chrom_task_names
                    ]
                else:
                    tsk_dat = [
                        (
                            args.fil_in,
                            chrom,
                            chrom_sizes,
                            args.usr_frg,
                            args.ref_fa,
                            start,
                            min(
                                start + args.indexed_window_size,
                                chrom_sizes[chrom]
                            )
                        )
                        for chrom in chrom_task_names
                        for start in range(
                            0,
                            chrom_sizes[chrom],
                            args.indexed_window_size
                        )
                    ]
                if profile is not None:
                    profile["tasks"] = [
                        {
                            "task_id": task_id,
                            "kind": "indexed_bed",
                            "bed_strategy": bed_strategy,
                            "chrom": task[1],
                            "start": task[5],
                            "end": task[6],
                            "chrom_size": int(chrom_sizes[task[1]]),
                        }
                        for task_id, task in enumerate(tsk_dat)
                    ]
                    profile["n_tasks"] = len(tsk_dat)
                record_phase(profile, "build_indexed_bed_tasks", time_phase)

                time_phase = time.perf_counter()
                try:
                    result_list = list(iter_task_results(
                        "indexed_bed",
                        collect_bed_indexed_task,
                        tsk_dat,
                        args.threads,
                        args.executor_mode,
                        profile["tasks"] if profile is not None else None
                    ))
                except OSError:
                    if args.prototype_bed_strategy != "auto":
                        raise
                    if profile is not None:
                        profile["prototype_bed_strategy_resolved"] = (
                            "serial_fallback"
                        )
                    frg_tup = parse_bam(
                        args.fil_in,
                        args.usr_frg,
                        args.ref_fa,
                        chrom_sizes=chrom_sizes
                    )
                    record_phase(profile, "parse_bed_fragments", time_phase)
                    time_phase = time.perf_counter()
                    row_count = write_bed_fragment_dict(frg_tup, fil_out)
                    record_phase(profile, "write_bed", time_phase)
                    if profile is not None:
                        profile["output_row_count"] = row_count
                        profile["phases_s"]["total"] = (
                            time.perf_counter() - time_total_start
                        )
                        write_profile(args.profile_json, profile)
                    return 0
                record_phase(
                    profile, "receive_indexed_bed_results", time_phase
                )
                summarize_task_profiles(profile)
                if profile is not None:
                    profile["result_payload_bytes"] = sum(
                        estimate_bed_payload_bytes(result)
                        for result in result_list
                    )

                time_phase = time.perf_counter()
                row_count = write_bed_array_results(result_list, fil_out)
                record_phase(
                    profile, "merge_sort_write_indexed_bed", time_phase
                )
                if profile is not None:
                    profile["output_row_count"] = row_count

            if profile is not None:
                profile["phases_s"]["total"] = (
                    time.perf_counter() - time_total_start
                )
                write_profile(args.profile_json, profile)
            return 0

        #  Otherwise, compute and write bedGraph signal.
        if args.engine in PUBLIC_ENGINE_STRATEGY:
            strategy = PUBLIC_ENGINE_STRATEGY[args.engine]
            is_prototype_strategy = strategy != "serial"

            if strategy == "serial":
                time_phase = time.perf_counter()
                frg_arr, frg_tot = collect_fragment_arrays(
                    args.fil_in,
                    chrom_sizes=chrom_sizes,
                    usr_frg=args.usr_frg,
                    ref_fa=args.ref_fa
                )
                record_phase(profile, "parse_fragment_arrays", time_phase)

                if is_norm and frg_tot == 0:
                    raise ValueError(
                        "Normalization requires non-zero total fragments. "
                        "Check the BAM or CRAM file."
                    )

                time_phase = time.perf_counter()
                tsk_dat = [
                    (
                        chrom,
                        arrays[0],
                        arrays[1],
                        arrays[2],
                        chrom_sizes[chrom],
                        frg_tot,
                        args.siz_bin,
                        is_len,
                        is_norm,
                        args.scl_fct
                    )
                    for chrom, arrays in frg_arr.items()
                ]
                if profile is not None:
                    profile["fragments_total"] = frg_tot
                    profile["tasks"] = [
                        {
                            "task_id": task_id,
                            "kind": "chrom",
                            "chrom": task[0],
                            "fragments": int(task[1].size),
                            "chrom_size": int(task[4]),
                            "n_bins": math.ceil(task[4] / args.siz_bin),
                            "payload_bytes": int(
                                task[1].nbytes
                                + task[2].nbytes
                                + task[3].nbytes
                            ),
                        }
                        for task_id, task in enumerate(tsk_dat)
                    ]
                    profile["n_tasks"] = len(tsk_dat)
                record_phase(profile, "build_tasks", time_phase)

                results = iter_task_results(
                    "chrom",
                    calc_sig_array_task,
                    tsk_dat,
                    args.threads,
                    args.executor_mode,
                    profile["tasks"] if profile is not None else None
                )

            else:
                time_phase = time.perf_counter()
                check_indexed_alignment(args.fil_in, args.ref_fa, strategy)
                record_phase(profile, "validate_indexed_alignment", time_phase)

                time_phase = time.perf_counter()
                chrom_task_names = [
                    chrom
                    for chrom in header_sizes
                    if chrom in chrom_sizes
                ]

                if strategy == "indexed_chrom":
                    tsk_dat = [
                        (
                            args.fil_in,
                            chrom,
                            chrom_sizes,
                            args.usr_frg,
                            args.ref_fa,
                            None,
                            None,
                            args.siz_bin,
                            is_len or is_norm,
                            args.prototype_result_format
                        )
                        for chrom in chrom_task_names
                    ]
                else:
                    tsk_dat = [
                        (
                            args.fil_in,
                            chrom,
                            chrom_sizes,
                            args.usr_frg,
                            args.ref_fa,
                            start,
                            min(
                                start + args.indexed_window_size,
                                chrom_sizes[chrom]
                            ),
                            args.siz_bin,
                            is_len or is_norm,
                            args.prototype_result_format
                        )
                        for chrom in chrom_task_names
                        for start in range(
                            0,
                            chrom_sizes[chrom],
                            args.indexed_window_size
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
                            "chrom_size": int(chrom_sizes[task[1]]),
                            "n_bins": math.ceil(
                                chrom_sizes[task[1]] / args.siz_bin
                            ),
                        }
                        for task_id, task in enumerate(tsk_dat)
                    ]
                    profile["n_tasks"] = len(tsk_dat)
                record_phase(profile, "build_tasks", time_phase)

                results = iter_task_results(
                    strategy,
                    calc_sig_indexed_fetch_task,
                    tsk_dat,
                    args.threads,
                    args.executor_mode,
                    profile["tasks"] if profile is not None else None
                )

        else:
            is_prototype_strategy = False
            if is_norm:
                time_phase = time.perf_counter()
                frg_tot = count_fragments(
                    args.fil_in,
                    chrom_sizes=chrom_sizes,
                    usr_frg=args.usr_frg,
                    ref_fa=args.ref_fa
                )
                record_phase(profile, "count_fragments", time_phase)
            else:
                frg_tot = 0

            if is_norm and frg_tot == 0:
                raise ValueError(
                    "Normalization requires non-zero total fragments. Check "
                    "the BAM or CRAM file."
                )

            #  Prepare fixed-size fragment chunks for signal computation
            time_phase = time.perf_counter()
            if profile is None:
                tsk_dat = (
                    (
                        chunk,
                        frg_tot,
                        args.siz_bin,
                        is_len,
                        is_norm,
                        args.scl_fct
                    )
                    for chunk in iter_fragment_chunks(
                        path_bam=args.fil_in,
                        chrom_sizes=chrom_sizes,
                        chunk_size=args.chunk_size,
                        usr_frg=args.usr_frg,
                        ref_fa=args.ref_fa
                    )
                )
            else:
                chunks = list(iter_fragment_chunks(
                    path_bam=args.fil_in,
                    chrom_sizes=chrom_sizes,
                    chunk_size=args.chunk_size,
                    usr_frg=args.usr_frg,
                    ref_fa=args.ref_fa
                ))
                tsk_dat = [
                    (
                        chunk,
                        frg_tot,
                        args.siz_bin,
                        is_len,
                        is_norm,
                        args.scl_fct
                    )
                    for chunk in chunks
                ]
                profile["fragments_total"] = (
                    frg_tot if is_norm else sum(len(chunk) for chunk in chunks)
                )
                profile["tasks"] = [
                    {
                        "task_id": task_id,
                        "kind": "chunk",
                        "fragments": len(chunk),
                        "payload_items": len(chunk),
                    }
                    for task_id, chunk in enumerate(chunks)
                ]
                profile["n_tasks"] = len(tsk_dat)
            record_phase(profile, "build_tasks", time_phase)

            #  Compute per-chunk signal, in parallel if '--threads > 1'
            results = iter_task_results(
                "chunk",
                calc_sig_chunk_task,
                tsk_dat,
                args.threads,
                args.executor_mode,
                profile["tasks"] if profile is not None else None
            )

        #  Receive worker results, then merge in the parent process.
        time_collect_merge = time.perf_counter()
        time_phase = time.perf_counter()
        result_list = list(results)
        record_phase(profile, "receive_worker_results", time_phase)
        summarize_task_profiles(profile)

        time_phase = time.perf_counter()
        sig_cmb, frg_tot_merged, n_result_bins, result_payload_bytes = (
            merge_signal_results(
                results=result_list,
                merge_strategy=args.prototype_merge_strategy,
                chrom_sizes=chrom_sizes,
                siz_bin=args.siz_bin,
                is_prototype_strategy=is_prototype_strategy,
                profile=profile
            )
        )
        record_phase(profile, "parent_merge_results", time_phase)
        if is_prototype_strategy:
            frg_tot = frg_tot_merged

        if is_prototype_strategy:
            time_phase = time.perf_counter()
            apply_signal_adjustments(
                sig_cmb,
                frg_tot,
                is_norm,
                args.scl_fct
            )
            record_phase(profile, "apply_signal_adjustments", time_phase)
        record_phase(profile, "collect_and_merge_results", time_collect_merge)
        if profile is not None:
            if is_prototype_strategy:
                profile["fragments_total"] = frg_tot
            profile["result_bins_before_merge"] = n_result_bins
            profile["result_bins_after_merge"] = signal_row_count(sig_cmb)
            profile["result_payload_bytes"] = result_payload_bytes

        #  Write, digest, or skip bedGraph output for private benchmark modes.
        time_phase = time.perf_counter()
        if args.prototype_write_mode == "profile_only":
            if profile is not None:
                profile["output_written"] = False
                profile["output_digest_mode"] = False
                profile["output_row_count"] = signal_row_count(sig_cmb)
                profile["output_sha256"] = None
        elif args.prototype_write_mode == "digest":
            row_count, output_sha256 = digest_bdg_rows(
                sig_cmb,
                args.siz_bin,
                args.dp,
                chrom_sizes
            )
            if profile is not None:
                profile["output_written"] = False
                profile["output_digest_mode"] = True
                profile["output_row_count"] = row_count
                profile["output_sha256"] = output_sha256
        elif args.prototype_writer_strategy == "parallel_ordered":
            write_bdg_parallel_ordered(
                sig_cmb,
                fil_out,
                args.siz_bin,
                args.dp,
                chrom_sizes,
                args.prototype_writer_workers
            )
        elif is_array_signal(sig_cmb):
            write_bdg_array_sparse(
                sig_cmb,
                fil_out,
                args.siz_bin,
                args.dp,
                chrom_sizes
            )
        else:
            write_bdg(
                sig_cmb,
                fil_out,
                args.siz_bin,
                args.dp,
                chrom_sizes=chrom_sizes
            )
        if profile is not None and args.prototype_write_mode == "full":
            profile["output_written"] = True
            profile["output_digest_mode"] = False
            profile["output_row_count"] = signal_row_count(sig_cmb)
            profile["output_sha256"] = None
        record_phase(profile, "write_bedgraph", time_phase)
        if profile is not None:
            profile["phases_s"]["total"] = (
                time.perf_counter() - time_total_start
            )
            write_profile(args.profile_json, profile)

        return 0

    except ValueError as e:
        raise SystemExit(str(e)) from None

    except OSError as e:
        raise SystemExit(f"I/O error: {e}") from None


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        with suppress(Exception):
            sys.stdout.close()
        with suppress(Exception):
            sys.stderr.close()
        raise SystemExit(0) from None
