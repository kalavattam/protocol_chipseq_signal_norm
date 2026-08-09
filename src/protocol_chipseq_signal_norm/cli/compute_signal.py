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
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


"""
Calculate binned signal or fragment-coordinate output from BAM/CRAM input.

The CLI accepts input, output, signal-mode, engine, scaling, and formatting
options. It writes bedGraph-like signal tracks or BED-like fragment-coordinate
records. When writing bedGraph, finite values are rounded to at most '--dp'
decimal places and trailing zeros are stripped.

Examples
--------
python -m protocol_chipseq_signal_norm.cli.compute_signal \\
    --fil_in <file> --fil_out <file> [options]
"""

from __future__ import annotations

import argparse
import hashlib
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

from protocol_chipseq_signal_norm.utilities.utils_bdg import (
    key_bin,
    load_chromosome_sizes,
    write_bdg,
)
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

# TODO: Reassess chunked processing, compiled-core options, and downstream
# bigWig compatibility after the current indexed engines are benchmarked.

# Map accepted `--method` values to canonical internal names.
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
    "indexed_window",
)
PROTOTYPE_BED_STRATEGY_CHOICES = (
    "auto",
    "serial",
    "indexed_chrom",
    "indexed_window",
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
    "event_np",
)
PROTOTYPE_MERGE_STRATEGY_CHOICES = (
    "dict_merge",
    "chrom_array_merge",
    "vectorized_merge",
    "array_sparse_merge_legacy",
    "array_sparse_merge",
    "array_dense_merge",
    "event_diff_merge",
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
    "direct_sparse_local_bincount_np",
)


def get_alignment_chrom_sizes(
    alignment_path: str,
    ref_fa: str | None = None,
) -> dict[str, int]:
    """
    Read chromosome sizes from a BAM or CRAM header.
    """

    kwargs = {}

    if ref_fa is not None:
        kwargs["reference_filename"] = ref_fa

    with pysam.AlignmentFile(alignment_path, "rb", **kwargs) as alignment_file:
        return {
            chrom: size
            for chrom, size in zip(
                alignment_file.references,
                alignment_file.lengths,
                strict=True,
            )
            if size is not None and size > 0
        }


def resolve_chrom_sizes(
    header_sizes: dict[str, int],
    chromosome_sizes_path: str | None = None,
) -> dict[str, int]:
    """
    Resolve chromosome sizes from BAM/CRAM headers plus an optional TSV.
    """

    chrom_sizes = dict(header_sizes)

    if chromosome_sizes_path is not None:
        file_sizes = load_chromosome_sizes(chromosome_sizes_path)
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
                f"but chr.sizes file has {file_size}.",
            )

        chrom_sizes.update(file_sizes)

    if not chrom_sizes:
        raise ValueError(
            "Chromosome sizes are required to trim bedGraph bins. Provide a "
            "BAM/CRAM with sequence lengths in the header or pass "
            "'--chr_sizes' with a UCSC-style two-column TSV.",
        )

    return chrom_sizes


def start_profile(
    args: argparse.Namespace,
    output_path: str,
    output_format: str,
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
        "fil_out": output_path,
        "out_fmt": output_format,
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
    siz_bin: int,
) -> object:
    """
    Convert a signal dict to one private benchmark result representation.

    Parameters
    ----------
    sig : dict[tuple[str, int], float]
        Sparse signal keyed by chromosome and bin start.
    result_format : str
        Requested private result representation.
    siz_bin : int
        Bin width used to derive array coordinates.

    Returns
    -------
    result : object
        Signal in the requested private benchmark representation.

    Raises
    ------
    ValueError
        If 'result_format' is not recognized.
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
                    np.asarray([value for _, value in rows], dtype=np.float64),
                )
                for chrom, rows in by_chrom.items()
            ],
        )

    if result_format == "dense_np":
        dense_parts = []

        for chrom, rows in by_chrom.items():
            starts = [start for start, _ in rows]
            first = min(starts)
            last = max(starts)
            values = np.zeros(
                ((last - first) // siz_bin) + 1,
                dtype=np.float64,
            )

            for start, value in rows:
                values[(start - first) // siz_bin] += value

            dense_parts.append((chrom, first, values))

        return "dense_np", dense_parts

    raise ValueError(f"Unknown prototype result format: {result_format!r}.")


def count_result_bins(result: Any) -> int:
    """
    Count signal bins for standard dict results and prototype tuple results.

    Parameters
    ----------
    result : Any
        Supported mapping- or tuple-backed signal result.

    Returns
    -------
    count : int
        Number of represented signal bins.
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
                    _,
                ) in result[1]
            )

        if tag == "dense_np":
            return sum(
                int(np.count_nonzero(values)) for _, _, values in result[1]
            )

        if isinstance(tag, dict):
            return len(tag)

    return len(result)


def estimate_result_payload_bytes(result: Any) -> int:
    """
    Estimate NumPy-array payload bytes for private benchmark result formats.

    Parameters
    ----------
    result : Any
        Supported mapping- or array-backed signal result.

    Returns
    -------
    byte_count : int
        Sum of owned NumPy array payload bytes.
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
                    + touch_values.nbytes,
                )
                for (
                    _,
                    _,
                    edge_bins,
                    edge_values,
                    diff_bins,
                    diff_values,
                    touch_bins,
                    touch_values,
                ) in result[1]
            )

        if tag == "dense_np":
            return sum(int(values.nbytes) for _, _, values in result[1])

        if isinstance(tag, dict):
            return 0

    return 0


def check_indexed_alignment(
    alignment_path: str,
    ref_fa: str | None,
    strategy: str,
) -> None:
    """
    Validate prerequisites for indexed signal engines.

    Parameters
    ----------
    alignment_path : str
        BAM or CRAM path that must have a corresponding index.
    ref_fa : str | None
        Reference FASTA required for indexed CRAM access.
    strategy : str
        Selected signal-engine strategy.

    Raises
    ------
    ValueError
        If the path, index, or required CRAM reference is absent.
    """

    if strategy == "serial":
        return

    if alignment_path == "-":
        raise ValueError(
            "Indexed signal engines require a named BAM or CRAM path.",
        )

    if alignment_path.lower().endswith(".cram") and ref_fa is None:
        raise ValueError(
            "Indexed CRAM signal engines require '--ref_fa' so every worker "
            "can open the CRAM deterministically.",
        )

    kwargs = {}

    if ref_fa is not None:
        kwargs["reference_filename"] = ref_fa

    with pysam.AlignmentFile(alignment_path, "rb", **kwargs) as alignment_file:
        if not alignment_file.has_index():
            raise ValueError(
                "Indexed signal engines require an alignment index "
                "(.bai for BAM or .crai for CRAM).",
            )


def has_indexed_alignment(
    alignment_path: str,
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
            alignment_path,
            "rb",
            **kwargs,
        ) as alignment_file:
            return bool(alignment_file.has_index())
    except (OSError, ValueError):
        return False


def keep_read(
    read: pysam.AlignedSegment,
    allow_secondary: bool = True,
    allow_supplementary: bool = False,
    allow_duplicates: bool = True,
) -> bool:
    """
    Return whether a read passes the shared alignment-filter policy.
    """

    if read.is_unmapped or read.reference_id < 0:
        return False

    if not allow_secondary and read.is_secondary:
        return False

    if not allow_supplementary and read.is_supplementary:
        return False

    return not (not allow_duplicates and read.is_duplicate)


def read_to_fragment(
    read: pysam.AlignedSegment,
    get_reference_name: Callable[[int], str],
    chrom_sizes: dict[str, int],
    user_fragment_length: int | None = None,
) -> tuple[str, int, int, int] | None:
    """
    Convert one accepted alignment record to a processed fragment interval.

    Parameters
    ----------
    read : pysam.AlignedSegment
        Accepted alignment record.
    get_reference_name : Callable[[int], str]
        Resolver for numeric reference identifiers.
    chrom_sizes : dict[str, int]
        Maximum coordinate for each accepted chromosome.
    user_fragment_length : int | None
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
    chrom_len = chrom_sizes.get(chrom)
    if chrom_len is None:
        raise ValueError(
            f"Chromosome {chrom!r} is missing from chromosome sizes. "
            "Provide '--chr_sizes' with a UCSC-style two-column TSV.",
        )

    # Handle paired-end alignments: one fragment is emitted per leftmost
    # anchor in a proper pair.
    is_leftmost_pe = (
        read.is_paired
        and read.is_proper_pair
        and read.reference_id == read.next_reference_id
        and read.template_length > 0
    )

    if is_leftmost_pe:
        start = read.reference_start
        tlen = read.template_length
        fragment_length = (
            user_fragment_length if user_fragment_length is not None else tlen
        )

        if user_fragment_length is not None and fragment_length <= 0:
            raise ValueError("usr_frg must be > 0 for paired-end extension.")

        fragment_start = start
        fragment_end = start + fragment_length

    elif not read.is_paired:
        fragment_length = (
            user_fragment_length
            if user_fragment_length is not None
            else read.query_alignment_length
        )

        if user_fragment_length is not None and fragment_length <= 0:
            raise ValueError("usr_frg must be > 0 for single-end extension.")

        if user_fragment_length is None and fragment_length <= 0:
            return None

        if read.is_reverse:
            ref_end = (
                read.reference_end
                if read.reference_end is not None
                else read.reference_start
            )
            fragment_end = ref_end
            fragment_start = fragment_end - fragment_length
        else:
            fragment_start = read.reference_start
            fragment_end = fragment_start + fragment_length
    else:
        return None

    if fragment_start < 0:
        fragment_start = 0

    if fragment_end > chrom_len:
        fragment_end = chrom_len

    if fragment_end <= fragment_start:
        return None

    return chrom, fragment_start, fragment_end, fragment_length


def iter_alignment_fragments(
    alignment_path: str,
    chrom_sizes: dict[str, int],
    user_fragment_length: int | None = None,
    ref_fa: str | None = None,
    allow_secondary: bool = True,
    allow_supplementary: bool = False,
    allow_duplicates: bool = True,
) -> Iterator[tuple[str, int, int, int]]:
    """
    Stream processed fragment intervals from a BAM or CRAM file.

    Parameters
    ----------
    alignment_path : str
        BAM or CRAM input path.
    chrom_sizes : dict[str, int]
        Maximum coordinate for each accepted chromosome.
    user_fragment_length : int | None
        Optional fragment length used to extend alignments.
    ref_fa : str | None
        Reference FASTA used for CRAM decoding.
    allow_secondary, allow_supplementary, allow_duplicates : bool
        Whether to retain secondary, supplementary, and duplicate records.

    Yields
    ------
    fragment : tuple[str, int, int, int]
        Chromosome, clipped start, clipped end, and fragment length.
    """

    kwargs = {}

    if ref_fa is not None:
        kwargs["reference_filename"] = ref_fa

    with pysam.AlignmentFile(alignment_path, "rb", **kwargs) as alignment_file:
        for read in alignment_file.fetch(until_eof=True):
            if not keep_read(
                read,
                allow_secondary,
                allow_supplementary,
                allow_duplicates,
            ):
                continue

            fragment = read_to_fragment(
                read,
                alignment_file.get_reference_name,
                chrom_sizes,
                user_fragment_length,
            )
            if fragment is not None:
                yield fragment


def iter_indexed_fragments(
    alignment_path: str,
    chrom: str,
    chrom_sizes: dict[str, int],
    user_fragment_length: int | None = None,
    ref_fa: str | None = None,
    start: int | None = None,
    end: int | None = None,
    allow_secondary: bool = True,
    allow_supplementary: bool = False,
    allow_duplicates: bool = True,
) -> Iterator[tuple[str, int, int, int]]:
    """
    Stream processed fragments from one indexed chromosome or window fetch.

    Parameters
    ----------
    alignment_path : str
        Indexed BAM or CRAM input path.
    chrom : str
        Chromosome to fetch.
    chrom_sizes : dict[str, int]
        Maximum coordinate for each accepted chromosome.
    user_fragment_length : int | None
        Optional fragment length used to extend alignments.
    ref_fa : str | None
        Reference FASTA used for CRAM decoding.
    start, end : int | None
        Optional half-open fetch bounds.
    allow_secondary, allow_supplementary, allow_duplicates : bool
        Whether to retain secondary, supplementary, and duplicate records.

    Yields
    ------
    fragment : tuple[str, int, int, int]
        Chromosome, clipped start, clipped end, and fragment length.
    """

    kwargs = {}

    if ref_fa is not None:
        kwargs["reference_filename"] = ref_fa

    with pysam.AlignmentFile(alignment_path, "rb", **kwargs) as alignment_file:
        if start is None and end is None:
            reads = alignment_file.fetch(chrom)
        else:
            reads = alignment_file.fetch(chrom, start, end)

        for read in reads:
            if not keep_read(
                read,
                allow_secondary,
                allow_supplementary,
                allow_duplicates,
            ):
                continue

            if (
                start is not None
                and end is not None
                and (
                    read.reference_start < start or read.reference_start >= end
                )
            ):
                # Partition by read anchor, not computed fragment start.
                # This prevents double-counting while preserving reverse-
                # strand single-end fragments that extend left of a window.
                continue

            fragment = read_to_fragment(
                read,
                alignment_file.get_reference_name,
                chrom_sizes,
                user_fragment_length,
            )
            if fragment is not None:
                yield fragment


def collect_fragment_arrays_from_iter(
    fragments: Iterator[tuple[str, int, int, int]],
) -> tuple[dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]], int]:
    """
    Collect processed fragments into per-chromosome NumPy arrays.
    """

    starts: dict[str, list[int]] = defaultdict(list)
    ends: dict[str, list[int]] = defaultdict(list)
    lengths: dict[str, list[int]] = defaultdict(list)
    fragment_count = 0

    for chrom, fragment_start, fragment_end, fragment_length in fragments:
        starts[chrom].append(fragment_start)
        ends[chrom].append(fragment_end)
        lengths[chrom].append(fragment_length)
        fragment_count += 1

    arrays: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    for chrom in starts:
        arrays[chrom] = (
            np.asarray(starts[chrom], dtype=np.int64),
            np.asarray(ends[chrom], dtype=np.int64),
            np.asarray(lengths[chrom], dtype=np.float64),
        )

    return arrays, fragment_count


def collect_bed_arrays_from_iter(
    fragments: Iterator[tuple[str, int, int, int]],
) -> tuple[str, list[tuple[str, np.ndarray, np.ndarray, np.ndarray]], int]:
    """
    Collect processed fragments into per-chromosome BED output arrays.
    """

    starts: dict[str, list[int]] = defaultdict(list)
    ends: dict[str, list[int]] = defaultdict(list)
    lengths: dict[str, list[int]] = defaultdict(list)
    n_fragments = 0

    for chrom, fragment_start, fragment_end, fragment_length in fragments:
        starts[chrom].append(fragment_start)
        ends[chrom].append(fragment_end)
        lengths[chrom].append(fragment_length)
        n_fragments += 1

    arrays = [
        (
            chrom,
            np.asarray(starts[chrom], dtype=np.int64),
            np.asarray(ends[chrom], dtype=np.int64),
            np.asarray(lengths[chrom], dtype=np.int64),
        )
        for chrom in starts
    ]

    return "bed_np", arrays, n_fragments


def iter_fragment_chunks(
    alignment_path: str,
    chrom_sizes: dict[str, int],
    chunk_size: int,
    user_fragment_length: int | None = None,
    ref_fa: str | None = None,
) -> Iterator[list[tuple[str, int, int, int]]]:
    """
    Stream processed fragments in fixed-size chunks.
    """

    validate_comparison(chunk_size, "gt", 0, "chunk_size", allow_none=False)

    chunk: list[tuple[str, int, int, int]] = []

    for fragment in iter_alignment_fragments(
        alignment_path=alignment_path,
        chrom_sizes=chrom_sizes,
        user_fragment_length=user_fragment_length,
        ref_fa=ref_fa,
    ):
        chunk.append(fragment)

        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []

    if chunk:
        yield chunk


def count_fragments(
    alignment_path: str,
    chrom_sizes: dict[str, int],
    user_fragment_length: int | None = None,
    ref_fa: str | None = None,
) -> int:
    """
    Count processed fragments after filtering and coordinate clamping.
    """

    return sum(
        1
        for _ in iter_alignment_fragments(
            alignment_path=alignment_path,
            chrom_sizes=chrom_sizes,
            user_fragment_length=user_fragment_length,
            ref_fa=ref_fa,
        )
    )


def parse_bam(
    alignment_path: str,
    user_fragment_length: int | None = None,
    ref_fa: str | None = None,
    chrom_sizes: dict[str, int] | None = None,
    allow_secondary: bool = True,
    allow_supplementary: bool = False,
    allow_duplicates: bool = True,
) -> dict[str, list[tuple[int, int, int]]]:
    """
    Parse a BAM or CRAM into chromosome-grouped fragment intervals.

    Parameters
    ----------
    alignment_path : str
        Input BAM or CRAM path.
    user_fragment_length : int | None
        Optional fixed fragment length override. If provided:
            - paired-end read alignments: overrides TLEN for leftmost anchors.
            - single-end read alignments: used for both strands from the 5’
                                          end.
    ref_fa : str | None
        Reference FASTA file for CRAM decoding.
    chrom_sizes : dict[str, int] | None
        Optional chromosome sizes used to filter and clamp fragments.
    allow_secondary : bool
        Include secondary (multi-mapping) alignments. Default True.
    allow_supplementary : bool
        Include supplementary alignments. Default False.
    allow_duplicates : bool
        Include duplicate-marked alignments. Default True.

    Returns
    -------
    fragment_records : dict[str, list[tuple[int, int, int]]]
        Chromosome-keyed lists of start, end, and fragment-length tuples. Every
        interval is half-open and uses zero-based coordinates.

    Raises
    ------
    FileNotFoundError
        If 'alignment_path' does not exist.
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
            + 'allow_secondary': Include or exclude secondary alignments.
            + 'allow_supplementary': Include or exclude supplementary
              alignments.
            + 'allow_duplicates': Include or exclude duplicate-marked
              alignments.

            These toggles apply uniformly to both paired-end and single-end
            alignments. Keeping 'allow_duplicates=True' retains
            duplicate-marked proper-pair alignments (e.g., FLAGs 1123 and 1187,
            which correspond to duplicate-marked versions of 99 and 163). If
            the current data lack secondary alignments (which is expected in
            alignment output from the Tsukiyama Lab Bio-protocol workflow),
            setting 'allow_secondary=True' does nothing.

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

    fragment_records = defaultdict(list)

    try:
        if chrom_sizes is None:
            header_sizes = get_alignment_chrom_sizes(alignment_path, ref_fa)
            chrom_sizes = resolve_chrom_sizes(header_sizes)

        for (
            chrom,
            fragment_start,
            fragment_end,
            fragment_length,
        ) in iter_alignment_fragments(
            alignment_path=alignment_path,
            chrom_sizes=chrom_sizes,
            user_fragment_length=user_fragment_length,
            ref_fa=ref_fa,
            allow_secondary=allow_secondary,
            allow_supplementary=allow_supplementary,
            allow_duplicates=allow_duplicates,
        ):
            fragment_records[chrom].append(
                (fragment_start, fragment_end, fragment_length),
            )

    except FileNotFoundError:
        print(
            f"Error: BAM or CRAM file '{alignment_path}' not found.",
            file=sys.stderr,
        )

        raise
    except ValueError:
        raise
    except Exception as error:
        print(
            f"Unexpected error with BAM or CRAM file '{alignment_path}': "
            f"{error}",
            file=sys.stderr,
        )

        raise

    return fragment_records


def collect_fragment_arrays(
    alignment_path: str,
    chrom_sizes: dict[str, int],
    user_fragment_length: int | None = None,
    ref_fa: str | None = None,
    allow_secondary: bool = True,
    allow_supplementary: bool = False,
    allow_duplicates: bool = True,
) -> tuple[dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]], int]:
    """
    Parse a BAM or CRAM once into per-chromosome NumPy fragment arrays.
    """

    return collect_fragment_arrays_from_iter(
        iter_alignment_fragments(
            alignment_path=alignment_path,
            chrom_sizes=chrom_sizes,
            user_fragment_length=user_fragment_length,
            ref_fa=ref_fa,
            allow_secondary=allow_secondary,
            allow_supplementary=allow_supplementary,
            allow_duplicates=allow_duplicates,
        ),
    )


def calc_sig_chrom(
    chrom: str,
    fragment_records: list[tuple[int, int, int]],
    fragment_count: int,
    siz_bin: int,
    is_len: bool,
    is_norm: bool,
    scl_fct: float | None = None,
) -> dict[tuple[str, int], float]:
    """
    Compute one chromosome of exact fragment-bin overlap signal.

    Function respects half-open '[start, end)' fragments and avoids per-base
    loops. Overlap is measured in bases, and output is per-bin sums of per-base
    contributions.

    If provided, 'scl_fct' is applied after any optional normalization.
    ('scl_fct > 0' is required to avoid silent zeroing.)

    Parameters
    ----------
    chrom : str
        Chromosome name.
    fragment_records : list[tuple[int, int, int]]
        Start, end, and fragment-length tuples.
    fragment_count : int
        Total number of fragments (used when 'is_norm=True').
    siz_bin : int
        Bin size in base pairs.
    is_len : bool
        If 'True', normalize by fragment length.
    is_norm : bool
        If 'True', normalize by both fragment length and total fragment count
        so that the genome-wide summed signal is approximately 1.
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
        - If 'is_len' or 'is_norm' is True and any fragment length is
          nonpositive.
        - If 'is_norm' is True and 'fragment_count' <= 0.
        - If a provided 'scl_fct' <= 0.

    Notes
    -----
    - If 'is_len=True' or 'is_norm=True', each fragment contributes one divided
      by the fragment length per covered base.
    - If 'is_norm=True', signal is additionally divided by 'fragment_count' so
      that genome-wide summed signal is approximately 1.
    - If 'scl_fct' is provided, signal values are scaled accordingly.
    - Signal is accumulated per bin according to the number of overlapping
      bases contributed by each fragment.
    """

    validate_comparison(siz_bin, "gt", 0, "siz_bin", allow_none=False)

    sig_bin = defaultdict(float)

    for fragment_start, fragment_end, fragment_length in fragment_records:
        if fragment_end <= fragment_start:
            continue

        if (is_len or is_norm) and fragment_length <= 0:
            raise ValueError("'frg_len' must be > 0 when using normalization.")

        per_bas = (1.0 / fragment_length) if (is_len or is_norm) else 1.0

        first_bin_start = (fragment_start // siz_bin) * siz_bin
        last_bin_start = ((fragment_end - 1) // siz_bin) * siz_bin

        for bin_start in range(first_bin_start, last_bin_start + 1, siz_bin):
            bin_end = bin_start + siz_bin
            overlap_start = max(fragment_start, bin_start)
            overlap_end = min(fragment_end, bin_end)
            overlap = overlap_end - overlap_start

            if overlap > 0:
                sig_bin[(chrom, bin_start)] += per_bas * overlap

    # Normalized signal sums to approximately 1 across the genome.
    if is_norm:
        if fragment_count <= 0:
            raise ValueError(
                "Normalization requires non-zero total fragments.",
            )

        for k in sig_bin:
            sig_bin[k] /= fragment_count

    if scl_fct is not None:
        validate_comparison(scl_fct, "gt", 0, "scl_fct")

        for k in sig_bin:
            sig_bin[k] *= scl_fct

    return sig_bin


def calc_sig_chrom_array(
    chrom: str,
    starts: np.ndarray,
    ends: np.ndarray,
    lengths: np.ndarray,
    chrom_size: int,
    fragment_count: int,
    siz_bin: int,
    is_len: bool,
    is_norm: bool,
    scl_fct: float | None = None,
) -> dict[tuple[str, int], float]:
    """
    Compute binned signal for one chromosome with vectorized range additions.

    Parameters
    ----------
    chrom : str
        Chromosome represented by the fragment arrays.
    starts, ends, lengths : np.ndarray
        Fragment starts, ends, and lengths in matching order.
    chrom_size : int
        Chromosome length used to clip bins.
    fragment_count : int
        Total accepted fragments used for optional normalization.
    siz_bin : int
        Signal-bin width.
    is_len : bool
        Whether fragments contribute inverse-length weights.
    is_norm : bool
        Whether to apply depth normalization.
    scl_fct : float | None
        Optional multiplicative scale factor.

    Returns
    -------
    signal : dict[tuple[str, int], float]
        Sparse nonzero signal keyed by chromosome and bin start.

    Raises
    ------
    ValueError
        If chromosome size, fragment lengths, or normalization inputs are
        invalid.
    """

    validate_comparison(siz_bin, "gt", 0, "siz_bin", allow_none=False)

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
            raise ValueError("'frg_len' must be > 0 when using normalization.")

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
            (ends[same] - starts[same]) * weights[same],
        )
        touched[start_bins[same]] = True

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
        touched[start_bins_subset] = True
        touched[end_bins_subset] = True

        has_interior = end_bins_subset > start_bins_subset + 1

        if np.any(has_interior):
            diff = np.zeros(n_bins + 1, dtype=np.float64)
            touched_diff = np.zeros(n_bins + 1, dtype=np.int64)
            first = start_bins_subset[has_interior] + 1
            last_exclusive = end_bins_subset[has_interior]
            value = siz_bin * weights_subset[has_interior]

            np.add.at(diff, first, value)
            np.add.at(diff, last_exclusive, -value)

            np.add.at(touched_diff, first, 1)
            np.add.at(touched_diff, last_exclusive, -1)

            sig += np.cumsum(diff[:-1])
            touched |= np.cumsum(touched_diff[:-1]) > 0

    if is_norm:
        if fragment_count <= 0:
            raise ValueError(
                "Normalization requires non-zero total fragments.",
            )

        sig /= fragment_count

    if scl_fct is not None:
        validate_comparison(scl_fct, "gt", 0, "scl_fct")

        sig *= scl_fct

    nonzero_indices = np.flatnonzero(touched & (sig != 0.0))

    return {
        (chrom, int(index) * siz_bin): float(sig[index])
        for index in nonzero_indices
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
    emit_touched_only: bool = False,
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
    return_bin_indices : bool
        Whether sparse results use bin indices instead of base coordinates.
    local_span : bool
        Whether allocation is bounded to the observed local span.
    use_bincount : bool
        Whether aggregation uses NumPy bincount.
    emit_touched_only : bool
        Whether sparse output omits untouched or zero-valued bins.

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

    genome_n_bins = math.ceil(chrom_size / siz_bin)
    if starts.size == 0:
        return (
            "direct_sparse_idx_np"
            if return_bin_indices
            else "direct_sparse_np",
            [],
        )

    valid = ends > starts

    if not np.all(valid):
        starts = starts[valid]
        ends = ends[valid]
        lengths = lengths[valid]

    if starts.size == 0:
        return (
            "direct_sparse_idx_np"
            if return_bin_indices
            else "direct_sparse_np",
            [],
        )

    if is_len:
        if np.any(lengths <= 0):
            raise ValueError("'frg_len' must be > 0 when using normalization.")

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
            "direct_sparse_idx_np"
            if return_bin_indices
            else "direct_sparse_np",
            [],
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
                minlength=n_bins,
            )[:n_bins]
        else:
            np.add.at(sig, start_bins[same], same_values)

        if touched is not None:
            touched[start_bins[same]] = True

    multi = ~same

    if np.any(multi):
        starts_subset = starts[multi]
        ends_subset = ends[multi]
        start_bins_subset = start_bins[multi]
        end_bins_subset = end_bins[multi]
        global_start_bins_subset = start_bins_global[multi]
        global_end_bins_subset = end_bins_global[multi]
        weights_subset = weights[multi]

        left_overlap = (
            (global_start_bins_subset + 1) * siz_bin - starts_subset
        ) * weights_subset
        right_overlap = (
            ends_subset - global_end_bins_subset * siz_bin
        ) * weights_subset

        if use_bincount:
            sig += np.bincount(
                start_bins_subset,
                weights=left_overlap,
                minlength=n_bins,
            )[:n_bins]
            sig += np.bincount(
                end_bins_subset,
                weights=right_overlap,
                minlength=n_bins,
            )[:n_bins]
        else:
            np.add.at(sig, start_bins_subset, left_overlap)
            np.add.at(sig, end_bins_subset, right_overlap)

        if touched is not None:
            touched[start_bins_subset] = True
            touched[end_bins_subset] = True

        has_interior = end_bins_subset > start_bins_subset + 1

        if np.any(has_interior):
            first = start_bins_subset[has_interior] + 1
            last_exclusive = end_bins_subset[has_interior]
            value = siz_bin * weights_subset[has_interior]
            touched_diff = None

            if use_bincount:
                diff = np.bincount(
                    first,
                    weights=value,
                    minlength=n_bins + 1,
                ).astype(np.float64, copy=False)
                diff -= np.bincount(
                    last_exclusive,
                    weights=value,
                    minlength=n_bins + 1,
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
            "direct_sparse_idx_np"
            if return_bin_indices
            else "direct_sparse_np",
            [],
        )

    if return_bin_indices:
        index_dtype = (
            np.int32 if n_bins <= np.iinfo(np.int32).max else np.int64
        )

        return (
            "direct_sparse_idx_np",
            [
                (
                    chrom,
                    (idx + first_bin).astype(index_dtype, copy=False),
                    sig[idx].astype(np.float64, copy=False),
                ),
            ],
        )

    return (
        "direct_sparse_np",
        [
            (
                chrom,
                ((idx + first_bin) * siz_bin).astype(np.int64, copy=False),
                sig[idx].astype(np.float64, copy=False),
            ),
        ],
    )


def calc_sig_chrom_direct_dense_np(
    chrom: str,
    starts: np.ndarray,
    ends: np.ndarray,
    lengths: np.ndarray,
    chrom_size: int,
    siz_bin: int,
    is_len: bool,
) -> object:
    """
    Compute binned signal and return one dense chromosome array.

    Parameters
    ----------
    chrom : str
        Chromosome represented by the fragment arrays.
    starts, ends, lengths : np.ndarray
        Fragment starts, ends, and lengths in matching order.
    chrom_size : int
        Chromosome length used to clip bins.
    siz_bin : int
        Signal-bin width.
    is_len : bool
        Whether fragments contribute inverse-length weights.

    Returns
    -------
    result : object
        Tagged chromosome result containing dense signal values.

    Raises
    ------
    ValueError
        If chromosome size or fragment lengths are invalid.
    """

    validate_comparison(siz_bin, "gt", 0, "siz_bin", allow_none=False)

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
            raise ValueError("'frg_len' must be > 0 when using normalization.")

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
            (ends[same] - starts[same]) * weights[same],
        )
        touched[start_bins[same]] = True

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
        np.add.at(values, start_bins_subset, left_overlap)
        np.add.at(values, end_bins_subset, right_overlap)
        touched[start_bins_subset] = True
        touched[end_bins_subset] = True

        has_interior = end_bins_subset > start_bins_subset + 1

        if np.any(has_interior):
            diff = np.zeros(n_bins + 1, dtype=np.float64)
            touched_diff = np.zeros(n_bins + 1, dtype=np.int64)
            first = start_bins_subset[has_interior] + 1
            last_exclusive = end_bins_subset[has_interior]
            value = siz_bin * weights_subset[has_interior]

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
    is_len: bool,
) -> object:
    """
    Return bin edge and range-add events for parent-side materialization.

    Parameters
    ----------
    chrom : str
        Chromosome represented by the fragment arrays.
    starts, ends, lengths : np.ndarray
        Fragment starts, ends, and lengths in matching order.
    chrom_size : int
        Chromosome length used to clip events.
    siz_bin : int
        Signal-bin width.
    is_len : bool
        Whether fragments contribute inverse-length weights.

    Returns
    -------
    result : object
        Tagged chromosome result containing range-add event arrays.

    Raises
    ------
    ValueError
        If chromosome size or fragment lengths are invalid.
    """

    validate_comparison(siz_bin, "gt", 0, "siz_bin", allow_none=False)

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
            raise ValueError("'frg_len' must be > 0 when using normalization.")

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
                copy=False,
            ),
        )

    multi = ~same

    if np.any(multi):
        starts_subset = starts[multi]
        ends_subset = ends[multi]
        start_bins_subset = start_bins[multi]
        end_bins_subset = end_bins[multi]
        weights_subset = weights[multi]

        edge_bin_parts.extend(
            [
                start_bins_subset.astype(np.int64, copy=False),
                end_bins_subset.astype(np.int64, copy=False),
            ],
        )
        edge_value_parts.extend(
            [
                (
                    ((start_bins_subset + 1) * siz_bin - starts_subset)
                    * weights_subset
                ).astype(np.float64, copy=False),
                (
                    (ends_subset - end_bins_subset * siz_bin) * weights_subset
                ).astype(np.float64, copy=False),
            ],
        )
        has_interior = end_bins_subset > start_bins_subset + 1

        if np.any(has_interior):
            first = (start_bins_subset[has_interior] + 1).astype(
                np.int64,
                copy=False,
            )
            last_exclusive = end_bins_subset[has_interior].astype(
                np.int64,
                copy=False,
            )
            value = (siz_bin * weights_subset[has_interior]).astype(
                np.float64,
                copy=False,
            )
            diff_bin_parts.extend([first, last_exclusive])
            diff_value_parts.extend([value, -value])
            touch_bin_parts.extend([first, last_exclusive])
            touch_value_parts.extend(
                [
                    np.ones(first.shape, dtype=np.int64),
                    -np.ones(last_exclusive.shape, dtype=np.int64),
                ],
            )

    edge_bins = (
        np.concatenate(edge_bin_parts)
        if edge_bin_parts
        else np.empty(0, dtype=np.int64)
    )
    edge_values = (
        np.concatenate(edge_value_parts)
        if edge_value_parts
        else np.empty(0, dtype=np.float64)
    )
    diff_bins = (
        np.concatenate(diff_bin_parts)
        if diff_bin_parts
        else np.empty(0, dtype=np.int64)
    )
    diff_values = (
        np.concatenate(diff_value_parts)
        if diff_value_parts
        else np.empty(0, dtype=np.float64)
    )
    touch_bins = (
        np.concatenate(touch_bin_parts)
        if touch_bin_parts
        else np.empty(0, dtype=np.int64)
    )
    touch_values = (
        np.concatenate(touch_value_parts)
        if touch_value_parts
        else np.empty(0, dtype=np.int64)
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
                touch_values,
            ),
        ],
    )


def calc_sig_task(data: tuple[object, ...]) -> object:
    """
    Unpack one worker-task tuple and dispatch to 'calc_sig_chrom()'.
    """

    return calc_sig_chrom(*data)


def calc_sig_array_task(data: tuple[object, ...]) -> object:
    """
    Unpack one worker-task tuple and dispatch to 'calc_sig_chrom_array()'.
    """

    return calc_sig_chrom_array(*data)


def calc_sig_indexed_fetch_task(data: tuple[object, ...]) -> object:
    """
    Fetch one indexed region, parse fragments, and compute unscaled signal.

    Parameters
    ----------
    data : tuple[object, ...]
        Serialized indexed-fetch task arguments.

    Returns
    -------
    result : object
        Tagged task result for parent-side signal assembly.
    """

    (
        alignment_path,
        chrom,
        chrom_sizes,
        user_fragment_length,
        ref_fa,
        start,
        end,
        siz_bin,
        is_len,
        result_format,
    ) = data

    fragment_arrays, fragment_count = collect_fragment_arrays_from_iter(
        iter_indexed_fragments(
            alignment_path=alignment_path,
            chrom=chrom,
            chrom_sizes=chrom_sizes,
            user_fragment_length=user_fragment_length,
            ref_fa=ref_fa,
            start=start,
            end=end,
        ),
    )

    arrays = fragment_arrays.get(chrom)

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
        local_span = result_format in (
            "direct_sparse_local_np",
            "direct_sparse_local_bincount_np",
        )
        use_bincount = result_format in (
            "direct_sparse_bincount_np",
            "direct_sparse_touched_bincount_np",
            "direct_sparse_local_bincount_np",
        )
        emit_touched_only = result_format in (
            "direct_sparse_touched_np",
            "direct_sparse_touched_bincount_np",
        )
        sig = calc_sig_chrom_direct_sparse_np(
            chrom=chrom,
            starts=arrays[0],
            ends=arrays[1],
            lengths=arrays[2],
            chrom_size=chrom_sizes[chrom],
            siz_bin=siz_bin,
            is_len=is_len,
            return_bin_indices=result_format == "direct_sparse_idx_np",
            local_span=local_span,
            use_bincount=use_bincount,
            emit_touched_only=emit_touched_only,
        )
    elif result_format == "direct_dense_np":
        sig = calc_sig_chrom_direct_dense_np(
            chrom=chrom,
            starts=arrays[0],
            ends=arrays[1],
            lengths=arrays[2],
            chrom_size=chrom_sizes[chrom],
            siz_bin=siz_bin,
            is_len=is_len,
        )
    elif result_format == "event_np":
        sig = calc_sig_chrom_event_np(
            chrom=chrom,
            starts=arrays[0],
            ends=arrays[1],
            lengths=arrays[2],
            chrom_size=chrom_sizes[chrom],
            siz_bin=siz_bin,
            is_len=is_len,
        )
    else:
        sig = calc_sig_chrom_array(
            chrom=chrom,
            starts=arrays[0],
            ends=arrays[1],
            lengths=arrays[2],
            chrom_size=chrom_sizes[chrom],
            fragment_count=1,
            siz_bin=siz_bin,
            is_len=is_len,
            is_norm=False,
            scl_fct=None,
        )

    if result_format in (
        *DIRECT_SPARSE_RESULT_FORMATS,
        "direct_dense_np",
        "event_np",
    ):
        return sig, fragment_count

    return signal_dict_to_result(sig, result_format, siz_bin), fragment_count


def collect_bed_indexed_task(data: tuple[object, ...]) -> object:
    """
    Fetch one indexed region and return compact BED row arrays.
    """

    (
        alignment_path,
        chrom,
        chrom_sizes,
        user_fragment_length,
        ref_fa,
        start,
        end,
    ) = data

    return collect_bed_arrays_from_iter(
        iter_indexed_fragments(
            alignment_path=alignment_path,
            chrom=chrom,
            chrom_sizes=chrom_sizes,
            user_fragment_length=user_fragment_length,
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


def estimate_bed_payload_bytes(result: Any) -> int:
    """
    Estimate NumPy payload bytes in one compact BED array result.
    """

    if not isinstance(result, tuple) or len(result) < 2:
        return 0

    return sum(
        int(starts.nbytes + ends.nbytes + lengths.nbytes)
        for _chrom, starts, ends, lengths in result[1]
    )


def calc_sig_chunk_task(data: tuple[object, ...]) -> object:
    """
    Unpack one worker-task tuple and compute signal for a fragment chunk.
    """

    fragments, fragment_count, siz_bin, is_len, is_norm, scl_fct = data

    by_chrom = defaultdict(list)

    for chrom, start, end, fragment_length in fragments:
        by_chrom[chrom].append((start, end, fragment_length))

    sig_chunk = defaultdict(float)

    for chrom, entries in by_chrom.items():
        sig_chrom = calc_sig_chrom(
            chrom=chrom,
            fragment_records=entries,
            fragment_count=fragment_count,
            siz_bin=siz_bin,
            is_len=is_len,
            is_norm=is_norm,
            scl_fct=scl_fct,
        )

        for key, value in sig_chrom.items():
            sig_chunk[key] += value

    return sig_chunk


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

    return (
        task_id,
        {
            "worker_start_s": start,
            "worker_end_s": end,
            "worker_elapsed_s": end - start,
            "result_bins": result_count,
        },
        result,
    )


def iter_task_results(
    kind: str,
    task_func: Callable[[object], object],
    task_data: Iterable[object],
    threads: int,
    executor_mode: str,
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
    executor_mode : str
        Ordered-map or completion-order execution mode.
    task_profiles : list[dict] | None
        Optional collection receiving worker timing records.

    Yields
    ------
    result : object
        Task results in the order defined by 'executor_mode'.
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
                        "worker_elapsed_s": end - start,
                        "parent_received_s": received,
                        "parent_receive_lag_s": received - end,
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
                profiled_data,
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
    fragment_count: int,
    is_norm: bool,
    scl_fct: float | None = None,
) -> None:
    """
    Apply global depth normalization and scaling to already merged signal.

    Parameters
    ----------
    signal_bins : dict[tuple[str, int], float]
        Mutable sparse signal values.
    fragment_count : int
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

    if is_array_signal(signal_bins):
        if is_norm:
            if fragment_count <= 0:
                raise ValueError(
                    "Normalization requires non-zero total fragments.",
                )

            for values in iter_signal_value_arrays(signal_bins):
                values /= fragment_count

        if scl_fct is not None:
            validate_comparison(scl_fct, "gt", 0, "scl_fct")

            for values in iter_signal_value_arrays(signal_bins):
                values *= scl_fct

        return

    if is_norm:
        if fragment_count <= 0:
            raise ValueError(
                "Normalization requires non-zero total fragments.",
            )

        for key in signal_bins:
            signal_bins[key] /= fragment_count

    if scl_fct is not None:
        validate_comparison(scl_fct, "gt", 0, "scl_fct")

        for key in signal_bins:
            signal_bins[key] *= scl_fct


def is_array_sparse_signal(signal_bins: object) -> bool:
    """
    Return True for the merged sparse-array signal container.
    """

    return (
        isinstance(signal_bins, tuple)
        and len(signal_bins) == 2
        and signal_bins[0] == "array_sparse_merged"
    )


def is_array_dense_signal(signal_bins: object) -> bool:
    """
    Return True for the merged dense-array signal container.
    """

    return (
        isinstance(signal_bins, tuple)
        and len(signal_bins) == 3
        and signal_bins[0] == "array_dense_merged"
    )


def is_array_signal(signal_bins: object) -> bool:
    """
    Return True for merged signal containers that avoid the final dict.
    """

    return is_array_sparse_signal(signal_bins) or is_array_dense_signal(
        signal_bins,
    )


def iter_signal_value_arrays(signal_bins: Any) -> Iterator[np.ndarray]:
    """
    Yield mutable value arrays from merged array-backed signal containers.
    """

    if is_array_sparse_signal(signal_bins):
        for _, _, values in signal_bins[1]:
            yield values
    elif is_array_dense_signal(signal_bins):
        for _, values, _ in signal_bins[2]:
            yield values


def signal_row_count(signal_bins: Any) -> int:
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


def sparse_parts_are_ordered(starts_parts: list[np.ndarray]) -> bool:
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


def coalesce_sparse_signal_parts(
    parts_by_chrom: dict[str, list[tuple[np.ndarray, np.ndarray]]],
    optimize: bool = True,
    stats: dict[str, int] | None = None,
) -> object:
    """
    Sort and sum sparse signal arrays by chromosome without building a dict.

    Parameters
    ----------
    parts_by_chrom : dict[str, list[tuple[np.ndarray, np.ndarray]]]
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
            merged.append((chrom, out_starts, out_values))

            if stats is not None:
                stats["rows_after"] += int(out_starts.size)

    return "array_sparse_merged", merged


def sparse_payload_to_starts(
    tag: str,
    coords: np.ndarray,
    siz_bin: int,
) -> np.ndarray:
    """
    Convert sparse payload coordinates to bedGraph bin starts.
    """

    if tag == "direct_sparse_idx_np":
        return coords.astype(np.int64, copy=False) * siz_bin

    return coords.astype(np.int64, copy=False)


def materialize_event_signal_parts(
    parts_by_chrom: dict[str, list[tuple]],
    siz_bin: int,
) -> object:
    """
    Combine event arrays and materialize final sparse signal arrays.

    Parameters
    ----------
    parts_by_chrom : dict[str, list[tuple]]
        Range-add event-array parts grouped by chromosome.
    siz_bin : int
        Signal-bin width used to derive output coordinates.

    Returns
    -------
    signal : object
        Ordered sparse signal arrays grouped by chromosome.
    """

    merged = []

    for chrom in sorted(parts_by_chrom, key=sort_chrom):
        n_bins = max(part[0] for part in parts_by_chrom[chrom])
        values = np.zeros(n_bins, dtype=np.float64)
        touched = np.zeros(n_bins, dtype=bool)

        edge_bins_parts = [part[1] for part in parts_by_chrom[chrom]]
        edge_values_parts = [part[2] for part in parts_by_chrom[chrom]]
        edge_bins = (
            np.concatenate(edge_bins_parts)
            if edge_bins_parts
            else (np.empty(0, dtype=np.int64))
        )
        edge_values = (
            np.concatenate(edge_values_parts)
            if edge_values_parts
            else np.empty(0, dtype=np.float64)
        )

        if edge_bins.size > 0:
            np.add.at(values, edge_bins, edge_values)
            touched[edge_bins] = True

        diff_bins_parts = [part[3] for part in parts_by_chrom[chrom]]
        diff_values_parts = [part[4] for part in parts_by_chrom[chrom]]
        diff_bins = (
            np.concatenate(diff_bins_parts)
            if diff_bins_parts
            else (np.empty(0, dtype=np.int64))
        )
        diff_values = (
            np.concatenate(diff_values_parts)
            if diff_values_parts
            else np.empty(0, dtype=np.float64)
        )

        if diff_bins.size > 0:
            diff = np.zeros(n_bins + 1, dtype=np.float64)
            np.add.at(diff, diff_bins, diff_values)
            values += np.cumsum(diff[:-1])

        touch_bins_parts = [part[5] for part in parts_by_chrom[chrom]]
        touch_values_parts = [part[6] for part in parts_by_chrom[chrom]]
        touch_bins = (
            np.concatenate(touch_bins_parts)
            if touch_bins_parts
            else (np.empty(0, dtype=np.int64))
        )
        touch_values = (
            np.concatenate(touch_values_parts)
            if touch_values_parts
            else np.empty(0, dtype=np.int64)
        )

        if touch_bins.size > 0:
            touch_diff = np.zeros(n_bins + 1, dtype=np.int64)
            np.add.at(touch_diff, touch_bins, touch_values)
            touched |= np.cumsum(touch_diff[:-1]) > 0

        idx = np.flatnonzero(touched & (values != 0.0))

        if idx.size > 0:
            merged.append(
                (
                    chrom,
                    (idx * siz_bin).astype(np.int64, copy=False),
                    values[idx].astype(np.float64, copy=False),
                ),
            )

    return "array_sparse_merged", merged


def _append_mapping_sparse_parts(
    signal_result: object,
    sparse_parts: dict[str, list[tuple[np.ndarray, np.ndarray]]],
) -> None:
    """
    Convert one mapping result into chromosome-keyed sparse arrays.

    Parameters
    ----------
    signal_result : object
        Mapping from chromosome/start keys to signal values.
    sparse_parts : dict[str, list[tuple[np.ndarray, np.ndarray]]]
        Mutable chromosome inventory receiving start and value arrays.
    """

    by_chrom: dict[str, list[tuple[int, float]]] = defaultdict(list)

    for (chrom, start), value in signal_result.items():
        by_chrom[chrom].append((start, value))

    for chrom, rows in by_chrom.items():
        sparse_parts[chrom].append(
            (
                np.asarray(
                    [start for start, _ in rows],
                    dtype=np.int64,
                ),
                np.asarray(
                    [value for _, value in rows],
                    dtype=np.float64,
                ),
            ),
        )


def _accumulate_sparse_merge_result(
    signal_result: object,
    sparse_parts: dict[str, list[tuple[np.ndarray, np.ndarray]]],
    siz_bin: int,
) -> None:
    """
    Append one worker result to sparse-merge chromosome parts.

    Mapping, sparse-array, and dense-array worker representations are
    normalized to chromosome-keyed start and value arrays.

    Parameters
    ----------
    signal_result : object
        Mapping- or tagged-tuple worker result.
    sparse_parts : dict[str, list[tuple[np.ndarray, np.ndarray]]]
        Mutable chromosome inventory receiving sparse arrays.
    siz_bin : int
        Signal-bin width used to convert bin indexes to starts.
    """

    if not isinstance(signal_result, tuple) or not signal_result:
        _append_mapping_sparse_parts(signal_result, sparse_parts)

        return

    tag = signal_result[0]

    if tag in (
        "sparse_np",
        "direct_sparse_np",
        "direct_sparse_idx_np",
    ):
        for chrom, coords, values in signal_result[1]:
            sparse_parts[chrom].append(
                (
                    sparse_payload_to_starts(tag, coords, siz_bin),
                    values,
                ),
            )

        return

    if tag == "dense_np":
        for chrom, first, values in signal_result[1]:
            idx = np.flatnonzero(values != 0.0)
            sparse_parts[chrom].append(
                (
                    first + idx.astype(np.int64) * siz_bin,
                    values[idx],
                ),
            )

        return

    _append_mapping_sparse_parts(signal_result, sparse_parts)


def _accumulate_dense_merge_result(
    signal_result: object,
    ensure_array: Callable[[str], np.ndarray],
    ensure_touched: Callable[[str], np.ndarray],
    siz_bin: int,
) -> None:
    """
    Accumulate one worker result into dense arrays and touch masks.

    Parameters
    ----------
    signal_result : object
        Mapping- or tagged-tuple worker result.
    ensure_array : Callable[[str], np.ndarray]
        Resolver for a chromosome's mutable dense signal array.
    ensure_touched : Callable[[str], np.ndarray]
        Resolver for a chromosome's mutable boolean touch mask.
    siz_bin : int
        Signal-bin width used to convert starts to array indexes.
    """

    if not isinstance(signal_result, tuple) or not signal_result:
        for (chrom, start), value in signal_result.items():
            ensure_array(chrom)[start // siz_bin] += value
            ensure_touched(chrom)[start // siz_bin] = True

        return

    tag = signal_result[0]

    if tag == "direct_dense_np":
        for chrom, values, touched in signal_result[1]:
            arr = ensure_array(chrom)
            arr[: values.size] += values
            ensure_touched(chrom)[: touched.size] |= touched

        return

    if tag == "dense_np":
        for chrom, first, values in signal_result[1]:
            arr = ensure_array(chrom)
            first_idx = first // siz_bin
            arr[first_idx : first_idx + values.size] += values
            ensure_touched(chrom)[first_idx : first_idx + values.size] |= (
                values != 0.0
            )

        return

    if tag in (
        "sparse_np",
        "direct_sparse_np",
        "direct_sparse_idx_np",
    ):
        for chrom, coords, values in signal_result[1]:
            if tag == "direct_sparse_idx_np":
                idx = coords.astype(np.int64, copy=False)
            else:
                idx = coords // siz_bin

            np.add.at(ensure_array(chrom), idx, values)
            ensure_touched(chrom)[idx] = True

        return

    for (chrom, start), value in signal_result.items():
        ensure_array(chrom)[start // siz_bin] += value
        ensure_touched(chrom)[start // siz_bin] = True


def _accumulate_event_merge_result(
    signal_result: object,
    event_parts: dict[str, list[tuple]],
    sparse_parts: dict[str, list[tuple[np.ndarray, np.ndarray]]],
    siz_bin: int,
) -> None:
    """
    Append one worker result to event or sparse fallback parts.

    Parameters
    ----------
    signal_result : object
        Mapping- or tagged-tuple worker result.
    event_parts : dict[str, list[tuple]]
        Mutable chromosome inventory receiving event-array payloads.
    sparse_parts : dict[str, list[tuple[np.ndarray, np.ndarray]]]
        Mutable chromosome inventory receiving fallback sparse arrays.
    siz_bin : int
        Signal-bin width used to convert bin indexes to starts.
    """

    if not isinstance(signal_result, tuple) or not signal_result:
        for (chrom, start), value in signal_result.items():
            sparse_parts[chrom].append(
                (
                    np.asarray([start], dtype=np.int64),
                    np.asarray([value], dtype=np.float64),
                ),
            )

        return

    tag = signal_result[0]

    if tag == "event_np":
        for (
            chrom,
            n_bins,
            edge_bins,
            edge_values,
            diff_bins,
            diff_values,
            touch_bins,
            touch_values,
        ) in signal_result[1]:
            event_parts[chrom].append(
                (
                    n_bins,
                    edge_bins,
                    edge_values,
                    diff_bins,
                    diff_values,
                    touch_bins,
                    touch_values,
                ),
            )

        return

    if tag in (
        "sparse_np",
        "direct_sparse_np",
        "direct_sparse_idx_np",
    ):
        for chrom, coords, values in signal_result[1]:
            sparse_parts[chrom].append(
                (
                    sparse_payload_to_starts(tag, coords, siz_bin),
                    values,
                ),
            )

        return

    for (chrom, start), value in signal_result.items():
        sparse_parts[chrom].append(
            (
                np.asarray([start], dtype=np.int64),
                np.asarray([value], dtype=np.float64),
            ),
        )


def _accumulate_dict_merge_result(
    signal_result: object,
    combined: defaultdict[tuple[str, int], float],
    siz_bin: int,
) -> None:
    """
    Convert and accumulate one worker result into a signal dictionary.

    Parameters
    ----------
    signal_result : object
        Mapping- or tagged-tuple worker result.
    combined : defaultdict[tuple[str, int], float]
        Mutable chromosome/start signal mapping.
    siz_bin : int
        Signal-bin width used to convert array indexes to starts.
    """

    if not isinstance(signal_result, tuple) or not signal_result:
        for key, value in signal_result.items():
            combined[key] += value

        return

    tag = signal_result[0]

    if tag in (
        "sparse_np",
        "direct_sparse_np",
        "direct_sparse_idx_np",
    ):
        for chrom, coords, values in signal_result[1]:
            starts = sparse_payload_to_starts(tag, coords, siz_bin)

            for start, value in zip(starts, values, strict=True):
                combined[(chrom, int(start))] += float(value)

        return

    if tag == "direct_dense_np":
        for chrom, values, touched in signal_result[1]:
            idx = np.flatnonzero(touched & (values != 0.0))

            for index in idx:
                combined[(chrom, int(index) * siz_bin)] += float(values[index])

        return

    if tag == "dense_np":
        for chrom, first, values in signal_result[1]:
            for offset, value in enumerate(values):
                if value != 0.0:
                    combined[(chrom, first + offset * siz_bin)] += float(value)

        return

    for key, value in signal_result.items():
        combined[key] += value


def _accumulate_array_merge_result(
    signal_result: object,
    ensure_array: Callable[[str], np.ndarray],
    siz_bin: int,
) -> None:
    """
    Accumulate one worker result into chromosome arrays.

    Parameters
    ----------
    signal_result : object
        Mapping- or tagged-tuple worker result.
    ensure_array : Callable[[str], np.ndarray]
        Resolver for a chromosome's mutable dense signal array.
    siz_bin : int
        Signal-bin width used to convert starts to array indexes.
    """

    if not isinstance(signal_result, tuple) or not signal_result:
        for (chrom, start), value in signal_result.items():
            ensure_array(chrom)[start // siz_bin] += value

        return

    tag = signal_result[0]

    if tag in (
        "sparse_np",
        "direct_sparse_np",
        "direct_sparse_idx_np",
    ):
        for chrom, coords, values in signal_result[1]:
            arr = ensure_array(chrom)

            if tag == "direct_sparse_idx_np":
                idx = coords.astype(np.int64, copy=False)
            else:
                idx = coords // siz_bin

            np.add.at(arr, idx, values)

        return

    if tag == "dense_np":
        for chrom, first, values in signal_result[1]:
            arr = ensure_array(chrom)
            first_idx = first // siz_bin
            arr[first_idx : first_idx + values.size] += values

        return

    for (chrom, start), value in signal_result.items():
        ensure_array(chrom)[start // siz_bin] += value


def merge_signal_results(
    results: Iterable[Any],
    merge_strategy: str,
    chrom_sizes: dict[str, int],
    siz_bin: int,
    is_prototype_strategy: bool,
    profile: dict[str, Any] | None = None,
) -> tuple[Any, int, int, int]:
    """
    Merge worker signal results and return signal, fragments, bins, and bytes.

    Parameters
    ----------
    results : Iterable[Any]
        Per-chromosome signal results to validate and merge.
    merge_strategy : str
        Named implementation used to combine worker results.
    chrom_sizes : dict[str, int]
        Chromosome lengths keyed by chromosome identifier.
    siz_bin : int
        Positive signal-bin width in base pairs.
    is_prototype_strategy : bool
        Whether the selected merge strategy is explicitly experimental.
    profile : dict[str, Any] | None
        Optional mutable timing and merge-statistics payload.

    Returns
    -------
    combined_signal, fragment_count, bin_count, payload_bytes : tuple[
        Any, int, int, int
    ]
        Mapping- or array-backed merged signal and execution statistics.

    Notes
    -----
    The first return value intentionally remains 'Any' because benchmark
    strategies use several tagged tuple representations in addition to the
    production signal mapping.
    """

    combined_signal: Any = defaultdict(float)
    combined_arrays: dict[str, np.ndarray] = {}
    array_touched: dict[str, np.ndarray] = {}
    sparse_parts: dict[str, list[tuple[np.ndarray, np.ndarray]]] = defaultdict(
        list,
    )
    event_parts: dict[str, list[tuple]] = defaultdict(list)

    fragment_count = 0
    bin_count = 0
    payload_bytes = 0

    metadata_seconds = 0.0
    accumulation_seconds = 0.0
    array_to_dict_seconds = 0.0
    array_coalesce_seconds = 0.0
    dense_finalize_seconds = 0.0
    event_materialize_seconds = 0.0

    sparse_coalesce_stats: dict[str, int] = {}

    def ensure_array(chrom: str) -> np.ndarray:
        """
        Ensure a dense chromosome array exists for one chromosomal key.

        Parameters
        ----------
        chrom : str
            Chromosome name used as the lookup key.

        Returns
        -------
        array : np.ndarray
            Mutable float64 array sized to the chromosome span under the
            current bin width.
        """

        if chrom not in combined_arrays:
            combined_arrays[chrom] = np.zeros(
                math.ceil(chrom_sizes[chrom] / siz_bin),
                dtype=np.float64,
            )

        return combined_arrays[chrom]

    def ensure_touched(chrom: str) -> np.ndarray:
        """
        Ensure a boolean touch mask exists for one chromosomal key.

        Parameters
        ----------
        chrom : str
            Chromosome name used as the lookup key.

        Returns
        -------
        mask : np.ndarray
            Mutable boolean mask sized to the chromosome span under the current
            bin width.
        """

        if chrom not in array_touched:
            array_touched[chrom] = np.zeros(
                math.ceil(chrom_sizes[chrom] / siz_bin),
                dtype=bool,
            )

        return array_touched[chrom]

    for result in results:
        time_meta = time.perf_counter()

        if is_prototype_strategy:
            signal_result, result_fragment_count = result
            fragment_count += result_fragment_count
        else:
            signal_result = result

        bin_count += count_result_bins(signal_result)
        payload_bytes += estimate_result_payload_bytes(signal_result)

        if profile is not None:
            metadata_seconds += time.perf_counter() - time_meta

        time_accumulate = time.perf_counter()

        if merge_strategy in (
            "array_sparse_merge",
            "array_sparse_merge_legacy",
        ):
            _accumulate_sparse_merge_result(
                signal_result,
                sparse_parts,
                siz_bin,
            )
        elif merge_strategy == "array_dense_merge":
            _accumulate_dense_merge_result(
                signal_result,
                ensure_array,
                ensure_touched,
                siz_bin,
            )
        elif merge_strategy == "event_diff_merge":
            _accumulate_event_merge_result(
                signal_result,
                event_parts,
                sparse_parts,
                siz_bin,
            )
        elif merge_strategy == "dict_merge":
            _accumulate_dict_merge_result(
                signal_result,
                combined_signal,
                siz_bin,
            )
        else:
            _accumulate_array_merge_result(
                signal_result,
                ensure_array,
                siz_bin,
            )

        if profile is not None:
            accumulation_seconds += time.perf_counter() - time_accumulate

    dict_converting_merges = ("chrom_array_merge", "vectorized_merge")

    if merge_strategy in dict_converting_merges:
        time_convert = time.perf_counter()

        for chrom, array in combined_arrays.items():
            nonzero_indexes = np.flatnonzero(array != 0.0)

            for index in nonzero_indexes:
                combined_signal[(chrom, int(index) * siz_bin)] = float(
                    array[index],
                )

        if profile is not None:
            array_to_dict_seconds = time.perf_counter() - time_convert

    if merge_strategy in ("array_sparse_merge", "array_sparse_merge_legacy"):
        time_coalesce = time.perf_counter()
        combined_signal = coalesce_sparse_signal_parts(
            sparse_parts,
            optimize=merge_strategy == "array_sparse_merge",
            stats=sparse_coalesce_stats,
        )

        if profile is not None:
            array_coalesce_seconds = time.perf_counter() - time_coalesce

    if merge_strategy == "array_dense_merge":
        time_dense = time.perf_counter()
        combined_signal = (
            "array_dense_merged",
            siz_bin,
            [
                (
                    chrom,
                    combined_arrays[chrom],
                    array_touched[chrom],
                )
                for chrom in sorted(combined_arrays, key=sort_chrom)
                if np.any(
                    array_touched[chrom] & (combined_arrays[chrom] != 0.0),
                )
            ],
        )

        if profile is not None:
            dense_finalize_seconds = time.perf_counter() - time_dense

    if merge_strategy == "event_diff_merge":
        time_event = time.perf_counter()
        combined_signal = materialize_event_signal_parts(
            event_parts,
            siz_bin,
        )

        if sparse_parts:
            _, sparse_merged = coalesce_sparse_signal_parts(sparse_parts)
            combined_parts: dict[str, list[tuple[np.ndarray, np.ndarray]]] = (
                defaultdict(list)
            )

            for chrom, starts, values in combined_signal[1] + sparse_merged:
                combined_parts[chrom].append((starts, values))

            combined_signal = coalesce_sparse_signal_parts(combined_parts)

        if profile is not None:
            event_materialize_seconds = time.perf_counter() - time_event

    if profile is not None:
        phases = profile["phases_s"]
        phases["merge_result_metadata"] = metadata_seconds
        phases["merge_accumulate_results"] = accumulation_seconds
        phases["merge_array_to_dict"] = array_to_dict_seconds
        phases["merge_array_coalesce"] = array_coalesce_seconds
        phases["merge_array_dense_finalize"] = dense_finalize_seconds
        phases["merge_event_materialize"] = event_materialize_seconds
        profile["sparse_coalesce_stats"] = sparse_coalesce_stats

    return (
        combined_signal,
        fragment_count,
        bin_count,
        payload_bytes,
    )


def format_bdg_rows(
    rows: list[tuple[tuple[str, int], float]],
    siz_bin: int,
    decimal_places: int,
    chrom_sizes: dict[str, int] | None,
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
                    f"Missing chromosome size for bedGraph row: {chrom!r}.",
                )

            if bin_start < 0 or bin_start >= chrom_size:
                raise ValueError(
                    "bedGraph bin start is outside chromosome bounds: "
                    f"{chrom}:{bin_start} (chromosome size {chrom_size}).",
                )

            bin_end = min(bin_end, chrom_size)

        lines.append(
            f"{chrom}\t{bin_start}\t{bin_end}\t"
            f"{format_bdg_value(value, decimal_places)}\n",
        )

    return "".join(lines)


def format_bdg_value(value: float, decimal_places: int) -> str:
    """
    Format one bedGraph value with the same rounding used by write_bdg().
    """

    rounded_value = round(float(value), decimal_places)

    if rounded_value == 0.0:
        rounded_value = 0.0

    rendered_value = f"{rounded_value:.{decimal_places}f}"

    if "." in rendered_value:
        rendered_value = rendered_value.rstrip("0").rstrip(".")

    if rendered_value == "-0":
        rendered_value = "0"

    return rendered_value


def iter_sorted_bed_array_rows(
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
            yield (chrom, int(starts[idx]), int(ends[idx]), int(lengths[idx]))


def write_bed_array_results(
    results: Iterable[object],
    fil_out: str,
) -> int:
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
    fragment_records: dict[str, list[tuple[int, int, int]]],
    fil_out: str,
) -> int:
    """
    Write chromosome-grouped fragment tuples in deterministic BED order.
    """

    row_count = 0

    with open_out(fil_out) as bed_file:
        for chrom in sorted(fragment_records.keys(), key=sort_chrom):
            for start, end, length in sorted(
                fragment_records[chrom],
                key=lambda t: t[0],
            ):
                bed_file.write(f"{chrom}\t{start}\t{end}\t{length}\n")
                row_count += 1

    return row_count


def format_bdg_array_rows(
    chrom: str,
    starts: np.ndarray,
    values: np.ndarray,
    siz_bin: int,
    decimal_places: int,
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
            f"{format_bdg_value(float(value), decimal_places)}\n",
        )

    return "".join(lines)


def format_bdg_rows_task(data: tuple[object, ...]) -> str:
    """
    ProcessPool-friendly wrapper for formatting one bedGraph shard.
    """

    return format_bdg_rows(*data)


def format_bdg_array_rows_task(data: tuple[object, ...]) -> str:
    """
    ProcessPool-friendly wrapper for sparse-array bedGraph formatting.
    """

    return format_bdg_array_rows(*data)


def iter_array_signal_chunks(
    signal_bins: Any,
    target_chunks: int,
) -> Iterator[tuple[str, np.ndarray, np.ndarray]]:
    """
    Yield ordered array-backed signal chunks for parallel formatting.

    Parameters
    ----------
    signal_bins : Any
        Supported array-backed signal representation.
    target_chunks : int
        Approximate number of ordered output chunks.

    Yields
    ------
    chunk : tuple[str, np.ndarray, np.ndarray]
        Chromosome with matching bin-index and value arrays.
    """

    total_rows = max(signal_row_count(signal_bins), 1)
    chunk_size = max(1, math.ceil(total_rows / max(target_chunks, 1)))

    if is_array_sparse_signal(signal_bins):
        for chrom, starts, values in signal_bins[1]:
            for i in range(0, starts.size, chunk_size):
                yield (
                    chrom,
                    starts[i : i + chunk_size],
                    values[i : i + chunk_size],
                )
    elif is_array_dense_signal(signal_bins):
        dense_siz_bin = signal_bins[1]

        for chrom, values, touched in signal_bins[2]:
            idx = np.flatnonzero(touched & (values != 0.0))
            starts = (idx * dense_siz_bin).astype(np.int64, copy=False)

            for i in range(0, idx.size, chunk_size):
                chunk_starts = starts[i : i + chunk_size]
                chunk_values = values[idx[i : i + chunk_size]]

                yield (chrom, chunk_starts, chunk_values)


def digest_bdg_rows(
    coverage: Any,
    siz_bin: int,
    decimal_places: int,
    chrom_sizes: dict[str, int] | None,
) -> tuple[int, str]:
    """
    Return row count and SHA-256 for deterministic bedGraph-formatted rows.

    Parameters
    ----------
    coverage : Any
        Supported signal representation.
    siz_bin : int
        Signal-bin width.
    decimal_places : int
        Decimal precision for rendered values.
    chrom_sizes : dict[str, int] | None
        Optional chromosome sizes used to clip final intervals.

    Returns
    -------
    row_count, digest : tuple[int, str]
        Rendered row count and hexadecimal SHA-256 digest.
    """

    if is_array_signal(coverage):
        digest = hashlib.sha256()
        row_count = 0

        for chrom, starts, values in iter_array_signal_chunks(coverage, 1):
            text = format_bdg_array_rows(
                chrom,
                starts,
                values,
                siz_bin,
                decimal_places,
                chrom_sizes.get(chrom) if chrom_sizes is not None else None,
            )
            row_count += text.count("\n")
            digest.update(text.encode("utf-8"))

        return row_count, digest.hexdigest()

    digest = hashlib.sha256()
    row_count = 0
    items = sorted(
        coverage.items(),
        key=lambda item: key_bin(item[0][0], item[0][1]),
    )
    chunk_size = 100000

    for index in range(0, len(items), chunk_size):
        text = format_bdg_rows(
            items[index : index + chunk_size],
            siz_bin,
            decimal_places,
            chrom_sizes,
        )
        row_count += text.count("\n")
        digest.update(text.encode("utf-8"))

    return row_count, digest.hexdigest()


def write_bdg_array_sparse(
    signal_bins: Any,
    fil_out: str,
    siz_bin: int,
    decimal_places: int,
    chrom_sizes: dict[str, int] | None,
) -> None:
    """
    Write a sparse-array signal directly to bedGraph.
    """

    with open_out(fil_out) as handle:
        for chrom, starts, values in iter_array_signal_chunks(signal_bins, 1):
            handle.write(
                format_bdg_array_rows(
                    chrom,
                    starts,
                    values,
                    siz_bin,
                    decimal_places,
                    chrom_sizes.get(chrom)
                    if chrom_sizes is not None
                    else None,
                ),
            )


def write_bdg_parallel_ordered(
    coverage: Any,
    fil_out: str,
    siz_bin: int,
    decimal_places: int,
    chrom_sizes: dict[str, int] | None,
    workers: int,
) -> None:
    """
    Format sorted bedGraph shards in worker processes and write in order.

    Parameters
    ----------
    coverage : Any
        Supported signal representation.
    fil_out : str
        Output bedGraph path.
    siz_bin : int
        Signal-bin width.
    decimal_places : int
        Decimal precision for rendered values.
    chrom_sizes : dict[str, int] | None
        Optional chromosome sizes used to clip final intervals.
    workers : int
        Number of formatting worker processes.
    """

    validate_comparison(
        workers,
        "ge",
        1,
        "prototype_writer_workers",
        allow_none=False,
    )

    if is_array_signal(coverage):
        if workers == 1:
            write_bdg_array_sparse(
                coverage,
                fil_out,
                siz_bin,
                decimal_places,
                chrom_sizes,
            )

            return

        task_data = [
            (
                chrom,
                starts,
                values,
                siz_bin,
                decimal_places,
                chrom_sizes.get(chrom) if chrom_sizes is not None else None,
            )
            for chrom, starts, values in iter_array_signal_chunks(
                coverage,
                workers,
            )
        ]

        try:
            with ProcessPoolExecutor(
                max_workers=min(workers, os.cpu_count() or 1),
            ) as ex:
                formatted = list(ex.map(format_bdg_array_rows_task, task_data))
        except OSError:
            write_bdg_array_sparse(
                coverage,
                fil_out,
                siz_bin,
                decimal_places,
                chrom_sizes,
            )

            return

        with open_out(fil_out) as handle:
            for shard in formatted:
                handle.write(shard)

        return

    if workers == 1:
        write_bdg(
            coverage,
            fil_out,
            siz_bin,
            decimal_places,
            chrom_sizes=chrom_sizes,
        )

        return

    items = sorted(
        coverage.items(),
        key=lambda item: key_bin(item[0][0], item[0][1]),
    )

    if not items:
        with open_out(fil_out):
            return

    chunk_size = math.ceil(len(items) / workers)
    chunks = [
        items[index : index + chunk_size]
        for index in range(0, len(items), chunk_size)
    ]
    task_data = [
        (chunk, siz_bin, decimal_places, chrom_sizes) for chunk in chunks
    ]

    try:
        with ProcessPoolExecutor(
            max_workers=min(workers, os.cpu_count() or 1),
        ) as ex:
            formatted = list(ex.map(format_bdg_rows_task, task_data))
    except OSError:
        write_bdg(
            coverage,
            fil_out,
            siz_bin,
            decimal_places,
            chrom_sizes=chrom_sizes,
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
            "(https://github.com/kalavattam/siQ-ChIP/tree/protocol)."
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
            "parallel.\n\n"
        ),
    )

    parser.add_argument(
        "-fi",
        "--fil_in",
        dest="fil_in",
        required=True,
        help="Input file path for the BAM or CRAM file.\n\n",
    )
    parser.add_argument(
        "-rf",
        "--ref_fa",
        dest="ref_fa",
        default=None,
        help=(
            "Reference FASTA file for CRAM decoding. Required for CRAM inputs "
            "unless the reference is otherwise available to htslib/pysam.\n\n"
        ),
    )
    parser.add_argument(
        "-cs",
        "--chr_sizes",
        "--chrom_sizes",
        dest="chr_sizes",
        default=None,
        help=(
            "Chromosome sizes file in UCSC-style TSV format with chromosome "
            "name and positive integer size columns. Header sizes from "
            "BAM/CRAM are used when available; this file can supplement "
            "missing header sizes, but conflicting sizes are rejected.\n\n"
        ),
    )
    parser.add_argument(
        "-fo",
        "--fil_out",
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
        ),
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
            "  - Normalized-coverage aliases: 'n', 'nrm', 'norm', "
            "'normalized'. Internally standardized to 'norm'.\n"
            "\n"
            "Note: 'norm' normalizes for both fragment length and total "
            "fragment count so that the genome-wide summed signal is "
            "approximately 1.\n\n"
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
            "%(default)s).\n\n"
        ),
    )
    parser.add_argument(
        "-eg",
        "--engine",
        dest="engine",
        choices=ENGINE_CHOICES,
        default="chrom",
        help=(
            "Processing engine for bedGraph output (default: %(default)s). "
            "'chrom' parallelizes indexed chromosome fetching and "
            "array-backed signal calculation. 'window' parallelizes indexed "
            "coordinate-window fetching for finer load balance on large BAM "
            "inputs.\n\n"
        ),
    )
    parser.add_argument(
        "-ck",
        "--chunk_size",
        dest="chunk_size",
        type=int,
        default=100000,
        help=(
            "Number of records to process per chunk, retained for "
            "compatibility with older workflows (default: %(default)s). "
            "Ignored by public engines.\n\n"
        ),
    )
    parser.add_argument(
        "--chunk-size",
        dest="chunk_size",
        type=int,
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--chnk_size",
        "--chnk-size",
        dest="chunk_size",
        type=int,
        default=argparse.SUPPRESS,
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
        "-uf",
        "--usr_frg",
        dest="user_fragment_length",
        type=int,
        default=None,
        help=(
            "Fixed fragment length to use instead of read lengths (single-end "
            "alignments) or template lengths (paired-end alignments; default: "
            "%(default)s).\n\n"
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
            "values (default: %(default)s). After rounding, non-informative "
            "trailing zeros are stripped.\n\n"
        ),
    )

    parser.add_argument(
        "-pj",
        "--profile_json",
        dest="profile_json",
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-em",
        "--executor_mode",
        dest="executor_mode",
        choices=EXECUTOR_MODE_CHOICES,
        default="map",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-pps",
        "--prototype_parse_strategy",
        dest="prototype_parse_strategy",
        choices=PROTOTYPE_PARSE_STRATEGY_CHOICES,
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-pbs",
        "--prototype_bed_strategy",
        dest="prototype_bed_strategy",
        choices=PROTOTYPE_BED_STRATEGY_CHOICES,
        default="auto",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-ws",
        "--indexed_window_size",
        dest="indexed_window_size",
        type=int,
        default=100000,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-prf",
        "--prototype_result_format",
        dest="prototype_result_format",
        choices=PROTOTYPE_RESULT_FORMAT_CHOICES,
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-pms",
        "--prototype_merge_strategy",
        dest="prototype_merge_strategy",
        choices=PROTOTYPE_MERGE_STRATEGY_CHOICES,
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-pws",
        "--prototype_writer_strategy",
        dest="prototype_writer_strategy",
        choices=PROTOTYPE_WRITER_STRATEGY_CHOICES,
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-pww",
        "--prototype_writer_workers",
        dest="prototype_writer_workers",
        type=int,
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-pwm",
        "--prototype_write_mode",
        dest="prototype_write_mode",
        choices=PROTOTYPE_WRITE_MODE_CHOICES,
        default="full",
        help=argparse.SUPPRESS,
    )

    argv_parse = sys.argv[1:] if argv is None else argv

    if not argv_parse:
        parser.print_help(sys.stderr)
        raise SystemExit(0)

    return parser.parse_args(argv_parse)


def _write_bed_output(
    args: argparse.Namespace,
    fil_out: str,
    header_sizes: dict[str, int],
    chrom_sizes: dict[str, int],
    profile: dict | None,
    time_total_start: float,
) -> int:
    """
    Compute and write fragment-coordinate BED output.

    Parameters
    ----------
    args : argparse.Namespace
        Validated command-line arguments.
    fil_out : str
        Writable BED output path.
    header_sizes : dict[str, int]
        Alignment-header chromosome lengths in source order.
    chrom_sizes : dict[str, int]
        Resolved chromosome lengths selected for output.
    profile : dict | None
        Optional mutable execution profile.
    time_total_start : float
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

    bed_strategy = args.prototype_bed_strategy

    if bed_strategy == "auto":
        if args.threads > 1 and has_indexed_alignment(
            args.fil_in,
            args.ref_fa,
        ):
            bed_strategy = "indexed_chrom"
        else:
            bed_strategy = "serial"

    if profile is not None:
        profile["prototype_bed_strategy_resolved"] = bed_strategy

    if bed_strategy == "serial":
        time_phase = time.perf_counter()
        fragments = parse_bam(
            args.fil_in,
            args.user_fragment_length,
            args.ref_fa,
            chrom_sizes=chrom_sizes,
        )
        record_phase(profile, "parse_bed_fragments", time_phase)

        time_phase = time.perf_counter()
        row_count = write_bed_fragment_dict(fragments, fil_out)
        record_phase(profile, "write_bed", time_phase)

        if profile is not None:
            profile["output_row_count"] = row_count
    else:
        time_phase = time.perf_counter()
        check_indexed_alignment(args.fil_in, args.ref_fa, bed_strategy)
        record_phase(
            profile,
            "validate_indexed_bed_alignment",
            time_phase,
        )

        time_phase = time.perf_counter()
        chrom_task_names = [
            chrom for chrom in header_sizes if chrom in chrom_sizes
        ]

        if bed_strategy == "indexed_chrom":
            tasks = [
                (
                    args.fil_in,
                    chrom,
                    chrom_sizes,
                    args.user_fragment_length,
                    args.ref_fa,
                    None,
                    None,
                )
                for chrom in chrom_task_names
            ]
        else:
            tasks = [
                (
                    args.fil_in,
                    chrom,
                    chrom_sizes,
                    args.user_fragment_length,
                    args.ref_fa,
                    start,
                    min(
                        start + args.indexed_window_size,
                        chrom_sizes[chrom],
                    ),
                )
                for chrom in chrom_task_names
                for start in range(
                    0,
                    chrom_sizes[chrom],
                    args.indexed_window_size,
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
                for task_id, task in enumerate(tasks)
            ]
            profile["n_tasks"] = len(tasks)

        record_phase(profile, "build_indexed_bed_tasks", time_phase)
        time_phase = time.perf_counter()

        try:
            results = list(
                iter_task_results(
                    "indexed_bed",
                    collect_bed_indexed_task,
                    tasks,
                    args.threads,
                    args.executor_mode,
                    profile["tasks"] if profile is not None else None,
                ),
            )
        except OSError:
            if args.prototype_bed_strategy != "auto":
                raise

            if profile is not None:
                profile["prototype_bed_strategy_resolved"] = "serial_fallback"

            fragments = parse_bam(
                args.fil_in,
                args.user_fragment_length,
                args.ref_fa,
                chrom_sizes=chrom_sizes,
            )
            record_phase(profile, "parse_bed_fragments", time_phase)

            time_phase = time.perf_counter()
            row_count = write_bed_fragment_dict(fragments, fil_out)
            record_phase(profile, "write_bed", time_phase)

            if profile is not None:
                profile["output_row_count"] = row_count
                profile["phases_s"]["total"] = (
                    time.perf_counter() - time_total_start
                )
                write_profile(args.profile_json, profile)

            return 0

        record_phase(
            profile,
            "receive_indexed_bed_results",
            time_phase,
        )
        summarize_task_profiles(profile)

        if profile is not None:
            profile["result_payload_bytes"] = sum(
                estimate_bed_payload_bytes(result) for result in results
            )

        time_phase = time.perf_counter()
        row_count = write_bed_array_results(results, fil_out)
        record_phase(
            profile,
            "merge_sort_write_indexed_bed",
            time_phase,
        )

        if profile is not None:
            profile["output_row_count"] = row_count

    if profile is not None:
        profile["phases_s"]["total"] = time.perf_counter() - time_total_start
        write_profile(args.profile_json, profile)

    return 0


def _write_bedgraph_output(
    args: argparse.Namespace,
    signal_data: Any,
    fil_out: str,
    chrom_sizes: dict[str, int],
    profile: dict | None,
    time_total_start: float,
) -> int:
    """
    Write, digest, or profile the merged bedGraph signal.

    Parameters
    ----------
    args : argparse.Namespace
        Validated command-line arguments and output strategy.
    signal_data : Any
        Merged sparse, dense, or dictionary signal representation.
    fil_out : str
        Writable bedGraph output path.
    chrom_sizes : dict[str, int]
        Resolved chromosome lengths keyed by chromosome.
    profile : dict | None
        Optional mutable execution profile.
    time_total_start : float
        Monotonic start time for total-duration reporting.

    Returns
    -------
    status : int
        Zero after output, digest, or profile-only handling succeeds.
    """

    time_phase = time.perf_counter()

    if args.prototype_write_mode == "profile_only":
        if profile is not None:
            profile["output_written"] = False
            profile["output_digest_mode"] = False
            profile["output_row_count"] = signal_row_count(signal_data)
            profile["output_sha256"] = None
    elif args.prototype_write_mode == "digest":
        row_count, output_sha256 = digest_bdg_rows(
            signal_data,
            args.siz_bin,
            args.dp,
            chrom_sizes,
        )

        if profile is not None:
            profile["output_written"] = False
            profile["output_digest_mode"] = True
            profile["output_row_count"] = row_count
            profile["output_sha256"] = output_sha256
    elif args.prototype_writer_strategy == "parallel_ordered":
        write_bdg_parallel_ordered(
            signal_data,
            fil_out,
            args.siz_bin,
            args.dp,
            chrom_sizes,
            args.prototype_writer_workers,
        )
    elif is_array_signal(signal_data):
        write_bdg_array_sparse(
            signal_data,
            fil_out,
            args.siz_bin,
            args.dp,
            chrom_sizes,
        )
    else:
        write_bdg(
            signal_data,
            fil_out,
            args.siz_bin,
            args.dp,
            chrom_sizes=chrom_sizes,
        )

    if profile is not None and args.prototype_write_mode == "full":
        profile["output_written"] = True
        profile["output_digest_mode"] = False
        profile["output_row_count"] = signal_row_count(signal_data)
        profile["output_sha256"] = None

    record_phase(profile, "write_bedgraph", time_phase)

    if profile is not None:
        profile["phases_s"]["total"] = time.perf_counter() - time_total_start
        write_profile(args.profile_json, profile)

    return 0


def _public_signal_results(
    args: argparse.Namespace,
    header_sizes: dict[str, int],
    chrom_sizes: dict[str, int],
    is_len: bool,
    is_norm: bool,
    profile: dict | None,
) -> tuple[Iterable[Any], bool, int]:
    """
    Build serial or indexed tasks for a public signal engine.

    Parameters
    ----------
    args : argparse.Namespace
        Validated command-line arguments and public engine selection.
    header_sizes : dict[str, int]
        Alignment-header chromosome lengths in source order.
    chrom_sizes : dict[str, int]
        Resolved chromosome lengths selected for signal computation.
    is_len : bool
        Whether signal weights use fragment length.
    is_norm : bool
        Whether signal is normalized by total fragments.
    profile : dict | None
        Optional mutable task and timing profile.

    Returns
    -------
    results, is_prototype, fragment_count : tuple[Iterable[Any], bool, int]
        Worker results, prototype-strategy status, and known fragment total.

    Raises
    ------
    ValueError
        If normalization is requested for an empty alignment.
    """

    strategy = PUBLIC_ENGINE_STRATEGY[args.engine]
    is_prototype_strategy = strategy != "serial"

    if strategy == "serial":
        time_phase = time.perf_counter()
        fragment_arrays, fragment_total = collect_fragment_arrays(
            args.fil_in,
            chrom_sizes=chrom_sizes,
            user_fragment_length=args.user_fragment_length,
            ref_fa=args.ref_fa,
        )
        record_phase(profile, "parse_fragment_arrays", time_phase)

        if is_norm and fragment_total == 0:
            raise ValueError(
                (
                    "Normalization requires non-zero total fragments. Check "
                    "the BAM or CRAM file."
                ),
            )

        time_phase = time.perf_counter()
        tasks = [
            (
                chrom,
                arrays[0],
                arrays[1],
                arrays[2],
                chrom_sizes[chrom],
                fragment_total,
                args.siz_bin,
                is_len,
                is_norm,
                args.scl_fct,
            )
            for chrom, arrays in fragment_arrays.items()
        ]

        if profile is not None:
            profile["fragments_total"] = fragment_total
            profile["tasks"] = [
                {
                    "task_id": task_id,
                    "kind": "chrom",
                    "chrom": task[0],
                    "fragments": int(task[1].size),
                    "chrom_size": int(task[4]),
                    "n_bins": math.ceil(task[4] / args.siz_bin),
                    "payload_bytes": int(
                        task[1].nbytes + task[2].nbytes + task[3].nbytes,
                    ),
                }
                for task_id, task in enumerate(tasks)
            ]
            profile["n_tasks"] = len(tasks)

        record_phase(profile, "build_tasks", time_phase)
        results = iter_task_results(
            "chrom",
            calc_sig_array_task,
            tasks,
            args.threads,
            args.executor_mode,
            profile["tasks"] if profile is not None else None,
        )

        return results, is_prototype_strategy, fragment_total

    time_phase = time.perf_counter()
    check_indexed_alignment(args.fil_in, args.ref_fa, strategy)
    record_phase(profile, "validate_indexed_alignment", time_phase)

    time_phase = time.perf_counter()
    chrom_task_names = [
        chrom for chrom in header_sizes if chrom in chrom_sizes
    ]

    if strategy == "indexed_chrom":
        tasks = [
            (
                args.fil_in,
                chrom,
                chrom_sizes,
                args.user_fragment_length,
                args.ref_fa,
                None,
                None,
                args.siz_bin,
                is_len or is_norm,
                args.prototype_result_format,
            )
            for chrom in chrom_task_names
        ]
    else:
        tasks = [
            (
                args.fil_in,
                chrom,
                chrom_sizes,
                args.user_fragment_length,
                args.ref_fa,
                start,
                min(
                    start + args.indexed_window_size,
                    chrom_sizes[chrom],
                ),
                args.siz_bin,
                is_len or is_norm,
                args.prototype_result_format,
            )
            for chrom in chrom_task_names
            for start in range(
                0,
                chrom_sizes[chrom],
                args.indexed_window_size,
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
                "n_bins": math.ceil(chrom_sizes[task[1]] / args.siz_bin),
            }
            for task_id, task in enumerate(tasks)
        ]
        profile["n_tasks"] = len(tasks)

    record_phase(profile, "build_tasks", time_phase)
    results = iter_task_results(
        strategy,
        calc_sig_indexed_fetch_task,
        tasks,
        args.threads,
        args.executor_mode,
        profile["tasks"] if profile is not None else None,
    )

    return results, is_prototype_strategy, 0


def _chunk_signal_results(
    args: argparse.Namespace,
    chrom_sizes: dict[str, int],
    is_len: bool,
    is_norm: bool,
    profile: dict | None,
) -> tuple[Iterable[Any], bool, int]:
    """
    Build fragment-chunk tasks for the compatibility signal engine.

    Parameters
    ----------
    args : argparse.Namespace
        Validated command-line arguments.
    chrom_sizes : dict[str, int]
        Resolved chromosome lengths selected for signal computation.
    is_len : bool
        Whether signal weights use fragment length.
    is_norm : bool
        Whether signal is normalized by total fragments.
    profile : dict | None
        Optional mutable task and timing profile.

    Returns
    -------
    results, is_prototype, fragment_count : tuple[Iterable[Any], bool, int]
        Worker results, false prototype status, and fragment total.

    Raises
    ------
    ValueError
        If normalization is requested for an empty alignment.
    """

    if is_norm:
        time_phase = time.perf_counter()
        fragment_total = count_fragments(
            args.fil_in,
            chrom_sizes=chrom_sizes,
            user_fragment_length=args.user_fragment_length,
            ref_fa=args.ref_fa,
        )
        record_phase(profile, "count_fragments", time_phase)
    else:
        fragment_total = 0

    if is_norm and fragment_total == 0:
        raise ValueError(
            (
                "Normalization requires non-zero total fragments. Check the "
                "BAM or CRAM file."
            ),
        )

    time_phase = time.perf_counter()

    if profile is None:
        tasks = (
            (
                chunk,
                fragment_total,
                args.siz_bin,
                is_len,
                is_norm,
                args.scl_fct,
            )
            for chunk in iter_fragment_chunks(
                alignment_path=args.fil_in,
                chrom_sizes=chrom_sizes,
                chunk_size=args.chunk_size,
                user_fragment_length=args.user_fragment_length,
                ref_fa=args.ref_fa,
            )
        )
    else:
        chunks = list(
            iter_fragment_chunks(
                alignment_path=args.fil_in,
                chrom_sizes=chrom_sizes,
                chunk_size=args.chunk_size,
                user_fragment_length=args.user_fragment_length,
                ref_fa=args.ref_fa,
            ),
        )
        tasks = [
            (
                chunk,
                fragment_total,
                args.siz_bin,
                is_len,
                is_norm,
                args.scl_fct,
            )
            for chunk in chunks
        ]
        profile["fragments_total"] = (
            fragment_total if is_norm else sum(len(chunk) for chunk in chunks)
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
        profile["n_tasks"] = len(tasks)

    record_phase(profile, "build_tasks", time_phase)
    results = iter_task_results(
        "chunk",
        calc_sig_chunk_task,
        tasks,
        args.threads,
        args.executor_mode,
        profile["tasks"] if profile is not None else None,
    )

    return results, False, fragment_total


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

    time_total_start = time.perf_counter()

    args = parse_args(argv)

    if args.fil_in == "-":
        raise SystemExit(
            "'--fil_in -' is no longer supported; provide a BAM or CRAM path.",
        )

    if args.fil_out == "-":
        raise SystemExit(
            "'--fil_out -' is no longer supported; provide an output path.",
        )

    try:
        check_exists(args.fil_in, kind="file", label="Alignment file")

        if args.ref_fa is not None:
            check_exists(args.ref_fa, kind="file", label="Reference FASTA")

        if args.chr_sizes is not None:
            check_exists(args.chr_sizes, kind="file", label="chr.sizes file")
    except FileNotFoundError as error:
        raise SystemExit(str(error)) from None

    try:
        output_path, output_format, _ = validate_output_path(
            args.fil_out,
            ALLOWED_OUTPUT_FORMATS,
        )

        check_writable(output_path, "file")

        if args.profile_json is not None:
            check_writable(args.profile_json, "file")
    except (
        ValueError,
        FileNotFoundError,
        PermissionError,
        IsADirectoryError,
    ) as error:
        raise SystemExit(str(error)) from None

    if output_format != "bed":
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

    try:
        validate_comparison(args.threads, "ge", 1, "threads", allow_none=False)

        if output_format != "bed":
            validate_comparison(
                args.siz_bin,
                "gt",
                0,
                "siz_bin",
                allow_none=False,
            )
            validate_comparison(
                args.chunk_size,
                "gt",
                0,
                "chunk_size",
                allow_none=False,
            )
            validate_comparison(
                args.indexed_window_size,
                "gt",
                0,
                "indexed_window_size",
                allow_none=False,
            )
            validate_comparison(
                args.prototype_writer_workers,
                "ge",
                1,
                "prototype_writer_workers",
                allow_none=False,
            )
            validate_comparison(
                args.prototype_writer_workers,
                "le",
                4,
                "prototype_writer_workers",
                allow_none=False,
            )
            validate_comparison(args.dp, "ge", 0, "dp", allow_none=False)

            if (
                args.prototype_writer_strategy == "serial"
                and args.prototype_writer_workers != 1
            ):
                raise ValueError(
                    "'--prototype_writer_workers' must be 1 when "
                    "'--prototype_writer_strategy serial'.",
                )
        else:
            validate_comparison(
                args.indexed_window_size,
                "gt",
                0,
                "indexed_window_size",
                allow_none=False,
            )

        validate_comparison(args.scl_fct, "gt", 0, "scl_fct", allow_none=True)
        validate_comparison(
            args.user_fragment_length,
            "gt",
            0,
            "usr_frg",
            allow_none=True,
        )

    except ValueError as error:
        raise SystemExit(str(error)) from None

    requested_method = args.method
    args.method = METHOD_CANON[args.method]

    is_len = args.method == "frag"
    is_norm = args.method == "norm"
    profile = start_profile(args, output_path, output_format)

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
            print(f"--fil_out {output_path}")

            if output_format == "bed":
                print(f"--usr_frg {args.user_fragment_length}")
                print(
                    "\n\n(BED output mode: signal computation arguments "
                    "ignored)\n",
                )
            else:
                if requested_method != args.method:
                    method_display = (
                        f"{requested_method}  (standardized internally to "
                        f"{args.method})"
                    )

                    print(f"--method  {method_display}")
                else:
                    print(f"--method  {args.method}")

                print(f"--siz_bin {args.siz_bin}")
                print(f"--engine  {args.engine}")
                print(f"--chunk_size {args.chunk_size}")
                print(f"--scl_fct {args.scl_fct}")
                print(f"--usr_frg {args.user_fragment_length}")
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
        raise SystemExit(f"Alignment file not found: {args.fil_in}") from None
    except ValueError as error:
        raise SystemExit(str(error)) from None
    except OSError as error:
        raise SystemExit(
            f"I/O error while reading BAM or CRAM: {error}",
        ) from None

    try:
        if output_format == "bed":
            _write_bed_output(
                args,
                output_path,
                header_sizes,
                chrom_sizes,
                profile,
                time_total_start,
            )

            return 0

        # Otherwise, compute and write bedGraph signal.
        if args.engine in PUBLIC_ENGINE_STRATEGY:
            (
                results,
                is_prototype_strategy,
                fragment_count,
            ) = _public_signal_results(
                args,
                header_sizes,
                chrom_sizes,
                is_len,
                is_norm,
                profile,
            )

        else:
            (
                results,
                is_prototype_strategy,
                fragment_count,
            ) = _chunk_signal_results(
                args,
                chrom_sizes,
                is_len,
                is_norm,
                profile,
            )

        # Receive worker results, then merge in the parent process.
        time_collect_merge = time.perf_counter()
        time_phase = time.perf_counter()
        result_list = list(results)
        record_phase(profile, "receive_worker_results", time_phase)
        summarize_task_profiles(profile)

        time_phase = time.perf_counter()
        (
            combined_signal,
            merged_fragment_count,
            merged_bin_count,
            result_payload_bytes,
        ) = merge_signal_results(
            results=result_list,
            merge_strategy=args.prototype_merge_strategy,
            chrom_sizes=chrom_sizes,
            siz_bin=args.siz_bin,
            is_prototype_strategy=is_prototype_strategy,
            profile=profile,
        )

        record_phase(profile, "parent_merge_results", time_phase)

        if is_prototype_strategy:
            fragment_count = merged_fragment_count
            time_phase = time.perf_counter()
            apply_signal_adjustments(
                combined_signal,
                fragment_count,
                is_norm,
                args.scl_fct,
            )
            record_phase(profile, "apply_signal_adjustments", time_phase)

        record_phase(profile, "collect_and_merge_results", time_collect_merge)

        if profile is not None:
            if is_prototype_strategy:
                profile["fragments_total"] = fragment_count

            profile["result_bins_before_merge"] = merged_bin_count
            profile["result_bins_after_merge"] = signal_row_count(
                combined_signal,
            )
            profile["result_payload_bytes"] = result_payload_bytes

        _write_bedgraph_output(
            args,
            combined_signal,
            output_path,
            chrom_sizes,
            profile,
            time_total_start,
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
