#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: compare_bam_blocks.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Compare the BGZF blocks that hold two BAM files' records.

The blocks after the BAM header hold the records. Two files written by the same
htslib from the same records at the same compression level hold the same bytes
there, whatever their headers, so comparing those blocks tells whether a file
was written at a given level: rewrite it at that level and compare. The header
must end on a block boundary, as htslib writes it. Empty blocks are left out,
since 'samtools reheader' can leave two end-of-file blocks. Only the header's
blocks are inflated; the rest are read as a stream and hashed.

Prints 'same <n>' or 'differ <n>', where <n> is the number of blocks with
records in the first file.
"""

from __future__ import annotations

import argparse
import hashlib
import struct
import sys
import zlib

# Bytes of the BGZF block header before its extra field, and of the gzip
# trailer (CRC32 and ISIZE).
LEN_HEAD = 12
LEN_TAIL = 8


def parse_args(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    """
    Parse the two BAM paths.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Compare the BGZF blocks that hold two BAM files' records, and "
            "print 'same <n>' or 'differ <n>'."
        ),
    )
    parser.add_argument("bam_a", help="First BAM file.")
    parser.add_argument("bam_b", help="Second BAM file.")

    return parser.parse_args(argv)


def read_block(fh) -> bytes | None:
    """
    Read one whole BGZF block from an open file.

    Parameters
    ----------
    fh : BinaryIO
        BAM file opened for binary reading, at the start of a block.

    Returns
    -------
    block : bytes or None
        The block's bytes, or None at the end of the file.
    """

    head = fh.read(LEN_HEAD)
    if not head:
        return None

    if len(head) < LEN_HEAD or head[:2] != b"\x1f\x8b":
        raise ValueError("not a BGZF block")

    len_xtra = struct.unpack_from("<H", head, 10)[0]
    xtra = fh.read(len_xtra)
    size = None
    idx = 0
    while idx < len_xtra:
        len_sub = struct.unpack_from("<H", xtra, idx + 2)[0]
        if xtra[idx : idx + 2] == b"BC":
            size = struct.unpack_from("<H", xtra, idx + 4)[0]
        idx += 4 + len_sub

    if size is None:
        raise ValueError("BGZF block lacks its size field ('BC')")

    rest = fh.read(size + 1 - LEN_HEAD - len_xtra)

    return head + xtra + rest


def measure_header(data: bytes) -> int | None:
    """
    Return the length of the BAM header, or None if more bytes are needed.

    Parameters
    ----------
    data : bytes
        Inflated bytes from the start of the BAM file.

    Returns
    -------
    length : int or None
        Bytes of magic, header text, and reference list.
    """

    if len(data) < 8:
        return None

    pos = 8 + struct.unpack_from("<i", data, 4)[0]
    if len(data) < pos + 4:
        return None

    n_ref = struct.unpack_from("<i", data, pos)[0]
    pos += 4
    for _ in range(n_ref):
        if len(data) < pos + 4:
            return None

        pos += 4 + struct.unpack_from("<i", data, pos)[0] + 4

    if len(data) < pos:
        return None

    return pos


def hash_record_blocks(path: str) -> tuple[str, int]:
    """
    Hash the non-empty BGZF blocks after a BAM file's header.

    Parameters
    ----------
    path : str
        BAM file.

    Returns
    -------
    digest : str
        SHA-256 of the blocks' bytes, in order.
    n_blk : int
        Number of blocks hashed.
    """

    sha = hashlib.sha256()
    n_blk = 0
    hdr = b""
    end = None

    with open(path, "rb") as fh:
        while (block := read_block(fh)) is not None:
            len_xtra = struct.unpack_from("<H", block, 10)[0]

            # Inflate blocks until the header is complete.
            if end is None or len(hdr) < end:
                body = block[LEN_HEAD + len_xtra : -LEN_TAIL]
                hdr += zlib.decompress(body, -15)
                end = measure_header(hdr)
                if end is not None and len(hdr) > end:
                    raise ValueError(
                        f"the header of '{path}' does not end on a block "
                        "boundary"
                    )

                continue

            # Leave out empty blocks, such as end-of-file markers.
            if struct.unpack_from("<I", block, len(block) - 4)[0] == 0:
                continue

            sha.update(block)
            n_blk += 1

    return sha.hexdigest(), n_blk


def main(argv: list[str] | None = None) -> int:
    """
    Print 'same <n>' or 'differ <n>' for the two files.
    """

    args = parse_args(argv)

    try:
        dig_a, n_a = hash_record_blocks(args.bam_a)
        dig_b, _ = hash_record_blocks(args.bam_b)
    except (OSError, ValueError, zlib.error) as err:
        print(f"error(compare_bam_blocks.py): {err}", file=sys.stderr)
        return 1

    print(f"{'same' if dig_a == dig_b else 'differ'} {n_a}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
