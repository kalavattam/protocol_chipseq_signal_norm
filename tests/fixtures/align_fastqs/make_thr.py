#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: make_thr.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Write the thread reference and reads for the align-fastqs determinism checks.

One chromosome of 100,000 random bases carries 120,000 read pairs of 50-base
mates: first 100,000 pairs from 190- to 210-base fragments, one in every 50 of
them replaced by a pair of random reads that do not align, then 20,000 pairs
from 590- to 610-base fragments. Each pair is on either strand.

The set is larger than the 10,000,000 bases that BWA-MEM reads per batch at one
thread without '-K', so its first batch holds only the short fragments and its
second only the long ones. Each batch then estimates insert sizes alone, and
the long pairs are proper. At three threads, one batch holds both, and the long
pairs are not. So the thread count moves pairs between files unless '-K' fixes
the batch size. Every random draw is via 'random.Random(0).random()', whose
stream Python guarantees across versions.
"""

from __future__ import annotations

import argparse
import random

BASES = "ACGT"
LEN_CHROM = 100_000
LEN_READ = 50
WIDTH_FASTA = 60

# Fragment blocks, in order: number of pairs and the fragment-length range.
BLOCKS = ((100_000, 190, 210), (20_000, 590, 610))

# In the first block, every 'EVERY_UNAL'-th pair is two random reads.
EVERY_UNAL = 50


def parse_args(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    """
    Parse the output paths.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Write the thread reference and reads for the align-fastqs "
            "determinism checks."
        ),
    )
    parser.add_argument(
        "--fil_ref",
        required=True,
        help="Reference FASTA to write.",
    )
    parser.add_argument(
        "--fil_r1",
        required=True,
        help="Paired-end R1 FASTQ.",
    )
    parser.add_argument(
        "--fil_r2",
        required=True,
        help="Paired-end R2 FASTQ.",
    )

    return parser.parse_args(argv)


def draw_seq(rng: random.Random, length: int) -> str:
    """
    Return a random sequence of the given length.

    Parameters
    ----------
    rng : random.Random
        Seeded generator.
    length : int
        Sequence length in bases.

    Returns
    -------
    seq : str
        Sequence over 'ACGT'.
    """

    return "".join(BASES[int(rng.random() * 4)] for _ in range(length))


def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a sequence.

    Parameters
    ----------
    seq : str
        Sequence over 'ACGT'.

    Returns
    -------
    seq_rc : str
        Reverse complement of 'seq'.
    """

    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def build_pairs(
    rng: random.Random,
    chrom: str,
) -> list[tuple[str, str, str]]:
    """
    Return the read pairs, in order, as name, R1, and R2.

    Parameters
    ----------
    rng : random.Random
        Seeded generator, after the chromosome's draws.
    chrom : str
        Chromosome sequence.

    Returns
    -------
    pairs : list[tuple[str, str, str]]
        Each pair's name, R1 sequence, and R2 sequence.
    """

    pairs = []

    for idx_blk, (n_pair, len_min, len_max) in enumerate(BLOCKS):
        for i in range(n_pair):
            len_frg = len_min + int(rng.random() * (len_max - len_min + 1))
            start = int(rng.random() * (len(chrom) - len_frg))
            r1 = chrom[start : start + LEN_READ]
            r2 = reverse_complement(
                chrom[start + len_frg - LEN_READ : start + len_frg]
            )

            if rng.random() < 0.5:
                r1, r2 = r2, r1

            if idx_blk == 0 and i % EVERY_UNAL == 0:
                r1 = draw_seq(rng, LEN_READ)
                r2 = draw_seq(rng, LEN_READ)

            pairs.append((f"thr_{len(pairs):06d}", r1, r2))

    return pairs


def write_fasta(path: str, name: str, seq: str) -> None:
    """
    Write one chromosome as FASTA, 'WIDTH_FASTA' bases per line.

    Parameters
    ----------
    path : str
        Output path.
    name : str
        Chromosome name.
    seq : str
        Chromosome sequence.
    """

    with open(path, "w", encoding="ascii") as out:
        out.write(f">{name}\n")

        for i in range(0, len(seq), WIDTH_FASTA):
            out.write(f"{seq[i : i + WIDTH_FASTA]}\n")


def write_fastq(path: str, reads: list[tuple[str, str]]) -> None:
    """
    Write named reads with uniform base qualities.

    Parameters
    ----------
    path : str
        Output path.
    reads : list[tuple[str, str]]
        Name and sequence of each read, in order.
    """

    with open(path, "w", encoding="ascii") as out:
        for name, read in reads:
            out.write(f"@{name}\n{read}\n+\n{'I' * len(read)}\n")


def main(argv: list[str] | None = None) -> int:
    """
    Write the reference, then the paired reads.
    """

    args = parse_args(argv)
    rng = random.Random(0)
    chrom = draw_seq(rng, LEN_CHROM)
    pairs = build_pairs(rng, chrom)

    write_fasta(args.fil_ref, "chrH", chrom)
    write_fastq(args.fil_r1, [(n, r1) for n, r1, _ in pairs])
    write_fastq(args.fil_r2, [(n, r2) for n, _, r2 in pairs])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
