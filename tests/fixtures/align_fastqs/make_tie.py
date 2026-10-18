#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: make_tie.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Write the tie reference and reads for the align-fastqs determinism checks.

Two chromosomes of random unique sequence each hold three exact copies of one
600-base repeat, so a read inside a copy has six equally good placements. The
same 1,000 read pairs are written under three header forms, plus two single-end
sets of the R1 reads. Every random draw is via 'random.Random(0).random()',
whose stream Python guarantees across versions.
"""

from __future__ import annotations

import argparse
import random

BASES = "ACGT"
LEN_UNIQ = 3000
LEN_REP = 600
LEN_READ = 50
LEN_FRAG_MIN = 150
LEN_FRAG_MAX = 400
N_PAIR = 1000
P_MISMATCH = 0.3
WIDTH_FASTA = 60


def parse_args(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    """
    Parse the output paths.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Write the tie reference and reads for the align-fastqs "
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
        help="Paired-end R1 FASTQ with bare names.",
    )
    parser.add_argument(
        "--fil_r2",
        required=True,
        help="Paired-end R2 FASTQ with bare names.",
    )
    parser.add_argument(
        "--fil_cmt_r1",
        required=True,
        help="Paired-end R1 FASTQ with CASAVA comments.",
    )
    parser.add_argument(
        "--fil_cmt_r2",
        required=True,
        help="Paired-end R2 FASTQ with CASAVA comments.",
    )
    parser.add_argument(
        "--fil_idx_r1",
        required=True,
        help="Paired-end R1 FASTQ with CASAVA comments and another index.",
    )
    parser.add_argument(
        "--fil_idx_r2",
        required=True,
        help="Paired-end R2 FASTQ with CASAVA comments and another index.",
    )
    parser.add_argument(
        "--fil_se",
        required=True,
        help="Single-end FASTQ of the R1 reads with bare names.",
    )
    parser.add_argument(
        "--fil_cmt_se",
        required=True,
        help="Single-end FASTQ of the R1 reads with CASAVA comments.",
    )

    return parser.parse_args(argv)


def draw_index(rng: random.Random, n: int) -> int:
    """
    Return a random index below 'n'.

    Parameters
    ----------
    rng : random.Random
        Seeded generator.
    n : int
        Number of choices.

    Returns
    -------
    index : int
        Index in the range from 0 to 'n' minus 1.
    """

    return int(rng.random() * n)


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

    return "".join(BASES[draw_index(rng, 4)] for _ in range(length))


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


def add_mismatch(rng: random.Random, read: str) -> str:
    """
    Return the read with one base changed, at probability 'P_MISMATCH'.

    Parameters
    ----------
    rng : random.Random
        Seeded generator.
    read : str
        Read sequence.

    Returns
    -------
    read_out : str
        The read, unchanged or with one substituted base.
    """

    if rng.random() >= P_MISMATCH:
        return read

    pos = draw_index(rng, len(read))
    alts = [base for base in BASES if base != read[pos]]
    alt = alts[draw_index(rng, len(alts))]

    return read[:pos] + alt + read[pos + 1 :]


def build_chroms(rng: random.Random) -> dict[str, str]:
    """
    Return the tie chromosomes by name.

    Each chromosome is three unique blocks, each followed by a copy of one
    shared repeat, then a last unique block a third as long.

    Parameters
    ----------
    rng : random.Random
        Seeded generator.

    Returns
    -------
    chroms : dict[str, str]
        Sequence of each chromosome, in reference order.
    """

    rep = draw_seq(rng, LEN_REP)
    chroms: dict[str, str] = {}

    for name in ("chrT1", "chrT2"):
        parts: list[str] = []

        for _ in range(3):
            parts.append(draw_seq(rng, LEN_UNIQ))
            parts.append(rep)

        parts.append(draw_seq(rng, LEN_UNIQ // 3))
        chroms[name] = "".join(parts)

    return chroms


def draw_pairs(
    rng: random.Random,
    chroms: dict[str, str],
) -> tuple[list[str], list[str]]:
    """
    Return R1 and R2 reads cut from random fragments on random strands.

    R1 reads the fragment's start and R2 the reverse complement of its end.

    Parameters
    ----------
    rng : random.Random
        Seeded generator.
    chroms : dict[str, str]
        Sequence of each chromosome.

    Returns
    -------
    reads_r1, reads_r2 : tuple[list[str], list[str]]
        R1 and R2 reads, in pair order.
    """

    names = list(chroms)
    reads_r1: list[str] = []
    reads_r2: list[str] = []

    for _ in range(N_PAIR):
        chrom = chroms[names[draw_index(rng, len(names))]]
        len_frag = LEN_FRAG_MIN + draw_index(
            rng,
            LEN_FRAG_MAX - LEN_FRAG_MIN + 1,
        )
        start = draw_index(rng, len(chrom) - len_frag + 1)
        frag = chrom[start : start + len_frag]

        if rng.random() < 0.5:
            frag = reverse_complement(frag)

        reads_r1.append(add_mismatch(rng, frag[:LEN_READ]))
        reads_r2.append(add_mismatch(rng, reverse_complement(frag)[:LEN_READ]))

    return reads_r1, reads_r2


def write_fasta(path: str, chroms: dict[str, str]) -> None:
    """
    Write chromosomes as FASTA, 'WIDTH_FASTA' bases per line.

    Parameters
    ----------
    path : str
        Output path.
    chroms : dict[str, str]
        Sequence of each chromosome.
    """

    with open(path, "w", encoding="ascii") as out:
        for name, seq in chroms.items():
            out.write(f">{name}\n")

            for i in range(0, len(seq), WIDTH_FASTA):
                out.write(f"{seq[i : i + WIDTH_FASTA]}\n")


def write_fastq(path: str, reads: list[str], suffix: str) -> None:
    """
    Write reads named 'tie_<n>' plus 'suffix', with uniform qualities.

    Parameters
    ----------
    path : str
        Output path.
    reads : list[str]
        Read sequences, in order.
    suffix : str
        Text appended to each name, such as a comment after a space.
    """

    with open(path, "w", encoding="ascii") as out:
        for i, read in enumerate(reads):
            out.write(f"@tie_{i}{suffix}\n{read}\n+\n{'I' * len(read)}\n")


def main(argv: list[str] | None = None) -> int:
    """
    Write the reference, then the same reads under each header form.
    """

    args = parse_args(argv)
    rng = random.Random(0)
    chroms = build_chroms(rng)
    reads_r1, reads_r2 = draw_pairs(rng, chroms)

    write_fasta(args.fil_ref, chroms)

    # CASAVA comments differ between mates in the read number; the second
    # index sequence is a change both mates share.
    write_fastq(args.fil_r1, reads_r1, "")
    write_fastq(args.fil_r2, reads_r2, "")
    write_fastq(args.fil_cmt_r1, reads_r1, " 1:N:0:AAAAAAAA")
    write_fastq(args.fil_cmt_r2, reads_r2, " 2:N:0:AAAAAAAA")
    write_fastq(args.fil_idx_r1, reads_r1, " 1:N:0:CCCCCCCC")
    write_fastq(args.fil_idx_r2, reads_r2, " 2:N:0:CCCCCCCC")
    write_fastq(args.fil_se, reads_r1, "")
    write_fastq(args.fil_cmt_se, reads_r1, " 1:N:0:AAAAAAAA")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
