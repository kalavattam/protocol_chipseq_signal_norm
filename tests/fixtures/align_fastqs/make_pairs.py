#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: make_pairs.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Write the pairs reference and reads for the align-fastqs unaligned-read checks.

Two chromosomes of random unique sequence, the first holding three exact copies
of one 300-base repeat, carry one read pair per pair class: concordant,
overlapping, contained, beyond 1,000 bases, dovetailed, outward-facing (RF),
same-strand (FF), on different chromosomes, with one mate unaligned, with both
mates unaligned, and with one mate in the repeat. A single-end set holds
unique, repeat, and unaligned reads. A duplicate set, whose names carry
flow-cell coordinates, holds a concordant pair and pairs with one mate
unaligned whose aligned mates duplicate that pair's first mate, each other, or
each other on the same tile. A filled set adds 300 concordant pairs to the pair
classes, so that BWA can estimate insert sizes and flag pairs proper. Every
random draw is via 'random.Random(0).random()', whose stream Python guarantees
across versions; the filled set draws nothing.
"""

from __future__ import annotations

import argparse
import random

BASES = "ACGT"
LEN_READ = 50
LEN_REP = 300
N_FILL = 300
WIDTH_FASTA = 60

# Read names for the duplicate set: '<instrument>:<run>:<flow cell>:<lane>:'
# then the tile and the x and y coordinates.
PFX_DUP = "A00123:8:HFLOWCELL:1"


def parse_args(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    """
    Parse the output paths.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Write the pairs reference and reads for the align-fastqs "
            "unaligned-read checks."
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
        help="Paired-end R1 FASTQ with one pair per pair class.",
    )
    parser.add_argument(
        "--fil_r2",
        required=True,
        help="Paired-end R2 FASTQ with one pair per pair class.",
    )
    parser.add_argument(
        "--fil_se",
        required=True,
        help="Single-end FASTQ of unique, repeat, and unaligned reads.",
    )
    parser.add_argument(
        "--fil_dup_r1",
        required=True,
        help="Paired-end R1 FASTQ of the duplicate set.",
    )
    parser.add_argument(
        "--fil_dup_r2",
        required=True,
        help="Paired-end R2 FASTQ of the duplicate set.",
    )
    parser.add_argument(
        "--fil_fill_r1",
        required=True,
        help="Paired-end R1 FASTQ of the pair classes and the filler pairs.",
    )
    parser.add_argument(
        "--fil_fill_r2",
        required=True,
        help="Paired-end R2 FASTQ of the pair classes and the filler pairs.",
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


def build_chroms(rng: random.Random) -> dict[str, str]:
    """
    Return the pairs chromosomes by name.

    'chrP1' is 6,000 bases: unique blocks with copies of one repeat at offsets
    1,000, 2,500, and 4,000 (zero-based). 'chrP2' is 3,000 unique bases.

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
    chr_p1 = "".join(
        [
            draw_seq(rng, 1000),
            rep,
            draw_seq(rng, 1200),
            rep,
            draw_seq(rng, 1200),
            rep,
            draw_seq(rng, 1700),
        ]
    )

    return {"chrP1": chr_p1, "chrP2": draw_seq(rng, 3000)}


def cut_fwd(chrom: str, start: int, length: int = LEN_READ) -> str:
    """
    Return the forward-strand read at a zero-based start.

    Parameters
    ----------
    chrom : str
        Chromosome sequence.
    start : int
        Zero-based start of the read.
    length : int
        Read length (default: 'LEN_READ').

    Returns
    -------
    read : str
        Read sequence.
    """

    return chrom[start : start + length]


def cut_rev(chrom: str, start: int, length: int = LEN_READ) -> str:
    """
    Return the reverse-strand read whose alignment starts at 'start'.

    Parameters
    ----------
    chrom : str
        Chromosome sequence.
    start : int
        Zero-based leftmost position of the read's alignment.
    length : int
        Read length (default: 'LEN_READ').

    Returns
    -------
    read : str
        Read sequence, the reverse complement of the reference.
    """

    return reverse_complement(chrom[start : start + length])


def build_pairs(
    rng: random.Random,
    chroms: dict[str, str],
) -> list[tuple[str, str, str]]:
    """
    Return one read pair per pair class as '(name, R1, R2)'.

    Parameters
    ----------
    rng : random.Random
        Seeded generator.
    chroms : dict[str, str]
        Sequence of each chromosome.

    Returns
    -------
    pairs : list[tuple[str, str, str]]
        Name, R1 sequence, and R2 sequence of each pair.
    """

    p1 = chroms["chrP1"]
    p2 = chroms["chrP2"]

    # Fragments start at the R1 offset and end where R2's alignment ends; the
    # repeat copies sit at 1,000, 2,500, and 4,000.
    return [
        ("conc", cut_fwd(p1, 100), cut_rev(p1, 350)),
        ("ovlp", cut_fwd(p1, 500), cut_rev(p1, 520)),
        ("cont", cut_fwd(p1, 700), cut_rev(p1, 720, 30)),
        ("far", cut_fwd(p1, 4400), cut_rev(p1, 5850)),
        ("dove", cut_fwd(p1, 1500), cut_rev(p1, 1490)),
        ("rf", cut_rev(p1, 1700), cut_fwd(p1, 1950)),
        ("ff", cut_fwd(p1, 2000), cut_fwd(p1, 2200)),
        ("chrs", cut_fwd(p1, 3000), cut_rev(p2, 500)),
        ("one", cut_fwd(p1, 3300), draw_seq(rng, LEN_READ)),
        ("none", draw_seq(rng, LEN_READ), draw_seq(rng, LEN_READ)),
        ("rpt", cut_fwd(p1, 2550), cut_rev(p1, 2850)),
        ("conc2", cut_fwd(p2, 1000), cut_rev(p2, 1200)),
    ]


def build_singles(
    rng: random.Random,
    chroms: dict[str, str],
) -> list[tuple[str, str]]:
    """
    Return the single-end reads as '(name, read)'.

    Parameters
    ----------
    rng : random.Random
        Seeded generator.
    chroms : dict[str, str]
        Sequence of each chromosome.

    Returns
    -------
    reads : list[tuple[str, str]]
        Name and sequence of each read.
    """

    p1 = chroms["chrP1"]
    p2 = chroms["chrP2"]

    return [
        ("uniq", cut_fwd(p1, 200)),
        ("uniq2", cut_rev(p2, 2000)),
        ("rpt", cut_fwd(p1, 1100)),
        ("none", draw_seq(rng, LEN_READ)),
    ]


def build_dups(
    rng: random.Random,
    chroms: dict[str, str],
) -> list[tuple[str, str, str]]:
    """
    Return the duplicate set as '(name, R1, R2)'.

    Each name carries a tile and x and y coordinates. 'dupB' is concordant. The
    rest have R2 unaligned: 'dupA1' has R1 at 'dupB''s R1; 'dupA2' and 'dupA3'
    share an R1 position far apart on the tile; 'dupA4' and 'dupA5' share an R1
    position within 2,500 pixels of each other.

    Parameters
    ----------
    rng : random.Random
        Seeded generator.
    chroms : dict[str, str]
        Sequence of each chromosome.

    Returns
    -------
    pairs : list[tuple[str, str, str]]
        Name, R1 sequence, and R2 sequence of each pair.
    """

    p2 = chroms["chrP2"]

    def name(x: int, y: int) -> str:
        return f"{PFX_DUP}:1101:{x}:{y}"

    return [
        (name(1000, 1000), cut_fwd(p2, 1500), cut_rev(p2, 1750)),
        (name(20000, 20000), cut_fwd(p2, 1500), draw_seq(rng, LEN_READ)),
        (name(1000, 30000), cut_fwd(p2, 2200), draw_seq(rng, LEN_READ)),
        (name(25000, 2000), cut_fwd(p2, 2200), draw_seq(rng, LEN_READ)),
        (name(5000, 5000), cut_fwd(p2, 2600), draw_seq(rng, LEN_READ)),
        (name(5100, 5100), cut_fwd(p2, 2600), draw_seq(rng, LEN_READ)),
    ]


def build_fills(chroms: dict[str, str]) -> list[tuple[str, str, str]]:
    """
    Return 'N_FILL' concordant filler pairs as '(name, R1, R2)'.

    The pairs lie in the unique blocks, away from the repeat copies, with R1
    forward and R2 reverse. Fragment lengths cycle from 200 to 300 bases in
    steps of 10, and each pair starts at its own position, ending in 5, so no
    filler pair duplicates another or a pair class. Nothing is drawn at random.

    Parameters
    ----------
    chroms : dict[str, str]
        Sequence of each chromosome.

    Returns
    -------
    pairs : list[tuple[str, str, str]]
        Name, R1 sequence, and R2 sequence of each pair.
    """

    # Unique blocks as '(chromosome, start, end)', zero-based and half-open.
    blocks = [
        ("chrP1", 0, 1000),
        ("chrP1", 1300, 2500),
        ("chrP1", 2800, 4000),
        ("chrP1", 4300, 6000),
        ("chrP2", 0, 3000),
    ]
    starts = []
    for nam, beg, end in blocks:
        for pos in range(beg + 5, end - 300, 10):
            starts.append((nam, pos))

    step = len(starts) // N_FILL
    pairs = []
    for idx in range(N_FILL):
        nam, pos = starts[idx * step]
        len_frg = 200 + 10 * (idx % 11)
        chrom = chroms[nam]
        pairs.append(
            (
                f"fill_{idx:03d}",
                cut_fwd(chrom, pos),
                cut_rev(chrom, pos + len_frg - LEN_READ),
            )
        )

    return pairs


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
    Write the reference, then the pairs, single-end, duplicate, and filled
    reads.
    """

    args = parse_args(argv)
    rng = random.Random(0)
    chroms = build_chroms(rng)
    pairs = build_pairs(rng, chroms)
    singles = build_singles(rng, chroms)
    dups = build_dups(rng, chroms)
    fills = [(f"pair_{n}", r1, r2) for n, r1, r2 in pairs]
    fills += build_fills(chroms)

    write_fasta(args.fil_ref, chroms)
    write_fastq(args.fil_r1, [(f"pair_{n}", r1) for n, r1, _ in pairs])
    write_fastq(args.fil_r2, [(f"pair_{n}", r2) for n, _, r2 in pairs])
    write_fastq(args.fil_se, [(f"se_{n}", rd) for n, rd in singles])
    write_fastq(args.fil_dup_r1, [(n, r1) for n, r1, _ in dups])
    write_fastq(args.fil_dup_r2, [(n, r2) for n, _, r2 in dups])
    write_fastq(args.fil_fill_r1, [(n, r1) for n, r1, _ in fills])
    write_fastq(args.fil_fill_r2, [(n, r2) for n, _, r2 in fills])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
