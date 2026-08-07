#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_chrom.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6) were
# used in design, development, and documentation, with all output reviewed,
# edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Chromosome sorting helpers.

Notes
-----
Sorting prioritizes yeast-style Roman numerals, then numeric names,
mitochondrial aliases, and lexical scaffold names.
"""

# Map Roman numerals used in S. cerevisiae and S. pombe chromosome names.
ROMAN_TO_INT: dict[str, int] = {
    "I": 1,
    "II": 2,
    "III": 3,
    "IV": 4,
    "V": 5,
    "VI": 6,
    "VII": 7,
    "VIII": 8,
    "IX": 9,
    "X": 10,
    "XI": 11,
    "XII": 12,
    "XIII": 13,
    "XIV": 14,
    "XV": 15,
    "XVI": 16,
}

# Recognize common mitochondrial-reference names.
MITO_NAMES = {"M", "MT", "MITO", "MITOCHONDRION"}

# Integer sentinel used when numeric tie-breaker is unused (tiers 2–3); safe
# for all typical yeast genomes.
#
# One quintillion bases also safely covers maintained yeast workflows.
INT_MAX = 10**18  # 10**12


def sort_chrom(chrom: str, int_max: int = INT_MAX) -> tuple[int, int, str]:
    """
    Return a sortable key for chromosome names with the following precedence:.

        1. Roman numerals I..XVI (case-insensitive; optional 'chr' prefix),
           ordered numerically by Arabic whole-number value (e.g., I < II <
           ... < XVI).
        2. Purely numeric names (with or without 'chr' prefix) ordered by
           integer value (e.g., 1 < 2 < ... < 10 < 11 ...).
        3. Common mitochondrial labels grouped after numerics.
        4. Everything else sorted lexically.

    Parameters
    ----------
    chrom : str
        Chromosome/contig name (e.g., 'X', 'chrII', 'chr10', '2', 'MT',
        'scaffold_42', etc.).
    int_max : int = INT_MAX
        Integer sentinel used as the numeric tie-breaker for tiers that
        do not have a real number (mitochondrial and “other” scaffolds).
        Larger values keep those tiers sorted after purely numeric
        chromosomes; override only if you need an even bigger sentinel for
        extreme contig ranges.

    Returns
    -------
    key : tuple[int, int, str]
        A 3-tuple sort key.

    Notes
    -----
    - tier: precedence bucket (lower sorts first).
    - num_key: numeric tie-breaker within tiers that use numbers; 'int_max'
               for tiers without a natural numeric key.
    - lex_key: lexical tie-breaker within the catch-all tier; "" otherwise.
    """

    key = chrom[3:] if chrom.lower().startswith("chr") else chrom
    up = key.upper()

    # First, sort Roman numerals I–XVI as whole numbers.
    r = ROMAN_TO_INT.get(up)
    if r is not None:
        return (0, r, "")

    # Next, sort decimal chromosome names numerically.
    if key.isdigit():
        return (1, int(key), "")

    # Then, place mitochondrial references after numeric chromosomes.
    if up in MITO_NAMES:
        return (2, int_max, "")  # The lexical tie-breaker is unused.

    # Finally, sort all remaining names lexically.
    return (3, int_max, key)
