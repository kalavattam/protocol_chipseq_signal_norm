#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_chrom.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6) were
# used in development and documentation.
#
# Distributed under the MIT license.


"""
Chromosome sorting helpers.

Notes
-----
Sorting prioritizes yeast-style Roman numerals, then numeric names,
mitochondrial aliases, and lexical scaffold names.
"""

#  Roman numerals used in S. cerevisiae and S. pombe chromosome names
ROMAN_TO_INT: dict[str, int] = {
    "I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6, "VII": 7,
    "VIII": 8, "IX": 9, "X": 10, "XI": 11, "XII": 12, "XIII": 13,
    "XIV": 14, "XV": 15, "XVI": 16
}

#  Common names for mitochondrial references
MITO_NAMES = {"M", "MT", "MITO", "MITOCHONDRION"}

#  Integer sentinel used when numeric tie-breaker is unused (tiers 2–3); safe
#  for all typical yeast genomes
#
#  1e12 bp is appropriate for yeast genomes; 1e18 bp is appropriate for
#  mammalian genomes (and is actually fine to use for yeast genomes too)
INT_MAX = 10**18  # 10**12


def sort_chrom(chrom: str, int_max: int = INT_MAX) -> tuple[int, int, str]:
    """
    Return a sortable key for chromosome names with the following precedence:
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
        (tier, num_key, lex_key) : tuple[int, int, str]
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

    #  First, sort Roman numerals I..XVI as if they were Arabic whole numbers
    r = ROMAN_TO_INT.get(up)
    if r is not None:
        return (0, r, "")

    #  Next, numerically sort whole numbers
    if key.isdigit():
        return (1, int(key), "")

    #  Organize “Mito” group after the numerics
    if up in MITO_NAMES:
        return (2, int_max, "")  # Tie-breaker unused

    #  Finally, sort everything else lexically
    return (3, int_max, key)
