#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_chrom.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5-series models) was used in development.
#
# Distributed under the MIT license.

"""
Script
------
utils_chrom.py


Description
-----------
Helper function for chromosome/contig naming with yeast-aware ordering.
Provides a stable sort key that ranks S. cerevisiae/S. pombe Roman-numeral
chromosome names (I..XVI) numerically, then purely numeric names, then common
mitochondrial aliases, and finally any remaining scaffolds lexically.

This keeps bedGraph or BED outputs in biologically familiar order without
re-writing names.


Variables
---------
ROMAN_TO_INT
    Mapping of Roman numerals I..XVI → 1..16 used to sort S. cerevisiae and S.
    pombe chromosomes numerically. Matching is case-insensitive and ignores an
    optional 'chr' prefix (e.g., 'chriv' == 'IV' == 4).

MITO_NAMES
    Set of labels treated as names of mitochondrial references: {'M', 'MT',
    'MITO', 'MITOCHONDRION'}. Matching is case-insensitive and ignores an
    optional 'chr' prefix, so 'chrM' and 'm' are equivalent.

INT_MAX
    Large integer sentinel (default 10**18) used as the numeric tie-breaker
    for tiers without a natural number (mitochondrial and “other” scaffolds or
    contigs). Ensures those tiers sort after true numeric chromosomes; override
    if contig naming scheme requires a larger sentinel.


Functions
---------
sort_chrom()
    Return a sortable 3-tuple '(tier, num_key, lex_key)' for a chromosome or
    contig name. Roman numerals sort by Arabic integer values, decimal names
    sort by their integer values, mitochondrial aliases are grouped after
    numerics, and everything else sorts lexically. Case-insensitive; optional
    'chr' prefix is ignored in sorting.


Example
-------
'''python
from scripts.functions.utils_chrom import sort_chrom

chroms = ["chrII", "chr10", "chrM", "scaffold_42", "V", "chrXVI", "2"]
for c in sorted(chroms, key=sort_chrom):
    print(c)
'''

'''ipython
In [4]: for c in sorted(chroms, key=sort_chrom):
   ...:     print(c)
chrII
V
chrXVI
2
chr10
chrM
scaffold_42
'''


Notes
-----
- Roman numerals I..XVI are recognized case-insensitively (with or without a
  ‘chr’ prefix) via 'ROMAN_TO_INT'.
- 'MITO_NAMES' collects common mitochondrial aliases; the 'chr' prefix is
  stripped before matching, so 'chrM' and 'M' are treated the same.
- 'INT_MAX' is an integer sentinel for tiers without a numeric tie-breaker
  (mitochondrial and “other”); the default 1e18 is safe for mammalian-scale
  genomes and fine for yeast. Override via the 'int_max' parameter only if
  an even larger sentinel is needed.
- The returned key is suitable for Python's default tuple ordering and is used
  by 'utils_bdg.key_bin()' to keep bins in chromosome order first, then by
  start.
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

    Args:
        chrom : str
            Chromosome/contig name (e.g., 'X', 'chrII', 'chr10', '2', 'MT',
            'scaffold_42', etc.).
        int_max : int = INT_MAX
            Integer sentinel used as the numeric tie-breaker for tiers that
            do not have a real number (mitochondrial and “other” scaffolds).
            Larger values keep those tiers sorted after purely numeric
            chromosomes; override only if you need an even bigger sentinel for
            extreme contig ranges.

    Returns:
        (tier, num_key, lex_key) : tuple[int, int, str]
            A 3-tuple sort key.

    Notes:
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
