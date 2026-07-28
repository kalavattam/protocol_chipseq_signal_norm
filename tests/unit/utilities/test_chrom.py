#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_chrom.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


from protocol_chipseq_signal_norm.utilities.utils_chrom import sort_chrom


def test_sort_chrom_orders_roman_numeric_mito_and_other() -> None:
    chroms = ["scaffold", "chr2", "chrX", "chrI", "MT", "chr10"]

    assert sorted(chroms, key=sort_chrom) == [
        "chrI",
        "chrX",
        "chr2",
        "chr10",
        "MT",
        "scaffold",
    ]
