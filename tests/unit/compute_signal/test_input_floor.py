#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_input_floor.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


import pytest

from protocol_chipseq_signal_norm.cli.compute_input_floor import (
    compute_input_floor,
    count_aln_bed,
    infer_fmt,
    parse_flag_csv,
)
from protocol_chipseq_signal_norm.utilities.utils_io import DEF_SKP_PFX


def test_parse_flag_csv_accepts_decimal_and_hex():
    assert parse_flag_csv("99, 0x400", "flags") == {99, 1024}


def test_parse_flag_csv_rejects_invalid_token():
    with pytest.raises(SystemExit):
        parse_flag_csv("99, bad", "flags")


def test_parse_flag_csv_rejects_empty_list():
    with pytest.raises(SystemExit):
        parse_flag_csv(",", "flags")


def test_infer_fmt_recognizes_supported_suffixes_and_stdin_hints():
    assert infer_fmt("x.bam") == "bam"
    assert infer_fmt("x.bed.gz") == "bed"
    assert infer_fmt("x.bdg") == "bedgraph"
    assert infer_fmt("-", "bedGraph") == "bedgraph"
    assert infer_fmt("x.txt") == "other"


def test_count_aln_bed_skips_headers_and_blanks(tmp_path):
    path = tmp_path / "alignments.bed"
    path.write_text(
        "# h\n"
        "track name=x\n"
        "chrI\t0\t10\n"
        "\n"
        "   \n"
        "chrI\t10\t20\n",
        encoding="utf-8",
    )

    assert count_aln_bed(str(path), DEF_SKP_PFX) == 2


def test_compute_input_floor_norm_and_dist_modes(tmp_path):
    bdg = tmp_path / "input.bdg"
    bdg.write_text(
        "chrI 0 10 0\n"
        "chrI 10 20 -1\n"
        "chrI 20 30 1\n"
        "chrI 30 40 2\n"
        "chrI 40 50 10\n",
        encoding="utf-8",
    )

    assert compute_input_floor("-", 10, 100, mode="norm") == pytest.approx(
        1 / 9
    )
    assert compute_input_floor(
        str(bdg),
        10,
        100,
        mode="dist",
        method="qntl_nz",
        qntl_nz=50,
    ) == 2.0


def test_compute_input_floor_rejects_bin_size_at_or_above_genome_size():
    with pytest.raises(SystemExit):
        compute_input_floor("-", 100, 100, mode="norm")
