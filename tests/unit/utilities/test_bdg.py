#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_bdg.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


import pytest

from protocol_chipseq_signal_norm.utilities.utils_bdg import (
    canon_nonfinite,
    check_grid_bin,
    check_size_bin,
    generate_name_track,
    iter_rows_bdg,
    load_chr_sizes,
    try_float,
    validate_bounds_bdg,
    write_bdg,
)
from protocol_chipseq_signal_norm.utilities.utils_io import (
    DEF_SKP_PFX,
    is_header,
)


def test_canon_nonfinite_and_try_float():
    assert canon_nonfinite("+Inf") == "inf"
    assert canon_nonfinite("-INF") == "-inf"
    assert canon_nonfinite("1") is None
    assert try_float("1.5") == 1.5
    assert try_float("bad") is None


def test_iter_rows_bdg_skips_headers_and_malformed_rows():
    rows = list(iter_rows_bdg(
        [
            "track name=x\n",
            "chrI 0 10 1\n",
            "chrI 10 10 2\n",
            "chrI a 20 3\n",
            "chrI 20 30 NaN\n",
        ],
        lambda line: is_header(line, DEF_SKP_PFX),
    ))

    assert rows == [
        ("chrI", 0, 10, "1", 1.0),
        ("chrI", 20, 30, "NaN", None),
    ]


def test_check_size_bin_rejects_mismatched_bins(tmp_path):
    file_a = tmp_path / "a.bdg"
    file_b = tmp_path / "b.bdg"
    file_a.write_text("chrI 0 10 1\n", encoding="utf-8")
    file_b.write_text("chrI 0 20 1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Mismatched bin sizes"):
        check_size_bin(str(file_a), str(file_b))


def test_write_bdg_sorts_and_formats_values(tmp_path):
    out = tmp_path / "signal.bdg"

    write_bdg({("chrII", 0): 1.2, ("chrI", 10): -0.0001}, str(out), 10, 2)

    assert out.read_text(encoding="utf-8") == (
        "chrI\t10\t20\t0\n"
        "chrII\t0\t10\t1.2\n"
    )


def test_write_bdg_clamps_final_bin_to_chromosome_size(tmp_path):
    out = tmp_path / "signal.bdg"

    write_bdg(
        {("chrI", 70): 1.0},
        str(out),
        10,
        2,
        chrom_sizes={"chrI": 75},
    )

    assert out.read_text(encoding="utf-8") == "chrI\t70\t75\t1\n"


def test_load_chr_sizes_parses_tsv_and_rejects_duplicates(tmp_path):
    sizes = tmp_path / "chr.sizes"
    sizes.write_text("# comment\nchrI\t80\nchrII 120\n", encoding="utf-8")

    assert load_chr_sizes(str(sizes)) == {"chrI": 80, "chrII": 120}

    sizes.write_text("chrI\t80\nchrI\t90\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate chromosome"):
        load_chr_sizes(str(sizes))


def test_validate_bounds_bdg_rejects_out_of_bounds_rows(tmp_path):
    bdg = tmp_path / "signal.bdg"
    bdg.write_text("chrI\t70\t81\t1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="extends beyond"):
        validate_bounds_bdg(str(bdg), {"chrI": 80})


def test_check_grid_bin_rejects_late_end_mismatch(tmp_path):
    file_a = tmp_path / "a.bdg"
    file_b = tmp_path / "b.bdg"
    file_a.write_text(
        "chrI 0 10 1\nchrI 10 20 1\nchrI 20 30 1\n",
        encoding="utf-8",
    )
    file_b.write_text(
        "chrI 0 10 1\nchrI 10 20 1\nchrI 20 31 1\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Mismatched bedGraph grids"):
        check_grid_bin(str(file_a), str(file_b))


def test_generate_name_track_preserves_gzip_suffix():
    assert generate_name_track("signal.bedGraph.gz") == (
        "signal.track.bedGraph.gz"
    )
