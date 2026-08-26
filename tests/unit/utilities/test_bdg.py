#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_bdg.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


from pathlib import Path

import pytest

from protocol_chipseq_signal_norm.utilities.utils_bdg import (
    canon_nonfinite,
    check_grid_bin,
    check_size_bin,
    generate_name_track,
    iter_rows_bdg,
    load_chromosome_sizes,
    sum_counts_bdg,
    try_float,
    validate_bounds_bdg,
    write_bdg,
)
from protocol_chipseq_signal_norm.utilities.utils_io import (
    DEF_SKP_PFX,
    is_header,
)


def test_canon_nonfinite_and_try_float() -> None:
    assert canon_nonfinite("+Inf") == "inf"
    assert canon_nonfinite("-INF") == "-inf"
    assert canon_nonfinite("1") is None
    assert try_float("1.5") == 1.5
    assert try_float("bad") is None


def test_iter_rows_bdg_skips_headers_and_malformed_rows() -> None:
    rows = list(
        iter_rows_bdg(
            [
                "track name=x\n",
                "chrI 0 10 1\n",
                "chrI 10 10 2\n",
                "chrI a 20 3\n",
                "chrI 20 30 NaN\n",
            ],
            lambda line: is_header(line, DEF_SKP_PFX),
        ),
    )

    assert rows == [
        ("chrI", 0, 10, "1", 1.0),
        ("chrI", 20, 30, "NaN", None),
    ]


def test_check_size_bin_rejects_mismatched_bins(tmp_path: Path) -> None:
    file_a = tmp_path / "a.bdg"
    file_b = tmp_path / "b.bdg"
    file_a.write_text("chrI 0 10 1\n", encoding="utf-8")
    file_b.write_text("chrI 0 20 1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Mismatched bin sizes"):
        check_size_bin(str(file_a), str(file_b))


def test_write_bdg_sorts_and_formats_values(tmp_path: Path) -> None:
    out = tmp_path / "signal.bdg"

    write_bdg({("chrII", 0): 1.2, ("chrI", 10): -0.0001}, str(out), 10, 2)

    assert out.read_text(encoding="utf-8") == (
        "chrI\t10\t20\t0\nchrII\t0\t10\t1.2\n"
    )


def test_write_bdg_clamps_final_bin_to_chromosome_size(tmp_path: Path) -> None:
    out = tmp_path / "signal.bdg"

    write_bdg(
        {("chrI", 70): 1.0},
        str(out),
        10,
        2,
        chrom_sizes={"chrI": 75},
    )

    assert out.read_text(encoding="utf-8") == "chrI\t70\t75\t1\n"


def test_load_chromosome_sizes_parses_tsv_and_rejects_duplicates(
    tmp_path: Path,
) -> None:
    sizes = tmp_path / "chr.sizes"
    sizes.write_text("# comment\nchrI\t80\nchrII 120\n", encoding="utf-8")

    assert load_chromosome_sizes(str(sizes)) == {"chrI": 80, "chrII": 120}

    sizes.write_text("chrI\t80\nchrI\t90\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Duplicate chromosome"):
        load_chromosome_sizes(str(sizes))


def test_validate_bounds_bdg_rejects_out_of_bounds_rows(
    tmp_path: Path,
) -> None:
    bdg = tmp_path / "signal.bdg"
    bdg.write_text("chrI\t70\t81\t1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="extends beyond"):
        validate_bounds_bdg(str(bdg), {"chrI": 80})


def test_check_grid_bin_rejects_late_end_mismatch(tmp_path: Path) -> None:
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


def test_generate_name_track_preserves_gzip_suffix() -> None:
    assert generate_name_track("signal.bedGraph.gz") == (
        "signal.track.bedGraph.gz"
    )


def test_sum_counts_bdg_expands_run_lengths_over_bins(
    tmp_path: Path,
) -> None:
    """
    An interval spanning 'k' bins contributes 'k' times its value, not once.
    """

    path = tmp_path / "runs.bdg"
    path.write_text(
        "chrI\t0\t10\t2\nchrI\t10\t40\t3\n",
        encoding="utf-8",
    )

    # One bin at 2, then three bins at 3, on a grid the track itself implies.
    result = sum_counts_bdg(str(path))

    assert result.siz_bin == 10
    assert result.total == 2.0 + 3.0 * 3


def test_sum_counts_bdg_counts_a_partial_terminal_bin_as_one_bin(
    tmp_path: Path,
) -> None:
    """
    A partial bin still holds a real count, so it rounds up rather than down.
    """

    path = tmp_path / "partial.bdg"
    path.write_text(
        "chrI\t0\t10\t1\nchrI\t10\t14\t5\n",
        encoding="utf-8",
    )

    result = sum_counts_bdg(str(path), 10)

    # Scoring the 4 bp remainder as 0.4 of a bin would give 3.0 instead.
    assert result.total == 1.0 + 5.0


def test_sum_counts_bdg_infers_the_bin_width_from_the_track(
    tmp_path: Path,
) -> None:
    path = tmp_path / "fifty.bdg"
    path.write_text(
        "chrI\t0\t50\t1\nchrI\t50\t150\t2\n",
        encoding="utf-8",
    )

    inferred = sum_counts_bdg(str(path))
    supplied = sum_counts_bdg(str(path), 50)

    assert inferred.siz_bin == 50
    assert inferred.total == supplied.total


def test_sum_counts_bdg_rejects_a_width_the_track_contradicts(
    tmp_path: Path,
) -> None:
    """
    A wrong width silently rescales every derived pseudocount, so refuse it.
    """

    path = tmp_path / "ten.bdg"
    path.write_text(
        "chrI\t0\t10\t1\nchrI\t10\t20\t2\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="imply a bin width of 10"):
        sum_counts_bdg(str(path), 20)


def test_sum_counts_bdg_rejects_a_nonpositive_width(tmp_path: Path) -> None:
    path = tmp_path / "any.bdg"
    path.write_text("chrI\t0\t10\t1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="must be a positive integer"):
        sum_counts_bdg(str(path), 0)


def test_sum_counts_bdg_rejects_a_track_with_no_usable_interval(
    tmp_path: Path,
) -> None:
    path = tmp_path / "none.bdg"
    path.write_text("track name=empty\n", encoding="utf-8")

    with pytest.raises(ValueError, match="no usable bedGraph interval"):
        sum_counts_bdg(str(path))


def test_sum_counts_bdg_rounds_up_once_per_chromosome(
    tmp_path: Path,
) -> None:
    """
    Only a chromosome's final interval can be partial, so only it rounds up.
    """

    path = tmp_path / "two_chrom.bdg"
    path.write_text(
        "chrI\t0\t10\t1\nchrI\t10\t15\t1\nchrII\t0\t10\t1\nchrII\t10\t15\t1\n",
        encoding="utf-8",
    )

    result = sum_counts_bdg(str(path), 10)

    # Two full bins and two rounded-up partials, one partial per chromosome.
    assert result.total == 4.0


def test_sum_counts_bdg_rejects_a_noninteger_width(tmp_path: Path) -> None:
    """
    A float width corrupts the floor and remainder arithmetic downstream.
    """

    path = tmp_path / "any.bdg"
    path.write_text("chrI\t0\t10\t1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="must be a positive integer"):
        sum_counts_bdg(str(path), 2.5)


def test_sum_counts_bdg_survives_a_clamped_terminal_bin(
    tmp_path: Path,
) -> None:
    """
    Track ends are clamped to chromosome size, so terminal runs are short.

    A 75 bp chromosome on a 10 bp grid ends in a 5 bp interval. Folding that
    width into the divisor would infer a 5 bp grid and inflate the library
    size, so terminal widths stay out of it.
    """

    path = tmp_path / "clamped.bdg"
    path.write_text(
        "chrI\t0\t10\t2\nchrI\t10\t70\t3\nchrI\t70\t75\t4\n",
        encoding="utf-8",
    )

    result = sum_counts_bdg(str(path))

    assert result.siz_bin == 10
    assert result.total == 2.0 + 3.0 * 6 + 4.0


def test_sum_counts_bdg_infers_a_multiple_when_runs_share_a_factor(
    tmp_path: Path,
) -> None:
    """
    Pin the documented limit of inference so it is not mistaken for exact.

    Every run here spans an even number of 10 bp bins, so the divisor lands
    on 20 and the library size halves. Real tracks carry runs of differing
    lengths and do not do this; supplying 'siz_bin' is the stated remedy.
    """

    path = tmp_path / "even_runs.bdg"
    path.write_text(
        "chrI\t0\t20\t1\nchrI\t20\t40\t2\nchrI\t40\t60\t3\n",
        encoding="utf-8",
    )

    inferred = sum_counts_bdg(str(path))
    supplied = sum_counts_bdg(str(path), 10)

    assert inferred.siz_bin == 20
    assert inferred.total == 6.0
    assert supplied.total == 12.0
