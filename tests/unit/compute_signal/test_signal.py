#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_signal.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


import shutil
from pathlib import Path

import numpy as np
import pytest

from protocol_chipseq_signal_norm.cli.compute_signal import (
    calc_sig_chrom,
    calc_sig_chrom_array,
    calc_sig_task,
    get_alignment_chrom_sizes,
    iter_alignment_fragments,
    iter_indexed_fragments,
    main,
    parse_args,
    resolve_chrom_sizes,
)

ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tests" / "fixtures" / "compute_signal"


@pytest.mark.parametrize(
    "alias",
    (
        "-ck",
        "--chunk_size",
        "--chunk-size",
        "--chnk_size",
        "--chnk-size",
    ),
)
def test_chunk_size_aliases_are_accepted(alias: str) -> None:
    args = parse_args(
        [
            "--fil_in",
            "input.bam",
            "--fil_out",
            "output.bdg",
            alias,
            "17",
        ],
    )

    assert args.chunk_size == 17


@pytest.mark.parametrize(
    "alias",
    (
        "--chunk_s",
        "--chnk_s",
        "--chunck_size",
        "--chnk__size",
        "--chunk_sizes",
    ),
)
def test_chunk_size_near_misses_are_rejected(alias: str) -> None:
    with pytest.raises(SystemExit) as raised:
        parse_args(
            [
                "--fil_in",
                "input.bam",
                "--fil_out",
                "output.bdg",
                alias,
                "17",
            ],
        )

    assert raised.value.code == 2


def test_chunk_size_hidden_aliases_are_absent_from_help(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as raised:
        parse_args(["--help"])

    assert raised.value.code == 0

    rendered = capsys.readouterr()
    help_text = rendered.out + rendered.err

    assert "-ck" in help_text
    assert "--chunk_size" in help_text
    assert "--chunk-size" not in help_text
    assert "--chnk_size" not in help_text
    assert "--chnk-size" not in help_text


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def collect_indexed_windows(
    fil_in: Path,
    chrom_sizes: dict[str, int],
    user_fragment_length: int | None,
    window_size: int,
) -> list[tuple[str, int, int, int]]:
    fragments = []

    for chrom, chrom_size in chrom_sizes.items():
        for start in range(0, chrom_size, window_size):
            fragments.extend(
                iter_indexed_fragments(
                    str(fil_in),
                    chrom=chrom,
                    chrom_sizes=chrom_sizes,
                    user_fragment_length=user_fragment_length,
                    start=start,
                    end=min(start + window_size, chrom_size),
                ),
            )

    return fragments


def test_calc_sig_chrom_accumulates_fragment_overlap_by_bin() -> None:
    observed = calc_sig_chrom(
        "chrI",
        [(0, 10, 10), (5, 15, 10)],
        fragment_count=2,
        siz_bin=10,
        is_len=False,
        is_norm=False,
    )

    assert observed[("chrI", 0)] == 15.0
    assert observed[("chrI", 10)] == 5.0


def test_calc_sig_chrom_length_and_count_normalization() -> None:
    observed = calc_sig_chrom(
        "chrI",
        [(0, 10, 10), (5, 15, 10)],
        fragment_count=2,
        siz_bin=10,
        is_len=True,
        is_norm=True,
        scl_fct=4.0,
    )

    assert observed[("chrI", 0)] == pytest.approx(3.0)
    assert observed[("chrI", 10)] == pytest.approx(1.0)


def test_calc_sig_task_dispatches_to_calc_sig_chrom() -> None:
    task = ("chrI", [(0, 5, 5)], 1, 10, False, False, None)

    assert calc_sig_task(task)[("chrI", 0)] == 5.0


def test_calc_sig_chrom_array_matches_python_kernel() -> None:
    fragments = [(0, 10, 10), (5, 25, 20), (12, 15, 3), (79, 80, 1)]
    expected = calc_sig_chrom(
        "chrI",
        fragments,
        fragment_count=len(fragments),
        siz_bin=10,
        is_len=True,
        is_norm=True,
        scl_fct=4.0,
    )

    observed = calc_sig_chrom_array(
        "chrI",
        np.asarray([frag[0] for frag in fragments]),
        np.asarray([frag[1] for frag in fragments]),
        np.asarray([frag[2] for frag in fragments], dtype=float),
        chrom_size=80,
        fragment_count=len(fragments),
        siz_bin=10,
        is_len=True,
        is_norm=True,
        scl_fct=4.0,
    )

    assert observed.keys() == expected.keys()

    for key in expected:
        assert observed[key] == pytest.approx(expected[key])


def test_calc_sig_chrom_array_emits_only_touched_bins() -> None:
    fragments = [
        (0, 26040, 3519),
        (26080, 40000, 3519),
    ]
    expected = calc_sig_chrom(
        "chrI",
        fragments,
        fragment_count=len(fragments),
        siz_bin=10,
        is_len=True,
        is_norm=True,
    )

    observed = calc_sig_chrom_array(
        "chrI",
        np.asarray([frag[0] for frag in fragments]),
        np.asarray([frag[1] for frag in fragments]),
        np.asarray([frag[2] for frag in fragments], dtype=float),
        chrom_size=50000,
        fragment_count=len(fragments),
        siz_bin=10,
        is_len=True,
        is_norm=True,
    )

    assert observed.keys() == expected.keys()
    assert ("chrI", 26040) not in observed
    assert ("chrI", 26050) not in observed
    assert ("chrI", 26060) not in observed
    assert ("chrI", 26070) not in observed

    for key in expected:
        assert observed[key] == pytest.approx(expected[key])


def test_calc_sig_chrom_rejects_invalid_bin_size() -> None:
    with pytest.raises(ValueError, match="siz_bin"):
        calc_sig_chrom("chrI", [], 0, 0, False, False)


def test_resolve_chrom_sizes_rejects_conflicts(tmp_path: Path) -> None:
    sizes = tmp_path / "chr.sizes"
    sizes.write_text("chrI\t81\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Conflicting chromosome sizes"):
        resolve_chrom_sizes({"chrI": 80}, str(sizes))


def test_compute_signal_engines_match_and_clamp_bedgraph_end(
    tmp_path: Path,
) -> None:
    fil_in = FIXTURES / "bam" / "se" / "tiny_se.bam"
    out_chrom = tmp_path / "chrom.bdg"
    out_window = tmp_path / "window.bdg"

    chrom_status = main(
        [
            "--fil_in",
            str(fil_in),
            "--fil_out",
            str(out_chrom),
            "--method",
            "unadj",
            "--siz_bin",
            "10",
        ],
    )
    window_status = main(
        [
            "--fil_in",
            str(fil_in),
            "--fil_out",
            str(out_window),
            "--method",
            "unadj",
            "--siz_bin",
            "10",
            "--engine",
            "window",
            "--indexed_window_size",
            "15",
        ],
    )

    assert chrom_status == 0
    assert window_status == 0
    assert out_chrom.read_text(encoding="utf-8") == out_window.read_text(
        encoding="utf-8",
    )

    for line in out_chrom.read_text(encoding="utf-8").splitlines():
        chrom, _start, end, _value = line.split("\t")

        assert chrom == "I"
        assert int(end) <= 80


@pytest.mark.parametrize(
    "subdir,user_fragment_length",
    [
        ("bam/se/tiny_se.bam", 30),
        ("bam/pe/tiny_pe.bam", None),
    ],
)
def test_indexed_chrom_and_window_fetch_match_serial_fragments(
    subdir: str,
    user_fragment_length: int | None,
) -> None:
    fil_in = FIXTURES / subdir
    chrom_sizes = get_alignment_chrom_sizes(str(fil_in))
    serial = sorted(
        iter_alignment_fragments(
            str(fil_in),
            chrom_sizes=chrom_sizes,
            user_fragment_length=user_fragment_length,
        ),
    )
    indexed_chrom = []

    for chrom in chrom_sizes:
        indexed_chrom.extend(
            iter_indexed_fragments(
                str(fil_in),
                chrom=chrom,
                chrom_sizes=chrom_sizes,
                user_fragment_length=user_fragment_length,
            ),
        )

    indexed_windows = collect_indexed_windows(
        fil_in,
        chrom_sizes,
        user_fragment_length=user_fragment_length,
        window_size=15,
    )

    assert sorted(indexed_chrom) == serial
    assert sorted(indexed_windows) == serial


@pytest.mark.parametrize(
    "subdir,ref_args,user_fragment_length",
    [
        ("bam/se/tiny_se.bam", [], "30"),
        (
            "cram/se/tiny_se.cram",
            ["--ref_fa", str(FIXTURES / "reference" / "tiny.fa")],
            "30",
        ),
        ("bam/pe/tiny_pe.bam", [], None),
        (
            "cram/pe/tiny_pe.cram",
            ["--ref_fa", str(FIXTURES / "reference" / "tiny.fa")],
            None,
        ),
    ],
)
def test_compute_signal_public_engines_match_default(
    tmp_path: Path,
    subdir: str,
    ref_args: list[str],
    user_fragment_length: int | None,
) -> None:
    fil_in = FIXTURES / subdir
    out_default = tmp_path / f"{fil_in.stem}.default.bdg"
    out_chrom = tmp_path / f"{fil_in.stem}.chrom.bdg"
    out_window = tmp_path / f"{fil_in.stem}.window.bdg"
    base_args = [
        "--fil_in",
        str(fil_in),
        "--method",
        "norm",
        "--scl_fct",
        "2.5",
        "--siz_bin",
        "10",
        "--dp",
        "12",
        *ref_args,
    ]

    if user_fragment_length is not None:
        base_args.extend(["--usr_frg", user_fragment_length])

    default_status = main([*base_args, "--fil_out", str(out_default)])
    chrom_status = main(
        [
            *base_args,
            "--fil_out",
            str(out_chrom),
            "--engine",
            "chrom",
        ],
    )
    window_status = main(
        [
            *base_args,
            "--fil_out",
            str(out_window),
            "--engine",
            "window",
            "--indexed_window_size",
            "15",
        ],
    )

    assert default_status == 0
    assert chrom_status == 0
    assert window_status == 0
    assert read_text(out_chrom) == read_text(out_default)
    assert read_text(out_window) == read_text(out_default)


@pytest.mark.parametrize(
    "result_format",
    [
        "direct_sparse_local_np",
        "direct_sparse_bincount_np",
        "direct_sparse_local_bincount_np",
    ],
)
def test_compute_signal_hidden_sparse_variants_match_public_baseline(
    tmp_path: Path,
    result_format: str,
) -> None:
    fil_in = FIXTURES / "bam" / "se" / "tiny_se.bam"
    out_default = tmp_path / "default.bdg"
    out_variant = tmp_path / f"{result_format}.bdg"
    base_args = [
        "--fil_in",
        str(fil_in),
        "--method",
        "norm",
        "--scl_fct",
        "2.5",
        "--siz_bin",
        "10",
        "--dp",
        "12",
        "--usr_frg",
        "30",
    ]

    default_status = main([*base_args, "--fil_out", str(out_default)])
    variant_status = main(
        [
            *base_args,
            "--fil_out",
            str(out_variant),
            "--engine",
            "window",
            "--indexed_window_size",
            "15",
            "--prototype_result_format",
            result_format,
        ],
    )

    assert default_status == 0
    assert variant_status == 0
    assert read_text(out_variant) == read_text(out_default)


@pytest.mark.parametrize(
    "bed_strategy",
    [
        "indexed_chrom",
        "indexed_window",
    ],
)
def test_compute_signal_hidden_indexed_bed_strategies_match_serial(
    tmp_path: Path,
    bed_strategy: str,
) -> None:
    fil_in = FIXTURES / "bam" / "pe" / "tiny_pe.bam"
    out_serial = tmp_path / "serial.bed"
    out_indexed = tmp_path / f"{bed_strategy}.bed"

    serial_status = main(
        [
            "--fil_in",
            str(fil_in),
            "--fil_out",
            str(out_serial),
            "--prototype_bed_strategy",
            "serial",
        ],
    )
    indexed_status = main(
        [
            "--fil_in",
            str(fil_in),
            "--fil_out",
            str(out_indexed),
            "--threads",
            "1",
            "--prototype_bed_strategy",
            bed_strategy,
            "--indexed_window_size",
            "15",
        ],
    )

    assert serial_status == 0
    assert indexed_status == 0
    assert read_text(out_indexed) == read_text(out_serial)


def test_compute_signal_auto_bed_falls_back_without_index(
    tmp_path: Path,
) -> None:
    fil_in = tmp_path / "tiny_no_index.bam"
    shutil.copyfile(FIXTURES / "bam" / "se" / "tiny_se.bam", fil_in)
    out_bed = tmp_path / "coord.bed"

    auto_status = main(
        [
            "--fil_in",
            str(fil_in),
            "--fil_out",
            str(out_bed),
            "--threads",
            "2",
            "--prototype_bed_strategy",
            "auto",
        ],
    )

    assert auto_status == 0
    assert out_bed.read_text(encoding="utf-8")


def test_compute_signal_public_engines_require_index(tmp_path: Path) -> None:
    fil_in = tmp_path / "tiny_no_index.bam"
    shutil.copyfile(FIXTURES / "bam" / "se" / "tiny_se.bam", fil_in)

    with pytest.raises(SystemExit, match="alignment index"):
        main(
            [
                "--fil_in",
                str(fil_in),
                "--fil_out",
                str(tmp_path / "out.bdg"),
                "--method",
                "unadj",
                "--engine",
                "chrom",
            ],
        )


def test_compute_signal_public_cram_requires_reference_fasta(
    tmp_path: Path,
) -> None:
    fil_in = FIXTURES / "cram" / "pe" / "tiny_pe.cram"

    with pytest.raises(SystemExit, match="CRAM signal engines require"):
        main(
            [
                "--fil_in",
                str(fil_in),
                "--fil_out",
                str(tmp_path / "out.bdg"),
                "--method",
                "unadj",
                "--engine",
                "chrom",
            ],
        )


def test_compute_signal_public_engines_and_bed_output(tmp_path: Path) -> None:
    fil_in = FIXTURES / "bam" / "pe" / "tiny_pe.bam"

    for engine in ("chunk", "indexed_chrom"):
        with pytest.raises(SystemExit):
            main(
                [
                    "--fil_in",
                    str(fil_in),
                    "--fil_out",
                    str(tmp_path / f"{engine}.bdg"),
                    "--method",
                    "unadj",
                    "--engine",
                    engine,
                ],
            )

    out_bed = tmp_path / "coord.bed"

    chrom_status = main(
        [
            "--fil_in",
            str(fil_in),
            "--fil_out",
            str(out_bed),
            "--engine",
            "chrom",
        ],
    )

    assert chrom_status == 0
    assert out_bed.read_text(encoding="utf-8")


def test_compute_signal_rejects_dash_io(tmp_path: Path) -> None:
    fil_in = FIXTURES / "bam" / "se" / "tiny_se.bam"

    with pytest.raises(SystemExit, match="--fil_in -"):
        main(["--fil_in", "-", "--fil_out", str(tmp_path / "x.bdg")])

    with pytest.raises(SystemExit, match="--fil_out -"):
        main(["--fil_in", str(fil_in), "--fil_out", "-"])
