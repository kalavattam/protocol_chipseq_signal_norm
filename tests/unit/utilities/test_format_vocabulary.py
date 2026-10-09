#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_format_vocabulary.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Tests for the bedGraph-only format vocabulary.

Notes
-----
Format names and file suffixes are matched in any letter case and resolve to
'bam', 'cram', 'bed', or 'bedGraph'; the short bedGraph spellings 'bdg' and
'bg' are refused, as names and as suffixes, for input and output alike.
"""

import warnings
from collections.abc import Callable
from pathlib import Path

import pytest

from protocol_chipseq_signal_norm.cli import (
    compute_input_floor,
    compute_pseudo,
    compute_pseudo_deeptools,
    compute_signal,
    compute_signal_ratio,
    merge_bins_bdg,
    sum_bdg,
)
from protocol_chipseq_signal_norm.utilities.utils_check import (
    canonicalize_format,
    check_bedgraph_path,
    format_from_path,
    validate_output_path,
)
from protocol_chipseq_signal_norm.utilities.utils_cli import (
    warn_as_note,
)

ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tests" / "fixtures"
BAM = str(FIXTURES / "compute_signal" / "bam" / "se" / "tiny_se.bam")
PAIR_A = str(FIXTURES / "compute_pseudo" / "bedgraph" / "pair_A.bedGraph")
RATIO_A = str(FIXTURES / "compute_signal" / "bedgraph" / "ratio_A.bedGraph")
RATIO_B = str(FIXTURES / "compute_signal" / "bedgraph" / "ratio_B.bedGraph")


@pytest.mark.parametrize(
    ("name", "fmt"),
    (
        ("bedGraph", "bedGraph"),
        ("bedgraph", "bedGraph"),
        ("bEdGrApH", "bedGraph"),
        ("BED", "bed"),
        ("CrAm", "cram"),
        ("bam", "bam"),
    ),
)
def test_format_names_resolve_in_any_letter_case(name: str, fmt: str) -> None:
    assert canonicalize_format(name) == fmt


@pytest.mark.parametrize("name", ("bdg", "BDG", "bg", "Bg"))
def test_short_bedgraph_names_are_refused(name: str) -> None:
    with pytest.raises(ValueError, match="is not accepted; use 'bedGraph'"):
        canonicalize_format(name)


def test_unknown_format_names_list_the_choices() -> None:
    with pytest.raises(
        ValueError,
        match="choose 'bam', 'cram', 'bed', or 'bedGraph'",
    ):
        canonicalize_format("txt")


@pytest.mark.parametrize(
    ("path", "fmt"),
    (
        ("x.bedGraph", "bedGraph"),
        ("x.BEDGRAPH.gz", "bedGraph"),
        ("x.bedGraph.GZ", "bedGraph"),
        ("x.CrAm", "cram"),
        ("x.bed.gz", "bed"),
        ("x.bam.gz", None),
        ("x.txt", None),
    ),
)
def test_suffixes_resolve_in_any_letter_case(
    path: str,
    fmt: str | None,
) -> None:
    assert format_from_path(path) == fmt


@pytest.mark.parametrize(
    ("path", "want"),
    (
        ("x.bdg", ".bedGraph"),
        ("x.BDG.gz", ".bedGraph.gz"),
        ("x.bg", ".bedGraph"),
        ("x.bg.gz", ".bedGraph.gz"),
    ),
)
def test_short_bedgraph_suffixes_are_refused(path: str, want: str) -> None:
    for check in (format_from_path, check_bedgraph_path, validate_output_path):
        with pytest.raises(ValueError, match=f"must end in '{want}'"):
            check(path)


@pytest.mark.parametrize("path", ("-", "/dev/fd/63", "x.txt", "noext"))
def test_other_input_names_are_read_as_bedgraph(path: str) -> None:
    check_bedgraph_path(path)


def test_canonical_output_suffix_draws_no_warning() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")

        assert validate_output_path("x.bedGraph.gz") == (
            "x.bedGraph.gz",
            "bedGraph",
            True,
        )


@pytest.mark.parametrize(
    ("path", "canonical"),
    (
        ("x.BEDGRAPH.gz", ".bedGraph.gz"),
        ("x.bedgraph", ".bedGraph"),
        ("x.bedGraph.GZ", ".bedGraph.gz"),
        ("x.BED", ".bed"),
    ),
)
def test_case_variant_output_is_written_as_given_with_a_warning(
    path: str,
    canonical: str,
) -> None:
    with pytest.warns(UserWarning) as record:
        output_path, _, _ = validate_output_path(path)

    assert output_path == path
    assert [str(item.message) for item in record] == [
        f"writing '{path}' as given; the canonical suffix is '{canonical}'.",
    ]


def test_warn_as_note_prints_each_warning(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with warn_as_note():
        warnings.warn("first", stacklevel=1)
        warnings.warn("second", stacklevel=1)

    assert capsys.readouterr().err == "Note: first\nNote: second\n"


def _run(main: Callable[[list[str]], int], argv: list[str]) -> str:
    """
    Run a CLI expected to refuse, and return its message.
    """

    with pytest.raises(SystemExit) as error:
        main(argv)

    return str(error.value.code)


@pytest.mark.parametrize(
    ("main", "argv"),
    (
        pytest.param(
            compute_signal_ratio.main,
            [
                "--fil_A",
                "a.bdg",
                "--fil_B",
                RATIO_B,
                "--fil_out",
                "r.bedGraph",
            ],
            id="compute_signal_ratio fil_A",
        ),
        pytest.param(
            compute_signal_ratio.main,
            [
                "--fil_A",
                RATIO_A,
                "--fil_B",
                "a.bdg",
                "--fil_out",
                "r.bedGraph",
            ],
            id="compute_signal_ratio fil_B",
        ),
        pytest.param(
            compute_pseudo.main,
            ["--fil_A", "a.bdg"],
            id="compute_pseudo fil_A",
        ),
        pytest.param(
            compute_pseudo_deeptools.main,
            ["--fil_A", "a.bdg", "--typ_sig", "CPM"],
            id="compute_pseudo_deeptools fil_A",
        ),
        pytest.param(
            merge_bins_bdg.main,
            ["--fil_in", "a.bdg", "--fil_out", "m.bedGraph"],
            id="merge_bins_bdg fil_in",
        ),
        pytest.param(sum_bdg.main, ["a.bdg"], id="sum_bdg path"),
    ),
)
def test_clis_refuse_short_input_suffixes(
    main: Callable[[list[str]], int],
    argv: list[str],
) -> None:
    assert _run(main, argv) == (
        "'a.bdg': '.bdg' is not accepted; the file name must end in "
        "'.bedGraph'."
    )


def test_compute_input_floor_refuses_a_short_input_suffix(
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = compute_input_floor.main(
        ["--mode", "dist", "--fil_in", "a.bdg"],
    )

    assert status == 1
    assert "'a.bdg': '.bdg' is not accepted" in capsys.readouterr().err


@pytest.mark.parametrize(
    ("main", "argv"),
    (
        pytest.param(
            compute_signal.main,
            ["--threads", "1", "--fil_in", BAM, "--fil_out", "x.bdg.gz"],
            id="compute_signal",
        ),
        pytest.param(
            compute_signal_ratio.main,
            ["--fil_A", RATIO_A, "--fil_B", RATIO_B, "--fil_out", "x.bdg.gz"],
            id="compute_signal_ratio",
        ),
        pytest.param(
            merge_bins_bdg.main,
            ["--fil_in", PAIR_A, "--fil_out", "x.bdg.gz"],
            id="merge_bins_bdg",
        ),
    ),
)
def test_clis_refuse_short_output_suffixes(
    main: Callable[[list[str]], int],
    argv: list[str],
) -> None:
    assert _run(main, argv) == (
        "'x.bdg.gz': '.bdg.gz' is not accepted; the file name must end in "
        "'.bedGraph.gz'."
    )


@pytest.mark.parametrize(
    ("main", "argv"),
    (
        pytest.param(
            compute_signal.main,
            ["--threads", "1", "--fil_in", BAM],
            id="compute_signal",
        ),
        pytest.param(
            compute_signal_ratio.main,
            ["--fil_A", RATIO_A, "--fil_B", RATIO_B],
            id="compute_signal_ratio",
        ),
        pytest.param(
            merge_bins_bdg.main,
            ["--fil_in", PAIR_A],
            id="merge_bins_bdg",
        ),
    ),
)
def test_clis_write_a_case_variant_output_with_a_note(
    main: Callable[[list[str]], int],
    argv: list[str],
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    fil_out = tmp_path / "x.BEDGRAPH"

    assert main([*argv, "--fil_out", str(fil_out)]) in (0, None)
    assert fil_out.is_file()
    assert (
        f"Note: writing '{fil_out}' as given; the canonical suffix is "
        "'.bedGraph'."
    ) in capsys.readouterr().err
