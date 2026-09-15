#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_report_counts.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5, Fable 5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


from pathlib import Path

import pytest

from protocol_chipseq_signal_norm.cli.compute_signal import (
    count_frgs_and_bins,
    get_siz_chr,
    iter_aln_frg,
    main,
    parse_args,
)

ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tests" / "fixtures" / "compute_signal"

# One fragment at a time, so each expected span is hand-checkable: a half-open
# fragment '[start, end)' touches bins 'start // siz_bin' through
# '(end - 1) // siz_bin' inclusive, and a partial terminal bin counts as one
# bin like any other.
SINGLE_FRAGMENT_CASES = [
    pytest.param(0, 10, 10, 1, id="bin-aligned-single-bin"),
    pytest.param(3, 9, 10, 1, id="interior-single-bin"),
    pytest.param(5, 15, 10, 2, id="boundary-straddling"),
    pytest.param(9, 11, 10, 2, id="two-bins-by-one-base-each"),
    pytest.param(10, 20, 10, 1, id="second-bin-aligned"),
    pytest.param(0, 21, 10, 3, id="three-bins"),
    pytest.param(60, 80, 30, 1, id="partial-terminal-bin"),
    pytest.param(30, 80, 30, 2, id="full-plus-partial-terminal"),
]


def read_count(path: Path) -> int:
    """
    Read one report file, asserting it holds exactly one integer line.
    """

    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()

    assert text.endswith("\n")
    assert len(lines) == 1

    return int(lines[0])


def run_report(
    tmp_path: Path,
    fil_in: Path,
    extra_args: list[str],
) -> tuple[int, int]:
    """
    Run report-only mode with both report flags and return '(N, L)'.
    """

    path_n = tmp_path / "report_n.txt"
    path_l = tmp_path / "report_l.txt"

    status = main(
        [
            "--fil_in",
            str(fil_in),
            "--report_n_frg",
            str(path_n),
            "--report_n_bin",
            str(path_l),
            *extra_args,
        ],
    )

    assert status == 0

    return read_count(path_n), read_count(path_l)


@pytest.mark.parametrize(
    ("start", "end", "siz_bin", "expected_bins"),
    SINGLE_FRAGMENT_CASES,
)
def test_count_frgs_and_bins_single_fragment_span(
    start: int,
    end: int,
    siz_bin: int,
    expected_bins: int,
) -> None:
    observed = count_frgs_and_bins(
        iter([("I", start, end, end - start)]),
        siz_bin,
    )

    assert observed == (1, expected_bins)


@pytest.mark.parametrize("offset", range(10))
def test_length_rule_fragment_span_is_offset_independent(offset: int) -> None:
    # 'F = (k - 1) * b + 1' is the one length whose bin span cannot depend on
    # placement: with 'F - 1' an exact multiple of the bin size, the last base
    # sits exactly 'k - 1' bins after the first at every offset.
    length = 2 * 10 + 1

    observed = count_frgs_and_bins(
        iter([("I", offset, offset + length, length)]),
        10,
    )

    assert observed == (1, 3)


def test_count_frgs_and_bins_accumulates_across_chromosomes() -> None:
    fragments = [
        ("I", 0, 25, 25),
        ("II", 95, 105, 10),
        ("II", 5, 8, 3),
    ]

    observed = count_frgs_and_bins(iter(fragments), 10)

    assert observed == (3, 6)


def test_count_frgs_and_bins_skips_nonpositive_spans() -> None:
    fragments = [
        ("I", 10, 10, 0),
        ("I", 20, 15, -5),
        ("I", 0, 5, 5),
    ]

    observed = count_frgs_and_bins(iter(fragments), 10)

    assert observed == (1, 1)


def test_count_frgs_and_bins_empty_input_is_zero_zero() -> None:
    assert count_frgs_and_bins(iter([]), 10) == (0, 0)


def test_count_frgs_and_bins_ignores_the_length_field() -> None:
    # 'L' derives from the clamped coordinates alone; the carried
    # fragment-length field weights signal, not spans.
    short_length = count_frgs_and_bins(iter([("I", 0, 25, 1)]), 10)
    long_length = count_frgs_and_bins(iter([("I", 0, 25, 400)]), 10)

    assert short_length == long_length == (1, 3)


def test_report_flag_aliases_leave_fil_out_optional() -> None:
    args = parse_args(
        [
            "--fil_in",
            "input.bam",
            "-rnf",
            "n.txt",
            "-rnb",
            "l.txt",
        ],
    )

    assert args.fil_out is None
    assert args.report_n_frg == "n.txt"
    assert args.report_n_bin == "l.txt"


def test_report_counts_pe_tlen_fragments(tmp_path: Path) -> None:
    fil_in = FIXTURES / "bam" / "pe" / "tiny_pe.bam"

    # The two leftmost proper-pair anchors give TLEN-derived fragments
    # '[10, 40)' and '[40, 60)', spanning three and two ten-base bins.
    observed = run_report(tmp_path, fil_in, ["--siz_bin", "10"])

    assert observed == (2, 5)


def test_report_counts_se_read_length_fragments(tmp_path: Path) -> None:
    fil_in = FIXTURES / "bam" / "se" / "tiny_se.bam"

    # Without '--usr_frg', both ten-base reads span exactly their aligned bin:
    # '[0, 10)' forward and '[20, 30)' reverse.
    observed = run_report(tmp_path, fil_in, ["--siz_bin", "10"])

    assert observed == (2, 2)


def test_report_counts_recover_fragment_length_via_bins(
    tmp_path: Path,
) -> None:
    fil_in = FIXTURES / "bam" / "se" / "tiny_se.bam"
    siz_bin = 10
    usr_frg = 2 * siz_bin + 1

    fragment_count, bin_count = run_report(
        tmp_path,
        fil_in,
        ["--siz_bin", str(siz_bin), "--usr_frg", str(usr_frg)],
    )
    spanned = bin_count // fragment_count

    assert (fragment_count, bin_count) == (2, 6)
    assert bin_count % fragment_count == 0
    assert (spanned - 1) * siz_bin + 1 == usr_frg


def test_report_counts_clamp_to_partial_terminal_bin(tmp_path: Path) -> None:
    fil_in = FIXTURES / "bam" / "se" / "tiny_se.bam"

    # Extension with '--usr_frg 100' overruns the 80-base chromosome: the
    # forward fragment clamps to '[0, 80)' and spans three 30-base bins, with
    # the partial terminal bin '[60, 80)' counting as one; the reverse fragment
    # clamps to '[0, 30)' and spans one.
    observed = run_report(
        tmp_path,
        fil_in,
        ["--siz_bin", "30", "--usr_frg", "100"],
    )

    assert observed == (2, 4)


def test_cli_counts_match_direct_iterator_counts(tmp_path: Path) -> None:
    fil_in = FIXTURES / "bam" / "pe" / "tiny_pe.bam"
    siz_chr = get_siz_chr(str(fil_in))
    expected = count_frgs_and_bins(
        iter_aln_frg(str(fil_in), siz_chr=siz_chr),
        10,
    )

    observed = run_report(tmp_path, fil_in, ["--siz_bin", "10"])

    assert observed == expected


@pytest.mark.parametrize(
    ("with_n", "with_l"),
    [
        pytest.param(True, False, id="n-only"),
        pytest.param(False, True, id="l-only"),
        pytest.param(True, True, id="both"),
    ],
)
def test_report_only_mode_counts_without_writing_a_track(
    tmp_path: Path,
    with_n: bool,
    with_l: bool,
) -> None:
    fil_in = FIXTURES / "bam" / "se" / "tiny_se.bam"
    path_n = tmp_path / "n.txt"
    path_l = tmp_path / "l.txt"
    arguments = ["--fil_in", str(fil_in), "--siz_bin", "10"]

    if with_n:
        arguments.extend(["--report_n_frg", str(path_n)])

    if with_l:
        arguments.extend(["--report_n_bin", str(path_l)])

    status = main(arguments)
    created = sorted(entry.name for entry in tmp_path.iterdir())
    expected = sorted(
        name
        for name, wanted in (("n.txt", with_n), ("l.txt", with_l))
        if wanted
    )

    assert status == 0
    assert created == expected

    if with_n:
        assert read_count(path_n) == 2

    if with_l:
        assert read_count(path_l) == 2


def test_missing_fil_out_without_report_flags_is_rejected() -> None:
    fil_in = FIXTURES / "bam" / "se" / "tiny_se.bam"

    with pytest.raises(SystemExit, match="'--fil_out' is required"):
        main(["--fil_in", str(fil_in)])


def test_reporting_leaves_the_written_track_byte_identical(
    tmp_path: Path,
) -> None:
    fil_in = FIXTURES / "bam" / "pe" / "tiny_pe.bam"
    out_plain = tmp_path / "plain.bdg"
    out_reported = tmp_path / "reported.bdg"
    base_args = [
        "--fil_in",
        str(fil_in),
        "--method",
        "norm",
        "--siz_bin",
        "10",
    ]

    plain_status = main([*base_args, "--fil_out", str(out_plain)])
    reported_status = main(
        [
            *base_args,
            "--fil_out",
            str(out_reported),
            "--report_n_frg",
            str(tmp_path / "n.txt"),
            "--report_n_bin",
            str(tmp_path / "l.txt"),
        ],
    )

    assert plain_status == 0
    assert reported_status == 0
    assert out_reported.read_bytes() == out_plain.read_bytes()


@pytest.mark.parametrize("engine", ["chrom", "window"])
@pytest.mark.parametrize("method", ["unadj", "frag", "norm"])
def test_emitted_track_carries_no_zero_valued_row(
    tmp_path: Path,
    engine: str,
    method: str,
) -> None:
    # 'ZR1': 'compute_signal' emits no zero-valued row. Lab-facing advice (that
    # '--skip_00' is redundant on our own tracks) depends on this, so the
    # property is pinned here against changes to the emission predicate.
    fil_in = FIXTURES / "bam" / "se" / "tiny_se.bam"
    fil_out = tmp_path / "track.bdg"

    status = main(
        [
            "--fil_in",
            str(fil_in),
            "--fil_out",
            str(fil_out),
            "--method",
            method,
            "--siz_bin",
            "10",
            "--usr_frg",
            "30",
            "--engine",
            engine,
        ],
    )
    rows = fil_out.read_text(encoding="utf-8").splitlines()

    assert status == 0
    assert rows

    for row in rows:
        value = row.split("\t")[3]

        assert float(value) != 0.0
        assert value != "-0"


@pytest.mark.parametrize(
    ("nam_out", "base"),
    [
        pytest.param("track.bedGraph.gz", "track", id="bedgraph-gz"),
        pytest.param("track.bedGraph", "track", id="bedgraph-plain"),
        pytest.param("track.bed", "track", id="bed"),
        pytest.param("a.b.track.bedGraph.gz", "a.b.track", id="dotted-stem"),
    ],
)
def test_bare_report_flags_derive_paths_from_fil_out(
    tmp_path: Path,
    nam_out: str,
    base: str,
) -> None:
    # The derivation strips one trailing '.gz' and then one extension, which is
    # the rule the execute wrapper applies in Bash. Pinning both spellings and
    # a dotted stem here keeps the two derivations from drifting apart.
    fil_in = FIXTURES / "bam" / "se" / "tiny_se.bam"
    fil_out = tmp_path / nam_out

    status = main(
        [
            "--fil_in",
            str(fil_in),
            "--fil_out",
            str(fil_out),
            "--siz_bin",
            "10",
            "--report_n_frg",
            "--report_n_bin",
        ],
    )

    assert status == 0
    assert read_count(tmp_path / f"{base}.n_frg.txt") > 0
    assert read_count(tmp_path / f"{base}.n_bin.txt") > 0


def test_bare_report_flags_match_explicit_paths(tmp_path: Path) -> None:
    fil_in = FIXTURES / "bam" / "se" / "tiny_se.bam"
    dir_bare = tmp_path / "bare"
    dir_expl = tmp_path / "expl"

    dir_bare.mkdir()
    dir_expl.mkdir()

    base_args = [
        "--fil_in",
        str(fil_in),
        "--siz_bin",
        "10",
        "--method",
        "norm",
    ]

    status_bare = main(
        [
            *base_args,
            "--fil_out",
            str(dir_bare / "track.bedGraph.gz"),
            "--report_n_frg",
            "--report_n_bin",
        ],
    )
    status_expl = main(
        [
            *base_args,
            "--fil_out",
            str(dir_expl / "track.bedGraph.gz"),
            "--report_n_frg",
            str(dir_expl / "n.txt"),
            "--report_n_bin",
            str(dir_expl / "l.txt"),
        ],
    )

    bare_n = read_count(dir_bare / "track.n_frg.txt")
    bare_l = read_count(dir_bare / "track.n_bin.txt")
    expl_n = read_count(dir_expl / "n.txt")
    expl_l = read_count(dir_expl / "l.txt")

    assert status_bare == 0
    assert status_expl == 0
    assert bare_n == expl_n
    assert bare_l == expl_l


def test_bare_report_flag_without_fil_out_is_rejected(tmp_path: Path) -> None:
    # Report-only mode has no output path to derive from. The bare form must
    # say so rather than guess a location beside the input alignment.
    fil_in = FIXTURES / "bam" / "se" / "tiny_se.bam"

    with pytest.raises(SystemExit) as excinfo:
        main(["--fil_in", str(fil_in), "--report_n_frg"])

    message = str(excinfo.value)

    assert "--report_n_frg" in message
    assert "--fil_out" in message


def test_explicit_report_path_is_not_overridden(tmp_path: Path) -> None:
    fil_in = FIXTURES / "bam" / "se" / "tiny_se.bam"
    fil_out = tmp_path / "track.bedGraph.gz"
    path_n = tmp_path / "elsewhere.n_frg.txt"

    status = main(
        [
            "--fil_in",
            str(fil_in),
            "--fil_out",
            str(fil_out),
            "--siz_bin",
            "10",
            "--report_n_frg",
            str(path_n),
        ],
    )

    assert status == 0
    assert read_count(path_n) > 0
    assert not (tmp_path / "track.n_frg.txt").exists()


def test_hidden_hyphen_spelling_derives_identically(tmp_path: Path) -> None:
    # The hyphen aliases are separate argparse actions sharing one dest, so an
    # unmirrored 'nargs' / 'const' would make the bare hyphen form behave
    # differently from the bare underscore form.
    fil_in = FIXTURES / "bam" / "se" / "tiny_se.bam"
    fil_out = tmp_path / "track.bedGraph.gz"

    status = main(
        [
            "--fil_in",
            str(fil_in),
            "--fil_out",
            str(fil_out),
            "--siz_bin",
            "10",
            "--report-n-frg",
            "--report-n-bin",
        ],
    )

    assert status == 0
    assert read_count(tmp_path / "track.n_frg.txt") > 0
    assert read_count(tmp_path / "track.n_bin.txt") > 0
