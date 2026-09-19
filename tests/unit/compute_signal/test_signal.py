#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_signal.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5, Fable 5).
#
# Distributed under the MIT license.


import argparse
import gzip
import shutil
from pathlib import Path

import numpy as np
import pysam
import pytest

from protocol_chipseq_signal_norm.cli.compute_signal import (
    METHOD_CANON,
    METHOD_DEPOSIT,
    calc_sig_chrom_direct_sparse_np,
    get_siz_chr,
    iter_aln_frg,
    iter_idx_frg,
    main,
    parse_args,
    read_to_frg,
    resolve_siz_chr,
)
from protocol_chipseq_signal_norm.utilities.utils_cli import (
    CapArgumentParser,
)

ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tests" / "fixtures" / "compute_signal"


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def collect_indexed_windows(
    fil_in: Path,
    siz_chr: dict[str, int],
    usr_frg: int | None,
    window_size: int,
) -> list[tuple[str, int, int, int]]:
    fragments = []

    for chrom, chrom_size in siz_chr.items():
        for start in range(0, chrom_size, window_size):
            fragments.extend(
                iter_idx_frg(
                    str(fil_in),
                    chrom=chrom,
                    siz_chr=siz_chr,
                    usr_frg=usr_frg,
                    start=start,
                    end=min(start + window_size, chrom_size),
                ),
            )

    return fragments


def test_resolve_siz_chr_rejects_an_empty_header() -> None:
    """
    Sizes come from the BAM/CRAM header alone, and an empty set is an error.
    """

    assert resolve_siz_chr({"chrI": 80}) == {"chrI": 80}

    with pytest.raises(ValueError, match="usable sequence lengths"):
        resolve_siz_chr({})


def test_chr_siz_is_rejected(tmp_path: Path) -> None:
    """
    The retired option is gone outright, with no compatibility retention.
    """

    sizes = tmp_path / "chr.sizes"
    sizes.write_text("chrI\t80\n", encoding="utf-8")

    for spelling in ("--chr_siz", "--chr-siz", "-cs"):
        with pytest.raises(SystemExit):
            parse_args(
                [
                    "--fil_in",
                    "x.bam",
                    "--fil_out",
                    "y.bdg",
                    spelling,
                    str(sizes),
                ],
            )


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
            "--siz_win",
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
    "subdir,usr_frg",
    [
        ("bam/se/tiny_se.bam", 30),
        ("bam/pe/tiny_pe.bam", None),
    ],
)
def test_indexed_chrom_and_window_fetch_match_serial_fragments(
    subdir: str,
    usr_frg: int | None,
) -> None:
    fil_in = FIXTURES / subdir
    siz_chr = get_siz_chr(str(fil_in))
    serial = sorted(
        iter_aln_frg(
            str(fil_in),
            siz_chr=siz_chr,
            usr_frg=usr_frg,
        ),
    )
    idx_chrom = []

    for chrom in siz_chr:
        idx_chrom.extend(
            iter_idx_frg(
                str(fil_in),
                chrom=chrom,
                siz_chr=siz_chr,
                usr_frg=usr_frg,
            ),
        )

    indexed_windows = collect_indexed_windows(
        fil_in,
        siz_chr,
        usr_frg=usr_frg,
        window_size=15,
    )

    assert sorted(idx_chrom) == serial
    assert sorted(indexed_windows) == serial


@pytest.mark.parametrize(
    "subdir,ref_args,usr_frg",
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
    usr_frg: int | None,
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

    if usr_frg is not None:
        base_args.extend(["--usr_frg", usr_frg])

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
            "--siz_win",
            "15",
        ],
    )

    assert default_status == 0
    assert chrom_status == 0
    assert window_status == 0
    assert read_text(out_chrom) == read_text(out_default)
    assert read_text(out_window) == read_text(out_default)


@pytest.mark.parametrize(
    "strat_bed",
    [
        "idx_chrom",
        "idx_win",
    ],
)
def test_compute_signal_hidden_indexed_bed_strategies_match_serial(
    tmp_path: Path,
    strat_bed: str,
) -> None:
    fil_in = FIXTURES / "bam" / "pe" / "tiny_pe.bam"
    out_serial = tmp_path / "serial.bed"
    out_indexed = tmp_path / f"{strat_bed}.bed"

    serial_status = main(
        [
            "--fil_in",
            str(fil_in),
            "--fil_out",
            str(out_serial),
            "--strat_bed",
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
            "--strat_bed",
            strat_bed,
            "--siz_win",
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
            "--strat_bed",
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

    for engine in ("chunk", "idx_chrom"):
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


# Hidden hyphen aliases are separate 'add_argument' calls that must restate the
# primary's 'type', 'choices', and action; argparse enforces none of that. Each
# pair is exercised below, and a completeness guard fails when an alias reaches
# the parser without a row here.
HYPHEN_ALIASES = (
    pytest.param("--fil_in",  "--fil-in",    "value", id="fil_in"),
    pytest.param("--fil_out", "--fil-out",   "value", id="fil_out"),
    pytest.param("--ref_fa",  "--ref-fa",    "value", id="ref_fa"),
    pytest.param(
        "--report_n_bin", "--report-n-bin", "value", id="report_n_bin"
    ),
    pytest.param(
        "--report_n_frg", "--report-n-frg", "value", id="report_n_frg"
    ),
    pytest.param("--scl_fct",  "--scl-fct",  "0.5",   id="scl_fct"),
    pytest.param("--siz_bin",  "--siz-bin",  "7",     id="siz_bin"),
    pytest.param("--siz_win",  "--siz-win",  "7",     id="siz_win"),
    pytest.param("--usr_frg",  "--usr-frg",  "7",     id="usr_frg"),
)


def _alias_parser() -> argparse.ArgumentParser:
    """
    Return the parser built by the module under test.
    """

    captured: dict[str, argparse.ArgumentParser] = {}
    original = CapArgumentParser.parse_args

    def capture(
        self: CapArgumentParser,
        *args: object,
        **kwargs: object,
    ) -> None:
        captured["parser"] = self

        raise SystemExit(0)

    CapArgumentParser.parse_args = capture

    try:
        parse_args(["-fi", "x.bam", "-fo", "y.bdg"])
    except SystemExit:
        pass
    finally:
        CapArgumentParser.parse_args = original

    return captured["parser"]


def _registered_hyphen_aliases() -> set[str]:
    """
    Return hidden hyphen spellings that alias a visible option.

    A hidden option may also carry a hyphen spelling; those are excluded, as
    they alias nothing the user can see.
    """

    parser = _alias_parser()
    visible = {
        action.dest
        for action in parser._actions
        if action.help is not argparse.SUPPRESS
    }

    return {
        option
        for action in parser._actions
        if action.help is argparse.SUPPRESS and action.dest in visible
        for option in action.option_strings
        if option.startswith("--") and "_" not in option
    }


@pytest.mark.parametrize(("primary", "alias", "value"), HYPHEN_ALIASES)
def test_hidden_hyphen_alias_matches_its_primary(
    primary: str,
    alias: str,
    value: str | None,
) -> None:
    """
    Each hidden hyphen alias parses to the primary's value and type.
    """

    base = list(["-fi", "x.bam", "-fo", "y.bdg"])
    supplied = [primary] if value is None else [primary, value]
    aliased = [alias] if value is None else [alias, value]
    destination = primary.lstrip("-")

    from_primary = getattr(parse_args(base + supplied), destination)
    from_alias = getattr(parse_args(base + aliased), destination)

    assert from_primary == from_alias
    assert type(from_primary) is type(from_alias)


def test_every_hidden_hyphen_alias_is_covered() -> None:
    """
    No hidden hyphen alias may exist without a row in 'HYPHEN_ALIASES'.
    """

    covered = {row.values[1] for row in HYPHEN_ALIASES}
    registered = _registered_hyphen_aliases()

    assert registered == covered


def test_hidden_hyphen_aliases_stay_out_of_rendered_help(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Hidden aliases remain usable but are never advertised in help.
    """

    with pytest.raises(SystemExit):
        parse_args(["--help"])

    rendered = capsys.readouterr().out
    registered = _registered_hyphen_aliases()

    for alias in registered:
        assert alias not in rendered


def test_hidden_aliases_satisfy_the_required_input_option() -> None:
    """
    A run supplying only hyphen spellings is accepted.

    'argparse' enforces 'required' per action, so a hidden alias cannot satisfy
    a required primary. The requirement is enforced on the parsed value in
    'main' instead, and this pins that the alias route works.
    """

    args = parse_args(["--fil-in", "input.bam", "--fil-out", "output.bdg"])

    assert args.fil_in == "input.bam"
    assert args.fil_out == "output.bdg"


def test_missing_input_is_still_rejected() -> None:
    """
    Dropping the requirement from argparse must not drop the requirement.
    """

    with pytest.raises(SystemExit) as error:
        main(["--fil_out", "output.bdg"])

    assert "'--fil_in' is required" in str(error.value)


def test_get_siz_chr_reads_header_lengths(tmp_path: Path) -> None:
    """
    Sizes come from the header's reference names and lengths.
    """

    fil_in = tmp_path / "hdr.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "I", "LN": 80}, {"SN": "II", "LN": 40}],
    }

    with pysam.AlignmentFile(str(fil_in), "wb", header=header):
        pass

    assert get_siz_chr(str(fil_in)) == {"I": 80, "II": 40}


def test_get_siz_chr_drops_nonpositive_lengths(tmp_path: Path) -> None:
    """
    A reference whose length is zero is unusable and is dropped.

    This is the case '--chr_siz' once rescued. Nothing rescues it now, so the
    reference simply does not appear, and a read on it fails loudly.
    """

    fil_in = tmp_path / "zero.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "I", "LN": 80}, {"SN": "II", "LN": 0}],
    }

    with pysam.AlignmentFile(str(fil_in), "wb", header=header):
        pass

    assert get_siz_chr(str(fil_in)) == {"I": 80}


def test_read_to_frg_rejects_a_chromosome_absent_from_sizes() -> None:
    """
    The error names the header, not a retired option.
    """

    read = pysam.AlignedSegment()
    read.reference_id = 0
    read.reference_start = 0
    read.query_name = "r1"
    read.query_sequence = "A" * 10
    read.cigarstring = "10M"

    with pytest.raises(ValueError, match="usable sequence lengths"):
        read_to_frg(read, lambda _: "absent", {"I": 80})


# A single fragment laid across a bin grid it does not divide evenly. At
# 'siz_bin' 3 the fragment '[0, 10)' fully covers bins 0, 1 and 2 and reaches
# one base into bin 3, which is the geometry that separates the two deposition
# families: the fractional methods weight that last bin by its overlap, and the
# whole-count methods give it the same 1 as any other touched bin.
PARTIAL_FRAGMENT = {
    "starts": np.array([0], dtype=np.int64),
    "ends": np.array([10], dtype=np.int64),
    "lengths": np.array([10.0], dtype=np.float64),
    "chrom_size": 30,
    "siz_bin": 3,
}


def read_bdg_rows(path: Path) -> list[tuple[int, int, float]]:
    """
    Read one bedGraph track as '(start, end, value)' rows.
    """

    if path.suffix == ".gz":
        text = gzip.decompress(path.read_bytes()).decode("utf-8")
    else:
        text = path.read_text(encoding="utf-8")

    rows = []

    for line in text.splitlines():
        _chrom, start, end, value = line.split("\t")
        rows.append((int(start), int(end), float(value)))

    return rows


def run_method(
    tmp_path: Path,
    method: str,
    extra_args: list[str] | None = None,
    fil_in: Path | None = None,
    name: str | None = None,
) -> list[tuple[int, int, float]]:
    """
    Run one bedGraph method at full precision and return its rows.
    """

    source = fil_in or (FIXTURES / "bam" / "se" / "tiny_se.bam")
    fil_out = tmp_path / f"{name or method}.bedGraph"

    status = main(
        [
            "--fil_in",
            str(source),
            "--fil_out",
            str(fil_out),
            "--method",
            method,
            "--dp",
            "24",
            *(extra_args or []),
        ],
    )

    assert status == 0

    return read_bdg_rows(fil_out)


def test_whole_count_deposits_one_in_a_partly_covered_bin() -> None:
    """
    A whole count lands in a bin the fragment only partly covers.

    This is the test that fails if the whole-count family were built by
    reweighting the fractional accumulator instead of given its own path. A
    per-bin tweak of the fractional engine cannot produce both columns below:
    the fractional deposit is the overlap divided by fragment length, so the
    one-base terminal bin is worth a tenth of a full bin, while the whole-count
    deposit does not know how much of the bin was covered.
    """

    _tag, whole = calc_sig_chrom_direct_sparse_np(
        chrom="I",
        **PARTIAL_FRAGMENT,
        is_len=False,
        is_whole=True,
    )
    _tag, fractional = calc_sig_chrom_direct_sparse_np(
        chrom="I",
        **PARTIAL_FRAGMENT,
        is_len=True,
    )

    np.testing.assert_array_equal(whole[0][1], np.array([0, 3, 6, 9]))
    np.testing.assert_array_equal(fractional[0][1], np.array([0, 3, 6, 9]))

    # Every touched bin takes exactly one, the partial terminal bin included.
    np.testing.assert_array_equal(
        whole[0][2],
        np.array([1.0, 1.0, 1.0, 1.0]),
    )

    # The same terminal bin under the fractional family is worth its one base
    # of overlap out of a fragment length of ten.
    np.testing.assert_allclose(
        fractional[0][2],
        np.array([0.3, 0.3, 0.3, 0.1]),
        rtol=1e-12,
    )

    assert whole[0][2][-1] != pytest.approx(fractional[0][2][-1])


def test_whole_count_refuses_fragment_length_weighting() -> None:
    """
    Asking for both deposition families at once is refused, not resolved.

    The two flags encode three valid states, which leaves a fourth that
    describes no method. Letting it resolve to whichever flag the
    implementation reads first would make a caller error look like a deliberate
    choice.
    """

    with pytest.raises(ValueError, match="no fragment-length weighting"):
        calc_sig_chrom_direct_sparse_np(
            chrom="I",
            **PARTIAL_FRAGMENT,
            is_len=True,
            is_whole=True,
        )


def test_count_and_frag_disagree_on_a_partly_covered_bin(
    tmp_path: Path,
) -> None:
    """
    The same separation survives the whole command-line path.
    """

    count_rows = run_method(tmp_path, "count", ["--siz_bin", "3"])
    frag_rows = run_method(tmp_path, "frag", ["--siz_bin", "3"])

    assert [row[:2] for row in count_rows] == [row[:2] for row in frag_rows]
    assert {row[2] for row in count_rows} == {1.0}

    partial = [row for row in frag_rows if row[0] == 9]

    assert len(partial) == 1
    assert partial[0][2] == pytest.approx(0.1)


def test_whole_count_methods_share_one_bin_support(tmp_path: Path) -> None:
    """
    Both whole-count methods touch the same bins as the fractional family.
    """

    bins_count = [row[:2] for row in run_method(tmp_path, "count")]
    bins_cpm = [row[:2] for row in run_method(tmp_path, "cpm")]
    bins_frag = [row[:2] for row in run_method(tmp_path, "frag")]

    assert bins_count == bins_frag
    assert bins_cpm == bins_frag


@pytest.mark.parametrize("siz_bin", ["3", "10", "7"])
def test_cpm_is_the_count_track_scaled_to_one_million(
    tmp_path: Path,
    siz_bin: str,
) -> None:
    """
    A 'cpm' track is its 'count' track divided by that track's own total.
    """

    count_rows = run_method(tmp_path, "count", ["--siz_bin", siz_bin])
    cpm_rows = run_method(tmp_path, "cpm", ["--siz_bin", siz_bin])

    assert [row[:2] for row in count_rows] == [row[:2] for row in cpm_rows]

    bin_tot = sum(row[2] for row in count_rows)

    assert bin_tot > 0

    expected = [row[2] * 1e6 / bin_tot for row in count_rows]

    np.testing.assert_allclose(
        [row[2] for row in cpm_rows],
        expected,
        rtol=1e-12,
    )

    assert sum(row[2] for row in cpm_rows) == pytest.approx(1e6, rel=1e-12)


def test_scl_fct_scales_the_whole_count_methods(tmp_path: Path) -> None:
    """
    A scale factor means the same thing for the whole-count methods.
    """

    for method in ("count", "cpm"):
        plain = run_method(
            tmp_path,
            method,
            ["--siz_bin", "3"],
            name=f"{method}_plain",
        )
        scaled = run_method(
            tmp_path,
            method,
            ["--siz_bin", "3", "--scl_fct", "2.5"],
            name=f"{method}_scaled",
        )

        np.testing.assert_allclose(
            [row[2] for row in scaled],
            [row[2] * 2.5 for row in plain],
            rtol=1e-12,
        )


@pytest.mark.parametrize("method", ["count", "cpm"])
def test_whole_count_engines_agree(tmp_path: Path, method: str) -> None:
    """
    Window fetch deposits each fragment once, as chromosome fetch does.
    """

    by_chrom = run_method(
        tmp_path,
        method,
        ["--siz_bin", "3"],
        name=f"{method}_chrom",
    )
    by_window = run_method(
        tmp_path,
        method,
        ["--siz_bin", "3", "--engine", "window", "--siz_win", "15"],
        name=f"{method}_window",
    )

    assert by_chrom == by_window


@pytest.mark.parametrize(
    ("alias", "canonical"),
    [
        ("count", "count"),
        ("cpm", "cpm"),
    ],
)
def test_whole_count_aliases_standardize(alias: str, canonical: str) -> None:
    args = parse_args(
        [
            "--fil_in",
            "in.bam",
            "--fil_out",
            "out.bedGraph",
            "--method",
            alias,
        ],
    )

    assert METHOD_CANON[args.method] == canonical


def test_every_canonical_method_has_a_deposition_rule() -> None:
    """
    The deposition map is total over the methods the parser accepts.

    A method reaching the accumulator without a deposition rule would take
    whichever family the lookup defaulted to, which is the silent wrong-object
    failure the two families exist to keep apart. Keeping the map total makes
    that a loud lookup error instead.
    """

    assert set(METHOD_DEPOSIT) == set(METHOD_CANON.values())
    assert set(METHOD_DEPOSIT.values()) == {"fractional", "whole"}


def test_the_method_vocabulary_is_exactly_the_ruled_set() -> None:
    """
    Every accepted spelling is one the vocabulary ruling admits.

    Each alias here is also a spelling a sibling tool mirrors, so an extra one
    is not free: it becomes a second file's edit in another workstream.
    """

    # A list comparison pins membership, count and sequence at once. Sequence
    # matters because argparse renders the choices display straight from this
    # mapping, so a reorder changes published help; the ruled order puts each
    # canonical name before its alias.
    assert list(METHOD_CANON) == [
        "unadj",
        "frag",
        "norm",
        "nc",
        "count",
        "cpm",
    ]

    for retired in ("ct", "cp", "cnt", "c"):
        assert retired not in METHOD_CANON


def test_a_bare_c_method_is_rejected() -> None:
    """
    'c' names no method, because it would name either whole-count member.

    The two differ by a per-sample scalar, so guessing one would hand back a
    plausible-looking track rather than an error.
    """

    assert "c" not in METHOD_CANON

    with pytest.raises(SystemExit):
        parse_args(
            [
                "--fil_in",
                "in.bam",
                "--fil_out",
                "out.bedGraph",
                "--method",
                "c",
            ],
        )


@pytest.mark.parametrize("method", ["count", "cpm"])
def test_whole_count_writes_a_gzipped_bedgraph(
    tmp_path: Path,
    method: str,
) -> None:
    fil_out = tmp_path / f"{method}.bedGraph.gz"

    status = main(
        [
            "--fil_in",
            str(FIXTURES / "bam" / "se" / "tiny_se.bam"),
            "--fil_out",
            str(fil_out),
            "--method",
            method,
            "--siz_bin",
            "3",
        ],
    )

    assert status == 0
    assert read_bdg_rows(fil_out)


def test_method_help_names_every_canonical_method(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Rendered help names every method the parser accepts.

    Deriving the expectation from 'METHOD_CANON' rather than listing methods
    here means a method added without a help row fails at this assertion,
    instead of shipping a method the help never mentions.
    """

    with pytest.raises(SystemExit):
        parse_args(["--help"])

    rendered = capsys.readouterr().out

    for canonical in sorted(set(METHOD_CANON.values())):
        assert f"'{canonical}'" in rendered

    for spelling in METHOD_CANON:
        assert f"'{spelling}'" in rendered


def test_method_help_names_both_deposition_families(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Rendered help separates the fractional and whole-count families.

    The families are what a reader has to understand before the individual
    methods mean anything, so their absence is a documentation defect even when
    every method is listed.
    """

    with pytest.raises(SystemExit):
        parse_args(["--help"])

    rendered = capsys.readouterr().out

    assert "Fractional methods deposit" in rendered
    assert "Whole-count methods deposit" in rendered
    assert "'--report_n_bin'" in rendered

