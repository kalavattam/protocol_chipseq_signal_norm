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


import hashlib
import inspect
import re
from itertools import pairwise
from pathlib import Path
from shutil import copyfile

import pytest

from protocol_chipseq_signal_norm.cli import (
    compute_input_floor as input_floor_module,
)
from protocol_chipseq_signal_norm.cli.compute_input_floor import (
    AlignmentReadError,
    FlagParseError,
    InputFloorValidationError,
    _count_alignment_records,
    _count_bed_records,
    _validate_data_arguments,
    compute_input_floor,
    infer_input_format,
    main,
    parse_args,
    parse_flag_csv,
)
from protocol_chipseq_signal_norm.utilities.utils_io import DEF_SKP_PFX

ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tests" / "fixtures" / "compute_signal"
BAM_SE = FIXTURES / "bam" / "se" / "tiny_se.bam"
CRAM_SE = FIXTURES / "cram" / "se" / "tiny_se.cram"
REFERENCE_FASTA = FIXTURES / "reference" / "tiny.fa"

MIXED_CASE_FORMAT_HINTS = (
    ("bAm", "bam"),
    ("CrAm", "cram"),
    ("bEd", "bed"),
    ("bedGraph", "bedgraph"),
    ("bedgraph", "bedgraph"),
    ("BeDgRaPh", "bedgraph"),
    ("BdG", "bedgraph"),
    ("bG", "bedgraph"),
)
DIST_DIMENSION_CASES = (
    (10, 100),
    (0, 100),
    (10, 0),
    (-1, 100),
    (10, -1),
    (100, 100),
    (101, 100),
    (51, 100),
)
INVALID_MODE_DIMENSION_CASES = (
    (0, 100, "Error: 'siz_bin' must be positive (got siz_bin=0).\n"),
    (10, 0, "Error: 'siz_gen' must be positive (got siz_gen=0).\n"),
    (-1, 100, "Error: 'siz_bin' must be positive (got siz_bin=-1).\n"),
    (10, -1, "Error: 'siz_gen' must be positive (got siz_gen=-1).\n"),
    (
        100,
        100,
        "Error: 'siz_bin' must be smaller than 'siz_gen' "
        "(got siz_bin=100, siz_gen=100).\n",
    ),
    (
        101,
        100,
        "Error: 'siz_bin' must be smaller than 'siz_gen' "
        "(got siz_bin=101, siz_gen=100).\n",
    ),
)


def test_parse_flag_csv_accepts_decimal_and_hex() -> None:
    assert parse_flag_csv("99, 0x400", "flags") == {99, 1024}


def test_parse_flag_csv_rejects_invalid_token() -> None:
    with pytest.raises(
        FlagParseError,
        match="Error: Invalid FLAG 'bad' in 'flags'",
    ):
        parse_flag_csv("99, bad", "flags")


def test_parse_flag_csv_rejects_empty_list() -> None:
    with pytest.raises(
        FlagParseError,
        match="Error: No valid FLAGs parsed for 'flags'",
    ):
        parse_flag_csv(",", "flags")


def test_infer_input_format_recognizes_suffixes_and_stdin_hints() -> None:
    assert infer_input_format("x.bam") == "bam"
    assert infer_input_format("x.cram") == "cram"
    assert infer_input_format("x.bed.gz") == "bed"
    assert infer_input_format("x.bdg") == "bedgraph"
    assert infer_input_format("-", "cram") == "cram"
    assert infer_input_format("-", "bedGraph") == "bedgraph"
    assert infer_input_format("x.txt") == "other"


@pytest.mark.parametrize(
    ("hint", "expected"),
    MIXED_CASE_FORMAT_HINTS,
)
def test_input_format_hints_are_case_insensitive_and_canonical(
    hint: str,
    expected: str,
) -> None:
    assert infer_input_format("-", hint) == expected

    args = parse_args(
        [
            "--mode",
            "dist",
            "--fil_in",
            "-",
            "--infmt",
            hint,
        ],
    )

    assert args.infmt == expected


def test_named_path_format_inference_ignores_explicit_hint() -> None:
    assert infer_input_format("x.BAM", "BeDgRaPh") == "bam"
    assert infer_input_format("x.BEDGRAPH.GZ", "CrAm") == "bedgraph"


def test_compute_input_floor_accepts_mixed_case_direct_hint(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        input_floor_module,
        "iter_vals_bdg",
        lambda *args, **kwargs: iter([2.0]),
    )

    result = compute_input_floor(
        "-",
        10,
        100,
        mode="dist",
        infmt="BeDgRaPh",
    )

    assert result == 2.0


def test_main_accepts_mixed_case_cli_hint(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr(
        input_floor_module,
        "iter_vals_bdg",
        lambda *args, **kwargs: iter([2.0]),
    )

    status = main(
        [
            "--mode",
            "dist",
            "--fil_in",
            "-",
            "--infmt",
            "BdG",
        ],
    )
    captured = capsys.readouterr()

    assert status == 0
    assert captured.out == "2\n"
    assert captured.err == ""


def test_help_uses_approved_applicability_wording(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as error:
        parse_args(["--help"])

    captured = capsys.readouterr()
    normalized_help = " ".join(captured.out.split())
    applicability = (
        "Use 'dist' for new analyses. Use 'frag' or 'norm' when reproducing "
        "the fragment-normalization or normalized-coverage floor calculations "
        "used in the Dickson/siQ-ChIP and *Bio-protocol* workflows."
    )

    assert error.value.code == 0
    assert applicability in normalized_help
    assert "- dist: Recommended for new analyses." in captured.out
    assert "Retained to reproduce the legacy workflow." not in captured.out
    assert "legacy-reproducible" not in captured.out
    assert "The command returns one scalar 'dep_min'." in captured.out
    assert "PMF" not in captured.out
    assert "CDF" not in captured.out
    assert "probability mass function" not in captured.out
    assert "cumulative distribution function" not in captured.out


def test_help_shows_one_public_bedgraph_spelling(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as error:
        parse_args(["--help"])

    captured = capsys.readouterr()
    public_choices = "{bam,cram,bed,bedGraph,bdg,bg}"
    public_choice_count = captured.out.count(public_choices)

    assert error.value.code == 0
    assert public_choice_count == 2
    assert "bedGraph,bedgraph" not in captured.out
    assert "bedgraph,bedGraph" not in captured.out


def test_help_examples_use_exact_bash_blocks(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as error:
        parse_args(["--help"])

    captured = capsys.readouterr()
    expected_examples = (
        "Examples\n"
        "--------\n"
        "  1. Compute a first-percentile floor from bedGraph values.\n"
        "    '''bash\n"
        "    compute_input_floor \\\n"
        "        --mode dist \\\n"
        "        --fil_in signal.bdg \\\n"
        "        --method qntl_nz \\\n"
        "        --qntl_nz 1\n"
        "    '''\n\n"
        "  2. Compute a normalized floor from explicit dimensions.\n"
        "    '''bash\n"
        "    compute_input_floor --mode norm --siz_bin 30 "
        "--siz_gen 12157105\n"
        "    '''\n"
    )

    assert error.value.code == 0
    assert captured.out.endswith(expected_examples)
    assert expected_examples.count("\\\n") == 4


def test_help_output_is_byte_exact(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as error:
        parse_args(["--help"])

    rendered = capsys.readouterr().out.encode()

    assert error.value.code == 0
    assert len(rendered) == 5719
    assert hashlib.sha256(rendered).hexdigest() == (
        "f905d67d248556fae5962375b151126e270dc66601dac7ec536662d5fb45836f"
    )


def _first_indivisible_unit(line: str) -> str:
    tokens = line.strip().split()
    unit = [tokens.pop(0)]
    quote_is_open = unit[0].count("'") % 2 == 1

    while quote_is_open:
        token = tokens.pop(0)
        unit.append(token)
        quote_is_open = token.count("'") % 2 == 0

    return " ".join(unit)


@pytest.mark.parametrize(
    ("option", "destination"),
    [
        ("--flags_pe", "flags_pe"),
        ("--flags-pe", "flags_pe"),
        ("--flags_se", "flags_se"),
        ("--flags-se", "flags_se"),
    ],
)
def test_flag_spellings_share_canonical_destinations(
    option: str,
    destination: str,
) -> None:
    args = parse_args([option, "99"])

    assert getattr(args, destination) == "99"


def test_help_displays_only_canonical_flag_spellings(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as error:
        parse_args(["--help"])

    help_text = capsys.readouterr().out

    assert error.value.code == 0
    assert "--flags_pe" in help_text
    assert "--flags_se" in help_text
    assert "--flags-pe" not in help_text
    assert "--flags-se" not in help_text


def test_help_uses_approved_semantic_cli_order(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as error:
        parse_args(["--help"])

    help_text = capsys.readouterr().out
    usage = help_text.split("\n\n", maxsplit=1)[0]
    options = (
        "\n"
        + help_text.split("Options\n-------\n", maxsplit=1)[1].split(
            "\nExamples\n",
            maxsplit=1,
        )[0]
    )
    usage_tokens = (
        "[--help]",
        "[--verbose]",
        "[--mode {dist,frag,norm}]",
        "[--fil_in FIL_IN]",
        "[--infmt {bam,cram,bed,bedGraph,bdg,bg}]",
        "[--ref_fa REF_FA]",
        "[--skp_pfx SKP_PFX]",
        "[--method {qntl_nz,frc_mdn_nz,frc_avg_nz,min_nz}]",
        "[--qntl_nz QNTL_NZ]",
        "[--coef COEF]",
        "[--eps EPS]",
        "[--mode_nz {closed,open,off}]",
        "[--floor FLOOR]",
        "[--siz_bin SIZ_BIN]",
        "[--siz_gen SIZ_GEN]",
        "[--flags_pe FLAGS_PE]",
        "[--flags_se FLAGS_SE]",
        "[--dp DP]",
    )
    option_tokens = (
        "\n  -h, --help",
        "\n  -v, --verbose",
        "\n  -md, --mode {dist,frag,norm}",
        "\n  -fi, --fil_in FIL_IN",
        "\n  -if, --infmt {bam,cram,bed,bedGraph,bdg,bg}",
        "\n  -rf, --ref_fa REF_FA",
        "\n  -sp, --skp_pfx SKP_PFX",
        "\n  -m, --method {qntl_nz,frc_mdn_nz,frc_avg_nz,min_nz}",
        "\n  -qn, --qntl_nz QNTL_NZ",
        "\n  -c, --coef COEF",
        "\n  -e, --eps EPS",
        "\n  -mn, --mode_nz {closed,open,off}",
        "\n  -f, --floor FLOOR",
        "\n  -sb, --siz_bin SIZ_BIN",
        "\n  -sg, --siz_gen SIZ_GEN",
        "\n  -fp, --flags_pe FLAGS_PE",
        "\n  -fs, --flags_se FLAGS_SE",
        "\n  -dp, --dp DP",
    )

    assert error.value.code == 0
    assert [usage.index(token) for token in usage_tokens] == sorted(
        usage.index(token) for token in usage_tokens
    )
    assert [options.index(token) for token in option_tokens] == sorted(
        options.index(token) for token in option_tokens
    )


def test_verbose_argument_report_uses_applicable_semantic_order(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    bedgraph = tmp_path / "input.bdg"
    bedgraph.write_text("chrI\t0\t10\t2\n", encoding="utf-8")

    status = main(
        [
            "--verbose",
            "--mode",
            "dist",
            "--fil_in",
            str(bedgraph),
            "--infmt",
            "BdG",
            "--ref_fa",
            "unused.fa",
            "--skp_pfx",
            "#",
        ],
    )
    captured = capsys.readouterr()
    reported_options = [
        line.split(maxsplit=1)[0]
        for line in captured.err.splitlines()
        if line.startswith("--")
    ]

    assert status == 0
    assert captured.out == "2\n"
    assert reported_options == [
        "--verbose",
        "--mode",
        "--fil_in",
        "--infmt",
        "--ref_fa",
        "--skp_pfx",
        "--method",
        "--qntl_nz",
        "--coef",
        "--eps",
        "--mode_nz",
        "--floor",
        "--siz_bin",
        "--siz_gen",
        "--dp",
    ]


def test_callable_docstring_summary_matches_signature() -> None:
    docstring = inspect.getdoc(compute_input_floor)
    source = inspect.getsource(compute_input_floor)

    assert docstring is not None

    extended_summary = " ".join(
        docstring.split(
            "\n\nParameters\n----------",
            maxsplit=1,
        )[0].splitlines(),
    )

    assert (
        "The 'mode' parameter selects 'dist', 'frag', or 'norm'. Use 'dist' "
        "for new analyses. Use 'frag' or 'norm' when reproducing the "
        "fragment-normalization or normalized-coverage floor calculations "
        "used in the Dickson/siQ-ChIP and *Bio-protocol* workflows."
    ) in extended_summary
    assert docstring.index("The 'mode' parameter") < docstring.index(
        "Parameters\n----------",
    )

    parameter_section = docstring.split(
        "Parameters\n----------\n",
        maxsplit=1,
    )[1].split("\n\nReturns\n-------", maxsplit=1)[0]
    documented_parameters = re.findall(
        r"^([a-z_]+) : ",
        parameter_section,
        flags=re.MULTILINE,
    )
    signature_parameters = list(
        inspect.signature(compute_input_floor).parameters,
    )
    docstring_start = source.index('    """')
    docstring_end = source.index('    """', docstring_start + 7) + 7
    source_docstring = source[docstring_start:docstring_end]
    source_lines = source_docstring.splitlines()
    prose_paragraphs: list[list[str]] = [[source_lines[1]]]
    extended_start = source_lines.index(
        "    The 'mode' parameter selects 'dist', 'frag', or 'norm'. Use "
        "'dist' for new",
    )
    parameters_start = source_lines.index("    Parameters")
    prose_paragraphs.append(
        source_lines[extended_start : parameters_start - 1],
    )
    description: list[str] = []

    for line in source_lines[parameters_start:]:
        if line.startswith("        "):
            description.append(line)
        elif description:
            prose_paragraphs.append(description)
            description = []

    assert not description
    assert len(prose_paragraphs) == 20

    for paragraph in prose_paragraphs:
        assert all(len(line) <= 79 for line in paragraph)

        for previous, continuation in pairwise(paragraph):
            first_unit = _first_indivisible_unit(continuation)

            assert len(previous) + 1 + len(first_unit) > 79

    assert documented_parameters == signature_parameters
    assert (
        "    fil_in : str\n"
        "        Input path. 'dist' accepts 'bedGraph', 'bdg', or 'bg', "
        "optionally with\n"
        "        '.gz'; 'frag' accepts 'bam', 'cram', 'bed', or 'bed.gz'; "
        "'norm' ignores\n"
        "        'fil_in'. For 'dist' and 'frag', '-' reads standard input "
        "and requires\n"
        "        'infmt'."
    ) in source_docstring
    assert len(docstring.encode()) == 3416
    assert hashlib.sha256(docstring.encode()).hexdigest() == (
        "644927a73d710242a05cbc3f89fe672f146de70f56613bd4e92310901cea34a0"
    )
    assert hashlib.sha256(
        " ".join(docstring.split()).encode(),
    ).hexdigest() == (
        "360deeb6da232e614a1404499eecddc68ae064db540ca2553fdcf7b64667427c"
    )
    assert hashlib.sha256((source_docstring + "\n").encode()).hexdigest() == (
        "d1e0bc008338922769443de1990f4bd8275f2aff74b96af2281c2d9184fbfdb8"
    )


def test_main_trims_noninformative_finite_output(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    bedgraph = tmp_path / "input.bdg"
    bedgraph.write_text("chrI\t0\t10\t1.23\n", encoding="utf-8")

    status = main(
        [
            "--mode",
            "dist",
            "--fil_in",
            str(bedgraph),
            "--dp",
            "4",
        ],
    )
    captured = capsys.readouterr()

    assert status == 0
    assert captured.out == "1.23\n"
    assert captured.err == ""


def test_main_renders_normalized_floor_with_bounded_precision(
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = main(
        [
            "--mode",
            "norm",
            "--siz_bin",
            "10",
            "--siz_gen",
            "100",
            "--dp",
            "4",
        ],
    )
    captured = capsys.readouterr()

    assert status == 0
    assert captured.out == "0.1111\n"
    assert captured.err == ""


def test_main_normalizes_negative_zero_output(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr(
        input_floor_module,
        "_compute_input_floor_from_args",
        lambda *args, **kwargs: -0.0001,
    )

    status = main(
        [
            "--mode",
            "norm",
            "--siz_bin",
            "10",
            "--siz_gen",
            "100",
            "--dp",
            "2",
        ],
    )
    captured = capsys.readouterr()

    assert status == 0
    assert captured.out == "0\n"
    assert captured.err == ""


def test_count_aln_bed_skips_headers_and_blanks(tmp_path: Path) -> None:
    path = tmp_path / "alignments.bed"
    path.write_text(
        "# h\ntrack name=x\nchrI\t0\t10\n\n   \nchrI\t10\t20\n",
        encoding="utf-8",
    )

    assert _count_bed_records(str(path), DEF_SKP_PFX) == 2


def test_alignment_helpers_raise_specific_reusable_errors(
    tmp_path: Path,
) -> None:
    missing_bed = tmp_path / "missing.bed"
    missing_alignment = tmp_path / "missing.bam"

    with pytest.raises(AlignmentReadError, match="Cannot open bed file"):
        _count_bed_records(str(missing_bed), DEF_SKP_PFX)

    with pytest.raises(
        AlignmentReadError,
        match="Cannot process bam alignment",
    ):
        _count_alignment_records(
            str(missing_alignment),
            alignment_format="bam",
        )


def test_alignment_counting_supports_cram_without_an_index(
    tmp_path: Path,
) -> None:
    copied_cram = tmp_path / "tiny_se.cram"
    copyfile(CRAM_SE, copied_cram)

    cram_count = _count_alignment_records(
        str(CRAM_SE),
        alignment_format="cram",
        ref_fa=str(REFERENCE_FASTA),
    )
    unindexed_cram_count = _count_alignment_records(
        str(copied_cram),
        alignment_format="cram",
        ref_fa=str(REFERENCE_FASTA),
    )
    bam_count = _count_alignment_records(
        str(BAM_SE),
        alignment_format="bam",
    )

    assert cram_count == 2
    assert unindexed_cram_count == 2
    assert bam_count == 2


def test_cram_reference_validation_preserves_anticipated_errors(
    tmp_path: Path,
) -> None:
    arguments = (
        str(CRAM_SE),
        10,
        100,
    )

    with pytest.raises(
        InputFloorValidationError,
        match="CRAM input requires '--ref_fa'",
    ):
        compute_input_floor(*arguments, mode="frag")

    missing_reference = tmp_path / "missing.fa"

    with pytest.raises(
        InputFloorValidationError,
        match="Error: Reference FASTA not found:",
    ):
        compute_input_floor(
            *arguments,
            mode="frag",
            ref_fa=str(missing_reference),
        )


def test_cram_stdin_hint_reaches_reader_as_resolved_format(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    observed: dict[str, object] = {}

    def count_records(
        alignment_path: str,
        paired_flags: set[int] | None = None,
        single_flags: set[int] | None = None,
        *,
        alignment_format: str,
        ref_fa: str | None = None,
    ) -> int:
        observed.update(
            {
                "alignment_path": alignment_path,
                "paired_flags": paired_flags,
                "single_flags": single_flags,
                "alignment_format": alignment_format,
                "ref_fa": ref_fa,
            },
        )

        return 2

    monkeypatch.setattr(
        input_floor_module,
        "_count_alignment_records",
        count_records,
    )

    result = compute_input_floor(
        "-",
        10,
        100,
        mode="frag",
        infmt="cram",
        ref_fa=str(REFERENCE_FASTA),
    )

    assert result == pytest.approx(2 / 9)
    assert observed == {
        "alignment_path": "-",
        "paired_flags": None,
        "single_flags": None,
        "alignment_format": "cram",
        "ref_fa": str(REFERENCE_FASTA),
    }


def test_data_validation_raises_without_owning_process_exit() -> None:
    args = parse_args(["--mode", "dist"])

    with pytest.raises(
        InputFloorValidationError,
        match="'--fil_in' is required",
    ):
        _validate_data_arguments(args)


def test_compute_input_floor_norm_and_dist_modes(tmp_path: Path) -> None:
    bdg = tmp_path / "input.bdg"
    bdg.write_text(
        "chrI 0 10 0\n"
        "chrI 10 20 -1\n"
        "chrI 20 30 1\n"
        "chrI 30 40 2\n"
        "chrI 40 50 10\n",
        encoding="utf-8",
    )

    normalized_floor = compute_input_floor("-", 10, 100, mode="norm")

    observed_floor = compute_input_floor(
        str(bdg),
        10,
        100,
        mode="dist",
        method="qntl_nz",
        qntl_nz=50,
    )

    assert normalized_floor == pytest.approx(1 / 9)
    assert observed_floor == 2.0


@pytest.mark.parametrize(
    ("siz_bin", "siz_gen"),
    DIST_DIMENSION_CASES,
)
def test_compute_input_floor_dist_ignores_all_dimension_cases(
    tmp_path: Path,
    siz_bin: int,
    siz_gen: int,
) -> None:
    bedgraph = tmp_path / "input.bdg"
    bedgraph.write_text("chrI\t0\t10\t2\n", encoding="utf-8")

    result = compute_input_floor(
        str(bedgraph),
        siz_bin,
        siz_gen,
        mode="dist",
    )

    assert result == 2.0


@pytest.mark.parametrize("mode", ("frag", "norm"))
@pytest.mark.parametrize(
    ("siz_bin", "siz_gen", "error_text"),
    INVALID_MODE_DIMENSION_CASES,
)
def test_compute_input_floor_frag_and_norm_reject_invalid_dimensions(
    mode: str,
    siz_bin: int,
    siz_gen: int,
    error_text: str,
) -> None:
    with pytest.raises(
        InputFloorValidationError,
        match=re.escape(error_text.rstrip()),
    ):
        compute_input_floor(
            "-",
            siz_bin,
            siz_gen,
            mode=mode,
        )


@pytest.mark.parametrize("mode", ("frag", "norm"))
def test_compute_input_floor_frag_and_norm_accept_large_fraction(
    tmp_path: Path,
    mode: str,
) -> None:
    bed = tmp_path / "input.bed"
    bed.write_text("chrI\t0\t10\n", encoding="utf-8")
    input_path = str(bed) if mode == "frag" else "-"

    result = compute_input_floor(
        input_path,
        51,
        100,
        mode=mode,
    )

    assert result == pytest.approx(51 / 49)


def test_main_preserves_flag_error_diagnostic_and_status(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    bed = tmp_path / "input.bed"
    bed.write_text("chrI\t0\t10\n", encoding="utf-8")

    status = main(
        [
            "--mode",
            "frag",
            "--fil_in",
            str(bed),
            "--flags_pe",
            "bad",
        ],
    )

    assert status == 1
    assert capsys.readouterr().err == (
        "Error: Invalid FLAG 'bad' in 'flags_pe'. Use decimal (e.g., 99) or "
        "hex (e.g., 0x63).\n"
    )


def test_main_translates_data_validation_error_to_diagnostic_and_status(
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = main(["--mode", "dist"])

    assert status == 1
    assert capsys.readouterr().err == (
        "Error: '--fil_in' is required when '--mode dist' or '--mode frag'.\n"
    )


def test_main_counts_reference_backed_cram(
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = main(
        [
            "--mode",
            "frag",
            "--fil_in",
            str(CRAM_SE),
            "--ref_fa",
            str(REFERENCE_FASTA),
            "--siz_bin",
            "10",
            "--siz_gen",
            "100",
        ],
    )
    captured = capsys.readouterr()

    assert status == 0
    assert float(captured.out) == pytest.approx(2 / 9)
    assert captured.err == ""


def test_main_reports_missing_cram_reference_as_anticipated_error(
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = main(
        [
            "--mode",
            "frag",
            "--fil_in",
            str(CRAM_SE),
        ],
    )

    assert status == 1
    assert capsys.readouterr().err == (
        "Error: CRAM input requires '--ref_fa' for reference-backed "
        "decoding.\n"
    )


@pytest.mark.parametrize(
    ("siz_bin", "siz_gen"),
    DIST_DIMENSION_CASES,
)
def test_main_dist_ignores_all_dimension_cases_without_warning(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    siz_bin: int,
    siz_gen: int,
) -> None:
    bedgraph = tmp_path / "input.bdg"
    bedgraph.write_text("chrI\t0\t10\t2\n", encoding="utf-8")

    status = main(
        [
            "--mode",
            "dist",
            "--fil_in",
            str(bedgraph),
            "--siz_bin",
            str(siz_bin),
            "--siz_gen",
            str(siz_gen),
            "--dp",
            "1",
        ],
    )

    captured = capsys.readouterr()

    assert status == 0
    assert captured.out == "2\n"
    assert captured.err == ""


@pytest.mark.parametrize("mode", ("frag", "norm"))
@pytest.mark.parametrize(
    ("siz_bin", "siz_gen", "error_text"),
    INVALID_MODE_DIMENSION_CASES,
)
def test_main_frag_and_norm_reject_invalid_dimensions(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    mode: str,
    siz_bin: int,
    siz_gen: int,
    error_text: str,
) -> None:
    bed = tmp_path / "input.bed"
    bed.write_text("chrI\t0\t10\n", encoding="utf-8")
    arguments = [
        "--mode",
        mode,
        "--siz_bin",
        str(siz_bin),
        "--siz_gen",
        str(siz_gen),
    ]

    if mode == "frag":
        arguments.extend(["--fil_in", str(bed)])

    status = main(arguments)
    captured = capsys.readouterr()

    assert status == 1
    assert captured.out == ""
    assert captured.err == error_text


@pytest.mark.parametrize("mode", ("frag", "norm"))
def test_main_frag_and_norm_warn_and_succeed_for_large_fraction(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    mode: str,
) -> None:
    bed = tmp_path / "input.bed"
    bed.write_text("chrI\t0\t10\n", encoding="utf-8")
    arguments = [
        "--mode",
        mode,
        "--siz_bin",
        "51",
        "--siz_gen",
        "100",
        "--dp",
        "4",
    ]

    if mode == "frag":
        arguments.extend(["--fil_in", str(bed)])

    status = main(arguments)
    captured = capsys.readouterr()

    assert status == 0
    assert captured.out == "1.0408\n"
    assert captured.err == (
        "Warning: 'siz_bin' is a large fraction of 'siz_gen'; 'dep_min' may "
        "be very large.\n"
    )
