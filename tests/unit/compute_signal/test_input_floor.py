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


from pathlib import Path
from shutil import copyfile

import pytest

import protocol_chipseq_signal_norm.cli.compute_input_floor as input_floor_module
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


def test_input_floor_rejects_bin_size_at_or_above_genome() -> None:
    with pytest.raises(InputFloorValidationError):
        compute_input_floor("-", 100, 100, mode="norm")


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
            "--flags-pe",
            "bad",
        ],
    )

    assert status == 1
    assert capsys.readouterr().err == (
        "Error: Invalid FLAG 'bad' in 'flags-pe'. Use decimal (e.g., 99) or "
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


def test_main_preserves_computation_error_diagnostic_and_status(
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = main(
        [
            "--mode",
            "norm",
            "--siz_bin",
            "100",
            "--siz_gen",
            "100",
        ],
    )

    assert status == 1
    assert capsys.readouterr().err == (
        "Warning: 'siz_bin' is a large fraction of 'siz_gen'; 'dep_min' may "
        "be very large.\n"
        "Error: 'siz_bin' must be smaller than 'siz_gen' "
        "(got siz_bin=100, siz_gen=100).\n"
    )
