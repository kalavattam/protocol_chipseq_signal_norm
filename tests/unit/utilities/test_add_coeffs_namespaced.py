#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_add_coeffs_namespaced.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


from pathlib import Path

import pytest

from protocol_chipseq_signal_norm.cli.add_coeffs_namespaced import (
    gmean,
    hmean,
    main,
    median,
    safe_refs,
)


def test_reference_summary_helpers() -> None:
    assert median([3, 1, 2]) == 2
    assert gmean([1, 4]) == 2
    assert hmean([1, 2]) == pytest.approx(4 / 3)


def test_safe_refs_omits_unavailable_means(
    capsys: pytest.CaptureFixture[str],
) -> None:
    refs = safe_refs([0.0, 2.0], "coef")

    assert refs["min"] == 0.0
    assert refs["max"] == 2.0
    assert "gmean" not in refs
    assert "hmean" not in refs
    assert "unavailable" in capsys.readouterr().err


def test_main_augments_namespaced_coefficients(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    path = tmp_path / "coefficients.tsv"
    path.write_text(
        "sample\tpair_s\tpair_alpha_rxinput\tip_alpha_ip\tin_alpha_in\n"
        "A\t1\t2\t4\t8\n"
        "B\t2\t4\t8\t16\n",
        encoding="utf-8",
    )

    assert main([str(path), "--dp", "2"]) == 0

    out = capsys.readouterr().out

    assert "pair_s_median" in out.splitlines()[0]
    assert out.splitlines()[1].startswith("A\t1.00\t1.00\t0.67")
