#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_io.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


import gzip
import io
from pathlib import Path

import pytest

from protocol_chipseq_signal_norm.utilities.utils_io import (
    DEF_SKP_PFX,
    ensure_single_stdin,
    is_header,
    open_in,
    open_out,
    parse_skp_pfx,
    read_data_line,
)


def test_parse_skp_pfx_special_cases() -> None:
    assert parse_skp_pfx(None) == DEF_SKP_PFX
    assert parse_skp_pfx("__default__") == DEF_SKP_PFX
    assert parse_skp_pfx("") == ()
    assert parse_skp_pfx("#,track, custom ") == ("#", "track", "custom")


def test_is_header_and_read_data_line_skip_metadata() -> None:
    lines = io.StringIO("# note\ntrack name=x\n\nchrI 0 10 1\n")

    assert is_header("  browser position chrI:1-10\n")
    assert read_data_line(lines, DEF_SKP_PFX) == "chrI 0 10 1"


def test_open_in_and_open_out_handle_gzip(tmp_path: Path) -> None:
    path = tmp_path / "out.txt.gz"

    with open_out(str(path)) as handle:
        handle.write("hello\n")

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        assert handle.read() == "hello\n"

    with open_in(str(path)) as handle:
        assert handle.read() == "hello\n"


def test_ensure_single_stdin_rejects_multiple_stdin_paths() -> None:
    with pytest.raises(ValueError, match="At most one"):
        ensure_single_stdin(["-", "file.bdg", "-"])
