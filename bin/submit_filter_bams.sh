#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: submit_filter_bams.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in design, development, and documentation, with all output
# reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


#  Compatibility wrapper for the renamed 'filter_alignments' chain
set -euo pipefail

dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"

exec "${BASH}" "${dir_scr}/submit_filter_alignments.sh" "$@"
