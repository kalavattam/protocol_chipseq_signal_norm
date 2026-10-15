#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: submit_filter_crams.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
#   GPT-5.6);
# - Anthropic Claude Code (Opus 5.5).
#
# Distributed under the MIT license.


# Compatibility wrapper for the renamed 'filter_alignments' chain.
set -euo pipefail

dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"

# Add this alias to the call path that alignment-file records give
# ('OUTPUT.ALIGNMENT.RECORD').
PIPELINE_CALLS="${PIPELINE_CALLS:+${PIPELINE_CALLS} > }submit_filter_crams.sh"
export PIPELINE_CALLS

exec "${BASH}" "${dir_scr}/submit_filter_alignments.sh" "$@"
