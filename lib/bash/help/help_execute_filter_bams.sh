#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_execute_filter_bams.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in design, development,
# and documentation, with all output reviewed, edited, and approved by the
# author.
#
# Distributed under the MIT license.


# Compatibility alias for the renamed filter_alignments help function.

# shellcheck source=lib/bash/help/help_execute_filter_alignments.sh
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
)/help_execute_filter_alignments.sh"


function help_execute_filter_bams() {
    help_execute_filter_alignments "$@"
}
