#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: stub_wget_ok.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5.5).
#
# Distributed under the MIT license.


set -euo pipefail

output=""
while (( $# > 0 )); do
    case "${1}" in
        -O)
            output="${2}"
            shift 2
            ;;

        *)
            url="${1}"
            shift
            ;;
    esac
done

printf '%s\n' "${url:-}" > "${output%/*}/wget_url.txt"
printf '%s\n' '@synthetic-fastq' 'ACGT' '+' 'IIII' | gzip -c > "${output}"
