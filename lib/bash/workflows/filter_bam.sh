#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: filter_bam.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


# Legacy forwarding functions for the renamed filter_alignment helper.

# shellcheck source=lib/bash/workflows/filter_alignment.sh
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
)/filter_alignment.sh"


function _parse_args_filter_bam() {
    _parse_args_filter_alignment "$@"
}


function _validate_args_filter_bam() {
    _validate_args_filter_alignment "$@"
}


function _check_chr_bam() {
    _check_chr_alignment "$@"
}


function _finalize_bam_filter() {
    _finalize_alignment_filter "$@"
}


function _cleanup_filter_bam_tmp() {
    _cleanup_filter_alignment_tmp "$@"
}


function filter_bam_sc() {
    filter_alignment_sc "$@"
}


function filter_bam_sp() {
    filter_alignment_sp "$@"
}
