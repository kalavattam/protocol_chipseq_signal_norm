#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Compatibility aliases for the renamed filter_alignment helper.

# shellcheck disable=SC1091
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)/filter_alignment.sh"

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
