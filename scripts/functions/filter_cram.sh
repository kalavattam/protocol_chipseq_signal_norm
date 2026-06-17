#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Compatibility aliases for the renamed filter_alignment helper.

# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
)/filter_alignment.sh"

function _parse_args_filter_cram() {
    _parse_args_filter_alignment "$@"
}


function _validate_args_filter_cram() {
    _validate_args_filter_alignment "$@"
}


function _check_chr_cram() {
    _check_chr_alignment "$@"
}


function _finalize_cram_filter() {
    _finalize_alignment_filter "$@"
}


function _cleanup_filter_cram_tmp() {
    _cleanup_filter_alignment_tmp "$@"
}


function filter_cram_sc() {
    filter_alignment_sc "$@"
}


function filter_cram_sp() {
    filter_alignment_sp "$@"
}
