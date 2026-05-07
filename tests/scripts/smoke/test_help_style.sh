#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_help_style.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


#TODO:
#+ - detect generic `Arguments:` where only `Keyword arguments:` or only `Positional arguments:` would clearly be more specific
#+   + the generic `Arguments:` check is hard to detect reliably without false positives; make it WARN-only forever
#+ - detect `Options:` as a stale synonym
#+ - detect missing two-blank-line separators before major headings
#+ - detect `--dry-run` documented when `--dry_run` should be canonical


TEST_NAME="help style"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

#  Scan top-level script help text and extracted help files for style drift
files=()
while IFS= read -r file; do
    files+=( "${file}" )
done < <(
    {
        find "${ROOT_REPO}/scripts" -maxdepth 1 -type f -name '*.sh' -print
        find "${ROOT_REPO}/scripts/functions/help" -type f -name '*.sh' -print
    } | sort
)

#  Restrict example line-continuation checks to extracted help files for now
help_files=()
while IFS= read -r file; do
    help_files+=( "${file}" )
done < <(
    find "${ROOT_REPO}/scripts/functions/help" \
        -type f -name '*.sh' -print \
        | sort
)


#  Scan a caller-supplied file list for a fixed pattern
function scan_fixed_files() {
    local label="${1:-}"
    local pattern="${2:-}"

    shift 2

    local scan_files=( "$@" )
    local found=0
    local file line

    for file in "${scan_files[@]}"; do
        while IFS= read -r line; do
            found=1
            rec_warn "${label}: $(rec_relpath "${file}"):${line}"
        done < <(grep -n -- "${pattern}" "${file}" || true)
    done

    if (( found == 0 )); then
        rec_pass "no ${label} findings"
    fi
}


#  Scan the default style-check file set for a fixed pattern
function scan_fixed() {
    local label="${1:-}"
    local pattern="${2:-}"

    scan_fixed_files "${label}" "${pattern}" "${files[@]}"
}


#  Warn on hard tabs in user-facing help text
function scan_tabs() {
    local found=0
    local file line

    for file in "${files[@]}"; do
        while IFS= read -r line; do
            found=1
            rec_warn "hard tab: $(rec_relpath "${file}"):${line}"
        done < <(grep -n $'\t' "${file}" || true)
    done

    if (( found == 0 )); then
        rec_pass "no hard tab findings"
    fi
}


#  Warn when '--rnd' appears as the canonical Usage spelling
function scan_usage_rnd() {
    local found=0
    local file line

    for file in "${files[@]}"; do
        while IFS= read -r line; do
            found=1
            rec_warn "--rnd in Usage block: $(rec_relpath "${file}"):${line}"
        done < <(
            awk '
                /^Usage:/ { in_usage = 1 }
                in_usage && /^$/ { in_usage = 0 }
                in_usage && /--rnd/ { print FNR ":" $0 }
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        rec_pass "no --rnd canonical Usage findings"
    fi
}


#  Run lightweight help-style scans ===========================================
scan_tabs
scan_fixed "literal <flg>" '<flg>'
scan_fixed "literal <mlt>" '<mlt>'
scan_fixed "Markdown triple-backtick fence" '```'
scan_fixed_files \
    "possible shell line continuation in sourced help text" \
    '\\$' \
    "${help_files[@]}"
rec_skip \
    "script-local heredoc/example backslash parsing is a future help-style" \
    "enhancement"
scan_usage_rnd
scan_fixed "stale helper name echo_error" 'echo_error'
scan_fixed "stale helper name echo_warning" 'echo_warning'
scan_fixed "stale helper dependency submit.sh" 'submit\.sh'

finish
