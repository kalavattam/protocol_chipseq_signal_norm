#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_shell_wrapper_style.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="shell wrapper style"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Fail when eligible top-level CLIs keep main help inline.
function scan_inline_top_level_help() {
    local found=0
    local file line

    for file in "${files_top[@]}"; do
        case "$(basename "${file}")" in
            install_envs_entrypoint.sh)
                continue
                ;;
        esac

        while IFS= read -r line; do
            found=1
            record_fail \
                "inline top-level main help:" \
                "$(print_relpath "${file}"):${line}"
        done < <(
            grep -nE \
                '^function (show_help_main|help_[A-Za-z0-9_]+)\(\) \{' \
                "${file}" || true
        )
    done

    if (( found == 0 )); then
        record_pass "no inline top-level main-help findings"
    fi
}


# Fail when top-level CLIs lack externalized main help.
function scan_external_main_help() {
    local found=0
    local file base stem fil_hlp fun_hlp

    for file in "${files_top[@]}"; do
        base="$(basename "${file}")"
        stem="${base%.sh}"

        case "${base}" in
            install_envs_entrypoint.sh|\
            execute_filter_bams.sh|execute_filter_crams.sh|\
            submit_filter_bams.sh|submit_filter_crams.sh)
                continue
                ;;
        esac

        fil_hlp="${ROOT_REPO}/lib/bash/help/help_${stem}.sh"
        fun_hlp="help_${stem}"

        if [[ ! -f "${fil_hlp}" ]]; then
            found=1
            record_fail \
                "top-level CLI missing external main help:" \
                "$(print_relpath "${file}") ->" \
                "lib/bash/help/help_${stem}.sh"
            continue
        fi

        if ! \
            grep -q "function ${fun_hlp}()" "${fil_hlp}"
        then
            found=1
            record_fail \
                "external main help missing expected function:" \
                "$(print_relpath "${fil_hlp}"):${fun_hlp}"
        fi

        if ! \
            grep -q "help_${stem}" "${file}"
        then
            found=1
            record_fail \
                "top-level CLI does not call external main help:" \
                "$(print_relpath "${file}")"
        fi
    done

    if (( found == 0 )); then
        record_pass "no external main-help findings"
    fi
}


# Fail when real submit wrappers drift from the required '--dir_scr' contract.
function scan_submit_dir_scr_policy() {
    local found=0
    local file base stem fil_hlp

    for file in "${files_submit[@]}"; do
        base="$(basename "${file}")"
        stem="${base%.sh}"

        case "${base}" in
            submit_filter_bams.sh|submit_filter_crams.sh)
                continue
                ;;
        esac

        fil_hlp="${ROOT_REPO}/lib/bash/help/help_${stem}.sh"

        if [[ ! -f "${fil_hlp}" ]]; then
            found=1
            record_fail \
                "submit wrapper missing external main help:" \
                "$(print_relpath "${file}")"
            continue
        fi

        if ! \
            grep -q "help_${stem}.sh" "${file}"
        then
            found=1
            record_fail \
                "submit wrapper missing bootstrap main-help source:" \
                "$(print_relpath "${file}")"
        fi

        if ! \
            grep -q "help/${stem/#submit_/help_submit_}" "${file}"
        then
            found=1
            record_fail \
                "submit wrapper missing helper-sourced main help:" \
                "$(print_relpath "${file}")"
        fi

        if ! \
            grep -q -- '--dir_scr' "${fil_hlp}"
        then
            found=1
            record_fail \
                "submit help missing --dir_scr:" \
                "$(print_relpath "${fil_hlp}")"
        fi

        # '--dir_scr' is optional: a worker resolves its own location from
        # 'BASH_SOURCE' when the flag is absent, which is correct for direct
        # invocation, GNU Parallel, serial execution, and 'sbatch --wrap'. Only
        # 'sbatch <script>' needs it, because Slurm runs a copy. The usage
        # block must therefore bracket it.
        if ! \
            awk '
                /^Usage:$/ { in_usage = 1; next }
                /^Usage$/ { pending_usage = 1; next }
                pending_usage && /^-{3,}$/ {
                    in_usage = 1
                    pending_usage = 0
                    next
                }
                pending_usage { pending_usage = 0 }
                in_usage && /^$/ { in_usage = 0 }
                in_usage && /\[[^]]*--dir_scr/ { found = 1 }
                END { exit found ? 0 : 1 }
            ' "${fil_hlp}"
        then
            found=1
            record_fail \
                "submit help must show --dir_scr as optional:" \
                "$(print_relpath "${fil_hlp}")"
        fi

        if ! \
            grep -q -- '-ds|--dir\[_-\]scr' "${file}"
        then
            found=1
            record_fail \
                "submit wrapper missing -ds|--dir_scr parser:" \
                "$(print_relpath "${file}")"
        fi
    done

    if (( found == 0 )); then
        record_pass "no submit --dir_scr policy findings"
    fi
}


# Fail when execute wrappers do not pass '--dir_scr' to submit wrappers.
function scan_execute_submit_dir_scr_policy() {
    local found=0
    local file base line

    for file in "${files_execute[@]}"; do
        base="$(basename "${file}")"

        case "${base}" in
            execute_filter_bams.sh|execute_filter_crams.sh)
                continue
                ;;
        esac

        if ! grep -q 'scr_sub=.*submit_.*\.sh' "${file}"; then
            continue
        fi

        # shellcheck disable=SC2016
        if ! grep -Fq 'BASH_SOURCE[0]' "${file}" \
            || ! grep -Fq 'dir_scr="$(cd "$(dirname' "${file}"
        then
            found=1
            record_fail \
                "execute wrapper missing physical dir_scr default:" \
                "$(print_relpath "${file}")"
        fi

        if ! grep -Eq \
            'scr_sub="\$\{dir_scr\}/submit_[A-Za-z0-9_]+\.sh"' \
            "${file}"
        then
            found=1
            record_fail \
                "execute wrapper submit path not derived from dir_scr:" \
                "$(print_relpath "${file}")"
        fi

        while IFS= read -r line; do
            found=1
            record_fail \
                "execute submit command missing --dir_scr \"\${dir_scr}\":" \
                "$(print_relpath "${file}"):${line}"
        done < <(
            awk '
                /^[[:space:]]*cmd_bld=\(/ {
                    in_cmd = 1
                    start = FNR
                    has_submit = 0
                    has_dir_scr_opt = 0
                    has_dir_scr_val = 0
                }

                in_cmd {
                    if (index($0, "\"${scr_sub}\"")) {
                        has_submit = 1
                    }
                    if (index($0, "--dir_scr")) {
                        has_dir_scr_opt = 1
                    }
                    if (index($0, "\"${dir_scr}\"")) {
                        has_dir_scr_val = 1
                    }
                    if ($0 ~ /^[[:space:]]*\)/) {
                        missing_dir_scr = !(has_dir_scr_opt && has_dir_scr_val)
                        if (has_submit && missing_dir_scr) {
                            print start ":" $0
                        }
                        in_cmd = 0
                    }
                }
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no execute submit --dir_scr policy findings"
    fi
}


print_section "${TEST_NAME}"

# Scan maintained top-level shell CLIs for wrapper structural drift.
files_top=()
while IFS= read -r file; do
    files_top+=( "${file}" )
done < <(
    {
        find "${ROOT_REPO}/bin" -maxdepth 1 -type f -name '*.sh' -print
        find "${ROOT_REPO}/install/scripts" \
            -maxdepth 1 -type f -name '*.sh' -print
    } | sort
)

files_submit=()
while IFS= read -r file; do
    files_submit+=( "${file}" )
done < <(
    find "${ROOT_REPO}/bin" -maxdepth 1 -type f -name 'submit_*.sh' -print \
        | sort
)

files_execute=()
while IFS= read -r file; do
    files_execute+=( "${file}" )
done < <(
    find "${ROOT_REPO}/bin" \
        -maxdepth 1 -type f -name 'execute_*.sh' -print \
        | sort
)

scan_inline_top_level_help
scan_external_main_help
scan_submit_dir_scr_policy
scan_execute_submit_dir_scr_policy

finish
