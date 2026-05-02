#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: write_header.sh
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


#  Require Bash >= 4.4 before doing any work
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be run under Bash >= 4.4." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

#  Run in safe mode, exiting on errors, unset variables, and pipe failures
set -euo pipefail

#  Set path to the 'scripts' directory
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"


#  Source and define functions ================================================
# arr_fnc=(  #TODO: record in help documentation before deleting
#     check_file_dir_exists  ## check_inputs ##
#     echo_err               ## format_outputs ##
#     echo_warn              ## format_outputs ##  ## NOTE: not here now ##
# )

dir_fnc="${dir_scr}/functions"
fnc_src="${dir_fnc}/source_helpers.sh"

if [[ ! -f "${fnc_src}" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "script not found: '${fnc_src}'." >&2
    exit 1
fi

# shellcheck disable=SC1090
source "${fnc_src}" || {
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to source '${fnc_src}'." >&2
    exit 1
}

source_helpers "${dir_fnc}" \
    check_args \
    check_inputs \
    format_outputs \
    || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source required helper scripts." >&2
        exit 1
    }


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
mode="alpha"
fil_out=""

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  write_header.sh [--help] [--verbose] [--dry-run] [--mode <str>] --fil_out <str>

Description:
  Write a predefined tab-delimited header to the specified output file.

Options:
   -h, --help     Display this help message and exit.
   -v, --verbose  Print the header before writing.
  -dr, --dry-run  Print the header but do not write to a file.
  -md, --mode     Type of header to write: 'alpha' or 'spike' (default: '${mode}').
  -fo, --fil_out  Output file where the header should be written.

Dependencies:
  - Bash >= 4.4

Example:
  '''bash
  write_header.sh -v -t -md alpha -o results/ChIP_samples_alpha_6nd.tsv
  '''

Note:
  - Script was implemented to pre-write the header before running Slurm jobs, preventing race conditions that occurred when header-writing logic was previously in 'submit' scripts.
EOM
)

#  Parse arguments
if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
    echo "${show_help}" >&2
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case "${1}" in
        -v|--verbose)
            verbose=true
            shift 1
            ;;

        -dr|--dry|--dry[_-]run)
            dry_run=true
            shift 1
            ;;

        -md|--mode)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                echo "${show_help}" >&2
                exit 1
            }
            mode="${2}"
            shift 2
            ;;

        -fo|--fil[_-]out)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                echo "${show_help}" >&2
                exit 1
            }
            fil_out="${2}"
            shift 2
            ;;

        *)
            echo_err "unknown option/parameter passed: '${1}'."
            echo >&2
            echo "${show_help}" >&2
            exit 1
            ;;
    esac
done

#  Check arguments
case "${mode}" in
    alpha|spike) : ;;
    *)
        echo_err \
            "header mode ('--mode') was assigned '${mode}' but must be" \
            "'alpha' or 'spike'."
        exit 1
        ;;
esac

validate_var "fil_out" "${fil_out}"
validate_var_dir "dir_out" "$(dirname "${fil_out}")"


#  Do the main work ===========================================================
#  Define the header column names as an array
case "${mode}" in
    alpha)
        nam_col=(
            "fil_ip" "fil_in" "alpha" "eqn"
            "mass_ip" "mass_in" "vol_all" "vol_in" "dep_ip" "dep_in"
            "len_ip" "len_in"
            "dm_fr_1" "dm_fr_5" "dm_fr_10" "dm_fr_20" "dm_fr_30" "dm_fr_40"
            "dm_fr_50"
            "dm_nm_1" "dm_nm_5" "dm_nm_10" "dm_nm_20" "dm_nm_30" "dm_nm_40"
            "dm_nm_50"
        )
        ;;

    spike) nam_col=(
            "main_ip" "spike_ip" "main_in" "spike_in" "spike"
            "num_mp" "num_sp" "num_mn" "num_sn"
            "dm_fr_1" "dm_fr_5" "dm_fr_10" "dm_fr_20" "dm_fr_30" "dm_fr_40"
            "dm_fr_50"
            "dm_nm_1" "dm_nm_5" "dm_nm_10" "dm_nm_20" "dm_nm_30" "dm_nm_40"
            "dm_nm_50"
        )
        ;;
esac

#  Generate 'printf' format string dynamically
fmt_str=$(printf "%s\t" "${nam_col[@]}")
fmt_str="${fmt_str%$'\t'}\n"  # Remove trailing tab and add newline

#  Print the formatted header line
# shellcheck disable=SC2059
header=$(printf "${fmt_str}" "${nam_col[@]}")

#  Print the header (if in dry-run or verbose modes)
if [[ "${dry_run}" == "true" ]] || [[ "${verbose}" == "true" ]]; then
    echo "##################"
    echo "## Table header ##"
    echo "##################"
    echo
    echo "${header}"
    echo
    echo
fi

#  Prepend the header (only if not already present); create file if missing
if [[ "${dry_run}" == "false" ]]; then
    tmp="${fil_out}.tmp.$$"
    if [[ -f "${fil_out}" ]]; then
        lin_fst="$(head -n1 "${fil_out}")"
        if [[ "${lin_fst}" == "${header%$'\n'}" ]]; then
            :  # Header already present, so do nothing
        else
            { printf '%s\n' "${header%$'\n'}"; cat "${fil_out}"; } > "${tmp}" \
                && mv "${tmp}" "${fil_out}"
        fi
    else
        printf '%s\n' "${header%$'\n'}" > "${tmp}" && mv "${tmp}" "${fil_out}"
    fi
fi
