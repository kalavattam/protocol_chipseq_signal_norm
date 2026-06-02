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
mode="siq"
fil_out=""

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  write_header.sh [--help] [--verbose] [--dry_run] [--mode <enum:siq,spike>] --fil_out <file>


Description:
  Write a predefined tab-delimited header to the specified output file.


Arguments:
  -h, --help  <flag>
    Display this help message and exit.

  -v, --verbose  <flag>
    Print the header before writing.

  -dr, --dry, --dry_run  <flag>
    Print the header and planned file action without creating or modifying a file.

  -md, --mode  <enum:siq,spike>
    Type of header to write: 'siq' or 'spike' (default: '${mode}').

  -o, -fo, --outfile, --fil_out  <file>
    Output file where the header should be written.


Dependencies:
  External programs:
    - Bash >= 4.4
    - cat
    - head
    - mv

  Sourced function scripts:
    - source_helpers.sh
      + source_helpers
    - check_args.sh
      + require_optarg
    - check_inputs.sh
      + validate_var
      + validate_var_dir
    - format_outputs.sh
      + echo_err


Notes:
  - The script writes or prepends the selected core scaling-factor header.
  - 'combine_parts_scaling_factor.sh' calls this script while assembling a final table from deterministic part files.
  - If the output file already begins with the expected header, the file is left unchanged.
  - If the output file exists but does not begin with the expected header, the header is prepended.
  - If the output file does not exist, it is created with the header only.
  - If '--dry_run' is enabled, the script validates arguments and reports the planned action without creating or modifying the output file.


Examples:
  1. Write a siQ-ChIP-mode (siq-mode) header
    '''bash
    bash write_header.sh
      --mode siq
      --fil_out results/ChIP_samples_siq_6nd.tsv
    '''

  2. Preview a spike-mode header without modifying the output file
    '''bash
    bash write_header.sh
      --dry_run
      --mode spike
      --fil_out results/ChIP_samples_spike.tsv
    '''
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

        -o|-fo|--outfile|--fil[_-]out)
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
    siq|spike) : ;;
    *)
        echo_err \
            "header mode ('--mode') was assigned '${mode}' but must be" \
            "'siq' or 'spike'."
        exit 1
        ;;
esac

validate_var "fil_out" "${fil_out}"
validate_var_dir "dir_out" "$(dirname "${fil_out}")"


#  Do the main work ===========================================================
#  Define the header column names as an array
case "${mode}" in
    siq)
        nam_col=(
            "fil_ip" "fil_in" "siq" "eqn"
            "mass_ip" "mass_in" "vol_all" "vol_in" "dep_ip" "dep_in"
            "len_ip" "len_in"
        )
        ;;

    spike)
        nam_col=(
            "main_ip" "spike_ip" "main_in" "spike_in" "spike" "coef"
            "num_mp" "num_sp" "num_mn" "num_sn"
        )
        ;;
esac

#  Generate 'printf' format string dynamically
fmt_str=$(printf "%s\t" "${nam_col[@]}")
fmt_str="${fmt_str%$'\t'}\n"  # Remove trailing tab and add newline

#  Print the formatted header line
# shellcheck disable=SC2059
header=$(printf "${fmt_str}" "${nam_col[@]}")

#  Report the planned file action if in dry-run mode
if [[ "${dry_run}" == "true" ]]; then
    if [[ -f "${fil_out}" ]]; then
        lin_fst="$(head -n1 "${fil_out}")"

        if [[ "${lin_fst}" == "${header%$'\n'}" ]]; then
            msg="Dry run: header already present; would not modify"
            msg+=" '${fil_out}'."
        else
            msg="Dry run: would prepend header to existing file '${fil_out}'."
        fi
    else
        msg="Dry run: would create '${fil_out}' and write the header."
    fi

    echo "${msg}" >&2
    echo >&2
fi

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

unset fmt_str header nam_col
unset lin_fst msg tmp
