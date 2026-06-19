#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: write_header.sh
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models) were used in
# development.
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
in_place=false
mode="siq"
fil_in=""
fil_out=""

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  write_header.sh
    [--help] [--verbose] [--dry_run] [--mode <enum:siq,spike>] [--fil_in <file>] [--fil_out <file>] [--in_place]

Description:
  Create a header-only scaling-factor table, write a headered copy of a data table, or add a header to a data table in place.

Arguments:
  -h, --hlp, --help  <flag>
    Display this help message and exit.

  -v, --verbose  <flag>
    Print the header before writing.

  -dr, --dry, --dry_run  <flag>
    Print the header and planned file action without creating or modifying a file.

  -md, --mode  <enum:siq,spike>
    Type of header to write: 'siq' (siQ-ChIP normalization) or 'spike' (normalization with a spike-in coefficient) (default: '${mode}').

  -i, -fi, --infile, --fil_in  <file>
    Input data table to header. If omitted, '--fil_out' creates a header-only utility table.

  -o, -fo, --outfile, --fil_out  <file>
    Output file to create. With '--fil_in', writes a headered copy. Without '--fil_in', writes a header-only table.

  --in_place, --in-place, --inplace  <flag>
    Modify '--fil_in' in place instead of writing a separate output file.

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
  - The script writes the selected core scaling-factor header.
  - 'execute_calculate_scaling_factor.sh' calls this script with '--fil_in <final_table> --in_place' after successful part-file combination unless '--no_header' is supplied there.
  - If '--fil_in' already begins with the expected header, the copy or in-place operation does not duplicate the header.
  - If '--fil_in' is data-only, the selected header is prepended.
  - If '--fil_in' is omitted, '--fil_out' creates a header-only utility table.
  - If '--dry_run' is enabled, the script validates arguments and reports the planned action without creating or modifying files.

Examples:
  1. Create a header-only siQ-ChIP-mode utility table
    '''bash
    bash write_header.sh
      --mode siq
      --fil_out results/header.siq.tsv
    '''

  2. Add a spike-mode header to a data-only table in place
    '''bash
    bash write_header.sh
      --mode spike
      --fil_in results/ChIP_samples_spike.tsv
      --in_place
    '''

  3. Preview writing a headered siQ-ChIP copy
    '''bash
    bash write_header.sh
      --dry_run
      --mode siq
      --fil_in results/ChIP_samples_siq_6nd.data.tsv
      --fil_out results/ChIP_samples_siq_6nd.tsv
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

        -i|-fi|--infile|--fil[_-]in)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                echo "${show_help}" >&2
                exit 1
            }
            fil_in="${2}"
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

        --in[_-]place|--inplace)
            in_place=true
            shift 1
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

if [[ -n "${fil_in}" ]]; then
    validate_var_file "fil_in" "${fil_in}"
fi

if [[ -n "${fil_out}" ]]; then
    validate_var_dir "dir_out" "$(dirname "${fil_out}")"
fi

if [[ -z "${fil_in}" && -n "${fil_out}" && -e "${fil_out}" ]]; then
    echo_err \
        "'--fil_out' already exists in header-only mode: '${fil_out}'. Use" \
        "'--fil_in' with '--fil_out' or '--in_place' to header an existing" \
        "table."
    exit 1
fi

if [[ "${in_place}" == "true" && -z "${fil_in}" ]]; then
    echo_err "'--in_place' requires '--fil_in'."
    exit 1
fi

if [[ -z "${fil_in}" && -z "${fil_out}" ]]; then
    echo_err \
        "either '--fil_out' or '--fil_in' with an output mode is required."
    exit 1
fi

if [[ "${in_place}" == "true" && -n "${fil_out}" ]]; then
    echo_err "use either '--fil_out' or '--in_place', not both."
    exit 1
fi

if [[ -n "${fil_in}" && -z "${fil_out}" && "${in_place}" == "false" ]]; then
    echo_err "'--fil_in' requires either '--fil_out' or '--in_place'."
    exit 1
fi

if [[ -n "${fil_in}" && -n "${fil_out}" ]]; then
    fil_in_abs="$(cd "$(dirname "${fil_in}")" > /dev/null 2>&1 && pwd)/$(
        basename "${fil_in}"
    )"
    fil_out_abs="$(cd "$(dirname "${fil_out}")" > /dev/null 2>&1 && pwd)/$(
        basename "${fil_out}"
    )"

    if [[ "${fil_in_abs}" == "${fil_out_abs}" ]]; then
        echo_err \
            "'--fil_in' and '--fil_out' refer to the same path; use" \
            "'--in_place'."
        exit 1
    fi
fi

if [[ "${in_place}" == "true" ]]; then
    fil_dst="${fil_in}"
elif [[ -n "${fil_out}" ]]; then
    fil_dst="${fil_out}"
else
    fil_dst=""
fi


#  Do the main work ===========================================================
#  Define the header column names as an array
case "${mode}" in
    siq)
        if [[ -n "${fil_in}" ]]; then
            num_fld="$(awk -F '\t' 'NF > 0 { print NF; exit }' "${fil_in}")"
        else
            num_fld=12
        fi

        case "${num_fld}" in
            12)
                nam_col=(
                    "fil_ip" "fil_in" "siq" "eqn"
                    "mass_ip" "mass_in" "vol_all" "vol_in"
                    "dep_ip" "dep_in" "len_ip" "len_in"
                )
                ;;

            14)
                nam_col=(
                    "fil_ip" "fil_in" "siq" "eqn"
                    "mass_ip" "mass_in" "vol_all" "vol_in"
                    "dep_ip" "dep_in" "len_ip" "len_in"
                    "lib_vol_ip" "lib_vol_in"
                )
                ;;

            *)
                echo_err \
                    "cannot choose siQ header for '${fil_in}': first data" \
                    "row has '${num_fld}' fields, but expected 12 or 14."
                exit 1
                ;;
        esac
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
    if [[ -n "${fil_in}" ]]; then
        lin_fst="$(head -n1 "${fil_in}")"

        if [[ "${lin_fst}" == "${header%$'\n'}" ]]; then
            if [[ "${in_place}" == "true" ]]; then
                msg="Dry run: header already present; would not modify"
                msg+=" '${fil_in}'."
            else
                msg="Dry run: header already present; would copy"
                msg+=" '${fil_in}' to '${fil_out}'."
            fi
        else
            if [[ "${in_place}" == "true" ]]; then
                msg="Dry run: would prepend header to '${fil_in}' in place."
            else
                msg="Dry run: would prepend header from '${fil_in}' to"
                msg+=" '${fil_out}'."
            fi
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

#  Prepend the header, copy a headered input, or create a header-only output
if [[ "${dry_run}" == "false" ]]; then
    tmp="${fil_dst}.tmp.$$"
    if [[ -n "${fil_in}" ]]; then
        lin_fst="$(head -n1 "${fil_in}")"
        if [[ "${lin_fst}" == "${header%$'\n'}" ]]; then
            if [[ "${in_place}" == "false" ]]; then
                cat "${fil_in}" > "${tmp}" && mv "${tmp}" "${fil_dst}"
            fi
        else
            { printf '%s\n' "${header%$'\n'}"; cat "${fil_in}"; } > "${tmp}" \
                && mv "${tmp}" "${fil_dst}"
        fi
    else
        printf '%s\n' "${header%$'\n'}" > "${tmp}" && mv "${tmp}" "${fil_dst}"
    fi
fi

unset fmt_str header nam_col
unset fil_dst fil_in_abs fil_out_abs lin_fst msg tmp
