#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: combine_parts_scaling_factor.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development and
# documentation.
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


#  Source shared helpers
function source_helpers_script() {
    local fnc_src

    dir_fnc="${dir_scr}/functions"
    fnc_src="${dir_fnc}/source_helpers.sh"

    if [[ ! -f "${fnc_src}" ]]; then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "script not found: '${fnc_src}'." >&2
        return 1
    fi

    # shellcheck disable=SC1090
    source "${fnc_src}" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${fnc_src}'." >&2
        return 1
    }

    source_helpers "${dir_fnc}" \
        check_args \
        check_inputs \
        format_outputs \
        || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper scripts." >&2
            return 1
        }
}


function init_arg_defs() {
    dry_run=false
    force=false
    no_parts=false
    mode=""
    csv_infile=""
    fil_out=""
    num_exp=0
    num_max_siq=0

    unset arr_infile arr_in_ord
    declare -ga arr_infile arr_in_ord
}


function show_help_main() {
    cat << EOM
Usage:
  combine_parts_scaling_factor.sh
    [--help] [--dry_run] [--force] [--no_parts] --mode <enum:siq,spike> --csv_infile <csv:file> --fil_out <file>

Description:
  Combine generated calculate_scaling_factor part files into one data-only tab-delimited output table.

Arguments:
  -h, --hlp, --help  <flag>
    Display this help message and exit.

  -dr, --dry, --dry_run  <flag>
    Validate inputs and print the ordered combination plan without writing output.

  -f, --force  <flag>
    Replace an existing output file.

  -np, --no_parts  <flag>
    Remove validated part files after successfully writing the final table.

  -md, --mode  <enum:siq,spike>
    Scaling-factor part-file mode. Required.

  -i, -fi, -ci, --infile, --infiles, --fil_in, --csv_infile, --csv_infiles  <csv:file>
    Comma-separated list of generated scaling-factor part files. Required.

  -o, -fo, --outfile, --fil_out  <file>
    Final tab-delimited output table. Required.

Dependencies:
  External programs:
    - awk
    - Bash >= 4.4
    - basename
    - cat
    - dirname
    - mktemp
    - mv
    - rm
    - sed
    - sort

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
  - Input basenames must end in '.part.<digits>'.
  - Each input must contain exactly one non-empty core scaling-factor row.
  - Core rows contain 12 or 14 fields for 'siq' and 10 fields for 'spike'.
  - Rows are ordered by numeric part-file index, not by input-list order.
  - This script does not write a header; use 'write_header.sh' separately or let 'execute_calculate_scaling_factor.sh' orchestrate header insertion.
  - Duplicate paths, duplicate indexes, and header rows are rejected.
  - Existing output files are replaced only when '--force' is supplied.
  - Validated part files are retained unless '--no_parts' is supplied.

Examples:
  1. Combine spike-in scaling-factor parts.
    '''bash
    bash combine_parts_scaling_factor.sh
      --mode spike
      --csv_infile results.spike.tsv.part.000001,results.spike.tsv.part.000000
      --fil_out results.tsv
    '''

  2. Preview siQ-ChIP part combination.
    '''bash
    bash combine_parts_scaling_factor.sh
      --dry_run
      --mode siq
      --csv_infile results.siq.tsv.part.000001,results.siq.tsv.part.000000
      --fil_out results.tsv
    '''
EOM
}


function parse_args() {
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -dr|--dry|--dry[_-]run)
                dry_run=true
                shift 1
                ;;

            -f|--force)
                force=true
                shift 1
                ;;

            -np|--no[_-]parts)
                no_parts=true
                shift 1
                ;;

            -md|--mode)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main >&2
                    return 1
                }
                mode="${2}"
                shift 2
                ;;

            -i|-fi|-ci|--infile|--infiles|--fil[_-]in|--csv[_-]infile|--csv[_-]infiles)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main >&2
                    return 1
                }
                csv_infile="${2}"
                shift 2
                ;;

            -o|-fo|--outfile|--fil[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main >&2
                    return 1
                }
                fil_out="${2}"
                shift 2
                ;;

            *)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                show_help_main >&2
                return 1
                ;;
        esac
    done
}


function validate_args() {
    validate_var "mode" "${mode}"
    validate_var "csv_infile" "${csv_infile}"
    validate_var "fil_out" "${fil_out}"
    validate_var_dir "fil_out parent directory" "$(dirname "${fil_out}")"

    case "${mode}" in
        spike)
            num_exp=10
            ;;
        siq)
            num_exp=12
            ;;
        *)
            echo_err "'--mode' must be 'siq' or 'spike', but got '${mode}'."
            return 1
            ;;
    esac

    if [[ -e "${fil_out}" && "${force}" == "false" ]]; then
        echo_err \
            "output file already exists: '${fil_out}'." \
            "Supply '--force' to replace it."
        return 1
    fi
}


function prepare_parts() {
    local base
    local fil_in
    local fil_in_abs
    local fil_out_abs
    local fst_fld
    local idx
    local idx_norm
    local num_fld
    local num_lin
    local -A fil_seen idx_seen
    local -a arr_ord

    IFS=',' read -r -a arr_infile <<< "${csv_infile}"
    unset IFS

    if (( ${#arr_infile[@]} == 0 )); then
        echo_err "'--csv_infile' resolved to an empty vector."
        return 1
    fi

    num_max_siq=0
    fil_out_abs="$(
        cd "$(dirname "${fil_out}")" > /dev/null 2>&1 && pwd
    )/$(basename "${fil_out}")"

    for fil_in in "${arr_infile[@]}"; do
        if [[ -z "${fil_in}" ]]; then
            echo_err "'--csv_infile' contains an empty element."
            return 1
        elif [[ ! -e "${fil_in}" ]]; then
            echo_err "input part file does not exist: '${fil_in}'."
            return 1
        elif [[ ! -r "${fil_in}" ]]; then
            echo_err "input part file is not readable: '${fil_in}'."
            return 1
        elif [[ ! -s "${fil_in}" ]]; then
            echo_err "input part file is empty: '${fil_in}'."
            return 1
        fi

        fil_in_abs="$(cd "$(dirname "${fil_in}")" > /dev/null 2>&1 && pwd)/$(
            basename "${fil_in}"
        )"

        if [[ "${fil_in_abs}" == "${fil_out_abs}" ]]; then
            echo_err "'--fil_out' is also an input part file: '${fil_out}'."
            return 1
        elif [[ -n "${fil_seen[${fil_in_abs}]:-}" ]]; then
            echo_err "duplicate input part file: '${fil_in}'."
            return 1
        fi
        fil_seen["${fil_in_abs}"]=1

        base="$(basename "${fil_in}")"
        if [[ "${base}" =~ \.part\.([0-9]+)$ ]]; then
            idx="${BASH_REMATCH[1]}"
        else
            echo_err \
                "input basename must end in '.part.<digits>':" \
                "'${base}'."
            return 1
        fi

        idx_norm="$(printf '%s\n' "${idx}" | sed 's/^0*//; s/^$/0/')"
        if [[ -n "${idx_seen[${idx_norm}]:-}" ]]; then
            echo_err "duplicate numeric part index '${idx_norm}'."
            return 1
        fi
        idx_seen["${idx_norm}"]=1

        num_lin="$(awk 'END { print NR + 0 }' "${fil_in}")"
        if [[ "${num_lin}" -ne 1 ]]; then
            echo_err \
                "input part file must contain exactly one row:" \
                "'${fil_in}' has '${num_lin}'."
            return 1
        fi

        fst_fld="$(awk -F '\t' 'NR == 1 { print $1 }' "${fil_in}")"
        case "${fst_fld}" in
            fil_ip|main_ip)
                echo_err \
                    "input part file appears to contain a header: '${fil_in}'."
                return 1
                ;;
        esac

        num_fld="$(awk -F '\t' 'NR == 1 { print NF }' "${fil_in}")"
        if [[ "${mode}" == "siq" ]]; then
            if [[ "${num_fld}" -ne 12 && "${num_fld}" -ne 14 ]]; then
                echo_err \
                    "input part file has '${num_fld}' fields; expected 12" \
                    "or 14 for mode '${mode}': '${fil_in}'."
                return 1
            fi
            if (( num_fld > num_max_siq )); then
                num_max_siq="${num_fld}"
            fi
        elif [[ "${num_fld}" -ne "${num_exp}" ]]; then
            echo_err \
                "input part file has '${num_fld}' fields; expected" \
                "'${num_exp}' for mode '${mode}': '${fil_in}'."
            return 1
        fi

        arr_ord+=( "${idx_norm}"$'\t'"${fil_in_abs}" )
    done

    while IFS=$'\t' read -r idx fil_in; do
        arr_in_ord+=( "${fil_in}" )
    done < <(printf '%s\n' "${arr_ord[@]}" | LC_ALL=C sort -n -k1,1)
}


function report_plan() {
    local fil_in

    echo "Ordered scaling-factor part files:"
    for fil_in in "${arr_in_ord[@]}"; do
        echo "  ${fil_in}"
    done

    if [[ "${dry_run}" == "false" ]]; then
        return 0
    fi

    echo
    echo "Dry run: would write combined scaling-factor table '${fil_out}'."

    if [[ "${no_parts}" == "true" ]]; then
        echo "Dry run: would remove validated part files after writing."
    fi
}


function run_jobs() {
    local fil_in
    local tmp=""

    if [[ "${dry_run}" == "true" ]]; then
        return 0
    fi

    tmp="$(mktemp "${fil_out}.tmp.XXXXXX")"
    trap 'rm -f "${tmp:-}"' RETURN

    for fil_in in "${arr_in_ord[@]}"; do
        if [[ "${mode}" == "siq" && "${num_max_siq}" -eq 14 ]]; then
            awk -F '\t' '
                BEGIN { OFS = "\t" }
                NF == 12 { print $0, "NA", "NA"; next }
                { print }
            ' \
                "${fil_in}" >> "${tmp}" || return 1
        else
            cat "${fil_in}" >> "${tmp}" || return 1
        fi
    done

    mv "${tmp}" "${fil_out}" || return 1
    trap - RETURN

    echo
    echo "Wrote combined scaling-factor table '${fil_out}'."

    if [[ "${no_parts}" == "true" ]]; then
        if ! \
            rm -f -- "${arr_in_ord[@]}"
        then
            echo_err \
                "wrote combined scaling-factor table but failed to remove" \
                "one or more validated part files requested by '--no_parts'."
            return 1
        fi

        echo "Removed validated scaling-factor part files."
    fi
}


function main() {
    init_arg_defs
    source_helpers_script || return 1

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        show_help_main >&2
        return 0
    fi

    parse_args "$@" || return 1
    validate_args   || return 1
    prepare_parts   || return 1
    report_plan     || return 1
    run_jobs        || return 1
}


main "$@"
