#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: write_header.sh
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
#   GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


# Require Bash >= 4.4 before doing any work.
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

# Run in safe mode, exiting on errors, unset variables, and pipe failures.
set -euo pipefail

# Resolve 'dir_scr' for 'sbatch', which copies this script elsewhere.
dir_scr=""
_arr_arg=( "$@" )

for (( _idx = 0; _idx < ${#_arr_arg[@]}; _idx++ )); do
    case "${_arr_arg[_idx]}" in
        -ds|--dir[_-]scr)
            if \
                   (( _idx + 1 < ${#_arr_arg[@]} )) \
                && [[ -n "${_arr_arg[_idx + 1]}" ]] \
                && [[ "${_arr_arg[_idx + 1]}" != -* ]]
            then
                dir_scr="${_arr_arg[_idx + 1]}"
            fi

            break
            ;;
    esac
done
unset _arr_arg _idx

if [[ -z "${dir_scr}" ]]; then
    dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
fi

# 'sbatch' runs a copy, so 'BASH_SOURCE' can resolve outside the repo; if so,
# fail here rather than at a subsequent 'source'.
if [[ ! -d "${dir_scr}/../lib/bash" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "cannot locate helper libraries from '${dir_scr}'; pass '--dir_scr'" \
        "when this script is run from a copy." >&2
    exit 1
fi


# Source shared helpers.
function source_helpers_script() {
    local fnc_src

    dir_fnc="${dir_scr}/../lib/bash"
    fnc_src="${dir_fnc}/core/source_helpers.sh"

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
        help/help_write_header \
        || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper scripts." >&2
            return 1
        }
}


function init_arg_defs() {
    verbose=false
    dry_run=false
    in_place=false
    mode="siq"
    fil_in=""
    fil_out=""
    fil_dst=""
    header=""
}


function parse_args() {
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        help_write_header >&2
        return 0
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

            -ds|--dir[_-]scr)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_write_header >&2
                    return 1
                }
                shift 2
                ;;

            -md|--mode)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_write_header >&2
                    return 1
                }
                mode="${2}"
                shift 2
                ;;

            -fi|--fil[_-]in)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_write_header >&2
                    return 1
                }
                fil_in="${2}"
                shift 2
                ;;

            -fo|--fil[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_write_header >&2
                    return 1
                }
                fil_out="${2}"
                shift 2
                ;;

            -ip|--in[_-]place)
                in_place=true
                shift 1
                ;;

            *)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                help_write_header >&2
                return 1
                ;;
        esac
    done
}


function validate_args() {
    local fil_in_abs
    local fil_out_abs

    case "${mode}" in
        siq|spike) : ;;
        *)
            echo_err \
                "header mode ('--mode') was assigned '${mode}' but must be" \
                "'siq' or 'spike'."
            return 1
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
            "'--fil_out' already exists in header-only mode: '${fil_out}'." \
            "Use '--fil_in' with '--fil_out' or '--in_place' to header an" \
            "existing table."
        return 1
    fi

    if [[ "${in_place}" == "true" && -z "${fil_in}" ]]; then
        echo_err "'--in_place' requires '--fil_in'."
        return 1
    fi

    if [[ -z "${fil_in}" && -z "${fil_out}" ]]; then
        echo_err \
            "either '--fil_out' or '--fil_in' with an output mode is required."
        return 1
    fi

    if [[ "${in_place}" == "true" && -n "${fil_out}" ]]; then
        echo_err "use either '--fil_out' or '--in_place', not both."
        return 1
    fi

    if [[ -n "${fil_in}" && -z "${fil_out}" && "${in_place}" == "false" ]]
    then
        echo_err "'--fil_in' requires either '--fil_out' or '--in_place'."
        return 1
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
            return 1
        fi
    fi

    if [[ "${in_place}" == "true" ]]; then
        fil_dst="${fil_in}"
    elif [[ -n "${fil_out}" ]]; then
        fil_dst="${fil_out}"
    else
        fil_dst=""
    fi
}


function build_header() {
    local fmt_str
    local num_fld
    local -a nam_col

    case "${mode}" in
        siq)
            if [[ -n "${fil_in}" ]]; then
                num_fld="$(awk -F '\t' 'NF > 0 { print NF; exit }' \
                    "${fil_in}")"
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
                        "cannot choose siQ header for '${fil_in}': first" \
                        "data row has '${num_fld}' fields, but expected 12" \
                        "or 14."
                    return 1
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

    fmt_str=$(printf "%s\t" "${nam_col[@]}")
    fmt_str="${fmt_str%$'\t'}\n"

    # shellcheck disable=SC2059
    header=$(printf "${fmt_str}" "${nam_col[@]}")
}


function report_plan() {
    local lin_fst
    local msg

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
                    msg="Dry run: would prepend header to '${fil_in}'"
                    msg+=" in place."
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

    if [[ "${dry_run}" == "true" ]] || [[ "${verbose}" == "true" ]]; then
        echo "##################"
        echo "## Table header ##"
        echo "##################"
        echo
        echo "${header}"
        echo
        echo
    fi
}


function run_jobs() {
    local lin_fst
    local tmp

    if [[ "${dry_run}" == "true" ]]; then
        return 0
    fi

    tmp="${fil_dst}.tmp.$$"
    if [[ -n "${fil_in}" ]]; then
        lin_fst="$(head -n1 "${fil_in}")"
        if [[ "${lin_fst}" == "${header%$'\n'}" ]]; then
            if [[ "${in_place}" == "false" ]]; then
                cat "${fil_in}" > "${tmp}" && mv "${tmp}" "${fil_dst}"
            fi
        else
            {
                printf '%s\n' "${header%$'\n'}"
                cat "${fil_in}"
            } > "${tmp}" \
                && mv "${tmp}" "${fil_dst}"
        fi
    else
        printf '%s\n' "${header%$'\n'}" > "${tmp}" \
            && mv "${tmp}" "${fil_dst}"
    fi
}


function main() {
    init_arg_defs
    source_helpers_script || return 1

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        help_write_header >&2
        return 0
    fi

    parse_args "$@" || return 1
    validate_args   || return 1
    build_header    || return 1
    report_plan     || return 1
    run_jobs        || return 1
}


main "$@"
