#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: execute_calculate_scaling_factor.sh
#
# Copyright 2024-2026 by Kris Alavattam
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

#  Set the path to the 'scripts' directory
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
    check_env \
    check_inputs \
    check_numbers \
    format_outputs \
    handle_env \
    help/help_execute_calculate_scaling_factor \
    manage_parallel \
    wrap_cmd \
    || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source required helper scripts." >&2
        exit 1
    }

unset fnc_src


function build_cmd() {
    local idx="${1:-}"
    local csv_mip_i csv_min_i csv_sip_i csv_sin_i
    local len_mip_i len_min_i
    local dep_mip_i dep_min_i dep_sip_i dep_sin_i
    local show_help

    unset cmd_bld && declare -ga cmd_bld

    show_help=$(cat << EOM
Usage:
  build_cmd [-h|--hlp|--help] [idx]

Description:
  Construct the command array 'cmd_bld' for one call to 'submit_calculate_scaling_factor.sh'.

Positional arguments:
  1  idx  <int|UNSET>  Optional zero-based sample index.

                       If omitted or set to 'UNSET', construct a non-indexed command using the current scalar/global argument values.

                       If an index is supplied, construct a per-sample command using indexed values from reconstructed input arrays and broadcastable override arrays.

Expected globals:
  Required scalar globals:
    scr_sub dir_scr env_nam threads mode fil_out dp err_out nam_job aln_typ

  Mode-specific scalar globals:
    'spike' mode:
      method csv_sip csv_sin
    'siq' mode:
      tbl_met cfg_met eqn

  Optional scalar override globals:
    len_def len_mip len_min dep_mip dep_min dep_sip dep_sin

  Required reconstructed arrays for indexed calls:
    arr_mip arr_min

  'spike'-only reconstructed arrays for indexed calls:
    arr_sip arr_sin

  Optional broadcastable override arrays for indexed calls:
    arr_len_mip arr_len_min
    arr_dep_mip arr_dep_min arr_dep_sip arr_dep_sin

Returns:
  0 if 'cmd_bld' is constructed successfully; 1 if argument parsing fails.

Notes:
  - 'cmd_bld' is written as a global indexed array.
  - Broadcastable override arrays may have length 0, 1, or n, where n is the number of samples.
EOM
    )

    if [[ "${idx}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${idx}" ]]; then
        idx="UNSET"
    fi

    if [[ "${idx}" != "UNSET" ]]; then
        check_int_nonneg "${idx}" "idx" || return 1
    fi

    csv_mip_i="${csv_mip}"
    csv_min_i="${csv_min}"
    csv_sip_i="${csv_sip}"
    csv_sin_i="${csv_sin}"

    len_mip_i="${len_mip}"
    len_min_i="${len_min}"
    dep_mip_i="${dep_mip}"
    dep_min_i="${dep_min}"
    dep_sip_i="${dep_sip}"
    dep_sin_i="${dep_sin}"

    if [[ "${idx}" != "UNSET" ]]; then
        csv_mip_i="${arr_mip[idx]}"
        csv_min_i="${arr_min[idx]}"

        if [[ "${mode}" == "spike" ]]; then
            csv_sip_i="${arr_sip[idx]}"
            csv_sin_i="${arr_sin[idx]}"
        fi

        if [[ ${#arr_len_mip[@]} -gt 0 ]]; then
            if [[ ${#arr_len_mip[@]} -eq 1 ]]; then
                len_mip_i="${arr_len_mip[0]}"
            else
                len_mip_i="${arr_len_mip[idx]}"
            fi
        fi

        if [[ ${#arr_len_min[@]} -gt 0 ]]; then
            if [[ ${#arr_len_min[@]} -eq 1 ]]; then
                len_min_i="${arr_len_min[0]}"
            else
                len_min_i="${arr_len_min[idx]}"
            fi
        fi

        if [[ ${#arr_dep_mip[@]} -gt 0 ]]; then
            if [[ ${#arr_dep_mip[@]} -eq 1 ]]; then
                dep_mip_i="${arr_dep_mip[0]}"
            else
                dep_mip_i="${arr_dep_mip[idx]}"
            fi
        fi

        if [[ ${#arr_dep_min[@]} -gt 0 ]]; then
            if [[ ${#arr_dep_min[@]} -eq 1 ]]; then
                dep_min_i="${arr_dep_min[0]}"
            else
                dep_min_i="${arr_dep_min[idx]}"
            fi
        fi

        if [[ ${#arr_dep_sip[@]} -gt 0 ]]; then
            if [[ ${#arr_dep_sip[@]} -eq 1 ]]; then
                dep_sip_i="${arr_dep_sip[0]}"
            else
                dep_sip_i="${arr_dep_sip[idx]}"
            fi
        fi

        if [[ ${#arr_dep_sin[@]} -gt 0 ]]; then
            if [[ ${#arr_dep_sin[@]} -eq 1 ]]; then
                dep_sin_i="${arr_dep_sin[0]}"
            else
                dep_sin_i="${arr_dep_sin[idx]}"
            fi
        fi
    fi

    cmd_bld=(
        "${scr_sub}"
        --env_nam "${env_nam}"
        --dir_scr "${dir_scr}"
        --threads "${threads}"
        --mode "${mode}"
    )

    if [[ "${mode}" == "spike" ]]; then
        cmd_bld+=( --method "${method}" )
    fi

    cmd_bld+=(
        --csv_mip "${csv_mip_i}"
        --csv_min "${csv_min_i}"
    )

    if [[ "${mode}" == "spike" ]]; then
        cmd_bld+=(
            --csv_sip "${csv_sip_i}"
            --csv_sin "${csv_sin_i}"
        )
    fi

    cmd_bld+=(
        --aln_typ "${aln_typ}"
        --fil_out "${fil_out}"
    )

    if [[ "${mode}" == "siq" ]]; then
        cmd_bld+=(
            --tbl_met "${tbl_met}"
            --cfg_met "${cfg_met}"
            --eqn "${eqn}"
        )
    fi

    if [[ -n "${len_def}" ]]; then
        cmd_bld+=( --len_def "${len_def}" )
    fi

    if [[ -n "${len_mip_i}" ]]; then
        cmd_bld+=( --len_mip "${len_mip_i}" )
    fi
    if [[ -n "${len_min_i}" ]]; then
        cmd_bld+=( --len_min "${len_min_i}" )
    fi
    if [[ -n "${dep_mip_i}" ]]; then
        cmd_bld+=( --dep_mip "${dep_mip_i}" )
    fi
    if [[ -n "${dep_min_i}" ]]; then
        cmd_bld+=( --dep_min "${dep_min_i}" )
    fi
    if [[ -n "${dep_sip_i}" ]]; then
        cmd_bld+=( --dep_sip "${dep_sip_i}" )
    fi
    if [[ -n "${dep_sin_i}" ]]; then
        cmd_bld+=( --dep_sin "${dep_sin_i}" )
    fi

    cmd_bld+=(
        --dp "${dp}"
        --err_out "${err_out}"
        --nam_job "${nam_job}"
    )
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
# shellcheck disable=SC2269
{
    env_nam="env_protocol"
    dir_scr="${dir_scr}"
    scr_hdr="${dir_scr}/write_header.sh"
    scr_sub="${dir_scr}/submit_calculate_scaling_factor.sh"
    par_job=""
}
#TODO: 'scr_hdr' is defined and validated, but not directly used in the top-level
#+     script:
#+       - remove 'scr_hdr' from this wrapper if unused or
#+       - implement its use


#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false

threads=1

mode="spike"
method=""

csv_mip=""
csv_min=""
csv_sip=""
csv_sin=""
aln_typ="auto"

fil_out=""

tbl_met=""
cfg_met=""
eqn="6nd"

len_def=""
len_mip=""
len_min=""
dep_mip=""
dep_min=""
dep_sip=""
dep_sin=""

dp=24

err_out=""
nam_job=""
max_job=6
slurm=false
time="0:30:00"


#  Parse arguments
if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
    help_execute_calculate_scaling_factor >&2
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

        -t|--thr|--threads)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            threads="${2}"
            shift 2
            ;;

        -md|--mode)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            mode="${2,,}"
            shift 2
            ;;

        -me|--method)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            method="${2,,}"
            shift 2
            ;;

        -mp|--csv[_-]mip)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            csv_mip="${2}"
            shift 2
            ;;

        -mn|--csv[_-]min)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            csv_min="${2}"
            shift 2
            ;;

        -sp|--csv[_-]sip)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            csv_sip="${2}"
            shift 2
            ;;

        -sn|--csv[_-]sin)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            csv_sin="${2}"
            shift 2
            ;;

        -at|--aln[_-]typ|--align[_-]typ)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            aln_typ="${2,,}"
            shift 2
            ;;

        -fo|--fil[_-]out)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            fil_out="${2}"
            shift 2
            ;;

        -tb|--tbl[_-]met)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            tbl_met="${2}"
            shift 2
            ;;

        -cm|--cfg[_-]met)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            cfg_met="${2}"
            shift 2
            ;;

        -eq|--eqn)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            eqn="${2,,}"
            shift 2
            ;;

        -ld|--len[_-]def)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            len_def="${2}"
            shift 2
            ;;

        -lmp|--len[_-]mip)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            len_mip="${2}"
            shift 2
            ;;

        -lmn|--len[_-]min)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            len_min="${2}"
            shift 2
            ;;

        -dmp|--dep[_-]mip)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            dep_mip="${2}"
            shift 2
            ;;

        -dmn|--dep[_-]min)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            dep_min="${2}"
            shift 2
            ;;

        -dsp|--dep[_-]sip)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            dep_sip="${2}"
            shift 2
            ;;

        -dsn|--dep[_-]sin)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            dep_sin="${2}"
            shift 2
            ;;

        -dp|--dp|--rnd|--round|--decimals|--digits)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            dp="${2}"
            shift 2
            ;;

        -eo|--err[_-]out)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            err_out="${2}"
            shift 2
            ;;

        -nj|--nam[_-]job)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            nam_job="${2}"
            shift 2
            ;;

        -mj|--max[_-]job)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            max_job="${2}"
            shift 2
            ;;

        -sl|--slurm)
            slurm=true
            shift 1
            ;;

        -tm|--time)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_calculate_scaling_factor >&2
                exit 1
            }
            time="${2}"
            shift 2
            ;;

        *)
            echo_err "unknown option/parameter passed: '${1}'."
            echo >&2
            help_execute_calculate_scaling_factor >&2
            exit 1
            ;;
    esac
done

#  Check arguments
validate_var "env_nam" "${env_nam}"
check_env_installed "${env_nam}"

validate_var_dir  "dir_scr" "${dir_scr}" 0 false

validate_var_file "scr_hdr" "${scr_hdr}"

validate_var_file "scr_sub" "${scr_sub}"

validate_var "threads" "${threads}"
check_int_pos "${threads}" "threads"

case "${mode}" in
    siq|alpha)
        #TODO: Phase out 'alpha' as an alias for 'siq'; 'alpha' is ambiguous
        #+     because it is used for both spike-in and siQ-ChIP coefficients
        #+     in the literature.
        mode="siq"

        if [[ -n "${method}" ]]; then
            echo_err \
                "'--method' may be used only when '--mode spike' is active."
            exit 1
        fi
        ;;

    spike|spk)
        mode="spike"

        if [[ -z "${method}" ]]; then
            method="fractional"
        fi

        case "${method}" in
            fractional|bioprotocol|bio_protocol|s)
                #TODO: Phase out hidden test alias 's' unless it remains useful
                #+     for manuscript/blog drafting.
                method="fractional"
                ;;
            chiprx_alpha_ratio|alpha_chiprx_ratio|chiprx_ratio|r)
                #TODO: Decide whether to keep hidden test alias 'r'.
                method="chiprx_alpha_ratio"
                ;;
            chiprx_alpha_ip|alpha_chiprx_ip|chiprx_ip)
                method="chiprx_alpha_ip"
                ;;
            chiprx_alpha_in|alpha_chiprx_in|chiprx_in)
                method="chiprx_alpha_in"
                ;;
            rxinput_alpha|alpha_rxinput|rxi_alpha|alpha_rxi|rxinput|rxi)
                method="rxinput_alpha"
                ;;
            *)
                echo_err \
                    "spike method ('--method') was assigned '${method}' but" \
                    "is not recognized."
                exit 1
                ;;
        esac
        ;;

    *)
        echo_err \
            "scaling-factor mode ('--mode') was assigned '${mode}' but must" \
            "be 'siq' or 'spike'."
        exit 1
        ;;
esac

validate_var "csv_mip" "${csv_mip}"
validate_var_dir "csv_mip parent directory" \
    "$(dirname "${csv_mip%%[,;]*}")" 0 false
check_str_delim "csv_mip" "${csv_mip}"

validate_var "csv_min" "${csv_min}"
validate_var_dir "csv_min parent directory" \
    "$(dirname "${csv_min%%[,;]*}")" 0 false
check_str_delim "csv_min" "${csv_min}"

if [[ "${mode}" == "spike" ]]; then
    validate_var "csv_sip" "${csv_sip}"
    validate_var_dir "csv_sip parent directory" \
        "$(dirname "${csv_sip%%[,;]*}")" 0 false
    check_str_delim "csv_sip" "${csv_sip}"

    validate_var "csv_sin" "${csv_sin}"
    validate_var_dir "csv_sin parent directory" \
        "$(dirname "${csv_sin%%[,;]*}")" 0 false
    check_str_delim "csv_sin" "${csv_sin}"
fi

case "${aln_typ}" in
    pe|paired) aln_typ="pe" ;;
    se|single) aln_typ="se" ;;
    auto) : ;;
    *)
        echo_err \
            "alignment type ('--aln_typ') was assigned '${aln_typ}' but" \
            "must be 'pe', 'se', or 'auto'."
        exit 1
        ;;
esac

validate_var "fil_out" "${fil_out}"

if [[ "${fil_out}" == "-" ]]; then
    echo_err "'-' is not allowed for '--fil_out' in this wrapper."
    exit 1
fi

validate_var_dir "fil_out parent directory" \
    "$(dirname "${fil_out}")"

if [[ "${mode}" == "siq" ]]; then
    validate_var_file "tbl_met" "${tbl_met}"

    validate_var_file "cfg_met" "${cfg_met}"

    validate_var "eqn" "${eqn}"

    case "${eqn}" in
        5|5nd|6|6nd) : ;;
        *)
            echo_err \
                "equation ('--eqn') was assigned '${eqn}' but must be '5'," \
                "'5nd', '6', or '6nd'."
            exit 1
            ;;
    esac

    if [[
             "${aln_typ}" == "se" \
        &&   -z "${len_def}" \
        && ( -z "${len_mip}" || -z "${len_min}" )
    ]]; then
        echo_err \
            "for '--mode siq' with SE data, supply '--len_def' or both" \
            "'--len_mip' and '--len_min'."
        exit 1
    fi

    unset csv_sip csv_sin dep_sip dep_sin
else
    unset tbl_met cfg_met eqn
fi

if [[ -n "${len_def}" ]]; then check_int_pos "${len_def}" "len_def"; fi

validate_var "dp" "${dp}"
check_int_pos "${dp}" "dp"

if [[ -z "${err_out}" ]]; then err_out="$(dirname "${fil_out}")/logs"; fi
validate_var_dir "err_out" "${err_out}"

if [[ -z "${nam_job}" ]]; then
    if [[ "${mode}" == "siq" ]]; then
        nam_job="calc_sf_siq_${eqn}"
    else
        nam_job="calc_sf_spike_${method}"
    fi
fi

validate_var "nam_job" "${nam_job}"


#  Parse and validate input vector elements -----------------------------------
unset arr_mip arr_min && declare -a arr_mip arr_min
IFS=',' read -r -a arr_mip <<< "${csv_mip}"
IFS=',' read -r -a arr_min <<< "${csv_min}"

check_arr_nonempty "arr_mip" "csv_mip"
check_arr_nonempty "arr_min" "csv_min"
check_arr_lengths  "arr_mip" "arr_min"

for file in "${arr_mip[@]}" "${arr_min[@]}"; do
    validate_var_file "file" "${file}"
done
unset file

if [[ "${mode}" == "spike" ]]; then
    unset arr_sip arr_sin && declare -a arr_sip arr_sin

    IFS=',' read -r -a arr_sip <<< "${csv_sip}"
    IFS=',' read -r -a arr_sin <<< "${csv_sin}"

    check_arr_nonempty "arr_sip" "csv_sip"
    check_arr_nonempty "arr_sin" "csv_sin"
    check_arr_lengths  "arr_sip" "arr_sin"
    check_arr_lengths  "arr_sip" "arr_mip"

    for file in "${arr_sip[@]}" "${arr_sin[@]}"; do
        validate_var_file "file" "${file}"
    done
    unset file
fi

#  Parse optional override vectors
unset \
    arr_len_mip arr_len_min \
    arr_dep_mip arr_dep_min arr_dep_sip arr_dep_sin
declare -a \
    arr_len_mip arr_len_min \
    arr_dep_mip arr_dep_min arr_dep_sip arr_dep_sin

if [[ -n "${len_mip}" ]]; then
    IFS=',' read -r -a arr_len_mip <<< "${len_mip}"
fi

if [[ -n "${len_min}" ]]; then
    IFS=',' read -r -a arr_len_min <<< "${len_min}"
fi

if [[ -n "${dep_mip}" ]]; then
    IFS=',' read -r -a arr_dep_mip <<< "${dep_mip}"
fi

if [[ -n "${dep_min}" ]]; then
    IFS=',' read -r -a arr_dep_min <<< "${dep_min}"
fi

if [[ "${mode}" == "spike" && -n "${dep_sip}" ]]; then
    IFS=',' read -r -a arr_dep_sip <<< "${dep_sip}"
fi

if [[ "${mode}" == "spike" && -n "${dep_sin}" ]]; then
    IFS=',' read -r -a arr_dep_sin <<< "${dep_sin}"
fi

check_arr_len_bcst \
    "${#arr_mip[@]}" \
    arr_len_mip arr_len_min arr_dep_mip arr_dep_min

if [[ "${mode}" == "spike" ]]; then
    check_arr_len_bcst \
        "${#arr_mip[@]}" \
        arr_dep_sip arr_dep_sin
fi

#  Check optional override-vector values after broadcast-compatible lengths are
#+ known
if (( ${#arr_len_mip[@]} > 0 )); then
    check_arr_num_pos arr_len_mip len_mip
fi

if (( ${#arr_len_min[@]} > 0 )); then
    check_arr_num_pos arr_len_min len_min
fi

if (( ${#arr_dep_mip[@]} > 0 )); then
    check_arr_int_pos arr_dep_mip dep_mip
fi

if (( ${#arr_dep_min[@]} > 0 )); then
    check_arr_int_pos arr_dep_min dep_min
fi

if [[ "${mode}" == "spike" ]]; then
    if (( ${#arr_dep_sip[@]} > 0 )); then
        check_arr_int_pos arr_dep_sip dep_sip
    fi

    if (( ${#arr_dep_sin[@]} > 0 )); then
        check_arr_int_pos arr_dep_sin dep_sin
    fi
fi


#  Parse job execution parameters ---------------------------------------------
validate_var "max_job" "${max_job}"
check_int_pos "${max_job}" "max_job"

if [[ "${slurm}" == "true" ]]; then
    max_job="$(reset_max_job "${max_job}" "${#arr_mip[@]}")"

    validate_var "time" "${time}"
    check_format_time "${time}"
else
    IFS=';' read -r threads par_job < <(
        set_params_parallel "${threads}" "${max_job}"
    ) || exit 1
    unset max_job time

    validate_var "par_job" "${par_job}"
    check_int_pos "${par_job}" "par_job"
fi

#  Debug parallelization information
if [[ "${mode}" == "spike" ]]; then
    unset extra && typeset -a extra
    extra=( "arr_sip" "arr_sin" )
fi

print_parallel_info \
    "${slurm}" "${max_job:-UNSET}" "${par_job}" "${threads}" \
    "arr_mip" "arr_min" "${extra[@]}"

if [[ "${mode}" == "spike" ]]; then unset extra; fi


#  Activate environment and check that dependencies are in PATH ---------------
env_msg=(
    "'handle_env' failed for 'env_nam=${env_nam}'. Check that Conda/Mamba are"
    "available and that the environment exists."
)

if [[ "${verbose}" == "true" ]]; then
    if out="$(handle_env "${env_nam}")"; then
        print_banner_pretty -tx "${out:-"${env_nam} already active."}" -w "%"
        echo
        echo
    else
        echo_err "${env_msg[*]}"
        exit 1
    fi
else
    if ! handle_env "${env_nam}" > /dev/null 2>&1; then
        echo_err "${env_msg[*]}"
        exit 1
    fi
fi

check_pgrm_path awk
check_pgrm_path python
check_pgrm_path samtools

if [[ "${slurm}" == "true" ]]; then
    check_pgrm_path sbatch
elif [[ "${par_job}" -gt 1 ]]; then
    check_pgrm_path parallel
fi


#  Do the main work ===========================================================
#  Report argument variable assignments if in "verbose mode"
if [[ "${verbose}" == "true" ]]; then
    print_banner_pretty "Hardcoded variable assignments"
    echo
    echo "env_nam=${env_nam}"
    echo "dir_scr=${dir_scr}"
    echo "scr_sub=${scr_sub}"
    echo "scr_hdr=${scr_hdr}"
    echo "par_job=${par_job:-UNSET}"
    echo
    echo

    print_banner_pretty "Argument variable assignments"
    echo
    echo "verbose=${verbose}"
    echo "dry_run=${dry_run}"
    echo "threads=${threads}"
    echo
    echo "mode=${mode}"
    echo "method=${method:-UNSET}"
    echo
    echo "csv_mip=${csv_mip}"
    echo "csv_min=${csv_min}"
    echo "csv_sip=${csv_sip:-UNSET}"
    echo "csv_sin=${csv_sin:-UNSET}"
    echo "aln_typ=${aln_typ}"
    echo "fil_out=${fil_out}"
    echo
    echo "tbl_met=${tbl_met:-UNSET}"
    echo "cfg_met=${cfg_met:-UNSET}"
    echo "eqn=${eqn:-UNSET}"
    echo
    echo "len_def=${len_def:-UNSET}"
    echo "len_mip=${len_mip:-UNSET}"
    echo "len_min=${len_min:-UNSET}"
    echo "dep_mip=${dep_mip:-UNSET}"
    echo "dep_min=${dep_min:-UNSET}"
    echo "dep_sip=${dep_sip:-UNSET}"
    echo "dep_sin=${dep_sin:-UNSET}"
    echo
    echo "dp=${dp}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo "max_job=${max_job:-UNSET}"
    echo "slurm=${slurm}"
    echo "time=${time:-UNSET}"
    echo
    echo

    print_banner_pretty "Arrays derived from variables"
    echo
    echo "arr_mip=( ${arr_mip[*]} )"
    echo
    echo "arr_min=( ${arr_min[*]} )"
    echo

    if [[ "${mode}" == "spike" ]]; then
        echo "arr_sip=( ${arr_sip[*]} )"
        echo
        echo "arr_sin=( ${arr_sin[*]} )"
        echo
    fi

    echo "arr_len_mip=( ${arr_len_mip[*]} )"
    echo
    echo "arr_len_min=( ${arr_len_min[*]} )"
    echo
    echo "arr_dep_mip=( ${arr_dep_mip[*]} )"
    echo
    echo "arr_dep_min=( ${arr_dep_min[*]} )"
    echo

    if [[ "${mode}" == "spike" ]]; then
        echo "arr_dep_sip=( ${arr_dep_sip[*]} )"
        echo
        echo "arr_dep_sin=( ${arr_dep_sin[*]} )"
        echo
    fi

    echo
fi

if [[ "${slurm}" == "true" ]]; then
    #  Slurm execution
    build_cmd "UNSET"

    unset cmd_slurm && declare -a cmd_slurm
    cmd_slurm=(
        sbatch
        --job-name="${nam_job}"
        --nodes=1
        --cpus-per-task="${threads}"
        --time="${time}"
        --output="${err_out}/${nam_job}.%A-%a.stdout.txt"
        --error="${err_out}/${nam_job}.%A-%a.stderr.txt"
        --array="1-${#arr_mip[@]}%${max_job}"
        "${cmd_bld[@]}"
    )

    if [[ "${dry_run}" == "true" || "${verbose}" == "true" ]]; then
        print_banner_pretty "Call to 'sbatch'"
        echo
        printf '%q ' "${cmd_slurm[@]}"
        echo
        echo
    fi

    if [[ "${dry_run}" == "false" ]]; then
        "${cmd_slurm[@]}"
    fi
else
    #  Non-Slurm execution: GNU Parallel ('par_job > 1') or serial
    #+ ('par_job == 1')
    if [[ "${par_job}" -gt 1 ]]; then
        config="${err_out}/${nam_job}.config_parallel.txt"

        if [[ -f "${config}" ]]; then rm "${config}"; fi

        for idx in "${!arr_mip[@]}"; do
            build_cmd "${idx}"
            IFS=';' read -r log_out log_err < <(
                get_submit_logs "${arr_mip[idx]}"
            )

            {
                print_built_cmd "${log_out}" "${log_err}"
            } >> "${config}" || {
                echo_err "failed to write command, index no. '${idx}'."
                exit 1
            }
        done

        if [[ "${dry_run}" == "true" || "${verbose}" == "true" ]]; then
            print_banner_pretty "GNU Parallel execution"
            echo
            parallel --jobs "${par_job}" --dryrun < "${config}"
            echo
            echo
        fi

        if [[ "${dry_run}" == "false" ]]; then
            parallel --jobs "${par_job}" < "${config}"
        fi
    else
        if [[ "${dry_run}" == "true" || "${verbose}" == "true" ]]; then
            print_banner_pretty "Serial execution"
            echo

            for idx in "${!arr_mip[@]}"; do
                build_cmd "${idx}"
                IFS=';' read -r log_out log_err < <(
                    get_submit_logs "${arr_mip[idx]}"
                )
                print_built_cmd "${log_out}" "${log_err}"
                echo
            done

            echo
        fi

        if [[ "${dry_run}" == "false" ]]; then
            for idx in "${!arr_mip[@]}"; do
                build_cmd "${idx}"
                IFS=';' read -r log_out log_err < <(
                    get_submit_logs "${arr_mip[idx]}"
                )
                "${cmd_bld[@]}" >> "${log_out}" 2>> "${log_err}"
            done
        fi
    fi
fi
