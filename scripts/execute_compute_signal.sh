#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: execute_compute_signal.sh
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
    help/help_execute_compute_signal \
    manage_parallel \
    populate_array_empty \
    wrap_cmd \
    || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source required helper scripts." >&2
        exit 1
    }

unset fnc_src


#  Auto-generate prefix for ratio track filenames based on 'method' and
#+ 'scl_fct'
function generate_prefix() {
    local method="${1:-}"
    local scl_fct="${2:-}"
    local pfx
    local scaled="false"
    local show_help

    show_help=$(cat << EOM
Usage:
  generate_prefix [-h|--hlp|--help] [method] [scl_fct]

Description:
  Generate a default output filename prefix for ratio tracks based on the resolved ratio method and whether scaling is in effect.

Positional arguments:
  1  method   <str|UNSET>
    Ratio method name.

    Recognized values:

    | Before        | After            |
    | :------------ | :--------------- |
    | 'log2'        | 'log2_rat'       |
    | 'log2_r'      | 'log2_recip_rat' |
    | 'unadj_r'     | 'recip_rat'      |
    | anything else | 'rat'            |

  2  scl_fct  <str|UNSET>
    Comma-delimited scaling-factor string.

    Scaling is considered active if at least one comma-delimited element is non-empty and not 'NA' after whitespace is removed.

Returns:
  Prints the resolved prefix to stdout.

Notes:
  - If scaling is active, the prefix is prefixed with 'scl_'.
  - If scaling is not active, the base prefix is returned unchanged.
  - On checking whether scaling is in effect:

    | Condition               | Assessment |
    | :---------------------- | :--------- |
    | Empty 'scl_fct'         | No scaling |
    | All 'NA' or blank       | No scaling |
    | Any non-NA or non-blank | Scaled     |
EOM
    )

    case "${method}" in
        log2)    pfx="log2_rat"       ;;
        log2_r)  pfx="log2_recip_rat" ;;
        unadj_r) pfx="recip_rat"      ;;
        *)       pfx="rat"            ;;
    esac

    if [[ -n "${scl_fct}" ]]; then
        local raw="${scl_fct}"
        local arr val

        IFS=',' read -r -a arr <<< "${raw}"
        for val in "${arr[@]}"; do
            #  Strip all whitespace
            val="${val//[[:space:]]/}"
            if [[ -n "${val}" && "${val}" != "NA" ]]; then
                scaled="true"
                break
            fi
        done
    fi

    if [[ "${scaled}" == "true" ]]; then
        echo "scl_${pfx}"
    else
        echo "${pfx}"
    fi
}


function build_cmd() {
    local idx="${1:-}"
    local infile fil_A fil_B outfile
    local scl_fct usr_frg dep_min pseudo
    local show_help

    unset cmd_bld && declare -ga cmd_bld

    show_help=$(cat << EOM
Usage:
  build_cmd [-h|--hlp|--help] [idx]

Description:
  Construct the command array 'cmd_bld' for one call to
  'submit_compute_signal.sh'.

Positional arguments:
  1  idx  <int|UNSET>
    Optional zero-based sample index.

    If omitted or set to 'UNSET', construct a non-indexed command using the current scalar/global argument values.

    If an index is supplied, construct a per-sample command using indexed values from reconstructed input arrays.

Expected globals:
  Required scalar globals:
    scr_sub dir_scr env_nam threads mode rnd err_out nam_job

  Additional scalar globals:
    method siz_bin eps skip_00 drp_nan skp_pfx track

  Required scalar serialized inputs:
    For 'signal' or 'coord':
      csv_infile csv_outfile csv_usr_frg
    For 'signal':
      csv_scl_fct
    For 'ratio':
      csv_fil_A csv_fil_B csv_outfile csv_scl_fct csv_dep_min csv_pseudo

  Required reconstructed arrays for indexed calls:
    For 'signal' or 'coord':
      arr_infile arr_outfile arr_usr_frg
    For 'signal':
      arr_scl_fct
    For 'ratio':
      arr_fil_A arr_fil_B arr_outfile arr_scl_fct arr_dep_min arr_pseudo

Returns:
  0 if 'cmd_bld' is constructed successfully; 1 if argument parsing fails.

Notes:
  - 'cmd_bld' is written as a global indexed array.
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

    #  Assign default local values from global variables or arrays, but only
    #+ for the active mode
    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        infile="${csv_infile}"
        outfile="${csv_outfile}"
        usr_frg="${csv_usr_frg}"

        if [[ "${mode}" == "signal" ]]; then
            scl_fct="${csv_scl_fct}"
        fi
    elif [[ "${mode}" == "ratio" ]]; then
        fil_A="${csv_fil_A}"
        fil_B="${csv_fil_B}"
        outfile="${csv_outfile}"
        scl_fct="${csv_scl_fct}"
        dep_min="${csv_dep_min}"
        pseudo="${csv_pseudo}"
    fi

    if [[ "${idx}" != "UNSET" ]]; then
        if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
            infile="${arr_infile[idx]}"
            outfile="${arr_outfile[idx]}"
            usr_frg="${arr_usr_frg[idx]}"

            if [[ "${mode}" == "signal" ]]; then
                scl_fct="${arr_scl_fct[idx]}"
            fi
        elif [[ "${mode}" == "ratio" ]]; then
            fil_A="${arr_fil_A[idx]}"
            fil_B="${arr_fil_B[idx]}"
            outfile="${arr_outfile[idx]}"
            scl_fct="${arr_scl_fct[idx]}"
            dep_min="${arr_dep_min[idx]}"
            pseudo="${arr_pseudo[idx]}"
        fi
    fi

    cmd_bld=(
        "${scr_sub}"
        --env_nam "${env_nam}"
        --dir_scr "${dir_scr}"
        --threads "${threads}"
        --mode "${mode}"
    )

    if [[ "${mode}" != "coord" ]]; then
        cmd_bld+=( --method "${method}" )
    fi

    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        cmd_bld+=( --csv_infile "${infile}" )
    else
        cmd_bld+=( --csv_fil_A "${fil_A}" --csv_fil_B "${fil_B}" )
    fi

    cmd_bld+=( --csv_outfile "${outfile}" )

    if [[ "${mode}" == "signal" ]]; then
        cmd_bld+=(
            --siz_bin "${siz_bin}"
            --csv_scl_fct "${scl_fct}"
            --csv_usr_frg "${usr_frg}"
        )
    elif [[ "${mode}" == "coord" ]]; then
        cmd_bld+=( --csv_usr_frg "${usr_frg}" )
    else
        cmd_bld+=(
            --csv_scl_fct "${scl_fct}"
            --csv_dep_min "${dep_min}"
            --csv_pseudo "${pseudo}"
        )

        if [[ -n "${eps}" ]]; then
            cmd_bld+=( --eps "${eps}" )
        fi

        if [[ -n "${skip_00}" ]]; then
            cmd_bld+=( --skip_00 "${skip_00}" )
        fi

        if [[ "${drp_nan}" == "true" ]]; then
            cmd_bld+=( --drp_nan )
        fi

        if [[ -n "${skp_pfx}" ]]; then
            cmd_bld+=( --skp_pfx "${skp_pfx}" )
        fi

        if [[ "${track}" == "true" ]]; then
            cmd_bld+=( --track )
        fi
    fi

    cmd_bld+=(
        --rnd "${rnd}"
        --err_out "${err_out}"
        --nam_job "${nam_job}"
    )
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_protocol"
scr_sub="${dir_scr}/submit_compute_signal.sh"
par_job=""

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=4
mode="signal"
method=""
csv_infile=""
csv_fil_A=""
csv_fil_B=""
dir_out=""
typ_out="bedGraph.gz"
prefix=""
track=false
siz_bin=""
csv_scl_fct=""
csv_usr_frg=""
csv_dep_min=""
csv_pseudo=""
eps=""
skip_00=""
drp_nan=false
skp_pfx=""
rnd=24
err_out=""
nam_job=""
max_job=6
slurm=false
time="0:30:00"

#  Parse arguments
if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
    help_execute_compute_signal >&2
    exit 0
else
    if [[ "${1}" =~ ^(-d|--det(ail)?(s)?)$ ]]; then
        detail_execute_compute_signal >&2
        exit 0
    elif [[ "${1}" =~ ^(-ah|--all(_|-)?h[e]?lp)$ ]]; then
        help_execute_compute_signal >&2
        echo >&2
        detail_execute_compute_signal "--no-usage" >&2
        exit 0
    fi
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
                help_execute_compute_signal >&2
                exit 1
            }
            threads="${2}"
            shift 2
            ;;

        -md|--mode)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            mode="${2,,}"
            shift 2
            ;;

        -me|--method)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            method="${2,,}"
            shift 2
            ;;

        -i|-fi|-ci|--infile|--infiles|--fil[_-]in|--csv[_-]infile|--csv[_-]infiles)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            csv_infile="${2}"
            shift 2
            ;;

        -fA|-f1|-cA|-c1|--fil[_-]A|--fil[_-]1|--csv[_-]A|--csv[_-]1|--csv[_-]fil[_-]A)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            csv_fil_A="${2}"
            shift 2
            ;;

        -fB|-f2|-cB|-c2|--fil[_-]B|--fil[_-]2|--csv[_-]B|--csv[_-]2|--csv[_-]fil[_-]B)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            csv_fil_B="${2}"
            shift 2
            ;;

        -do|--dir[_-]out)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            dir_out="${2}"
            shift 2
            ;;

        -to|--typ[_-]out)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            typ_out="${2}"
            shift 2
            ;;

        -px|-pr|--pfx|--prfx|--prefix)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            prefix="${2}"
            shift 2
            ;;

        -tr|--trk|--track)
            track=true
            shift 1
            ;;

        -sb|--siz[_-]bin)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            siz_bin="${2}"
            shift 2
            ;;

        -sf|--scale|--scl[_-]fct|--csv[_-]scl[_-]fct)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            csv_scl_fct="${2}"
            shift 2
            ;;

        -uf|--usr[_-]frg|--csv[_-]usr[_-]frg)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            csv_usr_frg="${2}"
            shift 2
            ;;

        -dm|--dep[_-]min|--depth[_-]min|--csv[_-]dep[_-]min)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            csv_dep_min="${2}"
            shift 2
            ;;

        -ps|--pseudo|--pseudocount|--csv[_-]pseudo)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            csv_pseudo="${2}"
            shift 2
            ;;

        -e|--eps)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            eps="${2}"
            shift 2
            ;;

        -s0|--skp_00|--skip[_-]00)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            skip_00="${2,,}"
            shift 2
            ;;

        -dn|--drp[_-]nan|--drop[_-]nan)
            drp_nan=true
            shift 1
            ;;

        -sk|--skp[_-]pfx|--skip[_-]pfx|--skip[_-]prefix)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            skp_pfx="${2}"
            shift 2
            ;;

        -dp|--dp|--rnd|--round|--decimals|--digits)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            rnd="${2}"
            shift 2
            ;;

        -eo|--err[_-]out)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            err_out="${2}"
            shift 2
            ;;

        -nj|--nam[_-]job|--name[_-]job)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
                exit 1
            }
            nam_job="${2}"
            shift 2
            ;;

        -mj|--max[_-]job)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_compute_signal >&2
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
                help_execute_compute_signal >&2
                exit 1
            }
            time="${2}"
            shift 2
            ;;

        *)
            echo_err "unknown option/parameter passed: '${1}'."
            echo >&2
            help_execute_compute_signal >&2
            exit 1
            ;;
    esac
done


#  Check arguments ------------------------------------------------------------
validate_var "env_nam" "${env_nam}"
check_env_installed "${env_nam}"

validate_var_dir  "dir_scr" "${dir_scr}" 0 false

validate_var_file "scr_sub" "${scr_sub}"

validate_var "threads" "${threads}"
check_int_pos "${threads}" "threads"

case "${mode}" in
    s|sig|signal)
        mode="signal"

        if [[ -z "${method}" ]]; then
            method="norm"
        fi

        case "${method}" in
            u|unadj|unadjusted|s|smp|simple|r|raw)
                method="unadj"
                ;;
            f|frg|frag|frg[_-]len|frag[_-]len|l|len|len[_-]frg|len[_-]frag)
                method="frag"
                ;;
            n|nrm|norm|normalized)
                method="norm"
                ;;
            *)
                echo_err \
                    "invalid value for '--method': '${method}'. Expected" \
                    "'u', 'unadj', 'unadjusted', 's', 'smp', 'simple', 'r'," \
                    "or 'raw' ('method=unadj'); 'f', 'frg', 'frag', 'l'," \
                    "'len', 'len_frg', or 'len_frag' ('method=frag'); 'n'," \
                    "'nrm', 'norm', or 'normalized' ('method=norm')."
                exit 1
                ;;
        esac

        validate_var "csv_infile" "${csv_infile}"
        validate_var_dir "csv_infile parent directory" \
            "$(dirname "${csv_infile%%[,;]*}")" 0 false
        check_str_delim "csv_infile" "${csv_infile}"
        ;;

    r|rat|ratio)
        mode="ratio"

        if [[ -z "${method}" ]]; then
            method="unadj"
        fi

        case "${method}" in
            u|unadj|unadjusted|s|smp|simple|r|raw)
                method="unadj"
                ;;
            2|l2|lg2|log2)
                method="log2"
                ;;
            ur|unadj[_-]r|unadjusted[_-]r|sr|smp[_-]r|simple[_-]r|rr|raw[_-]r)
                method="unadj_r"
                ;;
            2r|l2r|l2[_-]r|lg2[_-]r|log2[_-]r)
                method="log2_r"
                ;;
            *)
                echo_err \
                    "invalid value for '--method': '${method}'. Expected" \
                    "'u', 'unadj', 'unadjusted', 's', 'smp', 'simple', 'r'," \
                    "or 'raw' ('method=unadj'); '2', 'l2', 'lg2', or 'log2'" \
                    " ('method=log2'); 'ur', 'unadj_r', 'unadjusted_r', 'sr'," \
                    "'smp_r', 'simple_r', 'rr', or 'raw_r'" \
                    " ('method=unadj_r'); or '2r', 'l2r', 'l2_r', 'lg2_r'," \
                    " or 'log2_r' ('method=log2_r')."
                exit 1
                ;;

        esac

        validate_var "csv_fil_A" "${csv_fil_A}"
        validate_var_dir "csv_fil_A parent directory" \
            "$(dirname "${csv_fil_A%%[,;]*}")" 0 false
        check_str_delim "csv_fil_A" "${csv_fil_A}"

        validate_var "csv_fil_B" "${csv_fil_B}"
        validate_var_dir "csv_fil_B parent directory" \
            "$(dirname "${csv_fil_B%%[,;]*}")" 0 false
        check_str_delim "csv_fil_B" "${csv_fil_B}"
        ;;

    c|coord|coordinates)
        mode="coord"

        if [[ -n "${method}" ]]; then
            echo_warn \
                "argument '--method' is not applicable with '--mode" \
                "coord'. Ignoring/unsetting 'method=${method}'."
        fi
        unset method

        validate_var "csv_infile" "${csv_infile}"
        validate_var_dir "csv_infile parent directory" \
            "$(dirname "${csv_infile%%[,;]*}")" 0 false
        check_str_delim "csv_infile" "${csv_infile}"
        ;;

    *)
        echo_err \
            "invalid value for '--mode': '${mode}'. Expected 's', 'sig'," \
            "'signal', 'r', 'rat', 'ratio', 'c', 'coord', or 'coordinates'."
        exit 1
        ;;
esac

validate_var_dir "dir_out" "${dir_out}"

if [[ "${mode}" =~ ^(signal|ratio)$ ]]; then
    validate_var "typ_out" "${typ_out}"
    case "${typ_out}" in
        bedGraph|bedGraph.gz|bedgraph|bedgraph.gz|bdg|bdg.gz|bg|bg.gz) : ;;
        *)
            echo_warn \
                "unsupported value for '--typ_out' with '--mode ${mode}':" \
                "'${typ_out}'. Coercing '--typ_out' to 'bedGraph.gz'."
            typ_out="bedGraph.gz"
            ;;
    esac

    if [[ -n "${prefix}" && ! "${prefix}" =~ ^[a-zA-Z0-9._-]+$ ]]; then
        echo_warn \
            "user-supplied '--prefix' contains unusual characters:" \
            "'${prefix}'. Proceeding, but this may result in malformed" \
            "filenames or other issues."
    fi

    if [[ "${mode}" == "signal" ]]; then
        #  If user didn’t supply '--siz_bin', apply a hardcoded default: 10 bp
        if [[ -z "${siz_bin}" ]]; then siz_bin=10; fi

        check_int_pos "${siz_bin}" "siz_bin"
    else
        #  For 'mode=ratio', ignore '--siz_bin' if the user tried to supply it
        if [[ -n "${siz_bin}" ]]; then
            echo_warn \
                "argument '--siz_bin' is not applicable with '--mode ${mode}'." \
                "Ignoring/unsetting '--siz_bin ${siz_bin}'."
            unset siz_bin
        fi
    fi
else
    validate_var "typ_out" "${typ_out}"
    case "${typ_out}" in
        bedGraph|bedGraph.gz|bedgraph|bedgraph.gz|bdg|bdg.gz|bg|bg.gz)
            echo_warn \
                "unsupported value for '--typ_out' with '--mode ${mode}':" \
                "'${typ_out}'. Coercing '--typ_out' to 'bed.gz'."
            typ_out="bed.gz"
            ;;

        bed|bed.gz) : ;;

        *)
            echo_err \
                "invalid value for '--typ_out': '${typ_out}'. Expected" \
                "'bed' or 'bed.gz'."
            exit 1
            ;;
    esac

    #  '--siz_bin' is not applicable with 'mode=coord', so ignore if supplied
    if [[ -n "${siz_bin}" ]]; then
        echo_warn \
            "argument '--siz_bin' is not applicable with '--mode ${mode}'." \
            "Ignoring/unsetting '--siz_bin ${siz_bin}'."
        unset siz_bin
    fi
fi

if [[ -n "${csv_scl_fct}" ]]; then
    check_str_delim "csv_scl_fct" "${csv_scl_fct}"
fi

if [[ -n "${csv_usr_frg}" ]]; then
    check_str_delim "csv_usr_frg" "${csv_usr_frg}"
fi

if [[ -n "${csv_dep_min}" ]]; then
    check_str_delim "csv_dep_min" "${csv_dep_min}"
fi

if [[ -n "${csv_pseudo}" ]]; then
    check_str_delim "csv_pseudo" "${csv_pseudo}"
fi

if [[ "${mode}" == "ratio" ]]; then
    if [[ -n "${eps}" ]]; then check_flt_nonneg "${eps}" "eps"; fi

    if [[ -n "${skip_00}" ]]; then
        case "${skip_00}" in
            pre_scale|post_scale) : ;;
            *)
                echo_err \
                    "invalid value for '--skip_00': '${skip_00}'." \
                    "Expected 'pre_scale' or 'post_scale'."
                exit 1
                ;;
        esac
    fi

    if [[ -n "${csv_dep_min}" && -n "${csv_pseudo}" ]]; then
        echo_warn \
            "both '--csv_dep_min' and '--csv_pseudo' were supplied for" \
            "'--mode ratio'. This is allowed, but interpretability may be" \
            "reduced because both arguments stabilize low-depth ratio" \
            "behavior in different ways."
    fi
fi

validate_var "rnd" "${rnd}"
check_int_pos "${rnd}" "rnd"

if [[ -z "${err_out}" ]]; then err_out="${dir_out}/logs"; fi
validate_var_dir "err_out" "${err_out}"

if [[ -z "${nam_job}" ]]; then
    if [[ "${mode}" == "coord" ]]; then
        nam_job="compute_${mode}"
    else
        nam_job="compute_${mode}_${method}"
    fi
fi
validate_var "nam_job" "${nam_job}"


#  Parse and validate input vector elements -----------------------------------
if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
    if [[ -n "${csv_infile}" ]]; then
        IFS=',' read -r -a arr_infile <<< "${csv_infile}"
        check_arr_nonempty "arr_infile" "csv_infile"

        for infile in "${arr_infile[@]}"; do
            validate_var_file "infile" "${infile}"
        done
        unset infile
    else
        echo_err "variable 'csv_infile' unexpectedly empty."
        exit 1
    fi

    unset arr_outfile && declare -a arr_outfile
    for i in "${arr_infile[@]}"; do
        base="$(basename "${i}" .bam)"

        if [[ -n "${prefix}" ]]; then
            if [[ "${base}" =~ ^IP_ ]]; then
                echo_warn \
                    "filename '${base}' starts with 'IP_'. Stripping this" \
                    "prefix before applying custom '--prefix ${prefix}'."
                base="${base#IP_}"
            elif [[ "${base}" =~ ^in_ ]]; then
                echo_warn \
                    "filename '${base}' starts with 'in_'. Stripping this" \
                    "prefix before applying custom '--prefix ${prefix}'."
                base="${base#in_}"
            fi
            arr_outfile+=( "${dir_out}/${prefix}.${base}.${typ_out}" )
        else
            arr_outfile+=( "${dir_out}/${base}.${typ_out}" )
        fi
    done

    check_arr_lengths "arr_outfile" "arr_infile"

    #  Scaling factors are only used for '--mode signal'
    if [[ "${mode}" == "signal" ]]; then
        if [[ -z "${csv_scl_fct}" ]]; then
            unset arr_scl_fct && declare -a arr_scl_fct
            populate_array_empty arr_scl_fct "${#arr_infile[@]}"
        else
            IFS=',' read -r -a arr_scl_fct <<< "${csv_scl_fct}"
        fi

        for s in "${arr_scl_fct[@]}"; do
            if [[ "${s}" != "NA" ]]; then check_flt_pos "${s}" "csv_scl_fct"; fi
        done
        unset s

        check_arr_lengths "arr_scl_fct" "arr_infile"
    fi

    #  User-supplied fragment lengths are allowed for both 'signal' and 'coord'
    if [[ -z "${csv_usr_frg}" ]]; then
        unset arr_usr_frg && declare -a arr_usr_frg
        populate_array_empty arr_usr_frg "${#arr_infile[@]}"
    else
        IFS=',' read -r -a arr_usr_frg <<< "${csv_usr_frg}"
    fi

    for u in "${arr_usr_frg[@]}"; do
        if [[ "${u}" != "NA" ]]; then check_flt_pos "${u}" "csv_usr_frg"; fi
    done
    unset u

    check_arr_lengths "arr_usr_frg" "arr_infile"
else
    if [[ -n "${csv_fil_A}" ]]; then
        IFS=',' read -r -a arr_fil_A <<< "${csv_fil_A}"
        check_arr_nonempty "arr_fil_A" "csv_fil_A"

        for file in "${arr_fil_A[@]}"; do
            validate_var_file "csv_fil_A" "${file}"
        done
        unset file
    else
        echo_err "variable 'csv_fil_A' unexpectedly empty."
        exit 1
    fi

    if [[ -n "${csv_fil_B}" ]]; then
        IFS=',' read -r -a arr_fil_B <<< "${csv_fil_B}"
        check_arr_nonempty "arr_fil_B" "csv_fil_B"

        for file in "${arr_fil_B[@]}"; do
            validate_var_file "csv_fil_B" "${file}"
        done
        unset file
    else
        echo_err "variable 'csv_fil_B' unexpectedly empty."
        exit 1
    fi

    pfx_lcl="${prefix}"
    exts=( bedGraph bedGraph.gz bedgraph bedgraph.gz bdg bdg.gz bg bg.gz )
    unset arr_outfile && declare -a arr_outfile
    for i in "${arr_fil_A[@]}"; do
        base=$(basename "${i}")

        #  Strip 'IP_' prefix (if present)
        base="${base#IP_}"

        #  If not user-assigned, determine filename 'prefix' based on provided
        #+ arguments
        if [[ -z "${pfx_lcl}" ]]; then
            pfx_lcl="$(generate_prefix "${method}" "${csv_scl_fct}")"
        fi

        #  Remove file extensions
        for ext in "${exts[@]}"; do
            base="${base%."${ext}"}"
        done
        unset ext

        arr_outfile+=( "${dir_out}/${pfx_lcl}_${base}.${typ_out}" )
    done
    unset i base pfx_lcl exts

    if [[ -z "${csv_scl_fct}" ]]; then
        #  Produce per-sample sentinel entries ("NA")
        unset arr_scl_fct && declare -a arr_scl_fct
        populate_array_empty arr_scl_fct "${#arr_fil_A[@]}"
    else
        IFS=',' read -r -a arr_scl_fct <<< "${csv_scl_fct}"
    fi

    for s in "${arr_scl_fct[@]}"; do
        if [[ "${s}" != "NA" ]]; then check_flt_pos "${s}" "csv_scl_fct"; fi
    done
    unset s

    if [[ -z "${csv_dep_min}" ]]; then
        #  Produce per-sample sentinel entries ("NA")
        unset arr_dep_min && declare -a arr_dep_min
        populate_array_empty arr_dep_min "${#arr_fil_A[@]}"
    else
        IFS=',' read -r -a arr_dep_min <<< "${csv_dep_min}"
    fi

    if [[ -z "${csv_pseudo}" ]]; then
        #  Produce per-sample sentinel entries ("NA")
        unset arr_pseudo && declare -a arr_pseudo
        populate_array_empty arr_pseudo "${#arr_fil_A[@]}"
    else
        IFS=',' read -r -a arr_pseudo <<< "${csv_pseudo}"
    fi

    for p in "${arr_pseudo[@]}"; do
        if [[ "${p}" == "NA" ]]; then continue; fi

        if [[ "${p}" != *:* ]]; then
            check_flt_nonneg "${p}" "csv_pseudo"
            continue
        fi

        if [[ "${p}" == *:*:* ]]; then
            echo_err \
                "invalid pseudocount spec in '--csv_pseudo': '${p}'." \
                "Expected 'A' or 'A:B'."
            exit 1
        fi

        IFS=':' read -r pseudo_A pseudo_B <<< "${p}"

        if [[ -z "${pseudo_A}" || -z "${pseudo_B}" ]]; then
            echo_err \
                "invalid pseudocount spec in '--csv_pseudo': '${p}'." \
                "Expected 'A' or 'A:B'."
            exit 1
        fi

        check_flt_nonneg "${pseudo_A}" "csv_pseudo"
        check_flt_nonneg "${pseudo_B}" "csv_pseudo"
    done
    unset p pseudo_A pseudo_B

    for d in "${arr_dep_min[@]}"; do
        if [[ "${d}" != "NA" ]]; then check_flt_pos "${d}" "csv_dep_min"; fi
    done
    unset d

    check_arr_lengths "arr_fil_B"   "arr_fil_A"
    check_arr_lengths "arr_outfile" "arr_fil_A"
    check_arr_lengths "arr_scl_fct" "arr_fil_A"
    check_arr_lengths "arr_dep_min" "arr_fil_A"
    check_arr_lengths "arr_pseudo"  "arr_fil_A"
fi


#  Parse job execution parameters ---------------------------------------------
validate_var "max_job" "${max_job}"
check_int_pos "${max_job}" "max_job"

if [[ "${slurm}" == "true" ]]; then
    if [[ "${mode}" == "ratio" ]]; then
        max_job=$(reset_max_job "${max_job}" "${#arr_outfile[@]}")
    else
        max_job=$(reset_max_job "${max_job}" "${#arr_infile[@]}")
    fi

    validate_var "time" "${time}"
    check_format_time "${time}"
else
    IFS=';' read -r threads par_job < <(
        set_params_parallel "${threads}" "${max_job}"
    )
    unset time

    validate_var "par_job" "${par_job}"
    check_int_pos "${par_job}" "par_job"
fi

#  Debug parallelization information and summary output of resolved states
print_parallel_info \
    "${slurm}" "${max_job:-UNSET}" "${par_job:-UNSET}" "${threads}" \
    "arr_infile" "arr_fil_A" "arr_fil_B" "arr_outfile" \
    "arr_scl_fct" "arr_usr_frg" "arr_dep_min" "arr_pseudo"

if [[ "${mode}" == "signal" ]]; then
    summarize_sig_norm "${method}" "${csv_scl_fct}"
fi


#  Activate environment and check that dependencies are in PATH ---------------
#TODO: this env activation block is repeated verbatim many times: modularize?
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

check_pgrm_path python

if [[ "${slurm}" == "true" ]]; then
    check_pgrm_path sbatch
elif [[ "${par_job}" -gt 1 ]]; then
    check_pgrm_path parallel
fi


#  Do the main work ===========================================================
if [[ "${verbose}" == "true" ]]; then
    print_banner_pretty "Hardcoded variable assignments"
    echo
    echo "env_nam=${env_nam}"
    echo "scr_sub=${scr_sub}"
    echo "par_job=${par_job:-UNSET}"
    echo
    echo
    print_banner_pretty "Argument variable assignments"
    echo
    echo "verbose=${verbose}"
    echo "dry_run=${dry_run}"
    echo "threads=${threads}"
    echo "mode=${mode}"
    echo "method=${method:-UNSET}"
    echo "csv_infile=${csv_infile:-UNSET}"
    echo "csv_fil_A=${csv_fil_A:-UNSET}"
    echo "csv_fil_B=${csv_fil_B:-UNSET}"
    echo "dir_out=${dir_out}"
    echo "typ_out=${typ_out}"
    echo "prefix=${prefix:-UNSET}"
    echo "track=${track}"
    echo "siz_bin=${siz_bin:-UNSET}"
    echo "csv_scl_fct=${csv_scl_fct:-UNSET}"
    echo "csv_usr_frg=${csv_usr_frg:-UNSET}"
    echo "csv_dep_min=${csv_dep_min:-UNSET}"
    echo "csv_pseudo=${csv_pseudo:-UNSET}"
    echo "eps=${eps:-UNSET}"
    echo "skip_00=${skip_00:-UNSET}"
    echo "drp_nan=${drp_nan}"
    echo "skp_pfx=${skp_pfx:-UNSET}"
    echo "rnd=${rnd}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo "max_job=${max_job:-UNSET}"
    echo "slurm=${slurm}"
    echo "time=${time:-UNSET}"
    echo
    echo
    print_banner_pretty "Arrays derived from variables"
    echo

    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        echo "arr_infile=( ${arr_infile[*]} )"
        echo
    elif [[ "${mode}" == "ratio" ]]; then
        echo "arr_fil_A=( ${arr_fil_A[*]} )"
        echo
        echo "arr_fil_B=( ${arr_fil_B[*]} )"
        echo
    fi

    echo "arr_outfile=( ${arr_outfile[*]} )"
    echo

    if [[ "${mode}" == "signal" ]]; then
        echo "arr_scl_fct=( ${arr_scl_fct[*]} )"
        echo
        echo "arr_usr_frg=( ${arr_usr_frg[*]} )"
        echo
    elif [[ "${mode}" == "coord" ]]; then
        echo "arr_usr_frg=( ${arr_usr_frg[*]} )"
        echo
    elif [[ "${mode}" == "ratio" ]]; then
        echo "arr_scl_fct=( ${arr_scl_fct[*]} )"
        echo
        echo "arr_dep_min=( ${arr_dep_min[*]} )"
        echo
        echo "arr_pseudo=( ${arr_pseudo[*]} )"
        echo
    fi

    echo
fi

if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
    csv_infile=$(echo "${arr_infile[*]}"  | tr ' ' ',')
    csv_outfile=$(echo "${arr_outfile[*]}" | tr ' ' ',')
    csv_usr_frg=$(echo "${arr_usr_frg[*]}" | tr ' ' ',')

    if [[ "${mode}" == "signal" ]]; then
        csv_scl_fct=$(echo "${arr_scl_fct[*]}" | tr ' ' ',')
    fi
elif [[ "${mode}" == "ratio" ]]; then
    csv_fil_A=$(echo "${arr_fil_A[*]}"  | tr ' ' ',')
    csv_fil_B=$(echo "${arr_fil_B[*]}"  | tr ' ' ',')
    csv_outfile=$(echo "${arr_outfile[*]}" | tr ' ' ',')
    csv_scl_fct=$(echo "${arr_scl_fct[*]}" | tr ' ' ',')
    csv_dep_min=$(echo "${arr_dep_min[*]}" | tr ' ' ',')
    csv_pseudo=$(echo "${arr_pseudo[*]}" | tr ' ' ',')
fi

if [[ "${verbose}" == "true" ]]; then
    print_banner_pretty "Variable assignments constructed from arrays"
    echo

    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        echo "csv_infile=\"${csv_infile}\""
        echo
    elif [[ "${mode}" == "ratio" ]]; then
        echo "csv_fil_A=\"${csv_fil_A}\""
        echo
        echo "csv_fil_B=\"${csv_fil_B}\""
        echo
    fi

    echo "csv_outfile=\"${csv_outfile}\""
    echo

    if [[ "${mode}" != "coord" ]]; then
        echo "csv_scl_fct=\"${csv_scl_fct}\""
        echo
    fi

    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        echo "csv_usr_frg=\"${csv_usr_frg}\""
        echo
    elif [[ "${mode}" == "ratio" ]]; then
        echo "csv_dep_min=\"${csv_dep_min}\""
        echo
        echo "csv_pseudo=\"${csv_pseudo}\""
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
        --array="1-${#arr_outfile[@]}%${max_job}"
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

        for idx in "${!arr_outfile[@]}"; do
            build_cmd "${idx}"
            IFS=';' read -r log_out log_err < <(
                get_submit_logs "${arr_outfile[idx]}"
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

            for idx in "${!arr_outfile[@]}"; do
                build_cmd "${idx}"
                IFS=';' read -r log_out log_err < <(
                    get_submit_logs "${arr_outfile[idx]}"
                )
                print_built_cmd "${log_out}" "${log_err}"
                echo
            done

            echo
        fi

        if [[ "${dry_run}" == "false" ]]; then
            for idx in "${!arr_outfile[@]}"; do
                build_cmd "${idx}"
                IFS=';' read -r log_out log_err < <(
                    get_submit_logs "${arr_outfile[idx]}"
                )
                "${cmd_bld[@]}" >> "${log_out}" 2>> "${log_err}"
            done
        fi
    fi
fi
