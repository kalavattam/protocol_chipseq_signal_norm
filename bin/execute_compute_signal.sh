#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: execute_compute_signal.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in development and documentation.
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


#  Source shared helpers
function source_helpers_execute() {
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
            return 1
        }
}


function generate_pfx() {
    local method="${1:-}"
    local scl_fct="${2:-}"
    local pfx
    local scaled="false"
    local show_help

    show_help=$(cat << EOM
Usage
-----
  generate_pfx
    [--help] [method] [scl_fct]

  Generate a default output filename prefix for ratio tracks based on the resolved ratio method and whether scaling is in effect.

Parameters
----------
  1  method : str
    Workflow method for ratio-track prefix generation.

    Recognized values:

    | Before        | After            |
    | :----         | :----            |
    | 'log2'        | 'log2_rat'       |
    | 'log2_r'      | 'log2_recip_rat' |
    | 'unadj_r'     | 'recip_rat'      |
    | anything else | 'rat'            |

  2  scl_fct : list of structured string
    Scaling factor. Comma-delimited scaling-factor string.

    Scaling is considered active if at least one comma-delimited element is non-empty and not 'NA' after whitespace is removed.

Returns
-------
  Prints the resolved prefix to stdout.

Notes
-----
  Runtime requirements:
    bash >= 4.4

  - If scaling is active, the prefix is prefixed with 'scl_'.
  - If scaling is not active, the base prefix is returned unchanged.
  - On checking whether scaling is in effect:

    | Condition               | Assessment |
    | :----                   | :----      |
    | Empty 'scl_fct'         | No scaling |
    | All 'NA' or blank       | No scaling |
    | Any non-NA or non-blank | Scaled     |

Examples
--------
  1. Generate an unscaled logarithmic-ratio prefix.
    '''bash
    generate_pfx log2 NA
    '''

  2. Generate a scaled reciprocal-ratio prefix.
    '''bash
    generate_pfx unadj_r 2:1
    '''
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


function check_scl_fct_ratio() {
    local scl_fct="${1:-}"
    local scl_A scl_B
    local show_help

    show_help=$(cat << EOM
Usage
-----
  check_scl_fct_ratio
    [--help] scl_fct

  Validate one ratio-mode scaling-factor element.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  scl_fct : structured string
    Scaling factor. Either 'NA', a positive scalar float, or a positive 'A:B' scaling-factor spec.

Returns
-------
  0 if the scaling-factor element is valid; 1 otherwise.

Notes
-----
  Runtime requirements:
    bash >= 4.4

  - Scalar 'A' means scale file A by A and file B by 1.0.
  - Spec 'A:B' means scale file A by A and file B by B.

Examples
--------
  1. Accept the sentinel that disables ratio scaling.
    '''bash
    check_scl_fct_ratio NA
    '''

  2. Accept distinct positive scaling factors for files A and B.
    '''bash
    check_scl_fct_ratio 2:1
    '''
EOM
    )

    if [[ "${scl_fct}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${scl_fct}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "missing required scaling-factor element."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if [[ "${scl_fct}" == "NA" ]]; then
        return 0
    elif [[ "${scl_fct}" != *:* ]]; then
        check_flt_pos "${scl_fct}" "csv_scl_fct" || return 1
        return 0
    elif [[ "${scl_fct}" == *:*:* ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "invalid scaling-factor spec in '--csv_scl_fct': '${scl_fct}'." \
            "Expected 'A' or 'A:B'."
        return 1
    fi

    IFS=':' read -r scl_A scl_B <<< "${scl_fct}"

    if [[ -z "${scl_A}" || -z "${scl_B}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "invalid scaling-factor spec in '--csv_scl_fct': '${scl_fct}'." \
            "Expected 'A' or 'A:B'."
        return 1
    fi

    check_flt_pos "${scl_A}" "csv_scl_fct" || return 1
    check_flt_pos "${scl_B}" "csv_scl_fct" || return 1
}


function build_cmd() {
    local idx="${1:-}"
    local fil_in fil_A fil_B fil_out
    local scl_fct usr_frg dep_min pseudo
    local show_help

    unset cmd_bld && declare -ga cmd_bld

    show_help=$(cat << EOM
Usage
-----
  build_cmd
    [--help] [idx]

  Construct the command array 'cmd_bld' for one call to
  'submit_compute_signal.sh'.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  idx : int
    Optional zero-based sample index.

    If omitted or set to 'UNSET', construct a non-indexed command using the current scalar/global argument values.

    If an index is supplied, construct a per-sample command using indexed values from reconstructed input arrays.

Expected globals
----------------
  scr_sub, ref_fa : file
    Submission-script and optional reference-FASTA paths, respectively.

  dir_scr, dir_eo : dir
    Script and wrapper-log directories, respectively.

  env_nam, mode, method, engine, skip_00, skp_pfx, nam_job : str
    Environment, mode, method, engine, zero-bin policy, skipped prefix list, and job-name values, respectively.

  threads, siz_bin, chunk_size, dp : int
    Thread count, bin size, chunk size, and rounding precision, respectively.

  eps : num
    Optional numerical tolerance.

  drp_nan, track : bool
    Missing-value and track-line switches, respectively.

  csv_fil_in, csv_fil_out, csv_usr_frg, csv_scl_fct, csv_fil_A, csv_fil_B, csv_dep_min, csv_pseudo : str
    Mode-dependent serialized input, output, fragment-length, scaling-factor, numerator, denominator, minimum-depth, and pseudocount values, respectively.

  arr_fil_in, arr_fil_out, arr_usr_frg, arr_scl_fct, arr_fil_A, arr_fil_B, arr_dep_min, arr_pseudo : array
    Reconstructed input, output, fragment-length, scaling-factor, numerator, denominator, minimum-depth, and pseudocount arrays, respectively.

Returns
-------
  0 if 'cmd_bld' is constructed successfully; 1 if argument parsing fails.

Notes
-----
  Runtime requirements:
    bash >= 4.4

  'cmd_bld' is written as a global indexed array.

Examples
--------
  1. Build the non-indexed command for the active validated mode.
    '''bash
    build_cmd
    declare -p cmd_bld
    '''

  2. Build the first per-sample command from prepared mode-specific arrays.
    '''bash
    build_cmd 0
    declare -p cmd_bld
    '''
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
        fil_in="${csv_fil_in}"
        fil_out="${csv_fil_out}"
        usr_frg="${csv_usr_frg}"

        if [[ "${mode}" == "signal" ]]; then
            scl_fct="${csv_scl_fct}"
        fi
    elif [[ "${mode}" == "ratio" ]]; then
        fil_A="${csv_fil_A}"
        fil_B="${csv_fil_B}"
        fil_out="${csv_fil_out}"
        scl_fct="${csv_scl_fct}"
        dep_min="${csv_dep_min}"
        pseudo="${csv_pseudo}"
    fi

    if [[ "${idx}" != "UNSET" ]]; then
        if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
            fil_in="${arr_fil_in[idx]}"
            fil_out="${arr_fil_out[idx]}"
            usr_frg="${arr_usr_frg[idx]}"

            if [[ "${mode}" == "signal" ]]; then
                scl_fct="${arr_scl_fct[idx]}"
            fi
        elif [[ "${mode}" == "ratio" ]]; then
            fil_A="${arr_fil_A[idx]}"
            fil_B="${arr_fil_B[idx]}"
            fil_out="${arr_fil_out[idx]}"
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
        cmd_bld+=( --csv_fil_in "${fil_in}" )

        if [[ -n "${ref_fa}" ]]; then
            cmd_bld+=( --ref_fa "${ref_fa}" )
        fi
    else
        cmd_bld+=( --csv_fil_A "${fil_A}" --csv_fil_B "${fil_B}" )
    fi

    if [[ -n "${chr_sizes}" ]]; then
        cmd_bld+=( --chr_sizes "${chr_sizes}" )
    fi

    cmd_bld+=( --csv_fil_out "${fil_out}" )

    if [[ "${mode}" == "signal" ]]; then
        cmd_bld+=(
            --siz_bin "${siz_bin}"
            --engine "${engine}"
            --chunk_size "${chunk_size}"
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

        if [[ "${strict_bins}" == "true" ]]; then
            cmd_bld+=( --strict_bins )
        fi
    fi

    cmd_bld+=(
        --dp "${dp}"
        --dir_eo "${dir_eo}"
        --nam_job "${nam_job}"
    )
}


#  Initialize hardcoded argument variables
function init_args_hardcoded() {
    env_nam="env_protocol"
    scr_sub="${dir_scr}/submit_compute_signal.sh"
    par_job=""
}


#  Initialize argument variables, assigning default values where applicable
function init_arg_defs() {
    verbose=false
    dry_run=false
    threads=4
    mode="signal"
    method=""
    csv_fil_in=""
    ref_fa=""
    chr_sizes=""
    csv_fil_A=""
    csv_fil_B=""
    dir_out=""
    typ_out="bedGraph.gz"
    prefix=""
    track=false
    siz_bin=""
    engine="chrom"
    chunk_size=100000
    csv_scl_fct=""
    csv_usr_frg=""
    csv_dep_min=""
    csv_pseudo=""
    eps=""
    skip_00=""
    strict_bins=false
    drp_nan=false
    skp_pfx=""
    dp=24
    dir_eo=""
    nam_job=""
    max_job=6
    slurm=false
    time="0:30:00"
}


#  Initialize hardcoded arguments and user-facing argument defaults
function init_defs() {
    init_args_hardcoded
    init_arg_defs
}


#  Parse keyword arguments
function parse_args() {
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

            -en|--env|--env[_-]nam)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                env_nam="${2}"
                shift 2
                ;;

            -t|--thr|--threads)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                threads="${2}"
                shift 2
                ;;

            -md|--mode)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                mode="${2,,}"
                shift 2
                ;;

            -me|--method)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                method="${2,,}"
                shift 2
                ;;

            -ci|--csv[_-]fil_in)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                csv_fil_in="${2}"
                shift 2
                ;;

            -rf|--ref[_-]fa)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                ref_fa="${2}"
                shift 2
                ;;

            -cs|--chr[_-]sizes|--chrom[_-]sizes)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                chr_sizes="${2}"
                shift 2
                ;;

            -cA|--csv[_-]fil[_-]A)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                csv_fil_A="${2}"
                shift 2
                ;;

            -cB|--csv[_-]fil[_-]B)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                csv_fil_B="${2}"
                shift 2
                ;;

            -do|--dir[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                dir_out="${2}"
                shift 2
                ;;

            -to|--typ[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                typ_out="${2}"
                shift 2
                ;;

            -px|--pfx|--prefix)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
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
                    return 1
                }
                siz_bin="${2}"
                shift 2
                ;;

            -eg|--engine)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                engine="${2,,}"
                shift 2
                ;;

            -ck|--chunk[_-]size)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                chunk_size="${2}"
                shift 2
                ;;

            -csf|--csv[_-]scl[_-]fct)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                csv_scl_fct="${2}"
                shift 2
                ;;

            -cuf|--csv[_-]usr[_-]frg)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                csv_usr_frg="${2}"
                shift 2
                ;;

            -cdm|--csv[_-]dep[_-]min)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                csv_dep_min="${2}"
                shift 2
                ;;

            -cps|--csv[_-]pseudo)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                csv_pseudo="${2}"
                shift 2
                ;;

            -e|--eps)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                eps="${2}"
                shift 2
                ;;

            -s0|--skp_00|--skip[_-]00)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                skip_00="${2,,}"
                shift 2
                ;;

            --strict[_-]bins)
                strict_bins=true
                shift 1
                ;;

            -dn|--drp[_-]nan|--drop[_-]nan)
                drp_nan=true
                shift 1
                ;;

            -sp|--skp[_-]pfx)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                skp_pfx="${2}"
                shift 2
                ;;

            -dp|--dp)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                dp="${2}"
                shift 2
                ;;

            -deo|--dir[_-]eo)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                dir_eo="${2}"
                shift 2
                ;;

            -nj|--nam[_-]job)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
                }
                nam_job="${2}"
                shift 2
                ;;

            -mj|--max[_-]job)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_compute_signal >&2
                    return 1
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
                    return 1
                }
                time="${2}"
                shift 2
                ;;

            *)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                help_execute_compute_signal >&2
                return 1
                ;;
        esac
    done
}


#  Canonicalize mode and method aliases
function canonicalize_args() {
    case "${mode}" in
        s|sig|signal)
            mode="signal"

            if [[ -z "${method}" ]]; then method="norm"; fi

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
                        "'u', 'unadj', 'unadjusted', 's', 'smp', 'simple'," \
                        "'r', or 'raw' ('method=unadj'); 'f', 'frg', 'frag'," \
                        "'l', 'len', 'len_frg', or 'len_frag'" \
                        "('method=frag'); 'n', 'nrm', 'norm', or" \
                        "'normalized' ('method=norm')."
                    return 1
                    ;;
            esac

            validate_var "csv_fil_in" "${csv_fil_in}"
            validate_var_dir "csv_fil_in parent directory" \
                "$(dirname "${csv_fil_in%%[,;]*}")" 0 false
            check_str_delim "csv_fil_in" "${csv_fil_in}"
            ;;

        r|rat|ratio)
            mode="ratio"

            if [[ -z "${method}" ]]; then method="unadj"; fi

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
                        "'u', 'unadj', 'unadjusted', 's', 'smp', 'simple'," \
                        "'r', or 'raw' ('method=unadj'); '2', 'l2', 'lg2'," \
                        "or 'log2' ('method=log2'); 'ur', 'unadj_r'," \
                        "'unadjusted_r', 'sr', 'smp_r', 'simple_r', 'rr', or" \
                        "'raw_r' ('method=unadj_r'); or '2r', 'l2r', 'l2_r'," \
                        "'lg2_r', or 'log2_r' ('method=log2_r')."
                    return 1
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

            validate_var "csv_fil_in" "${csv_fil_in}"
            validate_var_dir "csv_fil_in parent directory" \
                "$(dirname "${csv_fil_in%%[,;]*}")" 0 false
            check_str_delim "csv_fil_in" "${csv_fil_in}"
            ;;

        *)
            echo_err \
                "invalid value for '--mode': '${mode}'. Expected 's', 'sig'," \
                "'signal', 'r', 'rat', 'ratio', 'c', 'coord', or" \
                "'coordinates'."
            return 1
            ;;
    esac
}


#  Validate scalar arguments and assign derived scalar defaults
function validate_args() {
    validate_var "env_nam" "${env_nam}"
    check_env_installed "${env_nam}"

    validate_var_dir  "dir_scr" "${dir_scr}" 0 false

    validate_var_file "scr_sub" "${scr_sub}"

    validate_var "threads" "${threads}"
    check_int_pos "${threads}" "threads"

    validate_var_dir "dir_out" "${dir_out}"

    if [[ "${mode}" =~ ^(signal|ratio)$ ]]; then
        validate_var "typ_out" "${typ_out}"
        case "${typ_out}" in
            bedGraph|bedGraph.gz|bedgraph|bedgraph.gz|bdg|bdg.gz|bg|bg.gz) : ;;
            *)
                echo_warn \
                    "unsupported value for '--typ_out' with" \
                    "'--mode ${mode}': '${typ_out}'. Coercing '--typ_out' to" \
                    "'bedGraph.gz'."
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
            #  If user didn’t supply '--siz_bin', apply a hardcoded default:
            #+ 10 bp
            if [[ -z "${siz_bin}" ]]; then siz_bin=10; fi

            check_int_pos "${siz_bin}" "siz_bin"

            case "${engine}" in
                chrom|window) : ;;
                *)
                    echo_err \
                        "'--engine' must be 'chrom' or 'window': '${engine}'."
                    return 1
                    ;;
            esac

            check_int_pos "${chunk_size}" "chunk_size"
        else
            #  For 'mode=ratio', ignore '--siz_bin' if the user supplied it
            if [[ -n "${siz_bin}" ]]; then
                echo_warn \
                    "argument '--siz_bin' is not applicable with '--mode" \
                    "${mode}'. Ignoring/unsetting '--siz_bin ${siz_bin}'."
                unset siz_bin
            fi
        fi
    else
        validate_var "typ_out" "${typ_out}"
        case "${typ_out}" in
            bedGraph|bedGraph.gz|bedgraph|bedgraph.gz|bdg|bdg.gz|bg|bg.gz)
                echo_warn \
                    "unsupported value for '--typ_out' with" \
                    "'--mode ${mode}': '${typ_out}'. Coercing '--typ_out' to" \
                    "'bed.gz'."
                typ_out="bed.gz"
                ;;

            bed|bed.gz) : ;;

            *)
                echo_err \
                    "invalid value for '--typ_out': '${typ_out}'. Expected" \
                    "'bed' or 'bed.gz'."
                return 1
                ;;
        esac

        #  '--siz_bin' is not applicable with 'mode=coord', so ignore if
        #+ supplied
        if [[ -n "${siz_bin}" ]]; then
            echo_warn \
                "argument '--siz_bin' is not applicable with" \
                "'--mode ${mode}'. Ignoring/unsetting '--siz_bin ${siz_bin}'."
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

    if [[ -n "${ref_fa}" ]]; then
        validate_var_file "ref_fa" "${ref_fa}"
    fi

    if [[ -n "${chr_sizes}" ]]; then
        validate_var_file "chr_sizes" "${chr_sizes}"
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
                    return 1
                    ;;
            esac
        fi

        if [[ -n "${csv_dep_min}" && -n "${csv_pseudo}" ]]; then
            echo_warn \
                "both '--csv_dep_min' and '--csv_pseudo' were supplied for" \
                "'--mode ratio'. This is allowed, but interpretability may" \
                "be reduced because both arguments stabilize low-depth ratio" \
                "behavior in different ways."
        fi
    fi

    validate_var "dp" "${dp}"
    check_int_pos "${dp}" "dp"

    if [[ -z "${dir_eo}" ]]; then dir_eo="${dir_out}/logs"; fi
    validate_var_dir "dir_eo" "${dir_eo}"

    if [[ -z "${nam_job}" ]]; then
        if [[ "${mode}" == "coord" ]]; then
            nam_job="compute_${mode}"
        else
            nam_job="compute_${mode}_${method}"
        fi
    fi
    validate_var "nam_job" "${nam_job}"
}


#  Parse input vectors, derive outputs, and validate per-element values
function prepare_vecs() {
    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        if [[ -n "${csv_fil_in}" ]]; then
            IFS=',' read -r -a arr_fil_in <<< "${csv_fil_in}"
            check_arr_nonempty "arr_fil_in" "csv_fil_in"

            for fil_in in "${arr_fil_in[@]}"; do
                validate_var_file "fil_in" "${fil_in}"
            done
            unset fil_in
        else
            echo_err "variable 'csv_fil_in' unexpectedly empty."
            return 1
        fi

        unset arr_fil_out && declare -ga arr_fil_out
        for i in "${arr_fil_in[@]}"; do
            base="$(basename "${i}")"
            base="${base%.bam}"
            base="${base%.cram}"

            if [[ -n "${prefix}" ]]; then
                if [[ "${base}" =~ ^IP_ ]]; then
                    echo_warn \
                        "filename '${base}' starts with 'IP_'. Stripping" \
                        "this prefix before applying custom '--prefix" \
                        "${prefix}'."
                    base="${base#IP_}"
                elif [[ "${base}" =~ ^in_ ]]; then
                    echo_warn \
                        "filename '${base}' starts with 'in_'. Stripping" \
                        "this prefix before applying custom '--prefix" \
                        "${prefix}'."
                    base="${base#in_}"
                fi
                arr_fil_out+=( "${dir_out}/${prefix}.${base}.${typ_out}" )
            else
                arr_fil_out+=( "${dir_out}/${base}.${typ_out}" )
            fi
        done

        check_arr_lengths "arr_fil_out" "arr_fil_in"

        for fil_in in "${arr_fil_in[@]}"; do
            if [[ "${fil_in,,}" == *.cram && -z "${ref_fa}" ]]; then
                echo_err \
                    "'--ref_fa' is required when '--csv_fil_in' contains" \
                    "CRAM input: '${fil_in}'."
                return 1
            fi
        done
        unset fil_in

        #  Scaling factors are only used for '--mode signal'
        if [[ "${mode}" == "signal" ]]; then
            if [[ -z "${csv_scl_fct}" ]]; then
                unset arr_scl_fct && declare -ga arr_scl_fct
                populate_array_empty arr_scl_fct "${#arr_fil_in[@]}"
            else
                IFS=',' read -r -a arr_scl_fct <<< "${csv_scl_fct}"
            fi

            for s in "${arr_scl_fct[@]}"; do
                if [[ "${s}" != "NA" ]]; then
                    check_flt_pos "${s}" "csv_scl_fct"
                fi
            done
            unset s

            check_arr_lengths "arr_scl_fct" "arr_fil_in"
        fi

        #  User-supplied fragment lengths are allowed for both 'signal' and
        #+ 'coord'
        if [[ -z "${csv_usr_frg}" ]]; then
            unset arr_usr_frg && declare -ga arr_usr_frg
            populate_array_empty arr_usr_frg "${#arr_fil_in[@]}"
        else
            IFS=',' read -r -a arr_usr_frg <<< "${csv_usr_frg}"
        fi

        for u in "${arr_usr_frg[@]}"; do
            if [[ "${u}" != "NA" ]]; then
                check_flt_pos "${u}" "csv_usr_frg"
            fi
        done
        unset u

        check_arr_lengths "arr_usr_frg" "arr_fil_in"
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
            return 1
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
            return 1
        fi

        pfx_lcl="${prefix}"
        exts=( bedGraph bedGraph.gz bedgraph bedgraph.gz bdg bdg.gz bg bg.gz )
        unset arr_fil_out && declare -ga arr_fil_out
        for i in "${arr_fil_A[@]}"; do
            base=$(basename "${i}")

            #  Strip 'IP_' prefix (if present)
            base="${base#IP_}"

            #  If not user-assigned, determine filename 'prefix' based on
            #+ provided arguments
            if [[ -z "${pfx_lcl}" ]]; then
                pfx_lcl="$(generate_pfx "${method}" "${csv_scl_fct}")"
            fi

            #  Remove file extensions
            for ext in "${exts[@]}"; do
                base="${base%."${ext}"}"
            done
            unset ext

            arr_fil_out+=( "${dir_out}/${pfx_lcl}_${base}.${typ_out}" )
        done
        unset i base pfx_lcl exts

        if [[ -z "${csv_scl_fct}" ]]; then
            #  Produce per-sample sentinel entries ("NA")
            unset arr_scl_fct && declare -ga arr_scl_fct
            populate_array_empty arr_scl_fct "${#arr_fil_A[@]}"
        else
            IFS=',' read -r -a arr_scl_fct <<< "${csv_scl_fct}"
        fi

        for s in "${arr_scl_fct[@]}"; do
            check_scl_fct_ratio "${s}"
        done
        unset s

        if [[ -z "${csv_dep_min}" ]]; then
            #  Produce per-sample sentinel entries ("NA")
            unset arr_dep_min && declare -ga arr_dep_min
            populate_array_empty arr_dep_min "${#arr_fil_A[@]}"
        else
            IFS=',' read -r -a arr_dep_min <<< "${csv_dep_min}"
        fi

        if [[ -z "${csv_pseudo}" ]]; then
            #  Produce per-sample sentinel entries ("NA")
            unset arr_pseudo && declare -ga arr_pseudo
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
                return 1
            fi

            IFS=':' read -r pseudo_A pseudo_B <<< "${p}"

            if [[ -z "${pseudo_A}" || -z "${pseudo_B}" ]]; then
                echo_err \
                    "invalid pseudocount spec in '--csv_pseudo': '${p}'." \
                    "Expected 'A' or 'A:B'."
                return 1
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
        check_arr_lengths "arr_fil_out" "arr_fil_A"
        check_arr_lengths "arr_scl_fct" "arr_fil_A"
        check_arr_lengths "arr_dep_min" "arr_fil_A"
        check_arr_lengths "arr_pseudo"  "arr_fil_A"
    fi
}


#  Validate final vector shapes after preparation and auto-population
function validate_vecs() {
    if [[ "${mode}" == "signal" ]]; then
        check_arr_lengths "arr_fil_in" "arr_fil_out" || return 1
        check_arr_lengths "arr_fil_in" "arr_scl_fct" || return 1
        check_arr_lengths "arr_fil_in" "arr_usr_frg" || return 1
    elif [[ "${mode}" == "coord" ]]; then
        check_arr_lengths "arr_fil_in" "arr_fil_out" || return 1
        check_arr_lengths "arr_fil_in" "arr_usr_frg" || return 1
    else
        check_arr_lengths "arr_fil_A" "arr_fil_B"    || return 1
        check_arr_lengths "arr_fil_A" "arr_fil_out"  || return 1
        check_arr_lengths "arr_fil_A" "arr_scl_fct"  || return 1
        check_arr_lengths "arr_fil_A" "arr_dep_min"  || return 1
        check_arr_lengths "arr_fil_A" "arr_pseudo"   || return 1
    fi
}


#  Configure Slurm, GNU Parallel, or serial execution parameters
function config_exec() {
    validate_var "max_job" "${max_job}"
    check_int_pos "${max_job}" "max_job"

    if [[ "${slurm}" == "true" ]]; then
        if [[ "${mode}" == "ratio" ]]; then
            max_job=$(reset_max_job "${max_job}" "${#arr_fil_out[@]}")
        else
            max_job=$(reset_max_job "${max_job}" "${#arr_fil_in[@]}")
        fi

        validate_var "time" "${time}"
        check_format_time "${time}"
    elif [[ "${max_job}" -le 1 ]]; then
        #  Serial local execution does not require parallel job detection
        par_job=1
        unset time

        validate_var "par_job" "${par_job}"
        check_int_pos "${par_job}" "par_job"
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
        "arr_fil_in" "arr_fil_A" "arr_fil_B" "arr_fil_out" \
        "arr_scl_fct" "arr_usr_frg" "arr_dep_min" "arr_pseudo"

    if [[ "${mode}" == "signal" ]]; then
        summarize_sig_norm "${method}" "${csv_scl_fct}"
    fi
}


#  Activate environment
function setup_env() {
    local out
    local -a env_msg

    #TODO: this environment activation block is repeated verbatim many across
    #+     the 'execute_*.sh' scripts: modularize
    env_msg=(
        "'handle_env' failed for 'env_nam=${env_nam}'. Check that Conda/Mamba"
        "are available and that the environment exists."
    )

    if [[ "${verbose}" == "true" ]]; then
        if \
            out="$(handle_env "${env_nam}")"
        then
            print_banner_pretty \
                -tx "${out:-"${env_nam} already active."}" \
                -w "%"
            echo
            echo
        else
            echo_err "${env_msg[*]}"
            return 1
        fi
    else
        if ! \
            handle_env "${env_nam}" > /dev/null 2>&1
        then
            echo_err "${env_msg[*]}"
            return 1
        fi
    fi
}


#  Check tools needed by the selected dispatch mode
function check_tools() {
    check_pgrm_path compute_signal || return 1
    check_pgrm_path compute_signal_ratio || return 1

    if [[ "${slurm}" == "true" ]]; then
        check_pgrm_path sbatch || return 1
    elif [[ "${par_job}" -gt 1 ]]; then
        check_pgrm_path parallel || return 1
    fi
}


#  Print hardcoded defaults, argument values, and prepared arrays
function print_state_debug() {
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
        echo "csv_fil_in=${csv_fil_in:-UNSET}"
        echo "ref_fa=${ref_fa:-UNSET}"
        echo "chr_sizes=${chr_sizes:-UNSET}"
        echo "csv_fil_A=${csv_fil_A:-UNSET}"
        echo "csv_fil_B=${csv_fil_B:-UNSET}"
        echo "dir_out=${dir_out}"
        echo "typ_out=${typ_out}"
        echo "prefix=${prefix:-UNSET}"
        echo "track=${track}"
        echo "siz_bin=${siz_bin:-UNSET}"
        echo "chunk_size=${chunk_size:-UNSET}"
        echo "engine=${engine:-UNSET}"
        echo "csv_scl_fct=${csv_scl_fct:-UNSET}"
        echo "csv_usr_frg=${csv_usr_frg:-UNSET}"
        echo "csv_dep_min=${csv_dep_min:-UNSET}"
        echo "csv_pseudo=${csv_pseudo:-UNSET}"
        echo "eps=${eps:-UNSET}"
        echo "skip_00=${skip_00:-UNSET}"
        echo "strict_bins=${strict_bins}"
        echo "drp_nan=${drp_nan}"
        echo "skp_pfx=${skp_pfx:-UNSET}"
        echo "dp=${dp}"
        echo "dir_eo=${dir_eo}"
        echo "nam_job=${nam_job}"
        echo "max_job=${max_job:-UNSET}"
        echo "slurm=${slurm}"
        echo "time=${time:-UNSET}"
        echo
        echo
        print_banner_pretty "Arrays derived from variables"
        echo

        if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
            echo "arr_fil_in=( ${arr_fil_in[*]} )"
            echo
        elif [[ "${mode}" == "ratio" ]]; then
            echo "arr_fil_A=( ${arr_fil_A[*]} )"
            echo
            echo "arr_fil_B=( ${arr_fil_B[*]} )"
            echo
        fi

        echo "arr_fil_out=( ${arr_fil_out[*]} )"
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
}


#  Serialize prepared arrays into comma-delimited submit-wrapper arguments
function serialize_vecs() {
    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        csv_fil_in=$(echo "${arr_fil_in[*]}"  | tr ' ' ',')
        csv_fil_out=$(echo "${arr_fil_out[*]}" | tr ' ' ',')
        csv_usr_frg=$(echo "${arr_usr_frg[*]}" | tr ' ' ',')

        if [[ "${mode}" == "signal" ]]; then
            csv_scl_fct=$(echo "${arr_scl_fct[*]}" | tr ' ' ',')
        fi
    elif [[ "${mode}" == "ratio" ]]; then
        csv_fil_A=$(echo "${arr_fil_A[*]}"  | tr ' ' ',')
        csv_fil_B=$(echo "${arr_fil_B[*]}"  | tr ' ' ',')
        csv_fil_out=$(echo "${arr_fil_out[*]}" | tr ' ' ',')
        csv_scl_fct=$(echo "${arr_scl_fct[*]}" | tr ' ' ',')
        csv_dep_min=$(echo "${arr_dep_min[*]}" | tr ' ' ',')
        csv_pseudo=$(echo "${arr_pseudo[*]}" | tr ' ' ',')
    fi
}


#  Print serialized variables used for submit-wrapper command construction
function print_vecs_serialized() {
    if [[ "${verbose}" == "true" ]]; then
        print_banner_pretty "Variable assignments constructed from arrays"
        echo

        if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
            echo "csv_fil_in=\"${csv_fil_in}\""
            echo
            echo "ref_fa=\"${ref_fa:-UNSET}\""
            echo
            echo "chr_sizes=\"${chr_sizes:-UNSET}\""
            echo
        elif [[ "${mode}" == "ratio" ]]; then
            echo "csv_fil_A=\"${csv_fil_A}\""
            echo
            echo "csv_fil_B=\"${csv_fil_B}\""
            echo
        fi

        echo "csv_fil_out=\"${csv_fil_out}\""
        echo
        echo "chr_sizes=\"${chr_sizes:-UNSET}\""
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
}


#  Dispatch Slurm, GNU Parallel, or serial work
function run_jobs() {
    local config idx log_out log_err
    local -a cmd_slurm

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
                --output="${dir_eo}/${nam_job}.%A-%a.stdout.txt"
                --error="${dir_eo}/${nam_job}.%A-%a.stderr.txt"
                --array="1-${#arr_fil_out[@]}%${max_job}"
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
            config="${dir_eo}/${nam_job}.config_parallel.txt"

            if [[ -f "${config}" ]]; then rm "${config}"; fi

            for idx in "${!arr_fil_out[@]}"; do
                build_cmd "${idx}"
                cmd_bld=( "${BASH}" "${cmd_bld[@]}" )

                IFS=';' read -r log_out log_err < <(
                    get_submit_logs "${arr_fil_out[idx]}"
                )

                {
                    print_built_cmd "${log_out}" "${log_err}"
                } >> "${config}" || {
                    echo_err "failed to write command, index no. '${idx}'."
                    return 1
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

                for idx in "${!arr_fil_out[@]}"; do
                    build_cmd "${idx}"
                    cmd_bld=( "${BASH}" "${cmd_bld[@]}" )

                    IFS=';' read -r log_out log_err < <(
                        get_submit_logs "${arr_fil_out[idx]}"
                    )

                    print_built_cmd "${log_out}" "${log_err}"
                    echo
                done
                echo
            fi

            if [[ "${dry_run}" == "false" ]]; then
                for idx in "${!arr_fil_out[@]}"; do
                    build_cmd "${idx}"
                    cmd_bld=( "${BASH}" "${cmd_bld[@]}" )

                    IFS=';' read -r log_out log_err < <(
                        get_submit_logs "${arr_fil_out[idx]}"
                    )

                    "${cmd_bld[@]}" >> "${log_out}" 2>> "${log_err}"
                done
            fi
        fi
    fi
}


#  Main script execution
function main() {
    init_defs
    source_helpers_execute || return 1

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        help_execute_compute_signal >&2
        return 0
    elif [[ "${1}" =~ ^(-d|--details)$ ]]; then
        detail_execute_compute_signal >&2
        return 0
    elif [[ "${1}" =~ ^(-ah|--all[_-]h[e]?lp)$ ]]; then
        help_execute_compute_signal >&2
        echo >&2
        detail_execute_compute_signal "--no-usage" >&2
        return 0
    fi

    parse_args "$@"       || return 1
    canonicalize_args     || return 1
    validate_args         || return 1
    prepare_vecs          || return 1
    validate_vecs         || return 1
    config_exec           || return 1
    setup_env             || return 1
    check_tools           || return 1
    print_state_debug     || return 1
    serialize_vecs        || return 1
    print_vecs_serialized || return 1
    run_jobs              || return 1
}


main "$@"
