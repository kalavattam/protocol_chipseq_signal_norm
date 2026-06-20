#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: execute_download_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
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

#  Set the path to the 'scripts' directory
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"


#  Source shared helpers
function source_helpers_execute() {
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
        check_env \
        check_inputs \
        check_numbers \
        format_outputs \
        handle_env \
        help/help_execute_download_fastqs \
        wrap_cmd \
        || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper scripts." >&2
            return 1
        }
}


function build_cmd() {
    local idx="${1:-}"
    local show_help

    unset cmd_bld && declare -ga cmd_bld

    show_help=$(cat << EOM
Usage:
  build_cmd [--help] idx

Description:
  Construct the command array 'cmd_bld' for one call to 'submit_download_fastqs.sh'.

Positional arguments:
  1  idx  <int>
    Zero-based index into the unique parsed download arrays.

Expected globals:
  scr_sub list_job_acc list_job_url_1 list_job_url_2 dir_out dir_sym list_job_cus dir_eo nam_job

Generated globals:
  cmd_bld
    Indexed array containing one positional call to 'submit_download_fastqs.sh'.

Returns:
  0 if 'cmd_bld' is constructed successfully; otherwise 1.

Notes:
  This preserves the positional interface expected by 'submit_download_fastqs.sh'.
EOM
    )

    if [[ "${idx}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${idx}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'idx', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    check_int_nonneg "${idx}" "idx" || return 1

    # shellcheck disable=SC2034
    cmd_bld=(
        "${scr_sub}"
            "${list_job_acc[idx]}"
            "${list_job_url_1[idx]}"
            "${list_job_url_2[idx]}"
            "${dir_out}"
            "${dir_sym}"
            "${list_job_cus[idx]}"
            "${dir_eo}"
            "${nam_job}"
    )
}


function link_name_custom() {
    local idx="${1:-}"
    local acc cus url_2
    local src dst

    check_int_nonneg "${idx}" "idx" || return 1

    acc="${list_acc[idx]}"
    cus="${list_cus[idx]}"
    url_2="${list_url_2[idx]}"

    if [[ "${url_2}" != "NA" ]]; then
        src="${dir_out}/${acc}_R1.fastq.gz"
        dst="${dir_sym}/${cus}_R1.fastq.gz"
        if [[ ! -s "${src}" ]]; then
            echo_err "downloaded FASTQ is missing or empty: '${src}'."
            return 1
        fi
        ln -sf "${src}" "${dst}" || {
            echo_err "failed to create symlink '${dst}' -> '${src}'."
            return 1
        }

        src="${dir_out}/${acc}_R2.fastq.gz"
        dst="${dir_sym}/${cus}_R2.fastq.gz"
        if [[ ! -s "${src}" ]]; then
            echo_err "downloaded FASTQ is missing or empty: '${src}'."
            return 1
        fi
        ln -sf "${src}" "${dst}" || {
            echo_err "failed to create symlink '${dst}' -> '${src}'."
            return 1
        }
    else
        src="${dir_out}/${acc}.fastq.gz"
        dst="${dir_sym}/${cus}.fastq.gz"
        if [[ ! -s "${src}" ]]; then
            echo_err "downloaded FASTQ is missing or empty: '${src}'."
            return 1
        fi
        ln -sf "${src}" "${dst}" || {
            echo_err "failed to create symlink '${dst}' -> '${src}'."
            return 1
        }
    fi
}


#  Initialize hardcoded argument variables
function init_args_hardcoded() {
    env_nam="env_protocol"
    scr_sub="${dir_scr}/submit_download_fastqs.sh"
}


#  Initialize argument variables, assigning default values where applicable
function init_arg_defs() {
    verbose=false
    dry_run=false
    threads=1
    infile=""
    dir_out=""
    dir_sym=""
    nam_job="download_fastqs"
    dir_eo=""
    slurm=false
    time="3:00:00"
    config=""
    has_dup_dl=false

    unset list_acc list_url_1 list_url_2 list_cus
    unset list_job_acc list_job_url_1 list_job_url_2 list_job_cus
    unset cmd_bld

    declare -ga list_acc list_url_1 list_url_2 list_cus
    declare -ga list_job_acc list_job_url_1 list_job_url_2 list_job_cus
    declare -ga cmd_bld
}


#  Initialize all script defaults
function init_defs() {
    init_args_hardcoded
    init_arg_defs
}


#  Parse command-line arguments
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

            -t|--thr|--threads)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_download_fastqs >&2
                    return 1
                }
                threads="${2}"
                shift 2
                ;;

            -i|-fi|--infile|--fil[_-]in)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_download_fastqs >&2
                    return 1
                }
                infile="${2}"
                shift 2
                ;;

            -do|--dir[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_download_fastqs >&2
                    return 1
                }
                dir_out="${2}"
                shift 2
                ;;

            -dy|--dir[_-]sym|--dir[_-]symlink)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_download_fastqs >&2
                    return 1
                }
                dir_sym="${2}"
                shift 2
                ;;

            -nj|--nam[_-]job)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_download_fastqs >&2
                    return 1
                }
                nam_job="${2}"
                shift 2
                ;;

            -eo|--dir[_-]eo)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_download_fastqs >&2
                    return 1
                }
                dir_eo="${2}"
                shift 2
                ;;

            -sl|--slurm)
                slurm=true
                shift 1
                ;;

            -tm|--time)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_download_fastqs >&2
                    return 1
                }
                time="${2}"
                shift 2
                ;;

            *)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                help_execute_download_fastqs >&2
                return 1
                ;;
        esac
    done
}


#  Check parsed argument values
function validate_args() {
    validate_var "env_nam" "${env_nam}" || return 1
    check_env_installed "${env_nam}" || return 1

    validate_var_dir "dir_scr" "${dir_scr}" 0 false || return 1

    validate_var_file "scr_sub" "${scr_sub}" || return 1

    validate_var "threads" "${threads}" || return 1
    check_int_pos "${threads}" "threads" || return 1

    validate_var_file "infile" "${infile}" || return 1

    validate_var_dir "dir_out" "${dir_out}" || return 1

    validate_var_dir "dir_sym" "${dir_sym}" || return 1

    validate_var "nam_job" "${nam_job}" || return 1

    if [[ -z "${dir_eo}" ]]; then dir_eo="${dir_out}/logs"; fi
    validate_var_dir "dir_eo" "${dir_eo}" || return 1

    if [[ "${slurm}" == "true" ]]; then
        validate_var "time" "${time}" || return 1
        check_format_time "${time}" || return 1
    fi
}


#  Activate environment
function setup_env() {
    local out
    local -a env_msg

    #TODO: this environment activation block is repeated verbatim many across
    #+     the 'execute_*.sh' scripts: modularize
    env_msg=(
        "'handle_env' failed for 'env_nam=${env_nam}'. Check that Conda/Mamba are"
        "available and that the environment exists."
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
    check_pgrm_path cut  || return 1
    check_pgrm_path wget || return 1

    if [[ "${slurm}" == "true" ]]; then
        check_pgrm_path sbatch || return 1
    fi

    if [[ "${threads}" -gt 1 ]]; then
        check_pgrm_path parallel || return 1
    fi
}


#  Parse the download metadata TSV into raw vectors
function parse_tbl_dl() {
    local run_acc_idx=""
    local custom_name_idx=""
    local url_col_idx=""
    local iter=0
    local line i
    local run_acc custom_name urls
    local -a arr_hdr arr_url_fq

    unset list_acc list_url_1 list_url_2 list_cus
    declare -ga list_acc list_url_1 list_url_2 list_cus

    while IFS=$'\t' read -r line || [[ -n "${line}" ]]; do
        (( iter++ )) || true

        if [[ "${verbose}" == "true" ]]; then
            echo "Processing line #${iter}: ${line}"
        fi

        if [[ "${iter}" -eq 1 ]]; then
            IFS=$'\t' read -r -a arr_hdr <<< "${line}"

            for i in "${!arr_hdr[@]}"; do
                case "${arr_hdr[i]}" in
                    "run_accession") run_acc_idx=${i}     ;;
                    "custom_name")   custom_name_idx=${i} ;;
                    "fastq_ftp")     url_col_idx=${i}     ;;
                    "fastq_https")   url_col_idx=${i}     ;;
                esac
            done

            if [[ -z "${run_acc_idx}" ]]; then
                echo_err "required column 'run_accession' was not found in header."
                return 1
            elif [[ -z "${custom_name_idx}" ]]; then
                echo_err "required column 'custom_name' was not found in header."
                return 1
            elif [[ -z "${url_col_idx}" ]]; then
                echo_err \
                    "no valid URL column found in header. Expected" \
                    "'fastq_ftp' or 'fastq_https'."
                return 1
            fi

            continue
        fi

        run_acc=$(echo "${line}" | cut -f $(( run_acc_idx + 1 )))
        custom_name=$(echo "${line}" | cut -f $(( custom_name_idx + 1 )))
        urls=$(echo "${line}" | cut -f $(( url_col_idx + 1 )))

        if [[ -z "${run_acc}" || "${run_acc}" == "NA" ]]; then
            run_acc="SRR_undefined_${iter}"
        fi

        if [[ -z "${custom_name}" || "${custom_name}" == "NA" ]]; then
            echo_err \
                "missing required custom_name for TSV line '${iter}'" \
                "('${run_acc}')."
            return 1
        fi

        if [[ -z "${urls}" || "${urls}" == "NA" ]]; then
            echo_err \
                "missing required FASTQ URL field for TSV line '${iter}'" \
                "('${run_acc}')."
            return 1
        fi

        IFS=';' read -r -a arr_url_fq <<< "${urls}"

        list_acc+=( "${run_acc}" )
        list_cus+=( "${custom_name}" )

        if [[ ${#arr_url_fq[@]} -eq 1 ]]; then
            if [[ -z "${arr_url_fq[0]}" || "${arr_url_fq[0]}" == "NA" ]]; then
                echo_err \
                    "missing FASTQ URL for TSV line '${iter}'" \
                    "('${run_acc}')."
                return 1
            fi

            list_url_1+=( "${arr_url_fq[0]}" )
            list_url_2+=( "NA" )
        elif [[ ${#arr_url_fq[@]} -eq 2 ]]; then
            if [[ -z "${arr_url_fq[0]}" || -z "${arr_url_fq[1]}" ]]; then
                echo_err \
                    "missing one or more paired-end FASTQ URLs for TSV line" \
                    "'${iter}' ('${run_acc}')."
                return 1
            fi

            list_url_1+=( "${arr_url_fq[0]}" )
            list_url_2+=( "${arr_url_fq[1]}" )
        else
            echo_err "unexpected number of FASTQ URLs for '${run_acc}'."
            return 1
        fi
    done < "${infile}"
}


#  Validate raw download vectors
function validate_vecs_dl() {
    check_arr_nonempty "list_acc"   "download entries"       || return 1
    check_arr_nonempty "list_url_1" "FASTQ URL entries"      || return 1
    check_arr_nonempty "list_url_2" "second FASTQ URL entries" || return 1
    check_arr_nonempty "list_cus"   "custom-name entries"    || return 1

    check_arr_lengths "list_acc" "list_url_1" || return 1
    check_arr_lengths "list_acc" "list_url_2" || return 1
    check_arr_lengths "list_acc" "list_cus"   || return 1
}


#  Build the unique raw-download job vectors
function parse_jobs_uniq() {
    local i key
    local -A seen_key_acc seen_cus

    unset list_job_acc list_job_url_1 list_job_url_2 list_job_cus

    declare -ga list_job_acc list_job_url_1 list_job_url_2 list_job_cus

    has_dup_dl=false

    for i in "${!list_acc[@]}"; do
        key="${list_url_1[i]}"$'\t'"${list_url_2[i]}"

        if [[ -n "${seen_cus[${list_cus[i]}]:-}" ]]; then
            echo_err "duplicate custom_name in input TSV: '${list_cus[i]}'."
            return 1
        fi
        seen_cus["${list_cus[i]}"]=1

        if [[ -n "${seen_key_acc[${list_acc[i]}]:-}" ]]; then
            if [[ "${seen_key_acc[${list_acc[i]}]}" != "${key}" ]]; then
                echo_err \
                    "run_accession '${list_acc[i]}' appears with conflicting" \
                    "FASTQ URL values."
                return 1
            fi
            has_dup_dl=true
            continue
        fi

        seen_key_acc["${list_acc[i]}"]="${key}"
        list_job_acc+=( "${list_acc[i]}" )
        list_job_url_1+=( "${list_url_1[i]}" )
        list_job_url_2+=( "${list_url_2[i]}" )
        list_job_cus+=( "${list_cus[i]}" )
    done

    check_arr_nonempty "list_job_acc" "unique download entries" || return 1
    check_arr_lengths  "list_job_acc" "list_job_url_1"          || return 1
    check_arr_lengths  "list_job_acc" "list_job_url_2"          || return 1
    check_arr_lengths  "list_job_acc" "list_job_cus"            || return 1

    if [[ "${slurm}" == "true" && "${has_dup_dl}" == "true" ]]; then
        echo_err \
            "duplicate run_accession rows require local execution so symlinks for" \
            "all custom_name values can be created after the unique download" \
            "finishes."
        return 1
    fi
}


#  Report current scalar and vector state in verbose mode
function print_state_debug() {
    local el

    if [[ "${verbose}" != "true" ]]; then
        return 0
    fi

    print_banner_pretty "Hardcoded variable assignments"
    echo
    echo "env_nam=${env_nam}"
    echo "dir_scr=${dir_scr}"
    echo "scr_sub=${scr_sub}"
    echo
    echo
    print_banner_pretty "Argument variable assignments"
    echo
    echo "verbose=${verbose}"
    echo "dry_run=${dry_run}"
    echo "threads=${threads}"
    echo "infile=${infile}"
    echo "dir_out=${dir_out}"
    echo "dir_sym=${dir_sym}"
    echo "nam_job=${nam_job}"
    echo "dir_eo=${dir_eo}"
    echo "slurm=${slurm}"
    echo "time=${time}"
    echo

    echo
    print_banner_pretty "list_acc elements"
    echo
    for el in "${list_acc[@]}"; do echo "${el}"; done
    echo
    echo

    print_banner_pretty "list_cus elements"
    echo
    for el in "${list_cus[@]}"; do echo "${el}"; done
    echo
    echo

    print_banner_pretty "list_url_1 elements"
    echo
    for el in "${list_url_1[@]}"; do echo "${el}"; done
    echo
    echo

    print_banner_pretty "list_url_2 elements"
    echo
    for el in "${list_url_2[@]}"; do echo "${el}"; done
    echo
    echo
}


#  Create a command file for GNU Parallel-backed execution if needed
function prepare_config_parallel() {
    local i

    if [[ "${threads}" -le 1 ]]; then
        return 0
    fi

    config="${dir_eo}/${nam_job}.config_parallel.txt"

    if [[ -f "${config}" ]]; then rm "${config}"; fi

    for i in "${!list_job_acc[@]}"; do
        build_cmd "${i}" || return 1
        cmd_bld=( "${BASH}" "${cmd_bld[@]}" )

        {
            print_built_cmd
        } >> "${config}" || {
            echo_err "failed to write command, index no. '${i}'."
            return 1
        }
    done
}


#  Create all requested custom-name symlinks after local unique downloads
function link_names_custom() {
    local i

    for i in "${!list_acc[@]}"; do
        link_name_custom "${i}" || return 1
    done
}


#  Dispatch Slurm, GNU Parallel, or serial download work
function run_jobs() {
    local i cmd_wrap
    local -a cmd_para cmd_slurm

    if [[ "${slurm}" == "true" ]]; then
        if [[ "${threads}" -gt 1 ]]; then
            cmd_para=(
                parallel
                    --jobs "${threads}"
            )

            cmd_wrap="$(
                printf '%q ' "${cmd_para[@]}"
                printf '< %q' "${config}"
            )"

            cmd_slurm=(
                sbatch
                    --job-name="${nam_job}"
                    --nodes=1
                    --cpus-per-task="${threads}"
                    --time="${time}"
                    --error="${config%.config_parallel.txt}.%A.stderr.txt"
                    --output="${config%.config_parallel.txt}.%A.stdout.txt"
                    --wrap="${cmd_wrap}"
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
            echo_err \
                "Slurm submissions require 'threads > 1'; current value:" \
                "threads=${threads}."
            return 1
        fi
    else
        if [[ "${threads}" -gt 1 ]]; then
            if [[ "${dry_run}" == "true" || "${verbose}" == "true" ]]; then
                print_banner_pretty "GNU Parallel execution"
                echo
                parallel --jobs "${threads}" --dryrun < "${config}"
                echo
                echo
            fi

            if [[ "${dry_run}" == "false" ]]; then
                parallel --jobs "${threads}" < "${config}"
                link_names_custom || return 1
            fi
        else
            for i in "${!list_job_acc[@]}"; do
                build_cmd "${i}" || return 1
                cmd_bld=( "${BASH}" "${cmd_bld[@]}" )

                if [[ "${dry_run}" == "true" || "${verbose}" == "true" ]]; then
                    print_banner_pretty "Serial execution"
                    echo
                    print_built_cmd
                    echo
                fi

                if [[ "${dry_run}" == "false" ]]; then
                    "${cmd_bld[@]}"
                fi
            done

            if [[ "${dry_run}" == "false" ]]; then
                link_names_custom || return 1
            fi
        fi
    fi
}


#  Main script execution
function main() {
    init_defs
    source_helpers_execute || return 1

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        help_execute_download_fastqs >&2
        return 0
    fi

    parse_args "$@"         || return 1
    validate_args           || return 1
    setup_env               || return 1
    check_tools             || return 1
    parse_tbl_dl            || return 1
    validate_vecs_dl        || return 1
    parse_jobs_uniq         || return 1
    print_state_debug       || return 1
    prepare_config_parallel || return 1
    run_jobs                || return 1
}


main "$@"
