#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: execute_download_fastqs.sh
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
    help/help_execute_download_fastqs \
    wrap_cmd \
    || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source required helper scripts." >&2
        exit 1
    }

unset fnc_src


function build_cmd() {
    local idx="${1:-}"
    local show_help

    unset cmd_bld && declare -ga cmd_bld

    show_help=$(cat << EOM
Usage:
  build_cmd [-h|--hlp|--help] idx

Description:
  Construct the command array 'cmd_bld' for one call to 'submit_download_fastqs.sh'.

Positional argument:
  1  idx  <int>  Zero-based index into the parsed download arrays.

Expected globals:
  scr_sub list_acc list_url_1 list_url_2 dir_out dir_sym list_cus err_out nam_job

Returns:
  0 if 'cmd_bld' is constructed successfully; otherwise 1.

Notes:
  - 'cmd_bld' is written as a global indexed array.
  - This preserves the positional interface expected by 'submit_download_fastqs.sh'.
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
            "${list_acc[idx]}"
            "${list_url_1[idx]}"
            "${list_url_2[idx]}"
            "${dir_out}"
            "${dir_sym}"
            "${list_cus[idx]}"
            "${err_out}"
            "${nam_job}"
    )
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_protocol"
scr_sub="${dir_scr}/submit_download_fastqs.sh"

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=1
infile=""
dir_out=""
dir_sym=""
nam_job="download_fastqs"
err_out=""
slurm=false
time="3:00:00"

#  Define help message
show_help=$(:)

#  Parse arguments
if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
    help_execute_download_fastqs >&2
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
                help_execute_download_fastqs >&2
                exit 1
            }
            threads="${2}"
            shift 2
            ;;

        -i|-fi|--infile|--fil[_-]in)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_download_fastqs >&2
                exit 1
            }
            infile="${2}"
            shift 2
            ;;

        -do|--dir[_-]out)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_download_fastqs >&2
                exit 1
            }
            dir_out="${2}"
            shift 2
            ;;

        -dy|--dir[_-]sym|--dir[_-]symlink)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_download_fastqs >&2
                exit 1
            }
            dir_sym="${2}"
            shift 2
            ;;

        -nj|--nam[_-]job)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_download_fastqs >&2
                exit 1
            }
            nam_job="${2}"
            shift 2
            ;;

        -eo|--err[_-]out)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_download_fastqs >&2
                exit 1
            }
            err_out="${2}"
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
                exit 1
            }
            time="${2}"
            shift 2
            ;;

        *)
            echo_err "unknown option/parameter passed: '${1}'."
            echo >&2
            help_execute_download_fastqs >&2
            exit 1
            ;;
    esac
done

#  Check arguments
validate_var "env_nam" "${env_nam}"
check_env_installed "${env_nam}"

validate_var_dir  "dir_scr" "${dir_scr}" 0 false

validate_var_file "scr_sub" "${scr_sub}"

validate_var "threads" "${threads}"
check_int_pos "${threads}" "threads"

validate_var_file "infile" "${infile}"

validate_var_dir "dir_out" "${dir_out}"

validate_var_dir "dir_sym" "${dir_sym}"

validate_var "nam_job" "${nam_job}"

if [[ -z "${err_out}" ]]; then err_out="${dir_out}/logs"; fi
validate_var_dir "err_out" "${err_out}"

if [[ "${slurm}" == "true" ]]; then
    validate_var "time" "${time}"
    check_format_time "${time}"
fi

#  Activate environment and check that dependencies are in PATH
env_msg=(
    "'handle_env' failed for 'env_nam=${env_nam}'. Check that Conda/Mamba are"
    "available and that the environment exists."
)

if [[ "${verbose}" == "true" ]]; then
    if ! handle_env "${env_nam}"; then
        echo_err "${env_msg[*]}"
        exit 1
    fi

    echo
else
    if ! handle_env "${env_nam}" > /dev/null 2>&1; then
        echo_err "${env_msg[*]}"
        exit 1
    fi
fi

check_pgrm_path cut
check_pgrm_path wget

if [[ "${slurm}" == "true" ]]; then
    check_pgrm_path sbatch
fi

if [[ "${threads}" -gt 1 ]]; then
    check_pgrm_path parallel
fi


#  Do the main work ===========================================================
#  Report argument variable assignments if in "verbose mode"
if [[ "${verbose}" == "true" ]]; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo
    echo "env_nam=${env_nam}"
    echo "dir_scr=${dir_scr}"
    echo "scr_sub=${scr_sub}"
    echo
    echo
    echo "###################################"
    echo "## Argument variable assignments ##"
    echo "###################################"
    echo
    echo "verbose=${verbose}"
    echo "dry_run=${dry_run}"
    echo "threads=${threads}"
    echo "infile=${infile}"
    echo "dir_out=${dir_out}"
    echo "dir_sym=${dir_sym}"
    echo "nam_job=${nam_job}"
    echo "err_out=${err_out}"
    echo "slurm=${slurm}"
    echo "time=${time}"
    echo
fi

#  Initialize TSV-parsing variables used under 'set -u'
run_acc_idx=""
custom_name_idx=""
url_col_idx=""
# url_col=""  #NOTE: no longer used

#  Create arrays to store SRR entries, URLs, and custom names
unset      list_acc list_url_1 list_url_2 list_cus
typeset -a list_acc list_url_1 list_url_2 list_cus

#  Read the TSV file, processing each line to extract SRR accessions (if
#+ available), URLs, and custom names
iter=0
while IFS=$'\t' read -r line; do
    (( iter++ )) || true  # (Prevent script exit if `set -e`)

    if [[ "${verbose}" == "true" ]]; then echo "Processing line #${iter}: ${line}"; fi

    #  Parse the header and detect available columns
    if [[ "${iter}" -eq 1 ]]; then
        IFS=$'\t' read -r -a headers <<< "${line}"

        #  Determine the index of the required columns dynamically
        for i in "${!headers[@]}"; do
            case "${headers[i]}" in
                "run_accession") run_acc_idx=${i}     ;;
                "custom_name")   custom_name_idx=${i} ;;
                "fastq_ftp")     url_col_idx=${i}     ;;
                "fastq_https")   url_col_idx=${i}     ;;

                #NOTE: 'url_col' is no longer used
                # "fastq_ftp")   url_col_idx=${i}; url_col='fastq_ftp'   ;;
                # "fastq_https") url_col_idx=${i}; url_col='fastq_https' ;;
            esac
        done

        #  Ensure required columns were found
        if [[ -z "${run_acc_idx}" ]]; then
            echo_err "required column 'run_accession' was not found in header."
            exit 1
        elif [[ -z "${custom_name_idx}" ]]; then
            echo_err "required column 'custom_name' was not found in header."
            exit 1
        elif [[ -z "${url_col_idx}" ]]; then
            echo_err \
                "no valid URL column found in header. Expected" \
                "'fastq_ftp' or 'fastq_https'."
            exit 1
        fi

        continue
    fi

    #  Read each column based on detected header indices
    run_acc=$(echo "${line}" | cut -f $(( run_acc_idx + 1 )))
    custom_name=$(echo "${line}" | cut -f $(( custom_name_idx + 1 )))
    urls=$(echo "${line}" | cut -f $(( url_col_idx + 1 )))

    #  Handle missing run accession entries
    if [[ -z "${run_acc}" || "${run_acc}" == "NA" ]]; then
        run_acc="SRR_undefined_${iter}"
    fi

    if [[ -z "${custom_name}" || "${custom_name}" == "NA" ]]; then
        echo_err \
            "missing required custom_name for TSV line '${iter}'" \
            "('${run_acc}')."
        exit 1
    fi

    if [[ -z "${urls}" || "${urls}" == "NA" ]]; then
        echo_err \
            "missing required FASTQ URL field for TSV line '${iter}'" \
            "('${run_acc}')."
        exit 1
    fi

    #  Parse the FASTQ URLs (paired-end or single-end)
    IFS=';' read -r -a fastq_urls <<< "${urls}"

    #  Add to arrays
    list_acc+=( "${run_acc}" )
    list_cus+=( "${custom_name}" )

    #  For single-end data
    if [[ ${#fastq_urls[@]} -eq 1 ]]; then
        if [[ -z "${fastq_urls[0]}" || "${fastq_urls[0]}" == "NA" ]]; then
            echo_err \
                "missing FASTQ URL for TSV line '${iter}'" \
                "('${run_acc}')."
            exit 1
        fi

        list_url_1+=( "${fastq_urls[0]}" )
        list_url_2+=( "NA" )  # No second URL for single-end data

    #  For paired-end data
    elif [[ ${#fastq_urls[@]} -eq 2 ]]; then
        if [[ -z "${fastq_urls[0]}" || -z "${fastq_urls[1]}" ]]; then
            echo_err \
                "missing one or more paired-end FASTQ URLs for TSV line" \
                "'${iter}' ('${run_acc}')."
            exit 1
        fi

        list_url_1+=( "${fastq_urls[0]}" )
        list_url_2+=( "${fastq_urls[1]}" )
    else
        echo_err "unexpected number of FASTQ URLs for '${run_acc}'."
        exit 1
    fi
done < "${infile}"

check_arr_nonempty "list_acc"   "download entries"
check_arr_nonempty "list_url_1" "FASTQ URL entries"
check_arr_nonempty "list_url_2" "second FASTQ URL entries"
check_arr_nonempty "list_cus"   "custom-name entries"

check_arr_lengths "list_acc" "list_url_1"
check_arr_lengths "list_acc" "list_url_2"
check_arr_lengths "list_acc" "list_cus"

#  Report array element assignments if in "verbose mode"
if [[ "${verbose}" == "true" ]]; then
    echo
    echo  #TODO: switch to using 'print_banner_pretty'
    echo "#######################"
    echo "## list_acc elements ##"
    echo "#######################"
    echo
    for el in "${list_acc[@]}"; do echo "${el}"; done
    echo
    echo

    echo "#######################"
    echo "## list_cus elements ##"
    echo "#######################"
    echo
    for el in "${list_cus[@]}"; do echo "${el}"; done
    echo
    echo

    echo "#########################"
    echo "## list_url_1 elements ##"
    echo "#########################"
    echo
    for el in "${list_url_1[@]}"; do echo "${el}"; done
    echo
    echo

    echo "#########################"
    echo "## list_url_2 elements ##"
    echo "#########################"
    echo
    for el in "${list_url_2[@]}"; do echo "${el}"; done
    echo
    echo

    unset el
fi

#  Create a command file for GNU Parallel-backed execution if needed
if [[ "${threads}" -gt 1 ]]; then
    config="${err_out}/${nam_job}.config_parallel.txt"

    if [[ -f "${config}" ]]; then rm "${config}"; fi

    for i in "${!list_acc[@]}"; do
        build_cmd "${i}"

        {
            print_built_cmd
        } >> "${config}" || {
            echo_err "failed to write command, index no. '${i}'."
            exit 1
        }
    done
fi

#  Execute download and symlink creation based on Slurm flag and thread count
if [[ "${slurm}" == "true" ]]; then
    if [[ "${threads}" -gt 1 ]]; then
        unset cmd_para  && declare -a cmd_para
        unset cmd_slurm && declare -a cmd_slurm
        cmd_wrap=""

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
        exit 1
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
        fi
    else
        for i in "${!list_acc[@]}"; do
            build_cmd "${i}"

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
    fi
fi
