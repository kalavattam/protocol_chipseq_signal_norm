#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: execute_filter_bams.sh
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
    help/help_execute_filter_bams \
    manage_parallel \
    wrap_cmd \
    || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source required helper scripts." >&2
        exit 1
    }

unset fnc_src


function build_cmd() {
    local idx="${1:-UNSET}"
    local infile_i=""
    local show_help

    unset cmd_bld && declare -ag cmd_bld

    show_help=$(cat << EOM
Usage:
  build_cmd [-h|--hlp|--help] [idx]

Description:
  Construct the command array 'cmd_bld' for one call to 'submit_filter_bams.sh'.

Positional argument:
  1  idx  <int|UNSET>  Optional zero-based index into 'arr_infile'.

                       If omitted or set to 'UNSET', construct a non-indexed command using the full serialized 'csv_infile' string.

                       If supplied, construct a per-entry command using 'arr_infile[idx]'.

Expected globals:
  scr_sub env_nam dir_scr retain threads csv_infile dir_out mito tg mtr chk_chr err_out nam_job

Returns:
  0 if 'cmd_bld' is constructed successfully; otherwise 1.

Notes:
  - 'cmd_bld' is written as a global indexed array.
  - Flag-only options are added as separate array elements.
  - On index handling:

    idx=UNSET  ->  --csv_infile "\${csv_infile}"          # Full serialized list
    idx=0..n   ->  --csv_infile "\${arr_infile[idx]}"  # One BAM entry
EOM
    )

    if [[ "${idx}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${idx}" ]]; then
        idx="UNSET"
    fi

    if [[ "${idx}" == "UNSET" ]]; then
        #  Use the full serialized input list for Slurm or whole-wrapper calls
        infile_i="${csv_infile}"
    else
        #  Use one parsed BAM path for per-sample local/parallel calls
        check_int_nonneg "${idx}" "idx" || return 1
        infile_i="${arr_infile[idx]}"
    fi

    # shellcheck disable=SC2034
    cmd_bld=(
        "${scr_sub}"
            --env_nam "${env_nam}"
            --dir_scr "${dir_scr}"
            --retain "${retain}"
            --threads "${threads}"
            --csv_infile "${infile_i}"
            --dir_out "${dir_out}"
            --err_out "${err_out}"
            --nam_job "${nam_job}"
    )

    if [[ "${mito}" == "true" ]]; then
        cmd_bld+=( --mito )
    fi

    if [[ "${tg}" == "true" ]]; then
        cmd_bld+=( --tg )
    fi

    if [[ "${mtr}" == "true" ]]; then
        cmd_bld+=( --mtr )
    fi

    if [[ "${chk_chr}" == "true" ]]; then
        cmd_bld+=( --chk_chr )
    fi
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_protocol"
scr_sub="${dir_scr}/submit_filter_bams.sh"
par_job=""

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=1
csv_infile=""
dir_out=""
retain="sc"
mito=false
tg=false
mtr=false
chk_chr=false
err_out=""
nam_job="filter_bams"
max_job=6
slurm=false
time="0:30:00"

#  Parse arguments
if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
    help_execute_filter_bams >&2
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
                help_execute_filter_bams >&2
                exit 1
            }
            threads="${2}"
            shift 2
            ;;

        -i|-fi|-ci|--infile|--infiles|--fil[_-]in|--csv[_-]infile|--csv[_-]infiles)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_filter_bams >&2
                exit 1
            }
            csv_infile="${2}"
            shift 2
            ;;

        -do|--dir[_-]out)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_filter_bams >&2
                exit 1
            }
            dir_out="${2}"
            shift 2
            ;;

        -rt|--retain)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_filter_bams >&2
                exit 1
            }
            retain="${2,,}"
            shift 2
            ;;

        -m|--mito)
            mito=true
            shift 1
            ;;

        -tg|--tg)
            tg=true
            shift 1
            ;;

        -mr|--mtr)
            mtr=true
            shift 1
            ;;

        -cc|--chk[_-]chr)
            chk_chr=true
            shift 1
            ;;

        -eo|--err[_-]out)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_filter_bams >&2
                exit 1
            }
            err_out="${2}"
            shift 2
            ;;

        -nj|--nam[_-]job)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_filter_bams >&2
                exit 1
            }
            nam_job="${2}"
            shift 2
            ;;

        -mj|--max[_-]job)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_execute_filter_bams >&2
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
                help_execute_filter_bams >&2
                exit 1
            }
            time="${2}"
            shift 2
            ;;

        *)
            echo_err "unknown option/parameter passed: '${1}'."
            echo >&2
            help_execute_filter_bams >&2
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

validate_var "csv_infile" "${csv_infile}"
validate_var_dir "csv_infile parent directory" \
    "$(dirname "${csv_infile%%,*}")" 0 false
check_str_delim "csv_infile" "${csv_infile}"

validate_var_dir "dir_out" "${dir_out}"

case "${retain}" in
    sc|sp) : ;;
    *)
        echo_err \
            "selection associated with '--retain' is not valid: '${retain}'." \
            "Selection must be 'sc' or 'sp'."
        exit 1
        ;;
esac

if [[ "${retain}" == "sc" ]]; then
    if [[ "${tg}" == "true" && "${mtr}" == "true" ]]; then
        echo_warn \
            "flags '--tg' and '--mtr' were supplied with '--retain sc' and" \
            "will be ignored."
        tg=false
        mtr=false
    elif [[ "${tg}" == "true" ]]; then
        echo_warn \
            "flag '--tg' was supplied with '--retain sc' and will be" \
            "ignored."
        tg=false
    elif [[ "${mtr}" == "true" ]]; then
        echo_warn \
            "flag '--mtr' was supplied with '--retain sc' and will be" \
            "ignored."
        mtr=false
    fi
fi

if [[ -z "${err_out}" ]]; then err_out="${dir_out}/logs"; fi
validate_var_dir "err_out" "${err_out}"

validate_var "nam_job" "${nam_job}"


#  Parse and validate infiles -------------------------------------------------
IFS=',' read -r -a arr_infile <<< "${csv_infile}"

check_arr_nonempty "arr_infile" "csv_infile"

for infile in "${arr_infile[@]}"; do
    validate_var_file "infile" "${infile}"
done
unset infile


#  Parse job execution parameters ---------------------------------------------
validate_var "max_job" "${max_job}"
check_int_pos "${max_job}" "max_job"

if [[ "${slurm}" == "true" ]]; then
    max_job="$(reset_max_job "${max_job}" "${#arr_infile[@]}")"

    validate_var "time" "${time}"
    check_format_time "${time}"
else
    IFS=';' read -r threads par_job < <(
        set_params_parallel "${threads}" "${max_job}"
    )
    unset max_job time

    validate_var "par_job" "${par_job}"
    check_int_pos "${par_job}" "par_job"
fi

#  Debug parallelization information
if [[ "${verbose}" == "true" ]]; then
    print_parallel_info \
        "${slurm}" "${max_job:-UNSET}" "${par_job:-UNSET}" "${threads}" \
        "arr_infile"
fi


#  Activate environment and check that dependencies are in PATH ---------------
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

check_pgrm_path awk
check_pgrm_path grep
check_pgrm_path mv
check_pgrm_path rm
check_pgrm_path samtools

if [[ "${slurm}" == "true" ]]; then
    check_pgrm_path sbatch
elif [[ ${par_job} -gt 1 ]]; then
    check_pgrm_path parallel
fi


#  Do the main work ===========================================================
if [[ "${verbose}" == "true" ]]; then
    #TODO: incorporate 'print_banner_pretty' where appropriate
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo
    echo "env_nam=${env_nam}"
    echo "dir_scr=${dir_scr}"
    echo "scr_sub=${scr_sub}"
    echo "par_job=${par_job:-UNSET}"
    echo
    echo
    echo "###################################"
    echo "## Argument variable assignments ##"
    echo "###################################"
    echo
    echo "verbose=${verbose}"
    echo "dry_run=${dry_run}"
    echo "threads=${threads}"
    echo "csv_infile=${csv_infile}"
    echo "dir_out=${dir_out}"
    echo "retain=${retain}"
    echo "mito=${mito}"
    echo "tg=${tg}"
    echo "mtr=${mtr}"
    echo "chk_chr=${chk_chr}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo "max_job=${max_job:-UNSET}"
    echo "slurm=${slurm}"
    echo "time=${time:-UNSET}"
    echo
    echo
    echo "#################################"
    echo "## Array derived from variable ##"
    echo "#################################"
    echo
    echo "arr_infile=( ${arr_infile[*]} )"
    echo
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
            --error="${err_out}/${nam_job}.%A-%a.stderr.txt"
            --output="${err_out}/${nam_job}.%A-%a.stdout.txt"
            --array="1-${#arr_infile[@]}%${max_job}"
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

        for idx in "${!arr_infile[@]}"; do
            build_cmd "${idx}"

            IFS=';' read -r log_out log_err < <(
                get_submit_logs "${arr_infile[idx]}"
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
        #  Serial execution
        build_cmd "UNSET"

        log_out="${err_out}/${nam_job}_ser.stdout.txt"
        log_err="${err_out}/${nam_job}_ser.stderr.txt"

        if [[ "${dry_run}" == "true" || "${verbose}" == "true" ]]; then
            print_banner_pretty "Serial execution"
            echo
            print_built_cmd "${log_out}" "${log_err}"
            echo
        fi

        if [[ "${dry_run}" == "false" ]]; then
            "${cmd_bld[@]}" >> "${log_out}" 2>> "${log_err}"
        fi
    fi
fi
