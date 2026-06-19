#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: execute_align_fastqs.sh
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
    help/help_execute_align_fastqs \
    manage_parallel \
    process_sequences \
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

    unset cmd_bld && declare -ga cmd_bld

    show_help=$(cat << EOM
Usage:
  build_cmd [-h|--hlp|--help] [idx]

Description:
  Construct the command array 'cmd_bld' for one call to 'submit_align_fastqs.sh'.

Positional argument:
  1  idx  <int|UNSET>
    Optional zero-based index into 'arr_infile'.

    If omitted or set to 'UNSET', construct a non-indexed command using the full serialized 'csv_infile' string.

    If supplied, construct a per-entry command using 'arr_infile[idx]'.

Expected globals:
  scr_sub env_nam dir_scr threads aligner bt2_aln bwa_alg mapq req_flg index ref_fa out_ext csv_infile dir_out qname sfx_se sfx_pe err_out nam_job

Returns:
  0 if 'cmd_bld' is constructed successfully; otherwise 1.

Notes:
  - 'cmd_bld' is written as a global indexed array.
  - The reference FASTA argument is included only when 'out_ext=cram'.
  - Flag-only options are added as separate array elements.
  - On index handling:

    idx=UNSET  ->  --csv_infile "\${csv_infile}"       # Full serialized list
    idx=0..n   ->  --csv_infile "\${arr_infile[idx]}"  # One FASTQ entry
EOM
    )

    if [[ "${idx}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${idx}" ]]; then
        idx="UNSET"
    fi

    if [[ "${idx}" == "UNSET" ]]; then
        #  Use the full serialized input list for Slurm or whole-wrapper calls.
        infile_i="${csv_infile}"
    else
        #  Use one parsed FASTQ entry for per-sample local/parallel calls.
        check_int_nonneg "${idx}" "idx" || return 1
        infile_i="${arr_infile[idx]}"
    fi

    cmd_bld=(
        "${scr_sub}"
            --env_nam "${env_nam}"
            --dir_scr "${dir_scr}"
            --threads "${threads}"
            --aligner "${aligner}"
            --mapq "${mapq}"
            --index "${index}"
            --csv_infile "${infile_i}"
            --dir_out "${dir_out}"
            --out_ext "${out_ext}"
            --sfx_se "${sfx_se}"
            --sfx_pe "${sfx_pe}"
            --err_out "${err_out}"
            --nam_job "${nam_job}"
    )

    case "${aligner}" in
        bowtie2) cmd_bld+=( --bt2_aln "${bt2_aln}" ) ;;
        bwa)     cmd_bld+=( --bwa_alg "${bwa_alg}" ) ;;
    esac

    if [[ "${req_flg}" == "true" ]]; then
        cmd_bld+=( --req_flg )
    fi

    if [[ "${out_ext}" == "cram" ]]; then
        cmd_bld+=( --ref "${ref_fa}" )
    fi

    if [[ "${qname}" == "true" ]]; then
        cmd_bld+=( --qname )
    fi
}


#  Initialize hardcoded argument variables
function init_args_hardcoded() {
    env_nam="env_protocol"
    scr_sub="${dir_scr}/submit_align_fastqs.sh"
    par_job=""
}


#  Initialize argument variables, assigning default values where applicable
function init_arg_defs() {
    verbose=false
    dry_run=false
    threads=1
    aligner="bowtie2"
    bt2_aln="global"
    bwa_alg="mem"
    ref_fa=""
    out_ext="bam"
    mapq=1
    req_flg=false
    index=""
    csv_infile=""
    dir_out=""
    qname=false
    sfx_se=".atria.fastq.gz"
    sfx_pe="_R1.atria.fastq.gz"
    err_out=""
    nam_job="align_fastqs"
    max_job=6
    slurm=false
    time="1:00:00"

    unset arr_infile cmd_bld
    declare -ga arr_infile cmd_bld
}


#  Initialize hardcoded arguments and user-facing argument defaults
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
                    help_execute_align_fastqs >&2
                    return 1
                }
                threads="${2}"
                shift 2
                ;;

            -a|--aln|--aligner)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_align_fastqs >&2
                    return 1
                }
                aligner="${2,,}"
                shift 2
                ;;

            #TODO: not fully aligned with 'submit_align_fastqs'
            -2a|-bn|--bt2[_-]aln)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_align_fastqs >&2
                    return 1
                }
                bt2_aln="${2,,}"
                shift 2
                ;;

            -ba|--bwa[_-]alg)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_align_fastqs >&2
                    return 1
                }
                bwa_alg="${2,,}"
                shift 2
                ;;

            -r|--ref)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_align_fastqs >&2
                    return 1
                }
                ref_fa="${2}"
                shift 2
                ;;

            -mq|--mapq)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_align_fastqs >&2
                    return 1
                }
                mapq="${2}"
                shift 2
                ;;

            -rf|--req[_-]flg)
                req_flg=true
                shift 1
                ;;

            -ix|--index)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_align_fastqs >&2
                    return 1
                }
                index="${2}"
                shift 2
                ;;

            -i|-fi|-ci|--infile|--infiles|--fil[_-]in|--csv[_-]infile|--csv[_-]infiles)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_align_fastqs >&2
                    return 1
                }
                csv_infile="${2}"
                shift 2
                ;;

            -do|--dir[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_align_fastqs >&2
                    return 1
                }
                dir_out="${2}"
                shift 2
                ;;

            -ox|--out[_-]ext)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_align_fastqs >&2
                    return 1
                }
                out_ext="${2,,}"
                shift 2
                ;;

            -qn|--qname)
                qname=true
                shift 1
                ;;

            -sxs|--sfx[_-]se|--suffix[_-]se)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_align_fastqs >&2
                    return 1
                }
                sfx_se="${2}"
                shift 2
                ;;

            -sxp|--sfx[_-]pe|--suffix[_-]pe)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_align_fastqs >&2
                    return 1
                }
                sfx_pe="${2}"
                shift 2
                ;;

            -eo|--err[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_align_fastqs >&2
                    return 1
                }
                err_out="${2}"
                shift 2
                ;;

            -nj|--nam[_-]job)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_align_fastqs >&2
                    return 1
                }
                nam_job="${2}"
                shift 2
                ;;

            -mj|--max[_-]job)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_execute_align_fastqs >&2
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
                    help_execute_align_fastqs >&2
                    return 1
                }
                time="${2}"
                shift 2
                ;;

            *)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                help_execute_align_fastqs >&2
                return 1
                ;;
        esac
    done
}


#  Validate scalar arguments and assign derived scalar defaults
function validate_args() {
    validate_var "env_nam" "${env_nam}" || return 1
    check_env_installed "${env_nam}" || return 1

    validate_var_dir  "dir_scr" "${dir_scr}" 0 false || return 1

    validate_var_file "scr_sub" "${scr_sub}" || return 1

    validate_var "threads" "${threads}" || return 1
    check_int_pos "${threads}" "threads" || return 1

    case "${aligner}" in
        bowtie2)
            case "${bt2_aln}" in
                local|global|end-to-end) : ;;
                *)
                    echo_err \
                        "selection associated with '--bt2_aln' is not valid:" \
                        "'${bt2_aln}'. Selection must be 'local', 'global'," \
                        "or 'end-to-end'."
                    return 1
                    ;;
            esac
            ;;

        bwa)
            case "${bwa_alg}" in
                mem|aln) : ;;
                *)
                    echo_err \
                        "selection associated with '--bwa_alg' is not valid:" \
                        "'${bwa_alg}'. Selection must be 'mem' or 'aln'."
                    return 1
                    ;;
            esac
            ;;

        bwa-mem2)
            if [[ "${bwa_alg}" != "mem" ]]; then
                echo_err \
                    "'--bwa_alg' must be 'mem' when '--aligner bwa-mem2'."
                return 1
            fi
            ;;

        *)
            echo_err \
                "selection associated with '--aligner' is not valid:" \
                "'${aligner}'. Selection must be 'bowtie2', 'bwa', or" \
                "'bwa-mem2'."
            return 1
            ;;
    esac

    validate_var "mapq" "${mapq}" || return 1
    check_int_nonneg "${mapq}" "mapq" || return 1

    case "${out_ext}" in
        bam|cram) : ;;
        *)
            echo_err \
                "'--out_ext' must be 'bam' or 'cram': '${out_ext}'."
            return 1
            ;;
    esac

    if [[ "${out_ext}" == "cram" ]]; then
        validate_var_file "ref_fa" "${ref_fa}" || return 1
    fi

    validate_var "index" "${index}" || return 1
    validate_var_dir "index parent directory" "$(dirname "${index}")" 0 false \
        || return 1

    if [[ "${aligner}" =~ ^(bwa|bwa-mem2)$ ]]; then
        validate_var_file "index" "${index}" || return 1
    fi

    validate_var "csv_infile" "${csv_infile}" || return 1
    validate_var_dir "csv_infile parent directory" \
        "$(dirname "${csv_infile%%[,;]*}")" 0 false || return 1

    validate_var_dir "dir_out" "${dir_out}" || return 1

    validate_var "sfx_se" "${sfx_se}" || return 1
    validate_var "sfx_pe" "${sfx_pe}" || return 1

    if [[ -z "${err_out}" ]]; then err_out="${dir_out}/logs"; fi
    validate_var_dir "err_out" "${err_out}" || return 1

    validate_var "nam_job" "${nam_job}" || return 1
}


#  Parse input FASTQ vector
function prepare_vecs() {
    unset arr_infile && declare -ga arr_infile
    IFS=';' read -r -a arr_infile <<< "${csv_infile}"
}


#  Validate parsed FASTQ vector
function validate_vecs() {
    check_arr_nonempty "arr_infile" "csv_infile" || return 1
    check_string_fastqs "${csv_infile}" "${sfx_se}" "${sfx_pe}" || return 1
}


#  Configure local, GNU Parallel, or Slurm execution
function config_exec() {
    validate_var "max_job" "${max_job}" || return 1
    check_int_pos "${max_job}" "max_job" || return 1

    if [[ "${slurm}" == "true" ]]; then
        max_job="$(reset_max_job "${max_job}" "${#arr_infile[@]}")"

        validate_var "time" "${time}" || return 1
        check_format_time "${time}" || return 1
    elif [[ "${max_job}" -le 1 ]]; then
        #  Serial local execution does not require parallel job detection
        par_job=1
        unset max_job time
    else
        IFS=';' read -r threads par_job < <(
            set_params_parallel "${threads}" "${max_job}"
        ) || return 1
        unset max_job time

        validate_var "par_job" "${par_job}" || return 1
        check_int_pos "${par_job}" "par_job" || return 1
    fi

    if [[ "${verbose}" == "true" ]]; then
        print_parallel_info \
            "${slurm}" "${max_job:-UNSET}" "${par_job:-UNSET}" "${threads}" \
            "arr_infile" \
            || return 1
    fi
}


#  Activate environment
function setup_env() {
    local -a env_msg

    env_msg=(
        "'handle_env' failed for 'env_nam=${env_nam}'. Check that"
        "Conda/Mamba are available and that the environment exists."
    )

    if [[ "${verbose}" == "true" ]]; then
        if ! handle_env "${env_nam}"; then
            echo_err "${env_msg[*]}"
            return 1
        fi

        echo
    else
        if ! handle_env "${env_nam}" > /dev/null 2>&1; then
            echo_err "${env_msg[*]}"
            return 1
        fi
    fi
}


#  Check tools needed by the selected dispatch mode
function check_tools() {
    case "${aligner}" in
        bowtie2)  check_pgrm_path bowtie2  || return 1 ;;
        bwa)      check_pgrm_path bwa      || return 1 ;;
        bwa-mem2) check_pgrm_path bwa-mem2 || return 1 ;;
    esac

    check_pgrm_path samtools || return 1

    if [[ "${slurm}" == "true" ]]; then
        check_pgrm_path sbatch || return 1
    elif [[ "${par_job}" -gt 1 ]]; then
        check_pgrm_path parallel || return 1
    fi
}


#  Report argument variable assignments if in verbose mode
function print_state_debug() {
    if [[ "${verbose}" == "true" ]]; then
        print_banner_pretty "Hardcoded variable assignments"
        echo
        echo "env_nam=${env_nam}"
        echo "dir_scr=${dir_scr}"
        echo "scr_sub=${scr_sub}"
        echo "par_job=${par_job:-UNSET}"
        echo
        echo

        print_banner_pretty "Argument variable assignments"
        echo
        echo "verbose=${verbose}"
        echo "dry_run=${dry_run}"
        echo "threads=${threads}"
        echo "aligner=${aligner}"
        echo "bt2_aln=${bt2_aln}"
        echo "bwa_alg=${bwa_alg}"
        echo "ref_fa=${ref_fa}"
        echo "out_ext=${out_ext}"
        echo "mapq=${mapq}"
        echo "req_flg=${req_flg}"
        echo "index=${index}"
        echo "csv_infile=${csv_infile}"
        echo "dir_out=${dir_out}"
        echo "qname=${qname}"
        echo "sfx_se=${sfx_se}"
        echo "sfx_pe=${sfx_pe}"
        echo "err_out=${err_out}"
        echo "nam_job=${nam_job}"
        echo "max_job=${max_job:-UNSET}"
        echo "slurm=${slurm}"
        echo "time=${time:-UNSET}"
        echo
        echo

        print_banner_pretty "Array derived from variable"
        echo
        echo "arr_infile=( ${arr_infile[*]} )"
        echo
        echo
    fi
}


#  Dispatch Slurm, GNU Parallel, or serial jobs
function run_jobs() {
    local config idx log_err log_out
    local -a cmd_slurm

    if [[ "${slurm}" == "true" ]]; then
        #  Slurm execution
        build_cmd "UNSET" || return 1

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
            "${cmd_slurm[@]}" || return 1
        fi
    else
        #  Non-Slurm execution: GNU Parallel ('par_job > 1') or serial
        #+ ('par_job == 1')
        if [[ "${par_job}" -gt 1 ]]; then
            config="${err_out}/${nam_job}.config_parallel.txt"

            if [[ -f "${config}" ]]; then rm "${config}"; fi

            for idx in "${!arr_infile[@]}"; do
                build_cmd "${idx}" || return 1
                cmd_bld=( "${BASH}" "${cmd_bld[@]}" )

                IFS=';' read -r log_out log_err < <(
                    get_submit_logs "${arr_infile[idx]}"
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
                parallel --jobs "${par_job}" --dryrun < "${config}" \
                    || return 1
                echo
                echo
            fi

            if [[ "${dry_run}" == "false" ]]; then
                parallel --jobs "${par_job}" < "${config}" || return 1
            fi
        else
            #  Serial execution
            build_cmd "UNSET" || return 1
            cmd_bld=( "${BASH}" "${cmd_bld[@]}" )

            log_out="${err_out}/${nam_job}_ser.stdout.txt"
            log_err="${err_out}/${nam_job}_ser.stderr.txt"

            if [[ "${dry_run}" == "true" || "${verbose}" == "true" ]]; then
                print_banner_pretty "Serial execution"
                echo
                print_built_cmd "${log_out}" "${log_err}"
                echo
            fi

            if [[ "${dry_run}" == "false" ]]; then
                "${cmd_bld[@]}" >> "${log_out}" 2>> "${log_err}" \
                    || return 1
            fi
        fi
    fi
}


#  Main script execution
function main() {
    init_defs

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        help_execute_align_fastqs >&2
        return 0
    fi

    parse_args "$@"   || return 1
    validate_args     || return 1
    prepare_vecs      || return 1
    validate_vecs     || return 1
    config_exec       || return 1
    setup_env         || return 1
    check_tools       || return 1
    print_state_debug || return 1
    run_jobs          || return 1
}


main "$@"
