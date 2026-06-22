#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: submit_convert_bam_bed.sh
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


#  Initialize hardcoded argument variables
function init_args_hardcoded() {
    #  If true, run script in debug mode
    debug=true
}


#  Initialize argument variables, assigning default values where applicable
function init_arg_defs() {
    env_nam="env_protocol"
    dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
    threads=1
    csv_infile=""
    pth_scr_py=""
    dir_out=""
    dir_eo=""
    nam_job="convert_bam_bed"
    use_awk=false
    pth_scr_py_set=false
    ref_fa=""

    dir_fnc=""
    dir_rep=""

    unset arr_infiles
    declare -ga arr_infiles
}


#  Initialize hardcoded arguments and user-facing argument defaults
function init_defs() {
    init_args_hardcoded
    init_arg_defs
}


#  Print main usage text
function show_help_main() {
    cat << EOM
Usage:
  submit_convert_bam_bed.sh
    [--help] [--env_nam <str>] [--dir_scr <dir>] [--threads <int>]
    --csv_infile <csv:file> [--ref_fa <file>]
    [(--pth_scr_py <file> | --use_awk)]
    --dir_out <dir> --dir_eo <dir> [--nam_job <str>]

Description:
  Convert BAM or CRAM input files to BED output files.

  When run inside a Slurm array task, this script processes the indexed input selected by 'SLURM_ARRAY_TASK_ID'. Otherwise, it processes every input file serially in the current shell.

Keyword arguments:
  -h, --hlp, --help  <flag>
    Display this help message and exit.

  -en, --env, --env_nam  <str>
    Conda/Mamba environment to activate (default: '${env_nam}').

  -ds, --dir_scr  <dir>
    Directory containing scripts and shared shell functions (default: '${dir_scr}').

  -t, --thr, --threads  <int>
    Number of threads to use per task (default: ${threads}).

  -i, -fi, --fil_in, --csv_infile, --csv_infiles  <csv:file>
    Comma-separated list of BAM or CRAM input files.

  -r, --ref, --ref_fa, --reference  <file>
    Reference FASTA. Required when any input file is CRAM.

  -pp, --pth_scr_py  <file>
    Path to a custom Python BAM/CRAM-to-BED script (default: '\${dir_scr}/compute_signal.py'). Mutually exclusive with '--use_awk'.

  -awk, --use_awk  <flag>
    Run AWK processing code rather than the Python script. Mutually exclusive with explicit '--pth_scr_py'. Do not use with single-end data.

  -do, --dir_out  <dir>
    Directory to save BED output files.

  -eo, --dir_eo, --err_out  <dir>
    Directory for stderr and stdout files.

  -nj, --nam, --name, --nam_job, --name_job  <str>
    Name of job (default: '${nam_job}').

Dependencies:
  Recommended environment:
    - env_protocol

  External programs:
    - samtools
    - awk
    - sort
    - gzip
    - python

  Python scripts:
    - compute_signal.py or compatible BAM/CRAM-to-BED script

Notes:
  - The AWK and Python branches are not equivalent:
    + The AWK branch assumes QNAME-sorted paired-end records in adjacent pairs and writes paired-fragment intervals.
    + The Python branch uses 'compute_signal.py', which does not require QNAME sorting, handles single-end data, and emits fragments according to its own 'parse_bam()' policy.
  - Python conversion is the default. Use '--pth_scr_py' only to supply a custom converter, or use '--use_awk' to select the AWK branch.
  - CRAM inputs require '--ref_fa' in both AWK and Python branches.
  - This is a submit-wrapper script: it supports serial/local iteration and Slurm-array task execution, but it does not submit Slurm jobs itself.
EOM
}


#  Resolve '--dir_scr' before sourced parser helpers are available
function resolve_dir_scr() {
    local scr="${1:-}"
    shift

    local -a args=( "$@" )
    local i=0

    if [[ -z "${scr}" ]]; then
        scr="unknown_script"
    fi

    for (( i = 0; i < ${#args[@]}; i++ )); do
        case "${args[i]}" in
            -ds|--dir[_-]scr)
                if (( i + 1 >= ${#args[@]} )) \
                    || [[ -z "${args[i + 1]:-}" || "${args[i + 1]}" == -* ]]
                then
                    echo "error(${scr}):" \
                        "option '${args[i]}' requires a value." >&2
                    echo >&2
                    show_help_main
                    return 1
                fi

                printf "%s\n" "${args[i + 1]}"
                return 0
                ;;
        esac
    done

    printf "%s\n" "${dir_scr}"
}


#  Source 'source_helpers.sh' and requested helper scripts from 'dir_scr'
function source_helpers_submit() {
    local scr="${1:-}"
    local dir_scr_arg="${2:-}"
    local fnc_src

    if (( $# < 2 )); then
        echo "error(${scr:-unknown_script}):" \
            "expected at least two arguments: 'scr' and 'dir_scr_arg'." >&2
        return 1
    fi

    shift 2

    if [[ -z "${scr}" ]]; then
        scr="unknown_script"
    fi

    if [[ -z "${dir_scr_arg}" ]]; then
        echo "error(${scr}):" \
            "positional argument 2, 'dir_scr_arg', is missing." >&2
        return 1
    elif [[ ! -d "${dir_scr_arg}" ]]; then
        echo "error(${scr}):" \
            "script directory not found: '${dir_scr_arg}'." >&2
        return 1
    elif (( $# < 1 )); then
        echo "error(${scr}):" \
            "at least one helper script name must be supplied." >&2
        return 1
    fi

    dir_fnc="${dir_scr_arg}/functions"
    fnc_src="${dir_fnc}/source_helpers.sh"

    if [[ ! -f "${fnc_src}" ]]; then
        echo "error(${scr}): script not found: '${fnc_src}'." >&2
        return 1
    fi

    # shellcheck disable=SC1090
    source "${fnc_src}" || {
        echo "error(${scr}): failed to source '${fnc_src}'." >&2
        return 1
    }

    source_helpers "${dir_fnc}" "$@" || {
        echo "error(${scr}): failed to source required helper scripts." >&2
        return 1
    }
}


#  Parse keyword arguments after helper scripts have been sourced
function parse_args() {
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -h|--hlp|--help)
                show_help_main
                echo >&2
                return 0
                ;;

            -en|--env|--env[_-]nam)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                env_nam="${2}"
                shift 2
                ;;

            -ds|--dir[_-]scr)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                dir_scr="${2}"
                shift 2
                ;;

            -t|--thr|--threads)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                threads="${2}"
                shift 2
                ;;

            -i|-fi|--fil[_-]in|--csv[_-]infile|--csv[_-]infiles)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                csv_infile="${2}"
                shift 2
                ;;

            -r|--ref|--ref[_-]fa|--reference)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                ref_fa="${2}"
                shift 2
                ;;

            -pp|--pth[_-]scr[_-]py)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                pth_scr_py="${2}"
                pth_scr_py_set=true
                shift 2
                ;;

            -do|--dir[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                dir_out="${2}"
                shift 2
                ;;

            -eo|--dir[_-]eo|--err[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                dir_eo="${2}"
                shift 2
                ;;

            -nj|--nam|--name|--nam[_-]job|--name_job)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                nam_job="${2}"
                shift 2
                ;;

            -awk|--use[_-]awk)
                use_awk=true
                shift 1
                ;;

            *)
                echo_err "unknown parameter passed: '${1}'."
                echo >&2
                show_help_main
                return 1
                ;;
        esac
    done
}


#  Resolve script paths after 'dir_scr' has been parsed
function resolve_paths_scrs() {
    validate_var_dir "dir_scr" "${dir_scr}" 0 false || return 1

    if [[ -z "${pth_scr_py}" ]]; then
        pth_scr_py="${dir_scr}/compute_signal.py"
    fi
}


#  Validate required arguments and scalar values
function validate_args() {
    validate_var     "env_nam"    "${env_nam}"         || return 1
    validate_var_dir "dir_scr"    "${dir_scr}" 0 false || return 1
    validate_var     "threads"    "${threads}"         || return 1
    validate_var     "csv_infile" "${csv_infile}"      || return 1
    validate_var_dir "dir_out"    "${dir_out}"         || return 1
    validate_var_dir "dir_eo"     "${dir_eo}"          || return 1
    validate_var     "nam_job"    "${nam_job}"         || return 1

    check_int_pos "${threads}" "threads" || return 1

    if [[ "${use_awk}" == "true" && "${pth_scr_py_set}" == "true" ]]; then
        echo_err \
            "arguments '--pth_scr_py' and '--use_awk' are mutually" \
            "exclusive."
        return 1
    fi

    if [[ "${use_awk}" != "true" ]]; then
        validate_var_file "pth_scr_py" "${pth_scr_py}" || return 1
    fi
}


#  Print parsed scalar arguments when debugging is enabled
function print_state_debug() {
    if [[ "${debug}" != "true" ]]; then
        return 0
    fi

    debug_var \
        "env_nam=${env_nam}" \
        "dir_scr=${dir_scr}" \
        "threads=${threads}" \
        "csv_infile=${csv_infile}" \
        "pth_scr_py=${pth_scr_py}" \
        "dir_out=${dir_out}" \
        "dir_eo=${dir_eo}" \
        "nam_job=${nam_job}" \
        "use_awk=${use_awk}" \
        "ref_fa=${ref_fa:-UNSET}"
}


#  Activate the requested Conda/Mamba environment
function setup_env() {
    handle_env "${env_nam}" || return 1

    dir_rep="$(cd "${dir_scr}/.." && pwd)"
    if [[ -z "${PYTHONPATH:-}" ]]; then
        export PYTHONPATH="${dir_rep}"
    else
        export PYTHONPATH="${dir_rep}:${PYTHONPATH}"
    fi
}


#  Check runtime tool dependencies
function check_tools() {
    if [[ "${use_awk}" == "true" ]]; then
        check_pgrm_path samtools || return 1
        check_pgrm_path awk      || return 1
        check_pgrm_path sort     || return 1
        check_pgrm_path gzip     || return 1
    else
        check_pgrm_path python   || return 1
    fi
}


#  Reconstruct input alignment vector
function prepare_vecs() {
    unset arr_infiles && declare -ga arr_infiles
    IFS=',' read -r -a arr_infiles <<< "${csv_infile}"
}


#  Validate reconstructed input alignment vector
function validate_vecs() {
    local idx infile need_ref=false

    check_arr_nonempty "arr_infiles" "csv_infile" || return 1

    for idx in "${!arr_infiles[@]}"; do
        infile="${arr_infiles[idx]}"
        validate_var_file "arr_infiles" "${infile}" "${idx}" || return 1

        case "${infile,,}" in
            *.bam) : ;;
            *.cram) need_ref=true ;;
            *)
                echo_err \
                    "input file must end in '.bam' or '.cram': '${infile}'."
                return 1
                ;;
        esac
    done

    if [[ "${need_ref}" == "true" && -z "${ref_fa}" ]]; then
        echo_err "'--ref_fa' is required when '--csv_infile' contains CRAM."
        return 1
    fi

    if [[ -n "${ref_fa}" ]]; then
        validate_var_file "ref_fa" "${ref_fa}" || return 1
    fi
}


#  Print reconstructed array values when debugging is enabled
function print_vecs_debug() {
    if [[ "${debug}" != "true" ]]; then
        return 0
    fi

    debug_var "\${#arr_infiles[@]}=${#arr_infiles[@]}"
    echo "arr_infiles=( ${arr_infiles[*]} )" >&2 && echo >&2
}


#  Strip known alignment suffixes before assigning a BED output path
function strip_sfx_aln() {
    local base="${1:-}"

    case "${base,,}" in
        *.qnam.bam)  printf '%s\n' "${base%.qnam.bam}"  ;;
        *.qnam.cram) printf '%s\n' "${base%.qnam.cram}" ;;
        *.bam)       printf '%s\n' "${base%.bam}"       ;;
        *.cram)      printf '%s\n' "${base%.cram}"      ;;
        *)           printf '%s\n' "${base}"            ;;
    esac
}


#  Derive the BED output path for one input alignment file
function derive_outfile() {
    local infile="${1:-}"
    local base stem

    base="$(basename "${infile}")"
    stem="$(strip_sfx_aln "${base}")"

    printf '%s\n' "${dir_out}/${stem}.bed.gz"
}


#  Convert one BAM or CRAM file to BED
function run_job() {
    local idx="${1:-}"
    local infile outfile samp log_out log_err

    if ! [[ "${idx}" =~ ^[0-9]+$ ]]; then
        echo_err_func "${FUNCNAME[0]}" "invalid task index: '${idx}'."
        return 1
    elif (( idx >= ${#arr_infiles[@]} )); then
        echo_err_func "${FUNCNAME[0]}" "task index is out of bounds: '${idx}'."
        return 1
    fi

    infile="${arr_infiles[idx]}"
    outfile="$(derive_outfile "${infile}")" || return 1
    samp="$(basename "${outfile}")"
    log_out="${dir_eo}/${nam_job}.${samp}.stdout.txt"
    log_err="${dir_eo}/${nam_job}.${samp}.stderr.txt"

    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "idx=${idx}" \
            "infile=${infile}" \
            "outfile=${outfile}" \
            "samp=${samp}"
    fi

    if [[ "${use_awk}" == "true" ]]; then
        if ! \
            convert_alignments_bed_awk \
                "${threads}" \
                "${infile}" \
                "${outfile}" \
                "${ref_fa}" \
                "${log_out}" \
                "${log_err}"
        then
            echo_err \
                "BED file conversion failed for infile '${infile}'." \
                "See log: '${log_err}'."
            return 1
        fi
    else
        if ! \
            convert_alignments_bed_python \
                "${pth_scr_py}" \
                "${infile}" \
                "${outfile}" \
                "${ref_fa}" \
                "${log_out}" \
                "${log_err}"
        then
            echo_err \
                "BED file conversion failed for infile '${infile}'." \
                "See log: '${log_err}'."
            return 1
        fi
    fi

    echo "Successfully converted '${infile}' to '${outfile}'." >&2
}


#  Parse Slurm task state and run one array task
function run_job_slurm() {
    local err_dsc err_ini id_job id_tsk idx out_dsc out_ini
    local infile outfile samp

    id_job="${SLURM_ARRAY_JOB_ID:-}"
    id_tsk="${SLURM_ARRAY_TASK_ID:-}"

    if [[ -z "${id_job}" ]]; then
        echo_err "Slurm array job ID is missing."
        return 1
    elif ! [[ "${id_tsk}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err "Slurm task ID is invalid: '${id_tsk}'."
        return 1
    elif (( id_tsk > ${#arr_infiles[@]} )); then
        echo_err \
            "Slurm task ID '${id_tsk}' exceeds number of input files:" \
            "'${#arr_infiles[@]}'."
        return 1
    fi

    idx=$(( id_tsk - 1 ))
    infile="${arr_infiles[idx]}"
    outfile="$(derive_outfile "${infile}")" || return 1
    samp="$(basename "${outfile}")"

    # shellcheck disable=SC2034
    IFS=',' read -r err_ini out_ini err_dsc out_dsc < <(
        set_logs_slurm \
            "${id_job}" "${id_tsk}" "${samp}" "${dir_eo}" "${nam_job}"
    ) || return 1

    run_job "${idx}" || return 1

    rm -f "${err_ini}" "${out_ini}" || {
        echo_warn \
            "failed to remove initial Slurm log file(s): '${err_ini}' and/or" \
            "'${out_ini}'."
    }
}


#  Dispatch one Slurm-array task or local loop over input files
function run_jobs() {
    local idx

    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        run_job_slurm || return 1
    else
        for idx in "${!arr_infiles[@]}"; do
            run_job "${idx}" || return 1
        done
    fi
}


#  Main script execution
function main() {
    init_defs

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        show_help_main
        echo >&2
        return 0
    fi

    dir_scr="$(resolve_dir_scr "${0##*/}" "$@")" || return 1

    source_helpers_submit "${0##*/}" "${dir_scr}" \
        check_args \
        check_env \
        check_inputs \
        check_numbers \
        format_outputs \
        handle_env \
        manage_slurm \
        process_alignments \
        || return 1

    parse_args "$@"    || return 1
    resolve_paths_scrs || return 1
    validate_args      || return 1
    print_state_debug  || return 1
    setup_env          || return 1
    check_tools        || return 1
    prepare_vecs       || return 1
    validate_vecs      || return 1
    print_vecs_debug   || return 1
    run_jobs           || return 1
}


main "$@"
