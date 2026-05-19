#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: submit_filter_alignments.sh
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

#  If true, run script in debug mode
debug=true


#  Define functions
#  Parse a BAM/CRAM input entry into 'samp', 'nam_fnc', and 'outfile'
function parse_filter_alignment_entry() {
    local infile="${1:-}"   # Input BAM/CRAM file
    local retain="${2:-}"   # Species selector
    local dir_out="${3:-}"  # Directory for output alignment files
    local out_ext="${4:-bam}" # Output extension
    local samp      # Sample name derived from infile
    local nam_fnc   # Function name derived from 'retain'
    local outfile   # Output alignment file
    local show_help # Help message

    show_help=$(cat << EOM
Usage:
  parse_filter_alignment_entry
    [-h|--hlp|--help] infile retain dir_out [out_ext]

Description:
  Parse one BAM or CRAM input entry into 'samp', 'nam_fnc', and 'outfile'.

  This helper derives the sample name from the BAM/CRAM filename, determines which downstream filtering function to use based on '--retain', and constructs the corresponding output path.

Positional arguments:
  1  infile   <str>  Input BAM or CRAM file.
  2  retain   <str>  Species selector; must be 'sc' or 'sp'.
  3  dir_out  <str>  Output directory for filtered alignment files.
  4  out_ext  <str>  Output extension: 'bam' or 'cram' (default: bam).

Returns:
  Prints a comma-delimited record to stdout:

    samp,nam_fnc,outfile

  where:
    - 'samp' is derived from the alignment basename without trailing '.bam' or '.cram'
    - 'nam_fnc' is 'filter_alignment_sc' or 'filter_alignment_sp'
    - 'outfile' is the filtered path in 'dir_out'

Notes:
  - This helper validates required inputs with 'validate_var'.
  - If 'retain=sc', the output path is '\${dir_out}/\${samp}.sc.\${out_ext}'.
  - If 'retain=sp', the output path is '\${dir_out}/\${samp}.sp.\${out_ext}'.
EOM
    )

    if [[ "${infile}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${infile}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'infile', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Validate input arguments
    validate_var "infile"  "${infile}"  || return 1
    validate_var "retain"  "${retain}"  || return 1
    validate_var "dir_out" "${dir_out}" || return 1
    validate_var "out_ext" "${out_ext}" || return 1

    case "${out_ext}" in
        bam|cram) : ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "'out_ext' must be 'bam' or 'cram': '${out_ext}'."
            return 1
            ;;
    esac

    #  Extract sample and function names from input values, and assign outfile
    #+ name based on species selector
    samp="$(basename "${infile}")"
    samp="${samp%.bam}"
    samp="${samp%.cram}"

    case "${retain}" in
        sc)
            nam_fnc="filter_alignment_sc"
            outfile="${dir_out}/${samp}.sc.${out_ext}"
            ;;
        sp)
            nam_fnc="filter_alignment_sp"
            outfile="${dir_out}/${samp}.sp.${out_ext}"
            ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "'--retain' must be 'sc' or 'sp'."
            return 1
            ;;
    esac

    #  Return values
    echo "${samp},${nam_fnc},${outfile}"
}


#  Execute filtering using the specified function
function run_filtering() {
    local nam_fnc="${1:-}"   # Name of function to run
    local threads="${2:-}"   # Number of threads
    local infile="${3:-}"    # Input BAM/CRAM file
    local outfile="${4:-}"   # Output alignment file
    local mito="${5:-}"      # Retain mito. chr. (true/false)
    local tg="${6:-}"        # Retain SP_II_TG chr. (true/false)
    local mtr="${7:-}"       # Retain SP_MTR chr. (true/false)
    local chk_chr="${8:-}"   # Check chr. in output (true/false)
    local err_out="${9:-}"   # Directory for stderr and stdout logs
    local nam_job="${10:-}"  # Job name for log file naming
    local samp="${11:-}"     # Sample name for log file naming
    local ref_fa="${12:-}"   # Reference FASTA for CRAM input/output
    local log_out log_err    # 'nam_fnc' stdout and stderr log files
    local -a cmd_filter      # Command array for filtering function
    local show_help          # Help message

    show_help=$(cat << EOM
Usage:
  run_filtering
    [-h|--hlp|--help] nam_fnc threads infile outfile mito tg mtr chk_chr err_out nam_job samp [ref_fa]

Description:
  Execute the specified alignment-filtering function and write stdout/stderr
  logs to sample-specific files.

  Log files are written to:

    \${err_out}/\${nam_job}.\${samp}.stdout.txt
    \${err_out}/\${nam_job}.\${samp}.stderr.txt

Positional arguments:
   1  nam_fnc   <str>  Name of downstream filtering function to run.
   2  threads   <int>  Number of threads.
   3  infile    <str>  Input BAM or CRAM file.
   4  outfile   <str>  Output BAM or CRAM file.
   5  mito      <flag>  If 'true', pass '--mito'.
   6  tg        <flag>  If 'true', pass '--tg'.
   7  mtr       <flag>  If 'true', pass '--mtr'.
   8  chk_chr   <flag>  If 'true', pass '--chk_chr'.
   9  err_out   <str>  Directory for stderr/stdout log files.
  10  nam_job   <str>  Job name used in log-file naming.
  11  samp      <str>  Sample name used in log-file naming.
  12  ref_fa    <str>  Reference FASTA for CRAM input/output, or empty string.

Notes:
  - This helper is a thin wrapper around either 'filter_alignment_sc' or 'filter_alignment_sp'.
  - '--mito', '--tg', '--mtr', and '--chk_chr' are passed only when their values are 'true'.
EOM
    )

    if [[ "${nam_fnc}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${nam_fnc}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'nam_fnc', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Define paths for log output files
    log_out="${err_out}/${nam_job}.${samp}.stdout.txt"
    log_err="${err_out}/${nam_job}.${samp}.stderr.txt"

    cmd_filter=(
        "${nam_fnc}"
            --threads "${threads}"
            --infile "${infile}"
            --outfile "${outfile}"
    )

    if [[ -n "${ref_fa}" ]]; then
        cmd_filter+=( --ref_fa "${ref_fa}" )
    fi

    if [[ "${mito}" == "true" ]]; then
        cmd_filter+=( --mito )
    fi

    if [[ "${tg}" == "true" ]]; then
        cmd_filter+=( --tg )
    fi

    if [[ "${mtr}" == "true" ]]; then
        cmd_filter+=( --mtr )
    fi

    if [[ "${chk_chr}" == "true" ]]; then
        cmd_filter+=( --chk_chr )
    fi

    #  Run the filtering function and capture logs
    if ! "${cmd_filter[@]}" > "${log_out}" 2> "${log_err}"; then
        echo_err_func "${FUNCNAME[0]}" \
            "filtering failed for sample '${samp}'. See log: '${log_err}'."
        return 1
    fi
}


#  Resolve '--dir_scr' before sourced parser helpers are available
function resolve_dir_scr() {
    local script="${1:-}"
    shift

    local -a args=( "$@" )
    local i=0

    if [[ -z "${script}" ]]; then
        script="unknown_script"
    fi

    for (( i = 0; i < ${#args[@]}; i++ )); do
        case "${args[i]}" in
            -ds|--dir[_-]scr)
                if (( i + 1 >= ${#args[@]} )) \
                    || [[ -z "${args[i + 1]:-}" || "${args[i + 1]}" == -* ]]
                then
                    echo "error(${script}):" \
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

    echo "error(${script}):" \
        "required option '--dir_scr' was not supplied." >&2
    echo >&2
    show_help_main
    return 1
}


#  Source 'source_helpers.sh' and requested helper scripts from 'dir_scr'
function source_submit_helpers() {
    local script="${1:-}"
    local dir_scr_arg="${2:-}"
    local fnc_src

    if (( $# < 2 )); then
        echo "error(${script:-unknown_script}):" \
            "expected at least two arguments: 'script' and 'dir_scr_arg'." >&2
        return 1
    fi

    shift 2

    if [[ -z "${script}" ]]; then
        script="unknown_script"
    fi

    if [[ -z "${dir_scr_arg}" ]]; then
        echo "error(${script}):" \
            "positional argument 2, 'dir_scr_arg', is missing." >&2
        return 1
    elif [[ ! -d "${dir_scr_arg}" ]]; then
        echo "error(${script}):" \
            "script directory not found: '${dir_scr_arg}'." >&2
        return 1
    elif (( $# < 1 )); then
        echo "error(${script}):" \
            "at least one helper script name must be supplied." >&2
        return 1
    fi

    dir_fnc="${dir_scr_arg}/functions"
    fnc_src="${dir_fnc}/source_helpers.sh"

    if [[ ! -f "${fnc_src}" ]]; then
        echo "error(${script}):" \
            "script not found: '${fnc_src}'." >&2
        return 1
    fi

    # shellcheck disable=SC1090
    source "${fnc_src}" || {
        echo "error(${script}):" \
            "failed to source '${fnc_src}'." >&2
        return 1
    }

    source_helpers "${dir_fnc}" "$@" || {
        echo "error(${script}):" \
            "failed to source required helper scripts." >&2
        return 1
    }
}


#  Parse keyword arguments after helper scripts have been sourced
function parse_args() {
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
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

            -i|-fi|-ci|--infile|--infiles|--fil[_-]in|--csv[_-]infile|--csv[_-]infiles)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                csv_infile="${2}"
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

            -rt|--retain)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                #TODO: take advantage of Bash 4.4, i.e., "${2,,}"
                retain="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
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

            -ox|--out[_-]ext)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                out_ext="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
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
                    show_help_main
                    return 1
                }
                err_out="${2}"
                shift 2
                ;;

            -nj|--nam[_-]job)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                nam_job="${2}"
                shift 2
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


#  Validate required arguments and paths
function validate_args() {
    validate_var     "env_nam"    "${env_nam}"         || return 1
    validate_var_dir "dir_scr"    "${dir_scr}" 0 false || return 1
    validate_var     "threads"    "${threads}"         || return 1
    validate_var     "csv_infile" "${csv_infile}"      || return 1
    validate_var_dir "dir_out"    "${dir_out}"         || return 1
    validate_var_dir "err_out"    "${err_out}"         || return 1
    validate_var     "nam_job"    "${nam_job}"         || return 1
    validate_var     "out_ext"    "${out_ext}"         || return 1

    if [[ -n "${ref_fa}" ]]; then
        validate_var_file "ref_fa" "${ref_fa}" || return 1
    fi

    case "${out_ext}" in
        bam|cram) : ;;
        *)
            echo_err "'--out_ext' must be 'bam' or 'cram': '${out_ext}'."
            return 1
            ;;
    esac

    case "${retain}" in
        sc|sp) : ;;
        *)
            echo_err "'--retain' must be 'sc' or 'sp': '${retain}'."
            return 1
            ;;
    esac
}


#  Print debug argument variable assignments
function print_debug_args() {
    if [[ "${debug}" == "true" ]]; then
        echo
        debug_var \
            "env_nam=${env_nam}" \
            "dir_scr=${dir_scr}" \
            "retain=${retain}" \
            "threads=${threads}" \
            "csv_infile=${csv_infile}" \
            "dir_out=${dir_out}" \
            "out_ext=${out_ext}" \
            "mito=${mito}" \
            "tg=${tg}" \
            "mtr=${mtr}" \
            "chk_chr=${chk_chr}" \
            "ref_fa=${ref_fa:-UNSET}" \
            "err_out=${err_out}" \
            "nam_job=${nam_job}"
    fi
}


#  Normalize species-specific optional flags
function normalize_flags() {
    if [[ "${retain}" == "sc" ]]; then
        if [[ "${tg}" == "true" && "${mtr}" == "true" ]]; then
            echo_warn \
                "'--tg' and '--mtr' were supplied with '--retain sc' and will be" \
                "ignored."
            tg=false
            mtr=false
        elif [[ "${tg}" == "true" ]]; then
            echo_warn \
                "'--tg' was supplied with '--retain sc' and will be ignored."
            tg=false
        elif [[ "${mtr}" == "true" ]]; then
            echo_warn \
                "'--mtr' was supplied with '--retain sc' and will be ignored."
            mtr=false
        fi
    fi
}


#  Initialize argument variables, assigning default values where applicable
env_nam="env_protocol"
dir_scr=""
retain="sc"
threads=4
csv_infile=""
dir_out=""
out_ext="bam"
mito=false
tg=false
mtr=false
chk_chr=false
ref_fa=""
err_out=""
nam_job="filter_alignments"


function show_help_main() {
    cat << EOM >&2
Usage:
  submit_filter_alignments.sh
    [-h|--hlp|--help] [-en|--env_nam <str>] -ds|--dir_scr <str> [-t|--threads <int>]
    -i|--csv_infile <str> -do|--dir_out <str> [-ox|--out_ext <str>]
    [-rt|--retain <str>] [-r|--ref_fa <file>] [-m|--mito] [-tg|--tg] [-mr|--mtr]
    [-cc|--chk_chr] -eo|--err_out <str> [-nj|--nam_job <str>]

Description:
  Submit or execute one or more alignment-filtering jobs by calling downstream functions 'filter_alignment_sc' or 'filter_alignment_sp'.

  This wrapper
    - parses a comma-delimited list of BAM or CRAM input files,
    - determines which downstream filtering function to run based on '--retain',
    - activates the requested Conda/Mamba environment, and then
    - runs filtering either under Slurm array execution or by serial/GNU-Parallel-style iteration, depending on how the script is invoked.

  For each input BAM or CRAM file, this script writes BAM or CRAM output and log files to:

    \${err_out}/\${nam_job}.\${samp}.stdout.txt
    \${err_out}/\${nam_job}.\${samp}.stderr.txt

Keyword arguments:
  -en, --env_nam  <str>
    Conda/Mamba environment to activate.

  -ds, --dir_scr  <str>
    Directory containing scripts and functions.

  -t, --threads  <int>
    Number of threads to use.

  -i, --csv_infile  <str>
    Comma-delimited list of input BAM or CRAM files.

    Compatibility aliases include '--infile', '--infiles', '--fil_in', and '--csv_infiles'.

  -do, --dir_out  <str>
    Directory in which filtered alignment files will be written.

  -ox, --out_ext  <str>
    Filtered output extension: 'bam' or 'cram' (default: '${out_ext}').

  -rt, --retain  <str>
    Species chromosomes to retain: 'sc' or 'sp'.

  -r, --ref_fa  <file>
    Reference FASTA required when any input file is CRAM or '--out_ext cram'.

  -m, --mito  <flag>
    If supplied, retain the mitochondrial chromosome.

  -tg, --tg  <flag>
    If supplied, retain chromosome 'SP_II_TG'.

  -mr, --mtr  <flag>
    If supplied, retain chromosome 'SP_MTR'.

  -cc, --chk_chr  <flag>
    If supplied, check chromosomes in output alignment files.

  -eo, --err_out  <str>
    Directory in which stderr/stdout log files will be written.

  -nj, --nam_job  <str>
    Job name used in log-file naming.

Notes:
  - All arguments are required with the following notes and exceptions:
    + '--env_nam' defaults to 'env_nam=${env_nam}' if not specified.
    + '--retain' defaults to 'retain=${retain}' if not specified.
    + '--threads' defaults to 'threads=${threads}' if not specified.
    + '--out_ext' defaults to 'out_ext=${out_ext}' if not specified.
    + '--ref_fa' is required when '--csv_infile' contains CRAM input or '--out_ext cram'.
    + '--mito', '--tg', '--mtr', and '--chk_chr' are optional flags.
    + '--tg' and '--mtr' are only meaningful when '--retain sp' is used; if supplied with '--retain sc', they are ignored with a warning.
    + '--nam_job' defaults to 'nam_job=${nam_job}' if not specified.
EOM
}


#  Main script execution
function main() {
    #  Display help message if a help option or no arguments are given
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        show_help_main
        echo >&2
        exit 0
    fi

    #  First-pass parse: resolve 'dir_scr' before using sourced parser helpers
    dir_scr="$(resolve_dir_scr "${0##*/}" "$@")" || exit 1

    source_submit_helpers "${0##*/}" "${dir_scr}" \
        check_args \
        check_inputs \
        filter_alignment \
        format_outputs \
        handle_env \
        manage_slurm \
        || exit 1

    parse_args "$@"  || exit 1
    validate_args    || exit 1
    normalize_flags  || exit 1
    print_debug_args || exit 1

    #  Activate environment
    handle_env "${env_nam}" || exit 1

    #  Reconstruct array from serialized string
    IFS=',' read -r -a arr_infile <<< "${csv_infile}"
    check_arr_nonempty "arr_infile" "csv_infile" || exit 1

    for infile in "${arr_infile[@]}"; do
        if [[ "${infile,,}" == *.cram && -z "${ref_fa}" ]]; then
            echo_err \
                "'--ref_fa' is required when '--csv_infile' contains CRAM" \
                "input: '${infile}'."
            exit 1
        fi
    done
    unset infile

    if [[ "${out_ext}" == "cram" && -z "${ref_fa}" ]]; then
        echo_err "'--ref_fa' is required when '--out_ext cram'."
        exit 1
    fi

    #  Debug output to check number of array elements and array element values
    if [[ "${debug}" == "true" ]]; then
        echo "\${#arr_infile[@]}=${#arr_infile[@]}" && echo
        echo "arr_infile=( ${arr_infile[*]} )"      && echo
    fi

    #  Determine and run mode: Slurm or GNU Parallel/serial
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        #  Mode: Slurm
        id_job="${SLURM_ARRAY_JOB_ID}"
        id_tsk="${SLURM_ARRAY_TASK_ID}"

        if ! [[ "${id_tsk}" =~ ^[1-9][0-9]*$ ]]; then
            echo_err "Slurm task ID is invalid: '${id_tsk}'."
            exit 1
        elif (( id_tsk > ${#arr_infile[@]} )); then
            echo_err \
                "Slurm task ID '${id_tsk}' exceeds number of BAM entries:" \
                "'${#arr_infile[@]}'."
            exit 1
        else
            idx=$(( id_tsk - 1 ))
        fi

        infile="${arr_infile[idx]}"

        if [[ "${debug}" == "true" ]]; then debug_var "infile=${infile}"; fi

        IFS=',' read -r samp nam_fnc outfile < <(
            parse_filter_alignment_entry \
                "${infile}" "${retain}" "${dir_out}" "${out_ext}"
        ) || exit 1

        if [[ "${debug}" == "true" ]]; then
            debug_var \
                "samp=${samp}" \
                "nam_fnc=${nam_fnc}" \
                "outfile=${outfile}"
        fi

        IFS=',' read -r err_ini out_ini err_dsc out_dsc < <(
            set_logs_slurm \
                "${id_job}" "${id_tsk}" "${samp}" "${err_out}" "${nam_job}"
        ) || exit 1

        if [[ "${debug}" == "true" ]]; then
            debug_var \
                "err_ini=${err_ini}" \
                "out_ini=${out_ini}" \
                "err_dsc=${err_dsc}" \
                "out_dsc=${out_dsc}"
        fi

        if ! \
            run_filtering \
                "${nam_fnc}" "${threads}" "${infile}" "${outfile}" "${mito}" \
                "${tg}" "${mtr}" "${chk_chr}" "${err_out}" "${nam_job}" \
                "${samp}" "${ref_fa}"
        then
            echo_err "failed to filter alignment file: '${infile}'."
            exit 1
        fi

        rm "${err_ini}" "${out_ini}"
    else
        #  Mode: GNU Parallel/serial
        for idx in "${!arr_infile[@]}"; do
            infile="${arr_infile[idx]}"

            if [[ "${debug}" == "true" ]]; then debug_var "infile=${infile}"; fi

            IFS=',' read -r samp nam_fnc outfile < <(
                parse_filter_alignment_entry \
                    "${infile}" "${retain}" "${dir_out}" "${out_ext}"
            ) || exit 1

            if [[ "${debug}" == "true" ]]; then
                debug_var \
                    "samp=${samp}" \
                    "nam_fnc=${nam_fnc}" \
                    "outfile=${outfile}"
            fi

            if ! \
                run_filtering \
                    "${nam_fnc}" "${threads}" "${infile}" "${outfile}" \
                    "${mito}" "${tg}" "${mtr}" "${chk_chr}" "${err_out}" \
                    "${nam_job}" "${samp}" "${ref_fa}"
            then
                echo_err "failed to filter alignment file: '${infile}'."
                exit 1
            fi
        done
    fi
}


main "$@"
