#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: submit_filter_alignments.sh
#
# Copyright 2024-2026 by Kris Alavattam
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

# 'sbatch' copies this script, so '--dir_scr' must outrank 'BASH_SOURCE'.
function resolve_dir_scr() {
    local -a args=( "$@" )
    local i

    for (( i = 0; i < ${#args[@]}; i++ )); do
        case "${args[i]}" in
            -ds|--dir[_-]scr)
                if \
                       (( i + 1 < ${#args[@]} )) \
                    && [[ -n "${args[i + 1]}" ]] \
                    && [[ "${args[i + 1]}" != -* ]]
                then
                    printf '%s\n' "${args[i + 1]}"
                    return 0
                fi
                ;;
        esac
    done

    cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
}


dir_scr="$(resolve_dir_scr "$@")" || {
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to resolve the script directory." >&2
    exit 1
}

# shellcheck source=lib/bash/help/help_submit_filter_alignments.sh
source "${dir_scr}/../lib/bash/help/help_submit_filter_alignments.sh" || {
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to source main help." >&2
    exit 1
}


# Define functions.
function parse_filter_alignment_entry() {
    local fil_in="${1:-}"     # Input BAM/CRAM file
    local retain="${2:-}"     # Species selector
    local dir_out="${3:-}"    # Directory for output alignment files
    local out_ext="${4:-bam}" # Final output extension. Output extension
    local samp                # Sample name derived from fil_in
    local nam_fnc             # Function name derived from 'retain'
    local fil_out             # Output file path. Output alignment file
    local show_help           # Help message

    show_help=$(cat << EOM
Usage
-----
  parse_filter_alignment_entry
    [--help] fil_in retain dir_out [out_ext]

  Parse one BAM or CRAM input entry into 'samp', 'nam_fnc', and 'fil_out'.

  This helper derives the sample name from the BAM/CRAM filename, determines which downstream filtering function to use based on '--retain', and constructs the corresponding output path.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fil_in : file
    Input file path. Input BAM or CRAM file.

  2  retain : {'sc', 'sp'}
    Species chromosomes to retain: 'sc' or 'sp'.

  3  dir_out : dir
    Output directory for filtered alignment files.

  4  out_ext : {'bam', 'cram'}
    Final output extension: 'bam' or 'cram' (default: 'bam').

Returns
-------
  Prints a comma-delimited record to stdout:

    samp,nam_fnc,fil_out

  where:
    - 'samp' is derived from the alignment basename without trailing '.bam' or '.cram'
    - 'nam_fnc' is 'filter_alignment_sc' or 'filter_alignment_sp'
    - 'fil_out' is the filtered path in 'dir_out'

Notes
-----
  Runtime requirements:
    - basename
    - bash >= 4.4

  - This helper validates required inputs with 'validate_var'.
  - If 'retain=sc', the output path is '\${dir_out}/\${samp}.sc.\${out_ext}'.
  - If 'retain=sp', the output path is '\${dir_out}/\${samp}.sp.\${out_ext}'.

Examples
--------
  1. Derive an S. cerevisiae BAM output using the default extension.
    '''bash
    parse_filter_alignment_entry work/sample.bam sc work/out
    '''

  2. Derive an S. pombe CRAM output from a CRAM input.
    '''bash
    parse_filter_alignment_entry work/sample.cram sp work/out cram
    '''
EOM
    )

    if [[ "${fil_in}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${fil_in}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'fil_in', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    # Validate input arguments.
    validate_var "fil_in"  "${fil_in}"  || return 1
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

    # Extract sample and function names from input values, and assign 'fil_out'
    # name based on species selector.
    samp="$(basename "${fil_in}")"
    samp="${samp%.bam}"
    samp="${samp%.cram}"

    case "${retain}" in
        sc)
            nam_fnc="filter_alignment_sc"
            fil_out="${dir_out}/${samp}.sc.${out_ext}"
            ;;
        sp)
            nam_fnc="filter_alignment_sp"
            fil_out="${dir_out}/${samp}.sp.${out_ext}"
            ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "'--retain' must be 'sc' or 'sp'."
            return 1
            ;;
    esac

    # Return values.
    echo "${samp},${nam_fnc},${fil_out}"
}


# Execute filtering using the specified function.
function run_filtering() {
    local nam_fnc="${1:-}"   # Name of function to run
    local threads="${2:-}"   # Number of threads
    local fil_in="${3:-}"    # Input BAM/CRAM file
    local fil_out="${4:-}"   # Output file path. Output alignment file
    local mito="${5:-}"      # Retain mito. chr. (true/false)
    local tg="${6:-}"        # Retain SP_II_TG chr. (true/false)
    local mtr="${7:-}"       # Retain SP_MTR chr. (true/false)
    local chk_chr="${8:-}"   # Check chr. in output (true/false)
    local dir_eo="${9:-}"   # Directory for stderr and stdout log files
    local nam_job="${10:-}"  # Job name for log file naming
    local samp="${11:-}"     # Sample name for log file naming
    local ref_fa="${12:-}"   # Reference FASTA file for CRAM input/output
    local log_out log_err    # 'nam_fnc' stdout and stderr log files
    local -a cmd_filter      # Command array for filtering function
    local show_help          # Help message

    show_help=$(cat << EOM
Usage
-----
  run_filtering
    [--help] nam_fnc threads fil_in fil_out mito tg mtr chk_chr dir_eo nam_job samp [ref_fa]

  Execute the specified alignment-filtering function and write stdout/stderr logs to sample-specific files.

  Log files are written to:

    \${dir_eo}/\${nam_job}.\${samp}.stdout.txt
    \${dir_eo}/\${nam_job}.\${samp}.stderr.txt

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  01  nam_fnc : str
    Name of downstream filtering function to run.

  02  threads : int
    Number of threads to use.

  03  fil_in : file
    Input file path. Input BAM or CRAM file.

  04  fil_out : file
    Output file path. Output BAM or CRAM file.

  05  mito : flag
    Retain mitochondrial chromosome. If 'true', pass '--mito'.

  06  tg : flag
    Retain SP_II_TG chromosome. If 'true', pass '--tg'.

  07  mtr : flag
    Retain SP_MTR chromosome. If 'true', pass '--mtr'.

  08  chk_chr : flag
    Check chromosomes in output alignment files. If 'true', pass '--chk_chr'.

  09  dir_eo : dir
    Directory for stderr and stdout log files.

  10  nam_job : str
    Job name used in log-file naming.

  11  samp : str
    Sample name used in log-file naming.

  12  ref_fa : file
    Reference FASTA file for CRAM input/output, or empty string.

Returns
-------
  Returns 0 when filtering finishes successfully; 1 otherwise.

Notes
-----
  Runtime requirements:
    - awk
    - basename
    - bash >= 4.4
    - dirname
    - Input BAM or CRAM index (when filtering S. cerevisiae alignments)
    - mv (when writing S. cerevisiae BAM output)
    - Reference FASTA and required index (when processing CRAM)
    - rm
    - samtools
    - sort (when 'chk_chr' is true)
    - uniq (when 'chk_chr' is true)

  - This helper is a thin wrapper around either 'filter_alignment_sc' or 'filter_alignment_sp'.
  - '--mito', '--tg', '--mtr', and '--chk_chr' are passed only when their values are 'true'.

Examples
--------
  1. Display the filtering dispatch contract without running a workflow.
    '''bash
    run_filtering --help
    '''

  2. Run S. cerevisiae BAM filtering with per-sample log capture.
    '''bash
    run_filtering filter_alignment_sc \\
        2 \\
        work/sample.bam \\
        work/sample.sc.bam \\
        false \\
        false \\
        false \\
        true \\
        work/logs \\
        filter_alignments sample \\
        ''
    '''
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

    # Define paths for log output files.
    log_out="${dir_eo}/${nam_job}.${samp}.stdout.txt"
    log_err="${dir_eo}/${nam_job}.${samp}.stderr.txt"

    cmd_filter=(
        "${nam_fnc}"
            --threads "${threads}"
            --fil_in "${fil_in}"
            --fil_out "${fil_out}"
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

    # Run the filtering function and capture logs.
    if ! "${cmd_filter[@]}" > "${log_out}" 2> "${log_err}"; then
        echo_err_func "${FUNCNAME[0]}" \
            "filtering failed for sample '${samp}'. See log: '${log_err}'."
        return 1
    fi
}


# Source 'source_helpers.sh' and requested helper scripts from 'dir_scr'.
function source_helpers_submit() {
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

    dir_fnc="${dir_scr_arg}/../lib/bash"
    fnc_src="${dir_fnc}/core/source_helpers.sh"

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


# Parse keyword arguments after helper scripts have been sourced.
function parse_args() {
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -en|--env|--env[_-]nam)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_filter_alignments
                    return 1
                }
                env_nam="${2}"
                shift 2
                ;;

            -ds|--dir[_-]scr)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_filter_alignments
                    return 1
                }
                dir_scr="${2}"
                shift 2
                ;;

            -t|--thr|--threads)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_filter_alignments
                    return 1
                }
                threads="${2}"
                shift 2
                ;;

            -ci|--csv[_-]fil_in)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_filter_alignments
                    return 1
                }
                csv_fil_in="${2}"
                shift 2
                ;;

            -do|--dir[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_filter_alignments
                    return 1
                }
                dir_out="${2}"
                shift 2
                ;;

            -rt|--retain)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_filter_alignments
                    return 1
                }
                retain="${2}"
                shift 2
                ;;

            -rf|--ref[_-]fa)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_filter_alignments
                    return 1
                }
                ref_fa="${2}"
                shift 2
                ;;

            -ox|--out[_-]ext)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_filter_alignments
                    return 1
                }
                out_ext="${2}"
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

            -deo|--dir[_-]eo)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_filter_alignments
                    return 1
                }
                dir_eo="${2}"
                shift 2
                ;;

            -nj|--nam[_-]job)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_filter_alignments
                    return 1
                }
                nam_job="${2}"
                shift 2
                ;;

            *)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                help_submit_filter_alignments
                return 1
                ;;
        esac
    done
}


# Canonicalize scalar argument aliases.
function canonicalize_args() {
    retain="${retain,,}"
    out_ext="${out_ext,,}"
}


# Validate required arguments and paths.
function validate_args() {
    validate_var     "env_nam"    "${env_nam}"         || return 1
    validate_var_dir "dir_scr"    "${dir_scr}" 0 false || return 1
    validate_var     "threads"    "${threads}"         || return 1
    validate_var     "csv_fil_in" "${csv_fil_in}"      || return 1
    validate_var_dir "dir_out"    "${dir_out}"         || return 1
    validate_var_dir "dir_eo"    "${dir_eo}"         || return 1
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


# Print debug argument variable assignments.
function print_state_debug() {
    if [[ "${debug}" == "true" ]]; then
        echo
        debug_var \
            "env_nam=${env_nam}" \
            "dir_scr=${dir_scr}" \
            "retain=${retain}" \
            "threads=${threads}" \
            "csv_fil_in=${csv_fil_in}" \
            "dir_out=${dir_out}" \
            "out_ext=${out_ext}" \
            "mito=${mito}" \
            "tg=${tg}" \
            "mtr=${mtr}" \
            "chk_chr=${chk_chr}" \
            "ref_fa=${ref_fa:-UNSET}" \
            "dir_eo=${dir_eo}" \
            "nam_job=${nam_job}"
    fi
}


# Normalize species-specific optional flags.
function normalize_flags() {
    if [[ "${retain}" == "sc" ]]; then
        if [[ "${tg}" == "true" && "${mtr}" == "true" ]]; then
            echo_warn \
                "'--tg' and '--mtr' were supplied with '--retain sc' and" \
                "will be ignored."
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


# Initialize hardcoded argument variables.
function init_args_hardcoded() {
    # If true, run script in debug mode.
    debug=true
}


# Initialize argument variables, assigning default values where applicable.
function init_arg_defs() {
    env_nam="env_protocol"
    retain="sc"
    threads=4
    csv_fil_in=""
    dir_out=""
    out_ext="bam"
    mito=false
    tg=false
    mtr=false
    chk_chr=false
    ref_fa=""
    dir_eo=""
    nam_job="filter_alignments"
}


# Initialize hardcoded arguments and user-facing argument defaults.
function init_defs() {
    init_args_hardcoded
    init_arg_defs
}


# Reconstruct array from serialized inputs.
function prepare_vecs() {
    unset arr_fil_in && declare -ga arr_fil_in
    IFS=',' read -r -a arr_fil_in <<< "${csv_fil_in}"
}


# Validate reconstructed input arrays.
function validate_vecs() {
    local fil_in

    check_arr_nonempty "arr_fil_in" "csv_fil_in" || return 1

    for fil_in in "${arr_fil_in[@]}"; do
        if [[ "${fil_in,,}" == *.cram && -z "${ref_fa}" ]]; then
            echo_err \
                "'--ref_fa' is required when '--csv_fil_in' contains CRAM" \
                "input: '${fil_in}'."
            return 1
        fi
    done
    unset fil_in

    if [[ "${out_ext}" == "cram" && -z "${ref_fa}" ]]; then
        echo_err "'--ref_fa' is required when '--out_ext cram'."
        return 1
    fi
}


# Activate environment.
function setup_env() {
    handle_env "${env_nam}" || return 1
}


# Print debug output for reconstructed input arrays.
function print_vecs_debug() {
    if [[ "${debug}" == "true" ]]; then
        echo "\${#arr_fil_in[@]}=${#arr_fil_in[@]}" && echo
        echo "arr_fil_in=( ${arr_fil_in[*]} )"      && echo
    fi
}


# Dispatch Slurm-array or local worker jobs.
function run_jobs() {
    local err_dsc err_ini id_job id_tsk idx fil_in nam_fnc out_dsc out_ini
    local fil_out samp

    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        # Workflow mode: Slurm.
        id_job="${SLURM_ARRAY_JOB_ID}"
        id_tsk="${SLURM_ARRAY_TASK_ID}"

        if ! [[ "${id_tsk}" =~ ^[1-9][0-9]*$ ]]; then
            echo_err "Slurm task ID is invalid: '${id_tsk}'."
            return 1
        elif (( id_tsk > ${#arr_fil_in[@]} )); then
            echo_err \
                "Slurm task ID '${id_tsk}' exceeds number of BAM entries:" \
                "'${#arr_fil_in[@]}'."
            return 1
        else
            idx=$(( id_tsk - 1 ))
        fi

        fil_in="${arr_fil_in[idx]}"

        if [[ "${debug}" == "true" ]]; then debug_var "fil_in=${fil_in}"; fi

        IFS=',' read -r samp nam_fnc fil_out < <(
            parse_filter_alignment_entry \
                "${fil_in}" "${retain}" "${dir_out}" "${out_ext}"
        ) || return 1

        if [[ "${debug}" == "true" ]]; then
            debug_var \
                "samp=${samp}" \
                "nam_fnc=${nam_fnc}" \
                "fil_out=${fil_out}"
        fi

        IFS=',' read -r err_ini out_ini err_dsc out_dsc < <(
            set_logs_slurm \
                "${id_job}" "${id_tsk}" "${samp}" "${dir_eo}" "${nam_job}"
        ) || return 1

        if [[ "${debug}" == "true" ]]; then
            debug_var \
                "err_ini=${err_ini}" \
                "out_ini=${out_ini}" \
                "err_dsc=${err_dsc}" \
                "out_dsc=${out_dsc}"
        fi

        if ! \
            run_filtering \
                "${nam_fnc}" "${threads}" "${fil_in}" "${fil_out}" "${mito}" \
                "${tg}" "${mtr}" "${chk_chr}" "${dir_eo}" "${nam_job}" \
                "${samp}" "${ref_fa}"
        then
            echo_err "failed to filter alignment file: '${fil_in}'."
            return 1
        fi

        rm "${err_ini}" "${out_ini}"
    else
        # Workflow mode: GNU Parallel/serial.
        for idx in "${!arr_fil_in[@]}"; do
            fil_in="${arr_fil_in[idx]}"

            if [[ "${debug}" == "true" ]]; then debug_var "fil_in=${fil_in}"; fi

            IFS=',' read -r samp nam_fnc fil_out < <(
                parse_filter_alignment_entry \
                    "${fil_in}" "${retain}" "${dir_out}" "${out_ext}"
            ) || return 1

            if [[ "${debug}" == "true" ]]; then
                debug_var \
                    "samp=${samp}" \
                    "nam_fnc=${nam_fnc}" \
                    "fil_out=${fil_out}"
            fi

            if ! \
                run_filtering \
                    "${nam_fnc}" "${threads}" "${fil_in}" "${fil_out}" \
                    "${mito}" "${tg}" "${mtr}" "${chk_chr}" "${dir_eo}" \
                    "${nam_job}" "${samp}" "${ref_fa}"
            then
                echo_err "failed to filter alignment file: '${fil_in}'."
                return 1
            fi
        done
    fi
}


# Main script execution.
function main() {
    init_defs

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        help_submit_filter_alignments
        echo >&2
        return 0
    fi

    # First-pass parse: resolve 'dir_scr' before using sourced parser helpers.
    source_helpers_submit "${0##*/}" "${dir_scr}" \
        check_args \
        check_inputs \
        filter_alignment \
        format_outputs \
        handle_env \
        manage_slurm \
        help/help_submit_filter_alignments \
        || return 1

    parse_args "$@"   || return 1
    canonicalize_args || return 1
    validate_args     || return 1
    normalize_flags   || return 1
    print_state_debug || return 1
    setup_env         || return 1
    prepare_vecs      || return 1
    validate_vecs     || return 1
    print_vecs_debug  || return 1
    run_jobs          || return 1
}


main "$@"
