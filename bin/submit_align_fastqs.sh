#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: submit_align_fastqs.sh
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

# 'sbatch' runs a copy, so 'BASH_SOURCE' can resolve outside the repo; if so,
# fail here rather than at a subsequent 'source'.
if [[ ! -d "${dir_scr}/../lib/bash" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "cannot locate helper libraries from '${dir_scr}'; pass '--dir_scr'" \
        "when this script is run from a copy." >&2
    exit 1
fi

# shellcheck source=lib/bash/help/help_submit_align_fastqs.sh
source "${dir_scr}/../lib/bash/help/help_submit_align_fastqs.sh" || {
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to source main help." >&2
    exit 1
}


# Parse one FASTQ entry into 'fq_1', 'fq_2', 'samp', and alignment 'fil_out'.
function parse_entry_align_fastq() {
    local fil_in="${1:-}"      # Input FASTQ file(s).
    local sfx_se="${2:-}"      # Suffix for SE FASTQ files.
    local sfx_pe="${3:-}"      # Suffix for PE FASTQ files (FASTQ #1).
    local dir_out="${4:-}"     # Directory for output files.
    local out_ext="${5:-bam}"  # Alignment output extension.
    local fq_1                 # FASTQ file #1.
    local fq_2                 # FASTQ file #2, or 'NA' for SE.
    local samp                 # Sample name.
    local fil_out              # Output file path. Output alignment file.
    local show_help            # Help message.

    # TODO: break parameter options under Usage across semantic paragraphs.
    show_help=$(cat << EOM
Usage
-----
  parse_entry_align_fastq
    [--help] fil_in sfx_se sfx_pe dir_out [out_ext]

  Parse one Input file path. FASTQ input entry into 'fq_1', 'fq_2', 'samp', and alignment 'fil_out'.

  The input entry may represent either:
    - single-end data: one FASTQ path
    - paired-end data: two comma-delimited FASTQ paths

    This helper calls 'process_sequences.sh::parse_fastq_entry' to parse the FASTQ input and then appends the alignment output path:

    \${dir_out}/\${samp}.\${out_ext}

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fil_in : file
    Input file path. FASTQ input entry. For single-end data, this is one FASTQ file. For paired-end data, this is a comma-delimited FASTQ pair.

  2  sfx_se : str
    Suffix to strip from single-end FASTQ filenames when deriving the sample name.

  3  sfx_pe : str
    Suffix to strip from paired-end FASTQ read-1 filenames when deriving the sample name.

  4  dir_out : dir
    Output directory used to construct the alignment fil_out path.

  5  out_ext : {'bam', 'cram'}
    Final output extension for alignment files: 'bam' or 'cram' (default: 'bam').

Returns
-------
  - Prints 'fq_1;fq_2;samp;fil_out' to stdout, where 'fq_2' is 'NA' for single-end data.
  - Returns 0 when the FASTQ entry is parsed successfully; 1 otherwise.

Notes
-----
  Runtime requirements:
    - awk
    - basename
    - bash >= 4.4

  - This helper validates required inputs with 'validate_var'.
  - Parsing of the FASTQ entry itself is delegated to 'parse_fastq_entry'.
  - The output alignment file path is always constructed as '\${dir_out}/\${samp}.\${out_ext}'.

Examples
--------
  1. Parse one single-end fixture and derive its BAM output path.
    '''bash
    parse_entry_align_fastq \\
        tests/fixtures/align_fastqs/fastq/se/tiny_se.atria.fastq.gz \\
        .atria.fastq.gz \\
        _R1.atria.fastq.gz \\
        work \\
        bam
    '''

  2. Parse paired-end fixtures and derive their CRAM output path.
    '''bash
    parse_entry_align_fastq \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R1.atria.fastq.gz,tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        .atria.fastq.gz \\
        _R1.atria.fastq.gz \\
        work \\
        cram
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

    validate_var "fil_in"  "${fil_in}"  || return 1
    validate_var "sfx_se"  "${sfx_se}"  || return 1
    validate_var "sfx_pe"  "${sfx_pe}"  || return 1
    validate_var "dir_out" "${dir_out}" || return 1

    IFS=';' read -r fq_1 fq_2 samp < <(
        parse_fastq_entry "${fil_in}" "${sfx_se}" "${sfx_pe}"
    ) || return 1

    #MAYBE: these 'fq_(1|2)' checks are a little redundant
    validate_var "fq_1" "${fq_1}" || return 1

    if [[ "${fq_2}" != "NA" ]]; then
        validate_var "fq_2" "${fq_2}" || return 1
    fi

    case "${out_ext}" in
        bam|cram) : ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "positional argument 5, 'out_ext', must be 'bam' or 'cram':" \
                " '${out_ext}'."
            return 1
            ;;
    esac

    fil_out="${dir_out}/${samp}.${out_ext}"

    echo "${fq_1};${fq_2};${samp};${fil_out}"
}
# MAYBE: 'parse_entry_align_fastq' still uses 'validate_var' rather than
# file-aware helpers such as 'validate_file' and 'validate_var_file'; revisit
# if stronger file validation is wanted here.


# Execute alignment using function 'align_fastqs.sh::align_fastqs'.
function run_alignment() {
    local threads="${1:-}"
    local aligner="${2:-}"
    local bt2_mode="${3:-}"
    local bwa_alg="${4:-}"
    local mapq="${5:-}"
    local req_flg="${6:-}"
    local index="${7:-}"
    local ref_fa="${8:-}"
    local fq_1="${9:-}"
    local fq_2="${10:-}"
    local fil_out="${11:-}"
    local qname="${12:-}"
    local dir_eo="${13:-}"
    local nam_job="${14:-}"
    local samp="${15:-}"
    local log_out log_err
    local -a cmd_aln
    local show_help

    # TODO: break parameter options under Usage across semantic paragraphs.
    show_help=$(cat << EOM
Usage
-----
  run_alignment
    [--help] threads aligner bt2_mode bwa_alg mapq req_flg index ref_fa fq_1 fq_2 fil_out qname dir_eo nam_job samp

  Construct log-file paths and then run 'align_fastqs.sh::align_fastqs', writing stdout to

    \${dir_eo}/\${nam_job}.\${samp}.stdout.txt

  and stderr to

    \${dir_eo}/\${nam_job}.\${samp}.stderr.txt

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  01  threads : int
    Number of threads to use.

  02  aligner : {'bowtie2', 'bwa', 'bwa-mem2'}
    Alignment program to use: 'bowtie2', 'bwa', or 'bwa-mem2'.

  03  bt2_mode : {'local', 'global', 'end-to-end'}
    Bowtie 2 alignment type when '--aligner bowtie2': 'local', 'global', or 'end-to-end'.

  04  bwa_alg : {'mem', 'aln'}
    BWA algorithm when '--aligner bwa': 'mem' or 'aln'.

  05  mapq : int
    MAPQ threshold for filtering alignment output files. If 'mapq > 0', '--mapq' is passed to 'align_fastqs'; otherwise, it is omitted.

  06  req_flg : flag
    Require SAM flag bit 2 for properly paired alignments. If 'true', pass '--req_flg'.

  07  index : path
    Path to the aligner index/reference.

  08  ref_fa : file
    Reference FASTA file for CRAM output; ignored for BAM output.

  09  fq_1 : file
    First FASTQ input file.

  10  fq_2 : file
    Second FASTQ input file, or 'NA' for single-end data.

  11  fil_out : file
    Output file path. Output alignment file; must end in '.bam' or '.cram'.

  12  qname : flag
    Retain queryname-sorted intermediate alignment files. If 'true', pass '--qname'.

  13  dir_eo : dir
    Directory for stderr and stdout log files.

  14  nam_job : str
    Job name used in log-file naming.

  15  samp : str
    Sample name used in log-file naming.

Returns
-------
  - Writes alignment stdout and stderr to per-sample log files.
  - Returns 0 when alignment completes successfully; 1 otherwise.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - bowtie2 (when 'aligner' is 'bowtie2')
    - bwa (when 'aligner' is 'bwa')
    - bwa-mem2 (when 'aligner' is 'bwa-mem2')
    - grep
    - Reference FASTA and required index (when writing CRAM)
    - samtools

  - This helper is a thin wrapper around 'align_fastqs'.
  - '--fq_2' is passed only when 'fq_2 != NA'.
  - '--mapq' is passed only when 'mapq > 0'.
  - '--req_flg' and '--qname' are passed only when their values are 'true'.
  - '--ref_fa' is passed only when 'fil_out' ends in '.cram'.

Examples
--------
  1. Display the execution helper's full positional contract.
    '''bash
    run_alignment --help
    '''

  2. Run one single-end Bowtie 2 fixture and capture per-sample logs.
    '''bash
    run_alignment \\
        1 \\
        bowtie2 \\
        global \\
        mem \\
        0 \\
        false \\
        tests/fixtures/align_fastqs/bowtie2/tiny \\
        '' \\
        tests/fixtures/align_fastqs/fastq/se/tiny_se.atria.fastq.gz NA \\
        work/tiny_se.bam \\
        false \\
        work/logs \\
        align_fastqs \\
        tiny_se
    '''
EOM
)

    if [[ "${threads}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${threads}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'threads', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    # Define paths for log output files.
    log_out="${dir_eo}/${nam_job}.${samp}.stdout.txt"
    log_err="${dir_eo}/${nam_job}.${samp}.stderr.txt"

    # Build the alignment command as an array so optional arguments can be
    # appended without word splitting paths or option values.
    cmd_aln=(
        align_fastqs
            --threads "${threads}"
            --aligner "${aligner}"
            --index "${index}"
            --fq_1 "${fq_1}"
            --fil_out "${fil_out}"
    )

    case "${aligner}" in
        bowtie2) cmd_aln+=( --bt2_mode "${bt2_mode}" ) ;;
        bwa)     cmd_aln+=( --bwa_alg "${bwa_alg}" ) ;;
    esac

    if [[ "${mapq}" -gt 0 ]]; then
        cmd_aln+=( --mapq "${mapq}" )
    fi

    if [[ "${req_flg}" == "true" ]]; then
        cmd_aln+=( --req_flg )
    fi

    if [[ "${fil_out}" == *.cram ]]; then
        cmd_aln+=( --ref_fa "${ref_fa}" )
    fi

    if [[ "${fq_2}" != "NA" ]]; then
        cmd_aln+=( --fq_2 "${fq_2}" )
    fi

    if [[ "${qname}" == "true" ]]; then
        cmd_aln+=( --qname )
    fi

    # Run alignment function with specified parameters.
    if ! \
        "${cmd_aln[@]}" > "${log_out}" 2> "${log_err}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "alignment failed for sample '${samp}'. See log: '${log_err}'."
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


# Initialize hardcoded argument variables.
function init_args_hardcoded() {
    # If true, run script in debug mode.
    debug=true
}


# Initialize argument variables, assigning default values where applicable.
function init_arg_defs() {
    env_nam="env_protocol"
    threads=4
    aligner="bowtie2"
    bt2_mode="end-to-end"
    bwa_alg="mem"
    ref_fa=""
    out_ext="bam"
    mapq=0
    req_flg=false
    index=""
    csv_fil_in=""
    dir_out=""
    qname=false
    sfx_se=""
    sfx_pe=""
    dir_eo=""
    nam_job="align_fastqs"

    unset arr_fil_in
    declare -ga arr_fil_in
}


# Initialize hardcoded arguments and user-facing argument defaults.
function init_defs() {
    init_args_hardcoded
    init_arg_defs
}


# Parse keyword arguments after helper scripts have been sourced.
function parse_args() {
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -en|--env|--env[_-]nam)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_align_fastqs
                    return 1
                }
                env_nam="${2}"
                shift 2
                ;;

            -ds|--dir[_-]scr)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_align_fastqs
                    return 1
                }
                dir_scr="${2}"
                shift 2
                ;;

            -t|--thr|--threads)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_align_fastqs
                    return 1
                }
                threads="${2}"
                shift 2
                ;;

            -a|--aln|--aligner)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_align_fastqs
                    return 1
                }
                aligner="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -2m|--bt2[_-]mode)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_align_fastqs
                    return 1
                }
                bt2_mode="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -ba|--bwa[_-]alg)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_align_fastqs
                    return 1
                }
                bwa_alg="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -mq|--mapq)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_align_fastqs
                    return 1
                }
                mapq="${2}"
                shift 2
                ;;

            -rq|--req[_-]flg)
                req_flg=true
                shift 1
                ;;

            -ix|--index)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_align_fastqs
                    return 1
                }
                index="${2}"
                shift 2
                ;;

            -ci|--csv[_-]fil[_-]in)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_align_fastqs
                    return 1
                }
                csv_fil_in="${2}"
                shift 2
                ;;

            -do|--dir[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_align_fastqs
                    return 1
                }
                dir_out="${2}"
                shift 2
                ;;

            -ox|--out[_-]ext)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_align_fastqs
                    return 1
                }
                out_ext="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -rf|--ref[_-]fa)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_align_fastqs
                    return 1
                }
                ref_fa="${2}"
                shift 2
                ;;

            -qn|--qnam|--qname)
                qname=true
                shift 1
                ;;

            -sxs|--sfx[_-]se|--suffix[_-]se)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_align_fastqs
                    return 1
                }
                sfx_se="${2}"
                shift 2
                ;;

            -sxp|--sfx[_-]pe|--suffix[_-]pe)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_align_fastqs
                    return 1
                }
                sfx_pe="${2}"
                shift 2
                ;;

            -deo|--dir[_-]eo)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_align_fastqs
                    return 1
                }
                dir_eo="${2}"
                shift 2
                ;;

            -nj|--nam[_-]job)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_align_fastqs
                    return 1
                }
                nam_job="${2}"
                shift 2
                ;;

            *)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                help_submit_align_fastqs
                return 1
                ;;
        esac
    done
}


# Validate required arguments and paths.
function validate_args() {
    validate_var     "env_nam"    "${env_nam}"         || return 1
    validate_var_dir "dir_scr"    "${dir_scr}" 0 false || return 1
    validate_var     "threads"    "${threads}"         || return 1
    validate_var     "aligner"    "${aligner}"         || return 1
    validate_var     "bt2_mode"   "${bt2_mode}"        || return 1
    validate_var     "bwa_alg"    "${bwa_alg}"         || return 1
    validate_var     "mapq"       "${mapq}"            || return 1
    validate_var     "index"      "${index}"           || return 1
    validate_var     "csv_fil_in" "${csv_fil_in}"      || return 1
    validate_var_dir "dir_out"    "${dir_out}"         || return 1
    validate_var     "sfx_se"     "${sfx_se}"          || return 1
    validate_var     "sfx_pe"     "${sfx_pe}"          || return 1
    validate_var_dir "dir_eo"     "${dir_eo}"          || return 1
    validate_var     "nam_job"    "${nam_job}"         || return 1

    check_int_pos    "${threads}" "threads" || return 1
    check_int_nonneg "${mapq}"    "mapq"    || return 1

    case "${aligner}" in
        bowtie2|bwa|bwa-mem2) : ;;
        *)
            echo_err \
                "'--aligner' must be 'bowtie2', 'bwa', or 'bwa-mem2':" \
                "'${aligner}'."
            return 1
            ;;
    esac

    case "${bt2_mode}" in
        local|global|end-to-end) : ;;
        *)
            echo_err \
                "'--bt2_mode' must be 'local', 'global', or 'end-to-end':" \
                "'${bt2_mode}'."
            return 1
            ;;
    esac

    case "${bwa_alg}" in
        mem|aln) : ;;
        *)
            echo_err \
                "'--bwa_alg' must be 'mem' or 'aln': '${bwa_alg}'."
            return 1
            ;;
    esac

    if [[ "${aligner}" == "bwa-mem2" && "${bwa_alg}" != "mem" ]]; then
        echo_err "'--aligner bwa-mem2' requires '--bwa_alg mem'."
        return 1
    fi

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
}


# Print debug argument variable assignments.
function print_state_debug() {
    if [[ "${debug}" == "true" ]]; then
        echo
        debug_var \
            "env_nam=${env_nam}" \
            "dir_scr=${dir_scr}" \
            "threads=${threads}" \
            "aligner=${aligner}" \
            "bt2_mode=${bt2_mode}" \
            "bwa_alg=${bwa_alg}" \
            "ref_fa=${ref_fa}" \
            "out_ext=${out_ext}" \
            "mapq=${mapq}" \
            "req_flg=${req_flg}" \
            "index=${index}" \
            "csv_fil_in=${csv_fil_in}" \
            "dir_out=${dir_out}" \
            "qname=${qname}" \
            "sfx_se=${sfx_se}" \
            "sfx_pe=${sfx_pe}" \
            "dir_eo=${dir_eo}" \
            "nam_job=${nam_job}"
    fi
}


# Activate environment.
function setup_env() {
    handle_env "${env_nam}" || return 1
}


# Reconstruct input FASTQ vector.
function prepare_vecs() {
    unset arr_fil_in && declare -ga arr_fil_in
    IFS=';' read -r -a arr_fil_in <<< "${csv_fil_in}"
}


# Validate reconstructed input FASTQ vector.
function validate_vecs() {
    check_arr_nonempty "arr_fil_in" "csv_fil_in" || return 1
}


# Print debug vector assignments.
function print_vecs_debug() {
    if [[ "${debug}" == "true" ]]; then
        echo "\${#arr_fil_in[@]}=${#arr_fil_in[@]}" && echo
        echo "arr_fil_in=( ${arr_fil_in[*]} )"      && echo
    fi
}


# Parse one input entry and run one alignment.
function run_job() {
    local fil_in="${1:-}"
    local fq_1 fq_2 samp fil_out

    if [[ -z "${fil_in}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'fil_in', is missing."
        return 1
    fi

    if [[ "${debug}" == "true" ]]; then debug_var "fil_in=${fil_in}"; fi

    IFS=';' read -r fq_1 fq_2 samp fil_out < <(
        parse_entry_align_fastq \
            "${fil_in}" "${sfx_se}" "${sfx_pe}" "${dir_out}" "${out_ext}"
    ) || return 1

    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "fq_1=${fq_1}" \
            "fq_2=${fq_2}" \
            "samp=${samp}" \
            "fil_out=${fil_out}"
    fi

    if ! \
        run_alignment \
            "${threads}" "${aligner}" "${bt2_mode}" "${bwa_alg}" \
            "${mapq}"    "${req_flg}" "${index}"   "${ref_fa}" \
            "${fq_1}"    "${fq_2}"    "${fil_out}" "${qname}" \
            "${dir_eo}" "${nam_job}" "${samp}"
    then
        echo_err "failed to perform alignment."
        return 1
    fi
}


# Parse Slurm task state and run one array task.
function run_job_slurm() {
    local err_dsc err_ini id_job id_tsk idx fil_in out_dsc out_ini
    local fq_1 fq_2 samp fil_out

    id_job="${SLURM_ARRAY_JOB_ID:-}"
    id_tsk="${SLURM_ARRAY_TASK_ID:-}"

    if [[ -z "${id_job}" ]]; then
        echo_err "Slurm array job ID is missing."
        return 1
    elif ! [[ "${id_tsk}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err "Slurm task ID is invalid: '${id_tsk}'."
        return 1
    elif (( id_tsk > ${#arr_fil_in[@]} )); then
        echo_err \
            "Slurm task ID '${id_tsk}' exceeds number of FASTQ entries:" \
            "'${#arr_fil_in[@]}'."
        return 1
    else
        idx=$(( id_tsk - 1 ))
    fi

    fil_in="${arr_fil_in[idx]}"

    if [[ "${debug}" == "true" ]]; then debug_var "fil_in=${fil_in}"; fi

    IFS=';' read -r fq_1 fq_2 samp fil_out < <(
        parse_entry_align_fastq \
            "${fil_in}" "${sfx_se}" "${sfx_pe}" "${dir_out}" "${out_ext}"
    ) || return 1

    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "fq_1=${fq_1}" \
            "fq_2=${fq_2}" \
            "samp=${samp}" \
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
        run_alignment \
            "${threads}" "${aligner}" "${bt2_mode}" "${bwa_alg}" \
            "${mapq}"    "${req_flg}" "${index}"   "${ref_fa}" \
            "${fq_1}"    "${fq_2}"    "${fil_out}" "${qname}" \
            "${dir_eo}" "${nam_job}" "${samp}"
    then
        echo_err "failed to perform alignment."
        return 1
    fi

    rm -f -- "${err_ini}" "${out_ini}" || {
        echo_warn \
            "failed to remove initial Slurm log files:" \
            "'${err_ini}', '${out_ini}'."
    }
}


# Dispatch one Slurm array task or local loop over input entries.
function run_jobs() {
    local idx fil_in

    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        run_job_slurm || return 1
    else
        for idx in "${!arr_fil_in[@]}"; do
            fil_in="${arr_fil_in[idx]}"
            run_job "${fil_in}" || return 1
        done
    fi
}


# Main script execution.
function main() {
    init_defs

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        help_submit_align_fastqs
        echo >&2
        return 0
    fi

    source_helpers_submit "${0##*/}" "${dir_scr}" \
        align_fastqs \
        check_args \
        check_inputs \
        check_numbers \
        format_outputs \
        handle_env \
        manage_slurm \
        process_sequences \
        help/help_submit_align_fastqs \
        || return 1

    parse_args "$@"   || return 1
    validate_args     || return 1
    print_state_debug || return 1
    setup_env         || return 1
    prepare_vecs      || return 1
    validate_vecs     || return 1
    print_vecs_debug  || return 1
    run_jobs          || return 1
}


main "$@"
