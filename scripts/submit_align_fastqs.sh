#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: submit_align_fastqs.sh
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
#  Parse one FASTQ entry into 'fq_1', 'fq_2', 'samp', and BAM 'outfile'
function parse_entry_align_fastq() {
    local infile="${1:-}"      # Input FASTQ file(s)
    local sfx_se="${2:-}"      # Suffix for SE FASTQ files
    local sfx_pe="${3:-}"      # Suffix for PE FASTQ files (FASTQ #1)
    local dir_out="${4:-}"     # Directory for output files
    local out_ext="${5:-bam}"  # Alignment output extension
    local fq_1                 # FASTQ file #1
    local fq_2                 # FASTQ file #2, or 'NA' for SE
    local samp                 # Sample name
    local outfile              # Output alignment file
    local show_help            # Help message

    show_help=$(cat << EOM
Usage:
  parse_entry_align_fastq
    [-h|--hlp|--help] infile sfx_se sfx_pe dir_out [out_ext]

Description:
  Parse one FASTQ input entry into 'fq_1', 'fq_2', 'samp', and alignment 'outfile'.

  The input entry may represent either:
    - single-end data: one FASTQ path
    - paired-end data: two comma-delimited FASTQ paths

  This helper calls 'process_sequences.sh::parse_fastq_entry' to parse the FASTQ input and then appends the alignment output path:

    \${dir_out}/\${samp}.\${out_ext}

Positional arguments:
  1  infile   <str>  FASTQ input entry. For SE data, this is one FASTQ file. For PE data, this is a comma-delimited FASTQ pair.
  2  sfx_se   <str>  Suffix to strip from SE FASTQ filenames when deriving the sample name.
  3  sfx_pe   <str>  Suffix to strip from PE FASTQ read-1 filenames when deriving the sample name.
  4  dir_out  <str>  Output directory used to construct the alignment outfile path.
  5  out_ext  <str>  Final alignment output extension: 'bam' or 'cram' (default: bam).

Returns:
  Prints a semicolon-delimited record to stdout:

    fq_1;fq_2;samp;outfile

  where 'fq_2' is 'NA' for SE data.

Notes:
  - This helper validates required inputs with 'validate_var'.
  - Parsing of the FASTQ entry itself is delegated to 'parse_fastq_entry'.
  - The output BAM path is always constructed as '\${dir_out}/\${samp}.\${out_ext}'.
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

    validate_var "infile"  "${infile}"  || return 1
    validate_var "sfx_se"  "${sfx_se}"  || return 1
    validate_var "sfx_pe"  "${sfx_pe}"  || return 1
    validate_var "dir_out" "${dir_out}" || return 1

    IFS=';' read -r fq_1 fq_2 samp < <(
        parse_fastq_entry "${infile}" "${sfx_se}" "${sfx_pe}"
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

    outfile="${dir_out}/${samp}.${out_ext}"

    echo "${fq_1};${fq_2};${samp};${outfile}"
}
#MAYBE: 'parse_entry_align_fastq' still uses 'validate_var' rather than file-
#      aware helpers such as 'validate_file' and 'validate_var_file'; revisit
#      if stronger file validation is wanted here


#  Execute alignment using function 'align_fastqs.sh::align_fastqs'
function run_alignment() {
    local threads="${1:-}"
    local aligner="${2:-}"
    local bt2_aln="${3:-}"
    local bwa_alg="${4:-}"
    local mapq="${5:-}"
    local req_flg="${6:-}"
    local index="${7:-}"
    local ref_fa="${8:-}"
    local fq_1="${9:-}"
    local fq_2="${10:-}"
    local outfile="${11:-}"
    local qname="${12:-}"
    local err_out="${13:-}"
    local nam_job="${14:-}"
    local samp="${15:-}"
    local log_out log_err
    local show_help

    show_help=$(cat << EOM
Usage:
  run_alignment
    [-h|--hlp|--help] threads aligner bt2_aln bwa_alg mapq req_flg index ref_fa
    fq_1 fq_2 outfile qname err_out nam_job samp

Description:
  Construct log-file paths and then run 'align_fastqs.sh::align_fastqs', writing stdout to

    \${err_out}/\${nam_job}.\${samp}.stdout.txt

  and stderr to

    \${err_out}/\${nam_job}.\${samp}.stderr.txt

Positional arguments:
   1  threads   <int>  Number of threads.
   2  aligner   <str>  Aligner program: 'bowtie2', 'bwa', or 'bwa-mem2'.
   3  bt2_aln   <str>  Bowtie 2 alignment type.
   4  bwa_alg   <str>  BWA algorithm: 'mem' or 'aln'.
   5  mapq      <int>  MAPQ threshold. If 'mapq > 0', '--mapq' is passed to 'align_fastqs'; otherwise it is omitted.
   6  req_flg   <flg>  If 'true', pass '--req_flg' to require flag bit 2.
   7  index     <str>  Path/prefix for the alignment index.
   8  ref_fa    <str>  Reference FASTA for CRAM output; ignored for BAM output.
   9  fq_1      <str>  FASTQ file #1.
  10  fq_2      <str>  FASTQ file #2 or, for SE data, 'NA'.
  11  outfile   <str>  Output alignment file; must end in '.bam' or '.cram'.
  12  qname     <flg>  If 'true', pass '--qname' to retain a queryname-sorted intermediate.
  13  err_out   <str>  Directory for stderr/stdout log files.
  14  nam_job   <str>  Job name used in log file naming.
  15  samp      <str>  Sample name used in log file naming.

Notes:
  - This helper is a thin wrapper around 'align_fastqs'.
  - '--fq_2' is passed only when 'fq_2 != NA'.
  - '--mapq' is passed only when 'mapq > 0'.
  - '--req_flg' and '--qname' are passed only when their values are 'true'.
  - '--ref' is passed only when 'outfile' ends in '.cram'.
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

    #  Define paths for log output files
    log_out="${err_out}/${nam_job}.${samp}.stdout.txt"
    log_err="${err_out}/${nam_job}.${samp}.stderr.txt"

    #  Run alignment function with specified parameters
    # shellcheck disable=SC2046
    if ! \
        align_fastqs \
            --threads "${threads}" \
            --aligner "${aligner}" \
            --bt2_aln "${bt2_aln}" \
            --bwa_alg "${bwa_alg}" \
            $(if [[ "${mapq}" -gt 0 ]]; then echo "--mapq ${mapq}"; fi) \
            $(if [[ "${req_flg}" == "true" ]]; then echo "--req_flg"; fi) \
            --index "${index}" \
            $(if [[ "${outfile}" == *.cram ]]; then echo "--ref ${ref_fa}"; fi) \
            --fq_1 "${fq_1}" \
            $(if [[ "${fq_2}" != "NA" ]]; then echo "--fq_2 ${fq_2}"; fi) \
            --outfile "${outfile}" \
            $(if [[ "${qname}" == "true" ]]; then echo "--qname"; fi) \
                 > "${log_out}" \
                2> "${log_err}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "alignment failed for sample '${samp}'. See log: '${log_err}'."
        return 1
    fi
    #TODO: 'fq_2' handling is vulnerable to word splitting if paths contain
    #+     spaces
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

            -a|--aln|--aligner)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                aligner="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -2a|-bn|--bt2[_-]aln|--bt2[_-]mode|--bowtie2[_-]aln|--bowtie2[_-]mode)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                bt2_aln="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -ba|--bwa[_-]alg)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                bwa_alg="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -r|--ref)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                ref_fa="${2}"
                shift 2
                ;;

            -mq|--mapq)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
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
                    show_help_main
                    return 1
                }
                index="${2}"
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

            -ox|--out[_-]ext)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                out_ext="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -qn|--qname)
                qname=true
                shift 1
                ;;

            -sxs|--sfx[_-]se|--suffix[_-]se)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                sfx_se="${2}"
                shift 2
                ;;

            -sxp|--sfx[_-]pe|--suffix[_-]pe)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                sfx_pe="${2}"
                shift 2
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
    validate_var     "aligner"    "${aligner}"         || return 1
    validate_var     "bt2_aln"    "${bt2_aln}"         || return 1
    validate_var     "bwa_alg"    "${bwa_alg}"         || return 1
    validate_var     "mapq"       "${mapq}"            || return 1
    validate_var     "index"      "${index}"           || return 1
    validate_var     "csv_infile" "${csv_infile}"      || return 1
    validate_var_dir "dir_out"    "${dir_out}"         || return 1
    validate_var     "sfx_se"     "${sfx_se}"          || return 1
    validate_var     "sfx_pe"     "${sfx_pe}"          || return 1
    validate_var_dir "err_out"    "${err_out}"         || return 1
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

    case "${bt2_aln}" in
        local|global|end-to-end) : ;;
        *)
            echo_err \
                "'--bt2_aln' must be 'local', 'global', or 'end-to-end':" \
                "'${bt2_aln}'."
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


#  Print debug argument variable assignments
function print_debug_args() {
    if [[ "${debug}" == "true" ]]; then
        echo
        debug_var \
            "env_nam=${env_nam}" \
            "dir_scr=${dir_scr}" \
            "threads=${threads}" \
            "aligner=${aligner}" \
            "bt2_aln=${bt2_aln}" \
            "bwa_alg=${bwa_alg}" \
            "ref_fa=${ref_fa}" \
            "out_ext=${out_ext}" \
            "mapq=${mapq}" \
            "req_flg=${req_flg}" \
            "index=${index}" \
            "csv_infile=${csv_infile}" \
            "dir_out=${dir_out}" \
            "qname=${qname}" \
            "sfx_se=${sfx_se}" \
            "sfx_pe=${sfx_pe}" \
            "err_out=${err_out}" \
            "nam_job=${nam_job}"
    fi
}


#  Initialize argument variables, assigning default values where applicable
env_nam="env_protocol"
dir_scr=""
threads=4
aligner="bowtie2"
bt2_aln="end-to-end"
bwa_alg="mem"
ref_fa=""
out_ext="bam"
mapq=0
req_flg=false
index=""
csv_infile=""
dir_out=""
qname=false
sfx_se=""
sfx_pe=""
err_out=""
nam_job="align_fastqs"


function show_help_main() {
    cat << EOM >&2
Usage:
  submit_align_fastqs.sh
    [-h|--hlp|--help] [-en|--env_nam <str>] -ds|--dir_scr <str>
    [-t|--threads <int>] [-a|--aligner <str>] [-2a|-bn|--bt2_aln <str>] [-ba|--bwa_alg <str>] [-mq|--mapq <int>] [-rf|--req_flg] -ix|--index <str> [-r|--ref <str>]
    -i|--csv_infile <str> -do|--dir_out <str> [-ox|--out_ext <str>] [-qn|--qname] -sxs|--sfx_se <str> -sxp|--sfx_pe <str> -eo|--err_out <str> [-nj|--nam_job <str>]

Description:
  Submit or execute one or more FASTQ-alignment jobs by calling the downstream function 'align_fastqs'.

  This wrapper
    - parses a semicolon-delimited list of FASTQ input entries,
    - derives sample names and BAM outfile paths,
    - activates the requested Conda/Mamba environment, and then
    - runs alignment either under Slurm array execution or by serial/GNU-Parallel-style iteration, depending on how the script is invoked.

  For each input entry, this script writes log files to:

    \${err_out}/\${nam_job}.\${samp}.stdout.txt
    \${err_out}/\${nam_job}.\${samp}.stderr.txt

Keyword arguments:
  -en, --env_nam  <str>
    Conda/Mamba environment to activate.

  -ds, --dir_scr  <str>
    Directory containing scripts and functions.

  -t, --threads  <int>
    Number of threads to use.

  -a, --aligner  <str>
    Alignment program to use: 'bowtie2', 'bwa', or 'bwa-mem2'.

  -2a, -bn, --bt2_aln, --bowtie2_aln, --bowtie2_mode  <str>
    Bowtie 2 alignment type: 'local', 'global', or 'end-to-end'.

  -ba, --bwa_alg  <str>
    BWA algorithm when '--aligner bwa': 'mem' or 'aln'.

  -mq, --mapq  <int>
    MAPQ threshold for filtering alignment outfiles.

  -rf, --req_flg  <flg>
    Require flag bit 2 (proper PE alignments) for filtering alignment outfiles.

  -ix, --index  <str>
    Path to the aligner index/reference.

  -r, --ref  <str>
    Reference FASTA path required when '--out_ext cram'.

  -i, --csv_infile  <str>
    Semicolon-delimited serialized string of FASTQ input entries. For single-end data, each entry is one FASTQ file. For paired-end data, each entry contains a comma-delimited FASTQ pair.

    Compatibility aliases include '--infile', '--infiles', '--fil_in', and '--csv_infiles'.

  -do, --dir_out  <str>
    Directory to write alignment outfiles.

  -ox, --out_ext  <str>
    Final alignment output extension: 'bam' or 'cram'.

  -qn, --qname  <flg>
    Retain queryname-sorted intermediate BAM files.

  -sxs, --sfx_se, --suffix_se  <str>
    Suffix to strip from SE FASTQ files.

  -sxp, --sfx_pe, --suffix_pe  <str>
    Suffix to strip from PE FASTQ files.

  -eo, --err_out  <str>
    Directory to store stderr and stdout outfiles.

  -nj, --nam_job  <str>
    Name of job.

Notes:
  - All arguments are required with the following notes and exceptions:
    + '--env_nam' defaults to 'env_nam=${env_nam}' if not specified.
    + '--threads' defaults to 'threads=${threads}' if not specified.
    + '--aligner' defaults to 'aligner=${aligner}' if not specified.
    + '--bt2_aln' defaults to 'bt2_aln=${bt2_aln}' if not specified.
    + '--bwa_alg' defaults to 'bwa_alg=${bwa_alg}' if not specified.
    + '--out_ext' defaults to 'out_ext=${out_ext}' if not specified.
    + '--ref' is required when '--out_ext cram'.
    + '--mapq' defaults to 'mapq=${mapq}' if not specified.
    + '--req_flg' and '--qname' are optional flags.
    + '--nam_job' defaults to 'nam_job=${nam_job}' if not specified.
EOM
}


#  Main script execution
function main() {
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        show_help_main
        echo >&2
        exit 0
    fi

    dir_scr="$(resolve_dir_scr "${0##*/}" "$@")" || exit 1

    source_submit_helpers "${0##*/}" "${dir_scr}" \
        align_fastqs \
        check_args \
        check_inputs \
        check_numbers \
        format_outputs \
        handle_env \
        manage_slurm \
        process_sequences \
        || exit 1

    parse_args "$@"  || exit 1
    validate_args    || exit 1
    print_debug_args || exit 1

    #  Activate environment
    handle_env "${env_nam}" || exit 1

    IFS=';' read -r -a arr_infile <<< "${csv_infile}"
    check_arr_nonempty "arr_infile" "csv_infile" || exit 1

    if [[ "${debug}" == "true" ]]; then
        echo "\${#arr_infile[@]}=${#arr_infile[@]}" && echo
        echo "arr_infile=( ${arr_infile[*]} )"      && echo
    fi

    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        id_job="${SLURM_ARRAY_JOB_ID}"
        id_tsk="${SLURM_ARRAY_TASK_ID}"

        if ! [[ "${id_tsk}" =~ ^[1-9][0-9]*$ ]]; then
            echo_err "Slurm task ID is invalid: '${id_tsk}'."
            exit 1
        elif (( id_tsk > ${#arr_infile[@]} )); then
            echo_err \
                "Slurm task ID '${id_tsk}' exceeds number of FASTQ entries:" \
                "'${#arr_infile[@]}'."
            exit 1
        else
            idx=$(( id_tsk - 1 ))
        fi

        infile="${arr_infile[idx]}"

        if [[ "${debug}" == "true" ]]; then debug_var "infile=${infile}"; fi

        IFS=';' read -r fq_1 fq_2 samp outfile < <(
            parse_entry_align_fastq \
                "${infile}" "${sfx_se}" "${sfx_pe}" "${dir_out}" "${out_ext}"
        ) || exit 1

        if [[ "${debug}" == "true" ]]; then
            debug_var \
                "fq_1=${fq_1}" \
                "fq_2=${fq_2}" \
                "samp=${samp}" \
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
            run_alignment \
                "${threads}" "${aligner}" "${bt2_aln}" "${bwa_alg}" \
                "${mapq}"    "${req_flg}" "${index}"   "${ref_fa}" \
                "${fq_1}"    "${fq_2}"    "${outfile}" "${qname}" \
                "${err_out}" "${nam_job}" "${samp}"
        then
            echo_err "failed to perform alignment."
            exit 1
        fi

        rm "${err_ini}" "${out_ini}"
    else
        for idx in "${!arr_infile[@]}"; do
            infile="${arr_infile[idx]}"

            if [[ "${debug}" == "true" ]]; then debug_var "infile=${infile}"; fi

            IFS=';' read -r fq_1 fq_2 samp outfile < <(
                parse_entry_align_fastq \
                    "${infile}" "${sfx_se}" "${sfx_pe}" "${dir_out}" \
                    "${out_ext}"
            ) || exit 1

            if [[ "${debug}" == "true" ]]; then
                debug_var \
                    "fq_1=${fq_1}" \
                    "fq_2=${fq_2}" \
                    "samp=${samp}" \
                    "outfile=${outfile}"
            fi

            if ! \
                run_alignment \
                    "${threads}" "${aligner}" "${bt2_aln}" "${bwa_alg}" \
                    "${mapq}"    "${req_flg}" "${index}"   "${ref_fa}" \
                    "${fq_1}"    "${fq_2}"    "${outfile}" "${qname}" \
                    "${err_out}" "${nam_job}" "${samp}"
            then
                echo_err "failed to perform alignment."
                exit 1
            fi
        done
    fi
}


main "$@"
