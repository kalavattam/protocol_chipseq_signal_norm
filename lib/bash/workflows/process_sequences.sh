#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: process_sequences.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in development and documentation.
#
# Distributed under the MIT license.


# check_seq_type
# check_string_fastqs
# get_paired_suffix
# parse_fastq_entry
# pair_fastqs
# pair_fqs
# trim_fastqs_atria


#  Require Bash >= 4.4 before defining functions
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be sourced or run under Bash >= 4.4." >&2

    if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
        return 1
    else
        exit 1
    fi
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2

    if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
        return 1
    else
        exit 1
    fi
fi

#  Source required helper functions if needed
{
    _dir_src_seq="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    # shellcheck source=lib/bash/core/source_helpers.sh
    source "${_dir_src_seq}/../core/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_seq}/../core/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_seq}/.." \
        check_inputs check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_seq
}


#TODO: audit current usage; keep this helper even if unused
#MAYBE: make this function "private"
function check_seq_type() {
    local bam="${1:-}"  # Input BAM file
    local line_pg       # Recognized '@PG' lines from BAM header
    local show_help     # Help message

show_help=$(cat << EOM
Usage
-----
  check_seq_type
    [--help] bam

  Determine whether a BAM file appears to represent single- or paired-end sequencing data based on recognized program-group ('@PG') entries in the BAM header.

  This helper looks for '@PG' lines containing 'ID:bowtie2' or 'ID:bwa', then classifies the BAM as paired-end if those program-group lines contain paired-read indicators such as '_R2', '_r2', or '_2'. Otherwise, it reports single-end.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  bam : file
    Full path to the BAM file.

Returns
-------
  Prints 'single' if the BAM file appears to contain single-end alignments, or 'paired' if it appears to contain paired-end alignments. Returns 1 if the BAM file cannot be checked or if no recognized aligner/program-group information is found.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - grep
    - samtools

  - Helper currently supports BAM input only.
  - The paired-end detection pattern is based on recognized '@PG' lines matching

      ^@PG.*ID:(bowtie2|bwa).*(_[Rr]2|_2)

Examples
--------
  1. Classify a BAM whose recognized aligner program group describes paired-end reads.
    '''bash
    check_seq_type "alignments/sample.paired.bam"
    '''

  2. Handle rejection of an existing BAM whose header has no recognized Bowtie 2 or BWA program group.
    '''bash
    if ! check_seq_type "alignments/no-aligner-pg.bam"; then
        printf '%s\n' "No supported sequencing type was found."
    fi
    '''
EOM
)

    #  Parse and check function argument
    if [[ "${bam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${bam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'bam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var_file "bam" "${bam}" || return 1

    #  Check that Samtools is available in PATH
    if ! command -v samtools > /dev/null 2>&1; then
        echo_err_func "${FUNCNAME[0]}" \
            "samtools is either not installed or not in PATH."
        return 1
    fi

    #  Extract @PG lines from BAM header and analyze sequencing type
    line_pg=$(
        samtools view -H "${bam}" | grep -E '^@PG.*ID:(bowtie2|bwa)'
    )

    if [[ -z "${line_pg}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "no recognized '@PG' entry with 'ID:bowtie2' or 'ID:bwa' in the" \
            "BAM header."
        return 1
    fi

    #  Determine sequencing type based on paired-end indicators '_R2', '_r2',
    #+ or '_2'
    if \
        printf '%s\n' "${line_pg}" | grep -E '(_[Rr]2|_2)' > /dev/null 2>&1
    then
        echo "paired"
    else
        echo "single"
    fi
}


#MAYBE: make this function "private"
function check_string_fastqs() {
    local csv_fil_in="${1:-}"  # Serialized string of FASTQ file paths
    local sfx_se="${2:-}"   # Expected suffix for SE FASTQ files
    local sfx_pe="${3:-}"   # Expected suffix for PE FASTQ read-1 files
    local -a arr_fq         # FASTQ entries split on semicolons
    local fq                # One FASTQ entry
    local fq_1              # First FASTQ file
    local fq_2              # Second FASTQ file (for PE reads)
    local sfx_pe_2          # Expected suffix for PE FASTQ read-2 files
    local num_prt           # Number of comma-delimited fields in one entry
    local show_help         # Help message

    show_help=$(cat << EOM
Usage
-----
  check_string_fastqs
    [--help] csv_fil_in sfx_se sfx_pe

  Validate a serialized FASTQ-entry string for expected semicolon/comma structure, file existence, and suffix consistency.

  Each semicolon-delimited entry must be either
    - one SE FASTQ path or
    - one comma-delimited PE FASTQ pair.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  csv_fil_in : list of structured string
    Comma-separated list of input file paths. Serialized FASTQ-entry string.

  2  sfx_se : str
    Suffix to strip from single-end FASTQ filenames.

  3  sfx_pe : str
    Suffix to strip from paired-end FASTQ read-1 filenames.

Returns
-------
  0 if all FASTQ entries are valid; otherwise 1.

Notes
-----
  Runtime requirements:
    - awk
    - bash >= 4.4

Examples
--------
  1. Validate one existing single-end FASTQ file.
    '''bash
    check_string_fastqs "reads/control.fastq.gz" ".fastq.gz" "_R1.fastq.gz"
    '''

  2. Validate a serialized mixture containing one existing paired-end entry and one existing single-end entry.
    '''bash
    csv_fil_in="reads/sample_R1.fastq.gz,reads/sample_R2.fastq.gz;reads/control.fastq.gz"
    check_string_fastqs "\${csv_fil_in}" ".fastq.gz" "_R1.fastq.gz"
    '''
EOM
    )

    if [[ "${csv_fil_in}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${csv_fil_in}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'csv_fil_in', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${sfx_se}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'sfx_se', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${sfx_pe}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'sfx_pe', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    IFS=';' read -r -a arr_fq <<< "${csv_fil_in}"
    unset IFS

    for fq in "${arr_fq[@]}"; do
        num_prt=$(awk -F',' '{ print NF }' <<< "${fq}")

        if (( num_prt > 2 )); then
            echo_err_func "${FUNCNAME[0]}" \
                "too many comma-delimited values in FASTQ entry: '${fq}'."
            return 1
        fi

        if (( num_prt == 2 )); then
            IFS=',' read -r fq_1 fq_2 <<< "${fq}"
            unset IFS

            validate_var_file "fq_1" "${fq_1}" || return 1
            validate_var_file "fq_2" "${fq_2}" || return 1

            if [[ "${fq_1}" != *"${sfx_pe}" ]]; then
                echo_err_func "${FUNCNAME[0]}" \
                    "PE FASTQ read-1 file '${fq_1}' does not match expected" \
                    "suffix '${sfx_pe}'."
                return 1
            fi

            sfx_pe_2="$(get_paired_suffix "${fq_1}" "${sfx_pe}")" || return 1

            if [[ "${fq_2}" != *"${sfx_pe_2}" ]]; then
                echo_err_func "${FUNCNAME[0]}" \
                    "PE FASTQ read-2 file '${fq_2}' does not match expected" \
                    "suffix '${sfx_pe_2}'."
                return 1
            fi
        else
            fq_1="${fq}"
            fq_2="NA"

            validate_var_file "fq_1" "${fq_1}" || return 1

            if [[ "${fq_1}" != *"${sfx_se}" ]]; then
                echo_err_func "${FUNCNAME[0]}" \
                    "SE FASTQ file '${fq_1}' does not match expected suffix" \
                    "'${sfx_se}'."
                return 1
            fi
        fi
    done

    return 0
}


#TODO: audit current usage; keep this helper even if unused
#MAYBE: make this function "private"
function get_paired_suffix() {
    local file="${1:-}"    # File to check existence
    local sfx_pe="${2:-}"  # User-supplied PE suffix (_R1.fastq.gz, etc.)
    local sfx_pe_2         # Modified suffix for second read
    local sfx_mod          # Candidate modified suffix
    local file_2           # Candidate read-2 FASTQ file
    local i                # Index for iterative replacement
    local show_help        # Help message

    show_help=$(cat << EOM
Usage
-----
  get_paired_suffix
    [--help] file sfx_pe

  Infer the expected read-2 suffix corresponding to a user-supplied read-1 paired-end FASTQ suffix.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  file : file
    Read-1 FASTQ file path.

  2  sfx_pe : str
    Suffix to strip from paired-end FASTQ read-1 filenames.

Returns
-------
  Prints the inferred read-2 suffix to stdout.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Convert an '_R1.fastq.gz' read-1 suffix to its '_R2.fastq.gz' mate suffix.
    '''bash
    get_paired_suffix "reads/sample_R1.fastq.gz" "_R1.fastq.gz"
    '''

  2. Convert a '.1.fastq.gz' read-1 suffix to its '.2.fastq.gz' mate suffix.
    '''bash
    fq_r1="reads/sample.1.fastq.gz"
    get_paired_suffix "\${fq_r1}" ".1.fastq.gz"
    '''
EOM
    )

    if [[ "${file}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${file}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'file', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${sfx_pe}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'sfx_pe', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var_file "file" "${file}" || return 1

    if [[ "${sfx_pe}" =~ ([rR])1 ]]; then
        sfx_pe_2="${sfx_pe/${BASH_REMATCH[1]}1/${BASH_REMATCH[1]}2}"
    elif [[ "${sfx_pe}" =~ 1(\.f(ast)?q) ]]; then
        sfx_pe_2="${sfx_pe/1${BASH_REMATCH[1]}/2${BASH_REMATCH[1]}}"
    elif [[ "${sfx_pe}" =~ 1(\.atria) ]]; then
        sfx_pe_2="${sfx_pe/1${BASH_REMATCH[1]}/2${BASH_REMATCH[1]}}"
    else
        i=0
        while (( i < ${#sfx_pe} )); do
            if [[ "${sfx_pe:${i}:1}" == "1" ]]; then
                sfx_mod="${sfx_pe:0:${i}}2${sfx_pe:$(( i + 1 ))}"
                file_2="${file/${sfx_pe}/${sfx_mod}}"

                if [[ -f "${file_2}" ]]; then
                    sfx_pe_2="${sfx_mod}"
                    break
                fi
            fi
            (( i++ )) || true
        done

        [[ -z "${sfx_pe_2}" ]] && sfx_pe_2="${sfx_pe/1/2}"
    fi

    echo "${sfx_pe_2}"
}


#TODO: audit current usage; keep this helper even if unused
#MAYBE: make this function "private"
function parse_fastq_entry() {
    local fil_in="${1:-}"  # Input file path. Input FASTQ entry
    local sfx_se="${2:-}"  # Expected suffix for SE FASTQ files
    local sfx_pe="${3:-}"  # Expected suffix for PE FASTQ read-1 files
    local fq_1             # FASTQ file #1
    local fq_2             # FASTQ file #2, or 'NA' for SE
    local samp             # Sample name derived from 'fq_1'
    local show_help        # Help message

    show_help=$(cat << EOM
Usage
-----
  parse_fastq_entry
    [--help] fil_in sfx_se sfx_pe

  Parse one FASTQ entry into 'fq_1', 'fq_2', and sample name.

  The input entry may represent either
    - one FASTQ path for single-end data or
    - two comma-delimited FASTQ paths for paired-end data.

  If a paired-end entry explicitly supplies FASTQ file #2, that path is preserved rather than being re-derived from 'fq_1'.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fil_in : file
    Input file path. Input FASTQ entry.

  2  sfx_se : str
    Suffix to strip from single-end FASTQ filenames.

  3  sfx_pe : str
    Suffix to strip from paired-end FASTQ read-1 filenames.

Returns
-------
  Prints a semicolon-delimited record to stdout:

    fq_1;fq_2;samp

  where 'fq_2' is 'NA' for SE data.

Notes
-----
  Runtime requirements:
    - awk
    - basename
    - bash >= 4.4

Examples
--------
  1. Parse one existing single-end FASTQ and derive the sample name by removing the SE suffix.
    '''bash
    parse_fastq_entry "reads/control.fastq.gz" ".fastq.gz" "_R1.fastq.gz"
    '''

  2. Parse an existing paired-end FASTQ entry while preserving its explicitly supplied read-2 path.
    '''bash
    fq_pair="reads/sample_R1.fastq.gz,reads/sample_R2.fastq.gz"
    parse_fastq_entry "\${fq_pair}" ".fastq.gz" "_R1.fastq.gz"
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
    elif [[ -z "${sfx_se}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'sfx_se', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${sfx_pe}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'sfx_pe', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    check_string_fastqs "${fil_in}" "${sfx_se}" "${sfx_pe}" || return 1

    if [[ "${fil_in}" == *,* ]]; then
        IFS=',' read -r fq_1 fq_2 <<< "${fil_in}"
        unset IFS
        samp="$(basename "${fq_1%%"${sfx_pe}"}")"
    else
        fq_1="${fil_in}"
        fq_2="NA"
        samp="$(basename "${fq_1%%"${sfx_se}"}")"
    fi

    validate_var "fq_1" "${fq_1}" || return 1
    validate_var "samp" "${samp}" || return 1

    echo "${fq_1};${fq_2};${samp}"
}


#  Trim one Input file path. FASTQ input entry with Atria
# TODO: Fill in any missing maintainer-facing runtime assumptions for
# 'trim_fastqs_atria()' after the Atria wrapper behavior is next reviewed.
function trim_fastqs_atria() {
    local threads="${1:-}"  # Number of threads to use
    local fq_1="${2:-}"     # FASTQ file #1
    local fq_2="${3:-}"     # FASTQ file #2, or 'NA' for SE data
    local dir_out="${4:-}"  # Output directory for trimmed FASTQs
    local log_out="${5:-}"  # Stdout log file
    local log_err="${6:-}"  # Stderr log file
    local samp="${7:-}"     # Sample name for error reporting
    local show_help         # Help message
    local -a arr_cmd_atria  # Atria command array

    show_help=$(cat << EOM
Usage
-----
  trim_fastqs_atria
    [--help] threads fq_1 fq_2 dir_out log_out log_err samp

  Trim one single- or paired-end Input file path. FASTQ input entry with Atria.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  threads : int
    Number of threads to use.

  2  fq_1 : file
    First FASTQ input file.

  3  fq_2 : file
    Second FASTQ input file, or 'NA' for single-end data.

  4  dir_out : dir
    Output directory for trimmed FASTQ files.

  5  log_out : file
    Stdout log file.

  6  log_err : file
    Stderr log file.

  7  samp : str
    Sample name used in error messages.

Returns
-------
  Runs Atria, writing stdout and stderr to the supplied log files. Returns 0 when trimming completes successfully; otherwise 1.

Notes
-----
  Runtime requirements:
    - atria
    - bash >= 4.4
    - dirname

  - '-R <fq_2>' is passed only when 'fq_2 != NA'.
  - Atria is always called with '--length-range 35:500'.

Examples
--------
  1. Trim one existing single-end FASTQ with Atria after creating the output and log directories.
    '''bash
    mkdir -p "trimmed" "logs"
    trim_fastqs_atria 4 "reads/control.fastq.gz" NA "trimmed" "logs/control.out" "logs/control.err" control
    '''

  2. Trim one existing paired-end FASTQ pair with Atria after creating the output and log directories.
    '''bash
    mkdir -p "trimmed" "logs"
    trim_fastqs_atria 4 "reads/sample_R1.fastq.gz" "reads/sample_R2.fastq.gz" "trimmed" "logs/sample.out" "logs/sample.err" sample
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
    elif [[ -z "${fq_1}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'fq_1', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${fq_2}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'fq_2', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${dir_out}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 4, 'dir_out', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${log_out}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 5, 'log_out', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${log_err}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 6, 'log_err', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${samp}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 7, 'samp', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var_file "fq_1" "${fq_1}" || return 1

    if [[ "${fq_2}" != "NA" ]]; then
        validate_var_file "fq_2" "${fq_2}" || return 1
    fi

    validate_var_dir \
        "dir_out" "${dir_out}" || return 1
    validate_var_dir \
        "log_out parent directory" "$(dirname "${log_out}")" || return 1
    validate_var_dir \
        "log_err parent directory" "$(dirname "${log_err}")" || return 1

    arr_cmd_atria=(
        atria
            -t "${threads}"
            -r "${fq_1}"
    )

    if [[ "${fq_2}" != "NA" ]]; then
        arr_cmd_atria+=( -R "${fq_2}" )
    fi

    arr_cmd_atria+=(
            -o "${dir_out}"
            --length-range 35:500
    )

    if ! "${arr_cmd_atria[@]}" > "${log_out}" 2> "${log_err}"; then
        echo_err_func "${FUNCNAME[0]}" \
            "read trimming failed for sample '${samp}'. See log:" \
            "'${log_err}'."
        return 1
    fi
}


function pair_fastqs() {
    local arg="${1:-}"  # Optional help flag
    local show_help     # Help message

    show_help=$(cat << EOM
Usage
-----
  pair_fastqs
    [--help]

  Read FASTQ paths from stdin and emit a serialized FASTQ-entry string in which sample-specific entries are separated by semicolons and, for paired-end samples, read-1/read-2 FASTQ paths are separated by commas.

  This helper is intended for use in pipelines where a stream of FASTQ paths is generated upstream, then converted here into the semicolon/comma-delimited format expected by downstream wrapper/helper functions such as 'check_string_fastqs' and 'parse_fastq_entry'.

  Paired-end detection is purely local and order-dependent: when a line is recognized as a read-1 FASTQ path, this helper assumes the very next input line is the corresponding read-2 FASTQ path for the same sample. It does not buffer, sort, or search ahead for matching mates.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

Returns
-------
  Writes serialized FASTQ entries to stdout. Returns 1 if an expected read-2 mate is missing after a detected read-1 entry.

Notes
-----
  Runtime requirements:
    - awk
    - bash >= 4.4

  - Input is read from stdin, one FASTQ path per line.
  - For paired-end data, upstream logic must already ensure that each read-1 FASTQ path is immediately followed by its matching read-2 FASTQ path for the same sample.
  - This helper does not sort input, group by sample, or search ahead for mates.
  - Read-1 files must end in '_R1' or '_r1', optionally followed by '.atria', then by '.fastq', '.fq', '.fastq.gz', or '.fq.gz'.
  - The corresponding paired read is expected to end in '_R2' or '_r2' with the same optional '.atria' and FASTQ/FQ extension pattern.
  - Lines not matching the read-1 pattern are treated as unpaired entries and emitted as standalone semicolon-delimited records.
  - This function has an alias: 'pair_fqs'.

Examples
--------
  1. Serialize one adjacent paired-end FASTQ pair followed by one single-end FASTQ.
    '''bash
    printf '%s\n' "sampleA_R1.fastq.gz" "sampleA_R2.fastq.gz" "sampleB.fq.gz" | pair_fastqs
    '''

  2. Handle rejection when a detected read-1 FASTQ is not followed by its read-2 mate.
    '''bash
    if ! printf '%s\n' "sampleA_R1.fastq.gz" | pair_fastqs; then
        printf '%s\n' "The read-2 mate is missing."
    fi
    '''
EOM
    )

    if [[ "${arg}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif (( $# > 0 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "this function does not accept options other than '-h', '--hlp'," \
            "or '--help': '${1}'."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if [[ -t 0 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "stdin is a terminal. Pipe FASTQ paths into '${FUNCNAME[0]}' or" \
            "redirect them from a file."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    awk '
        BEGIN { OFS = "" }
        {
            #  Check for lines ending with "_R1" or "_r1", optionally followed
            #+ by ".atria", then by ".fastq", ".fq", ".fastq.gz", or ".fq.gz"
            if ($0 ~ /_[Rr]1(\.atria)?\.(fastq|fq)(\.gz)?$/) {
                r1 = $0
                getline
                #  Ensure the next line is a matching "_R2" or "_r2" file
                if ($0 ~ /_[Rr]2(\.atria)?\.(fastq|fq)(\.gz)?$/) {
                    r2 = $0
                    print r1 ",", r2, ";"
                } else {
                    #  Print an error message and exit with a failure code if
                    #+ no matching "_R2" or "_r2" file is found
                    print \
                        "error(pair_fastqs): missing [Rr]2 file for " r1 \
                        > "/dev/stderr"
                    exit 1
                }
            } else {
                #  Handle unpaired files
                print $0 ";"
            }
        }
    '
}


#  Alias for 'pair_fastqs'
function pair_fqs() { pair_fastqs "$@"; }


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
