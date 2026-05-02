#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: process_sequences.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# check_seq_type
# check_string_fastqs
# get_paired_suffix
# parse_fastq_entry
# pair_fastqs
# pair_fqs


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
# shellcheck disable=SC1091
{
    _dir_src_seq="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    source "${_dir_src_seq}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_seq}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_seq}" \
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
Usage:
  check_seq_type [-h|--hlp|--help] bam

Description:
  Determine whether a BAM file appears to represent single- or paired-end sequencing data based on recognized program-group ('@PG') entries in the BAM header.

  This helper looks for '@PG' lines containing 'ID:bowtie2' or 'ID:bwa', then classifies the BAM as paired-end if those program-group lines contain paired-read indicators such as '_R2', '_r2', or '_2'. Otherwise, it reports single-end.

Positional argument:
  1  bam  <str>  Full path to the BAM file.

Returns:
  Prints 'single' if the BAM file appears to contain single-end alignments, or 'paired' if it appears to contain paired-end alignments. Returns 1 if the BAM file cannot be checked or if no recognized aligner/program-group information is found.

Dependencies:
  - Bash
  - Samtools

Notes:
  - Helper currently supports BAM input only.
  - The paired-end detection pattern is based on recognized '@PG' lines matching

      ^@PG.*ID:(bowtie2|bwa).*(_[Rr]2|_2)

Example:
  '''bash
  seq_type=\$(check_seq_type "\${HOME}/path/to/file.bam")
  echo "The sequencing type is '\${seq_type}'."
  '''

  '''txt
  The sequencing type is 'paired'.
  '''

#TODO:
  - Extend input support beyond BAM to SAM and CRAM.
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
    local infiles="${1:-}"  # Serialized string of FASTQ file paths
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
Usage:
  check_string_fastqs [-h|--hlp|--help] infiles sfx_se sfx_pe

Description:
  Validate a serialized FASTQ-entry string for expected semicolon/comma structure, file existence, and suffix consistency.

  Each semicolon-delimited entry must be either
    - one SE FASTQ path or
    - one comma-delimited PE FASTQ pair.

Positional arguments:
  1  infiles  <str>  Serialized FASTQ-entry string.
  2  sfx_se   <str>  Expected suffix for SE FASTQ files.
  3  sfx_pe   <str>  Expected suffix for PE FASTQ read-1 files.

Returns:
  0 if all FASTQ entries are valid; otherwise 1.
EOM
    )

    if [[ "${infiles}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${infiles}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'infiles', is missing."
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

    IFS=';' read -r -a arr_fq <<< "${infiles}"
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
Usage:
  get_paired_suffix [-h|--hlp|--help] file sfx_pe

Description:
  Infer the expected read-2 suffix corresponding to a user-supplied read-1 paired-end FASTQ suffix.

Positional arguments:
  1  file    <str>  Read-1 FASTQ file path.
  2  sfx_pe  <str>  User-supplied read-1 suffix.

Returns:
  Prints the inferred read-2 suffix to stdout.
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
    local infile="${1:-}"  # Input FASTQ entry
    local sfx_se="${2:-}"  # Expected suffix for SE FASTQ files
    local sfx_pe="${3:-}"  # Expected suffix for PE FASTQ read-1 files
    local fq_1             # FASTQ file #1
    local fq_2             # FASTQ file #2, or 'NA' for SE
    local samp             # Sample name derived from 'fq_1'
    local show_help        # Help message

    show_help=$(cat << EOM
Usage:
  parse_fastq_entry [-h|--hlp|--help] infile sfx_se sfx_pe

Description:
  Parse one FASTQ entry into 'fq_1', 'fq_2', and sample name.

  The input entry may represent either
    - one FASTQ path for single-end data or
    - two comma-delimited FASTQ paths for paired-end data.

  If a paired-end entry explicitly supplies FASTQ file #2, that path is preserved rather than being re-derived from 'fq_1'.

Positional arguments:
  1  infile  <str>  Input FASTQ entry.
  2  sfx_se  <str>  Expected suffix for SE FASTQ files.
  3  sfx_pe  <str>  Expected suffix for PE FASTQ read-1 files.

Returns:
  Prints a semicolon-delimited record to stdout:

    fq_1;fq_2;samp

  where 'fq_2' is 'NA' for SE data.
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

    check_string_fastqs "${infile}" "${sfx_se}" "${sfx_pe}" || return 1

    if [[ "${infile}" == *,* ]]; then
        IFS=',' read -r fq_1 fq_2 <<< "${infile}"
        unset IFS
        samp="$(basename "${fq_1%%"${sfx_pe}"}")"
    else
        fq_1="${infile}"
        fq_2="NA"
        samp="$(basename "${fq_1%%"${sfx_se}"}")"
    fi

    validate_var "fq_1" "${fq_1}" || return 1
    validate_var "samp" "${samp}" || return 1

    echo "${fq_1};${fq_2};${samp}"
}


function pair_fastqs() {
    local arg="${1:-}"  # Optional help flag
    local show_help     # Help message

    show_help=$(cat << EOM
Usage:
  pair_fastqs [-h|--hlp|--help]

Description:
  Read FASTQ paths from stdin and emit a serialized FASTQ-entry string in which sample-specific entries are separated by semicolons and, for paired-end samples, read-1/read-2 FASTQ paths are separated by commas.

  This helper is intended for use in pipelines where a stream of FASTQ paths is generated upstream, then converted here into the semicolon/comma-delimited format expected by downstream wrapper/helper functions such as 'check_string_fastqs' and 'parse_fastq_entry'.

  Paired-end detection is purely local and order-dependent: when a line is recognized as a read-1 FASTQ path, this helper assumes the very next input line is the corresponding read-2 FASTQ path for the same sample. It does not buffer, sort, or search ahead for matching mates.

Returns:
  Writes serialized FASTQ entries to stdout. Returns 1 if an expected read-2 mate is missing after a detected read-1 entry.

Notes:
  - Input is read from stdin, one FASTQ path per line.
  - For paired-end data, upstream logic must already ensure that each read-1 FASTQ path is immediately followed by its matching read-2 FASTQ path for the same sample.
  - This helper does not sort input, group by sample, or search ahead for mates.
  - Read-1 files must end in '_R1' or '_r1', optionally followed by '.atria', then by '.fastq', '.fq', '.fastq.gz', or '.fq.gz'.
  - The corresponding paired read is expected to end in '_R2' or '_r2' with the same optional '.atria' and FASTQ/FQ extension pattern.
  - Lines not matching the read-1 pattern are treated as unpaired entries and emitted as standalone semicolon-delimited records.
  - This function has an alias: 'pair_fqs'.

Example:
  '''bash
  printf "%s\n"
      "sampleA_R1.fastq.gz"
      "sampleA_R2.fastq.gz"
      "sampleB.fq.gz"
      | pair_fastqs
  '''

  '''txt
  sampleA_R1.fastq.gz,sampleA_R2.fastq.gz;sampleB.fq.gz;
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
