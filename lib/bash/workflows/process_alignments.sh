#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: process_alignments.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6) were
# used in design, development, and documentation, with all output reviewed,
# edited, and approved by the author.
#
# Distributed under the MIT license.


#MAYBE: further modularization of primitives?


# _validate_process_alignments
# qsort_file_alignments
# convert_alignments_bed_awk
# convert_alignments_bed_python


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
    _dir_src_pro_aln="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    # shellcheck source=lib/bash/core/source_helpers.sh
    source "${_dir_src_pro_aln}/../core/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_pro_aln}/../core/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_pro_aln}/.." \
        check_inputs check_numbers check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_pro_aln
}


#  Validate shared alignment-processing paths
function _validate_process_alignments() {
    local func="${1:-}"
    local fil_in="${2:-}"
    local fil_out="${3:-}"
    local ref_fa="${4:-}"
    local log_out="${5:-}"
    local log_err="${6:-}"
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _validate_process_alignments
    [--help] func fil_in fil_out ref_fa log_out log_err

  Validate the common input, output, and log paths used by alignment processing helpers.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  fil_in : file
    Input file path. Input BAM or CRAM file.

  3  fil_out : file
    Output file path. Output BAM or CRAM file.

  4  ref_fa : file
    Reference FASTA file. Reference FASTA required for CRAM input.

  5  log_out : file
    Stdout log file.

  6  log_err : file
    Stderr log file.

Returns
-------
  Returns 0 when validation succeeds; returns 1 if an input, output, or log path is missing or invalid.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - dirname

Examples
--------
  1. Validate an existing BAM input without a reference FASTA after creating output and log directories.
    '''bash
    mkdir -p "results" "logs"
    _validate_process_alignments \\
        qsort_file_alignments \\
        "alignments/sample.bam" \\
        "results/sample.qnam.bam" \\
        "" \\
        "logs/qsort.out" \\
        "logs/qsort.err"
    '''

  2. Validate an existing CRAM input with its required reference FASTA.
    '''bash
    mkdir -p "results" "logs"
    _validate_process_alignments \\
        qsort_file_alignments \\
        "alignments/sample.cram" \\
        "results/sample.qnam.cram" \\
        "reference/genome.fa" \\
        "logs/qsort.out" \\
        "logs/qsort.err"
    '''
EOM
    )

    #MAYBE: [[ -z "${func}" || "${threads}" =~ ^(-h|--h[e]?lp)$ ]]
    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    #MAYBE: check parameter assignments

    validate_var_file "fil_in"  "${fil_in}"  || return 1
    validate_var      "fil_out" "${fil_out}" || return 1
    validate_var      "log_out" "${log_out}" || return 1
    validate_var      "log_err" "${log_err}" || return 1

    validate_var_dir \
        "fil_out parent directory" \
        "$(dirname "${fil_out}")" \
        0 \
        true \
        || return 1

    validate_var_dir \
        "log_out parent directory" \
        "$(dirname "${log_out}")" \
        0 \
        true \
        || return 1

    validate_var_dir \
        "log_err parent directory" \
        "$(dirname "${log_err}")" \
        0 \
        true \
        || return 1

    if [[ "${fil_in,,}" == *.cram ]]; then
        if [[ -z "${ref_fa}" ]]; then
            echo_err_func "${func}" \
                "positional argument 4, 'ref_fa', is required for CRAM."
            return 1
        fi

        validate_var_file "ref_fa" "${ref_fa}" || return 1
    fi
}


#  Sort one BAM or CRAM file by query name
function qsort_file_alignments() {
    local threads="${1:-}"
    local fil_in="${2:-}"
    local fil_out="${3:-}"
    local ref_fa="${4:-}"
    local log_out="${5:-}"
    local log_err="${6:-}"
    local tmp_bam
    local show_help

    show_help=$(cat << EOM
Usage
-----
  qsort_file_alignments
    [--help] threads fil_in fil_out ref_fa log_out log_err

  Sort one BAM or CRAM alignment file by query name.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  threads : int
    Number of threads to use for samtools.

  2  fil_in : file
    Input file path. Input BAM or CRAM file.

  3  fil_out : file
    Output file path. Output queryname-sorted BAM or CRAM file.

  4  ref_fa : file
    Reference FASTA file. Reference FASTA required for CRAM input.

  5  log_out : file
    Stdout log file for the processing command.

  6  log_err : file
    Stderr log file for the processing command.

Returns
-------
  Returns 0 if the output alignment file is written successfully; otherwise 1.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - Reference FASTA and required index (when processing CRAM)
    - rm (when processing CRAM)
    - samtools

  - CRAM input is converted through a temporary BAM before being written as CRAM.

Examples
--------
  1. Queryname-sort an existing BAM without a reference FASTA.
    '''bash
    mkdir -p "results" "logs"
    qsort_file_alignments \\
        4 \\
        "alignments/sample.bam" \\
        "results/sample.qnam.bam" \\
        "" \\
        "logs/qsort.out" \\
        "logs/qsort.err"
    '''

  2. Queryname-sort an existing CRAM using its reference FASTA.
    '''bash
    mkdir -p "results" "logs"
    qsort_file_alignments \\
        4 \\
        "alignments/sample.cram" \\
        "results/sample.qnam.cram" \\
        "reference/genome.fa" \\
        "logs/qsort.out" \\
        "logs/qsort.err"
    '''
EOM
    )

    #MAYBE: [[ -z "${threads}" || "${threads}" =~ ^(-h|--h[e]?lp)$ ]]
    if [[ "${threads}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    #TODO: check parameter assignments

    check_int_pos "${threads}" "threads" || return 1
    _validate_process_alignments \
        "${FUNCNAME[0]}" \
        "${fil_in}" \
        "${fil_out}" \
        "${ref_fa}" \
        "${log_out}" \
        "${log_err}" \
        || return 1

    if [[ "${fil_in,,}" == *.cram ]]; then
        tmp_bam="${fil_out%.cram}.tmp.$$.bam"

        if ! {
            samtools view \
                -@ "${threads}" \
                -T "${ref_fa}" \
                -u \
                "${fil_in}" \
                | samtools sort \
                    -@ "${threads}" \
                    -n \
                    -o "${tmp_bam}" \
                    -

            samtools view \
                -@ "${threads}" \
                -T "${ref_fa}" \
                -O cram \
                -o "${fil_out}" \
                "${tmp_bam}"
            } > "${log_out}" 2> "${log_err}"
        then
            rm -f "${tmp_bam}"
            return 1
        fi

        rm -f "${tmp_bam}"
    else
        if ! \
            samtools sort \
                -@ "${threads}" \
                -n \
                -o "${fil_out}" \
                "${fil_in}" \
                > "${log_out}" \
                2> "${log_err}"
        then
            return 1
        fi
    fi
}


#  Convert one BAM or CRAM file to BED.GZ with AWK-based PE fragment handling
function convert_alignments_bed_awk() {
    local threads="${1:-}"
    local fil_in="${2:-}"
    local fil_out="${3:-}"
    local ref_fa="${4:-}"
    local log_out="${5:-}"
    local log_err="${6:-}"
    local show_help
    local -a arr_ref_arg=()

    show_help=$(cat << EOM
Usage
-----
  convert_alignments_bed_awk
    [--help] threads fil_in fil_out ref_fa log_out log_err

  Convert one queryname-sorted paired-end BAM or CRAM file to BED.GZ of alignment records using samtools, awk, sort, and gzip.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  threads : int
    Number of threads to use for samtools.

  2  fil_in : file
    Input file path. Input queryname-sorted BAM or CRAM file.

  3  fil_out : file
    Output file path. Output BED.GZ file.

  4  ref_fa : file
    Reference FASTA file. Reference FASTA required for CRAM input.

  5  log_out : file
    Stdout log file for the processing command.

  6  log_err : file
    Stderr log file for the processing command.

Returns
-------
  Returns 0 if the output BED.GZ file is written successfully; otherwise 1.

Notes
-----
  Runtime requirements:
    - awk
    - bash >= 4.4
    - gzip
    - Reference FASTA and required index (when processing CRAM)
    - samtools
    - sort

  - This helper assumes adjacent queryname-sorted paired-end records.
  - Use the Python converter helper for single-end data.

Examples
--------
  1. Convert an existing queryname-sorted paired-end BAM to compressed BED without a reference FASTA.
    '''bash
    mkdir -p "results" "logs"
    convert_alignments_bed_awk 4 \\
        "alignments/sample.qnam.bam" \\
        "results/sample.bed.gz" \\
        "" \\
        "logs/convert.out" \\
        "logs/convert.err"
    '''

  2. Convert an existing queryname-sorted paired-end CRAM to compressed BED using its reference FASTA.
    '''bash
    mkdir -p "results" "logs"
    convert_alignments_bed_awk 4 \\
        "alignments/sample.qnam.cram" \\
        "results/sample.bed.gz" \\
        "reference/genome.fa" \\
        "logs/convert.out" \\
        "logs/convert.err"
    '''
EOM
    )

    #MAYBE: [[ -z "${threads}" || "${threads}" =~ ^(-h|--h[e]?lp)$ ]]
    if [[ "${threads}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    #TODO: check parameter assignments

    check_int_pos "${threads}" "threads" || return 1
    _validate_process_alignments \
        "${FUNCNAME[0]}" \
        "${fil_in}" \
        "${fil_out}" \
        "${ref_fa}" \
        "${log_out}" \
        "${log_err}" \
        || return 1

    if [[ "${fil_in,,}" == *.cram ]]; then
        arr_ref_arg=( -T "${ref_fa}" )
    fi

    {
        samtools view \
            -@ "${threads}" \
            "${arr_ref_arg[@]}" \
            "${fil_in}" \
            | awk '{
                if (NR % 2 == 1) {
                    c1 = $3
                    s1 = $4
                    l1 = length($10)
                } else {
                    c2 = $3
                    s2 = $4
                    l2 = length($10)

                    if (c1 == c2) {
                        st = (s1 < s2) ? s1 : s2
                        en = (s1 < s2) ? s2 + l2 - 1 : s1 + l1 - 1
                        lf = en - st + 1

                        print c1, st, en, lf
                    }
                }
            }' OFS='\t' \
            | sort -k1,1 -k2,2n \
            | gzip \
                > "${fil_out}"
    } > "${log_out}" 2> "${log_err}"
}


#  Convert one BAM or CRAM file to BED.GZ with a Python converter script
function convert_alignments_bed_python() {
    local pth_scr_py="${1:-}"
    local fil_in="${2:-}"
    local fil_out="${3:-}"
    local ref_fa="${4:-}"
    local log_out="${5:-}"
    local log_err="${6:-}"
    local show_help
    local -a arr_ref_arg=()

    show_help=$(cat << EOM
Usage
-----
  convert_alignments_bed_python
    [--help] pth_scr_py fil_in fil_out ref_fa log_out log_err

  Convert one BAM or CRAM file to BED.GZ of alignment records using a Python converter script.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  pth_scr_py : file
    Python converter script.

  2  fil_in : file
    Input file path. Input BAM or CRAM file.

  3  fil_out : file
    Output file path. Output BED.GZ file.

  4  ref_fa : file
    Reference FASTA file. Reference FASTA required for CRAM input.

  5  log_out : file
    Stdout log file for the processing command.

  6  log_err : file
    Stderr log file for the processing command.

Returns
-------
  Returns 0 if the output BED.GZ file is written successfully; otherwise 1.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - python >= 3.11
    - Reference FASTA and required index (when processing CRAM)

  - CRAM input is passed to the Python converter with '--ref_fa'.

Examples
--------
  1. Convert an existing BAM with the repository Python converter and no reference FASTA.
    '''bash
    mkdir -p "results" "logs"
    convert_alignments_bed_python \\
        "compute_signal" \\
        "alignments/sample.bam" \\
        "results/sample.bed.gz" \\
        "" \\
        "logs/convert.out" \\
        "logs/convert.err"
    '''

  2. Convert an existing CRAM with the repository Python converter and its reference FASTA.
    '''bash
    mkdir -p "results" "logs"
    convert_alignments_bed_python \\
        "compute_signal" \\
        "alignments/sample.cram" \\
        "results/sample.bed.gz" \\
        "reference/genome.fa" \\
        "logs/convert.out" \\
        "logs/convert.err"
    '''
EOM
    )

    #MAYBE: [[ -z "${pth_scr_py}" || "${pth_scr_py}" =~ ^(-h|--h[e]?lp)$ ]]
    if [[ "${pth_scr_py}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    #TODO: check parameter assignments

    if [[ "${pth_scr_py}" == */* || "${pth_scr_py}" == *.py ]]; then
        validate_var_file "pth_scr_py" "${pth_scr_py}" || return 1
    else
        check_pgrm_path "${pth_scr_py}" || return 1
    fi
    _validate_process_alignments \
        "${FUNCNAME[0]}" \
        "${fil_in}" \
        "${fil_out}" \
        "${ref_fa}" \
        "${log_out}" \
        "${log_err}" \
        || return 1

    if [[ "${fil_in,,}" == *.cram ]]; then
        arr_ref_arg=( --ref_fa "${ref_fa}" )
    fi

    if [[ "${pth_scr_py}" == */* || "${pth_scr_py}" == *.py ]]; then
        python "${pth_scr_py}" \
            --fil_in "${fil_in}" \
            --fil_out "${fil_out}" \
            "${arr_ref_arg[@]}" \
            > "${log_out}" \
            2> "${log_err}"
    else
        "${pth_scr_py}" \
            --fil_in "${fil_in}" \
            --fil_out "${fil_out}" \
            "${arr_ref_arg[@]}" \
            > "${log_out}" \
            2> "${log_err}"
    fi
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
