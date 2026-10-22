#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: process_alignments.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6);
# - Anthropic Claude Code (Opus 5.5).
#
# Distributed under the MIT license.


# MAYBE: further modularization of primitives?


# _validate_process_alignments
# qsort_file_alignments
# convert_alignments_bed_awk
# convert_alignments_bed_python


# Require Bash >= 4.4 before defining functions.
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

# Source required helper functions if needed.
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
        check_inputs \
        check_numbers \
        check_source \
        format_outputs \
        record_alignment || {
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


# Validate shared alignment-processing paths.
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
    Input BAM or CRAM file.

  3  fil_out : file
    Output BAM or CRAM file.

  4  ref_fa : file
    Reference FASTA file. Required for CRAM input; otherwise, ignored with a warning.

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

    # MAYBE: [[ -z "${func}" || "${threads}" =~ ^(-h|--h[e]?lp)$ ]]
    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    # MAYBE: check parameter assignments.

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
            echo_err_func "${func}" "'ref_fa' is required for CRAM input."
            return 1
        fi

        validate_var_file "ref_fa" "${ref_fa}" || return 1
    elif [[ -n "${ref_fa}" ]]; then
        echo_warn_func "${func}" \
            "'ref_fa' has no effect without CRAM input and is ignored."
    fi
}


# Build the pipeline record ('OUTPUT.ALIGNMENT.RECORD') of one 'qsort_bam' sort
# as '@CO' header lines, content lines first.
function _build_qsort_record() {
    local rec_nam="${1:-}"
    local threads="${2:-}"
    local fil_in="${3:-}"
    local fil_out="${4:-}"
    local ref_fa="${5:-}"
    local -a arr_cnt arr_run arr_pth
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _build_qsort_record
    [--help] rec_nam threads fil_in fil_out ref_fa

  Fill the named array with the '@CO' lines of the pipeline record of one 'qsort_bam' sort: 'pipeline:' lines for what determines the file's content, then 'pipeline-run:' lines for the run's context.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  rec_nam : str
    Name of the indexed array to fill, one header line per element.

  2  threads : int
    Number of threads to use.

  3  fil_in : file
    Input BAM or CRAM file.

  4  fil_out : file
    Output queryname-sorted BAM or CRAM file.

  5  ref_fa : file
    Reference FASTA, recorded only for CRAM input, or empty.

Returns
-------
  Returns 0 after filling the array; prints an error and returns 1 if 'rec_nam' is missing.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - date
    - dirname
    - git (optional; without it, the commit is 'unknown')
    - hostname

Examples
--------
  1. Build the record for a BAM sort and print it.
    '''bash
    _build_qsort_record arr_rec 4 \\
        alignments/sample.bam results/sample.qnam.bam ''
    printf '%s\\n' "\${arr_rec[@]}"
    '''

  2. Display the helper's argument contract.
    '''bash
    _build_qsort_record --help
    '''
EOM
    )

    if [[ "${rec_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${rec_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'rec_nam', is missing."
        return 1
    fi

    local -n rec_ref="${rec_nam}"

    # Content lines: the chain, the commit, and the input by file name.
    arr_cnt=(
        $'@CO\tpipeline: chain=qsort_bam'
        $'@CO\tpipeline: commit='"$(print_record_commit)"
        $'@CO\tpipeline: --fil_in='"${fil_in##*/}"
    )

    arr_pth=( fil_in "${fil_in}" )
    if [[ "${fil_in,,}" == *.cram ]]; then
        arr_cnt+=( $'@CO\tpipeline: --ref_fa='"${ref_fa##*/}" )
        arr_pth+=( ref_fa "${ref_fa}" )
    fi
    arr_pth+=( fil_out "${fil_out}" )

    # Run lines: the thread count, the start time, the host, the caller, the
    # environment, any Slurm job, and the full paths.
    arr_run=()
    add_record_run_lines \
        arr_run "${BASH_SOURCE[0]##*/}::qsort_file_alignments" "${threads}" \
        "${arr_pth[@]}" \
        || return 1

    # shellcheck disable=SC2034  # Read by name in the caller.
    rec_ref=( "${arr_cnt[@]}" "${arr_run[@]}" )
}


# Write a pipeline record into the header of a BAM file.
function _write_qsort_record() {
    local fil_bam="${1:-}"
    local rec_nam="${2:-}"
    local hdr_rec bam_rec
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _write_qsort_record
    [--help] fil_bam rec_nam

  Write the pipeline record in the named array into the header of a BAM file as '@CO' lines, after the header's existing lines.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fil_bam : file
    BAM file to rewrite in place.

  2  rec_nam : str
    Name of the indexed array holding the record, one header line per element.

Returns
-------
  Returns 0 when the record is written; otherwise 1.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - mv
    - rm
    - samtools

  - The header is read with '--no-PG', so reading it adds no '@PG' line; 'samtools reheader' adds its own.

Examples
--------
  1. Write a record built by '_build_qsort_record()' into a sorted BAM.
    '''bash
    _write_qsort_record results/sample.qnam.bam arr_rec
    '''

  2. Display the helper's argument contract.
    '''bash
    _write_qsort_record --help
    '''
EOM
    )

    if [[ "${fil_bam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    local -n rec_wrt_ref="${rec_nam}"
    hdr_rec="${fil_bam%.bam}.header.sam"
    bam_rec="${fil_bam%.bam}.record.bam"

    if ! \
        {
            samtools view -H --no-PG "${fil_bam}" \
                && printf '%s\n' "${rec_wrt_ref[@]}"
        } > "${hdr_rec}"
    then
        rm -f "${hdr_rec}"
        return 1
    fi

    if ! \
        samtools reheader "${hdr_rec}" "${fil_bam}" > "${bam_rec}"
    then
        rm -f "${hdr_rec}" "${bam_rec}"
        return 1
    fi

    rm -f "${hdr_rec}"
    mv -f "${bam_rec}" "${fil_bam}"
}


# Sort one BAM or CRAM file by query name.
function qsort_file_alignments() {
    local threads="${1:-}"
    local fil_in="${2:-}"
    local fil_out="${3:-}"
    local ref_fa="${4:-}"
    local log_out="${5:-}"
    local log_err="${6:-}"
    local tmp_bam

    # shellcheck disable=SC2034  # Filled and read by name in the helpers.
    local -a arr_rec
    local show_help

    show_help=$(cat << EOM
Usage
-----
  qsort_file_alignments
    [--help] threads fil_in fil_out ref_fa log_out log_err

  Sort one BAM or CRAM alignment file by query name, and write the chain's pipeline record into the output's header.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  threads : int
    Number of threads to use for samtools.

  2  fil_in : file
    Input BAM or CRAM file.

  3  fil_out : file
    Output queryname-sorted BAM or CRAM file.

  4  ref_fa : file
    Reference FASTA file. Required for CRAM input; otherwise, ignored with a warning.

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
  - The output keeps the input's header, any earlier pipeline record included, and adds this chain's record as '@CO' lines:
    + 'pipeline:' lines give the commit and the input file name;
    + 'pipeline-run:' lines give the thread count, the start time, the host, the scripts and functions that ran, and the full paths.

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

    # MAYBE: [[ -z "${threads}" || "${threads}" =~ ^(-h|--h[e]?lp)$ ]]
    if [[ "${threads}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    # TODO: check parameter assignments.

    check_int_pos "${threads}" "threads" || return 1
    _validate_process_alignments \
        "${FUNCNAME[0]}" \
        "${fil_in}" \
        "${fil_out}" \
        "${ref_fa}" \
        "${log_out}" \
        "${log_err}" \
        || return 1

    _build_qsort_record \
        arr_rec "${threads}" "${fil_in}" "${fil_out}" "${ref_fa}" \
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
                    - \
                && _write_qsort_record "${tmp_bam}" arr_rec \
                && samtools view \
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
        if ! {
            samtools sort \
                -@ "${threads}" \
                -n \
                -o "${fil_out}" \
                "${fil_in}" \
                && _write_qsort_record "${fil_out}" arr_rec
            } > "${log_out}" 2> "${log_err}"
        then
            return 1
        fi
    fi
}


# Convert one BAM or CRAM file to BED.GZ with AWK-based PE fragment handling.
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
    Input queryname-sorted BAM or CRAM file.

  3  fil_out : file
    Output BED.GZ file.

  4  ref_fa : file
    Reference FASTA file. Required for CRAM input; otherwise, ignored with a warning.

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
  - Records flagged unaligned (SAM flag 4) or with an unaligned mate (SAM flag 8) are dropped first ('samtools view -F 12'), so a file that keeps unaligned reads gives no fragment for them.
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

    # MAYBE: [[ -z "${threads}" || "${threads}" =~ ^(-h|--h[e]?lp)$ ]]
    if [[ "${threads}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    # TODO: check parameter assignments.

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

    # Drop unaligned reads and their mates, which have no fragment.
    {
        samtools view \
            -@ "${threads}" \
            -F 12 \
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


# Convert one BAM or CRAM file to BED.GZ with a Python converter script.
function convert_alignments_bed_python() {
    local threads="${1:-}"
    local pth_scr_py="${2:-}"
    local fil_in="${3:-}"
    local fil_out="${4:-}"
    local ref_fa="${5:-}"
    local log_out="${6:-}"
    local log_err="${7:-}"
    local show_help
    local -a arr_ref_arg=()

    show_help=$(cat << EOM
Usage
-----
  convert_alignments_bed_python
    [--help] threads pth_scr_py fil_in fil_out ref_fa log_out log_err

  Convert one BAM or CRAM file to BED.GZ of alignment records using a Python converter script.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  threads : int
    Number of threads to pass to the Python converter.

  2  pth_scr_py : file
    Python converter script.

  3  fil_in : file
    Input BAM or CRAM file.

  4  fil_out : file
    Output BED.GZ file.

  5  ref_fa : file
    Reference FASTA file. Required for CRAM input; otherwise, ignored with a warning.

  6  log_out : file
    Stdout log file for the processing command.

  7  log_err : file
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

  - The Python converter is called with '--fil_in', '--fil_out', and '--threads', and with '--ref_fa' for CRAM input.

Examples
--------
  1. Convert an existing BAM with the repository Python converter and no reference FASTA.
    '''bash
    mkdir -p "results" "logs"
    convert_alignments_bed_python \\
        4 \\
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
        4 \\
        "compute_signal" \\
        "alignments/sample.cram" \\
        "results/sample.bed.gz" \\
        "reference/genome.fa" \\
        "logs/convert.out" \\
        "logs/convert.err"
    '''
EOM
    )

    # MAYBE: [[ -z "${threads}" || "${threads}" =~ ^(-h|--h[e]?lp)$ ]]
    if [[ "${threads}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    # TODO: check parameter assignments.

    check_int_pos "${threads}" "threads" || return 1
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
            --threads "${threads}" \
            "${arr_ref_arg[@]}" \
            > "${log_out}" \
            2> "${log_err}"
    else
        "${pth_scr_py}" \
            --fil_in "${fil_in}" \
            --fil_out "${fil_out}" \
            --threads "${threads}" \
            "${arr_ref_arg[@]}" \
            > "${log_out}" \
            2> "${log_err}"
    fi
}


# Print an error message when function script is executed directly.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
