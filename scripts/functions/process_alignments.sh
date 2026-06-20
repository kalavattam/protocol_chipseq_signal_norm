#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: process_alignments.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5-series models) were used in development.
#
# Distributed under the MIT license.


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
# shellcheck disable=SC1091
{
    _dir_src_pro_aln="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    source "${_dir_src_pro_aln}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_pro_aln}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_pro_aln}" \
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
Usage:
  qsort_file_alignments [-h|--hlp|--help] threads fil_in fil_out ref_fa log_out log_err

Description:
  Sort one BAM or CRAM alignment file by query name.

Positional arguments:
  1  threads  <int>
    Number of threads to pass to samtools.

  2  fil_in  <file>
    Input BAM or CRAM file.

  3  fil_out  <file>
    Output queryname-sorted BAM or CRAM file.

  4  ref_fa  <file>
    Reference FASTA. Required for CRAM input.

  5  log_out  <file>
    File to receive stdout from the processing command.

  6  log_err  <file>
    File to receive stderr from the processing command.

Dependencies:
  External programs:
    - samtools

Returns:
  Returns 0 if the output alignment file is written successfully; otherwise 1.

Notes:
  - CRAM input is converted through a temporary BAM before being written as CRAM.
EOM
    )

    if [[ "${threads}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
Usage:
  convert_alignments_bed_awk [-h|--hlp|--help] threads fil_in fil_out ref_fa log_out log_err

Description:
  Convert one queryname-sorted paired-end BAM or CRAM file to BED.GZ using samtools, AWK, sort, and gzip.

Positional arguments:
  1  threads  <int>
    Number of threads to pass to samtools.

  2  fil_in  <file>
    Input queryname-sorted BAM or CRAM file.

  3  fil_out  <file>
    Output BED.GZ file.

  4  ref_fa  <file>
    Reference FASTA. Required for CRAM input.

  5  log_out  <file>
    File to receive stdout from the processing command.

  6  log_err  <file>
    File to receive stderr from the processing command.

Dependencies:
  External programs:
    - samtools
    - awk
    - sort
    - gzip

Returns:
  Returns 0 if the output BED.GZ file is written successfully; otherwise 1.

Notes:
  - This helper assumes adjacent queryname-sorted paired-end records.
  - Use the Python converter helper for single-end data.
EOM
    )

    if [[ "${threads}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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
Usage:
  convert_alignments_bed_python [-h|--hlp|--help] pth_scr_py fil_in fil_out ref_fa log_out log_err

Description:
  Convert one BAM or CRAM file to BED.GZ using a Python converter script.

Positional arguments:
  1  pth_scr_py  <file>
    Python converter script.

  2  fil_in  <file>
    Input BAM or CRAM file.

  3  fil_out  <file>
    Output BED.GZ file.

  4  ref_fa  <file>
    Reference FASTA. Required for CRAM input.

  5  log_out  <file>
    File to receive stdout from the processing command.

  6  log_err  <file>
    File to receive stderr from the processing command.

Dependencies:
  External programs:
    - python

  Python scripts:
    - pth_scr_py

Returns:
  Returns 0 if the output BED.GZ file is written successfully; otherwise 1.

Notes:
  - CRAM input is passed to the Python converter with '--ref_fa'.
EOM
    )

    if [[ "${pth_scr_py}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    validate_var_file "pth_scr_py" "${pth_scr_py}" || return 1
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

    python "${pth_scr_py}" \
        -i "${fil_in}" \
        -o "${fil_out}" \
        "${arr_ref_arg[@]}" \
        > "${log_out}" \
        2> "${log_err}"
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
