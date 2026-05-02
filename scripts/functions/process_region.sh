#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: process_region.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# check_region
# check_region_bam
# check_region_bdg


#  Source required helper functions if needed
# shellcheck disable=SC1091
{
    _dir_src_region="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    source "${_dir_src_region}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_region}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_region}" \
        check_inputs check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_region
}


#TODO: audit current usage; keep this helper even if unused
#MAYBE: make function "private"?
function check_region() {
    local region="${1:-}"
    local pat_rmn="I|II|III|IV|V|VI|VII|VIII|IX|X|XI|XII|XIII|XIV|XV|XVI"
    local pat_int="[0-9]+"
    local pat_ro pat_pi pat_cr pat_ci pat_ra
    local chr start end
    local show_help

    show_help=$(cat << EOM
Usage:
  check_region [-h|--hlp|--help] region

Description:
  Validate a genomic region string.

Positional argument:
  1  region  <str>  Region string to validate.

Returns:
  0 if 'region' is valid; otherwise 1.

Accepted formats:
  - Roman-numeral chromosome name: 'I' ... 'XVI'
  - Positive-integer chromosome name: '1', '2', ...
  - 'chr' followed by a Roman numeral: 'chrI' ... 'chrXVI'
  - 'chr' followed by a positive integer: 'chr1', 'chr2', ...
  - Genomic range:

      <chromosome>:<start>-<end>

Notes:
  - '<start>' and '<end>' must be positive integers.
  - If a range is supplied, '<end>' must be greater than or equal to '<start>'.

Examples:
  '''bash
  check_region chrI
  check_region XII:100-500
  check_region chr2:1-1000
  '''
EOM
    )

    if [[ "${region}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${region}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'region', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    pat_ro="^(${pat_rmn})$"
    pat_pi="^${pat_int}$"
    pat_cr="^chr(${pat_rmn})$"
    pat_ci="^chr${pat_int}$"
    pat_ra="^(${pat_rmn}|chr(${pat_rmn})|${pat_int}|chr${pat_int}):(${pat_int})-(${pat_int})$"

    if [[ "${region}" =~ ${pat_ro} ]]; then
        return 0
    elif [[ "${region}" =~ ${pat_pi} ]]; then
        return 0
    elif [[ "${region}" =~ ${pat_cr} ]]; then
        return 0
    elif [[ "${region}" =~ ${pat_ci} ]]; then
        return 0
    elif [[ "${region}" =~ ${pat_ra} ]]; then
        chr="${BASH_REMATCH[1]}"
        start="${BASH_REMATCH[3]}"
        end="${BASH_REMATCH[4]}"

        if (( end < start )); then
            echo_err_func "${FUNCNAME[0]}" \
                "range end must be greater than or equal to range start for" \
                "chromosome '${chr}' in region '${region}'."
            return 1
        fi

        return 0
    else
        echo_err_func "${FUNCNAME[0]}" \
            "invalid region format '${region}'. Expected a chromosome name" \
            "or range of the form '<chromosome>:<start>-<end>'."
        return 1
    fi
}


function check_region_bam() {
    #TODO: handle (S|CR)AM too
    #TODO: handle bedGraph as well, or create a separate function for that
    local bam="${1:-}"
    local region="${2:-}"
    local chr start end
    local ref_chr ref_siz
    local siz_chr=0
    local i
    local -a arr_chr arr_siz
    local show_help

    show_help=$(cat << EOM
Usage:
  check_region_bam [-h|--hlp|--help] bam region

Description:
  Check that a genomic region is valid and lies within the chromosome bounds encoded in a BAM file.

Positional arguments:
  1  bam     <str>  Input BAM file.
  2  region  <str>  Region string to validate against BAM chromosome bounds.

Returns:
  0 if 'region' is valid and in bounds for 'bam'; otherwise 1.

Dependency:
  - Samtools

Notes:
  - Chromosome sizes are derived from 'samtools idxstats'.
  - If 'region' is a chromosome name only, the full chromosome extent is assumed.
EOM
    )

    if [[ "${bam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${bam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'bam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${region}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'region', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var_file "bam" "${bam}" || return 1

    if ! command -v samtools > /dev/null 2>&1; then
        echo_err_func "${FUNCNAME[0]}" \
            "'samtools' is not installed or not in PATH."
        return 1
    fi

    check_region "${region}" || return 1

    while IFS=$'\t' read -r ref_chr ref_siz _; do
        if [[ "${ref_chr}" == "*" ]]; then
            continue
        fi

        arr_chr+=( "${ref_chr}" )
        arr_siz+=( "${ref_siz}" )
    done < <(samtools idxstats "${bam}")

    if (( ${#arr_chr[@]} == 0 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to retrieve chromosome sizes from BAM '${bam}'."
        return 1
    fi

    if [[ "${region}" =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]]; then
        chr="${BASH_REMATCH[1]}"
        start="${BASH_REMATCH[2]}"
        end="${BASH_REMATCH[3]}"
    else
        chr="${region}"
        start=1
        end=0
    fi

    for i in "${!arr_chr[@]}"; do
        if [[ "${arr_chr[i]}" == "${chr}" ]]; then
            siz_chr="${arr_siz[i]}"
            break
        fi
    done

    if (( siz_chr == 0 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "chromosome '${chr}' was not found in BAM '${bam}'."
        return 1
    fi

    if (( end == 0 )); then
        end="${siz_chr}"
    fi

    if (( start < 1 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "region start must be at least 1: '${region}'."
        return 1
    fi

    if (( end < start )); then
        echo_err_func "${FUNCNAME[0]}" \
            "region end must be greater than or equal to region start:" \
            "'${region}'."
        return 1
    fi

    if (( start > siz_chr || end > siz_chr )); then
        echo_err_func "${FUNCNAME[0]}" \
            "region '${region}' is out of bounds for chromosome '${chr}' in" \
            "BAM '${bam}' (chromosome size: ${siz_chr})."
        return 1
    fi

    return 0
}


#TODO: audit current usage; keep this helper even if unused
function check_region_bdg() {
    local fil_bdg="${1:-}"
    local region="${2:-}"
    local chr start end
    local ref_chr ref_end
    local siz_chr=0
    local i
    local cmd_cat
    local -a arr_chr arr_siz
    local show_help

    show_help=$(cat << EOM
Usage:
  check_region_bdg [-h|--hlp|--help] fil_bdg region

Description:
  Check that a genomic region is valid and lies within the chromosome bounds implied by a bedGraph file.

Positional arguments:
  1  fil_bdg  <str>  Input bedGraph file; may be plain text or '.gz'-compressed.
  2  region   <str>  Region string to validate against bedGraph chromosome bounds.

Returns:
  0 if 'region' is valid and in bounds for 'fil_bdg'; otherwise 1.

Notes:
  - Chromosome bounds are inferred from the maximum end coordinate observed for each chromosome in the bedGraph file.
  - If 'region' is a chromosome name only, the full observed chromosome extent is assumed.
  - This function assumes standard bedGraph coordinates: chromosome, start, end, value.
EOM
    )

    if [[ "${fil_bdg}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${fil_bdg}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'fil_bdg', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${region}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'region', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var_file "fil_bdg" "${fil_bdg}" || return 1
    check_region "${region}" || return 1

    if [[ "${fil_bdg}" =~ \.gz$ ]]; then
        if command -v gzip > /dev/null 2>&1; then
            cmd_cat=( gzip -cd "${fil_bdg}" )
        elif command -v zcat > /dev/null 2>&1; then
            cmd_cat=( zcat "${fil_bdg}" )
        else
            echo_err_func "${FUNCNAME[0]}" \
                "cannot read compressed bedGraph '${fil_bdg}': neither" \
                "'gzip' nor 'zcat' is available."
            return 1
        fi
    else
        cmd_cat=( cat "${fil_bdg}" )
    fi

    while IFS=$'\t' read -r ref_chr _ ref_end _; do
        [[ -n "${ref_chr}" ]] || continue
        [[ "${ref_chr}" == track* ]] && continue
        [[ "${ref_chr}" == browser* ]] && continue
        [[ "${ref_end}" =~ ^[0-9]+$ ]] || continue

        siz_chr=0
        for i in "${!arr_chr[@]}"; do
            if [[ "${arr_chr[i]}" == "${ref_chr}" ]]; then
                siz_chr="${arr_siz[i]}"
                break
            fi
        done

        if (( siz_chr == 0 )); then
            arr_chr+=( "${ref_chr}" )
            arr_siz+=( "${ref_end}" )
        elif (( ref_end > siz_chr )); then
            arr_siz[i]="${ref_end}"
        fi
    done < <("${cmd_cat[@]}")

    if (( ${#arr_chr[@]} == 0 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to retrieve chromosome bounds from bedGraph '${fil_bdg}'."
        return 1
    fi

    if [[ "${region}" =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]]; then
        chr="${BASH_REMATCH[1]}"
        start="${BASH_REMATCH[2]}"
        end="${BASH_REMATCH[3]}"
    else
        chr="${region}"
        start=1
        end=0
    fi

    siz_chr=0
    for i in "${!arr_chr[@]}"; do
        if [[ "${arr_chr[i]}" == "${chr}" ]]; then
            siz_chr="${arr_siz[i]}"
            break
        fi
    done

    if (( siz_chr == 0 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "chromosome '${chr}' was not found in bedGraph '${fil_bdg}'."
        return 1
    fi

    if (( end == 0 )); then
        end="${siz_chr}"
    fi

    if (( start < 1 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "region start must be at least 1: '${region}'."
        return 1
    fi

    if (( end < start )); then
        echo_err_func "${FUNCNAME[0]}" \
            "region end must be greater than or equal to region start:" \
            "'${region}'."
        return 1
    fi

    if (( start > siz_chr || end > siz_chr )); then
        echo_err_func "${FUNCNAME[0]}" \
            "region '${region}' is out of bounds for chromosome '${chr}' in" \
            "bedGraph '${fil_bdg}' (maximum observed end: ${siz_chr})."
        return 1
    fi

    return 0
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
