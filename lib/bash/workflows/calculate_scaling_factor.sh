#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: calculate_scaling_factor.sh
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in development and documentation.
#
# Distributed under the MIT license.


# _get_len_idx
# _get_dep_idx
# _set_ref_arg_cram
# _detect_typ_aln
# _resolve_typ_fil
# _get_expr_filter
# _count_alignments
# _calculate_frag_avg
# _compute_scl_fct
# _import_shell_asgmt
# _parse_metadata
# _calculate_dep_fct   #LEGACY
# _calculate_dep_arr   #LEGACY
# _compute_dep_all     #LEGACY
# _generate_fmt_str    #LEGACY
# _get_fil_out_part
# process_samp_siq
# process_samp_spike


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

#  Source required helper functions
_dir_src_sf="$(
    cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
)"

# shellcheck source=lib/bash/core/source_helpers.sh
source "${_dir_src_sf}/../core/source_helpers.sh" || {
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to source '${_dir_src_sf}/../core/source_helpers.sh'." >&2

    if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
        return 1
    else
        exit 1
    fi
}

source_helpers "${_dir_src_sf}/.." \
    check_args check_inputs check_source format_outputs run_python || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source required helper dependencies." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

unset _dir_src_sf


function _get_len_idx() {
    local key="${1:-}"
    local idx="${2:-}"
    local val=""
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _get_len_idx
    [--help] key idx

  Return the per-sample fragment-length override for a given key and sample index.

  If the corresponding override array contains exactly one value, that value is broadcast to all sample indices.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  key : str
    Key for which to retrieve the override; must be 'mip' or 'min'.

  2  idx : int
    Zero-based sample index.

Returns
-------
  Prints the matching override value if present and valid.

  If no usable override is available for the requested key/index combination, prints an empty line and returns 0.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Broadcast one main-IP fragment-length override to sample index 3.
    '''bash
    arr_len_mip=( 150 )
    _get_len_idx mip 3
    '''

  2. Confirm that an unsupported override key is rejected.
    '''bash
    if ! _get_len_idx unknown 0; then
        printf '%s\n' 'unsupported key rejected as expected'
    fi
    '''
EOM
    )

    if [[ "${key}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${key}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'key', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${idx}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'idx', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! [[ "${idx}" =~ ^[0-9]+$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'idx', must be a non-negative integer:" \
            "'${idx}'."
        return 1
    fi

    case "${key}" in
        mip)
            val="${arr_len_mip[idx]:-}"
            if [[ -z "${val}" ]] && (( ${#arr_len_mip[@]} == 1 )); then
                val="${arr_len_mip[0]}"
            fi
        ;;
        min)
            val="${arr_len_min[idx]:-}"
            if [[ -z "${val}" ]] && (( ${#arr_len_min[@]} == 1 )); then
                val="${arr_len_min[0]}"
            fi
        ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "positional argument 1, 'key', is '${key}' but must be" \
                "'mip' or 'min'."
            echo >&2
            echo "${show_help}" >&2
            return 1
        ;;
    esac

    if [[ "${val}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        echo "${val}"
        return 0
    fi

    echo
}


function _get_dep_idx() {
    local key="${1:-}"
    local idx="${2:-}"
    local val=""
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _get_dep_idx
    [--help] key idx

  Return the per-sample alignment-depth override for a given key and sample index.

  If the corresponding override array contains exactly one value, that value is broadcast to all sample indices.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  key : str
    Key for which to retrieve the override; must be 'mip', 'min', 'sip', or 'sin'.

  2  idx : int
    Zero-based sample index.

Returns
-------
  Prints the matching override value if present and valid.

  If no usable override is available for the requested key/index combination, prints an empty line and returns 0.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Select the third spike-input depth override.
    '''bash
    arr_dep_sin=( 120000 135000 142000 )
    _get_dep_idx sin 2
    '''

  2. Confirm that a nonnumeric sample index is rejected.
    '''bash
    if ! _get_dep_idx mip sample_one; then
        printf '%s\n' 'nonnumeric index rejected as expected'
    fi
    '''
EOM
    )

    if [[ "${key}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${key}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'key', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${idx}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'idx', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! [[ "${idx}" =~ ^[0-9]+$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'idx', must be a non-negative integer:" \
            "'${idx}'."
        return 1
    fi

    case "${key}" in
        mip)
            val="${arr_dep_mip[idx]:-}"
            if [[ -z "${val}" ]] && (( ${#arr_dep_mip[@]} == 1 )); then
                val="${arr_dep_mip[0]}"
            fi
        ;;

        min)
            val="${arr_dep_min[idx]:-}"
            if [[ -z "${val}" ]] && (( ${#arr_dep_min[@]} == 1 )); then
                val="${arr_dep_min[0]}"
            fi
        ;;

        sip)
            val="${arr_dep_sip[idx]:-}"
            if [[ -z "${val}" ]] && (( ${#arr_dep_sip[@]} == 1 )); then
                val="${arr_dep_sip[0]}"
            fi
        ;;

        sin)
            val="${arr_dep_sin[idx]:-}"
            if [[ -z "${val}" ]] && (( ${#arr_dep_sin[@]} == 1 )); then
                val="${arr_dep_sin[0]}"
            fi
        ;;

        *)
            echo_err_func "${FUNCNAME[0]}" \
                "positional argument 1, 'key', is '${key}' but must be" \
                "'mip', 'min', 'sip', or 'sin'."
            echo >&2
            echo "${show_help}" >&2
            return 1
        ;;
    esac

    if [[ "${val}" =~ ^[1-9][0-9]*$ ]]; then
        echo "${val}"
        return 0
    fi

    echo
}


#  Set Samtools reference arguments for one BAM or CRAM input
function _set_ref_arg_cram() {
    local fil_aln="${1:-}"
    local arr_ref_nam="${2:-}"
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _set_ref_arg_cram
    [--help] fil_aln arr_ref_nam

  Populate a named array with Samtools reference arguments for CRAM input.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fil_aln : file
    Input BAM or CRAM alignment file.

  2  arr_ref_nam : str
    Name of the caller-owned array to populate.

Returns
-------
  Returns 0 for non-CRAM input or after populating the named array; returns 1 if CRAM input requires 'ref_fa' and validation fails.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - Reference FASTA and required index (when processing CRAM)

Examples
--------
  1. Leave the reference-argument array empty for BAM input.
    '''bash
    ref_args=()
    _set_ref_arg_cram sample.bam ref_args
    declare -p ref_args
    '''

  2. Confirm that CRAM input without global 'ref_fa' is rejected.
    '''bash
    unset ref_fa
    ref_args=()
    if ! _set_ref_arg_cram sample.cram ref_args; then
        printf '%s\n' 'missing CRAM reference rejected as expected'
    fi
    '''
EOM
    )

    if [[ "${fil_aln}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    local -n arr_ref="${arr_ref_nam}"

    arr_ref=()

    if [[ "${fil_aln,,}" != *.cram ]]; then
        return 0
    elif [[ -z "${ref_fa:-}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "global variable 'ref_fa' is required for CRAM input:" \
            "'${fil_aln}'."
        return 1
    fi

    validate_var_file "ref_fa" "${ref_fa}" || return 1
    # shellcheck disable=SC2034
    arr_ref=( -T "${ref_fa}" )
}


function _detect_typ_aln() {
    local fil_aln="${1:-}"
    local n="${2:-200000}"
    local flag
    local seen=false
    local -a arr_ref_arg=()
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _detect_typ_aln
    [--help] fil_aln [n]

  Detect whether a BAM or CRAM alignment file appears to contain paired-end ("pe") or single-end ("se") alignments by sampling up to 'n' FLAG values from the file.

  If any sampled alignment has bit 0x1 set, the file is treated as paired-end; otherwise it is treated as single-end.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fil_aln : file
    Input BAM or CRAM alignment file.

  2  n : int
    Maximum number of alignments to sample (default: 200000).

Returns
-------
  Prints 'pe' or 'se'.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - cut
    - head
    - Reference FASTA and required index (when processing CRAM)
    - samtools

Examples
--------
  1. Detect the library layout from the default BAM sample size.
    '''bash
    _detect_typ_aln sample.bam
    '''

  2. Detect a CRAM layout from at most 50000 alignments.
    '''bash
    ref_fa=reference.fa
    _detect_typ_aln sample.cram 50000
    '''
EOM
    )

    if [[ "${fil_aln}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${fil_aln}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'fil_aln', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var_file "fil_aln" "${fil_aln}" || return 1

    if ! command -v samtools > /dev/null 2>&1; then
        echo_err_func "${FUNCNAME[0]}" \
            "'samtools' is not installed or not in PATH."
        return 1
    fi

    if ! [[ "${n}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'n', must be a positive integer: '${n}'."
        return 1
    fi

    _set_ref_arg_cram "${fil_aln}" arr_ref_arg || return 1

    while IFS= read -r flag; do
        [[ "${flag}" =~ ^[0-9]+$ ]] || continue
        seen=true

        if (( flag & 1 )); then echo "pe" && return 0; fi
    done < <(
        samtools view "${arr_ref_arg[@]}" "${fil_aln}" 2>/dev/null \
            | head -n "${n}" \
            | cut -f 2
    )

    if [[ "${seen}" != "true" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to read usable FLAG values from alignment file" \
            "'${fil_aln}'."
        return 1
    fi

    echo "se"
}


function _resolve_typ_fil() {
    local fil_aln="${1:-}"
    local pref="${2:-${aln_typ:-auto}}"
    local typ
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _resolve_typ_fil
    [--help] fil_aln [pref]

  Resolve the desired library end type for a BAM or CRAM file. If 'pref' is 'pe'/'paired' or 'se'/'single', that choice is returned directly. If 'pref' is 'auto' (or empty), the function calls '_detect_typ_aln' to infer the type from the alignment file.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fil_aln : file
    Input BAM or CRAM alignment file.

  2  pref : str
    Preferred library type; must be 'pe', 'paired', 'se', 'single', 'auto', or empty (default: \${aln_typ:-auto}).

Returns
-------
  Prints 'pe' or 'se'.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - cut (when 'pref' is 'auto' or empty)
    - head (when 'pref' is 'auto' or empty)
    - Reference FASTA and required index (when 'pref' is 'auto' or empty and processing CRAM)
    - samtools (when 'pref' is 'auto' or empty)

Examples
--------
  1. Use an explicit paired-end preference without sampling the BAM.
    '''bash
    _resolve_typ_fil sample.bam pe
    '''

  2. Auto-detect the layout of a CRAM using global 'ref_fa'.
    '''bash
    ref_fa=reference.fa
    _resolve_typ_fil sample.cram
    '''
EOM
    )

    if [[ "${fil_aln}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${fil_aln}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'fil_aln', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var_file "fil_aln" "${fil_aln}" || return 1

    case "${pref}" in
        pe|paired) echo "pe" ;;
        se|single) echo "se" ;;
        auto|"")
            typ="$(_detect_typ_aln "${fil_aln}")" || {
                echo_err_func "${FUNCNAME[0]}" \
                    "failed to auto-detect library type for alignment file" \
                    "'${fil_aln}'."
                return 1
            }

            case "${typ}" in
                pe|se) echo "${typ}" ;;
                *)
                    echo_err_func "${FUNCNAME[0]}" \
                        "unexpected detected library type '${typ}' for" \
                        "alignment file" \
                        "'${fil_aln}'."
                    return 1
                ;;
            esac
        ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "positional argument 2, 'pref', must be 'pe', 'paired'," \
                "'se', 'single', 'auto', or empty: '${pref}'."
            return 1
        ;;
    esac
}


function _get_expr_filter() {
    local aln_typ="${1:-pe}"  # Alignment type for file
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _get_expr_filter
    [--help] [aln_typ]

  Return the 'samtools view --expr' filtering expression corresponding to the requested alignment type.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  aln_typ : {'paired', 'pe', 'single', 'se'}
    Alignment layout type for input alignment files: 'paired', 'pe', 'single', or 'se' (default: 'pe').

Returns
-------
  Prints the filtering expression for the requested alignment type.

Notes
-----
  Runtime requirements:
    bash >= 4.4

  - PE data filtering expression: (flag == 99) || (flag == 1123) || (flag == 163) || (flag == 1187)
  - SE data filtering expression: (flag == 0) || (flag == 1024) || (flag == 16) || (flag == 1040)

Examples
--------
  1. Print the default paired-end filtering expression.
    '''bash
    _get_expr_filter
    '''

  2. Print the single-end filtering expression.
    '''bash
    _get_expr_filter se
    '''
EOM
    )

    if [[ "${aln_typ}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    fi

    case "${aln_typ}" in
        paired|pe) \
            echo \
                "(flag == 99) || (flag == 1123) || (flag == 163) ||" \
                "(flag == 1187)"
        ;;

        single|se) \
            echo \
                "(flag == 0) || (flag == 1024) || (flag == 16) ||" \
                "(flag == 1040)"
        ;;

        *)
            echo_err_func "${FUNCNAME[0]}" \
                "positional argument 1, 'aln_typ', must be 'paired', 'pe'," \
                "'single', or 'se': '${aln_typ}'."
            return 1
        ;;
    esac
}


function _count_alignments() {
    local threads="${1:-}"    # No. threads for parallelization
    local fil_aln="${2:-}"    # Input BAM or CRAM alignment file
    local aln_typ="${3:-pe}"  # "paired", "pe", "single", "se" (default: "pe")
    local expr                # Samtools filtration expression
    local -a arr_ref_arg=()   # Samtools CRAM reference arguments
    local show_help           # Help message/documentation

    show_help=$(cat << EOM
Usage
-----
  _count_alignments
    [--help] threads fil_aln [aln_typ]

  Count alignments in a BAM or CRAM file based on whether the data is made up of paired- ("paired" or "pe") or single-end ("single" or "se") sequenced read alignments. Uses 'samtools view' with filtering expressions to count specific alignment flags.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  threads : int
    Number of threads to use for 'samtools view' decompression.

  2  fil_aln : file
    Input BAM or CRAM alignment file for which to count alignments.

  3  aln_typ : {'paired', 'pe', 'single', 'se'}
    Alignment layout type for input alignment files: 'paired', 'pe', 'single', or 'se' (default: '${aln_typ}').

Returns
-------
  An integer representing the count of alignments matching the given type.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - Reference FASTA and required index (when processing CRAM)
    - samtools

Examples
--------
  1. Count paired-end alignments in a BAM using the default layout.
    '''bash
    _count_alignments 8 sample.bam
    '''

  2. Count single-end alignments in a CRAM using its reference FASTA.
    '''bash
    ref_fa=reference.fa
    _count_alignments 4 sample.cram single
    '''
EOM
    )

    if [[ "${threads}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${threads}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'threads', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${fil_aln}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'fil_aln', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! [[ "${threads}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'threads', must be a positive integer:" \
            "'${threads}'."
        return 1
    fi

    validate_var_file "fil_aln" "${fil_aln}" || return 1

    if ! command -v samtools > /dev/null 2>&1; then
        echo_err_func "${FUNCNAME[0]}" \
            "'samtools' is not installed or not in PATH."
        return 1
    fi

    case "${aln_typ}" in
        single|se|paired|pe) : ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "positional argument 3, 'aln_typ', must be 'paired', 'pe'," \
                "'single', or 'se': '${aln_typ}'."
            return 1
        ;;
    esac

    #  Determine filtering flags based on alignment type
    expr="$(_get_expr_filter "${aln_typ}")" || return 1
    _set_ref_arg_cram "${fil_aln}" arr_ref_arg || return 1

    #  Count alignments based on alignment type
    samtools view \
        -@ "${threads}" \
        "${arr_ref_arg[@]}" \
        -c \
        --expr "${expr}" \
        "${fil_aln}" \
        || {
            echo_err_func "${FUNCNAME[0]}" \
                "'samtools view' failed for '${fil_aln}' with type" \
                "'${aln_typ}' and expression '${expr}'."
            return 1
        }
}


function _calculate_frag_avg() {
    local threads="${1:-}"    # No. threads for parallelization
    local fil_aln="${2:-}"    # Input BAM or CRAM alignment file
    local aln_typ="${3:-pe}"  # "paired", "pe", "single", "se"
    local len_lcl="${4:-}"    # Optional default length for SE libraries
    local expr=""             # Samtools filtration expression
    local -a arr_ref_arg=()   # Samtools CRAM reference arguments
    local show_help           # Help message/documentation

    show_help=$(cat << EOM
Usage
-----
  _calculate_frag_avg
    [--help] threads fil_aln [aln_typ] [len_lcl]

  Computes the average fragment length from a BAM or CRAM file based on whether the data is paired-end ("paired" or "pe") or single-end ("single" or "se"). Uses 'samtools view' with filtering expressions and 'awk' to process fragment lengths.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  threads : int
    Number of threads to use for 'samtools view' decompression.

  2  fil_aln : file
    Input BAM or CRAM alignment file for which to compute fragment lengths.

  3  aln_typ : {'paired', 'pe', 'single', 'se'}
    Alignment layout type for input alignment files: 'paired', 'pe', 'single', or 'se' (default: 'pe').

  4  len_lcl : int
    Default fragment length to use for single-end libraries when TLEN is not meaningful.

Expected globals
----------------
  len_def : int
    Default fragment length used for single-end libraries when 'len_lcl' is not supplied.

Returns
-------
  A floating-point value representing the average fragment length.

Notes
-----
  Runtime requirements:
    - awk (when deriving paired-end fragment length)
    - bash >= 4.4
    - Reference FASTA and required index (when deriving fragment length from CRAM)
    - samtools

Examples
--------
  1. Compute the mean positive TLEN for paired-end BAM alignments.
    '''bash
    _calculate_frag_avg 8 sample.bam paired
    '''

  2. Use a fixed fragment length for a single-end library.
    '''bash
    _calculate_frag_avg 4 sample.bam single 150
    '''
EOM
    )

    if [[ "${threads}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${threads}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'threads', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${fil_aln}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'fil_aln', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! [[ "${threads}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'threads', must be a positive integer:" \
            "'${threads}'."
        return 1
    fi

    validate_var_file "fil_aln" "${fil_aln}" || return 1

    if ! command -v samtools > /dev/null 2>&1; then
        echo_err_func "${FUNCNAME[0]}" \
            "'samtools' is not installed or not in PATH."
        return 1
    fi

    case "${aln_typ}" in
        single|se) aln_typ="single" ;;
        paired|pe) aln_typ="paired" ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "positional argument 3, 'aln_typ', must be 'paired', 'pe'," \
                "'single', or 'se': '${aln_typ}'."
            return 1
        ;;
    esac

    #  If single-end, TLEN is not meaningful, so use provided default (per SAM
    #+ spec, TLEN is 0 for single-end reads)
    if [[ "${aln_typ}" == "single" ]]; then
        if [[ -n "${len_lcl}" ]]; then
            if ! [[ "${len_lcl}" =~ ^[1-9][0-9]*([.][0-9]+)?$ ]]; then
                echo_err_func "${FUNCNAME[0]}" \
                    "positional argument 4, 'len_lcl', must be a positive" \
                    "number: '${len_lcl}'."
                return 1
            fi
            echo "${len_lcl}"
            return 0
        elif [[ -n "${len_def:-}" ]]; then
            if ! [[ "${len_def}" =~ ^[1-9][0-9]*([.][0-9]+)?$ ]]; then
                echo_err_func "${FUNCNAME[0]}" \
                    "global variable 'len_def' must be a positive number:" \
                    "'${len_def}'."
                return 1
            fi
            echo "${len_def}"
            return 0
        else
            echo_err_func "${FUNCNAME[0]}" \
                "single-end library detected but no default fragment length" \
                "provided. Pass a 4th positional argument ('len_lcl') or set" \
                "global variable 'len_def'."
            return 1
        fi
    fi

    #  Determine filtering flags based on alignment type
    expr="$(_get_expr_filter "${aln_typ}")" || return 1
    _set_ref_arg_cram "${fil_aln}" arr_ref_arg || return 1

    #  Compute average fragment length using samtools and awk
    samtools view -@ "${threads}" "${arr_ref_arg[@]}" --expr "${expr}" "${fil_aln}" \
        | awk '{
            if ($9 > 0) { sum += $9; count++ }
        } END {
            if (count > 0) { print sum / count }
            else {
                print \
                    "error(_calculate_frag_avg): "
                        "no valid fragment lengths found." > "/dev/stderr"
                exit 1
            }
        }' || {
            echo_err_func "${FUNCNAME[0]}" \
            "failed for alignment file '${fil_aln}' with type '${aln_typ}'."
            return 1
        }
}


function _compute_scl_fct() {
    local mode="${1:-}"     # Workflow mode. Coefficient family: "siq" or "spike"
    local scr_siq="${2:-}"  # Entry point for siQ-ChIP scaling factor
    local scr_spk="${3:-}"  # Entry point for spike-in scaling factor
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _compute_scl_fct
    [--help] mode scr_siq scr_spk [args...]

  Compute a scaling coefficient by dispatching to the appropriate Python entry point.

  If 'mode' is 'siq', the function runs 'scr_siq'. If 'mode' is 'spike', the function runs 'scr_spk'. Any additional arguments are passed through unchanged to the selected Python script.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  mode : {'siq', 'spike'}
    Workflow mode. Coefficient family; must be 'siq' or 'spike'.

  2  scr_siq : structured string
    Python entry point for siQ-ChIP coefficient calculation.

  3  scr_spk : structured string
    Python entry point for spike-in coefficient calculation.

  4+ args : structured string
    Additional arguments passed to the selected Python entry point.

Returns
-------
  Prints the output produced by the selected Python script.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - python >= 3.11

  - In console mode, 'scr_siq' and 'scr_spk' are installed command names.
  - In module mode, they identify CLIs in the bundled package source tree.

Examples
--------
  1. Dispatch an equation-6 calculation to the siQ-ChIP entry point.
    '''bash
    _compute_scl_fct \\
        siq \\
        calculate_scaling_factor_siqchip \\
        calculate_scaling_factor_spike \\
        --eqn 6 \\
        --mass_ip 5 --mass_in 25 --vol_all 50 --vol_in 10 \\
        --dep_ip 120000 --dep_in 110000 --len_ip 160 --len_in 150
    '''

  2. Dispatch a fractional-coefficient calculation to the spike-in entry point.
    '''bash
    _compute_scl_fct \\
        spike \\
        calculate_scaling_factor_siqchip \\
        calculate_scaling_factor_spike \\
        --coef fractional \\
        --main_ip 120000 --spike_ip 30000 \\
        --main_in 110000 --spike_in 25000
    '''
EOM
    )

    if [[ "${mode}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${mode}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'mode', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${scr_siq}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'scr_siq', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${scr_spk}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'scr_spk', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    shift 3
    local params=( "$@" )

    case "${mode}" in
        siq)   run_py "${scr_siq}" "${params[@]}" ;;
        spike) run_py "${scr_spk}" "${params[@]}" ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "positional argument 1, 'mode', is '${mode}' but must be" \
                "either 'siq' or 'spike'."
            return 1
            ;;
    esac
}


function _import_shell_asgmt() {
    local func="${1:-}"
    local arr_ref_nam="${2:-}"
    shift 2

    local -A arr_allowed=()
    local line nam rhs stmt
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _import_shell_asgmt
    [--help] func arr_ref allowed_name [allowed_name ...]

  Import validated shell assignments from an array of Python-emitted lines.

  Each non-empty input line must have the form:

      export name=value

  Only variable names explicitly listed in the allowed-name argument list are accepted.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Calling-function name for error reporting.

  2  arr_ref : str
    Name of array variable containing shell-assignment lines.

  3+ allowed_name : str
    One or more allowed variable names.

Returns
-------
  Assigns validated values into the current shell and returns 0 on success.

Notes
-----
  Runtime requirements:
    bash >= 4.4

  - This helper is intentionally restrictive.
  - It rejects unexpected variable names and selected unsafe shell constructs.

Examples
--------
  1. Import two explicitly allowed metadata assignments.
    '''bash
    arr_shell=( 'export mass_ip=5' 'export mass_in=25' )
    _import_shell_asgmt _parse_metadata arr_shell mass_ip mass_in
    '''

  2. Confirm that an unexpected assignment name is rejected.
    '''bash
    arr_shell=( 'export unexpected=5' )
    if ! _import_shell_asgmt _parse_metadata arr_shell mass_ip; then
        printf '%s\n' 'unexpected assignment rejected as expected'
    fi
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${func}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'func', is missing."
        return 1
    elif [[ -z "${arr_ref_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'arr_ref', is missing."
        return 1
    elif [[ ! "${arr_ref_nam}" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'arr_ref', is not a valid shell" \
            "identifier: '${arr_ref_nam}'."
        return 1
    elif [[ "$#" -eq 0 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "at least one allowed variable name must be supplied."
        return 1
    fi

    local -n arr_lines_ref="${arr_ref_nam}"

    for nam in "$@"; do
        if [[ ! "${nam}" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "allowed variable name is not a valid shell identifier:" \
                "'${nam}'."
            return 1
        fi
        arr_allowed["${nam}"]=1
    done

    for line in "${arr_lines_ref[@]}"; do
        [[ -z "${line}" ]] && continue

        if [[
            ! "${line}" =~ ^export[[:space:]]+([A-Za-z_][A-Za-z0-9_]*)=(.*)$
        ]]; then
            echo_err_func "${func}" \
                "metadata parser emitted unsupported shell assignment:" \
                "'${line}'."
            return 1
        fi

        nam="${BASH_REMATCH[1]}"
        rhs="${BASH_REMATCH[2]}"

        if [[ -z "${arr_allowed[${nam}]+x}" ]]; then
            echo_err_func "${func}" \
                "metadata parser attempted to assign unexpected variable" \
                "'${nam}'."
            return 1
        fi

        # shellcheck disable=SC2016
        if [[
            "${rhs}" == *'$('* || "${rhs}" == *'`'* || "${rhs}" == *';'*
        ]]; then
            echo_err_func "${func}" \
                "metadata parser emitted unsafe shell value for variable" \
                "'${nam}': '${rhs}'."
            return 1
        fi

        stmt="${nam}=${rhs}"

        if ! \
            eval "${stmt}"
        then
            echo_err_func "${func}" \
                "failed to import parsed metadata assignment for variable" \
                "'${nam}'."
            return 1
        fi
    done
}


function _parse_metadata() {
    local scr_met="${1:-}"  # Python entry point for parsing metadata
    local fil_aln="${2:-}"  # Alignment file to process
    local tbl_met="${3:-}"  # siQ-ChIP metadata table
    local cfg_met="${4:-}"  # YAML configuration for metadata parsing
    local -a arr_shell      # Python-emitted shell assignment lines
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _parse_metadata
    [--help] scr_met fil_aln tbl_met cfg_met

  Use the siQ-ChIP metadata Python parser to extract metadata values for an alignment file and assign them in the current shell.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  scr_met : structured string
    Python entry point for parsing siQ-ChIP metadata.

  2  fil_aln : file
    Input BAM or CRAM alignment file for which metadata should be retrieved.

  3  tbl_met : file
    siQ-ChIP metadata table.

  4  cfg_met : file
    YAML configuration file for metadata parsing.

Returns
-------
  Assigns parsed metadata variables in the current shell environment from validated 'export key=value' lines emitted by the Python helper.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - python >= 3.11

Examples
--------
  1. Import metadata for one BAM through the configured table schema.
    '''bash
    _parse_metadata \\
        parse_metadata_siqchip \\
        sample.bam \\
        metadata.tsv \\
        metadata.yaml
    '''

  2. Confirm that omitting the metadata configuration is rejected.
    '''bash
    if ! _parse_metadata \
        parse_metadata_siqchip sample.bam metadata.tsv
    then
        printf '%s\n' 'missing configuration rejected as expected'
    fi
    '''
EOM
    )

    if [[ "${scr_met}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${scr_met}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'scr_met', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${fil_aln}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'fil_aln', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${tbl_met}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'tbl_met', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${cfg_met}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 4, 'cfg_met', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var_file "fil_aln" "${fil_aln}" || return 1
    validate_var_file "tbl_met" "${tbl_met}" || return 1
    validate_var_file "cfg_met" "${cfg_met}" || return 1

    # shellcheck disable=SC2034
    if ! \
        mapfile -t arr_shell < <(
            run_py "${scr_met}" \
                --verbose \
                --alignment "${fil_aln}" \
                --tbl_met "${tbl_met}" \
                --cfg "${cfg_met}" \
                --shell
        )
    then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to use script '${scr_met}' with configuration" \
            "'${cfg_met}' to parse siQ-ChIP metadata in '${tbl_met}' for" \
            "file '${fil_aln}'."
        return 1
    fi

    _import_shell_asgmt \
        "${FUNCNAME[0]}" \
        arr_shell \
        vol_in vol_all mass_in mass_ip \
        conc_in conc_ip len_in len_ip lib_vol_in lib_vol_ip dep_in dep_ip \
        || return 1
}


function _get_fil_out_part() {
    local fil_out="${1:-}"
    local idx="${2:-}"
    local idx_pad
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _get_fil_out_part
    [--help] fil_out idx

  Construct the per-sample part-file path derived from a base output-table path.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fil_out : file
    Output file path. Base output-table path.

  2  idx : int
    Zero-based, zero-padded sample index.

Returns
-------
  Prints a path of the form:

      <fil_out>.part.<idx>

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Derive the first zero-padded part-file name.
    '''bash
    _get_fil_out_part results.tsv 0
    '''

  2. Confirm that a nonnumeric part index is rejected.
    '''bash
    if ! _get_fil_out_part results.tsv sample_one; then
        printf '%s\n' 'nonnumeric index rejected as expected'
    fi
    '''
EOM
    )

    if [[ "${fil_out}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${fil_out}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'fil_out', is missing."
        return 1
    elif [[ -z "${idx}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'idx', is missing."
        return 1
    elif [[ ! "${idx}" =~ ^[0-9]+$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'idx', must be a non-negative integer:" \
            "'${idx}'."
        return 1
    fi

    printf -v idx_pad '%06d' "${idx}"
    printf '%s.part.%s\n' "${fil_out}" "${idx_pad}"
}


# ------------------------------ Begin #LEGACY ------------------------------ #
#  The following global variable and four functions are retained for reference
#+ and, e.g., future comparisons with Python helpers; they are no longer called
#+ by the production scaling-factor row writers
DEP_BINS_DFLT="${DEP_BINS_DFLT:-1,5,10,20,30,40,50}"


function _calculate_dep_fct() {
    local n_in="${1:-}"             # Alignment count for input (not IP) BAM file
    local siz_bin="${2:-10}"        # Bin size (in bp)
    local siz_gen="${3:-12157105}"  # Effective genome size for model organism (in bp)
    local mode="${4:-norm}"         # "frag" or "norm"
    local dp="${5:-24}"            # Number of decimal points for rounding
    local fct_dep                   # Variable for calculations
    local show_help                 # Help message/documentation

    show_help=$(cat << EOM
Usage
-----
  _calculate_dep_fct
    [--help] n_in [siz_bin] [siz_gen] [mode] [dp]

  Computes a minimum input depth factor from an input-alignment count, bin size, and effective genome size.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  n_in : int
    Input-alignment count.

  2  siz_bin : int
    Bin size in base pairs (default: ${siz_bin}).

  3  siz_gen : int
    Effective genome size (in base pairs; default: ${siz_gen} [appropriate for S. cerevisiae]).

  4  mode : {'frag', 'norm'}
    Workflow mode. Mode of calculation; options: "frag" or "norm" (default: '${mode}').

  5  dp : int
    Maximum number of decimal places retained for finite emitted values. Number of decimal places for rounding (default: ${dp}).

Returns
-------
  The computed depth factor.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - bc

Examples
--------
  1. Compute normalized-coverage depth with default precision.
    '''bash
    _calculate_dep_fct 12851824 30 12157105 norm
    '''

  2. Compute fragment-adjusted depth with 12 decimal places.
    '''bash
    _calculate_dep_fct 12851824 20 12157105 frag 12
    '''
EOM
    )

    if [[ "${n_in}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${n_in}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'n_in', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! command -v bc > /dev/null 2>&1; then
        echo_err_func "${FUNCNAME[0]}" \
            "'bc' is not installed or not in PATH."
        return 1
    fi

    for var in n_in siz_bin siz_gen dp; do
        if ! [[ "${!var}" =~ ^[1-9][0-9]*$ ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "'${var}' must be a positive integer: '${!var}'."
            return 1
        fi
    done

    case "${mode}" in
        frag|norm) : ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "positional argument 4, 'mode', must be 'frag' or 'norm':" \
                "'${mode}'."
            return 1
        ;;
    esac

    #  Compute depth factor via 'bc' operation based on 'mode'
    if [[ "${mode}" == "norm" ]]; then
        #  For "normalized coverage"
        fct_dep=$(bc -l <<< "
            scale=${dp};
            (${siz_bin}) / (${siz_gen} * (1 - (${siz_bin} / ${siz_gen})))
        ")
    else
        #  For fragment-length-adjusted signal
        fct_dep=$(bc -l <<< "
            scale=${dp};
            (${n_in} * ${siz_bin}) / (${siz_gen} * (1 - (${siz_bin} / ${siz_gen})))
        ")
    fi

    #  Add leading zero if bc returns .ddd or -.ddd
    if [[ "${fct_dep}" =~ ^\.[0-9] ]]; then
        fct_dep="0${fct_dep}"
    elif [[ "${fct_dep}" =~ ^-\.[0-9] ]]; then
        fct_dep="-0${fct_dep:1}"
    fi

    echo "${fct_dep}"
}


function _calculate_dep_arr() {
    local dep="${1:-}"                      # Number of mapped reads in sample
    local mod="${2:-norm}"                  # Data transformation mode
    local egs="${3:-12157105}"              # Effective genome size
    local dp="${4:-24}"                    # Rounding precision for output
    local csv_bin="${5:-${DEP_BINS_DFLT}}"  # Comma-delimited bin sizes
    local bin
    local -a arr_dep arr_bin
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _calculate_dep_arr
    [--help] dep [mod] [egs] [dp] [csv_bin]

  Compute a comma-delimited series of minimum input depth values across one or more bin sizes for either fragment-length-adjusted signal ('frag') or "normalized coverage" ('norm'), following Dickson et al., Sci Rep 2023.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  dep : int
    Number of mapped reads/alignments in the sample.

  2  mod : str
    Data transformation mode; must be 'frag' or 'norm' (default: ${mod}).

  3  egs : int
    Effective genome size (default: ${egs} [appropriate for S. cerevisiae]).

  4  dp : int
    Maximum number of decimal places retained for finite emitted values. Number of decimal places to round to (default: ${dp}).

  5  csv_bin : str
    Comma-separated list of bin sizes to use (default: ${csv_bin}).

Returns
-------
  Prints a comma-delimited list of depth values, one per bin size, in the order given by 'csv_bin'.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - bc

Examples
--------
  1. Compute normalized depths for the default bin series.
    '''bash
    _calculate_dep_arr 12851824
    '''

  2. Compute fragment-adjusted depths for a custom bin series.
    '''bash
    _calculate_dep_arr 12851824 frag 12157105 24 1,10,25,50
    '''
EOM
    )

    if [[ "${dep}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${dep}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'dep', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${mod}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'mod', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${egs}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'egs', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${dp}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 4, 'dp', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${csv_bin}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 5, 'csv_bin', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    for var in dep egs dp; do
        if ! [[ "${!var}" =~ ^[1-9][0-9]*$ ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "'${var}' is '${!var}' but must be a positive integer."
            return 1
        fi
    done

    IFS=',' read -r -a arr_bin <<< "${csv_bin}"

    if [[ ${#arr_bin[@]} -eq 0 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "no usable bin sizes were parsed from positional argument 5," \
            "'csv_bin': '${csv_bin}'."
        return 1
    fi

    for bin in "${arr_bin[@]}"; do
        if ! [[ "${bin}" =~ ^[1-9][0-9]*$ ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "bin size '${bin}' in positional argument 5, 'csv_bin', must" \
                "be a positive integer."
            return 1
        fi
    done

    case "${mod}" in
        frag|norm) : ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "positional argument 2, 'mod', is '${mod}' but must be" \
                "'frag' or 'norm'."
            return 1
            ;;
    esac

    for bin in "${arr_bin[@]}"; do
        arr_dep+=( "$(
            _calculate_dep_fct "${dep}" "${bin}" "${egs}" "${mod}" "${dp}"
        )" )
    done

    IFS=',' printf "%s\n" "${arr_dep[*]}"
}


function _compute_dep_all() {
    local dep="${1:-}"                      # Number of alignments in sample BAM
    local dp="${2:-24}"                    # Number of decimals to round to
    local csv_bin="${3:-${DEP_BINS_DFLT}}"  # Comma-delimited bin sizes
    local -a arr_dm_fr arr_dm_nm
    local output
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _compute_dep_all
    [--help] dep [dp] [csv_bin]

  Compute minimum input depth values needed downstream for both fragment-length-adjusted signal ('frag') and "normalized coverage" ('norm'), following Dickson et al., Sci Rep 2023.

  Internally, this function calls '_calculate_dep_arr' twice and concatenates the two comma-delimited outputs.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  dep : int
    Number of mapped reads/alignments in the sample.

  2  dp : int
    Maximum number of decimal places retained for finite emitted values. Number of decimal places to round to (default: ${dp}).

  3  csv_bin : str
    Comma-separated list of bin sizes to use (default: ${csv_bin}).

Returns
-------
  Prints one comma-delimited string containing the 'frag' depth values followed by the 'norm' depth values, each in the order given by 'csv_bin'.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - bc

Examples
--------
  1. Compute both depth families using the default precision and bins.
    '''bash
    _compute_dep_all 12851824
    '''

  2. Compute both depth families for a custom bin series.
    '''bash
    _compute_dep_all 12851824 12 1,10,25,50
    '''
EOM
    )

    if [[ "${dep}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${dep}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'dep', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${dp}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'dp', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${csv_bin}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'csv_bin', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! [[ "${dep}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'dep', must be a positive integer:" \
            "'${dep}'."
        return 1
    fi

    if ! [[ "${dp}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'dp', must be a positive integer:" \
            "'${dp}'."
        return 1
    fi

    #  Note: detailed value validation is delegated to '_calculate_dep_arr'
    IFS=',' read -r -a arr_dm_fr < <(
        _calculate_dep_arr "${dep}" "frag" "12157105" "${dp}" "${csv_bin}"
    ) || return 1
    IFS=',' read -r -a arr_dm_nm < <(
        _calculate_dep_arr "${dep}" "norm" "12157105" "${dp}" "${csv_bin}"
    ) || return 1

    output=$(IFS=','; echo "${arr_dm_fr[*]},${arr_dm_nm[*]}")
    echo "${output}"
}


function _generate_fmt_str() {
    local num_fld="${1:-}"  # Number of fields in the output row
    local fmt_str=""        # Variable for format string
    local i
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _generate_fmt_str
    [--help] num_fld

  Construct a tab-delimited 'printf' format string containing 'num_fld' '%s' fields and no trailing tab.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  num_fld : int
    Number of fields in the output row.

Returns
-------
  Prints a format string suitable for tab-delimited row output.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Generate a format string for three tab-separated fields.
    '''bash
    _generate_fmt_str 3
    '''

  2. Confirm that a zero-field format is rejected.
    '''bash
    if ! _generate_fmt_str 0; then
        printf '%s\n' 'zero-field format rejected as expected'
    fi
    '''
EOM
    )

    if [[ "${num_fld}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${num_fld}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'num_fld', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! [[ "${num_fld}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'num_fld', must be a positive integer:" \
            "'${num_fld}'."
        return 1
    fi

    for ((i = 1; i <= num_fld; i++)); do fmt_str+="%s\t"; done
    printf "%s\n" "${fmt_str%$'\t'}"
}
# ------------------------------- End #LEGACY ------------------------------- #


#  Compute siQ-ChIP alpha scaling factor and related values for a sample
#+
#+ Workflow function that processes a sample using global variables; extracts
#+ siQ-ChIP metadata from a TSV table, computes alignment counts, and computes
#+ average fragment lengths
#+
# shellcheck disable=SC2154
function process_samp_siq() {
    local idx="${1:-}"  # Array sample index

    #  Declare local variables
    local fil_ip fil_in siq mass_ip mass_in vol_all vol_in
    local lib_vol_ip lib_vol_in
    local typ_ip typ_in dep_ip dep_in len_ip len_in v
    local dep_ip_met dep_in_met len_ip_met len_in_met
    # local fmt_str  # Reserved in case formatted output generation is revived
    local len_ip_ovrd len_in_ovrd fil_out_part
    local -a arr_fields arr_arg_siq
    local show_help

    show_help=$(cat << EOM
Usage
-----
  process_samp_siq
    [--help] idx

  Processes a single sample (array index 'idx') to compute the siQ-ChIP alpha scaling factor and related QC values.

  For the IP/input alignment-file pair at 'idx', this function
    - parses siQ-ChIP metadata,
    - detects PE/SE for each alignment file,
    - counts alignments,
    - computes average fragment lengths (PE only; SE uses a provided default),
    - writes a per-sample tab-separated results row.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  idx : int
    Zero-based array index into arr_mip/arr_min.

Expected globals
----------------
  arr_mip, arr_min : array
    Sample alignment-file arrays for IP and input, respectively.

  scr_met, scr_siq, scr_spk : structured string
    Metadata-parser, siQ-ChIP, and spike-in Python entry points, respectively.

  cfg_met, tbl_met, ref_fa, fil_out : file
    Metadata-configuration, metadata-table, reference-FASTA, and output-table paths, respectively.

  eqn, aln_typ : str
    siQ-ChIP equation and optional alignment-type override, respectively.

  threads, len_def, dp, idx_out : int
    Thread count, default fragment length, rounding precision, and optional output-part index, respectively.

  debug : bool
    If 'true', prints 'debug_var' lines.

Returns
-------
  Returns 0 after writing the sample row; otherwise 1.

Output
------
  Prints these row fields in order:
    fil_ip, fil_in, siq, eqn, mass_ip, mass_in, vol_all, vol_in, dep_ip, dep_in, len_ip, len_in

    If metadata supplies a library-volume correction, appends lib_vol_ip and lib_vol_in.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - awk (when deriving paired-end fragment length)
    - bash >= 4.4
    - cut (when auto-detecting alignment layout)
    - dirname
    - head (when auto-detecting alignment layout)
    - python >= 3.11
    - Reference FASTA and required index (when processing CRAM)
    - samtools (when alignment depth, fragment length, or layout must be derived)

  - For single-end libraries, TLEN is not meaningful; this function requires 'len_def' (or an optional 4th positional argument, 'len_lcl') to supply a fixed fragment length.
  - End-type detection is per-file in case of mixed inputs.
  - Individual alignment files are expected to be a single input type: either entirely PE or entirely SE.

Examples
--------
  1. Process the first fixture-backed siQ-ChIP IP/input pair after initializing wrapper globals.
    '''bash
    fix=tests/fixtures/calculate_scaling_factor
    arr_mip=( "\${fix}/bam/pe/IP_WT_G1_Hho1_6336.sc.bam" )
    arr_min=( "\${fix}/bam/pe/in_WT_G1_Hho1_6336.sc.bam" )
    arr_len_mip=()
    arr_len_min=()
    arr_dep_mip=()
    arr_dep_min=()
    scr_met=parse_metadata_siqchip
    cfg_met="\${fix}/config/parse_metadata_siqchip_field_to_column.yml"
    tbl_met="\${fix}/metadata/measurements_siqchip.tsv"
    ref_fa="\${fix}/reference/tiny.fa"
    scr_siq=calculate_scaling_factor_siqchip
    scr_spk=calculate_scaling_factor_spike
    fil_out="\${TMPDIR:-/tmp}/scaling.siq.tsv"
    eqn=6 aln_typ=pe threads=1 len_def=150 dp=12 idx_out=0 debug=false
    process_samp_siq 0
    '''

  2. Confirm that a nonnumeric sample index is rejected before processing.
    '''bash
    if ! process_samp_siq sample_one; then
        printf '%s\n' 'nonnumeric index rejected as expected'
    fi
    '''
EOM
    )

    if [[ "${idx}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${idx}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'idx', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! [[ "${idx}" =~ ^[0-9]+$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'idx', must be a non-negative integer:" \
            "'${idx}'."
        return 1
    fi

    #  Assign alignment files based on sample index
    fil_ip="${arr_mip[idx]}"
    fil_in="${arr_min[idx]}"

    #  Check that alignment files exist
    validate_var_file "fil_ip" "${fil_ip}" "${idx}" || return 1
    validate_var_file "fil_in" "${fil_in}" "${idx}" || return 1

    if [[ "${debug}" == "true" ]]; then
        debug_var "idx=${idx}" "fil_ip=${fil_ip}" "fil_in=${fil_in}"
    fi

    #  Parse siQ-ChIP metadata, assigning global variables
    _parse_metadata \
        "${scr_met}" "${fil_ip}" "${tbl_met}" "${cfg_met}" || {
            echo_err_func "${FUNCNAME[0]}" \
                "failed while running '_parse_metadata' for IP alignment" \
                "file '${fil_ip}'."
            return 1
        }

    dep_ip_met="${dep_ip:-}"
    dep_in_met="${dep_in:-}"
    len_ip_met="${len_ip:-}"
    len_in_met="${len_in:-}"
    lib_vol_ip="${lib_vol_ip:-NA}"
    lib_vol_in="${lib_vol_in:-NA}"

    #  Determine end type per alignment file (robust if inputs ever mix)
    typ_ip="$(_resolve_typ_fil "${fil_ip}")" || {
        echo_err_func "${FUNCNAME[0]}" \
            "failed while resolving 'typ_ip'."
        return 1
    }
    typ_in="$(_resolve_typ_fil "${fil_in}")" || {
        echo_err_func "${FUNCNAME[0]}" \
            "failed while resolving 'typ_in'."
        return 1
    }

    #  Count alignments in alignment files
    dep_ip="$(_get_dep_idx mip "${idx}")"
    if [[ -z "${dep_ip}" ]]; then
        if [[ -n "${dep_ip_met}" && "${dep_ip_met}" != "NA" ]]; then
            dep_ip="${dep_ip_met}"
        else
            dep_ip="$(_count_alignments "${threads}" "${fil_ip}" "${typ_ip}")" || {
                echo_err_func "${FUNCNAME[0]}" \
                    "failed while counting alignments for IP file '${fil_ip}'" \
                    "with type '${typ_ip}'."
                return 1
            }
        fi
    fi

    dep_in="$(_get_dep_idx min "${idx}")"
    if [[ -z "${dep_in}" ]]; then
        if [[ -n "${dep_in_met}" && "${dep_in_met}" != "NA" ]]; then
            dep_in="${dep_in_met}"
        else
            dep_in="$(_count_alignments "${threads}" "${fil_in}" "${typ_in}")" || {
                echo_err_func "${FUNCNAME[0]}" \
                    "failed while counting alignments for input file '${fil_in}'" \
                    "with type '${typ_in}'."
                return 1
            }
        fi
    fi

    #  Compute average fragment lengths for alignment files; overrides take
    #+ precedence over metadata, TLEN, or 'len_def'
    len_ip_ovrd="$(_get_len_idx mip "${idx}")"
    len_in_ovrd="$(_get_len_idx min "${idx}")"

    if [[ -n "${len_ip_ovrd}" ]]; then
        len_ip="${len_ip_ovrd}"
    elif [[ -n "${len_ip_met}" && "${len_ip_met}" != "NA" ]]; then
        len_ip="${len_ip_met}"
    else
        len_ip="$(
            _calculate_frag_avg \
                "${threads}" "${fil_ip}" "${typ_ip}" "${len_def:-}"
        )" || {
            echo_err_func "${FUNCNAME[0]}" \
                "failed while computing average fragment length for IP file" \
                "'${fil_ip}' with type '${typ_ip}'."
            return 1
        }
    fi

    if [[ -n "${len_in_ovrd}" ]]; then
        len_in="${len_in_ovrd}"
    elif [[ -n "${len_in_met}" && "${len_in_met}" != "NA" ]]; then
        len_in="${len_in_met}"
    else
        len_in="$(
            _calculate_frag_avg \
                "${threads}" "${fil_in}" "${typ_in}" "${len_def:-}"
        )" || {
            echo_err_func "${FUNCNAME[0]}" \
                "failed while computing average fragment length for input" \
                "file '${fil_in}' with type '${typ_in}'."
            return 1
        }
    fi

    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "eqn=${eqn}"         "dp=${dp}" \
            "mass_ip=${mass_ip}" "mass_in=${mass_in}" \
            "vol_all=${vol_all}" "vol_in=${vol_in}" \
            "typ_ip=${typ_ip}"   "typ_in=${typ_in}" \
            "dep_ip=${dep_ip}"   "dep_in=${dep_in}" \
            "len_ip=${len_ip}"   "len_in=${len_in}" \
            "lib_vol_ip=${lib_vol_ip}" "lib_vol_in=${lib_vol_in}"
    fi

    #  Compute siQ-ChIP alpha scaling factor
    arr_arg_siq=(
        --eqn     "${eqn}"     --dp      "${dp}"
        --mass_ip "${mass_ip}" --mass_in "${mass_in}"
        --vol_all "${vol_all}" --vol_in  "${vol_in}"
        --dep_ip  "${dep_ip}"  --dep_in  "${dep_in}"
        --len_ip  "${len_ip}"  --len_in  "${len_in}"
    )

    if [[ "${lib_vol_ip}" != "NA" && "${lib_vol_in}" != "NA" ]]; then
        arr_arg_siq+=(
            --lib_vol_ip "${lib_vol_ip}"
            --lib_vol_in "${lib_vol_in}"
        )
    fi

    siq=$(
        _compute_scl_fct \
            "siq" "${scr_siq}" "${scr_spk}" \
            "${arr_arg_siq[@]}"
    ) || return 1

    #  Build a row of results, printing them tab-separated with no trailing tab
    if [[ -z "${fil_out:-}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "global variable 'fil_out' is empty or unset."
        return 1
    fi

    if ! [[ -d "$(dirname "${fil_out}")" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "parent directory of 'fil_out' does not exist:" \
            "'$(dirname "${fil_out}")'."
        return 1
    fi

    fil_out_part="$(
        _get_fil_out_part "${fil_out}" "${idx_out:-${idx}}"
    )" || return 1

    if ! {
            arr_fields=(
                "${fil_ip}" "${fil_in}" "${siq}" "${eqn}"
                "${mass_ip}" "${mass_in}" "${vol_all}" "${vol_in}"
                "${dep_ip}" "${dep_in}" "${len_ip}" "${len_in}"
            )

            if [[ "${lib_vol_ip}" != "NA" && "${lib_vol_in}" != "NA" ]]; then
                arr_fields+=( "${lib_vol_ip}" "${lib_vol_in}" )
            fi

            printf '%s' "${arr_fields[0]}"
            for v in "${arr_fields[@]:1}"; do printf '\t%s' "${v}"; done
            printf '\n'
        } > "${fil_out_part}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to write per-sample results file '${fil_out_part}'."
        return 1
    fi
    #  The execute-layer driver combines per-sample part files deterministically
    #+ after all workers finish successfully.
}


#  Compute spike-in scaling factor and related values for a sample
#+
#+ Workflow function that processes a sample using global variables; computes
#+ alignment counts and a requested scaling-factor coefficient
#+
# shellcheck disable=SC2154
function process_samp_spike() {
    local idx="${1:-}"  # Array sample index

    #  Declare local variables
    local mp mn sp sn typ_mp typ_sp typ_mn typ_sn
    local num_mp num_sp num_mn num_sn coef_lcl coef_val v fil_out_part
    # local fmt_str  # Reserved in case formatted output generation is revived
    local -a arr_fields
    local show_help

    show_help=$(cat << EOM
Usage
-----
  process_samp_spike
    [--help] idx

  Processes a single sample (array index 'idx') to compute a spike-in scaling factor and related QC values.

  For the main/spike IP and input alignment files at 'idx', this function detects paired- or single-end sequenced read alignment status per file, counts alignments, computes the requested scaling-factor coefficient, and writes a per-sample tab-separated results row.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  idx : int
    Zero-based array index into arrays: arr_mip (main IP), arr_sip (spike IP), arr_min (main input), arr_sin (spike input).

Expected globals
----------------
  arr_mip, arr_sip, arr_min, arr_sin : array
    Main-IP, spike-IP, main-input, and spike-input alignment arrays, respectively.

  scr_siq, scr_spk : structured string
    siQ-ChIP and spike-in Python entry points, respectively.

  ref_fa, fil_out : file
    Reference-FASTA and output-table paths, respectively.

  aln_typ, coef_spk : str
    Alignment-type override and requested spike-in coefficient, respectively.

  threads, dp, idx_out : int
    Thread count, rounding precision, and optional output-part index, respectively.

  debug : bool
    If 'true', prints 'debug_var' lines.

Returns
-------
  Returns 0 after writing the sample row; otherwise 1.

Output
------
  Prints these row fields in order:
    mp, sp, mn, sn, coef_val, coef_lcl, num_mp, num_sp, num_mn, num_sn

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - bash >= 4.4
    - cut (when auto-detecting alignment layout)
    - dirname
    - head (when auto-detecting alignment layout)
    - python >= 3.11
    - Reference FASTA and required index (when processing CRAM)
    - samtools (when alignment depth or layout must be derived)

  - End-type detection is per-file in case of mixed inputs.
  - Individual alignment files are expected to be a single input type: either entirely PE or entirely SE.

Examples
--------
  1. Process the first fixture-backed main/spike IP/input quartet using depth overrides.
    '''bash
    fix=tests/fixtures/calculate_scaling_factor/bam/pe
    arr_mip=( "\${fix}/IP_WT_G1_Hho1_6336.sc.bam" )
    arr_sip=( "\${fix}/IP_WT_G1_Hho1_6336.sp.bam" )
    arr_min=( "\${fix}/in_WT_G1_Hho1_6336.sc.bam" )
    arr_sin=( "\${fix}/in_WT_G1_Hho1_6336.sp.bam" )
    arr_dep_mip=( 120000 )
    arr_dep_sip=( 30000 )
    arr_dep_min=( 110000 )
    arr_dep_sin=( 25000 )
    ref_fa=tests/fixtures/calculate_scaling_factor/reference/tiny.fa
    scr_siq=calculate_scaling_factor_siqchip
    scr_spk=calculate_scaling_factor_spike
    fil_out="\${TMPDIR:-/tmp}/scaling.spike.tsv"
    aln_typ=pe coef_spk=chiprx_alpha_ratio threads=1 dp=12 idx_out=0 debug=false
    process_samp_spike 0
    '''

  2. Confirm that a nonnumeric sample index is rejected before processing.
    '''bash
    if ! process_samp_spike sample_one; then
        printf '%s\n' 'nonnumeric index rejected as expected'
    fi
    '''
EOM
    )

    if [[ "${idx}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${idx}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'idx', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! [[ "${idx}" =~ ^[0-9]+$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'idx', must be a non-negative integer:" \
            "'${idx}'."
        return 1
    fi

    #  Assign alignment files based on sample index
    mp="${arr_mip[idx]}"
    mn="${arr_min[idx]}"
    sp="${arr_sip[idx]}"
    sn="${arr_sin[idx]}"

    #  Check that alignment files exist
    validate_var_file "mp" "${mp}" "${idx}" || return 1
    validate_var_file "sp" "${sp}" "${idx}" || return 1
    validate_var_file "mn" "${mn}" "${idx}" || return 1
    validate_var_file "sn" "${sn}" "${idx}" || return 1

    if [[ "${debug}" == "true" ]]; then
        debug_var "idx=${idx}" "mp=${mp}" "mn=${mn}" "sp=${sp}" "sn=${sn}"
    fi

    #  Determine end type per alignment file
    typ_mp="$(_resolve_typ_fil "${mp}")" || {
        echo_err_func "${FUNCNAME[0]}" \
            "failed while resolving 'typ_mp'."
        return 1
    }
    typ_sp="$(_resolve_typ_fil "${sp}")" || {
        echo_err_func "${FUNCNAME[0]}" \
            "failed while resolving 'typ_sp'."
        return 1
    }
    typ_mn="$(_resolve_typ_fil "${mn}")" || {
        echo_err_func "${FUNCNAME[0]}" \
            "failed while resolving 'typ_mn'."
        return 1
    }
    typ_sn="$(_resolve_typ_fil "${sn}")" || {
        echo_err_func "${FUNCNAME[0]}" \
            "failed while resolving 'typ_sn'."
        return 1
    }

    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "typ_mp=${typ_mp}" "typ_sp=${typ_sp}" \
            "typ_mn=${typ_mn}" "typ_sn=${typ_sn}"
    fi

    #  Count alignments in alignment files
    num_mp="$(_get_dep_idx mip "${idx}")"
    if [[ -z "${num_mp}" ]]; then
        num_mp="$(_count_alignments "${threads}" "${mp}" "${typ_mp}")" || {
            echo_err_func "${FUNCNAME[0]}" \
                "failed while counting alignments for main IP file '${mp}'" \
                "with type '${typ_mp}'."
            return 1
        }
    fi

    num_mn="$(_get_dep_idx min "${idx}")"
    if [[ -z "${num_mn}" ]]; then
        num_mn="$(_count_alignments "${threads}" "${mn}" "${typ_mn}")" || {
            echo_err_func "${FUNCNAME[0]}" \
                "failed while counting alignments for main input file '${mn}'" \
                "with type '${typ_mn}'."
            return 1
        }
    fi

    num_sp="$(_get_dep_idx sip "${idx}")"
    if [[ -z "${num_sp}" ]]; then
        num_sp="$(_count_alignments "${threads}" "${sp}" "${typ_sp}")" || {
            echo_err_func "${FUNCNAME[0]}" \
                "failed while counting alignments for spike IP file '${sp}'" \
                "with type '${typ_sp}'."
            return 1
        }
    fi

    num_sn="$(_get_dep_idx sin "${idx}")"
    if [[ -z "${num_sn}" ]]; then
        num_sn="$(_count_alignments "${threads}" "${sn}" "${typ_sn}")" || {
            echo_err_func "${FUNCNAME[0]}" \
                "failed while counting alignments for spike input file" \
                "'${sn}' with type '${typ_sn}'."
            return 1
        }
    fi

    coef_lcl="${coef_spk:-chiprx_alpha_ratio}"

    if [[ "${coef_lcl}" =~ ^([Aa][Ll][Ll])$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'coef_spk=all' is not supported here because this function" \
            "expects a single coefficient value per sample row."
        return 1
    fi

    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "num_mp=${num_mp}" "num_sp=${num_sp}" \
            "num_mn=${num_mn}" "num_sn=${num_sn}" \
            "coef_spk=${coef_lcl}" "dp=${dp}"
    fi

    #  Compute requested spike-in coefficient
    coef_val="$(
        _compute_scl_fct \
            "spike" "${scr_siq}" "${scr_spk}" \
            --coef    "${coef_lcl}" \
            --format  "plain" \
            --main_ip "${num_mp}" --spike_ip "${num_sp}" \
            --main_in "${num_mn}" --spike_in "${num_sn}" \
            --dp      "${dp}"
    )" || {
        echo_err_func "${FUNCNAME[0]}" \
            "failed while computing spike-in coefficient '${coef_lcl}' for" \
            "sample idx '${idx}'."
        return 1
    }

    #  Build a row of results, printing them tab-separated with no trailing tab
    if [[ -z "${fil_out:-}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "global variable 'fil_out' is empty or unset."
        return 1
    fi

    if ! [[ -d "$(dirname "${fil_out}")" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "parent directory of 'fil_out' does not exist:" \
            "'$(dirname "${fil_out}")'."
        return 1
    fi

    fil_out_part="$(
        _get_fil_out_part "${fil_out}" "${idx_out:-${idx}}"
    )" || return 1

    if ! {
            arr_fields=(
                "${mp}" "${sp}" "${mn}" "${sn}" "${coef_val}" "${coef_lcl}"
                "${num_mp}" "${num_sp}" "${num_mn}" "${num_sn}"
            )

            printf '%s' "${arr_fields[0]}"
            for v in "${arr_fields[@]:1}"; do printf '\t%s' "${v}"; done
            printf '\n'
        } > "${fil_out_part}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to write per-sample results file '${fil_out_part}'."
        return 1
    fi
    #  The execute-layer driver combines per-sample part files deterministically
    #+ after all workers finish successfully.
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
