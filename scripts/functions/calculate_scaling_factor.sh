#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: calculate_scaling_factor.sh
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# _get_len_idx
# _get_dep_idx
# _detect_typ_bam
# _resolve_typ_fil
# _get_expr_filter
# _count_align_bam
# _calculate_frag_avg
# _compute_scl_fct
# _import_shell_asgmt
# _parse_metadata
# _calculate_dep_fct
# _calculate_dep_arr
# _compute_dep_all
# _generate_fmt_str
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

#  Set default bin sizes for minimum input depth calculations
DEP_BINS_DFLT="${DEP_BINS_DFLT:-1,5,10,20,30,40,50}"


#  Source required helper functions
_dir_src_sf="$(
    cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
)"

# shellcheck disable=SC1091
source "${_dir_src_sf}/source_helpers.sh" || {
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to source '${_dir_src_sf}/source_helpers.sh'." >&2

    if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
        return 1
    else
        exit 1
    fi
}

source_helpers "${_dir_src_sf}" \
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
Usage:
  _get_len_idx [-h|--hlp|--help] key idx

Description:
  Return the per-sample fragment-length override for a given key and sample index. If the corresponding override array contains exactly one value, that value is broadcast to all sample indices.

Positional arguments:
  1  key  <str>  Key for which to retrieve the override; must be 'mip' or 'min'.
  2  idx  <int>  Zero-based sample index.

Returns:
  Prints the matching override value if present and valid.

  If no usable override is available for the requested key/index combination, prints an empty line and returns 0.

Examples:
  '''bash
  _get_len_idx mip 0
  _get_len_idx min 3
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
Usage:
  _get_dep_idx [-h|--hlp|--help] key idx

Description:
  Return the per-sample alignment-depth override for a given key and sample index. If the corresponding override array contains exactly one value, that value is broadcast to all sample indices.

Positional arguments:
  1  key  <str>  Key for which to retrieve the override; must be 'mip', 'min', 'sip', or 'sin'.
  2  idx  <int>  Zero-based sample index.

Returns:
  Prints the matching override value if present and valid.

  If no usable override is available for the requested key/index combination, prints an empty line and returns 0.

Examples:
  '''bash
  _get_dep_idx mip 0
  _get_dep_idx sin 2
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


function _detect_typ_bam() {
    local bam="${1:-}"
    local n="${2:-200000}"
    local flag
    local seen=false
    local show_help

    show_help=$(cat << EOM
Usage:
  _detect_typ_bam [-h|--hlp|--help] bam [n]

Description:
  Detect whether a BAM file appears to contain paired-end ("pe") or single-end ("se") alignments by sampling up to 'n' FLAG values from the file.

  If any sampled alignment has bit 0x1 set, the file is treated as paired-end; otherwise it is treated as single-end.

Positional arguments:
  1  bam  <str>  Input BAM file.
  2  n    <int>  Maximum number of alignments to sample (default: 200000).

Returns:
  Prints 'pe' or 'se'.

Dependency:
  - Samtools

Examples:
  '''bash
  _detect_typ_bam sample.bam
  _detect_typ_bam sample.bam 50000
  '''
EOM
    )

    if [[ "${bam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${bam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'bam', is missing."
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

    if ! [[ "${n}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'n', must be a positive integer: '${n}'."
        return 1
    fi

    while IFS= read -r flag; do
        [[ "${flag}" =~ ^[0-9]+$ ]] || continue
        seen=true

        if (( flag & 1 )); then echo "pe" && return 0; fi
    done < <(
        samtools view "${bam}" 2>/dev/null \
            | head -n "${n}" \
            | cut -f 2
    )

    if [[ "${seen}" != "true" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to read usable FLAG values from BAM '${bam}'."
        return 1
    fi

    echo "se"
}


function _resolve_typ_fil() {
    local bam="${1:-}"
    local pref="${2:-${aln_typ:-auto}}"
    local typ
    local show_help

    show_help=$(cat << EOM
Usage:
  _resolve_typ_fil [-h|--hlp|--help] bam [pref]

Description:
  Resolve the desired library end type for a BAM file. If 'pref' is 'pe'/'paired' or 'se'/'single', that choice is returned directly. If 'pref' is 'auto' (or empty), the function calls '_detect_typ_bam' to infer the type from the BAM file.

Positional arguments:
  1  bam   <str>  Input BAM file.
  2  pref  <str>  Preferred library type; must be 'pe', 'paired', 'se', 'single', 'auto', or empty (default: \${aln_typ:-auto}).

Returns:
  Prints 'pe' or 'se'.

Examples:
  '''bash
  _resolve_typ_fil sample.bam
  _resolve_typ_fil sample.bam pe
  _resolve_typ_fil sample.bam auto
  '''
EOM
    )

    if [[ "${bam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${bam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'bam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var_file "bam" "${bam}" || return 1

    case "${pref}" in
        pe|paired) echo "pe" ;;
        se|single) echo "se" ;;
        auto|"")
            typ="$(_detect_typ_bam "${bam}")" || {
                echo_err_func "${FUNCNAME[0]}" \
                    "failed to auto-detect library type for BAM '${bam}'."
                return 1
            }

            case "${typ}" in
                pe|se) echo "${typ}" ;;
                *)
                    echo_err_func "${FUNCNAME[0]}" \
                        "unexpected detected library type '${typ}' for BAM" \
                        "'${bam}'."
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
Usage:
  _get_expr_filter [-h|--hlp|--help] [aln_typ]

Description:
  Return the 'samtools view --expr' filtering expression corresponding to the requested alignment type.

Positional argument:
  1  aln_typ  <str>  Alignment type; use 'paired', 'pe', 'single', or 'se' (default: 'pe').

Returns:
  Prints the filtering expression for the requested alignment type.

Notes:
  - PE data filtering expression: (flag == 99) || (flag == 1123) || (flag == 163) || (flag == 1187)
  - SE data filtering expression: (flag == 0) || (flag == 1024) || (flag == 16) || (flag == 1040)

Examples:
  '''bash
  _get_expr_filter pe
  _get_expr_filter paired
  _get_expr_filter se
  _get_expr_filter single
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


function _count_align_bam() {
    local threads="${1:-}"    # No. threads for parallelization
    local fil_bam="${2:-}"    # Input (not IP) BAM file
    local aln_typ="${3:-pe}"  # "paired", "pe", "single", "se" (default: "pe")
    local expr                # Samtools filtration expression
    local show_help           # Help message/documentation

    show_help=$(cat << EOM
Usage:
  _count_align_bam [-h|--hlp|--help] threads fil_bam [aln_typ]

Description:
  Counts the number of alignments in a BAM file based on whether the data is made up of paired- ("paired" or "pe") or single-end ("single" or "se") sequenced read alignments. Uses 'samtools view' with filtering expressions to count specific alignment flags.

Positional arguments:
  1  threads  <int>  Number of threads for 'samtools view' decompression.
  2  fil_bam  <str>  BAM infile for which to count alignments.
  3  aln_typ  <str>  Alignment type; options: 'paired', 'pe', 'single', or 'se' (default: ${aln_typ}).

Returns:
  An integer representing the count of alignments matching the given type.

Dependency:
  - Samtools

Examples:
  '''bash
  #  Count alignments in a BAM file of paired-end alignments using 8 threads
  _count_align_bam 8 sample.bam paired

  #  Count alignments in a BAM file of single-end alignments using 4 threads
  _count_align_bam 4 sample.bam single
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
    elif [[ -z "${fil_bam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'fil_bam', is missing."
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

    validate_var_file "fil_bam" "${fil_bam}" || return 1

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

    #  Count alignments based on alignment type
    samtools view -@ "${threads}" -c --expr "${expr}" "${fil_bam}" || {
        echo_err_func "${FUNCNAME[0]}" \
            "'samtools view' failed for '${fil_bam}' with type '${aln_typ}'" \
            "and expression '${expr}'."
        return 1
    }
}


function _calculate_frag_avg() {
    local threads="${1:-}"    # No. threads for parallelization
    local fil_bam="${2:-}"    # BAM infile
    local aln_typ="${3:-pe}"  # "paired", "pe", "single", "se"
    local len_lcl="${4:-}"    # Optional default length for SE libraries
    local expr=""             # Samtools filtration expression
    local show_help           # Help message/documentation

    show_help=$(cat << EOM
Usage:
  _calculate_frag_avg [-h|--hlp|--help] threads fil_bam [aln_typ] [len_lcl]

Description:
  Computes the average fragment length from a BAM file based on whether the data is paired-end ("paired" or "pe") or single-end ("single" or "se"). Uses 'samtools view' with filtering expressions and 'awk' to process fragment lengths.

Positional arguments:
  1  threads  <int>  Number of threads for 'samtools view' decompression.
  2  fil_bam  <str>  BAM infile for which to compute fragment lengths.
  3  aln_typ  <str>  Alignment type; options: 'paired', 'pe', 'single', or 'se' (default: 'pe').
  4  len_lcl  <int>  Default fragment length to use for single-end libraries when TLEN is not meaningful.

Optional global variable:
  Read:
    len_def  <int>  Default fragment length used for single-end libraries when 'len_lcl' is not supplied.

Returns:
  A floating-point value representing the average fragment length.

Dependency:
  - Samtools

Examples:
  '''bash
  #  Compute average fragment length for paired-end alignments using 8 threads
  _calculate_frag_avg 8 sample.bam paired

  #  Compute average fragment length for single-end alignments using 4 threads
  _calculate_frag_avg 4 sample.bam single
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
    elif [[ -z "${fil_bam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'fil_bam', is missing."
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

    validate_var_file "fil_bam" "${fil_bam}" || return 1

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

    #  Compute average fragment length using samtools and awk
    samtools view -@ "${threads}" --expr "${expr}" "${fil_bam}" \
        | awk '{
            if ($9 > 0) { sum += $9; count++ }
        } END {
            if (count > 0) { print sum / count }
            else {
                print \
                    "error(_calculate_frag_avg): no valid fragment " \
                    "lengths found." > "/dev/stderr"
                exit 1
            }
        }' || {
            echo_err_func "${FUNCNAME[0]}" \
                "failed for BAM '${fil_bam}' with type '${aln_typ}'."
            return 1
        }
}


function _compute_scl_fct() {
    local mode="${1:-}"     # Coefficient family: "siq" or "spike"
    local scr_siq="${2:-}"  # Script to compute siQ-ChIP scaling factor
    local scr_spk="${3:-}"  # Script to compute spike-in scaling factor
    local show_help

    show_help=$(cat << EOM
Usage:
  _compute_scl_fct [-h|--hlp|--help] mode scr_siq scr_spk [args...]

Description:
  Compute a scaling coefficient by dispatching to the appropriate Python entry point.

  If 'mode' is 'siq', the function runs 'scr_siq'. If 'mode' is 'spike', the function runs 'scr_spk'. Any additional arguments are passed through unchanged to the selected Python script.

Positional arguments:
  1  mode     <str>  Coefficient family; must be 'siq' or 'spike'.
  2  scr_siq  <str>  Python entry point for siQ-ChIP coefficient calculation.
  3  scr_spk  <str>  Python entry point for spike-in coefficient calculation.
  4+ args     <...>  Additional arguments passed to the selected Python entry point.

Returns:
  Prints the output produced by the selected Python script.

Note:
  - 'scr_siq' and 'scr_spk' must exist as readable Python entry-point files.

Examples (do not run):
  '''bash
  _compute_scl_fct
      siq
      path/to/calculate_scaling_factor_siq_chip.py
      path/to/calculate_scaling_factor_spike.py
      --eqn 6

  _compute_scl_fct
      spike
      path/to/calculate_scaling_factor_siq_chip.py
      path/to/calculate_scaling_factor_spike.py
      --coef fractional
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

    validate_var_file "scr_siq" "${scr_siq}" || return 1
    validate_var_file "scr_spk" "${scr_spk}" || return 1

    if [[ ! -r "${scr_siq}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'scr_siq' is not readable: '${scr_siq}'."
        return 1
    fi

    if [[ ! -r "${scr_spk}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'scr_spk' is not readable: '${scr_spk}'."
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
    local arr_ref="${2:-}"
    shift 2

    local -n arr_lines="${arr_ref}"
    local -A allowed=()
    local line nam rhs stmt
    local show_help

    show_help=$(cat << EOM
Usage:
  _import_shell_asgmt [-h|--hlp|--help] func arr_ref allowed_name [allowed_name ...]

Description:
  Import validated shell assignments from an array of Python-emitted lines.

  Each non-empty input line must have the form:

      export name=value

  Only variable names explicitly listed in the allowed-name argument list are accepted.

Positional arguments:
  1  func          <str>  Calling-function name for error reporting.
  2  arr_ref       <str>  Name of array variable containing shell-assignment lines.
  3+ allowed_name  <str>  One or more allowed variable names.

Returns:
  Assigns validated values into the current shell and returns 0 on success.

Notes:
  - This helper is intentionally restrictive.
  - It rejects unexpected variable names and selected unsafe shell constructs.

Example:
  '''bash
  _import_shell_asgmt "_parse_metadata" arr_shell mass_ip mass_in vol_all vol_in
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
    elif [[ -z "${arr_ref}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'arr_ref', is missing."
        return 1
    elif [[ "$#" -eq 0 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "at least one allowed variable name must be supplied."
        return 1
    fi

    for nam in "$@"; do
        if [[ ! "${nam}" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "allowed variable name is not a valid shell identifier:" \
                "'${nam}'."
            return 1
        fi
        allowed["${nam}"]=1
    done

    for line in "${arr_lines[@]}"; do
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

        if [[ -z "${allowed[${nam}]+x}" ]]; then
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

        stmt="declare -- ${nam}=${rhs}"

        if ! eval "${stmt}"; then
            echo_err_func "${func}" \
                "failed to import parsed metadata assignment for variable" \
                "'${nam}'."
            return 1
        fi
    done
}


function _parse_metadata() {
    local scr_met="${1:-}"  # Python script for parsing metadata
    local fil_bam="${2:-}"  # BAM file to process
    local tbl_met="${3:-}"  # siQ-ChIP metadata table
    local cfg_met="${4:-}"  # YAML configuration for metadata parsing
    local eqn="${5:-}"      # Equation included with parsing
    local -a arr_shell      # Python-emitted shell assignment lines
    local show_help

    show_help=$(cat << EOM
Usage:
  _parse_metadata [-h|--hlp|--help] scr_met fil_bam tbl_met cfg_met eqn

Description:
  Use the siQ-ChIP metadata Python parser to extract metadata values for a BAM file and assign them in the current shell.

Positional arguments:
  1  scr_met  <str>  Python entry point for parsing siQ-ChIP metadata.
  2  fil_bam  <str>  BAM file for which metadata should be retrieved.
  3  tbl_met  <str>  siQ-ChIP metadata table.
  4  cfg_met  <str>  YAML configuration file for metadata parsing.
  5  eqn      <str>  Equation identifier included with parsed metadata.

Returns:
  Assigns parsed metadata variables in the current shell environment from validated 'export key=value' lines emitted by the Python helper.

Example:
  '''bash
  _parse_metadata "\${scr_met}" "\${fil_bam}" "\${tbl_met}" "\${cfg_met}" "\${eqn}"
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
    elif [[ -z "${fil_bam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'fil_bam', is missing."
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
    elif [[ -z "${eqn}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 5, 'eqn', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var_file "scr_met" "${scr_met}" || return 1
    validate_var_file "fil_bam" "${fil_bam}" || return 1
    validate_var_file "tbl_met" "${tbl_met}" || return 1
    validate_var_file "cfg_met" "${cfg_met}" || return 1

    # if ! \
    #     eval "$(
    #         run_py "${scr_met}" \
    #             --verbose \
    #             --bam "${fil_bam}" \
    #             --tbl_met "${tbl_met}" \
    #             --cfg "${cfg_met}" \
    #             --eqn "${eqn}" \
    #             --shell
    #     )"
    # then
    #     echo_err_func "${FUNCNAME[0]}" \
    #         "failed to use script '${scr_met}' with configuration" \
    #         "'${cfg_met}' to parse siQ-ChIP metadata in '${tbl_met}' for" \
    #         "file '${fil_bam}'."
    #     return 1
    # fi

    # shellcheck disable=SC2034
    if ! mapfile -t arr_shell < <(
        run_py "${scr_met}" \
            --verbose \
            --bam "${fil_bam}" \
            --tbl_met "${tbl_met}" \
            --cfg "${cfg_met}" \
            --eqn "${eqn}" \
            --shell
    ); then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to use script '${scr_met}' with configuration" \
            "'${cfg_met}' to parse siQ-ChIP metadata in '${tbl_met}' for" \
            "file '${fil_bam}'."
        return 1
    fi

    _import_shell_asgmt \
        "${FUNCNAME[0]}" \
        arr_shell \
        eqn \
        vol_in vol_all mass_in mass_ip \
        conc_in conc_ip len_in len_ip dep_in dep_ip \
        || return 1
}


function _calculate_dep_fct() {
    local n_in="${1:-}"             # Alignment count for input (not IP) BAM file
    local siz_bin="${2:-10}"        # Bin size (in bp)
    local siz_gen="${3:-12157105}"  # Effective genome size for model organism (in bp)
    local mode="${4:-norm}"         # "frag" or "norm"
    local rnd="${5:-24}"            # Number of decimal points for rounding
    local fct_dep                   # Variable for calculations
    local show_help                 # Help message/documentation

    show_help=$(cat << EOM
Usage:
  _calculate_dep_fct [-h|--hlp|--help] n_in [siz_bin] [siz_gen] [mode] [rnd]

Description:
  Computes a minimum input depth factor from an input-alignment count, bin size, and effective genome size.

Positional arguments:
  1  n_in     <int>  Input-alignment count.
  2  siz_bin  <int>  Bin size (in base pairs; default: ${siz_bin}).
  3  siz_gen  <int>  Effective genome size (in base pairs; default: ${siz_gen} [appropriate for S. cerevisiae]).
  4  mode     <str>  Mode of calculation; options: "frag" or "norm" (default: ${mode}).
  5  rnd      <int>  Number of decimal places for rounding (default: ${rnd}).

Returns:
  The computed depth factor.

Dependency:
  - bc

Examples:
  '''bash
  #  Compute depth factor for fragment-length-adjusted signal
  _calculate_dep_fct 12851824 20 12157105 "frag" 12

  #  Compute depth factor for "normalized coverage"
  _calculate_dep_fct 12851824 30 12157105 "norm"
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

    for var in n_in siz_bin siz_gen rnd; do
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
            scale=${rnd};
            (${siz_bin}) / (${siz_gen} * (1 - (${siz_bin} / ${siz_gen})))
        ")
    else
        #  For fragment-length-adjusted signal
        fct_dep=$(bc -l <<< "
            scale=${rnd};
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
    local rnd="${4:-24}"                    # Rounding precision for output
    local csv_bin="${5:-${DEP_BINS_DFLT}}"  # Comma-delimited bin sizes
    local bin
    local -a arr_dep arr_bin
    local show_help

    show_help=$(cat << EOM
Usage:
  _calculate_dep_arr dep [mod] [egs] [rnd] [csv_bin]

Description:
  Compute a comma-delimited series of minimum input depth values across one or more bin sizes for either fragment-length-adjusted signal ('frag') or "normalized coverage" ('norm'), following Dickson et al., Sci Rep 2023.

Positional arguments:
  1  dep      <int>  Number of mapped reads/alignments in the sample.
  2  mod      <str>  Data transformation mode; must be 'frag' or 'norm' (default: ${mod}).
  3  egs      <int>  Effective genome size (default: ${egs} [appropriate for S. cerevisiae]).
  4  rnd      <int>  Number of decimal places to round to (default: ${rnd}).
  5  csv_bin  <str>  Comma-separated list of bin sizes to use (default: ${csv_bin}).

Returns:
  Prints a comma-delimited list of depth values, one per bin size, in the order given by 'csv_bin'.

Examples:
  '''bash
  _calculate_dep_arr 12851824
  _calculate_dep_arr 12851824 frag
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
    elif [[ -z "${rnd}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 4, 'rnd', is missing."
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

    for var in dep egs rnd; do
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
            _calculate_dep_fct "${dep}" "${bin}" "${egs}" "${mod}" "${rnd}"
        )" )
    done

    IFS=',' printf "%s\n" "${arr_dep[*]}"
}


function _compute_dep_all() {
    local dep="${1:-}"                      # Number of alignments in sample BAM
    local rnd="${2:-24}"                    # Number of decimals to round to
    local csv_bin="${3:-${DEP_BINS_DFLT}}"  # Comma-delimited bin sizes
    local -a arr_dm_fr arr_dm_nm
    local output
    local show_help

    show_help=$(cat << EOM
Usage:
  _compute_dep_all dep [rnd] [csv_bin]

Description:
  Compute minimum input depth values needed downstream for both fragment-length-adjusted signal ('frag') and "normalized coverage" ('norm'), following Dickson et al., Sci Rep 2023.

  Internally, this function calls '_calculate_dep_arr' twice and concatenates the two comma-delimited outputs.

Positional arguments:
  1  dep      <int>  Number of mapped reads/alignments in the sample.
  2  rnd      <int>  Number of decimal places to round to (default: ${rnd}).
  3  csv_bin  <str>  Comma-separated list of bin sizes to use (default: ${csv_bin}).

Returns:
  Prints one comma-delimited string containing the 'frag' depth values followed by the 'norm' depth values, each in the order given by 'csv_bin'.

Examples:
  '''bash
  _compute_dep_all 12851824 24
  _compute_dep_all 12851824 24 1,10,25,50
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
    elif [[ -z "${rnd}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'rnd', is missing."
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

    if ! [[ "${rnd}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'rnd', must be a positive integer:" \
            "'${rnd}'."
        return 1
    fi

    #  Note: detailed value validation is delegated to '_calculate_dep_arr'
    IFS=',' read -r -a arr_dm_fr < <(
        _calculate_dep_arr "${dep}" "frag" "12157105" "${rnd}" "${csv_bin}"
    ) || return 1
    IFS=',' read -r -a arr_dm_nm < <(
        _calculate_dep_arr "${dep}" "norm" "12157105" "${rnd}" "${csv_bin}"
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
Usage:
  _generate_fmt_str num_fld

Description:
  Construct a tab-delimited 'printf' format string containing 'num_fld' '%s' fields and no trailing tab.

Positional argument:
  1  num_fld  <int>  Number of fields in the output row.

Returns:
  Prints a format string suitable for tab-delimited row output.

Examples:
  '''bash
  _generate_fmt_str 3
  _generate_fmt_str 14
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


function _get_fil_out_part() {
    local fil_out="${1:-}"
    local idx="${2:-}"
    local tag="${3:-part}"
    local idx_pad
    local show_help

    show_help=$(cat << EOM
Usage:
  _get_fil_out_part [-h|--hlp|--help] fil_out idx [tag]

Description:
  Construct the per-sample part-file path derived from a base output-table path.

Positional arguments:
  1  fil_out  <str>  Base output-table path.
  2  idx      <int>  Zero-based sample index.
  3  tag      <str>  Label inserted into the part-file name (default: ${tag}).

Returns:
  Prints a path of the form:

      <fil_out>.part.<tag>.<zero-padded idx>

Examples:
  '''bash
  _get_fil_out_part results.tsv 0
  _get_fil_out_part results.tsv 12 siq
  _get_fil_out_part results.tsv 12 spike
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
    printf '%s.part.%s.%s\n' "${fil_out}" "${tag}" "${idx_pad}"
}


#  Compute siQ-ChIP alpha scaling factor and related values for a sample
#+
#+ Workflow function that processes a sample using global variables; extracts
#+ siQ-ChIP metadata from a TSV table, computes alignment counts and average
#+ fragment lengths, and minimum input depth values
#+
# shellcheck disable=SC2154
function process_samp_siq() {
    local idx="${1:-}"  # Array sample index

    #  Declare local variables
    local fil_ip fil_in siq mass_ip mass_in vol_all vol_in
    local typ_ip typ_in dep_ip dep_in len_ip len_in v
    # local fmt_str  # Reserved in case formatted output generation is revived
    local len_ip_ovrd len_in_ovrd fil_out_part
    local -a fields arr_dm
    local show_help

    show_help=$(cat << 'EOM'
Usage:
  process_samp_siq [-h|--hlp|--help] idx

Description:
  Processes a single sample (array index 'idx') to compute the siQ-ChIP alpha scaling factor and related QC values.

  For the IP/input BAM pair at 'idx', this function
    - parses siQ-ChIP metadata,
    - detects PE/SE for each BAM,
    - counts alignments,
    - computes average fragment lengths (PE only; SE uses a provided default),
    - computes minimum input depth values across the currently configured default bin sizes, and
    - writes a per-sample tab-separated results row.

Positional argument:
  1  idx  <int>  Zero-based array index into arr_mip/arr_min.

Expected global variables:
  Read:
    arr_mip, arr_min                Sample BAM arrays for IP and input
    scr_met, cfg_met, tbl_met, eqn  siQ-ChIP metadata parser, YAML configuration, table, and equation
    threads                         Samtools decompression threads
    len_def                         Default fragment length to use for SE libraries
    aln_typ                         Optional override: 'pe', 'se', or 'auto' (default: auto)
    scr_siq, scr_spk                Python entry points for scaling-factor calculation
    rnd                             Rounding precision for depth-factor outputs
    debug                           If 'true', prints 'debug_var' lines

  Write:
    fil_out                         Base path used to derive the per-sample results file

Output row fields (in order):
  fil_ip, fil_in, siq, eqn, mass_ip, mass_in, vol_all, vol_in, dep_ip, dep_in, len_ip, len_in, <depth factors for bins 1,5,10,20,30,40,50 for frag and norm modes>

Dependencies:
  This function expects 'validate_var_file' and 'debug_var' to be available in the current shell via the helper-layer sourcing block above.

Notes:
  - For single-end libraries, TLEN is not meaningful; this function requires 'len_def' (or an optional 4th positional argument, 'len_lcl') to supply a fixed fragment length.
  - End-type detection is per-file in case of mixed inputs.
  - Individual alignment files are expected to be a single input type: either entirely PE or entirely SE.
  - Minimum input depth values are currently computed only for the default bin sizes 1,5,10,20,30,40,50. The workflow helper does not yet expose a way to pass custom bin sizes through to '_compute_dep_all'.
EOM
    )
    #TODO: if custom bin sizes are threaded through to this function, the
    #+     output-column semantics will no longer necessarily correspond to
    #+     the fixed default bin sequence 1,5,10,20,30,40,50

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

    #  Assign BAM files based on sample index
    fil_ip="${arr_mip[idx]}"
    fil_in="${arr_min[idx]}"

    #  Check that BAM files exist
    validate_var_file "fil_ip" "${fil_ip}" "${idx}" || return 1
    validate_var_file "fil_in" "${fil_in}" "${idx}" || return 1

    if [[ "${debug}" == "true" ]]; then
        debug_var "idx=${idx}" "fil_ip=${fil_ip}" "fil_in=${fil_in}"
    fi

    #  Parse siQ-ChIP metadata, assigning global variables
    _parse_metadata \
        "${scr_met}" "${fil_ip}" "${tbl_met}" "${cfg_met}" "${eqn}" || {
            echo_err_func "${FUNCNAME[0]}" \
                "failed while running '_parse_metadata' for IP BAM '${fil_ip}'."
            return 1
        }

    #  Determine end type per BAM (robust if inputs ever mix)
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

    #  Count alignments in BAM files
    dep_ip="$(_get_dep_idx mip "${idx}")"
    if [[ -z "${dep_ip}" ]]; then
        dep_ip="$(_count_align_bam "${threads}" "${fil_ip}" "${typ_ip}")" || {
            echo_err_func "${FUNCNAME[0]}" \
                "failed while counting alignments for IP BAM '${fil_ip}'" \
                "with type '${typ_ip}'."
            return 1
        }
    fi

    dep_in="$(_get_dep_idx min "${idx}")"
    if [[ -z "${dep_in}" ]]; then
        dep_in="$(_count_align_bam "${threads}" "${fil_in}" "${typ_in}")" || {
            echo_err_func "${FUNCNAME[0]}" \
                "failed while counting alignments for input BAM '${fil_in}'" \
                "with type '${typ_in}'."
            return 1
        }
    fi

    #  Compute average fragment lengths for BAM files; overrides take
    #+ precedence over TLEN or 'len_def'
    len_ip_ovrd="$(_get_len_idx mip "${idx}")"
    len_in_ovrd="$(_get_len_idx min "${idx}")"

    if [[ -n "${len_ip_ovrd}" ]]; then
        len_ip="${len_ip_ovrd}"
    else
        len_ip="$(
            _calculate_frag_avg \
                "${threads}" "${fil_ip}" "${typ_ip}" "${len_def:-}"
        )" || {
            echo_err_func "${FUNCNAME[0]}" \
                "failed while computing average fragment length for IP BAM" \
                "'${fil_ip}' with type '${typ_ip}'."
            return 1
        }
    fi

    if [[ -n "${len_in_ovrd}" ]]; then
        len_in="${len_in_ovrd}"
    else
        len_in="$(
            _calculate_frag_avg \
                "${threads}" "${fil_in}" "${typ_in}" "${len_def:-}"
        )" || {
            echo_err_func "${FUNCNAME[0]}" \
                "failed while computing average fragment length for input" \
                "BAM '${fil_in}' with type '${typ_in}'."
            return 1
        }
    fi

    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "eqn=${eqn}"         "rnd=${rnd}" \
            "mass_ip=${mass_ip}" "mass_in=${mass_in}" \
            "vol_all=${vol_all}" "vol_in=${vol_in}" \
            "dep_ip=${dep_ip}"   "dep_in=${dep_in}" \
            "len_ip=${len_ip}"   "len_in=${len_in}"
    fi

    #  Compute siQ-ChIP alpha scaling factor
    siq=$(
        _compute_scl_fct \
            "siq" "${scr_siq}" "${scr_spk}" \
            --eqn     "${eqn}"     --rnd     "${rnd}" \
            --mass_ip "${mass_ip}" --mass_in "${mass_in}" \
            --vol_all "${vol_all}" --vol_in  "${vol_in}" \
            --dep_ip  "${dep_ip}"  --dep_in  "${dep_in}" \
            --len_ip  "${len_ip}"  --len_in  "${len_in}"
    ) || return 1

    #  Compute minimum input depth values
    #TODO: thread custom bin sizes through this call when the higher-level
    #+     interface exposes them
    IFS="," read -r -a arr_dm < <(
        _compute_dep_all "${dep_in}" "${rnd}"
    ) || return 1
    #TODO: add regression tests for process-substitution patterns of the form
    #+
    #+         read ... < <(producer_cmd ...) || return 1
    #+
    #+     across the codebase; confirm, for representative producer failures,
    #+     that stderr/error messages still appear and that failure status is
    #+     handled the way the surrounding code expects

    # #  Dynamically construct string for output formatting
    # #TODO: 'fmt_str' is computed but never used; is this dead code now? Can I
    # #+     still make use of this? Do I even need to?
    # fmt_str=$( _generate_fmt_str $((12 + ${#arr_dm[@]})) )

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
        _get_fil_out_part "${fil_out}" "${idx}" "siq"
    )" || return 1

    if ! {
            fields=(
                "${fil_ip}" "${fil_in}" "${siq}" "${eqn}"
                "${mass_ip}" "${mass_in}" "${vol_all}" "${vol_in}"
                "${dep_ip}" "${dep_in}" "${len_ip}" "${len_in}"
            )
            fields+=( "${arr_dm[@]}" )

            printf '%s' "${fields[0]}"
            for v in "${fields[@]:1}"; do printf '\t%s' "${v}"; done
            printf '\n'
        } > "${fil_out_part}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to write per-sample results file '${fil_out_part}'."
        return 1
    fi
    #TODO: concatenate per-sample part files deterministically in the higher-
    #+     level driver instead of writing directly to shared 'fil_out'
}


#  Compute spike-in scaling factor and related values for a sample
#+
#+ Workflow function that processes a sample using global variables; computes
#+ alignment counts and minimum input depth values
#+
# shellcheck disable=SC2154
function process_samp_spike() {
    local idx="${1:-}"  # Array sample index

    #  Declare local variables
    local mp mn sp sn typ_mp typ_sp typ_mn typ_sn
    local num_mp num_sp num_mn num_sn coef_lcl coef_val v fil_out_part
    # local fmt_str  # Reserved in case formatted output generation is revived
    local -a fields arr_dm
    local show_help

    show_help=$(cat << 'EOM'
Usage:
  process_samp_spike [-h|--hlp|--help] idx

Description:
  Processes a single sample (array index 'idx') to compute a spike-in scaling factor and related QC values.

  For the main/spike IP and input BAMs at 'idx', this function detects paired- or single-end sequenced read alignment status per BAM, counts alignments, computes minimum input depth values across common bin sizes, and writes a per-sample tab-separated results row.

Positional argument:
  1  idx  <int>  Zero-based array index into arrays: arr_mip (main IP), arr_sip (spike IP), arr_min (main input), arr_sin (spike input)

Expected global variables:
  Read:
    arr_mip, arr_sip, arr_min, arr_sin   Sample BAM arrays
    threads                              Samtools decompression threads
    aln_typ                              Optional override: 'pe', 'se', or 'auto' (default: auto)
    scr_siq, scr_spk                     Python entry points for scaling-factor calculation
    coef_spk                             Spike-in coefficient to request from 'calculate_scaling_factor_spike.py' (default: fractional)
    rnd                                  Rounding precision for depth-factor outputs
    debug                                If 'true', prints 'debug_var' lines

  Write:
    fil_out                              Base path used to derive the per-sample results file

Output row fields (in order):
  mp, sp, mn, sn, coef_val, num_mp, num_sp, num_mn, num_sn, <depth factors for bins 1,5,10,20,30,40,50 for frag and norm modes>

Dependencies:
  This function expects 'validate_var_file' and 'debug_var' to be available in the current shell via the helper-layer sourcing block above.

Notes:
  - End-type detection is per-file in case of mixed inputs.
  - Individual alignment files are expected to be a single input type: either entirely PE or entirely SE.
  - Minimum input depth values are currently computed only for the default bin sizes 1,5,10,20,30,40,50. The workflow helper does not yet expose a way to pass custom bin sizes through to '_compute_dep_all'.
EOM
    )
    #TODO: if custom bin sizes are threaded through to this function, the
    #+     output-column semantics will no longer necessarily correspond to
    #+     the fixed default bin sequence 1,5,10,20,30,40,50

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

    #  Assign BAM files based on sample index
    mp="${arr_mip[idx]}"
    mn="${arr_min[idx]}"
    sp="${arr_sip[idx]}"
    sn="${arr_sin[idx]}"

    #  Check that BAM files exist
    validate_var_file "mp" "${mp}" "${idx}" || return 1
    validate_var_file "sp" "${sp}" "${idx}" || return 1
    validate_var_file "mn" "${mn}" "${idx}" || return 1
    validate_var_file "sn" "${sn}" "${idx}" || return 1

    if [[ "${debug}" == "true" ]]; then
        debug_var "idx=${idx}" "mp=${mp}" "mn=${mn}" "sp=${sp}" "sn=${sn}"
    fi

    #  Determine end type per BAM
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

    #  Count alignments in BAM files
    num_mp="$(_get_dep_idx mip "${idx}")"
    if [[ -z "${num_mp}" ]]; then
        num_mp="$(_count_align_bam "${threads}" "${mp}" "${typ_mp}")" || {
            echo_err_func "${FUNCNAME[0]}" \
                "failed while counting alignments for main IP BAM '${mp}'" \
                "with type '${typ_mp}'."
            return 1
        }
    fi

    num_mn="$(_get_dep_idx min "${idx}")"
    if [[ -z "${num_mn}" ]]; then
        num_mn="$(_count_align_bam "${threads}" "${mn}" "${typ_mn}")" || {
            echo_err_func "${FUNCNAME[0]}" \
                "failed while counting alignments for main input BAM '${mn}'" \
                "with type '${typ_mn}'."
            return 1
        }
    fi

    num_sp="$(_get_dep_idx sip "${idx}")"
    if [[ -z "${num_sp}" ]]; then
        num_sp="$(_count_align_bam "${threads}" "${sp}" "${typ_sp}")" || {
            echo_err_func "${FUNCNAME[0]}" \
                "failed while counting alignments for spike IP BAM '${sp}'" \
                "with type '${typ_sp}'."
            return 1
        }
    fi

    num_sn="$(_get_dep_idx sin "${idx}")"
    if [[ -z "${num_sn}" ]]; then
        num_sn="$(_count_align_bam "${threads}" "${sn}" "${typ_sn}")" || {
            echo_err_func "${FUNCNAME[0]}" \
                "failed while counting alignments for spike input BAM" \
                "'${sn}' with type '${typ_sn}'."
            return 1
        }
    fi

    coef_lcl="${coef_spk:-fractional}"

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
            "coef_spk=${coef_lcl}" "rnd=${rnd}"
    fi

    #  Compute requested spike-in coefficient
    coef_val="$(
        _compute_scl_fct \
            "spike" "${scr_siq}" "${scr_spk}" \
            --coef    "${coef_lcl}" \
            --format  "plain" \
            --main_ip "${num_mp}" --spike_ip "${num_sp}" \
            --main_in "${num_mn}" --spike_in "${num_sn}" \
            --rnd     "${rnd}"
    )" || {
        echo_err_func "${FUNCNAME[0]}" \
            "failed while computing spike-in coefficient '${coef_lcl}' for" \
            "sample idx '${idx}'."
        return 1
    }

    #  Compute minimum input depth values
    #TODO: thread custom bin sizes through this call when the higher-level
    #+     interface exposes them
    IFS="," read -r -a arr_dm < <(
        _compute_dep_all "${num_mn}" "${rnd}"
    ) || return 1
    #TODO: add regression tests for process-substitution patterns of the form
    #+
    #+         read ... < <(producer_cmd ...) || return 1
    #+
    #+     across the codebase; confirm, for representative producer failures,
    #+     that stderr/error messages still appear and that failure status is
    #+     handled the way the surrounding code expects

    # #  Dynamically construct string for output formatting
    # #TODO: 'fmt_str' is computed but never used; is this dead code now? Can I
    # #+     still make use of this? Do I even need to?
    # fmt_str=$( _generate_fmt_str $((9 + ${#arr_dm[@]})) )

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
        _get_fil_out_part "${fil_out}" "${idx}" "spike"
    )" || return 1

    if ! {
            fields=(
                "${mp}" "${sp}" "${mn}" "${sn}" "${coef_val}"
                "${num_mp}" "${num_sp}" "${num_mn}" "${num_sn}"
            )
            fields+=( "${arr_dm[@]}" )

            printf '%s' "${fields[0]}"
            for v in "${fields[@]:1}"; do printf '\t%s' "${v}"; done
            printf '\n'
        } > "${fil_out_part}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to write per-sample results file '${fil_out_part}'."
        return 1
    fi
    #TODO: concatenate per-sample part files deterministically in the higher-
    #+     level driver instead of writing directly to shared 'fil_out'
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
