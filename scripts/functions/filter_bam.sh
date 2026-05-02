#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: filter_bam.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# _parse_args_filter_bam
# _validate_args_filter_bam
# _check_chr_bam
# _finalize_bam_filter
# _filter_sam_chr
# _cleanup_filter_bam_tmp
# filter_bam_sc
# filter_bam_sp


#  Require Bash >= 4.4 before defining functions
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be sourced or run under Bash >= 4.4." >&2

    if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
        return 1
    else
        exit 1
    fi
elif (( BASH_VERSINFO[0] < 5 )); then
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
    _dir_src_bam="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    source "${_dir_src_bam}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_bam}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_bam}" \
        check_args check_inputs check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_bam
}


#  Parse keyword arguments for BAM-filter helpers
#+
#+ - Assign parsed values back to caller variables named by the caller
#+ - Support shared arguments for S. cerevisiae and S. pombe filtering
#+ - Argument support can be extended for additional organisms
#+ - Restrict organism-specific flags such as '--tg' and '--mtr' to the
#+   appropriate caller
#+ - Assumes caller-level help handling so this helper only parses arguments
#+   after public entry-point validation has begun
function _parse_args_filter_bam() {
    local func="${1:-}"
    local chr_nam="${2:-}"
    shift 2

    local threads_ref="${1:-}"
    local infile_ref="${2:-}"
    local outfile_ref="${3:-}"
    local mito_ref="${4:-}"
    local tg_ref="${5:-}"
    local mtr_ref="${6:-}"
    local chk_chr_ref="${7:-}"
    local show_help="${8:-}"
    shift 8

    #  Parse arguments, assigning parsed values back to caller variables whose
    #+ names are passed in
    #+
    #+ Use scalar '--infile'/'--outfile' aliases as the canonical helper API;
    #+ deprecated list-style long aliases remain accepted for compatibility
    #+ with older direct helper calls, but should not be used in new callers
    while (( $# > 0 )); do
        case "${1}" in
            -t|--thr|--threads)
                require_optarg "${1}" "${2:-}" "${func}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                printf -v "${threads_ref}" '%s' "${2}"
                shift 2
                ;;

            #TODO: old plural/list-style aliases remain accepted temporarily
            -i|-fi|--infile|--fil[_-]in|--infiles|--csv[_-]infile|--csv[_-]infiles)
                require_optarg "${1}" "${2:-}" "${func}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                printf -v "${infile_ref}" '%s' "${2}"
                shift 2
                ;;

            #TODO: old plural/list-style aliases remain accepted temporarily
            -o|-fo|--outfile|--fil[_-]out|--outfiles|--csv[_-]outfile|--csv[_-]outfiles)
                require_optarg "${1}" "${2:-}" "${func}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                printf -v "${outfile_ref}" '%s' "${2}"
                shift 2
                ;;

            -m|--mit|--mito)
                printf -v "${mito_ref}" '%s' true
                shift 1
                ;;

            -cc|--chk[_-]chr)
                printf -v "${chk_chr_ref}" '%s' true
                shift 1
                ;;

            -tg|--tg)
                if [[ "${chr_nam}" != "sp" ]]; then
                    echo_err_func "${func}" \
                        "option '${1}' is not supported here."
                    return 1
                fi
                printf -v "${tg_ref}" '%s' true
                shift 1
                ;;

            -mr|--mr|--mtr)
                if [[ "${chr_nam}" != "sp" ]]; then
                    echo_err_func "${func}" \
                        "option '${1}' is not supported here."
                    return 1
                fi
                printf -v "${mtr_ref}" '%s' true
                shift 1
                ;;

            *)
                #TODO: update to current error messaging
                echo "## Unknown argument passed: '${1}' ##" >&2
                echo >&2
                echo "${show_help}" >&2
                return 1
                ;;
        esac
    done
}


#  Validate common required arguments for BAM-filter helpers
#+ - Check thread count, input-file existence, and output-directory existence
#+ - Shared by organism-specific BAM-filter entry-point functions
function _validate_args_filter_bam() {
    local func="${1:-}"
    local threads="${2:-}"
    local infile="${3:-}"
    local outfile="${4:-}"
    local outdir

    if [[ -z "${threads}" ]]; then
        echo_err_func "${func}" \
            "'--threads' is required."
        return 1
    fi

    if [[ ! "${threads}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err_func "${func}" \
            "'--threads' is '${threads}' but must be a positive integer."
        return 1
    fi

    if [[ -z "${infile}" ]]; then
        echo_err_func "${func}" \
            "'--infile' is required."
        return 1
    fi

    validate_var_file "infile" "${infile}" || return 1

    if [[ -z "${outfile}" ]]; then
        echo_err_func "${func}" \
            "'--outfile' is required."
        return 1
    fi

    outdir="$(dirname "${outfile}")"
    if [[ ! -d "${outdir}" ]]; then
        echo_err_func "${func}" \
            "directory associated with '--outfile' does not exist:" \
            "'${outdir}'."
        return 1
    fi

    return 0
}


#  Print unique reference-sequence names present in a BAM file
#+ - Used for optional post-filter chromosome checking
function _check_chr_bam() {
    local outfile="${1:-}"

    validate_var_file "outfile" "${outfile}" || return 1

    samtools view -h "${outfile}" \
        | awk '!/^@/ { print $3 }' \
        | sort \
        | uniq
}
#MAYBE: move to a shared helper script later if reused elsewhere


#  "Finalize" a filtered BAM file
#+ - Index the BAM and optionally print retained chromosome names
function _finalize_bam_filter() {
    local threads="${1:-}"
    local outfile="${2:-}"
    local chk_chr="${3:-false}"

    samtools index -@ "${threads}" "${outfile}" || return 1

    if [[ "${chk_chr}" == "true" ]]; then
        _check_chr_bam "${outfile}" || return 1
    fi
}


#  Filter a SAM file by retained reference-sequence names
#+ - Keep non-'@SQ' header lines, retain only matching '@SQ' lines, and keep
#+   only alignments whose reference sequence is in the supplied chromosome set
#+ - Used when chromosome filtering is done by rewriting SAM header/body
#+   content rather than by BAM filtering plus reheadering
function _filter_sam_chr() {
    local func="${1:-}"
    local infile="${2:-}"
    local outfile="${3:-}"
    local chrs="${4:-}"

    if [[ -z "${func}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'func', is missing."
        return 1
    elif [[ -z "${infile}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'infile', is missing."
        return 1
    elif [[ -z "${outfile}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'outfile', is missing."
        return 1
    elif [[ -z "${chrs}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 4, 'chrs', is missing."
        return 1
    fi

    validate_var_file "infile" "${infile}" || return 1

    #TODO: explicitly use 'gawk'?
    if ! \
        awk -v chrs="${chrs}" '
            BEGIN {
                split(chrs, arr_chrs, " ")
                for (i in arr_chrs) {
                    chrom_map[arr_chrs[i]] = 1
                }
            }

            function keep_sq_line(line, a, b, sn) {
                if (line !~ /SN:/) {
                    return 0
                }

                split(line, a, /SN:/)
                split(a[2], b, /[ \t]/)
                sn = b[1]

                return (sn in chrom_map)
            }

            /^@SQ/ {
                if (keep_sq_line($0)) {
                    print
                }
                next
            }

            /^@/ {
                print
                next
            }

            ($3 in chrom_map) {
                print
            }
        ' "${infile}" \
            > "${outfile}"
    then
        echo_err_func "${func}" \
            "failed to generate '${outfile}'."
        return 1
    fi
}
#TODO: benchmark GNU Awk-based filtering approach against alternatives if
#+     performance becomes limiting for larger genomes


#  Remove temporary files used during BAM filtering
#+ - Best-effort cleanup helper for explicit call sites
function _cleanup_filter_bam_tmp() {
    local pth_1="${1:-}"
    local pth_2="${2:-}"

    [[ -n "${pth_1}" ]] && rm -f "${pth_1}"
    [[ -n "${pth_2}" ]] && rm -f "${pth_2}"
}


#  Filter and reheader a BAM file for S. cerevisiae chromosomes
#+ - Public entry-point function for main-organism yeast BAM filtering
#+ - Uses direct BAM filtering followed by BAM reheadering
function filter_bam_sc() {
    local threads=1
    local infile=""
    local outfile=""
    local mito=false
    local tg=false
    local mtr=false
    local chk_chr=false
    local chrs pattern
    local outdir outbam bam_rh_init bam_rh_sort
    local show_help

    show_help=$(cat << EOM
Usage:
  filter_bam_sc
    [--help] [--threads <int>] --infile <str> --outfile <str> [--mito] [--chk_chr]

Description:
  Filter and reheader a BAM file for S. cerevisiae chromosomes.

Keyword arguments:
   -t, --threads  <int>  Number of threads to use (default: ${threads}).
   -i, --infile   <str>  Coordinate-sorted BAM infile.
   -o, --outfile  <str>  Filtered BAM outfile.
   -m, --mito     <flg>  Retain mitochondrial chromosome (optional).
  -cc, --chk_chr  <flg>  Check chromosomes in filtered BAM outfile (optional).

Returns:
  Creates a BAM outfile filtered and reheadered for S. cerevisiae chromosomes at the specified path.

Dependencies:
  - Programs
    + Bash >= 4.4
    + grep
    + Samtools >= 1.21

Examples:
  '''bash
  filter_bam_sc
      --threads 4
      --infile "sample.bam"
      --outfile "sample.sc.bam"

  filter_bam_sc
      --threads 4
      --infile "sample.bam"
      --outfile "sample.sc.bam"
      --mito
      --chk_chr
  '''

#TODO:
  - Somewhere and somehow, need to handle more than S. cerevisiae as "main" organism.
  - The filter activity needs to be recorded in the BAM header.
  - Support SAM and CRAM too.
EOM
    )

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    _parse_args_filter_bam \
        "${FUNCNAME[0]}" \
        "sc" \
        threads infile outfile mito tg mtr chk_chr \
        "${show_help}" \
        "$@" \
        || return $?

    _validate_args_filter_bam \
        "${FUNCNAME[0]}" \
        "${threads}" \
        "${infile}" \
        "${outfile}" \
        || return 1

    chrs="I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI"
    if [[ "${mito}" == "true" ]]; then
        chrs="${chrs} Mito"
    fi

    pattern="^@HD|^@SQ.*SN:($(echo "${chrs}" | tr ' ' '|'))|^@PG"

    outdir="$(dirname "${outfile}")"
    outbam="$(basename "${outfile}")"
    bam_rh_init="${outdir}/rehead.${outbam}"
    bam_rh_sort="${outdir}/txt_rh_sort.${outbam}"

    # shellcheck disable=SC2086
    if ! \
        samtools view \
            -@ "${threads}" \
            -b \
            -o "${outfile}" \
            "${infile}" \
            ${chrs}
    then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to filter '${infile}'."
        return 1
    fi

    if ! \
        samtools reheader \
            -c "grep -E '${pattern}'" \
            "${outfile}" \
                > "${bam_rh_init}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to reheader '${outfile}'."
        rm -f "${bam_rh_init}" "${bam_rh_sort}"
        return 1
    fi

    if ! \
        samtools sort \
            -@ "${threads}" \
            -o "${bam_rh_sort}" \
            "${bam_rh_init}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to sort reheadered BAM intermediate."
        rm -f "${bam_rh_init}" "${bam_rh_sort}"
        return 1
    fi

    if ! mv -f "${bam_rh_sort}" "${outfile}"; then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to rename sorted reheadered BAM file."
        rm -f "${bam_rh_init}" "${bam_rh_sort}"
        return 1
    fi

    if ! rm -f "${bam_rh_init}"; then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to delete reheader BAM intermediate."
        return 1
    fi

    _finalize_bam_filter "${threads}" "${outfile}" "${chk_chr}" || return 1
}


#  Filter and reheader a BAM file for S. pombe chromosomes
#+ - Public entry-point function for spike-in-organism yeast BAM filtering
#+ - Uses SAM-level rewriting because retained references may include optional
#+   contigs that are easier to handle through header/body filtering
function filter_bam_sp() {
    local threads=1
    local infile=""
    local outfile=""
    local mito=false
    local tg=false
    local mtr=false
    local chk_chr=false
    local chrs
    local outdir insam outsam pth_in pth_out
    local show_help

    show_help=$(cat << EOM
Usage:
  filter_bam_sp
    [--help] [--threads <int>] --infile <str> --outfile <str> [--mito] [--tg] [--mtr] [--chk_chr]

Description:
  Filter and reheader a BAM file for S. pombe chromosomes.

Keyword arguments:
   -t, --threads  <int>  Number of threads to use (default: ${threads}).
   -i, --infile   <str>  Coordinate-sorted BAM infile.
   -o, --outfile  <str>  Filtered BAM outfile.
   -m, --mito     <flg>  Retain SP_Mito chromosome.
  -tg, --tg       <flg>  Retain SP_II_TG chromosome.
  -mr, --mtr      <flg>  Retain SP_MTR chromosome.
  -cc, --chk_chr  <flg>  Check chromosomes in filtered BAM outfile.

Returns:
  Creates a BAM outfile filtered and reheadered for S. pombe chromosomes at the specified path.

Dependencies:
  - Programs
    + AWK >= 5
    + Bash >= 4.4
    + Samtools >= 1.21

Examples:
  '''bash
  filter_bam_sp
      --threads 4
      --infile "sample.bam"
      --outfile "sample.sp.bam"

  filter_bam_sp
      --threads 4
      --infile "sample.bam"
      --outfile "sample.sp.bam"
      --mito
      --tg
      --mtr
      --chk_chr
  '''

#TODO:
  - Somewhere and somehow, need to handle more than S. pombe as "spike-in" organism.
  - #IMPORTANT: The filter activity needs to be recorded in the BAM header.
  - Support SAM and CRAM too.
EOM
    )

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    _parse_args_filter_bam \
        "${FUNCNAME[0]}" \
        "sp" \
        threads infile outfile mito tg mtr chk_chr \
        "${show_help}" \
        "$@" \
        || return $?

    _validate_args_filter_bam \
        "${FUNCNAME[0]}" \
        "${threads}" \
        "${infile}" \
        "${outfile}" \
        || return 1

    chrs="SP_I SP_II SP_III"
    if [[ "${tg}" == "true" ]]; then
        chrs="SP_II_TG ${chrs}"
    fi

    if [[ "${mtr}" == "true" ]]; then
        chrs="${chrs} SP_MTR"
    fi

    if [[ "${mito}" == "true" ]]; then
        chrs="${chrs} SP_Mito"
    fi

    outdir="$(dirname "${outfile}")"
    insam="$(basename "${infile/.bam/.sam}")"
    outsam="$(basename "${outfile/.bam/.sam}")"
    pth_in="${outdir}/${insam}"
    pth_out="${outdir}/${outsam}"

    # shellcheck disable=SC2086
    if ! \
        samtools view \
            -@ "${threads}" \
            -h \
            -o "${pth_in}" \
            "${infile}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to generate '${pth_in}'."
        _cleanup_filter_bam_tmp "${pth_in}" "${pth_out}"
        return 1
    fi

    if ! \
        _filter_sam_chr \
            "${FUNCNAME[0]}" \
            "${pth_in}" \
            "${pth_out}" \
            "${chrs}"
    then
        _cleanup_filter_bam_tmp "${pth_in}" "${pth_out}"
        return 1
    fi

    if ! \
        samtools view -b "${pth_out}" > "${outfile}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to generate '${outfile}'."
        _cleanup_filter_bam_tmp "${pth_in}" "${pth_out}"
        return 1
    fi

    #  Remove intermediates immediately on success
    _cleanup_filter_bam_tmp "${pth_in}" "${pth_out}"

    _finalize_bam_filter "${threads}" "${outfile}" "${chk_chr}" || return 1
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
