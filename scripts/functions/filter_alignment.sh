#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: filter_alignment.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# _parse_args_filter_alignment
# _validate_args_filter_alignment
# _check_chr_alignment
# _finalize_alignment_filter
# _sanitize_pg_value
# _build_filter_pg_cl
# _filter_sam_header_chr
# _filter_sam_chr
# _cleanup_filter_alignment_tmp
# filter_alignment_sc
# filter_alignment_sp


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
    _dir_src_alignment="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    source "${_dir_src_alignment}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_alignment}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_alignment}" \
        check_args check_inputs check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_alignment
}


#  Parse keyword arguments for alignment-filter helpers
#+
#+ - Assign parsed values back to caller variables named by the caller
#+ - Support shared arguments for S. cerevisiae and S. pombe filtering
#+ - Argument support can be extended for additional organisms
#+ - Restrict organism-specific flags such as '--tg' and '--mtr' to the
#+   appropriate caller
#+ - Assumes caller-level help handling so this helper only parses arguments
#+   after public entry-point validation has begun
function _parse_args_filter_alignment() {
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
    local ref_fa_ref="${8:-}"
    local show_help="${9:-}"
    shift 9

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

            -r|--ref|--ref[_-]fa|--reference)
                require_optarg "${1}" "${2:-}" "${func}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                printf -v "${ref_fa_ref}" '%s' "${2}"
                shift 2
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


#  Validate common required arguments for alignment-filter helpers
#+ - Check thread count, input-file existence, and output-directory existence
#+ - Shared by organism-specific alignment-filter entry-point functions
function _validate_args_filter_alignment() {
    local func="${1:-}"
    local threads="${2:-}"
    local infile="${3:-}"
    local outfile="${4:-}"
    local ref_fa="${5:-}"
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

    case "${outfile,,}" in
        *.bam|*.cram) : ;;
        *)
            echo_err_func "${func}" \
                "'--outfile' must end in '.bam' or '.cram': '${outfile}'."
            return 1
            ;;
    esac

    if [[ "${infile,,}" == *.cram || "${outfile,,}" == *.cram ]]; then
        if [[ -z "${ref_fa}" ]]; then
            echo_err_func "${func}" \
                "'--ref_fa' is required when input or output is CRAM."
            return 1
        fi

        validate_var_file "ref_fa" "${ref_fa}" || return 1
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


#  Print unique reference-sequence names present in an alignment file
#+ - Used for optional post-filter chromosome checking
function _check_chr_alignment() {
    local outfile="${1:-}"
    local ref_fa="${2:-}"
    local -a ref_arg=()

    validate_var_file "outfile" "${outfile}" || return 1

    if [[ "${outfile,,}" == *.cram ]]; then
        validate_var_file "ref_fa" "${ref_fa}" || return 1
        ref_arg=( -T "${ref_fa}" )
    fi

    samtools view -h "${ref_arg[@]}" "${outfile}" \
        | awk '!/^@/ { print $3 }' \
        | sort \
        | uniq
}
#MAYBE: move to a shared helper script later if reused elsewhere


#  "Finalize" a filtered alignment file
#+ - Index the output and optionally print retained chromosome names
function _finalize_alignment_filter() {
    local threads="${1:-}"
    local outfile="${2:-}"
    local chk_chr="${3:-false}"
    local ref_fa="${4:-}"

    samtools index -@ "${threads}" "${outfile}" || return 1

    if [[ "${chk_chr}" == "true" ]]; then
        _check_chr_alignment "${outfile}" "${ref_fa}" || return 1
    fi
}


#  Sanitize one value for inclusion in a SAM @PG field
function _sanitize_pg_value() {
    local value="${1:-}"

    value="${value//$'\t'/ }"
    value="${value//$'\r'/ }"
    value="${value//$'\n'/ }"

    printf '%s' "${value}"
}


#  Build the @PG CL field describing one filter_alignments operation
function _build_filter_pg_cl() {
    local func="${1:-}"
    local retain="${2:-}"
    local infile="${3:-}"
    local outfile="${4:-}"
    local mito="${5:-false}"
    local tg="${6:-false}"
    local mtr="${7:-false}"
    local chk_chr="${8:-false}"
    local ref_fa="${9:-}"
    local out_ext="${outfile##*.}"
    local cl

    cl="$(_sanitize_pg_value "${func}")"
    cl+=" retain=$(_sanitize_pg_value "${retain}")"
    cl+=" infile=$(_sanitize_pg_value "${infile}")"
    cl+=" outfile=$(_sanitize_pg_value "${outfile}")"
    cl+=" mito=$(_sanitize_pg_value "${mito}")"

    if [[ "${retain}" == "sp" ]]; then
        cl+=" tg=$(_sanitize_pg_value "${tg}")"
        cl+=" mtr=$(_sanitize_pg_value "${mtr}")"
    fi

    cl+=" chk_chr=$(_sanitize_pg_value "${chk_chr}")"
    cl+=" out_ext=$(_sanitize_pg_value "${out_ext}")"

    if [[ -n "${ref_fa}" ]]; then
        cl+=" ref_fa=$(_sanitize_pg_value "${ref_fa}")"
    fi

    printf '%s' "${cl}"
}


#  Filter SAM header lines and append a valid filter_alignments @PG record
function _filter_sam_header_chr() {
    local func="${1:-}"
    local chrs="${2:-}"
    local pg_id="${3:-}"
    local pg_pn="${4:-filter_alignments}"
    local pg_cl="${5:-}"

    if [[ -z "${func}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'func', is missing."
        return 1
    elif [[ -z "${chrs}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'chrs', is missing."
        return 1
    fi

    awk \
        -v chrs="${chrs}" \
        -v pg_id="${pg_id}" \
        -v pg_pn="${pg_pn}" \
        -v pg_cl="${pg_cl}" \
        '
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

            function pg_tag_value(line, tag, fields, i, n) {
                n = split(line, fields, "\t")
                for (i = 1; i <= n; i++) {
                    if (fields[i] ~ "^" tag ":") {
                        return substr(fields[i], length(tag) + 2)
                    }
                }

                return ""
            }

            function emit_pg(    id, i) {
                if (pg_id == "" || pg_done) {
                    return
                }

                id = pg_id
                i = 1
                while (id in pg_ids) {
                    id = pg_id "." i
                    i++
                }

                print "@PG\tID:" id "\tPN:" pg_pn "\tCL:" pg_cl
                pg_done = 1
            }

            /^@HD/ {
                print
                next
            }

            /^@SQ/ {
                if (keep_sq_line($0)) {
                    print
                }
                next
            }

            /^@PG/ {
                id = pg_tag_value($0, "ID")
                if (id != "") {
                    pg_ids[id] = 1
                }
                print
                next
            }

            /^@/ {
                print
                next
            }

            END {
                emit_pg()
            }
        '
}


#  Filter a SAM file by retained reference-sequence names
#+ - Keep non-'@SQ' header lines, retain only matching '@SQ' lines, and keep
#+   only alignments whose reference sequence is in the supplied chromosome set
#+ - Used when chromosome filtering is done by rewriting SAM header/body
#+   content rather than by alignment filtering plus reheadering
function _filter_sam_chr() {
    local func="${1:-}"
    local infile="${2:-}"
    local outfile="${3:-}"
    local chrs="${4:-}"
    local pg_id="${5:-}"
    local pg_pn="${6:-filter_alignments}"
    local pg_cl="${7:-}"

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
        awk \
            -v chrs="${chrs}" \
            -v pg_id="${pg_id}" \
            -v pg_pn="${pg_pn}" \
            -v pg_cl="${pg_cl}" \
            '
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

            function pg_tag_value(line, tag, fields, i, n) {
                n = split(line, fields, "\t")
                for (i = 1; i <= n; i++) {
                    if (fields[i] ~ "^" tag ":") {
                        return substr(fields[i], length(tag) + 2)
                    }
                }

                return ""
            }

            function emit_pg(    id, i) {
                if (pg_id == "" || pg_done) {
                    return
                }

                id = pg_id
                i = 1
                while (id in pg_ids) {
                    id = pg_id "." i
                    i++
                }

                print "@PG\tID:" id "\tPN:" pg_pn "\tCL:" pg_cl
                pg_done = 1
            }

            /^@SQ/ {
                if (keep_sq_line($0)) {
                    print
                }
                next
            }

            /^@PG/ {
                id = pg_tag_value($0, "ID")
                if (id != "") {
                    pg_ids[id] = 1
                }
                print
                next
            }

            /^@/ {
                print
                next
            }

            ($3 in chrom_map) {
                emit_pg()
                print
            }

            END {
                emit_pg()
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


#  Remove temporary files used during alignment filtering
#+ - Best-effort cleanup helper for explicit call sites
function _cleanup_filter_alignment_tmp() {
    local pth_1="${1:-}"
    local pth_2="${2:-}"

    [[ -n "${pth_1}" ]] && rm -f "${pth_1}"
    [[ -n "${pth_2}" ]] && rm -f "${pth_2}"
}


#  Filter and reheader a BAM or CRAM file for S. cerevisiae chromosomes
#+ - Public entry-point function for main-organism yeast alignment filtering
#+ - Uses direct filtering followed by reheadering
function filter_alignment_sc() {
    local threads=1
    local infile=""
    local outfile=""
    local mito=false
    local tg=false
    local mtr=false
    local chk_chr=false
    local ref_fa=""
    local chrs
    local outdir outbam hdr_rh bam_filter bam_rh_init bam_rh_sort pg_cl
    local -a ref_arg=()
    local show_help

    show_help=$(cat << EOM
Usage:
  filter_alignment_sc
    [--help] [--threads <int>] --infile <str> --outfile <str> [--ref_fa <str>] [--mito] [--chk_chr]

Description:
  Filter and reheader a BAM or CRAM file for S. cerevisiae chromosomes, writing BAM or CRAM output.

Keyword arguments:
   -t, --threads  <int>  Number of threads to use (default: ${threads}).
   -i, --infile   <str>  Coordinate-sorted BAM or CRAM infile.
   -o, --outfile  <str>  Filtered BAM or CRAM outfile.
   -r, --ref_fa   <str>  Reference FASTA required when input or output is CRAM.
   -m, --mito     <flg>  Retain mitochondrial chromosome (optional).
  -cc, --chk_chr  <flg>  Check chromosomes in filtered outfile (optional).

Returns:
  Creates a BAM or CRAM outfile filtered and reheadered for S. cerevisiae chromosomes at the specified path.

Dependencies:
  - Programs
    + Bash >= 4.4
    + grep
    + Samtools >= 1.21

Examples:
  '''bash
  filter_alignment_sc
      --threads 4
      --infile "sample.bam"
      --outfile "sample.sc.bam"

  filter_alignment_sc
      --threads 4
      --infile "sample.cram"
      --outfile "sample.sc.bam"
      --ref_fa "reference.fa"

  filter_alignment_sc
      --threads 4
      --infile "sample.bam"
      --outfile "sample.sc.cram"
      --ref_fa "reference.fa"

  filter_alignment_sc
      --threads 4
      --infile "sample.bam"
      --outfile "sample.sc.bam"
      --mito
      --chk_chr
  '''

#TODO:
  - Somewhere and somehow, need to handle more than S. cerevisiae as "main" organism.
  - Support SAM input too.
EOM
    )

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    _parse_args_filter_alignment \
        "${FUNCNAME[0]}" \
        "sc" \
        threads infile outfile mito tg mtr chk_chr ref_fa \
        "${show_help}" \
        "$@" \
        || return $?

    _validate_args_filter_alignment \
        "${FUNCNAME[0]}" \
        "${threads}" \
        "${infile}" \
        "${outfile}" \
        "${ref_fa}" \
        || return 1

    if [[ "${infile,,}" == *.cram ]]; then
        ref_arg=( -T "${ref_fa}" )
    fi

    chrs="I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI"
    if [[ "${mito}" == "true" ]]; then
        chrs="${chrs} Mito"
    fi

    pg_cl="$(
        _build_filter_pg_cl \
            "${FUNCNAME[0]}" \
            sc \
            "${infile}" \
            "${outfile}" \
            "${mito}" \
            false \
            false \
            "${chk_chr}" \
            "${ref_fa}"
    )"

    outdir="$(dirname "${outfile}")"
    outbam="$(basename "${outfile}")"
    if [[ "${outfile,,}" == *.cram ]]; then
        hdr_rh="${outdir}/tmp.${outbam%.cram}.header.sam"
        bam_filter="${outdir}/tmp.${outbam%.cram}.filter.bam"
        bam_rh_init="${outdir}/tmp.${outbam%.cram}.rehead.bam"
        bam_rh_sort="${outdir}/tmp.${outbam%.cram}.sort.bam"
    else
        hdr_rh="${outdir}/tmp.${outbam%.bam}.header.sam"
        bam_filter="${outfile}"
        bam_rh_init="${outdir}/rehead.${outbam}"
        bam_rh_sort="${outdir}/txt_rh_sort.${outbam}"
    fi

    # shellcheck disable=SC2086
    if ! \
        samtools view \
            -@ "${threads}" \
            -b \
            -o "${bam_filter}" \
            "${ref_arg[@]}" \
            "${infile}" \
            ${chrs}
    then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to filter '${infile}'."
        return 1
    fi

    if ! \
        samtools view -H "${bam_filter}" \
            | _filter_sam_header_chr \
                "${FUNCNAME[0]}" \
                "${chrs}" \
                "${FUNCNAME[0]}" \
                "filter_alignments" \
                "${pg_cl}" \
                > "${hdr_rh}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to build filtered header for '${bam_filter}'."
        if [[ "${outfile,,}" == *.cram ]]; then
            rm -f "${hdr_rh}" "${bam_filter}" "${bam_rh_init}" "${bam_rh_sort}"
        else
            rm -f "${hdr_rh}" "${bam_rh_init}" "${bam_rh_sort}"
        fi
        return 1
    fi

    if ! \
        samtools reheader "${hdr_rh}" "${bam_filter}" > "${bam_rh_init}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to reheader '${bam_filter}'."
        if [[ "${outfile,,}" == *.cram ]]; then
            rm -f "${hdr_rh}" "${bam_filter}" "${bam_rh_init}" "${bam_rh_sort}"
        else
            rm -f "${hdr_rh}" "${bam_rh_init}" "${bam_rh_sort}"
        fi
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
        if [[ "${outfile,,}" == *.cram ]]; then
            rm -f "${hdr_rh}" "${bam_filter}" "${bam_rh_init}" "${bam_rh_sort}"
        else
            rm -f "${hdr_rh}" "${bam_rh_init}" "${bam_rh_sort}"
        fi
        return 1
    fi

    if [[ "${outfile,,}" == *.cram ]]; then
        if ! \
            samtools view \
                -@ "${threads}" \
                -C \
                -T "${ref_fa}" \
                -o "${outfile}" \
                "${bam_rh_sort}"
        then
            echo_err_func "${FUNCNAME[0]}" \
                "failed to convert filtered BAM intermediate to CRAM."
            rm -f "${hdr_rh}" "${bam_filter}" "${bam_rh_init}" "${bam_rh_sort}"
            return 1
        fi
    else
        if ! mv -f "${bam_rh_sort}" "${outfile}"; then
            echo_err_func "${FUNCNAME[0]}" \
                "failed to rename sorted reheadered BAM file."
            rm -f "${hdr_rh}" "${bam_rh_init}" "${bam_rh_sort}"
            return 1
        fi
    fi

    if [[ "${outfile,,}" == *.cram ]]; then
        if ! rm -f "${hdr_rh}" "${bam_filter}" "${bam_rh_init}" "${bam_rh_sort}"; then
            echo_err_func "${FUNCNAME[0]}" \
                "failed to delete filter BAM intermediates."
            return 1
        fi
    elif ! rm -f "${hdr_rh}" "${bam_rh_init}"; then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to delete reheader BAM intermediate."
        return 1
    fi

    _finalize_alignment_filter "${threads}" "${outfile}" "${chk_chr}" "${ref_fa}" \
        || return 1
}


#  Filter and reheader a BAM or CRAM file for S. pombe chromosomes
#+ - Public entry-point function for spike-in-organism yeast alignment filtering
#+ - Uses SAM-level rewriting because retained references may include optional
#+   contigs that are easier to handle through header/body filtering
function filter_alignment_sp() {
    local threads=1
    local infile=""
    local outfile=""
    local mito=false
    local tg=false
    local mtr=false
    local chk_chr=false
    local ref_fa=""
    local chrs
    local outdir outbase pth_in pth_out pg_cl
    local -a ref_arg=()
    local show_help

    show_help=$(cat << EOM
Usage:
  filter_alignment_sp
    [--help] [--threads <int>] --infile <str> --outfile <str> [--ref_fa <str>] [--mito] [--tg] [--mtr] [--chk_chr]

Description:
  Filter and reheader a BAM or CRAM file for S. pombe chromosomes, writing BAM or CRAM output.

Keyword arguments:
   -t, --threads  <int>  Number of threads to use (default: ${threads}).
   -i, --infile   <str>  Coordinate-sorted BAM or CRAM infile.
   -o, --outfile  <str>  Filtered BAM or CRAM outfile.
   -r, --ref_fa   <str>  Reference FASTA required when input or output is CRAM.
   -m, --mito     <flg>  Retain SP_Mito chromosome.
  -tg, --tg       <flg>  Retain SP_II_TG chromosome.
  -mr, --mtr      <flg>  Retain SP_MTR chromosome.
  -cc, --chk_chr  <flg>  Check chromosomes in filtered outfile.

Returns:
  Creates a BAM or CRAM outfile filtered and reheadered for S. pombe chromosomes at the specified path.

Dependencies:
  - Programs
    + AWK >= 5
    + Bash >= 4.4
    + Samtools >= 1.21

Examples:
  '''bash
  filter_alignment_sp
      --threads 4
      --infile "sample.bam"
      --outfile "sample.sp.bam"

  filter_alignment_sp
      --threads 4
      --infile "sample.cram"
      --outfile "sample.sp.bam"
      --ref_fa "reference.fa"

  filter_alignment_sp
      --threads 4
      --infile "sample.bam"
      --outfile "sample.sp.cram"
      --ref_fa "reference.fa"

  filter_alignment_sp
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
  - Support SAM input too.
EOM
    )

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    _parse_args_filter_alignment \
        "${FUNCNAME[0]}" \
        "sp" \
        threads infile outfile mito tg mtr chk_chr ref_fa \
        "${show_help}" \
        "$@" \
        || return $?

    _validate_args_filter_alignment \
        "${FUNCNAME[0]}" \
        "${threads}" \
        "${infile}" \
        "${outfile}" \
        "${ref_fa}" \
        || return 1

    if [[ "${infile,,}" == *.cram ]]; then
        ref_arg=( -T "${ref_fa}" )
    fi

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

    pg_cl="$(
        _build_filter_pg_cl \
            "${FUNCNAME[0]}" \
            sp \
            "${infile}" \
            "${outfile}" \
            "${mito}" \
            "${tg}" \
            "${mtr}" \
            "${chk_chr}" \
            "${ref_fa}"
    )"

    outdir="$(dirname "${outfile}")"
    outbase="$(basename "${outfile}")"
    outbase="${outbase%.bam}"
    outbase="${outbase%.cram}"
    pth_in="${outdir}/tmp.${outbase}.in.sam"
    pth_out="${outdir}/tmp.${outbase}.out.sam"

    # shellcheck disable=SC2086
    if ! \
        samtools view \
            -@ "${threads}" \
            -h \
            -o "${pth_in}" \
            "${ref_arg[@]}" \
            "${infile}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to generate '${pth_in}'."
        _cleanup_filter_alignment_tmp "${pth_in}" "${pth_out}"
        return 1
    fi

    if ! \
        _filter_sam_chr \
            "${FUNCNAME[0]}" \
            "${pth_in}" \
            "${pth_out}" \
            "${chrs}" \
            "${FUNCNAME[0]}" \
            "filter_alignments" \
            "${pg_cl}"
    then
        _cleanup_filter_alignment_tmp "${pth_in}" "${pth_out}"
        return 1
    fi

    if [[ "${outfile,,}" == *.cram ]]; then
        if ! \
            samtools view \
                -C \
                -T "${ref_fa}" \
                -o "${outfile}" \
                "${pth_out}"
        then
            echo_err_func "${FUNCNAME[0]}" \
                "failed to generate '${outfile}'."
            _cleanup_filter_alignment_tmp "${pth_in}" "${pth_out}"
            return 1
        fi
    else
        if ! \
            samtools view -b "${pth_out}" > "${outfile}"
        then
            echo_err_func "${FUNCNAME[0]}" \
                "failed to generate '${outfile}'."
            _cleanup_filter_alignment_tmp "${pth_in}" "${pth_out}"
            return 1
        fi
    fi

    #  Remove intermediates immediately on success
    _cleanup_filter_alignment_tmp "${pth_in}" "${pth_out}"

    _finalize_alignment_filter "${threads}" "${outfile}" "${chk_chr}" "${ref_fa}" \
        || return 1
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
