#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: filter_alignment.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in design, development, and documentation, with all output
# reviewed, edited, and approved by the author.
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
    _dir_src_alignment="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    # shellcheck source=lib/bash/core/source_helpers.sh
    source "${_dir_src_alignment}/../core/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_alignment}/../core/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_alignment}/.." \
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
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _parse_args_filter_alignment
    [--help] func chr_nam threads_ref fil_in_ref fil_out_ref mito_ref tg_ref mtr_ref chk_chr_ref ref_fa_ref -- [--threads <int>] [--fil_in <file>] [--fil_out <file>] [--mito] [--chk_chr] [--ref_fa <file>] [--tg] [--mtr]

  Parse shared alignment-filter arguments into caller-owned variables.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  chr_nam : str
    Organism code that controls organism-specific options.

  3  threads_ref : str
    Name of the caller variable that receives parsed threads.

  4  fil_in_ref : str
    Name of the caller variable that receives the input path.

  5  fil_out_ref : str
    Name of the caller variable that receives the output path.

  6  mito_ref : str
    Name of the caller variable that receives the mitochondrial flag.

  7  tg_ref : str
    Name of the caller variable that receives the TG flag.

  8  mtr_ref : str
    Name of the caller variable that receives the MTR flag.

  9  chk_chr_ref : str
    Name of the caller variable that receives the chromosome-check flag.

  10  ref_fa_ref : str
    Name of the caller variable that receives the reference FASTA path.

  -t, --thr, --threads : int
    Number of threads to use. Keyword option accepted after the '--' delimiter.

  -fi, --fil_in : file
    Input file path. Keyword option accepted after the '--' delimiter.

  -fo, --fil_out : file
    Output file path. Keyword option accepted after the '--' delimiter.

  -m, --mit, --mito : flag
    Retain mitochondrial chromosome. Keyword option accepted after the '--' delimiter.

  -cc, --chk_chr : flag
    Check chromosomes in output alignment files. Keyword option accepted after the '--' delimiter.

  -rf, --ref_fa : file
    Reference FASTA file. Keyword option accepted after the '--' delimiter.

  -tg, --tg : flag
    Retain SP_II_TG chromosome. Keyword option accepted after the '--' delimiter only when 'chr_nam=sp'.

  -mr, --mr, --mtr : flag
    Retain SP_MTR chromosome. Keyword option accepted after the '--' delimiter only when 'chr_nam=sp'.

Returns
-------
  Returns 0 when parsing succeeds; returns 1 on invalid input or unsupported organism-specific options.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Display the parser contract used by both organism-specific callers.
    '''bash
    _parse_args_filter_alignment --help
    '''

  2. Initialize caller variables without supplying optional filter arguments.
    '''bash
    threads=1; fil_in=; fil_out=; mito=false; tg=false; mtr=false
    chk_chr=false; ref_fa=; help_text='filter help'
    _parse_args_filter_alignment filter_alignment_sc sc
      threads fil_in fil_out mito tg mtr chk_chr ref_fa "\${help_text}" --
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    shift 2

    local threads_ref="${1:-}"
    local fil_in_ref="${2:-}"
    local fil_out_ref="${3:-}"
    local mito_ref="${4:-}"
    local tg_ref="${5:-}"
    local mtr_ref="${6:-}"
    local chk_chr_ref="${7:-}"
    local ref_fa_ref="${8:-}"
    local show_help="${9:-}"
    shift 9

    #  Consume the documented parser-wiring delimiter when supplied. Keep it
    #+ optional for compatibility with older internal callers.
    if [[ "${1:-}" == "--" ]]; then
        shift 1
    fi

    #  Parse arguments, assigning parsed values back to caller variables whose
    #+ names are passed in
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

            -fi|--fil_in)
                require_optarg "${1}" "${2:-}" "${func}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                printf -v "${fil_in_ref}" '%s' "${2}"
                shift 2
                ;;

            -fo|--fil_out)
                require_optarg "${1}" "${2:-}" "${func}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                printf -v "${fil_out_ref}" '%s' "${2}"
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

            -rf|--ref[_-]fa)
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
                echo_err_func "${FUNCNAME[0]}" \
                    "unknown option/parameter passed: '${1}'."
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
    local fil_in="${3:-}"
    local fil_out="${4:-}"
    local ref_fa="${5:-}"
    local outdir
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _validate_args_filter_alignment
    [--help] func threads fil_in fil_out ref_fa

  Validate shared alignment-filter arguments before dispatch.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  threads : str
    Number of threads to use. Requested thread count.

  3  fil_in : file
    Input file path. Input BAM or CRAM file.

  4  fil_out : file
    Output file path. Output BAM or CRAM file.

  5  ref_fa : file
    Reference FASTA file. Reference FASTA required for CRAM.

Returns
-------
  Returns 0 when validation succeeds; returns 1 when a required argument is missing or invalid.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - dirname

Examples
--------
  1. Validate a prepared BAM input and an output path in an existing directory.
    '''bash
    tmp="\$(mktemp -d)"
    trap 'rm -r -- "\${tmp}"' EXIT
    touch "\${tmp}/input.bam"
    _validate_args_filter_alignment filter_alignment_sc 1
      "\${tmp}/input.bam" "\${tmp}/output.bam" ''
    '''

  2. Confirm that a nonpositive thread count is rejected before dispatch.
    '''bash
    if ! \\
        _validate_args_filter_alignment filter_alignment_sc 0 in.bam out.bam ''
    then
        printf '%s\n' 'invalid thread count rejected as expected'
    fi
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

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

    if [[ -z "${fil_in}" ]]; then
        echo_err_func "${func}" \
            "'--fil_in' is required."
        return 1
    fi

    validate_var_file "fil_in" "${fil_in}" || return 1

    if [[ -z "${fil_out}" ]]; then
        echo_err_func "${func}" \
            "'--fil_out' is required."
        return 1
    fi

    case "${fil_out,,}" in
        *.bam|*.cram) : ;;
        *)
            echo_err_func "${func}" \
                "'--fil_out' must end in '.bam' or '.cram': '${fil_out}'."
            return 1
            ;;
    esac

    if [[ "${fil_in,,}" == *.cram || "${fil_out,,}" == *.cram ]]; then
        if [[ -z "${ref_fa}" ]]; then
            echo_err_func "${func}" \
                "'--ref_fa' is required when input or output is CRAM."
            return 1
        fi

        validate_var_file "ref_fa" "${ref_fa}" || return 1
    fi

    outdir="$(dirname "${fil_out}")"
    if [[ ! -d "${outdir}" ]]; then
        echo_err_func "${func}" \
            "directory associated with '--fil_out' does not exist:" \
            "'${outdir}'."
        return 1
    fi

    return 0
}


#  Print unique reference-sequence names present in an alignment file
#+ - Used for optional post-filter chromosome checking
function _check_chr_alignment() {
    local fil_out="${1:-}"
    local ref_fa="${2:-}"
    local -a arr_ref_arg=()
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _check_chr_alignment
    [--help] fil_out ref_fa

  Print the unique reference-sequence names observed in an alignment file.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fil_out : file
    Output file path. BAM or CRAM file to inspect.

  2  ref_fa : file
    Reference FASTA file. Reference FASTA required for CRAM input.

Returns
-------
  Prints chromosome names to stdout. Returns 1 when validation or Samtools access fails.

Notes
-----
  Runtime requirements:
    - awk
    - bash >= 4.4
    - Reference FASTA and required index (when processing CRAM)
    - samtools
    - sort
    - uniq

Examples
--------
  1. Print chromosomes observed in a filtered BAM file.
    '''bash
    _check_chr_alignment work/sample.sc.bam ''
    '''

  2. Print chromosomes from CRAM using its reference FASTA.
    '''bash
    _check_chr_alignment work/sample.sp.cram \\
        tests/fixtures/filter_alignments/reference/filter_sc_sp.fa
    '''
EOM
    )

    if [[ "${fil_out}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    validate_var_file "fil_out" "${fil_out}" || return 1

    if [[ "${fil_out,,}" == *.cram ]]; then
        validate_var_file "ref_fa" "${ref_fa}" || return 1
        arr_ref_arg=( -T "${ref_fa}" )
    fi

    samtools view -h "${arr_ref_arg[@]}" "${fil_out}" \
        | awk '!/^@/ { print $3 }' \
        | sort \
        | uniq
}
#MAYBE: move to a shared helper script later if reused elsewhere


#  "Finalize" a filtered alignment file
#+ - Index the output and optionally print retained chromosome names
function _finalize_alignment_filter() {
    local threads="${1:-}"
    local fil_out="${2:-}"
    local chk_chr="${3:-false}"
    local ref_fa="${4:-}"
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _finalize_alignment_filter
    [--help] threads fil_out [chk_chr] [ref_fa]

  Index a filtered alignment file and optionally print retained chromosomes.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  threads : int
    Number of threads to use for Samtools indexing.

  2  fil_out : file
    Output file path. Filtered alignment file to index.

  3  chk_chr : flag
    Check chromosomes in output alignment files. If true, print retained chromosomes after indexing.

  4  ref_fa : file
    Reference FASTA file. Reference FASTA required when chromosome checking a CRAM file.

Returns
-------
  Returns 0 when indexing succeeds; returns 1 if indexing or chromosome inspection fails.

Notes
-----
  Runtime requirements:
    - awk (when 'chk_chr' is true)
    - bash >= 4.4
    - Reference FASTA and required index (when 'chk_chr' is true and processing CRAM)
    - samtools
    - sort (when 'chk_chr' is true)
    - uniq (when 'chk_chr' is true)

Examples
--------
  1. Index a filtered BAM without printing chromosome names.
    '''bash
    _finalize_alignment_filter 2 work/sample.sc.bam false
    '''

  2. Index a CRAM and print the retained chromosomes using its reference.
    '''bash
    _finalize_alignment_filter 4 work/sample.sp.cram true \\
        tests/fixtures/filter_alignments/reference/filter_sc_sp.fa
    '''
EOM
    )

    if [[ "${threads}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    samtools index -@ "${threads}" "${fil_out}" || return 1

    if [[ "${chk_chr}" == "true" ]]; then
        _check_chr_alignment "${fil_out}" "${ref_fa}" || return 1
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
    local fil_in="${3:-}"
    local fil_out="${4:-}"
    local mito="${5:-false}"
    local tg="${6:-false}"
    local mtr="${7:-false}"
    local chk_chr="${8:-false}"
    local ref_fa="${9:-}"
    local out_ext="${fil_out##*.}"
    local cl

    cl="$(_sanitize_pg_value "${func}")"
    cl+=" retain=$(_sanitize_pg_value "${retain}")"
    cl+=" fil_in=$(_sanitize_pg_value "${fil_in}")"
    cl+=" fil_out=$(_sanitize_pg_value "${fil_out}")"
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
    local fil_in="${2:-}"
    local fil_out="${3:-}"
    local chrs="${4:-}"
    local pg_id="${5:-}"
    local pg_pn="${6:-filter_alignments}"
    local pg_cl="${7:-}"
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _filter_sam_chr
    [--help] func fil_in fil_out chrs [pg_id] [pg_pn] [pg_cl]

  Rewrite a SAM file so that only selected chromosomes and their matching
  alignments remain.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  fil_in : file
    Input file path. Input SAM, BAM, or CRAM path.

  3  fil_out : file
    Output file path. Output SAM path.

  4  chrs : str
    Space-delimited chromosome set to retain.

  5  pg_id : str
    Program-group identifier to append when emitting a new @PG line.

  6  pg_pn : str
    Program-group name to append when emitting a new @PG line.

  7  pg_cl : str
    Sanitized program-group command line.

Returns
-------
  Returns 0 when the filtered SAM is written successfully; returns 1 on validation or AWK failure.

Notes
-----
  Runtime requirements:
    - awk
    - bash >= 4.4

Examples
--------
  1. Retain canonical S. cerevisiae chromosomes from the committed SAM fixture.
    '''bash
    tmp="\$(mktemp -d)"
    trap 'rm -r -- "\${tmp}"' EXIT
    _filter_sam_chr filter_alignment_sc
      tests/fixtures/filter_alignments/sam/filter_sc_sp.sam
      "\${tmp}/sample.sc.sam" 'I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI'
    '''

  2. Retain selected S. pombe chromosomes and append a program-group record.
    '''bash
    tmp="\$(mktemp -d)"
    trap 'rm -r -- "\${tmp}"' EXIT
    _filter_sam_chr \\
        filter_alignment_sp \\
        tests/fixtures/filter_alignments/sam/filter_sc_sp.sam \\
        "\${tmp}/sample.sp.sam" \\
        'SP_I SP_II SP_III SP_Mito' \\
        filter_alignment_sp \\
        filter_alignments \\
        'retain=sp mito=true'
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ -z "${func}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'func', is missing."
        return 1
    elif [[ -z "${fil_in}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'fil_in', is missing."
        return 1
    elif [[ -z "${fil_out}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'fil_out', is missing."
        return 1
    elif [[ -z "${chrs}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 4, 'chrs', is missing."
        return 1
    fi

    validate_var_file "fil_in" "${fil_in}" || return 1

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
        ' "${fil_in}" \
            > "${fil_out}"
    then
        echo_err_func "${func}" \
            "failed to generate '${fil_out}'."
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
    local fil_in=""
    local fil_out=""
    local mito=false
    local tg=false
    local mtr=false
    local chk_chr=false
    local ref_fa=""
    local chrs
    local outdir outbam hdr_rh bam_filter bam_rh_init bam_rh_sort pg_cl
    local -a arr_ref_arg=()
    local show_help

    show_help=$(cat << EOM
Usage
-----
  filter_alignment_sc
    [--help] [--threads <int>] --fil_in <file> --fil_out <file> [--ref_fa <file>] [--mito] [--chk_chr]

  Filter and reheader a BAM or CRAM file for S. cerevisiae chromosomes, writing BAM or CRAM output.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -t, --thr, --threads : int
    Number of threads to use (default: ${threads}).

  -fi, --fil_in : file
    Input file path. Coordinate-sorted BAM or CRAM fil_in.

  -fo, --fil_out : file
    Output file path. Filtered BAM or CRAM fil_out.

  -rf, --ref_fa : file
    Reference FASTA file. Reference FASTA required when input or output is CRAM.

  -m, --mit, --mito : flag
    Retain mitochondrial chromosome.

  -cc, --chk_chr : flag
    Check chromosomes in output alignment files.

Returns
-------
  Creates a BAM or CRAM fil_out filtered and reheadered for S. cerevisiae chromosomes at the specified path.

Notes
-----
  Runtime requirements:
    - awk
    - basename
    - bash >= 4.4
    - dirname
    - Input BAM or CRAM index
    - mv (when writing BAM)
    - Reference FASTA and required index (when processing CRAM)
    - rm
    - samtools
    - sort (when '--chk_chr' is specified)
    - uniq (when '--chk_chr' is specified)

Examples
--------
  1. Retain canonical S. cerevisiae chromosomes in BAM output.
    '''bash
    filter_alignment_sc \\
        --threads 2 \\
        --fil_in work/filter_sc_sp.bam \\
        --fil_out work/filter_sc_sp.sc.bam
    '''

  2. Convert CRAM input to filtered CRAM while retaining mitochondria.
    '''bash
    filter_alignment_sc \\
        --threads 4 \\
        --fil_in work/filter_sc_sp.cram \\
        --fil_out work/filter_sc_sp.sc.cram \\
        --ref_fa tests/fixtures/filter_alignments/reference/filter_sc_sp.fa \\
        --mito \\
        --chk_chr
    '''
EOM
    )

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    _parse_args_filter_alignment \
        "${FUNCNAME[0]}" \
        "sc" \
        threads fil_in fil_out mito tg mtr chk_chr ref_fa \
        "${show_help}" \
        -- \
        "$@" \
        || return $?

    _validate_args_filter_alignment \
        "${FUNCNAME[0]}" \
        "${threads}" \
        "${fil_in}" \
        "${fil_out}" \
        "${ref_fa}" \
        || return 1

    if [[ "${fil_in,,}" == *.cram ]]; then
        arr_ref_arg=( -T "${ref_fa}" )
    fi

    chrs="I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI"
    if [[ "${mito}" == "true" ]]; then
        chrs="${chrs} Mito"
    fi

    pg_cl="$(
        _build_filter_pg_cl \
            "${FUNCNAME[0]}" \
            sc \
            "${fil_in}" \
            "${fil_out}" \
            "${mito}" \
            false \
            false \
            "${chk_chr}" \
            "${ref_fa}"
    )"

    outdir="$(dirname "${fil_out}")"
    outbam="$(basename "${fil_out}")"
    if [[ "${fil_out,,}" == *.cram ]]; then
        hdr_rh="${outdir}/tmp.${outbam%.cram}.header.sam"
        bam_filter="${outdir}/tmp.${outbam%.cram}.filter.bam"
        bam_rh_init="${outdir}/tmp.${outbam%.cram}.rehead.bam"
        bam_rh_sort="${outdir}/tmp.${outbam%.cram}.sort.bam"
    else
        hdr_rh="${outdir}/tmp.${outbam%.bam}.header.sam"
        bam_filter="${fil_out}"
        bam_rh_init="${outdir}/rehead.${outbam}"
        bam_rh_sort="${outdir}/txt_rh_sort.${outbam}"
    fi

    # shellcheck disable=SC2086
    if ! \
        samtools view \
            -@ "${threads}" \
            -b \
            -o "${bam_filter}" \
            "${arr_ref_arg[@]}" \
            "${fil_in}" \
            ${chrs}
    then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to filter '${fil_in}'."
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
        if [[ "${fil_out,,}" == *.cram ]]; then
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
        if [[ "${fil_out,,}" == *.cram ]]; then
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
        if [[ "${fil_out,,}" == *.cram ]]; then
            rm -f "${hdr_rh}" "${bam_filter}" "${bam_rh_init}" "${bam_rh_sort}"
        else
            rm -f "${hdr_rh}" "${bam_rh_init}" "${bam_rh_sort}"
        fi
        return 1
    fi

    if [[ "${fil_out,,}" == *.cram ]]; then
        if ! \
            samtools view \
                -@ "${threads}" \
                -C \
                -T "${ref_fa}" \
                -o "${fil_out}" \
                "${bam_rh_sort}"
        then
            echo_err_func "${FUNCNAME[0]}" \
                "failed to convert filtered BAM intermediate to CRAM."
            rm -f "${hdr_rh}" "${bam_filter}" "${bam_rh_init}" "${bam_rh_sort}"
            return 1
        fi
    else
        if ! mv -f "${bam_rh_sort}" "${fil_out}"; then
            echo_err_func "${FUNCNAME[0]}" \
                "failed to rename sorted reheadered BAM file."
            rm -f "${hdr_rh}" "${bam_rh_init}" "${bam_rh_sort}"
            return 1
        fi
    fi

    if [[ "${fil_out,,}" == *.cram ]]; then
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

    _finalize_alignment_filter "${threads}" "${fil_out}" "${chk_chr}" "${ref_fa}" \
        || return 1
}


#  Filter and reheader a BAM or CRAM file for S. pombe chromosomes
#+ - Public entry-point function for spike-in-organism yeast alignment filtering
#+ - Uses SAM-level rewriting because retained references may include optional
#+   contigs that are easier to handle through header/body filtering
function filter_alignment_sp() {
    local threads=1
    local fil_in=""
    local fil_out=""
    local mito=false
    local tg=false
    local mtr=false
    local chk_chr=false
    local ref_fa=""
    local chrs
    local outdir outbase pth_in pth_out pg_cl
    local -a arr_ref_arg=()
    local show_help

    show_help=$(cat << EOM
Usage
-----
  filter_alignment_sp
    [--help] [--threads <int>] --fil_in <file> --fil_out <file> [--ref_fa <file>] [--mito] [--tg] [--mtr] [--chk_chr]

  Filter and reheader a BAM or CRAM file for S. pombe chromosomes, writing BAM or CRAM output.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -t, --thr, --threads : int
    Number of threads to use (default: ${threads}).

  -fi, --fil_in : file
    Input file path. Coordinate-sorted BAM or CRAM fil_in.

  -fo, --fil_out : file
    Output file path. Filtered BAM or CRAM fil_out.

  -rf, --ref_fa : file
    Reference FASTA file. Reference FASTA required when input or output is CRAM.

  -m, --mit, --mito : flag
    Retain mitochondrial chromosome. Uses chromosome 'SP_Mito'.

  -tg, --tg : flag
    Retain SP_II_TG chromosome.

  -mr, --mr, --mtr : flag
    Retain SP_MTR chromosome.

  -cc, --chk_chr : flag
    Check chromosomes in output alignment files.

Returns
-------
  Creates a BAM or CRAM fil_out filtered and reheadered for S. pombe chromosomes at the specified path.

Notes
-----
  Runtime requirements:
    - awk
    - basename
    - bash >= 4.4
    - dirname
    - Reference FASTA and required index (when processing CRAM)
    - rm
    - samtools
    - sort (when '--chk_chr' is specified)
    - uniq (when '--chk_chr' is specified)

Examples
--------
  1. Retain canonical S. pombe chromosomes in BAM output.
    '''bash
    filter_alignment_sp \\
        --threads 2 \\
        --fil_in work/filter_sc_sp.bam \\
        --fil_out work/filter_sc_sp.sp.bam
    '''

  2. Write filtered CRAM while retaining optional S. pombe contigs.
    '''bash
    filter_alignment_sp \\
        --threads 4 \\
        --fil_in work/filter_sc_sp.bam \\
        --fil_out work/filter_sc_sp.sp.cram \\
        --ref_fa tests/fixtures/filter_alignments/reference/filter_sc_sp.fa \\
        --mito \\
        --tg \\
        --mtr \\
        --chk_chr
    '''
EOM
    )

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    _parse_args_filter_alignment \
        "${FUNCNAME[0]}" \
        "sp" \
        threads fil_in fil_out mito tg mtr chk_chr ref_fa \
        "${show_help}" \
        -- \
        "$@" \
        || return $?

    _validate_args_filter_alignment \
        "${FUNCNAME[0]}" \
        "${threads}" \
        "${fil_in}" \
        "${fil_out}" \
        "${ref_fa}" \
        || return 1

    if [[ "${fil_in,,}" == *.cram ]]; then
        arr_ref_arg=( -T "${ref_fa}" )
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
            "${fil_in}" \
            "${fil_out}" \
            "${mito}" \
            "${tg}" \
            "${mtr}" \
            "${chk_chr}" \
            "${ref_fa}"
    )"

    outdir="$(dirname "${fil_out}")"
    outbase="$(basename "${fil_out}")"
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
            "${arr_ref_arg[@]}" \
            "${fil_in}"
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

    if [[ "${fil_out,,}" == *.cram ]]; then
        if ! \
            samtools view \
                -C \
                -T "${ref_fa}" \
                -o "${fil_out}" \
                "${pth_out}"
        then
            echo_err_func "${FUNCNAME[0]}" \
                "failed to generate '${fil_out}'."
            _cleanup_filter_alignment_tmp "${pth_in}" "${pth_out}"
            return 1
        fi
    else
        if ! \
            samtools view -b "${pth_out}" > "${fil_out}"
        then
            echo_err_func "${FUNCNAME[0]}" \
                "failed to generate '${fil_out}'."
            _cleanup_filter_alignment_tmp "${pth_in}" "${pth_out}"
            return 1
        fi
    fi

    #  Remove intermediates immediately on success
    _cleanup_filter_alignment_tmp "${pth_in}" "${pth_out}"

    _finalize_alignment_filter "${threads}" "${fil_out}" "${chk_chr}" "${ref_fa}" \
        || return 1
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
