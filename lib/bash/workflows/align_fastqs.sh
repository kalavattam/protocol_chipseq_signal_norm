#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: align_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
#   GPT-5.6);
# - Anthropic Claude Code (Opus 5.5).
#
# Distributed under the MIT license.


# _print_dry_run_cmd
# _run_or_print_cmd
# _filter_bam_mapq_template
# _validate_args_align_fastqs
# _align_bowtie2
# _align_bwa_mem
# _align_bwa_aln
# _align_reads
# _sort_qname_bam
# _sort_coord_bam
# _check_optical_names
# _mark_dup_bam
# _build_align_record
# _print_optical_line
# _write_align_record
# _finalize_align_output
# align_fastqs


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
    _dir_src_aln="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    # shellcheck source=lib/bash/core/source_helpers.sh
    source "${_dir_src_aln}/../core/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_aln}/../core/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_aln}/.." \
        check_args \
        check_source \
        format_outputs \
        record_alignment \
        workflows/align_fastqs_separate || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_aln
}


# Print a command for a dry run: the named command arrays as a pipeline, one
# option per line with continuations, then an optional redirect of the last
# command's stdout, then a blank line.
function _print_dry_run_cmd() {
    local fil_red=""
    local idx_cmd=0
    local decl idx ind nam_cmd opt_val pfx tok
    local -a arr_lin
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _print_dry_run_cmd
    [--help] [--out <file>] arr_nam [arr_nam ...]

  Print a command to stdout for a dry run of 'align_fastqs()': each named command array in turn as a pipeline, then ' > file' when '--out' is given, then a blank line.

  The output is valid shell: words are shell-escaped, and each option or option-value pair is on its own continuation line.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  --out : file
    File that the last command's stdout is redirected to.

  1+  arr_nam : str
    Name of an indexed command array; several names print a pipeline.

Returns
-------
  Returns 0 after printing; prints an error and returns 1 if no array name is given or a name is not a nonempty indexed array.

Notes
-----
  Runtime requirements:
    bash >= 4.4

  - The program name, with its subcommand for 'samtools', 'bwa', and 'bwa-mem2', starts the command; 'mv', 'rm', 'cut', 'uniq', 'awk', and 'gzip' also keep their leading options there, as in 'mv -f'.
  - An option shares a line with its value only when it is in this function's table of value-taking options, kept per program in the source, since shape alone cannot tell 'samtools markdup -t in.bam' (a flag, then a file) from 'bwa aln -t 4' (an option and its value).
  - A new value-taking option belongs in that table; without it, the option and its value print on separate lines, which changes the layout but not the command.

Examples
--------
  1. Print one command whose path contains a space.
    '''bash
    arr_cmd=( samtools index -@ 2 'work dir/tiny_pe.bam' )
    _print_dry_run_cmd arr_cmd
    '''

  2. Print a pipeline that writes a file.
    '''bash
    arr_view=( samtools view -F 0x900 work/tiny_se.qnam.bam )
    arr_cut=( cut -f 1 )
    _print_dry_run_cmd --out work/tiny_se.names.txt arr_view arr_cut
    '''
EOM
    )

    if [[ "${1:-}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ "${1:-}" == "--out" ]]; then
        fil_red="${2:-}"
        shift 2
    fi

    if (( $# == 0 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'arr_nam', is missing."
        return 1
    fi

    for nam_cmd in "$@"; do
        decl="$(declare -p "${nam_cmd}" 2> /dev/null)" || decl=""

        if [[ "${decl}" != "declare -a"* ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "'${nam_cmd}' is not an indexed array."
            return 1
        fi

        local -n arr_prn_ref="${nam_cmd}"

        if (( ${#arr_prn_ref[@]} == 0 )); then
            echo_err_func "${FUNCNAME[0]}" \
                "array '${nam_cmd}' is empty."
            return 1
        fi

        # A piped command starts with '| ' and is indented one level more.
        ind=""
        pfx=""
        if (( idx_cmd > 0 )); then ind="    "; pfx="| "; fi

        # Start with the program and, where it has one, its subcommand.
        idx=1
        case "${arr_prn_ref[0]}" in
            samtools|bwa|bwa-mem2) idx=2 ;;
        esac

        # Options that take a value, by program.
        case "${arr_prn_ref[*]:0:idx}" in
            bowtie2)                  opt_val="-p --seed -x -I -X -1 -2 -U" ;;
            "samtools view")          opt_val="-@ -F -f -q -o -U -T -O -N" ;;
            "samtools markdup")       opt_val="-@ -d -f -O" ;;
            "samtools merge")         opt_val="-@ -o -O -s" ;;
            samtools\ *)              opt_val="-@ -o -O" ;;
            "bwa mem"|"bwa-mem2 mem") opt_val="-t -K" ;;
            "bwa aln")                opt_val="-t" ;;
            awk)                      opt_val="-v" ;;
            cut)                      opt_val="-f" ;;
            *)                        opt_val="" ;;
        esac

        # A short utility keeps its leading options, with their values, on the
        # first line, as in 'mv -f' or 'uniq -c'; only its operands go on
        # continuation lines.
        case "${arr_prn_ref[0]}" in
            mv|rm|cut|uniq|awk|gzip)
                while [[ "${arr_prn_ref[idx]:-}" == -* ]]; do
                    if [[ " ${opt_val} " == *" ${arr_prn_ref[idx]} "* ]]; then
                        idx=$(( idx + 1 ))
                    fi

                    idx=$(( idx + 1 ))
                done
                ;;
        esac

        printf -v tok '%q ' "${arr_prn_ref[@]:0:idx}"
        arr_lin+=( "${ind}${pfx}${tok% }" )

        while (( idx < ${#arr_prn_ref[@]} )); do
            printf -v tok '%q' "${arr_prn_ref[idx]}"

            if \
                [[ " ${opt_val} " == *" ${arr_prn_ref[idx]} "* ]] \
                    && (( idx + 1 < ${#arr_prn_ref[@]} ))
            then
                printf -v tok '%s %q' "${tok}" "${arr_prn_ref[idx+1]}"
                idx=$(( idx + 1 ))
            fi

            # An unquoted comma is literal in Bash, and '%q' escapes the braces
            # that could make it special, so print it bare, as in
            # '-O bam,level=1'; a quoted word is left as it is.
            if [[ "${tok}" != *"'"* ]]; then tok="${tok//\\,/,}"; fi

            arr_lin+=( "${ind}    ${tok}" )
            idx=$(( idx + 1 ))
        done

        unset -n arr_prn_ref
        idx_cmd=$(( idx_cmd + 1 ))
    done

    if [[ -n "${fil_red}" ]]; then
        printf -v tok '%q' "${fil_red}"
        arr_lin+=( "${ind}    > ${tok}" )
    fi

    # Continue every line but the last, then leave a blank line.
    for (( idx = 0; idx < ${#arr_lin[@]} - 1; idx++ )); do
        printf '%s \\\n' "${arr_lin[idx]}"
    done

    printf '%s\n\n' "${arr_lin[-1]}"
}


# Run one command array or, for a dry run, print it instead.
function _run_or_print_cmd() {
    local dry_run="${1:-}"
    local nam_cmd="${2:-}"
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _run_or_print_cmd
    [--help] dry_run arr_nam

  Run the command held in a named array or, when 'dry_run' is 'true', print it to stdout with '_print_dry_run_cmd' instead, so a dry run shows the same argument vector that a real run executes.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  dry_run : bool
    Whether to print the command instead of running it.

  2  arr_nam : str
    Name of the indexed command array.

Returns
-------
  Returns the command's exit status or 0 after printing it; prints an error and returns 1 if 'arr_nam' is missing.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Print an index command instead of running it.
    '''bash
    arr_cmd=( samtools index -@ 2 work/tiny_pe.bam )
    _run_or_print_cmd true arr_cmd
    '''

  2. Display the helper's argument contract.
    '''bash
    _run_or_print_cmd --help
    '''
EOM
    )

    if [[ "${dry_run}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${nam_cmd}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'arr_nam', is missing."
        return 1
    fi

    if [[ "${dry_run}" == "true" ]]; then
        _print_dry_run_cmd "${nam_cmd}"
        return
    fi

    local -n arr_run_ref="${nam_cmd}"
    "${arr_run_ref[@]}"
}


# Keep or drop each template (all records sharing a read name) of a queryname-
# grouped SAM or BAM file as a whole. A template is kept when its primary
# records pass a MAPQ threshold and, optionally, carry the proper-pair flag.
# Write the kept records to a BAM file and, optionally, the dropped ones to
# another.
function _filter_bam_mapq_template() {
    local bam_in="${1:-}"
    local bam_out="${2:-}"
    local mapq="${3:-}"
    local n_pri="${4:-}"
    local threads="${5:-1}"
    local func="${6:-${FUNCNAME[1]:-${FUNCNAME[0]}}}"
    local dry_run="${7:-false}"
    local req_flg="${8:-false}"
    local bam_rej="${9:-}"
    local txt_nam="${bam_out%.bam}.names.txt"
    local -a arr_view arr_cut arr_uniq arr_awk arr_cmd
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _filter_bam_mapq_template
    [--help] bam_in bam_out mapq n_pri [threads] [func] [dry_run] [req_flg] [bam_rej]

  Keep or drop each template (all records sharing a read name) of a queryname-grouped SAM or BAM file as a whole, by whether its primary records pass a MAPQ threshold. With 'req_flg', they must also carry the proper-pair flag (SAM flag 2). Write the kept records to a BAM file.

  Secondary and supplementary records follow their template, whatever their own MAPQ and flags.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  bam_in : file
    Queryname-grouped SAM or BAM file to filter.

  2  bam_out : file
    BAM file to write the kept records to, at compression level 1.

  3  mapq : int
    MAPQ threshold that every primary record of a kept template must meet.

  4  n_pri : int
    Number of primary records per template: 2 for paired-end data, 1 for single-end data.

  5  threads : int
    Number of threads to use for 'samtools view' (default: '${threads}').

  6  func : str
    Name of the calling function for diagnostics (default: the caller of this helper).

  7  dry_run : bool
    Whether to print the commands to stdout instead of running them (default: '${dry_run}').

  8  req_flg : bool
    Whether every primary record of a kept template must also carry the proper-pair flag (default: '${req_flg}').

  9  bam_rej : file
    BAM file to write the dropped templates to, at compression level 1, or '' to discard them (default: '').

Returns
-------
  Writes 'bam_out' (and 'bam_rej', if given) and returns 0 on success or prints the commands and returns 0 for a dry run; prints an error and returns 1 on failure.

Notes
-----
  Runtime requirements:
    - awk
    - bash >= 4.4
    - samtools

Examples
--------
  1. Keep paired-end templates whose two primary records reach MAPQ 30.
    '''bash
    _filter_bam_mapq_template \\
        sample.qnam.bam sample.mapq.bam 30 2 4
    '''

  2. Print, without running, the commands that keep single-end reads whose primary record reaches MAPQ 1.
    '''bash
    _filter_bam_mapq_template \\
        sample.qnam.bam sample.mapq.bam 1 1 1 align_fastqs true
    '''
EOM
    )

    if [[ "${bam_in}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        echo >&2
        return 0
    elif [[ -z "${bam_in}" || -z "${bam_out}" || -z "${mapq}" ]]; then
        echo_err_func "${func}" \
            "positional arguments 1 to 3, 'bam_in', 'bam_out', and 'mapq'," \
            "are required."
        return 1
    elif [[ -z "${n_pri}" ]]; then
        echo_err_func "${func}" \
            "positional argument 4, 'n_pri', is missing."
        return 1
    fi

    # List each read name that has exactly 'n_pri' aligned primary records at
    # or above the MAPQ threshold: 2 for paired-end data, 1 for single-end.
    #
    # '-F 2308' (4 + 256 + 2048) leaves out unaligned, secondary, and
    # supplementary records. So an unaligned record is never counted, even if
    # one reached this filter.
    #
    # With 'req_flg', '-f 2' counts only records that still carry the
    # proper-pair flag, so a pair that 'samtools fixmate' unflagged is dropped
    # whole. One such pair is a cross-reference pair that 'bwa sampe' flags
    # proper (https://github.com/lh3/bwa/issues/460).
    arr_view=( samtools view -@ "${threads}" -F 2308 )
    if [[ "${req_flg}" == "true" ]]; then arr_view+=( -f 2 ); fi
    arr_view+=( -q "${mapq}" "${bam_in}" )
    arr_cut=( cut -f 1 )
    arr_uniq=( uniq -c )
    # shellcheck disable=SC2016  # The program is for awk, not the shell.
    arr_awk=( awk -v "n=${n_pri}" '$1 == n { print $2 }' )

    if [[ "${dry_run}" == "true" ]]; then
        _print_dry_run_cmd \
            --out "${txt_nam}" arr_view arr_cut arr_uniq arr_awk
    elif ! (
        set -o pipefail
        "${arr_view[@]}" \
            | "${arr_cut[@]}" \
            | "${arr_uniq[@]}" \
            | "${arr_awk[@]}" \
            > "${txt_nam}"
    ); then
        rm -f "${txt_nam}"
        echo_err_func "${func}" \
            "failed to list read names whose primary records pass MAPQ" \
            "'${mapq}' in '${bam_in}'."
        return 1
    fi

    # Write every record of the listed names, secondary and supplementary
    # records included, at compression level 1, since 'align_fastqs()' sorts
    # the file again in Step 3. With 'bam_rej', write the other records there.
    arr_cmd=(
        samtools view
            -@ "${threads}"
            -O "bam,level=1"
            -N "${txt_nam}"
            -o "${bam_out}"
    )
    if [[ -n "${bam_rej}" ]]; then arr_cmd+=( -U "${bam_rej}" ); fi
    arr_cmd+=( "${bam_in}" )

    if ! \
        _run_or_print_cmd "${dry_run}" arr_cmd
    then
        rm -f "${txt_nam}" "${bam_out}" ${bam_rej:+"${bam_rej}"}
        echo_err_func "${func}" \
            "failed to filter '${bam_in}' by read name."
        return 1
    fi

    arr_cmd=( rm -f "${txt_nam}" )
    _run_or_print_cmd "${dry_run}" arr_cmd
}


# Validate the parsed arguments of 'align_fastqs()', in the order the checks
# have always run, refusing or warning per 'HELP.PARAMETER.APPLICABILITY'.
# TODO: "semantic paragraphs" for options under Usage in 'show_help'.
function _validate_args_align_fastqs() {
    local func="${1:-}"
    local threads="${2:-}"
    local aligner="${3:-}"
    local bt2_mode="${4:-}"
    local bt2_mode_set="${5:-}"
    local bt2_X="${6:-}"
    local bt2_X_set="${7:-}"
    local bt2_inherited="${8:-}"
    local bt2_relaxed="${9:-}"
    local bwa_alg="${10:-}"
    local bwa_alg_set="${11:-}"
    local mapq="${12:-}"
    local req_flg="${13:-}"
    local unal_mode="${14:-}"
    local index="${15:-}"
    local ref_fa="${16:-}"
    local fq_1="${17:-}"
    local fq_2="${18:-}"
    local fil_out="${19:-}"
    local out_fmt
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _validate_args_align_fastqs
    [--help] func threads aligner bt2_mode bt2_mode_set bt2_X bt2_X_set bt2_inherited bt2_relaxed bwa_alg bwa_alg_set mapq req_flg unal_mode index ref_fa fq_1 fq_2 fil_out

  Validate the parsed arguments of 'align_fastqs()' before alignment.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  01  func : str
    Name of the calling function for diagnostics.

  02  threads : int
    Number of threads to use.

  03  aligner : str
    Alignment program: 'bowtie2', 'bwa', or 'bwa-mem2'.

  04  bt2_mode : str
    Bowtie 2 alignment type: 'local', 'global', or 'end-to-end'.

  05  bt2_mode_set : bool
    Whether the caller supplied '--bt2_mode'.

  06  bt2_X : int
    Maximum fragment length for a Bowtie 2 paired-end alignment.

  07  bt2_X_set : bool
    Whether the caller supplied '--bt2_X'.

  08  bt2_inherited : bool
    Whether the caller asked for the earlier "inherited" Bowtie 2 paired-end call.

  09  bt2_relaxed : bool
    Whether the caller asked for the relaxed Bowtie 2 paired-end call.

  10  bwa_alg : str
    BWA algorithm: 'mem' or 'aln'.

  11  bwa_alg_set : bool
    Whether the caller supplied '--bwa_alg'.

  12  mapq : int
    MAPQ threshold.

  13  req_flg : bool
    Whether properly paired alignments are required.

  14  unal_mode : str
    What to do with unaligned reads and their mates: 'remove' or 'separate'.

  15  index : path
    Path to the aligner index/reference.

  16  ref_fa : file
    Reference FASTA path; required for CRAM output, otherwise ignored with a warning.

  17  fq_1 : file
    First FASTQ input file.

  18  fq_2 : file
    Second FASTQ input file, or empty for single-end data.

  19  fil_out : file
    Path to the final alignment file.

Returns
-------
  Returns 0 when validation succeeds; prints an error and returns 1 when an argument is missing, invalid, or refused.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - dirname

Examples
--------
  1. Validate a single-end Bowtie 2 call on the committed fixtures.
    '''bash
    _validate_args_align_fastqs \\
        align_fastqs \\
        1 \\
        bowtie2 \\
        end-to-end \\
        false \\
        1000 \\
        false \\
        false \\
        false \\
        mem \\
        false \\
        1 \\
        false \\
        remove \\
        tests/fixtures/align_fastqs/bowtie2/tiny \\
        '' \\
        tests/fixtures/align_fastqs/fastq/se/tiny_se.atria.fastq.gz \\
        '' \\
        work/tiny_se.bam
    '''

  2. Validate a paired-end BWA-MEM2 call that writes CRAM.
    '''bash
    _validate_args_align_fastqs \\
        align_fastqs \\
        4 \\
        bwa-mem2 \\
        end-to-end \\
        false \\
        1000 \\
        false \\
        false \\
        false \\
        mem \\
        false \\
        30 \\
        true \\
        remove \\
        tests/fixtures/align_fastqs/bwa-mem2/tiny.fa \\
        tests/fixtures/align_fastqs/reference/tiny.fa \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R1.atria.fastq.gz \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        work/tiny_pe.cram
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
            "'--threads' was assigned '${threads}' but must be a positive" \
            "integer greater than or equal to 1."
        return 1
    fi

    case "${aligner}" in
        bowtie2)
            case "${bt2_mode}" in
                local|global|end-to-end) : ;;
                *)
                    echo_err_func "${func}" \
                        "'--bt2_mode' must be 'local', 'global', or" \
                        "'end-to-end': '${bt2_mode}'."
                    return 1
                    ;;
            esac

            if [[ ! "${bt2_X}" =~ ^[1-9][0-9]*$ ]]; then
                echo_err_func "${func}" \
                    "'--bt2_X' was assigned '${bt2_X}' but must be a positive" \
                    "integer."
                return 1
            fi
            ;;

        bwa)
            case "${bwa_alg}" in
                mem|aln) : ;;
                *)
                    echo_err_func "${func}" \
                        "'--bwa_alg' must be 'mem' or 'aln': '${bwa_alg}'."
                    return 1
                    ;;
            esac
            ;;

        bwa-mem2) : ;;

        *)
            echo_err_func "${func}" \
                "'--aligner' must be 'bowtie2', 'bwa', or 'bwa-mem2':" \
                "'${aligner}'."
            return 1
            ;;
    esac

    if [[ ! "${mapq}" =~ ^[0-9]+$ ]]; then
        echo_err_func "${func}" \
            "--mapq was assigned '${mapq}' but must be an integer greater" \
            "than or equal to 0."
        return 1
    fi

    case "${unal_mode}" in
        remove|separate) : ;;
        *)
            echo_err_func "${func}" \
                "'--unal_mode' must be 'remove' or 'separate': '${unal_mode}'."
            return 1
            ;;
    esac

    if [[ -z "${index}" ]]; then
        echo_err_func "${func}" \
            "'--index' is a required argument."
        return 1
    fi

    if [[ ! -d "$(dirname "${index}")" ]]; then
        echo_err_func "${func}" \
            "directory associated with '--index' does not exist:" \
            "'$(dirname "${index}")'."
        return 1
    fi

    if [[ -z "${fil_out}" ]]; then
        echo_err_func "${func}" \
            "'--fil_out' is a required argument."
        return 1
    fi

    if [[ "${aligner}" =~ ^(bwa|bwa-mem2)$ && ! -f "${index}" ]]; then
        echo_err_func "${func}" \
            "file associated with '--index' does not exist: '${index}'."
        return 1
    fi

    if [[ ! -d "$(dirname "${fil_out}")" ]]; then
        echo_err_func "${func}" \
            "directory associated with '--fil_out' does not exist:" \
            "'$(dirname "${fil_out}")'"
        return 1
    fi

    case "${fil_out}" in
        *.bam)
            out_fmt="bam"
            ;;

        *.sam)
            echo_err_func "${func}" \
                "SAM output is not supported by 'align_fastqs()'. Please use" \
                "'.bam' or '.cram' for '--fil_out'."
            return 1
            ;;

        *.cram)
            out_fmt="cram"
            ;;

        *)
            echo_err_func "${func}" \
                "'--fil_out' must end in '.bam' or '.cram', but was assigned" \
                "'${fil_out}'."
            return 1
            ;;
    esac

    if [[ "${out_fmt}" == "cram" ]]; then
        if [[ -z "${ref_fa}" ]]; then
            echo_err_func "${func}" \
                "'--ref_fa' is required when '--fil_out' ends in '.cram'."
            return 1
        fi

        if [[ ! -f "${ref_fa}" ]]; then
            echo_err_func "${func}" \
                "file associated with '--ref_fa' does not exist: '${ref_fa}'."
            return 1
        fi
    elif [[ -n "${ref_fa}" ]]; then
        echo_warn_func "${func}" \
            "'--ref_fa' has no effect without CRAM output and is ignored."
    fi

    if [[ -z "${fq_1}" ]]; then
        echo_err_func "${func}" \
            "'--fq_1' is a required argument."
        return 1
    fi

    if [[ ! -f "${fq_1}" ]]; then
        echo_err_func "${func}" \
            "file associated with '--fq_1' does not exist: '${fq_1}'."
        return 1
    fi

    if [[ -n "${fq_2}" ]]; then
        if [[ ! -f "${fq_2}" ]]; then
            echo_err_func "${func}" \
                "file associated with '--fq_2' does not exist: '${fq_2}'."
            return 1
        fi
    fi

    # Per 'HELP.PARAMETER.APPLICABILITY', refuse the other aligner's option and
    # warn about a pair requirement with no pairs.
    if [[ "${aligner}" != "bowtie2" && "${bt2_mode_set}" == "true" ]]; then
        echo_err_func "${func}" \
            "'--bt2_mode' is for '--aligner bowtie2'; it has no effect with" \
            "'--aligner ${aligner}'."
        return 1
    fi

    if [[ "${aligner}" != "bowtie2" && "${bt2_X_set}" == "true" ]]; then
        echo_err_func "${func}" \
            "'--bt2_X' is for '--aligner bowtie2'; it has no effect with" \
            "'--aligner ${aligner}'."
        return 1
    fi

    if [[ "${aligner}" != "bowtie2" && "${bt2_inherited}" == "true" ]]; then
        echo_err_func "${func}" \
            "'--bt2_inherited' is for '--aligner bowtie2'; it has no effect" \
            "with '--aligner ${aligner}'."
        return 1
    fi

    if [[ "${aligner}" != "bowtie2" && "${bt2_relaxed}" == "true" ]]; then
        echo_err_func "${func}" \
            "'--bt2_relaxed' is for '--aligner bowtie2'; it has no effect" \
            "with '--aligner ${aligner}'."
        return 1
    fi

    # The inherited call passes no '-X', so '--bt2_X' would change nothing.
    if [[ "${bt2_X_set}" == "true" && "${bt2_inherited}" == "true" ]]; then
        echo_err_func "${func}" \
            "'--bt2_X' has no effect with '--bt2_inherited', which passes no" \
            "'-X'."
        return 1
    fi

    # Each makes its own paired-end call, so they cannot be combined.
    if [[ "${bt2_relaxed}" == "true" && "${bt2_inherited}" == "true" ]]; then
        echo_err_func "${func}" \
            "'--bt2_relaxed' cannot be combined with '--bt2_inherited'; each" \
            "makes its own Bowtie 2 paired-end call."
        return 1
    fi

    # The inherited call keeps '--no-unal', so Bowtie 2 writes no unaligned
    # reads to separate.
    if [[ "${unal_mode}" == "separate" && "${bt2_inherited}" == "true" ]]
    then
        echo_err_func "${func}" \
            "'--unal_mode separate' has no effect with '--bt2_inherited'," \
            "whose call writes no unaligned reads."
        return 1
    fi

    # Relaxed pairing is a paired-end call; single-end input has no pairs.
    if [[ -z "${fq_2}" && "${bt2_relaxed}" == "true" ]]; then
        echo_err_func "${func}" \
            "'--bt2_relaxed' is for paired-end input; it has no effect with" \
            "single-end input."
        return 1
    fi

    # 'bwa-mem2' has only 'mem', so asking for it changes nothing.
    if [[ "${aligner}" != "bwa" && "${bwa_alg_set}" == "true" ]]; then
        if [[ "${aligner}" == "bwa-mem2" && "${bwa_alg}" == "mem" ]]; then
            echo_warn_func "${func}" \
                "'--bwa_alg' has no effect with '--aligner bwa-mem2' and is" \
                "ignored."
        else
            echo_err_func "${func}" \
                "'--bwa_alg' is for '--aligner bwa'; it has no effect with" \
                "'--aligner ${aligner}'."
            return 1
        fi
    fi

    if [[ -z "${fq_2}" && "${req_flg}" == "true" ]]; then
        echo_warn_func "${func}" \
            "'--req_flg' has no effect with single-end input and is ignored."
    fi

    if [[ -z "${fq_2}" && "${bt2_X_set}" == "true" ]]; then
        echo_warn_func "${func}" \
            "'--bt2_X' has no effect with single-end input and is ignored."
    fi

    if [[ -z "${fq_2}" && "${bt2_inherited}" == "true" ]]; then
        echo_warn_func "${func}" \
            "'--bt2_inherited' has no effect with single-end input and is" \
            "ignored."
    fi

    # Per 'HELP.ALIAS.PUBLIC', a hidden measurement setting warns on every run,
    # and so does a side file that is not for downstream use.
    if [[ "${bt2_relaxed}" == "true" ]]; then
        echo_warn_func "${func}" \
            "'--bt2_relaxed' makes the relaxed Bowtie 2 paired-end call; its" \
            "output measures pairing and is not for signal or scaling."
    fi

    if [[ "${unal_mode}" == "separate" ]]; then
        echo_warn_func "${func}" \
            "'--unal_mode separate' writes the records that the alignment" \
            "file leaves out to diagnostic files, which are not for signal" \
            "or scaling."
    fi
}


# Step 1 with Bowtie 2: align paired- or single-end reads and pipe them to the
# Samtools call that the caller built.
function _align_bowtie2() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local bt2_mode="${4:-}"
    local bt2_X="${5:-}"
    local bt2_inherited="${6:-}"
    local bt2_relaxed="${7:-}"
    local unal_mode="${8:-}"
    local index="${9:-}"
    local fq_1="${10:-}"
    local fq_2="${11:-}"
    local -a arr_bt2 arr_sam
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _align_bowtie2
    [--help] func dry_run threads bt2_mode bt2_X bt2_inherited bt2_relaxed unal_mode index fq_1 fq_2 [samtools_arg ...]

  Align paired- or single-end reads with Bowtie 2 and pipe the alignments to 'samtools view' with the caller's arguments (Step 1 of 'align_fastqs()').

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  01  func : str
    Name of the calling function for diagnostics.

  02  dry_run : bool
    Whether to print the command to stdout instead of running it.

  03  threads : int
    Number of threads to use.

  04  bt2_mode : str
    Bowtie 2 alignment type: 'local', 'global', or 'end-to-end'.

  05  bt2_X : int
    Maximum fragment length for a paired-end alignment, passed as '-X'.

  06  bt2_inherited : bool
    Whether to make the earlier paired-end call ('--no-overlap --no-dovetail', no '-X') instead.

  07  bt2_relaxed : bool
    Whether to make the relaxed paired-end call ('-I 0 -X bt2_X', without '--no-mixed' and without '--no-discordant') instead.

  08  unal_mode : str
    What happens to unaligned reads: 'remove' passes '--no-unal', so Bowtie 2 writes none; 'separate' does not.

  09  index : path
    Bowtie 2 index stem.

  10  fq_1 : file
    First FASTQ input file.

  11  fq_2 : file
    Second FASTQ input file, or empty for single-end data.

  12+  samtools_arg : str
    Arguments for 'samtools view', ending with its output path.

Returns
-------
  Returns 0 when alignment and filtering succeed or after printing the command for a dry run; prints an error and returns 1 otherwise.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - bowtie2
    - samtools

Examples
--------
  1. Align the single-end fixture and write the filtered BAM work file.
    '''bash
    _align_bowtie2 \\
        align_fastqs \\
        false \\
        1 \\
        end-to-end \\
        1000 \\
        false \\
        false \\
        remove \\
        tests/fixtures/align_fastqs/bowtie2/tiny \\
        tests/fixtures/align_fastqs/fastq/se/tiny_se.atria.fastq.gz \\
        '' \\
        -@ 1 -F 4 -o work/tiny_se.bam
    '''

  2. Print, without running, the paired-end call in local mode that requires proper pairs.
    '''bash
    _align_bowtie2 \\
        align_fastqs \\
        true \\
        2 \\
        local \\
        1000 \\
        false \\
        false \\
        remove \\
        tests/fixtures/align_fastqs/bowtie2/tiny \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R1.atria.fastq.gz \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        -@ 2 -F 12 -f 2 -o work/tiny_pe.bam
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    # Assign Bowtie 2 arguments. '--reorder' keeps input order at any thread
    # count. '--seed 0' is the project's seed and Bowtie 2's default, so new
    # files have the same records as earlier ones.
    arr_bt2=( bowtie2 -p "${threads}" --reorder --seed 0 -x "${index}" )

    if [[ "${bt2_mode}" == "local" ]]; then
        arr_bt2+=( --very-sensitive-local )
    elif [[ "${bt2_mode}" =~ ^(end-to-end|global)$ ]]; then
        arr_bt2+=( --very-sensitive )
    fi

    # Exclude unaligned reads from output, unless they are kept in a side file.
    if [[ "${unal_mode}" == "remove" ]]; then arr_bt2+=( --no-unal ); fi
    arr_bt2+=( --phred33 )

    # Run Bowtie 2 with relaxed paired-end (top), paired-end (middle), or
    # single-end (bottom) sequenced reads.
    if [[ -n "${fq_2}" && "${bt2_relaxed}" == "true" ]]; then
        # Hidden per 'HELP.ALIAS.PUBLIC'. The relaxed call drops the pairing
        # options/filters, so non-standard pairs are written rather than lost,
        # and what the options would otherwise remove can now be measured. '-X'
        # sets the boundary between concordant and discordant pairs.
        arr_bt2+=( -I 0 -X "${bt2_X}" )
        arr_bt2+=( -1 "${fq_1}" -2 "${fq_2}" )
    elif [[ -n "${fq_2}" ]]; then
        arr_bt2+=( --no-mixed --no-discordant )

        # Keep overlapping mates and the bound fragment length unless asked for
        # for the earlier "inherited" call, which drops overlapping mates and
        # leaves '-X' at Bowtie 2's default of 500.
        if [[ "${bt2_inherited}" == "true" ]]; then
            arr_bt2+=( --no-overlap --no-dovetail )
        else
            arr_bt2+=( -X "${bt2_X}" )
        fi

        arr_bt2+=( -1 "${fq_1}" -2 "${fq_2}" )
    else
        arr_bt2+=( -U "${fq_1}" )
    fi

    arr_sam=( samtools view "${@:12}" )

    if [[ "${dry_run}" == "true" ]]; then
        _print_dry_run_cmd arr_bt2 arr_sam
        return
    fi

    if ! (
        set -o pipefail
        "${arr_bt2[@]}" | "${arr_sam[@]}"
    ); then
        if [[ -n "${fq_2}" ]]; then
            echo_err_func "${func}" \
                "Step #1: failed to align paired-end data with 'bowtie2'" \
                "('${bt2_mode}'), and/or failed to process the alignments" \
                "with 'samtools view'."
            return 1
        else
            echo_err_func "${func}" \
                "Step #1: failed to align single-end data with 'bowtie2'" \
                "('${bt2_mode}'), and/or failed to process the alignments" \
                "with 'samtools view'."
            return 1
        fi
    fi
}


# Step 1 with BWA-MEM or BWA-MEM2: align paired- or single-end reads and pipe
# them to the Samtools call that the caller built, which excludes unaligned
# reads.
function _align_bwa_mem() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local prog="${4:-}"
    local index="${5:-}"
    local fq_1="${6:-}"
    local fq_2="${7:-}"
    local typ_dat
    local -a arr_bwa arr_sam
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _align_bwa_mem
    [--help] func dry_run threads prog index fq_1 fq_2 [samtools_arg ...]

  Align paired- or single-end reads with 'bwa mem' or 'bwa-mem2 mem' and pipe the alignments to 'samtools view' with the caller's arguments (Step 1 of 'align_fastqs()').

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  dry_run : bool
    Whether to print the command to stdout instead of running it.

  3  threads : int
    Number of threads to use.

  4  prog : str
    Aligner program: 'bwa' or 'bwa-mem2'.

  5  index : file
    Indexed reference FASTA path.

  6  fq_1 : file
    First FASTQ input file.

  7  fq_2 : file
    Second FASTQ input file, or empty for single-end data.

  8+  samtools_arg : str
    Arguments for 'samtools view', ending with its output path.

Returns
-------
  Returns 0 when alignment and filtering succeed or after printing the command for a dry run; prints an error and returns 1 otherwise.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - bwa (when 'prog' is 'bwa')
    - bwa-mem2 (when 'prog' is 'bwa-mem2')
    - samtools

Examples
--------
  1. Align the paired-end fixtures with BWA-MEM2 and write the BAM work file.
    '''bash
    _align_bwa_mem \\
        align_fastqs \\
        false \\
        2 \\
        bwa-mem2 \\
        tests/fixtures/align_fastqs/bwa-mem2/tiny.fa \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R1.atria.fastq.gz \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        -@ 2 -F 12 -o work/tiny_pe.bam
    '''

  2. Print, without running, the single-end call with BWA-MEM.
    '''bash
    _align_bwa_mem \\
        align_fastqs \\
        true \\
        1 \\
        bwa \\
        tests/fixtures/align_fastqs/bwa/tiny.fa \\
        tests/fixtures/align_fastqs/fastq/se/tiny_se.atria.fastq.gz \\
        '' \\
        -@ 1 -F 4 -o work/tiny_se.bam
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    # Assign BWA or BWA-MEM2 arguments. By default, BWA estimates PE insert
    # sizes per batch of input bases, and the batch grows with the thread
    # count, so output changes with '-t'. A fixed '-K' sets the batch, making
    # outputs the same at any thread count.
    arr_bwa=(
        "${prog}" mem -t "${threads}" -K 100000000 "${index}" "${fq_1}"
    )

    if [[ -n "${fq_2}" ]]; then
        arr_bwa+=( "${fq_2}" )
    fi

    arr_sam=( samtools view "${@:8}" )

    if [[ "${dry_run}" == "true" ]]; then
        _print_dry_run_cmd arr_bwa arr_sam
        return
    fi

    if [[ -n "${fq_2}" ]]; then
        typ_dat="paired-end"
    else
        typ_dat="single-end"
    fi

    if ! (
        set -o pipefail
        "${arr_bwa[@]}" | "${arr_sam[@]}"
    ); then
        echo_err_func "${func}" \
            "Step #1: failed to align ${typ_dat} data with '${prog} mem'," \
            "and/or failed to process the alignments with 'samtools view'."
        return 1
    fi
}


# Step 1 with BWA-backtrack: write '.sai' files with 'bwa aln', pair them with
# 'bwa sampe' or 'bwa samse', and pipe the alignments to the Samtools call that
# the caller built, which excludes unaligned reads.
function _align_bwa_aln() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local index="${4:-}"
    local fq_1="${5:-}"
    local fq_2="${6:-}"
    local fil_wrk="${7:-}"
    local fil_sai_1 fil_sai_2
    local -a arr_aln arr_sam arr_bwa_sam arr_cmd
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _align_bwa_aln
    [--help] func dry_run threads index fq_1 fq_2 fil_wrk [samtools_arg ...]

  Align paired- or single-end reads with 'bwa aln' and 'bwa sampe' or 'bwa samse', and pipe the alignments to 'samtools view' with the caller's arguments (Step 1 of 'align_fastqs()').

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  dry_run : bool
    Whether to print the commands to stdout instead of running them.

  3  threads : int
    Number of threads to use for 'bwa aln'.

  4  index : file
    Indexed reference FASTA path.

  5  fq_1 : file
    First FASTQ input file.

  6  fq_2 : file
    Second FASTQ input file, or empty for single-end data.

  7  fil_wrk : file
    BAM work file; the temporary '.sai' files are written beside it and deleted.

  8+  samtools_arg : str
    Arguments for 'samtools view', ending with its output path.

Returns
-------
  Returns 0 when alignment and filtering succeed or after printing the commands for a dry run; prints an error and returns 1 otherwise.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - bwa
    - rm
    - samtools

Examples
--------
  1. Align the single-end fixture with BWA-backtrack and write the BAM work file.
    '''bash
    _align_bwa_aln \\
        align_fastqs \\
        false \\
        1 \\
        tests/fixtures/align_fastqs/bwa/tiny.fa \\
        tests/fixtures/align_fastqs/fastq/se/tiny_se.atria.fastq.gz \\
        '' \\
        work/tiny_se.bam \\
        -@ 1 -F 4 -o work/tiny_se.bam
    '''

  2. Print, without running, the paired-end calls with BWA-backtrack.
    '''bash
    _align_bwa_aln \\
        align_fastqs \\
        true \\
        2 \\
        tests/fixtures/align_fastqs/bwa/tiny.fa \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R1.atria.fastq.gz \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        work/tiny_pe.bam \\
        -@ 2 -F 12 -o work/tiny_pe.bam
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    # Assign BWA-backtrack arguments.
    fil_sai_1="${fil_wrk%.bam}.R1.sai"
    fil_sai_2=""
    arr_sam=( samtools view "${@:8}" )

    arr_aln=( bwa aln -t "${threads}" "${index}" "${fq_1}" )

    if [[ "${dry_run}" == "true" ]]; then
        _print_dry_run_cmd --out "${fil_sai_1}" arr_aln
    elif ! \
        "${arr_aln[@]}" > "${fil_sai_1}"
    then
        echo_err_func "${func}" \
            "Step #1: failed to generate first '.sai' file with 'bwa aln'."
        return 1
    fi

    if [[ -n "${fq_2}" ]]; then
        fil_sai_2="${fil_wrk%.bam}.R2.sai"
        arr_aln=( bwa aln -t "${threads}" "${index}" "${fq_2}" )

        if [[ "${dry_run}" == "true" ]]; then
            _print_dry_run_cmd --out "${fil_sai_2}" arr_aln
        elif ! \
            "${arr_aln[@]}" > "${fil_sai_2}"
        then
            echo_err_func "${func}" \
                "Step #1: failed to generate second '.sai' file with" \
                "'bwa aln'."
            rm -f "${fil_sai_1}"
            return 1
        fi

        arr_bwa_sam=(
            bwa sampe
                "${index}" "${fil_sai_1}" "${fil_sai_2}"
                "${fq_1}" "${fq_2}"
        )

        if [[ "${dry_run}" == "true" ]]; then
            _print_dry_run_cmd arr_bwa_sam arr_sam
        elif ! (
            set -o pipefail
            "${arr_bwa_sam[@]}" | "${arr_sam[@]}"
        ); then
            rm -f "${fil_sai_1}" "${fil_sai_2}"
            echo_err_func "${func}" \
                "Step #1: failed to align paired-end data with 'bwa aln' or" \
                "'bwa sampe', and/or failed to process the alignments with" \
                "'samtools view'."
            return 1
        fi

        arr_cmd=( rm -f "${fil_sai_1}" "${fil_sai_2}" )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #1: failed to delete temporary '.sai' files."
            return 1
        fi

    else
        arr_bwa_sam=( bwa samse "${index}" "${fil_sai_1}" "${fq_1}" )

        if [[ "${dry_run}" == "true" ]]; then
            _print_dry_run_cmd arr_bwa_sam arr_sam
        elif ! (
            set -o pipefail
            "${arr_bwa_sam[@]}" | "${arr_sam[@]}"
        ); then
            rm -f "${fil_sai_1}"
            echo_err_func "${func}" \
                "Step #1: failed to align single-end data with 'bwa aln' or" \
                "'bwa samse', and/or failed to process the alignments with" \
                "'samtools view'."
            return 1
        fi

        arr_cmd=( rm -f "${fil_sai_1}" )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #1: failed to delete temporary '.sai' file."
            return 1
        fi
    fi
}


# Step 1: build the Samtools call that drops unaligned reads, then align with
# the selected aligner.
# TODO: "semantic paragraphs" for options under Usage in 'show_help'.
function _align_reads() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local aligner="${4:-}"
    local bt2_mode="${5:-}"
    local bt2_X="${6:-}"
    local bt2_inherited="${7:-}"
    local bt2_relaxed="${8:-}"
    local bwa_alg="${9:-}"
    local req_flg="${10:-}"
    local unal_mode="${11:-}"
    local index="${12:-}"
    local fq_1="${13:-}"
    local fq_2="${14:-}"
    local fil_wrk="${15:-}"
    local -a arr_sam
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _align_reads
    [--help] func dry_run threads aligner bt2_mode bt2_X bt2_inherited bt2_relaxed bwa_alg req_flg unal_mode index fq_1 fq_2 fil_wrk

  Build the 'samtools view' arguments that drop unaligned reads, then align with the selected aligner and write the BAM work file (Step 1 of 'align_fastqs()').

  With 'unal_mode' 'separate', the dropped records are written instead to the path of 'fil_wrk' with '.setU.bam' for '.bam'.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  01  func : str
    Name of the calling function for diagnostics.

  02  dry_run : bool
    Whether to print the commands to stdout instead of running them.

  03  threads : int
    Number of threads to use.

  04  aligner : str
    Alignment program: 'bowtie2', 'bwa', or 'bwa-mem2'.

  05  bt2_mode : str
    Bowtie 2 alignment type: 'local', 'global', or 'end-to-end'.

  06  bt2_X : int
    Maximum fragment length for a Bowtie 2 paired-end alignment.

  07  bt2_inherited : bool
    Whether to make the earlier "inherited" Bowtie 2 paired-end call.

  08  bt2_relaxed : bool
    Whether to make the relaxed Bowtie 2 paired-end call.

  09  bwa_alg : str
    BWA algorithm: 'mem' or 'aln'.

  10  req_flg : bool
    Whether properly paired alignments are required.

  11  unal_mode : str
    What to do with unaligned reads and their mates: 'remove' or 'separate'.

  12  index : path
    Path to the aligner index/reference.

  13  fq_1 : file
    First FASTQ input file.

  14  fq_2 : file
    Second FASTQ input file, or empty for single-end data.

  15  fil_wrk : file
    BAM work file to write.

Returns
-------
  Returns 0 when alignment succeeds or after printing the commands for a dry run; prints an error and returns 1 otherwise.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - bowtie2 (when 'aligner' is 'bowtie2')
    - bwa (when 'aligner' is 'bwa')
    - bwa-mem2 (when 'aligner' is 'bwa-mem2')
    - samtools

Examples
--------
  1. Align the paired-end fixtures with Bowtie 2, requiring proper pairs.
    '''bash
    _align_reads \\
        align_fastqs \\
        false \\
        2 \\
        bowtie2 \\
        end-to-end \\
        1000 \\
        false \\
        false \\
        mem \\
        true \\
        remove \\
        tests/fixtures/align_fastqs/bowtie2/tiny \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R1.atria.fastq.gz \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        work/tiny_pe.bam
    '''

  2. Print without running the single-end calls with BWA-backtrack, writing the unaligned reads to a separate work file.
    '''bash
    _align_reads \\
        align_fastqs \\
        true \\
        1 \\
        bwa \\
        end-to-end \\
        1000 \\
        false \\
        false \\
        aln \\
        false \\
        separate \\
        tests/fixtures/align_fastqs/bwa/tiny.fa \\
        tests/fixtures/align_fastqs/fastq/se/tiny_se.atria.fastq.gz \\
        '' \\
        work/tiny_se.bam
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ "${dry_run}" == "true" ]]; then
        printf '\n%s\n' \
            "# Step 1: align the reads into the BAM work file."
    fi

    # Drop unaligned reads with Samtools.
    #
    # In paired-end data, '-F 12' (0xC) drops a read whose mate is flagged
    # unaligned. Samtools tests each record alone, so this keeps pairs whole
    # only when both records' flags agree. BWA-backtrack, however, can flag one
    # mate unaligned and leave its partner's 0x8 clear. So, Step 2 drops any
    # template left with one primary record.
    arr_sam=( -@ "${threads}" )
    if [[ -n "${fq_2}" ]]; then
        arr_sam+=( -F 12 )
    else
        arr_sam+=( -F 4 )
    fi

    # With '--unal_mode separate', write the dropped records ("Set U") to their
    # own work file for the side file. Proper pairs are then required only in
    # Step 2, after the templates split between the two files move whole into
    # Set U. Otherwise, require them here, as before.
    if [[ "${unal_mode}" == "separate" ]]; then
        arr_sam+=( -U "${fil_wrk%.bam}.setU.bam" )
    elif [[ -n "${fq_2}" && "${req_flg}" == "true" ]]; then
        arr_sam+=( -f 2 )
    fi

    # The MAPQ threshold is not applied here, as a per-record '-q' can split a
    # template whenever its records' MAPQs differ. Instead, Step 2 applies it
    # by template for every aligner.
    #
    # Step 2 rewrites the work file and, with '--unal_mode separate', Set U.
    # Both are written at compression level 1, which is faster to write than
    # the default level, 6.
    arr_sam+=( -O "bam,level=1" -o "${fil_wrk}" )

    # Based on parsed arguments, run Bowtie 2, BWA, or BWA-MEM2 with output
    # piped to the above-constructed Samtools call.
    if [[ "${aligner}" == "bowtie2" ]]; then
        _align_bowtie2 \
            "${func}" "${dry_run}" "${threads}" "${bt2_mode}" "${bt2_X}" \
            "${bt2_inherited}" "${bt2_relaxed}" "${unal_mode}" "${index}" \
            "${fq_1}" "${fq_2}" "${arr_sam[@]}" \
            || return 1
    elif [[ "${aligner}" == "bwa" ]]; then
        if [[ "${bwa_alg}" == "mem" ]]; then
            _align_bwa_mem \
                "${func}" "${dry_run}" "${threads}" bwa "${index}" \
                "${fq_1}" "${fq_2}" "${arr_sam[@]}" \
                || return 1
        elif [[ "${bwa_alg}" == "aln" ]]; then
            _align_bwa_aln \
                "${func}" "${dry_run}" "${threads}" "${index}" \
                "${fq_1}" "${fq_2}" "${fil_wrk}" "${arr_sam[@]}" \
                || return 1
        fi
    elif [[ "${aligner}" == "bwa-mem2" ]]; then
        _align_bwa_mem \
            "${func}" "${dry_run}" "${threads}" bwa-mem2 "${index}" \
            "${fq_1}" "${fq_2}" "${arr_sam[@]}" \
            || return 1
    fi
}


# Step 2: sort the BAM work file by queryname, fix mates if the file holds
# paired-end alignments, and keep or drop each template whole by its primary
# records' MAPQ values.
function _sort_qname_bam() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local aligner="${4:-}"
    local mapq="${5:-}"
    local req_flg="${6:-}"
    local unal_mode="${7:-}"
    local fq_2="${8:-}"
    local fil_wrk="${9:-}"
    local bam_qnam="${10:-}"
    local bam_mate bam_mapq bam_set_u bam_rej bam_rej0 n_pri req_tpl
    local -a arr_cmd arr_srt arr_fix
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _sort_qname_bam
    [--help] func dry_run threads aligner mapq req_flg unal_mode fq_2 fil_wrk bam_qnam

  Sort the BAM work file by queryname into 'bam_qnam' and fix mates for paired-end data (Step 2 of 'align_fastqs()').

  Then, for every aligner and layout, keep or drop each template whole by whether its primary records pass 'mapq'. For paired-end data with 'req_flg', they must also still carry the proper-pair flag after 'fixmate'. A paired-end template that lost a mate in Step 1 is dropped too.

  The sort streams into 'fixmate' without a file between them. Every file written here is at compression level 1, making it faster for Step 3 to sort the results again.

  With 'unal_mode=separate', first move every template in the work file ("Set M") that has records in Step 1's unaligned work file ("Set U") wholly into Set U. Fix mates in Set U on its own. Then, with 'req_flg', require proper pairs in the rest as Step 1 does under 'remove', and keep the records it rejects.

  Last, sort those records and the templates the template filter drops into an improper and a low-MAPQ work file, so nothing is discarded. Steps 3 and 4 then write them out beside 'bam_qnam'.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  01  func : str
    Name of the calling function for diagnostics.

  02  dry_run : bool
    Whether to print the commands to stdout instead of running them.

  03  threads : int
    Number of threads to use.

  04  aligner : str
    Alignment program that wrote the BAM work file.

  05  mapq : int
    MAPQ threshold that a template's primary records must pass.

  06  req_flg : bool
    Whether properly paired alignments are required. The template filter applies it after 'fixmate'. Step 1 applies it before 'fixmate' under 'remove'; with 'unal_mode=separate', this step does so instead.

  07  unal_mode : str
    What to do with unaligned reads and their mates: 'remove' or 'separate'.

  08  fq_2 : file
    Second FASTQ input file, or empty for single-end data; selects mate fixing.

  09  fil_wrk : file
    BAM work file from Step 1.

  10  bam_qnam : file
    Queryname-sorted BAM file to write.

Returns
-------
  Returns 0 when sorting, mate fixing, and filtering succeed or after printing the commands for a dry run; prints an error and returns 1 otherwise.

Notes
-----
  Runtime requirements:
    - awk
    - bash >= 4.4
    - cut
    - mv
    - samtools
    - uniq

Examples
--------
  1. Queryname-sort and mate-fix a paired-end BWA work file, filtering templates at MAPQ 1.
    '''bash
    _sort_qname_bam \\
        align_fastqs \\
        false \\
        2 \\
        bwa \\
        1 \\
        false \\
        remove \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        work/tiny_pe.bam \\
        work/tiny_pe.qnam.bam
    '''

  2. Print, without running, the queryname sort and template filter of a single-end Bowtie 2 work file whose unaligned reads go to a side file.
    '''bash
    _sort_qname_bam \\
        align_fastqs \\
        true \\
        1 \\
        bowtie2 \\
        1 \\
        false \\
        separate \\
        '' \\
        work/tiny_se.bam \\
        work/tiny_se.qnam.bam
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ "${dry_run}" == "true" ]]; then
        printf '\n%s %s\n%s\n' \
            "# Step 2: sort the BAM work file by queryname," \
            "fix mates in paired-end data," \
            "# and filter templates."
    fi

    # A dry run writes no Step 1 file, so it prints the commands as they would
    # run.
    if [[ "${dry_run}" == "true" || -f "${fil_wrk}" ]]; then
        # With '--unal_mode separate', move templates split between the two
        # work files into Set U; for PE data, this also sorts Set M.
        if [[ "${unal_mode}" == "separate" ]]; then
            bam_set_u="${fil_wrk%.bam}.setU.bam"

            if ! \
                _route_unal_bam \
                    "${func}" \
                    "${dry_run}" \
                    "${threads}" \
                    "${fq_2}" \
                    "${fil_wrk}" \
                    "${bam_qnam}" \
                    "${bam_set_u}"
            then
                return 1
            fi
        fi

        # Sort by queryname and, for PE data, fix mates in the same stream. The
        # below template filter reads the result twice, so it is written at
        # compression level 1.
        arr_srt=(
            samtools sort
                -@ "${threads}"
                -n
        )
        arr_fix=(
            samtools fixmate
                -@ "${threads}"
                -c
                -m
                -O "bam,level=1"
        )

        if [[ -z "${fq_2}" ]]; then
            arr_srt+=( -O "bam,level=1" -o "${bam_qnam}" "${fil_wrk}" )

            if ! \
                _run_or_print_cmd "${dry_run}" arr_srt
            then
                echo_err_func "${func}" \
                    "Step #2: failed to queryname-sort ${aligner}-aligned" \
                    "BAM file."
                return 1
            fi
        elif [[ "${unal_mode}" == "remove" ]]; then
            arr_srt+=( -u -o - "${fil_wrk}" )
            arr_fix+=( - "${bam_qnam}" )

            if [[ "${dry_run}" == "true" ]]; then
                _print_dry_run_cmd arr_srt arr_fix
            elif ! (
                set -o pipefail
                "${arr_srt[@]}" | "${arr_fix[@]}"
            ); then
                echo_err_func "${func}" \
                    "Step #2: failed to queryname-sort ${aligner}-aligned" \
                    "BAM file and fix its mate pairs."
                return 1
            fi
        else
            # With '--unal_mode separate', Set M is already sorted. With '--
            # req_flg', keep the records the aligner flagged proper, then fix
            # mates. '-f 2' runs before 'fixmate', as in Step 1, so the
            # aligner's own 0x2 decides. The records it rejects go to their own
            # work file ("Set I0") and get their mates fixed there, so Markdup
            # 2 can judge them.
            bam_mate="${fil_wrk%.bam}.mate.bam"

            if [[ "${req_flg}" == "true" ]]; then
                bam_rej0="${fil_wrk%.bam}.setI0.bam"
                arr_cmd=(
                    samtools view
                        -@ "${threads}"
                        -u
                        -f 2
                        -U "${bam_rej0}"
                        "${bam_qnam}"
                )
                arr_fix+=( - "${bam_mate}" )

                if [[ "${dry_run}" == "true" ]]; then
                    _print_dry_run_cmd arr_cmd arr_fix
                elif ! (
                    set -o pipefail
                    "${arr_cmd[@]}" | "${arr_fix[@]}"
                ); then
                    echo_err_func "${func}" \
                        "Step #2: failed to keep properly paired alignments" \
                        "and fix their mate pairs."
                    return 1
                fi
            else
                arr_fix+=( "${bam_qnam}" "${bam_mate}" )

                if ! \
                    _run_or_print_cmd "${dry_run}" arr_fix
                then
                    echo_err_func "${func}" \
                        "Step #2: failed to fix mate pairs in" \
                        "queryname-sorted BAM file."
                    return 1
                fi
            fi

            arr_cmd=( mv -f "${bam_mate}" "${bam_qnam}" )

            if ! \
                _run_or_print_cmd "${dry_run}" arr_cmd
            then
                echo_err_func "${func}" \
                    "Step #2: failed to rename mate pair-fixed," \
                    "queryname-sorted BAM file."
                return 1
            fi
        fi

        # Keep or drop all records of each template (i.e., those sharing its
        # read name) by whether its primary records pass '--mapq'. In PE data,
        # both mates must pass. A per-record '-q' would split a template whose
        # records' MAPQs differ, such as BWA mates, supplementary records, and
        # Bowtie 2 secondaries under '-k'. Step 1 can also leave a PE template
        # with one primary record (#75). This filter runs for every call.
        #
        # With '--req_flg', the primary records must also still carry 0x2,
        # which 'fixmate' clears on a pair whose mates map to different
        # references or align in the same orientation (both forward or both
        # reverse). With '--unal_mode separate', the dropped templates go to a
        # work file of their own ("Set R").
        bam_mapq="${fil_wrk%.bam}.mapq.bam"
        if [[ -n "${fq_2}" ]]; then n_pri=2; else n_pri=1; fi
        req_tpl=false
        if [[ -n "${fq_2}" && "${req_flg}" == "true" ]]; then req_tpl=true; fi
        bam_rej=""
        if [[ "${unal_mode}" == "separate" ]]; then
            bam_rej="${fil_wrk%.bam}.setR.bam"
        fi

        if ! \
            _filter_bam_mapq_template \
                "${bam_qnam}" \
                "${bam_mapq}" \
                "${mapq}" \
                "${n_pri}" \
                "${threads}" \
                "${func}" \
                "${dry_run}" \
                "${req_tpl}" \
                "${bam_rej}"
        then
            echo_err_func "${func}" \
                "Step #2: failed to filter queryname-sorted BAM file by" \
                "template MAPQ."
            return 1
        fi

        arr_cmd=( mv -f "${bam_mapq}" "${bam_qnam}" )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #2: failed to rename MAPQ-filtered, queryname-sorted" \
                "BAM file."
            return 1
        fi

        # With '--unal_mode separate', sort the dropped templates into the
        # improper pairs and the low-MAPQ templates, then count what each work
        # file holds for '.diag.tsv.gz'.
        if [[ "${unal_mode}" == "separate" ]]; then
            if ! \
                _route_filt_bam \
                    "${func}" \
                    "${dry_run}" \
                    "${threads}" \
                    "${req_tpl}" \
                    "${fq_2}" \
                    "${fil_wrk}"
            then
                return 1
            fi
        fi

    else
        echo_err_func "${func}" \
            "Step #2: BAM work file from ${aligner} alignment does not exist."
        return 1
    fi
}


# Step 3: sort the queryname-sorted, fixmate-adjusted BAM file by coordinates
# into the BAM work file.
function _sort_coord_bam() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local fil_wrk="${4:-}"
    local bam_qnam="${5:-}"
    local bam_coor
    local -a arr_cmd
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _sort_coord_bam
    [--help] func dry_run threads fil_wrk bam_qnam

  Sort the queryname-sorted BAM file by coordinates into the BAM work file, deleting the queryname-sorted file (Step 3 of 'align_fastqs()').

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  dry_run : bool
    Whether to print the commands to stdout instead of running them.

  3  threads : int
    Number of threads to use.

  4  fil_wrk : file
    BAM work file to replace with the coordinate-sorted alignments.

  5  bam_qnam : file
    Queryname-sorted BAM file from Step 2.

Returns
-------
  Returns 0 when sorting succeeds or after printing the commands for a dry run; otherwise, prints an error and returns 1.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - mv
    - rm
    - samtools

  The work file is written at compression level 1 and is not indexed. Step 4 reads it whole and rewrites it, and Step 5 indexes the final file.

Examples
--------
  1. Coordinate-sort a work file and delete the queryname-sorted file.
    '''bash
    _sort_coord_bam \\
        align_fastqs \\
        false \\
        2 \\
        work/tiny_pe.bam \\
        work/tiny_pe.qnam.bam
    '''

  2. Display the helper's argument contract.
    '''bash
    _sort_coord_bam --help
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ "${dry_run}" == "true" ]]; then
        printf '\n%s\n' \
            "# Step 3: sort the BAM work file by coordinates."
    fi

    if [[ "${dry_run}" == "true" || -f "${bam_qnam}" ]]; then
        bam_coor="${fil_wrk%.bam}.coor.bam"
        arr_cmd=(
            samtools sort
                -@ "${threads}"
                -O "bam,level=1"
                -o "${bam_coor}"
                "${bam_qnam}"
        )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #3: failed to coordinate-sort queryname-sorted BAM file."
            return 1
        fi

        # With '--qname', Step 6 name-sorts the final, duplicate-marked BAM, so
        # this file is not kept.
        arr_cmd=( rm -f "${bam_qnam}" )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #3: failed to delete intermediate queryname-sorted BAM" \
                "file."
            return 1
        fi

        arr_cmd=( mv -f "${bam_coor}" "${fil_wrk}" )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #3: failed to rename coordinate-sorted BAM file."
            return 1
        fi
    else
        echo_err_func "${func}" \
            "Step #3: queryname-sorted BAM file from Step 2 does not exist."
        return 1
    fi
}


# For Step 4: decide whether 'samtools markdup' can tell optical duplicates
# apart based on whether every read name carries flow-cell coordinates and all
# tile prefixes have one length.
# TODO: revisit once pinned Samtools includes the fix for Samtools Issue #2399.
function _check_optical_names() {
    local func="${1:-}"
    local threads="${2:-}"
    local fil_wrk="${3:-}"
    local dist_nam="${4:-}"
    local cnt n_nam n_bad n_len fil
    local -a arr_fil
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _check_optical_names
    [--help] func threads fil_wrk dist_nam [fil_add ...]

  Decide whether 'samtools markdup -d' can tell optical duplicates apart in a BAM file, or in several marked together (for Step 4 of 'align_fastqs()').

  Every read name must carry Illumina flow-cell coordinates, and all tile prefixes ('machine:run:flowcell:lane:tile') must have one length.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  threads : int
    Number of threads to use for 'samtools view'.

  3  fil_wrk : file
    BAM file whose read names to check.

  4  dist_nam : str
    Name of the variable that receives the 'OPTICAL DISTANCE' value: '2500' or 'not used' and the reason.

  5+  fil_add : file
    Further BAM files whose read names are checked with those of 'fil_wrk', for files that 'samtools markdup' reads together.

Returns
-------
  Sets 'dist_nam' and returns 0 when optical duplicates can be told apart.

  Sets it, prints a warning, and returns 1 when a read name has no coordinates or the tile prefixes differ in length.

  Sets it and returns 1 silently when the file has no reads.

  Prints an error and returns 2 when the names cannot be read.

Notes
-----
  Runtime requirements:
    - awk
    - bash >= 4.4
    - cut
    - samtools

  - A name is parsed as Samtools 1.24 parses it: it must have 3, 4, 6, or 7 colons, with digits where the x and y coordinates go; the tile prefix is everything before x.
  - One name per template is read (primary records that are not the second mate).
  - Samtools 1.24 calls duplicates on different tiles optical when tile prefixes in a duplicate set differ in length, which includes names it cannot parse (Samtools Issue #2399), so '-d' is not passed then.

Examples
--------
  1. Check whether a coordinate-sorted work file allows optical calls.
    '''bash
    _check_optical_names align_fastqs 2 work/tiny_pe.bam dist_opt
    '''

  2. Display the helper's argument contract.
    '''bash
    _check_optical_names --help
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    arr_fil=( "${fil_wrk}" "${@:5}" )

    # Count names, names without coordinates, and distinct prefix lengths,
    # following the colon count that 'get_coordinates_colons()' in
    # 'bam_markdup.c' uses.
    # shellcheck disable=SC2016  # The program is for awk, not the shell.
    if ! \
        cnt="$(
            set -o pipefail
            for fil in "${arr_fil[@]}"; do
                samtools view -@ "${threads}" -F 0x980 "${fil}" || exit 1
            done \
                | cut -f 1 \
                | awk -F ':' '
                    {
                        n_nam++
                        sep = NF - 1

                        if (sep == 3) { ix = 3; iy = 4 }
                        else if (sep == 4) { ix = 4; iy = 5 }
                        else if (sep == 6 || sep == 7) { ix = 6; iy = 7 }
                        else { n_bad++; next }

                        if ($ix !~ /^[-+]?[0-9]/ || $iy !~ /^[-+]?[0-9]/) {
                            n_bad++
                            next
                        }

                        len = 0
                        for (i = 1; i < ix; i++) len += length($i) + 1
                        lens[len] = 1
                    }

                    END {
                        n_len = 0
                        for (k in lens) n_len++
                        print n_nam + 0, n_bad + 0, n_len
                    }
                '
        )"
    then
        echo_err_func "${func}" \
            "Step #4: failed to read the read names of coordinate-sorted BAM" \
            "file."
        return 2
    fi

    local -n dist_ref="${dist_nam}"
    read -r n_nam n_bad n_len <<< "${cnt}"

    if (( n_nam == 0 )); then
        dist_ref="not used; no reads"
        return 1
    elif (( n_bad > 0 )); then
        dist_ref="not used; ${n_bad} of ${n_nam} read names carry no"
        dist_ref+=" flow-cell coordinates"
        echo_warn_func "${func}" \
            "Step #4: ${n_bad} of ${n_nam} read names carry no flow-cell" \
            "coordinates, so optical duplicates are not told apart."
        return 1
    elif (( n_len > 1 )); then
        dist_ref="not used; read names have ${n_len} tile-prefix lengths"
        dist_ref+=" (Samtools Issue #2399)"
        echo_warn_func "${func}" \
            "Step #4: read names have ${n_len} tile-prefix lengths, so" \
            "optical duplicates are not told apart (Samtools Issue #2399)."
        return 1
    fi

    dist_ref="2500"
    return 0
}


# Step 4: mark duplicate alignments in the coordinate-sorted BAM work file.
function _mark_dup_bam() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local fil_wrk="${4:-}"
    local txt_mrk="${5:-}"
    local dist_nam="${6:-}"
    local fin="${7:-true}"
    local bam_mrk dist_opt rc_opt
    local -a arr_cmd arr_opt
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _mark_dup_bam
    [--help] func dry_run threads fil_wrk txt_mrk dist_nam [fin]

  Mark duplicate alignments in the coordinate-sorted BAM work file and write the duplicate statistics, skipping files already marked (Step 4 of 'align_fastqs()'). The file is not indexed, as Step 5 rewrites its header and indexes it.

  Optical duplicates are told apart ('-d 2500') only when '_check_optical_names()' allows it.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  dry_run : bool
    Whether to print the commands to stdout instead of running them; a dry run does not check for earlier marking.

  3  threads : int
    Number of threads to use.

  4  fil_wrk : file
    Coordinate-sorted BAM work file from Step 3.

  5  txt_mrk : file
    Path for the duplicate statistics, without '.gz'; the file is written gzip-compressed, as 'txt_mrk.gz'.

  6  dist_nam : str
    Name of the variable to set to the optical duplicate distance, '2500', or to 'not used' and the reason; set to empty when the file was already marked.

  7  fin : bool
    Whether the marked file ships as it is, at the default compression level (6). 'false' writes it at level 1, since a later step rewrites it: the split of '--unal_mode separate', or Step 6's CRAM conversion (default: '${fin}').

Returns
-------
  Returns 0 when marking succeeds or is skipped or after printing the commands for a dry run; prints an error and returns 1 otherwise.

Notes
-----
  Runtime requirements:
    - awk
    - bash >= 4.4
    - cut
    - grep
    - gzip
    - mv
    - samtools

  A dry run cannot read the names, so it prints the '-d 2500' form after a comment saying when it applies.

Examples
--------
  1. Mark duplicates in a coordinate-sorted work file and write its statistics.
    '''bash
    _mark_dup_bam \\
        align_fastqs \\
        false \\
        2 \\
        work/tiny_pe.bam \\
        work/tiny_pe.markdup.txt \\
        dist_mrk
    '''

  2. Display the helper's argument contract.
    '''bash
    _mark_dup_bam --help
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ "${dry_run}" == "true" ]]; then
        # Note: the weird line break is necessary for output text formatting.
        printf '\n%s %s\n%s\n' \
            "# Step 4: mark duplicates in the BAM work file." \
            "'-d 2500' applies only if all" \
            "# read names carry flow-cell coordinates of one format."
    fi

    local -n dist_out="${dist_nam}"
    dist_out=""

    if [[ "${dry_run}" == "true" || -f "${fil_wrk}" ]]; then
        if \
            [[ "${dry_run}" == "false" ]] \
                && samtools view -H "${fil_wrk}" \
                | grep -q "@PG.*CL:samtools markdup"
        then
            echo_err_func "${func}" \
                "Step #4: duplicate alignments are already marked in" \
                "coordinate-sorted BAM file. Skipping markdup operations."
        else
            bam_mrk="${fil_wrk%.bam}.mark.bam"

            # Tell optical duplicates apart only when the read names allow it.
            arr_opt=()
            dist_opt="2500"

            if [[ "${dry_run}" == "false" ]]; then
                rc_opt=0
                _check_optical_names \
                    "${func}" \
                    "${threads}" \
                    "${fil_wrk}" \
                    dist_opt \
                    || rc_opt=$?

                if (( rc_opt == 2 )); then
                    return 1
                fi
            fi

            # shellcheck disable=SC2034  # Read by name in the caller.
            dist_out="${dist_opt}"
            if [[ "${dist_opt}" == "2500" ]]; then
                arr_opt=( -d 2500 )
            fi

            # A work file that a later step rewrites is written at compression
            # level 1. ('markdup -l' is the maximum read length, not a level.)
            if [[ "${fin}" == "false" ]]; then
                arr_opt+=( -O "bam,level=1" )
            fi

            arr_cmd=(
                samtools markdup
                    -@ "${threads}"
                    -t
                    "${arr_opt[@]}"
                    -f "${txt_mrk}"
                    "${fil_wrk}"
                    "${bam_mrk}"
            )

            if ! \
                _run_or_print_cmd "${dry_run}" arr_cmd
            then
                echo_err_func "${func}" \
                    "Step #4: failed to mark duplicates in coordinate-sorted" \
                    "BAM file."
                return 1
            fi

            # As the statistics' last line, record whether optical duplicates
            # were told apart and, if not, why.
            if [[ "${dry_run}" == "true" ]]; then
                printf '%s %q\n\n' \
                    "printf 'OPTICAL DISTANCE: %s\\n' 2500 >>" \
                    "${txt_mrk}"
            else
                printf 'OPTICAL DISTANCE: %s\n' "${dist_opt}" >> "${txt_mrk}"
            fi

            arr_cmd=( gzip -f "${txt_mrk}" )

            if ! \
                _run_or_print_cmd "${dry_run}" arr_cmd
            then
                echo_err_func "${func}" \
                    "Step #4: failed to compress duplicate statistics."
                return 1
            fi

            # Replace the original coordinate-sorted BAM with one in which
            # duplicate alignments are marked.
            arr_cmd=( mv -f "${bam_mrk}" "${fil_wrk}" )

            if ! \
                _run_or_print_cmd "${dry_run}" arr_cmd
            then
                echo_err_func "${func}" \
                    "Step #4: failed to rename duplicate-marked," \
                    "coordinate-sorted BAM file."
                return 1
            fi

        fi
    else
        echo_err_func "${func}" \
            "Step #4: coordinate-sorted BAM file does not exist."
        return 1
    fi
}


# Build the pipeline record ('OUTPUT.ALIGNMENT.RECORD') that Step 5 writes into
# the header of the BAM work file as '@CO' lines.
function _build_align_record() {
    local func="${1:-}"
    local cnt_nam="${2:-}"
    local run_nam="${3:-}"
    local threads="${4:-}"
    local aligner="${5:-}"
    local bt2_mode="${6:-}"
    local bt2_X="${7:-}"
    local bt2_inherited="${8:-}"
    local bt2_relaxed="${9:-}"
    local bwa_alg="${10:-}"
    local mapq="${11:-}"
    local req_flg="${12:-}"
    local unal_mode="${13:-}"
    local index="${14:-}"
    local ref_fa="${15:-}"
    local fq_1="${16:-}"
    local fq_2="${17:-}"
    local fil_out="${18:-}"
    local qname="${19:-}"
    local md5_idx md5_ref md5_fq_1 md5_fq_2 prg_aln pkg_aln pkg_sam
    local -a arr_pth
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _build_align_record
    [--help] func cnt_nam run_nam threads aligner bt2_mode bt2_X bt2_inherited bt2_relaxed bwa_alg mapq req_flg unal_mode index ref_fa fq_1 fq_2 fil_out qname

  Fill two named arrays with the '@CO' lines of the pipeline record for 'align_fastqs()': one with the 'pipeline:' lines for what determines the file's content, the other with the 'pipeline-run:' lines for the run's context.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  01  func : str
    Name of the calling function for diagnostics.

  02  cnt_nam : str
    Name of the indexed array to fill with the content lines, one header line per element.

  03  run_nam : str
    Name of the indexed array to fill with the run lines, one header line per element.

  04  threads : int
    Number of threads to use.

  05  aligner : str
    Alignment program: 'bowtie2', 'bwa', or 'bwa-mem2'.

  06  bt2_mode : str
    Bowtie 2 alignment type: 'local', 'global', or 'end-to-end'.

  07  bt2_X : int
    Maximum fragment length for a Bowtie 2 paired-end alignment.

  08  bt2_inherited : bool
    Whether the earlier "inherited" Bowtie 2 paired-end call is made.

  09  bt2_relaxed : bool
    Whether the relaxed Bowtie 2 paired-end call is made.

  10  bwa_alg : str
    BWA algorithm: 'mem' or 'aln'.

  11  mapq : int
    MAPQ threshold.

  12  req_flg : bool
    Whether properly paired alignments are required.

  13  unal_mode : str
    What is done with unaligned reads and their mates: 'remove' or 'separate'.

  14  index : path
    Path to the aligner index/reference.

  15  ref_fa : file
    Reference FASTA, recorded only for CRAM output, or empty.

  16  fq_1 : file
    First FASTQ input file.

  17  fq_2 : file
    Second FASTQ input file, or empty for single-end data.

  18  fil_out : file
    Path to the final alignment file.

  19  qname : bool
    Whether a queryname-sorted copy is also written.

Returns
-------
  Returns 0 after filling the arrays; prints an error and returns 1 if 'cnt_nam' or 'run_nam' is missing or an input cannot be checksummed.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - date
    - dirname
    - git (optional; without it, the commit is 'unknown')
    - hostname
    - md5sum

  - Options are recorded only where they are in effect; for example, '--bt2_X' is left out for single-end data and with '--bt2_inherited'.
  - Inputs are named without their directories in the 'pipeline:' lines; their full paths, the thread count, the start time, the host, the call path, the environment, and any Slurm job are in the 'pipeline-run:' lines.
  - The commit is marked '-dirty' when tracked files have uncommitted changes; untracked files do not count.
  - Each input's name is followed by the MD5 checksum of its bytes ('index_md5', 'ref_fa_md5', 'fq_1_md5', 'fq_2_md5'), from 'print_record_md5' and 'print_record_index_md5'.
  - 'pkg_<aligner>' and 'pkg_samtools' give each program's Conda package version and build, from 'print_record_conda_pkg', or 'unknown' outside Conda.

Examples
--------
  1. Build the record for a single-end Bowtie 2 run and print it.
    '''bash
    _build_align_record \\
        align_fastqs \\
        arr_rec_cnt \\
        arr_rec_run \\
        1 \\
        bowtie2 \\
        end-to-end \\
        1000 \\
        false \\
        false \\
        mem \\
        1 \\
        false \\
        remove \\
        tests/fixtures/align_fastqs/bowtie2/tiny \\
        '' \\
        tests/fixtures/align_fastqs/fastq/se/tiny_se.atria.fastq.gz \\
        '' \\
        work/tiny_se.bam \\
        false
    printf '%s\\n' "\${arr_rec_cnt[@]}" "\${arr_rec_run[@]}"
    '''

  2. Display the helper's argument contract.
    '''bash
    _build_align_record --help
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${cnt_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'cnt_nam', is missing."
        return 1
    elif [[ -z "${run_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'run_nam', is missing."
        return 1
    fi

    local -n cnt_ref="${cnt_nam}"
    local -n run_ref="${run_nam}"
    cnt_ref=()
    # shellcheck disable=SC2034  # Filled by name in 'add_record_run_lines'.
    run_ref=()

    # Content lines: the chain, the commit, the packages that ran, the options
    # in effect, the inputs by file name, and the steps skipped from the start.
    # Package versions come from Conda, not from each program's own '@PG' line,
    # which can be wrong (BWA-MEM2 2.3 reports 2.2.1, bwa-mem2#276).
    prg_aln="${aligner}"
    pkg_aln="$(print_record_conda_pkg "${prg_aln}")" || return 1
    pkg_sam="$(print_record_conda_pkg samtools)" || return 1
    cnt_ref+=(
        $'@CO\tpipeline: chain=align_fastqs'
        $'@CO\tpipeline: commit='"$(print_record_commit)"
        $'@CO\tpipeline: pkg_'"${prg_aln}=${pkg_aln}"
        $'@CO\tpipeline: pkg_samtools='"${pkg_sam}"
        $'@CO\tpipeline: --aligner='"${aligner}"
    )

    case "${aligner}" in
        bowtie2)
            cnt_ref+=(
                $'@CO\tpipeline: --bt2_mode='"${bt2_mode}"
                $'@CO\tpipeline: bt2_seed=0'
            )
            if [[ -n "${fq_2}" && "${bt2_inherited}" == "true" ]]; then
                cnt_ref+=( $'@CO\tpipeline: --bt2_inherited=true' )
            elif [[ -n "${fq_2}" ]]; then
                cnt_ref+=( $'@CO\tpipeline: --bt2_X='"${bt2_X}" )
            fi

            # Per 'HELP.ALIAS.PUBLIC', the hidden measurement setting records
            # itself and its warning.
            if [[ -n "${fq_2}" && "${bt2_relaxed}" == "true" ]]; then
                cnt_ref+=(
                    $'@CO\tpipeline: --bt2_relaxed=true'
                    $'@CO\tpipeline: warning=not for signal or scaling'
                )
            fi
            ;;

        bwa)
            cnt_ref+=( $'@CO\tpipeline: --bwa_alg='"${bwa_alg}" )
            if [[ "${bwa_alg}" == "mem" ]]; then
                cnt_ref+=( $'@CO\tpipeline: bwa_K=100000000' )
            fi
            ;;

        bwa-mem2)
            cnt_ref+=( $'@CO\tpipeline: bwa_K=100000000' )
            ;;
    esac

    cnt_ref+=( $'@CO\tpipeline: --mapq='"${mapq}" )
    if [[ -n "${fq_2}" ]]; then
        cnt_ref+=( $'@CO\tpipeline: --req_flg='"${req_flg}" )
    fi
    cnt_ref+=( $'@CO\tpipeline: --unal_mode='"${unal_mode}" )

    # Each input is named and then identified by the checksum of its bytes, so
    # a rebuild can confirm it has the same files rather than the same names.
    md5_idx="$(print_record_index_md5 "${aligner}" "${index}")" || return 1
    cnt_ref+=(
        $'@CO\tpipeline: --qname='"${qname}"
        $'@CO\tpipeline: --index='"${index##*/}"
        $'@CO\tpipeline: index_md5='"${md5_idx}"
    )
    if [[ "${fil_out}" == *.cram ]]; then
        md5_ref="$(print_record_md5 "${ref_fa}")" || return 1
        cnt_ref+=(
            $'@CO\tpipeline: --ref_fa='"${ref_fa##*/}"
            $'@CO\tpipeline: ref_fa_md5='"${md5_ref}"
        )
    fi
    md5_fq_1="$(print_record_md5 "${fq_1}")" || return 1
    cnt_ref+=(
        $'@CO\tpipeline: --fq_1='"${fq_1##*/}"
        $'@CO\tpipeline: fq_1_md5='"${md5_fq_1}"
    )
    if [[ -n "${fq_2}" ]]; then
        md5_fq_2="$(print_record_md5 "${fq_2}")" || return 1
        cnt_ref+=(
            $'@CO\tpipeline: --fq_2='"${fq_2##*/}"
            $'@CO\tpipeline: fq_2_md5='"${md5_fq_2}"
        )
    else
        cnt_ref+=( $'@CO\tpipeline: skipped=fixmate (single-end)' )
    fi

    # Run lines: the thread count, the start time, the host, the caller, the
    # environment, any Slurm job, and the full paths.
    arr_pth=( index "${index}" fq_1 "${fq_1}" fq_2 "${fq_2}" )
    if [[ "${fil_out}" == *.cram ]]; then arr_pth+=( ref_fa "${ref_fa}" ); fi
    arr_pth+=( fil_out "${fil_out}" )

    add_record_run_lines \
        "${run_nam}" "${BASH_SOURCE[0]##*/}::align_fastqs" "${threads}" \
        "${arr_pth[@]}" \
        || return 1
}


# For Step 5: print the record's optical-duplicate line for an optical distance
# that '_mark_dup_bam()' or '_mark_dup_all_bam()' set.
function _print_optical_line() {
    local dist="${1:-}"

    if [[ "${dist}" == "2500" ]]; then
        printf '%s\n' $'@CO\tpipeline: optical_distance=2500'
    elif [[ -n "${dist}" ]]; then
        printf '%s\n' $'@CO\tpipeline: skipped=optical duplicates ('"${dist})"
    else
        printf '%s\n' \
            $'@CO\tpipeline: optical_distance=unknown (marked earlier)'
    fi
}


# Step 5: write the pipeline record into the header of the BAM work file,
# skipping files that already carry a record, then index the file.
function _write_align_record() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local fil_wrk="${4:-}"
    local cnt_nam="${5:-}"
    local run_nam="${6:-}"
    local idx="${7:-true}"
    local skp=false
    local hdr_rec bam_rec
    local -a arr_cmd
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _write_align_record
    [--help] func dry_run threads fil_wrk cnt_nam run_nam [idx]

  Write the pipeline record, content lines first, into the header of the duplicate-marked BAM work file as '@CO' lines, skipping files that already carry a record, then index the file (Step 5 of 'align_fastqs()').

  This is the only index the chain builds for a BAM file, so a skipped file is indexed too.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  dry_run : bool
    Whether to print the commands to stdout instead of running them; a dry run does not check for an earlier record.

  3  threads : int
    Number of threads to use.

  4  fil_wrk : file
    Duplicate-marked, coordinate-sorted BAM work file from Step 4.

  5  cnt_nam : str
    Name of the indexed array holding the record's content lines.

  6  run_nam : str
    Name of the indexed array holding the record's run lines.

  7  idx : bool
    Whether to index the file; 'false' for CRAM output, as Step 6 converts the file and indexes the CRAM file instead (default: '${idx}').

Returns
-------
  Returns 0 when writing and indexing succeed, when writing is skipped and indexing succeeds, or after printing the commands for a dry run; prints an error and returns 1 otherwise.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - grep
    - mv
    - rm
    - samtools

  - The header is read with '--no-PG', so reading it adds no '@PG' line; 'samtools reheader' adds its own.

Examples
--------
  1. Write a record built by '_build_align_record()' into a work file.
    '''bash
    _write_align_record \\
        align_fastqs \\
        false \\
        2 \\
        work/tiny_pe.bam \\
        arr_rec_cnt \\
        arr_rec_run
    '''

  2. Display the helper's argument contract.
    '''bash
    _write_align_record --help
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ "${dry_run}" == "true" ]]; then
        printf '\n%s\n' \
            "# Step 5: write the pipeline record into the BAM work file."
    fi

    if [[ "${dry_run}" == "false" && ! -f "${fil_wrk}" ]]; then
        echo_err_func "${func}" \
            "Step #5: duplicate-marked BAM file does not exist."
        return 1
    elif \
        [[ "${dry_run}" == "false" ]] \
            && samtools view -H --no-PG "${fil_wrk}" \
                | grep -q $'^@CO\tpipeline: chain=align_fastqs$'
    then
        echo_err_func "${func}" \
            "Step #5: duplicate-marked BAM file already carries a pipeline" \
            "record. Skipping record operations."
        skp=true
    fi

    if [[ "${skp}" == "false" ]]; then
        local -n rec_cnt_ref="${cnt_nam}"
        local -n rec_run_ref="${run_nam}"
        hdr_rec="${fil_wrk%.bam}.header.sam"
        bam_rec="${fil_wrk%.bam}.record.bam"

        # Write the header, then the record after it.
        arr_cmd=( samtools view -H --no-PG "${fil_wrk}" )

        if [[ "${dry_run}" == "true" ]]; then
            _print_dry_run_cmd --out "${hdr_rec}" arr_cmd
            printf 'printf %q' '%s\n'
            printf ' \\\n    %q' "${rec_cnt_ref[@]}" "${rec_run_ref[@]}"
            printf ' \\\n    >> %q\n\n' "${hdr_rec}"
        elif ! \
            {
                "${arr_cmd[@]}" \
                    && printf '%s\n' "${rec_cnt_ref[@]}" "${rec_run_ref[@]}"
            } > "${hdr_rec}"
        then
            echo_err_func "${func}" \
                "Step #5: failed to write the header with the pipeline record."
            return 1
        fi

        # Replace the header of the BAM work file with the new one.
        arr_cmd=( samtools reheader "${hdr_rec}" "${fil_wrk}" )

        if [[ "${dry_run}" == "true" ]]; then
            _print_dry_run_cmd --out "${bam_rec}" arr_cmd
        elif ! \
            "${arr_cmd[@]}" > "${bam_rec}"
        then
            echo_err_func "${func}" \
                "Step #5: failed to write the pipeline record into the" \
                "duplicate-marked BAM file."
            return 1
        fi

        arr_cmd=( mv -f "${bam_rec}" "${fil_wrk}" )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #5: failed to rename the BAM file with the pipeline" \
                "record."
            return 1
        fi

        arr_cmd=( rm -f "${hdr_rec}" )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #5: failed to remove the header file."
            return 1
        fi
    fi

    # Index the file here, after the new header moves the records, and not in
    # Steps 3 and 4, which rewrite the file.
    if [[ "${idx}" == "false" ]]; then return 0; fi

    arr_cmd=( samtools index -@ "${threads}" "${fil_wrk}" )

    if ! \
        _run_or_print_cmd "${dry_run}" arr_cmd
    then
        echo_err_func "${func}" \
            "Step #5: failed to index the BAM file with the pipeline record."
        return 1
    fi
}


# Step 6: with '--qname', name-sort the final BAM work file into a copy; then,
# for CRAM output, convert the work file and any copy to CRAM.
function _finalize_align_output() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local out_fmt="${4:-}"
    local fil_out="${5:-}"
    local fil_wrk="${6:-}"
    local ref_fa="${7:-}"
    local qname="${8:-}"
    local bam_qnam="${9:-}"
    local fil_qnam_out="${10:-}"
    local hdr="${11:-true}"
    local -a arr_cmd arr_lvl
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _finalize_align_output
    [--help] func dry_run threads out_fmt fil_out fil_wrk ref_fa qname bam_qnam fil_qnam_out [hdr]

  With 'qname', name-sort the final, duplicate-marked BAM work file into a queryname-sorted copy; then convert the work file and the copy to CRAM when requested, or rename the copy into place (Step 6 of 'align_fastqs()').

  For CRAM output, the copy is a work file and is written at compression level 1.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  01  func : str
    Name of the calling function for diagnostics.

  02  dry_run : bool
    Whether to print the commands to stdout instead of running them.

  03  threads : int
    Number of threads to use for sorting and indexing.

  04  out_fmt : str
    Final output format: 'bam' or 'cram'.

  05  fil_out : file
    Path to the final alignment file.

  06  fil_wrk : file
    BAM work file from Step 5.

  07  ref_fa : file
    Reference FASTA path, used for CRAM output.

  08  qname : bool
    Whether to write a queryname-sorted copy of the final alignments.

  09  bam_qnam : file
    Path for the queryname-sorted BAM copy of the work file.

  10  fil_qnam_out : file
    Path for the queryname-sorted output.

  11  hdr : bool
    Whether a dry run prints the step's heading; 'false' for the side file, which follows the main output (default: '${hdr}').

Returns
-------
  Returns 0 when conversion and cleanup succeed or after printing the commands for a dry run; prints an error and returns 1 otherwise.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - mv
    - Reference FASTA and required index (when 'out_fmt' is 'cram')
    - rm
    - samtools

Examples
--------
  1. Name-sort a paired-end work file into a copy, then convert both to CRAM.
    '''bash
    _finalize_align_output \\
        align_fastqs \\
        false \\
        2 \\
        cram \\
        work/tiny_pe.cram \\
        work/tiny_pe.bam \\
        tests/fixtures/align_fastqs/reference/tiny.fa \\
        true \\
        work/tiny_pe.qnam.bam \\
        work/tiny_pe.qnam.cram
    '''

  2. Print, without running, the step for BAM output with a queryname-sorted copy.
    '''bash
    _finalize_align_output \\
        align_fastqs \\
        true \\
        1 \\
        bam \\
        work/tiny_se.bam \\
        work/tiny_se.bam \\
        '' \\
        true \\
        work/tiny_se.qnam.bam \\
        work/tiny_se.qnam.bam
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ "${dry_run}" == "true" && "${hdr}" == "true" ]]; then
        printf '\n%s\n' \
            "# Step 6: write the final output and the retained files."

        # BAM output without '--qname' is the BAM work file itself, so nothing
        # runs.
        if [[ "${out_fmt}" == "bam" && "${qname}" == "false" ]]; then
            printf '%s\n\n' \
                "# Nothing to do: the BAM work file is the final output."
        fi
    fi

    # Queryname-sort the final, duplicate-marked work file. For CRAM output,
    # the copy is converted below, so it is written at compression level 1.
    if [[ "${qname}" == "true" ]]; then
        arr_lvl=()
        if [[ "${out_fmt}" == "cram" ]]; then arr_lvl=( -O "bam,level=1" ); fi

        arr_cmd=(
            samtools sort
                -@ "${threads}"
                -n
                "${arr_lvl[@]}"
                -o "${bam_qnam}"
                "${fil_wrk}"
        )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #6: failed to queryname-sort the final BAM work file."
            return 1
        fi
    fi

    if [[ "${out_fmt}" == "cram" ]]; then
        arr_cmd=(
            samtools view
                -T "${ref_fa}"
                -O cram
                -o "${fil_out}"
                "${fil_wrk}"
        )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #6: failed to convert BAM work file to final CRAM file."
            return 1
        fi

        arr_cmd=( samtools index -@ "${threads}" "${fil_out}" )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #6: failed to index final CRAM file."
            return 1
        fi

        if [[ "${qname}" == "true" ]]; then
            arr_cmd=(
                samtools view
                    -T "${ref_fa}"
                    -O cram
                    -o "${fil_qnam_out}"
                    "${bam_qnam}"
            )

            if ! \
                _run_or_print_cmd "${dry_run}" arr_cmd
            then
                echo_err_func "${func}" \
                    "Step #6: failed to convert queryname-sorted BAM" \
                    "intermediate to retained CRAM file."
                return 1
            fi
        fi

        arr_cmd=( rm -f "${fil_wrk}" "${fil_wrk}.bai" )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #6: failed to delete BAM work file and/or BAM index" \
                "after CRAM conversion."
            return 1
        fi

        if [[ "${qname}" == "true" ]]; then
            arr_cmd=( rm -f "${bam_qnam}" )

            if ! \
                _run_or_print_cmd "${dry_run}" arr_cmd
            then
                echo_err_func "${func}" \
                    "Step #6: failed to delete queryname-sorted BAM work" \
                    "intermediate after CRAM conversion."
                return 1
            fi
        fi
    elif [[ "${out_fmt}" == "bam" && "${qname}" == "true" ]]; then
        if [[ "${bam_qnam}" != "${fil_qnam_out}" ]]; then
            # shellcheck disable=SC2034  # Read by name in '_run_or_print_cmd'.
            arr_cmd=( mv -f "${bam_qnam}" "${fil_qnam_out}" )

            if ! \
                _run_or_print_cmd "${dry_run}" arr_cmd
            then
                echo_err_func "${func}" \
                    "Step #6: failed to rename retained queryname-sorted BAM" \
                    "intermediate."
                return 1
            fi
        fi
    fi
}


# TODO: "semantic paragraphs" for options under Usage.
function align_fastqs() {
    local dry_run threads aligner bt2_mode bt2_X bt2_inherited bt2_relaxed
    local bwa_alg mapq req_flg unal_mode index ref_fa
    local bt2_mode_set bt2_X_set bwa_alg_set
    local fq_1 fq_2 fil_out qname out_fmt
    local fil_wrk fil_qnam_out bam_qnam txt_mrk dist_mrk
    local txt_all dist_all fil_tsv lin_opt lin_opt_all kind
    local fin_mrk fin_spl idx_rec
    local -A arr_wrk arr_out arr_qnam_out arr_bam_qnam arr_cnt
    local -a arr_kind

    # shellcheck disable=SC2034  # Filled and read by name in the helpers.
    local -a arr_rec_cnt arr_rec_run arr_rec_unal arr_rec_improp arr_rec_mapq
    local show_help

    # Assign default argument values.
    dry_run=false
    threads=1
    aligner="bowtie2"
    bt2_mode="end-to-end"
    bt2_mode_set=false
    bt2_X=1000
    bt2_X_set=false
    bt2_inherited=false
    bt2_relaxed=false
    bwa_alg="mem"
    bwa_alg_set=false
    ref_fa=""
    mapq=1
    req_flg=false
    unal_mode="remove"
    fq_1=""
    fq_2=""
    fil_out=""
    qname=false

    show_help=$(cat << EOM
Usage
-----
  align_fastqs
    [--help] [--dry_run] [--threads <int>] [--aligner <aligner>] [--bt2_mode <mode>] [--bt2_X <int>] [--bwa_alg <algorithm>] [--mapq <int>] [--req_flg] [--unal_mode <mode>] --index <path> [--ref_fa <file>] --fq_1 <file> [--fq_2 <file>] --fil_out <file> [--qname]

  Align single- or paired-end Illumina short reads using Bowtie 2, BWA, or BWA-MEM2, and convert the output to a sorted, duplicate-marked, mate-fixed (if working with paired-end reads) alignment file using Samtools.

  Final output format is inferred from '--fil_out': '.bam' or '.cram'. Internal post-processing is performed in BAM format. If CRAM output is requested, the final BAM work file and any retained queryname-sorted BAM intermediate are converted to CRAM in the final step.

  When '--aligner bwa' is selected, '--bwa_alg' may be used to choose either 'mem' or the older backtrack workflow ('aln' with downstream 'samse' or 'sampe').

  Optionally, '--qname' also writes a queryname-sorted copy of the final alignment file.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -dr, --dry_run : flag
    Run script in dry-run mode. Validate the arguments, then print each step's commands to stdout, in order, without running them or writing files. Commands that read a file an earlier step would write are printed as they would run.

  -t, --thr, --threads : int
    Number of threads to use (default: ${threads}).

  -a, --aln, --aligner : {'bowtie2', 'bwa', 'bwa-mem2'}
    Alignment program to use: 'bowtie2', 'bwa', or 'bwa-mem2' (default: '${aligner}').

  -2m, --bt2_mode : {'local', 'global', 'end-to-end'}
    Bowtie 2 alignment type when '--aligner bowtie2': 'local', 'global', or 'end-to-end' (default: '${bt2_mode}'); refused otherwise.

  -2X, --bt2_X : int
    Maximum fragment length for a Bowtie 2 paired-end alignment, passed as 'bowtie2 -X' (default: ${bt2_X}); refused unless '--aligner bowtie2'; ignored with a warning for single-end data.

  -ba, --bwa_alg : {'mem', 'aln'}
    BWA algorithm when '--aligner bwa': 'mem' or 'aln' (default: '${bwa_alg}').

    With '--aligner bwa-mem2', which has only 'mem', 'mem' is ignored with a warning and 'aln' is refused; with '--aligner bowtie2', it is refused.

  -mq, --mpq, --mapq : int
    MAPQ threshold for filtering alignment output files (default: ${mapq}).

    To disable MAPQ-based filtering, specify 0.

    With paired-end data, a pair is kept only when both mates pass. Secondary (SAM flag 256) and supplementary (SAM flag 2048) alignments are kept or dropped with their associated primary read alignments, whatever their own MAPQ. With Bowtie 2, a MAPQ of 1 or more does not mean a unique alignment, as a read alignment with multiple equally good placements can get MAPQ 1, and its placement is chosen at random.

  -rq, --req_flg : flag
    Require SAM flag bit 2 for properly paired alignments.

    The flag is checked as the aligner set it and again after 'samtools fixmate', which clears it on a pair whose mates map to different references or align in the same orientation (both forward or both reverse).

    Ignored with a warning for single-end data.

  -um, --unal_mode : {'remove', 'separate'}
    What to do with the records that the alignment file leaves out: 'remove' or 'separate' (default: '${unal_mode}'). These are the unaligned reads and, in paired-end data, their mates (aligned or not), plus the templates that '--req_flg' and '--mapq' reject.

    With 'remove', they are dropped. With 'separate', nothing is dropped: each goes to a diagnostic file beside the alignment file, and the alignment file holds the same records as with 'remove' (see Notes).

    Refused as 'separate' with '--bt2_inherited', whose call writes no unaligned reads.

  -ix, --idx, --index : path
    Path to the aligner index/reference.

  -rf, --ref_fa : file
    Reference FASTA path required when '--fil_out' ends in '.cram'; otherwise, ignored with a warning.

  -f1, --fq_1, --fastq_1 : file
    First FASTQ input file. For paired-end data, this is read 1.

  -f2, --fq_2, --fastq_2 : file
    Second FASTQ input file. Required for paired-end data.

  -fo, --fil_out : file
    Path to the final alignment file (must end in '.bam' or '.cram').

  -qn, --qnam, --qname : flag
    Also write a queryname-sorted copy of the final alignment file.

Returns
-------
  Returns 0 when alignment and post-processing complete successfully; otherwise 1.

Output
------
  Creates a final '.bam' or '.cram' file at '--fil_out'. A corresponding index ('.bai' or '.crai') is also created.

  With '--unal_mode separate', also creates the diagnostic files '<stem>.unal.bam', '<stem>.improp.bam' (paired-end data only), and '<stem>.mapq.bam', or the same with '.cram'.

  Each comes with its index and, with '--qname', a '.qnam' copy (for example, '<stem>.unal.qnam.bam').

  It also writes Markdup 2's duplicate statistics ('<stem>.all.markdup.txt.gz') and the counts ('<stem>.diag.tsv.gz').

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - bowtie2 (when '--aligner bowtie2' is specified)
    - bwa (when '--aligner bwa' is specified)
    - bwa-mem2 (when '--aligner bwa-mem2' is specified)
    - grep
    - md5sum
    - Reference FASTA and required index (when writing CRAM)
    - samtools

  - For '--index' when using Bowtie 2, the path should end with the index stem: 'path/to/dir/stem'; when using 'bwa' or 'bwa-mem2', the path should be the indexed reference FASTA path (for example, '.fa').
  - '--bwa_alg' is used only when '--aligner bwa'. It is refused with '--aligner bowtie2'; with '--aligner bwa-mem2', which always runs 'bwa-mem2 mem', 'mem' is ignored with a warning and 'aln' is refused.
  - Unaligned reads are excluded from the output; with paired-end data, so are the mates of unaligned reads. With '--unal_mode separate', they are written to the side file instead.
  - With '--unal_mode separate', every record the aligner writes lands in exactly one file ('docs/design/align_fastqs_separate.md'):
    + '<stem>.unal.bam' (or '.cram'), the side file, holds the unaligned reads and their mates:
      - With single-end data, these are the reads flagged unaligned (SAM flag 4).
      - With Bowtie 2 and paired-end data, they are the pairs that did not align as concordant pairs. Bowtie 2 writes both mates of such a pair unaligned even when one mate could align alone.
      - With BWA or BWA-MEM2, they are the records flagged unaligned or with an unaligned mate.
      - When BWA-backtrack flags one mate unaligned but leaves its partner's mate flag clear, the partner goes beside its mate.
    + The templates that '--req_flg' and '--mapq' reject go to their own files:
      - '<stem>.improp.bam' holds the pairs that '--req_flg' rejects. Single-end data has no improper file, and without '--req_flg' the improper file is empty.
      - '<stem>.mapq.bam' holds the templates that '--mapq' rejects.
      - A template that fails both is improper.
      - '--mapq' and '--req_flg' do not apply to the side file, as both judge aligned pairs and would remove every pair it holds.
    + Each template is kept whole: all of a read pair's records go to the same file.
    + For each diagnostic file:
      - It is coordinate-sorted and indexed.
      - Its duplicates are marked against all the files together.
      - It carries the pipeline record and follows the alignment file into CRAM and into a '--qname' copy.
    + The diagnostic files are not for signal or scaling.
    + '<stem>.diag.tsv.gz' counts the templates and records each file holds, by category.
  - With Bowtie 2 and paired-end data, mates may overlap or contain one another, so fragments shorter than the two reads together are kept; mates that extend past each other ("dovetailed") are not.
  - Duplicates are marked (SAM flag 1024) in Step #4, after MAPQ filtering, and kept.
    + Like any other read, '--mapq' filters them by the MAPQ the aligner gave them.
    + When every read name carries Illumina flow-cell coordinates ('...:tile:x:y') in one format, duplicates on the same tile as their original and within 2500 pixels of it in x and y are tagged as optical ('dt:Z:SQ'); the rest are tagged as library duplicates ('dt:Z:LB').
    + Otherwise, as with most SRA-dumped names, no 'dt' tags are written and a warning says why.
    + Duplicate statistics are written to '<stem>.markdup.txt.gz' beside the alignment file; Samtools estimates library size for paired-end data only. The file's last line, 'OPTICAL DISTANCE:', gives 2500 when optical duplicates were told apart, or 'not used' and the reason.
    + With '--unal_mode separate', duplicates are marked twice:
      - Markdup 1 runs on the alignment file alone, as with 'remove', so its flags and '<stem>.markdup.txt.gz' are the same in both modes.
      - Markdup 2 runs on all files together, so a read in a diagnostic file can be marked as a duplicate of a read in the alignment file. Its statistics are written to '<stem>.all.markdup.txt.gz'.
  - Bowtie 2 runs with '--reorder', so its output order does not depend on the thread count. It also runs with '--seed 0', its default, which the record gives as 'bt2_seed=0'.
    + Which of several equally good alignments a read gets depends on the seed and on the read's sequence, qualities, and whole FASTQ header line, so renaming reads can move them.
    + For paired-end data, the two mates' seeds are combined, which cancels the seed and any header text identical in both mates; only header differences between the mates move pairs.
  - BWA-MEM and BWA-MEM2 run with '-K 100000000', which fixes the batch of read bases over which paired-end insert sizes are estimated. Their output therefore does not depend on the thread count; the record gives the setting as 'bwa_K=100000000'.
  - Each alignment file's header, the '.qnam' copy's included, carries a pipeline record as '@CO' lines:
    + 'pipeline:' lines give the commit, the options in effect, the input file names, and any skipped steps;
    + 'pipeline-run:' lines give the thread count, the start time, the host, the scripts and functions that ran, and the full paths.
  - The record names each input and gives the MD5 checksum of its bytes, as stored: 'index_md5' for the index files, 'fq_1_md5' and 'fq_2_md5' for the FASTQs, and 'ref_fa_md5' for the reference with CRAM output.
  - It also gives the Conda package version and build of the aligner and Samtools ('pkg_bowtie2=2.5.5-h9e91881_0', for example), since a program's own '@PG' version can be wrong: BWA-MEM2 2.3 reports 2.2.1 (bwa-mem2 issue #276).
  - For '--qname', the retained queryname-sorted output will have the same path and stem assigned to '--fil_out', except '.qnam' will be inserted before the final extension (for example, '.qnam.bam' or '.qnam.cram'). The copy is made in Step #6 from the final, duplicate-marked alignments, so it holds the same records and flags.
  - '--ref_fa' is required when '--fil_out' ends in '.cram', since CRAM writing requires a reference FASTA.

Examples
--------
  1. Align one single-end fixture with Bowtie 2 and write BAM output.
    '''bash
    align_fastqs \\
        --threads 1 \\
        --aligner bowtie2 \\
        --bt2_mode global \\
        --index tests/fixtures/align_fastqs/bowtie2/tiny \\
        --fq_1 tests/fixtures/align_fastqs/fastq/se/tiny_se.atria.fastq.gz \\
        --fil_out work/tiny_se.bam
    '''

  2. Align paired-end fixtures with BWA-MEM2 and retain CRAM queryname output.
    '''bash
    align_fastqs \\
        --threads 4 \\
        --aligner bwa-mem2 \\
        --index tests/fixtures/align_fastqs/bwa-mem2/tiny.fa \\
        --ref_fa tests/fixtures/align_fastqs/reference/tiny.fa \\
        --fq_1 tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R1.atria.fastq.gz \\
        --fq_2 tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        --fil_out work/tiny_pe.cram \\
        --qname
    '''
EOM
    )

    # Step 0: describe, parse, and validate keyword arguments.
    if [[ -z "${1:-}" || "${1:-}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    # Parse keyword arguments.
    while (( $# > 0 )); do
        case "${1}" in
            -dr|--dry[_-]run)
                dry_run=true
                shift 1
                ;;

            -t|--thr|--threads)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                threads="${2:-}"
                shift 2
                ;;

            -a|--aln|--aligner)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                aligner="${2,,}"
                shift 2
                ;;

            -2m|--bt2[_-]mode)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                bt2_mode="${2,,}"
                bt2_mode_set=true
                shift 2
                ;;

            -2X|--bt2[_-]X)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                bt2_X="${2}"
                bt2_X_set=true
                shift 2
                ;;

            # Hidden per 'HELP.ALIAS.PUBLIC': reproduces the Bowtie 2 PE calls
            # made before 2026-10 ('--no-overlap --no-dovetail', no '-X') and
            # their alignments.
            -2i|--bt2[_-]inherited)
                bt2_inherited=true
                shift 1
                ;;

            # Hidden per 'HELP.ALIAS.PUBLIC': a measurement setting that makes
            # the relaxed Bowtie 2 PE call ('-I 0 -X <bt2_X>', without
            # '--no-mixed' or '--no-discordant') to measure how many pairs the
            # default call's pairing options remove.
            -2r|--bt2[_-]relaxed)
                bt2_relaxed=true
                shift 1
                ;;

            -ba|--bwa[_-]alg)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                bwa_alg="${2,,}"
                bwa_alg_set=true
                shift 2
                ;;

            -mq|--mpq|--mapq)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                mapq="${2:-}"
                shift 2
                ;;

            -rq|--req[_-]flg)
                req_flg=true
                shift 1
                ;;

            -um|--unal[_-]mode)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                unal_mode="${2,,}"
                shift 2
                ;;

            -ix|--idx|--index)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                index="${2:-}"
                shift 2
                ;;

            -rf|--ref[_-]fa)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                ref_fa="${2:-}"
                shift 2
                ;;

            -f1|--fq[_-]1|--fastq[_-]1)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                fq_1="${2:-}"
                shift 2
                ;;

            -f2|--fq[_-]2|--fastq[_-]2)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                fq_2="${2:-}"
                shift 2
                ;;

            -fo|--fil[_-]out)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                fil_out="${2:-}"
                shift 2
                ;;

            -qn|--qnam|--qname)
                qname=true
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

    # Validate keyword arguments.
    _validate_args_align_fastqs \
        "${FUNCNAME[0]}" \
        "${threads}" \
        "${aligner}" \
        "${bt2_mode}" \
        "${bt2_mode_set}" \
        "${bt2_X}" \
        "${bt2_X_set}" \
        "${bt2_inherited}" \
        "${bt2_relaxed}" \
        "${bwa_alg}" \
        "${bwa_alg_set}" \
        "${mapq}" \
        "${req_flg}" \
        "${unal_mode}" \
        "${index}" \
        "${ref_fa}" \
        "${fq_1}" \
        "${fq_2}" \
        "${fil_out}" \
        || return 1

    # Derive the work and output paths from the validated '--fil_out'.
    case "${fil_out}" in
        *.bam)
            out_fmt="bam"
            fil_wrk="${fil_out}"
            fil_qnam_out="${fil_out%.bam}.qnam.bam"
            ;;

        *.cram)
            out_fmt="cram"
            fil_wrk="${fil_out%.cram}.bam"
            fil_qnam_out="${fil_out%.cram}.qnam.cram"
            ;;
    esac

    bam_qnam="${fil_wrk%.bam}.qnam.bam"
    txt_mrk="${fil_wrk%.bam}.markdup.txt"

    # The diagnostic files of '--unal_mode separate' share the stem with
    # '.unal', '.improp' (paired-end data only), and/or '.mapq', and so do
    # Markdup 2's statistics ('.all.markdup.txt.gz') and the table of counts
    # ('.diag.tsv.gz').
    arr_kind=( unal )
    if [[ -n "${fq_2}" ]]; then arr_kind+=( improp ); fi
    arr_kind+=( mapq )

    for kind in "${arr_kind[@]}"; do
        arr_wrk[${kind}]="${fil_wrk%.bam}.${kind}.bam"
        arr_out[${kind}]="${arr_wrk[${kind}]}"
        arr_qnam_out[${kind}]="${fil_wrk%.bam}.${kind}.qnam.bam"
        if [[ "${out_fmt}" == "cram" ]]; then
            arr_out[${kind}]="${fil_out%.cram}.${kind}.cram"
            arr_qnam_out[${kind}]="${fil_out%.cram}.${kind}.qnam.cram"
        fi
        arr_bam_qnam[${kind}]="${fil_wrk%.bam}.${kind}.qnam.bam"
    done

    txt_all="${fil_wrk%.bam}.all.markdup.txt"
    fil_tsv="${fil_wrk%.bam}.diag.tsv.gz"

    # Build the pipeline record now, so that its start time is the run's.
    _build_align_record \
        "${FUNCNAME[0]}" \
        arr_rec_cnt \
        arr_rec_run \
        "${threads}" \
        "${aligner}" \
        "${bt2_mode}" \
        "${bt2_X}" \
        "${bt2_inherited}" \
        "${bt2_relaxed}" \
        "${bwa_alg}" \
        "${mapq}" \
        "${req_flg}" \
        "${unal_mode}" \
        "${index}" \
        "${ref_fa}" \
        "${fq_1}" \
        "${fq_2}" \
        "${fil_out}" \
        "${qname}" \
        || return 1


    # With '--unal_mode separate', nothing is discarded: every record the
    # aligner writes lands in exactly one file (for more details, see
    # 'docs/design/align_fastqs_separate.md'):
    # - Step 1 writes the unaligned reads and their mates ("Set U") to one file
    #   and the aligned reads ("Set M"), the main file's candidates, to
    #   another.
    # - Step 2 moves every template with records in both sets into Set U and
    #   fixes mates in each set. It then filters Set M as under 'remove', but
    #   keeps what the filters reject: the improper pairs ("Set I") and the
    #   low-MAPQ templates ("Set Q").
    # - Step 3 sorts the main file, as under 'remove'. It also sorts Sets U, I,
    #   and Q together.
    # - Step 4 marks duplicates in the main file alone (Markdup 1), as under
    #   'remove'. It then marks them in all files together (Markdup 2), keeps
    #   the marked records of Sets U, I, and Q, and splits them into their
    #   files.
    # - Step 5 writes the pipeline record into each file, noting what the file
    #   holds. It indexes each file and counts what each holds in
    #   '.diag.tsv.gz'.
    # - Step 6 does two optional things to each file: (a) converts it to CRAM
    #   and (b) makes a '--qname' copy.
    #
    # With '--unal_mode remove', Set U is dropped in Step 1 and Sets I and Q in
    # Step 2, so the main file is the same in both modes.

    # Step 1: align the reads into the BAM work file, dropping unaligned reads
    # or, with '--unal_mode separate', setting them aside.
    _align_reads \
        "${FUNCNAME[0]}" \
        "${dry_run}" \
        "${threads}" \
        "${aligner}" \
        "${bt2_mode}" \
        "${bt2_X}" \
        "${bt2_inherited}" \
        "${bt2_relaxed}" \
        "${bwa_alg}" \
        "${req_flg}" \
        "${unal_mode}" \
        "${index}" \
        "${fq_1}" \
        "${fq_2}" \
        "${fil_wrk}" \
        || return 1


    # Step 2: sort by queryname, fix mates in paired-end data, and keep or drop
    # each template whole. A template is kept by its primary records' MAPQ and,
    # with '--req_flg', their proper-pair flags. With '--unal_mode separate',
    # the dropped templates are set aside.
    _sort_qname_bam \
        "${FUNCNAME[0]}" \
        "${dry_run}" \
        "${threads}" \
        "${aligner}" \
        "${mapq}" \
        "${req_flg}" \
        "${unal_mode}" \
        "${fq_2}" \
        "${fil_wrk}" \
        "${bam_qnam}" \
        || return 1


    # Steps 3 and 4 do not index the BAM work file, as no step reads it by
    # region and each step rewrites it. Step 5 indexes the final file. Work
    # files that a later step rewrites are written at compression level 1;
    # shipped files are written at the default compression level, 6.
    #
    # Step 3: sort the BAM work file by coordinates. With '--unal_mode
    # separate', also sort the diagnostic work files into their own.
    _sort_coord_bam \
        "${FUNCNAME[0]}" \
        "${dry_run}" \
        "${threads}" \
        "${fil_wrk}" \
        "${bam_qnam}" \
        || return 1

    if [[ "${unal_mode}" == "separate" ]]; then
        _sort_coord_diag_bam \
            "${FUNCNAME[0]}" \
            "${dry_run}" \
            "${threads}" \
            "${fq_2}" \
            "${fil_wrk}" \
            || return 1
    fi


    # Step 4: mark duplicates in the main work file alone (Markdup 1), as in
    # every mode. With '--unal_mode separate', also mark them in all files
    # together (Markdup 2), keeping only the diagnostic files' records. The
    # marked files ship as they are unless Step 6 converts them to CRAM.
    fin_mrk=true
    if [[ "${out_fmt}" == "cram" ]]; then fin_mrk=false; fi

    _mark_dup_bam \
        "${FUNCNAME[0]}" \
        "${dry_run}" \
        "${threads}" \
        "${fil_wrk}" \
        "${txt_mrk}" \
        dist_mrk \
        "${fin_mrk}" \
        || return 1

    if [[ "${unal_mode}" == "separate" ]]; then
        _mark_dup_all_bam \
            "${FUNCNAME[0]}" \
            "${dry_run}" \
            "${threads}" \
            "${fil_wrk}" \
            "${txt_all}" \
            dist_all \
            || return 1
    fi

    # With '--unal_mode separate', split the marked diagnostic records into
    # their files.
    if [[ "${unal_mode}" == "separate" ]]; then
        fin_spl="${fin_mrk}"

        _split_diag_bam \
            "${FUNCNAME[0]}" \
            "${dry_run}" \
            "${threads}" \
            "${fq_2}" \
            "${fil_wrk}" \
            "${arr_wrk[unal]}" \
            "${arr_wrk[improp]:-}" \
            "${arr_wrk[mapq]}" \
            "${fin_spl}" \
            || return 1
    fi


    # Step 5: write the pipeline record into the header, including whether
    # optical duplicates were told apart, and then index the BAM work file,
    # unless Step 6 converts it to CRAM and indexes that instead. Under
    # '--unal_mode separate', each diagnostic file carries the same record,
    # plus what it holds and that Markdup 2 marked it. In each, the optical
    # line ends the content lines.
    if [[ "${unal_mode}" == "separate" ]]; then
        lin_opt_all="$(_print_optical_line "${dist_all}")"
        arr_cnt=(
            [unal]="unaligned reads and their mates, aligned or not"
            [improp]="pairs with both mates aligned that --req_flg rejects"
            [mapq]="templates that --mapq rejects and --req_flg does not"
        )

        for kind in "${arr_kind[@]}"; do
            local -n arr_rec_kind="arr_rec_${kind}"
            # shellcheck disable=SC2034  # Read in '_write_align_record'.
            arr_rec_kind=(
                "${arr_rec_cnt[@]}"
                $'@CO\tpipeline: content='"${arr_cnt[${kind}]}"
                $'@CO\tpipeline: warning=not for signal or scaling'
                $'@CO\tpipeline: markdup=all files (Markdup 2)'
                "${lin_opt_all}"
            )
            unset -n arr_rec_kind
        done
    fi

    lin_opt="$(_print_optical_line "${dist_mrk}")"
    arr_rec_cnt+=( "${lin_opt}" )

    idx_rec=true
    if [[ "${out_fmt}" == "cram" ]]; then idx_rec=false; fi

    _write_align_record \
        "${FUNCNAME[0]}" \
        "${dry_run}" \
        "${threads}" \
        "${fil_wrk}" \
        arr_rec_cnt \
        arr_rec_run \
        "${idx_rec}" \
        || return 1

    if [[ "${unal_mode}" == "separate" ]]; then
        for kind in "${arr_kind[@]}"; do
            _write_align_record \
                "${FUNCNAME[0]}" \
                "${dry_run}" \
                "${threads}" \
                "${arr_wrk[${kind}]}" \
                "arr_rec_${kind}" \
                arr_rec_run \
                "${idx_rec}" \
                || return 1
        done

        _write_diag_tsv \
            "${FUNCNAME[0]}" \
            "${dry_run}" \
            "${threads}" \
            "${aligner}" \
            "${bt2_X}" \
            "${fq_2}" \
            "${mapq}" \
            "${fil_wrk}" \
            "${arr_wrk[unal]}" \
            "${arr_wrk[improp]:-}" \
            "${arr_wrk[mapq]}" \
            "${txt_all}.gz" \
            "${fil_tsv}" \
            || return 1
    fi


    # Step 6: convert to CRAM if requested; with '--qname', name-sort the final
    # BAM into a copy.
    _finalize_align_output \
        "${FUNCNAME[0]}" \
        "${dry_run}" \
        "${threads}" \
        "${out_fmt}" \
        "${fil_out}" \
        "${fil_wrk}" \
        "${ref_fa}" \
        "${qname}" \
        "${bam_qnam}" \
        "${fil_qnam_out}" \
        || return 1

    if [[ "${unal_mode}" == "separate" ]]; then
        for kind in "${arr_kind[@]}"; do
            _finalize_align_output \
                "${FUNCNAME[0]}" \
                "${dry_run}" \
                "${threads}" \
                "${out_fmt}" \
                "${arr_out[${kind}]}" \
                "${arr_wrk[${kind}]}" \
                "${ref_fa}" \
                "${qname}" \
                "${arr_bam_qnam[${kind}]}" \
                "${arr_qnam_out[${kind}]}" \
                false \
                || return 1
        done
    fi
}


# Print an error message when function script is executed directly.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
