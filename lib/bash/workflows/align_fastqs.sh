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
# _mark_dup_bam
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
        check_args check_source format_outputs || {
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

  Print a command to stdout for a dry run of 'align_fastqs()': each named command array in turn as a pipeline, then ' > file' when '--out' is given, then a blank line. The output is valid shell: words are shell-escaped, and each option or option-value pair is on its own continuation line.

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

  - The program name, with its subcommand for 'samtools', 'bwa', and 'bwa-mem2', starts the command; 'mv', 'rm', 'cut', 'uniq', and 'awk' also keep their leading options there, as in 'mv -f'. An option shares a line with its value only when the option is listed here as taking one, since shape alone cannot tell 'samtools markdup -t in.bam' (a flag, then a file) from 'bwa aln -t 4' (an option and its value). A new value-taking option belongs in this list; without it, the option and its value print on separate lines, which changes the layout but not the command.

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
            bowtie2)                  opt_val="-p -x -1 -2 -U" ;;
            "samtools view")          opt_val="-@ -F -f -q -o -T -O -N" ;;
            samtools\ *)              opt_val="-@ -o" ;;
            "bwa mem"|"bwa-mem2 mem") opt_val="-t" ;;
            "bwa aln")                opt_val="-t" ;;
            awk)                      opt_val="-v" ;;
            cut)                      opt_val="-f" ;;
            *)                        opt_val="" ;;
        esac

        # A short utility keeps its leading options, with their values, on the
        # first line, as in 'mv -f' or 'uniq -c'; only its operands go on
        # continuation lines.
        case "${arr_prn_ref[0]}" in
            mv|rm|cut|uniq|awk)
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


# Keep or drop each template (all records sharing a read name) of a SAM or BAM
# file queryname-grouped as a whole, by whether its primary records pass a MAPQ
# threshold; write the kept records to a BAM file.
function _filter_bam_mapq_template() {
    local bam_in="${1:-}"
    local bam_out="${2:-}"
    local mapq="${3:-}"
    local n_pri="${4:-}"
    local threads="${5:-1}"
    local func="${6:-${FUNCNAME[1]:-${FUNCNAME[0]}}}"
    local dry_run="${7:-false}"
    local txt_nam="${bam_out%.bam}.names.txt"
    local -a arr_view arr_cut arr_uniq arr_awk arr_cmd
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _filter_bam_mapq_template
    [--help] bam_in bam_out mapq n_pri [threads] [func] [dry_run]

  Keep or drop each template (all records sharing a read name) of a queryname-grouped SAM or BAM file as a whole, by whether its primary records pass a MAPQ threshold, and write the kept records to a BAM file.

  Secondary and supplementary records follow their template, whatever their own MAPQ.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  bam_in : file
    Queryname-grouped SAM or BAM file to filter.

  2  bam_out : file
    BAM file to write the kept records to.

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

Returns
-------
  Writes 'bam_out' and returns 0 on success or prints the commands and returns 0 for a dry run; prints an error and returns 1 on failure.

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

    # List the names whose primary records, 'n_pri' of them, all pass.
    arr_view=(
        samtools view -@ "${threads}" -F 0x900 -q "${mapq}" "${bam_in}"
    )
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
    # records included.
    arr_cmd=(
        samtools view
            -@ "${threads}"
            -b
            -N "${txt_nam}"
            -o "${bam_out}"
            "${bam_in}"
    )

    if ! \
        _run_or_print_cmd "${dry_run}" arr_cmd
    then
        rm -f "${txt_nam}" "${bam_out}"
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
    local bwa_alg="${6:-}"
    local bwa_alg_set="${7:-}"
    local mapq="${8:-}"
    local req_flg="${9:-}"
    local index="${10:-}"
    local ref_fa="${11:-}"
    local fq_1="${12:-}"
    local fq_2="${13:-}"
    local fil_out="${14:-}"
    local out_fmt
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _validate_args_align_fastqs
    [--help] func threads aligner bt2_mode bt2_mode_set bwa_alg bwa_alg_set mapq req_flg index ref_fa fq_1 fq_2 fil_out

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

  06  bwa_alg : str
    BWA algorithm: 'mem' or 'aln'.

  07  bwa_alg_set : bool
    Whether the caller supplied '--bwa_alg'.

  08  mapq : int
    MAPQ threshold.

  09  req_flg : bool
    Whether properly paired alignments are required.

  10  index : path
    Path to the aligner index/reference.

  11  ref_fa : file
    Reference FASTA path; required for CRAM output, otherwise ignored with a warning.

  12  fq_1 : file
    First FASTQ input file.

  13  fq_2 : file
    Second FASTQ input file, or empty for single-end data.

  14  fil_out : file
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
        mem \\
        false \\
        1 \\
        false \\
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
        mem \\
        false \\
        30 \\
        true \\
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
}


# Step 1 with Bowtie 2: align paired- or single-end reads and pipe them to the
# Samtools call that the caller built.
function _align_bowtie2() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local bt2_mode="${4:-}"
    local index="${5:-}"
    local fq_1="${6:-}"
    local fq_2="${7:-}"
    local -a arr_bt2 arr_sam
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _align_bowtie2
    [--help] func dry_run threads bt2_mode index fq_1 fq_2 [samtools_arg ...]

  Align paired- or single-end reads with Bowtie 2 and pipe the alignments to 'samtools view' with the caller's arguments (Step 1 of 'align_fastqs()').

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

  4  bt2_mode : str
    Bowtie 2 alignment type: 'local', 'global', or 'end-to-end'.

  5  index : path
    Bowtie 2 index stem.

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
        tests/fixtures/align_fastqs/bowtie2/tiny \\
        tests/fixtures/align_fastqs/fastq/se/tiny_se.atria.fastq.gz \\
        '' \\
        -@ 1 -F 4 -q 1 -o work/tiny_se.bam
    '''

  2. Print, without running, the paired-end call in local mode that requires proper pairs.
    '''bash
    _align_bowtie2 \\
        align_fastqs \\
        true \\
        2 \\
        local \\
        tests/fixtures/align_fastqs/bowtie2/tiny \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R1.atria.fastq.gz \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        -@ 2 -F 12 -f 2 -q 1 -o work/tiny_pe.bam
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    # Assign Bowtie 2 arguments.
    arr_bt2=( bowtie2 -p "${threads}" -x "${index}" )

    if [[ "${bt2_mode}" == "local" ]]; then
        arr_bt2+=( --very-sensitive-local )
    elif [[ "${bt2_mode}" =~ ^(end-to-end|global)$ ]]; then
        arr_bt2+=( --very-sensitive )
    fi

    # Exclude unaligned reads from output.
    arr_bt2+=( --no-unal --phred33 )

    # Run Bowtie 2 with paired-end (top) or single-end (bottom) sequenced reads.
    if [[ -n "${fq_2}" ]]; then
        arr_bt2+=( --no-mixed --no-discordant )
        arr_bt2+=( --no-overlap --no-dovetail )
        arr_bt2+=( -1 "${fq_1}" -2 "${fq_2}" )
    else
        arr_bt2+=( -U "${fq_1}" )
    fi

    arr_sam=( samtools view "${@:8}" )

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

    # Assign BWA or BWA-MEM2 arguments.
    arr_bwa=( "${prog}" mem -t "${threads}" "${index}" "${fq_1}" )

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


# Step 1: build the Samtools call that drops unaligned reads (and applies the
# MAPQ filter for Bowtie 2), then align with the selected aligner.
function _align_reads() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local aligner="${4:-}"
    local bt2_mode="${5:-}"
    local bwa_alg="${6:-}"
    local mapq="${7:-}"
    local req_flg="${8:-}"
    local index="${9:-}"
    local fq_1="${10:-}"
    local fq_2="${11:-}"
    local fil_wrk="${12:-}"
    local -a arr_sam
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _align_reads
    [--help] func dry_run threads aligner bt2_mode bwa_alg mapq req_flg index fq_1 fq_2 fil_wrk

  Build the 'samtools view' arguments that drop unaligned reads (and, for Bowtie 2, apply the MAPQ threshold), then align with the selected aligner and write the BAM work file (Step 1 of 'align_fastqs()').

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

  06  bwa_alg : str
    BWA algorithm: 'mem' or 'aln'.

  07  mapq : int
    MAPQ threshold; applied here only for Bowtie 2.

  08  req_flg : bool
    Whether properly paired alignments are required.

  09  index : path
    Path to the aligner index/reference.

  10  fq_1 : file
    First FASTQ input file.

  11  fq_2 : file
    Second FASTQ input file, or empty for single-end data.

  12  fil_wrk : file
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
        mem \\
        1 \\
        true \\
        tests/fixtures/align_fastqs/bowtie2/tiny \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R1.atria.fastq.gz \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        work/tiny_pe.bam
    '''

  2. Print, without running, the single-end calls with BWA-backtrack; MAPQ is filtered later.
    '''bash
    _align_reads \\
        align_fastqs \\
        true \\
        1 \\
        bwa \\
        end-to-end \\
        aln \\
        30 \\
        false \\
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

    # Based on parsed arguments, construct the call to Samtools. BWA and
    # BWA-MEM2 write unaligned reads, so use Samtools to drop them, including
    # the mates of unaligned reads in PE data (i.e., keep pairs whole).
    arr_sam=( -@ "${threads}" )
    if [[ -n "${fq_2}" ]]; then
        arr_sam+=( -F 12 )
        if [[ "${req_flg}" == "true" ]]; then arr_sam+=( -f 2 ); fi
    else
        arr_sam+=( -F 4 )
    fi

    # Bowtie 2 gives both mates of a pair one MAPQ and, without '-k', writes no
    # supplementary or secondary records, so a per-record MAPQ filter keeps or
    # drops whole templates. BWA and BWA-MEM2 do neither, so their MAPQ filter
    # runs via read alignment name in Step 2.
    if [[ "${aligner}" == "bowtie2" ]]; then arr_sam+=( -q "${mapq}" ); fi
    arr_sam+=( -o "${fil_wrk}" )

    # Based on parsed arguments, run Bowtie 2, BWA, or BWA-MEM2 with output
    # piped to the above-constructed Samtools call.
    if [[ "${aligner}" == "bowtie2" ]]; then
        _align_bowtie2 \
            "${func}" "${dry_run}" "${threads}" "${bt2_mode}" "${index}" \
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
# paired-end alignments, and, for BWA and BWA-MEM2, filter templates by MAPQ.
function _sort_qname_bam() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local aligner="${4:-}"
    local mapq="${5:-}"
    local fq_2="${6:-}"
    local fil_wrk="${7:-}"
    local bam_qnam="${8:-}"
    local bam_mate bam_mapq n_pri
    local -a arr_cmd
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _sort_qname_bam
    [--help] func dry_run threads aligner mapq fq_2 fil_wrk bam_qnam

  Sort the BAM work file by queryname into 'bam_qnam', fix mates for paired-end data, and, for BWA and BWA-MEM2 with 'mapq' above 0, keep or drop each template whole by its primary records' MAPQ (Step 2 of 'align_fastqs()').

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

  4  aligner : str
    Alignment program that wrote the BAM work file.

  5  mapq : int
    MAPQ threshold; applied here only for BWA and BWA-MEM2.

  6  fq_2 : file
    Second FASTQ input file, or empty for single-end data; selects mate fixing.

  7  fil_wrk : file
    BAM work file from Step 1.

  8  bam_qnam : file
    Queryname-sorted BAM file to write.

Returns
-------
  Returns 0 when sorting, mate fixing, and filtering succeed or after printing the commands for a dry run; prints an error and returns 1 otherwise.

Notes
-----
  Runtime requirements:
    - awk (when 'aligner' is 'bwa' or 'bwa-mem2' and 'mapq' is above 0)
    - bash >= 4.4
    - mv
    - samtools

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
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        work/tiny_pe.bam \\
        work/tiny_pe.qnam.bam
    '''

  2. Print, without running, the queryname sort of a single-end Bowtie 2 work file; MAPQ was filtered in Step 1.
    '''bash
    _sort_qname_bam \\
        align_fastqs \\
        true \\
        1 \\
        bowtie2 \\
        1 \\
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
        printf '\n%s\n' \
            "# Step 2: sort the BAM work file by queryname."
    fi

    # A dry run writes no Step 1 file, so it prints the commands as they would
    # run.
    if [[ "${dry_run}" == "true" || -f "${fil_wrk}" ]]; then
        arr_cmd=(
            samtools sort
                -@ "${threads}"
                -n
                -o "${bam_qnam}"
                "${fil_wrk}"
        )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #2: failed to queryname-sort ${aligner}-aligned BAM" \
                "file."
            return 1
        fi

        if [[ -n "${fq_2}" ]]; then
            bam_mate="${fil_wrk%.bam}.mate.bam"

            arr_cmd=(
                samtools fixmate
                    -@ "${threads}"
                    -c
                    -m
                    "${bam_qnam}"
                    "${bam_mate}"
            )

            if ! \
                _run_or_print_cmd "${dry_run}" arr_cmd
            then
                echo_err_func "${func}" \
                    "Step #2: failed to fix mate pairs in queryname-sorted" \
                    "BAM file."
                return 1
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

        # BWA and BWA-MEM2 can give two mates different MAPQ values and write
        # supplementary records with their own MAPQ values, so keep or drop
        # each template (all records sharing a read name, secondary records
        # included) as a whole, by whether its primary records (both mates in
        # paired-end data) pass the '--mapq' gate.
        if [[ "${aligner}" != "bowtie2" && "${mapq}" -gt 0 ]]; then
            bam_mapq="${fil_wrk%.bam}.mapq.bam"
            if [[ -n "${fq_2}" ]]; then n_pri=2; else n_pri=1; fi

            if ! \
                _filter_bam_mapq_template \
                    "${bam_qnam}" \
                    "${bam_mapq}" \
                    "${mapq}" \
                    "${n_pri}" \
                    "${threads}" \
                    "${func}" \
                    "${dry_run}"
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
                    "Step #2: failed to rename MAPQ-filtered," \
                    "queryname-sorted BAM file."
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
# into the BAM work file, then index it.
function _sort_coord_bam() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local qname="${4:-}"
    local fq_2="${5:-}"
    local fil_wrk="${6:-}"
    local bam_qnam="${7:-}"
    local bam_coor
    local -a arr_cmd
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _sort_coord_bam
    [--help] func dry_run threads qname fq_2 fil_wrk bam_qnam

  Sort the queryname-sorted BAM file by coordinates into the BAM work file and index it, deleting the queryname-sorted file unless it is to be retained (Step 3 of 'align_fastqs()').

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

  4  qname : bool
    Whether to retain the queryname-sorted file.

  5  fq_2 : file
    Second FASTQ input file, or empty for single-end data; selects the error message.

  6  fil_wrk : file
    BAM work file to replace with the coordinate-sorted alignments.

  7  bam_qnam : file
    Queryname-sorted BAM file from Step 2.

Returns
-------
  Returns 0 when sorting and indexing succeed or after printing the commands for a dry run; prints an error and returns 1 otherwise.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - mv
    - rm
    - samtools

Examples
--------
  1. Coordinate-sort a single-end work file and delete the queryname-sorted file.
    '''bash
    _sort_coord_bam \\
        align_fastqs \\
        false \\
        2 \\
        false \\
        '' \\
        work/tiny_se.bam \\
        work/tiny_se.qnam.bam
    '''

  2. Print, without running, the coordinate sort of a paired-end work file that retains the queryname-sorted file.
    '''bash
    _sort_coord_bam \\
        align_fastqs \\
        true \\
        2 \\
        true \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        work/tiny_pe.bam \\
        work/tiny_pe.qnam.bam
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ "${dry_run}" == "true" ]]; then
        printf '\n%s\n' \
            "# Step 3: sort the BAM work file by coordinates and index it."
    fi

    if [[ "${dry_run}" == "true" || -f "${bam_qnam}" ]]; then
        bam_coor="${fil_wrk%.bam}.coor.bam"
        arr_cmd=(
            samtools sort
                -@ "${threads}"
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

        if [[ "${qname}" == "false" ]]; then
            arr_cmd=( rm -f "${bam_qnam}" )

            if ! \
                _run_or_print_cmd "${dry_run}" arr_cmd
            then
                echo_err_func "${func}" \
                    "Step #3: failed to delete intermediate queryname-sorted" \
                    "BAM file."
                return 1
            fi
        fi

        arr_cmd=( mv -f "${bam_coor}" "${fil_wrk}" )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #3: failed to rename coordinate-sorted BAM file."
            return 1
        fi

        arr_cmd=( samtools index -@ "${threads}" "${fil_wrk}" )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #3: failed to index renamed, coordinate-sorted BAM file."
            return 1
        fi
    else
        if [[ -n "${fq_2}" ]]; then
            echo_err_func "${func}" \
                "Step #3: mate-fixed, queryname-sorted BAM file does not" \
                "exist."
        else
            echo_err_func "${func}" \
                "Step #3: queryname-sorted BAM file does not exist."
        fi
        return 1
    fi
}


# Step 4: mark duplicate alignments in the coordinate-sorted BAM work file,
# then index the file again.
function _mark_dup_bam() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local fil_wrk="${4:-}"
    local bam_mrk
    local -a arr_cmd
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _mark_dup_bam
    [--help] func dry_run threads fil_wrk

  Mark duplicate alignments in the coordinate-sorted BAM work file and index it again, skipping files already marked (Step 4 of 'align_fastqs()').

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

Returns
-------
  Returns 0 when marking and indexing succeed or are skipped or after printing the commands for a dry run; prints an error and returns 1 otherwise.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - grep
    - mv
    - samtools

Examples
--------
  1. Mark duplicates in a coordinate-sorted work file.
    '''bash
    _mark_dup_bam align_fastqs false 2 work/tiny_pe.bam
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
        printf '\n%s\n' \
            "# Step 4: mark duplicates and index the BAM work file again."
    fi

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

            arr_cmd=(
                samtools markdup
                    -@ "${threads}"
                    -t
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

            # Index the duplicate-marked coordinate-sorted BAM file.
            arr_cmd=( samtools index -@ "${threads}" "${fil_wrk}" )

            if ! \
                _run_or_print_cmd "${dry_run}" arr_cmd
            then
                echo_err_func "${func}" \
                    "Step #4: failed to index duplicate-marked," \
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


# Step 5: convert the final BAM work file to CRAM if requested; otherwise
# retain BAM output as is, renaming any retained queryname-sorted file.
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
    local -a arr_cmd
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _finalize_align_output
    [--help] func dry_run threads out_fmt fil_out fil_wrk ref_fa qname bam_qnam fil_qnam_out

  Convert the BAM work file and any retained queryname-sorted file to CRAM when requested, or rename the retained queryname-sorted BAM file into place (Step 5 of 'align_fastqs()').

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  01  func : str
    Name of the calling function for diagnostics.

  02  dry_run : bool
    Whether to print the commands to stdout instead of running them.

  03  threads : int
    Number of threads to use for indexing.

  04  out_fmt : str
    Final output format: 'bam' or 'cram'.

  05  fil_out : file
    Path to the final alignment file.

  06  fil_wrk : file
    BAM work file from Step 4.

  07  ref_fa : file
    Reference FASTA path, used for CRAM output.

  08  qname : bool
    Whether the queryname-sorted file is retained.

  09  bam_qnam : file
    Queryname-sorted BAM work file.

  10  fil_qnam_out : file
    Path for the retained queryname-sorted output.

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
  1. Convert a paired-end work file and its retained queryname-sorted file to CRAM.
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

  2. Print, without running, the step for BAM output that keeps the queryname-sorted file in place.
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

    if [[ "${dry_run}" == "true" ]]; then
        printf '\n%s\n' \
            "# Step 5: write the final output and the retained files."

        # BAM output is the BAM work file itself; with '--qname', so is the
        # retained file, so nothing runs.
        if [[
            "${out_fmt}" == "bam"
                && ( "${qname}" == "false" || "${bam_qnam}" == "${fil_qnam_out}" )
        ]]; then
            printf '%s\n\n' \
                "# Nothing to do: the BAM work file is the final output."
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
                "Step #5: failed to convert BAM work file to final CRAM file."
            return 1
        fi

        arr_cmd=( samtools index -@ "${threads}" "${fil_out}" )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #5: failed to index final CRAM file."
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
                    "Step #5: failed to convert queryname-sorted BAM" \
                    "intermediate to retained CRAM file."
                return 1
            fi
        fi

        arr_cmd=( rm -f "${fil_wrk}" "${fil_wrk}.bai" )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #5: failed to delete BAM work file and/or BAM index" \
                "after CRAM conversion."
            return 1
        fi

        if [[ "${qname}" == "true" ]]; then
            arr_cmd=( rm -f "${bam_qnam}" )

            if ! \
                _run_or_print_cmd "${dry_run}" arr_cmd
            then
                echo_err_func "${func}" \
                    "Step #5: failed to delete queryname-sorted BAM work" \
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
                    "Step #5: failed to rename retained queryname-sorted BAM" \
                    "intermediate."
                return 1
            fi
        fi
    fi
}


function align_fastqs() {
    local dry_run threads aligner bt2_mode bwa_alg mapq req_flg index ref_fa
    local bt2_mode_set bwa_alg_set
    local fq_1 fq_2 fil_out qname out_fmt
    local fil_wrk fil_qnam_out bam_qnam
    local show_help

    # Assign default argument values.
    dry_run=false
    threads=1
    aligner="bowtie2"
    bt2_mode="end-to-end"
    bt2_mode_set=false
    bwa_alg="mem"
    bwa_alg_set=false
    ref_fa=""
    mapq=1
    req_flg=false
    fq_1=""
    fq_2=""
    fil_out=""
    qname=false

    show_help=$(cat << EOM
Usage
-----
  align_fastqs
    [--help] [--dry_run] [--threads <int>] [--aligner <aligner>] [--bt2_mode <mode>] [--bwa_alg <algorithm>] [--mapq <int>] [--req_flg] --index <path> [--ref_fa <file>] --fq_1 <file> [--fq_2 <file>] --fil_out <file> [--qname]

  Align single- or paired-end Illumina short reads using Bowtie 2, BWA, or BWA-MEM2, and convert the output to a sorted, duplicate-marked, mate-fixed (if working with paired-end reads) alignment file using Samtools.

  Final output format is inferred from '--fil_out': '.bam' or '.cram'. Internal post-processing is performed in BAM format. If CRAM output is requested, the final BAM work file and any retained queryname-sorted BAM intermediate are converted to CRAM in the final step.

  When '--aligner bwa' is selected, '--bwa_alg' may be used to choose either 'mem' or the older backtrack workflow ('aln' with downstream 'samse' or 'sampe').

  Optionally, a queryname-sorted output can be retained with the '--qname' flag; otherwise, the queryname-sorted intermediate/work file is deleted.

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

  -ba, --bwa_alg : {'mem', 'aln'}
    BWA algorithm when '--aligner bwa': 'mem' or 'aln' (default: '${bwa_alg}').

    With '--aligner bwa-mem2', which has only 'mem', 'mem' is ignored with a warning and 'aln' is refused; with '--aligner bowtie2', it is refused.

  -mq, --mpq, --mapq : int
    MAPQ threshold for filtering alignment output files (default: ${mapq}).

    To disable MAPQ-based filtering, specify 0.

    With paired-end data, a pair is kept only when both mates pass. Secondary (SAM flag 256) and supplementary (SAM flag 2048) alignments are kept or dropped with their associated primary read alignments, whatever their own MAPQ. With Bowtie 2, a MAPQ of 1 or more does not mean a unique alignment, as a read alignment with multiple equally good placements can get MAPQ 1, and its placement is chosen at random.

  -rq, --req_flg : flag
    Require SAM flag bit 2 for properly paired alignments. Ignored with a warning for single-end data.

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
    Retain queryname-sorted intermediate alignment files.

Returns
-------
  Returns 0 when alignment and post-processing complete successfully; otherwise 1.

Output
------
  Creates a final '.bam' or '.cram' file at '--fil_out'. A corresponding index ('.bai' or '.crai') is also created.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - bowtie2 (when '--aligner bowtie2' is specified)
    - bwa (when '--aligner bwa' is specified)
    - bwa-mem2 (when '--aligner bwa-mem2' is specified)
    - grep
    - Reference FASTA and required index (when writing CRAM)
    - samtools

  - For '--index' when using Bowtie 2, the path should end with the index stem: 'path/to/dir/stem'; when using 'bwa' or 'bwa-mem2', the path should be the indexed reference FASTA path (for example, '.fa').
  - '--bwa_alg' is used only when '--aligner bwa'. It is refused with '--aligner bowtie2'; with '--aligner bwa-mem2', which always runs 'bwa-mem2 mem', 'mem' is ignored with a warning and 'aln' is refused.
  - Unaligned reads are excluded from the output; with paired-end data, so are the mates of unaligned reads.
  - Duplicates are marked (SAM flag 1024) in Step #4, after MAPQ filtering, and kept. Like any other read, '--mapq' filters them by the MAPQ the aligner gave them.
  - For '--qname', the retained queryname-sorted output will have the same path and stem assigned to '--fil_out', except '.qnam' will be inserted before the final extension (for example, '.qnam.bam' or '.qnam.cram'). For CRAM output, the retained queryname-sorted BAM work file is converted in Step #5.
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
        "${bwa_alg}" \
        "${bwa_alg_set}" \
        "${mapq}" \
        "${req_flg}" \
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


    # Step 1: align the reads into the BAM work file, dropping unaligned reads
    # and, for Bowtie 2, applying the MAPQ threshold.
    _align_reads \
        "${FUNCNAME[0]}" \
        "${dry_run}" \
        "${threads}" \
        "${aligner}" \
        "${bt2_mode}" \
        "${bwa_alg}" \
        "${mapq}" \
        "${req_flg}" \
        "${index}" \
        "${fq_1}" \
        "${fq_2}" \
        "${fil_wrk}" \
        || return 1


    # Step 2: sort by queryname, fix mates in paired-end data, and, for BWA and
    # BWA-MEM2, keep or drop each template whole by its primary records' MAPQ.
    _sort_qname_bam \
        "${FUNCNAME[0]}" \
        "${dry_run}" \
        "${threads}" \
        "${aligner}" \
        "${mapq}" \
        "${fq_2}" \
        "${fil_wrk}" \
        "${bam_qnam}" \
        || return 1


    # Step 3: sort by coordinates into the BAM work file and index it.
    _sort_coord_bam \
        "${FUNCNAME[0]}" \
        "${dry_run}" \
        "${threads}" \
        "${qname}" \
        "${fq_2}" \
        "${fil_wrk}" \
        "${bam_qnam}" \
        || return 1


    # Step 4: mark duplicates and index the BAM work file again.
    _mark_dup_bam \
        "${FUNCNAME[0]}" \
        "${dry_run}" \
        "${threads}" \
        "${fil_wrk}" \
        || return 1


    # Step 5: convert to CRAM if requested, and keep or delete the
    # queryname-sorted file.
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
}


# Print an error message when function script is executed directly.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
