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


# _check_path_whitespace
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


# MAYBE: move to 'check_inputs.sh'?
# MAYBE: make "public"?
function _check_path_whitespace() {
    local path="${1:-}"
    local nam_arg="${2:-path}"
    local func="${3:-${FUNCNAME[1]:-${FUNCNAME[0]}}}"

    if [[ "${path}" =~ [[:space:]] ]]; then
        echo_err_func "${func}" \
            "'${nam_arg}' contains whitespace, which is not supported by" \
            "'align_fastqs()' because Step #1 currently constructs" \
            "string-built argument bundles."
        return 1
    fi
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
    local txt_nam="${bam_out%.bam}.names.txt"
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _filter_bam_mapq_template
    [--help] bam_in bam_out mapq n_pri [threads] [func]

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

Returns
-------
  Writes 'bam_out' and returns 0 on success; prints an error and returns 1 on failure.

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
    _filter_bam_mapq_template sample.qnam.bam sample.mapq.bam 30 2 4
    '''

  2. Keep single-end reads whose primary record reaches MAPQ 1.
    '''bash
    _filter_bam_mapq_template sample.qnam.bam sample.mapq.bam 1 1
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
    if ! (
        set -o pipefail
        samtools view -@ "${threads}" -F 0x900 -q "${mapq}" "${bam_in}" \
            | cut -f 1 \
            | uniq -c \
            | awk -v n="${n_pri}" '$1 == n { print $2 }' \
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
    if ! \
        samtools view \
            -@ "${threads}" \
            -b \
            -N "${txt_nam}" \
            -o "${bam_out}" \
            "${bam_in}"
    then
        rm -f "${txt_nam}" "${bam_out}"
        echo_err_func "${func}" \
            "failed to filter '${bam_in}' by read name."
        return 1
    fi

    rm -f "${txt_nam}"
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

  1  func : str
    Name of the calling function for diagnostics.

  2  threads : int
    Number of threads to use.

  3  aligner : str
    Alignment program: 'bowtie2', 'bwa', or 'bwa-mem2'.

  4  bt2_mode : str
    Bowtie 2 alignment type: 'local', 'global', or 'end-to-end'.

  5  bt2_mode_set : bool
    Whether the caller supplied '--bt2_mode'.

  6  bwa_alg : str
    BWA algorithm: 'mem' or 'aln'.

  7  bwa_alg_set : bool
    Whether the caller supplied '--bwa_alg'.

  8  mapq : int
    MAPQ threshold.

  9  req_flg : bool
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
    _validate_args_align_fastqs align_fastqs 1 bowtie2 end-to-end false \\
        mem false 1 false tests/fixtures/align_fastqs/bowtie2/tiny '' \\
        tests/fixtures/align_fastqs/fastq/se/tiny_se.atria.fastq.gz '' \\
        work/tiny_se.bam
    '''

  2. Validate a paired-end BWA-MEM2 call that writes CRAM.
    '''bash
    _validate_args_align_fastqs align_fastqs 4 bwa-mem2 end-to-end false \\
        mem false 30 true tests/fixtures/align_fastqs/bwa-mem2/tiny.fa \\
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

    _check_path_whitespace \
        "${index}" "--index" "${func}" \
        || return 1

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

    _check_path_whitespace \
        "${fil_out}" "--fil_out" "${func}" \
        || return 1

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
                "'--fil_out' must end in '.bam' or '.cram', but was" \
                "assigned '${fil_out}'."
            return 1
            ;;
    esac

    if [[ "${out_fmt}" == "cram" ]]; then
        if [[ -z "${ref_fa}" ]]; then
            echo_err_func "${func}" \
                "'--ref_fa' is required when '--fil_out' ends in '.cram'."
            return 1
        fi

        _check_path_whitespace \
            "${ref_fa}" "--ref_fa" "${func}" \
            || return 1

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

    _check_path_whitespace \
        "${fq_1}" "--fq_1" "${func}" \
        || return 1

    if [[ -n "${fq_2}" ]]; then
        if [[ ! -f "${fq_2}" ]]; then
            echo_err_func "${func}" \
                "file associated with '--fq_2' does not exist: '${fq_2}'."
            return 1
        fi

        _check_path_whitespace \
            "${fq_2}" "--fq_2" "${func}" \
            || return 1
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
    local threads="${2:-}"
    local bt2_mode="${3:-}"
    local index="${4:-}"
    local fq_1="${5:-}"
    local fq_2="${6:-}"
    local args_sam="${7:-}"
    local args_bt2
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _align_bowtie2
    [--help] func threads bt2_mode index fq_1 fq_2 args_sam

  Align paired- or single-end reads with Bowtie 2 and pipe the alignments to 'samtools view' with the caller's arguments (Step 1 of 'align_fastqs()').

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  threads : int
    Number of threads to use.

  3  bt2_mode : str
    Bowtie 2 alignment type: 'local', 'global', or 'end-to-end'.

  4  index : path
    Bowtie 2 index stem.

  5  fq_1 : file
    First FASTQ input file.

  6  fq_2 : file
    Second FASTQ input file, or empty for single-end data.

  7  args_sam : str
    Arguments for 'samtools view', ending with its output path.

Returns
-------
  Returns 0 when alignment and filtering succeed; prints an error and returns 1 otherwise.

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
    _align_bowtie2 align_fastqs 1 end-to-end \\
        tests/fixtures/align_fastqs/bowtie2/tiny \\
        tests/fixtures/align_fastqs/fastq/se/tiny_se.atria.fastq.gz '' \\
        '-@ 1 -F 4 -q 1 -o work/tiny_se.bam'
    '''

  2. Align the paired-end fixtures in local mode, requiring proper pairs.
    '''bash
    _align_bowtie2 align_fastqs 2 local \\
        tests/fixtures/align_fastqs/bowtie2/tiny \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R1.atria.fastq.gz \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        '-@ 2 -F 12 -f 2 -q 1 -o work/tiny_pe.bam'
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    # Assign Bowtie 2 arguments.
    args_bt2="-p ${threads} -x ${index}"

    if [[ "${bt2_mode}" == "local" ]]; then
        args_bt2+=" --very-sensitive-local"
    elif [[ "${bt2_mode}" =~ ^(end-to-end|global)$ ]]; then
        args_bt2+=" --very-sensitive"
    fi

    # Exclude unaligned reads from output.
    args_bt2+=" --no-unal --phred33"

    # Run Bowtie 2 with paired-end (top) or single-end (bottom) sequenced reads.
    if [[ -n "${fq_2}" ]]; then
        args_bt2+=" --no-mixed --no-discordant"
        args_bt2+=" --no-overlap --no-dovetail"
        args_bt2+=" -1 ${fq_1} -2 ${fq_2}"
    else
        args_bt2+=" -U ${fq_1}"
    fi

    # shellcheck disable=SC2086,SC2046
    if ! (
        set -o pipefail
        bowtie2 ${args_bt2} | samtools view ${args_sam}
    ); then
        if [[ -n "${fq_2}" ]]; then
            echo_err_func "${func}" \
                "Step #1: failed to align paired-end data with" \
                "'bowtie2' ('${bt2_mode}'), and/or failed to process" \
                "the alignments with 'samtools view'."
            return 1
        else
            echo_err_func "${func}" \
                "Step #1: failed to align single-end data with" \
                "'bowtie2' ('${bt2_mode}'), and/or failed to process" \
                "the alignments with 'samtools view'."
            return 1
        fi
    fi
}


# Step 1 with BWA-MEM or BWA-MEM2: align paired- or single-end reads and pipe
# them to the Samtools call that the caller built, which excludes unaligned
# reads.
function _align_bwa_mem() {
    local func="${1:-}"
    local threads="${2:-}"
    local prog="${3:-}"
    local index="${4:-}"
    local fq_1="${5:-}"
    local fq_2="${6:-}"
    local args_sam="${7:-}"
    local args_bwa
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _align_bwa_mem
    [--help] func threads prog index fq_1 fq_2 args_sam

  Align paired- or single-end reads with 'bwa mem' or 'bwa-mem2 mem' and pipe the alignments to 'samtools view' with the caller's arguments (Step 1 of 'align_fastqs()').

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  threads : int
    Number of threads to use.

  3  prog : str
    Aligner program: 'bwa' or 'bwa-mem2'.

  4  index : file
    Indexed reference FASTA path.

  5  fq_1 : file
    First FASTQ input file.

  6  fq_2 : file
    Second FASTQ input file, or empty for single-end data.

  7  args_sam : str
    Arguments for 'samtools view', ending with its output path.

Returns
-------
  Returns 0 when alignment and filtering succeed; prints an error and returns 1 otherwise.

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
    _align_bwa_mem align_fastqs 2 bwa-mem2 \\
        tests/fixtures/align_fastqs/bwa-mem2/tiny.fa \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R1.atria.fastq.gz \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        '-@ 2 -F 12 -o work/tiny_pe.bam'
    '''

  2. Align the single-end fixture with BWA-MEM.
    '''bash
    _align_bwa_mem align_fastqs 1 bwa tests/fixtures/align_fastqs/bwa/tiny.fa \\
        tests/fixtures/align_fastqs/fastq/se/tiny_se.atria.fastq.gz '' \\
        '-@ 1 -F 4 -o work/tiny_se.bam'
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    # Assign BWA or BWA-MEM2 arguments.
    args_bwa="-t ${threads} ${index}"

    # shellcheck disable=SC2086,SC2046
    if [[ -n "${fq_2}" ]]; then
        if ! (
            set -o pipefail
            "${prog}" mem ${args_bwa} "${fq_1}" "${fq_2}" \
                | samtools view ${args_sam}
        ); then
            echo_err_func "${func}" \
                "Step #1: failed to align paired-end data with" \
                "'${prog} mem', and/or failed to process the alignments" \
                "with 'samtools view'."
            return 1
        fi
    else
        if ! (
            set -o pipefail
            "${prog}" mem ${args_bwa} "${fq_1}" | samtools view ${args_sam}
        ); then
            echo_err_func "${func}" \
                "Step #1: failed to align single-end data with" \
                "'${prog} mem', and/or failed to process the alignments" \
                "with 'samtools view'."
            return 1
        fi
    fi
}


# Step 1 with BWA-backtrack: write '.sai' files with 'bwa aln', pair them with
# 'bwa sampe' or 'bwa samse', and pipe the alignments to the Samtools call that
# the caller built, which excludes unaligned reads.
function _align_bwa_aln() {
    local func="${1:-}"
    local threads="${2:-}"
    local index="${3:-}"
    local fq_1="${4:-}"
    local fq_2="${5:-}"
    local fil_wrk="${6:-}"
    local args_sam="${7:-}"
    local args_bwa_aln fil_sai_1 fil_sai_2
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _align_bwa_aln
    [--help] func threads index fq_1 fq_2 fil_wrk args_sam

  Align paired- or single-end reads with 'bwa aln' and 'bwa sampe' or 'bwa samse', and pipe the alignments to 'samtools view' with the caller's arguments (Step 1 of 'align_fastqs()').

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  threads : int
    Number of threads to use for 'bwa aln'.

  3  index : file
    Indexed reference FASTA path.

  4  fq_1 : file
    First FASTQ input file.

  5  fq_2 : file
    Second FASTQ input file, or empty for single-end data.

  6  fil_wrk : file
    BAM work file; the temporary '.sai' files are written beside it and deleted.

  7  args_sam : str
    Arguments for 'samtools view', ending with its output path.

Returns
-------
  Returns 0 when alignment and filtering succeed; prints an error and returns 1 otherwise.

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
    _align_bwa_aln align_fastqs 1 tests/fixtures/align_fastqs/bwa/tiny.fa \\
        tests/fixtures/align_fastqs/fastq/se/tiny_se.atria.fastq.gz '' \\
        work/tiny_se.bam '-@ 1 -F 4 -o work/tiny_se.bam'
    '''

  2. Align the paired-end fixtures with BWA-backtrack.
    '''bash
    _align_bwa_aln align_fastqs 2 tests/fixtures/align_fastqs/bwa/tiny.fa \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R1.atria.fastq.gz \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        work/tiny_pe.bam '-@ 2 -F 12 -o work/tiny_pe.bam'
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    # Assign BWA-backtrack arguments.
    args_bwa_aln="-t ${threads}"

    fil_sai_1="${fil_wrk%.bam}.R1.sai"
    fil_sai_2=""

    # shellcheck disable=SC2086,SC2046
    if ! \
        bwa aln ${args_bwa_aln} "${index}" "${fq_1}" > "${fil_sai_1}"
    then
        echo_err_func "${func}" \
            "Step #1: failed to generate first '.sai' file with" \
            "'bwa aln'."
        return 1
    fi

    # shellcheck disable=SC2086,SC2046
    if [[ -n "${fq_2}" ]]; then
        fil_sai_2="${fil_wrk%.bam}.R2.sai"

        if ! \
            bwa aln ${args_bwa_aln} "${index}" "${fq_2}" \
                > "${fil_sai_2}"
        then
            echo_err_func "${func}" \
                "Step #1: failed to generate second '.sai' file with" \
                "'bwa aln'."
            rm -f "${fil_sai_1}"
            return 1
        fi

        if ! (
            set -o pipefail
            bwa sampe \
                "${index}" "${fil_sai_1}" "${fil_sai_2}" \
                "${fq_1}" "${fq_2}" \
                | samtools view ${args_sam}
        ); then
            rm -f "${fil_sai_1}" "${fil_sai_2}"
            echo_err_func "${func}" \
                "Step #1: failed to align paired-end data with" \
                "'bwa aln'/'bwa sampe', and/or failed to process the" \
                "alignments with 'samtools view'."
            return 1
        fi

        if ! \
            rm -f "${fil_sai_1}" "${fil_sai_2}"
        then
            echo_err_func "${func}" \
                "Step #1: failed to delete temporary '.sai' files."
            return 1
        fi

    else
        if ! (
            set -o pipefail
            bwa samse "${index}" "${fil_sai_1}" "${fq_1}" \
                | samtools view ${args_sam}
        ); then
            rm -f "${fil_sai_1}"
            echo_err_func "${func}" \
                "Step #1: failed to align single-end data with" \
                "'bwa aln'/'bwa samse', and/or failed to process the" \
                "alignments with 'samtools view'."
            return 1
        fi

        if ! \
            rm -f "${fil_sai_1}"
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
    local threads="${2:-}"
    local aligner="${3:-}"
    local bt2_mode="${4:-}"
    local bwa_alg="${5:-}"
    local mapq="${6:-}"
    local req_flg="${7:-}"
    local index="${8:-}"
    local fq_1="${9:-}"
    local fq_2="${10:-}"
    local fil_wrk="${11:-}"
    local args_sam
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _align_reads
    [--help] func threads aligner bt2_mode bwa_alg mapq req_flg index fq_1 fq_2 fil_wrk

  Build the 'samtools view' arguments that drop unaligned reads (and, for Bowtie 2, apply the MAPQ threshold), then align with the selected aligner and write the BAM work file (Step 1 of 'align_fastqs()').

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  threads : int
    Number of threads to use.

  3  aligner : str
    Alignment program: 'bowtie2', 'bwa', or 'bwa-mem2'.

  4  bt2_mode : str
    Bowtie 2 alignment type: 'local', 'global', or 'end-to-end'.

  5  bwa_alg : str
    BWA algorithm: 'mem' or 'aln'.

  6  mapq : int
    MAPQ threshold; applied here only for Bowtie 2.

  7  req_flg : bool
    Whether properly paired alignments are required.

  8  index : path
    Path to the aligner index/reference.

  9  fq_1 : file
    First FASTQ input file.

  10  fq_2 : file
    Second FASTQ input file, or empty for single-end data.

  11  fil_wrk : file
    BAM work file to write.

Returns
-------
  Returns 0 when alignment succeeds; prints an error and returns 1 otherwise.

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
    _align_reads align_fastqs 2 bowtie2 end-to-end mem 1 true \\
        tests/fixtures/align_fastqs/bowtie2/tiny \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R1.atria.fastq.gz \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        work/tiny_pe.bam
    '''

  2. Align the single-end fixture with BWA-backtrack; MAPQ is filtered later.
    '''bash
    _align_reads align_fastqs 1 bwa end-to-end aln 30 false \\
        tests/fixtures/align_fastqs/bwa/tiny.fa \\
        tests/fixtures/align_fastqs/fastq/se/tiny_se.atria.fastq.gz '' \\
        work/tiny_se.bam
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    # Based on parsed arguments, construct the call to Samtools. BWA and
    # BWA-MEM2 write unaligned reads, so use Samtools to drop them, including
    # the mates of unaligned reads in PE data (i.e., keep pairs whole).
    args_sam="-@ ${threads}"
    if [[ -n "${fq_2}" ]]; then
        args_sam+=" -F 12"
        if [[ "${req_flg}" == "true" ]]; then args_sam+=" -f 2"; fi
    else
        args_sam+=" -F 4"
    fi

    # Bowtie 2 gives both mates of a pair one MAPQ and, without '-k', writes no
    # supplementary or secondary records, so a per-record MAPQ filter keeps or
    # drops whole templates. BWA and BWA-MEM2 do neither, so their MAPQ filter
    # runs via read alignment name in Step 2.
    if [[ "${aligner}" == "bowtie2" ]]; then args_sam+=" -q ${mapq}"; fi
    args_sam+=" -o ${fil_wrk}"

    # Based on parsed arguments, run Bowtie 2, BWA, or BWA-MEM2 with output
    # piped to the above-constructed Samtools call.
    if [[ "${aligner}" == "bowtie2" ]]; then
        _align_bowtie2 \
            "${func}" "${threads}" "${bt2_mode}" "${index}" \
            "${fq_1}" "${fq_2}" "${args_sam}" \
            || return 1
    elif [[ "${aligner}" == "bwa" ]]; then
        if [[ "${bwa_alg}" == "mem" ]]; then
            _align_bwa_mem \
                "${func}" "${threads}" bwa "${index}" \
                "${fq_1}" "${fq_2}" "${args_sam}" \
                || return 1
        elif [[ "${bwa_alg}" == "aln" ]]; then
            _align_bwa_aln \
                "${func}" "${threads}" "${index}" \
                "${fq_1}" "${fq_2}" "${fil_wrk}" "${args_sam}" \
                || return 1
        fi
    elif [[ "${aligner}" == "bwa-mem2" ]]; then
        _align_bwa_mem \
            "${func}" "${threads}" bwa-mem2 "${index}" \
            "${fq_1}" "${fq_2}" "${args_sam}" \
            || return 1
    fi
}


# Step 2: sort the BAM work file by queryname, fix mates if the file holds
# paired-end alignments, and, for BWA and BWA-MEM2, filter templates by MAPQ.
function _sort_qname_bam() {
    local func="${1:-}"
    local threads="${2:-}"
    local aligner="${3:-}"
    local mapq="${4:-}"
    local fq_2="${5:-}"
    local fil_wrk="${6:-}"
    local bam_qnam="${7:-}"
    local bam_mate bam_mapq n_pri
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _sort_qname_bam
    [--help] func threads aligner mapq fq_2 fil_wrk bam_qnam

  Sort the BAM work file by queryname into 'bam_qnam', fix mates for paired-end data, and, for BWA and BWA-MEM2 with 'mapq' above 0, keep or drop each template whole by its primary records' MAPQ (Step 2 of 'align_fastqs()').

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  threads : int
    Number of threads to use.

  3  aligner : str
    Alignment program that wrote the BAM work file.

  4  mapq : int
    MAPQ threshold; applied here only for BWA and BWA-MEM2.

  5  fq_2 : file
    Second FASTQ input file, or empty for single-end data; selects mate fixing.

  6  fil_wrk : file
    BAM work file from Step 1.

  7  bam_qnam : file
    Queryname-sorted BAM file to write.

Returns
-------
  Returns 0 when sorting, mate fixing, and filtering succeed; prints an error and returns 1 otherwise.

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
    _sort_qname_bam align_fastqs 2 bwa 1 \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        work/tiny_pe.bam work/tiny_pe.qnam.bam
    '''

  2. Queryname-sort a single-end Bowtie 2 work file; MAPQ was filtered in Step 1.
    '''bash
    _sort_qname_bam align_fastqs 1 bowtie2 1 '' \\
        work/tiny_se.bam work/tiny_se.qnam.bam
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ -f "${fil_wrk}" ]]; then
        if ! \
            samtools sort \
                -@ "${threads}" \
                -n \
                -o "${bam_qnam}" \
                "${fil_wrk}"
        then
            echo_err_func "${func}" \
                "Step #2: failed to queryname-sort ${aligner}-aligned BAM" \
                "file."
            return 1
        fi

        if [[ -n "${fq_2}" ]]; then
            bam_mate="${fil_wrk%.bam}.mate.bam"

            if ! \
                samtools fixmate \
                    -@ "${threads}" \
                    -c \
                    -m \
                    "${bam_qnam}" \
                    "${bam_mate}"
            then
                echo_err_func "${func}" \
                    "Step #2: failed to fix mate pairs in queryname-sorted" \
                    "BAM file."
                return 1
            fi

            if ! \
                mv -f "${bam_mate}" "${bam_qnam}"
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
                    "${func}"
            then
                echo_err_func "${func}" \
                    "Step #2: failed to filter queryname-sorted BAM file by" \
                    "template MAPQ."
                return 1
            fi

            if ! \
                mv -f "${bam_mapq}" "${bam_qnam}"
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
    local threads="${2:-}"
    local qname="${3:-}"
    local fq_2="${4:-}"
    local fil_wrk="${5:-}"
    local bam_qnam="${6:-}"
    local bam_coor
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _sort_coord_bam
    [--help] func threads qname fq_2 fil_wrk bam_qnam

  Sort the queryname-sorted BAM file by coordinates into the BAM work file and index it, deleting the queryname-sorted file unless it is to be retained (Step 3 of 'align_fastqs()').

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  threads : int
    Number of threads to use.

  3  qname : bool
    Whether to retain the queryname-sorted file.

  4  fq_2 : file
    Second FASTQ input file, or empty for single-end data; selects the error message.

  5  fil_wrk : file
    BAM work file to replace with the coordinate-sorted alignments.

  6  bam_qnam : file
    Queryname-sorted BAM file from Step 2.

Returns
-------
  Returns 0 when sorting and indexing succeed; prints an error and returns 1 otherwise.

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
    _sort_coord_bam align_fastqs 2 false '' \\
        work/tiny_se.bam work/tiny_se.qnam.bam
    '''

  2. Coordinate-sort a paired-end work file and retain the queryname-sorted file.
    '''bash
    _sort_coord_bam align_fastqs 2 true \\
        tests/fixtures/align_fastqs/fastq/pe/tiny_pe_R2.atria.fastq.gz \\
        work/tiny_pe.bam work/tiny_pe.qnam.bam
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    # shellcheck disable=SC2086
    if [[ -f "${bam_qnam}" ]]; then
        bam_coor="${fil_wrk%.bam}.coor.bam"

        if ! \
            samtools sort \
                -@ ${threads} \
                -o "${bam_coor}" \
                "${bam_qnam}"
        then
            echo_err_func "${func}" \
                "Step #3: failed to coordinate-sort queryname-sorted BAM file."
            return 1
        fi

        if [[ "${qname}" == "false" ]]; then
            if ! \
                rm -f "${bam_qnam}"
            then
                echo_err_func "${func}" \
                    "Step #3: failed to delete intermediate queryname-sorted" \
                    "BAM file."
                return 1
            fi
        fi

        if ! \
            mv -f "${bam_coor}" "${fil_wrk}"
        then
            echo_err_func "${func}" \
                "Step #3: failed to rename coordinate-sorted BAM file."
            return 1
        fi

        if ! \
            samtools index -@ ${threads} "${fil_wrk}"
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
    local threads="${2:-}"
    local fil_wrk="${3:-}"
    local bam_mrk
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _mark_dup_bam
    [--help] func threads fil_wrk

  Mark duplicate alignments in the coordinate-sorted BAM work file and index it again, skipping files already marked (Step 4 of 'align_fastqs()').

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  threads : int
    Number of threads to use.

  3  fil_wrk : file
    Coordinate-sorted BAM work file from Step 3.

Returns
-------
  Returns 0 when marking and indexing succeed or are skipped; prints an error and returns 1 otherwise.

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
    _mark_dup_bam align_fastqs 2 work/tiny_pe.bam
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

    if [[ -f "${fil_wrk}" ]]; then
        if \
            samtools view -H "${fil_wrk}" | grep -q "@PG.*CL:samtools markdup"
        then
            echo_err_func "${func}" \
                "Step #4: duplicate alignments are already marked in" \
                "coordinate-sorted BAM file. Skipping markdup operations."
        else
            bam_mrk="${fil_wrk%.bam}.mark.bam"

            if ! \
                samtools markdup \
                    -@ "${threads}" \
                    -t \
                    "${fil_wrk}" \
                    "${bam_mrk}"
            then
                echo_err_func "${func}" \
                    "Step #4: failed to mark duplicates in coordinate-sorted" \
                    "BAM file."
                return 1
            fi

            # Replace the original coordinate-sorted BAM with one in which
            # duplicate alignments are marked.
            if ! \
                mv -f "${bam_mrk}" "${fil_wrk}"
            then
                echo_err_func "${func}" \
                    "Step #4: failed to rename duplicate-marked," \
                    "coordinate-sorted BAM file."
                return 1
            fi

            # Index the duplicate-marked coordinate-sorted BAM file.
            if ! \
                samtools index -@ "${threads}" "${fil_wrk}"
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
    local threads="${2:-}"
    local out_fmt="${3:-}"
    local fil_out="${4:-}"
    local fil_wrk="${5:-}"
    local ref_fa="${6:-}"
    local qname="${7:-}"
    local bam_qnam="${8:-}"
    local fil_qnam_out="${9:-}"
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _finalize_align_output
    [--help] func threads out_fmt fil_out fil_wrk ref_fa qname bam_qnam fil_qnam_out

  Convert the BAM work file and any retained queryname-sorted file to CRAM when requested, or rename the retained queryname-sorted BAM file into place (Step 5 of 'align_fastqs()').

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  threads : int
    Number of threads to use for indexing.

  3  out_fmt : str
    Final output format: 'bam' or 'cram'.

  4  fil_out : file
    Path to the final alignment file.

  5  fil_wrk : file
    BAM work file from Step 4.

  6  ref_fa : file
    Reference FASTA path, used for CRAM output.

  7  qname : bool
    Whether the queryname-sorted file is retained.

  8  bam_qnam : file
    Queryname-sorted BAM work file.

  9  fil_qnam_out : file
    Path for the retained queryname-sorted output.

Returns
-------
  Returns 0 when conversion and cleanup succeed; prints an error and returns 1 otherwise.

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
    _finalize_align_output align_fastqs 2 cram work/tiny_pe.cram \\
        work/tiny_pe.bam tests/fixtures/align_fastqs/reference/tiny.fa \\
        true work/tiny_pe.qnam.bam work/tiny_pe.qnam.cram
    '''

  2. Keep BAM output and move the retained queryname-sorted file into place.
    '''bash
    _finalize_align_output align_fastqs 1 bam work/tiny_se.bam \\
        work/tiny_se.bam '' true work/tiny_se.qnam.bam work/tiny_se.qnam.bam
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ "${out_fmt}" == "cram" ]]; then
        if ! \
            samtools view \
                -T "${ref_fa}" \
                -O cram \
                -o "${fil_out}" \
                "${fil_wrk}"
        then
            echo_err_func "${func}" \
                "Step #5: failed to convert BAM work file to final CRAM file."
            return 1
        fi

        if ! \
            samtools index -@ "${threads}" "${fil_out}"
        then
            echo_err_func "${func}" \
                "Step #5: failed to index final CRAM file."
            return 1
        fi

        if [[ "${qname}" == "true" ]]; then
            if ! \
                samtools view \
                    -T "${ref_fa}" \
                    -O cram \
                    -o "${fil_qnam_out}" \
                    "${bam_qnam}"
            then
                echo_err_func "${func}" \
                    "Step #5: failed to convert queryname-sorted BAM" \
                    "intermediate to retained CRAM file."
                return 1
            fi
        fi

        if ! \
            rm -f "${fil_wrk}" "${fil_wrk}.bai"
        then
            echo_err_func "${func}" \
                "Step #5: failed to delete BAM work file and/or BAM index" \
                "after CRAM conversion."
            return 1
        fi

        if [[ "${qname}" == "true" ]]; then
            if ! \
                rm -f "${bam_qnam}"
            then
                echo_err_func "${func}" \
                    "Step #5: failed to delete queryname-sorted BAM work" \
                    "intermediate after CRAM conversion."
                return 1
            fi
        fi
    elif [[ "${out_fmt}" == "bam" && "${qname}" == "true" ]]; then
        if [[ "${bam_qnam}" != "${fil_qnam_out}" ]]; then
            if ! \
                mv -f "${bam_qnam}" "${fil_qnam_out}"
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
    local threads aligner bt2_mode bwa_alg mapq req_flg index ref_fa
    local bt2_mode_set bwa_alg_set
    local fq_1 fq_2 fil_out qname out_fmt
    local fil_wrk fil_qnam_out bam_qnam
    local show_help

    # Assign default argument values.
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
    [--help] [--threads <int>] [--aligner <aligner>] [--bt2_mode <mode>] [--bwa_alg <algorithm>] [--mapq <int>] [--req_flg] --index <path> [--ref_fa <file>] --fq_1 <file> [--fq_2 <file>] --fil_out <file> [--qname]

  Align single- or paired-end Illumina short reads using Bowtie 2, BWA, or BWA-MEM2, and convert the output to a sorted, duplicate-marked, mate-fixed (if working with paired-end reads) alignment file using Samtools.

  Final output format is inferred from '--fil_out': '.bam' or '.cram'. Internal post-processing is performed in BAM format. If CRAM output is requested, the final BAM work file and any retained queryname-sorted BAM intermediate are converted to CRAM in the final step.

  When '--aligner bwa' is selected, '--bwa_alg' may be used to choose either 'mem' or the older backtrack workflow ('aln' with downstream 'samse' or 'sampe').

  Optionally, a queryname-sorted output can be retained with the '--qname' flag; otherwise, the queryname-sorted intermediate/work file is deleted.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

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
  - Because Step #1 currently builds aligner and Samtools argument strings, whitespace is not supported in '--index', '--ref_fa', '--fq_1', '--fq_2', or '--fil_out'.

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
    # TODO: add a '--dry_run' flag that shows what *would* happen


    # Step 0: describe, parse, and validate keyword arguments.
    if [[ -z "${1:-}" || "${1:-}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    # Parse keyword arguments.
    while (( $# > 0 )); do
        case "${1}" in
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
        "${threads}" \
        "${qname}" \
        "${fq_2}" \
        "${fil_wrk}" \
        "${bam_qnam}" \
        || return 1


    # Step 4: mark duplicates and index the BAM work file again.
    _mark_dup_bam \
        "${FUNCNAME[0]}" \
        "${threads}" \
        "${fil_wrk}" \
        || return 1


    # Step 5: convert to CRAM if requested, and keep or delete the
    # queryname-sorted file.
    _finalize_align_output \
        "${FUNCNAME[0]}" \
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
