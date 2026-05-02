#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: align_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# _check_path_whitespace
# align_fastqs


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
    _dir_src_aln="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    source "${_dir_src_aln}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_aln}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_aln}" \
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


#MAYBE: move to 'check_inputs.sh' and make "public"?
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


function align_fastqs() {
    local threads aligner bt2_aln bwa_alg mapq req_flg index ref_fa
    local fq_1 fq_2 outfile qname out_fmt
    local args_sam args_bt2 args_bwa args_bwa_aln
    local fil_wrk fil_qnam_out bam_qnam bam_mate bam_coor bam_mrk
    local fil_sai_1 fil_sai_2
    local show_help

    #  Assign default argument values
    threads=1
    aligner="bowtie2"
    bt2_aln="end-to-end"
    bwa_alg="mem"
    ref_fa=""
    mapq=0
    req_flg=false
    qname=false

    show_help=$(cat << EOM
Description:
  Align single- or paired-end Illumina short reads using Bowtie 2, BWA, or BWA-MEM2, and convert the output to a sorted, duplicate-marked, mate-fixed (if working with paired-end reads) alignment file using Samtools.

  Final output format is inferred from '--outfile': '.bam' or '.cram'. Internal post-processing is performed in BAM format. If CRAM output is requested, the final BAM work file and any retained queryname-sorted BAM intermediate are converted to CRAM in the final step.

  When '--aligner bwa' is selected, '--bwa_alg' may be used to choose either 'mem' or the older backtrack workflow ('aln' with downstream 'samse' or 'sampe').

  Optionally, a queryname-sorted output can be retained with the '--qname' flag; otherwise, the corresponding intermediate/work file is deleted.

Keyword arguments:
   -t, --threads  <int>  Number of threads to use (default: ${threads}).
   -a, --aligner  <str>  Alignment program to use: 'bowtie2', 'bt2' (alias of 'bowtie2'), 'bwa', or 'bwa-mem2' (default: '${aligner}').
  -2a, --bt2_aln  <str>  Bowtie 2 alignment type when '--aligner bowtie2': 'local', 'global' (alias of 'end-to-end'), or 'end-to-end' (default: '${bt2_aln}'; ignored if not '--aligner bowtie2').
  -ba, --bwa_alg  <str>  BWA algorithm to use when '--aligner bwa' is selected: 'mem' or 'aln' (default: '${bwa_alg}'; ignored if not '--aligner bwa').
  -mq, --mapq     <int>  MAPQ threshold for filtering alignments during BAM work-file generation (default: ${mapq}).
  -rf, --req_flg  <flg>  Require flag bit 2, signifying that paired-end alignments are properly paired; used for filtering alignments during BAM work-file generation (ignored if working with single-end sequenced reads).
  -ix, --index    <str>  Path to the aligner index/reference.
   -r, --ref      <str>  Reference FASTA path required when '--outfile' ends in '.cram'.
  -f1, --fq_1     <str>  Path to a FASTQ file; if working with paired-end sequenced reads, path to the first of two FASTQ files.
  -f2, --fq_2     <str>  If working with paired-end sequenced reads, path to the second of two FASTQ files.
   -o, --outfile  <str>  Path to the final alignment outfile (must end in '.bam' or '.cram').
  -qn, --qname    <flg>  Retain a queryname-sorted output using the same final extension as '--outfile'.

Returns:
  Creates a final '.bam' or '.cram' outfile at the specified path. A corresponding index ('.bai' or '.crai') is also created.

Dependencies:
  - Bash >= 4.4
  - Bowtie 2 == 2.5.4, BWA == 0.7.19, or BWA-MEM2 == 2.2.1
  - Samtools >= 1.21

Notes:
  - For '-ix' / '--index' when using Bowtie 2, the path should end with the index stem: "path/to/dir/stem"; when using 'bwa' or 'bwa-mem2', the path should be the indexed reference FASTA path (for example, ".fa").
  - '--bwa_alg' is used only when '--aligner bwa'. It is ignored when '--aligner bowtie2' and must remain 'mem' when '--aligner bwa-mem2'.
  - For '-qn' / '--qname', the retained queryname-sorted output will have the same path and stem assigned to '--outfile', except '.qnam' will be inserted before the final extension (for example, '.qnam.bam' or '.qnam.cram'). For CRAM output, the retained queryname-sorted BAM work file is converted in Step #5.
  - '--ref' is required when '--outfile' ends in '.cram', since CRAM writing requires a reference FASTA.
  - Because Step #1 currently builds aligner and Samtools argument strings, whitespace is not supported in '--index', '--ref', '--fq_1', '--fq_2', or '--outfile'.

Examples:
  '''bash
  align_fastqs
      --threads 4
      --aligner "bowtie2"
      --bt2_aln "global"
      --req_flg
      --index "\${HOME}/path/index"
      --ref "\${HOME}/path/reference.fa"
      --fq_1 "\${HOME}/path/infile_R1.fastq.gz"
      --fq_2 "\${HOME}/path/infile_R2.fastq.gz"
      --outfile "\${HOME}/path/outfile.cram"

  align_fastqs
      --threads \${threads}
      --aligner \${aligner}
      \$(if [[ \${aligner} =~ ^(bowtie2|bt2)$ ]]; then echo "--bt2_aln \${bt2_aln}"; fi)
      \$(if [[ \${aligner} == "bwa" ]]; then echo "--bwa_alg \${bwa_alg}"; fi)
      \$(if [[ \${mapq} -gt 0 ]]; then echo "--mapq \${mapq}"; fi)
      \$(if \${req_flg}; then echo "--req_flg"; fi)
      --index \${index}
      \$(if [[ \${outfile} == *.cram ]]; then echo "--ref \${ref_fa}"; fi)
      --fq_1 \${fq_1}
      \$(if [[ -n \${fq_2} ]]; then echo "--fq_2 \${fq_2}"; fi)
      --outfile \${outfile}
      \$(if \${qname}; then echo "--qname"; fi)
           > \${err_out}/\${nam_job}.\${samp}.stdout.txt
          2> \${err_out}/\${nam_job}.\${samp}.stderr.txt
  '''

#TODO:
  - Consider adding support for additional BWA-family tuning options.
  - Consider whether the BWA-backtrack branch warrants a small helper-layer split to reduce branching in Step 1.
EOM
    )


    #  Step 0 -----------------------------------------------------------------
    #  Describe keyword arguments
    if [[ -z "${1:-}" || "${1:-}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    #  Parse keyword arguments
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
                aligner="${2:-,,}"
                shift 2
                ;;
            -2a|-bn|--bt2[_-]aln)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                bt2_aln="${2:-,,}"
                shift 2
                ;;
            -ba|--bwa[_-]alg)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                bwa_alg="${2:-,,}"
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
            -rf|--req[_-]flg)
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
            -r|--ref|--reference|--fa[_-]ref|--fasta[_-]ref)
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
             -o|-fo|--outfile|--fil[_-]out)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                outfile="${2:-}"
                shift 2
                ;;
            -qn|--qnam|--qname)
                qname=true
                shift 1
                ;;
            *)
                echo "## Unknown argument passed: '${1}'. ##" >&2
                echo >&2
                echo "${show_help}" >&2
                return 1
                ;;
        esac
    done

    #  Validate keyword arguments
    if [[ -z "${threads}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'--threads' is required."
        return 1
    fi

    if [[ ! "${threads}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'--threads' was assigned '${threads}' but must be a positive" \
            "integer greater than or equal to 1."
        return 1
    fi

    case "${aligner}" in
        bt2) aligner="bowtie2" ;;
    esac

    case "${aligner}" in
        bowtie2)
            case "${bt2_aln}" in
                local|global|end-to-end) : ;;
                na)
                    echo_err_func "${FUNCNAME[0]}" \
                        "'--bt2_aln' cannot be 'NA' when using 'bowtie2'."
                    return 1
                    ;;
                *)
                    echo_err_func "${FUNCNAME[0]}" \
                        "selection associated with '--bt2_aln' is not valid:" \
                        "'${bt2_aln}'; must be 'local', 'global', or" \
                        "'end-to-end'."
                    return 1
                    ;;
            esac
            ;;

        bwa)
            case "${bwa_alg}" in
                mem|aln) : ;;
                *)
                    echo_err_func "${FUNCNAME[0]}" \
                        "selection associated with '--bwa_alg' is not" \
                        "valid: '${bwa_alg}'; must be 'mem' or 'aln'."
                    return 1
                    ;;
            esac
            ;;

        bwa-mem2)
            if [[ "${bwa_alg}" != "mem" ]]; then
                echo_err_func "${FUNCNAME[0]}" \
                    "'--bwa_alg' must be 'mem' when using 'bwa-mem2', but" \
                    "was assigned '${bwa_alg}'."
                return 1
            fi
            ;;

        *)
            echo_err_func "${FUNCNAME[0]}" \
                "selection associated with '--aligner' is not valid:" \
                "'${aligner}'; must be 'bowtie2', 'bwa', or 'bwa-mem2'."
            return 1
            ;;
    esac

    if [[ ! "${mapq}" =~ ^[0-9]+$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "--mapq was assigned '${mapq}' but must be an integer greater" \
            "than or equal to 0."
        return 1
    fi

    if [[ -z "${index}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'--index' is a required argument."
        return 1
    fi

    if [[ ! -d "$(dirname "${index}")" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "directory associated with '--index' does not exist:" \
            "'$(dirname "${index}")'."
        return 1
    fi

    _check_path_whitespace \
        "${index}" "--index" "${FUNCNAME[0]}" \
        || return 1

    if [[ -z "${outfile}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'--outfile' is a required argument."
        return 1
    fi

    if [[ "${aligner}" =~ ^(bwa|bwa-mem2)$ && ! -f "${index}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "file associated with '--index' does not exist: '${index}'."
        return 1
    fi

    _check_path_whitespace \
        "${outfile}" "--outfile" "${FUNCNAME[0]}" \
        || return 1

    if [[ ! -d "$(dirname "${outfile}")" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "directory associated with '--outfile' does not exist:" \
            "'$(dirname "${outfile}")'"
        return 1
    fi

    case "${outfile}" in
        *.bam)
            out_fmt="bam"
            fil_wrk="${outfile}"
            ;;

        *.sam)
            echo_err_func "${FUNCNAME[0]}" \
                "SAM output is not supported by 'align_fastqs()'. Please use" \
                "'.bam' or '.cram' for '--outfile'."
            return 1
            ;;

        *.cram)
            out_fmt="cram"
            fil_wrk="${outfile%.cram}.bam"
            ;;

        *)
            echo_err_func "${FUNCNAME[0]}" \
                "'--outfile' must end in '.bam' or '.cram', but was" \
                "assigned '${outfile}'."
            return 1
            ;;
    esac

    case "${out_fmt}" in
        bam)  fil_qnam_out="${outfile%.bam}.qnam.bam"   ;;
        cram) fil_qnam_out="${outfile%.cram}.qnam.cram" ;;
    esac

    if [[ "${out_fmt}" == "cram" ]]; then
        if [[ -z "${ref_fa}" ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "'--ref' is required when '--outfile' ends in '.cram'."
            return 1
        fi

        _check_path_whitespace \
            "${ref_fa}" "--ref" "${FUNCNAME[0]}" \
            || return 1

        if [[ ! -f "${ref_fa}" ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "file associated with '--ref' does not exist: '${ref_fa}'."
            return 1
        fi
    elif [[ -n "${ref_fa}" ]]; then
        _check_path_whitespace \
            "${ref_fa}" "--ref" "${FUNCNAME[0]}" \
            || return 1

        echo_warn_func "${FUNCNAME[0]}" \
            "'--ref' was supplied but will be ignored because '--outfile'" \
            "does not end in '.cram'."
    fi

    if [[ -z "${fq_1}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'--fq_1' is a required argument."
        return 1
    fi

    if [[ ! -f "${fq_1}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "file associated with '--fq_1' does not exist: '${fq_1}'."
        return 1
    fi

    _check_path_whitespace \
        "${fq_1}" "--fq_1" "${FUNCNAME[0]}" \
        || return 1

    if [[ -n "${fq_2}" ]]; then
        if [[ ! -f "${fq_2}" ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "file associated with '--fq_2' does not exist: '${fq_2}'."
            return 1
        fi

        _check_path_whitespace \
            "${fq_2}" "--fq_2" "${FUNCNAME[0]}" \
            || return 1
    fi

    #  Warn when supplied arguments will be ignored based on other selections
    if [[ "${aligner}" != "bowtie2" && "${bt2_aln}" != "end-to-end" ]]; then
        echo_warn_func "${FUNCNAME[0]}" \
            "'--bt2_aln' was supplied but will be ignored because '--aligner'" \
            "is not 'bowtie2'."
    fi

    if [[ "${aligner}" != "bwa" && "${bwa_alg}" != "mem" ]]; then
        echo_warn_func "${FUNCNAME[0]}" \
            "'--bwa_alg' was supplied but will be ignored because '--aligner'" \
            "is not 'bwa'."
    fi

    if [[ -z "${fq_2}" && "${req_flg}" == "true" ]]; then
        echo_warn_func "${FUNCNAME[0]}" \
            "'--req_flg' was supplied but will be ignored because only one" \
            "FASTQ file was provided."
    fi


    #  Step 1 -----------------------------------------------------------------
    #  Based on parsed arguments, construct the call to Samtools
    args_sam="-@ ${threads}"
    args_sam+="$(
        if [[ "${req_flg}" == "true" ]] && [[ -n "${fq_2}" ]]; then
            echo " -f 2"
        fi
    )"
    args_sam+=" -q ${mapq} -o ${fil_wrk}"

    #  Based on parsed arguments, construct and run the call to Bowtie 2,
    #+ BWA, or BWA-MEM2 with output piped to the above-constructed Samtools
    #+ call where appropriate
    # shellcheck disable=SC2086
    if [[ "${aligner}" == "bowtie2" ]]; then
        #  Assign Bowtie 2 arguments
        args_bt2="-p ${threads} -x ${index}"

        if [[ "${bt2_aln}" == "local" ]]; then
            args_bt2+=" --very-sensitive-local"
        elif [[ "${bt2_aln}" =~ ^(end-to-end|global)$ ]]; then
            args_bt2+=" --very-sensitive"
        fi

        #  Exclude unaligned reads from output
        args_bt2+=" --no-unal --phred33"

        #  Run Bowtie 2 with paired-end (top) or single-end (bottom) sequenced
        #+ reads
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
                echo_err_func "${FUNCNAME[0]}" \
                    "Step #1: failed to align paired-end data with" \
                    "'${aligner}' ('${bt2_aln}'), and/or failed to process" \
                    "the alignments with 'samtools view'."
                return 1
            else
                echo_err_func "${FUNCNAME[0]}" \
                    "Step #1: failed to align single-end data with" \
                    "'${aligner}' ('${bt2_aln}'), and/or failed to process" \
                    "the alignments with 'samtools view'."
                return 1
            fi
        fi
    elif [[ "${aligner}" == "bwa" ]]; then
        #  Assign BWA arguments
        args_bwa="-t ${threads} ${index}"

        if [[ "${bwa_alg}" == "mem" ]]; then
            #  Run BWA-MEM with paired-end (top) or single-end (bottom)
            #+ sequenced reads, excluding unaligned reads from output
            # shellcheck disable=SC2086,SC2046
            if [[ -n "${fq_2}" ]]; then
                if ! (
                    set -o pipefail
                    bwa mem ${args_bwa} "${fq_1}" "${fq_2}" \
                        | samtools view ${args_sam}
                ); then
                    echo_err_func "${FUNCNAME[0]}" \
                        "Step #1: failed to align paired-end data with" \
                        "'bwa mem', and/or failed to process the alignments" \
                        "with 'samtools view'."
                    return 1
                fi
            else
                if ! (
                    set -o pipefail
                    bwa mem ${args_bwa} "${fq_1}" | samtools view ${args_sam}
                ); then
                    echo_err_func "${FUNCNAME[0]}" \
                        "Step #1: failed to align single-end data with" \
                        "'bwa mem', and/or failed to process the alignments" \
                        "with 'samtools view'."
                    return 1
                fi
            fi
        elif [[ "${bwa_alg}" == "aln" ]]; then
            #  Assign BWA-backtrack arguments
            args_bwa_aln="-t ${threads}"

            fil_sai_1="${fil_wrk%.bam}.R1.sai"
            fil_sai_2=""

            if ! \
                bwa aln ${args_bwa_aln} "${index}" "${fq_1}" > "${fil_sai_1}"
            then
                echo_err_func "${FUNCNAME[0]}" \
                    "Step #1: failed to generate first '.sai' file with" \
                    "'bwa aln'."
                return 1
            fi

            if [[ -n "${fq_2}" ]]; then
                fil_sai_2="${fil_wrk%.bam}.R2.sai"

                if ! \
                    bwa aln ${args_bwa_aln} "${index}" "${fq_2}" \
                        > "${fil_sai_2}"
                then
                    echo_err_func "${FUNCNAME[0]}" \
                        "Step #1: failed to generate second '.sai' file" \
                        "with 'bwa aln'."
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
                    echo_err_func "${FUNCNAME[0]}" \
                        "Step #1: failed to align paired-end data with" \
                        "'bwa aln'/'bwa sampe', and/or failed to process the" \
                        "alignments with 'samtools view'."
                    return 1
                fi

                if ! \
                    rm -f "${fil_sai_1}" "${fil_sai_2}"
                then
                    echo_err_func "${FUNCNAME[0]}" \
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
                    echo_err_func "${FUNCNAME[0]}" \
                        "Step #1: failed to align single-end data with" \
                        "'bwa aln'/'bwa samse', and/or failed to process the" \
                        "alignments with 'samtools view'."
                    return 1
                fi

                if ! \
                    rm -f "${fil_sai_1}"
                then
                    echo_err_func "${FUNCNAME[0]}" \
                        "Step #1: failed to delete temporary '.sai' file."
                    return 1
                fi
            fi
        fi
    elif [[ "${aligner}" == "bwa-mem2" ]]; then
        #  Assign BWA-MEM2 arguments
        args_bwa="-t ${threads} ${index}"

        #  Run BWA-MEM2 with paired-end (top) or single-end (bottom)
        #+ sequenced reads, excluding unaligned reads from output
        # shellcheck disable=SC2086,SC2046
        if [[ -n "${fq_2}" ]]; then
            if ! (
                set -o pipefail
                bwa-mem2 mem ${args_bwa} "${fq_1}" "${fq_2}" \
                    | samtools view ${args_sam}
            ); then
                echo_err_func "${FUNCNAME[0]}" \
                    "Step #1: failed to align paired-end data with" \
                    "'bwa-mem2 mem', and/or failed to process the" \
                    "alignments with 'samtools view'."
                return 1
            fi
        else
            if ! (
                set -o pipefail
                bwa-mem2 mem ${args_bwa} "${fq_1}" | samtools view ${args_sam}
            ); then
                echo_err_func "${FUNCNAME[0]}" \
                    "Step #1: failed to align single-end data with" \
                    "'bwa-mem2 mem', and/or failed to process the" \
                    "alignments with 'samtools view'."
                return 1
            fi
        fi
    fi


    #  Step 2 -----------------------------------------------------------------
    #  Sort the BAM file by queryname, and fix mates if the BAM file is made up
    #+ of paired-end alignments
    if [[ -f "${fil_wrk}" ]]; then
        bam_qnam="${fil_wrk%.bam}.qnam.bam"

        if ! \
            samtools sort \
                -@ "${threads}" \
                -n \
                -o "${bam_qnam}" \
                "${fil_wrk}"
        then
            echo_err_func "${FUNCNAME[0]}" \
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
                echo_err_func "${FUNCNAME[0]}" \
                    "Step #2: failed to fix mate pairs in queryname-sorted" \
                    "BAM file."
                return 1
            fi

            if ! \
                mv -f "${bam_mate}" "${bam_qnam}"
            then
                echo_err_func "${FUNCNAME[0]}" \
                    "Step #2: failed to rename mate pair-fixed," \
                    "queryname-sorted BAM file."
                return 1
            fi
        fi

    else
        echo_err_func "${FUNCNAME[0]}" \
            "Step #2: BAM work file from ${aligner} alignment does not exist."
        return 1
    fi


    #  Step 3 -----------------------------------------------------------------
    #  Sort the queryname-sorted, fixmate-adjusted BAM file by coordinates,
    #+ then index the coordinate-sorted BAM file
    # shellcheck disable=SC2086
    if [[ -f "${bam_qnam}" ]]; then
        bam_coor="${fil_wrk%.bam}.coor.bam"

        if ! \
            samtools sort \
                -@ ${threads} \
                -o "${bam_coor}" \
                "${bam_qnam}"
        then
            echo_err_func "${FUNCNAME[0]}" \
                "Step #3: failed to coordinate-sort queryname-sorted BAM file."
            return 1
        fi

        if [[ "${qname}" == "false" ]]; then
            if ! \
                rm -f "${bam_qnam}"
            then
                echo_err_func "${FUNCNAME[0]}" \
                    "Step #3: failed to delete intermediate queryname-sorted" \
                    "BAM file."
                return 1
            fi
        fi

        if ! \
            mv -f "${bam_coor}" "${fil_wrk}"
        then
            echo_err_func "${FUNCNAME[0]}" \
                "Step #3: failed to rename coordinate-sorted BAM file."
            return 1
        fi

        if ! \
            samtools index -@ ${threads} "${fil_wrk}"
        then
            echo_err_func "${FUNCNAME[0]}" \
                "Step #3: failed to index renamed, coordinate-sorted BAM file."
            return 1
        fi
    else
        if [[ -n "${fq_2}" ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "Step #3: mate-fixed, queryname-sorted BAM file does not" \
                "exist."
        else
            echo_err_func "${FUNCNAME[0]}" \
                "Step #3: queryname-sorted BAM file does not exist."
        fi
        return 1
    fi


    #  Step #4 ----------------------------------------------------------------
    #  Mark duplicate alignments in the coordinate-sorted BAM file; then,
    #+ index the file again
    if [[ -f "${fil_wrk}" ]]; then
        if \
            samtools view -H "${fil_wrk}" | grep -q "@PG.*CL:samtools markdup"
        then
            echo_err_func "${FUNCNAME[0]}" \
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
                echo_err_func "${FUNCNAME[0]}" \
                    "Step #4: failed to mark duplicates in coordinate-sorted" \
                    "BAM file."
                return 1
            fi

            #  Replace the original coordinate-sorted BAM with one in which
            #+ duplicates alignments are marked
            if ! \
                mv -f "${bam_mrk}" "${fil_wrk}"
            then
                echo_err_func "${FUNCNAME[0]}" \
                    "Step #4: failed to rename duplicate-marked," \
                    "coordinate-sorted BAM file."
                return 1
            fi

            #  Index the duplicate-marked coordinate-sorted BAM file
            if ! \
                samtools index -@ "${threads}" "${fil_wrk}"
            then
                echo_err_func "${FUNCNAME[0]}" \
                    "Step #4: failed to index duplicate-marked," \
                    "coordinate-sorted BAM file."
                return 1
            fi
        fi
    else
        echo_err_func "${FUNCNAME[0]}" \
            "Step #4: coordinate-sorted BAM file does not exist."
        return 1
    fi


    #  Step 5 -----------------------------------------------------------------
    #  Convert the final BAM work file to CRAM if requested; otherwise retain
    #+ BAM output as-is
    if [[ "${out_fmt}" == "cram" ]]; then
        if ! \
            samtools view \
                -T "${ref_fa}" \
                -O cram \
                -o "${outfile}" \
                "${fil_wrk}"
        then
            echo_err_func "${FUNCNAME[0]}" \
                "Step #5: failed to convert BAM work file to final CRAM file."
            return 1
        fi

        if ! \
            samtools index -@ "${threads}" "${outfile}"
        then
            echo_err_func "${FUNCNAME[0]}" \
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
                echo_err_func "${FUNCNAME[0]}" \
                    "Step #5: failed to convert queryname-sorted BAM" \
                    "intermediate to retained CRAM file."
                return 1
            fi
        fi

        if ! \
            rm -f "${fil_wrk}" "${fil_wrk}.bai"
        then
            echo_err_func "${FUNCNAME[0]}" \
                "Step #5: failed to delete BAM work file and/or BAM index" \
                "after CRAM conversion."
            return 1
        fi

        if [[ "${qname}" == "true" ]]; then
            if ! \
                rm -f "${bam_qnam}"
            then
                echo_err_func "${FUNCNAME[0]}" \
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
                echo_err_func "${FUNCNAME[0]}" \
                    "Step #5: failed to rename retained queryname-sorted BAM" \
                    "intermediate."
                return 1
            fi
        fi
    fi
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
