#!/bin/bash

set -eo pipefail

#  If true, run script in "debug mode"
debug=true

#  If true, run script in "dry-run mode"
dry_run=false

#  Define helper function to validate and log errors for variables
function validate_var() {
    local var_nam="${1}"  # Name of the variable
    local var_val="${2}"  # Value of the variable
    local arr_ref="${3}"  # Reference to the array (as a string)
    local idx="${4}"      # Index of the array element
    local id_tsk="${5}"   # Task ID

    if [[ -z "${var_val}" ]]; then
        echo \
            "Error: Failed to retrieve ${var_nam} for id_tsk=${id_tsk}:" \
            "\${${arr_ref}[${idx}]}." >&2
        return 1
    fi
}


#  Define helper function to determine BAM infile's sequenced-read status
function check_seq_type() {
    local bam="${1}"
    local header

    if [[ -z "${bam}" ]]; then
        echo "Error: BAM infile is required." >&2
        return 1
    fi

    if [[ ! -f "${bam}" ]]; then
        echo "Error: BAM infile does not exist." >&2
        return 1
    fi

    #  Parse BAM header once
    header=$(samtools view -H "${bam}")

    #  Check for aligners Bowtie 2 or BWA
    if echo "${header}" | grep -Eq '^@PG.*ID:(bowtie2|bwa)'; then
        #  Check for presence of read #2
        if echo "${header}" | grep -Eq '^@PG.*ID:(bowtie2|bwa).*(_R2|_2)'; then
            echo "paired"
        else
            echo "single"
        fi
    else
        echo "Error: Unrecognized aligner in BAM file header." >&2
        return 1
    fi
}


#  Define helper function to validate matching values between first
#+ file/numerator and second file/denominator
function validate_match() {
    local var_dsc="${1}"  # Description of variable being compared
    local val_num="${2}"  # Numerator value
    local val_den="${3}"  # Denominator value
    local out_nam="${4}"  # Name of the output variable to assign if matched

    if [[ "${val_num}" == "${val_den}" ]]; then
        #  Assign the matched value to the output variable if specified
        if [[ -n "${out_nam}" ]]; then eval "${out_nam}=\"${val_num}\""; fi
    else
        echo \
            "Error: ${var_dsc} does not match between numerator" \
            "(${val_num}) and denominator (${val_den})." >&2
        return 1
    fi
}


#  Parse keyword arguments, assigning them to variables; most of the argument
#+ inputs are not checked, or not checked thoroughly, as this is performed by
#+ execute_*.sh and again by the deepTools program submitted to SLURM
verbose=false
threads=8
str_fil_num=""
str_fil_den=""
str_fil_stm=""
typ_out="bigwig"
oper="log2"
bin_siz=1
region=""
str_scl_fct=""
norm=""
scl_pre="None"
exact=false
str_usr_frg=""
err_out=""
nam_job="deeptools_compare"
env_nam="env_align"

show_help=$(cat << EOM
$(basename "${0}")
  --verbose <bol> --threads <int> --str_fil_num <str> --str_fil_den <str>
  --str_typ_fil <str> --str_fil_stm <str> --typ_out <str> --bin_siz <int>
  --region <str> --str_scl_fct <str> --norm <str> --exact <bol>
  --str_usr_frg <str> --err_out <str> --nam_job <str> --env_nam <str>

$(basename "${0}") takes the following keyword arguments:
   -v, --verbose      Run in 'verbose mode' (default: ${verbose}).
   -t, --threads      Number of threads to use (default: ${threads}).
  -sn, --str_fil_num  Comma-separated string of BAM, BIGWIG, or BEDGRAPH infiles
                      to be used as the first file/numerator in the comparisons.
  -sd, --str_fil_den  Comma-separated string of BAM, BIGWIG, or BEDGRAPH infiles
                      to be used as the second file/denominator in the
                      comparisons.
  -st, --str_typ_fil  Comma-separated string of numerator-denominator pair file
                      types.
  -ss, --str_fil_stm  Comma-separated string of outstems sans extensions.
  -to, --typ_out      Outfile type: 'bedgraph' or 'bigwig' (default: '${typ_out}').
  -op, --oper         Operation to perform for comparisons of two BAM infiles:
                      'log2', 'ratio', 'subtract', 'add', 'mean',
                      'reciprocal_ratio', 'first', 'second' (default: '${oper}').
                      For more details, see the deepTools documentation.
  -bs, --bin_siz      Bin size in base pairs (default: ${bin_siz}).
   -r, --region       Specify a genomic region to limit the operation. The
                      format can be either: 'chr' (entire chromosome) or
                      'chr:start-stop' (specific range within a chromosome).
  -sf, --str_scl_fct  Comma-separated string of scaling factors for numerator
                      and denominator files. Each scaling factor pair should be
                      in the format 'num:den', where 'num' is the scaling
                      factor for the numerator and 'den' is the scaling factor
                      for the denominator. Multiple pairs should be separated
                      by commas. For example, '0.7:1,0.3:0.6' scales the first
                      numerator file by 0.7 and the first denominator file by
                      1.0, and the second numerator file by 0.3 and the second
                      denominator file by 0.6. This option cannot be used with
                      '--norm <str>'.
  -no, --norm         Use one of the following normalization methods when
                      computing comparisons: 'None', 'RPKM', 'CPM', 'BPM', or
                      'RPGC'. Cannot be used with '--str_scl_fct <str>'.
  -sp, --scl_pre       For comparisons, method to compensate for sequencing
                      depth differences between the samples: 'readCount',
                      'SES', or 'None' (default: '${scl_pre}'). If using
                      '--norm <str>' or '--str_scl_fct <str>', will be set to
                      'None'.
  -ex, --exact        Compute scaling factors based on all alignments
                      (default: ${exact}). Only applicable if '--norm <str>' is
                      specified; otherwise, it is ignored. Significantly slows
                      comparison computation.
  -su, --str_usr_frg  Comma-separated string of fragment lengths for alignment
                      extension.
  -eo, --err_out      Directory to store stderr and stdout outstems.
  -nj, --nam_job      Name of job (default: '${nam_job}').
  -en, --env_nam      Mamba environment to activate (default: '${env_nam}').
EOM
)

if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case "${1}" in
         -v|--verbose)     verbose="${2}";     shift 2 ;;
         -t|--threads)     threads="${2}";     shift 2 ;;
        -sn|--str_fil_num) str_fil_num="${2}"; shift 2 ;;
        -sd|--str_fil_den) str_fil_den="${2}"; shift 2 ;;
        -st|--str_typ_fil) str_typ_fil="${2}"; shift 2 ;;
        -ss|--str_fil_stm) str_fil_stm="${2}"; shift 2 ;;
        -to|--typ_out)     typ_out="${2}";     shift 2 ;;
        -op|--oper)        oper="${2}";        shift 2 ;;
        -bs|--bin_siz)     bin_siz="${2}";     shift 2 ;;
         -r|--region)      region="${2}";      shift 2 ;;
        -sf|--str_scl_fct) str_scl_fct="${2}"; shift 2 ;;
        -no|--norm)        norm="${2}";        shift 2 ;;
        -sp|--scl_pre)     scl_pre="${2}";     shift 2 ;;
        -ex|--exact)       exact="${2}";       shift 2 ;;
        -su|--str_usr_frg) str_usr_frg="${2}"; shift 2 ;;
        -eo|--err_out)     err_out="${2}";     shift 2 ;;
        -nj|--nam_job)     nam_job="${2}";     shift 2 ;;
        -en|--env_nam)     env_nam="${2}";     shift 2 ;;
        *)
            echo "## Unknown parameter passed: ${1} ##" >&2
            echo "" >&2
            echo "${show_help}" >&2
            exit 1
            ;;
    esac
done

#  Validate required arguments
arr_arg_req=(
    "verbose" "threads" "str_fil_num" "str_fil_den" "str_typ_fil" "str_fil_stm"
    "typ_out" "bin_siz" "region" "str_scl_fct" "norm" "scl_pre" "exact" "str_usr_frg"
    "err_out" "nam_job" "env_nam"
)
for var in "${arr_arg_req[@]}"; do
    if [[ -z "${!var}" ]]; then
        echo "Error: '--${var}' is required." >&2
        echo "" >&2
        echo "${show_help}" >&2
        exit 1
    fi
done

#  Validate and standardize output file type
case "${typ_out}" in
    bedgraph|bdg|bg) typ_out="bedgraph" ;;
    bigwig|bw)       typ_out="bigwig"   ;;
    *)
        echo \
            "Error: Unsupported output type: '${typ_out}'. --typ_out must be" \
            "'bigwig', 'bw', 'bedgraph', 'bdg', or 'bg'." >&2
        exit 1
        ;;
esac

#  Debug argument assignments
if ${debug}; then
    echo "verbose=${verbose}"
    echo ""
    echo "threads=${threads}"
    echo ""
    echo "str_fil_num=${str_fil_num}"
    echo ""
    echo "str_fil_den=${str_fil_den}"
    echo ""
    echo "str_typ_fil=${str_typ_fil}"
    echo ""
    echo "str_fil_stm=${str_fil_stm}"
    echo ""
    echo "typ_out=${typ_out}"
    echo ""
    echo "bin_siz=${bin_siz}"
    echo ""
    echo "region=${region}"
    echo ""
    echo "str_scl_fct=${str_scl_fct}"
    echo ""
    echo "norm=${norm}"
    echo ""
    echo "scl_pre=${scl_pre}"
    echo ""
    echo "exact=${exact}"
    echo ""
    echo "str_usr_frg=${str_usr_frg}"
    echo ""
    echo "err_out=${err_out}"
    echo ""
    echo "nam_job=${nam_job}"
    echo ""
    echo "env_nam=${env_nam}"
    echo ""
fi

#  Activate environment
if [[ "${CONDA_DEFAULT_ENV}" != "${env_nam}" ]]; then
    eval "$(conda shell.bash hook)"
    conda activate "${env_nam}" ||
        {
            echo "Error: Failed to activate environment '${env_nam}'." >&2
            exit 1
        }
fi

#  Check that SLURM environment variables are set
if [[ -z "${SLURM_ARRAY_JOB_ID}" ]]; then
    echo "Error: SLURM_ARRAY_JOB_ID is not set." >&2
    exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID}" ]]; then
    echo "Error: SLURM_ARRAY_TASK_ID is not set." >&2
    exit 1
fi

#  Give important SLURM environmental variables shorter names
id_job=${SLURM_ARRAY_JOB_ID}
id_tsk=${SLURM_ARRAY_TASK_ID}

#  Reconstruct arrays from serialized strings
IFS=',' read -r -a arr_fil_num  <<< "${str_fil_num}"
IFS=',' read -r -a arr_fil_den  <<< "${str_fil_den}"
IFS=',' read -r -a arr_typ_fil  <<< "${str_typ_fil}"
IFS=',' read -r -a arr_fil_stm <<< "${str_fil_stm}"
IFS=',' read -r -a arr_scl_fct  <<< "${str_scl_fct}"
IFS=',' read -r -a arr_usr_frg  <<< "${str_usr_frg}"

#  Debug output to check the task ID/array index, number of array elements, and
#+ array element values
if ${debug}; then
    echo "SLURM_ARRAY_TASK_ID=${id_tsk}"
    echo ""
    echo "\${#arr_fil_num[@]}=${#arr_fil_num[@]}"
    echo ""
    echo "arr_fil_num=( ${arr_fil_num[*]} )"
    echo ""
    echo "\${#arr_fil_den[@]}=${#arr_fil_den[@]}"
    echo ""
    echo "arr_fil_den=( ${arr_fil_den[*]} )"
    echo ""
    echo "\${#arr_typ_fil[@]}=${#arr_typ_fil[@]}"
    echo ""
    echo "arr_typ_fil=( ${arr_typ_fil[*]} )"
    echo ""
    echo "\${#arr_fil_stm[@]}=${#arr_fil_stm[@]}"
    echo ""
    echo "arr_fil_stm=( ${arr_fil_stm[*]} )"
    echo ""
    echo "\${#arr_scl_fct[@]}=${#arr_scl_fct[@]}"
    echo ""
    echo "arr_scl_fct=( ${arr_scl_fct[*]} )"
    echo ""
    echo "\${#arr_usr_frg[@]}=${#arr_usr_frg[@]}"
    echo ""
    echo "arr_usr_frg=( ${arr_usr_frg[*]} )"
    echo ""
fi

#  Determine array index based on SLURM_ARRAY_TASK_ID
idx=$(( id_tsk - 1 ))

#  Assign variables from reconstructed arrays based on index
fil_num="${arr_fil_num[idx]}"
fil_den="${arr_fil_den[idx]}"
typ_fil="${arr_typ_fil[idx]}"
fil_stm="${arr_fil_stm[idx]}"
scl_fct="${arr_scl_fct[idx]}"
usr_frg="${arr_usr_frg[idx]}"

#  Debug variable assignments from reconstructed arrays
if ${debug}; then
    echo "fil_num=${fil_num}"
    echo ""
    echo "fil_den=${fil_den}"
    echo ""
    echo "typ_fil=${typ_fil}"
    echo ""
    echo "fil_stm=${fil_stm}"
    echo ""
    echo "scl_fct=${scl_fct}"
    echo ""
    echo "usr_frg=${usr_frg}"
    echo ""
fi

#  Exit if any variable assignment is empty
validate_var "fil_num" "${fil_num}" "arr_fil_num"  "${idx}" "${id_tsk}"
validate_var "fil_den" "${fil_den}" "arr_fil_den"  "${idx}" "${id_tsk}"
validate_var "typ_fil" "${fil_den}" "arr_typ_fil"  "${idx}" "${id_tsk}"
validate_var "fil_stm" "${fil_stm}" "arr_fil_stm" "${idx}" "${id_tsk}"
validate_var "scl_fct" "${scl_fct}" "arr_scl_fct"  "${idx}" "${id_tsk}"
validate_var "usr_frg" "${usr_frg}" "arr_usr_frg"  "${idx}" "${id_tsk}"

#  Determine BAM infiles' sequenced-read status and file type, ensuring
#+ consistency
if [[ "${typ_fil}" == "bam" ]]; then
    validate_match \
        "Sequenced-read statuses" \
        "$(check_seq_type "${fil_num}")" \
        "$(check_seq_type "${fil_den}")" \
        typ_seq
else
    typ_seq="#N/A"
fi

# validate_match \
#     "File extensions" \
#     "${fil_num##*.}" \
#     "${fil_den##*.}" \
#     typ_fil

#  Debug BAM infiles sequenced-read status and file type
# shellcheck disable=SC2154
if ${debug}; then
    echo "typ_seq=${typ_seq}"
    echo ""
    # echo "typ_fil=${typ_fil}"
    # echo ""
fi

#  Derive sample name from fil_stm assignment
samp="${fil_stm##*/}"

#  Debug sample name
if ${debug}; then
    echo "samp=${samp}"
    echo ""
fi

if ${debug}; then
    # shellcheck disable=SC2005
    case "${typ_fil}" in
        bedgraph|bdg|bg|bigwig|bw)
            echo "bigwigCompare \\"
            echo "$(if ${verbose}; then echo "    --verbose"; fi) \\"
            echo "    --numberOfProcessors ${threads} \\"
            echo "    --bigwig1 ${fil_num} \\"
            echo "    --bigwig2 ${fil_den} \\"
            echo "    --outFileName ${fil_stm}.${typ_out} \\"
            echo "    --outFileFormat ${typ_out} \\"
            echo "    --operation ${oper} \\"
            echo "$(
                if [[ "${oper}" =~ ^(log2|ratio)$ ]]; then
                    echo "    --pseudocount 1 1 \\"
                else
                    echo "    \\"
                fi
            )"
            echo "    --binSize ${bin_siz} \\"
            echo "$(
                if [[ "${region}" != "#N/A" ]]; then
                    echo "    --region ${region} \\"
                else
                    echo "    \\"
                fi
            )"
            echo "$(
                if [[ "${scl_fct}" != "#N/A" ]]; then
                    echo "    --scaleFactors ${scl_fct} \\"
                else
                    echo "    \\"
                fi
            )"
            echo "    --skipZeroOverZero \\"
            echo "    --skipNonCoveredRegions"
            echo ""
            echo ""
            ;;
        bam)
            echo "bamCompare \\"
            echo "$(if ${verbose}; then echo "    --verbose"; fi) \\"
            echo "    --numberOfProcessors ${threads} \\"
            echo "    --bamfile1 ${fil_num} \\"
            echo "    --bamfile2 ${fil_den} \\"
            echo "    --outFileName ${fil_stm}.${typ_out} \\"
            echo "    --outFileFormat ${typ_out} \\"
            echo "    --operation ${oper} \\"
            echo "$(
                if [[ "${oper}" =~ ^(log2|ratio)$ ]]; then
                    echo "    --pseudocount 1 1 \\"
                else
                    echo "    \\"
                fi
            )"
            echo "    --binSize ${bin_siz} \\"
            echo "$(
                if [[ "${region}" != "#N/A" ]]; then
                    echo "    --region ${region} \\"
                else
                    echo "    \\"
                fi
            )"
            echo "    --skipZeroOverZero \\"
            echo "    --skipNonCoveredRegions \\"
            echo "$(
                if [[ "${scl_fct}" != "#N/A" ]]; then
                    echo "    --scaleFactors ${scl_fct} \\"
                    echo "    --scaleFactorsMethod None \\"
                elif [[ "${norm}" != "#N/A" ]]; then
                    echo "    --scaleFactorsMethod None \\"
                    echo "    --normalizeUsing ${norm} \\"
                    if [[ "${norm}" == "RPGC" ]]; then
                        echo "    --effectiveGenomeSize ${gen_siz:-11624332} \\"
                    fi
                elif [[ "${scl_pre}" != "#N/A" ]]; then
                    echo "    --scaleFactorsMethod ${scl_pre:-readCount} \\"
                else
                    echo "    \\"
                fi
            )"
            echo "$(
                if ${exact}; then
                    echo "    --exactScaling \\"
                fi
            )"
            echo "$(
                if [[ "${typ_seq}" == "paired" ]]; then
                    echo "    --samFlagInclude 64 \\"
                fi
            )"
            echo "$(
                if [[
                       "${typ_seq}" == "paired" && "${usr_frg}" == "#N/A"
                ]]; then
                    echo "    --extendReads"
                elif [[ "${usr_frg}" != "#N/A" ]]; then
                    echo "    --extendReads ${usr_frg}"
                fi
            )"
            echo ""
            echo ""
            ;;
    esac
fi

#  Execute deepTools bamCompare
if ! ${dry_run}; then
    #  Assign stdout and stderr outfiles to variables
    err_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stderr.txt"
    out_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stdout.txt"
    err_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stderr.txt"
    out_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stdout.txt"

    #  Give the initial stderr and stdout TXT outfiles more descriptive names
    ln -f "${err_ini}" "${err_dsc}"
    ln -f "${out_ini}" "${out_dsc}"

    # shellcheck disable=SC2046
    case "${typ_fil}" in
        bedgraph|bdg|bg|bigwig|bw)
            bigwigCompare \
                $(if ${verbose}; then echo "--verbose"; fi) \
                --numberOfProcessors "${threads}" \
                --bigwig1 "${fil_num}" \
                --bigwig2 "${fil_den}" \
                --outFileName "${fil_stm}.${typ_out}" \
                --outFileFormat "${typ_out}" \
                --operation "${oper}" \
                $(
                    if [[ "${oper}" =~ ^(log2|ratio)$ ]]; then
                        echo "--pseudocount 1 1"
                    fi
                ) \
                --binSize "${bin_siz}" \
                $(
                    if [[ "${region}" != "#N/A" ]]; then
                        echo "--region ${region}"
                    fi
                ) \
                $(
                    if [[ "${scl_fct}" != "#N/A" ]]; then
                        echo "--scaleFactors ${scl_fct}"
                    fi
                ) \
                --skipZeroOverZero \
                --skipNonCoveredRegions
            ;;
        bam)
            bamCompare \
                $(if ${verbose}; then echo "--verbose"; fi) \
                --numberOfProcessors "${threads}" \
                --bamfile1 "${fil_num}" \
                --bamfile2 "${fil_den}" \
                --outFileName "${fil_stm}.${typ_out}" \
                --outFileFormat "${typ_out}" \
                --operation "${oper}" \
                $(
                    if [[ "${oper}" =~ ^(log2|ratio)$ ]]; then
                        echo "--pseudocount 1 1"
                    fi
                ) \
                --binSize "${bin_siz}" \
                $(
                    if [[ "${region}" != "#N/A" ]]; then
                        echo "--region ${region}"
                    fi
                ) \
                --skipZeroOverZero \
                --skipNonCoveredRegions \
                $(
                    if [[ "${scl_fct}" != "#N/A" ]]; then
                        echo "--scaleFactors ${scl_fct}"
                        echo "--scaleFactorsMethod None"
                    elif [[ "${norm}" != "#N/A" ]]; then
                        echo "--scaleFactorsMethod None"
                        echo "--normalizeUsing ${norm}"
                        if [[ "${norm}" == "RPGC" ]]; then
                            echo "--effectiveGenomeSize ${gen_siz:-11624332}"
                        fi
                    else
                        echo "--scaleFactorsMethod ${scl_pre:-readCount}"
                    fi
                ) \
                $(if ${exact}; then echo "--exactScaling"; fi) \
                $(
                    if [[ "${typ_seq}" == "paired" ]]; then
                        echo "--samFlagInclude 64"
                    fi
                ) \
                $(
                    if [[
                           "${typ_seq}" == "paired" && "${usr_frg}" == "#N/A"
                    ]]; then
                        echo "--extendReads"
                    elif [[ "${usr_frg}" != "#N/A" ]]; then
                        echo "--extendReads ${usr_frg}"
                    fi
                )
            ;;
    esac

    #  Remove the initial SLURM stderr and stdout TXT outfiles
    rm "${err_ini}" "${out_ini}"
fi
