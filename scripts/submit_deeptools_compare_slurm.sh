#!/bin/bash

#  If true, run script in debug mode
debug=true

#  Define helper function to determine BAM infile's sequenced-read status
function check_typ_seq() {
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


#  Parse keyword arguments, assigning them to variables; most of the argument
#+ inputs are not checked, or not checked thoroughly, as this is performed by
#+ execute_*.sh and again by the function submitted to SLURM
verbose=false
threads=8
str_fil_num=""
str_fil_den=""
str_outstem=""
typ_out="bigwig"
operation="log2"
bin_siz=1
region=""
str_scl_fct=""
norm=""
method="None"
exact=false
str_usr_frg=""
err_out=""
nam_job="deeptools_compare"
env_nam="env_align"

show_help=$(cat << EOM
$(basename "${0}")
  [--verbose] --threads <int> --str_fil_num <str> --str_fil_den <str>
  --str_outstem <str> --typ_out <str> --bin_siz <int> [--region <str>]
  --str_scl_fct <str> --norm <str> [--exact] [--str_usr_frg <str>]
  --err_out <str> --nam_job <str> --env_nam <str>

$(basename "${0}") takes the following keyword arguments:
   -v, --verbose      Run in 'verbose mode' (default: ${verbose}).
   -t, --threads      Number of threads to use (default: ${threads}).
  -sn, --str_fil_num  Comma-separated list of BAM, BIGWIG, or BEDGRAPH infiles
                      to be used as the numerator in the comparisons.
  -sd, --str_fil_den  Comma-separated list of BAM, BIGWIG, or BEDGRAPH infiles
                      to be used as the denominator in the comparisons.
  -so, --str_outstem  Comma-separated string of outfile stems.
  -to, --typ_out      Outfile type: 'bedgraph' or 'bigwig' (default: ${typ_out}).
  -op, --operation    Operation to perform for comparisons of two BAM infiles
                      (default: ${operation}).
  -bs, --bin_siz      Bin size in base pairs (default: ${bin_siz}).
   -r, --region       Region in 'chr' or 'chr:start-stop' format.
  -sf, --str_scl_fct  Comma-separated string of scaling factors. Cannot be
                      used with --norm.
  -no, --norm         Use one of the following normalization methods when
                      computing coverage: 'None', 'RPKM', 'CPM', 'BPM', or
                      'RPGC'. Cannot be used with --str_scl_fct.
   -m, --method       For comparisons, method to compensate for sequencing
                      depth differences between the samples (default: ${norm}).
                      If using --norm, will be set to 'None'.
  -ex, --exact        Compute scaling factors based on all alignments. Only
                      applicable if '--norm <str>' is specified. Significantly
                      slows coverage computation (default: ${exact}).
  -su, --str_usr_frg  Comma-separated string of fragment lengths for alignment
                      extension.
  -eo, --err_out      Directory to store stderr and stdout outfiles.
  -nj, --nam_job      Name of job (default ${nam_job}).
  -en, --env_nam      Mamba environment to activate (default: ${env_nam}).
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
        -so|--str_outstem) str_outstem="${2}"; shift 2 ;;
        -to|--typ_out)     typ_out="${2}";     shift 2 ;;
        -bs|--bin_siz)     bin_siz="${2}";     shift 2 ;;
         -r|--region)      region="${2}";      shift 2 ;;
        -sf|--str_scl_fct) str_scl_fct="${2}"; shift 2 ;;
        -no|--norm)        norm="${2}";        shift 2 ;;
         -m|--method)      method="${2}";      shift 2 ;;
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
    "verbose" "threads" "str_fil_num" "str_fil_den" "str_outstem" "typ_out"
    "bin_siz" "region" "str_scl_fct" "norm" "method" "exact" "str_usr_frg"
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

#  Validate specification of '--norm', '--str_scl_fct', and '--method'
if [[ -z "${norm}" || "${norm}" == "#N/A" ]]; then
    #  Ensure '--str_scl_fct' is specified and valid when '--norm' is not used
    if [[ -z "${str_scl_fct}" || "${str_scl_fct}" =~ (^|,)?'#N/A'(,|$) ]]; then
        echo \
            "Error: If '--norm <str>' is not specified, '--str_scl_fct' must" \
            "be provided and cannot contain '#N/A' elements." >&2
        exit 1
    fi

    #  If '--method' is provided, ensure it is compatible with missing '--norm'
    if [[ "${method}" != "None" ]]; then
        echo \
            "Warning: '--method' is ignored because '--norm <str>' is not" \
            "specified." >&2
        method="None"
    fi
else
    #  Ensure '--str_scl_fct' is not specified or is invalid when '--norm' is
    #+ used
    if [[
        -n "${str_scl_fct}" && ! "${str_scl_fct}" =~ (^|,)?'#N/A'(,|$)
    ]]; then
        echo \
            "Error: If '--norm <str>' is specified, '--str_scl_fct' must not" \
            "be provided or must only contain '#N/A' elements." >&2
        exit 1
    fi

    #  Automatically set '--method None' if '--norm' is specified
    if [[ "${method}" != "None" ]]; then
        echo \
            "Warning: '--method' is overridden to 'None' because '--norm' is" \
            "specified." >&2
        method="None"
    fi
fi

#  Validate and standardize output file type
case "${typ_out}" in
    bigwig|bw) typ_out="bigwig" ;;
    bedgraph|bg) typ_out="bedgraph" ;;
    *)
        echo \
            "Error: Unsupported output type: '${typ_out}'. Supported types:" \
            "'bigwig', 'bw', 'bedgraph', or 'bg'." >&2 
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
    echo "str_outstem=${str_outstem}"
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
    echo "method=${method}"
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
IFS=',' read -r -a arr_outstems <<< "${str_outstem}"
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
    echo "\${#arr_outstems[@]}=${#arr_outstems[@]}"
    echo ""
    echo "arr_outstems=( ${arr_outstems[*]} )"
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
outstem="${arr_outstems[idx]}"
scl_fct="${arr_scl_fct[idx]}"
usr_frg="${arr_usr_frg[idx]}"

#  Debug variable assignments from reconstructed arrays
if ${debug}; then
    echo "fil_num=${fil_num}"
    echo ""
    echo "fil_den=${fil_den}"
    echo ""
    echo "outstem=${outstem}"
    echo ""
    echo "scl_fct=${scl_fct}"
    echo ""
    echo "usr_frg=${usr_frg}"
    echo ""
fi

#  Exit if any variable assignment is empty
if [[ -z "${fil_num}" ]]; then
    echo \
        "Error: Failed to retrieve fil_num for id_tsk=${id_tsk}:" \
        "\${arr_fil_num[${idx}]}." >&2
    exit 1
fi

if [[ -z "${fil_den}" ]]; then
    echo \
        "Error: Failed to retrieve fil_den for id_tsk=${id_tsk}:" \
        "\${arr_fil_den[${idx}]}." >&2
    exit 1
fi

if [[ -z "${outstem}" ]]; then
    echo \
        "Error: Failed to retrieve outstem for id_tsk=${id_tsk}:" \
        "\${arr_outstems[${idx}]}." >&2
    exit 1
fi

if [[ -z "${scl_fct}" ]]; then
    echo \
        "Error: Failed to retrieve scl_fct for id_tsk=${id_tsk}:" \
        "\${arr_scl_fct[${idx}]}." >&2
    exit 1
fi

if [[ -z "${usr_frg}" ]]; then
    echo \
        "Error: Failed to retrieve usr_frg for id_tsk=${id_tsk}:" \
        "\${arr_usr_frg[${idx}]}." >&2
    exit 1
fi

#  Determine BAM infiles' sequenced-read status and file type, ensuring
#+ consistency
typ_seq_num=$(check_typ_seq "${fil_num}")
typ_seq_den=$(check_typ_seq "${fil_den}")
if [[ "${typ_seq_num}" == "${typ_seq_den}" ]]; then
    typ_seq="${typ_seq_num}"
else
    echo \
        "Error: Sequenced-read statuses do not match between numerator" \
        "(${typ_seq_num}) and denominator (${typ_seq_den})." >&2
    exit 1
fi

typ_fil_num="${fil_num##*.}"
typ_fil_den="${fil_den##*.}"
if [[ "${typ_fil_num}" == "${typ_fil_den}" ]]; then
    typ_fil="${typ_fil_num}"
else
    echo \
        "Error: File extensions do not match between numerator" \
        "(${typ_fil_num}) and denominator (${typ_fil_den})." >&2
    exit 1
fi

#  Debug BAM infiles sequenced-read status and file type
if ${debug}; then
    echo "typ_seq=${typ_seq}"
    echo ""
    echo "typ_fil=${typ_fil}"
    echo ""
fi


