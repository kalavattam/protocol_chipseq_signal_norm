#!/bin/bash

set -eo pipefail

#  If true, run script in debug mode
debug=true

#  Define helper function to validate and log errors for variables
function validate_var() {
    local var_nam="${1}"  # Name of variable
    local var_val="${2}"  # Value of variable
    local arr_ref="${3}"  # Reference to array (as string)
    local idx="${4}"      # Index of array element
    local id_tsk="${5}"   # Task ID

    if [[ -z "${var_val}" ]]; then
        echo \
            "Task ${id_tsk}: Error retrieving ${var_nam}. Value at" \
            "\${${arr_ref}[${idx}]} is empty." >&2
        echo "" >&2
        echo \
            "Hint: Verify that inputs to '${arr_ref}' are correctly" \
            "formatted and align with the expected index (${idx})." >&2
        return 1
    fi
}


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
#+ execute_*.sh and again by the deepTools program submitted to SLURM
verbose=false
threads=8
str_infile=""
str_outstem=""
typ_out="bedgraph"
siz_bin=1
siz_gen=12157105
region=""
str_scl_fct=""
typ_cvg="None"
exact=false
str_usr_frg=""
err_out=""
nam_job="deeptools_coverage"
env_nam="env_align"

show_help=$(cat << EOM
$(basename "${0}")
  [--verbose] --threads <int> --str_infile <str> --str_outstem <str>
  --typ_out <str> --siz_bin <int> [--siz_gen <int>] [--region <str>]
  --str_scl_fct <str> --typ_cvg <str> [--exact] [--str_usr_frg <str>]
  --err_out <str> --nam_job <str> --env_nam <str>

$(basename "${0}") takes the following keyword arguments:
   -v, --verbose      Run in 'verbose mode' (optional).
   -t, --threads      Number of threads to use (optional).
  -si, --str_infile   Comma-separated string of BAM infiles.
  -so, --str_outstem  Comma-separated string of outfile stems.
  -to, --typ_out      Outfile type: 'bedgraph' or 'bigwig' (default:
                      "${typ_out}").
  -sb, --siz_bin      Bin size in base pairs (default: ${siz_bin}).
  -sg, --siz_gen      Effective genome size of model organism (default:
                      ${siz_gen}).
   -r, --region       Region in 'chr' or 'chr:start-stop' format (optional).
  -sf, --str_scl_fct  Comma-separated string of scaling factors.
  -tv, --typ_cvg      Use one of the following normalization methods when
                      computing coverage: 'None', 'RPKM', 'CPM', 'BPM', or
                      'RPGC' (default: "${typ_cvg}").
  -ex, --exact        Compute scaling factors based on all alignments. Only
                      applicable if '--typ_cvg <str != None>' is specified.
                      Significantly slows coverage computation (optional).
  -su, --str_usr_frg  Comma-separated string of fragment lengths for alignment
                      extension (optional).
  -eo, --err_out      Directory to store stderr and stdout outfiles.
  -nj, --nam_job      Name of job (default "${nam_job}").
  -en, --env_nam      Mamba environment to activate (default: "${env_nam}").
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
        -si|--str_infile)  str_infile="${2}";  shift 2 ;;
        -so|--str_outstem) str_outstem="${2}"; shift 2 ;;
        -to|--typ_out)     typ_out="${2}";     shift 2 ;;
        -sb|--siz_bin)     siz_bin="${2}";     shift 2 ;;
        -sg|--siz_gen)     siz_gen="${2}";     shift 2 ;;
         -r|--region)      region="${2}";      shift 2 ;;
        -sf|--str_scl_fct) str_scl_fct="${2}"; shift 2 ;;
        -tv|--typ_cvg)     typ_cvg="${2}";     shift 2 ;;
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
    "verbose" "threads" "str_infile" "str_outstem" "typ_out" "siz_bin"
    "siz_gen" "region" "str_scl_fct" "typ_cvg" "exact" "str_usr_frg" "err_out"
    "nam_job" "env_nam"
)
for var in "${arr_arg_req[@]}"; do
    if [[ -z "${!var}" ]]; then
        echo "Error: '--${var}' is required." >&2
        echo "" >&2
        echo "${show_help}" >&2
        exit 1
    fi
done

#  Validate '--str_scl_fct' if specified
if [[ -n "${str_scl_fct}" ]]; then
    IFS=',' read -r -a arr_scl_fct <<< "${str_scl_fct}"
    for s in "${arr_scl_fct[@]}"; do
        if [[ "${s}" == "#N/A" ]]; then continue; fi
        if ! [[ "${s}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
            echo \
                "Error: Invalid scaling factor in '--str_scl_fct':" \
                "'${s}'." >&2
            exit 1
        fi
    done
fi

#  Validate '--typ_cvg' if specified
if [[ -n "${typ_cvg}" ]]; then
    case "${typ_cvg}" in
        None|RPKM|CPM|BPM|RPGC) : ;;
        "#N/A") typ_cvg="None" ;;
        *)
            echo \
                "Error: Unsupported normalization method in '--typ_cvg':" \
                "'${typ_cvg}'." >&2
            exit 1
            ;;
    esac
fi

#  Validate and standardize output file type
case "${typ_out}" in
          bigwig|bw) typ_out="bigwig"   ;;
    bedgraph|bdg|bg) typ_out="bedgraph" ;;
    *)
        echo \
            "Error: Unsupported output type: '${typ_out}'. Supported types:" \
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
    echo "str_infile=${str_infile}"
    echo ""
    echo "str_outstem=${str_outstem}"
    echo ""
    echo "typ_out=${typ_out}"
    echo ""
    echo "siz_bin=${siz_bin}"
    echo ""
    echo "region=${region}"
    echo ""
    echo "str_scl_fct=${str_scl_fct}"
    echo ""
    echo "typ_cvg=${typ_cvg}"
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
IFS=',' read -r -a arr_infiles  <<< "${str_infile}"
IFS=',' read -r -a arr_outstems <<< "${str_outstem}"
IFS=',' read -r -a arr_scl_fct  <<< "${str_scl_fct}"
IFS=',' read -r -a arr_usr_frg  <<< "${str_usr_frg}"

#  Debug output to check the task ID/array index, number of array elements, and
#+ array element values
if ${debug}; then
    echo "SLURM_ARRAY_TASK_ID=${id_tsk}"
    echo ""
    echo "\${#arr_infiles[@]}=${#arr_infiles[@]}"
    echo ""
    echo "arr_infiles=( ${arr_infiles[*]} )"
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
infile="${arr_infiles[idx]}"
outstem="${arr_outstems[idx]}"
scl_fct="${arr_scl_fct[idx]}"
usr_frg="${arr_usr_frg[idx]}"

#  Validate array index is not out of bounds
if (( idx >= ${#arr_infiles[@]} )); then
    echo "Error: SLURM task ID ${id_tsk} exceeds number of input files." >&2
    exit 1
fi

#  Debug variable assignments from reconstructed arrays
if ${debug}; then
    echo "infile=${infile}"
    echo ""
    echo "outstem=${outstem}"
    echo ""
    echo "scl_fct=${scl_fct}"
    echo ""
    echo "usr_frg=${usr_frg}"
    echo ""
fi

#  Exit if any variable assignment is empty
validate_var "infile"  "${infile}"  "arr_infiles"  "${idx}" "${id_tsk}"
validate_var "outstem" "${outstem}" "arr_outstems" "${idx}" "${id_tsk}"
validate_var "scl_fct" "${scl_fct}" "arr_scl_fct"  "${idx}" "${id_tsk}"
validate_var "usr_frg" "${usr_frg}" "arr_usr_frg"  "${idx}" "${id_tsk}"

#  Determine BAM infile's sequenced-read status: 'single' or 'paired'
typ_seq=$(check_typ_seq "${infile}")

#  Debug BAM infile's sequenced-read status
if ${debug}; then
    echo "typ_seq=${typ_seq}"
    echo ""
fi

#  Derive sample name from outstem assignment
samp="${outstem##*/}"
samp="${samp%.bam}"

#  Debug sample name
if ${debug}; then
    echo "samp=${samp}"
    echo ""
fi

#  Assign stdout and stderr outfiles to variables
err_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stderr.txt"
out_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stdout.txt"
err_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stderr.txt"
out_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stdout.txt"

#  Give the initial stderr and stdout TXT outfiles more descriptive names
ln -f "${err_ini}" "${err_dsc}"
ln -f "${out_ini}" "${out_dsc}"

#  Debug execution of deepTools bamCoverage
if ${debug}; then
    echo "bamCoverage \\"
    echo "    $(if ${verbose}; then echo "--verbose \\"; fi)"
    echo "    --numberOfProcessors ${threads} \\"
    echo "    --bam ${infile} \\"
    echo "    --outFileName ${outstem}.${typ_out} \\"
    echo "    --outFileFormat ${typ_out} \\"
    echo "    --binSize ${siz_bin} \\"
    echo "    $(
        if [[ -n "${region}" && "${region}" != "#N/A" ]]; then
            echo "--region ${region} \\"
        fi
    )"
    echo "    --skipNonCoveredRegions \\"
    echo "    $(
        if [[ -n "${scl_fct}" && "${scl_fct}" != "#N/A" ]]; then
            echo "--scaleFactor ${scl_fct} \\"
        fi
    )"
    echo "    $(
        if [[ -n "${typ_cvg}" ]]; then
            echo "--normalizeUsing ${typ_cvg} \\"
        fi
    )"
    echo "    $(
        if [[ "${typ_cvg}" == "RPGC" ]]; then
            echo "--effectiveGenomeSize ${siz_gen} \\"
        fi
    )"
    echo "    $(if ${exact}; then echo "--exactScaling \\"; fi)"
    echo "    $(
        if [[ "${typ_seq}" == "paired" ]]; then
            echo "--samFlagInclude 64 \\"
        fi
    )"
    echo "    $(
        if [[ "${typ_seq}" == "paired" && "${usr_frg}" == "#N/A" ]]; then
            echo "--extendReads"
        elif [[ "${usr_frg}" != "#N/A" ]]; then
            echo "--extendReads ${usr_frg}"
        fi
    )"
fi

#  Execute deepTools bamCoverage
# shellcheck disable=SC2046,SC2154
bamCoverage \
    $(if ${verbose}; then echo "--verbose"; fi) \
    --numberOfProcessors "${threads}" \
    --bam "${infile}" \
    --outFileName "${outstem}.${typ_out}" \
    --outFileFormat "${typ_out}" \
    --binSize "${siz_bin}" \
    $(
        if [[ -n "${region}" && "${region}" != "#N/A" ]]; then
            echo "--region ${region}"
        fi
    ) \
    --skipNonCoveredRegions \
    $(
        if [[ -n "${scl_fct}" && "${scl_fct}" != "#N/A" ]]; then
            echo "--scaleFactor ${scl_fct}"
        fi
    ) \
    $(
        if [[ -n "${typ_cvg}" ]]; then
            echo "--normalizeUsing ${typ_cvg}"
        fi
    ) \
    $(
        if [[ "${typ_cvg}" == "RPGC" ]]; then
            echo "--effectiveGenomeSize ${siz_gen}"
        fi
    ) \
    $(if ${exact}; then echo "--exactScaling"; fi) \
    $(
        if [[ "${typ_seq}" == "paired" ]]; then
            echo "--samFlagInclude 64"
        fi
    ) \
    $(
        if [[ "${typ_seq}" == "paired" && "${usr_frg}" == "#N/A" ]]; then
            echo "--extendReads"
        elif [[ "${usr_frg}" != "#N/A" ]]; then
            echo "--extendReads ${usr_frg}"
        fi
    )

#  Remove the initial SLURM stderr and stdout TXT outfiles
rm "${err_ini}" "${out_ini}"

#  - S. cerevisiae effective genome size if multimapping alignments are
#+   retained: 12157105 (via faCount)
#+ - S. cerevisiae effective genome size if MAPQ filtering is performed,
#+   i.e., no multimapping alignments are retained: 11624332 (via khmer
#+   unique-kmers.py for 50mers)

#TODO: Implement GNU Parallel and serial job submissions
