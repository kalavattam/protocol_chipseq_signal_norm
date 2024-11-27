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
str_infile=""
str_outfile=""
typ_out="bigwig"
bin_siz=1
region=""
str_scl_fct=""
norm=""
exact=false
str_usr_frg=""
err_out=""
nam_job="bamCoverage"
env_nam="env_align"

show_help=$(cat << EOM
$(basename "${0}")
  [--verbose] --threads <int> --str_infile <str> --str_outfile <str>
  --typ_out <str> --bin_siz <int> [--region <str>] --str_scl_fct <str>
  --norm <str> [--exact] [--str_usr_frg <str>] --err_out <str> --nam_job <str>
  --env_nam <str>

$(basename "${0}") takes the following keyword arguments:
   -v, --verbose      Run in 'verbose mode' (default: ${verbose}).
   -t, --threads      Number of threads to use (default: ${threads}).
  -si, --str_infile   Comma-separated string of BAM infiles.
  -so, --str_outfile  Comma-separated string of outfile stems.
  -to, --typ_out      Outfile type: 'bedgraph' or 'bigwig' (default: ${typ_out}).
  -bs, --bin_siz      Bin size in base pairs (default: ${bin_siz}).
   -r, --region       Region in 'chr' or 'chr:start-stop' format.
  -sf, --str_scl_fct  Comma-separated string of scaling factors. Cannot be
                      used with --norm.
  -no, --norm         Use one of the following normalization methods when
                      computing coverage: 'None', 'RPKM', 'CPM', 'BPM', or
                      'RPGC'. Cannot be used with --str_scl_fct.
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
        -si|--str_infile)  str_infile="${2}";  shift 2 ;;
        -so|--str_outfile) str_outfile="${2}"; shift 2 ;;
        -to|--typ_out)     typ_out="${2}";     shift 2 ;;
        -bs|--bin_siz)     bin_siz="${2}";     shift 2 ;;
         -r|--region)      region="${2}";      shift 2 ;;
        -sf|--str_scl_fct) str_scl_fct="${2}"; shift 2 ;;
        -no|--norm)        norm="${2}";        shift 2 ;;
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
    "verbose" "threads" "str_infile" "str_outfile" "typ_out" "bin_siz" "region"
    "str_scl_fct" "norm" "exact" "str_usr_frg" "err_out" "nam_job" "env_nam"
)
for var in "${arr_arg_req[@]}"; do
    if [[ -z "${!var}" ]]; then
        echo "Error: ${var} is required." >&2
        echo "" >&2
        echo "${show_help}" >&2
        exit 1
    fi
done

#  Validate specification of '--norm' or '--str_scl_fct'
if [[ "${norm}" == "#N/A" ]]; then
    #  If '--norm' is not specified, ensure '--str_scl_fct' is valid
    if [[
        -z "${str_scl_fct}" || "${str_scl_fct}" =~ (^|,)?'#N/A'(,|$)
    ]]; then
        echo \
            "Error: When '--norm <str>' is not specified, '--str_scl_fct'" \
            "must be provided and cannot contain '#N/A' elements." >&2
        exit 1
    fi
else
    #  If '--norm' is specified, ensure '--str_scl_fct' is not specified or
    #+ invalid
    if [[
        -n "${str_scl_fct}" && ! "${str_scl_fct}" =~ (^|,)?'#N/A'(,|$)
    ]]; then
        if [[ "${norm}" != "#N/A" ]]; then
            echo \
                "Error: When '--norm <str>' is specified, '--str_scl_fct' must" \
                "not be provided or must only contain '#N/A' elements." >&2
            exit 1
        fi
    fi
fi

#  Debug argument assignments
if ${debug}; then
    echo "verbose=${verbose}"
    echo ""
    echo "threads=${threads}"
    echo ""
    echo "str_infile=${str_infile}"
    echo ""
    echo "str_outfile=${str_outfile}"
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
IFS=',' read -r -a arr_outfiles <<< "${str_outfile}"
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
    echo "\${#arr_outfiles[@]}=${#arr_outfiles[@]}"
    echo ""
    echo "arr_outfiles=( ${arr_outfiles[*]} )"
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
outfile="${arr_outfiles[idx]}"
scl_fct="${arr_scl_fct[idx]}"
usr_frg="${arr_usr_frg[idx]}"

#  Debug variable assignments from reconstructed arrays
if ${debug}; then
    echo "infile=${infile}"
    echo ""
    echo "outfile=${outfile}"
    echo ""
    echo "scl_fct=${scl_fct}"
    echo ""
    echo "usr_frg=${usr_frg}"
    echo ""
fi

#  Exit if any variable assignment is empty
if [[ -z "${infile}" ]]; then
    echo \
        "Error: Failed to retrieve infile for id_tsk=${id_tsk}:" \
        "\${arr_infiles[${idx}]}." >&2
    exit 1
fi

if [[ -z "${outfile}" ]]; then
    echo \
        "Error: Failed to retrieve outfile for id_tsk=${id_tsk}:" \
        "\${arr_outfiles[${idx}]}." >&2
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

#  Determine BAM infile's sequenced-read status: 'single' or 'paired'
typ_seq=$(check_typ_seq "${infile}")

#  Debug BAM infile's sequenced-read status
if ${debug}; then
    echo "typ_seq=${typ_seq}"
    echo ""
fi

#  Derive sample name from infile assignment
samp="${outfile##*/}"
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
    echo "    --outFileName ${outfile}.${typ_out} \\"
    echo "    --outFileFormat ${typ_out} \\"
    echo "    --binSize ${bin_siz} \\"
    echo "    $(
        if [[ -n "${region}" && "${norm}" != "#N/A" ]]; then
            echo "--region ${region} \\"
        fi
    )"
    echo "    --skipNonCoveredRegions \\"
    echo "    $(
        if [[ "${scl_fct}" != "#N/A" && "${norm}" == "#N/A" ]]; then
            echo "--scaleFactor ${scl_fct} \\"
        fi
    )"
    echo "    $(
        if [[ "${scl_fct}" == "#N/A" && "${norm}" != "#N/A" ]]; then
            echo "--normalizeUsing ${norm} \\"
        fi
    )"
    echo "    $(
        if [[ "${norm}" == "RPGC" ]]; then
            echo "--effectiveGenomeSize ${gen_siz:-11624332} \\"
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
    --outFileName "${outfile}.${typ_out}" \
    --outFileFormat "${typ_out}" \
    --binSize "${bin_siz}" \
    $(
        if [[ -n "${region}" && "${region}" != "#N/A" ]]; then
            echo "--region ${region}"
        fi
    ) \
    --skipNonCoveredRegions \
    $(
        if [[ "${scl_fct}" != "#N/A" && "${norm}" == "#N/A" ]]; then
            echo "--scaleFactor ${scl_fct}"
        fi
    ) \
    $(
        if [[ "${scl_fct}" == "#N/A" && "${norm}" != "#N/A" ]]; then
            echo "--normalizeUsing ${norm}"
        fi
    ) \
    $(
        if [[ "${norm}" == "RPGC" ]]; then
            echo "--effectiveGenomeSize ${gen_siz:-11624332}"
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

#  - S. cerevisiae effective genome size if no MAPQ filtering is performed,
#+   i.e., all multimapping alignments are retained: 12157105 (via faCount)
#+ - S. cerevisiae effective genome size if any MAPQ filtering is performed,
#+   i.e., no or only a subset of multimapping alignments are retained:
#+   11624332 (via khmer unique-kmers.py for 50mers)
