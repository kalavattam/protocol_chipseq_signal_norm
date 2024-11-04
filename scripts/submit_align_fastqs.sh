#!/bin/bash

#  If true, run script in debug mode
debug=true

#  Parse keyword arguments, assigning them to variables; most of the argument
#+ inputs are not checked, as this is performed by execute_*.sh and again by
#+ the function submitted to SLURM
env_nam="env_align"
sc_func=""
threads=4
aligner="bowtie2"
a_type="end-to-end"
mapq=""
req_flg=false
index=""
infiles=""
dir_out=""
qname=false
sfx_se=""
sfx_pe=""
err_out=""
nam_job="align_fastqs"

show_help=$(cat << EOM
$(basename "${0}") takes the following keyword arguments:
  -en, --env_nam  Mamba environment to activate.
  -sf, --sc_func  Path and name of function to submit to SLURM.
   -t, --threads  Number of threads to use.
   -a, --aligner  Alignment program to use.
  -at, --a_type   Alignment type.
  -mq, --mapq     MAPQ threshold for filtering BAM outfiles.
  -rf, --req_flg  Require flag bit 2 for filtering BAM outfiles.
  -ix, --index    Path to directory containing aligner index.
   -i, --infiles  Semicolon-delimited serialized string of FASTQ infiles.
  -do, --dir_out  Directory to write BAM outfiles.
  -qn, --qname    Retain queryname-sorted intermediate BAM files.
  -ss, --sfx_se   Suffix to strip from SE FASTQ files.
  -sp, --sfx_pe   Suffix to strip from PE FASTQ files.
  -eo, --err_out  Directory to store stderr and stdout outfiles.
  -nj, --nam_job  Name of job.

All arguments are required with the following notes and exceptions:
  - --env_nam defaults to env_nam=${env_nam} if not specified.
  - --threads default to threads=${threads} if not specified.
  - --aligner defaults to aligner=${aligner} if not specified.
  - --a_type defaults to a_type=${a_type} if not specified.
  - --a_type is only required if aligner=bowtie2; if aligner=bwa, then the
    --a_type assignment is ignored.
  - --req_flg is optional.
  - --qname is optional.
  - --nam_job defaults to nam_job=${nam_job} if not specified.
EOM
)

if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case "${1}" in
        -en|--env_nam) env_nam="${2}";   shift 2 ;;
        -sf|--sc_func) sc_func="${2}";   shift 2 ;;
         -t|--threads) threads="${2}";   shift 2 ;;
         -a|--aligner) aligner="${2,,}"; shift 2 ;;
        -at|--a_type)  a_type="${2,,}";  shift 2 ;;
        -mq|--mapq)    mapq="${2}";      shift 2 ;;
        -rf|--req_flg) req_flg=true;     shift 1 ;;
        -ix|--index)   index="${2}";     shift 2 ;;
         -i|--infiles) infiles="${2}";   shift 2 ;;
        -do|--dir_out) dir_out="${2}";   shift 2 ;;
        -qn|--qname)   qname=true;       shift 1 ;;
        -ss|--sfx_se)  sfx_se="${2}";    shift 2 ;;
        -sp|--sfx_pe)  sfx_pe="${2}";    shift 2 ;;
        -eo|--err_out) err_out="${2}";   shift 2 ;;
        -nj|--nam_job) nam_job="${2}";   shift 2 ;;
        *)
            echo "## Unknown parameter passed: ${1} ##" >&2
            echo "" >&2
            echo "${show_help}" >&2
            exit 1
            ;;
    esac
done

#  Activate environment
if [[ "${CONDA_DEFAULT_ENV}" != "${env_nam}" ]]; then
    eval "$(conda shell.bash hook)"
    conda activate "${env_nam}" ||
        {
            echo "Error: Failed to activate environment '${env_nam}'." >&2
            exit 1
        }
fi

#  Source function to submit to SLURM
if [[ ! -f "${sc_func}" ]]; then
    echo "Error: Function script not found: ${sc_func}." >&2
    exit 1
else
    # shellcheck disable=SC1090
    source "${sc_func}" ||
        {
            echo "Error: Failed to source function: ${sc_func}." >&2
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

#  Reconstruct arr_infiles from the serialized string assigned to infiles
IFS=';' read -r -a arr_infiles <<< "${infiles}"

#  Debug output to check the task ID/array index, number of array elements, and
#+ array element values
if ${debug}; then
    echo "SLURM_ARRAY_TASK_ID=${id_tsk}"
    echo ""
    echo "\${#arr_infiles[@]}=${#arr_infiles[@]}"
    echo ""
    echo "arr_infiles=( ${arr_infiles[*]} )"
    echo ""
fi

#  Determine infile assignment based on SLURM_ARRAY_TASK_ID
infile="${arr_infiles[$(( id_tsk - 1 ))]}"

#  Debug output to check which file or file pair is being processed
if ${debug}; then
    echo "infile=${infile}"
    echo ""
fi

#  Exit if variable infile has an empty assignment
if [[ -z "${infile}" ]]; then
    echo \
        "Error: Failed to retrieve infile for id_tsk=${id_tsk}:" \
        "\${arr_infiles[$(( id_tsk - 1 ))]}." >&2
    exit 1
fi

#  Process the infile assignment
if [[ "${infile}" == *,* ]]; then
    #  If infile is a file pair, split the pair to variables fq_1 and fq_2, and
    #+ assign a sample stem variable too
    fq_1="${infile%%,*}"
    fq_2="${infile#*,}"
    samp="$(basename "${fq_1%%"${sfx_pe}"}")"
else
    #  Assign infile to fq_1 if the infile assignment is a single file; assign
    #+ a sample stem variable too
    fq_1="${infile}"
    unset fq_2
    samp="$(basename "${fq_1%%"${sfx_se}"}")"
fi

#  Assign the outfile based on the dir_out and infile assignments
outfile="${dir_out}/${samp}.bam"

#  Debug output to verify FASTQ file paths and sample name
if ${debug}; then
    echo "fq_1=${fq_1}"
    echo ""
    echo "fq_2=${fq_2}"
    echo ""
    echo "samp=${samp}"
    echo ""
    echo "outfile=${outfile}"
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

#  Run the function
# shellcheck disable=SC2046,2086
$(basename ${sc_func} .sh) \
    --threads ${threads} \
    --aligner ${aligner} \
    $(if [[ -n ${a_type} ]]; then echo "--a_type ${a_type}"; fi) \
    $(if [[ ${mapq} -gt 0 ]]; then echo "--mapq ${mapq}"; fi) \
    $(if ${req_flg}; then echo "--req_flg"; fi) \
    --index ${index} \
    --fq_1 ${fq_1} \
    $(if [[ -n ${fq_2} ]]; then echo "--fq_2 ${fq_2}"; fi) \
    --outfile ${outfile} \
    $(if ${qname}; then echo "--qname"; fi)

#  Remove the initial SLURM stderr and stdout TXT outfiles
rm "${err_ini}" "${out_ini}"
