#!/bin/bash

#  If true, run script in debug mode
debug=true

#  Parse keyword arguments, assigning them to variables; most of the argument
#+ inputs are not checked, as this is performed by execute_*.sh and again by
#+ the function submitted to SLURM
env_nam="env_align"
threads=4
infiles=""
dir_out=""
mito=false
tg=false
mtr=false
chk_chr=false
err_out=""
nam_job=""
scr_fnc=""

show_help=$(cat << EOM
$(basename "${0}") takes the following keyword arguments:
  -en, --env_nam  Mamba environment to activate.
   -t, --threads  Number of threads to use.
   -i, --infiles  Comma-separated serialized string of FASTQ infiles.
  -do, --dir_out  Directory to write BAM outfiles.
   -m, --mito     Retain mitochondrial chromosome in BAM outfiles.
  -tg, --tg       Retain SP_II_TG chromosome in BAM outfiles.
  -mr, --mtr      Retain SP_MTR chromosome in BAM outfiles.
  -cc, --chk_chr  Check chromosomes in BAM outfiles.
  -eo, --err_out  Directory to store stderr and stdout outfiles.
  -nj, --nam_job  Name of job.
  -sf, --scr_fnc  Path and name of function to submit to SLURM.

All arguments are required with the following notes and exceptions:
  - --env_nam defaults to env_nam=${env_nam} if not specified.
  - --threads default to threads=${threads} if not specified.
  - --mito, --tg, --mtr, and --chk_chr are optional.
EOM
)

if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case "${1}" in
        -en|--env_nam) env_nam="${2}";   shift 2 ;;
         -t|--threads) threads="${2}";   shift 2 ;;
         -i|--infiles) infiles="${2}";   shift 2 ;;
        -do|--dir_out) dir_out="${2}";   shift 2 ;;
         -m|--mito)    mito=true;        shift 1 ;;
        -tg|--tg)      tg=true;          shift 1 ;;
        -mr|--mtr)     mtr=true;         shift 1 ;;
        -cc|--chk_chr) chk_chr=true;     shift 1 ;;
        -eo|--err_out) err_out="${2}";   shift 2 ;;
        -nj|--nam_job) nam_job="${2}";   shift 2 ;;
        -sf|--scr_fnc) scr_fnc="${2}";   shift 2 ;;
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
if [[ ! -f "${scr_fnc}" ]]; then
    echo "Error: Function script not found: ${scr_fnc}." >&2
    exit 1
else
    # shellcheck disable=SC1090
    source "${scr_fnc}" ||
        {
            echo "Error: Failed to source function: ${scr_fnc}." >&2
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
IFS=',' read -r -a arr_infiles <<< "${infiles}"

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

#  Derive sample name from infile assignment
samp="${infile##*/}"
samp="${samp%.bam}"

#  Derive function name from function script
nam_fnc="$(basename "${scr_fnc}" ".sh")"

#  Perform pattern matching to assign the outfile name
# shellcheck disable=SC2086
if [[ ${nam_fnc} =~ "sc" ]]; then
    outfile="${dir_out}/${samp}.sc.bam"
elif [[ ${nam_fnc} =~ "sp" ]]; then
    outfile="${dir_out}/${samp}.sp.bam"
else
    echo "Error: Sample name could not be processed." >&2
    exit 1
fi

#  Debug output to verify sample name and outfile path
if ${debug}; then
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
${nam_fnc} \
    --threads ${threads} \
    --infile ${infile} \
    --outfile ${outfile} \
    $(if ${mito}; then echo "--mito"; fi) \
    $(if ${tg}; then echo "--tg"; fi) \
    $(if ${mtr}; then echo "--mtr"; fi) \
    $(if ${chk_chr}; then echo "--chk_chr"; fi)

#  Remove the initial SLURM stderr and stdout TXT outfiles
rm "${err_ini}" "${out_ini}"
