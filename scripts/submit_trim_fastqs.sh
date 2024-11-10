#!/bin/bash

#  Enable debug mode if true (optional)
debug=true

#  Display help message if no arguments or help option is given
show_help=$(cat << EOM
\${1}=env_nam  # Mamba environment to activate; use "default" for default env.
\${2}=threads  # Number of threads to use.
\${3}=infiles  # Semicolon-separated serialized string of FASTQ infiles.
\${4}=dir_out  # Directory for output FASTQ files.
\${5}=sfx_se   # Suffix to strip from SE FASTQ files.
\${6}=sfx_pe   # Suffix to strip from PE FASTQ files.
\${7}=err_out  # Directory to store stderr and stdout output files.
\${8}=nam_job  # Job name.
EOM
)

if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    cat << EOM
$(basename "${0}") requires 8 positional arguments:
${show_help}
EOM
    exit 0
fi

#  Check for exactly 8 arguments
if [[ $# -ne 8 ]]; then
    msg="but $# were supplied."
    [[ $# -eq 1 ]] && msg="but only $# was supplied."
    cat << EOM
Error: $(basename "${0}") requires 8 positional arguments, ${msg}

The necessary positional arguments:
${show_help}
EOM
    exit 1
fi

#  Parse positional arguments, assigning them to variables; most of the
#+ argument inputs are not checked, as this is performed by execute_*.sh and to
#+ a certain extent by the scripts submitted to SLURM
env_nam="${1}"
threads="${2}"
infiles="${3}"
dir_out="${4}"
sfx_se="${5}"
sfx_pe="${6}"
err_out="${7}"
nam_job="${8}"

#  Activate environment
if [[ "${env_nam,,}" == "default" ]]; then
    env_nam="env_analyze"
fi

if [[ "${CONDA_DEFAULT_ENV}" != "${env_nam}" ]]; then
    eval "$(conda shell.bash hook)"
    conda activate "${env_nam}" ||
        {
            echo "Error: Failed to activate environment '${env_nam}'." >&2
            exit 1
        }
fi

#  Check that SLURM environment variables are set
if [[ -z "${SLURM_ARRAY_JOB_ID}" || -z "${SLURM_ARRAY_TASK_ID}" ]]; then
    echo "Error: SLURM_ARRAY_JOB_ID and SLURM_ARRAY_TASK_ID are not set." >&2
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

#  Debug output to verify FASTQ file paths and sample name
if ${debug}; then
    echo "fq_1=${fq_1}"
    echo ""
    echo "fq_2=${fq_2}"
    echo ""
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

#  Run Atria
# shellcheck disable=SC2046,2086
atria \
    -t ${threads} \
    -r "${fq_1}" \
    $(if [[ -n ${fq_2} ]]; then echo "-R ${fq_2}"; fi) \
    -o "${dir_out}" \
    --length-range 35:500

#  Remove the initial SLURM stderr and stdout TXT outfiles
rm "${err_ini}" "${out_ini}"
