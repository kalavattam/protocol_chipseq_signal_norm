#!/bin/bash

#  submit_qsort_bam_slurm.sh
#  KA

#  If true, run script in debug mode
debug=true

#  Display help message if no arguments or help option is given
show_help=$(cat << EOM
\${1}=env_nam     # str: Name of Conda/Mamba environment to activate
\${2}=threads     # int: Number of threads to use
\${3}=str_infile  # str: Comma-separated list of BAM infiles
\${4}=dir_out     # str: Directory to save QNAME-sorted BAM outfiles
\${5}=err_out     # str: Directory for stdout and stderr files
\${6}=nam_job     # str: Name of job
EOM
)

if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    cat << EOM
$(basename "${0}") requires 6 positional arguments:
${show_help}
EOM
    exit 0
fi

#  Check for exactly 6 positional arguments
if [[ $# -ne 6 ]]; then
    msg="but $# were supplied."
    [[ $# -eq 1 ]] && msg="but only $# was supplied."
    cat << EOM
Error: $(basename "${0}") requires 6 positional arguments, ${msg}

The necessary positional arguments:
${show_help}
EOM
    exit 1
fi

#  Parse positional arguments, assigning them to variables
env_nam="${1}"
threads="${2}"
str_infile="${3}"
dir_out="${4}"
err_out="${5}"
nam_job="${6}"

#  Debug positional argument assignments
if ${debug}; then
    echo "\${env_nam}=${env_nam}"
    echo ""
    echo "\${threads}=${threads}"
    echo ""
    echo "\${str_infile}=${str_infile}"
    echo ""
    echo "\${dir_out}=${dir_out}"
    echo ""
    echo "\${err_out}=${err_out}"
    echo ""
    echo "\${nam_job}=${nam_job}"
    echo ""
fi

#  Validate specified directories
for dir in "$(dirname "${str_infile%%,*}")" "${dir_out}" "${err_out}"; do
    if [[ ! -d "${dir}" ]]; then
        echo "Error: Directory does not exist: '${dir}'" >&2
        exit 1
    fi
done

#  Validate number of threads
if ! [[ "${threads}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: Threads argument must be a positive integer: '${threads}'" >&2
    exit 1
fi

#  Activate environment
if [[ "${CONDA_DEFAULT_ENV}" != "${env_nam}" ]]; then
    eval "$(conda shell.bash hook)"
    conda activate "${env_nam}" ||
        {
            echo "Error: Failed to activate environment: '${env_nam}'" >&2
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

#  Reconstruct array from serialized string
IFS=',' read -r -a arr_infiles  <<< "${str_infile}"

#  Validate supplied infiles
for file in "${arr_infiles[@]}"; do
    if [[ ! -f "${file}" ]]; then
        echo "Error: File does not exist: '${file}'" >&2
        exit 1
    fi
done

if [[ "${id_tsk}" -gt "${#arr_infiles[@]}" || "${id_tsk}" -lt 1 ]]; then
    echo "Error: SLURM_ARRAY_TASK_ID is out of bounds: '${id_tsk}'" >&2
    exit 1
fi

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

#  Determine array index based on SLURM_ARRAY_TASK_ID
idx=$(( id_tsk - 1 ))

#  Assign variables from reconstructed array based on index
infile="${arr_infiles[idx]}"
outfile="${dir_out}/$(basename "${infile}" ".bam").qnam.bam"

#  Debug variable assignments from reconstructed array
if ${debug}; then
    echo "infile=${infile}"
    echo ""
    echo "outfile=${outfile}"
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

#  Derive sample name from outfile assignment
samp="$(basename "${outfile}" ".qnam.bam")_qnam"

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

#  Sort the BAM file by QNAME
samtools sort -@ "${threads}" -n -o "${outfile}" "${infile}" ||
    {
        echo "Error: Failed to sort file: '${infile}'." >&2
        exit 1
    }

echo "## Successfully sorted '${infile}' to '${outfile}' ##"

#  Remove the initial SLURM stderr and stdout TXT outfiles
rm "${err_ini}" "${out_ini}"
