#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: submit_qsort_bam_slurm.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


#  Require Bash >= 4.4 before doing any work
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be run under Bash >= 4.4." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

#  Run in safe mode, exiting on errors, unset variables, and pipe failures
set -euo pipefail

#  If true, run script in debug mode
debug=true

#  Display help message if no arguments or help option is given
show_help=$(cat << EOM
\${1}=env_nam     <str>  Name of Conda/Mamba environment to activate.
\${2}=threads     <int>  Number of threads to use.
\${3}=csv_infile  <str>  Comma-separated list of BAM infiles.
\${4}=dir_out     <str>  Directory to save QNAME-sorted BAM outfiles.
\${5}=err_out     <str>  Directory for stdout and stderr files.
\${6}=nam_job     <str>  Name of job.
EOM
)

if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
    cat >&2 << EOM
'$(basename "${0}")' requires 6 positional arguments:
${show_help}
EOM
    exit 0
fi

#  Check for exactly 6 positional arguments
if [[ $# -ne 6 ]]; then
    msg="but $# were supplied."
    [[ $# -eq 1 ]] && msg="but only $# was supplied."
    cat >&2 << EOM
error: '$(basename "${0}")' requires 6 positional arguments, ${msg}

The necessary positional arguments:
${show_help}
EOM
    exit 1
fi

#  Parse positional arguments, assigning them to variables
env_nam="${1}"
threads="${2}"
csv_infile="${3}"
dir_out="${4}"
err_out="${5}"
nam_job="${6}"

#  Debug positional argument assignments
if [[ "${debug:-false}" == "true" ]]; then
    printf '%s\n\n' \
        "\${env_nam}=${env_nam}" \
        "\${threads}=${threads}" \
        "\${csv_infile}=${csv_infile}" \
        "\${dir_out}=${dir_out}" \
        "\${err_out}=${err_out}" \
        "\${nam_job}=${nam_job}" \
        >&2
fi

#  Validate specified directories
for dir in "$(dirname "${csv_infile%%,*}")" "${dir_out}" "${err_out}"; do
    if [[ ! -d "${dir}" ]]; then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "directory does not exist: '${dir}'" >&2
        exit 1
    fi
done

#  Validate number of threads
if ! [[ "${threads}" =~ ^[1-9][0-9]*$ ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "threads argument must be a positive integer: '${threads}'" >&2
    exit 1
fi

#  Activate environment
if [[ "${CONDA_DEFAULT_ENV:-}" != "${env_nam}" ]]; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${env_nam}" ||
        {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to activate environment: '${env_nam}'" >&2
            exit 1
        }
fi

#  Check for necessary program
if ! \
    command -v samtools > /dev/null 2>&1
then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "required program is not in PATH: 'samtools'." >&2
    exit 1
fi

#  Check that SLURM environment variables are set
if [[ -z "${SLURM_ARRAY_JOB_ID:-}" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "'SLURM_ARRAY_JOB_ID' is not set." >&2
    exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "'SLURM_ARRAY_TASK_ID' is not set." >&2
    exit 1
fi

#  Give important SLURM environmental variables shorter names
id_job=${SLURM_ARRAY_JOB_ID}
id_tsk=${SLURM_ARRAY_TASK_ID}

#  Reconstruct array from serialized string
IFS=',' read -r -a arr_infiles  <<< "${csv_infile}"

#  Validate supplied infiles
for file in "${arr_infiles[@]}"; do
    if [[ ! -f "${file}" ]]; then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "file does not exist: '${file}'" >&2
        exit 1
    fi
done

#  Validates IDs
if ! [[ "${id_tsk}" =~ ^[1-9][0-9]*$ ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "'SLURM_ARRAY_TASK_ID' must be a positive integer: '${id_tsk}'" >&2
    exit 1
elif [[ "${id_tsk}" -gt "${#arr_infiles[@]}" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "'SLURM_ARRAY_TASK_ID' is out of bounds: '${id_tsk}'" >&2
    exit 1
fi

#  Debug output to check the task ID/array index, number of array elements, and
#+ array element values
if [[ "${debug:-false}" == "true" ]]; then
    printf '%s\n\n' \
        "SLURM_ARRAY_TASK_ID=${id_tsk}" \
        "\${#arr_infiles[@]}=${#arr_infiles[@]}" \
        "arr_infiles=( ${arr_infiles[*]} )" \
        >&2
fi

#  Determine array index based on SLURM_ARRAY_TASK_ID
idx=$(( id_tsk - 1 ))

#  Assign variables from reconstructed array based on index
infile="${arr_infiles[idx]}"
outfile="${dir_out}/$(basename "${infile}" ".bam").qnam.bam"

#  Debug variable assignments from reconstructed array
if [[ "${debug:-false}" == "true" ]]; then
    printf '%s\n\n' \
        "infile=${infile}" \
        "outfile=${outfile}" \
        >&2
fi

#  Exit if any variable assignment is empty
if [[ -z "${infile}" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to derive infile for 'id_tsk=${id_tsk}':" \
        "'\${arr_infiles[${idx}]}'." >&2
    exit 1
fi

if [[ -z "${outfile}" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to derive outfile for 'id_tsk=${id_tsk}':" \
        "'\${arr_infiles[${idx}]}'." >&2
    exit 1
fi

#  Derive sample name from outfile assignment
samp="$(basename "${outfile}" ".qnam.bam")_qnam"

#  Debug sample name
if [[ "${debug:-false}" == "true" ]]; then
    printf '%s\n\n' "samp=${samp}" >&2
fi

#  Assign stdout and stderr outfiles to variables
err_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stderr.txt"
out_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stdout.txt"
err_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stderr.txt"
out_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stdout.txt"

#  Give the initial stderr and stdout TXT outfiles more descriptive names
ln -f "${err_ini}" "${err_dsc}" || {
    echo "warning($(basename "${BASH_SOURCE[0]}")):" \
        "failed to link initial stderr log '${err_ini}' to '${err_dsc}'." >&2
}
ln -f "${out_ini}" "${out_dsc}" || {
    echo "warning($(basename "${BASH_SOURCE[0]}")):" \
        "failed to link initial stdout log '${out_ini}' to '${out_dsc}'." >&2
}

#  Sort the BAM file by QNAME
samtools sort -@ "${threads}" -n -o "${outfile}" "${infile}" ||
    {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to sort file: '${infile}'." >&2
        exit 1
    }

echo "Successfully sorted '${infile}' to '${outfile}'." >&2

#  Remove the initial SLURM stderr and stdout TXT outfiles
rm -f "${err_ini}" "${out_ini}" || {
    echo "warning($(basename "${BASH_SOURCE[0]}")):" \
        "failed to remove initial Slurm log file(s)." >&2
}
