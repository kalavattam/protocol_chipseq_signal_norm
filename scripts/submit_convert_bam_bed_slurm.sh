#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: submit_convert_bam_bed_slurm.sh
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
\${3}=csv_infile  <str>  Comma-separated list of QNAME-sorted BAM infiles (assumed PE for AWK processing code; use Python if data are SE).
\${4}=pth_scr     <str>  Path to 'compute_signal.py' or compatible BAM-to-BED script.
\${5}=dir_out     <str>  Directory to save BED outfiles.
\${6}=err_out     <str>  Directory for stdout and stderr files.
\${7}=nam_job     <str>  Name of job.
\${8}=use_awk     <bol>  Run the AWK processing code rather than the Python script (do not use with SE data). Boolean-like strings accepted: 'true' or 'false'.

Notes:
- The AWK and Python branches are not equivalent.
  + The AWK branch assumes QNAME-sorted PE records in adjacent pairs and writes paired-fragment intervals.
  + The Python branch uses 'compute_signal.py', which does not require QNAME sorting, handles SE, and emits fragments according to its own 'parse_bam()' policy.
- For tested PE fixtures, outputs are the same when the same data are passed to the AWK branch in QNAME-sorted format and to the Python branch in QNAME- or coordinate-sorted format.
EOM
)

if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
    cat >&2 << EOM
'$(basename "${0}")' requires 8 positional arguments:
${show_help}
EOM
    exit 0
fi

#  Check for exactly 8 positional arguments
if [[ $# -ne 8 ]]; then
    msg="but $# were supplied."
    [[ $# -eq 1 ]] && msg="but only $# was supplied."
    cat >&2 << EOM
error: '$(basename "${0}")' requires 8 positional arguments, ${msg}

The necessary positional arguments:
${show_help}
EOM
    exit 1
fi

#  Parse positional arguments, assigning them to variables
env_nam="${1}"
threads="${2}"
csv_infile="${3}"
pth_scr="${4}"
dir_out="${5}"
err_out="${6}"
nam_job="${7}"
use_awk="${8}"

#  Debug positional argument assignments
if [[ "${debug:-false}" == "true" ]]; then
    printf '%s\n\n' \
        "\${env_nam}=${env_nam}" \
        "\${threads}=${threads}" \
        "\${csv_infile}=${csv_infile}" \
        "\${pth_scr}=${pth_scr}" \
        "\${dir_out}=${dir_out}" \
        "\${err_out}=${err_out}" \
        "\${nam_job}=${nam_job}" \
        "\${use_awk}=${use_awk}" \
        >&2
fi

#  Validate number of threads
if ! [[ "${threads}" =~ ^[1-9][0-9]*$ ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "'threads' argument must be a positive integer: '${threads}'." >&2
    exit 1
fi

#  Validate AWK/Python selector
case "${use_awk}" in
    true|false) : ;;
    *)
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "'use_awk' argument must be 'true' or 'false': '${use_awk}'." >&2
        exit 1
        ;;
esac

#  Validate Python script for BAM-to-BED conversion
if [[ "${use_awk}" == "false" && ! -f "${pth_scr}" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "file does not exist: '${pth_scr}'." >&2
    exit 1
fi

#  Validate specified directories
for dir in "$(dirname "${csv_infile%%,*}")" "${dir_out}" "${err_out}"; do
    if [[ ! -d "${dir}" ]]; then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "directory does not exist: '${dir}'." >&2
        exit 1
    fi
done
unset dir

#  Activate environment
if [[ "${CONDA_DEFAULT_ENV:-}" != "${env_nam}" ]]; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${env_nam}" ||
        {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to activate environment: '${env_nam}'." >&2
            exit 1
        }
fi

#  Check for necessary programs
if [[ "${use_awk}" == "true" ]]; then
    for pgrm in samtools awk sort gzip; do
        if ! \
            command -v "${pgrm}" > /dev/null 2>&1
        then
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "required program is not in PATH: '${pgrm}'." >&2
            exit 1
        fi
    done
else
    if ! \
        command -v python > /dev/null 2>&1
    then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "required program is not in PATH: 'python'." >&2
        exit 1
    fi
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
            "file does not exist: '${file}'." >&2
        exit 1
    fi
done

#  Validates IDs
if ! [[ "${id_tsk}" =~ ^[1-9][0-9]*$ ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "'SLURM_ARRAY_TASK_ID' must be a positive integer: '${id_tsk}'." >&2
    exit 1
elif [[ "${id_tsk}" -gt "${#arr_infiles[@]}" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "'SLURM_ARRAY_TASK_ID' is out of bounds: '${id_tsk}'." >&2
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
outfile="${dir_out}/$(basename "${infile}" ".bam").bed.gz"

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
samp="$(basename "${outfile}")"

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

if [[ "${use_awk}" == "true" ]]; then
    #  Generate BED outfiles from QNAME-sorted BAM infiles using 'awk' (for
    #+ field formatting and fragment computation), 'sort' (to sort by
    #+ chromosome, column 1, and numerically by start position, column 2), and
    #+ 'gzip' (for compression)
    if ! \
        samtools view -@ "${threads}" "${infile}" \
            | awk '{
                #  Read paired-end BAM lines
                if (NR % 2 == 1) {
                    #  When NR is odd, process the first read in a pair
                    c1 = $3
                    s1 = $4
                    l1 = length($10)
                } else {
                    #  When NR is even, process the second read in a pair
                    c2 = $3
                    s2 = $4
                    l2 = length($10)

                    #  Process a read pair only if both are aligned to the same
                    #+ chromosome
                    if (c1 == c2) {
                        #  Compute fragment start, end, and length:
                        #+ - st: the smaller of the two read start positions
                        #+ - en: the larger of the two read end positions
                        #+ - lf: the difference between end and start plus one
                        st = (s1 < s2) ? s1 : s2
                        en = (s1 < s2) ? s2 + l2 - 1 : s1 + l1 - 1
                        lf = en - st + 1

                        #  Output fragment in BED format:
                        #+     chr, start, end, length
                        print c1, st, en, lf
                    }
                }
            }' OFS='\t' \
            | sort -k1,1 -k2,2n \
            | gzip \
                > "${outfile}"
    then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "BED file conversion failed for BAM infile: '${infile}'." >&2
        exit 1
    fi
else
    if ! \
        python "${pth_scr}" -i "${infile}" -o "${outfile}"
    then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "BED file conversion failed for BAM infile: '${infile}'." >&2
        exit 1
    fi
fi

echo "Successfully converted '${infile}' to '${outfile}'." >&2

#  Remove the initial SLURM stderr and stdout TXT outfiles
rm -f "${err_ini}" "${out_ini}" || {
    echo "warning($(basename "${BASH_SOURCE[0]}")):" \
        "failed to remove initial Slurm log file(s)." >&2
}
