#!/bin/bash

#  If true, run script in debug mode
debug=true

#  Parse keyword arguments, assigning them to variables; most of the argument
#+ inputs are not checked, as this is performed by execute_*.sh and to a
#+ certain extent by the scripts submitted to SLURM
threads=1
infiles=""
table=""
outfile=""
flg_dep=false
flg_len=false
flg_in=false
flg_mc=false
err_out=""
nam_job="calc_sf_alpha"
env_nam="env_analyze"
scr_par=""  #TODO
scr_alf=""  # calculate_scaling_factor_alpha.py

show_help=$(cat << EOM
$(basename "${0}") takes the following keyword arguments:
   -t, --threads  Number of threads to use.
   -i, --infiles  Comma-separated serialized string of sample IP BAM infiles.
  -tb, --table    Comma- or tab-separated table of siQ-ChIP metrics.
   -o, --outfile  Outfile of sample siQ-ChIP alpha values.
  -fd, --flg_dep  Calculate the number of alignments.
  -fl, --flg_len  Calculate the mean fragment length.
  -fi, --flg_in   Include sample input alpha values in outfile.
  -fm, --flg_mc   Include additional measurements, calculations in outfile.
  -eo, --err_out  Directory to store stderr and stdout outfiles.
  -nj, --nam_job  Name of job.
  -en, --env_nam  Name of Mamba environment to activate.
  -sp, --scr_par  Path to script that parses siQ-ChIP metadata.
  -sa, --scr_alf  Path to script that calculates siQ-ChIP alpha values.

All arguments are required. If not specified, --threads, --nam_job, and
--env_nam default to, respectively, threads=${threads}, nam_job=${nam_job}, and
env_nam=${env_nam}. Also, --flg_dep, --flg_len, --flg_in, and --flg_mc are
flags.
EOM
)

if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    exit 0
fi

# shellcheck disable=SC2034,SC2154
while [[ "$#" -gt 0 ]]; do
    case "${1}" in
         -t|--threads) threads="${2}"; shift 2 ;;
         -i|--infiles) infiles="${2}"; shift 2 ;;
        -tb|--table)   table="${2}";   shift 2 ;;
         -o|--outfile) outfile="${2}"; shift 2 ;;
        -fd|--flg_dep) flg_dep=true;   shift 1 ;;
        -fl|--flg_len) flg_len=true;   shift 1 ;;
        -fi|--flg_in)  flg_in=true;    shift 1 ;;
        -fm|--flg_mc)  flg_mc=true;    shift 1 ;;
        -eo|--err_out) err_out="${2}"; shift 2 ;;
        -nj|--nam_job) nam_job="${2}"; shift 2 ;;
        -en|--env_nam) env_nam="${2}"; shift 2 ;;
        -sp|--scr_par) scr_par="${2}"; shift 2 ;;
        -sa|--scr_alf) scr_alf="${2}"; shift 2 ;;
        *)
            echo "## Unknown argument passed: ${1} ##" >&2
            echo "" >&2
            echo "${show_help}" >&2
            exit 1
            ;;
    esac
done

#  Debug output to check argument assignments
if ${debug}; then
    echo "threads=${threads}"
    echo ""
    echo "infiles=${infiles}"
    echo ""
    echo "table=${table}"
    echo ""
    echo "outfile=${outfile}"
    echo ""
    echo "flg_dep=${flg_dep}"
    echo ""
    echo "flg_len=${flg_len}"
    echo ""
    echo "flg_in=${flg_in}"
    echo ""
    echo "flg_mc=${flg_in}"
    echo ""
    echo "err_out=${err_out}"
    echo ""
    echo "nam_job=${nam_job}"
    echo ""
    echo "env_nam=${env_nam}"
    echo ""
    echo "scr_par=${scr_par}"
    echo ""
    echo "scr_alf=${scr_alf}"
    echo ""
fi

#  Activate environment
if [[ "${CONDA_DEFAULT_ENV}" != "${env_nam}" ]]; then
    # eval "$(conda shell.bash hook)"
    # shellcheck disable=SC1091
    source "${HOME}/miniconda3/etc/profile.d/conda.sh"
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

#  Reconstruct arr_infiles from the serialized string assigned to infiles
IFS=',' read -r -a arr_infiles <<< "${infiles}"

#  Debug output to check the task ID/array index, number of array elements, and
#+ array element values
if ${debug}; then
    echo "SLURM_ARRAY_JOB_ID=${id_job}"
    echo ""
    echo "SLURM_ARRAY_TASK_ID=${id_tsk}"
    echo ""
    echo "\${#arr_infiles[@]}=${#arr_infiles[@]}"
    echo ""
    echo "arr_infiles=( ${arr_infiles[*]} )"
    echo ""
fi

#  Define IP and input infile assignments based on SLURM_ARRAY_TASK_ID
# shellcheck disable=SC2001
{    
    file_ip="${arr_infiles[$(( id_tsk - 1 ))]}"
    file_in="$(echo "${file_ip}" | sed 's:\/IP_:\/in_:g')"
}

#  Debug output to check which IP and input infiles are being processed
if ${debug}; then
    echo "file_ip=${file_ip}"
    echo ""
    echo "file_in=${file_in}"
    echo ""
fi

#  Exit if variables file_ip or file_in have empty assignments
if [[ -z "${file_ip}" ]]; then
    echo \
        "Error: Failed to retrieve file_ip for id_tsk=${id_tsk}:" \
        "\${arr_infiles[$(( id_tsk - 1 ))]}." >&2
    exit 1
elif [[ -z "${file_in}" ]]; then
    echo \
        "Error: Failed to retrieve file_in for id_tsk=${id_tsk}:" \
        "\$(echo \${arr_infiles[$(( id_tsk - 1 ))]} | sed" \
        "'s:\/IP_:\/in_:g')." >&2
    exit 1
fi

#  Exit if files assigned to variables file_ip or file_in do not exist
if [[ ! -f "${file_ip}" ]]; then
    echo \
        "Error: File not found for file_ip from id_tsk=${id_tsk}:" \
        "\${arr_infiles[$(( id_tsk - 1 ))]}." >&2
    exit 1
elif [[ ! -f "${file_in}" ]]; then
    echo \
        "Error: File not found for file_in from id_tsk=${id_tsk}:" \
        "\$(echo \${arr_infiles[$(( id_tsk - 1 ))]} | sed" \
        "'s:\/IP_:\/in_:g')." >&2
    exit 1
fi

#  Check call to sourced Python script
if ${debug}; then
    echo "source <("
    echo "    python \"${scr_par}\" \\"
    echo "        --text \"${table}\" \\"
    echo "        --bam \"${file_ip}\" \\"
    echo "        --shell"
    echo ")"
    echo ""
fi

#  Assign variables for siQ-ChIP measurements: volume, mass, concentration,
#+ length
# shellcheck disable=SC1090
source <(
    python "${scr_par}" \
        --text "${table}" \
        --bam "${file_ip}" \
        --shell
)

#  If --flg_dep, calculate sequencing depth (number of alignments) for IP and
#+ input samples
if ${flg_dep}; then
    depth_ip=$(samtools view -@ "${threads}" -c "${file_ip}")
    depth_in=$(samtools view -@ "${threads}" -c "${file_in}")
fi

#  If --flg_len, calculate mean fragment length for IP and input samples
if ${flg_len}; then
    length_ip="$(
        samtools view -@ "${threads}" "${file_ip}" \
            | awk '{
                if ($9 > 0) { sum += $9; count++ }
            } END {
                if (count > 0) { print sum / count }
            }'
    )"
    length_in="$(
        samtools view -@ "${threads}" "${file_in}" \
            | awk '{
                if ($9 > 0) { sum += $9; count++ }
            } END {
                if (count > 0) { print sum / count }
            }'
    )"
fi

#  Debug output to check IP and input volume, mass, concentration, and flg_len
#+ values
# shellcheck disable=SC2154
if ${debug}; then
    echo "mass_ip=${mass_ip}"
    echo ""
    echo "mass_in=${mass_in}"
    echo ""
    echo "volume_ip=${volume_ip}"
    echo ""
    echo "volume_in=${volume_in}"
    echo ""
    echo "depth_ip=${depth_ip}"
    echo ""
    echo "depth_in=${depth_in}"
    echo ""
    echo "length_ip=${length_ip}"
    echo ""
    echo "length_in=${length_in}"
    echo ""
fi

#  Derive sample name from file_ip assignment
# shellcheck disable=SC2001
{    
    samp="${file_ip##*/IP_}"
    samp=$(echo "${samp%.bam}" | sed 's:\.:_:g')
}

#  Debug output to verify sample name
if ${debug}; then
    echo "samp=${samp}"
    echo ""
fi

#  Assign stdout and stderr outstems to variables
err_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stderr.txt"
out_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stdout.txt"
err_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stderr.txt"
out_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stdout.txt"

#  Give the initial stderr and stdout TXT outfiles more descriptive names
ln -f "${err_ini}" "${err_dsc}"
ln -f "${out_ini}" "${out_dsc}"

#  Check call to calculate_scaling_factor_alpha.py
if ${debug}; then
    echo "python \"${scr_alf}\" \\"
    echo "    --mass_ip ${mass_ip} \\"
    echo "    --mass_in ${mass_in} \\"
    echo "    --volume_ip ${volume_ip} \\"
    echo "    --volume_in ${volume_in} \\"
    echo "    --depth_ip ${depth_ip} \\"
    echo "    --depth_in ${depth_in} \\"
    echo "    --length_ip ${length_ip} \\"
    echo "    --length_in ${length_in}"
    echo ""
fi

#  Run calculate_scaling_factor_alpha.py
# shellcheck disable=SC2046,2086
alpha=$(
    python "${scr_alf}" \
        --mass_ip ${mass_ip} \
        --mass_in ${mass_in} \
        --volume_ip ${volume_ip} \
        --volume_in ${volume_in} \
        --depth_ip ${depth_ip} \
        --depth_in ${depth_in} \
        --length_ip ${length_ip} \
        --length_in ${length_in}
)

#  Debug output to verify siQ-ChIP alpha value
if ${debug}; then
    echo "alpha=${alpha}"
    echo ""
fi

#  Write header with flock
{
    flock -n 200 || exit 1  # Only proceed if lock is acquired

    #  Check if the header is already present in the file
    if ! grep -q "^sample"$'\t'"alpha" "${outfile}" 2> /dev/null; then
        #  If the header is not present, write it
        if ! ${flg_mc}; then
            echo "sample"$'\t'"alpha" >> "${outfile}"
        else
            echo "sample"$'\t'"alpha"$'\t'"mass_ip"$'\t'"mass_in"$'\t'"volume_ip"$'\t'"volume_in"$'\t'"depth_ip"$'\t'"depth_in"$'\t'"length_ip"$'\t'"length_in" \
                >> "${outfile}"
        fi
    fi
} 200> "${outfile}.lock"  # Use a lock file

#  Print the IP sample and alpha value to the outfile
if ! ${flg_mc}; then
    echo "${file_ip##*/}"$'\t'"${alpha}" >> "${outfile}"

    #  If --flg_in, print the input sample and alpha value (#N/A) to the
    #+ outfile
    if ${flg_in}; then
        echo "${file_in##*/}"$'\t'"#N/A" >> "${outfile}"
    fi
else
    #  If --flg_mc, include mass, volume, depth, and length values in the outfile
    echo "${file_ip##*/}"$'\t'"${alpha}"$'\t'"${mass_ip}"$'\t'"${mass_in}"$'\t'"${volume_ip}"$'\t'"${volume_in}"$'\t'"${depth_ip}"$'\t'"${depth_in}"$'\t'"${length_ip}"$'\t'"${length_in}" \
        >> "${outfile}"

    #  If --flg_in, print the input sample and alpha value (#N/A), as well as
    #+ mass, volume, depth, and length values, to the outfile
    if ${flg_in}; then
        echo "${file_in##*/}"$'\t'"#N/A"$'\t'"${mass_ip}"$'\t'"${mass_in}"$'\t'"${volume_ip}"$'\t'"${volume_in}"$'\t'"${depth_ip}"$'\t'"${depth_in}"$'\t'"${length_ip}"$'\t'"${length_in}" \
            >> "${outfile}"
    fi
fi

#  Remove the initial SLURM stderr and stdout TXT outfiles
rm "${err_ini}" "${out_ini}"
