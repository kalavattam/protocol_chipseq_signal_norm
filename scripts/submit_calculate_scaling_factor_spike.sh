#!/bin/bash

#  If true, run script in debug mode
debug=true

#  Parse keyword arguments, assigning them to variables; most of the argument
#+ inputs are not checked, as this is performed by execute_*.sh and to a
#+ certain extent by the scripts submitted to SLURM
threads=1
infiles=""
outfile=""
flg_in=false
flg_mc=false
err_out=""
nam_job="calc_sf_spike"
scr_mng=""
fnc_env=""
env_nam="env_analyze"
scr_spk=""

show_help=$(cat << EOM
$(basename "${0}") takes the following keyword arguments:
   -t, --threads  Number of threads to use.
   -i, --infiles  Comma-separated serialized string of sample IP S. cerevisiae
                  BAM infiles.
   -o, --outfile  Outfile of sample spike-in-derived scaling factors.
  -fi, --flg_in   Include input spike-in-derived scaling factors in outfile.
  -fm, --flg_mc   Include additional measurements, calculations in outfile.
  -eo, --err_out  Directory to store stderr and stdout outfiles.
  -nj, --nam_job  Name of job.
  -sm, --scr_mng  Conda package manager shell script, conda.sh.
  -fe, --fnc_env  handle_env.sh utility function script.
  -en, --env_nam  Name of Conda/Mamba environment to activate.
  -ss, --scr_spk  Script that calculates spike-in-derived scaling factors,
                  calculate_scaling_factor_spike.py.

All arguments are required. If not specified, --threads, --nam_job, and
--env_nam default to, respectively, threads=${threads}, nam_job=${nam_job}, and
env_nam=${env_nam}. Also, --flg_in and --flg_mc are flags.
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
         -o|--outfile) outfile="${2}"; shift 2 ;;
        -fi|--flg_in)  flg_in=true;    shift 1 ;;
        -fm|--flg_mc)  flg_mc=true;    shift 1 ;;
        -eo|--err_out) err_out="${2}"; shift 2 ;;
        -nj|--nam_job) nam_job="${2}"; shift 2 ;;
        -sm|--scr_mng) scr_mng="${2}"; shift 2 ;;
        -fe|--fnc_env) fnc_env="${2}"; shift 2 ;;
        -en|--env_nam) env_nam="${2}"; shift 2 ;;
        -ss|--scr_spk) scr_spk="${2}"; shift 2 ;;
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
    echo "outfile=${outfile}"
    echo ""
    echo "flg_in=${flg_in}"
    echo ""
    echo "flg_mc=${flg_mc}"
    echo ""
    echo "err_out=${err_out}"
    echo ""
    echo "nam_job=${nam_job}"
    echo ""
    echo "scr_mng=${scr_mng}"
    echo ""
    echo "fnc_env=${fnc_env}"
    echo ""
    echo "env_nam=${env_nam}"
    echo ""
    echo "scr_spk=${scr_spk}"
    echo ""
fi

#  Activate environment
# shellcheck disable=SC1090,SC1091
if [[ "${CONDA_DEFAULT_ENV}" != "${env_nam}" ]]; then
    source "${scr_mng}"
    source "${fnc_env}"
    handle_env "${env_nam}"
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

#  Define main and spike-in IP and input infile assignments based on
#+ SLURM_ARRAY_TASK_ID
# shellcheck disable=SC2001
{    
    mp="${arr_infiles[$(( id_tsk - 1 ))]}"
    sp=$(echo "${mp}" | sed 's:\/sc\/:\/sp\/:g; s:\.sc\.:\.sp\.:g')
    mn=$(echo "${mp}" | sed 's:\/IP_:\/in_:g')
    sn=$(echo "${sp}" | sed 's:\/IP_:\/in_:g')
}

#  Debug assignments to variables mp, sp, mn, and mn
if ${debug}; then
    echo "mp=${mp}"
    echo ""
    echo "sp=${sp}"
    echo ""
    echo "mn=${mn}"
    echo ""
    echo "sn=${sn}"
    echo ""
fi

#  Exit if variables mp, sp, mn, or sn have empty assignments
if [[ -z "${mp}" ]]; then
    echo \
        "Error: Failed to retrieve mp for id_tsk=${id_tsk}:" \
        "\${arr_infiles[$(( id_tsk - 1 ))]}." >&2
    exit 1
elif [[ -z "${sp}" ]]; then
    echo \
        "Error: Failed to retrieve sp for id_tsk=${id_tsk}:" \
        "\$(echo \${arr_infiles[$(( id_tsk - 1 ))]} | sed s:\/sc\/:\/sp\/:g;" \
        "s:\.sc\.:\.sp\.:g')." >&2
    exit 1
elif [[ -z "${mn}" ]]; then
    echo \
        "Error: Failed to retrieve mn for id_tsk=${id_tsk}:" \
        "\$(echo \${arr_infiles[$(( id_tsk - 1 ))]} | sed" \
        "'s:\/IP_:\/in_:g')." >&2
elif [[ -z "${sn}" ]]; then
    echo \
        "Error: Failed to retrieve sn for id_tsk=${id_tsk}:" \
        "\$(echo \${arr_infiles[$(( id_tsk - 1 ))]} | sed 's:\/sc\/:\/sp\/:g;" \
        "s:\.sc\.:\.sp\.:g; s:\/IP_:\/in_:g')." >&2
fi

#  Exit if the files assigned to mp, sp, mn, or sn do not exist
if [[ ! -f "${mp}" ]]; then
    echo \
        "Error: File not found for mp from id_tsk=${id_tsk}:" \
        "\${arr_infiles[$(( id_tsk - 1 ))]}." >&2
    exit 1
elif [[ ! -f "${sp}" ]]; then
    echo \
        "Error: File not found for sp from id_tsk=${id_tsk}:" \
        "\$(echo \${arr_infiles[$(( id_tsk - 1 ))]} | sed s:\/sc\/:\/sp\/:g;" \
        "s:\.sc\.:\.sp\.:g')." >&2
    exit 1
elif [[ ! -f "${mn}" ]]; then
    echo \
        "Error: File not found for mn from id_tsk=${id_tsk}:" \
        "\$(echo \${arr_infiles[$(( id_tsk - 1 ))]} | sed" \
        "'s:\/IP_:\/in_:g')." >&2
elif [[ ! -f "${sn}" ]]; then
    echo \
        "Error: File not found for sn from id_tsk=${id_tsk}:" \
        "\$(echo \${arr_infiles[$(( id_tsk - 1 ))]} | sed 's:\/sc\/:\/sp\/:g;" \
        "s:\.sc\.:\.sp\.:g; s:\/IP_:\/in_:g')." >&2
fi

#  Check calls to samtools view -c for variables mp, sp, mn, and sn
if ${debug}; then
    echo "num_mp=\$(samtools view -c -@ ${threads} \"${mp}\")"
    echo "num_sp=\$(samtools view -c -@ ${threads} \"${sp}\")"
    echo "num_mn=\$(samtools view -c -@ ${threads} \"${mn}\")"
    echo "num_sn=\$(samtools view -c -@ ${threads} \"${sn}\")"
    echo ""
fi

#  Tally alignment numbers for mp, sp, mn, and sn
# shellcheck disable=SC2086
{
    num_mp=$(samtools view -c -@ ${threads} "${mp}")
    num_sp=$(samtools view -c -@ ${threads} "${sp}")
    num_mn=$(samtools view -c -@ ${threads} "${mn}")
    num_sn=$(samtools view -c -@ ${threads} "${sn}")
}

#  Debug values assigned to the num_{m|s}{p|n} variables
if ${debug}; then
    echo "num_mp=${num_mp}"
    echo ""
    echo "num_sp=${num_sp}"
    echo ""
    echo "num_mn=${num_mn}"
    echo ""
    echo "num_sn=${num_sn}"
    echo ""
fi

#  Derive sample name from mp assignment
# shellcheck disable=SC2001
{    
    samp="${mp##*/IP_}"
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

#  Calculate the spike-in derived scaling factor, which is assigned to variable
#+ sf, with script calculate_scaling_factor_spike.py
# shellcheck disable=SC2154
sf=$(
    python "${scr_spk}" \
        --main_ip  "${num_mp}" \
        --spike_ip "${num_sp}" \
        --main_in  "${num_mn}" \
        --spike_in "${num_sn}"
)

#  Debug output to verify spike-in-derived scaling factor
if ${debug}; then
    echo "sf=${sf}"
    echo ""
fi

#  Write header with flock
{
    flock -n 200 || exit 1  # Only proceed if lock is acquired

    #  Check if the header is already present in the file
    if ! grep -q "^sample"$'\t'"sf" "${outfile}" 2> /dev/null; then
        #  If the header is not present, write it
        if ! ${flg_mc}; then
            echo "sample"$'\t'"sf" >> "${outfile}"
        else
            echo "sample"$'\t'"sf"$'\t'"main_ip"$'\t'"spike_ip"$'\t'"main_in"$'\t'"spike_in" \
                >> "${outfile}"
        fi
    fi
} 200> "${outfile}.lock"  # Use a lock file

#  Print the IP sample and scaling factor to the outfile
if ! ${flg_mc}; then
    echo "${mp##*/}"$'\t'"${sf}" >> "${outfile}"

    #  If --flg_in, print the input sample and scaling factor (1) to the
    #+ outfile
    if ${flg_in}; then
        echo "${mn##*/}"$'\t'"1" >> "${outfile}"
    fi
else
    #  If --flg_mc, include num_{m|s}{p|n} values in the outfile
    echo "${mp##*/}"$'\t'"${sf}"$'\t'"${num_mp}"$'\t'"${num_sp}"$'\t'"${num_mn}"$'\t'"${num_sn}" \
        >> "${outfile}"

    #  If --flg_in, print the input sample and scaling factor (1), as well as
    #+ num_{m|s}{p|n} values, to the outfile
    if ${flg_in}; then
        echo "${mn##*/}"$'\t'"1"$'\t'"${num_mp}"$'\t'"${num_sp}"$'\t'"${num_mn}"$'\t'"${num_sn}" \
            >> "${outfile}"
    fi
fi

#  Remove the initial SLURM stderr and stdout TXT outfiles
rm "${err_ini}" "${out_ini}"
