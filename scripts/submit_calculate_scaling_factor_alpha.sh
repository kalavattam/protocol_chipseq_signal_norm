#!/bin/bash

#  submit_calculate_scaling_factor_alpha.sh
#  KA


#  If true, run script in debug mode
debug=true

#  Parse keyword arguments, assigning them to variables; most of the argument
#+ inputs are not checked, as this is performed by execute_*.sh and to a
#+ certain extent by the scripts submitted to SLURM
threads=1
ser_ip=""
ser_in=""
table=""
eqn="6nd"
outfile=""
flg_dep=false
flg_len=false
flg_mc=false
err_out=""
nam_job="calc_sf_alpha_${eqn}"
env_nam="env_analyze"
scr_met=""
scr_alf=""

show_help=$(cat << EOM
$(basename "${0}") takes the following keyword arguments:
   -t, --threads  Number of threads to use (default: ${threads}).
  -sp, --ser_ip   Comma-separated serialized string of sample IP BAM infiles.
  -sn, --ser_in   Comma-separated serialized string of corresponding sample
                  input BAM infiles.
  -tb, --table    Comma- or tab-separated table of siQ-ChIP metrics.
  -eq, --eqn      Alpha equation to compute; options: '5', '5nd', '6', '6nd'
                  (default: '${eqn}').
   -o, --outfile  Outfile of sample siQ-ChIP alpha values.
  -fd, --flg_dep  Calculate the number of alignments.
  -fl, --flg_len  Calculate the mean fragment length.
  -fm, --flg_mc   Include additional measurements, calculations in outfile.
  -eo, --err_out  Directory to store stderr and stdout outfiles.
  -nj, --nam_job  Name of job (default: '${nam_job}').
  -en, --env_nam  Name of Conda/Mamba environment to activate (default:
                  '${env_nam}').
  -sm, --scr_met  Script that parses siQ-ChIP metadata,
                  'parse_metadata_siq_chip.py'.
  -sa, --scr_alf  Script that calculates siQ-ChIP alpha values,
                  'calculate_scaling_factor_alpha.py'.

All arguments are required with the following notes and exceptions:
  - '--threads ${threads}'
  - '--eqn ${eqn}'
  - '--nam_job ${nam_job}'
  - '--env_nam ${env_nam}'
  - Also, '--flg_dep', '--flg_len', and '--flg_mc' are flags.
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
        -sp|--ser_ip)  ser_ip="${2}";  shift 2 ;;
        -sn|--ser_in)  ser_in="${2}";  shift 2 ;;
        -tb|--table)   table="${2}";   shift 2 ;;
        -eq|--eqn)     eqn="${2}";     shift 2 ;;
         -o|--outfile) outfile="${2}"; shift 2 ;;
        -fd|--flg_dep) flg_dep=true;   shift 1 ;;
        -fl|--flg_len) flg_len=true;   shift 1 ;;
        -fm|--flg_mc)  flg_mc=true;    shift 1 ;;
        -eo|--err_out) err_out="${2}"; shift 2 ;;
        -nj|--nam_job) nam_job="${2}"; shift 2 ;;
        -en|--env_nam) env_nam="${2}"; shift 2 ;;
        -sm|--scr_met) scr_met="${2}"; shift 2 ;;
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
    echo "ser_ip=${ser_ip}"
    echo ""
    echo "ser_in=${ser_in}"
    echo ""
    echo "table=${table}"
    echo ""
    echo "eqn=${eqn}"
    echo ""
    echo "outfile=${outfile}"
    echo ""
    echo "flg_dep=${flg_dep}"
    echo ""
    echo "flg_len=${flg_len}"
    echo ""
    echo "flg_mc=${flg_mc}"
    echo ""
    echo "err_out=${err_out}"
    echo ""
    echo "nam_job=${nam_job}"
    echo ""
    echo "env_nam=${env_nam}"
    echo ""
    echo "scr_met=${scr_met}"
    echo ""
    echo "scr_alf=${scr_alf}"
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

#  Construct 'arr_ip' and 'arr_in' from the serialized strings assigned to,
#+ respectively, 'ser_ip' and 'ser_in'
IFS=',' read -r -a arr_ip <<< "${ser_ip}"
IFS=',' read -r -a arr_in <<< "${ser_in}"

#  Debug output to check the task ID/array index, number of array elements, and
#+ array element values
if ${debug}; then
    echo "SLURM_ARRAY_JOB_ID=${id_job}"
    echo ""
    echo "SLURM_ARRAY_TASK_ID=${id_tsk}"
    echo ""
    echo "\${#arr_ip[@]}=${#arr_ip[@]}"
    echo ""
    echo "arr_ip=( ${arr_ip[*]} )"
    echo ""
    echo "\${#arr_in[@]}=${#arr_in[@]}"
    echo ""
    echo "arr_in=( ${arr_in[*]} )"
    echo ""
fi

#  Determine array index based on SLURM_ARRAY_TASK_ID
idx=$(( id_tsk - 1 ))

#  Define IP and input infile assignments based on index of reconstructed array
fil_ip="${arr_ip[idx]}"
fil_in="${arr_in[idx]}"

#  Debug output, checking which IP and input infiles are in use
if ${debug}; then
    echo "fil_ip=${fil_ip}"
    echo ""
    echo "fil_in=${fil_in}"
    echo ""
fi

#  Exit if variables 'fil_ip' or 'fil_in' have empty assignments
if [[ -z "${fil_ip}" ]]; then
    echo \
        "Error: Failed to retrieve 'fil_ip' for 'id_tsk=${id_tsk}':" \
        "'\${arr_ip[${idx}]}'." >&2
    exit 1
elif [[ -z "${fil_in}" ]]; then
    echo \
        "Error: Failed to retrieve 'fil_in' for 'id_tsk=${id_tsk}':" \
        "'\${arr_in[${idx}]}'." >&2
    exit 1
fi

#  Exit if files assigned to variables 'fil_ip' or 'fil_in' do not exist
if [[ ! -f "${fil_ip}" ]]; then
    echo \
        "Error: File not found for 'fil_ip' from 'id_tsk=${id_tsk}':" \
        "'\${arr_ip[${idx}]}'." >&2
    exit 1
elif [[ ! -f "${fil_in}" ]]; then
    echo \
        "Error: File not found for 'fil_in' from 'id_tsk=${id_tsk}':" \
        "'\${arr_in[${idx}]}'." >&2
    exit 1
fi

#  Debug call to sourced Python script
if ${debug}; then
    echo "source <("
    echo "    python \"${scr_met}\" \\"
    echo "        --text \"${table}\" \\"
    echo "        --eqn \"${eqn}\" \\"
    echo "        --bam \"${fil_ip}\" \\"
    echo "        --shell"
    echo ")"
    echo ""
fi

#  Assign variables for siQ-ChIP measurements: volume, mass, concentration,
#+ length
# shellcheck disable=SC1090
source <(
    python "${scr_met}" \
        --text "${table}" \
        --eqn "${eqn}" \
        --bam "${fil_ip}" \
        --shell
)

#  If '--flg_dep', calculate sequencing depth (number of alignments) for IP and
#+ input samples
if ${flg_dep}; then
    if ${debug}; then
        echo "{"
        echo "    dep_ip=\$(samtools view -@ ${threads} -c \"${fil_ip}\")"
        echo "    dep_in=\$(samtools view -@ ${threads} -c \"${fil_in}\")"
        echo "}"
        echo ""
    fi

    # shellcheck disable=SC2086
    {
        dep_ip=$(samtools view -@ ${threads} -c "${fil_ip}")
        dep_in=$(samtools view -@ ${threads} -c "${fil_in}")
    }
fi

#  If '--flg_len', calculate mean fragment length for IP and input samples
if ${flg_len}; then
    if ${debug}; then
        echo "{"
        echo "    len_ip=\"\$("
        echo "        samtools view -@ ${threads} \"${fil_ip}\" "
        echo "            | awk '{"
        echo "                if (\$9 > 0) { sum += \$9; count++ }"
        echo "            } END {"
        echo "                if (count > 0) { print sum / count }"
        echo "            }'"
        echo "    )\""
        echo "    len_in=\"\$("
        echo "        samtools view -@ ${threads} \"${fil_in}\" "
        echo "            | awk '{"
        echo "                if (\$9 > 0) { sum += \$9; count++ }"
        echo "            } END {"
        echo "                if (count > 0) { print sum / count }"
        echo "            }'"
        echo "    )\""
        echo "}"
        echo ""
    fi

    # shellcheck disable=SC2086
    {
        len_ip="$(
            samtools view -@ ${threads} "${fil_ip}" \
                | awk '{
                    if ($9 > 0) { sum += $9; count++ }
                } END {
                    if (count > 0) { print sum / count }
                }'
        )"
        len_in="$(
            samtools view -@ ${threads} "${fil_in}" \
                | awk '{
                    if ($9 > 0) { sum += $9; count++ }
                } END {
                    if (count > 0) { print sum / count }
                }'
        )"
    }
fi

#  Debug output, checking IP and input volume, mass, concentration, and
#+ 'flg_len' values
# shellcheck disable=SC2154
if ${debug}; then
    echo "mass_ip=${mass_ip}"
    echo ""
    echo "mass_in=${mass_in}"
    echo ""
    echo "vol_all=${vol_all}"
    echo ""
    echo "vol_in=${vol_in}"
    echo ""
    echo "dep_ip=${dep_ip}"
    echo ""
    echo "dep_in=${dep_in}"
    echo ""
    echo "len_ip=${len_ip}"
    echo ""
    echo "len_in=${len_in}"
    echo ""
fi

#  Derive sample name from fil_ip assignment
# shellcheck disable=SC2001
{    
    samp="${fil_ip##*/IP_}"
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
    echo "    --eqn \"${eqn}\" \\"
    echo "    --mass_ip ${mass_ip} \\"
    echo "    --mass_in ${mass_in} \\"
    echo "    --vol_all ${vol_all} \\"
    echo "    --vol_in ${vol_in} \\"
    echo "    --dep_ip ${dep_ip} \\"
    echo "    --dep_in ${dep_in} \\"
    echo "    --len_ip ${len_ip} \\"
    echo "    --len_in ${len_in}"
    echo ""
fi

#  Run calculate_scaling_factor_alpha.py
# shellcheck disable=SC2046,2086
alpha=$(
    python "${scr_alf}" \
        --eqn "${eqn}" \
        --mass_ip ${mass_ip} \
        --mass_in ${mass_in} \
        --vol_all ${vol_all} \
        --vol_in ${vol_in} \
        --dep_ip ${dep_ip} \
        --dep_in ${dep_in} \
        --len_ip ${len_ip} \
        --len_in ${len_in}
)

#  Debug output to verify siQ-ChIP alpha value
if ${debug}; then
    echo "alpha=${alpha}"
    echo ""
fi

#  Print the IP sample, input sample, and alpha value to the outfile
if ${flg_mc}; then
    #  If '--flg_mc', include mass, volume, depth, and length values in the
    #+ outfile
    echo -e \
        "${fil_ip}\t${fil_in}\t${alpha}\t${eqn}\t${mass_ip}\t${mass_in}\t${vol_all}\t${vol_in}\t${dep_ip}\t${dep_in}\t${len_ip}\t${len_in}" \
            >> "${outfile}"
else
    echo -e "${fil_ip}\t${fil_in}\t${alpha}\t${eqn}" >> "${outfile}"
fi

#  Remove the initial SLURM stderr and stdout TXT outfiles
rm "${err_ini}" "${out_ini}"
