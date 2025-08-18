#!/bin/bash

#  submit_calculate_scaling_factor.sh
#  KA


#  Run script in interactive mode (true) or command-line mode (false)
interactive=false

#  If true, run script in debug mode
debug=true


#  Set paths, values, arguments, etc. for interactive mode
function set_interactive() {
    #  Set base repository paths
    dir_bas="${HOME}/repos"
    dir_rep="${dir_bas}/protocol_chipseq_signal_norm"

    #  Set alignment parameters
    aligner="bowtie2"
    a_type="global"
    req_flg=true
    flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
    mapq=1
    str_det="${aligner}_${a_type}_flag-${flg}_mapq-${mapq}"

    #  Define data directories
    dir_dat="${dir_rep}/data"
    dir_pro="${dir_dat}/processed"
    dir_aln="${dir_pro}/align_fastqs"
    dir_det="${dir_aln}/${str_det}"
    dir_bam="${dir_det}/sc"

    #  Set output directories
    dir_sig="${dir_pro}/compute_signal/${str_det}"
    dir_out="${dir_sig}/tables"

    #  Set hardcoded argument assignments
    dir_scr="${HOME}/repos/protocol_chipseq_signal_norm/scripts"
    env_nam="env_protocol"
    threads=8
    mode="alpha"
    ser_mip="${dir_bam}/IP_WT_Q_Hmo1_7750.sc.bam,${dir_bam}/IP_WT_Q_Hmo1_7751.sc.bam"
    ser_min="${ser_mip//IP_/in_}"
    ser_sip="${ser_mip//sc/sp}"
    ser_sin="${ser_sip//IP_/in_}"
    tbl_met="${dir_dat}/raw/docs/measurements_siqchip.tsv"
    eqn="6nd"
    fil_out="${dir_out}/${mode}_test.tsv"
    rnd=24
    err_out="${dir_out}/logs"
    if [[ "${mode}" == "alpha" ]]; then
        nam_job="calc_sf_${mode}_${eqn}"
    else
        nam_job="calc_sf_${mode}"
    fi
}


#  Parse keyword arguments, assigning them to variables
#+ 
#+ Most of the argument inputs are not checked, as this is performed by
#+ 'execute_*.sh' and, to a certain extent, the script(s) submitted to SLURM
dir_scr="${HOME}/repos/protocol_chipseq_signal_norm/scripts"
env_nam="env_protocol"
threads=1
mode="alpha"
ser_mip=""
ser_min=""
ser_sip=""
ser_sin=""
tbl_met="$(dirname "${dir_scr}")/data/raw/docs/measurements_siqchip.tsv"
eqn="6nd"
fil_out=""
rnd=24
err_out=""
nam_job="calc_sf_${mode}_${eqn}"

#  Define the help message
show_help=$(cat << EOM

'$(basename "${0}")' takes the following keyword arguments:
  -ds, --dir_scr  Directory containing workflow scripts
  -en, --env_nam  Mamba environment to activate (default: '${env_nam}')
   -t, --threads  Number of threads to use (default: '${threads}')
   -m, --mode     Scaling factor mode: 'alpha' or 'spike' (default: '${mode}')
  -mp, --ser_mip  Comma-separated list of sample "main" IP BAM files
  -mn, --ser_min  Comma-separated list of sample "main" input BAM files
  -sp, --ser_sip  Comma-separated list of sample "spike-in" IP BAM files
                  (required if '--mode spike', ignored if not)
  -sn, --ser_sin  Comma-separated list of sample "spike-in" input BAM files
                  (required if '--mode spike', ignored if not)
  -tm, --tbl_met  Tab-delimited input file of siQ-ChIP metadata metrics
                  (required if '--mode alpha', ignored if not)
  -eq, --eqn      Alpha equation to compute: '5', '5nd', '6', '6nd' (required
                  if '--mode alpha', ignored if not; default: '${eqn}')
  -fo, --fil_out  Tab-delimited text output file in which scaling factors,
                  minimum input depth factors, and other values are written
  -fl, --flg_len  #TODO: Implement this
  -fd, --flg_dep  #TODO: Implement this
   -r, --rnd      Number of decimal places for rounding scaling factors and
                  minimum input depth factors (default: '${rnd}')
  -eo, --err_out  Directory for stderr and stdout output files
  -nj, --nam_job  Name of job (default: '${nam_job}')

Dependencies:
  - GNU AWK
  - GNU Parallel
  - Samtools
  - SLURM

Example:
\`\`\`
bash \${HOME}/repos/protocol_chipseq_signal_norm/scripts/submit_calculate_scaling_factor_alpha.sh
    -ds \${HOME}/repos/protocol_chipseq_signal_norm/scripts
    -en env_protocol
    -t 3
    -m spike
    -mp \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/align_reads/bowtie2_global_flag-2_mapq-1/sc/IP_WT_Q_Hmo1_7751.sc.bam
    -mn \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/align_reads/bowtie2_global_flag-2_mapq-1/sc/in_WT_Q_Hmo1_7751.sc.bam
    -sp \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/align_reads/bowtie2_global_flag-2_mapq-1/sp/IP_WT_Q_Hmo1_7751.sp.bam
    -sn \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/align_reads/bowtie2_global_flag-2_mapq-1/sp/in_WT_Q_Hmo1_7751.sp.bam
    -fo \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/compute_signal/bowtie2_global_flag-2_mapq-1/tables/test_spike.tsv
    -r 24
    -eo \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/compute_signal/bowtie2_global_flag-2_mapq-1/tables/logs
    -nj calc_sf_spike
         > \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/compute_signal/bowtie2_global_flag-2_mapq-1/tables/logs/calc_sf_spike_par.IP_WT_Q_Hmo1_7751.sc.stdout.txt
        2> \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/compute_signal/bowtie2_global_flag-2_mapq-1/tables/logs/calc_sf_spike_par.IP_WT_Q_Hmo1_7751.sc.stderr.txt
\`\`\`
EOM
)

#  Parse arguments
if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    echo ""
    if ! ${interactive:-false}; then
        if [[ -z "${1}" ]]; then exit 1; else exit 0; fi
    fi
fi

if ${interactive:-false}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -ds|--dir_scr) dir_scr="${2}"; shift 2 ;;
            -en|--env_nam) env_nam="${2}"; shift 2 ;;
             -t|--threads) threads="${2}"; shift 2 ;;
             -m|--mode)    mode="${2}";    shift 2 ;;
            -mp|--ser_mip) ser_mip="${2}"; shift 2 ;;
            -mn|--ser_min) ser_min="${2}"; shift 2 ;;
            -sp|--ser_sip) ser_sip="${2}"; shift 2 ;;
            -sn|--ser_sin) ser_sin="${2}"; shift 2 ;;
            -tm|--tbl_met) tbl_met="${2}"; shift 2 ;;
            -eq|--eqn)     eqn="${2}";     shift 2 ;;
            -fo|--fil_out) fil_out="${2}"; shift 2 ;;
             -r|--rnd)     rnd="${2}";     shift 2 ;;
            -eo|--err_out) err_out="${2}"; shift 2 ;;
            -nj|--nam_job) nam_job="${2}"; shift 2 ;;
            *)
                echo "## Unknown argument passed: ${1} ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                exit 1
                ;;
        esac
    done
fi

#  Validate 'mode'
case "${mode}" in
    alpha|spike) : ;;
    *)
        echo \
            "Error: Scaling factor computation mode ('--mode') was assigned" \
            "'${mode}' but must be 'alpha' or 'spike'." >&2
        if ! ${interactive:-false}; then exit 1; fi
esac

#  Assign and validate variables for scripts and functions
scr_met="${dir_scr}/parse_metadata_siq_chip.py"
scr_alf="${dir_scr}/calculate_scaling_factor_alpha.py"
scr_spk="${dir_scr}/calculate_scaling_factor_spike.py"
scr_sub="${dir_scr}/functions/submit.sh"
scr_clc="${dir_scr}/functions/calculate_scaling_factor.sh"

for fil in "${scr_met}" "${scr_alf}" "${scr_spk}" "${scr_sub}" "${scr_clc}"; do
    if [[ ! -f "${fil}" ]]; then
        echo "Error: Script not found: '${fil}'."
        if ! ${interactive:-false}; then exit 1; fi
    fi
done

#  Source necessary functions
# shellcheck disable=SC1090
for scr in "${scr_sub}" "${scr_clc}"; do
    if ! \
        source "${scr}"
    then
        echo \
            "Error: Failed to source '${scr}', which contains functions" \
            "necessary to run '${0}'." >&2
        if ! ${interactive:-false}; then exit 1; fi
    fi
done
unset scr

#  Debug argument assignments
# shellcheck disable=SC2046
if ${debug:-false}; then
    debug_var \
        "dir_scr=${dir_scr}" "env_nam=${env_nam}" "threads=${threads}" \
        "mode=${mode}"       "ser_mip=${ser_mip}" "ser_min=${ser_min}" \
        $(
            if [[ "${mode}" == "spike" ]]; then
                echo "ser_sip=${ser_sip}"; echo "ser_sin=${ser_sin}"
            fi
        ) \
        "tbl_met=${tbl_met}" "eqn=${eqn}"         "fil_out=${fil_out}" \
        "rnd=${rnd}"         "err_out=${err_out}" "nam_job=${nam_job}"
fi

#  Construct arrays from serialized strings
IFS=","
read -r -a arr_mip <<< "${ser_mip}"
read -r -a arr_min <<< "${ser_min}"

if [[ "${mode}" == "spike" ]]; then
    read -r -a arr_sip <<< "${ser_sip}"
    read -r -a arr_sin <<< "${ser_sin}"
fi
unset IFS

#  Debug output to check number of array elements and array element values
if ${debug:-false}; then
    echo "\${#arr_mip[@]}=${#arr_mip[@]}" && echo ""
    echo "arr_mip=( ${arr_mip[*]} )"      && echo ""
    echo "\${#arr_min[@]}=${#arr_min[@]}" && echo ""
    echo "arr_min=( ${arr_min[*]} )"      && echo ""
    
    if [[ "${mode}" == "spike" ]]; then
        echo "\${#arr_sip[@]}=${#arr_sip[@]}" && echo ""
        echo "arr_sip=( ${arr_sip[*]} )"      && echo ""
        echo "\${#arr_sin[@]}=${#arr_sin[@]}" && echo ""
        echo "arr_sin=( ${arr_sin[*]} )"      && echo ""
    fi
fi

#  Activate environment
activate_env "${env_nam}" || exit 1

#  Determine and run mode: SLURM or GNU Parallel/serial
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    #  Mode: SLURM
    id_job=${SLURM_ARRAY_JOB_ID}
    id_tsk=${SLURM_ARRAY_TASK_ID}

    if [[ ! "${id_tsk}" =~ ^[0-9]+$ || "${id_tsk}" -lt 1 ]]; then
        echo "Error: SLURM task ID is invalid: '${id_tsk}'." >&2
        exit 1
    else
        idx=$(( id_tsk - 1 ))
    fi

    #  Debug short names of environmental variables
    if ${debug:-false}; then
        debug_var "id_job=${id_job}" "id_tsk=${id_tsk}" "idx=${idx}"
    fi

    #  Derive sample name, which is needed for 'set_logs'
    # shellcheck disable=SC2001
    {
        samp="${arr_mip[idx]}"
        samp="${samp##*/IP_}"
        samp=$(echo "${samp%.bam}" | sed 's:\.:_:g')
    }

    #  Debug output to verify sample name
    if ${debug:-false}; then debug_var "samp=${samp}"; fi

    #  Run function to set SLURM and symlinked log files
    IFS=';' read -r err_ini out_ini err_dsc out_dsc < <(
        set_logs_slurm \
            "${id_job}" "${id_tsk}" "${samp}" "${err_out}" "${nam_job}"
    ) || exit 1
    unset IFS

    if ${debug:-false}; then
        debug_var \
            "err_ini=${err_ini}" "out_ini=${out_ini}" \
            "err_dsc=${err_dsc}" "out_dsc=${out_dsc}"
    fi

    #  Run processing function for scaling factor
    case "${mode}" in
        alpha) process_samp_alpha "${idx}" ;;
        spike) process_samp_spike "${idx}" ;;
    esac

    #  Remove the initial SLURM stderr and stdout TXT outfiles
    rm "${err_ini}" "${out_ini}"
else
    #  Mode: GNU Parallel/serial
    for idx in "${!arr_mip[@]}"; do
        #  Run processing function for scaling factor
        case "${mode}" in
            alpha) process_samp_alpha "${idx}" ;;
            spike) process_samp_spike "${idx}" ;;
        esac
    done
fi
