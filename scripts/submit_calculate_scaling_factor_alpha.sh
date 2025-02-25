#!/bin/bash

#  submit_calculate_scaling_factor_alpha.sh
#  KA


#  If true, run script in debug mode
debug=true

#  Run script in interactive/test mode (true) or command-line mode (false)
interactive=false


#  Define functions
function debug_var() { for var in "$@"; do echo "${var}" && echo ""; done; }


function validate_var() {
    local var_nam=${1}
    local var_val=${2}
    local idx=${3}
    if [[ -z "${var_val}" ]]; then
        echo \
            "Error: '${var_nam}' is empty or unset for array index" \
            "'${idx}'." >&2
        return 1
    fi
}


function exists_var() {
    local var_nam=${1}
    local var_val=${2}
    local idx=${3}
    if [[ ! -f "${var_val}" ]]; then
        echo \
            "Error: '${var_nam}' does not exist for array index '${idx}'." >&2
        return 1
    fi
}


function set_logs() {
    local job_id=${1}
    local task_id=${2}
    local nam_smp=${3}
    local dir_log=${4}
    err_ini="${dir_log}/${nam_job}.${job_id}-${task_id}.stderr.txt"
    out_ini="${dir_log}/${nam_job}.${job_id}-${task_id}.stdout.txt"
    err_dsc="${dir_log}/${nam_job}.${nam_smp}.${job_id}-${task_id}.stderr.txt"
    out_dsc="${dir_log}/${nam_job}.${nam_smp}.${job_id}-${task_id}.stdout.txt"
}


#  Define subroutine to process a sample
function process_sample() {
    local idx="${1}"
    local fil_ip fil_in
    local dm_fr_1 dm_fr_10 dm_fr_20 dm_fr_30 dm_fr_40 dm_fr_50
    local dm_nm_1 dm_nm_10 dm_nm_20 dm_nm_30 dm_nm_40 dm_nm_50
    local prt_1 prt_2 prt_3 prt_4 prt_5 prt_6 prt_7

    #  Assign variables based on idx
    fil_ip="${arr_ip[idx]}"
    fil_in="${arr_in[idx]}"

    #  Debug and validate assignments to variables 'fil_ip' and 'fil_in'
    if ${debug}; then
        debug_var \
            "idx=${idx} ('idx' passed to 'process_sample()')" \
            "fil_ip=${fil_ip}" "fil_in=${fil_in}"
    fi

    #  Exit if variables 'fil_ip' and 'fil_in' have empty assignments or do not
    #+ exist
    # shellcheck disable=SC2086
    {
        validate_var "fil_ip" "${fil_ip}" ${idx} || exit 1
        exists_var   "fil_ip" "${fil_ip}" ${idx} || exit 1
        validate_var "fil_in" "${fil_in}" ${idx} || exit 1
        exists_var   "fil_in" "${fil_in}" ${idx} || exit 1
    }

    #  Debug call to sourced 'parse_metadata_siq_chip.py'
    if ${debug}; then
        echo "source <("
        echo "    python \"${scr_met}\" \\"
        echo "        --tbl_met \"${tbl_met}\" \\"
        echo "        --eqn \"${eqn}\" \\"
        echo "        --bam \"${fil_ip}\" \\"
        echo "        --shell"
        echo ")"
        echo ""
    fi

    #  Run 'parse_metadata_siq_chip.py' to assign variables for siQ-ChIP
    #+ metadata: volume, mass, concentration, length
    # shellcheck disable=SC1090
    source <(
        python "${scr_met}" \
            --tbl_met "${tbl_met}" \
            --eqn "${eqn}" \
            --bam "${fil_ip}" \
            --shell
    )

    #  Check and calculate fragment depth for IP and input samples
    if ${debug}; then
        echo "{"
        echo "    ## WARNING: Assumes BAM files contain paired-end alignments ##"
        echo "    dep_ip=\$(count_alignments_bam ${threads} \"${fil_ip}\")"
        echo "    dep_in=\$(count_alignments_bam ${threads} \"${fil_in}\")"
        echo "}"
        echo ""
    fi

    # shellcheck disable=SC2086
    {
        ## WARNING: Assumes BAM files contain paired-end alignments ##
        dep_ip=$(count_alignments_bam ${threads} "${fil_ip}")
        dep_in=$(count_alignments_bam ${threads} "${fil_in}")
    }

    #  Debug and calculate mean fragment lengths for IP and input samples
    if ${debug}; then
        echo "{"
        echo "    ## WARNING: Assumes BAM files contain paired-end alignments ##"
        echo "    ## WARNING: Overwrites length observations in 'tbl_met' ##"
        echo "    len_ip=\"\$(calculate_frag_avg ${threads} \"${fil_ip}\")"
        echo "    len_in=\"\$(calculate_frag_avg ${threads} \"${fil_in}\")"
        echo "}"
        echo ""
    fi

    # shellcheck disable=SC2086
    {
        ## WARNING: Assumes BAM files contain paired-end alignments ##
        ## WARNING: Overwrites length observations in 'tbl_met' ##
        len_ip=$(calculate_frag_avg ${threads} "${fil_ip}")
        len_in=$(calculate_frag_avg ${threads} "${fil_in}")
    }

    #  Debug output, checking IP and input volume, mass, concentration, and
    #+ length values
    # shellcheck disable=SC2154
    if ${debug}; then
        debug_var \
            "mass_ip=${mass_ip}" "mass_in=${mass_in}" \
            "vol_all=${vol_all}" "vol_in=${vol_in}" \
            "dep_ip=${dep_ip}"   "dep_in=${dep_in}" \
            "len_ip=${len_ip}"   "len_in=${len_in}"
    fi

    #  Check call to 'calculate_scaling_factor_alpha.py'
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
        echo "    --len_in ${len_in} \\"
        echo "    --rnd ${rnd} \\"
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
            --len_in ${len_in} \
            --rnd ${rnd}
    )

    #  Debug output to verify siQ-ChIP alpha value
    if ${debug}; then debug_var "alpha=${alpha}"; fi

    if ${debug}; then
        echo "{"
        echo "     dm_fr_1=\$(calculate_factor_depth ${dep_in} 1  12157105 \"frag\" ${rnd})"
        echo "     dm_fr_5=\$(calculate_factor_depth ${dep_in} 5  12157105 \"frag\" ${rnd})"
        echo "    dm_fr_10=\$(calculate_factor_depth ${dep_in} 10 12157105 \"frag\" ${rnd})"
        echo "    dm_fr_20=\$(calculate_factor_depth ${dep_in} 20 12157105 \"frag\" ${rnd})"
        echo "    dm_fr_30=\$(calculate_factor_depth ${dep_in} 30 12157105 \"frag\" ${rnd})"
        echo "    dm_fr_40=\$(calculate_factor_depth ${dep_in} 40 12157105 \"frag\" ${rnd})"
        echo "    dm_fr_50=\$(calculate_factor_depth ${dep_in} 50 12157105 \"frag\" ${rnd})"
        echo ""
        echo "     dm_nm_1=\$(calculate_factor_depth ${dep_in} 1  12157105 \"norm\" ${rnd})"
        echo "     dm_nm_5=\$(calculate_factor_depth ${dep_in} 5  12157105 \"norm\" ${rnd})"
        echo "    dm_nm_10=\$(calculate_factor_depth ${dep_in} 10 12157105 \"norm\" ${rnd})"
        echo "    dm_nm_20=\$(calculate_factor_depth ${dep_in} 20 12157105 \"norm\" ${rnd})"
        echo "    dm_nm_30=\$(calculate_factor_depth ${dep_in} 30 12157105 \"norm\" ${rnd})"
        echo "    dm_nm_40=\$(calculate_factor_depth ${dep_in} 40 12157105 \"norm\" ${rnd})"
        echo "    dm_nm_50=\$(calculate_factor_depth ${dep_in} 50 12157105 \"norm\" ${rnd})"
        echo "}"
        echo ""
    fi

    #  Calculate input minimum depth values for common bin sizes
    # shellcheck disable=SC2086
    {
         dm_fr_1=$(calculate_factor_depth ${dep_in} 1  12157105 "frag" ${rnd})
         dm_fr_5=$(calculate_factor_depth ${dep_in} 5  12157105 "frag" ${rnd})
        dm_fr_10=$(calculate_factor_depth ${dep_in} 10 12157105 "frag" ${rnd})
        dm_fr_20=$(calculate_factor_depth ${dep_in} 20 12157105 "frag" ${rnd})
        dm_fr_30=$(calculate_factor_depth ${dep_in} 30 12157105 "frag" ${rnd})
        dm_fr_40=$(calculate_factor_depth ${dep_in} 40 12157105 "frag" ${rnd})
        dm_fr_50=$(calculate_factor_depth ${dep_in} 50 12157105 "frag" ${rnd})

         dm_nm_1=$(calculate_factor_depth ${dep_in} 1  12157105 "norm" ${rnd})
         dm_nm_5=$(calculate_factor_depth ${dep_in} 5  12157105 "norm" ${rnd})
        dm_nm_10=$(calculate_factor_depth ${dep_in} 10 12157105 "norm" ${rnd})
        dm_nm_20=$(calculate_factor_depth ${dep_in} 20 12157105 "norm" ${rnd})
        dm_nm_30=$(calculate_factor_depth ${dep_in} 30 12157105 "norm" ${rnd})
        dm_nm_40=$(calculate_factor_depth ${dep_in} 40 12157105 "norm" ${rnd})
        dm_nm_50=$(calculate_factor_depth ${dep_in} 50 12157105 "norm" ${rnd})
    }

    #  Debug output to verify input minimum depth values
    if ${debug}; then
        debug_var \
            "dm_fr_1=${dm_fr_1}"   "dm_fr_5=${dm_fr_5}"   "dm_fr_10=${dm_fr_10}" \
            "dm_fr_20=${dm_fr_20}" "dm_fr_30=${dm_fr_30}" "dm_fr_40=${dm_fr_40}" \
            "dm_fr_50=${dm_fr_50}" \
            "dm_nm_1=${dm_nm_1}"   "dm_nm_5=${dm_nm_5}"   "dm_nm_10=${dm_nm_10}" \
            "dm_nm_20=${dm_nm_20}" "dm_nm_30=${dm_nm_30}" "dm_nm_40=${dm_nm_40}" \
            "dm_nm_50=${dm_nm_50}"
    fi

    #  Print the IP sample, input sample, and alpha value to the outfile
    prt_1="${fil_ip}\t${fil_in}\t${alpha}\t${eqn}\t"
    prt_2="${mass_ip}\t${mass_in}\t${vol_all}\t${vol_in}\t"
    prt_3="${dep_ip}\t${dep_in}\t${len_ip}\t${len_in}\t"
    prt_4="${dm_fr_1}\t${dm_fr_5}\t${dm_fr_10}\t${dm_fr_20}\t${dm_fr_30}\t"
    prt_5="${dm_fr_40}\t${dm_fr_50}\t"
    prt_6="${dm_nm_1}\t${dm_nm_5}\t${dm_nm_10}\t${dm_nm_20}\t${dm_nm_30}\t"
    prt_7="${dm_nm_40}\t${dm_nm_50}"

    echo -e \
        "${prt_1}${prt_2}${prt_3}${prt_4}${prt_5}${prt_6}${prt_7}" \
            >> "${fil_out}"
}


function set_interactive() {
    #  Set hardcoded paths, values, etc.
    ## WARNING: If interactive=true, change values as needed ##
    dir_bas="${HOME}/repos"
    dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
    dir_dat="${dir_rep}/data"
    dir_pro="${dir_dat}/processed"

    aligner="bowtie2"
    a_type="global"
    req_flg=true
    flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
    mapq=1

    dir_aln="${dir_pro}/align_${aligner}_${a_type}/flag-${flg}_mapq-${mapq}"
    dir_bam="${dir_aln}/sc"
    dir_cvg="${dir_pro}/compute_coverage"
    dir_out="${dir_cvg}/${aligner}_${a_type}_flag-${flg}_mapq-${mapq}/tables"

    #  Set hardcoded argument assignments
    threads=8
    ser_ip="${dir_bam}/IP_WT_Q_Hho1_6336.sc.bam,${dir_bam}/IP_WT_Q_Hho1_6337.sc.bam"
    ser_in="${dir_bam}/in_WT_Q_Hho1_6336.sc.bam,${dir_bam}/in_WT_Q_Hho1_6337.sc.bam"
    tbl_met="${dir_dat}/raw/docs/measurements_siqchip.tsv"
    eqn="6nd"
    fil_out="${dir_out}/alpha_test.tsv"
    rnd=24
    err_out="${dir_out}/logs"
    nam_job="calc_sf_alpha_${eqn}"
    env_nam="env_protocol"
    dir_scr="${HOME}/repos/protocol_chipseq_signal_norm/scripts"
}


#  Parse keyword arguments, assigning them to variables; most of the argument
#+ inputs are not checked, as this is performed by execute_*.sh and to a
#+ certain extent by the scripts submitted to SLURM
threads=1
ser_ip=""
ser_in=""
tbl_met=""
eqn="6nd"
fil_out=""
rnd=24
err_out=""
nam_job="calc_sf_alpha_${eqn}"
env_nam="env_protocol"
dir_scr="${HOME}/repos/protocol_chipseq_signal_norm/scripts"

show_help=$(cat << EOM
$(basename "${0}") takes the following keyword arguments:
   -t, --threads  Number of threads to use (default: ${threads}).
  -sp, --ser_ip   Comma-separated serialized string of sample IP BAM infiles.
  -sn, --ser_in   Comma-separated serialized string of corresponding sample
                  input BAM infiles.
  -tm, --tbl_met  TSV table of siQ-ChIP metrics.
  -eq, --eqn      Alpha equation to compute; options: '5', '5nd', '6', '6nd'
                  (default: '${eqn}').
  -fo, --fil_out  Outfile of sample siQ-ChIP alpha values.
  -fl, --flg_len  ...
  -fd, --flg_dep  ...
   -r, --rnd      Number of decimal places for rounding the siQ-ChIP alpha
                  scaling and minimum input depth factors (default: ${rnd}).
  -eo, --err_out  Directory to store stderr and stdout outfiles.
  -nj, --nam_job  Name of job (default: '${nam_job}').
  -en, --env_nam  Name of Conda/Mamba environment to activate (default:
                  '${env_nam}').
  -ds, --dir_scr  Path to directory containing workflow scripts, functions,
                  etc. (default: '${dir_scr}').
EOM
)

if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    if ! ${interactive}; then exit 0; fi
fi

# shellcheck disable=SC2034,SC2154
if ${interactive}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -t|--threads) threads="${2}"; shift 2 ;;
            -sp|--ser_ip)  ser_ip="${2}";  shift 2 ;;
            -sn|--ser_in)  ser_in="${2}";  shift 2 ;;
            -tm|--tbl_met) tbl_met="${2}"; shift 2 ;;
            -eq|--eqn)     eqn="${2}";     shift 2 ;;
            -fo|--fil_out) fil_out="${2}"; shift 2 ;;
             -r|--rnd)     rnd="${2}";     shift 2 ;;
            -eo|--err_out) err_out="${2}"; shift 2 ;;
            -nj|--nam_job) nam_job="${2}"; shift 2 ;;
            -en|--env_nam) env_nam="${2}"; shift 2 ;;
            -ds|--dir_scr) dir_scr="${2}"; shift 2 ;;
            *)
                echo "## Unknown argument passed: ${1} ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                exit 1
                ;;
        esac
    done
fi

#  Debug argument variable assignments
if ${debug}; then
    debug_var \
        "threads=${threads}" "ser_ip=${ser_ip}"   "ser_in=${ser_in}" \
        "tbl_met=${tbl_met}" "eqn=${eqn}"         "fil_out=${fil_out}" \
        "rnd=${rnd}"         "err_out=${err_out}" "nam_job=${nam_job}" \
        "env_nam=${env_nam}" "dir_scr=${dir_scr}"
fi

#  Activate environment
if [[ "${CONDA_DEFAULT_ENV}" != "${env_nam}" ]]; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${env_nam}" ||
        {
            echo "Error: Failed to activate environment '${env_nam}'." >&2
            exit 1
        }
fi

#  Assign and validate variables for scripts, functions, etc.
scr_alf="${dir_scr}/calculate_scaling_factor_alpha.py"
scr_aln="${dir_scr}/functions/count_alignments_bam.sh"
scr_frg="${dir_scr}/functions/calculate_frag_avg.sh"
scr_met="${dir_scr}/parse_metadata_siq_chip.py"
scr_min="${dir_scr}/functions/calculate_factor_depth.sh"

if ${debug}; then
    debug_var \
        "scr_alf=${scr_alf}" "scr_aln=${scr_aln}" "scr_frg=${scr_frg}" \
        "scr_met=${scr_met}" "scr_min=${scr_min}"
fi

for fil in "${scr_alf}" "${scr_aln}" "${scr_frg}" "${scr_met}" "${scr_min}"; do
    if [[ ! -f "${fil}" ]]; then
        echo "Error: Script, function, or subroutine not found: '${fil}'."
        if ! ${interactive}; then exit 1; fi
    fi
done

#  Source necessary subroutines
# shellcheck disable=SC1090
{
    source "${scr_aln}"
    source "${scr_frg}"
    source "${scr_min}"
}

#  Construct arrays from serialized strings
IFS=',' read -r -a arr_ip <<< "${ser_ip}"
IFS=',' read -r -a arr_in <<< "${ser_in}"

#  Debug output to check number of array elements and array element values
if ${debug}; then
    echo "\${#arr_ip[@]}=${#arr_ip[@]}" && echo ""
    echo "arr_ip=( ${arr_ip[*]} )"      && echo ""
    echo "\${#arr_in[@]}=${#arr_in[@]}" && echo ""
    echo "arr_in=( ${arr_in[*]} )"      && echo ""
fi

#  Determine and run mode: SLURM or GNU Parallel/serial
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    #  Mode: SLURM
    id_job=${SLURM_ARRAY_JOB_ID}
    id_tsk=${SLURM_ARRAY_TASK_ID}
    idx=$(( id_tsk - 1 ))

    #  Debug short names of environmental variables
    if ${debug}; then
        debug_var "id_job=${id_job}" "id_tsk=${id_tsk}" "idx=${idx}"
    fi

    #  Derive sample name, which is needed for 'set_logs'
    # shellcheck disable=SC2001
    {
        samp="${arr_ip[idx]}"
        samp="${samp##*/IP_}"
        samp=$(echo "${samp%.bam}" | sed 's:\.:_:g')
    }

    #  Debug output to verify sample name
    if ${debug}; then debug_var "samp=${samp}"; fi

    #  Run subroutine to set SLURM and symlinked/better-named log files
    set_logs "${id_job}" "${id_tsk}" "${samp}" "${err_out}"
    set_logs "${id_job}" "${id_tsk}" "${samp}" "${err_out}"
    ln -f "${err_ini}" "${err_dsc}"
    ln -f "${out_ini}" "${out_dsc}"

    #  Run processing subroutine
    # shellcheck disable=SC2086
    process_sample ${idx}

    #  Remove the initial SLURM stderr and stdout TXT outfiles
    rm "${err_ini}" "${out_ini}"
else
    #  Mode: GNU Parallel/serial
    for idx in "${!arr_ip[@]}"; do
        #  Run processing subroutine
        # shellcheck disable=SC2086
        process_sample ${idx}
    done
fi
