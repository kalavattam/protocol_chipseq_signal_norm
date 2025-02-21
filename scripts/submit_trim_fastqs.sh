#!/bin/bash

#  If true, run script in debug mode
debug=true

#  Define functions
#  Print variable assignment
function debug_var() { printf "%s\n\n" "$@"; }


#  Validate that a variable is not empty or unset
function validate_var() {
    local var_nam="${1}"
    local var_val="${2}"
    local idx="${3}"

    if [[ -z "${var_val}" ]]; then
        echo \
            "Error: '${var_nam}' is empty or unset for array index" \
            "'${idx}'." >&2
        return 1
    fi
}


#  Set log file paths for SLURM job output, including symlinked and descriptive
#+ names
function set_logs() {
    local job_id="${1}"
    local task_id="${2}"
    local nam_smp="${3}"
    local dir_log="${4}"
    local nam_job="${5}"

    err_ini="${dir_log}/${nam_job}.${job_id}-${task_id}.stderr.txt"
    out_ini="${dir_log}/${nam_job}.${job_id}-${task_id}.stdout.txt"
    err_dsc="${dir_log}/${nam_job}.${nam_smp}.${job_id}-${task_id}.stderr.txt"
    out_dsc="${dir_log}/${nam_job}.${nam_smp}.${job_id}-${task_id}.stdout.txt"
}


#  Parse 'infile' into 'fq_1', 'fq_2', and 'samp' (sample name)
function parse_infile() {
    local infile="${1}"
    local sfx_se="${2}"
    local sfx_pe="${3}"

    if [[ "${infile}" == *,* ]]; then
        fq_1="${infile%%,*}"
        fq_2="${infile#*,}"
        samp="$(basename "${fq_1%%"${sfx_pe}"}")"
    else
        fq_1="${infile}"
        unset fq_2
        samp="$(basename "${fq_1%%"${sfx_se}"}")"
    fi
}


#  Subroutine to debug, validate, and parse input file
function process_infile() {
    local infile="${1}"
    local idx="${2}"
    local sfx_se="${3}"
    local sfx_pe="${4}"

    if ${debug}; then debug_var "infile=${infile}"; fi
    validate_var "infile"  "${infile}"  "${idx}" || return 1

    #  Parse infile into 'fq_1', 'fq_2', and 'samp'
    parse_infile "${infile}" "${sfx_se}" "${sfx_pe}" || return 1

    if ${debug}; then
        debug_var \
            "fq_1=${fq_1}" \
            "fq_2=${fq_2}" \
            "samp=${samp}"
    fi
}


#  Run Atria with global argument variable values
function run_atria() {
    # shellcheck disable=SC2046,2086
    atria \
        -t ${threads} \
        -r "${fq_1}" \
        $(if [[ -n ${fq_2} ]]; then echo "-R ${fq_2}"; fi) \
        -o "${dir_out}" \
        --length-range 35:500 \
             > ${err_out}/${nam_job}.${samp}.stdout.txt \
            2> ${err_out}/${nam_job}.${samp}.stderr.txt
}


#  Display help message if no arguments or help option is given
show_help=$(cat << EOM
\${1}=env_nam     # Mamba environment to activate.
\${2}=threads     # Number of threads to use.
\${3}=str_infile  # Semicolon-separated serialized string of FASTQ files.
\${4}=dir_out     # Directory for FASTQ output files.
\${5}=sfx_se      # Suffix to strip from SE FASTQ files.
\${6}=sfx_pe      # Suffix to strip from PE FASTQ files.
\${7}=err_out     # Directory to store stderr and stdout output files.
\${8}=nam_job     # Job name.
EOM
)

if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    cat << EOM
'$(basename "${0}")' requires 8 positional arguments:
${show_help}
EOM
    exit 0
fi

#  Check for exactly 8 arguments
if [[ $# -ne 8 ]]; then
    msg="but $# were supplied."
    [[ $# -eq 1 ]] && msg="but only $# was supplied."
    cat << EOM
Error: '$(basename "${0}")' requires 8 positional arguments, ${msg}

The necessary positional arguments:
${show_help}
EOM
    exit 1
fi

#  Parse positional arguments, assigning them to variables
#+ 
#+ Most of the argument inputs are not checked, as this is performed by
#+ execute_*.sh
env_nam="${1}"
threads="${2}"
str_infile="${3}"
dir_out="${4}"
sfx_se="${5}"
sfx_pe="${6}"
err_out="${7}"
nam_job="${8}"

##########
## TEST ##
##########
# env_nam=env_protocol
# threads=4
# str_infile='/Users/kalavattam/repos/protocol_chipseq_signal_norm/data/symlinked/IP_WT_Q_Hmo1_7750_R1.fastq.gz,/Users/kalavattam/repos/protocol_chipseq_signal_norm/data/symlinked/IP_WT_Q_Hmo1_7750_R2.fastq.gz'
# dir_out=/Users/kalavattam/repos/protocol_chipseq_signal_norm/data/processed/trim_fastqs
# sfx_se=.fastq.gz
# sfx_pe=_R1.fastq.gz
# err_out=/Users/kalavattam/repos/protocol_chipseq_signal_norm/data/processed/trim_fastqs/logs
# nam_job=trim_fastqs
##########
## TEST ##
##########

#  Debug argument variable assignments
if ${debug}; then
    debug_var \
        "env_nam=${env_nam}" \
        "threads=${threads}" \
        "str_infile=${str_infile}" \
        "dir_out=${dir_out}" \
        "sfx_se=${sfx_se}" \
        "sfx_pe=${sfx_pe}" \
        "err_out=${err_out}" \
        "nam_job=${nam_job}"
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

#  Reconstruct array from serialized string
IFS=';' read -r -a arr_infile <<< "${str_infile}"
# for infile in "${arr_infile[@]}"; do echo "${infile}"; done  # unset infile

#  Debug output to check number of array elements and array element values
if ${debug}; then
    echo "\${#arr_infile[@]}=${#arr_infile[@]}" && echo ""
    echo "arr_infile=( ${arr_infile[*]} )"      && echo ""
fi

#  Determine and run mode: SLURM or GNU Parallel/serial
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    #  Mode: SLURM
    id_job="${SLURM_ARRAY_JOB_ID}"
    id_tsk="${SLURM_ARRAY_TASK_ID}"
    idx=$(( id_tsk - 1 ))

    infile="${arr_infile[idx]}"

    #  Run subroutine to debug and validate 'infile', and then parse it into
    #+ 'fq_1', 'fq_2', and 'samp'
    process_infile "${infile}"  "${idx}" "${sfx_se}" "${sfx_pe}" || exit 1

    #  Run subroutine to set SLURM and symlinked/better-named log files
    set_logs "${id_job}" "${id_tsk}" "${samp}" "${err_out}" "${nam_job}"
    ln -f "${err_ini}" "${err_dsc}"
    ln -f "${out_ini}" "${out_dsc}"

    #  Run Atria
    run_atria "${fq_1}" "${fq_2}" "${dir_out}"

    #  Remove the initial SLURM stderr and stdout TXT outfiles
    rm "${err_ini}" "${out_ini}"
else
    #  Mode: GNU Parallel/serial
    for idx in "${!arr_infile[@]}"; do
        # idx=0
        infile="${arr_infile[idx]}"

        #  Run subroutine to debug and validate 'infile', and then parse it
        #+ into 'fq_1', 'fq_2', and 'samp'
        process_infile "${infile}"  "${idx}" "${sfx_se}" "${sfx_pe}" || exit 1

        #  Run Atria
        run_atria "${fq_1}" "${fq_2}" "${dir_out}"
    done
fi
