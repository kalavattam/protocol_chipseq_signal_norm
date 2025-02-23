#!/bin/bash

#  submit_trim_fastqs.sh
#  KA


#  If true, run script in debug mode
debug=true


#  Define functions
#  Function to print variable assignment(s)
function debug_var() { printf "%s\n\n" "$@"; }


#  Function to activate user-supplied Conda/Mamba environment
function activate_env() {
    local env_nam="${1}"
    if [[ "${CONDA_DEFAULT_ENV}" != "${env_nam}" ]]; then
        eval "$(conda shell.bash hook)"
        conda activate "${env_nam}" || 
            {
                echo "Error: Failed to activate environment '${env_nam}'." >&2
                return 1
            }
    fi
}


#  Function to validate that a variable is not empty or unset
function validate_var() {
    local var_nam="${1}"
    local var_val="${2}"
    local idx="${3:-0}"  # Default to '0' if no index is given

    if [[ -z "${var_val}" ]]; then
        if [[ "${idx}" -ne 0 ]]; then
            echo \
                "Error: '${var_nam}' is empty or unset for array index" \
                "'${idx}'." >&2
        else
            echo "Error: '${var_nam}' is empty or unset." >&2
        fi
        return 1
    fi
}


#  Function to set up SLURM log file paths and create symlinked log files
function setup_slurm_logs() {
    local id_job="${1}"
    local id_tsk="${2}"
    local samp="${3}"
    local err_out="${4}"
    local nam_job="${5}"
    local err_ini out_ini err_dsc out_dsc

    #  Set log paths
    err_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stderr.txt"
    out_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stdout.txt"
    err_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stderr.txt"
    out_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stdout.txt"

    #  Create symlinked log files
    ln -f "${err_ini}" "${err_dsc}"
    ln -f "${out_ini}" "${out_dsc}"

    #  Return values
    echo "${err_ini};${out_ini};${err_dsc};${out_dsc}"
}


#  Function to debug, validate, and parse 'infile' to output variable
#+ assignments: 'fq_1', 'fq_2', and 'samp' (sample name)
function process_infile() {
    local infile="${1}"
    local sfx_se="${2}"
    local sfx_pe="${3}"
    local fq_1 fq_2 samp

    validate_var "infile"  "${infile}" || return 1

    if [[ "${infile}" == *,* ]]; then
        fq_1="${infile%%,*}"
        fq_2="${infile#*,}"
        samp="$(basename "${fq_1%%"${sfx_pe}"}")"
    else
        fq_1="${infile}"
        unset fq_2
        samp="$(basename "${fq_1%%"${sfx_se}"}")"
    fi

    validate_var "fq_1" "${fq_1}" || return 1
    if [[ -n "${fq_2}" ]]; then validate_var "fq_2" "${fq_2}" || return 1; fi

    #  Return values
    echo "${fq_1};${fq_2};${samp}"
}


#  Function to run Atria with explicitly passed argument values
function run_atria() {
    local threads="${1}"
    local fq_1="${2}"
    local fq_2="${3}"
    local dir_out="${4}"
    local err_out="${5}"
    local nam_job="${6}"
    local samp="${7}"

    # shellcheck disable=SC2046
    atria \
        -t "${threads}" \
        -r "${fq_1}" \
        $(if [[ -n ${fq_2} ]]; then echo "-R ${fq_2}"; fi) \
        -o "${dir_out}" \
        --length-range 35:500 \
             > "${err_out}/${nam_job}.${samp}.stdout.txt" \
            2> "${err_out}/${nam_job}.${samp}.stderr.txt"
}


#  Define the help message
show_help=$(cat << EOM
\${1}=env_nam     # Conda/Mamba environment to activate.
\${2}=threads     # Number of threads to use.
\${3}=str_infile  # Semicolon-separated serialized string of FASTQ files.
\${4}=dir_out     # Directory for FASTQ output files.
\${5}=sfx_se      # Suffix to strip from SE FASTQ files.
\${6}=sfx_pe      # Suffix to strip from PE FASTQ files.
\${7}=err_out     # Directory to store stderr and stdout output files.
\${8}=nam_job     # Job name.
EOM
)

#  Display help message if no arguments or help option is given
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
activate_env "${env_nam}" || exit 1

#  Reconstruct array from serialized string
validate_var "str_infile" "${str_infile}" || exit 1
IFS=';' read -r -a arr_infile <<< "${str_infile}"
unset IFS

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

    #  Retrieve the input file by indexing into the reconstructed input file
    #+ array
    infile="${arr_infile[idx]}"

    if ${debug}; then debug_var "infile=${infile}"; fi

    #  Run function to validate 'infile', using it to assign values to
    #+ variables 'fq_1', 'fq_2', and 'samp'
    IFS=';' read -r fq_1 fq_2 samp < <(
        process_infile "${infile}" "${sfx_se}" "${sfx_pe}"
    ) || exit 1
    unset IFS

    if ${debug}; then
        debug_var \
            "fq_1=${fq_1}" \
            "fq_2=${fq_2}" \
            "samp=${samp}"
    fi

    #  Run function to set SLURM and symlinked log files
    IFS=';' read -r err_ini out_ini err_dsc out_dsc < <(
        setup_slurm_logs \
        "${id_job}" "${id_tsk}" "${samp}" "${err_out}" "${nam_job}"
    ) || exit 1
    unset IFS

    if ${debug}; then
        debug_var \
            "err_ini=${err_ini}" \
            "out_ini=${out_ini}" \
            "err_dsc=${err_dsc}" \
            "out_dsc=${out_dsc}"
    fi

    #  Run function that calls Atria
    if ! \
        run_atria \
            "${threads}" "${fq_1}" "${fq_2}" "${dir_out}" "${err_out}" \
            "${nam_job}" "${samp}"
    then
        echo "Error: Failed to perform read trimming." >&2
        exit 1
    fi

    #  Remove the initial SLURM stderr and stdout TXT outfiles
    rm "${err_ini}" "${out_ini}"
else
    #  Mode: GNU Parallel/serial
    for idx in "${!arr_infile[@]}"; do
        #  Retrieve the input file by indexing into the reconstructed input
        #+ file array
        infile="${arr_infile[idx]}"

        if ${debug}; then debug_var "infile=${infile}"; fi

        #  Run function to validate 'infile', using it to assign values to
        #+ variables 'fq_1', 'fq_2', and 'samp'
        IFS=';' read -r fq_1 fq_2 samp < <(
            process_infile "${infile}" "${sfx_se}" "${sfx_pe}"
        ) || exit 1
        unset IFS

        if ${debug}; then
            debug_var \
                "fq_1=${fq_1}" \
                "fq_2=${fq_2}" \
                "samp=${samp}"
        fi

        #  Run function that calls Atria
        if ! \
            run_atria \
                "${threads}" "${fq_1}" "${fq_2}" "${dir_out}" "${err_out}" \
                "${nam_job}" "${samp}"
        then
            echo "Error: Failed to perform read trimming." >&2
            exit 1
        fi
    done
fi
