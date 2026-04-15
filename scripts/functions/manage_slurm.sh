#!/bin/bash

#  submit.sh
#  KA

#  Support functions for various 'submit' scripts


#  Print variable assignment(s)
function debug_var() { printf "%s\n\n" "$@" >&2; }


#  Activate user-supplied Conda/Mamba environment
function activate_env() {
    local env_nam="${1}"  # Name of environment to activate
    if [[ "${CONDA_DEFAULT_ENV}" != "${env_nam}" ]]; then
        # shellcheck disable=SC1091
        if ! \
            source "$(conda info --base)/etc/profile.d/conda.sh"
        then
            echo "Error: Failed to source 'conda.sh'." >&2
            return 1
        fi

        if ! \
            conda activate "${env_nam}"
        then
            echo "Error: Failed to activate environment: '${env_nam}'." >&2
            return 1
        fi
    fi
}


#  Set up SLURM log file paths and create hard-linked log files
function set_logs_slurm() {
    local id_job="${1}"    # SLURM job ID
    local id_tsk="${2}"    # SLURM task ID within the job array
    local samp="${3}"      # Sample name associated with log files
    local err_out="${4}"   # Directory for stderr and stdout log files
    local nam_job="${5}"   # Name of the job (used in log file naming)
    local err_ini out_ini  # SLURM initial stderr and stdout log files
    local err_dsc out_dsc  # Hard-linked stderr and stdout with descriptive names

    #  Check required inputs
    for var in id_job id_tsk samp err_out nam_job; do
        if [[ -z "${!var}" ]]; then
            echo "Error: '${var}' is required but is unset or empty." >&2
            return 1
        fi
    done

    #  Check that 'err_out' directory exists and is writable
    if [[ ! -d "${err_out}" || ! -w "${err_out}" ]]; then
        echo \
            "Error: Log directory does not exist or is not writable:" \
            "'${err_out}'." >&2
        return 1
    fi

    #  Set log paths
    err_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stderr.txt"
    out_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stdout.txt"
    err_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stderr.txt"
    out_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stdout.txt"

    #  Create hard-linked log files
    if ! \
        ln -f "${err_ini}" "${err_dsc}"
    then
        echo \
            "Error: Failed to create hard link: '${err_ini}' to" \
            "'${err_dsc}'." >&2
        return 1
    fi

    if ! \
        ln -f "${out_ini}" "${out_dsc}"
    then
        echo \
            "Error: Failed to create hard link: '${out_ini}' to" \
            "'${out_dsc}'." >&2
        return 1
    fi

    #  Return values
    echo "${err_ini},${out_ini},${err_dsc},${out_dsc}"
}


#  Check that a variable is not empty or unset
function validate_var() {
    local var_nam="${1}"   # Variable name (for error messages)
    local var_val="${2-}"  # Value to check for emptiness/unset state
    local idx="${3:-0}"    # Optional index (for arrays); defaults to '0'

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


#  Check that a file assigned to a varibale exists
function validate_file() {
    local var_nam="${1}"  # Variable name (for error messages)
    local var_val="${2}"  # Value (file) to check for existence
    local idx="${3:-0}"   # Optional index (for arrays); defaults to '0'

    if [[ ! -f "${var_val}" ]]; then
        if [[ "${idx}" -ne 0 ]]; then
            echo \
                "Error: '${var_nam}' does not exist for array index" \
                "'${idx}'." >&2
        else
            echo \
                "Error: '${var_nam}' does not exist." >&2
        fi
        return 1
    fi
}


#  Check that a file exists
function validate_var_file() {
    local var_nam="${1}"  # Name of variable to validate
    local pth_fil="${2}"  # File to check for existence
    local idx="${3:-0}"   # Optional index (for arrays); defaults to '0'

    if ! \
        validate_var "${var_nam}" "${pth_fil}" "${idx}"
    then
        if [[ "${idx}" -ne 0 ]]; then
            echo "Error: '${var_nam}' is unset or empty at index '${idx}'." >&2
        else
            echo "Error: '${var_nam}' is unset or empty." >&2
        fi
        return 1
    fi

    if ! \
        validate_file "${var_nam}" "${pth_fil}" "${idx}"
    then
        if [[ "${idx}" -ne 0 ]]; then
            echo \
                "Error: '${pth_fil}' does not exist for variable" \
                "'${var_nam}', index '${idx}'." >&2
        else
            echo \
                "Error: '${pth_fil}' does not exist for variable" \
                "'${var_nam}'." >&2
        fi
        return 1
    fi
}
