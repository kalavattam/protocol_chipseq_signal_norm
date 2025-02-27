#!/bin/bash

#  Function to debug array contents
function debug_array_contents() {
    for arr_nam in "$@"; do
        #  Access the array indirectly using eval
        eval "arr=( \"\${${arr_nam}[@]}\" )"
        if [[ -n "${arr[*]}" ]]; then
            echo "  - ${arr_nam}=( ${arr[*]} )"
        fi
    done
}


#  Function to print parsed vector(s) and parameters for parallelization
function print_parallel_info() {
    local slurm="${1}"    # Boolean flag for SLURM: 'true' or 'false'
    local max_job="${2}"  # No. concurrent jobs: SLURM
    local par_job="${3}"  # No. concurrent jobs: GNU Parallel or serial
    local threads="${4}"  # No. threads per job
    shift 4

    echo "#########################################################"
    echo "## Parsed vector(s) and parameters for parallelization ##"
    echo "#########################################################"
    echo ""

    #  Pass remaining arguments as array names
    debug_array_contents "$@"

    #  Return relevant info
    if ${slurm}; then
        echo "  - Max concurrent jobs (SLURM): ${max_job}"
    elif [[ "${par_job}" -gt 1 ]]; then
        echo "  - Max concurrent jobs (GNU Parallel): ${par_job}"
    else
        echo "  - Jobs running in serial mode: ${par_job}"
    fi

    echo "  - Threads per job: ${threads}"
    echo ""
    echo ""
}
