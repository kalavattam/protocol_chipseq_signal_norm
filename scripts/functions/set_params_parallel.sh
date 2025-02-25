#!/bin/bash

#  Function to determine the number of available CPU cores
function determine_cores() {
    local cores

    if command -v nproc &> /dev/null; then
        cores=$(nproc)
    elif command -v sysctl &> /dev/null; then
        cores=$(sysctl -n hw.ncpu 2> /dev/null)
    else
        echo "Error: Unable to determine the number of CPU cores." >&2
        return 1
    fi

    #  Check that output is a valid integer
    if [[ ! "${cores}" =~ ^[0-9]+$ ]]; then
        echo "Error: Failed to retrieve a valid core count." >&2
        return 1
    fi

    echo "${cores}"
}


#  Function to determine parallelization parameters for running GNU Parallel
function set_params_parallel() {
    local threads="${1}"
    local max_job="${2}"
    local par_job="${3}"
    local n_cores

    #  Get system CPU core count
    n_cores=$(determine_cores) || return 1

    #  Assign parallel jobs based on 'max_job', 'n_cores', and 'threads' 
    if [[ ${max_job} -gt 1 ]]; then
        par_job=${max_job}

        #  Avoid exceeding available CPU cores
        if [[ ${threads} -gt ${n_cores} ]]; then
            echo \
                "Warning: Requested threads (${threads}) exceed available" \
                "cores (${n_cores}). Using ${n_cores}." >&2
            threads=${n_cores}
        fi

        #  Ensure 'par_job' does not exceed available cores
        if [[ ${par_job} -gt ${n_cores} ]]; then
            echo \
                "Warning: Requested parallel jobs (${par_job}) exceed" \
                "available cores (${n_cores}). Using ${n_cores}." >&2
            par_job=${n_cores}
        fi

        #  Ensure that the total thread allocation does not exceed 'par_job'
        if [[ ${threads} -lt ${par_job} ]]; then
            echo \
                "Warning: Threads (${threads}) per job is too low relative" \
                "to parallel jobs (${par_job}). Adjusting to 1 thread per" \
                "job." >&2
            threads=1
        else
            threads=$(( threads / par_job ))  # Perform floor division
        fi
    else
        par_job=1
    fi

    #  Output modified values
    echo "${threads};${par_job}"
}
