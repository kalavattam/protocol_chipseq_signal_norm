#!/bin/bash

#  Function to reset max_job, the maximum number of jobs to be run by SLURM at
#+ one time, if it exceeds the number of input files
reset_max_job() {
    local max_job="${1}"
    local num_fil="${2}"

    if [[ "${max_job}" -gt "${num_fil}" ]]; then
        echo_warning \
            "The maximum number of SLURM jobs to run at a time, ${max_job}," \
            "is greater than the number of infiles, ${num_fil}. Adjusting" \
            "max_job to ${num_fil}."
        echo "${num_fil}"  # Return adjusted max_job
    else
        echo "${max_job}"  # Return original max_job
    fi
}
