#!/bin/bash

#  Function to validate serialized string of FASTQ files
function validate_fastqs() {
    local str_fq="${1}"  # Serialized string containing FASTQ file paths
    local sfx_se="${2}"  # Expected suffix for SE FASTQ files
    local sfx_pe="${3}"  # Expected suffix for PE FASTQ files
    local fq_1
    local fq_2

    #  Reconstruct array from serialized string
    IFS=';' read -r -a arr_fq <<< "${str_fq}"
    
    #  Check that each FASTQ infile exists; if not, exit
    for fq in "${arr_fq[@]}"; do
        # fq="${arr_fq[0]}"

        #  Check that the FASTQ infile contains a comma (indicating PE)
        if [[ "${fq}" == *,* ]]; then
            fq_1="${fq%%,*}"
            fq_2="${fq#*,}"

            #  Ensure there is only one comma
            if [[ "${fq_2}" == *,* ]]; then
                echo \
                    "Error: It appears there is a problem with the" \
                    "serialized string output from running find_files.sh in" \
                    "'--fastqs' mode. More than one comma was found in a" \
                    "specific substring: '${fq}'."
                return 1
            fi

            #  Validate sample-specific paired-end FASTQ files exist
            if [[ ! -f "${fq_1}" ]]; then
                echo "Error: File does not exist: '${fq_1}'." >&2
                return 1
            fi

            if [[ ! -f "${fq_2}" ]]; then
                echo "Error: File does not exist: '${fq_2}'." >&2
                return 1
            fi

            #  Validate presence of file suffix in fq_1 assignment
            if [[ "${fq_1}" != *"${sfx_pe}" ]]; then
                echo \
                    "Error: Suffix '${sfx_pe}' not found in file name" \
                    "'${fq_1}'. Check '--sfx_pe' assignment." >&2
                return 1
            fi
        else
            fq_1="${fq}"
            unset fq_2  # Only affects local scope

            #  Validate sample-specific single-end FASTQ file exists
            if [[ ! -f "${fq_1}" ]]; then
                echo "Error: File does not exist: '${fq_1}'." >&2
                return 1
            fi

            #  Validate presence of file suffix in fq_1 assignment
            if [[ "${fq_1}" != *"${sfx_se}" ]]; then
                echo \
                    "Error: Suffix '${sfx_se}' not found in file name" \
                    "'${fq_1}'. Check '--sfx_se' assignment." >&2
                return 1
            fi
        fi
    done
}
