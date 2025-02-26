#!/bin/bash


#  Function to determine the suffix for the second PE FASTQ file
#+ based on the user-supplied PE suffix pattern
function get_paired_suffix() {
    local file="${1}"    # File to check existence
    local sfx_pe="${2}"  # User-supplied PE suffix (_R1.fastq.gz, etc.)
    local sfx_pe_2       # Modified suffix for second read

    #  First, try replacing '1' next to 'r' or 'R' (case insensitive)
    if [[ "${sfx_pe}" =~ ([rR])1 ]]; then
        sfx_pe_2="${sfx_pe/${BASH_REMATCH[1]}1/${BASH_REMATCH[1]}2}"

    #  If not found, try replacing '1' before '.fastq' or similar extensions
    elif [[ "${sfx_pe}" =~ 1(\.f(ast)?q) ]]; then
        sfx_pe_2="${sfx_pe/1${BASH_REMATCH[1]}/2${BASH_REMATCH[1]}}"

    #  If still not found, try replacing '1' before '.atria'
    elif [[ "${sfx_pe}" =~ 1(\.atria) ]]; then
        sfx_pe_2="${sfx_pe/1${BASH_REMATCH[1]}/2${BASH_REMATCH[1]}}"

    #  If neither pattern is found, iterate through '1' replacements
    else
        local sfx_mod
        local i

        i=0
        while [[ ${i} -lt ${#sfx_pe} ]]; do
            if [[ "${sfx_pe:${i}:1}" == "1" ]]; then
                sfx_mod="${sfx_pe:0:${i}}2${sfx_pe:$((i + 1))}"
                file_2="${file/${sfx_pe}/${sfx_mod}}"

                if [[ -f "${file_2}" ]]; then
                    sfx_pe_2="${sfx_mod}"
                    break
                fi
            fi
            (( i++ )) || true
        done

        #  If no valid file is found, fallback to a simple first-occurrence
        #+ replacement
        [[ -z "${sfx_pe_2}" ]] && sfx_pe_2="${sfx_pe/1/2}"
    fi

    #  Return suffix determined for second PE FASTQ file
    echo "${sfx_pe_2}"
}


#  Function to validate and parse a serialized string of FASTQ files, ensuring
#+ correct formatting, file existence, and suffix validity for both single-end
#+ (SE) and paired-end (PE) reads
function check_string_fastqs() {
    local infiles="${1}"  # Serialized string of FASTQ file paths (semicolon-separated)
    local sfx_se="${2}"   # Expected suffix for SE FASTQ files
    local sfx_pe="${3}"   # Expected suffix for PE FASTQ files
    local fq_1  # First FASTQ file
    local fq_2  # Second FASTQ file (for PE reads)

    #  Convert semicolon-separated string into an array
    IFS=';' read -r -a arr_fq <<< "${infiles}"
    unset IFS

    #  Validate each FASTQ entry
    for fq in "${arr_fq[@]}"; do
        # fq="${arr_fq[0]}"

        #  Run AWK to check number of fields when splitting by comma
        num_prt=$(awk -F',' '{ print NF }' <<< "${fq}")

        if [[ "${num_prt}" -gt 2 ]]; then
            echo \
                "Error: Too many comma-separated values in FASTQ entry:" \
                "'${fq}'." >&2
            return 1
        fi

        #  If PE, expect exactly 2 fields
        if [[ "${num_prt}" -eq 2 ]]; then
            IFS=',' read -r fq_1 fq_2 <<< "${fq}"
            unset IFS
        else
            fq_1="${fq}"
            fq_2="#N/A"
        fi

        #  Ensure FASTQ file(s) exist
        if [[ ! -f "${fq_1}" ]]; then
            echo "Error: FASTQ file does not exist: '${fq_1}'." >&2
            return 1
        fi

        #  Validate suffix for SE FASTQ file
        if [[ "${fq_2}" == "#N/A" && "${fq_1}" != *"${sfx_se}" ]]; then
            echo \
                "Error: SE file '${fq_1}' does not match expected suffix:" \
                "'${sfx_se}'." >&2
            return 1
        fi

        #  Validate suffix for FASTQ #1 of PE FASTQ file pair
        if [[ "${fq_2}" != "#N/A" && "${fq_1}" != *"${sfx_pe}" ]]; then
            echo \
                "Error: PE file '${fq_1}' does not match expected suffix:" \
                "'${sfx_pe}'." >&2
            return 1
        fi

        #  If PE, dynamically derive 'fq_2' and check file existence
        if [[ "${fq_2}" != "#N/A" ]]; then
            der_2="$(get_paired_suffix "${fq_1}" "${sfx_pe}")"
            fq_2="${fq_1/${sfx_pe}/${der_2}}"
            
            if [[ ! -f "${fq_2}" ]]; then
                echo \
                    "Error: Expected PE FASTQ file does not exist:" \
                    "'${fq_2}'." >&2
                return 1
            fi
        fi
    done

    #  Return success
    return 0
}
