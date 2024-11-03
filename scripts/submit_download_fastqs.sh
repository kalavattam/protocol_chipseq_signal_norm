#!/bin/bash

#  submit_download_fastqs.sh
#  KA


#  Parse arguments ------------------------------------------------------------
if [[ -z "${SLURM_ARRAY_TASK_ID}" ]]; then
    #  Use arguments directly serial execution
    srr="${1}"
    url_1="${2}"
    url_2="${3}"  # Will be #N/A for single-end data
    dir_out="${4}"
    dir_sym="${5}"
    nam_cus="${6}"
else
    #  For SLURM job arrays, retrieve index from SLURM_ARRAY_TASK_ID
    idx=$(( SLURM_ARRAY_TASK_ID - 1 ))

    #  Convert environment variable strings back to arrays
    # shellcheck disable=SC2154
    {
        IFS=',' read -ra list_acc <<< "${str_acc}"
        IFS=',' read -ra list_url_1 <<< "${str_url_1}"
        IFS=',' read -ra list_url_2 <<< "${str_url_2}"
        IFS=',' read -ra list_cus <<< "${str_cus}"
    }

    #  Retrieve SRR, URLs, directory paths, and custom name based on imported
    #+ variables and array index
    srr="${list_acc[idx]}"
    url_1="${list_url_1[idx]}"
    url_2="${list_url_2[idx]}"
    # shellcheck disable=SC2269
    {
        dir_out="${dir_out}"
        dir_sym="${dir_sym}"
    }
    nam_cus="${list_cus[idx]}"
fi


#  Do the main work -----------------------------------------------------------
#  Download the first FASTQ file (always present)
echo "Downloading ${srr} from ${url_1}."
# curl -s -L -o "${dir_out}/${srr}_R1.fastq.gz" "${url_1}" ||
wget --quiet --no-verbose -O "${dir_out}/${srr}_R1.fastq.gz" "${url_1}" ||
    {
        echo "Error: Failed to download ${url_1}."
        exit 1
    }

#  If paired-end, download the second FASTQ file
if [[ "${url_2}" != "#N/A" ]]; then
    echo "Downloading ${srr} from ${url_2}."
    # curl -s -L -o "${dir_out}/${srr}_R2.fastq.gz" "${url_2}" ||
    wget --quiet --no-verbose -O "${dir_out}/${srr}_R2.fastq.gz" "${url_2}" ||
        {
            echo "Error: Failed to download ${url_2}."
            exit 1
        }
fi

#  Create symlinks to the downloaded file(s) using the custom name
ln -sf \
    "${dir_out}/${srr}_R1.fastq.gz" \
    "${dir_sym}/${nam_cus}_R1.fastq.gz" ||
        {
            echo \
                "Error: Failed to create symlink for" \
                "${dir_out}/${srr}_R1.fastq.gz"
            exit 1
        }

if [[ "${url_2}" != "#N/A" ]]; then
    ln -sf \
        "${dir_out}/${srr}_R2.fastq.gz" \
        "${dir_sym}/${nam_cus}_R2.fastq.gz" ||
        {
            echo \
                "Error: Failed to create symlink for" \
                "${dir_out}/${srr}_R2.fastq.gz"
            exit 1
        }
fi
