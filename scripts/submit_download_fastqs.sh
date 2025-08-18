#!/bin/bash

#  submit_download_fastqs.sh
#  KA


#  Define the help message
show_help=$(cat << EOM
\${1}=srr      # str: NCBI SRA database run accession code
\${2}=url_1    # str: URL (FTP or HTTPS) for FASTQ file
\${3}=url_2    # str: Second FASTQ URL for PE data ("NA" for SE)
\${4}=dir_out  # str: Directory to save FASTQ file(s)
\${5}=dir_sym  # str: Directory for symlink(s) to FASTQ file(s)
\${6}=nam_cus  # str: Custom name for symlink(s)
\${7}=err_out  # str: Directory for stderr and stdout files
\${8}=nam_job  # str: Job name
EOM
)

#  Display help message if a help option or no arguments are given
if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    cat << EOM
'$(basename "${0}")' requires 8 positional arguments:

The necessary positional arguments:
${show_help}

EOM
    exit 0
fi

#  Check for exactly 8 positional arguments
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

#  Parse positional arguments, assigning them to variables; most of the
#+ argument inputs are not checked, as this is performed by execute_*.sh and,
#+ if applicable, to a certain extent by the scripts submitted to SLURM
srr="${1}"
url_1="${2}"
url_2="${3}"
dir_out="${4}"
dir_sym="${5}"
nam_cus="${6}"
err_out="${7}"
nam_job="${8}"

#  Download FASTQ file(s)
echo "Downloading ${srr} from ${url_1}."
if [[ "${url_2}" != "NA" ]]; then
    wget --progress=dot:mega -O "${dir_out}/${srr}_R1.fastq.gz" "${url_1}" \
         > "${err_out}/${nam_job}.${srr}_R1.stdout.txt" \
        2> "${err_out}/${nam_job}.${srr}_R1.stderr.txt" ||
        {
            echo "Error: Failed to download ${url_1}."
            exit 1
        }

    echo "Downloading ${srr} from ${url_2}."
    wget --progress=dot:mega -O "${dir_out}/${srr}_R2.fastq.gz" "${url_2}" \
         > "${err_out}/${nam_job}.${srr}_R2.stdout.txt" \
        2> "${err_out}/${nam_job}.${srr}_R2.stderr.txt" ||
        {
            echo "Error: Failed to download ${url_2}."
            exit 1
        }
else
    wget --progress=dot:mega -O "${dir_out}/${srr}.fastq.gz" "${url_1}" \
         > "${err_out}/${nam_job}.${srr}.stdout.txt" \
        2> "${err_out}/${nam_job}.${srr}.stderr.txt" ||
        {
            echo "Error: Failed to download ${url_1}."
            exit 1
        }
fi

#  Create symlinks to the downloaded file(s) using the custom name
echo "Symlinking ${srr} to ${nam_cus}."
if [[ "${url_2}" != "NA" ]]; then
    ln -sf \
        "${dir_out}/${srr}_R1.fastq.gz" \
        "${dir_sym}/${nam_cus}_R1.fastq.gz" ||
        {
            echo \
                "Error: Failed to create symlink for" \
                "${dir_out}/${srr}_R1.fastq.gz"
            exit 1
        }

    ln -sf \
        "${dir_out}/${srr}_R2.fastq.gz" \
        "${dir_sym}/${nam_cus}_R2.fastq.gz" ||
        {
            echo \
                "Error: Failed to create symlink for" \
                "${dir_out}/${srr}_R2.fastq.gz"
            exit 1
        }
else
    ln -sf \
        "${dir_out}/${srr}.fastq.gz" \
        "${dir_sym}/${nam_cus}.fastq.gz" ||
        {
            echo \
                "Error: Failed to create symlink for" \
                "${dir_out}/${srr}.fastq.gz"
            exit 1
        }
fi
