#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: submit_download_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


#  Require Bash >= 4.4 before doing any work
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be run under Bash >= 4.4." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

#  Run in safe mode, exiting on errors, unset variables, and pipe failures
set -euo pipefail

#  Define the help message
show_help=$(cat << EOM
\${1}=srr      <str>  NCBI SRA database run accession code.
\${2}=url_1    <str>  URL (FTP or HTTPS) for FASTQ file.
\${3}=url_2    <str>  Second FASTQ URL for PE data ("NA" for SE).
\${4}=dir_out  <str>  Directory to save FASTQ file(s).
\${5}=dir_sym  <str>  Directory for symlink(s) to FASTQ file(s).
\${6}=nam_cus  <str>  Custom name for symlink(s).
\${7}=err_out  <str>  Directory for stderr and stdout files.
\${8}=nam_job  <str>  Job name.
EOM
)

#  Display help message if a help option or no arguments are given
if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
    cat >&2 << EOM
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
    cat >&2 << EOM
error: '$(basename "${0}")' requires 8 positional arguments, ${msg}

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

#  Validate directories
for dir in "${dir_out}" "${dir_sym}" "${err_out}"; do
    if [[ ! -d "${dir}" ]]; then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "directory does not exist: '${dir}'" >&2
        exit 1
    fi
done
unset dir

#  Check for necessary program
if ! \
    command -v wget > /dev/null 2>&1
then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "required program is not in PATH: 'wget'." >&2
    exit 1
fi

#  Download FASTQ file(s)
echo "Downloading ${srr} from ${url_1}." >&2
if [[ "${url_2}" != "NA" ]]; then
    wget --progress=dot:mega -O "${dir_out}/${srr}_R1.fastq.gz" "${url_1}" \
         > "${err_out}/${nam_job}.${srr}_R1.stdout.txt" \
        2> "${err_out}/${nam_job}.${srr}_R1.stderr.txt" ||
        {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to download '${url_1}'." >&2
            exit 1
        }

    echo "Downloading ${srr} from ${url_2}." >&2
    wget --progress=dot:mega -O "${dir_out}/${srr}_R2.fastq.gz" "${url_2}" \
         > "${err_out}/${nam_job}.${srr}_R2.stdout.txt" \
        2> "${err_out}/${nam_job}.${srr}_R2.stderr.txt" ||
        {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to download '${url_2}'." >&2
            exit 1
        }
else
    wget --progress=dot:mega -O "${dir_out}/${srr}.fastq.gz" "${url_1}" \
         > "${err_out}/${nam_job}.${srr}.stdout.txt" \
        2> "${err_out}/${nam_job}.${srr}.stderr.txt" ||
        {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to download '${url_1}'." >&2
            exit 1
        }
fi

#  Create symlinks to the downloaded file(s) using the custom name
echo "Symlinking ${srr} to ${nam_cus}." >&2
if [[ "${url_2}" != "NA" ]]; then
    ln -sf \
        "${dir_out}/${srr}_R1.fastq.gz" \
        "${dir_sym}/${nam_cus}_R1.fastq.gz" ||
        {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to create symlink for" \
                "'${dir_out}/${srr}_R1.fastq.gz'" >&2
            exit 1
        }

    ln -sf \
        "${dir_out}/${srr}_R2.fastq.gz" \
        "${dir_sym}/${nam_cus}_R2.fastq.gz" ||
        {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to create symlink for" \
                "'${dir_out}/${srr}_R2.fastq.gz'" >&2
            exit 1
        }
else
    ln -sf \
        "${dir_out}/${srr}.fastq.gz" \
        "${dir_sym}/${nam_cus}.fastq.gz" ||
        {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to create symlink for" \
                "'${dir_out}/${srr}.fastq.gz'" >&2
            exit 1
        }
fi
