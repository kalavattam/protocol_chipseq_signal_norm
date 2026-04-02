#!/bin/bash
# -*- coding: utf-8 -*-
# 
# Script: process_sequences.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
# 
# Distributed under the MIT license.


function check_seq_type() {
    local bam="${1}"
    local line_pg
    local show_help

    show_help=$(cat << EOM
Usage:
  check_seq_type <bam>

Description:
  Determine whether a BAM file contains single- or paired-end alignments based on the program group (@PG) information in the BAM header.

Positional argument:
  1  bam  <str>  Full path to the BAM file.

Returns:
  String 'single' if the BAM file contains single-end alignments, 'paired' if it contains paired-end alignments.

Dependencies:
  - Bash or Zsh
  - Samtools

Note:
  - RegEx to detect presence of paired-end alignments in the BAM header: '^@PG.*ID:(bowtie2|bwa).*(_R2|_2)'.

Example:
  Determine the sequencing type of a BAM file:
  '''bash
  seq_type=\$(check_seq_type "\${HOME}/path/to/file.bam")
  echo "The sequencing type is '\${seq_type}'."
  '''

  '''txt
  0 ❯ The sequencing type is 'paired'.
  '''
EOM
    )

    #  Parse and check function argument
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    #  Check that BAM file is provided and exists
    if [[ -z "${bam}" ]]; then
        echo "Error: Positional argument 1, 'bam', is required." >&2
        return 1
    fi

    if [[ ! -f "${bam}" ]]; then
        echo "Error: Positional argument 1, 'bam', does not exist." >&2
        return 1
    fi

    #  Check that Samtools is available in PATH
    if ! command -v samtools > /dev/null; then
        echo "Error: Samtools is either not installed or not in PATH." >&2
        return 1
    fi

    #  Extract @PG lines from BAM header and analyze sequencing type
    line_pg=$(
        samtools view -H "${bam}" | grep -E '^@PG.*ID:(bowtie2|bwa)'
    )

    if [[ -z "${line_pg}" ]]; then
        echo \
            "Error: No recognized program group (@PG) with 'ID:bowtie2' or" \
            "'ID:bwa' in the BAM header." >&2
        return 1
    fi

    #  Determine sequencing type based on paired-end indicators '_R2' or '_2'
    if echo "${line_pg}" | grep -E '(_R2|_2)' > /dev/null; then
        echo "paired"
    else
        echo "single"
    fi

    # #  Use Samtools to view the BAM header and grep for 'ID:bowtie2' or
    # #+ 'ID:bwa' to determine the aligner used
    # if \
    #     samtools view -H "${bam}" \
    #         | grep -E '^@PG.*ID:(bowtie2|bwa)' \
    #             >/dev/null
    # then
    #     #  Check for '_R2' or '_2' in the program group ID to identify
    #     #+ paired-end reads
    #     if \
    #         samtools view -H "${bam}" \
    #             | grep -E '^@PG.*ID:(bowtie2|bwa).*(_R2|_2)' \
    #                 >/dev/null
    #     then
    #         echo "paired"
    #     else
    #         echo "single"
    #     fi
    # else
    #     echo \
    #         "Error: BAM file does not contain recognized program group" \
    #         "(@PG) information with 'ID:bowtie2' or 'ID:bwa'." >&2
    #     return 1
    # fi
}


#  Validate and parse a serialized string of FASTQ files, ensuring correct
#+ formatting, file existence, and suffix validity for both single-end (SE) and
#+ paired-end (PE) reads
function check_string_fastqs() {
    local infiles="${1}"  # Serialized string of FASTQ file paths (semicolon-separated)
    local sfx_se="${2}"   # Expected suffix for SE FASTQ files
    local sfx_pe="${3}"   # Expected suffix for PE FASTQ files
    local fq_1            # First FASTQ file
    local fq_2            # Second FASTQ file (for PE reads)

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
            fq_2="NA"
        fi

        #  Ensure FASTQ file(s) exist
        if [[ ! -f "${fq_1}" ]]; then
            echo "Error: FASTQ file does not exist: '${fq_1}'." >&2
            return 1
        fi

        #  Validate suffix for SE FASTQ file
        if [[ "${fq_2}" == "NA" && "${fq_1}" != *"${sfx_se}" ]]; then
            echo \
                "Error: SE file '${fq_1}' does not match expected suffix:" \
                "'${sfx_se}'." >&2
            return 1
        fi

        #  Validate suffix for FASTQ #1 of PE FASTQ file pair
        if [[ "${fq_2}" != "NA" && "${fq_1}" != *"${sfx_pe}" ]]; then
            echo \
                "Error: PE file '${fq_1}' does not match expected suffix:" \
                "'${sfx_pe}'." >&2
            return 1
        fi

        #  If PE, dynamically derive 'fq_2' and check file existence
        if [[ "${fq_2}" != "NA" ]]; then
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

    return 0
}


#  Determine the suffix for the second PE FASTQ file based on the user-supplied
#+ PE suffix pattern
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


#  Output a string of FASTQ files such that sample-specific ones are separated
#+ by semicolons and, within semicolon-delimited substrings, sample-specific
#+ FASTQ pairs are separated by commas
function pair_fastqs() {
    awk '
        BEGIN { OFS = "" }
        {
            #  Check for lines ending with "_R1.atria" or simply "_R1" followed
            #+ by ".fastq.gz" or ".fq.gz"
            if ($0 ~ /_R1(\.atria)?\.(fastq|fq)(\.gz)?$/) {
                r1 = $0
                getline
                #  Ensure the next line is a matching "_R2" file
                if ($0 ~ /_R2(\.atria)?\.(fastq|fq)(\.gz)?$/) {
                    r2 = $0
                    print r1 ",", r2, ";"
                } else {
                    #  Print an error message and exit with a failure code if
                    #+ no matching "_R2" file is found
                    print "Error: Missing R2 file for " r1 > "/dev/stderr"
                    exit 1
                }
            } else {
                #  Handle unpaired files
                print $0 ";"
            }
        }
    ' 
}
