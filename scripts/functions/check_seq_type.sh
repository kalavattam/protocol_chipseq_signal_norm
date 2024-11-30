#!/bin/bash

function check_seq_type() {
    local bam="${1}"
    local line_pg
    local show_help

    show_help=$(cat << EOM
--------------
check_seq_type
--------------

Description:
  Determines whether a BAM file contains single- or paired-end alignments 
  based on the program group (@PG) information in the BAM header.

Usage:
  check_seq_type <bam>

Positional parameters:
  1, bam (str): The full path of the BAM file.

Returns:
  str: 'single' if the BAM file contains single-end alignments, 'paired' if it
       contains paired-end alignments.

Dependencies:
  - Bash or Zsh
  - Samtools

Note:
  - RegEx to detect presence of paired-end alignments in the BAM header:
    '^@PG.*ID:(bowtie2|bwa).*(_R2|_2)'.

Example:
  Determine the sequencing type of a BAM file.
  \`\`\`
  ❯ seq_type=\$(check_seq_type "\${HOME}/path/to/file.bam")
  ❯ echo "The sequencing type is '\${seq_type}'."
  The sequencing type is 'paired'.
  \`\`\`
EOM
    )

    #  Parse and check function parameters
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    #  Check that BAM file is provided and exists
    if [[ -z "${bam}" ]]; then
        echo "Error: Positional parameter 1, bam, is required." >&2
        return 1
    fi

    if [[ ! -f "${bam}" ]]; then
        echo "Error: Positional parameter 1, bam, does not exist." >&2
        return 1
    fi

    #  Check that Samtools is available in PATH
    if ! command -v samtools > /dev/null; then
        echo "Error: Samtools is not installed or not in PATH." >&2
        return 1
    fi

    #  Extract @PG lines from BAM header and analyze sequencing type
    line_pg=$(samtools view -H "${bam}" | grep -E '^@PG.*ID:(bowtie2|bwa)')

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

    # #  Use Samtools to view the BAM header and grep for 'ID:bowtie2' or 'ID:bwa'
    # #+ to determine the aligner used
    # if \
    #     samtools view -H "${bam}" \
    #         | grep -E '^@PG.*ID:(bowtie2|bwa)' \
    #             >/dev/null
    # then
    #     #  Check for '_R2' or '_2' in the program group ID to identify paired-
    #     #+ end reads
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
    #         "Error: BAM file does not contain recognized program group (@PG)" \
    #         "information with 'ID:bowtie2' or 'ID:bwa'." >&2
    #     return 1
    # fi
}
