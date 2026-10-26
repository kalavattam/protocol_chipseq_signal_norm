#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: make.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5, Opus 5.5).
#
# Distributed under the MIT license.


# Require Bash >= 4.4 before doing any work.
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

# Run in safe mode, exiting on errors, unset variables, and pipe failures.
set -euo pipefail


# Resolve paths relative to 'tests/fixtures'.
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_fix="${dir_scr}"

# Source shared fixture-generation helpers.
# shellcheck source=tests/support/fixture_helpers.sh
source "${dir_scr}/../../support/fixture_helpers.sh"


# Define fixture paths, contig names, and the shared reference sequence.
dir_sam="${dir_fix}/sam"
dir_ref="${dir_fix}/reference"
dir_bam="${dir_fix}/bam"
dir_cram="${dir_fix}/cram"

sam="${dir_sam}/filter_sc_sp.sam"
ref="${dir_ref}/filter_sc_sp.fa"
fai="${ref}.fai"
bam="${dir_bam}/filter_sc_sp.bam"
cram="${dir_cram}/filter_sc_sp.cram"
cram_par_1="${dir_cram}/filter_sc_sp_parallel_1.cram"
cram_par_2="${dir_cram}/filter_sc_sp_parallel_2.cram"

contigs=(
    I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI Mito
    SP_I SP_II SP_III SP_II_TG SP_MTR SP_Mito
    chrUn
)

seq_100="ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"


# Remove stale generated fixture outputs. Every output below is rewritten
# unconditionally, so this sweep is not what makes regeneration correct today;
# it is what keeps regeneration correct after a later revision stops writing
# one of them.
rm_files \
    "${dir_fix}" \
    "${sam}" \
    "${ref}" \
    "${fai}" \
    "${bam}" \
    "${bam}.bai" \
    "${cram_par_1}" \
    "${cram_par_1}.crai" \
    "${cram_par_2}" \
    "${cram_par_2}.crai" \
    "${cram}" \
    "${cram}.crai"

# Create fixture output directories.
mkdirs "${dir_sam}" "${dir_ref}" "${dir_bam}" "${dir_cram}"


# Write the tiny reference FASTA and matching FASTA index.
: > "${ref}"
: > "${fai}"
offset=0
for contig in "${contigs[@]}"; do
    header=">${contig}"
    printf '%s\n%s\n' "${header}" "${seq_100}" >> "${ref}"
    offset=$(( offset + ${#header} + 1 ))
    printf '%s\t100\t%d\t100\t101\n' "${contig}" "${offset}" >> "${fai}"
    offset=$(( offset + 100 + 1 ))
done

# Write SAM headers and one representative read per retained/dropped class.
{
    write_sam_line '@HD' 'VN:1.6' 'SO:coordinate'

    for contig in "${contigs[@]}"; do
        write_sam_line '@SQ' "SN:${contig}" 'LN:100'
    done

    write_sam_line \
        'r_sc_I' '0' 'I' '1' '60' '10M' \
        '*' '0' '0' 'ACGTACGTAA' 'FFFFFFFFFF'

    write_sam_line \
        'r_sc_mito' '0' 'Mito' '1' '60' '10M' \
        '*' '0' '0' 'ACGTACGTAA' 'FFFFFFFFFF'

    write_sam_line \
        'r_sp_i' '0' 'SP_I' '1' '60' '10M' \
        '*' '0' '0' 'ACGTACGTAA' 'FFFFFFFFFF'

    write_sam_line \
        'r_sp_tg' '0' 'SP_II_TG' '1' '60' '10M' \
        '*' '0' '0' 'ACGTACGTAA' 'FFFFFFFFFF'

    write_sam_line \
        'r_sp_mtr' '0' 'SP_MTR' '1' '60' '10M' \
        '*' '0' '0' 'ACGTACGTAA' 'FFFFFFFFFF'

    write_sam_line \
        'r_sp_mito' '0' 'SP_Mito' '1' '60' '10M' \
        '*' '0' '0' 'ACGTACGTAA' 'FFFFFFFFFF'

    write_sam_line \
        'r_other' '0' 'chrUn' '1' '60' '10M' \
        '*' '0' '0' 'ACGTACGTAA' 'FFFFFFFFFF'
} > "${sam}"

# Build the BAM and CRAM inputs from the SAM, with their indexes; the CRAM
# index is the runner's sentinel, so it is written last.
require_env env_protocol "for filter-alignments BAM and CRAM fixtures."
require_cmd samtools "in 'env_protocol' to generate BAM/CRAM fixtures."

samtools view -b -o "${bam}" "${sam}"
samtools index "${bam}"
samtools view -C -T "${ref}" -o "${cram}" "${sam}"

# Copy the CRAM into two samples, so that the GNU Parallel test runs two jobs.
cp "${cram}" "${cram_par_1}"
cp "${cram}" "${cram_par_2}"
samtools index "${cram_par_1}"
samtools index "${cram_par_2}"
samtools index "${cram}"


succeed "generated filter-alignments fixtures under ${dir_fix}"
