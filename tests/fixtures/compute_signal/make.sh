#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: make.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in design, development,
# and documentation, with all output reviewed, edited, and approved by the
# author.
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


#  Resolve paths relative to 'tests/fixtures'
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_fix="${dir_scr}"

#  Source shared fixture-generation helpers
# shellcheck source=tests/support/fixture_helpers.sh
source "${dir_scr}/../../support/fixture_helpers.sh"


#  Define fixture directories, bedGraph paths, and alignment paths
dir_bdg="${dir_fix}/bedgraph"
dir_ref="${dir_fix}/reference"
dir_sam="${dir_fix}/sam"
dir_sam_se="${dir_sam}/se"
dir_sam_pe="${dir_sam}/pe"
dir_bam="${dir_fix}/bam"
dir_bam_se="${dir_bam}/se"
dir_bam_pe="${dir_bam}/pe"
dir_cram="${dir_fix}/cram"
dir_cram_se="${dir_cram}/se"
dir_cram_pe="${dir_cram}/pe"

fil_bg_A="${dir_bdg}/ratio_A.bdg"
fil_bg_B="${dir_bdg}/ratio_B.bdg"
fil_bg_hdr_A="${dir_bdg}/ratio_headers_A.bdg"
fil_bg_hdr_B="${dir_bdg}/ratio_headers_B.bdg"

fil_ref="${dir_ref}/tiny.fa"
fil_sam_se="${dir_sam_se}/tiny_se.sam"
fil_sam_pe="${dir_sam_pe}/tiny_pe.sam"

fil_bam_se="${dir_bam_se}/tiny_se.bam"
fil_bam_pe="${dir_bam_pe}/tiny_pe.bam"

fil_cram_se="${dir_cram_se}/tiny_se.cram"
fil_cram_pe="${dir_cram_pe}/tiny_pe.cram"

env_req="env_protocol"


#  Require the project environment for samtools-backed fixtures
require_env "${env_req}" "for compute-signal fixtures."

#  Create fixture output directories
mkdirs \
    "${dir_bdg}" \
    "${dir_ref}" \
    "${dir_sam}" \
    "${dir_sam_se}" \
    "${dir_sam_pe}" \
    "${dir_bam}" \
    "${dir_bam_se}" \
    "${dir_bam_pe}" \
    "${dir_cram}" \
    "${dir_cram_se}" \
    "${dir_cram_pe}"

#  Write tiny bedGraph fixtures for ratio-mode edge cases
{
    printf "I\t0\t10\t4\n"
    printf "I\t10\t20\t0\n"
    printf "I\t20\t30\t5\n"
    printf "I\t30\t40\t0\n"
    printf "I\t40\t50\t2\n"
    printf "I\t50\t60\t1\n"
    printf "I\t60\t70\t1\n"
    printf "I\t70\t80\t1\n"
} > "${fil_bg_A}"

{
    printf "I\t0\t10\t2\n"
    printf "I\t10\t20\t2\n"
    printf "I\t20\t30\t0\n"
    printf "I\t30\t40\t0\n"
    printf "I\t40\t50\t0.5\n"
    printf "I\t50\t60\t0.04\n"
    printf "I\t60\t70\t3\n"
    printf "I\t70\t80\t1\n"
} > "${fil_bg_B}"

#  Write bedGraph fixtures with skippable header/comment lines
{
    printf "track type=bedGraph name=\"ratio_A\"\n"
    printf "browser position I:0-80\n"
    printf "# comment in A\n"
    printf "customHeader sample=A\n"
    printf "I\t0\t10\t4\n"
    printf "I\t10\t20\t0\n"
    printf "  customHeader indented=A\n"
    printf "I\t20\t30\t5\n"
    printf "I\t30\t40\t0\n"
    printf "I\t40\t50\t2\n"
    printf "I\t50\t60\t1\n"
    printf "I\t60\t70\t1\n"
    printf "I\t70\t80\t1\n"
} > "${fil_bg_hdr_A}"

{
    printf "track type=bedGraph name=\"ratio_B\"\n"
    printf "browser position I:0-80\n"
    printf "# comment in B\n"
    printf "customHeader sample=B\n"
    printf "I\t0\t10\t2\n"
    printf "I\t10\t20\t2\n"
    printf "  customHeader indented=B\n"
    printf "I\t20\t30\t0\n"
    printf "I\t30\t40\t0\n"
    printf "I\t40\t50\t0.5\n"
    printf "I\t50\t60\t0.04\n"
    printf "I\t60\t70\t3\n"
    printf "I\t70\t80\t1\n"
} > "${fil_bg_hdr_B}"

#  Write tiny reference FASTA used for BAM/CRAM fixture generation
cat > "${fil_ref}" << EOM
>I
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOM

#  Write single-end SAM records for forward and reverse fragments
{
    write_sam_line '@HD' 'VN:1.6' 'SO:coordinate'
    write_sam_line '@SQ' 'SN:I' 'LN:80'

    write_sam_line \
        'se_fwd' '0' 'I' '1' '60' '10M' \
        '*' '0' '0' 'ACGTACGTAC' 'FFFFFFFFFF'

    write_sam_line \
        'se_rev' '16' 'I' '21' '60' '10M' \
        '*' '0' '0' 'ACGTACGTAC' 'FFFFFFFFFF'
} > "${fil_sam_se}"

#  Write paired-end SAM records for two fragment positions
{
    write_sam_line '@HD' 'VN:1.6' 'SO:coordinate'
    write_sam_line '@SQ' 'SN:I' 'LN:80'

    write_sam_line \
        'pe_1' '99' 'I' '11' '60' '10M' \
        '=' '31' '30' 'ACGTACGTAC' 'FFFFFFFFFF'

    write_sam_line \
        'pe_1' '147' 'I' '31' '60' '10M' \
        '=' '11' '-30' 'ACGTACGTAC' 'FFFFFFFFFF'

    write_sam_line \
        'pe_2' '99' 'I' '41' '60' '10M' \
        '=' '51' '20' 'ACGTACGTAC' 'FFFFFFFFFF'

    write_sam_line \
        'pe_2' '147' 'I' '51' '60' '10M' \
        '=' '41' '-20' 'ACGTACGTAC' 'FFFFFFFFFF'
} > "${fil_sam_pe}"

#  Require samtools for indexed FASTA, BAM, and CRAM generation
require_cmd samtools "in '${env_req}' to generate BAM/CRAM fixtures."

#  Index the tiny reference for CRAM generation and reading
samtools faidx "${fil_ref}"

#  Generate sorted and indexed BAM fixtures from SAM provenance
samtools view -bS "${fil_sam_se}" | samtools sort -o "${fil_bam_se}"
samtools index "${fil_bam_se}"

samtools view -bS "${fil_sam_pe}" | samtools sort -o "${fil_bam_pe}"
samtools index "${fil_bam_pe}"

#  Generate indexed CRAM fixtures from the BAM fixtures
samtools view -C -T "${fil_ref}" -o "${fil_cram_se}" "${fil_bam_se}"
samtools index "${fil_cram_se}"

samtools view -C -T "${fil_ref}" -o "${fil_cram_pe}" "${fil_bam_pe}"
samtools index "${fil_cram_pe}"

#  Validate the generated alignment fixtures
samtools quickcheck \
    "${fil_bam_se}" \
    "${fil_bam_pe}" \
    "${fil_cram_se}" \
    "${fil_cram_pe}"


succeed "generated compute-signal fixtures under '${dir_fix}'."
