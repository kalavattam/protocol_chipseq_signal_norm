#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: make.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
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


#  Remove temporary FASTQ and BWA ALN intermediates on exit or failure
function cleanup_tmp_fastqs() {
    rm_files \
        "${dir_fix}" \
        "${fq_se_tmp}" \
        "${fq_r1_tmp}" \
        "${fq_r2_tmp}" \
        "${tmp_bln_se}" \
        "${tmp_bln_r1}" \
        "${tmp_bln_r2}"
}


#  Define fixture directories, reference paths, FASTQ paths, and index prefixes
dir_ref="${dir_fix}/reference"
dir_fq="${dir_fix}/fastq"
dir_fq_se="${dir_fq}/se"
dir_fq_pe="${dir_fq}/pe"
dir_bt2="${dir_fix}/bowtie2"
dir_bwa="${dir_fix}/bwa"
dir_bm2="${dir_fix}/bwa-mem2"

ref="${dir_ref}/tiny.fa"
ref_bwa="${dir_bwa}/tiny.fa"
ref_bm2="${dir_bm2}/tiny.fa"
fq_se_tmp="${dir_fq_se}/tiny_se.atria.fastq.tmp"
fq_r1_tmp="${dir_fq_pe}/tiny_pe_R1.atria.fastq.tmp"
fq_r2_tmp="${dir_fq_pe}/tiny_pe_R2.atria.fastq.tmp"

fq_se_gz="${dir_fq_se}/tiny_se.atria.fastq.gz"
fq_r1_gz="${dir_fq_pe}/tiny_pe_R1.atria.fastq.gz"
fq_r2_gz="${dir_fq_pe}/tiny_pe_R2.atria.fastq.gz"

idx_bt2="${dir_bt2}/tiny"
idx_bwa="${ref_bwa}"
idx_bm2="${ref_bm2}"

tmp_bln_se="${dir_bwa}/tiny_se.sai"
tmp_bln_r1="${dir_bwa}/tiny_pe_R1.sai"
tmp_bln_r2="${dir_bwa}/tiny_pe_R2.sai"

env_req="env_protocol"


#  Register cleanup of temporary FASTQ and BWA ALN intermediates on exit
register_cleanup cleanup_tmp_fastqs

#  Require the project environment for aligner-backed fixtures
require_env "${env_req}" "for align-fastqs fixtures."

#  Require alignment tools used to generate and validate fixtures
require_cmds \
    "in '${env_req}' to generate align-fastqs fixtures." \
    bowtie2 \
    bowtie2-build \
    bowtie2-inspect \
    bwa \
    bwa-mem2 \
    gzip \
    samtools

#  Create fixture output directories
mkdirs \
    "${dir_ref}" \
    "${dir_fq}" \
    "${dir_fq_se}" \
    "${dir_fq_pe}" \
    "${dir_bt2}" \
    "${dir_bwa}" \
    "${dir_bm2}"

#  Remove stale temporary intermediates
cleanup_tmp_fastqs


#  Write tiny reference FASTA shared by all aligners
cat > "${ref}" << EOM
>I
GATCGTACCTAGGCTAACGTTGACCGTTAACGATCGTAGCTAGGATCCGTTACGATCGATGCTAGCTTACCGGATCAAGCTTAGGCTAATCGGCTAAGGTTCCGATTA
EOM

#  Copy the reference into aligner-specific index directories
cp "${ref}" "${ref_bwa}"
cp "${ref}" "${ref_bm2}"


#  Write tiny single-end FASTQ provenance and compressed input fixture
cat > "${fq_se_tmp}" << EOM
@tiny_se_read_1
ACGTTGACCGTTAACGATCGTAGCTAGGAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

#  Compress the single-end FASTQ fixture with deterministic gzip metadata
gzip_n "${fq_se_tmp}" "${fq_se_gz}"

#  Remove the temporary single-end FASTQ after compression
rm_file "${dir_fix}" "${fq_se_tmp}"


#  Write tiny paired-end FASTQ provenance and compressed input fixtures
cat > "${fq_r1_tmp}" << EOM
@tiny_pe_pair_1
ACGTTGACCGTTAACGATCGTAGCTAGGAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

cat > "${fq_r2_tmp}" << EOM
@tiny_pe_pair_1
CCTTAGCCGATTAGCCTAAGCTTGATCCGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

#  Compress paired-end FASTQ fixtures with deterministic gzip metadata
gzip_n "${fq_r1_tmp}" "${fq_r1_gz}"
gzip_n "${fq_r2_tmp}" "${fq_r2_gz}"

#  Remove temporary paired-end FASTQ intermediates after compression
rm_file "${dir_fix}" "${fq_r1_tmp}"
rm_file "${dir_fix}" "${fq_r2_tmp}"


#  Remove stale Bowtie2 index files before regeneration
rm_files \
    "${dir_fix}" \
    "${idx_bt2}.1.bt2" \
    "${idx_bt2}.2.bt2" \
    "${idx_bt2}.3.bt2" \
    "${idx_bt2}.4.bt2" \
    "${idx_bt2}.rev.1.bt2" \
    "${idx_bt2}.rev.2.bt2" \
    "${idx_bt2}.1.bt2l" \
    "${idx_bt2}.2.bt2l" \
    "${idx_bt2}.3.bt2l" \
    "${idx_bt2}.4.bt2l" \
    "${idx_bt2}.rev.1.bt2l" \
    "${idx_bt2}.rev.2.bt2l"

#  Build and validate the Bowtie2 index
log_bt2="${dir_bt2}/bowtie2-build.log"
if ! \
    bowtie2-build "${ref}" "${idx_bt2}" > "${log_bt2}" 2>&1
then
    cat "${log_bt2}" >&2
    exit 1
fi

rm_file "${dir_fix}" "${log_bt2}"

bowtie2-inspect -n "${idx_bt2}" > /dev/null

#  Validate Bowtie2 single-end and paired-end alignment paths
bowtie2 \
    -x "${idx_bt2}" \
    -x "${idx_bt2}" \
    -U "${fq_se_gz}" \
    -S /dev/null \
        > /dev/null \
        2> /dev/null

bowtie2 \
    -x "${idx_bt2}" \
    --very-sensitive \
    --no-mixed \
    --no-discordant \
    --no-overlap \
    --no-dovetail \
    -1 "${fq_r1_gz}" \
    -2 "${fq_r2_gz}" \
    -S /dev/null \
        > /dev/null \
        2> /dev/null


#  Remove stale BWA index files before regeneration
rm_files \
    "${dir_fix}" \
    "${idx_bwa}.amb" \
    "${idx_bwa}.ann" \
    "${idx_bwa}.bwt" \
    "${idx_bwa}.pac" \
    "${idx_bwa}.sa"

#  Build and validate the BWA index
log_bwa="${dir_bwa}/bwa-index.log"
if ! \
    bwa index "${idx_bwa}" > "${log_bwa}" 2>&1
then
    cat "${log_bwa}" >&2
    exit 1
fi

rm_file "${dir_fix}" "${log_bwa}"

#  Validate BWA MEM single-end and paired-end alignment paths
bwa mem "${idx_bwa}" "${fq_se_gz}" \
    > /dev/null 2> /dev/null

bwa mem "${idx_bwa}" "${fq_r1_gz}" "${fq_r2_gz}" \
    > /dev/null 2> /dev/null

#  Validate BWA ALN single-end alignment path
bwa aln -t 1 "${idx_bwa}" "${fq_se_gz}" \
    > "${tmp_bln_se}" 2> /dev/null

bwa samse "${idx_bwa}" "${tmp_bln_se}" "${fq_se_gz}" \
    > /dev/null 2> /dev/null

rm_file "${dir_fix}" "${tmp_bln_se}"

#  Validate BWA ALN paired-end alignment path
bwa aln -t 1 "${idx_bwa}" "${fq_r1_gz}" \
    > "${tmp_bln_r1}" 2> /dev/null

bwa aln -t 1 "${idx_bwa}" "${fq_r2_gz}" \
    > "${tmp_bln_r2}" 2> /dev/null

bwa sampe \
    "${idx_bwa}" "${tmp_bln_r1}" "${tmp_bln_r2}" \
    "${fq_r1_gz}" "${fq_r2_gz}" \
        > /dev/null 2> /dev/null

rm_file "${dir_fix}" "${tmp_bln_r1}"
rm_file "${dir_fix}" "${tmp_bln_r2}"


#  Remove stale BWA-MEM2 index files before regeneration
rm_files \
    "${dir_fix}" \
    "${idx_bm2}.0123" \
    "${idx_bm2}.amb" \
    "${idx_bm2}.ann" \
    "${idx_bm2}.bwt.2bit.64" \
    "${idx_bm2}.pac"

#  Build and validate the BWA-MEM2 index
log_bm2="${dir_bm2}/bwa-mem2-index.log"
if ! \
    bwa-mem2 index "${idx_bm2}" > "${log_bm2}" 2>&1
then
    cat "${log_bm2}" >&2
    exit 1
fi

rm_file "${dir_fix}" "${log_bm2}"

#  Validate BWA-MEM2 single-end and paired-end alignment paths
bwa-mem2 mem "${idx_bm2}" "${fq_se_gz}" \
    > /dev/null 2> /dev/null

bwa-mem2 mem "${idx_bm2}" "${fq_r1_gz}" "${fq_r2_gz}" \
    > /dev/null 2> /dev/null


succeed "generated align-fastqs fixtures under '${dir_fix}'."
