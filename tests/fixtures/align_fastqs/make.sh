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


# Remove temporary FASTQ and BWA ALN intermediates on exit or failure.
function cleanup_tmp_fastqs() {
    rm_files \
        "${dir_fix}" \
        "${fq_se_tmp}" \
        "${fq_r1_tmp}" \
        "${fq_r2_tmp}" \
        "${fq_rpt_se_tmp}" \
        "${fq_rpt_r1_tmp}" \
        "${fq_rpt_r2_tmp}" \
        "${tmp_bln_se}" \
        "${tmp_bln_r1}" \
        "${tmp_bln_r2}"
}


# Print the reverse complement of a sequence.
function print_rc() {
    rev <<< "${1}" | tr ACGT TGCA
}


# Print one FASTQ record with uniform base qualities.
function print_fq() {
    printf '@%s\n%s\n+\n%s\n' "${1}" "${2}" "${2//?/I}"
}


# Define fixture directories, reference paths, FASTQ paths, and index prefixes.
dir_ref="${dir_fix}/reference"
dir_fq="${dir_fix}/fastq"
dir_fq_se="${dir_fq}/se"
dir_fq_pe="${dir_fq}/pe"
dir_bt2="${dir_fix}/bowtie2"
dir_bwa="${dir_fix}/bwa"
dir_bm2="${dir_fix}/bwa-mem2"
dir_sam="${dir_fix}/sam"

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

ref_rpt="${dir_ref}/repeat.fa"
ref_rpt_bwa="${dir_bwa}/repeat.fa"
ref_rpt_bm2="${dir_bm2}/repeat.fa"
fq_rpt_se_tmp="${dir_fq_se}/repeat_se.atria.fastq.tmp"
fq_rpt_r1_tmp="${dir_fq_pe}/repeat_pe_R1.atria.fastq.tmp"
fq_rpt_r2_tmp="${dir_fq_pe}/repeat_pe_R2.atria.fastq.tmp"

fq_rpt_se_gz="${dir_fq_se}/repeat_se.atria.fastq.gz"
fq_rpt_r1_gz="${dir_fq_pe}/repeat_pe_R1.atria.fastq.gz"
fq_rpt_r2_gz="${dir_fq_pe}/repeat_pe_R2.atria.fastq.gz"

idx_rpt_bt2="${dir_bt2}/repeat"
idx_rpt_bwa="${ref_rpt_bwa}"
idx_rpt_bm2="${ref_rpt_bm2}"

sam_sec_pe="${dir_sam}/repeat_secondary_pe.sam"
sam_sec_se="${dir_sam}/repeat_secondary_se.sam"

env_req="env_protocol"


# Register cleanup of temporary FASTQ and BWA ALN intermediates on exit.
register_cleanup cleanup_tmp_fastqs

# Require the project environment for aligner-backed fixtures.
require_env "${env_req}" "for align-fastqs fixtures."

# Require alignment tools used to generate and validate fixtures.
require_cmds \
    "in '${env_req}' to generate align-fastqs fixtures." \
    bowtie2 \
    bowtie2-build \
    bowtie2-inspect \
    bwa \
    bwa-mem2 \
    fold \
    gzip \
    rev \
    samtools \
    tr


# Create fixture output directories.
mkdirs \
    "${dir_ref}" \
    "${dir_fq}" \
    "${dir_fq_se}" \
    "${dir_fq_pe}" \
    "${dir_bt2}" \
    "${dir_bwa}" \
    "${dir_bm2}" \
    "${dir_sam}"

# Remove stale temporary intermediates.
cleanup_tmp_fastqs


# Write tiny reference FASTA shared by all aligners.
cat > "${ref}" << EOM
>I
GATCGTACCTAGGCTAACGTTGACCGTTAACGATCGTAGCTAGGATCCGTTACGATCGATGCTAGCTTACCGGATCAAGCTTAGGCTAATCGGCTAAGGTTCCGATTA
EOM

# Copy the reference into aligner-specific index directories.
cp "${ref}" "${ref_bwa}"
cp "${ref}" "${ref_bm2}"


# Write tiny single-end FASTQ provenance and compressed input fixture.
cat > "${fq_se_tmp}" << EOM
@tiny_se_read_1
ACGTTGACCGTTAACGATCGTAGCTAGGAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

# Compress the single-end FASTQ fixture with deterministic gzip metadata.
gzip_n "${fq_se_tmp}" "${fq_se_gz}"

# Remove the temporary single-end FASTQ after compression.
rm_file "${dir_fix}" "${fq_se_tmp}"


# Write tiny paired-end FASTQ provenance and compressed input fixtures.
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

# Compress paired-end FASTQ fixtures with deterministic gzip metadata.
gzip_n "${fq_r1_tmp}" "${fq_r1_gz}"
gzip_n "${fq_r2_tmp}" "${fq_r2_gz}"

# Remove temporary paired-end FASTQ intermediates after compression.
rm_file "${dir_fix}" "${fq_r1_tmp}"
rm_file "${dir_fix}" "${fq_r2_tmp}"


# Remove stale Bowtie2 index files before regeneration.
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

# Build and validate the Bowtie2 index.
log_bt2="${dir_bt2}/bowtie2-build.log"
if ! \
    bowtie2-build "${ref}" "${idx_bt2}" > "${log_bt2}" 2>&1
then
    cat "${log_bt2}" >&2
    exit 1
fi

rm_file "${dir_fix}" "${log_bt2}"

bowtie2-inspect -n "${idx_bt2}" > /dev/null

# Validate Bowtie2 single-end and paired-end alignment paths.
bowtie2 \
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


# Remove stale BWA index files before regeneration.
rm_files \
    "${dir_fix}" \
    "${idx_bwa}.amb" \
    "${idx_bwa}.ann" \
    "${idx_bwa}.bwt" \
    "${idx_bwa}.pac" \
    "${idx_bwa}.sa"

# Build and validate the BWA index.
log_bwa="${dir_bwa}/bwa-index.log"
if ! \
    bwa index "${idx_bwa}" > "${log_bwa}" 2>&1
then
    cat "${log_bwa}" >&2
    exit 1
fi

rm_file "${dir_fix}" "${log_bwa}"

# Validate BWA MEM single-end and paired-end alignment paths.
bwa mem "${idx_bwa}" "${fq_se_gz}" \
    > /dev/null 2> /dev/null

bwa mem "${idx_bwa}" "${fq_r1_gz}" "${fq_r2_gz}" \
    > /dev/null 2> /dev/null

# Validate BWA ALN single-end alignment path.
bwa aln -t 1 "${idx_bwa}" "${fq_se_gz}" \
    > "${tmp_bln_se}" 2> /dev/null

bwa samse "${idx_bwa}" "${tmp_bln_se}" "${fq_se_gz}" \
    > /dev/null 2> /dev/null

rm_file "${dir_fix}" "${tmp_bln_se}"

# Validate BWA ALN paired-end alignment path.
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


# Remove stale BWA-MEM2 index files before regeneration.
rm_files \
    "${dir_fix}" \
    "${idx_bm2}.0123" \
    "${idx_bm2}.amb" \
    "${idx_bm2}.ann" \
    "${idx_bm2}.bwt.2bit.64" \
    "${idx_bm2}.pac"

# Build and validate the BWA-MEM2 index.
log_bm2="${dir_bm2}/bwa-mem2-index.log"
if ! \
    bwa-mem2 index "${idx_bm2}" > "${log_bm2}" 2>&1
then
    cat "${log_bm2}" >&2
    exit 1
fi

rm_file "${dir_fix}" "${log_bm2}"

# Validate BWA-MEM2 single-end and paired-end alignment paths.
bwa-mem2 mem "${idx_bm2}" "${fq_se_gz}" \
    > /dev/null 2> /dev/null

bwa-mem2 mem "${idx_bm2}" "${fq_r1_gz}" "${fq_r2_gz}" \
    > /dev/null 2> /dev/null


# Write the repeat reference's sequence blocks, 50 bases per line: one 150-base
# repeat and four 400-base unique blocks.
blk_rep="$(tr -d '\n' << 'EOM'
TAGGCCAAACAAGTTGCAACTCCCCTAAATACGATTCCGCATACCAAATT
CATGAGAGGTCTCACATAATCTGCAAATATCCGAATCATATAGATATCGT
CATTTGAGTACTAATGTTATCTATATCTTGACGATATAGGTCACTACTGA
EOM
)"

blk_u_0="$(tr -d '\n' << 'EOM'
GTCTAGCTCTTCCGTAGTGATGTACGCCCCTGCTCCTTAATATCAGTAAC
CACCGTCCTAACCAATATCCATCGTGGCTCCGTCCAGACACCTCGCAACA
GATTTATATGCGAGACATAATAGCATAGCTAGACTTGGGCAGAGACACAC
TAGGACTGGCGTGCAGTCAATTTCTCAGTGCTGAACCGAGAGTAGCAAAC
ACCGAAACTCCAGCAGCATTCCCATAAAAGGGGGTTACCTCGGCAAGTTG
TAGCTAGAGATAACACCATCTGCTTTAGTGTCGCTGTTTCCTTCCAAACG
CGCTGCTTTTCCATCGTTCGACAAGACAGCCACTTACCCTTATCTCTTAT
TTCCTTATACCCCCGGTATAGAGAATCAAGCGACATATCGTATTAAAACC
EOM
)"

blk_u_1="$(tr -d '\n' << 'EOM'
CGAGAATAATTCGCGACCAAAATCGAGGGGATGCTTCCAAGCGGCCCCGA
CGTGCCCCATATTCGTGCCTAATCATAACACCGGGACCTGGTCTTTGCGC
GTTTGGAATGATCTGCGCGGCTTATATCTAGACTGACATTTGTAAGTGCT
ATAGTTACCGCAACCGTAATAACTCCCACCTGAAAAGGTCGTCTGCGATT
GTACGCACAGCTTTTTGAGAAGCTGCACCGTGTACGACAAACATCCGTGG
AGCGACCCCCCACATTTATAGCGATGGTGCATGTTCCCTCTGGTTTGATG
TAAATAATTGTGGGGCAGTTACTCGGATCGATGCACATGGAGCCACAACC
AGGTTTAGCGGGCCCGTGGCCTCACAATCATTACAAGTCGCGTTGGATAT
EOM
)"

blk_u_2="$(tr -d '\n' << 'EOM'
ACACCCAGGAATCAGGCAGACCTGCCGGTTAGCCCCGCGTGATTTCTATG
TGCTCTAAGACTCAGGGGTTTCGTTAACGCCCGCGTCGATCTTCAGCATT
ACTATTGAATCTAGCATCTGGTCGCTAGTTGGATGTCCCGGCGTATCTTG
CACCCTGATCCCGATGTTCCAAAATTTCGGCATGGAGATCTTCAAGCGCC
ACGCTACACACAGAGTTGGCCTAAGATGATCTACAAGCACAAACCGCGCT
TTCAATTGCAGATAGAGTGAAGAGACCGGCGTATAATTCGGGCATTGGAG
ACGCAGACCGTGACTATGGTTCGGCTCGGGCACGACAAGGAGCTGTCCCA
CTCCGGAGTCAACCGTTTTAGGTTCGAAATGGAACAATGGCACACCTGTA
EOM
)"

blk_u_3="$(tr -d '\n' << 'EOM'
TGGTTGTGTCCTAAGTACCCGGAAAAGTATGACTTCAAAAACGGGGTAGA
GGAATGTTGAAAGTTGCTTTGTCGCTTTATCCTGCGGGCGGCATAGTAAT
GAGTAGTTGCGGTCGAAGGTATACCTCTGCAATGAAGGGAGTGCGTCACC
GCCACGGACTCTCGTATAGCGGGCCATGTCTGTTACACTAAGTGAATCTA
CAAGAACTAACATATGCGACGTTTTGGAGTTCAAGCTGTAATTTTACCTC
TGATAAGTCCATGGTGTAAGGTACTGGGGTCCACGCACGGAAGGTAAAGC
TGCCGGAGGATCGCAGTTCAAAATTTGTTTCGTGAGTCTAGACCATCGAA
ACTGAGCGAGCATGCGGAGGCGGCAACGTAAGCCCGCGGGTATGCGAGTG
EOM
)"

# Write the repeat reference: the unique blocks separated by three copies of
# the repeat, so a read in the repeat has three equally good placements.
seq_rpt="${blk_u_0}${blk_rep}${blk_u_1}${blk_rep}"
seq_rpt+="${blk_u_2}${blk_rep}${blk_u_3}"

{
    echo ">chrR"
    fold -w 60 <<< "${seq_rpt}"
} > "${ref_rpt}"

cp "${ref_rpt}" "${ref_rpt_bwa}"
cp "${ref_rpt}" "${ref_rpt_bm2}"


# Write the repeat reads, 100 bases each, cut from the blocks (offsets are
# zero-based):
# - 'chim': R1 is 60 bases of repeat then 40 bases of unique block 3, a split
#   read whose primary lies in the repeat; R2 is unique.
# - 'rchm': R1 is 60 bases of unique block 3 then 40 bases of repeat, a split
#   read whose primary is unique; R2 is unique.
# - 'half': R1 lies wholly in the repeat; R2 is unique.
# - 'ctrl': both mates unique.
# Single-end input holds the R1 of each read but 'half'.
{
    print_fq chim "${blk_rep:20:60}${blk_u_3:100:40}"
    print_fq rchm "${blk_u_3:100:60}${blk_rep:20:40}"
    print_fq half "${blk_rep:25:100}"
    print_fq ctrl "${blk_u_1:150:100}"
} > "${fq_rpt_r1_tmp}"

{
    print_fq chim "$(print_rc "${blk_u_1:100:100}")"
    print_fq rchm "$(print_rc "${blk_u_3:250:100}")"
    print_fq half "$(print_rc "${blk_u_1:50:100}")"
    print_fq ctrl "$(print_rc "${blk_u_1:280:100}")"
} > "${fq_rpt_r2_tmp}"

{
    print_fq chim "${blk_rep:20:60}${blk_u_3:100:40}"
    print_fq rchm "${blk_u_3:100:60}${blk_rep:20:40}"
    print_fq ctrl "${blk_u_1:150:100}"
} > "${fq_rpt_se_tmp}"

# Compress the repeat FASTQ fixtures with deterministic gzip metadata, then
# remove the temporary FASTQs.
gzip_n "${fq_rpt_se_tmp}" "${fq_rpt_se_gz}"
gzip_n "${fq_rpt_r1_tmp}" "${fq_rpt_r1_gz}"
gzip_n "${fq_rpt_r2_tmp}" "${fq_rpt_r2_gz}"

rm_file "${dir_fix}" "${fq_rpt_se_tmp}"
rm_file "${dir_fix}" "${fq_rpt_r1_tmp}"
rm_file "${dir_fix}" "${fq_rpt_r2_tmp}"


# Remove stale repeat index files before regeneration.
rm_files \
    "${dir_fix}" \
    "${idx_rpt_bt2}.1.bt2" \
    "${idx_rpt_bt2}.2.bt2" \
    "${idx_rpt_bt2}.3.bt2" \
    "${idx_rpt_bt2}.4.bt2" \
    "${idx_rpt_bt2}.rev.1.bt2" \
    "${idx_rpt_bt2}.rev.2.bt2" \
    "${idx_rpt_bwa}.amb" \
    "${idx_rpt_bwa}.ann" \
    "${idx_rpt_bwa}.bwt" \
    "${idx_rpt_bwa}.pac" \
    "${idx_rpt_bwa}.sa" \
    "${idx_rpt_bm2}.0123" \
    "${idx_rpt_bm2}.amb" \
    "${idx_rpt_bm2}.ann" \
    "${idx_rpt_bm2}.bwt.2bit.64" \
    "${idx_rpt_bm2}.pac"

# Build the repeat indexes for Bowtie2, BWA, and BWA-MEM2.
log_rpt="${dir_ref}/repeat-index.log"
if ! \
    {
        bowtie2-build "${ref_rpt}" "${idx_rpt_bt2}" \
            && bwa index "${idx_rpt_bwa}" \
            && bwa-mem2 index "${idx_rpt_bm2}"
    } > "${log_rpt}" 2>&1
then
    cat "${log_rpt}" >&2
    exit 1
fi

rm_file "${dir_fix}" "${log_rpt}"

# Validate the repeat paired-end alignment path for each aligner.
bowtie2 \
    -x "${idx_rpt_bt2}" \
    -1 "${fq_rpt_r1_gz}" \
    -2 "${fq_rpt_r2_gz}" \
    -S /dev/null \
        > /dev/null \
        2> /dev/null

bwa mem "${idx_rpt_bwa}" "${fq_rpt_r1_gz}" "${fq_rpt_r2_gz}" \
    > /dev/null 2> /dev/null

bwa-mem2 mem "${idx_rpt_bm2}" "${fq_rpt_r1_gz}" "${fq_rpt_r2_gz}" \
    > /dev/null 2> /dev/null


# Write name-grouped SAM records with secondary records for the MAPQ template
# filter: at MAPQ 30, 'sec_keep' has passing primaries and a failing secondary,
# and 'sec_drop' a failing primary and a passing secondary.
{
    write_sam_line '@HD' 'VN:1.6' 'SO:queryname'
    write_sam_line '@SQ' 'SN:chrR' "LN:${#seq_rpt}"

    write_sam_line \
        'sec_drop' '99' 'chrR' '1' '0' '10M' \
        '=' '101' '110' 'ACGTACGTAA' 'FFFFFFFFFF'

    write_sam_line \
        'sec_drop' '147' 'chrR' '101' '60' '10M' \
        '=' '1' '-110' 'ACGTACGTAA' 'FFFFFFFFFF'

    write_sam_line \
        'sec_drop' '355' 'chrR' '501' '255' '10M' \
        '=' '101' '0' 'ACGTACGTAA' 'FFFFFFFFFF'

    write_sam_line \
        'sec_keep' '99' 'chrR' '201' '60' '10M' \
        '=' '301' '110' 'ACGTACGTAA' 'FFFFFFFFFF'

    write_sam_line \
        'sec_keep' '147' 'chrR' '301' '60' '10M' \
        '=' '201' '-110' 'ACGTACGTAA' 'FFFFFFFFFF'

    write_sam_line \
        'sec_keep' '355' 'chrR' '601' '0' '10M' \
        '=' '301' '0' 'ACGTACGTAA' 'FFFFFFFFFF'
} > "${sam_sec_pe}"

{
    write_sam_line '@HD' 'VN:1.6' 'SO:queryname'
    write_sam_line '@SQ' 'SN:chrR' "LN:${#seq_rpt}"

    write_sam_line \
        'sec_drop' '0' 'chrR' '1' '0' '10M' \
        '*' '0' '0' 'ACGTACGTAA' 'FFFFFFFFFF'

    write_sam_line \
        'sec_drop' '256' 'chrR' '501' '255' '10M' \
        '*' '0' '0' 'ACGTACGTAA' 'FFFFFFFFFF'

    write_sam_line \
        'sec_keep' '0' 'chrR' '201' '60' '10M' \
        '*' '0' '0' 'ACGTACGTAA' 'FFFFFFFFFF'

    write_sam_line \
        'sec_keep' '256' 'chrR' '601' '0' '10M' \
        '*' '0' '0' 'ACGTACGTAA' 'FFFFFFFFFF'
} > "${sam_sec_se}"

# Validate the secondary-record SAM fixtures.
samtools quickcheck "${sam_sec_pe}" "${sam_sec_se}"


succeed "generated align-fastqs fixtures under '${dir_fix}'."
