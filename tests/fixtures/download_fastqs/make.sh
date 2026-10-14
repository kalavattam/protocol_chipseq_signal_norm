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
# - Anthropic Claude Code (Opus 5.5).
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


# Remove temporary FASTQ intermediates on exit or failure.
function cleanup_fastqs_tmp() {
    rm_files \
        "${dir_fix}" \
        "${fq_se_tmp}" \
        "${fq_r1_tmp}" \
        "${fq_r2_tmp}" \
        "${fq_nam_tmp}"
}


# Write one gzipped FASTQ of 8-base records whose headers are the arguments.
function write_fastq_names() {
    local fil_out="${1}"
    local hdr
    shift 1

    for hdr in "$@"; do
        printf '%s\n' "${hdr}" ACGTACGT + IIIIIIII
    done > "${fq_nam_tmp}"

    gzip_n "${fq_nam_tmp}" "${fil_out}"
    rm_file "${dir_fix}" "${fq_nam_tmp}"
}


# Define fixture directories, FASTQ paths, and metadata templates.
dir_src="${dir_fix}/source"
dir_src_se="${dir_src}/se"
dir_src_pe="${dir_src}/pe"
dir_mtd="${dir_fix}/metadata"
dir_src_nam="${dir_src}/names"
dir_exp="${dir_fix}/expected"

fq_se_tmp="${dir_src_se}/tiny_download_se.fastq.tmp"
fq_r1_tmp="${dir_src_pe}/tiny_download_pe_R1.fastq.tmp"
fq_r2_tmp="${dir_src_pe}/tiny_download_pe_R2.fastq.tmp"
fq_nam_tmp="${dir_src_nam}/names.fastq.tmp"

fq_se_gz="${dir_src_se}/tiny_download_se.fastq.gz"
fq_r1_gz="${dir_src_pe}/tiny_download_pe_R1.fastq.gz"
fq_r2_gz="${dir_src_pe}/tiny_download_pe_R2.fastq.gz"
fq_se_trn="${dir_src_se}/tiny_download_se_truncated.fastq.gz"

tpl_se="${dir_mtd}/local_se.template.tsv"
tpl_pe="${dir_mtd}/local_pe.template.tsv"
tpl_pe_dup="${dir_mtd}/local_pe_duplicate.template.tsv"
tpl_pe_dup_cus="${dir_mtd}/local_pe_duplicate_custom.template.tsv"
tpl_pe_cnfl="${dir_mtd}/local_pe_conflicting_accession.template.tsv"
tpl_mix="${dir_mtd}/local_mixed.template.tsv"
tpl_se_vrf="${dir_mtd}/local_se_verified.template.tsv"
tpl_pe_vrf="${dir_mtd}/local_pe_verified.template.tsv"
tpl_pe_md5="${dir_mtd}/local_pe_md5_only.template.tsv"
tpl_se_bad_md5="${dir_mtd}/local_se_wrong_md5.template.tsv"
tpl_se_bad_siz="${dir_mtd}/local_se_wrong_size.template.tsv"
tpl_pe_bad_cnt="${dir_mtd}/local_pe_md5_count.template.tsv"
tpl_pe_bad_r2="${dir_mtd}/local_pe_wrong_md5_r2.template.tsv"
tpl_pe_cnfl_md5="${dir_mtd}/local_pe_conflicting_md5.template.tsv"
tpl_se_trn="${dir_mtd}/local_se_truncated.template.tsv"
tpl_se_404="${dir_mtd}/local_se_missing.template.tsv"
tpl_nam_pe="${dir_mtd}/local_names_pe.template.tsv"
tpl_nam_pe_r2="${dir_mtd}/local_names_pe_r2_name.template.tsv"
tpl_nam_pe_ord="${dir_mtd}/local_names_pe_order.template.tsv"
tpl_nam_pe_mat="${dir_mtd}/local_names_pe_mate.template.tsv"
tpl_nam_pe_404="${dir_mtd}/local_names_pe_missing.template.tsv"
tpl_nam_se="${dir_mtd}/local_names_se.template.tsv"
tpl_nam_se_nmt="${dir_mtd}/local_names_se_no_mate.template.tsv"
tpl_nam_se_lst="${dir_mtd}/local_names_se_lost.template.tsv"
tpl_nam_seq="${dir_mtd}/local_names_sequencer_pe.template.tsv"
tpl_nam_seq_se="${dir_mtd}/local_names_sequencer_se.template.tsv"


# Register cleanup of temporary FASTQ intermediates on exit.
register_cleanup cleanup_fastqs_tmp

# Require gzip for deterministic compressed workflow inputs.
require_cmd gzip "to generate download-fastqs fixtures."
require_cmd md5sum "to write expected MD5 values into metadata templates."

# Create fixture output directories.
mkdirs \
    "${dir_src}" \
    "${dir_src_se}" \
    "${dir_src_pe}" \
    "${dir_src_nam}" \
    "${dir_mtd}" \
    "${dir_exp}"

# Remove stale temporary intermediates.
cleanup_fastqs_tmp

# Write tiny single-end FASTQ provenance and compressed source fixture.
cat > "${fq_se_tmp}" << EOM
@tiny_download_se_read_1
GATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

# Compress single-end FASTQ fixture with deterministic gzip metadata.
gzip_n "${fq_se_tmp}" "${fq_se_gz}"

# Remove temporary single-end FASTQ intermediate after compression.
rm_file "${dir_fix}" "${fq_se_tmp}"

# Write tiny paired-end FASTQ provenance and compressed source fixtures.
cat > "${fq_r1_tmp}" << EOM
@tiny_download_pe_pair_1/1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

cat > "${fq_r2_tmp}" << EOM
@tiny_download_pe_pair_1/2
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

# Compress paired-end FASTQ fixtures with deterministic gzip metadata.
gzip_n "${fq_r1_tmp}" "${fq_r1_gz}"
gzip_n "${fq_r2_tmp}" "${fq_r2_gz}"

# Remove temporary paired-end FASTQ intermediates after compression.
rm_file "${dir_fix}" "${fq_r1_tmp}"
rm_file "${dir_fix}" "${fq_r2_tmp}"

# Write metadata templates with placeholders for runtime HTTP URLs; the smoke
# test replaces '__BASE_URL__' with a loopback server URL.
{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https'

    write_tsv_row \
        'SRR_LOCAL_SE' \
        'tiny_download_se' \
        '__BASE_URL__/se/tiny_download_se.fastq.gz'
} > "${tpl_se}"

{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https'

    write_tsv_row \
        'SRR_LOCAL_PE' \
        'tiny_download_pe' \
        '__BASE_URL__/pe/tiny_download_pe_R1.fastq.gz;__BASE_URL__/pe/tiny_download_pe_R2.fastq.gz'
} > "${tpl_pe}"

{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https'

    write_tsv_row \
        'SRR_LOCAL_PE' \
        'tiny_download_pe_h3k18' \
        '__BASE_URL__/pe/tiny_download_pe_R1.fastq.gz;__BASE_URL__/pe/tiny_download_pe_R2.fastq.gz'

    write_tsv_row \
        'SRR_LOCAL_PE' \
        'tiny_download_pe_h3k27' \
        '__BASE_URL__/pe/tiny_download_pe_R1.fastq.gz;__BASE_URL__/pe/tiny_download_pe_R2.fastq.gz'
} > "${tpl_pe_dup}"

{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https'

    write_tsv_row \
        'SRR_LOCAL_PE_DUP_CUS_A' \
        'tiny_download_pe_duplicate' \
        '__BASE_URL__/pe/tiny_download_pe_R1.fastq.gz;__BASE_URL__/pe/tiny_download_pe_R2.fastq.gz'

    write_tsv_row \
        'SRR_LOCAL_PE_DUP_CUS_B' \
        'tiny_download_pe_duplicate' \
        '__BASE_URL__/pe/tiny_download_pe_R1.fastq.gz;__BASE_URL__/pe/tiny_download_pe_R2.fastq.gz'
} > "${tpl_pe_dup_cus}"

{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https'

    write_tsv_row \
        'SRR_LOCAL_PE_CONFLICT' \
        'tiny_download_pe_conflict_a' \
        '__BASE_URL__/pe/tiny_download_pe_R1.fastq.gz;__BASE_URL__/pe/tiny_download_pe_R2.fastq.gz'

    write_tsv_row \
        'SRR_LOCAL_PE_CONFLICT' \
        'tiny_download_pe_conflict_b' \
        '__BASE_URL__/pe/tiny_download_pe_R1.fastq.gz;__BASE_URL__/pe/missing_R2.fastq.gz'
} > "${tpl_pe_cnfl}"

{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https'

    write_tsv_row \
        'SRR_LOCAL_SE' \
        'tiny_download_se' \
        '__BASE_URL__/se/tiny_download_se.fastq.gz'

    write_tsv_row \
        'SRR_LOCAL_PE' \
        'tiny_download_pe' \
        '__BASE_URL__/pe/tiny_download_pe_R1.fastq.gz;__BASE_URL__/pe/tiny_download_pe_R2.fastq.gz'
} > "${tpl_mix}"


# Write a single-end source whose gzip trailer (CRC and size, 8 bytes) is cut
# off, so it fails 'gzip -t'.
siz_se="$(wc -c < "${fq_se_gz}")"
siz_se="${siz_se//[[:space:]]/}"
head -c "$(( siz_se - 8 ))" "${fq_se_gz}" > "${fq_se_trn}"

# Record each source file's MD5 and size for the verification templates.
md5_se="$(md5sum < "${fq_se_gz}")"
md5_se="${md5_se%% *}"
md5_r1="$(md5sum < "${fq_r1_gz}")"
md5_r1="${md5_r1%% *}"
md5_r2="$(md5sum < "${fq_r2_gz}")"
md5_r2="${md5_r2%% *}"
siz_r1="$(wc -c < "${fq_r1_gz}")"
siz_r1="${siz_r1//[[:space:]]/}"
siz_r2="$(wc -c < "${fq_r2_gz}")"
siz_r2="${siz_r2//[[:space:]]/}"

url_se='__BASE_URL__/se/tiny_download_se.fastq.gz'
url_pe='__BASE_URL__/pe/tiny_download_pe_R1.fastq.gz'
url_pe+=';__BASE_URL__/pe/tiny_download_pe_R2.fastq.gz'

# Write templates with expected MD5 values and sizes: right values, an MD5 sans
# sizes, a wrong MD5, a wrong size, one MD5 for two URLs, a wrong R2 MD5, and
# one accession listed twice with different MD5 values.
{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https' \
        'fastq_md5' \
        'fastq_bytes'

    write_tsv_row \
        'SRR_LOCAL_SE' \
        'tiny_download_se' \
        "${url_se}" \
        "${md5_se}" \
        "${siz_se}"
} > "${tpl_se_vrf}"

{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https' \
        'fastq_md5' \
        'fastq_bytes'

    write_tsv_row \
        'SRR_LOCAL_PE' \
        'tiny_download_pe' \
        "${url_pe}" \
        "${md5_r1};${md5_r2}" \
        "${siz_r1};${siz_r2}"
} > "${tpl_pe_vrf}"

{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https' \
        'fastq_md5'

    write_tsv_row \
        'SRR_LOCAL_PE' \
        'tiny_download_pe' \
        "${url_pe}" \
        "${md5_r1};${md5_r2}"
} > "${tpl_pe_md5}"

{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https' \
        'fastq_md5' \
        'fastq_bytes'

    write_tsv_row \
        'SRR_LOCAL_SE' \
        'tiny_download_se' \
        "${url_se}" \
        '0123456789abcdef0123456789abcdef' \
        "${siz_se}"
} > "${tpl_se_bad_md5}"

{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https' \
        'fastq_md5' \
        'fastq_bytes'

    write_tsv_row \
        'SRR_LOCAL_SE' \
        'tiny_download_se' \
        "${url_se}" \
        "${md5_se}" \
        "$(( siz_se + 1 ))"
} > "${tpl_se_bad_siz}"

{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https' \
        'fastq_md5'

    write_tsv_row \
        'SRR_LOCAL_PE' \
        'tiny_download_pe' \
        "${url_pe}" \
        "${md5_r1}"
} > "${tpl_pe_bad_cnt}"

{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https' \
        'fastq_md5'

    write_tsv_row \
        'SRR_LOCAL_PE' \
        'tiny_download_pe' \
        "${url_pe}" \
        "${md5_r1};0123456789abcdef0123456789abcdef"
} > "${tpl_pe_bad_r2}"

{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https' \
        'fastq_md5'

    write_tsv_row \
        'SRR_LOCAL_PE' \
        'tiny_download_pe_a' \
        "${url_pe}" \
        "${md5_r1};${md5_r2}"

    write_tsv_row \
        'SRR_LOCAL_PE' \
        'tiny_download_pe_b' \
        "${url_pe}" \
        "${md5_r1};0123456789abcdef0123456789abcdef"
} > "${tpl_pe_cnfl_md5}"

# Write templates without MD5 values: a truncated source and a URL the server
# does not have.
{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https'

    write_tsv_row \
        'SRR_LOCAL_SE' \
        'tiny_download_se' \
        '__BASE_URL__/se/tiny_download_se_truncated.fastq.gz'
} > "${tpl_se_trn}"

{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https'

    write_tsv_row \
        'SRR_LOCAL_SE' \
        'tiny_download_se' \
        '__BASE_URL__/se/missing.fastq.gz'
} > "${tpl_se_404}"


# Write paired-end sources in ENA's form, '@<run>.<n> <Illumina name>/<mate>',
# plus three bad pairs that each fail one check: an R2 name unlike R1's, a
# skipped spot number, '/2' in R1. Also write single-end sources with '/1',
# without it, and with names the archive did not store.
pfx_nam='VH00699:369:AACWKJJM5:1:1101'
nam_1="${pfx_nam}:30009:1000"
nam_2="${pfx_nam}:30236:1000"

write_fastq_names \
    "${dir_src_nam}/ena_R1.fastq.gz" \
    "@SRR_NAMES_PE.1 ${nam_1}/1" \
    "@SRR_NAMES_PE.2 ${nam_2}/1"

write_fastq_names \
    "${dir_src_nam}/ena_R2.fastq.gz" \
    "@SRR_NAMES_PE.1 ${nam_1}/2" \
    "@SRR_NAMES_PE.2 ${nam_2}/2"

write_fastq_names \
    "${dir_src_nam}/r2_name_R2.fastq.gz" \
    "@SRR_NAMES_PE.1 ${nam_1}/2" \
    "@SRR_NAMES_PE.2 ${pfx_nam}:30237:1000/2"

write_fastq_names \
    "${dir_src_nam}/order_R1.fastq.gz" \
    "@SRR_NAMES_PE.1 ${nam_1}/1" \
    "@SRR_NAMES_PE.3 ${nam_2}/1"

write_fastq_names \
    "${dir_src_nam}/order_R2.fastq.gz" \
    "@SRR_NAMES_PE.1 ${nam_1}/2" \
    "@SRR_NAMES_PE.3 ${nam_2}/2"

write_fastq_names \
    "${dir_src_nam}/mate_R1.fastq.gz" \
    "@SRR_NAMES_PE.1 ${nam_1}/2" \
    "@SRR_NAMES_PE.2 ${nam_2}/1"

write_fastq_names \
    "${dir_src_nam}/ena_se.fastq.gz" \
    "@SRR_NAMES_SE.1 ${nam_1}/1" \
    "@SRR_NAMES_SE.2 ${nam_2}/1"

write_fastq_names \
    "${dir_src_nam}/ena_se_no_mate.fastq.gz" \
    "@SRR_NAMES_SE.1 ${nam_1}" \
    "@SRR_NAMES_SE.2 ${nam_2}"

write_fastq_names \
    "${dir_src_nam}/lost_se.fastq.gz" \
    "@SRR_NAMES_LOST.1 1/1" \
    "@SRR_NAMES_LOST.2 2/1"

# Write sources as they come from the sequencer, whose names need nothing: a
# pair with the bcl-convert comment and a single-end file without a comment.
write_fastq_names \
    "${dir_src_nam}/sequencer_R1.fastq.gz" \
    "@${nam_1} 1:N:0:ACGTACGT" \
    "@${nam_2} 1:N:0:ACGTACGT"

write_fastq_names \
    "${dir_src_nam}/sequencer_R2.fastq.gz" \
    "@${nam_1} 2:N:0:ACGTACGT" \
    "@${nam_2} 2:N:0:ACGTACGT"

write_fastq_names \
    "${dir_src_nam}/sequencer_se.fastq.gz" \
    "@${nam_1}" \
    "@${nam_2}"

# Write the expected restored records, uncompressed; both mates of the pair
# restore to the same headers.
for nam in pe se; do
    if [[ "${nam}" == "pe" ]]; then
        acc=SRR_NAMES_PE
    else
        acc=SRR_NAMES_SE
    fi

    printf '%s\n' \
        "@${nam_1} ${acc}.1" ACGTACGT + IIIIIIII \
        "@${nam_2} ${acc}.2" ACGTACGT + IIIIIIII \
        > "${dir_exp}/names_${nam}.fastq"
done

# Write the read-name templates, one source pair or file per template; the
# missing-file template has the right MD5 values and sizes but URLs the server
# does not have, so only files already in place can satisfy it.
md5_nam_r1="$(md5sum < "${dir_src_nam}/ena_R1.fastq.gz")"
md5_nam_r1="${md5_nam_r1%% *}"
md5_nam_r2="$(md5sum < "${dir_src_nam}/ena_R2.fastq.gz")"
md5_nam_r2="${md5_nam_r2%% *}"
siz_nam_r1="$(wc -c < "${dir_src_nam}/ena_R1.fastq.gz")"
siz_nam_r1="${siz_nam_r1//[[:space:]]/}"
siz_nam_r2="$(wc -c < "${dir_src_nam}/ena_R2.fastq.gz")"
siz_nam_r2="${siz_nam_r2//[[:space:]]/}"

url_nam='__BASE_URL__/names'

arr_rows=(
    "${tpl_nam_pe}|SRR_NAMES_PE|ena_R1;ena_R2"
    "${tpl_nam_pe_r2}|SRR_NAMES_PE|ena_R1;r2_name_R2"
    "${tpl_nam_pe_ord}|SRR_NAMES_PE|order_R1;order_R2"
    "${tpl_nam_pe_mat}|SRR_NAMES_PE|mate_R1;ena_R2"
    "${tpl_nam_se}|SRR_NAMES_SE|ena_se"
    "${tpl_nam_se_nmt}|SRR_NAMES_SE|ena_se_no_mate"
    "${tpl_nam_se_lst}|SRR_NAMES_LOST|lost_se"
    "${tpl_nam_seq_se}|SEQ_SE|sequencer_se"
)

for row in "${arr_rows[@]}"; do
    IFS='|' read -r tpl acc stems <<< "${row}"
    url="${url_nam}/${stems%%;*}.fastq.gz"

    if [[ "${stems}" == *";"* ]]; then
        url+=";${url_nam}/${stems#*;}.fastq.gz"
    fi

    {
        write_tsv_row \
            'run_accession' \
            'custom_name' \
            'fastq_https'

        write_tsv_row \
            "${acc}" \
            "names_${acc,,}" \
            "${url}"
    } > "${tpl}"
done

{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https' \
        'fastq_md5' \
        'fastq_bytes'

    write_tsv_row \
        'SRR_NAMES_PE' \
        'names_srr_names_pe' \
        "${url_nam}/missing_R1.fastq.gz;${url_nam}/missing_R2.fastq.gz" \
        "${md5_nam_r1};${md5_nam_r2}" \
        "${siz_nam_r1};${siz_nam_r2}"
} > "${tpl_nam_pe_404}"

md5_seq_r1="$(md5sum < "${dir_src_nam}/sequencer_R1.fastq.gz")"
md5_seq_r1="${md5_seq_r1%% *}"
md5_seq_r2="$(md5sum < "${dir_src_nam}/sequencer_R2.fastq.gz")"
md5_seq_r2="${md5_seq_r2%% *}"

{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https' \
        'fastq_md5'

    write_tsv_row \
        'SEQ_PE' \
        'names_seq_pe' \
        "${url_nam}/sequencer_R1.fastq.gz;${url_nam}/sequencer_R2.fastq.gz" \
        "${md5_seq_r1};${md5_seq_r2}"
} > "${tpl_nam_seq}"


succeed "generated download-fastqs fixtures under ${dir_fix}"
