#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: tests/download_fastqs/scripts/make_fixtures.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
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


#  Resolve paths relative to 'tests/download_fastqs/scripts'
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_dwn="$(cd "${dir_scr}/.." > /dev/null 2>&1 && pwd)"
dir_fix="${dir_dwn}/fixtures"

#  Source shared fixture-generation helpers
# shellcheck disable=SC1091
source "${dir_scr}/../../scripts/lib/fixture_helpers.sh"


#  Define fixture directories, FASTQ paths, and metadata templates
dir_src="${dir_fix}/source"
dir_mtd="${dir_fix}/metadata"

fq_se_tmp="${dir_src}/tiny_download_se.fastq.tmp"
fq_r1_tmp="${dir_src}/tiny_download_pe_R1.fastq.tmp"
fq_r2_tmp="${dir_src}/tiny_download_pe_R2.fastq.tmp"

fq_se_gz="${dir_src}/tiny_download_se.fastq.gz"
fq_r1_gz="${dir_src}/tiny_download_pe_R1.fastq.gz"
fq_r2_gz="${dir_src}/tiny_download_pe_R2.fastq.gz"

tpl_se="${dir_mtd}/local_se.template.tsv"
tpl_pe="${dir_mtd}/local_pe.template.tsv"
tpl_mix="${dir_mtd}/local_mixed.template.tsv"


#  Remove temporary FASTQ intermediates on exit or failure
function cleanup_tmp_fastqs() {
    rm_files \
        "${dir_fix}" \
        "${fq_se_tmp}" \
        "${fq_r1_tmp}" \
        "${fq_r2_tmp}"
}


#  Register cleanup of temporary FASTQ intermediates on exit
register_cleanup cleanup_tmp_fastqs

#  Require gzip for deterministic compressed workflow inputs
require_cmd gzip "to generate download-fastqs fixtures."

#  Create fixture output directories
mkdirs "${dir_src}" "${dir_mtd}"

#  Remove stale temporary intermediates
cleanup_tmp_fastqs

#  Write tiny single-end FASTQ provenance and compressed source fixture
cat > "${fq_se_tmp}" << EOM
@tiny_download_se_read_1
GATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

#  Compress single-end FASTQ fixture with deterministic gzip metadata
gzip_n "${fq_se_tmp}" "${fq_se_gz}"

#  Remove temporary single-end FASTQ intermediate after compression
rm_file "${dir_fix}" "${fq_se_tmp}"

#  Write tiny paired-end FASTQ provenance and compressed source fixtures
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

#  Compress paired-end FASTQ fixtures with deterministic gzip metadata
gzip_n "${fq_r1_tmp}" "${fq_r1_gz}"
gzip_n "${fq_r2_tmp}" "${fq_r2_gz}"

#  Remove temporary paired-end FASTQ intermediates after compression
rm_file "${dir_fix}" "${fq_r1_tmp}"
rm_file "${dir_fix}" "${fq_r2_tmp}"

#  Write metadata templates with placeholders for runtime HTTP URLs; the smoke
#+ test replaces '__BASE_URL__' with a loopback server URL
{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https'

    write_tsv_row \
        'SRR_LOCAL_SE' \
        'tiny_download_se' \
        '__BASE_URL__/tiny_download_se.fastq.gz'
} > "${tpl_se}"

{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https'

    write_tsv_row \
        'SRR_LOCAL_PE' \
        'tiny_download_pe' \
        '__BASE_URL__/tiny_download_pe_R1.fastq.gz;__BASE_URL__/tiny_download_pe_R2.fastq.gz'
} > "${tpl_pe}"

{
    write_tsv_row \
        'run_accession' \
        'custom_name' \
        'fastq_https'

    write_tsv_row \
        'SRR_LOCAL_SE' \
        'tiny_download_se' \
        '__BASE_URL__/tiny_download_se.fastq.gz'

    write_tsv_row \
        'SRR_LOCAL_PE' \
        'tiny_download_pe' \
        '__BASE_URL__/tiny_download_pe_R1.fastq.gz;__BASE_URL__/tiny_download_pe_R2.fastq.gz'
} > "${tpl_mix}"


succeed "generated download-fastqs fixtures under ${dir_fix}"
