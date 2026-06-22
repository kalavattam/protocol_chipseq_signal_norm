#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: tests/trim_fastqs/scripts/make_fixtures.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development and
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


#  Resolve paths relative to 'tests/trim_fastqs/scripts'
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_trm="$(cd "${dir_scr}/.." > /dev/null 2>&1 && pwd)"
dir_fix="${dir_trm}/fixtures"

#  Source shared fixture-generation helpers
# shellcheck disable=SC1091
source "${dir_scr}/../../scripts/lib/fixture_helpers.sh"


#  Remove temporary FASTQ intermediates on exit or failure
function cleanup_tmp_fastqs() {
    rm_files \
        "${dir_fix}" \
        "${fq_se_tmp}" \
        "${fq_r1_tmp}" \
        "${fq_r2_tmp}"
}


#  Define fixture directories and FASTQ paths
dir_fq="${dir_fix}/fastq"
dir_fq_se="${dir_fq}/se"
dir_fq_pe="${dir_fq}/pe"

fq_se_tmp="${dir_fq_se}/tiny_se.fastq.tmp"
fq_r1_tmp="${dir_fq_pe}/tiny_pe_R1.fastq.tmp"
fq_r2_tmp="${dir_fq_pe}/tiny_pe_R2.fastq.tmp"

fq_se_gz="${dir_fq_se}/tiny_se.fastq.gz"
fq_r1_gz="${dir_fq_pe}/tiny_pe_R1.fastq.gz"
fq_r2_gz="${dir_fq_pe}/tiny_pe_R2.fastq.gz"


#  Register cleanup of temporary FASTQ intermediates on exit
register_cleanup cleanup_tmp_fastqs

#  Require gzip for deterministic compressed workflow inputs
require_cmd gzip "to generate trim-fastqs fixtures."

#  Create fixture output directories
mkdirs "${dir_fq}" "${dir_fq_se}" "${dir_fq_pe}"

#  Remove stale temporary intermediates
cleanup_tmp_fastqs


#  Write tiny single-end FASTQ provenance and compressed input fixture
cat > "${fq_se_tmp}" << EOM
@tiny_trim_se_read_1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

#  Compress single-end FASTQ fixture with deterministic gzip metadata
gzip_n "${fq_se_tmp}" "${fq_se_gz}"

#  Remove temporary single-end FASTQ intermediate after compression
rm_file "${dir_fix}" "${fq_se_tmp}"


#  Write tiny paired-end FASTQ provenance and compressed input fixtures
cat > "${fq_r1_tmp}" << EOM
@tiny_trim_pe_pair_1/1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

cat > "${fq_r2_tmp}" << EOM
@tiny_trim_pe_pair_1/2
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


succeed "generated trim-fastqs fixtures under ${dir_fix}"
