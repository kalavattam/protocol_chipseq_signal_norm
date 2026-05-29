#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: make_fixtures.sh
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


#  Resolve paths relative to 'tests/filter_alignments/scripts'
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_flt="$(cd "${dir_scr}/.." > /dev/null 2>&1 && pwd)"
dir_fix="${dir_flt}/fixtures"

#  Source shared fixture-generation helpers
# shellcheck disable=SC1091
source "${dir_scr}/../../scripts/lib/fixture_helpers.sh"


#  Define fixture paths, contig names, and the shared reference sequence
dir_sam="${dir_fix}/sam"
dir_ref="${dir_fix}/reference"

sam="${dir_sam}/filter_sc_sp.sam"
ref="${dir_ref}/filter_sc_sp.fa"
fai="${ref}.fai"

contigs=(
    I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI Mito
    SP_I SP_II SP_III SP_II_TG SP_MTR SP_Mito
    chrUn
)

seq_100="ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"


#  Create fixture output directories
mkdirs "${dir_sam}" "${dir_ref}"

#  Write the tiny reference FASTA and matching FASTA index
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

#  Write SAM headers and one representative read per retained/dropped class
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


succeed "generated filter-alignments fixtures under ${dir_fix}"
