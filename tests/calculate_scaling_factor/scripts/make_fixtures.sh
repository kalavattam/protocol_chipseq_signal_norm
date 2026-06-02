#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: tests/calculate_scaling_factor/scripts/make_fixtures.sh
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


#  Resolve paths relative to 'tests/calculate_scaling_factor/scripts'
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_scl="$(cd "${dir_scr}/.." > /dev/null 2>&1 && pwd)"
dir_fix="${dir_scl}/fixtures"

#  Source shared fixture-generation helpers
# shellcheck disable=SC1091
source "${dir_scr}/../../scripts/lib/fixture_helpers.sh"


#  Generate one sorted SE BAM fixture with a known alignment count
function write_bam_fixture_se() {
    local sam="${1:-}"
    local bam="${2:-}"
    local pfx="${3:-read}"
    local num="${4:-}"
    local idx pos

    {
        write_sam_line '@HD' 'VN:1.6' 'SO:coordinate'
        write_sam_line '@SQ' 'SN:I' 'LN:100'

        for (( idx = 1; idx <= num; idx++ )); do
            pos=$(( 1 + ( idx - 1 ) * 10 ))

            write_sam_line \
                "${pfx}_${idx}" '0' 'I' "${pos}" '60' '10M' \
                '*' '0' '0' 'ACGTACGTAA' 'FFFFFFFFFF'
        done
    } > "${sam}"

    samtools view -bS "${sam}" | samtools sort -o "${bam}"
}


#  Generate one sorted PE BAM fixture with a known fragment count
function write_bam_fixture_pe() {
    local sam="${1:-}"
    local bam="${2:-}"
    local pfx="${3:-read}"
    local num="${4:-}"
    local idx pos mate_pos tlen

    {
        write_sam_line '@HD' 'VN:1.6' 'SO:coordinate'
        write_sam_line '@SQ' 'SN:I' 'LN:100'

        for (( idx = 1; idx <= num; idx++ )); do
            pos=$(( 1 + ( idx - 1 ) * 20 ))
            mate_pos=$(( pos + 10 ))
            tlen=20

            write_sam_line \
                "${pfx}_${idx}" '99' 'I' "${pos}" '60' '10M' \
                '=' "${mate_pos}" "${tlen}" 'ACGTACGTAA' 'FFFFFFFFFF'

            write_sam_line \
                "${pfx}_${idx}" '147' 'I' "${mate_pos}" '60' '10M' \
                '=' "${pos}" "-${tlen}" 'ACGTACGTAA' 'FFFFFFFFFF'
        done
    } > "${sam}"

    samtools view -bS "${sam}" | samtools sort -o "${bam}"
}


#  Generate one reference-backed CRAM fixture from a BAM fixture
function write_cram_fixture() {
    local bam="${1:-}"
    local cram="${2:-}"
    local ref_fa="${3:-}"

    samtools view -C -T "${ref_fa}" -o "${cram}" "${bam}"
    samtools index "${cram}"
}


#  Define fixture directories and scaling-factor part-file paths
dir_ref="${dir_fix}/reference"
dir_sam="${dir_fix}/sam"
dir_bam="${dir_fix}/bam"
dir_sam_pe="${dir_sam}/pe"
dir_bam_pe="${dir_bam}/pe"
dir_cram="${dir_fix}/cram"
dir_cram_se="${dir_cram}/se"
dir_cram_pe="${dir_cram}/pe"
dir_prt="${dir_fix}/parts"

ref_fa="${dir_ref}/tiny.fa"

sam_mip_0="${dir_sam}/IP_A.sc.sam"
sam_mip_1="${dir_sam}/IP_B.sc.sam"
sam_min_0="${dir_sam}/in_A.sc.sam"
sam_min_1="${dir_sam}/in_B.sc.sam"
sam_sip_0="${dir_sam}/IP_A.sp.sam"
sam_sip_1="${dir_sam}/IP_B.sp.sam"
sam_sin_0="${dir_sam}/in_A.sp.sam"
sam_sin_1="${dir_sam}/in_B.sp.sam"

bam_mip_0="${dir_bam}/IP_A.sc.bam"
bam_mip_1="${dir_bam}/IP_B.sc.bam"
bam_min_0="${dir_bam}/in_A.sc.bam"
bam_min_1="${dir_bam}/in_B.sc.bam"
bam_sip_0="${dir_bam}/IP_A.sp.bam"
bam_sip_1="${dir_bam}/IP_B.sp.bam"
bam_sin_0="${dir_bam}/in_A.sp.bam"
bam_sin_1="${dir_bam}/in_B.sp.bam"

sam_pe_mip_0="${dir_sam_pe}/IP_A.sc.sam"
sam_pe_mip_1="${dir_sam_pe}/IP_B.sc.sam"
sam_pe_min_0="${dir_sam_pe}/in_A.sc.sam"
sam_pe_min_1="${dir_sam_pe}/in_B.sc.sam"
sam_pe_sip_0="${dir_sam_pe}/IP_A.sp.sam"
sam_pe_sip_1="${dir_sam_pe}/IP_B.sp.sam"
sam_pe_sin_0="${dir_sam_pe}/in_A.sp.sam"
sam_pe_sin_1="${dir_sam_pe}/in_B.sp.sam"

bam_pe_mip_0="${dir_bam_pe}/IP_A.sc.bam"
bam_pe_mip_1="${dir_bam_pe}/IP_B.sc.bam"
bam_pe_min_0="${dir_bam_pe}/in_A.sc.bam"
bam_pe_min_1="${dir_bam_pe}/in_B.sc.bam"
bam_pe_sip_0="${dir_bam_pe}/IP_A.sp.bam"
bam_pe_sip_1="${dir_bam_pe}/IP_B.sp.bam"
bam_pe_sin_0="${dir_bam_pe}/in_A.sp.bam"
bam_pe_sin_1="${dir_bam_pe}/in_B.sp.bam"

cram_se_mip_0="${dir_cram_se}/IP_A.sc.cram"
cram_se_mip_1="${dir_cram_se}/IP_B.sc.cram"
cram_se_min_0="${dir_cram_se}/in_A.sc.cram"
cram_se_min_1="${dir_cram_se}/in_B.sc.cram"
cram_se_sip_0="${dir_cram_se}/IP_A.sp.cram"
cram_se_sip_1="${dir_cram_se}/IP_B.sp.cram"
cram_se_sin_0="${dir_cram_se}/in_A.sp.cram"
cram_se_sin_1="${dir_cram_se}/in_B.sp.cram"

cram_pe_mip_0="${dir_cram_pe}/IP_A.sc.cram"
cram_pe_mip_1="${dir_cram_pe}/IP_B.sc.cram"
cram_pe_min_0="${dir_cram_pe}/in_A.sc.cram"
cram_pe_min_1="${dir_cram_pe}/in_B.sc.cram"
cram_pe_sip_0="${dir_cram_pe}/IP_A.sp.cram"
cram_pe_sip_1="${dir_cram_pe}/IP_B.sp.cram"
cram_pe_sin_0="${dir_cram_pe}/in_A.sp.cram"
cram_pe_sin_1="${dir_cram_pe}/in_B.sp.cram"

spk_0="${dir_prt}/example_scaling_factors.spike.tsv.part.000000"
spk_2="${dir_prt}/example_scaling_factors.spike.tsv.part.000002"
siq_0="${dir_prt}/example_scaling_factors.siq.tsv.part.000000"
siq_2="${dir_prt}/example_scaling_factors.siq.tsv.part.000002"
bad_fld="${dir_prt}/malformed_scaling_factors.spike.tsv.part.000003"

env_req="env_protocol"


#  Require the project environment and Samtools for alignment fixtures
require_env "${env_req}" "for calculate-scaling-factor fixtures."
require_cmd samtools "in '${env_req}' to generate BAM/CRAM fixtures."

#  Create fixture output directories
mkdirs \
    "${dir_ref}" \
    "${dir_sam}" \
    "${dir_bam}" \
    "${dir_sam_pe}" \
    "${dir_bam_pe}" \
    "${dir_cram_se}" \
    "${dir_cram_pe}" \
    "${dir_prt}"

#  Write and index the tiny reference used for CRAM fixtures
cat > "${ref_fa}" << 'EOM'
>I
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOM

samtools faidx "${ref_fa}"

#  Generate role-specific SE BAM fixtures with known alignment counts
write_bam_fixture_se "${sam_mip_0}" "${bam_mip_0}" mip_A 3
write_bam_fixture_se "${sam_mip_1}" "${bam_mip_1}" mip_B 2
write_bam_fixture_se "${sam_min_0}" "${bam_min_0}" min_A 2
write_bam_fixture_se "${sam_min_1}" "${bam_min_1}" min_B 3
write_bam_fixture_se "${sam_sip_0}" "${bam_sip_0}" sip_A 1
write_bam_fixture_se "${sam_sip_1}" "${bam_sip_1}" sip_B 2
write_bam_fixture_se "${sam_sin_0}" "${bam_sin_0}" sin_A 2
write_bam_fixture_se "${sam_sin_1}" "${bam_sin_1}" sin_B 1

#  Generate role-specific PE BAM fixtures with matching fragment counts
write_bam_fixture_pe "${sam_pe_mip_0}" "${bam_pe_mip_0}" mip_A 3
write_bam_fixture_pe "${sam_pe_mip_1}" "${bam_pe_mip_1}" mip_B 2
write_bam_fixture_pe "${sam_pe_min_0}" "${bam_pe_min_0}" min_A 2
write_bam_fixture_pe "${sam_pe_min_1}" "${bam_pe_min_1}" min_B 3
write_bam_fixture_pe "${sam_pe_sip_0}" "${bam_pe_sip_0}" sip_A 1
write_bam_fixture_pe "${sam_pe_sip_1}" "${bam_pe_sip_1}" sip_B 2
write_bam_fixture_pe "${sam_pe_sin_0}" "${bam_pe_sin_0}" sin_A 2
write_bam_fixture_pe "${sam_pe_sin_1}" "${bam_pe_sin_1}" sin_B 1

#  Generate SE and PE CRAM fixtures from the role-specific BAM fixtures
write_cram_fixture "${bam_mip_0}" "${cram_se_mip_0}" "${ref_fa}"
write_cram_fixture "${bam_mip_1}" "${cram_se_mip_1}" "${ref_fa}"
write_cram_fixture "${bam_min_0}" "${cram_se_min_0}" "${ref_fa}"
write_cram_fixture "${bam_min_1}" "${cram_se_min_1}" "${ref_fa}"
write_cram_fixture "${bam_sip_0}" "${cram_se_sip_0}" "${ref_fa}"
write_cram_fixture "${bam_sip_1}" "${cram_se_sip_1}" "${ref_fa}"
write_cram_fixture "${bam_sin_0}" "${cram_se_sin_0}" "${ref_fa}"
write_cram_fixture "${bam_sin_1}" "${cram_se_sin_1}" "${ref_fa}"

write_cram_fixture "${bam_pe_mip_0}" "${cram_pe_mip_0}" "${ref_fa}"
write_cram_fixture "${bam_pe_mip_1}" "${cram_pe_mip_1}" "${ref_fa}"
write_cram_fixture "${bam_pe_min_0}" "${cram_pe_min_0}" "${ref_fa}"
write_cram_fixture "${bam_pe_min_1}" "${cram_pe_min_1}" "${ref_fa}"
write_cram_fixture "${bam_pe_sip_0}" "${cram_pe_sip_0}" "${ref_fa}"
write_cram_fixture "${bam_pe_sip_1}" "${cram_pe_sip_1}" "${ref_fa}"
write_cram_fixture "${bam_pe_sin_0}" "${cram_pe_sin_0}" "${ref_fa}"
write_cram_fixture "${bam_pe_sin_1}" "${cram_pe_sin_1}" "${ref_fa}"

samtools quickcheck \
    "${bam_mip_0}" \
    "${bam_mip_1}" \
    "${bam_min_0}" \
    "${bam_min_1}" \
    "${bam_sip_0}" \
    "${bam_sip_1}" \
    "${bam_sin_0}" \
    "${bam_sin_1}" \
    "${bam_pe_mip_0}" \
    "${bam_pe_mip_1}" \
    "${bam_pe_min_0}" \
    "${bam_pe_min_1}" \
    "${bam_pe_sip_0}" \
    "${bam_pe_sip_1}" \
    "${bam_pe_sin_0}" \
    "${bam_pe_sin_1}" \
    "${cram_se_mip_0}" \
    "${cram_se_mip_1}" \
    "${cram_se_min_0}" \
    "${cram_se_min_1}" \
    "${cram_se_sip_0}" \
    "${cram_se_sip_1}" \
    "${cram_se_sin_0}" \
    "${cram_se_sin_1}" \
    "${cram_pe_mip_0}" \
    "${cram_pe_mip_1}" \
    "${cram_pe_min_0}" \
    "${cram_pe_min_1}" \
    "${cram_pe_sip_0}" \
    "${cram_pe_sip_1}" \
    "${cram_pe_sin_0}" \
    "${cram_pe_sin_1}"

#  Remove stale generated fixture outputs
rm_files \
    "${dir_fix}" \
    "${spk_0}" \
    "${spk_2}" \
    "${siq_0}" \
    "${siq_2}" \
    "${bad_fld}"

#  Write realistic spike-in scaling-factor part files
write_tsv_row \
    '/path/to/IP_WT_G1_Hho1_6336.sc.bam' \
    '/path/to/IP_WT_G1_Hho1_6336.sp.bam' \
    '/path/to/in_WT_G1_Hho1_6336.sc.bam' \
    '/path/to/in_WT_G1_Hho1_6336.sp.bam' \
    '2.145085661466039' \
    'fractional' \
    '13492920' '217340' '12851824' '452406' \
    > "${spk_0}"

write_tsv_row \
    '/path/to/IP_WT_G1_Hho1_6337.sc.bam' \
    '/path/to/IP_WT_G1_Hho1_6337.sp.bam' \
    '/path/to/in_WT_G1_Hho1_6337.sc.bam' \
    '/path/to/in_WT_G1_Hho1_6337.sp.bam' \
    '3.228015767161274' \
    'fractional' \
    '13655994' '116947' '12030091' '339029' \
    > "${spk_2}"

#  Write realistic siQ-ChIP scaling-factor part files
write_tsv_row \
    '/path/to/IP_WT_G1_Hho1_6336.sc.bam' \
    '/path/to/in_WT_G1_Hho1_6336.sc.bam' \
    '0.0023106316806475436' \
    '6nd' \
    '2.7' '72.5' '300' '20' \
    '13492920' '12851824' '227.009' '197.186' \
    > "${siq_0}"

write_tsv_row \
    '/path/to/IP_WT_G1_Hho1_6337.sc.bam' \
    '/path/to/in_WT_G1_Hho1_6337.sc.bam' \
    '0.003816492180174805' \
    '6nd' \
    '5' '81.1' '300' '20' \
    '13655994' '12030091' '230.767' '199.994' \
    > "${siq_2}"

#  Write one malformed part file for negative validation
write_tsv_row \
    '/path/to/IP_WT_G1_bad.sc.bam' \
    '/path/to/IP_WT_G1_bad.sp.bam' \
    > "${bad_fld}"


succeed "generated calculate-scaling-factor fixtures under ${dir_fix}"
