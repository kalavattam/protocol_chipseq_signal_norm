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


#  Generate one sorted SE BAM fixture with a known alignment count
function write_fixture_bam_se() {
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
    samtools index "${bam}"
}


#  Generate one sorted PE BAM fixture with a known fragment count
function write_fixture_bam_pe() {
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
    samtools index "${bam}"
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
dir_sam_se="${dir_sam}/se"
dir_bam="${dir_fix}/bam"
dir_bam_se="${dir_bam}/se"
dir_sam_pe="${dir_sam}/pe"
dir_bam_pe="${dir_bam}/pe"
dir_cram="${dir_fix}/cram"
dir_cram_se="${dir_cram}/se"
dir_cram_pe="${dir_cram}/pe"
dir_prt="${dir_fix}/parts"
dir_met="${dir_fix}/metadata"
dir_cfg="${dir_fix}/config"

ref_fa="${dir_ref}/tiny.fa"

sam_se_mip_0="${dir_sam_se}/IP_WT_log_Brn1_rep1.sc.sam"
sam_se_mip_1="${dir_sam_se}/IP_WT_log_Brn1_rep2.sc.sam"
sam_se_min_0="${dir_sam_se}/in_WT_log_Brn1_rep1.sc.sam"
sam_se_min_1="${dir_sam_se}/in_WT_log_Brn1_rep2.sc.sam"
sam_se_sip_0="${dir_sam_se}/IP_WT_log_Brn1_rep1.sp.sam"
sam_se_sip_1="${dir_sam_se}/IP_WT_log_Brn1_rep2.sp.sam"
sam_se_sin_0="${dir_sam_se}/in_WT_log_Brn1_rep1.sp.sam"
sam_se_sin_1="${dir_sam_se}/in_WT_log_Brn1_rep2.sp.sam"

bam_se_mip_0="${dir_bam_se}/IP_WT_log_Brn1_rep1.sc.bam"
bam_se_mip_1="${dir_bam_se}/IP_WT_log_Brn1_rep2.sc.bam"
bam_se_min_0="${dir_bam_se}/in_WT_log_Brn1_rep1.sc.bam"
bam_se_min_1="${dir_bam_se}/in_WT_log_Brn1_rep2.sc.bam"
bam_se_sip_0="${dir_bam_se}/IP_WT_log_Brn1_rep1.sp.bam"
bam_se_sip_1="${dir_bam_se}/IP_WT_log_Brn1_rep2.sp.bam"
bam_se_sin_0="${dir_bam_se}/in_WT_log_Brn1_rep1.sp.bam"
bam_se_sin_1="${dir_bam_se}/in_WT_log_Brn1_rep2.sp.bam"

sam_pe_mip_0="${dir_sam_pe}/IP_WT_G1_Hho1_6336.sc.sam"
sam_pe_mip_1="${dir_sam_pe}/IP_WT_G1_Hho1_6337.sc.sam"
sam_pe_min_0="${dir_sam_pe}/in_WT_G1_Hho1_6336.sc.sam"
sam_pe_min_1="${dir_sam_pe}/in_WT_G1_Hho1_6337.sc.sam"
sam_pe_sip_0="${dir_sam_pe}/IP_WT_G1_Hho1_6336.sp.sam"
sam_pe_sip_1="${dir_sam_pe}/IP_WT_G1_Hho1_6337.sp.sam"
sam_pe_sin_0="${dir_sam_pe}/in_WT_G1_Hho1_6336.sp.sam"
sam_pe_sin_1="${dir_sam_pe}/in_WT_G1_Hho1_6337.sp.sam"

bam_pe_mip_0="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sc.bam"
bam_pe_mip_1="${dir_bam_pe}/IP_WT_G1_Hho1_6337.sc.bam"
bam_pe_min_0="${dir_bam_pe}/in_WT_G1_Hho1_6336.sc.bam"
bam_pe_min_1="${dir_bam_pe}/in_WT_G1_Hho1_6337.sc.bam"
bam_pe_mip_hu_0="${dir_bam_pe}/IP_WT_G1_HU_Hho1_6336.sc.bam"
bam_pe_mip_hu_1="${dir_bam_pe}/IP_WT_G1_HU_Hho1_6337.sc.bam"
bam_pe_min_hu_0="${dir_bam_pe}/in_WT_G1_HU_Hho1_6336.sc.bam"
bam_pe_min_hu_1="${dir_bam_pe}/in_WT_G1_HU_Hho1_6337.sc.bam"
bam_pe_sip_0="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sp.bam"
bam_pe_sip_1="${dir_bam_pe}/IP_WT_G1_Hho1_6337.sp.bam"
bam_pe_sin_0="${dir_bam_pe}/in_WT_G1_Hho1_6336.sp.bam"
bam_pe_sin_1="${dir_bam_pe}/in_WT_G1_Hho1_6337.sp.bam"

cram_se_mip_0="${dir_cram_se}/IP_WT_log_Brn1_rep1.sc.cram"
cram_se_mip_1="${dir_cram_se}/IP_WT_log_Brn1_rep2.sc.cram"
cram_se_min_0="${dir_cram_se}/in_WT_log_Brn1_rep1.sc.cram"
cram_se_min_1="${dir_cram_se}/in_WT_log_Brn1_rep2.sc.cram"
cram_se_sip_0="${dir_cram_se}/IP_WT_log_Brn1_rep1.sp.cram"
cram_se_sip_1="${dir_cram_se}/IP_WT_log_Brn1_rep2.sp.cram"
cram_se_sin_0="${dir_cram_se}/in_WT_log_Brn1_rep1.sp.cram"
cram_se_sin_1="${dir_cram_se}/in_WT_log_Brn1_rep2.sp.cram"

cram_pe_mip_0="${dir_cram_pe}/IP_WT_G1_Hho1_6336.sc.cram"
cram_pe_mip_1="${dir_cram_pe}/IP_WT_G1_Hho1_6337.sc.cram"
cram_pe_min_0="${dir_cram_pe}/in_WT_G1_Hho1_6336.sc.cram"
cram_pe_min_1="${dir_cram_pe}/in_WT_G1_Hho1_6337.sc.cram"
cram_pe_mip_hu_0="${dir_cram_pe}/IP_WT_G1_HU_Hho1_6336.sc.cram"
cram_pe_mip_hu_1="${dir_cram_pe}/IP_WT_G1_HU_Hho1_6337.sc.cram"
cram_pe_min_hu_0="${dir_cram_pe}/in_WT_G1_HU_Hho1_6336.sc.cram"
cram_pe_min_hu_1="${dir_cram_pe}/in_WT_G1_HU_Hho1_6337.sc.cram"
cram_pe_sip_0="${dir_cram_pe}/IP_WT_G1_Hho1_6336.sp.cram"
cram_pe_sip_1="${dir_cram_pe}/IP_WT_G1_Hho1_6337.sp.cram"
cram_pe_sin_0="${dir_cram_pe}/in_WT_G1_Hho1_6336.sp.cram"
cram_pe_sin_1="${dir_cram_pe}/in_WT_G1_Hho1_6337.sp.cram"

spk_0="${dir_prt}/example_scaling_factors.spike.tsv.part.000000"
spk_2="${dir_prt}/example_scaling_factors.spike.tsv.part.000002"
siq_0="${dir_prt}/example_scaling_factors.siq.tsv.part.000000"
siq_2="${dir_prt}/example_scaling_factors.siq.tsv.part.000002"
bad_fld="${dir_prt}/malformed_scaling_factors.spike.tsv.part.000003"
bad_hdr="${dir_prt}/header_scaling_factors.spike.tsv.part.000004"
dup_idx_a="${dir_prt}/duplicate_index_A.spike.tsv.part.000005"
dup_idx_b="${dir_prt}/duplicate_index_B.spike.tsv.part.000005"
met_siq="${dir_met}/measurements_siqchip.tsv"
met_siq_gz="${dir_met}/measurements_siqchip.tsv.gz"
met_siq_lib="${dir_met}/measurements_siqchip_lib_volume.tsv"
met_siq_lib_one="${dir_met}/measurements_siqchip_lib_volume_one_sided.tsv"
met_siq_lib_zero="${dir_met}/measurements_siqchip_lib_volume_zero.tsv"
met_siq_mis="${dir_met}/measurements_siqchip_missing_required.tsv"
met_siq_dup="${dir_met}/measurements_siqchip_duplicate_match.tsv"
met_siq_pre="${dir_met}/measurements_siqchip_precomputed.tsv"
cfg_map="${dir_cfg}/parse_metadata_siqchip_field_to_column.yml"

env_req="env_protocol"


#  Require the project environment and Samtools for alignment fixtures
require_env "${env_req}" "for calculate-scaling-factor fixtures."
require_cmd samtools "in '${env_req}' to generate BAM/CRAM fixtures."
require_cmd gzip "in '${env_req}' to generate compressed metadata fixtures."

#  Create fixture output directories
mkdirs \
    "${dir_ref}" \
    "${dir_sam}" \
    "${dir_sam_se}" \
    "${dir_bam}" \
    "${dir_bam_se}" \
    "${dir_sam_pe}" \
    "${dir_bam_pe}" \
    "${dir_cram_se}" \
    "${dir_cram_pe}" \
    "${dir_prt}" \
    "${dir_met}" \
    "${dir_cfg}"

#  Write and index the tiny reference used for CRAM fixtures
cat > "${ref_fa}" << EOM
>I
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOM

samtools faidx "${ref_fa}"

#  Generate role-specific SE BAM fixtures with known alignment counts
write_fixture_bam_se "${sam_se_mip_0}" "${bam_se_mip_0}" mip_A 3
write_fixture_bam_se "${sam_se_mip_1}" "${bam_se_mip_1}" mip_B 2
write_fixture_bam_se "${sam_se_min_0}" "${bam_se_min_0}" min_A 2
write_fixture_bam_se "${sam_se_min_1}" "${bam_se_min_1}" min_B 3
write_fixture_bam_se "${sam_se_sip_0}" "${bam_se_sip_0}" sip_A 1
write_fixture_bam_se "${sam_se_sip_1}" "${bam_se_sip_1}" sip_B 2
write_fixture_bam_se "${sam_se_sin_0}" "${bam_se_sin_0}" sin_A 2
write_fixture_bam_se "${sam_se_sin_1}" "${bam_se_sin_1}" sin_B 1

#  Generate role-specific PE BAM fixtures with matching fragment counts
write_fixture_bam_pe "${sam_pe_mip_0}" "${bam_pe_mip_0}" mip_A 3
write_fixture_bam_pe "${sam_pe_mip_1}" "${bam_pe_mip_1}" mip_B 2
write_fixture_bam_pe "${sam_pe_min_0}" "${bam_pe_min_0}" min_A 2
write_fixture_bam_pe "${sam_pe_min_1}" "${bam_pe_min_1}" min_B 3
write_fixture_bam_pe "${sam_pe_sip_0}" "${bam_pe_sip_0}" sip_A 1
write_fixture_bam_pe "${sam_pe_sip_1}" "${bam_pe_sip_1}" sip_B 2
write_fixture_bam_pe "${sam_pe_sin_0}" "${bam_pe_sin_0}" sin_A 2
write_fixture_bam_pe "${sam_pe_sin_1}" "${bam_pe_sin_1}" sin_B 1

#  Duplicate PE main-alignment fixtures under treatment-aware siQ basenames
cp "${bam_pe_mip_0}"     "${bam_pe_mip_hu_0}"
cp "${bam_pe_mip_0}.bai" "${bam_pe_mip_hu_0}.bai"
cp "${bam_pe_mip_1}"     "${bam_pe_mip_hu_1}"
cp "${bam_pe_mip_1}.bai" "${bam_pe_mip_hu_1}.bai"
cp "${bam_pe_min_0}"     "${bam_pe_min_hu_0}"
cp "${bam_pe_min_0}.bai" "${bam_pe_min_hu_0}.bai"
cp "${bam_pe_min_1}"     "${bam_pe_min_hu_1}"
cp "${bam_pe_min_1}.bai" "${bam_pe_min_hu_1}.bai"

#  Generate SE and PE CRAM fixtures from the role-specific BAM fixtures
write_cram_fixture "${bam_se_mip_0}" "${cram_se_mip_0}" "${ref_fa}"
write_cram_fixture "${bam_se_mip_1}" "${cram_se_mip_1}" "${ref_fa}"
write_cram_fixture "${bam_se_min_0}" "${cram_se_min_0}" "${ref_fa}"
write_cram_fixture "${bam_se_min_1}" "${cram_se_min_1}" "${ref_fa}"
write_cram_fixture "${bam_se_sip_0}" "${cram_se_sip_0}" "${ref_fa}"
write_cram_fixture "${bam_se_sip_1}" "${cram_se_sip_1}" "${ref_fa}"
write_cram_fixture "${bam_se_sin_0}" "${cram_se_sin_0}" "${ref_fa}"
write_cram_fixture "${bam_se_sin_1}" "${cram_se_sin_1}" "${ref_fa}"

write_cram_fixture "${bam_pe_mip_0}" "${cram_pe_mip_0}" "${ref_fa}"
write_cram_fixture "${bam_pe_mip_1}" "${cram_pe_mip_1}" "${ref_fa}"
write_cram_fixture "${bam_pe_min_0}" "${cram_pe_min_0}" "${ref_fa}"
write_cram_fixture "${bam_pe_min_1}" "${cram_pe_min_1}" "${ref_fa}"
write_cram_fixture "${bam_pe_mip_hu_0}" "${cram_pe_mip_hu_0}" "${ref_fa}"
write_cram_fixture "${bam_pe_mip_hu_1}" "${cram_pe_mip_hu_1}" "${ref_fa}"
write_cram_fixture "${bam_pe_min_hu_0}" "${cram_pe_min_hu_0}" "${ref_fa}"
write_cram_fixture "${bam_pe_min_hu_1}" "${cram_pe_min_hu_1}" "${ref_fa}"
write_cram_fixture "${bam_pe_sip_0}" "${cram_pe_sip_0}" "${ref_fa}"
write_cram_fixture "${bam_pe_sip_1}" "${cram_pe_sip_1}" "${ref_fa}"
write_cram_fixture "${bam_pe_sin_0}" "${cram_pe_sin_0}" "${ref_fa}"
write_cram_fixture "${bam_pe_sin_1}" "${cram_pe_sin_1}" "${ref_fa}"

samtools quickcheck \
    "${bam_se_mip_0}" \
    "${bam_se_mip_1}" \
    "${bam_se_min_0}" \
    "${bam_se_min_1}" \
    "${bam_se_sip_0}" \
    "${bam_se_sip_1}" \
    "${bam_se_sin_0}" \
    "${bam_se_sin_1}" \
    "${bam_pe_mip_0}" \
    "${bam_pe_mip_1}" \
    "${bam_pe_min_0}" \
    "${bam_pe_min_1}" \
    "${bam_pe_mip_hu_0}" \
    "${bam_pe_mip_hu_1}" \
    "${bam_pe_min_hu_0}" \
    "${bam_pe_min_hu_1}" \
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
    "${cram_pe_mip_hu_0}" \
    "${cram_pe_mip_hu_1}" \
    "${cram_pe_min_hu_0}" \
    "${cram_pe_min_hu_1}" \
    "${cram_pe_sip_0}" \
    "${cram_pe_sip_1}" \
    "${cram_pe_sin_0}" \
    "${cram_pe_sin_1}"

for idx in \
    "${bam_se_mip_0}.bai" \
    "${bam_se_mip_1}.bai" \
    "${bam_se_min_0}.bai" \
    "${bam_se_min_1}.bai" \
    "${bam_se_sip_0}.bai" \
    "${bam_se_sip_1}.bai" \
    "${bam_se_sin_0}.bai" \
    "${bam_se_sin_1}.bai" \
    "${bam_pe_mip_0}.bai" \
    "${bam_pe_mip_1}.bai" \
    "${bam_pe_min_0}.bai" \
    "${bam_pe_min_1}.bai" \
    "${bam_pe_mip_hu_0}.bai" \
    "${bam_pe_mip_hu_1}.bai" \
    "${bam_pe_min_hu_0}.bai" \
    "${bam_pe_min_hu_1}.bai" \
    "${bam_pe_sip_0}.bai" \
    "${bam_pe_sip_1}.bai" \
    "${bam_pe_sin_0}.bai" \
    "${bam_pe_sin_1}.bai"
do
    [[ -s "${idx}" ]] || die "BAM index missing or empty: '${idx}'."
done

for idx in \
    "${cram_se_mip_0}.crai" \
    "${cram_se_mip_1}.crai" \
    "${cram_se_min_0}.crai" \
    "${cram_se_min_1}.crai" \
    "${cram_se_sip_0}.crai" \
    "${cram_se_sip_1}.crai" \
    "${cram_se_sin_0}.crai" \
    "${cram_se_sin_1}.crai" \
    "${cram_pe_mip_0}.crai" \
    "${cram_pe_mip_1}.crai" \
    "${cram_pe_min_0}.crai" \
    "${cram_pe_min_1}.crai" \
    "${cram_pe_mip_hu_0}.crai" \
    "${cram_pe_mip_hu_1}.crai" \
    "${cram_pe_min_hu_0}.crai" \
    "${cram_pe_min_hu_1}.crai" \
    "${cram_pe_sip_0}.crai" \
    "${cram_pe_sip_1}.crai" \
    "${cram_pe_sin_0}.crai" \
    "${cram_pe_sin_1}.crai"
do
    [[ -s "${idx}" ]] || die "CRAM index missing or empty: '${idx}'."
done

#  Remove stale generated fixture outputs. Broad metadata/config cleanup clears
#+ ignored legacy parser fixtures from earlier generator versions.
rm_files \
    "${dir_fix}" \
    "${spk_0}" \
    "${spk_2}" \
    "${siq_0}" \
    "${siq_2}" \
    "${bad_fld}" \
    "${bad_hdr}" \
    "${dup_idx_a}" \
    "${dup_idx_b}"

rm -f -- \
    "${dir_met}"/measurements_siqchip*.tsv \
    "${dir_met}"/measurements_siqchip.tsv.gz \
    "${dir_cfg}"/parse_metadata_siqchip*.yml

#  Write realistic spike-in scaling-factor part files
write_tsv_row \
    '/path/to/IP_WT_G1_Hho1_6336.sc.bam' \
    '/path/to/IP_WT_G1_Hho1_6336.sp.bam' \
    '/path/to/in_WT_G1_Hho1_6336.sc.bam' \
    '/path/to/in_WT_G1_Hho1_6336.sp.bam' \
    '2.081558847888101748679901' \
    'chiprx_alpha_ratio' \
    '13492920' '217340' '12851824' '452406' \
    > "${spk_0}"

write_tsv_row \
    '/path/to/IP_WT_G1_Hho1_6337.sc.bam' \
    '/path/to/IP_WT_G1_Hho1_6337.sp.bam' \
    '/path/to/in_WT_G1_Hho1_6337.sc.bam' \
    '/path/to/in_WT_G1_Hho1_6337.sp.bam' \
    '2.898996981538645822951139' \
    'chiprx_alpha_ratio' \
    '13655994' '116947' '12030091' '339029' \
    > "${spk_2}"

#  Write realistic siQ-ChIP scaling-factor part files
write_tsv_row \
    '/path/to/IP_WT_G1_Hho1_6336.sc.bam' \
    '/path/to/in_WT_G1_Hho1_6336.sc.bam' \
    '0.0023106316806475436' \
    '6nd' \
    '2.7' '72.5' '300' '20' \
    '13492920' '12851824' '227.009' '197.186' 'NA' 'NA' \
    > "${siq_0}"

write_tsv_row \
    '/path/to/IP_WT_G1_Hho1_6337.sc.bam' \
    '/path/to/in_WT_G1_Hho1_6337.sc.bam' \
    '0.003816492180174805' \
    '6nd' \
    '5' '81.1' '300' '20' \
    '13655994' '12030091' '230.767' '199.994' 'NA' 'NA' \
    > "${siq_2}"

#  Write one malformed part file for negative validation
write_tsv_row \
    '/path/to/IP_WT_G1_bad.sc.bam' \
    '/path/to/IP_WT_G1_bad.sp.bam' \
    > "${bad_fld}"

#  Write one header-looking part file for negative validation
write_tsv_row \
    'main_ip' \
    'spike_ip' \
    'main_in' \
    'spike_in' \
    'spike' \
    'coef' \
    'num_mp' \
    'num_sp' \
    'num_mn' \
    'num_sn' \
    > "${bad_hdr}"

#  Write duplicate-index part files for negative validation
write_tsv_row \
    '/path/to/IP_WT_G1_dupA.sc.bam' \
    '/path/to/IP_WT_G1_dupA.sp.bam' \
    '/path/to/in_WT_G1_dupA.sc.bam' \
    '/path/to/in_WT_G1_dupA.sp.bam' \
    '1.25' \
    'chiprx_alpha_ratio' \
    '100' '80' '200' '64' \
    > "${dup_idx_a}"

write_tsv_row \
    '/path/to/IP_WT_G1_dupB.sc.bam' \
    '/path/to/IP_WT_G1_dupB.sp.bam' \
    '/path/to/in_WT_G1_dupB.sc.bam' \
    '/path/to/in_WT_G1_dupB.sp.bam' \
    '0.75' \
    'chiprx_alpha_ratio' \
    '100' '60' '200' '80' \
    > "${dup_idx_b}"

#  Write minimal siQ-ChIP metadata fixtures for production-YAML parsing tests
{
    write_tsv_row \
        'genotype_full' 'genotype' 'state' 'factor' 'strain_full' 'strain' \
        'volume_in' 'volume_all' 'mass_in' 'mass_ip' \
        'conc_in' 'conc_ip' 'length_in' 'length_ip'

    write_tsv_row \
        'HHO1-2L-3FLAG' 'WT' 'G1' 'Hho1' 'yTT6336' '6336' \
        '20' '300' '72.5' '2.7' \
        '6.04' '0.224' '450' '626'

    write_tsv_row \
        'HHO1-2L-3FLAG' 'WT' 'G1' 'Hho1' 'yTT6337' '6337' \
        '20' '300' '81.1' '5' \
        '6.76' '0.42' '437' '663'

    write_tsv_row \
        'HHO1-2L-3FLAG' 'WT' 'G2M' 'Hho1' 'yTT6336' '6336' \
        '20' '300' '104.9' '6.6' \
        '8.74' '0.55' '450' '667'

    write_tsv_row \
        'HMO1-2L-3FLAG' 'WT' 'G1' 'Hmo1' 'yTT7750' '7750' \
        '20' '300' '79.9' '8.4' \
        '6.66' '0.704' '450' '499'

    write_tsv_row \
        'BRN1-2L-3FLAG' 'WT' 'log' 'Brn1' 'rep1' 'rep1' \
        '20' '300' '60' '4' \
        '5' '0.333' '150' '150'

    write_tsv_row \
        'BRN1-2L-3FLAG' 'WT' 'log' 'Brn1' 'rep2' 'rep2' \
        '20' '300' '75' '5' \
        '6.25' '0.417' '150' '150'
} > "${met_siq}"

gzip_n "${met_siq}" "${met_siq_gz}"

{
    write_tsv_row \
        'factor' 'genotype' 'state' 'strain' \
        'volume_in' 'volume_all' 'mass_in' 'mass_ip' \
        'length_in' 'length_ip' 'lib_vol_in' 'lib_vol_ip'

    write_tsv_row \
        'Hho1' 'WT' 'G1' '6336' \
        '20' '300' '72.5' '2.7' \
        '450' '626' '4' '2'

    write_tsv_row \
        'Hho1' 'WT' 'G1' '6337' \
        '20' '300' '81.1' '5' \
        '437' '663' '4' '2'
} > "${met_siq_lib}"

{
    write_tsv_row \
        'factor' 'genotype' 'state' 'strain' \
        'volume_in' 'volume_all' 'mass_in' 'mass_ip' \
        'length_in' 'length_ip' 'lib_vol_in'

    write_tsv_row \
        'Hho1' 'WT' 'G1' '6336' \
        '20' '300' '72.5' '2.7' \
        '450' '626' '4'
} > "${met_siq_lib_one}"

{
    write_tsv_row \
        'factor' 'genotype' 'state' 'strain' \
        'volume_in' 'volume_all' 'mass_in' 'mass_ip' \
        'length_in' 'length_ip' 'lib_vol_in' 'lib_vol_ip'

    write_tsv_row \
        'Hho1' 'WT' 'G1' '6336' \
        '20' '300' '72.5' '2.7' \
        '450' '626' '0' '2'
} > "${met_siq_lib_zero}"

{
    write_tsv_row \
        'genotype' 'state' 'factor' 'strain' \
        'volume_in' 'volume_all' 'mass_in' 'length_in' 'length_ip'

    write_tsv_row \
        'WT' 'G1' 'Hho1' '6336' \
        '20' '300' '72.5' '450' '626'
} > "${met_siq_mis}"

{
    write_tsv_row \
        'genotype' 'state' 'factor' 'strain' \
        'volume_in' 'volume_all' 'mass_in' 'mass_ip' \
        'length_in' 'length_ip'

    write_tsv_row \
        'WT' 'G1' 'Hho1' '6336' \
        '20' '300' '72.5' '2.7' \
        '450' '626'

    write_tsv_row \
        'WT' 'G1' 'Hho1' '6336' \
        '20' '300' '73.5' '3.7' \
        '450' '626'
} > "${met_siq_dup}"

{
    write_tsv_row \
        'genotype' 'state' 'factor' 'strain' \
        'volume_in' 'volume_all' 'mass_in' 'mass_ip' \
        'length_in' 'length_ip' 'dep_in' 'dep_ip'

    write_tsv_row \
        'WT' 'G1' 'Hho1' '6336' \
        '20' '300' '72.5' '2.7' \
        '123' '456' '2222' '3333'
} > "${met_siq_pre}"

cat > "${cfg_map}" << EOM
version: 2
profile: "fixture-field-to-column"

filename:
    delimiter: "_"
    strip_extensions:
        - ".bam"
        - ".cram"
        - ".sam"
    strip_suffixes:
        - ".sc"
        - ".sp"
    fields:
        - assay
        - genotype
        - state
        - factor
        - id

matching:
    fields:
        - genotype
        - state
        - factor
        - id

field_to_column:
    id: strain

table:
    skip_prefixes:
        - "#"
        - "//"

calculator_inputs:
    siqchip:
        required:
            vol_in: volume_in
            vol_all: volume_all
            mass_in: mass_in
            mass_ip: mass_ip
        optional:
            len_in: length_in
            len_ip: length_ip
            dep_in: dep_in
            dep_ip: dep_ip
            lib_vol_in: lib_vol_in
            lib_vol_ip: lib_vol_ip
EOM

succeed "generated calculate-scaling-factor fixtures under ${dir_fix}"
