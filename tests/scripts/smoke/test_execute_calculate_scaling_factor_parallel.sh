#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_calculate_scaling_factor_parallel.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute calculate-scaling-factor GNU Parallel"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

if ! \
    is_parallel_enabled
then
    record_skip \
        "GNU Parallel calculate-scaling-factor check disabled; set" \
        "RUN_PARALLEL=1 to enable"
    finish
    exit $?
fi

#  Define fixture, output, and worker-input paths
scr_exe="${ROOT_REPO}/scripts/execute_calculate_scaling_factor.sh"
cfg_met="${ROOT_REPO}/data/raw/docs/parse_metadata_siqchip.yml"

dir_fix="${ROOT_REPO}/tests/calculate_scaling_factor/fixtures"
dir_bam="${dir_fix}/bam"
dir_bam_pe="${dir_bam}/pe"
dir_met="${dir_fix}/metadata"

tbl_met="${dir_met}/measurements_siqchip.tsv"

tmp="${TEST_DIR_TMP}/execute_calculate_scaling_factor_parallel"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/calculate_scaling_factor"

bam_pe_mip_0="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sc.bam"
bam_pe_mip_1="${dir_bam_pe}/IP_WT_G1_Hho1_6337.sc.bam"
bam_pe_min_0="${dir_bam_pe}/in_WT_G1_Hho1_6336.sc.bam"
bam_pe_min_1="${dir_bam_pe}/in_WT_G1_Hho1_6337.sc.bam"
bam_pe_sip_0="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sp.bam"
bam_pe_sip_1="${dir_bam_pe}/IP_WT_G1_Hho1_6337.sp.bam"
bam_pe_sin_0="${dir_bam_pe}/in_WT_G1_Hho1_6336.sp.bam"
bam_pe_sin_1="${dir_bam_pe}/in_WT_G1_Hho1_6337.sp.bam"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

# shellcheck disable=SC2154
if ! \
    require_env_parallel \
        "${env_nam}" \
        "${dir_log}/execute_calculate_scaling_factor_parallel_env.log"
then
    finish
    exit $?
fi

if ! \
    require_files_nonempty \
        "${scr_exe}" \
        "${cfg_met}" \
        "${tbl_met}" \
        "${bam_pe_mip_0}" \
        "${bam_pe_mip_0}.bai" \
        "${bam_pe_mip_1}" \
        "${bam_pe_mip_1}.bai" \
        "${bam_pe_min_0}" \
        "${bam_pe_min_0}.bai" \
        "${bam_pe_min_1}" \
        "${bam_pe_min_1}.bai" \
        "${bam_pe_sip_0}" \
        "${bam_pe_sip_0}.bai" \
        "${bam_pe_sip_1}" \
        "${bam_pe_sip_1}.bai" \
        "${bam_pe_sin_0}" \
        "${bam_pe_sin_0}.bai" \
        "${bam_pe_sin_1}" \
        "${bam_pe_sin_1}.bai"
then
    finish
    exit $?
fi

hdr_spk=$'^main_ip\tspike_ip\tmain_in\tspike_in\tspike\tcoef\tnum_mp\tnum_sp\tnum_mn\tnum_sn$'
hdr_siq=$'^fil_ip\tfil_in\tsiq\teqn\tmass_ip\tmass_in\tvol_all\tvol_in\tdep_ip\tdep_in\tlen_ip\tlen_in$'

row_spike_0="${bam_pe_mip_0}"$'\t'"${bam_pe_sip_0}"$'\t'"${bam_pe_min_0}"$'\t'"${bam_pe_sin_0}"
row_spike_1="${bam_pe_mip_1}"$'\t'"${bam_pe_sip_1}"$'\t'"${bam_pe_min_1}"$'\t'"${bam_pe_sin_1}"
tail_spike_0=$'2\tchiprx_alpha_ratio\t3\t1\t2\t2'
tail_spike_1=$'0.5\tchiprx_alpha_ratio\t2\t2\t3\t1'

row_siq_0="${bam_pe_mip_0}"$'\t'"${bam_pe_min_0}"
row_siq_1="${bam_pe_mip_1}"$'\t'"${bam_pe_min_1}"
tail_siq_0=$'0.002660098522167487697376\t6nd\t2.7\t72.5\t300\t20\t3\t2\t20\t20'
tail_siq_1=$'0.00440373436674299824356\t6nd\t5\t81.1\t300\t20\t2\t3\t20\t20'


#  GNU Parallel spike execution should write config, worker parts, and table
fil_out="${dir_out}/scaling.parallel.spike.tsv"
nam_job="test_execute_calculate_scaling_factor_parallel_spike"
cfg="${dir_err}/${nam_job}.config_parallel.txt"
log="${dir_log}/execute_parallel_spike.log"

if \
    run_capture \
        "execute calculate-scaling-factor spike GNU Parallel" \
        "${log}" \
        "${TEST_BASH}" "${scr_exe}" \
            --threads 2 \
            --mode spike \
            --csv_mip "${bam_pe_mip_0},${bam_pe_mip_1}" \
            --csv_min "${bam_pe_min_0},${bam_pe_min_1}" \
            --csv_sip "${bam_pe_sip_0},${bam_pe_sip_1}" \
            --csv_sin "${bam_pe_sin_0},${bam_pe_sin_1}" \
            --aln_typ auto \
            --fil_out "${fil_out}" \
            --err_out "${dir_err}" \
            --nam_job "${nam_job}" \
            --max_job 2
then
    record_pass "execute_calculate_scaling_factor.sh spike GNU Parallel exits 0"
else
    record_fail \
        "execute_calculate_scaling_factor.sh spike GNU Parallel failed; see" \
        "$(print_relpath "${log}")"
fi

assert_file_nonempty \
    "${cfg}" \
    "execute scaling-factor spike GNU Parallel config"

assert_pattern_found \
    "${cfg}" \
    "${TEST_BASH} ${ROOT_REPO}/scripts/submit_calculate_scaling_factor.sh" \
    "execute scaling-factor spike config uses Bash-prefixed submit command"

assert_file_nonempty \
    "${fil_out}" \
    "execute scaling-factor spike GNU Parallel final TSV"

assert_file_nonempty \
    "${fil_out}.part.000000" \
    "execute scaling-factor spike GNU Parallel first retained part"

assert_file_nonempty \
    "${fil_out}.part.000001" \
    "execute scaling-factor spike GNU Parallel second retained part"

assert_scaling_factor_header \
    "${fil_out}" \
    "${hdr_spk}" \
    true \
    "execute scaling-factor spike GNU Parallel final TSV"

assert_pattern_found \
    "${fil_out}" \
    "^${row_spike_0}"$'\t'"${tail_spike_0}"'$' \
    "execute scaling-factor spike GNU Parallel final TSV has first row"

assert_pattern_found \
    "${fil_out}" \
    "^${row_spike_1}"$'\t'"${tail_spike_1}"'$' \
    "execute scaling-factor spike GNU Parallel final TSV has second row"


#  GNU Parallel siQ execution should write config, worker parts, and table
fil_out="${dir_out}/scaling.parallel.siq.tsv"
nam_job="test_execute_calculate_scaling_factor_parallel_siq"
cfg="${dir_err}/${nam_job}.config_parallel.txt"
log="${dir_log}/execute_parallel_siq.log"

if \
    run_capture \
        "execute calculate-scaling-factor siQ GNU Parallel" \
        "${log}" \
        "${TEST_BASH}" "${scr_exe}" \
            --threads 2 \
            --mode siq \
            --csv_mip "${bam_pe_mip_0},${bam_pe_mip_1}" \
            --csv_min "${bam_pe_min_0},${bam_pe_min_1}" \
            --aln_typ auto \
            --tbl_met "${tbl_met}" \
            --cfg_met "${cfg_met}" \
            --eqn 6nd \
            --fil_out "${fil_out}" \
            --err_out "${dir_err}" \
            --nam_job "${nam_job}" \
            --max_job 2
then
    record_pass "execute_calculate_scaling_factor.sh siQ GNU Parallel exits 0"
else
    record_fail \
        "execute_calculate_scaling_factor.sh siQ GNU Parallel failed; see" \
        "$(print_relpath "${log}")"
fi

assert_file_nonempty \
    "${cfg}" \
    "execute scaling-factor siQ GNU Parallel config"

assert_pattern_found \
    "${cfg}" \
    "${TEST_BASH} ${ROOT_REPO}/scripts/submit_calculate_scaling_factor.sh" \
    "execute scaling-factor siQ config uses Bash-prefixed submit command"

assert_file_nonempty \
    "${fil_out}" \
    "execute scaling-factor siQ GNU Parallel final TSV"

assert_file_nonempty \
    "${fil_out}.part.000000" \
    "execute scaling-factor siQ GNU Parallel first retained part"

assert_file_nonempty \
    "${fil_out}.part.000001" \
    "execute scaling-factor siQ GNU Parallel second retained part"

assert_scaling_factor_header \
    "${fil_out}" \
    "${hdr_siq}" \
    true \
    "execute scaling-factor siQ GNU Parallel final TSV"

assert_pattern_found \
    "${fil_out}" \
    "^${row_siq_0}"$'\t'"${tail_siq_0}"'$' \
    "execute scaling-factor siQ GNU Parallel final TSV has first row"

assert_pattern_found \
    "${fil_out}" \
    "^${row_siq_1}"$'\t'"${tail_siq_1}"'$' \
    "execute scaling-factor siQ GNU Parallel final TSV has second row"

finish
