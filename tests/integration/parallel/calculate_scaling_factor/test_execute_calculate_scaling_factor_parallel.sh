#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_calculate_scaling_factor_parallel.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in design, development,
# and documentation, with all output reviewed, edited, and approved by the
# author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="execute calculate-scaling-factor GNU Parallel"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


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
scr_exe="${ROOT_REPO}/bin/execute_calculate_scaling_factor.sh"
cfg_met="${ROOT_REPO}/data/raw/docs/parse_metadata_siqchip.yml"

dir_fix="${ROOT_REPO}/tests/fixtures/calculate_scaling_factor"
dir_bam="${dir_fix}/bam"
dir_bam_pe="${dir_bam}/pe"
dir_met="${dir_fix}/metadata"

tbl_met="${dir_met}/measurements_siqchip.tsv"

tmp="${TEST_DIR_TMP}/execute_calculate_scaling_factor_parallel"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/calculate_scaling_factor"

fil_out_spk="${dir_out}/scaling.parallel.spike.tsv"
fil_out_siq="${dir_out}/scaling.parallel.siq.tsv"

prt_spk_0="${fil_out_spk}.part.000000"
prt_spk_1="${fil_out_spk}.part.000001"
prt_siq_0="${fil_out_siq}.part.000000"
prt_siq_1="${fil_out_siq}.part.000001"

nam_job_spk="test_execute_calculate_scaling_factor_parallel_spike"
nam_job_siq="test_execute_calculate_scaling_factor_parallel_siq"

cfg_spk="${dir_err}/${nam_job_spk}.config_parallel.txt"
cfg_siq="${dir_err}/${nam_job_siq}.config_parallel.txt"

log_env_parallel="${dir_log}/execute_calculate_scaling_factor_parallel_env.log"
log_spk="${dir_log}/execute_parallel_spike.log"
log_siq="${dir_log}/execute_parallel_siq.log"

bam_pe_mip_0="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sc.bam"
bam_pe_mip_1="${dir_bam_pe}/IP_WT_G1_Hho1_6337.sc.bam"
bam_pe_min_0="${dir_bam_pe}/in_WT_G1_Hho1_6336.sc.bam"
bam_pe_min_1="${dir_bam_pe}/in_WT_G1_Hho1_6337.sc.bam"
bam_pe_sip_0="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sp.bam"
bam_pe_sip_1="${dir_bam_pe}/IP_WT_G1_Hho1_6337.sp.bam"
bam_pe_sin_0="${dir_bam_pe}/in_WT_G1_Hho1_6336.sp.bam"
bam_pe_sin_1="${dir_bam_pe}/in_WT_G1_Hho1_6337.sp.bam"

bai_pe_mip_0="${bam_pe_mip_0}.bai"
bai_pe_mip_1="${bam_pe_mip_1}.bai"
bai_pe_min_0="${bam_pe_min_0}.bai"
bai_pe_min_1="${bam_pe_min_1}.bai"
bai_pe_sip_0="${bam_pe_sip_0}.bai"
bai_pe_sip_1="${bam_pe_sip_1}.bai"
bai_pe_sin_0="${bam_pe_sin_0}.bai"
bai_pe_sin_1="${bam_pe_sin_1}.bai"

hdr_spk=$'^main_ip\tspike_ip\tmain_in\tspike_in\tspike\tcoef\tnum_mp\tnum_sp\tnum_mn\tnum_sn$'
hdr_siq=$'^fil_ip\tfil_in\tsiq\teqn\tmass_ip\tmass_in\tvol_all\tvol_in\tdep_ip\tdep_in\tlen_ip\tlen_in$'

row_spk_0="${bam_pe_mip_0}"$'\t'"${bam_pe_sip_0}"$'\t'"${bam_pe_min_0}"$'\t'"${bam_pe_sin_0}"
row_spk_1="${bam_pe_mip_1}"$'\t'"${bam_pe_sip_1}"$'\t'"${bam_pe_min_1}"$'\t'"${bam_pe_sin_1}"
tail_spk_0=$'2\tchiprx_alpha_ratio\t3\t1\t2\t2'
tail_spk_1=$'0.5\tchiprx_alpha_ratio\t2\t2\t3\t1'

row_siq_0="${bam_pe_mip_0}"$'\t'"${bam_pe_min_0}"
row_siq_1="${bam_pe_mip_1}"$'\t'"${bam_pe_min_1}"
tail_siq_0=$'0.001912211397724232486706\t6nd\t2.7\t72.5\t300\t20\t3\t2\t626\t450'
tail_siq_1=$'0.002902612244746139314594\t6nd\t5\t81.1\t300\t20\t2\t3\t663\t437'


print_section "${TEST_NAME}"

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
        "${log_env_parallel}"
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
        "${bai_pe_mip_0}" \
        "${bam_pe_mip_1}" \
        "${bai_pe_mip_1}" \
        "${bam_pe_min_0}" \
        "${bai_pe_min_0}" \
        "${bam_pe_min_1}" \
        "${bai_pe_min_1}" \
        "${bam_pe_sip_0}" \
        "${bai_pe_sip_0}" \
        "${bam_pe_sip_1}" \
        "${bai_pe_sip_1}" \
        "${bam_pe_sin_0}" \
        "${bai_pe_sin_0}" \
        "${bam_pe_sin_1}" \
        "${bai_pe_sin_1}"
then
    finish
    exit $?
fi


#  GNU Parallel spike execution should write config, worker parts, and table
if \
    run_capture \
        "execute calculate-scaling-factor spike GNU Parallel" \
        "${log_spk}" \
        "${TEST_BASH}" "${scr_exe}" \
            --env_nam "${env_nam}" \
            --threads 2 \
            --mode spike \
            --aln_typ auto \
            --csv_mip "${bam_pe_mip_0},${bam_pe_mip_1}" \
            --csv_min "${bam_pe_min_0},${bam_pe_min_1}" \
            --csv_sip "${bam_pe_sip_0},${bam_pe_sip_1}" \
            --csv_sin "${bam_pe_sin_0},${bam_pe_sin_1}" \
            --fil_out "${fil_out_spk}" \
            --dir_eo "${dir_err}" \
            --nam_job "${nam_job_spk}" \
            --max_job 2
then
    record_pass "execute_calculate_scaling_factor.sh spike GNU Parallel exits 0"
else
    record_fail \
        "execute_calculate_scaling_factor.sh spike GNU Parallel failed; see" \
        "$(print_relpath "${log_spk}")"
fi

assert_file_nonempty \
    "${cfg_spk}" \
    "execute scaling-factor spike GNU Parallel config"

assert_pattern_found \
    "${cfg_spk}" \
    "${TEST_BASH} ${ROOT_REPO}/bin/submit_calculate_scaling_factor.sh" \
    "execute scaling-factor spike config uses Bash-prefixed submit command"

assert_file_nonempty \
    "${fil_out_spk}" \
    "execute scaling-factor spike GNU Parallel final TSV"

assert_file_nonempty \
    "${prt_spk_0}" \
    "execute scaling-factor spike GNU Parallel first retained part"

assert_file_nonempty \
    "${prt_spk_1}" \
    "execute scaling-factor spike GNU Parallel second retained part"

assert_scaling_factor_header \
    "${fil_out_spk}" \
    "${hdr_spk}" \
    true \
    "execute scaling-factor spike GNU Parallel final TSV"

assert_pattern_found \
    "${fil_out_spk}" \
    "^${row_spk_0}"$'\t'"${tail_spk_0}"'$' \
    "execute scaling-factor spike GNU Parallel final TSV has first row"

assert_pattern_found \
    "${fil_out_spk}" \
    "^${row_spk_1}"$'\t'"${tail_spk_1}"'$' \
    "execute scaling-factor spike GNU Parallel final TSV has second row"


#  GNU Parallel siQ execution should write config, worker parts, and table
if \
    run_capture \
        "execute calculate-scaling-factor siQ GNU Parallel" \
        "${log_siq}" \
        "${TEST_BASH}" "${scr_exe}" \
            --env_nam "${env_nam}" \
            --threads 2 \
            --mode siq \
            --aln_typ auto \
            --csv_mip "${bam_pe_mip_0},${bam_pe_mip_1}" \
            --csv_min "${bam_pe_min_0},${bam_pe_min_1}" \
            --tbl_met "${tbl_met}" \
            --cfg_met "${cfg_met}" \
            --eqn 6nd \
            --fil_out "${fil_out_siq}" \
            --dir_eo "${dir_err}" \
            --nam_job "${nam_job_siq}" \
            --max_job 2
then
    record_pass "execute_calculate_scaling_factor.sh siQ GNU Parallel exits 0"
else
    record_fail \
        "execute_calculate_scaling_factor.sh siQ GNU Parallel failed; see" \
        "$(print_relpath "${log_siq}")"
fi

assert_file_nonempty \
    "${cfg_siq}" \
    "execute scaling-factor siQ GNU Parallel config"

assert_pattern_found \
    "${cfg_siq}" \
    "${TEST_BASH} ${ROOT_REPO}/bin/submit_calculate_scaling_factor.sh" \
    "execute scaling-factor siQ config uses Bash-prefixed submit command"

assert_file_nonempty \
    "${fil_out_siq}" \
    "execute scaling-factor siQ GNU Parallel final TSV"

assert_file_nonempty \
    "${prt_siq_0}" \
    "execute scaling-factor siQ GNU Parallel first retained part"

assert_file_nonempty \
    "${prt_siq_1}" \
    "execute scaling-factor siQ GNU Parallel second retained part"

assert_scaling_factor_header \
    "${fil_out_siq}" \
    "${hdr_siq}" \
    true \
    "execute scaling-factor siQ GNU Parallel final TSV"

assert_pattern_found \
    "${fil_out_siq}" \
    "^${row_siq_0}"$'\t'"${tail_siq_0}"'$' \
    "execute scaling-factor siQ GNU Parallel final TSV has first row"

assert_pattern_found \
    "${fil_out_siq}" \
    "^${row_siq_1}"$'\t'"${tail_siq_1}"'$' \
    "execute scaling-factor siQ GNU Parallel final TSV has second row"

finish
