#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_calculate_scaling_factor_opt_applicability.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="calculate-scaling-factor option applicability"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. Per 'HELP.PARAMETER.APPLICABILITY', each
# mode refuses the other mode's options, and '--ref_fa' without CRAM input is
# ignored with a warning and not checked, so a missing path is not refused.
dir_fix="${ROOT_REPO}/tests/fixtures/calculate_scaling_factor"
dir_pe="${dir_fix}/bam/pe"
dir_cram="${dir_fix}/cram/pe"
ref_fa="${dir_fix}/reference/tiny.fa"
tbl_met="${dir_fix}/metadata/measurements_siqchip.tsv"
cfg_met="${ROOT_REPO}/data/raw/docs/parse_metadata_siqchip.yml"

smp_0="WT_G1_Hho1_6336"
smp_1="WT_G1_Hho1_6337"

tmp="${TEST_DIR_TMP}/calculate_scaling_factor_opt_applicability"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${dir_pe}/IP_${smp_0}.sc.bam" \
    "${dir_cram}/IP_${smp_0}.sc.cram" \
    "${ref_fa}" \
    "${tbl_met}" \
    "${cfg_met}" || {
    finish
    exit $?
}


# Print the alignment files of samples 0 and 1 for one input type ('bam' or
# 'mix' with sample 0 as CRAM), file prefix, and genome, comma-separated.
function csv_aln() {
    local typ="${1}" pfx="${2}" gen="${3}"
    local dir_0="${dir_pe}" ext_0="bam"

    if [[ "${typ}" == "mix" ]]; then dir_0="${dir_cram}"; ext_0="cram"; fi
    printf '%s,%s' \
        "${dir_0}/${pfx}_${smp_0}.${gen}.${ext_0}" \
        "${dir_pe}/${pfx}_${smp_1}.${gen}.bam"
}


# Print the input-list arguments for one mode and input type ('bam' or 'mix',
# where sample 0 is CRAM).
function args_mode() {
    local mode="${1}" typ="${2}"

    printf '%s\n' \
        --mode "${mode}" \
        --csv_mip "$(csv_aln "${typ}" IP sc)" \
        --csv_min "$(csv_aln "${typ}" in sc)"

    if [[ "${mode}" == "spike" ]]; then
        printf '%s\n' \
            --csv_sip "$(csv_aln "${typ}" IP sp)" \
            --csv_sin "$(csv_aln "${typ}" in sp)"
    else
        printf '%s\n' --tbl_met "${tbl_met}" --cfg_met "${cfg_met}"
    fi
}


# Both per-sample fragment-length lists, which leave '--len_def' no file.
len_both="--csv_len_mip 200 --csv_len_min 180"

# Every count each mode reads from a file, which leaves no file to read: the
# layout, the threads (in a serial local run), and the reference then warn.
dep_spk="--csv_dep_mip 10 --csv_dep_min 20 --csv_dep_sip 2 --csv_dep_sin 4"
dep_siq="--csv_dep_mip 10 --csv_dep_min 20 ${len_both}"
dep_3="--csv_dep_mip 10 --csv_dep_min 20 --csv_dep_sip 2"

# 'execute' rows are dry runs: label|mode|input|arguments|act option. Refused
# rows name the option; warned rows warn and drop it; control rows do neither
# and pass on the option they name.
rows_exe=(
    "--tbl_met in spike|spike|bam|--tbl_met ${tbl_met}|refuse --tbl_met"
    "--cfg_met in spike|spike|bam|--cfg_met ${cfg_met}|refuse --cfg_met"
    "default --eqn in spike|spike|bam|--eqn 6nd|refuse --eqn"
    "invalid --eqn in spike|spike|bam|--eqn 7|refuse --eqn"
    "--len_def in spike|spike|bam|--len_def 200|refuse --len_def"
    "--csv_len_mip in spike|spike|bam|--csv_len_mip 200|refuse --csv_len_mip"
    "--csv_len_min in spike|spike|bam|--csv_len_min 200|refuse --csv_len_min"
    "--method in siq|siq|bam|--method fractional|refuse --method"
    "--csv_dep_sip in siq|siq|bam|--csv_dep_sip 100|refuse --csv_dep_sip"
    "--csv_dep_sin in siq|siq|bam|--csv_dep_sin 100|refuse --csv_dep_sin"
    "--ref_fa with BAM only|spike|bam|--ref_fa ${ref_fa}|warn --ref_fa"
    "--time without --slurm|spike|bam|--time 1:00:00|warn --time"
    "missing --ref_fa, BAM only|spike|bam|--ref_fa ${tmp}/no.fa|warn --ref_fa"
    "--ref_fa with a CRAM entry|spike|mix|--ref_fa ${ref_fa}|pass --ref_fa"
    "--eqn in siq|siq|bam|--eqn 5|pass --eqn"
    "--csv_len_mip in siq|siq|bam|--csv_len_mip 200|pass --csv_len_mip"
    "--len_def with PE|siq|bam|--aln_typ pe --len_def 200|warn --len_def"
    "--len_def with lengths|siq|bam|${len_both} --len_def 200|warn --len_def"
    "--len_def in siq|siq|bam|--len_def 200|pass --len_def"
    "--len_def, one list|siq|bam|--csv_len_mip 200 --len_def 9|pass --len_def"
    "--csv_dep_sip in spike|spike|bam|--csv_dep_sip 100|pass --csv_dep_sip"
    "--aln_typ, every count|spike|bam|${dep_spk} --aln_typ pe|warn --aln_typ"
    "bad --aln_typ, all counts|spike|bam|${dep_spk} --aln_typ x|warn --aln_typ"
    "--aln_typ, every siq count|siq|bam|${dep_siq} --aln_typ pe|warn --aln_typ"
    "--threads, every count|spike|bam|${dep_spk}|warn --threads"
    "--ref_fa, all|spike|mix|${dep_spk} --ref_fa ${ref_fa}|warn --ref_fa"
    "--aln_typ, three counts|spike|bam|${dep_3} --aln_typ pe|pass --aln_typ"
)

for idx in "${!rows_exe[@]}"; do
    IFS='|' read -r lbl mode typ arg act <<< "${rows_exe[idx]}"
    mapfile -t args_in < <(args_mode "${mode}" "${typ}")
    read -r -a args <<< "${arg}"
    read -r -a acts <<< "${act}"
    dir_out="${tmp}/exe_${idx}"

    mkdir -p "${dir_out}/logs"

    rc=0
    out="$(
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_calculate_scaling_factor.sh" \
            --dry_run \
            --env_nam "${env_nam}" \
            --threads 1 \
            "${args_in[@]}" \
            --fil_out "${dir_out}/sf.tsv" \
            --max_job 1 \
            "${args[@]}" 2>&1
    )" || rc=$?

    cmds="$(
        grep 'submit_calculate_scaling_factor.sh --env' <<< "${out}" || true
    )"
    lbl="execute ${lbl}"
    opt="${acts[1]}"

    case "${acts[0]}" in
        refuse)
            if [[
                "${rc}" -ne 0
                && "${out}" == *"error("*"'${opt}' is for '--mode "*
            ]]; then
                record_pass "${lbl} is refused by name"
            else
                record_fail "${lbl} was not refused by name (exit ${rc})"
            fi
            ;;

        warn)
            if [[
                "${rc}" -eq 0
                && "${out}" == *"warning("*"'${opt}' has no effect"*
                && -n "${cmds}"
                && "${cmds}" != *" ${opt} "*
            ]]; then
                record_pass "${lbl} warns and is not passed on"
            else
                record_fail \
                    "${lbl} did not warn, or was passed on (exit ${rc})"
            fi
            ;;

        pass)
            # With a mixed list, only the CRAM sample gets the reference.
            n_exp=2
            if [[ "${typ}" == "mix" ]]; then n_exp=1; fi
            n_opt="$(grep -c -- " ${opt} " <<< "${cmds}" || true)"

            if [[
                "${rc}" -eq 0
                && "${n_opt}" -eq "${n_exp}"
                && "${out}" != *"has no effect"*
            ]]; then
                record_pass "${lbl} is accepted with no message"
            else
                record_fail \
                    "${lbl} was dropped or flagged (exit ${rc}; passed" \
                    "${n_opt} of ${n_exp} times)"
            fi
            ;;
    esac
done


# With every count supplied, a CRAM run needs no reference.
dir_out="${tmp}/exe_cram_no_ref"
mkdir -p "${dir_out}/logs"
mapfile -t args_in < <(args_mode spike mix)
read -r -a args <<< "${dep_spk}"

rc=0
out="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_calculate_scaling_factor.sh" \
        --dry_run \
        --env_nam "${env_nam}" \
        "${args_in[@]}" \
        --fil_out "${dir_out}/sf.tsv" \
        --max_job 1 \
        "${args[@]}" 2>&1
)" || rc=$?

if [[ "${rc}" -eq 0 && "${out}" != *"'--ref_fa' is required"* ]]; then
    record_pass "execute takes CRAM with every count and no '--ref_fa'"
else
    record_fail "execute required '--ref_fa' with every count (exit ${rc})"
fi


# Run 'submit' directly on sample 0, writing to its own output directory.
function run_submit() {
    local dir_out="${1}" mode="${2}"
    shift 2

    mkdir -p "${dir_out}/logs"
    "${TEST_BASH}" "${ROOT_REPO}/bin/submit_calculate_scaling_factor.sh" \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --threads 1 \
        --mode "${mode}" \
        --csv_mip "${dir_pe}/IP_${smp_0}.sc.bam" \
        --csv_min "${dir_pe}/in_${smp_0}.sc.bam" \
        --fil_out "${dir_out}/sf.tsv" \
        --idx_out 0 \
        --dir_eo "${dir_out}/logs" \
        "$@" 2>&1
}

spk=(
    --csv_sip "${dir_pe}/IP_${smp_0}.sp.bam"
    --csv_sin "${dir_pe}/in_${smp_0}.sp.bam"
)
siq=( --tbl_met "${tbl_met}" --cfg_met "${cfg_met}" )

# A refusal in 'submit' stops before any work: label|mode|arguments|option.
rows_sub=(
    "default --eqn in spike|spike|--eqn 6nd|--eqn"
    "--len_def in spike|spike|--len_def 200|--len_def"
    "--csv_len_min in spike|spike|--csv_len_min 200|--csv_len_min"
    "--tbl_met in spike|spike|--tbl_met ${tbl_met}|--tbl_met"
    "--method in siq|siq|--method fractional|--method"
    "--csv_dep_sip in siq|siq|--csv_dep_sip 100|--csv_dep_sip"
    "--csv_sin in siq|siq|--csv_sin x.bam|--csv_sin"
)

for idx in "${!rows_sub[@]}"; do
    IFS='|' read -r lbl mode arg opt <<< "${rows_sub[idx]}"
    read -r -a args <<< "${arg}"
    dir_out="${tmp}/sub_${idx}"
    lbl="submit ${lbl}"

    if [[ "${mode}" == "spike" ]]; then
        args=( "${spk[@]}" "${args[@]}" )
    else
        args=( "${siq[@]}" "${args[@]}" )
    fi

    rc=0
    out="$(run_submit "${dir_out}" "${mode}" "${args[@]}")" || rc=$?

    if [[
        "${rc}" -ne 0
        && "${out}" == *"'${opt}' is for '--mode "*
        && ! -e "${dir_out}/sf.tsv.part.000000"
    ]]; then
        record_pass "${lbl} is refused by name before any work"
    else
        record_fail "${lbl} was not refused (exit ${rc})"
    fi
done


# A reference with BAM input alone warns and is not checked, so a missing path
# is not refused, and the row is still written: label|reference.
rows_ref=(
    "--ref_fa|${ref_fa}"
    "missing --ref_fa|${tmp}/no.fa"
)

for idx in "${!rows_ref[@]}"; do
    IFS='|' read -r lbl fil_ref <<< "${rows_ref[idx]}"
    dir_out="${tmp}/sub_ref_fa_${idx}"
    lbl="submit ${lbl} with BAM only"

    rc=0
    out="$(
        run_submit "${dir_out}" spike "${spk[@]}" --ref_fa "${fil_ref}"
    )" || rc=$?

    # The helper passes the calculator only the counts the coefficient uses,
    # so no other note appears.
    if [[
        "${rc}" -eq 0
        && "${out}" == *"'--ref_fa' has no effect without CRAM input"*
        && "${out}" != *"has no effect with '--coef"*
        && -s "${dir_out}/sf.tsv.part.000000"
    ]]; then
        record_pass "${lbl} warns and still runs"
    else
        record_fail "${lbl} did not warn, warned twice, or failed (exit ${rc})"
    fi
done


# With every count supplied, 'submit' reads no file: CRAM needs no reference,
# and the layout and threads warn without being checked. Past its four paths,
# the row matches a BAM run of the same counts.
read -r -a args <<< "${dep_spk}"
fil_cram=(
    --csv_mip "${dir_cram}/IP_${smp_0}.sc.cram"
    --csv_min "${dir_cram}/in_${smp_0}.sc.cram"
    --csv_sip "${dir_cram}/IP_${smp_0}.sp.cram"
    --csv_sin "${dir_cram}/in_${smp_0}.sp.cram"
)
dir_bam="${tmp}/sub_cnt_bam"
dir_cnt="${tmp}/sub_cnt_cram"
mkdir -p "${dir_bam}/logs" "${dir_cnt}/logs"

rc_bam=0
"${TEST_BASH}" "${ROOT_REPO}/bin/submit_calculate_scaling_factor.sh" \
    --env_nam "${env_nam}" \
    --dir_scr "${ROOT_REPO}/bin" \
    --mode spike \
    --csv_mip "${dir_pe}/IP_${smp_0}.sc.bam" \
    --csv_min "${dir_pe}/in_${smp_0}.sc.bam" \
    "${spk[@]}" \
    --fil_out "${dir_bam}/sf.tsv" \
    --idx_out 0 \
    --dir_eo "${dir_bam}/logs" \
    "${args[@]}" > /dev/null 2>&1 || rc_bam=$?

rc=0
out="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/submit_calculate_scaling_factor.sh" \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --threads x \
        --mode spike \
        --aln_typ x \
        "${fil_cram[@]}" \
        --fil_out "${dir_cnt}/sf.tsv" \
        --idx_out 0 \
        --dir_eo "${dir_cnt}/logs" \
        "${args[@]}" 2>&1
)" || rc=$?

if [[
    "${rc_bam}" -eq 0
    && "${rc}" -eq 0
    && "${out}" == *"'--aln_typ' has no effect with all four"*
    && "${out}" == *"'--threads' has no effect with all four"*
    && -s "${dir_cnt}/sf.tsv.part.000000"
]] && [[
    "$(cut -f 5- "${dir_bam}/sf.tsv.part.000000")" \
        == "$(cut -f 5- "${dir_cnt}/sf.tsv.part.000000")"
]]; then
    record_pass "submit with every count reads no CRAM and warns"
else
    record_fail "submit with every count failed or did not warn (exit ${rc})"
fi

rc=0
out="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/submit_calculate_scaling_factor.sh" \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --mode spike \
        "${fil_cram[@]}" \
        --ref_fa "${tmp}/no.fa" \
        --fil_out "${dir_cnt}/sf_ref.tsv" \
        --idx_out 0 \
        --dir_eo "${dir_cnt}/logs" \
        "${args[@]}" 2>&1
)" || rc=$?

if [[
    "${rc}" -eq 0
    && "${out}" == *"'--ref_fa' has no effect with all four"*
    && -s "${dir_cnt}/sf_ref.tsv.part.000000"
]]; then
    record_pass "submit '--ref_fa' with every count warns and is not checked"
else
    record_fail \
        "submit '--ref_fa' with every count was not ignored (exit ${rc})"
fi


# siQ mode reads no CRAM either once both depths and both lengths are given.
read -r -a args <<< "${dep_siq}"
dir_siq="${tmp}/sub_cnt_siq_cram"
mkdir -p "${dir_siq}/logs"

rc=0
out="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/submit_calculate_scaling_factor.sh" \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --mode siq \
        --csv_mip "${dir_cram}/IP_${smp_0}.sc.cram" \
        --csv_min "${dir_cram}/in_${smp_0}.sc.cram" \
        "${siq[@]}" \
        --fil_out "${dir_siq}/sf.tsv" \
        --idx_out 0 \
        --dir_eo "${dir_siq}/logs" \
        "${args[@]}" 2>&1
)" || rc=$?

if [[ "${rc}" -eq 0 && -s "${dir_siq}/sf.tsv.part.000000" ]]; then
    record_pass "submit siq with every count reads no CRAM"
else
    record_fail "submit siq with every count failed on CRAM (exit ${rc})"
fi


# A default fragment length that no file can use warns and is not checked, so
# the run still writes its row: label|arguments.
rows_len=(
    "--len_def with PE|--aln_typ pe --len_def 0"
    "--len_def with paired|--aln_typ paired --len_def 0"
    "--len_def with lengths|${len_both} --len_def 0"
)

for idx in "${!rows_len[@]}"; do
    IFS='|' read -r lbl arg <<< "${rows_len[idx]}"
    read -r -a args <<< "${arg}"
    dir_out="${tmp}/sub_len_${idx}"
    lbl="submit ${lbl}"

    rc=0
    out="$(run_submit "${dir_out}" siq "${siq[@]}" "${args[@]}")" || rc=$?

    if [[
        "${rc}" -eq 0
        && "${out}" == *"'--len_def' has no effect"*
        && -s "${dir_out}/sf.tsv.part.000000"
    ]]; then
        record_pass "${lbl} warns and still runs"
    else
        record_fail "${lbl} did not warn or failed (exit ${rc})"
    fi
done


# Each method gets only the counts it uses, so the calculator never notes an
# unused one.
for method in \
    fractional \
    main_per_spike \
    chiprx_alpha_ip \
    chiprx_alpha_in \
    chiprx_alpha_ratio \
    rxinput_alpha
do
    dir_out="${tmp}/sub_method_${method}"
    lbl="submit --method ${method}"

    rc=0
    out="$(
        run_submit "${dir_out}" spike "${spk[@]}" --method "${method}"
    )" || rc=$?

    if [[
        "${rc}" -eq 0
        && "${out}" != *"has no effect"*
        && -s "${dir_out}/sf.tsv.part.000000"
    ]]; then
        record_pass "${lbl} passes only the counts it uses"
    else
        record_fail "${lbl} passed an unused count or failed (exit ${rc})"
    fi
done

finish
