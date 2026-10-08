#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_align_fastqs_dry_run.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="align-fastqs dry run"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. A dry run validates, then prints each step's
# commands from the arrays a real run executes, and writes nothing; 'execute'
# forwards it through 'submit' to 'align_fastqs'.
dir_fx="${ROOT_REPO}/tests/fixtures/align_fastqs"
in_se="${dir_fx}/fastq/se/tiny_se.atria.fastq.gz"
in_pe_1="${dir_fx}/fastq/pe/tiny_pe_R1.atria.fastq.gz"
in_pe_2="${dir_fx}/fastq/pe/tiny_pe_R2.atria.fastq.gz"
ref_fa="${dir_fx}/reference/tiny.fa"
idx_bt2="${dir_fx}/bowtie2/tiny"
idx_bwa="${dir_fx}/bwa/tiny.fa"
idx_bm2="${dir_fx}/bwa-mem2/tiny.fa"

tmp="${TEST_DIR_TMP}/align_fastqs_dry_run"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${in_se}" \
    "${in_pe_1}" \
    "${in_pe_2}" \
    "${ref_fa}" \
    "${idx_bt2}.1.bt2" \
    "${idx_bwa}.bwt" \
    "${idx_bm2}.0123" || {
    finish
    exit $?
}

# The expected lines below spell paths as written, so the test paths must not
# need shell escaping themselves.
for pth in "${dir_fx}" "${tmp}"; do
    if [[ "$(printf '%q' "${pth}")" != "${pth}" ]]; then
        record_fail "path needs shell escaping, so rows cannot match: ${pth}"
        finish
        exit $?
    fi
done


# Run the helper.
function run_helper() {
    # shellcheck disable=SC2016  # Expand in the child shell, not this one.
    "${TEST_BASH}" -c '
        # shellcheck source=lib/bash/core/source_helpers.sh
        source "${1}/lib/bash/core/source_helpers.sh"

        source_helpers "${1}/lib/bash" \
            core/check_args \
            core/check_inputs \
            core/format_outputs \
            workflows/align_fastqs
        shift 1
        align_fastqs "$@"
    ' _ "${ROOT_REPO}" "$@"
}


# Print the lines of a dry run from one step header up to the next.
function print_step() {
    awk -v s="# Step ${1}:" '
        index($0, "# Step ") == 1 { on = (index($0, s) == 1) }
        on
    ' <<< "${2}"
}


# Compare one block of dry-run output with the lines expected for it.
function check_block() {
    local lbl="${1}"
    local exp="${2}"
    local got="${3}"

    if [[ "${got}" == "${exp}" ]]; then
        record_pass "${lbl}"
    else
        record_fail "${lbl}"
        diff <(printf '%s\n' "${exp}") <(printf '%s\n' "${got}") >&2 || true
    fi
}


# Bowtie 2, paired-end, into a directory and file name with spaces: every step
# prints in order, and nothing is written.
dir_sp="${tmp}/out dir"
esc="${tmp}/out\\ dir/tiny\\ pe"
msg_nil="# Nothing to do: the BAM work file is the final output."
mkdir -p "${dir_sp}"

rc=0
out="$(
    run_helper --dry_run --threads 2 --aligner bowtie2 --index "${idx_bt2}" \
        --fq_1 "${in_pe_1}" --fq_2 "${in_pe_2}" --req_flg --qname \
        --fil_out "${dir_sp}/tiny pe.bam" 2> /dev/null
)" || rc=$?

exp="
# Step 1: align the reads into the BAM work file.
bowtie2 \\
    -p 2 \\
    -x ${idx_bt2} \\
    --very-sensitive \\
    --no-unal \\
    --phred33 \\
    --no-mixed \\
    --no-discordant \\
    --no-overlap \\
    --no-dovetail \\
    -1 ${in_pe_1} \\
    -2 ${in_pe_2} \\
    | samtools view \\
        -@ 2 \\
        -F 12 \\
        -f 2 \\
        -q 1 \\
        -o ${esc}.bam


# Step 2: sort the BAM work file by queryname.
samtools sort \\
    -@ 2 \\
    -n \\
    -o ${esc}.qnam.bam \\
    ${esc}.bam

samtools fixmate \\
    -@ 2 \\
    -c \\
    -m \\
    ${esc}.qnam.bam \\
    ${esc}.mate.bam

mv -f \\
    ${esc}.mate.bam \\
    ${esc}.qnam.bam


# Step 3: sort the BAM work file by coordinates and index it.
samtools sort \\
    -@ 2 \\
    -o ${esc}.coor.bam \\
    ${esc}.qnam.bam

mv -f \\
    ${esc}.coor.bam \\
    ${esc}.bam

samtools index \\
    -@ 2 \\
    ${esc}.bam


# Step 4: mark duplicates and index the BAM work file again.
samtools markdup \\
    -@ 2 \\
    -t \\
    ${esc}.bam \\
    ${esc}.mark.bam

mv -f \\
    ${esc}.mark.bam \\
    ${esc}.bam

samtools index \\
    -@ 2 \\
    ${esc}.bam


# Step 5: write the final output and the retained files.
${msg_nil}"

if [[ "${rc}" -eq 0 ]]; then
    record_pass "helper dry run returns 0"
else
    record_fail "helper dry run returned ${rc}"
fi

check_block \
    "Bowtie 2 PE dry run prints every step, spaces escaped" "${exp}" "${out}"

if [[ -z "$(find "${dir_sp}" -mindepth 1 -print -quit)" ]]; then
    record_pass "helper dry run writes no file"
else
    record_fail "helper dry run wrote files in '${dir_sp}'"
fi


# BWA-backtrack, paired-end, '--mapq 20': the '.sai' redirects, their removal,
# and the template filter's pipeline.
t="${tmp}/aln"
out="$(
    run_helper --dry_run --threads 1 --aligner bwa --bwa_alg aln \
        --mapq 20 --index "${idx_bwa}" --fq_1 "${in_pe_1}" \
        --fq_2 "${in_pe_2}" --fil_out "${t}.bam" 2> /dev/null
)" || true

exp="# Step 1: align the reads into the BAM work file.
bwa aln \\
    -t 1 \\
    ${idx_bwa} \\
    ${in_pe_1} \\
    > ${t}.R1.sai

bwa aln \\
    -t 1 \\
    ${idx_bwa} \\
    ${in_pe_2} \\
    > ${t}.R2.sai

bwa sampe \\
    ${idx_bwa} \\
    ${t}.R1.sai \\
    ${t}.R2.sai \\
    ${in_pe_1} \\
    ${in_pe_2} \\
    | samtools view \\
        -@ 1 \\
        -F 12 \\
        -o ${t}.bam

rm -f \\
    ${t}.R1.sai \\
    ${t}.R2.sai"
check_block "BWA-backtrack PE dry run prints Step 1" \
    "${exp}" "$(print_step 1 "${out}")"

exp="# Step 2: sort the BAM work file by queryname.
samtools sort \\
    -@ 1 \\
    -n \\
    -o ${t}.qnam.bam \\
    ${t}.bam

samtools fixmate \\
    -@ 1 \\
    -c \\
    -m \\
    ${t}.qnam.bam \\
    ${t}.mate.bam

mv -f \\
    ${t}.mate.bam \\
    ${t}.qnam.bam

samtools view \\
    -@ 1 \\
    -F 0x900 \\
    -q 20 \\
    ${t}.qnam.bam \\
    | cut -f 1 \\
    | uniq -c \\
    | awk -v n=2 \\
        \\\$1\\ ==\\ n\\ \\{\\ print\\ \\\$2\\ \\} \\
        > ${t}.mapq.names.txt

samtools view \\
    -@ 1 \\
    -b \\
    -N ${t}.mapq.names.txt \\
    -o ${t}.mapq.bam \\
    ${t}.qnam.bam

rm -f \\
    ${t}.mapq.names.txt

mv -f \\
    ${t}.mapq.bam \\
    ${t}.qnam.bam"
check_block "BWA-backtrack PE dry run prints the template filter" \
    "${exp}" "$(print_step 2 "${out}")"

# The printed commands are the ones a real run executes: run as a script, they
# write the same files, holding the same records.
dir_par="${tmp}/parity"
mkdir -p "${dir_par}/printed" "${dir_par}/real"
out="$(
    run_helper --dry_run --threads 1 --aligner bwa --bwa_alg aln \
        --mapq 20 --index "${idx_bwa}" --fq_1 "${in_pe_1}" \
        --fq_2 "${in_pe_2}" --fil_out "${dir_par}/printed/aln.bam" \
        2> /dev/null
)" || true

rc=0
"${TEST_BASH}" -e -c "${out}" > /dev/null 2>&1 || rc=$?
run_helper --threads 1 --aligner bwa --bwa_alg aln --mapq 20 \
    --index "${idx_bwa}" --fq_1 "${in_pe_1}" --fq_2 "${in_pe_2}" \
    --fil_out "${dir_par}/real/aln.bam" > /dev/null 2>&1 || rc=$?

if [[
    "${rc}" -eq 0
    && "$(cd "${dir_par}/printed" && find . -type f | LC_ALL=C sort)" \
        == "$(cd "${dir_par}/real" && find . -type f | LC_ALL=C sort)"
]] && cmp -s \
    <(run_samtools view "${dir_par}/printed/aln.bam" | cut -f 1-11) \
    <(run_samtools view "${dir_par}/real/aln.bam" | cut -f 1-11)
then
    record_pass "printed commands, run as a script, match a real run"
else
    record_fail "printed commands, run as a script, differ from a real run"
fi


# BWA-MEM2, single-end, CRAM with '--qname': Step 5 converts both files,
# indexes the output, and removes the BAM work files.
t="${tmp}/cram"
out="$(
    run_helper --dry_run --threads 1 --aligner bwa-mem2 --index "${idx_bm2}" \
        --ref_fa "${ref_fa}" --fq_1 "${in_se}" --qname \
        --fil_out "${t}.cram" 2> /dev/null
)" || true

exp="# Step 1: align the reads into the BAM work file.
bwa-mem2 mem \\
    -t 1 \\
    ${idx_bm2} \\
    ${in_se} \\
    | samtools view \\
        -@ 1 \\
        -F 4 \\
        -o ${t}.bam"
check_block "BWA-MEM2 SE dry run prints Step 1 without '-q'" \
    "${exp}" "$(print_step 1 "${out}")"

exp="# Step 5: write the final output and the retained files.
samtools view \\
    -T ${ref_fa} \\
    -O cram \\
    -o ${t}.cram \\
    ${t}.bam

samtools index \\
    -@ 1 \\
    ${t}.cram

samtools view \\
    -T ${ref_fa} \\
    -O cram \\
    -o ${t}.qnam.cram \\
    ${t}.qnam.bam

rm -f \\
    ${t}.bam \\
    ${t}.bam.bai

rm -f \\
    ${t}.qnam.bam"
check_block "CRAM dry run prints Step 5" \
    "${exp}" "$(print_step 5 "${out}")"


# A dry run still validates: a missing FASTQ is refused before any step.
rc=0
out="$(
    run_helper --dry_run --aligner bowtie2 --index "${idx_bt2}" \
        --fq_1 "${tmp}/absent.fastq.gz" --fil_out "${tmp}/absent.bam" 2>&1
)" || rc=$?

if [[
    "${rc}" -ne 0
    && "${out}" == *"align_fastqs): file associated with '--fq_1' does not"*
    && "${out}" != *"# Step "*
]]; then
    record_pass "helper dry run still refuses a missing '--fq_1'"
else
    record_fail "helper dry run did not refuse a missing '--fq_1' (rc ${rc})"
fi


# 'submit' prints each call with its log redirections, then the steps, and
# writes no logs or outputs.
dir_sub="${tmp}/sub"
mkdir -p "${dir_sub}/logs"
log="${dir_sub}/logs/align_fastqs.tiny_se"

rc=0
out="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/submit_align_fastqs.sh" \
        --dry_run \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --threads 1 \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --csv_fil_in "${in_se}" \
        --dir_out "${dir_sub}" \
        --dir_eo "${dir_sub}/logs" \
        --sfx_se ".atria.fastq.gz" \
        --sfx_pe "_R1.atria.fastq.gz" 2> /dev/null
)" || rc=$?

if [[
    "${rc}" -eq 0
    && "${out}" == *"# Sample 'tiny_se': call to 'align_fastqs'."*
    && "${out}" == *"    > ${log}.stdout.txt \\"$'\n'"    2> ${log}.stderr.txt"*
    && "${out}" == *"# Step 5: "*
]]; then
    record_pass "submit dry run prints the call, its logs, and the steps"
else
    record_fail "submit dry run output is incomplete (rc ${rc})"
fi

if [[
    -z "$(find "${dir_sub}" -mindepth 1 ! -path "${dir_sub}/logs" -print)"
]]; then
    record_pass "submit dry run writes no log or output"
else
    record_fail "submit dry run wrote files in '${dir_sub}'"
fi


# 'execute' prints its dispatch, then forwards the dry run through 'submit',
# for every entry of a mixed list, and writes nothing.
dir_exe="${tmp}/exe"
mkdir -p "${dir_exe}/logs"

rc=0
out="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_align_fastqs.sh" \
        --dry_run \
        --threads 1 \
        --csv_fil_in "${in_se};${in_pe_1},${in_pe_2}" \
        --dir_out "${dir_exe}" \
        --dir_eo "${dir_exe}/logs" \
        --max_job 1 \
        --aligner bowtie2 \
        --index "${idx_bt2}" 2> /dev/null
)" || rc=$?

if [[
    "${rc}" -eq 0
    && "${out}" == *"submit_align_fastqs.sh --env_nam "*
    && "${out}" == *"# Sample 'tiny_se': call to 'align_fastqs'."*
    && "${out}" == *"# Sample 'tiny_pe': call to 'align_fastqs'."*
    && "${out}" == *"-U ${in_se} \\"$'\n'"    | samtools view \\"*
    && "${out}" == *"-2 ${in_pe_2} \\"$'\n'"    | samtools view \\"*
]]; then
    record_pass "execute dry run prints the dispatch and each entry's steps"
else
    record_fail "execute dry run output is incomplete (rc ${rc})"
fi

if [[
    -z "$(find "${dir_exe}" -mindepth 1 ! -path "${dir_exe}/logs" -print)"
]]; then
    record_pass "execute dry run writes no log or output"
else
    record_fail "execute dry run wrote files in '${dir_exe}'"
fi


# Every layer accepts '--dry-run' and refuses the retired '--dry'.
for spl in --dry-run --dry; do
    declare -A out_lyr=()
    declare -A rc_lyr=( [helper]=0 [submit]=0 [execute]=0 )

    out_lyr[helper]="$(
        run_helper "${spl}" --aligner bowtie2 --index "${idx_bt2}" \
            --fq_1 "${in_se}" --fil_out "${tmp}/spell.bam" 2>&1
    )" || rc_lyr[helper]=$?
    out_lyr[submit]="$(
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_align_fastqs.sh" \
            "${spl}" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --index "${idx_bt2}" \
            --csv_fil_in "${in_se}" \
            --dir_out "${dir_sub}" \
            --dir_eo "${dir_sub}/logs" \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" 2>&1
    )" || rc_lyr[submit]=$?
    out_lyr[execute]="$(
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_align_fastqs.sh" \
            "${spl}" \
            --threads 1 \
            --csv_fil_in "${in_se}" \
            --dir_out "${dir_exe}" \
            --dir_eo "${dir_exe}/logs" \
            --max_job 1 \
            --index "${idx_bt2}" 2>&1
    )" || rc_lyr[execute]=$?

    for lyr in helper submit execute; do
        if [[ "${spl}" == "--dry-run" ]]; then
            if [[
                "${rc_lyr[${lyr}]}" -eq 0
                && "${out_lyr[${lyr}]}" == *"# Step 1: "*
            ]]; then
                record_pass "${lyr} accepts '--dry-run'"
            else
                record_fail "${lyr} did not accept '--dry-run'"
            fi
        elif [[
            "${rc_lyr[${lyr}]}" -ne 0
            && "${out_lyr[${lyr}]}" == *"'--dry'"*
            && "${out_lyr[${lyr}]}" != *"# Step 1: "*
        ]]; then
            record_pass "${lyr} refuses '--dry'"
        else
            record_fail "${lyr} did not refuse '--dry'"
        fi
    done
done


finish
