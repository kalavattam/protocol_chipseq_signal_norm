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


# Drop the pipeline record's lines from a dry run: some, such as the start time
# and the commit, change between runs, so 'test_align_fastqs_record.sh' checks
# them instead.
function drop_record_lines() {
    grep -v -F "    \$'@CO" <<< "${1}" || true
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
    run_helper \
        --dry_run \
        --threads 2 \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --fq_1 "${in_pe_1}" \
        --fq_2 "${in_pe_2}" \
        --req_flg \
        --qname \
        --fil_out "${dir_sp}/tiny pe.bam" 2> /dev/null
)" || rc=$?

exp="
# Step 1: align the reads into the BAM work file.
bowtie2 \\
    -p 2 \\
    --reorder \\
    --seed 0 \\
    -x ${idx_bt2} \\
    --very-sensitive \\
    --no-unal \\
    --phred33 \\
    --no-mixed \\
    --no-discordant \\
    -X 1000 \\
    -1 ${in_pe_1} \\
    -2 ${in_pe_2} \\
    | samtools view \\
        -@ 2 \\
        -F 12 \\
        -f 2 \\
        -O bam,level=1 \\
        -o ${esc}.bam


# Step 2: sort the BAM work file by queryname, fix mates in paired-end data,
# and filter templates.
samtools sort \\
    -@ 2 \\
    -n \\
    -u \\
    -o - \\
    ${esc}.bam \\
    | samtools fixmate \\
        -@ 2 \\
        -c \\
        -m \\
        -O bam,level=1 \\
        - \\
        ${esc}.qnam.bam

samtools view \\
    -@ 2 \\
    -F 2308 \\
    -f 2 \\
    -q 1 \\
    ${esc}.qnam.bam \\
    | cut -f 1 \\
    | uniq -c \\
    | awk -v n=2 \\
        \\\$1\\ ==\\ n\\ \\{\\ print\\ \\\$2\\ \\} \\
        > ${esc}.mapq.names.txt

samtools view \\
    -@ 2 \\
    -O bam,level=1 \\
    -N ${esc}.mapq.names.txt \\
    -o ${esc}.mapq.bam \\
    ${esc}.qnam.bam

rm -f \\
    ${esc}.mapq.names.txt

mv -f \\
    ${esc}.mapq.bam \\
    ${esc}.qnam.bam


# Step 3: sort the BAM work file by coordinates.
samtools sort \\
    -@ 2 \\
    -O bam,level=1 \\
    -o ${esc}.coor.bam \\
    ${esc}.qnam.bam

rm -f \\
    ${esc}.qnam.bam

mv -f \\
    ${esc}.coor.bam \\
    ${esc}.bam


# Step 4: mark duplicates in the BAM work file. '-d 2500' applies only if all
# read names carry flow-cell coordinates of one format.
samtools markdup \\
    -@ 2 \\
    -t \\
    -d 2500 \\
    -f ${esc}.markdup.txt \\
    ${esc}.bam \\
    ${esc}.mark.bam

printf 'OPTICAL DISTANCE: %s\\n' 2500 >> ${esc}.markdup.txt

gzip -f \\
    ${esc}.markdup.txt

mv -f \\
    ${esc}.mark.bam \\
    ${esc}.bam


# Step 5: write the pipeline record into the BAM work file.
samtools view \\
    -H \\
    --no-PG \\
    ${esc}.bam \\
    > ${esc}.header.sam

printf %s\\\\n \\
    >> ${esc}.header.sam

samtools reheader \\
    ${esc}.header.sam \\
    ${esc}.bam \\
    > ${esc}.record.bam

mv -f \\
    ${esc}.record.bam \\
    ${esc}.bam

rm -f \\
    ${esc}.header.sam

samtools index \\
    -@ 2 \\
    ${esc}.bam


# Step 6: write the final output and the retained files.
samtools sort \\
    -@ 2 \\
    -n \\
    -o ${esc}.qnam.bam \\
    ${esc}.bam"

if [[ "${rc}" -eq 0 ]]; then
    record_pass "helper dry run returns 0"
else
    record_fail "helper dry run returned ${rc}"
fi

check_block \
    "Bowtie 2 PE dry run prints every step, spaces escaped" \
    "${exp}" \
    "$(drop_record_lines "${out}")"

if [[ -z "$(find "${dir_sp}" -mindepth 1 -print -quit)" ]]; then
    record_pass "helper dry run writes no file"
else
    record_fail "helper dry run wrote files in '${dir_sp}'"
fi


# Bowtie 2, paired-end: '--bt2_X' sets '-X', and the hidden '--bt2_inherited'
# restores the earlier call token for token, with no '-X'.
t="${tmp}/bt2"
for arm in "--bt2_X 800|-X 800" "-2i|--no-overlap;--no-dovetail"; do
    IFS='|' read -r arg tok <<< "${arm}"
    read -r -a args <<< "${arg}"
    IFS=';' read -r -a toks <<< "${tok}"
    printf -v lin '    %s \\\n' "${toks[@]}"

    out="$(
        run_helper \
            --dry_run \
            --threads 1 \
            --aligner bowtie2 \
            --index "${idx_bt2}" \
            --fq_1 "${in_pe_1}" \
            --fq_2 "${in_pe_2}" \
            "${args[@]}" \
            --fil_out "${t}.bam" 2> /dev/null
    )" || true

    exp="# Step 1: align the reads into the BAM work file.
bowtie2 \\
    -p 1 \\
    --reorder \\
    --seed 0 \\
    -x ${idx_bt2} \\
    --very-sensitive \\
    --no-unal \\
    --phred33 \\
    --no-mixed \\
    --no-discordant \\
${lin}    -1 ${in_pe_1} \\
    -2 ${in_pe_2} \\
    | samtools view \\
        -@ 1 \\
        -F 12 \\
        -O bam,level=1 \\
        -o ${t}.bam"

    check_block \
        "Bowtie 2 PE dry run with '${arg}' prints its call" \
        "${exp}" \
        "$(print_step 1 "${out}")"
done


# Bowtie 2, paired-end, the hidden '--bt2_relaxed' ('-2r'): '-I 0 -X <bt2_X>'
# without the pairing options, with the preset from '--bt2_mode'. Bowtie 2
# keeps '--no-unal' unless '--unal_mode separate' ('-um') keeps the unaligned
# reads, which Step 1 then writes with '-U'. The record gives each setting.
t="${tmp}/rlx"
for arm in \
    "-2r|--very-sensitive|--no-unal|1000|" \
    "-2r --bt2_X 2000|--very-sensitive|--no-unal|2000|" \
    "-2r --bt2_mode local|--very-sensitive-local|--no-unal|1000|" \
    "-2r -um separate|--very-sensitive||1000|-U"
do
    IFS='|' read -r arg pst unl x_len set_a <<< "${arm}"
    read -r -a args <<< "${arg}"
    lin_unl=""
    if [[ -n "${unl}" ]]; then lin_unl="    ${unl} \\"$'\n'; fi
    lin_u=""
    if [[ -n "${set_a}" ]]; then
        lin_u="        -U ${t}.setU.bam \\"$'\n'
    fi

    out="$(
        run_helper \
            --dry_run \
            --threads 1 \
            --aligner bowtie2 \
            --index "${idx_bt2}" \
            --fq_1 "${in_pe_1}" \
            --fq_2 "${in_pe_2}" \
            "${args[@]}" \
            --fil_out "${t}.bam" 2> /dev/null
    )" || true

    exp="# Step 1: align the reads into the BAM work file.
bowtie2 \\
    -p 1 \\
    --reorder \\
    --seed 0 \\
    -x ${idx_bt2} \\
    ${pst} \\
${lin_unl}    --phred33 \\
    -I 0 \\
    -X ${x_len} \\
    -1 ${in_pe_1} \\
    -2 ${in_pe_2} \\
    | samtools view \\
        -@ 1 \\
        -F 12 \\
${lin_u}        -O bam,level=1 \\
        -o ${t}.bam"

    check_block \
        "Bowtie 2 PE dry run with '${arg}' prints the relaxed call" \
        "${exp}" \
        "$(print_step 1 "${out}")"

    if [[
        "${out}" == *"pipeline: --bt2_X=${x_len}'"*
            && "${out}" == *"pipeline: --bt2_relaxed=true'"*
            && "${out}" == *"pipeline: warning=not for signal or scaling'"*
    ]]; then
        record_pass "Bowtie 2 PE dry run with '${arg}' records the settings"
    else
        record_fail "Bowtie 2 PE dry run with '${arg}' lacks its record lines"
    fi
done

# With '--req_flg', '--unal_mode remove' requires proper pairs in Step 1, as
# before, and again by template after 'fixmate' in Step 2. 'separate' requires
# them only in Step 2, after the routing, and keeps what it rejects. Its Step 4
# splits the diagnostic files out without indexing any file.
t="${tmp}/req"
for mod in remove separate; do
    out="$(
        run_helper \
            --dry_run \
            --threads 1 \
            --aligner bowtie2 \
            --index "${idx_bt2}" \
            --fq_1 "${in_pe_1}" \
            --fq_2 "${in_pe_2}" \
            --req_flg \
            --unal_mode "${mod}" \
            --fil_out "${t}.${mod}.bam" 2> /dev/null
    )" || true
    st1="$(print_step 1 "${out}")"
    st2="$(print_step 2 "${out}")"
    st4="$(print_step 4 "${out}")"

    if [[ "${mod}" == "remove" ]]; then
        n_f2="$(grep -c -- '-f 2' <<< "${st2}" || true)"
        if [[
            "${st1}" == *$'    --no-unal \\\n'*
                && "${st1}" == *$'        -F 12 \\\n        -f 2 \\\n'*
                && "${n_f2}" == 1
                && "${st2}" == *$'    -F 2308 \\\n    -f 2 \\\n    -q 1'*
                && "${out}" != *"setU"*
                && "${out}" != *".unal."*
                && "${out}" != *".diag."*
        ]]; then
            record_pass "dry run with '--req_flg' and 'remove' is unchanged"
        else
            record_fail "dry run with '--req_flg' and 'remove' changed"
        fi
    elif [[
        "${st1}" != *"--no-unal"*
            && "${st1}" != *"-f 2"*
            && "${st1}" == *"-U ${t}.separate.setU.bam"*
            && "${st2}" == *"-N ${t}.separate.unal.names.txt"*
            && "${st2}" == *$'    -f 2 \\\n    -U '"${t}.separate.setI0.bam"*
            && "${st4}" == *"-o ${t}.separate.unal.bam"*
            && "${st4}" == *"-o ${t}.separate.improp.bam"*
            && "${st4}" != *"samtools index"*
    ]]; then
        record_pass "dry run with '--req_flg' and 'separate' moves '-f 2'"
    else
        record_fail "dry run with '--req_flg' and 'separate' is wrong"
    fi
done


# With '--unal_mode separate' and '--req_flg', Step 2 sorts Set U to disk and
# lists its names. It then streams the sorted Set M into the routing, the
# seeded merge into Set U into 'fixmate', and the '-f 2' view into Set M's
# 'fixmate'. The view keeps the records it rejects.
#
# Next, Step 2 filters Set M by template and keeps the templates it drops. It
# fixes mates in the rejected records and sorts the rejects into the improper
# and low-MAPQ work files, with a seeded merge again. Last, it lists every
# diagnostic work file's names. Every file it writes is at level 1.
t="${tmp}/sep"
out="$(
    run_helper \
        --dry_run \
        --threads 1 \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --fq_1 "${in_pe_1}" \
        --fq_2 "${in_pe_2}" \
        --req_flg \
        --unal_mode separate \
        --fil_out "${t}.bam" 2> /dev/null
)" || true

exp="# Step 2: sort the BAM work file by queryname, fix mates in paired-end data,
# and filter templates.
samtools sort \\
    -@ 1 \\
    -n \\
    -O bam,level=1 \\
    -o ${t}.setU.qnam.bam \\
    ${t}.setU.bam

mv -f \\
    ${t}.setU.qnam.bam \\
    ${t}.setU.bam

samtools view \\
    -@ 1 \\
    ${t}.setU.bam \\
    | cut -f 1 \\
    | uniq \\
        > ${t}.unal.names.txt

samtools sort \\
    -@ 1 \\
    -n \\
    -u \\
    -o - \\
    ${t}.bam \\
    | samtools view \\
        -@ 1 \\
        -O bam,level=1 \\
        -N ${t}.unal.names.txt \\
        -o ${t}.moved.bam \\
        -U ${t}.qnam.bam \\
        -

samtools merge \\
    -@ 1 \\
    -n \\
    -s 0 \\
    -u \\
    -o - \\
    ${t}.setU.bam \\
    ${t}.moved.bam \\
    | samtools fixmate \\
        -@ 1 \\
        -c \\
        -m \\
        -O bam,level=1 \\
        - \\
        ${t}.setU.mate.bam

mv -f \\
    ${t}.setU.mate.bam \\
    ${t}.setU.bam

rm -f \\
    ${t}.moved.bam

samtools view \\
    -@ 1 \\
    -u \\
    -f 2 \\
    -U ${t}.setI0.bam \\
    ${t}.qnam.bam \\
    | samtools fixmate \\
        -@ 1 \\
        -c \\
        -m \\
        -O bam,level=1 \\
        - \\
        ${t}.mate.bam

mv -f \\
    ${t}.mate.bam \\
    ${t}.qnam.bam

samtools view \\
    -@ 1 \\
    -F 2308 \\
    -f 2 \\
    -q 1 \\
    ${t}.qnam.bam \\
    | cut -f 1 \\
    | uniq -c \\
    | awk -v n=2 \\
        \\\$1\\ ==\\ n\\ \\{\\ print\\ \\\$2\\ \\} \\
        > ${t}.mapq.names.txt

samtools view \\
    -@ 1 \\
    -O bam,level=1 \\
    -N ${t}.mapq.names.txt \\
    -o ${t}.mapq.bam \\
    -U ${t}.setR.bam \\
    ${t}.qnam.bam

rm -f \\
    ${t}.mapq.names.txt

mv -f \\
    ${t}.mapq.bam \\
    ${t}.qnam.bam

samtools fixmate \\
    -@ 1 \\
    -c \\
    -m \\
    -O bam,level=1 \\
    ${t}.setI0.bam \\
    ${t}.setI0.mate.bam

mv -f \\
    ${t}.setI0.mate.bam \\
    ${t}.setI0.bam

samtools view \\
    -@ 1 \\
    -F 2308 \\
    -f 2 \\
    ${t}.setR.bam \\
    | cut -f 1 \\
    | uniq -c \\
    | awk \\
        \\\$1\\ ==\\ 2\\ \\{\\ print\\ \\\$2\\ \\} \\
        > ${t}.setQ.names.txt

samtools view \\
    -@ 1 \\
    -O bam,level=1 \\
    -N ${t}.setQ.names.txt \\
    -o ${t}.setQ.bam \\
    -U ${t}.setI1.bam \\
    ${t}.setR.bam

samtools merge \\
    -@ 1 \\
    -n \\
    -s 0 \\
    -O bam,level=1 \\
    -o ${t}.setI.bam \\
    ${t}.setI0.bam \\
    ${t}.setI1.bam

rm -f \\
    ${t}.setR.bam \\
    ${t}.setI0.bam \\
    ${t}.setI1.bam

samtools view \\
    -@ 1 \\
    ${t}.setI.bam \\
    | cut -f 1 \\
    | uniq \\
        > ${t}.setI.names.txt

cat \\
    ${t}.unal.names.txt \\
    ${t}.setQ.names.txt \\
    ${t}.setI.names.txt \\
    > ${t}.diag.names.txt"

check_block \
    "PE dry run with 'separate' and '--req_flg' prints Step 2" \
    "${exp}" \
    "$(print_step 2 "${out}")"

# Without '--req_flg', Set M's 'fixmate' reads the sorted file itself. In
# single-end data, Step 2 lists Set U's names, then sorts Set M to disk.
out="$(
    run_helper \
        --dry_run \
        --threads 1 \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --fq_1 "${in_pe_1}" \
        --fq_2 "${in_pe_2}" \
        --unal_mode separate \
        --fil_out "${t}.bam" 2> /dev/null
)" || true
exp="samtools fixmate \\
    -@ 1 \\
    -c \\
    -m \\
    -O bam,level=1 \\
    ${t}.qnam.bam \\
    ${t}.mate.bam"

if [[
    "$(print_step 2 "${out}")" == *"${exp}"*
        && "$(print_step 2 "${out}")" != *"-f 2"*
]]; then
    record_pass "PE dry run with 'separate' fixes Set M's mates from its file"
else
    record_fail "PE dry run with 'separate' fixes Set M's mates wrongly"
fi

out="$(
    run_helper \
        --dry_run \
        --threads 1 \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --fq_1 "${in_se}" \
        --unal_mode separate \
        --fil_out "${t}_se.bam" 2> /dev/null
)" || true
exp="> ${t}_se.unal.names.txt

samtools sort \\
    -@ 1 \\
    -n \\
    -O bam,level=1 \\
    -o ${t}_se.qnam.bam \\
    ${t}_se.bam

samtools view \\
    -@ 1 \\
    -F 2308 \\
    -q 1 \\"
st2="$(print_step 2 "${out}")"
lin_mov="mv -f \\
    ${t}_se.setR.bam \\
    ${t}_se.setQ.bam"

if [[
    "${st2}" == *"${exp}"*
        && "${st2}" != *"fixmate"*
        && "${st2}" != *"moved"*
        && "${st2}" != *"setI"*
        && "${st2}" == *"${lin_mov}"*
]]; then
    record_pass "SE dry run with 'separate' sorts Set M to disk, no routing"
else
    record_fail "SE dry run with 'separate' prints Step 2 wrongly"
fi


# Only Step 5 indexes a BAM file, once per file it writes, and not for CRAM
# output, which Step 6 indexes instead. A file that a later step rewrites is
# written at level 1; a file that ships is written at the default level. Under
# 'separate', Step 3 also sorts the diagnostic work files. Step 4 marks them
# with the main file (Markdup 2) and splits them into three files.
t="${tmp}/mtx"
for mod in remove separate; do
    for fmt in bam cram; do
        lbl="'${mod}' with ${fmt^^} output"
        out="$(
            run_helper \
                --dry_run \
                --threads 1 \
                --aligner bowtie2 \
                --index "${idx_bt2}" \
                --ref_fa "${ref_fa}" \
                --fq_1 "${in_pe_1}" \
                --fq_2 "${in_pe_2}" \
                --unal_mode "${mod}" \
                --qname \
                --fil_out "${t}.${mod}.${fmt}" 2> /dev/null
        )" || true
        st3="$(print_step 3 "${out}")"
        st4="$(print_step 4 "${out}")"
        st5="$(print_step 5 "${out}")"
        st6="$(print_step 6 "${out}")"
        n_fil=1
        if [[ "${mod}" == "separate" ]]; then n_fil=4; fi
        n_idx5="$(grep -c '^samtools index' <<< "${st5}" || true)"
        n_idx6="$(grep -c '^samtools index' <<< "${st6}" || true)"
        n_lvl4="$(grep -c 'level=1' <<< "${st4}" || true)"
        n_lvl6="$(grep -c 'level=1' <<< "${st6}" || true)"
        n_b4="$(grep -c '^    -b \\$' <<< "${st4}" || true)"

        exp_idx5="${n_fil}"
        exp_lvl3=1
        exp_lvl4=0
        exp_lvl6=0
        exp_b4=0
        if [[ "${fmt}" == "cram" ]]; then
            exp_idx5=0
            exp_lvl6="${n_fil}"
        fi
        if [[ "${mod}" == "separate" && "${fmt}" == "bam" ]]; then
            exp_lvl3=2
            exp_lvl4=1
            exp_b4=2
        elif [[ "${mod}" == "separate" ]]; then
            exp_lvl3=2
            exp_lvl4=4
        elif [[ "${fmt}" == "cram" ]]; then
            exp_lvl4=1
        fi

        if [[
            -n "${st3}" && -n "${st4}"
                && "${st3}${st4}" != *"samtools index"*
                && "${n_idx5}" == "${exp_idx5}"
        ]]; then
            record_pass "${lbl}: only Step 5 indexes, ${exp_idx5} file(s)"
        else
            record_fail "${lbl}: indexes in Steps 3 to 5 are wrong (${n_idx5})"
        fi

        if [[
            "${fmt}" == "bam" && "${n_idx6}" == 0
                || "${fmt}" == "cram" && "${n_idx6}" == "${n_fil}"
        ]]; then
            record_pass "${lbl}: Step 6 indexes ${n_idx6} file(s)"
        else
            record_fail "${lbl}: Step 6 indexes ${n_idx6} file(s)"
        fi

        if [[
            "${n_lvl4}" == "${exp_lvl4}" && "${n_b4}" == "${exp_b4}"
                && "${n_lvl6}" == "${exp_lvl6}"
                && "$(grep -c 'level=1' <<< "${st3}")" == "${exp_lvl3}"
        ]]; then
            record_pass "${lbl}: work files at level 1, shipped files not"
        else
            record_fail "${lbl}: levels wrong: Step 4 ${n_lvl4}," \
                "'-b' ${n_b4}, Step 6 ${n_lvl6}"
        fi
    done
done

# BWA-backtrack, paired-end, '--mapq 20': the '.sai' redirects, their removal,
# and the template filter's pipeline.
t="${tmp}/aln"
out="$(
    run_helper \
        --dry_run \
        --threads 1 \
        --aligner bwa \
        --bwa_alg aln \
        --mapq 20 \
        --index "${idx_bwa}" \
        --fq_1 "${in_pe_1}" \
        --fq_2 "${in_pe_2}" \
        --fil_out "${t}.bam" 2> /dev/null
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
        -O bam,level=1 \\
        -o ${t}.bam

rm -f \\
    ${t}.R1.sai \\
    ${t}.R2.sai"

check_block \
    "BWA-backtrack PE dry run prints Step 1" \
    "${exp}" \
    "$(print_step 1 "${out}")"

exp="# Step 2: sort the BAM work file by queryname, fix mates in paired-end data,
# and filter templates.
samtools sort \\
    -@ 1 \\
    -n \\
    -u \\
    -o - \\
    ${t}.bam \\
    | samtools fixmate \\
        -@ 1 \\
        -c \\
        -m \\
        -O bam,level=1 \\
        - \\
        ${t}.qnam.bam

samtools view \\
    -@ 1 \\
    -F 2308 \\
    -q 20 \\
    ${t}.qnam.bam \\
    | cut -f 1 \\
    | uniq -c \\
    | awk -v n=2 \\
        \\\$1\\ ==\\ n\\ \\{\\ print\\ \\\$2\\ \\} \\
        > ${t}.mapq.names.txt

samtools view \\
    -@ 1 \\
    -O bam,level=1 \\
    -N ${t}.mapq.names.txt \\
    -o ${t}.mapq.bam \\
    ${t}.qnam.bam

rm -f \\
    ${t}.mapq.names.txt

mv -f \\
    ${t}.mapq.bam \\
    ${t}.qnam.bam"

check_block \
    "BWA-backtrack PE dry run prints the template filter" \
    "${exp}" \
    "$(print_step 2 "${out}")"

# Without '--qname', BAM output is the work file itself: Step 6 runs nothing.
exp="# Step 6: write the final output and the retained files.
${msg_nil}"

check_block \
    "BAM dry run without '--qname' runs nothing in Step 6" \
    "${exp}" \
    "$(print_step 6 "${out}")"

# The printed commands are the ones a real run executes: run as a script, they
# write the same files, holding the same records.
dir_par="${tmp}/parity"
mkdir -p "${dir_par}/printed" "${dir_par}/real"
out="$(
    run_helper \
        --dry_run \
        --threads 1 \
        --aligner bwa \
        --bwa_alg aln \
        --mapq 20 \
        --index "${idx_bwa}" \
        --fq_1 "${in_pe_1}" \
        --fq_2 "${in_pe_2}" \
        --fil_out "${dir_par}/printed/aln.bam" \
        2> /dev/null
)" || true

rc=0
"${TEST_BASH}" -e -c "${out}" > /dev/null 2>&1 || rc=$?
run_helper \
    --threads 1 \
    --aligner bwa \
    --bwa_alg aln \
    --mapq 20 \
    --index "${idx_bwa}" \
    --fq_1 "${in_pe_1}" \
    --fq_2 "${in_pe_2}" \
    --fil_out "${dir_par}/real/aln.bam" > /dev/null 2>&1 || rc=$?

if [[
    "${rc}" -eq 0
        && "$(
            cd "${dir_par}/printed" && find . -type f | LC_ALL=C sort
        )" == "$(cd "${dir_par}/real" && find . -type f | LC_ALL=C sort)"
]] && cmp -s \
    <(run_samtools view "${dir_par}/printed/aln.bam" | cut -f 1-11) \
    <(run_samtools view "${dir_par}/real/aln.bam" | cut -f 1-11)
then
    record_pass "printed commands, run as a script, match a real run"
else
    record_fail "printed commands, run as a script, differ from a real run"
fi


# BWA-MEM2, single-end, CRAM with '--qname': Step 6 name-sorts the final BAM,
# converts both files, indexes the output, and removes the BAM work files.
t="${tmp}/cram"
out="$(
    run_helper \
        --dry_run \
        --threads 1 \
        --aligner bwa-mem2 \
        --index "${idx_bm2}" \
        --ref_fa "${ref_fa}" \
        --fq_1 "${in_se}" \
        --qname \
        --fil_out "${t}.cram" 2> /dev/null
)" || true

exp="# Step 1: align the reads into the BAM work file.
bwa-mem2 mem \\
    -t 1 \\
    -K 100000000 \\
    ${idx_bm2} \\
    ${in_se} \\
    | samtools view \\
        -@ 1 \\
        -F 4 \\
        -O bam,level=1 \\
        -o ${t}.bam"

check_block \
    "BWA-MEM2 SE dry run prints Step 1 without '-q'" \
    "${exp}" \
    "$(print_step 1 "${out}")"

exp="# Step 6: write the final output and the retained files.
samtools sort \\
    -@ 1 \\
    -n \\
    -O bam,level=1 \\
    -o ${t}.qnam.bam \\
    ${t}.bam

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

check_block \
    "CRAM dry run prints Step 6" \
    "${exp}" \
    "$(print_step 6 "${out}")"


# A dry run still validates: a missing FASTQ is refused before any step.
rc=0
out="$(
    run_helper \
        --dry_run \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --fq_1 "${tmp}/absent.fastq.gz" \
        --fil_out "${tmp}/absent.bam" 2>&1
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
        && "${out}" \
            == *"    > ${log}.stdout.txt \\"$'\n'"    2> ${log}.stderr.txt"*
        && "${out}" == *"# Step 6: "*
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
        run_helper \
            "${spl}" \
            --aligner bowtie2 \
            --index "${idx_bt2}" \
            --fq_1 "${in_se}" \
            --fil_out "${tmp}/spell.bam" 2>&1
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
