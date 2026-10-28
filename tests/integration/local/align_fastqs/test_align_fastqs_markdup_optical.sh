#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_align_fastqs_markdup_optical.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="align-fastqs optical duplicates"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths.
#
# Step 4 passes '-d 2500' to 'samtools markdup' only when every read name
# carries flow-cell coordinates. Samtools then tags a duplicate 'dt:Z:SQ'
# (optical) when it lies on the same tile as another read of its duplicate set
# and within 2500 pixels of it. Every other duplicate is tagged 'dt:Z:LB'
# (library). 'LB' covers PCR copies of one molecule and separate molecules that
# happen to share both ends (Samtools cannot tell them apart).
#
# The statistics go to '<stem>.markdup.txt.gz'. The reads are the 'opt_*' sets
# that 'make.sh' writes, and the fixture README lists their names.
dir_fx="${ROOT_REPO}/tests/fixtures/align_fastqs"
ref_rpt="${dir_fx}/reference/repeat.fa"
idx_bt2="${dir_fx}/bowtie2/repeat"
idx_bwa="${dir_fx}/bwa/repeat.fa"
dir_fq_pe="${dir_fx}/fastq/pe"
dir_fq_se="${dir_fx}/fastq/se"

tmp="${TEST_DIR_TMP}/align_fastqs_markdup_optical"

# Tile prefixes of the expected names: one tile, another tile.
pfx="LH00001:168:22ABCDEF:1:1101"
pfx_tile="LH00001:168:22ABCDEF:1:1102"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${ref_rpt}" \
    "${ref_rpt}.fai" \
    "${idx_bt2}.1.bt2" \
    "${idx_bwa}.bwt" \
    "${dir_fx}/sam/optical_names.sam" \
    "${dir_fq_se}/opt_none_se.atria.fastq.gz" || {
    finish
    exit $?
}


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
    ' _ "${ROOT_REPO}" --threads 1 "$@" 2>&1
}


# Print 'name<TAB>dt' for each primary first-mate or single-end record, with
# '-' for a record without a 'dt' tag, sorted by name.
function print_dt() {
    run_samtools view -F 0x980 -T "${ref_rpt}" "${1}" \
        | awk -F '\t' '{
            d = "-"
            for (i = 12; i <= NF; i++) if ($i ~ /^dt:Z:/) d = substr($i, 6)
            print $1 "\t" d
        }' \
        | LC_ALL=C sort
}


# Print the '@PG' command line of 'samtools markdup'.
function print_pg_markdup() {
    run_samtools view -H -T "${ref_rpt}" "${1}" \
        | grep '^@PG' \
        | grep -o 'CL:samtools markdup[^'$'\t'']*' \
        || true
}


# Align one case into its own directory and check the exit status; print the
# helper's combined output on failure. Sets 'out' and 'dir_out'.
function run_case() {
    local lbl="${1}"
    local nam="${2}"
    shift 2
    local rc=0

    dir_out="${tmp}/${nam}"
    mkdir -p "${dir_out}"
    out="$(run_helper "$@")" || rc=$?

    if [[ "${rc}" -eq 0 ]]; then
        record_pass "${lbl}: exit 0"
    else
        record_fail "${lbl}: exit ${rc}"
        printf '%s\n' "${out}" >&2
        return 1
    fi
}


# Check that a file's 'dt' tags match the expected 'name<TAB>dt' lines.
function check_dt() {
    local lbl="${1}"
    local fil="${2}"
    local exp="${3}"
    local got

    got="$(print_dt "${fil}")"
    exp="$(printf '%s\n' "${exp}" | LC_ALL=C sort)"

    if [[ "${got}" == "${exp}" ]]; then
        record_pass "${lbl}: 'dt' tags as expected"
    else
        record_fail "${lbl}: 'dt' tags differ"
        diff <(printf '%s\n' "${exp}") <(printf '%s\n' "${got}") >&2 || true
    fi
}


# Check that a case ran without '-d': no '-d' in '@PG', no 'dt' tag, and the
# given warning printed once.
function check_no_optical() {
    local lbl="${1}"
    local fil="${2}"
    local wrn="${3}"
    local n_wrn

    if \
        [[ -n "$(print_pg_markdup "${fil}")" ]] \
        && ! print_pg_markdup "${fil}" | grep -q -- ' -d '
    then
        record_pass "${lbl}: markdup ran without '-d'"
    else
        record_fail "${lbl}: markdup '@PG': '$(print_pg_markdup "${fil}")'"
    fi

    if ! \
        print_dt "${fil}" | cut -f 2 | grep -qv '^-$'
    then
        record_pass "${lbl}: no 'dt' tags"
    else
        record_fail "${lbl}: found 'dt' tags"
    fi

    n_wrn="$(grep -c -- "${wrn}" <<< "${out}" || true)"

    if [[ "${n_wrn}" -eq 1 ]]; then
        record_pass "${lbl}: warns once"
    else
        record_fail "${lbl}: warning printed ${n_wrn} times, expected once"
        printf '%s\n' "${out}" >&2
    fi
}


# Print a value from a gzip-compressed 'samtools markdup' statistics file or
# nothing when the file is missing, so the check that reads it fails and the
# test goes on.
function stat_of() {
    { gzip -cd "${1}" 2> /dev/null || true; } \
        | awk -F ': ' -v key="${2}" '$1 == key { print $2 }'
}


# Check the statistics' last line: 'OPTICAL DISTANCE: ' and the expected value.
function check_dist() {
    local lbl="${1}"
    local fil="${2}"
    local exp="OPTICAL DISTANCE: ${3}"
    local got

    got="$({ gzip -cd "${fil}" 2> /dev/null || true; } | tail -n 1)"

    if [[ "${got}" == "${exp}" ]]; then
        record_pass "${lbl}: statistics end with '${exp}'"
    else
        record_fail "${lbl}: statistics end with '${got}', expected '${exp}'"
    fi
}


# Paired-end, coordinate names. Against the original on tile 1101, reads 30 and
# 2000 pixels away on its tile are optical; reads 6000 away or on another tile
# are library duplicates.
lbl="PE coordinate names"
fq_1="${dir_fq_pe}/opt_pe_pe_R1.atria.fastq.gz"
fq_2="${dir_fq_pe}/opt_pe_pe_R2.atria.fastq.gz"
exp_pe="${pfx}:10000:10000	-
${pfx}:10030:10000	SQ
${pfx}:12000:10000	SQ
${pfx}:16000:10000	LB
${pfx_tile}:10000:10000	LB"

if \
    run_case \
        "${lbl}" \
        pe \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --fq_1 "${fq_1}" \
        --fq_2 "${fq_2}" \
        --qname \
        --fil_out "${tmp}/pe/x.bam"
then
    check_dt "${lbl}" "${dir_out}/x.bam" "${exp_pe}"

    if \
        print_pg_markdup "${dir_out}/x.bam" | grep -q -- ' -d 2500 '
    then
        record_pass "${lbl}: markdup ran with '-d 2500'"
    else
        record_fail "${lbl}: markdup '@PG' lacks '-d 2500'"
    fi

    # The statistics: gzip-compressed, two optical pairs.
    txt_mrk="${dir_out}/x.markdup.txt.gz"

    if \
        gzip -t "${txt_mrk}" 2> /dev/null
    then
        record_pass "${lbl}: statistics written gzip-compressed"
    else
        record_fail "${lbl}: no gzip-compressed '${txt_mrk}'"
    fi

    n_opt="$(stat_of "${txt_mrk}" "DUPLICATE PAIR OPTICAL")"
    n_dup="$(stat_of "${txt_mrk}" "DUPLICATE PAIR")"

    if [[ "${n_opt}" == "4" && "${n_dup}" == "8" ]]; then
        record_pass "${lbl}: statistics count 4 optical of 8 duplicate reads"
    else
        record_fail \
            "${lbl}: statistics count ${n_opt:-?} optical of ${n_dup:-?}" \
            "duplicate reads, expected 4 of 8"
    fi

    check_dist \
        "${lbl}" \
        "${txt_mrk}" \
        "2500"

    # The '--qname' copy carries the same 'dt' tags.
    check_dt "${lbl}, '--qname' copy" "${dir_out}/x.qnam.bam" "${exp_pe}"
fi


# Single-end, the same rule.
lbl="SE coordinate names"
fq_1="${dir_fq_se}/opt_se_se.atria.fastq.gz"

if \
    run_case \
        "${lbl}" \
        se \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --fq_1 "${fq_1}" \
        --fil_out "${tmp}/se/x.bam"
then
    check_dt "${lbl}" "${dir_out}/x.bam" "${pfx}:20000:20000	-
${pfx}:20030:20000	SQ
${pfx}:26000:20000	LB"

    n_opt="$(
        stat_of "${dir_out}/x.markdup.txt.gz" "DUPLICATE SINGLE OPTICAL"
    )"

    if [[ "${n_opt}" == "1" ]]; then
        record_pass "${lbl}: statistics count 1 optical read"
    else
        record_fail "${lbl}: statistics count ${n_opt:-?} optical reads"
    fi

    check_dist \
        "${lbl}" \
        "${dir_out}/x.markdup.txt.gz" \
        "2500"
fi


# SRA-style names: no coordinates, so no '-d'; flags as with coordinates.
lbl="PE SRA-style names"
fq_1="${dir_fq_pe}/opt_sra_pe_R1.atria.fastq.gz"
fq_2="${dir_fq_pe}/opt_sra_pe_R2.atria.fastq.gz"

if \
    run_case \
        "${lbl}" \
        sra \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --fq_1 "${fq_1}" \
        --fq_2 "${fq_2}" \
        --fil_out "${tmp}/sra/x.bam"
then
    check_no_optical \
        "${lbl}" \
        "${dir_out}/x.bam" \
        "5 of 5 read names carry no flow-cell coordinates"

    check_dist \
        "${lbl}" \
        "${dir_out}/x.markdup.txt.gz" \
        "not used; 5 of 5 read names carry no flow-cell coordinates"

    n_dup="$(run_samtools view -c -f 1024 "${dir_out}/x.bam")"

    if [[ "${n_dup}" -eq 8 ]]; then
        record_pass "${lbl}: the same 8 duplicate reads as with coordinates"
    else
        record_fail "${lbl}: ${n_dup} duplicate reads, expected 8"
    fi

    if \
        gzip -t "${dir_out}/x.markdup.txt.gz" 2> /dev/null
    then
        record_pass "${lbl}: statistics still written"
    else
        record_fail "${lbl}: no statistics file"
    fi
fi


# Older and longer Illumina names parse as Samtools parses them: four colons
# ('machine:lane:tile:x:y') and seven with a UMI after y.
fq_1="${dir_fq_se}/opt_hiseq_se.atria.fastq.gz"
lbl="SE 4-colon names"

if \
    run_case \
        "${lbl}" \
        hiseq \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --fq_1 "${fq_1}" \
        --fil_out "${tmp}/hiseq/x.bam"
then
    check_dt "${lbl}" "${dir_out}/x.bam" "HWI-ST100:8:1101:20000:20000	-
HWI-ST100:8:1101:20030:20000	SQ"

    check_dist \
        "${lbl}" \
        "${dir_out}/x.markdup.txt.gz" \
        "2500"
fi

fq_1="${dir_fq_se}/opt_umi_se.atria.fastq.gz"
lbl="SE 7-colon names (UMI)"

if \
    run_case \
        "${lbl}" \
        umi \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --fq_1 "${fq_1}" \
        --fil_out "${tmp}/umi/x.bam"
then
    check_dt "${lbl}" "${dir_out}/x.bam" "${pfx}:20000:20000:ACGTACGT	-
${pfx}:20030:20000:ACGTACGT	SQ"

    check_dist \
        "${lbl}" \
        "${dir_out}/x.markdup.txt.gz" \
        "2500"
fi


# Coordinate and SRA-style names mixed, in both orders: no '-d'.
for ord in first last; do
    lbl="PE mixed names, coordinates ${ord}"
    fq_1="${dir_fq_pe}/opt_mix_${ord}_pe_R1.atria.fastq.gz"
    fq_2="${dir_fq_pe}/opt_mix_${ord}_pe_R2.atria.fastq.gz"

    if \
        run_case \
            "${lbl}" \
            "mix_${ord}" \
            --aligner bowtie2 \
            --index "${idx_bt2}" \
            --fq_1 "${fq_1}" \
            --fq_2 "${fq_2}" \
            --fil_out "${tmp}/mix_${ord}/x.bam"
    then
        check_no_optical \
            "${lbl}" \
            "${dir_out}/x.bam" \
            "1 of 3 read names carry no flow-cell coordinates"

        check_dist \
            "${lbl}" \
            "${dir_out}/x.markdup.txt.gz" \
            "not used; 1 of 3 read names carry no flow-cell coordinates"
    fi
done


# Two reads whose names carry coordinates but whose tile prefixes differ in
# length. In the 'tile' set, the second read has a 5-digit tile. In the 'run'
# set, it has another run number. Before Samtools 1.25, 'markdup' compared such
# names wrongly (Samtools Issue #2399), so Step 4 dropped '-d' for them. Now
# '-d' applies. The two reads are not on the same tile of the same run, so the
# duplicate is tagged 'LB' and not 'SQ'. Samtools picks which read is the
# original, so the test checks the two tags ('-' and 'LB') and not which read
# carries each.
for nam in tile run; do
    lbl="PE mixed prefix lengths (${nam})"
    fq_1="${dir_fq_pe}/opt_len_${nam}_pe_R1.atria.fastq.gz"
    fq_2="${dir_fq_pe}/opt_len_${nam}_pe_R2.atria.fastq.gz"

    if \
        run_case \
            "${lbl}" \
            "len_${nam}" \
            --aligner bowtie2 \
            --index "${idx_bt2}" \
            --fq_1 "${fq_1}" \
            --fq_2 "${fq_2}" \
            --fil_out "${tmp}/len_${nam}/x.bam"
    then
        tag="$(print_dt "${dir_out}/x.bam" | cut -f 2 | LC_ALL=C sort)"

        if \
            print_pg_markdup "${dir_out}/x.bam" \
                | grep -q -- ' -d 2500' \
                && [[ "${tag}" == "$(printf '%s\n' - LB)" ]]
        then
            record_pass "${lbl}: '-d 2500', one library duplicate"
        else
            record_fail \
                "${lbl}: markdup '@PG' is" \
                "'$(print_pg_markdup "${dir_out}/x.bam")' and the tags are" \
                "'${tag//$'\n'/ }', expected '-d 2500' and '- LB'"
        fi

        check_dist \
            "${lbl}" \
            "${dir_out}/x.markdup.txt.gz" \
            "2500"
    fi
done


# CRAM output: the statistics take the CRAM file's stem.
lbl="PE CRAM output"
fq_1="${dir_fq_pe}/opt_pe_pe_R1.atria.fastq.gz"
fq_2="${dir_fq_pe}/opt_pe_pe_R2.atria.fastq.gz"

if \
    run_case \
        "${lbl}" \
        cram \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --fq_1 "${fq_1}" \
        --fq_2 "${fq_2}" \
        --ref_fa "${ref_rpt}" \
        --fil_out "${tmp}/cram/x.cram"
then
    got="$(
        cd "${dir_out}" && find . -type f | LC_ALL=C sort | paste -sd ' ' -
    )"
    exp="./x.cram ./x.cram.crai ./x.markdup.txt.gz"

    if [[ "${got}" == "${exp}" ]]; then
        record_pass "${lbl}: writes the CRAM, its index, and the statistics"
    else
        record_fail "${lbl}: wrote '${got}', expected '${exp}'"
    fi

    check_dt "${lbl}" "${dir_out}/x.cram" "${exp_pe}"

    check_dist \
        "${lbl}" \
        "${dir_out}/x.markdup.txt.gz" \
        "2500"
fi


# A path with a space.
lbl="PE path with a space"

if \
    run_case \
        "${lbl}" \
        "spaced dir" \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --fq_1 "${fq_1}" \
        --fq_2 "${fq_2}" \
        --fil_out "${tmp}/spaced dir/x y.bam"
then
    if \
        gzip -t "${dir_out}/x y.markdup.txt.gz" 2> /dev/null
    then
        record_pass "${lbl}: statistics written beside the output"
    else
        record_fail "${lbl}: no '${dir_out}/x y.markdup.txt.gz'"
    fi

    check_dt "${lbl}" "${dir_out}/x y.bam" "${exp_pe}"

    check_dist \
        "${lbl}" \
        "${dir_out}/x y.markdup.txt.gz" \
        "2500"
fi


# No aligned reads: no '-d' and no warning, statistics still written.
lbl="SE no aligned reads"
fq_1="${dir_fq_se}/opt_none_se.atria.fastq.gz"

if \
    run_case \
        "${lbl}" \
        none \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --fq_1 "${fq_1}" \
        --fil_out "${tmp}/none/x.bam"
then
    n_rec="$(run_samtools view -c "${dir_out}/x.bam")"

    if \
        [[ "${n_rec}" -eq 0 ]] \
        && ! grep -q 'optical duplicates are not told apart' <<< "${out}"
    then
        record_pass "${lbl}: empty output, no optical warning"
    else
        record_fail "${lbl}: ${n_rec} records; output: ${out}"
    fi

    if \
        gzip -t "${dir_out}/x.markdup.txt.gz" 2> /dev/null
    then
        record_pass "${lbl}: statistics still written"
    else
        record_fail "${lbl}: no statistics file"
    fi

    check_dist \
        "${lbl}" \
        "${dir_out}/x.markdup.txt.gz" \
        "not used; no reads"
fi


# A work file already marked: Step 4 is skipped and writes no statistics.
lbl="Step 4 on a marked file"
dir_out="${tmp}/marked"

mkdir -p "${dir_out}"
cp "${tmp}/pe/x.bam" "${dir_out}/x.bam"

rc=0

# shellcheck disable=SC2016  # Expand in the child shell, not this one.
out="$(
    "${TEST_BASH}" -c '
        # shellcheck source=lib/bash/core/source_helpers.sh
        source "${1}/lib/bash/core/source_helpers.sh"

        source_helpers "${1}/lib/bash" \
            core/check_args \
            core/check_inputs \
            core/format_outputs \
            workflows/align_fastqs
        _mark_dup_bam \
            align_fastqs \
            false \
            1 \
            "${2}/x.bam" \
            "${2}/x.markdup.txt"
    ' _ "${ROOT_REPO}" "${dir_out}" 2>&1
)" || rc=$?

if \
    [[ "${rc}" -eq 0 ]] \
    && grep -q 'Skipping markdup operations' <<< "${out}" \
    && [[ ! -e "${dir_out}/x.markdup.txt.gz" ]] \
    && [[ ! -e "${dir_out}/x.markdup.txt" ]]
then
    record_pass "${lbl}: skipped, no statistics written"
else
    record_fail "${lbl}: exit ${rc}; output: ${out}"
fi


# The names check on its own: a name with the colons of a flow-cell name but a
# letter where x goes carries no coordinates Samtools can read.
lbl="names with a letter where x goes"
rc=0

# shellcheck disable=SC2016  # Expand in the child shell, not this one.
out="$(
    "${TEST_BASH}" -c '
        # shellcheck source=lib/bash/core/source_helpers.sh
        source "${1}/lib/bash/core/source_helpers.sh"

        source_helpers "${1}/lib/bash" workflows/align_fastqs
        rc_chk=0
        _check_optical_names align_fastqs 1 "${2}" dist 2> /dev/null \
            || rc_chk=$?
        printf "%s|%s\n" "${rc_chk}" "${dist}"
    ' _ "${ROOT_REPO}" "${dir_fx}/sam/optical_names.sam" 2>&1
)" || rc=$?
exp="1|not used; 1 of 2 read names carry no flow-cell coordinates"

if [[ "${rc}" -eq 0 && "${out}" == "${exp}" ]]; then
    record_pass "${lbl}: no '-d', and the reason"
else
    record_fail "${lbl}: exit ${rc}; got '${out}', expected '${exp}'"
fi


# BWA-MEM keeps the coordinates in the names too.
lbl="PE BWA-MEM"
fq_1="${dir_fq_pe}/opt_pe_pe_R1.atria.fastq.gz"
fq_2="${dir_fq_pe}/opt_pe_pe_R2.atria.fastq.gz"

if \
    run_case \
        "${lbl}" \
        bwa \
        --aligner bwa \
        --index "${idx_bwa}" \
        --fq_1 "${fq_1}" \
        --fq_2 "${fq_2}" \
        --fil_out "${tmp}/bwa/x.bam"
then
    check_dt "${lbl}" "${dir_out}/x.bam" "${exp_pe}"

    check_dist \
        "${lbl}" \
        "${dir_out}/x.markdup.txt.gz" \
        "2500"
fi


finish
