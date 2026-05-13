#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: tests/align_fastqs/scripts/make_fixtures.sh
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


#  Resolve paths relative to 'tests/align_fastqs/scripts'
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_aln="$(cd "${dir_scr}/.." > /dev/null 2>&1 && pwd)"
dir_fix="${dir_aln}/fixtures"

dir_ref="${dir_fix}/reference"
dir_fq="${dir_fix}/fastq"
dir_bt2="${dir_fix}/bowtie2"
dir_bwa="${dir_fix}/bwa"
dir_bm2="${dir_fix}/bwa-mem2"

ref="${dir_ref}/tiny.fa"
ref_bwa="${dir_bwa}/tiny.fa"
ref_bm2="${dir_bm2}/tiny.fa"
tmp_fq_se="${dir_fq}/tiny_se.atria.fastq.tmp"
tmp_fq_pe_1="${dir_fq}/tiny_pe_R1.atria.fastq.tmp"
tmp_fq_pe_2="${dir_fq}/tiny_pe_R2.atria.fastq.tmp"

gz_fq_se="${dir_fq}/tiny_se.atria.fastq.gz"
gz_fq_pe_1="${dir_fq}/tiny_pe_R1.atria.fastq.gz"
gz_fq_pe_2="${dir_fq}/tiny_pe_R2.atria.fastq.gz"

idx_bt2="${dir_bt2}/tiny"
idx_bwa="${ref_bwa}"
idx_bm2="${ref_bm2}"

fil_read="${dir_fix}/README.md"

env_req="env_protocol"


#  Remove a generated fixture file only if it is inside the fixture directory
function rm_fixture_file() {
    local file="${1:-}"

    if [[ -z "${file}" ]]; then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "refusing to remove an empty file path." >&2
        exit 1
    elif [[ "${file}" != "${dir_fix}/"* ]]; then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "refusing to remove path outside fixture directory:" \
            "'${file}'." >&2
        exit 1
    elif [[ -d "${file}" ]]; then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "refusing to remove directory:" \
            "'${file}'." >&2
        exit 1
    fi

    rm -f -- "${file}"
}


#  Remove temporary FASTQ intermediates on normal exit or failure
function cleanup_tmp_fastqs() {
    rm_fixture_file "${tmp_fq_se}"
    rm_fixture_file "${tmp_fq_pe_1}"
    rm_fixture_file "${tmp_fq_pe_2}"
}


trap cleanup_tmp_fastqs EXIT


#  Check that the project environment is active before writing fixtures
if [[ "${CONDA_DEFAULT_ENV:-}" != "${env_req}" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "activate '${env_req}' before generating align-fastqs fixtures;" \
        "current environment: '${CONDA_DEFAULT_ENV:-none}'." >&2
    exit 1
fi


#  Require alignment tools used to generate and validate fixtures
for cmd in bowtie2 bowtie2-build bowtie2-inspect bwa bwa-mem2 gzip samtools; do
    if ! command -v "${cmd}" > /dev/null 2>&1; then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "'${cmd}' must be available in '${env_req}' to generate" \
            "align-fastqs fixtures." >&2
        exit 1
    fi
done
unset cmd


#  Create fixture output directories
mkdir -p \
    "${dir_ref}" \
    "${dir_fq}" \
    "${dir_bt2}" \
    "${dir_bwa}" \
    "${dir_bm2}"

#  Remove obsolete uncompressed FASTQ fixtures and stale temp intermediates
rm_fixture_file "${dir_fq}/tiny_se.atria.fastq"
rm_fixture_file "${dir_fq}/tiny_pe_R1.atria.fastq"
rm_fixture_file "${dir_fq}/tiny_pe_R2.atria.fastq"
cleanup_tmp_fastqs


#  Write tiny reference FASTA used for Bowtie2 alignment tests
cat > "${ref}" << 'EOM'
>I
GATCGTACCTAGGCTAACGTTGACCGTTAACGATCGTAGCTAGGATCCGTTACGATCGATGCTAGCTTACCGGATCAAGCTTAGGCTAATCGGCTAAGGTTCCGATTA
EOM

cp "${ref}" "${ref_bwa}"
cp "${ref}" "${ref_bm2}"


#  Write tiny single-end FASTQ provenance and compressed input fixture
cat > "${tmp_fq_se}" << 'EOM'
@tiny_se_read_1
ACGTTGACCGTTAACGATCGTAGCTAGGAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

gzip -n -c "${tmp_fq_se}" > "${gz_fq_se}"
rm_fixture_file "${tmp_fq_se}"


#  Write tiny paired-end FASTQ provenance and compressed input fixtures
cat > "${tmp_fq_pe_1}" << 'EOM'
@tiny_pe_pair_1
ACGTTGACCGTTAACGATCGTAGCTAGGAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

cat > "${tmp_fq_pe_2}" << 'EOM'
@tiny_pe_pair_1
CCTTAGCCGATTAGCCTAAGCTTGATCCGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

gzip -n -c "${tmp_fq_pe_1}" > "${gz_fq_pe_1}"
gzip -n -c "${tmp_fq_pe_2}" > "${gz_fq_pe_2}"

rm_fixture_file "${tmp_fq_pe_1}"
rm_fixture_file "${tmp_fq_pe_2}"


#  Generate and validate Bowtie2 index files
rm_fixture_file "${idx_bt2}.1.bt2"
rm_fixture_file "${idx_bt2}.2.bt2"
rm_fixture_file "${idx_bt2}.3.bt2"
rm_fixture_file "${idx_bt2}.4.bt2"
rm_fixture_file "${idx_bt2}.rev.1.bt2"
rm_fixture_file "${idx_bt2}.rev.2.bt2"
rm_fixture_file "${idx_bt2}.1.bt2l"
rm_fixture_file "${idx_bt2}.2.bt2l"
rm_fixture_file "${idx_bt2}.3.bt2l"
rm_fixture_file "${idx_bt2}.4.bt2l"
rm_fixture_file "${idx_bt2}.rev.1.bt2l"
rm_fixture_file "${idx_bt2}.rev.2.bt2l"

log_bt2="${dir_bt2}/bowtie2-build.log"

if ! \
    bowtie2-build "${ref}" "${idx_bt2}" > "${log_bt2}" 2>&1
then
    cat "${log_bt2}" >&2
    exit 1
fi

rm_fixture_file "${log_bt2}"

bowtie2-inspect -n "${idx_bt2}" > /dev/null

bowtie2 \
    -x "${idx_bt2}" \
    -U "${gz_fq_se}" \
    -S /dev/null \
        > /dev/null \
        2> /dev/null

bowtie2 \
    -x "${idx_bt2}" \
    --very-sensitive \
    --no-mixed \
    --no-discordant \
    --no-overlap \
    --no-dovetail \
    -1 "${gz_fq_pe_1}" \
    -2 "${gz_fq_pe_2}" \
    -S /dev/null \
        > /dev/null \
        2> /dev/null


#  Generate and validate BWA index files
rm_fixture_file "${idx_bwa}.amb"
rm_fixture_file "${idx_bwa}.ann"
rm_fixture_file "${idx_bwa}.bwt"
rm_fixture_file "${idx_bwa}.pac"
rm_fixture_file "${idx_bwa}.sa"

log_bwa="${dir_bwa}/bwa-index.log"

if ! \
    bwa index "${idx_bwa}" > "${log_bwa}" 2>&1
then
    cat "${log_bwa}" >&2
    exit 1
fi

rm_fixture_file "${log_bwa}"

bwa mem "${idx_bwa}" "${gz_fq_se}" \
        > /dev/null \
        2> /dev/null

bwa mem "${idx_bwa}" "${gz_fq_pe_1}" "${gz_fq_pe_2}" \
        > /dev/null \
        2> /dev/null


#  Generate and validate BWA-MEM2 index files
rm_fixture_file "${idx_bm2}.0123"
rm_fixture_file "${idx_bm2}.amb"
rm_fixture_file "${idx_bm2}.ann"
rm_fixture_file "${idx_bm2}.bwt.2bit.64"
rm_fixture_file "${idx_bm2}.pac"

log_bm2="${dir_bm2}/bwa-mem2-index.log"

if ! \
    bwa-mem2 index "${idx_bm2}" > "${log_bm2}" 2>&1
then
    cat "${log_bm2}" >&2
    exit 1
fi

rm_fixture_file "${log_bm2}"

bwa-mem2 mem "${idx_bm2}" "${gz_fq_se}" \
    > /dev/null 2> /dev/null

bwa-mem2 mem "${idx_bm2}" "${gz_fq_pe_1}" "${gz_fq_pe_2}" \
    > /dev/null 2> /dev/null


#  Write fixture documentation
cat > "${fil_read}" << 'EOM'
# Align-fastqs test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic tests of the align-fastqs workflow.

They are intentionally small, hand-checkable, and version-controlled directly in Git. Running `scripts/make_fixtures.sh` from `env_protocol` regenerates the fixture set deterministically.

The fixture set focuses on single-end and paired-end Bowtie2 alignment to BAM and CRAM output, plus BWA MEM and bwa-mem2 MEM alignment to BAM output. These fixtures are not derived from real sequencing data. The reference sequence and FASTQ records are encoded directly in `scripts/make_fixtures.sh` so that the expected read alignments are easy to inspect.

<br />

## Files
Readable provenance:
- `reference/tiny.fa`
- `scripts/make_fixtures.sh`

The FASTQ records are generated deterministically by `scripts/make_fixtures.sh` and committed only as compressed workflow inputs. The compressed FASTQ files can be inspected with `gzip -cd`.

Compressed FASTQ input:
- `fastq/tiny_se.atria.fastq.gz`
- `fastq/tiny_pe_R1.atria.fastq.gz`
- `fastq/tiny_pe_R2.atria.fastq.gz`

Bowtie2 index:
- `bowtie2/tiny.1.bt2`
- `bowtie2/tiny.2.bt2`
- `bowtie2/tiny.3.bt2`
- `bowtie2/tiny.4.bt2`
- `bowtie2/tiny.rev.1.bt2`
- `bowtie2/tiny.rev.2.bt2`

BWA MEM index:
- `bwa/tiny.fa`
- `bwa/tiny.fa.amb`
- `bwa/tiny.fa.ann`
- `bwa/tiny.fa.bwt`
- `bwa/tiny.fa.pac`
- `bwa/tiny.fa.sa`

bwa-mem2 MEM index:
- `bwa-mem2/tiny.fa`
- `bwa-mem2/tiny.fa.0123`
- `bwa-mem2/tiny.fa.amb`
- `bwa-mem2/tiny.fa.ann`
- `bwa-mem2/tiny.fa.bwt.2bit.64`
- `bwa-mem2/tiny.fa.pac`

The compressed FASTQ, Bowtie2 index, BWA index, and bwa-mem2 index files are generated by `scripts/make_fixtures.sh`. They are committed because align/trim/download workflows commonly operate on compressed FASTQ inputs and real aligner indexes. The uncompressed FASTQ files are not committed; use `gzip -cd` to inspect the compressed inputs.

<br />

## Expected alignment behavior
The single-end FASTQ fixture contains one read:

| Read name        | Sequence                         | Expected reference |
|:---              |:---                              |:---                |
| `tiny_se_read_1` | `ACGTTGACCGTTAACGATCGTAGCTAGGAT` | `I`                |

Bowtie2, BWA MEM, and bwa-mem2 MEM should align this read to chromosome `I` in the tiny reference. Smoke tests should verify that wrapper-generated BAM/CRAM files and indexes exist, pass `samtools quickcheck`, report or count one read alignment on `I`, and contain the expected read name. CRAM checks should read with `samtools view -T reference/tiny.fa`.

The paired-end FASTQ fixtures contain one read pair:

| Read name        | Mate | Sequence                         | Expected reference | Expected start |
|:---              |:---  |:---                              |:---                |:---            |
| `tiny_pe_pair_1` | R1   | `ACGTTGACCGTTAACGATCGTAGCTAGGAT` | `I`                | 17             |
| `tiny_pe_pair_1` | R2   | `CCTTAGCCGATTAGCCTAAGCTTGATCCGG` | `I`                | 70             |

R1 matches chromosome `I` in forward orientation. R2 is the reverse complement of positions 70-99 on chromosome `I`, so Bowtie2, BWA MEM, and bwa-mem2 MEM should report the pair in FR orientation. Smoke tests should verify that wrapper-generated BAM/CRAM files and indexes exist, pass `samtools quickcheck`, report or count two read alignments on `I`, contain the expected paired read name, and include proper-pair flag records when the aligner reports stable proper-pair flags. CRAM checks should read with `samtools view -T reference/tiny.fa`.

<br />

## Deferred fixture batches
Later align-fastqs batches should add BWA ALN indexes and local GNU Parallel coverage.
EOM

echo "success($(basename "${BASH_SOURCE[0]}")):" \
    "generated align-fastqs fixtures under '${dir_fix}'."
