#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: tests/compute_signal/scripts/make_fixtures.sh
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


#  Resolve paths relative to 'tests/compute_signal/scripts'
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_sig="$(cd "${dir_scr}/.." > /dev/null 2>&1 && pwd)"
dir_fix="${dir_sig}/fixtures"

dir_bdg="${dir_fix}/bedgraph"
dir_ref="${dir_fix}/reference"
dir_sam="${dir_fix}/sam"
dir_bam="${dir_fix}/bam"
dir_cram="${dir_fix}/cram"

fil_bg_A="${dir_bdg}/ratio_A.bdg"
fil_bg_B="${dir_bdg}/ratio_B.bdg"
fil_bg_hdr_A="${dir_bdg}/ratio_headers_A.bdg"
fil_bg_hdr_B="${dir_bdg}/ratio_headers_B.bdg"

fil_ref="${dir_ref}/tiny.fa"
fil_sam_se="${dir_sam}/tiny_se.sam"
fil_sam_pe="${dir_sam}/tiny_pe.sam"

fil_bam_se="${dir_bam}/tiny_se.bam"
fil_bam_pe="${dir_bam}/tiny_pe.bam"

fil_cram_se="${dir_cram}/tiny_se.cram"
fil_cram_pe="${dir_cram}/tiny_pe.cram"

fil_read="${dir_fix}/README.md"

env_req="env_protocol"


#  Check that the project environment is active before writing fixtures
if [[ "${CONDA_DEFAULT_ENV:-}" != "${env_req}" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "activate '${env_req}' before generating compute-signal fixtures;" \
        "current environment: '${CONDA_DEFAULT_ENV:-none}'." >&2
    exit 1
fi


#  Create fixture output directories
mkdir -p \
    "${dir_bdg}" \
    "${dir_ref}" \
    "${dir_sam}" \
    "${dir_bam}" \
    "${dir_cram}"


#  Write tiny bedGraph fixtures for ratio-mode tests
{
    printf "I\t0\t10\t4\n"
    printf "I\t10\t20\t0\n"
    printf "I\t20\t30\t5\n"
    printf "I\t30\t40\t0\n"
    printf "I\t40\t50\t2\n"
    printf "I\t50\t60\t1\n"
    printf "I\t60\t70\t1\n"
    printf "I\t70\t80\t1\n"
} > "${fil_bg_A}"

{
    printf "I\t0\t10\t2\n"
    printf "I\t10\t20\t2\n"
    printf "I\t20\t30\t0\n"
    printf "I\t30\t40\t0\n"
    printf "I\t40\t50\t0.5\n"
    printf "I\t50\t60\t0.04\n"
    printf "I\t60\t70\t3\n"
    printf "I\t70\t80\t1\n"
} > "${fil_bg_B}"

{
    printf "track type=bedGraph name=\"ratio_A\"\n"
    printf "browser position I:0-80\n"
    printf "# comment in A\n"
    printf "customHeader sample=A\n"
    printf "I\t0\t10\t4\n"
    printf "I\t10\t20\t0\n"
    printf "  customHeader indented=A\n"
    printf "I\t20\t30\t5\n"
    printf "I\t30\t40\t0\n"
    printf "I\t40\t50\t2\n"
    printf "I\t50\t60\t1\n"
    printf "I\t60\t70\t1\n"
    printf "I\t70\t80\t1\n"
} > "${fil_bg_hdr_A}"

{
    printf "track type=bedGraph name=\"ratio_B\"\n"
    printf "browser position I:0-80\n"
    printf "# comment in B\n"
    printf "customHeader sample=B\n"
    printf "I\t0\t10\t2\n"
    printf "I\t10\t20\t2\n"
    printf "  customHeader indented=B\n"
    printf "I\t20\t30\t0\n"
    printf "I\t30\t40\t0\n"
    printf "I\t40\t50\t0.5\n"
    printf "I\t50\t60\t0.04\n"
    printf "I\t60\t70\t3\n"
    printf "I\t70\t80\t1\n"
} > "${fil_bg_hdr_B}"


#  Write tiny reference FASTA used for BAM/CRAM fixture generation
cat > "${fil_ref}" << 'EOM'
>I
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOM

#  Write tiny single-end SAM provenance fixture
{
    printf "@HD\tVN:1.6\tSO:coordinate\n"
    printf "@SQ\tSN:I\tLN:80\n"
    printf "se_fwd\t0\tI\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tFFFFFFFFFF\n"
    printf "se_rev\t16\tI\t21\t60\t10M\t*\t0\t0\tACGTACGTAC\tFFFFFFFFFF\n"
} > "${fil_sam_se}"

#  Write tiny paired-end SAM provenance fixture
{
    printf "@HD\tVN:1.6\tSO:coordinate\n"
    printf "@SQ\tSN:I\tLN:80\n"
    printf "pe_1\t99\tI\t11\t60\t10M\t=\t31\t30\tACGTACGTAC\tFFFFFFFFFF\n"
    printf "pe_1\t147\tI\t31\t60\t10M\t=\t11\t-30\tACGTACGTAC\tFFFFFFFFFF\n"
    printf "pe_2\t99\tI\t41\t60\t10M\t=\t51\t20\tACGTACGTAC\tFFFFFFFFFF\n"
    printf "pe_2\t147\tI\t51\t60\t10M\t=\t41\t-20\tACGTACGTAC\tFFFFFFFFFF\n"
} > "${fil_sam_pe}"


#  Generate FASTA index, BAM, BAM indexes, CRAM, and CRAM indexes
if ! \
    command -v samtools > /dev/null 2>&1
then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "'samtools' must be available in '${env_req}' to generate BAM/CRAM" \
        "fixtures." >&2
    exit 1
fi

samtools faidx "${fil_ref}"

samtools view -bS "${fil_sam_se}" | samtools sort -o "${fil_bam_se}"
samtools index "${fil_bam_se}"

samtools view -bS "${fil_sam_pe}" | samtools sort -o "${fil_bam_pe}"
samtools index "${fil_bam_pe}"

samtools view -C -T "${fil_ref}" -o "${fil_cram_se}" "${fil_bam_se}"
samtools index "${fil_cram_se}"

samtools view -C -T "${fil_ref}" -o "${fil_cram_pe}" "${fil_bam_pe}"
samtools index "${fil_cram_pe}"

samtools quickcheck \
    "${fil_bam_se}" \
    "${fil_bam_pe}" \
    "${fil_cram_se}" \
    "${fil_cram_pe}"


#  Write fixture documentation
cat > "${fil_read}" << 'EOM'
# Compute-signal test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic tests of the compute-signal workflow.

They are intentionally small, hand-checkable, and version-controlled directly in Git. Running `scripts/make_fixtures.sh` from `env_protocol` regenerates the fixture set deterministically.

The first fixture batch focuses on ratio-mode tests using plain text bedGraph files. These fixtures are not derived from real sequencing data. They were written by hand so that the expected behavior of each row is easy to inspect and reason about.

<br />

## Ratio bedGraph fixtures
Files:
- `bedgraph/ratio_A.bdg`
- `bedgraph/ratio_B.bdg`
- `bedgraph/ratio_headers_A.bdg`
- `bedgraph/ratio_headers_B.bdg`

The `ratio_A.bdg` and `ratio_B.bdg` files are plain 4-column bedGraph files with matching bins and no header lines.

The `ratio_headers_A.bdg` and `ratio_headers_B.bdg` files preserve the same data rows, but include header-like lines for `--skp_pfx` coverage. They include default skipped prefixes (`track`, `browser`, and `#`) plus the custom prefix `customHeader`.

<br />

## Row-level expected behavior
| Bin       | A  | B    | Purpose |
|:---       |:---|:---  |:---     |
| `I:0-10`  | 4  | 2    | Simple unadjusted ratio: `A / B = 2`.                                               |
| `I:10-20` | 0  | 2    | Zero numerator: unadjusted ratio is `0`; with pseudocounts `1:1`, ratio is `1 / 3`. |
| `I:20-30` | 5  | 0    | Zero denominator: test non-finite, stabilized, or filtered behavior.                |
| `I:30-40` | 0  | 0    | Zero-zero bin for `skip_00` behavior before or after scaling.                       |
| `I:40-50` | 2  | 0.5  | Scaling-sensitive finite bin: baseline ratio is `4`.                                |
| `I:50-60` | 1  | 0.04 | Denominator-floor case: with `dep_min = 0.1`, ratio is `10`.                        |
| `I:60-70` | 1  | 3    | Decimal-rounding case: with `dp = 3`, ratio is `0.333`.                             |
| `I:70-80` | 1  | 1    | Stable finite row: ratio is `1`.                                                    |

<br />

## Alignment fixtures
Signal-mode and coord-mode tests will use tiny SAM/BAM/CRAM fixtures generated from synthetic SAM and FASTA provenance files.

Files:
- `reference/tiny.fa`
- `reference/tiny.fa.fai`
- `sam/tiny_se.sam`
- `sam/tiny_pe.sam`
- `bam/tiny_se.bam`
- `bam/tiny_se.bam.bai`
- `bam/tiny_pe.bam`
- `bam/tiny_pe.bam.bai`
- `cram/tiny_se.cram`
- `cram/tiny_se.cram.crai`
- `cram/tiny_pe.cram`
- `cram/tiny_pe.cram.crai`

The SAM and FASTA files are committed as readable provenance. The FASTA index, BAM files, BAM indexes, CRAM files, and CRAM indexes are generated with `samtools` by `scripts/make_fixtures.sh`.

These binary fixtures should be small enough to commit directly unless they unexpectedly grow large enough to justify Git LFS.

CRAM support should not be treated as a distant optional feature. If the current compute-signal path cannot pass a reference FASTA cleanly to every layer that needs it, that gap should be exposed by the tiny CRAM fixture tests and addressed as part of the compute-signal testing/refactor work.

<br />

## First intended test
The first wet wrapper-to-Python ratio test should call `submit_compute_signal.sh` in local serial ratio mode and verify that a bedGraph output is created with expected finite rows.

Golden output files are intentionally deferred until the wrapper command shape and decimal formatting are stable.

EOM

echo "success($(basename "${BASH_SOURCE[0]}")):" \
    "generated compute-signal fixtures under '${dir_fix}'."
