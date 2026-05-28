#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: tests/trim_fastqs/scripts/make_fixtures.sh
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


#  Resolve paths relative to 'tests/trim_fastqs/scripts'
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_trm="$(cd "${dir_scr}/.." > /dev/null 2>&1 && pwd)"
dir_fix="${dir_trm}/fixtures"
dir_fq="${dir_fix}/fastq"

tmp_fq_se="${dir_fq}/tiny_se.fastq.tmp"
tmp_fq_p1="${dir_fq}/tiny_pe_R1.fastq.tmp"
tmp_fq_p2="${dir_fq}/tiny_pe_R2.fastq.tmp"
gz_fq_se="${dir_fq}/tiny_se.fastq.gz"
gz_fq_p1="${dir_fq}/tiny_pe_R1.fastq.gz"
gz_fq_p2="${dir_fq}/tiny_pe_R2.fastq.gz"
fil_rdm="${dir_fix}/README.md"


#  Remove a generated fixture file only if it is inside the fixture directory
function rm_fixture_file() {
    local file="${1:-}"

    if [[ -z "${file}" ]]; then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "refusing to remove an empty file path." >&2
        exit 1
    elif [[ "${file}" != "${dir_fix}/"* ]]; then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "refusing to remove path outside fixture directory: '${file}'." >&2
        exit 1
    elif [[ -d "${file}" ]]; then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "refusing to remove directory: '${file}'." >&2
        exit 1
    fi

    rm -f -- "${file}"
}


#  Remove temporary FASTQ intermediates on normal exit or failure
function cleanup_tmp_fastqs() {
    rm_fixture_file "${tmp_fq_se}"
    rm_fixture_file "${tmp_fq_p1}"
    rm_fixture_file "${tmp_fq_p2}"
}


trap cleanup_tmp_fastqs EXIT


#  Require gzip for deterministic compressed workflow inputs
if ! \
    command -v gzip > /dev/null 2>&1
then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "'gzip' must be available to generate trim-fastqs fixtures." >&2
    exit 1
fi


#  Create fixture output directories
mkdir -p "${dir_fq}"

#  Remove stale temporary intermediates
cleanup_tmp_fastqs


#  Write tiny single-end FASTQ provenance and compressed input fixture
cat > "${tmp_fq_se}" << EOM
@tiny_trim_se_read_1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

gzip -n -c "${tmp_fq_se}" > "${gz_fq_se}"
rm_fixture_file "${tmp_fq_se}"


#  Write tiny paired-end FASTQ provenance and compressed input fixtures
cat > "${tmp_fq_p1}" << EOM
@tiny_trim_pe_pair_1/1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

cat > "${tmp_fq_p2}" << EOM
@tiny_trim_pe_pair_1/2
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

gzip -n -c "${tmp_fq_p1}" > "${gz_fq_p1}"
gzip -n -c "${tmp_fq_p2}" > "${gz_fq_p2}"

rm_fixture_file "${tmp_fq_p1}"
rm_fixture_file "${tmp_fq_p2}"


#  Write fixture documentation from the same source as the generated files
cat > "${fil_rdm}" << EOM
# Trim-fastqs test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic tests of the trim-fastqs workflow.

They are intentionally small, hand-checkable, and version-controlled directly in Git. Running \`scripts/make_fixtures.sh\` regenerates the fixture set deterministically.

The fixture set currently focuses on single-end and paired-end Atria trimming smoke tests. The FASTQ records are encoded directly in \`scripts/make_fixtures.sh\` so that the expected read names and sequences are easy to inspect.

<br />

## Files
Readable provenance:
- \`scripts/make_fixtures.sh\`

Regenerate fixtures from the repository root with:

\`\`\`bash
bash tests/trim_fastqs/scripts/make_fixtures.sh
\`\`\`

The FASTQ records are generated deterministically by \`scripts/make_fixtures.sh\` and committed only as compressed workflow inputs. The compressed FASTQ files can be inspected with \`gzip -cd\`.

Compressed FASTQ input:
- \`fastq/tiny_se.fastq.gz\`
- \`fastq/tiny_pe_R1.fastq.gz\`
- \`fastq/tiny_pe_R2.fastq.gz\`

The compressed FASTQ files are generated with \`gzip -n -c\` so gzip headers do not encode timestamps or source filenames. The uncompressed FASTQ files are not committed.

<br />

## Expected trimming behavior
The single-end FASTQ fixture contains one clean 64 bp read:

| Read name               | Sequence                                                             |
|:---                     |:---                                                                  |
| \`tiny_trim_se_read_1\` | \`ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\` |

\`submit_trim_fastqs.sh\` currently calls Atria with \`--length-range 35:500\`, so this read is intentionally longer than the minimum length filter. The submit-level smoke test should verify that Atria-backed trimming produces a non-empty gzip FASTQ output and preserves the expected read name.

The paired-end FASTQ fixtures contain one clean 64 bp read pair:

| Read name                 | Mate | Sequence                                                             |
|:---                       |:---  |:---                                                                  |
| \`tiny_trim_pe_pair_1/1\` | R1   | \`ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\` |
| \`tiny_trim_pe_pair_1/2\` | R2   | \`TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\` |

The paired-end submit-level smoke test should pass the two compressed inputs as one comma-delimited \`--csv_infile\` entry, verify that Atria-backed trimming produces non-empty gzip FASTQ outputs for R1 and R2, and confirm that each output contains the expected paired read name prefix.

<br />

## Deferred fixture batches
Later trim-fastqs batches should add execute-level SE/PE coverage, local GNU Parallel coverage, and optional Slurm coverage.
EOM


echo "Wrote trim-fastqs fixtures under: ${dir_fix}"
