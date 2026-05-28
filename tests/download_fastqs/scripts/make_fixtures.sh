#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: tests/download_fastqs/scripts/make_fixtures.sh
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


#  Resolve paths relative to 'tests/download_fastqs/scripts'
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_dwn="$(cd "${dir_scr}/.." > /dev/null 2>&1 && pwd)"
dir_fix="${dir_dwn}/fixtures"
dir_src="${dir_fix}/source"
dir_mtd="${dir_fix}/metadata"

tmp_fq_se="${dir_src}/tiny_download_se.fastq.tmp"
tmp_fq_p1="${dir_src}/tiny_download_pe_R1.fastq.tmp"
tmp_fq_p2="${dir_src}/tiny_download_pe_R2.fastq.tmp"
gz_fq_se="${dir_src}/tiny_download_se.fastq.gz"
gz_fq_p1="${dir_src}/tiny_download_pe_R1.fastq.gz"
gz_fq_p2="${dir_src}/tiny_download_pe_R2.fastq.gz"
tpl_mtd_se="${dir_mtd}/local_se.template.tsv"
tpl_mtd_pe="${dir_mtd}/local_pe.template.tsv"
tpl_mtd_mix="${dir_mtd}/local_mixed.template.tsv"
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
        "'gzip' must be available to generate download-fastqs fixtures." >&2
    exit 1
fi


#  Create fixture output directories
mkdir -p "${dir_src}" "${dir_mtd}"

#  Remove stale temporary intermediates
cleanup_tmp_fastqs


#  Write tiny single-end FASTQ provenance and compressed source fixture
cat > "${tmp_fq_se}" << EOM
@tiny_download_se_read_1
GATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

gzip -n -c "${tmp_fq_se}" > "${gz_fq_se}"
rm_fixture_file "${tmp_fq_se}"


#  Write tiny paired-end FASTQ provenance and compressed source fixtures
cat > "${tmp_fq_p1}" << EOM
@tiny_download_pe_pair_1/1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

cat > "${tmp_fq_p2}" << EOM
@tiny_download_pe_pair_1/2
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

gzip -n -c "${tmp_fq_p1}" > "${gz_fq_p1}"
gzip -n -c "${tmp_fq_p2}" > "${gz_fq_p2}"

rm_fixture_file "${tmp_fq_p1}"
rm_fixture_file "${tmp_fq_p2}"


#  Write metadata template. The smoke test replaces '__BASE_URL__' at runtime
#+ with a loopback HTTP server URL.
cat > "${tpl_mtd_se}" << EOM
run_accession	custom_name	fastq_https
SRR_LOCAL_SE	tiny_download_se	__BASE_URL__/tiny_download_se.fastq.gz
EOM

cat > "${tpl_mtd_pe}" << EOM
run_accession	custom_name	fastq_https
SRR_LOCAL_PE	tiny_download_pe	__BASE_URL__/tiny_download_pe_R1.fastq.gz;__BASE_URL__/tiny_download_pe_R2.fastq.gz
EOM

cat > "${tpl_mtd_mix}" << EOM
run_accession	custom_name	fastq_https
SRR_LOCAL_SE	tiny_download_se	__BASE_URL__/tiny_download_se.fastq.gz
SRR_LOCAL_PE	tiny_download_pe	__BASE_URL__/tiny_download_pe_R1.fastq.gz;__BASE_URL__/tiny_download_pe_R2.fastq.gz
EOM


#  Write fixture documentation from the same source as the generated files
cat > "${fil_rdm}" << EOM
# Download-fastqs test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic no-network tests of the download-fastqs workflow.

They are intentionally small, hand-checkable, and version-controlled directly in Git. Running \`scripts/make_fixtures.sh\` regenerates the fixture set deterministically.

The fixture set currently focuses on single-end and paired-end local-download smoke tests. The FASTQ records and metadata templates are encoded directly in \`scripts/make_fixtures.sh\` so that the expected read names, sequences, and table shape are easy to inspect.

<br />

## Files
Readable provenance:
- \`scripts/make_fixtures.sh\`

Regenerate fixtures from the repository root with:

\`\`\`bash
bash tests/download_fastqs/scripts/make_fixtures.sh
\`\`\`

Compressed source FASTQ:
- \`source/tiny_download_se.fastq.gz\`
- \`source/tiny_download_pe_R1.fastq.gz\`
- \`source/tiny_download_pe_R2.fastq.gz\`

Metadata template:
- \`metadata/local_se.template.tsv\`
- \`metadata/local_pe.template.tsv\`
- \`metadata/local_mixed.template.tsv\`

The compressed FASTQ files are generated with \`gzip -n -c\` so gzip headers do not encode timestamps or source filenames. The uncompressed FASTQ files are not committed.

The metadata template uses the columns required by \`execute_download_fastqs.sh\`:

\`\`\`text
run_accession	custom_name	fastq_https
SRR_LOCAL_SE	tiny_download_se	__BASE_URL__/tiny_download_se.fastq.gz
SRR_LOCAL_PE	tiny_download_pe	__BASE_URL__/tiny_download_pe_R1.fastq.gz;__BASE_URL__/tiny_download_pe_R2.fastq.gz
\`\`\`

The smoke test replaces \`__BASE_URL__\` at runtime with a loopback HTTP URL served from \`127.0.0.1\`. This keeps the wrapper-backed download test no-network while still exercising real \`wget\` download behavior.

\`metadata/local_mixed.template.tsv\` contains both rows in one table to exercise the common case where users provide one mixed SE/PE metadata TSV.

<br />

## Expected download behavior
The single-end FASTQ source contains one clean 64 bp read:

| Read name                   | Sequence                                                             |
| :---                        | :---                                                                 |
| \`tiny_download_se_read_1\` | \`GATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAG\` |

\`execute_download_fastqs.sh\` should download the source URL to \`SRR_LOCAL_SE.fastq.gz\`, create the custom symlink \`tiny_download_se.fastq.gz\`, preserve gzip integrity, and preserve the expected read name.

The paired-end FASTQ sources contain one clean 64 bp read pair:

| Read name                    | Mate | Sequence                                                         |
| :---                         | :--- | :---                                                             |
| \`tiny_download_pe_pair_1/1\` | R1   | \`ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\` |
| \`tiny_download_pe_pair_1/2\` | R2   | \`TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\` |

\`execute_download_fastqs.sh\` should parse the semicolon-delimited PE URL field, download the sources to \`SRR_LOCAL_PE_R1.fastq.gz\` and \`SRR_LOCAL_PE_R2.fastq.gz\`, create custom symlinks \`tiny_download_pe_R1.fastq.gz\` and \`tiny_download_pe_R2.fastq.gz\`, preserve gzip integrity, and preserve the expected mate-specific read names.

<br />

## Deferred fixture batches
Later download-fastqs batches should add execute-level multiple-record metadata coverage, local GNU Parallel coverage if applicable, optional Slurm coverage, and optional external-network coverage gated by \`RUN_NETWORK=1\`.
EOM


echo "Wrote download-fastqs fixtures under: ${dir_fix}"
