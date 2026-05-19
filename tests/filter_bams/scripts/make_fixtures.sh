#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: make_fixtures.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


set -euo pipefail

dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_fx="$(cd "${dir_scr}/../fixtures" > /dev/null 2>&1 && pwd)"
dir_sam="${dir_fx}/sam"
dir_ref="${dir_fx}/reference"

mkdir -p "${dir_sam}" "${dir_ref}"

sam="${dir_sam}/filter_sc_sp.sam"
ref="${dir_ref}/filter_sc_sp.fa"
fai="${ref}.fai"
readme="${dir_fx}/README.md"

contigs=(
    I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI
    Mito
    SP_I SP_II SP_III SP_II_TG SP_MTR SP_Mito
    chrUn
)

seq_100="ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"

: > "${ref}"
: > "${fai}"
offset=0
for contig in "${contigs[@]}"; do
    header=">${contig}"
    printf '%s\n%s\n' "${header}" "${seq_100}" >> "${ref}"
    offset=$(( offset + ${#header} + 1 ))
    printf '%s\t100\t%d\t100\t101\n' "${contig}" "${offset}" >> "${fai}"
    offset=$(( offset + 100 + 1 ))
done

cat > "${sam}" << 'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:I	LN:100
@SQ	SN:II	LN:100
@SQ	SN:III	LN:100
@SQ	SN:IV	LN:100
@SQ	SN:V	LN:100
@SQ	SN:VI	LN:100
@SQ	SN:VII	LN:100
@SQ	SN:VIII	LN:100
@SQ	SN:IX	LN:100
@SQ	SN:X	LN:100
@SQ	SN:XI	LN:100
@SQ	SN:XII	LN:100
@SQ	SN:XIII	LN:100
@SQ	SN:XIV	LN:100
@SQ	SN:XV	LN:100
@SQ	SN:XVI	LN:100
@SQ	SN:Mito	LN:100
@SQ	SN:SP_I	LN:100
@SQ	SN:SP_II	LN:100
@SQ	SN:SP_III	LN:100
@SQ	SN:SP_II_TG	LN:100
@SQ	SN:SP_MTR	LN:100
@SQ	SN:SP_Mito	LN:100
@SQ	SN:chrUn	LN:100
r_sc_I	0	I	1	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF
r_sc_mito	0	Mito	1	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF
r_sp_i	0	SP_I	1	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF
r_sp_tg	0	SP_II_TG	1	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF
r_sp_mtr	0	SP_MTR	1	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF
r_sp_mito	0	SP_Mito	1	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF
r_other	0	chrUn	1	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF
SAM

cat > "${readme}" << 'README'
# filter_bams fixtures

These fixtures exercise `submit_filter_bams.sh` and `execute_filter_bams.sh`
with tiny deterministic alignment inputs.

Regenerate them from the repository root with:

```bash
bash tests/filter_bams/scripts/make_fixtures.sh
```

## Files

- `sam/filter_sc_sp.sam`
  - Coordinate-sorted SAM with one read each on `I`, `Mito`, `SP_I`,
    `SP_II_TG`, `SP_MTR`, `SP_Mito`, and `chrUn`.
  - Used by BAM smoke tests to build BAM input at test runtime.
  - Used by CRAM smoke tests to build CRAM input at test runtime.

- `reference/filter_sc_sp.fa`
  - Tiny reference FASTA with all contigs present in the SAM header.
  - Every contig is 100 bp.
  - Required for deterministic CRAM input generation, CRAM reading, and CRAM
    output writing.

- `reference/filter_sc_sp.fa.fai`
  - Deterministic FASTA index for `reference/filter_sc_sp.fa`.
  - Generated directly by this script because the reference uses fixed
    one-line contigs.

## Expected smoke-test behavior

- `retain=sc` keeps chromosome `I` and drops `Mito` unless `--mito` is used.
- `retain=sc` drops S. pombe contigs such as `SP_I`.
- `retain=sp --tg --mtr --mito` keeps `SP_I`, `SP_II_TG`, `SP_MTR`, and
  `SP_Mito`.
- `retain=sp` drops S. cerevisiae contigs such as `I` and unrelated contigs
  such as `chrUn`.
- CRAM input tests must pass `--ref_fa reference/filter_sc_sp.fa` and still
  produce BAM output.
- CRAM output tests must pass `--out_ext cram --ref_fa reference/filter_sc_sp.fa`
  and produce indexed CRAM output.
README

printf 'Wrote filter_bams fixtures under %s\n' "${dir_fx}"
