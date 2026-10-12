
# Download-fastqs test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic loopback-only tests of the download-fastqs workflow.

They are intentionally small and hand-checkable. Running `make.sh` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/download_fastqs/make.sh
```

Generated fixture outputs are ignored by Git. When required inputs are missing, `tests/run_tests.sh` regenerates this fixture set automatically.

The fixture set currently supports single-end, paired-end, and mixed-metadata local-download tests. The FASTQ records and metadata templates are encoded directly in `make.sh` so that the expected read names, sequences, and table shape are easy to inspect.

<br />

## Files
Readable provenance:
- `make.sh`

Compressed source FASTQ:
- `source/se/tiny_download_se.fastq.gz`
- `source/pe/tiny_download_pe_R1.fastq.gz`
- `source/pe/tiny_download_pe_R2.fastq.gz`
- `source/se/tiny_download_se_truncated.fastq.gz` (the single-end file without its 8-byte gzip trailer)

Metadata template:
- `metadata/local_se.template.tsv`
- `metadata/local_pe.template.tsv`
- `metadata/local_pe_duplicate.template.tsv`
- `metadata/local_pe_duplicate_custom.template.tsv`
- `metadata/local_pe_conflicting_accession.template.tsv`
- `metadata/local_mixed.template.tsv`
- `metadata/local_se_verified.template.tsv`
- `metadata/local_pe_verified.template.tsv`
- `metadata/local_pe_md5_only.template.tsv`
- `metadata/local_se_wrong_md5.template.tsv`
- `metadata/local_se_wrong_size.template.tsv`
- `metadata/local_pe_md5_count.template.tsv`
- `metadata/local_pe_wrong_md5_r2.template.tsv`
- `metadata/local_pe_conflicting_md5.template.tsv`
- `metadata/local_se_truncated.template.tsv`
- `metadata/local_se_missing.template.tsv`

The compressed FASTQ files are generated with `gzip -n -c` so gzip headers do not encode timestamps or source filenames.

The metadata template uses the columns required by `execute_download_fastqs.sh`:
```text
run_accession	custom_name	fastq_https
SRR_LOCAL_SE	tiny_download_se	__BASE_URL__/se/tiny_download_se.fastq.gz
SRR_LOCAL_PE	tiny_download_pe	__BASE_URL__/pe/tiny_download_pe_R1.fastq.gz;__BASE_URL__/pe/tiny_download_pe_R2.fastq.gz
```

The test replaces `__BASE_URL__` at runtime with a loopback HTTP URL served from `127.0.0.1`. This avoids external network access while still exercising real `wget` download behavior.

`metadata/local_mixed.template.tsv` contains both rows in one table to exercise the case where users provide a metadata TSV containing SE and PE data.

`metadata/local_pe_duplicate.template.tsv` repeats one accession/URL pair with distinct custom names. It exercises duplicate-aware download planning: one raw accession download is reused for multiple logical custom-name symlinks. The malformed duplicate templates exercise pre-download rejection for duplicate custom names and one accession paired with conflicting URLs.

<br />

## Expected download behavior
The single-end FASTQ source contains one clean 64-bp read:

| Read name                 | Sequence                                                           |
| :---                      | :---                                                               |
| `tiny_download_se_read_1` | `GATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAG` |

`execute_download_fastqs.sh` should download the source URL to `SRR_LOCAL_SE.fastq.gz`, create the custom symlink `tiny_download_se.fastq.gz`, preserve `gzip` integrity, and preserve the expected read name.

The paired-end FASTQ sources contain one clean 64-bp read pair:

| Read name                   | Mate | Sequence                                                           |
| :---                        | :--- | :---                                                               |
| `tiny_download_pe_pair_1/1` | R1   | `ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT` |
| `tiny_download_pe_pair_1/2` | R2   | `TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA` |

`execute_download_fastqs.sh` should
- parse the semicolon-delimited PE URL field,
- download the sources to `SRR_LOCAL_PE_R1.fastq.gz` and `SRR_LOCAL_PE_R2.fastq.gz`,
- create custom symlinks `tiny_download_pe_R1.fastq.gz` and `tiny_download_pe_R2.fastq.gz`,
- preserve `gzip` integrity, and
- preserve the expected mate-specific read names.

When the same accession/URL pair appears more than once with different custom names, it should download the raw accession FASTQs once and create every requested custom-name symlink.

<br />

## Download verification
`make.sh` computes each source file's MD5 and size when it generates the fixtures and writes them into the verification templates, which add the optional `fastq_md5` and `fastq_bytes` columns:

| Template                   | Content                                  | Expected behavior                                        |
| :---                       | :---                                     | :---                                                     |
| `local_se_verified`        | right MD5 and size                       | downloads; no unverified-download warning                |
| `local_pe_verified`        | right MD5 and size for both mates        | downloads both mates                                     |
| `local_pe_md5_only`        | right MD5s, no `fastq_bytes` column      | downloads both mates                                     |
| `local_se_wrong_md5`       | wrong MD5                                | fails after 3 attempts; no file left                     |
| `local_se_wrong_size`      | size one byte too large                  | fails after 3 attempts; no file left                     |
| `local_pe_md5_count`       | one MD5 for two URLs                     | rejected while reading the TSV                           |
| `local_pe_wrong_md5_r2`    | right R1 MD5, wrong R2 MD5               | R1 downloads; R2 fails after 3 attempts; no R2 file left |
| `local_pe_conflicting_md5` | one accession twice, with different MD5s | rejected while planning the downloads                    |
| `local_se_truncated`       | no MD5; the truncated source             | fails `gzip -t` after 3 attempts; no file left           |
| `local_se_missing`         | no MD5; a URL the server does not have   | fails at once (server 404); no file left                 |

`local_se.template.tsv`, with no `fastq_md5` column, gives one warning that its download is checked only with `gzip -t`.

<br />

## Current and deferred test coverage
The test suite covers local single-end, paired-end, and mixed-metadata downloads, as well as download verification, over loopback HTTP when gated by `RUN_DOWNLOAD=1`. It also covers local GNU Parallel dispatch for mixed metadata when gated by both `RUN_DOWNLOAD=1` and `RUN_PARALLEL=1`.

`#TODO`: Add Slurm coverage, and consider an optional external-network test gated by a future `RUN_NETWORK=1` flag. Handle issues when running tests while `env_protocol` is active.
