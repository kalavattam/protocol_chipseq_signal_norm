# Download-fastqs test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic no-network tests of the download-fastqs workflow.

They are intentionally small, hand-checkable, and version-controlled directly in Git. Running `scripts/make_fixtures.sh` regenerates the fixture set deterministically.

The fixture set currently focuses on single-end and paired-end local-download smoke tests. The FASTQ records and metadata templates are encoded directly in `scripts/make_fixtures.sh` so that the expected read names, sequences, and table shape are easy to inspect.

<br />

## Files
Readable provenance:
- `scripts/make_fixtures.sh`

Regenerate fixtures from the repository root with:

```bash
bash tests/download_fastqs/scripts/make_fixtures.sh
```

Compressed source FASTQ:
- `source/tiny_download_se.fastq.gz`
- `source/tiny_download_pe_R1.fastq.gz`
- `source/tiny_download_pe_R2.fastq.gz`

Metadata template:
- `metadata/local_se.template.tsv`
- `metadata/local_pe.template.tsv`
- `metadata/local_mixed.template.tsv`

The compressed FASTQ files are generated with `gzip -n -c` so gzip headers do not encode timestamps or source filenames. The uncompressed FASTQ files are not committed.

The metadata template uses the columns required by `execute_download_fastqs.sh`:

```text
run_accession	custom_name	fastq_https
SRR_LOCAL_SE	tiny_download_se	__BASE_URL__/tiny_download_se.fastq.gz
SRR_LOCAL_PE	tiny_download_pe	__BASE_URL__/tiny_download_pe_R1.fastq.gz;__BASE_URL__/tiny_download_pe_R2.fastq.gz
```

The smoke test replaces `__BASE_URL__` at runtime with a loopback HTTP URL served from `127.0.0.1`. This keeps the wrapper-backed download test no-network while still exercising real `wget` download behavior.

`metadata/local_mixed.template.tsv` contains both rows in one table to exercise the common case where users provide one mixed SE/PE metadata TSV.

<br />

## Expected download behavior
The single-end FASTQ source contains one clean 64 bp read:

| Read name                   | Sequence                                                             |
| :---                        | :---                                                                 |
| `tiny_download_se_read_1` | `GATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAGGATTACAG` |

`execute_download_fastqs.sh` should download the source URL to `SRR_LOCAL_SE.fastq.gz`, create the custom symlink `tiny_download_se.fastq.gz`, preserve gzip integrity, and preserve the expected read name.

The paired-end FASTQ sources contain one clean 64 bp read pair:

| Read name                    | Mate | Sequence                                                         |
| :---                         | :--- | :---                                                             |
| `tiny_download_pe_pair_1/1` | R1   | `ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT` |
| `tiny_download_pe_pair_1/2` | R2   | `TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA` |

`execute_download_fastqs.sh` should parse the semicolon-delimited PE URL field, download the sources to `SRR_LOCAL_PE_R1.fastq.gz` and `SRR_LOCAL_PE_R2.fastq.gz`, create custom symlinks `tiny_download_pe_R1.fastq.gz` and `tiny_download_pe_R2.fastq.gz`, preserve gzip integrity, and preserve the expected mate-specific read names.

<br />

## Deferred fixture batches
Later download-fastqs batches should add execute-level multiple-record metadata coverage, local GNU Parallel coverage if applicable, optional Slurm coverage, and optional external-network coverage gated by `RUN_NETWORK=1`.
