
# `align_fastqs`: the files `--unal_mode separate` writes
`align_fastqs` aligns reads and keeps the paired alignments (or, for single-end data, individual alignments) that pass its filters. `--unal_mode` decides what happens to the rest. Under `remove`, the default mode, they are dropped. Under `separate`, nothing is dropped: every record the aligner writes lands in exactly one file, so the files together describe the whole library.

This page documents the files, the rules that decide which file a record lands in, the two duplicate markings, and the counts in `<stem>.diag.tsv`. The helper's help (`align_fastqs --help`) and the drivers' help provide the short form and point here.

<br />

## Files
With `--out_ext cram` (or a `--fil_out` ending in `.cram`), every alignment file below is CRAM instead of BAM (never both). Each alignment file has an index (`.bai` or `.crai`) and, with `--qname`, a queryname-sorted copy (with `.qnam.bam` or `.qnam.cram` in place of the extension).

Paired-end data:

| File                        | `remove` | `separate` | What it holds                                                                        |
| :---                        | :---:    | :---:      | :---                                                                                 |
| `<stem>.bam`                | yes      | yes        | The main file: aligned pairs that pass every filter, identical in both modes         |
| `<stem>.unal.bam`           | no       | yes        | Unaligned reads and their mates: pairs with both mates unaligned, and pairs with one |
| `<stem>.improp.bam`         | no       | yes        | Pairs with both mates aligned that `--req_flg` rejects; empty without `--req_flg`    |
| `<stem>.mapq.bam`           | no       | yes        | Templates that the `--mapq` filter rejects and `--req_flg` does not                  |
| `<stem>.markdup.txt.gz`     | yes      | yes        | Duplicate statistics for the main file alone (Markdup 1)                             |
| `<stem>.all.markdup.txt.gz` | no       | yes        | Duplicate statistics for all four alignment files together (Markdup 2)               |
| `<stem>.diag.tsv`           | no       | yes        | The counts, one row per file and category                                            |

Single-end data has no pairs, so it has no `.improp` file. `--req_flg` is ignored with a warning when every entry is single-end. `<stem>.unal.bam` holds the reads flagged unaligned (SAM flag 4), `<stem>.mapq.bam` holds the reads the `--mapq` filter rejects, and Markdup 2 covers the three alignment files.

In a submission that mixes layouts, each entry follows its own layout.

The `.unal`, `.improp`, and `.mapq` files are diagnostics. They carry the pipeline record, with a line saying what each holds and a warning that it is not for signal or scaling.

<br />

## Which file a record lands in
A template is all the records that share a read name: both mates of a pair, with any secondary or supplementary records. Every rule below moves whole templates, so no template has records in two files.

1. **Unaligned reads and their mates** go to `.unal`. In paired-end data, this is any template with a record flagged unaligned (0x4) or with an unaligned mate (0x8). When BWA-backtrack flags one mate unaligned but leaves its partner's mate flag clear ([#75](https://github.com/kalavattam/protocol_chipseq_signal_norm/issues/75)), the partner joins its mate in `.unal`.
2. **Improper pairs**, with `--req_flg` and paired-end data, go to `.improp`. A pair is improper when the aligner did not flag it proper (SAM flag 2), or when `samtools fixmate` clears the flag. `fixmate` clears it on a pair whose mates map to different references and on a pair whose mates align in the same orientation (both forward or both reverse). BWA-backtrack's `bwa sampe` judges pairs by distance in the concatenated genome and never checks the reference, so it can flag a pair proper whose mates sit on either side of a chromosome boundary ([lh3/bwa#460](https://github.com/lh3/bwa/issues/460)).
3. **Low-MAPQ templates** go to `.mapq`. A template passes `--mapq` when each of its primary records (two for a pair, one for a single-end read) reaches the threshold. Its secondary and supplementary records follow it.
4. **Everything else** is the main file.

The rules run in this order. A template that fails both `--req_flg` and `--mapq` goes to `.improp`, never to both; `.diag.tsv` counts how many of the improper templates have a primary record below `--mapq` (the `low_mapq` row).

Under `remove`, rules 1 to 3 drop the templates instead. Both modes apply the same filters in the same order, so the main file holds the same records in both.

`--req_flg` filters twice. It keeps the records the aligner flagged proper before `fixmate`, so the aligner's own judgment decides, and then keeps the templates whose primary records still carry the flag after `fixmate`.

<br />

## Duplicate marking
`separate` marks duplicates twice.

- **Markdup 1** runs on the main file alone, exactly as under `remove`. The main file's duplicate flags and `<stem>.markdup.txt.gz` are therefore identical in both modes.
- **Markdup 2** runs on all the alignment files together. The `.unal`, `.improp`, and `.mapq` files keep its flags, so a read in one of them can be marked as a duplicate of a read in the main file. `<stem>.all.markdup.txt.gz` holds its statistics, which describe the whole library. The main file's records are marked in this pass too, but only their counts are kept, in `.diag.tsv`.

Optical duplicates are told apart in Markdup 1 and Markdup 2 only when every read name it reads carries flow-cell coordinates of one format. Each statistics file's last line, `OPTICAL DISTANCE:`, says whether they were.

`samtools markdup` judges duplicates by alignment position. It cannot judge an unaligned read, so a pair with both mates unaligned is never marked, optical duplicates among such pairs included. In a pair with one mate unaligned, only the aligned mate can be marked. `.diag.tsv` counts the records `markdup` cannot judge (the `unjudged` row). A duplication estimate for such reads has to come from the reads themselves, not from their alignments.

<br />

## `<stem>.diag.tsv`
A tab-separated table with a header line and one row per file and category. The table has the same rows for every run of a layout, whether or not a count is zero, so tables from different runs line up.

| Column      | Meaning                                                               |
| :---        | :---                                                                  |
| `file`      | The alignment file: `main`, `unal`, `improp`, or `mapq`               |
| `category`  | What is counted (rows below)                                          |
| `templates` | Templates with at least one record in the category                    |
| `records`   | Records in the category, secondary and supplementary records included |

Rows, in order, for paired-end data:

| `file`   | `category`       | Counts                                                                                   |
| :---     | :---             | :---                                                                                     |
| `main`   | `all`            | Every template and record in the file                                                    |
| `main`   | `duplicate_all`  | The file's templates and records that Markdup 2 marks as duplicates                      |
| `unal`   | `all`            | Every template and record in the file                                                    |
| `unal`   | `both_unaligned` | Pairs with both mates unaligned                                                          |
| `unal`   | `one_unaligned`  | Pairs with one mate aligned                                                              |
| `unal`   | `unjudged`       | Unaligned records (flag 4), which `markdup` cannot judge, and the templates holding them |
| `unal`   | `duplicate_all`  | Records marked as duplicates by Markdup 2, and the templates holding them                |
| `improp` | `all`            | Every template and record in the file                                                    |
| `improp` | `diff_ref`       | Pairs whose mates map to different references                                            |
| `improp` | `same_strand`    | Other pairs whose mates align in the same orientation (both forward or both reverse)     |
| `improp` | `other`          | The remaining improper pairs (for example, too far apart or facing outward)              |
| `improp` | `low_mapq`       | Improper pairs with a primary record below `--mapq`, which would also have failed it     |
| `improp` | `duplicate_all`  | As for `unal`                                                                            |
| `mapq`   | `all`            | Every template and record in the file                                                    |
| `mapq`   | `duplicate_all`  | As for `unal`                                                                            |

For single-end data, the rows are `main all`, `main duplicate_all`, `unal all`, `unal unjudged`, `unal duplicate_all`, `mapq all`, and `mapq duplicate_all`.

The `all` rows' records sum to the records the aligner wrote. `both_unaligned` and `one_unaligned` add up to `unal all`, and so do `diff_ref`, `same_strand`, and `other` to `improp all`. The reasons are judged in that order, from the primary records: a pair whose mates map to different references is counted there, whatever their strands.

The main file's `duplicate_all` row comes from Markdup 2's statistics: its records are the duplicate total less the other files' duplicate records, and its templates are half of its duplicate primary records (`markdup` marks both mates of a duplicate pair) or, for single-end data, all of them.
