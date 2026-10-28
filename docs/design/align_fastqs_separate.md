
# `align_fastqs`: the files `--unal_mode separate` writes
`align_fastqs` aligns reads and keeps the paired alignments (or, for single-end data, individual alignments) that pass its filters. `--unal_mode` decides what happens to the rest. Under `remove`, the default mode, they are dropped. Under `separate`, nothing is dropped: every record the aligner writes lands in exactly one file, so the files together describe the whole library.

This page documents the files, the rules that decide which file a record lands in, the two duplicate markings, the counts in `<stem>.diag.tsv.gz`, and how each aligner judges a pair. The helper's help (`align_fastqs --help`) and the drivers' help provide the short form and point here.

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
| `<stem>.diag.tsv.gz`        | no       | yes        | The counts, one row per file and category                                            |

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

The rules run in this order: a template that fails both `--req_flg` and `--mapq` goes to `.improp`, never to both; `.diag.tsv.gz` counts how many of the improper templates have a primary record below `--mapq` (the `templates_mapq_low` and `records_mapq_low` columns).

Under `remove`, rules 1 to 3 drop the templates instead. Both modes apply the same filters in the same order, so the main file holds the same records in both.

`--req_flg` filters twice. It keeps the records the aligner flagged proper before `fixmate`, so the aligner's own judgment decides, and then keeps the templates whose primary records still carry the flag after `fixmate`.

<br />

## Duplicate marking
`separate` marks duplicates twice.

- **Markdup 1** runs on the main file alone, exactly as under `remove`. The main file's duplicate flags and `<stem>.markdup.txt.gz` are therefore identical in both modes.
- **Markdup 2** runs on all the alignment files together. The `.unal`, `.improp`, and `.mapq` files keep its flags, so a read in one of them can be marked as a duplicate of a read in the main file. `<stem>.all.markdup.txt.gz` holds its statistics, which describe the whole library. The main file's records are marked in this pass too, but only their counts are kept, in `.diag.tsv.gz`.

Optical duplicates are told apart in Markdup 1 and Markdup 2 only when every read name it reads carries flow-cell coordinates of one format. Each statistics file's last line, `OPTICAL DISTANCE:`, says whether they were.

`samtools markdup` judges duplicates by alignment position. It cannot judge an unaligned read, so a pair with both mates unaligned is never marked, optical duplicates among such pairs included. In a pair with one mate unaligned, only the aligned mate can be marked. `.diag.tsv.gz` counts, in every file, the records `markdup` skips: these unaligned records, and secondary and supplementary records (the `unjudged` rows). A duplication estimate for such reads has to come from the reads themselves, not from their alignments.

<br />

## `<stem>.diag.tsv.gz`
A gzip-compressed, tab-separated table with a header line and one row per file and category. The table has the same rows for every run of a layout, whether or not a count is zero, so tables from different runs line up. A cell is NA when the run cannot measure it.

| Column               | Meaning                                                                                                                 |
| :---                 | :---                                                                                                                    |
| `file`               | The alignment file: `main`, `unal`, `improp`, or `mapq`                                                                 |
| `category`           | What is counted (rows below)                                                                                            |
| `templates`          | Templates with at least one record in the category                                                                      |
| `records`            | Records in the category, secondary and supplementary records included                                                   |
| `templates_mapq_low` | Of those templates, the ones with a primary record below `--mapq`, which would also have failed it; NA outside `improp` |
| `records_mapq_low`   | Their records in the category; NA outside `improp`                                                                      |

Rows, in order, for paired-end data:

| `file`   | `category`       | Counts                                                                                                                             |
| :---     | :---             | :---                                                                                                                               |
| `main`   | `all`            | Every template and record in the file                                                                                              |
| `main`   | `unjudged`       | Records that `markdup` skips (unaligned, 0x4; secondary, 0x100; supplementary, 0x800), and the templates holding them              |
| `main`   | `duplicate_all`  | The file's templates and records that Markdup 2 marks as duplicates                                                                |
| `unal`   | `all`            | Every template and record in the file                                                                                              |
| `unal`   | `both_unaligned` | Pairs with both mates unaligned                                                                                                    |
| `unal`   | `one_unaligned`  | Pairs with one mate aligned                                                                                                        |
| `unal`   | `unjudged`       | As for `main`                                                                                                                      |
| `unal`   | `duplicate_all`  | Records marked as duplicates by Markdup 2, and the templates holding them                                                          |
| `improp` | `all`            | Every template and record in the file                                                                                              |
| `improp` | `diff_ref`       | Pairs whose mates map to different references                                                                                      |
| `improp` | `orient_ff`      | Pairs with both mates forward                                                                                                      |
| `improp` | `orient_rr`      | Pairs with both mates reverse                                                                                                      |
| `improp` | `orient_rf`      | Pairs facing outward: the reverse mate ends before the forward mate starts                                                         |
| `improp` | `dovetailed`     | Pairs in which one mate extends past the other: the reverse mate starts before the forward mate, or the forward mate ends after it |
| `improp` | `too_far`        | Pairs that span more than Bowtie 2's `-X` (`--bt2_X`); NA for BWA                                                                  |
| `improp` | `too_close`      | Pairs that span less than Bowtie 2's `-I`, which is 0 in every call, so the count is always 0; NA for BWA                          |
| `improp` | `failed_est`     | Pairs improper because the aligner's insert-size estimate failed; NA for every aligner for now                                     |
| `improp` | `other`          | The remaining improper pairs                                                                                                       |
| `improp` | `unjudged`       | As for `main`                                                                                                                      |
| `improp` | `duplicate_all`  | As for `unal`                                                                                                                      |
| `mapq`   | `all`            | Every template and record in the file                                                                                              |
| `mapq`   | `unjudged`       | As for `main`                                                                                                                      |
| `mapq`   | `duplicate_all`  | As for `unal`                                                                                                                      |

For single-end data, the rows are `main all`, `main unjudged`, `main duplicate_all`, `unal all`, `unal unjudged`, `unal duplicate_all`, `mapq all`, `mapq unjudged`, and `mapq duplicate_all`.

The `all` rows' records sum to the records the aligner wrote. `both_unaligned` and `one_unaligned` add up to `unal all`. The categories from `diff_ref` to `other` add up to `improp all` in every column, with NA counted as 0, since each improper pair is counted in exactly one of them.

An improper pair's category is judged from its first mate's primary record, which gives the mate's position and, through the mate CIGAR (`MC`, written by `samtools fixmate -m`), the mate's end. The categories are tried in the order of the rows, and the first that applies counts: a pair whose mates map to different references is counted there, whatever their orientation. A pair's span runs from the leftmost to the rightmost aligned base of its two mates, soft clips excluded, as Bowtie 2 measures it. As Bowtie 2 does, `too_far` raises `-X` to a mate's length when the mate is longer. One `gawk` pass (`lib/awk/classify_improp_pairs.awk`) counts every `improp` row.

The categories describe each pair, not the aligner's verdict, and each aligner judges pairs its own way (next section). So `other` holds the pairs that are improper for a reason the table does not measure:
- With BWA, mostly pairs outside the aligner's insert-size bounds, which it sets per batch ([#85](https://github.com/kalavattam/protocol_chipseq_signal_norm/issues/85)).
- With Bowtie 2's relaxed call, also pairs that Bowtie 2 writes as two unpaired alignments (`YT:Z:UP`) whose positions happen to look concordant.
- With the default Bowtie 2 call, nothing: `--no-mixed --no-discordant` writes only concordant pairs, so the improper file is empty.

The main file's `duplicate_all` row comes from Markdup 2's statistics. Its records are the duplicate total less the other files' duplicate records. Its templates are half of its duplicate primary records (`markdup` marks both mates of a duplicate pair) or, for single-end data, all of them. If the subtraction comes out negative, or leaves an odd count of primary duplicates in paired-end data, it cannot be right: the row is NA, and the helper warns. Only a counting defect or a change in `samtools markdup` can cause this, so every other count stands.

<br />

## How each aligner judges a pair
Each aligner decides for itself whether a pair is proper (SAM flag 2) from its own source:

|                              | Bowtie 2 2.5.5 (`pe.cpp`)                                       | BWA-MEM 0.7.19 and BWA-MEM2 2.3 (`bwamem_pair`)                | BWA-backtrack 0.7.19 (`bwape.c`)                                   |
| :---                         | :---                                                            | :---                                                           | :---                                                               |
| What the span is measured on | Aligned extents only; soft clips are excluded                   | 5′ start of the forward mate to the 3′ end of the reverse mate | 5′ start of the forward mate to the 3′ end of the reverse mate     |
| FF and RR                    | Improper                                                        | Improper unless that orientation's estimate succeeded          | Improper                                                           |
| RF                           | Improper when the mates don't overlap                           | Its own orientation, usually failed, so improper               | Improper                                                           |
| Dovetailed                   | Improper (neither call passes `--dovetail`)                     | Counts as FR with a shorter span, or as RF                     | Improper (the reverse mate must start after the forward mate)      |
| Too far                      | Span > `-X` (`--bt2_X`)                                         | Outside its own `[low, high]`, set per batch                   | Above `high_bayesian`, or above `-a` (500) when the estimate fails |
| Too close                    | Span < `-I`, which is 0 in every call, so the count is always 0 | Per batch                                                      | Shorter than the longest read                                      |
| Failed estimate              | No estimate, so NA                                              | Per orientation and per batch                                  | Per batch, after which it falls back to `-a`                       |

`fixmate` then clears the flag on pairs whose mates map to different references or align in the same orientation (see [Which file a record lands in](#which-file-a-record-lands-in)).
