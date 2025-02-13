# protocol_chipseq_signal_norm

**Code and documentation for the *Bio-protocol* manuscript, "ChIP-seq data processing and relative and quantitative signal normalization for *Saccharomyces cerevisiae*."**

Version 1 of the manuscript is available as a [preprint](https://www.bio-protocol.org/exchange/preprintdetail?type=3&id=2770); version 2 is under review.

Workflow details are documented in [`workflow.md`](./workflow.md). Wrangling genome files is covered in [`download_process_fasta_gff3.md`](./download_process_fasta_gff3.md). Validation of the new Python implementation of siQ-ChIP is in [`validate_siqchip.md`](./validate_siqchip.md).

***Note:*** We are actively testing, extending, and refining this codebase. If you're using this repository, please run `git pull` regularly (e.g., daily) to stay up to date with the latest improvements and fixes.

***Important:*** If you cloned this repo prior to 2024-0212, you must clone it again. Large files were purged from the commit history, altering commit hashes and preventing `git pull` from resolving.
