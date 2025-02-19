# protocol_chipseq_signal_norm

**Code and documentation for the *Bio-protocol* manuscript,  
"ChIP-seq data processing and relative and quantitative signal normalization for *Saccharomyces cerevisiae*."**

Version 1 of the manuscript is available as a [preprint](https://www.bio-protocol.org/exchange/preprintdetail?type=3&id=2770); version 2 is under review.

## Workflow documentation
- Workflow details (regularly updated): [`workflow.md`](./workflow.md).
- Genome file processing: [`download_process_fasta_gff3.md`](./download_process_fasta_gff3.md).
- Validation of the Python implementation of siQ-ChIP (notebook in progress): [`validate_siqchip.md`](./validate_siqchip.md).

***Note:*** SLURM job execution has undergone the most testing and is considered the most stable. Local and remote execution of parallelized (via [GNU Parallel](https://www.gnu.org/software/parallel/)) and serial jobs are still being refactored and tested.

## Repository status
***This codebase is actively being tested, extended, and refined.***  

If you are using this repository, please run `git pull` regularly (e.g., daily) to stay up to date with the latest improvements and fixes.

## Important notice
***If you cloned this repo before <mark>2025-0219</mark>, you must clone it again.***  

This is because large and otherwise unnecessary files were purged from the commit history, altering commit hashes. These changes prevent `git pull` from resolving updates correctly.
