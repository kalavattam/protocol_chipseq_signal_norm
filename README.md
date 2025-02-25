# protocol_chipseq_signal_norm

**Code and documentation for the *Bio-protocol* manuscript "ChIP-seq data processing and relative and quantitative signal normalization for *Saccharomyces cerevisiae*."**

Version 1: Available as a [preprint](https://www.bio-protocol.org/exchange/preprintdetail?type=3&id=2770). Version 2: Currently under review.

***Note:** The protocol has changed substantially since version 1, and some of the code no longer implements what is described (e.g., deepTools has been removed for generating siQ-ChIP and spike-in scaled signal tracks).*
<br />
<br />

## Workflow documentation
- Regularly updated workflow details: [`workflow.md`](./workflow.md).
- Genome file processing: [`download_process_fasta_gff3.md`](./download_process_fasta_gff3.md).
- Validation of the Python implementation of siQ-ChIP (in progress): [`validate_siq_chip.md`](./validate_siq_chip.md).
    + Notebook is in progress&mdash;needs cleanup and better documentation.
    + Rough figures still need to be added.

***Note:** SLURM job execution has undergone the most testing and is considered the most stable. Local and remote execution of parallelized (via [GNU Parallel](https://www.gnu.org/software/parallel/)) and serial jobs are still being refactored and tested.*
<br />
<br />

## Repository status
***This codebase is actively being tested, extended, and refined.***  

If you are using this repository, please run `git pull` regularly (e.g., daily) to stay up to date with the latest improvements and fixes.
<br />
<br />

## Important notice
***If you cloned this repo before <mark>2025-0219</mark>, you must clone it again.***  

This is because large and otherwise unnecessary files were purged from the commit history, altering commit hashes. These changes prevent `git pull` from resolving updates correctly.
