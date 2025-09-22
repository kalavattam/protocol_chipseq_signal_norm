# protocol_chipseq_signal_norm

**Code and documentation for the *Bio-protocol* manuscript “ChIP-seq data processing and relative and quantitative signal normalization for *Saccharomyces cerevisiae*.”**

- Version 1: Available as a [preprint](https://www.bio-protocol.org/exchange/preprintdetail?type=3&id=2770) 2024-12-16.
- Version 2: [Published](https://doi.org/10.21769/BioProtoc.5299) 2025-05-05 in [*Bio-protocol* (volume 15, issue 9)](https://bio-protocol.org/en/archive?vol=15&issid=1370). A PDF with live internal links is available [here](./docs/protocol_chipseq_signal_norm.pdf).

***Note:** The protocol has changed substantially since version 1; much of the code no longer implements what is described there.*
<br />
<br />

## Repository status
***This codebase is being tested, extended, and refined.*** If you use this repository, run `git pull` regularly to stay up to date with improvements and fixes.

***If you cloned this repo before 2025-02-19, please clone it again.*** This is because large and otherwise unnecessary files were purged from the commit history, altering commit hashes. These changes prevent `git pull` from resolving updates correctly.
<br />
<br />

## Workflow documentation
- Regularly updated workflow details: [`workflow.md`](./workflow.md).
- Genome file processing: [`download_process_fasta_gff3.md`](./download_process_fasta_gff3.md).
- Validation of the Python implementation of siQ-ChIP (in progress): [`validate_siq_chip.md`](./validate_siq_chip.md).
    + [ ] Notebook needs cleanup and clearer documentation.
    + [ ] Rough figures will be added.
    + [x] Reference implementations:
        - [Original implementation](https://github.com/BradleyDickson/siQ-ChIP).
        - [Adapted for *S. cerevisiae*](https://github.com/kalavattam/siQ-ChIP/tree/protocol).

***Note:** SLURM job execution has undergone the most testing and is currently the most stable. Local and remote execution of parallelized (via [GNU Parallel](https://www.gnu.org/software/parallel/)) and serial jobs are still being refactored and tested.*
<br />
<br />

## Citation and licenses
If you use this repository or the associated protocol in your work, please cite the following:

Alavattam KG, Dickson BM, Hirano R, Dell R, Tsukiyama T. ChIP-seq Data Processing and Relative and Quantitative Signal Normalization for Saccharomyces cerevisiae. *Bio Protoc.* [2025 May 5;15(9):e5299](https://bio-protocol.org/en/archive?vol=15&issid=1370). doi: [10.21769/BioProtoc.5299](https://doi.org/10.21769/BioProtoc.5299). PMID: [40364978](https://pubmed.ncbi.nlm.nih.gov/40364978/); PMCID: [PMC12067309](https://pmc.ncbi.nlm.nih.gov/articles/PMC12067309/).

The preprint and published protocols are available under the [CC BY-NC 4.0 license](https://creativecommons.org/licenses/by-nc/4.0/). All code, scripts, and documentation in this repository are provided under the [MIT License](./LICENSE).
<br />
<br />

## Support
### Emails
- Dry-lab workflow (this repository): Kris Alavattam at <i>kalavat&#8203;tam (at) g&#8203;mail (dot) c&#8203;om</i> or <i>kal&#8203;avatt (at) fre&#8203;dhutch (dot) o&#8203;rg</i>.
- General siQ-ChIP information and the [original implementation](https://github.com/BradleyDickson/siQ-ChIP): Brad Dickson at <i>br&#8203;adley (dot) dick&#8203;son (at) va&#8203;i (dot) or&#8203;g</i>.
- Benchwork, yeast strains, and other materials: Toshi Tsukiyama at <i>tts&#8203;ukiya (at) fredhut&#8203;ch (dot) o&#8203;rg</i>.

### Issues
If you encounter an issue (bugs, broken code, broken links, unexpected behavior, unclear writing, etc.), please open a [GitHub Issue](https://github.com/kalavattam/protocol_chipseq_signal_norm/issues).

<details>
<summary><i>Details and tips for filing an issue</i></summary>
<br />

**Before filing**
- Make sure you've pulled the latest code: `git pull`.
- Try again with `--dry-run` and/or `--verbose`.
- Check [`workflow.md`](./workflow.md) and [existing issues](https://github.com/kalavattam/protocol_chipseq_signal_norm/issues) for known problems.
<br />

**When filing, it’s good to include/do the following:**
- What you ran (full command(s) and options)
- What you observed and what you expected (including error messages, logs, unexpected output, etc.)
- Your setup:
    + OS with version info
    + Conda and/or Mamba version info (include both if you have them installed)
    + Output of `conda list` or `mamba list` (or equivalent) showing installed packages and versions
    + Tools (`samtools`, `python`, etc.), including version info (particularly if not covered above)
    + Whether you ran locally [serial or with GNU Parallel (with version info)] or on SLURM (with version info)
    + Etc.
- For readability, please format code (or any plain text such as command outputs) with a Markdown or HTML method so that it renders in a [fixed-width font](https://en.wikipedia.org/wiki/Monospaced_font). For example:
    + Format [inline code](https://www.markdownguide.org/basic-syntax/#code) with single backticks (\`).
    + For longer snippets, use triple backticks (\`\`\`) for [fenced code blocks](https://www.markdownguide.org/extended-syntax/#fenced-code-blocks).
<br />

**Optional (but helpful)**
- Minimal test data or recipe to reproduce *(but please don’t upload large, sensitive, and/or proprietary data; use downsampled, public, and/or synthetic data instead)*
- `--dry-run` output for the same command
<br />

**If it’s not about code**
- Broken link: Mention the file/section where you found it, and provide the correct link (if you know).
- Unclear writing: Point to the section, sentence, or phrase, with a quick note about what is unclear.
- Unexpected behavior/output: Screenshots or copy-paste snippets are helpful.
</details>
