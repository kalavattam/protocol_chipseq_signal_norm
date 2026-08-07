
# `protocol_chipseq_signal_norm`
**Code and documentation accompanying the *Bio-protocol* manuscript “ChIP-seq data processing and relative and quantitative signal normalization for *Saccharomyces cerevisiae*.”**
- [Preprint version](https://www.bio-protocol.org/exchange/preprintdetail?type=3&id=2770) posted on 2024-12-16.
- [Published version](https://doi.org/10.21769/BioProtoc.5299) available in [*Bio-protocol* (volume 15, issue 9)](https://bio-protocol.org/en/archive?vol=15&issid=1370) on 2025-05-05. A [PDF with live internal links](./docs/protocol_chipseq_signal_norm.pdf) is also available.

`protocol_chipseq_signal_norm` combines a Python package with command-line workflows implemented primarily in Bash ≥ 4.4. It is designed for routine ChIP-seq data processing and analysis in the [Tsukiyama Lab](https://research.fredhutch.org/tsukiyama/en.html) at [Fred Hutch Cancer Center](https://www.fredhutch.org). A major focus is the implementation of methods for signal normalization within and across samples, including siQ-ChIP ([Dickson et al., *J Biol Chem* 2020](https://pubmed.ncbi.nlm.nih.gov/32994221/); [Dickson et al., *Sci Rep* 2023](https://pubmed.ncbi.nlm.nih.gov/37160995/)) and ChIP-Rx ([Orlando et al., *Cell Rep* 2014](https://pubmed.ncbi.nlm.nih.gov/25437568/)). `protocol_chipseq_signal_norm` supports semi-automated parallel or serial execution of methods described in the *Bio-protocol* manuscript, along with methods developed following publication.

***Note:** The protocol has changed substantially since the preprint; much of the code no longer implements what is described there.*

***Repository status:** This codebase is being tested, extended, and refined.*

<br />

## Workflow documentation
- Running the workflow: [`notebooks/workflows/workflow.md`](./notebooks/workflows/workflow.md).
- [Input-floor command reference and design](./docs/design/compute_input_floor.md).
- Genome file processing: [`notebooks/workflows/download_process_fasta_gff3.md`](./notebooks/workflows/download_process_fasta_gff3.md).
- Validation of the Python implementation of siQ-ChIP (in progress): [`notebooks/validation/validate_siqchip.md`](./notebooks/validation/validate_siqchip.md).
    + [ ] Notebook needs cleanup and clearer documentation.
    + [ ] Rough figures will be added.
    + [x] Reference implementations:
        - [Original implementation by Brad Dickson](https://github.com/BradleyDickson/siQ-ChIP).
        - [Implementation adapted for use with *S. cerevisiae* data (and more)](https://github.com/kalavattam/siQ-ChIP/tree/protocol).
- Validation of the ChIP-Rx implementation is planned.

<br />

## Installation
Installation support files are organized under [`install/`](./install/). Environment YAML files are stored in [`install/envs/`](./install/envs/), and installation scripts are stored in [`install/scripts/`](./install/scripts/).

To set up the main repository environment, start with the following:
```sh
#  After cloning/fetching the repo and cd'ing into it
sh install/scripts/install_envs_entrypoint.sh --env_nam env_protocol --yes
```

<details>
<summary><i>(Click to view <code>install_envs_entrypoint.sh</code> details.)</i></summary>
<br />

> This POSIX-compatible entrypoint is intended to be runnable from various shell environments, including cases where the required Bash ≥ 4.4 is not yet available in `PATH`.
>
> If Bash ≥ 4.4 is already available in `PATH` and Conda or Mamba is also available in `PATH`, `install_envs_entrypoint.sh` will hand off to `install_envs.sh`.
>
> If Bash ≥ 4.4 is not yet available, or if Conda or Mamba is not available in `PATH`, `install_envs_entrypoint.sh` will print guidance for the next setup step, including installing [Miniforge](https://github.com/conda-forge/miniforge) when needed.
>
> After Bash ≥ 4.4 has been installed and is available in `PATH`, along with Conda or Mamba, repository environments can also be installed directly with:
> ```bash
> bash install/scripts/install_envs.sh --env_nam env_protocol --yes
> ```
>
> User-facing environments currently supported by `install_envs.sh` are `env_protocol`, `env_analyze`, and `env_siqchip`. By default, these are created from the corresponding YAML files in `install/envs/`.
>
> For institutional or site-specific channel configurations, such as systems that use mirrored Conda channels, pass channels directly to the installer:
> ```bash
> bash install/scripts/install_envs.sh \
>     --env_nam env_protocol \
>     --channels "mirror-conda-forge,mirror-bioconda" \
>     --override_channels \
>     --yes
> ```
>
> The `--if_exists reuse` option skips environment creation and refreshes the managed editable package installation in an existing `env_protocol`:
> ```bash
> bash install/scripts/install_envs.sh \
>     --env_nam env_protocol \
>     --if_exists reuse
> ```
</details>
<br />

The installer exposes the maintained Python CLIs as underscore-named commands. The installed command and explicit module form share one implementation:
```bash
compute_pseudo --help
python -m protocol_chipseq_signal_norm.cli.compute_pseudo --help
```

Maintained shell entrypoints live in `bin/`, sourced Bash in `lib/bash/`, and installable Python in `src/protocol_chipseq_signal_norm/`. See [`docs/standards/README.md`](./docs/standards/README.md) for the governed layout and source standards.

<br />

After `env_protocol` has been created, Julia and Atria can be installed with:
```bash
bash install/scripts/install_atria.sh --if_exists reuse
```

<details>
<summary><i>(Click to view <code>install_atria.sh</code> details.)</i></summary>
<br />

> By default, `install_atria.sh` installs Julia 1.8.5 and Atria 4.1.5 under a user-controlled installation directory. It checks the active project environment for Atria runtime dependencies such as `pigz`, `pbzip2`, and `Rscript`, verifies Julia archive checksums, and can reuse matching existing Julia or Atria installations when `--if_exists reuse` is specified.
>
> To print the resolved installation plan without downloading, cloning, building, or writing `PATH` snippets, run:
> ```bash
> bash install/scripts/install_atria.sh --dry_run --if_exists reuse
> ```
>
> To append Julia and Atria `PATH` lines to a shell configuration or snippet file, use `--pth_snp`:
> ```bash
> bash install/scripts/install_atria.sh \
>     --if_exists reuse \
>     --pth_snp "${HOME}/.bash_profile"
> ```
>
> Use `--help` with either installer for the full list of supported options:
> ```bash
> bash install/scripts/install_envs.sh --help
> bash install/scripts/install_atria.sh --help
> ```
</details>
<br />

## Testing
Run the default test suite from the repository root:
```bash
bash tests/run_tests.sh
conda run -n env_protocol ruff check --no-fix src tests dev
conda run -n env_protocol ruff format --check src tests dev
```

Many workflow integration tests expect `env_protocol`, the canonical project environment. If that environment or a required dependency is unavailable, the affected tests skip with an explicit message. To run through the environment without relying on shell activation:
```bash
#  With Conda
conda run -n env_protocol bash tests/run_tests.sh

#  With Mamba
mamba run -n env_protocol bash tests/run_tests.sh
```

Interactive shell activation also works:
```bash
conda activate env_protocol
bash tests/run_tests.sh
```

The runner
- bootstraps missing generated fixtures,
- writes logs and temporary products under `artifacts/tests/`, and
- leaves generated fixture outputs ignored by Git.

See [`tests/README.md`](./tests/README.md) for
- fixture policy,
- cleanup commands,
- optional gates, and
- coverage details.

Optional dependency classes are enabled explicitly:
```bash
RUN_ATRIA=1 conda run -n env_protocol bash tests/run_tests.sh
RUN_DOWNLOAD=1 conda run -n env_protocol bash tests/run_tests.sh
RUN_PARALLEL=1 conda run -n env_protocol bash tests/run_tests.sh
RUN_SLURM=1 conda run -n env_protocol bash tests/run_tests.sh
```

For a broader local run with common non-Slurm optional gates:
```bash
RUN_ATRIA=1 RUN_DOWNLOAD=1 RUN_PARALLEL=1 conda run -n env_protocol bash tests/run_tests.sh
```

For example, to time a broader local run without the GNU Parallel-gated tests while running through the environment with Mamba:
```bash
time RUN_DOWNLOAD=1 RUN_ATRIA=1 mamba run -n env_protocol bash tests/run_tests.sh
```

<br />

## Notes on compute-signal engine selection
For normal signal-track generation, the default `chrom` engine is recommended. That said, `window` tends to perform similarly well and, under certain conditions, beats `chrom`. The following holds for *S. cerevisiae* test data:

| Input or goal                       | Recommended engine |
| :---                                | :---               |
| General use or uncertain input type | `chrom`            |
| CRAM input                          | `chrom`            |
| Large BAM input                     | `window`           |

<br />

### How the engines work
- `chrom` parallelizes indexed chromosome fetches and array-backed signal computation.
- `window` parallelizes indexed coordinate-window fetches for finer load balance and can be faster on large BAM inputs.
- Fragment-coordinate BED output ignores signal-engine settings.

<br />

## Citation and licenses
If you use this repository or the associated protocol in your work, please cite the following:

Alavattam KG, Dickson BM, Hirano R, Dell R, Tsukiyama T. ChIP-seq Data Processing and Relative and Quantitative Signal Normalization for Saccharomyces cerevisiae. *Bio Protoc.* [2025 May 5;15(9):e5299](https://bio-protocol.org/en/archive?vol=15&issid=1370). doi: [10.21769/BioProtoc.5299](https://doi.org/10.21769/BioProtoc.5299). PMID: [40364978](https://pubmed.ncbi.nlm.nih.gov/40364978/); PMCID: [PMC12067309](https://pmc.ncbi.nlm.nih.gov/articles/PMC12067309/).

The preprint and published protocols are available under the [CC BY-NC 4.0 license](https://creativecommons.org/licenses/by-nc/4.0/). All code, scripts, documentation, etc. in this repository are provided under the [MIT License](./LICENSE).

<br />

## Support
### Emails
- Dry-lab workflow (this repository): Kris Alavattam at <i>kalavat&#8203;tam (at) g&#8203;mail (dot) c&#8203;om</i> or <i>kal&#8203;avatt (at) fre&#8203;dhutch (dot) o&#8203;rg</i>.
- General siQ-ChIP information and the [original implementation](https://github.com/BradleyDickson/siQ-ChIP): Brad Dickson at <i>br&#8203;adley (dot) dick&#8203;son (at) va&#8203;i (dot) or&#8203;g</i>.
- Benchwork, yeast strains, and other materials: Toshi Tsukiyama at <i>tts&#8203;ukiya (at) fredhut&#8203;ch (dot) o&#8203;rg</i>.

<br />

### Issues
If you encounter an issue (bugs, broken code, broken links, unexpected behavior, unclear writing, etc.), please open a [GitHub Issue](https://github.com/kalavattam/protocol_chipseq_signal_norm/issues).

<details>
<summary><i>(Click to see details and tips for filing an issue.)</i></summary>
<br />

> **Before filing**
> - Make sure you've pulled the latest code: `git pull`.
> - Try again in dry-run mode and/or with `--verbose`.
> - Check [`notebooks/workflows/workflow.md`](./notebooks/workflows/workflow.md) and [existing issues](https://github.com/kalavattam/protocol_chipseq_signal_norm/issues) for known problems.
>
> <br />
>
> **When filing, it’s good to include/do the following:**
> - What you ran (full command(s) and options).
> - What you observed and what you expected (including error messages, logs, unexpected output, etc.).
> - Your setup:
>     + OS with version info.
>     + Conda and/or Mamba version info (include both if you have them installed).
>     + Output of `conda list` or `mamba list` (or equivalent) showing installed packages and versions.
>     + Tools (`samtools`, `python`, etc.), including version info (particularly if not covered above).
>     + Whether you ran locally [serial or with GNU Parallel (with version info)] or on Slurm (with version info).
>     + Etc.
> - For readability, please format code (or any plain text such as command outputs) with a Markdown or HTML method so that it renders in a [fixed-width font](https://en.wikipedia.org/wiki/Monospaced_font). For example:
>     + Format [inline code](https://www.markdownguide.org/basic-syntax/#code) with single backticks (\`...\`) or appropriate HTML tags (\<code\>...\</code\>).
>     + For longer snippets, use triple (or more) backticks (e.g., \`\`\`) for [fenced code blocks](https://www.markdownguide.org/extended-syntax/#fenced-code-blocks).
>
> <br />
>
> **Optional (but helpful)**
> - Minimal test data or recipe to reproduce. *(But please don’t upload large, sensitive, and/or proprietary data; use public and/or synthetic data instead. Smaller data are preferable.)*
> - Dry-run mode output for the same command.
>
> <br />
>
> **If it’s not about code**
> - Broken link: mention the file/section where you found it, and provide the correct link (if you know).
> - Unclear writing: point to the section, sentence, or phrase, with a quick note about what is unclear.
> - Unexpected behavior/output: screenshots and/or copy-paste snippets are helpful.
</details>
<br />

## Authorship and AI assistance in code, documentation, and related materials
- Director, responsible author, and maintainer: Kris Alavattam.
- Code, documentation, and related materials are planned, designed, developed, reviewed, and revised with AI assistance from ChatGPT and Codex (OpenAI), and Claude Cowork and Claude Code (Anthropic). All AI-assisted output is reviewed, edited, and approved by the author.
    + ChatGPT: GPT-4- and GPT-5-series models (most recent: GPT-5.6 Terra and Sol).
    + Codex: GPT-5.4, GPT-5.5, and GPT-5.6 Terra and Sol models.
    + Claude Cowork: Opus 4.8 model.
    + Claude Code: Opus 5 model.
- Ongoing maintenance is performed with AI assistance from Codex (most recent: GPT-5.6 Terra and Sol) and Claude Code (most recent: Opus 5).
