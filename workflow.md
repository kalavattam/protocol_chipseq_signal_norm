
ChIP-seq Protocol Workflow
==========================

**Supporting code and documentation for the [manuscript](https://www.researchgate.net/publication/390849439_ChIP-seq_Data_Processing_and_Relative_and_Quantitative_Signal_Normalization_for_Saccharomyces_cerevisiae) "ChIP-seq Data Processing and Relative and Quantitative Signal Normalization for *Saccharomyces cerevisiae*," published 2025-05-05 in [*Bio-protocol* (volume 15, issue 9)](https://bio-protocol.org/en/archive?vol=15&issid=1370).**

**Author:** *Kris Alavattam*

This notebook provides a guide to the ChIP-seq data processing workflow detailed in the manuscript, including code snippets, explanations, and step-by-step instructions.

*Note: If using a high-performance computing cluster (HPCC), request an interactive node to ensure adequate resources for running code in the below chunks. The specific command (e.g., `grabnode` at Fred Hutch Cancer Center) will depend on the job scheduler setup. (This step is unnecessary if running the code on a local machine.)*

*Note: For detailed instructions on keeping your local version of the [`protocol_chipseq_signal_norm`](https://github.com/kalavattam/protocol_chipseq_signal_norm) repository up-to-date, please see [this GitHub gist](https://gist.github.com/kalavattam/76f123011e8dcd77b445a72d23a64036).*

---
<br />

## Table of contents
<details>
<summary><i>Table of contents</i></summary>
<br />
<!-- MarkdownTOC -->

1. [Procedures](#procedures)
    1. [A. Install and configure Miniforge.](#a-install-and-configure-miniforge)
    1. [B. Clone the protocol repository and install the project environment.](#b-clone-the-protocol-repository-and-install-the-project-environment)
    1. [C. Install and configure Atria for FASTQ adapter and quality trimming.](#c-install-and-configure-atria-for-fastq-adapter-and-quality-trimming)
    1. [D. Install and configure Integrative Genomics Viewer \(IGV\).](#d-install-and-configure-integrative-genomics-viewer-igv)
1. [Data analysis](#data-analysis)
    1. [A. Prepare and concatenate FASTA and GFF3 files for model and spike-in organisms.](#a-prepare-and-concatenate-fasta-and-gff3-files-for-model-and-spike-in-organisms)
        1. [Text](#text)
    1. [B. Generate Bowtie 2 index files from the concatenated FASTA file.](#b-generate-bowtie-2-index-files-from-the-concatenated-fasta-file)
        1. [Text](#text-1)
        1. [Bash code #1](#bash-code-1)
        1. [Bash code #2](#bash-code-2)
    1. [C. Obtain and organize ChIP-seq FASTQ files.](#c-obtain-and-organize-chip-seq-fastq-files)
        1. [Text](#text-2)
        1. [Bash code](#bash-code)
    1. [D. Use Atria to perform adapter and quality trimming of sequenced reads.](#d-use-atria-to-perform-adapter-and-quality-trimming-of-sequenced-reads)
    1. [E. Align sequenced reads with Bowtie 2 and process the read alignments.](#e-align-sequenced-reads-with-bowtie-2-and-process-the-read-alignments)
    1. [F. Compute normalized coverage.](#f-compute-normalized-coverage)
        1. [Text](#text-3)
        1. [Bash code](#bash-code-1)
    1. [G. Construct sample tables recording computed scaling factors for normalization.](#g-construct-sample-tables-recording-computed-scaling-factors-for-normalization)
        1. [1. Construct sample table recording siQ-ChIP $\alpha$ scaling factors.](#1-construct-sample-table-recording-siq-chip-%24alpha%24-scaling-factors)
            1. [Text](#text-4)
            1. [Bash code](#bash-code-2)
        1. [2. Construct sample table recording spike-in scaling factors.](#2-construct-sample-table-recording-spike-in-scaling-factors)
            1. [Text](#text-5)
            1. [Bash code](#bash-code-3)
    1. [H. Compute log2 ratios of IP to input normalized coverage.](#h-compute-log2-ratios-of-ip-to-input-normalized-coverage)
        1. [Text](#text-6)
        1. [Bash code](#bash-code-4)
    1. [I. Compute signal with the siQ-ChIP method.](#i-compute-signal-with-the-siq-chip-method)
        1. [Text](#text-7)
            1. [Overview](#overview)
            1. [Steps performed in code chunk](#steps-performed-in-code-chunk)
            1. [Important note](#important-note)
        1. [Bash code](#bash-code-5)
    1. [J. Compute signal with the spike-in method.](#j-compute-signal-with-the-spike-in-method)
        1. [Text](#text-8)
        1. [Bash code](#bash-code-6)

<!-- /MarkdownTOC -->
</details>
<br />

<a id="procedures"></a>
## Procedures
<a id="a-install-and-configure-miniforge"></a>
### A. Install and configure Miniforge.
<details>
<summary><i>Text: Install and configure Miniforge.</i></summary>
<br />

[Miniforge](https://github.com/conda-forge/miniforge) is an open source tool for managing bioinformatics software in isolated environments&mdash;self-contained workspaces that prevent software conflicts. This protocol uses Miniforge to create and manage a single environment, `env_protocol`, which contains programs for read alignment, alignment processing, signal computation, and signal visualization. It is recommended to [uninstall other software managers (e.g., Anaconda)](https://stackoverflow.com/questions/42182706/how-to-uninstall-anaconda-completely-from-macos) before installing Miniforge.
</details>
<br />

<details>
<summary><i>Bash code: Install and configure Miniforge.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
debug=true

#  Define variable for base directory
dir_bas="${HOME}"  ## WARNING: Change as needed ##

#  Determine the appropriate Miniforge installation script to use based on
#+ operating system (OS) and system architecture
case $(uname -s) in
    Darwin) os="MacOSX" ;;
    Linux)  os="Linux"  ;;
    *) echo "Error: Unsupported operating system: '$(uname -s)'." >&2 ;;
esac

ar=$(uname -m)  # e.g., "x86_64" for Intel/AMD, "arm64" for ARM

#  Set Miniforge installation script URL and name
https="https://github.com/conda-forge/miniforge/releases/latest/download"
script="Miniforge3-${os}-${ar}.sh"

if ${debug:-false}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "debug=${debug}"
    echo ""
    echo "dir_bas=${dir_bas}"
    echo ""
    echo "os=${os}"
    echo "ar=${ar}"
    echo ""
    echo "https=${https}"
    echo "script=${script}"
    echo ""
    echo ""
fi


#  Download and install Miniforge ---------------------------------------------
#  Move to base directory
cd "${dir_bas}" || echo "Error: Failed to change directory: '${dir_bas}'." >&2

#  Download Miniforge installation script
if [[ ! -f "${script}" ]]; then curl -L -O "${https}/${script}"; fi

#  Run Miniforge installation script
bash "${script}"
```
</details>
<br />

<a id="b-clone-the-protocol-repository-and-install-the-project-environment"></a>
### B. Clone the protocol repository and install the project environment.
<details>
<summary><i>Text: Clone the protocol repository and install the project environment.</i></summary>
<br />

Use the script [`install_envs.sh`](./scripts/install_envs.sh) to set up the project environment, `env_protocol`.
</details>
<br />

<details>
<summary><i>Bash code: Clone the protocol repository and install project the environment.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
debug=true

dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
repo="protocol_chipseq_signal_norm"
https="https://github.com/kalavattam/${repo}.git"

env="env_protocol"

if ${debug:-false}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "debug=${debug}"
    echo ""
    echo "dir_bas=${dir_bas}"
    echo ""
    echo "repo=${repo}"
    echo "https=${https}"
    echo ""
    echo "env=${env}"
    echo ""
    echo ""
fi


#  Do the main work -----------------------------------------------------------
#  Make the base directory
mkdir -p "${dir_bas}"

#  Move to base directory
cd "${dir_bas}" || echo "Error: Failed to change directory: '${dir_bas}'." >&2

#  Clone the protocol repository
if [[ ! -d "${repo}" ]]; then git clone "${https}"; fi

#  Move to repository directory
cd "${repo}" || echo "Error: Failed to change directory: '${repo}'." >&2

#  Install Conda/Mamba environment with install_envs.sh
bash "scripts/install_envs.sh" --env_nam "${env}" --yes
```
</details>
<br />

<a id="c-install-and-configure-atria-for-fastq-adapter-and-quality-trimming"></a>
### C. Install and configure Atria for FASTQ adapter and quality trimming.
<details>
<summary><i>Text: Install and configure Atria for FASTQ adapter and quality trimming.</i></summary>
<br />

To promote accurate alignment of sequenced reads, it is important to remove adapter sequences and low-quality base calls from FASTQ files. For this, we use [Atria](https://github.com/cihga39871/Atria), a tool written in [Julia](https://julialang.org/) that excels in adapter and quality trimming. See the following code chunk to install and configure Atria.
</details>
<br />

<details>
<summary><i>Bash code: Install and configure Atria for FASTQ adapter and quality trimming.</i></summary>
<br />

*Note: The following code chunks should be run in order, as some later chunks depend on variables defined in earlier ones.*

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##
```

First, if necessary, download Julia 1.8.5 for the appropriate OS and system architecture from the [official page](https://julialang.org/downloads/oldreleases/).

```bash
#  Define variables -----------------------------------------------------------
debug=true

#  Define OS and architecture variables
os="linux"     # "mac"         ## WARNING: Change as needed ##
ar_s="x64"     # "aarch64"     ## WARNING: Change as needed ##
ar_l="x86_64"  # "macaarch64"  ## WARNING: Change as needed ##

#  Set variables for Julia version
ver_s="1.8"
ver_l="${ver_s}.5"

#  Define URL and tarball filename
https="https://julialang-s3.julialang.org/bin/${os}/${ar_s}/${ver_s}"
case "${os}" in
    linux) tarball="julia-${ver_l}-${os}-${ar_l}.tar.gz" ;;
      mac) tarball="julia-${ver_l}-${ar_l}.tar.gz" ;;
esac

if ${debug:-false}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "debug=${debug}"
    echo ""
    echo "os=${os}"
    echo "ar_s=${ar_s}"
    echo "ar_l=${ar_l}"
    echo ""
    echo "ver_s=${ver_s}"
    echo "ver_l=${ver_l}"
    echo ""
    echo "https=${https}"
    echo "tarball=${tarball}"
    echo ""
    echo ""
fi


#  Do the main work -----------------------------------------------------------
if [[ "$(pwd)" != "${HOME}" ]]; then
    cd "${HOME}" || echo "Error: Failed to change directory: '${HOME}'." >&2
fi

#  Download the tarball
curl -L -O "${https}/${tarball}"
```

Unpack the file in the `HOME` directory.

```bash
#  Unpack the tarball
tar -xzf "${tarball}"
```

Add Julia to the `PATH` variable by appending the following line to the shell configuration file (e.g., `.bashrc`, `.bash_profile`, or `.zshrc`), and then source the file:

```bash
#  Define variable for shell configuration file
config="${HOME}/.bashrc"  # "${HOME}/.zshrc"  ## WARNING: Change as needed ##

#  Add Julia to PATH and apply changes
# shellcheck disable=SC1090,SC2016
{
    echo 'export PATH="${PATH}:${HOME}/julia-1.8.5/bin"' >> "${config}"
    source "${config}"
}
```

Next, clone the [Atria repository](https://github.com/cihga39871/Atria), activate [`env_protocol`](b-clone-the-protocol-repository-and-install-the-project-environment) (which contains its dependencies), and build Atria using Julia:

```bash
#  Navigate to the repository directory
cd "${HOME}/repos" || echo "Error: Failed to change directory: '${repo}'." >&2

#  Clone Atria repository from GitHub
if [[ ! -d Atria ]]; then
    git clone "https://github.com/cihga39871/Atria.git"
fi

#  Change into newly cloned Atria directory
cd Atria || echo "Error: Failed to change directory: 'Atria'." >&2

#  Activate environment containing Atria dependencies 
mamba activate env_protocol

#  Build Atria using the provided Julia build script
julia build_atria.jl
```

Locate the Atria binary and add its path to the shell configuration file. It will be in a directory like this: `atria-version/bin`, e.g., `atria-4.0.3/bin`. For example,

```bash
#  Add Atria to PATH and apply changes  ## WARNING: Change version as needed ##
# shellcheck disable=SC1090,SC2016
{
    echo 'export PATH="${PATH}:${HOME}/repos/Atria/atria-version/bin"' \
        >> "${config}"
    source "${config}"
}
```

Ensure the `env_protocol` environment is active when running Atria.
</details>
<br />

<a id="d-install-and-configure-integrative-genomics-viewer-igv"></a>
### D. Install and configure Integrative Genomics Viewer (IGV).
<details>
<summary><i>Text: Install and configure Integrative Genomics Viewer (IGV).</i></summary>
<br />

Integrative Genomics Viewer (IGV) is a graphical tool for the interactive exploration of ChIP-seq and other genomic data (Robinson et al., 2011, 2017, 2023; Thorvaldsdóttir et al., 2013). To install IGV, visit the [IGV download page](https://igv.org/doc/desktop/#DownloadPage/), select the appropriate bundle for the OS (e.g., "With Java Included"), unzip the file, and move the application to a preferred directory.
</details>
<br />
<br />

<a id="data-analysis"></a>
## Data analysis
<a id="a-prepare-and-concatenate-fasta-and-gff3-files-for-model-and-spike-in-organisms"></a>
### A. Prepare and concatenate FASTA and GFF3 files for model and spike-in organisms.
<a id="text"></a>
#### Text
<details>
<summary><i>Text: Prepare and concatenate FASTA and GFF3 files for model and spike-in organisms.</i></summary>
<br />

For detailed instructions and an example implementation of preparing and concatenating FASTA and GFF3 files for both the model and spike-in organisms, refer to the supplementary notebook [`download_process_fasta_gff3.md`](./download_process_fasta_gff3.md), which provides an annotated, step-by-step implementation of the following:
1. Download FASTA and GFF3 files from the [Saccharomyces Genome Database](http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases) (*S. cerevisiae*) and [Pombase](https://www.pombase.org/data/releases/) (*S. pombe*).
2. Process the files by standardizing chromosome names and removing incompatible formatting.
    + Chromosome names in *S. cerevisiae* and *S. pombe* FASTA files are standardized and simplified, with *S. pombe* chromosome names prefixed with "SP_" to enable the downstream separation of *S. pombe* alignments from *S. cerevisiae* alignments.
    + Chromosome names in GFF3 files are standardized in the same way. In the *S. cerevisiae* GFF3 file, gene and autonomously replicating sequence (ARS) `Name` fields are reassigned from their systematic names (e.g., "YEL021W") to their standard, more interpretable names (e.g., "URA3"). Additionally, URL-encoded characters and HTML entities are converted to readable equivalents. These steps are unnecessary for the *S. pombe* file, as its `Name` field uses standard names, and the file contains no special characters.
3. Concatenate the processed files for alignment and visualization.
</details>
<br />

<a id="b-generate-bowtie-2-index-files-from-the-concatenated-fasta-file"></a>
### B. Generate Bowtie 2 index files from the concatenated FASTA file.
<a id="text-1"></a>
#### Text
<details>
<summary><i>Text: Generate Bowtie 2 index files from the concatenated FASTA file.</i></summary>
<br />

To align ChIP-seq reads against both *S. cerevisiae* and *S. pombe*, e.g., for downstream spike-in normalization, Bowtie 2 index files are generated from the concatenated FASTA file. Additionally, a separate code chunk is provided for generating Bowtie 2 index files using only the *S. cerevisiae* FASTA file, for cases where spike-in normalization is not required.

**Steps overview:**
1. *Set up directories and file paths* for inputs and outputs.
2. *Activate the required environment* and check for dependencies.
3. *Run Bowtie 2 index generation* using `bowtie2-build`, logging output.
4. *Optional cleanup:* Compress logs and remove temporary files to save space.
</details>
<br />

<a id="bash-code-1"></a>
#### Bash code #1
<details>
<summary><i>Bash code #1: Generate Bowtie 2 index files from the concatenated FASTA file.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
debug=true

#  Set base repository paths
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"

#  Define data directories
dir_dat="${dir_rep}/data"
dir_gen="${dir_dat}/genomes"
dir_cat="${dir_gen}/concat"
dir_fas="${dir_cat}/fasta/proc"
dir_idx="${dir_cat}/index/bowtie2"
dir_log="${dir_idx}/logs"

#  Define data file
fil_fas="${dir_fas}/sc_sp_proc.fasta"

#  Set execution parameters
env_nam="env_protocol"
day="$(date '+%Y-%m%d')"
log_exc="${dir_log}/${day}.execute_bowtie2_build"

if ${debug:-false}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "debug=${debug}"
    echo ""
    echo "#  Set base repository paths"
    echo "dir_bas=${dir_bas}"
    echo "dir_rep=${dir_rep}"
    echo "dir_scr=${dir_scr}"
    echo "dir_fnc=${dir_fnc}"
    echo ""
    echo "#  Define data directories"
    echo "dir_dat=${dir_dat}"
    echo "dir_gen=${dir_gen}"
    echo "dir_cat=${dir_cat}"
    echo "dir_fas=${dir_fas}"
    echo "dir_idx=${dir_idx}"
    echo "dir_log=${dir_log}"
    echo ""
    echo "#  Define data file"
    echo "fil_fas=${fil_fas}"
    echo ""
    echo "#  Set execution parameters"
    echo "env_nam=${env_nam}"
    echo "day=${day}"
    echo "log_exc=${log_exc}"
    echo ""
    echo ""
fi


#  Create required directories if necessary -----------------------------------
# shellcheck disable=SC2086
mkdir -p ${dir_idx}/logs

if ${debug:-false}; then
    echo "#################################################"
    echo "## In- and output directory paths and contents ##"
    echo "#################################################"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_fas %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_fas}"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_idx %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_idx}"
    echo ""
    ls -lhaFG "${dir_idx}/logs"
    echo ""
    echo ""
fi


#  Activate the environment and check dependencies ----------------------------
#  Source utility functions
# shellcheck disable=SC1091
{
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/handle_env.sh"
}

#  Activate the required environment
handle_env "${env_nam}"

#  Ensure access to bowtie2-build
check_program_path "bowtie2-build"


#  Do the main work -----------------------------------------------------------
#  Generate Bowtie 2 index files for the concatenated genome
#  If necessary, decompress the FASTA file
if [[ ! -f "${fil_fas}" && -f "${fil_fas}.gz" ]]; then
    gunzip -c "${fil_fas}.gz" > "${fil_fas}"
fi

#  "Build" the Bowtie 2 index files using the decompressed FASTA file
if ${debug:-false}; then
    echo "###########################"
    echo "## Call to bowtie2-build ##"
    echo "###########################"
    echo ""
    echo "bowtie2-build \\"
    echo "    ${fil_fas} \\"
    echo "    ${dir_idx}/$(basename ${fil_fas} ".fasta") \\"
    echo "         >> >(tee -a ${log_exc}.stdout.txt) \\"
    echo "        2>> >(tee -a ${log_exc}.stderr.txt)"
    echo ""
    echo ""
fi

bowtie2-build \
    "${fil_fas}" \
    "${dir_idx}/$(basename ${fil_fas} ".fasta")" \
         > >(tee -a "${log_exc}.stdout.txt") \
        2> >(tee -a "${log_exc}.stderr.txt")


#  Cleanup: Compress logs and remove empty files ------------------------------
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_idx}/logs"
# ls -lhaFG "${dir_idx}/logs"  ## Uncomment to check log directory ##
```
</details>
<br />

<a id="bash-code-2"></a>
#### Bash code #2
<details>
<summary><i>Bash code #2: Generate Bowtie 2 index files from the </i>S. cerevisiae<i> FASTA file.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
debug=true

#  Set base repository paths
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"

#  Define data directories
dir_dat="${dir_rep}/data"
dir_gen="${dir_dat}/genomes"
dir_cer="${dir_gen}/cerevisiae"
dir_fas="${dir_cer}/fasta/proc"
dir_idx="${dir_cer}/index/bowtie2"
dir_log="${dir_idx}/logs"

#  Define data file
fil_fas="${dir_fas}/S288C_R64-5-1_proc.fasta"

#  Set execution parameters
env_nam="env_protocol"
day="$(date '+%Y-%m%d')"
log_exc="${dir_log}/${day}.execute_bowtie2_build"

if ${debug:-false}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "debug=${debug}"
    echo ""
    echo "#  Set base repository paths"
    echo "dir_bas=${dir_bas}"
    echo "dir_rep=${dir_rep}"
    echo "dir_scr=${dir_scr}"
    echo "dir_fnc=${dir_fnc}"
    echo ""
    echo "#  Define data directories"
    echo "dir_dat=${dir_dat}"
    echo "dir_gen=${dir_gen}"
    echo "dir_cer=${dir_cer}"
    echo "dir_fas=${dir_fas}"
    echo "dir_idx=${dir_idx}"
    echo "dir_log=${dir_log}"
    echo ""
    echo "#  Define data file"
    echo "fil_fas=${fil_fas}"
    echo ""
    echo "#  Set execution parameters"
    echo "env_nam=${env_nam}"
    echo "day=${day}"
    echo "log_exc=${log_exc}"
    echo ""
    echo ""
fi


#  Create required directories if necessary -----------------------------------
# shellcheck disable=SC2086
mkdir -p ${dir_idx}/logs

if ${debug:-false}; then
    echo "#################################################"
    echo "## In- and output directory paths and contents ##"
    echo "#################################################"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_fas %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_fas}"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_idx %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_idx}"
    echo ""
    ls -lhaFG "${dir_idx}/logs"
    echo ""
    echo ""
fi


#  Activate the environment and check dependencies ----------------------------
#  Source utility functions
# shellcheck disable=SC1091
{
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/handle_env.sh"
}

#  Activate the required environment
handle_env "${env_nam}"

#  Ensure access to bowtie2-build
check_program_path "bowtie2-build"


#  Generate Bowtie 2 index files for the concatenated genome ------------------
#  If necessary, decompress the FASTA file
if [[ ! -f "${fil_fas}" && -f "${fil_fas}.gz" ]]; then
    gunzip -c "${fil_fas}.gz" > "${fil_fas}"
fi

#  "Build" the Bowtie 2 index using the decompressed FASTA file
if ${debug:-false}; then
    echo "###########################"
    echo "## Call to bowtie2-build ##"
    echo "###########################"
    echo ""
    echo "bowtie2-build \\"
    echo "    ${fil_fas} \\"
    echo "    ${dir_idx}/$(basename ${fil_fas} ".fasta") \\"
    echo "         >> >(tee -a ${log_exc}.stdout.txt) \\"
    echo "        2>> >(tee -a ${log_exc}.stderr.txt)"
    echo ""
    echo ""
fi

bowtie2-build \
    "${fil_fas}" \
    "${dir_idx}/$(basename ${fil_fas} ".fasta")" \
         > >(tee -a "${log_exc}.stdout.txt") \
        2> >(tee -a "${log_exc}.stderr.txt")

#  Optional: Once the index is built, delete the decompressed FASTA file
if [[ -f "${fil_fas}" ]]; then rm "${fil_fas}"; fi


#  Cleanup: Compress logs and remove empty files ------------------------------
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_idx}/logs"
# ls -lhaFG "${dir_idx}/logs"  ## Uncomment to check log directory ##
```
</details>
<br />

<a id="c-obtain-and-organize-chip-seq-fastq-files"></a>
### C. Obtain and organize ChIP-seq FASTQ files.
<a id="text-2"></a>
#### Text
<details>
<summary><i>Text: Obtain and organize ChIP-seq FASTQ files.</i></summary>
<br />

`#TODO`
</details>
<br />

<a id="bash-code"></a>
#### Bash code
<details>
<summary><i>Bash code: Obtain and organize ChIP-seq FASTQ files.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
debug=true

#  Set base repository paths
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"

#  Define data directories
dir_dat="${dir_rep}/data"
dir_raw="${dir_dat}/raw"
dir_doc="${dir_raw}/docs"
dir_log="${dir_raw}/logs"
dir_pro="${dir_dat}/processed"
dir_sym="${dir_dat}/symlinked"

#  Define data files
fil_tbl="datasets_dropbox.tsv"  ## WARNING: Change as needed ##
pth_tsv="${dir_doc}/${fil_tbl}"

#  Set execution parameters
env_nam="env_protocol"
threads=6
slurm=true
time="6:00:00"
day="$(date '+%Y-%m%d')"

if ${debug:-false}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "debug=${debug}"
    echo ""
    echo "#  Set base repository paths"
    echo "dir_bas=${dir_bas}"
    echo "dir_rep=${dir_rep}"
    echo "dir_scr=${dir_scr}"
    echo "dir_fnc=${dir_fnc}"
    echo ""
    echo "#  Define data directories"
    echo "dir_dat=${dir_dat}"
    echo "dir_raw=${dir_raw}"
    echo "dir_doc=${dir_doc}"
    echo "dir_log=${dir_log}"
    echo "dir_pro=${dir_pro}"
    echo "dir_sym=${dir_sym}"
    echo ""
    echo "#  Define data files"
    echo "fil_tbl=${fil_tbl}"
    echo "pth_tsv=${pth_tsv}"
    echo ""
    echo "#  Set execution parameters"
    echo "env_nam=${env_nam}"
    echo "threads=${threads}"
    echo "slurm=${slurm}"
    echo "time=${time}"
    echo "day=${day}"
    echo ""
    echo ""
fi


#  Activate the environment and check dependencies ----------------------------
#  Source utility functions
# shellcheck disable=SC1091
{
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/echo_warning.sh"
    source "${dir_fnc}/handle_env.sh"
}

#  Activate the required environment
handle_env "${env_nam}"

#  Ensure access to necessary dependencies such as GNU Parallel, SLURM sbatch
check_program_path parallel
check_program_path wget

if ${slurm:-false}; then
    check_program_path sbatch &> /dev/null ||
        echo_warning \
            "SLURM is not available on this system. Do not use the '--slurm'" \
            "flag with the driver script."
fi



#  Create required directories if necessary -----------------------------------
#  If needed, create directory structure for raw and symlinked FASTQ files and
#+ logs
mkdir -p {${dir_raw},${dir_sym}}/logs

if ${debug:-false}; then
    echo "#########################################"
    echo "## Output directory paths and contents ##"
    echo "#########################################"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_raw %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_raw}"
    echo ""
    ls -lhaFG "${dir_raw}/logs"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_sym %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_sym}"
    echo ""
    ls -lhaFG "${dir_sym}/logs"
    echo ""
    echo ""
fi


#  Do the main work -----------------------------------------------------------
#  Download and symlink FASTQ files
if ${debug:-false}; then
    echo "###########################"
    echo "## Call to driver script ##"
    echo "###########################"
    echo ""
    echo "bash ${dir_scr}/execute_download_fastqs.sh \\"
    echo "    --threads ${threads} \\"
    echo "    --infile ${pth_tsv} \\"
    echo "    --dir_out ${dir_raw} \\"
    echo "    --dir_sym ${dir_sym} \\"
    echo "    --err_out ${dir_raw}/logs \\"
    echo "$(if ${slurm:-false}; then echo "$(if ${slurm:-false}; then echo "    --slurm \\"; fi)"; fi)"
    echo "    --time ${time} \\"
    echo "         > >(tee -a ${dir_raw}/logs/${day}.execute.stdout.txt) \\"
    echo "        2> >(tee -a ${dir_raw}/logs/${day}.execute.stderr.txt)"
fi

bash "${dir_scr}/execute_download_fastqs.sh" \
    --threads "${threads}" \
    --infile "${pth_tsv}" \
    --dir_out "${dir_raw}" \
    --dir_sym "${dir_sym}" \
    --err_out "${dir_raw}/logs" \
    $(if ${slurm:-false}; then echo "--slurm"; fi) \
    --time "${time}" \
         > >(tee -a "${dir_raw}/logs/${day}.execute.stdout.txt") \
        2> >(tee -a "${dir_raw}/logs/${day}.execute.stderr.txt")


#  Cleanup: Compress logs and remove empty files ------------------------------
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_raw}/logs"
# ls -lhaFG "${dir_raw}/logs"  ## Uncomment to check log directory ##
```
</details>
<br />

<a id="d-use-atria-to-perform-adapter-and-quality-trimming-of-sequenced-reads"></a>
### D. Use Atria to perform adapter and quality trimming of sequenced reads.
<details>
<summary><i>Bash code: Use Atria to perform adapter and quality trimming of sequenced reads.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
debug=true

#  Set base repository paths
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"

#  Define data directories
dir_dat="${dir_rep}/data"
dir_sym="${dir_dat}/symlinked"
dir_pro="${dir_dat}/processed"
dir_trm="${dir_pro}/trim_fastqs"

#  Set execution parameters
env_nam="env_protocol"
threads=6
max_job=2
slurm=false
day="$(date '+%Y-%m%d')"

#  Locate input files
infiles="$(  ## WARNING: Change search parameters as needed ##
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_sym}" \
        --pattern "*.fastq.gz" \
        --include "*Q_Hmo1*" \
        --depth 1 \
        --follow \
        --fastqs
)"
# --include "*H?o1*" \

if ${debug:-false}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "debug=${debug}"
    echo ""
    echo "#  Set base repository paths"
    echo "dir_bas=${dir_bas}"
    echo "dir_rep=${dir_rep}"
    echo "dir_scr=${dir_scr}"
    echo "dir_fnc=${dir_fnc}"
    echo ""
    echo "#  Define data directories"
    echo "dir_dat=${dir_dat}"
    echo "dir_sym=${dir_sym}"
    echo "dir_pro=${dir_pro}"
    echo "dir_trm=${dir_trm}"
    echo ""
    echo "#  Set execution parameters"
    echo "env_nam=${env_nam}"
    echo "threads=${threads}"
    echo "max_job=${max_job}"
    echo "slurm=${slurm}"
    echo "day=${day}"
    echo ""
    echo "#  Locate input files"
    echo "infiles=${infiles}"
    echo ""
    echo ""
fi


#  Activate the environment and check dependencies ----------------------------
#  Source utility functions
# shellcheck disable=SC1091
{
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/echo_warning.sh"
    source "${dir_fnc}/handle_env.sh"
}

#  Activate the required environment
handle_env "${env_nam}"

#  Check availability of Atria and other necessary tools
check_program_path atria
check_program_path pbzip2
check_program_path pigz

if ${slurm:-false}; then
    check_program_path sbatch &> /dev/null ||
        echo_warning \
            "SLURM is not available on this system. Do not use the '--slurm'" \
            "flag with the driver script."
elif [[ "${threads}" -gt 1 ]]; then
    check_program_path parallel
fi


#  Create required directories if necessary -----------------------------------
# shellcheck disable=SC2086
mkdir -p ${dir_trm}/logs

if ${debug:-false}; then
    echo "#################################################"
    echo "## In- and output directory paths and contents ##"
    echo "#################################################"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_sym %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_sym}"
    echo ""
    ls -lhaFG "${dir_sym}/logs"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_trm %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_trm}"
    echo ""
    ls -lhaFG "${dir_trm}/logs"
    echo ""
    echo ""
fi


#  Do the main work -----------------------------------------------------------
#  Run the driver script to trim FASTQ files with Atria
if ${debug:-false}; then
    echo "###########################"
    echo "## Call to driver script ##"
    echo "###########################"
    echo ""
    echo "bash ${dir_scr}/execute_trim_fastqs.sh \\"
    echo "    --verbose \\"
    echo "    --threads ${threads} \\"
    echo "    --infiles ${infiles} \\"
    echo "    --dir_out ${dir_trm} \\"
    echo "    --err_out ${dir_trm}/logs \\"
    echo "    --max_job ${max_job} \\"
    echo "$(if ${slurm:-false}; then echo "    --slurm \\"; fi)"
    echo "         > >(tee -a ${dir_trm}/logs/${day}.execute.stdout.txt) \\"
    echo "        2> >(tee -a ${dir_trm}/logs/${day}.execute.stderr.txt)"
fi

bash "${dir_scr}/execute_trim_fastqs.sh" \
    --verbose \
    --threads "${threads}" \
    --infiles "${infiles}" \
    --dir_out "${dir_trm}" \
    --err_out "${dir_trm}/logs" \
    --max_job ${max_job} \
    $(if ${slurm:-false}; then echo "--slurm"; fi) \
         >> >(tee -a "${dir_trm}/logs/${day}.execute.stdout.txt") \
        2>> >(tee -a "${dir_trm}/logs/${day}.execute.stderr.txt")


#  Cleanup: Compress logs and remove empty files ------------------------------
#  Move Atria LOG and JSON files to the logs directory
mv ${dir_trm}/*.{log,json} "${dir_trm}/logs"

#  Compress large stdout, stderr, LOG, and JSON files, and remove files with
#+ size 0
bash "${dir_scr}/compress_remove_files.sh" \
    --dir_fnd "${dir_trm}/logs"

bash "${dir_scr}/compress_remove_files.sh" \
    --dir_fnd "${dir_trm}/logs" \
    --pattern "*.log"

bash "${dir_scr}/compress_remove_files.sh" \
    --dir_fnd "${dir_trm}/logs" \
    --pattern "*.json"

# ls -lhaFG "${dir_trm}/logs"  ## Uncomment to check log directory ##
```
</details>
<br />

<a id="e-align-sequenced-reads-with-bowtie-2-and-process-the-read-alignments"></a>
### E. Align sequenced reads with Bowtie 2 and process the read alignments.
<details>
<summary><i>Bash code: Align sequenced reads with Bowtie 2 and process the read alignments.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
debug=true

#  Set base repository paths
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"

#  Define data directories
dir_dat="${dir_rep}/data"
dir_idx="${dir_dat}/genomes/concat/index"
dir_pro="${dir_dat}/processed"
dir_trm="${dir_pro}/trim_fastqs"

#  Set execution parameters and related variables
env_nam="env_protocol"
threads=6
aligner="bowtie2"
a_type="global"
req_flg=true
flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
mapq=1
nam_job="align_fastqs"
max_job=2
slurm=false
time="1:00:00"
day="$(date '+%Y-%m%d')"

#  Locate input files
infiles="$(  ## WARNING: Change search parameters as needed ##
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_trm}" \
        --pattern "*.atria.fastq.gz" \
        --depth 1 \
        --fastqs
)"

#  Set index path
str_idx="sc_sp_proc"
if [[ ${aligner} == "bwa" ]]; then
    stm_idx="${dir_idx}/${aligner}/${str_idx}.fa"
else
    stm_idx="${dir_idx}/${aligner}/${str_idx}"
fi

#  Set output directories
dir_aln="${dir_pro}/align_fastqs"
dir_out="${dir_aln}/${aligner}_${a_type}_flag-${flg}_mapq-${mapq}"

if ${debug:-false}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "debug=${debug}"
    echo ""
    echo "#  Set base repository paths"
    echo "dir_bas=${dir_bas}"
    echo "dir_rep=${dir_rep}"
    echo "dir_scr=${dir_scr}"
    echo "dir_fnc=${dir_fnc}"
    echo ""
    echo "#  Define data directories"
    echo "dir_dat=${dir_dat}"
    echo "dir_idx=${dir_idx}"
    echo "dir_pro=${dir_pro}"
    echo "dir_trm=${dir_trm}"
    echo ""
    echo "#  Set execution parameters and related variables"
    echo "env_nam=${env_nam}"
    echo "threads=${threads}"
    echo "aligner=${aligner}"
    echo "a_type=${a_type}"
    echo "req_flg=${req_flg}"
    echo "flg=${flg}"
    echo "mapq=${mapq}"
    echo "nam_job=${nam_job}"
    echo "max_job=${max_job}"
    echo "slurm=${slurm}"
    echo "time=${time}"
    echo "day=${day}"
    echo ""
    echo "#  Locate input files"
    echo "infiles=${infiles}"
    echo ""
    echo "#  Set index path"
    echo "str_idx=${str_idx}"
    echo "stm_idx=${stm_idx}"
    echo ""
    echo "#  Set output directories"
    echo "dir_aln=${dir_aln}"
    echo "dir_out=${dir_out}"
    echo ""
    echo ""
fi


#  Activate the environment and check dependencies ----------------------------
#  Source utility functions
# shellcheck disable=SC1091
{
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/echo_warning.sh"
    source "${dir_fnc}/handle_env.sh"
}

#  Activate the required environment
handle_env "${env_nam}"

#  Check the availability of Bowtie 2, Samtools, and SLURM sbatch
check_program_path bowtie2
check_program_path samtools

if ${slurm:-false}; then
    check_program_path sbatch &> /dev/null ||
        echo_warning \
            "SLURM is not available on this system. Do not use the '--slurm'" \
            "flag with the driver script."
elif [[ "${threads}" -gt 1 ]]; then
    check_program_path parallel
fi


#  Create required directories if necessary -----------------------------------
#  Create output directory structure for trimmed FASTQ files and logs
mkdir -p ${dir_out}/{init,sc,sp}/logs

if ${debug:-false}; then
    echo "#################################################"
    echo "## In- and output directory paths and contents ##"
    echo "#################################################"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_trm %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_trm}"
    echo ""
    ls -lhaFG "${dir_trm}/logs"
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%%"
    echo "%% ${dir_out}/init %%"
    echo "%%%%%%%%%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_out}/init"
    echo ""
    ls -lhaFG "${dir_out}/init/logs"
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%"
    echo "%% ${dir_out}/sc %%"
    echo "%%%%%%%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_out}/sc"
    echo ""
    ls -lhaFG "${dir_out}/sc/logs"
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%"
    echo "%% ${dir_out}/sp %%"
    echo "%%%%%%%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_out}/sp"
    echo ""
    ls -lhaFG "${dir_out}/sp/logs"
    echo ""
    echo ""
fi


#  Do the main work -----------------------------------------------------------
#  Run the driver script to align and post-process FASTQ files
if ${debug:-false}; then
    echo "###########################"
    echo "## Call to driver script ##"
    echo "###########################"
    echo ""
    echo "bash "${dir_scr}/execute_align_fastqs.sh" \\"
    echo "    --verbose \\"
    echo "    --threads ${threads} \\"
    echo "    --aligner ${aligner} \\"
    echo "    --a_type ${a_type} \\"
    echo "    --mapq ${mapq} \\"
    echo "    $(if ${req_flg:-false}; then echo --req_flg; fi) \\"
    echo "    --index ${stm_idx} \\"
    echo "    --infiles ${infiles} \\"
    echo "    --dir_out ${dir_out}/init \\"
    echo "    --err_out ${dir_out}/init/logs \\"
    echo "    --max_job ${max_job} \\"
    echo "    $(if ${slurm:-false}; then echo --slurm; fi) \\"
    echo "         >> >(tee -a ${dir_out}/init/logs/${day}.execute.stdout.txt) \\"
    echo "        2>> >(tee -a ${dir_out}/init/logs/${day}.execute.stderr.txt)"
fi

bash "${dir_scr}/execute_align_fastqs.sh" \
    --verbose \
    --threads "${threads}" \
    --aligner "${aligner}" \
    --a_type "${a_type}" \
    --mapq "${mapq}" \
    $(if ${req_flg:-false}; then echo "--req_flg"; fi) \
    --index "${stm_idx}" \
    --infiles "${infiles}" \
    --dir_out "${dir_out}/init" \
    --err_out "${dir_out}/init/logs" \
    --max_job "${max_job}" \
    $(if ${slurm:-false}; then echo "--slurm"; fi) \
         >> >(tee -a "${dir_out}/init/logs/${day}.execute.stdout.txt") \
        2>> >(tee -a "${dir_out}/init/logs/${day}.execute.stderr.txt")

#  Adjust variable 'infiles' assignment for filtering BAM files
infiles="$(  ## WARNING: Change search parameters as needed ##
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_out}/init" \
        --pattern "*.bam" \
        --depth 1
)"

#  Run the driver script to filter BAM files for S. cerevisiae ("sc")
#+ alignments (i.e., the "main" alignments)
if ${debug:-false}; then
    echo "###########################"
    echo "## Call to driver script ##"
    echo "###########################"
    echo ""
    echo "bash ${dir_scr}/execute_filter_bams.sh \\"
    echo "    --verbose \\"
    echo "    --threads ${threads} \\"
    echo "    --infiles ${infiles} \\"
    echo "    --dir_out ${dir_out}/sc \\"
    echo "    --retain sc \\"
    echo "    --err_out ${dir_out}/sc/logs \\"
    echo "    --max_job ${max_job} \\"
    echo "    $(if ${slurm:-false}; then echo --slurm; fi) \\"
    echo "         >> >(tee -a ${dir_out}/sc/logs/${day}.execute.stdout.txt) \\"
    echo "        2>> >(tee -a ${dir_out}/sc/logs/${day}.execute.stderr.txt)"
fi

bash "${dir_scr}/execute_filter_bams.sh" \
    --verbose \
    --threads "${threads}" \
    --infiles "${infiles}" \
    --dir_out "${dir_out}/sc" \
    --retain "sc" \
    --err_out "${dir_out}/sc/logs" \
    --max_job "${max_job}" \
    $(if ${slurm:-false}; then echo "--slurm"; fi) \
         >> >(tee -a "${dir_out}/sc/logs/${day}.execute.stdout.txt") \
        2>> >(tee -a "${dir_out}/sc/logs/${day}.execute.stderr.txt")

#  Run the driver script to filter BAM files for S. pombe ("sp") alignments
#+ (i.e., the "spike-in" alignments)
if ${debug:-false}; then
    echo "###########################"
    echo "## Call to driver script ##"
    echo "###########################"
    echo ""
    echo "bash ${dir_scr}/execute_filter_bams.sh \\"
    echo "    --verbose \\"
    echo "    --threads ${threads} \\"
    echo "    --infiles ${infiles} \\"
    echo "    --dir_out ${dir_out}/sp \\"
    echo "    --retain sp \\"
    echo "    --err_out ${dir_out}/sp/logs \\"
    echo "    --max_job ${max_job} \\"
    echo "    $(if ${slurm:-false}; then echo --slurm; fi) \\"
    echo "         >> >(tee -a ${dir_out}/sp/logs/${day}.execute.stdout.txt) \\"
    echo "        2>> >(tee -a ${dir_out}/sp/logs/${day}.execute.stderr.txt)"
fi

bash "${dir_scr}/execute_filter_bams.sh" \
    --verbose \
    --threads "${threads}" \
    --infiles "${infiles}" \
    --dir_out "${dir_out}/sp" \
    --retain "sp" \
    --err_out "${dir_out}/sp/logs" \
    --max_job "${max_job}" \
    $(if ${slurm:-false}; then echo "--slurm"; fi) \
         >> >(tee -a "${dir_out}/sp/logs/${day}.execute.stdout.txt") \
        2>> >(tee -a "${dir_out}/sp/logs/${day}.execute.stderr.txt")


#  Cleanup: Compress logs and remove empty files ------------------------------
#  Compress large stdout and stderr files, and remove files with size 0
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_out}/init/logs"
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_out}/sc/logs"
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_out}/sp/logs"

# ls -lhaFG "${dir_out}/init/logs"  ## Uncomment to check log directory ##
# ls -lhaFG "${dir_out}/sc/logs"    ## Uncomment to check log directory ##
# ls -lhaFG "${dir_out}/sp/logs"    ## Uncomment to check log directory ##
```
</details>
<br />

<a id="f-compute-normalized-coverage"></a>
### F. Compute normalized coverage.
<a id="text-3"></a>
#### Text
<details>
<summary><i>Text: Compute normalized coverage.</i></summary>
<br />

**Overview**  
This section demonstrates how to compute ChIP-seq normalized coverage ([Dickson et al., *Sci Rep*, 2023](https://www.nature.com/articles/s41598-023-34430-2))&mdash;i.e., binned coverage adjusted for fragment length that sums to unity.

In the corresponding [Bash code chunk](#bash-code), the signal type is specified by setting `typ_sig` to `"norm"` (normalized coverage). Other options are `"unadj"` (unadjusted binned coverage) and `"frag"` (fragment length-adjusted binned coverage; e.g., see [Dickson et al., *Sci Rep*, 2023](https://www.nature.com/articles/s41598-023-34430-2)). For details, see below and the [`compute_signal.py` documentation](https://github.com/kalavattam/protocol_chipseq_signal_norm/blob/main/scripts/compute_signal.py). BEDGRAPH and log output files are organized by signal type.

---

**Steps performed in code chunk**  
1. *Define paths and environment:* Define input BAM paths, output directories, and dependencies.
2. *Activate environment and check dependencies:* Ensure required software (GNU Parallel, Python, etc.) is available.
3. *Compute signal:* Call `execute_compute_signal.sh` to generate normalized coverage BEDGRAPH files.
4. *Verify normalization:* Optionally confirm the normalized coverage sums to ~1.
5. *Cleanup:* Compress logs and remove empty files.

---

**Output files**  
1. *BEDGRAPH coverage tracks* (`*.bdg.gz`) organized by signal type.
2. *Log files* tracking execution steps, stored separately.
</details>
<br />

<a id="bash-code-1"></a>
#### Bash code
<details>
<summary><i>Bash code: Compute normalized coverage.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
debug=true

#  Set base repository paths
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"

#  Set alignment parameters
aligner="bowtie2"
a_type="global"
req_flg=true
flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
mapq=1
prm_aln="${aligner}_${a_type}_flag-${flg}_mapq-${mapq}"

#  Define data directories
dir_dat="${dir_rep}/data"
dir_pro="${dir_dat}/processed"
dir_aln="${dir_pro}/align_fastqs"
dir_det="${dir_aln}/${prm_aln}"
dir_bam="${dir_det}/sc"

#  Set output directories
dir_cvg="${dir_pro}/compute_signal/${prm_aln}"
dir_trk="${dir_cvg}/norm"
dir_log="${dir_trk}/logs"

#  Set execution parameters and related variables
env_nam="env_protocol"
scr_exc="${dir_scr}/execute_compute_signal.sh"
log_exc="${dir_log}/$(date '+%Y-%m%d').$(basename "${scr_exc}" ".sh")"

threads=6
infiles="$(  ## WARNING: Change search parameters as needed ##
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_bam}" \
        --pattern "*.bam"
)"
dir_out="${dir_trk}"
typ_out="bdg.gz"
siz_bin=30
typ_sig="norm"
err_out="${dir_log}"
nam_job="compute_signal_${typ_sig}"
max_job=2
slurm=false

if ${debug:-false}; then
    echo "##########################"
    echo "## Variable assignments ##"
    echo "##########################"
    echo ""
    echo "debug=${debug}"
    echo ""
    echo "#  Set base repository paths"
    echo "dir_bas=${dir_bas}"
    echo "dir_rep=${dir_rep}"
    echo "dir_scr=${dir_scr}"
    echo "dir_fnc=${dir_fnc}"
    echo ""
    echo "#  Set alignment parameters"
    echo "aligner=${aligner}"
    echo "a_type=${a_type}"
    echo "req_flg=${req_flg}"
    echo "flg=${flg}"
    echo "mapq=${mapq}"
    echo "prm_aln=${prm_aln}"
    echo ""
    echo "#  Define data directories"
    echo "dir_dat=${dir_dat}"
    echo "dir_pro=${dir_pro}"
    echo "dir_aln=${dir_aln}"
    echo "dir_det=${dir_det}"
    echo "dir_bam=${dir_bam}"
    echo ""
    echo "#  Set output directories"
    echo "dir_cvg=${dir_cvg}"
    echo "dir_trk=${dir_trk}"
    echo "dir_log=${dir_log}"
    echo ""
    echo "#  Set execution parameters and related variables"
    echo "env_nam=${env_nam}"
    echo "scr_exc=${scr_exc}"
    echo "log_exc=${log_exc}"
    echo ""
    echo "threads=${threads}"
    echo "infiles=${infiles}"
    echo "dir_out=${dir_out}"
    echo "typ_out=${typ_out}"
    echo "siz_bin=${siz_bin}"
    echo "typ_sig=${typ_sig}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo "max_job=${max_job}"
    echo "slurm=${slurm}"
    echo ""
    echo ""
fi


#  Create required directories if necessary -----------------------------------
mkdir -p "${err_out}"

if ${debug:-false}; then
    echo "#################################################"
    echo "## In- and output directory paths and contents ##"
    echo "#################################################"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_bam %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_bam}"
    echo ""
    ls -lhaFG "${dir_bam}/logs"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_out %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_out}"
    echo ""
    ls -lhaFG "${err_out}"
    echo ""
    echo ""
fi


#  Activate the environment and check dependencies ----------------------------
#  Source utility functions
# shellcheck disable=SC1091
{
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/check_unity.sh"
    source "${dir_fnc}/echo_warning.sh"
    source "${dir_fnc}/handle_env.sh"
}

#  Activate the required environment
handle_env "${env_nam}"

#  Check the availability of necessary dependencies such as GNU Parallel,
#+ Python, and SLURM sbatch
check_program_path awk
check_program_path python

if ${slurm:-false}; then
    check_program_path sbatch &> /dev/null ||
        echo_warning \
            "SLURM is not available on this system. Do not use the '--slurm'" \
            "flag with the driver script."
elif [[ "${threads}" -gt 1 ]]; then
    check_program_path parallel
fi


#  Compute normalized coverage ------------------------------------------------
if ${debug:-false}; then
    echo "###########################"
    echo "## Call to driver script ##"
    echo "###########################"
    echo ""
    echo "bash ${scr_exc} \\"
    echo "    --verbose \\"
    echo "    --threads ${threads} \\"
    echo "    --infiles ${infiles} \\"
    echo "    --dir_out ${dir_trk} \\"
    echo "    --typ_out ${typ_out} \\"
    echo "    --siz_bin ${siz_bin} \\"
    echo "    --typ_sig ${typ_sig} \\"
    echo "    --err_out ${err_out} \\"
    echo "    --nam_job ${nam_job} \\"
    echo "    --max_job ${max_job} \\"
    echo "$(if ${slurm:-false}; then echo "    --slurm \\"; fi)"
    echo "         >> >(tee -a ${log_exc}.stdout.txt) \\"
    echo "        2>> >(tee -a ${log_exc}.stderr.txt)"
    echo ""
    echo ""
fi

#  Run the driver script: Include '--dry_run' to perform a test run without
#+ executing commands
bash "${scr_exc}" \
    --verbose \
    --threads "${threads}" \
    --infiles "${infiles}" \
    --dir_out "${dir_trk}" \
    --typ_out "${typ_out}" \
    --siz_bin "${siz_bin}" \
    --typ_sig "${typ_sig}" \
    --err_out "${err_out}" \
    --nam_job "${nam_job}" \
    --max_job "${max_job}" \
    $(if ${slurm:-false}; then echo "--slurm"; fi) \
         >> >(tee -a "${log_exc}.stdout.txt") \
        2>> >(tee -a "${log_exc}.stderr.txt")

#  Check that each normalized coverage file sums to unity
if ${debug:-false} && [[ "${typ_sig}" == "norm" ]]; then
    echo "#######################################################"
    echo "## Check that normalized coverage files sum to unity ##"
    echo "#######################################################"
    echo ""
    for bdg in "${dir_trk}/"*".${typ_out}"; do
        echo -n "${bdg##*/}: "
        check_unity "${bdg}"
        echo ""
    done
    echo ""
fi


#  Cleanup: Compress logs and remove empty files ------------------------------
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${err_out}"
# ls -lhaFG "${err_out}"  ## Uncomment to check directory for logs ##
```
</details>
<br />

<a id="g-construct-sample-tables-recording-computed-scaling-factors-for-normalization"></a>
### G. Construct sample tables recording computed scaling factors for normalization.
<a id="1-construct-sample-table-recording-siq-chip-%24alpha%24-scaling-factors"></a>
#### 1. Construct sample table recording siQ-ChIP $\alpha$ scaling factors.
<a id="text-4"></a>
##### Text
<details>
<summary><i>Text: Construct sample table recording siQ-ChIP $\alpha$ scaling factors.</i></summary>
<br />

`#TODO`
</details>
<br />

<a id="bash-code-2"></a>
##### Bash code
<details>
<summary><i>Bash code: Construct sample table recording siQ-ChIP $\alpha$ scaling factors.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
debug=true

#  Set base repository paths
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"

#  Set alignment parameters
aligner="bowtie2"
a_type="global"
req_flg=true
flag="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
mapq=1
prm_aln="${aligner}_${a_type}_flag-${flag}_mapq-${mapq}"

#  Define data directories
dir_dat="${dir_rep}/data"
dir_pro="${dir_dat}/processed"
dir_aln="${dir_pro}/align_fastqs/${prm_aln}"
dir_bam="${dir_aln}/sc"
dir_cvg="${dir_pro}/compute_signal/${prm_aln}"

#  Set output directories
dir_tbl="${dir_cvg}/tables"
dir_log="${dir_tbl}/logs"

#  Set execution parameters and related variables
typ_sig="alpha"
env_nam="env_protocol"
scr_exc="${dir_scr}/execute_calculate_scaling_factor_${typ_sig}.sh"
log_exc="${dir_log}/$(date '+%Y-%m%d').$(basename "${scr_exc}" ".sh")"

threads=6
ser_ip="$(  ## WARNING: Change search parameters as needed ##
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_bam}" \
        --pattern "*.bam" \
        --include "IP*" \
        --exclude "*Brn1*"
)"
ser_in="$(sed 's:\/IP:\/in:g' <<< "${ser_ip}")"
tbl_met="${dir_dat}/raw/docs/measurements_siqchip.tsv"
eqn="6nd"
fil_out="${dir_tbl}/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_${typ_sig}.tsv"
err_out="${dir_log}"
max_job=2
slurm=false

if ${debug:-false}; then
    echo "##########################"
    echo "## Variable assignments ##"
    echo "##########################"
    echo ""
    echo "debug=${debug}"
    echo ""
    echo "#  Set base repository paths"
    echo "dir_bas=${dir_bas}"
    echo "dir_rep=${dir_rep}"
    echo "dir_scr=${dir_scr}"
    echo "dir_fnc=${dir_fnc}"
    echo ""
    echo "#  Set alignment parameters"
    echo "aligner=${aligner}"
    echo "a_type=${a_type}"
    echo "req_flg=${req_flg}"
    echo "flag=${flag}"
    echo "mapq=${mapq}"
    echo "prm_aln=${prm_aln}"
    echo ""
    echo "#  Define data directories"
    echo "dir_dat=${dir_dat}"
    echo "dir_pro=${dir_pro}"
    echo "dir_aln=${dir_aln}"
    echo "dir_bam=${dir_bam}"
    echo "dir_cvg=${dir_cvg}"
    echo ""
    echo "#  Set output directories"
    echo "dir_tbl=${dir_tbl}"
    echo "dir_log=${dir_log}"
    echo ""
    echo "#  Set execution parameters and related variables"
    echo "typ_sig=${typ_sig}"
    echo "env_nam=${env_nam}"
    echo "scr_exc=${scr_exc}"
    echo "log_exc=${log_exc}"
    echo ""
    echo "threads=${threads}"
    echo "ser_ip=${ser_ip}"
    echo "ser_in=${ser_in}"
    echo "tbl_met=${tbl_met}"
    echo "eqn=${eqn}"
    echo "fil_out=${fil_out}"
    echo "err_out=${err_out}"
    echo "max_job=${max_job}"
    echo "slurm=${slurm}"
    echo ""
    echo ""
fi


#  Create required directories if necessary -----------------------------------
mkdir -p "${err_out}"

if ${debug:-false}; then
    echo "#################################################"
    echo "## In- and output directory paths and contents ##"
    echo "#################################################"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_bam %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_bam}"
    echo ""
    ls -lhaFG "${dir_bam}/logs"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_tbl %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_tbl}"
    echo ""
    ls -lhaFG "${err_out}"
    echo ""
    echo ""
fi


#  Activate the environment and check dependencies ----------------------------
#  Source utility functions
# shellcheck disable=SC1091
{
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/echo_warning.sh"
    source "${dir_fnc}/handle_env.sh"
}

#  Activate the required environment
handle_env "${env_nam}"

#  Check the availability of necessary dependencies such as GNU Parallel,
#+ Python, Samtools, and SLURM
check_program_path awk
check_program_path python
check_program_path samtools

if ${slurm:-false}; then
    check_program_path sbatch &> /dev/null ||
        echo_warning \
            "SLURM is not available on this system. Do not use the '--slurm'" \
            "flag with the driver script."
elif [[ "${threads}" -gt 1 ]]; then
    check_program_path parallel
fi


#  Construct table ------------------------------------------------------------
if ${debug:-false}; then
    echo "###########################"
    echo "## Call to driver script ##"
    echo "###########################"
    echo ""
    echo "bash \"${scr_exc}\" \\"
    echo "    --verbose \\"
    echo "    --threads ${threads} \\"
    echo "    --ser_ip \"${ser_ip}\" \\"
    echo "    --ser_in \"${ser_in}\" \\"
    echo "    --tbl_met \"${tbl_met}\" \\"
    echo "    --eqn \"${eqn}\" \\"
    echo "    --fil_out \"${fil_out}\" \\"
    echo "    --err_out \"${err_out}\" \\"
    echo "    --slurm"
    echo "         >> >(tee -a ${log_exc}.stdout.txt) \\"
    echo "        2>> >(tee -a ${log_exc}.stderr.txt)"
    echo ""
    echo ""
fi

#  Run the driver script: Include '--dry_run' to perform a test run without
#+ executing commands
bash "${scr_exc}" \
    --verbose \
    --threads ${threads} \
    --ser_ip "${ser_ip}" \
    --ser_in "${ser_in}" \
    --tbl_met "${tbl_met}" \
    --eqn "${eqn}" \
    --fil_out "${fil_out}" \
    --err_out "${err_out}" \
    $(if ${slurm:-false}; then echo "--slurm"; fi) \
         >> >(tee -a "${log_exc}.stdout.txt") \
        2>> >(tee -a "${log_exc}.stderr.txt")

if [[ -f "${fil_out}" ]]; then
    #  Sort the table by rows
    awk 'NR == 1; NR > 1 && NF { print | "sort" }' "${fil_out}" \
        > "${dir_tbl}/tmp.tsv"

    #  Replace the original table with the sorted version
    mv -f "${dir_tbl}/tmp.tsv" "${fil_out}"
fi
# cat "${fil_out}"  ## Uncomment to check table ##


#  Cleanup: Compress logs and remove empty files ------------------------------
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${err_out}"
# ls -lhaFG "${err_out}"  ## Uncomment to check directory for logs ##
```
</details>
<br />

<a id="2-construct-sample-table-recording-spike-in-scaling-factors"></a>
#### 2. Construct sample table recording spike-in scaling factors.
<a id="text-5"></a>
##### Text
<details>
<summary><i>Text: Construct sample table recording spike-in scaling factors.</i></summary>
<br />

`#TODO`
</details>
<br />

<a id="bash-code-3"></a>
##### Bash code
<details>
<summary><i>Bash code: Construct sample table recording spike-in scaling factors.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
debug=true

#  Set hardcoded paths, values, etc.
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"
dir_dat="${dir_rep}/data"
dir_pro="${dir_dat}/processed"

aligner="bowtie2"
a_type="global"
req_flg=true
flag="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
mapq=1
details="${aligner}_${a_type}_flag-${flag}_mapq-${mapq}"

dir_aln="${dir_pro}/align_fastqs/${details}"
dir_bam="${dir_aln}/sc"
dir_cvg="${dir_pro}/compute_signal/${details}"
dir_tbl="${dir_cvg}/tables"
dir_log="${dir_tbl}/logs"
typ_sig="spike"

pattern="*.bam"
include="IP*"
exclude="*Brn1*"

env_nam="env_protocol"
scr_exc="${dir_scr}/execute_calculate_scaling_factor_${typ_sig}.sh"
day="$(date '+%Y-%m%d')"
log_exc="${dir_log}/${day}.$(basename "${scr_exc}" ".sh")"

#  Set hardcoded argument assignments
threads=8
ser_mip="$(  ## WARNING: Change search parameters as needed ##
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_bam}" \
        --pattern "${pattern}" \
        --include "${include}" \
        --exclude "${exclude}"
)"
ser_sip="$(sed 's:sc:sp:g' <<< "${ser_mip}")"
ser_min="$(sed 's:\/IP:\/in:g' <<< "${ser_mip}")"
ser_sin="$(sed 's:\/IP:\/in:g' <<< "${ser_sip}")"
fil_out="${dir_tbl}/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_${typ_sig}.tsv"
err_out="${dir_log}"

if ${debug:-false}; then
    echo "##########################"
    echo "## Variable assignments ##"
    echo "##########################"
    echo ""
    echo "debug=${debug}"
    echo ""
    echo "dir_bas=${dir_bas}"
    echo "dir_rep=${dir_rep}"
    echo "dir_scr=${dir_scr}"
    echo "dir_fnc=${dir_fnc}"
    echo "dir_dat=${dir_dat}"
    echo "dir_pro=${dir_pro}"
    echo ""
    echo "aligner=${aligner}"
    echo "a_type=${a_type}"
    echo "req_flg=${req_flg}"
    echo "flag=${flag}"
    echo "mapq=${mapq}"
    echo "details=${details}"
    echo ""
    echo "dir_aln=${dir_aln}"
    echo "dir_bam=${dir_bam}"
    echo "dir_cvg=${dir_cvg}"
    echo "dir_tbl=${dir_tbl}"
    echo "dir_log=${dir_log}"
    echo "typ_sig=${typ_sig}"
    echo ""
    echo "pattern=${pattern}"
    echo "include=${include}"
    echo "exclude=${exclude}"
    echo ""
    echo "env_nam=${env_nam}"
    echo "scr_exc=${scr_exc}"
    echo "day=${day}"
    echo "log_exc=${log_exc}"
    echo ""
    echo "threads=${threads}"
    echo "ser_mip=${ser_mip}"
    echo "ser_sip=${ser_sip}"
    echo "ser_min=${ser_min}"
    echo "ser_sin=${ser_sin}"
    echo "fil_out=${fil_out}"
    echo "err_out=${err_out}"
    echo ""
    echo ""
fi


#  Create required directories if necessary -----------------------------------
mkdir -p "${err_out}"

if ${debug:-false}; then
    echo "#################################################"
    echo "## In- and output directory paths and contents ##"
    echo "#################################################"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_bam %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_bam}"
    echo ""
    ls -lhaFG "${dir_bam}/logs"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_tbl %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_tbl}"
    echo ""
    ls -lhaFG "${err_out}"
    echo ""
    echo ""
fi


#  Activate the environment and check dependencies ----------------------------
#  Source utility functions
# shellcheck disable=SC1091
{
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/echo_warning.sh"
    source "${dir_fnc}/handle_env.sh"
}

#  Activate the required environment
handle_env "${env_nam}"

#  Check the availability of necessary dependencies such as GNU Parallel,
#+ Python, Samtools, and SLURM
check_program_path awk
check_program_path python
check_program_path samtools

if ${slurm:-false}; then
    check_program_path sbatch &> /dev/null ||
        echo_warning \
            "SLURM is not available on this system. Do not use the '--slurm'" \
            "flag with the driver script."
elif [[ "${threads}" -gt 1 ]]; then
    check_program_path parallel
fi


#  Construct table ------------------------------------------------------------
if ${debug:-false}; then
    echo "###########################"
    echo "## Call to driver script ##"
    echo "###########################"
    echo ""
    echo "bash \"${scr_exc}\" \\"
    echo "    --verbose \\"
    echo "    --threads ${threads} \\"
    echo "    --ser_mip \"${ser_mip}\" \\"
    echo "    --ser_sip \"${ser_sip}\" \\"
    echo "    --ser_min \"${ser_min}\" \\"
    echo "    --ser_sin \"${ser_sin}\" \\"
    echo "    --fil_out \"${fil_out}\" \\"
    echo "    --err_out \"${err_out}\" \\"
    echo "    --slurm"
    echo "         >> >(tee -a ${log_exc}.stdout.txt) \\"
    echo "        2>> >(tee -a ${log_exc}.stderr.txt)"
    echo ""
    echo ""
fi

#  Run the driver script: Include '--dry_run' to perform a test run without
#+ executing commands
bash "${scr_exc}" \
    --verbose \
    --threads ${threads} \
    --ser_mip "${ser_mip}" \
    --ser_sip "${ser_sip}" \
    --ser_min "${ser_min}" \
    --ser_sin "${ser_sin}" \
    --fil_out "${fil_out}" \
    --err_out "${err_out}" \
    $(if ${slurm:-false}; then echo "--slurm"; fi) \
         >> >(tee -a "${log_exc}.stdout.txt") \
        2>> >(tee -a "${log_exc}.stderr.txt")

if [[ -f "${fil_out}" ]]; then
    #  Sort the table by rows
    awk 'NR == 1; NR > 1 && NF { print | "sort" }' "${fil_out}" \
        > "${dir_tbl}/tmp.tsv"

    #  Replace the original table with the sorted version
    mv -f "${dir_tbl}/tmp.tsv" "${fil_out}"
fi
# cat "${fil_out}"  ## Uncomment to check table ##


#  Cleanup: Compress logs and remove empty files ------------------------------
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${err_out}"
# ls -lhaFG "${err_out}"  ## Uncomment to check directory for logs ##
```
</details>
<br />

<a id="h-compute-log2-ratios-of-ip-to-input-normalized-coverage"></a>
### H. Compute log2 ratios of IP to input normalized coverage.
<a id="text-6"></a>
#### Text
<details>
<summary><i>Text: Compute log2 ratios of IP to input normalized coverage.</i></summary>
<br />

**Overview**  
This section describes how to compute $\log_{2}$-transformed ratios of IP to input normalized coverage tracks, a common method for visualizing ChIP-seq enrichment relative to background. Normalized coverage must first be computed as described [above](#f-compute-normalized-coverage).


In the corresponding Bash code chunk, the signal type is specified by setting `typ_sig` to `"log2"`, which ensures that $\log_2 \left( \text{IP} / \text{input} \right)$ transformations are applied to normalized coverage tracks. The required input files are identified from a sample metadata table, and BEDGRAPH outputs are organized accordingly.

---

**Steps performed in code chunk**  
1. *Define paths and environment:* Specify input normalized coverage BEDGRAPH tracks, output directories, and dependencies.
2.  *Activate environment and check dependencies:* Ensure required software (GNU Parallel, Python, etc.) is available.
3.  *Compute $\log_{2}$ ratios:* Call `execute_compute_signal_ratio.sh` to generate $\log_2 \left( \text{IP} / \text{input} \right)$ BEDGRAPH tracks.
4.  *Perform optional depth filtering:* Generate $\log_{2}$ ratios with and without applying a minimum depth filter (`dep_min`).
5.  *Cleanup:* Compress logs and remove empty files.

---

**Output files**  
1. *BEDGRAPH $\log_2 \left( \text{IP} / \text{input} \right)$ tracks* (`*.bdg.gz`) organized by depth-filtering status.
2. *Log files* tracking execution steps, stored separately.
</details>
<br />

<a id="bash-code-4"></a>
#### Bash code
<details>
<summary><i>Bash code: Compute log2 ratios of IP to input normalized coverage.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
debug=true

#  Define base directory for repository
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"
dir_dat="${dir_rep}/data"
dir_pro="${dir_dat}/processed"

aligner="bowtie2"
a_type="global"
req_flg=true
flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
mapq=1
details="${aligner}_${a_type}_flag-${flg}_mapq-${mapq}"

dir_aln="${dir_pro}/align_fastqs/${details}"
dir_bam="${dir_aln}/sc"
dir_cvg="${dir_pro}/compute_signal/${details}"
dir_nrm="${dir_cvg}/norm"
dir_tbl="${dir_cvg}/tables"
dir_lg2="${dir_cvg}/log2"
dir_log="${dir_lg2}/logs"

env_nam="env_protocol"
scr_exc="execute_compute_signal_ratio.sh"
day="$(date '+%Y-%m%d')"
log_exc="${dir_log}/${day}.$(basename "${scr_exc}" ".sh")"

typ_sig="log2"  ## NOTE: Use samples from either the alpha or spike table ##
tbl_lg2="${dir_tbl}/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_alpha.tsv"

#  Set hardcoded argument assignments
# shellcheck disable=SC1091
source "${dir_fnc}/extract_fld_str.sh"

typ_out="bdg.gz"
fil_ip=$(
    sed \
        -e "s:${dir_bam}:${dir_nrm}:g" \
        -e "s:.bam:.${typ_out}:g" \
        < <(extract_fld_str 1 "${tbl_lg2}")
)
fil_in=$(
    sed \
        -e "s:${dir_bam}:${dir_nrm}:g" \
        -e "s:.bam:.${typ_out}:g" \
        < <(extract_fld_str 2 "${tbl_lg2}")
)
dir_out="${dir_lg2}"
dep_min="$(extract_fld_str 24 "${tbl_lg2}")"  ## WARNING: See description ##
err_out="${dir_log}"
nam_job="compute_signal_ratio_${typ_sig}"

if ${debug:-false}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "debug=${debug}"
    echo ""
    echo "dir_bas=${dir_bas}"
    echo "dir_rep=${dir_rep}"
    echo "dir_scr=${dir_scr}"
    echo "dir_fnc=${dir_fnc}"
    echo "dir_dat=${dir_dat}"
    echo "dir_pro=${dir_pro}"
    echo ""
    echo "aligner=${aligner}"
    echo "a_type=${a_type}"
    echo "req_flg=${req_flg}"
    echo "flg=${flg}"
    echo "mapq=${mapq}"
    echo "details=${details}"
    echo ""
    echo "dir_aln=${dir_aln}"
    echo "dir_bam=${dir_bam}"
    echo "dir_cvg=${dir_cvg}"
    echo "dir_nrm=${dir_nrm}"
    echo "dir_tbl=${dir_tbl}"
    echo "dir_lg2=${dir_lg2}"
    echo "dir_log=${dir_log}"
    echo ""
    echo "env_nam=${env_nam}"
    echo "scr_exc=${scr_exc}"
    echo "day=${day}"
    echo "log_exc=${log_exc}"
    echo ""
    echo "typ_sig=${typ_sig}"
    echo "tbl_lg2=${tbl_lg2}"
    echo ""
    echo "typ_out=${typ_out}"
    echo "fil_ip=${fil_ip}"
    echo "fil_in=${fil_in}"
    echo "dir_out=${dir_out}"
    echo "dep_min=${dep_min}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo ""
    echo ""
fi


#  Create required directories if necessary -----------------------------------
# shellcheck disable=SC2086
mkdir -p ${dir_lg2}/{dm_y,dm_n,logs}

if ${debug:-false}; then
    echo "#################################################"
    echo "## In- and output directory paths and contents ##"
    echo "#################################################"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_nrm %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_nrm}"
    echo ""
    ls -lhaFG "${dir_nrm}/logs"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_tbl %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_tbl}"
    echo ""
    ls -lhaFG "${dir_tbl}/logs"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_lg2 %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_lg2}/dm_y"
    echo ""
    ls -lhaFG "${dir_lg2}/dm_n"
    echo ""
    ls -lhaFG "${err_out}"
    echo ""
    echo ""
fi


#  Activate the environment and check dependencies ----------------------------
#  Source utility functions
# shellcheck disable=SC1091
{
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/echo_warning.sh"
    source "${dir_fnc}/handle_env.sh"
}

#  Activate the required environment
handle_env "${env_nam}"

#  Check the availability of necessary dependencies such as GNU Parallel and
#+ SLURM sbatch
check_program_path awk
check_program_path python
check_program_path samtools

if ${slurm:-false}; then
    check_program_path sbatch &> /dev/null ||
        echo_warning \
            "SLURM is not available on this system. Do not use the '--slurm'" \
            "flag with the driver script."
elif [[ "${threads}" -gt 1 ]]; then
    check_program_path parallel
fi


#  Generate log2(IP/input) normalized coverage tracks -------------------------
if ${debug:-false}; then
    echo "#####################################################"
    echo "## Call to driver script with '--dep_min' argument ##"
    echo "#####################################################"
    echo ""
    echo "bash ${dir_scr}/${scr_exc} \\"
    echo "    --verbose \\"
    echo "    --fil_ip ${fil_ip} \\"
    echo "    --fil_in ${fil_in} \\"
    echo "    --dir_out ${dir_lg2}/dm_y \\"
    echo "    --typ_out ${typ_out} \\"
    echo "    --track \\"
    echo "    --dep_min ${dep_min} \\"
    echo "    --log2 \\"
    echo "    --err_out ${dir_log} \\"
    echo "$(if ${slurm:-false}; then echo "    --slurm \\"; fi)"
    echo "         >> >(tee -a ${log_exc}.stdout.txt) \\"
    echo "        2>> >(tee -a ${log_exc}.stderr.txt)"
    echo ""
    echo ""
fi

#  Run driver script to generate tracks with 'dep_min': Add '--dry_run' to view
#+ commands without execution
bash "${dir_scr}/${scr_exc}" \
    --verbose \
    --fil_ip "${fil_ip}" \
    --fil_in "${fil_in}" \
    --dir_out "${dir_lg2}/dm_y" \
    --typ_out "${typ_out}" \
    --track \
    --dep_min "${dep_min}" \
    --log2 \
    --err_out "${dir_log}" \
    $(if ${slurm:-false}; then echo "--slurm"; fi) \
         >> >(tee -a "${log_exc}.dm_y.stdout.txt") \
        2>> >(tee -a "${log_exc}.dm_y.stderr.txt")

#  Run driver script to generate tracks without 'dep_min': Add '--dry_run' to
#+ view commands without execution
bash "${dir_scr}/${scr_exc}" \
    --verbose \
    --fil_ip "${fil_ip}" \
    --fil_in "${fil_in}" \
    --dir_out "${dir_lg2}/dm_n" \
    --typ_out "${typ_out}" \
    --track \
    --log2 \
    --err_out "${dir_log}" \
    $(if ${slurm:-false}; then echo "--slurm"; fi) \
         >> >(tee -a "${log_exc}.dm_n.stdout.txt") \
        2>> >(tee -a "${log_exc}.dm_n.stderr.txt")


#  Cleanup: Compress logs and remove empty files ------------------------------
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_log}"
# ls -lhaFG "${dir_log}"  ## Uncomment to check directory for track logs ##
# find . -maxdepth 1 -type f -size 0 -delete -o -type f -not -name "*.gz" -exec gzip {} +
```
</details>
<br />

<a id="i-compute-signal-with-the-siq-chip-method"></a>
### I. Compute signal with the siQ-ChIP method.
<a id="text-7"></a>
#### Text
<details>
<summary><i>Text: Compute signal with the siQ-ChIP method.</i></summary>

<a id="overview"></a>
##### Overview
This section describes the steps to compute ChIP-seq signal normalized using the siQ-ChIP method. For more details, see the introduction to [Section G](#g-construct-sample-tables-recording-computed-scaling-factors-for-normalization). The procedure makes use of utility scripts and functions, environment handling, and parallel processing where applicable.

---

<a id="steps-performed-in-code-chunk"></a>
##### Steps performed in code chunk
1. *Set up directories and paths:* Define variables for key directories, data locations, and output locations.
2. *Activate environment and check dependencies:* Load the necessary computational environment and ensure dependencies are available.
3. *Calculate alpha scaling factors:* Use the driver script to compute siQ-ChIP alpha scaling factors and save the sample-specific values to a TSV file. The script can utilize SLURM for job scheduling if available; otherwise, it will use GNU Parallel for parallel processing.
4. *Sort and update output:* Sort the generated output file, replacing it with the sorted version.
5. *Optional cleanup:* Compress large log files and remove empty log files.

---

<a id="important-note"></a>
##### Important note
- The [`execute_calculate_scaling_factor.sh`](https://github.com/kalavattam/protocol_chipseq_signal_norm/blob/main/scripts/execute_calculate_scaling_factor.sh) script in this code chunk requires that *S. cerevisiae* IP BAM files follow a specific naming convention as described in the protocol manuscript. The expected filename format:
    ```txt
    assay_genotype_state_treatment_factor_strain/replicate.
    ```
- Required name components:
    + *assay:* Must be 'IP' or 'in', and must be followed by an underscore.
    + *factor:* A required component preceded by an underscore.
    + *strain/replicate:* A required component preceded by an underscore; it marks the end of the pattern.
- Optional name components:
    + *genotype:* If present, must be preceded by an underscore.
    + *state:* An optional component with preferred values (e.g., 'G1', 'G2M', 'log', or 'Q') but can also be flexible; if present, it must be preceded by an underscore.
    + *treatment:* If present, must be preceded by an underscore.
- Failure to adhere to this naming convention may cause the script to fail.
</details>
<br />

<a id="bash-code-5"></a>
#### Bash code
<details>
<summary><i>Bash code: Compute signal with the siQ-ChIP method.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
debug=true

#  Set hardcoded paths, values, etc.
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"
dir_dat="${dir_rep}/data"
dir_pro="${dir_dat}/processed"

aligner="bowtie2"
a_type="global"
req_flg=true
flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
mapq=1
details="${aligner}_${a_type}_flag-${flg}_mapq-${mapq}"

dir_aln="${dir_pro}/align_fastqs/${details}"
dir_bam="${dir_aln}/sc"
dir_cvg="${dir_pro}/compute_signal/${details}"
dir_nrm="${dir_cvg}/norm"
dir_tbl="${dir_cvg}/tables"
dir_alf="${dir_cvg}/alpha"
dir_log="${dir_alf}/logs"

env_nam="env_protocol"
scr_exc="execute_compute_signal_ratio.sh"
day="$(date '+%Y-%m%d')"
log_exc="${dir_log}/${day}.$(basename "${scr_exc}" ".sh")"

typ_sig="alpha"
tbl_alf="${dir_tbl}/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_${typ_sig}.tsv"

#  Set hardcoded argument assignments
# shellcheck disable=SC1091
source "${dir_fnc}/extract_fld_str.sh"  ## NOTE: Use to parse sample table ##

typ_out="bdg.gz"
fil_ip=$(
    sed \
        -e "s:${dir_bam}:${dir_nrm}:g" \
        -e "s:.bam:.${typ_out}:g" \
        < <(extract_fld_str 1 "${tbl_alf}")
)
fil_in=$(
    sed \
        -e "s:${dir_bam}:${dir_nrm}:g" \
        -e "s:.bam:.${typ_out}:g" \
        < <(extract_fld_str 2 "${tbl_alf}")
)
dir_out="${dir_alf}"
scl_fct="$(extract_fld_str  3 "${tbl_alf}")"
dep_min="$(extract_fld_str 24 "${tbl_alf}")"  ## WARNING: See description ##
err_out="${dir_log}"
nam_job="compute_signal_ratio_${typ_sig}"

if ${debug:-false}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "debug=${debug}"
    echo ""
    echo "dir_bas=${dir_bas}"
    echo "dir_rep=${dir_rep}"
    echo "dir_scr=${dir_scr}"
    echo "dir_fnc=${dir_fnc}"
    echo "dir_dat=${dir_dat}"
    echo "dir_pro=${dir_pro}"
    echo ""
    echo "aligner=${aligner}"
    echo "a_type=${a_type}"
    echo "req_flg=${req_flg}"
    echo "flg=${flg}"
    echo "mapq=${mapq}"
    echo "details=${details}"
    echo ""
    echo "dir_aln=${dir_aln}"
    echo "dir_bam=${dir_bam}"
    echo "dir_cvg=${dir_cvg}"
    echo "dir_nrm=${dir_nrm}"
    echo "dir_tbl=${dir_tbl}"
    echo "dir_alf=${dir_alf}"
    echo "dir_log=${dir_log}"
    echo ""
    echo "env_nam=${env_nam}"
    echo "scr_exc=${scr_exc}"
    echo "day=${day}"
    echo "log_exc=${log_exc}"
    echo ""
    echo "typ_sig=${typ_sig}"
    echo "tbl_alf=${tbl_alf}"
    echo ""
    echo "typ_out=${typ_out}"
    echo "fil_ip=${fil_ip}"
    echo "fil_in=${fil_in}"
    echo "dir_out=${dir_out}"
    echo "scl_fct=${scl_fct}"
    echo "dep_min=${dep_min}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo ""
    echo ""
fi


#  Create required directories if necessary -----------------------------------
# shellcheck disable=SC2086
mkdir -p ${dir_alf}/{dm_y,dm_n,logs}

if ${debug:-false}; then
    echo "#################################################"
    echo "## In- and output directory paths and contents ##"
    echo "#################################################"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_nrm %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_nrm}"
    echo ""
    ls -lhaFG "${dir_nrm}/logs"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_tbl %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_tbl}"
    echo ""
    ls -lhaFG "${dir_tbl}/logs"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_alf %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_alf}/dm_y"
    echo ""
    ls -lhaFG "${dir_alf}/dm_n"
    echo ""
    ls -lhaFG "${err_out}"
    echo ""
    echo ""
fi


#  Activate the environment and check dependencies ----------------------------
#  Source utility functions
# shellcheck disable=SC1091
{
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/echo_warning.sh"
    source "${dir_fnc}/handle_env.sh"
}

#  Activate the required environment
handle_env "${env_nam}"

#  Check the availability of necessary dependencies such as GNU Parallel and
#+ SLURM sbatch
check_program_path awk
check_program_path python
check_program_path samtools

if ${slurm:-false}; then
    check_program_path sbatch &> /dev/null ||
        echo_warning \
            "SLURM is not available on this system. Do not use the '--slurm'" \
            "flag with the driver script."
elif [[ "${threads}" -gt 1 ]]; then
    check_program_path parallel
fi


#  Generate siQ-ChIP alpha-scaled signal tracks -----------------------------
if ${debug:-false}; then
    echo "#####################################################"
    echo "## Call to driver script with '--dep_min' argument ##"
    echo "#####################################################"
    echo ""
    echo "bash ${dir_scr}/${scr_exc} \\"
    echo "    --verbose \\"
    echo "    --fil_ip ${fil_ip} \\"
    echo "    --fil_in ${fil_in} \\"
    echo "    --dir_out ${dir_alf}/dm_y \\"
    echo "    --typ_out ${typ_out} \\"
    echo "    --track \\"
    echo "    --scl_fct ${scl_fct} \\"
    echo "    --dep_min ${dep_min} \\"
    echo "    --err_out ${dir_log} \\"
    echo "$(if ${slurm:-false}; then echo "    --slurm \\"; fi)"
    echo "         >> >(tee -a ${log_exc}.stdout.txt) \\"
    echo "        2>> >(tee -a ${log_exc}.stderr.txt)"
    echo ""
    echo ""
fi

#  Run driver script to generate tracks with 'dep_min': Add '--dry_run' to view
#+ commands without execution
bash "${dir_scr}/${scr_exc}" \
    --verbose \
    --fil_ip "${fil_ip}" \
    --fil_in "${fil_in}" \
    --dir_out "${dir_alf}/dm_y" \
    --typ_out "${typ_out}" \
    --track \
    --scl_fct "${scl_fct}" \
    --dep_min "${dep_min}" \
    --err_out "${dir_log}" \
    $(if ${slurm:-false}; then echo "--slurm"; fi) \
         >> >(tee -a "${log_exc}.dm_y.stdout.txt") \
        2>> >(tee -a "${log_exc}.dm_y.stderr.txt")

#  Run driver script to generate tracks without 'dep_min': Add '--dry_run' to
#+ view commands without execution
bash "${dir_scr}/${scr_exc}" \
    --verbose \
    --fil_ip "${fil_ip}" \
    --fil_in "${fil_in}" \
    --dir_out "${dir_alf}/dm_n" \
    --typ_out "${typ_out}" \
    --track \
    --scl_fct "${scl_fct}" \
    --err_out "${dir_log}" \
    $(if ${slurm:-false}; then echo "--slurm"; fi) \
         >> >(tee -a "${log_exc}.dm_n.stdout.txt") \
        2>> >(tee -a "${log_exc}.dm_n.stderr.txt")


#  Cleanup: Compress logs and remove empty files ------------------------------
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_log}"
# ls -lhaFG "${dir_log}"  ## Uncomment to check directory for track logs ##
```
</details>
<br />

<a id="j-compute-signal-with-the-spike-in-method"></a>
### J. Compute signal with the spike-in method.
<a id="text-8"></a>
#### Text
<details>
<summary><i>Text: Compute signal with the spike-in method.</i></summary>
<br />

This section describes the steps to calculate spike-in normalized ChIP-seq signal. For more details, see the introduction to [Section G](#g-construct-sample-tables-recording-computed-scaling-factors-for-normalization). The procedure makes use of utility scripts and functions, environment handling, and parallel processing where applicable.
</details>
<br />

<a id="bash-code-6"></a>
#### Bash code
<details>
<summary><i>Bash code: Compute signal with the spike-in method.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
debug=true

#  Set hardcoded paths, values, etc.
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"
dir_dat="${dir_rep}/data"
dir_pro="${dir_dat}/processed"

aligner="bowtie2"
a_type="global"
req_flg=true
flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
mapq=1
details="${aligner}_${a_type}_flag-${flg}_mapq-${mapq}"

dir_aln="${dir_pro}/align_fastqs/${details}"
dir_bam="${dir_aln}/sc"
dir_cvg="${dir_pro}/compute_signal/${details}"
dir_nrm="${dir_cvg}/norm"
dir_tbl="${dir_cvg}/tables"
dir_spk="${dir_cvg}/spike"
dir_log="${dir_spk}/logs"

env_nam="env_protocol"
scr_exc="execute_compute_signal_ratio.sh"
day="$(date '+%Y-%m%d')"
log_exc="${dir_log}/${day}.$(basename "${scr_exc}" ".sh")"

typ_sig="spike"
tbl_spk="${dir_tbl}/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_${typ_sig}.tsv"

#  Set hardcoded argument assignments
# shellcheck disable=SC1091
source "${dir_fnc}/extract_fld_str.sh"  ## NOTE: Use to parse sample table ##

typ_out="bdg.gz"
fil_ip=$(
    sed \
        -e "s:${dir_bam}:${dir_nrm}:g" \
        -e "s:.bam:.${typ_out}:g" \
        < <(extract_fld_str 1 "${tbl_spk}")
)
fil_in=$(
    sed \
        -e "s:${dir_bam}:${dir_nrm}:g" \
        -e "s:.bam:.${typ_out}:g" \
        < <(extract_fld_str 3 "${tbl_spk}")
)
dir_out="${dir_spk}"
scl_fct="$(extract_fld_str  5 "${tbl_spk}")"
dep_min="$(extract_fld_str 21 "${tbl_spk}")"  ## WARNING: See description ##
err_out="${dir_log}"
nam_job="compute_signal_ratio_${typ_sig}"

if ${debug:-false}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "debug=${debug}"
    echo ""
    echo "dir_bas=${dir_bas}"
    echo "dir_rep=${dir_rep}"
    echo "dir_scr=${dir_scr}"
    echo "dir_fnc=${dir_fnc}"
    echo "dir_dat=${dir_dat}"
    echo "dir_pro=${dir_pro}"
    echo ""
    echo "aligner=${aligner}"
    echo "a_type=${a_type}"
    echo "req_flg=${req_flg}"
    echo "flg=${flg}"
    echo "mapq=${mapq}"
    echo "details=${details}"
    echo ""
    echo "dir_aln=${dir_aln}"
    echo "dir_bam=${dir_bam}"
    echo "dir_cvg=${dir_cvg}"
    echo "dir_nrm=${dir_nrm}"
    echo "dir_tbl=${dir_tbl}"
    echo "dir_spk=${dir_spk}"
    echo "dir_log=${dir_log}"
    echo ""
    echo "env_nam=${env_nam}"
    echo "scr_exc=${scr_exc}"
    echo "day=${day}"
    echo "log_exc=${log_exc}"
    echo ""
    echo "typ_sig=${typ_sig}"
    echo "tbl_spk=${tbl_spk}"
    echo ""
    echo "typ_out=${typ_out}"
    echo "fil_ip=${fil_ip}"
    echo "fil_in=${fil_in}"
    echo "dir_out=${dir_out}"
    echo "scl_fct=${scl_fct}"
    echo "dep_min=${dep_min}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo ""
    echo ""
fi


#  Create required directories if necessary -----------------------------------
# shellcheck disable=SC2086
mkdir -p ${dir_spk}/{dm_y,dm_n,logs}

if ${debug:-false}; then
    echo "#################################################"
    echo "## In- and output directory paths and contents ##"
    echo "#################################################"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_nrm %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_nrm}"
    echo ""
    ls -lhaFG "${dir_nrm}/logs"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_tbl %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_tbl}"
    echo ""
    ls -lhaFG "${dir_tbl}/logs"
    echo ""
    echo "%%%%%%%%%%%%%"
    echo "%% dir_spk %%"
    echo "%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_spk}/dm_y"
    echo ""
    ls -lhaFG "${dir_spk}/dm_n"
    echo ""
    ls -lhaFG "${err_out}"
    echo ""
    echo ""
fi


#  Activate the environment and check dependencies ----------------------------
#  Source utility functions
# shellcheck disable=SC1091
{
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/echo_warning.sh"
    source "${dir_fnc}/handle_env.sh"
}

#  Activate the required environment
handle_env "${env_nam}"

#  Check the availability of necessary dependencies such as GNU Parallel and
#+ SLURM sbatch
check_program_path awk
check_program_path python
check_program_path samtools

if ${slurm:-false}; then
    check_program_path sbatch &> /dev/null ||
        echo_warning \
            "SLURM is not available on this system. Do not use the '--slurm'" \
            "flag with the driver script."
elif [[ "${threads}" -gt 1 ]]; then
    check_program_path parallel
fi


#  Generate spike-in scaled signal tracks -------------------------------------
if ${debug:-false}; then
    echo "#####################################################"
    echo "## Call to driver script with '--dep_min' argument ##"
    echo "#####################################################"
    echo ""
    echo "bash ${dir_scr}/${scr_exc} \\"
    echo "    --verbose \\"
    echo "    --fil_ip ${fil_ip} \\"
    echo "    --fil_in ${fil_in} \\"
    echo "    --dir_out ${dir_spk}/dm_y \\"
    echo "    --typ_out ${typ_out} \\"
    echo "    --track \\"
    echo "    --scl_fct ${scl_fct} \\"
    echo "    --dep_min ${dep_min} \\"
    echo "    --err_out ${dir_log} \\"
    echo "$(if ${slurm:-false}; then echo "    --slurm \\"; fi)"
    echo "         >> >(tee -a ${log_exc}.stdout.txt) \\"
    echo "        2>> >(tee -a ${log_exc}.stderr.txt)"
    echo ""
    echo ""
fi

#  Run driver script to generate tracks with 'dep_min': Add '--dry_run' to view
#+ commands without execution
bash "${dir_scr}/${scr_exc}" \
    --verbose \
    --fil_ip "${fil_ip}" \
    --fil_in "${fil_in}" \
    --dir_out "${dir_spk}/dm_y" \
    --typ_out "${typ_out}" \
    --track \
    --scl_fct "${scl_fct}" \
    --dep_min "${dep_min}" \
    --err_out "${dir_log}" \
    $(if ${slurm:-false}; then echo "--slurm"; fi) \
         >> >(tee -a "${log_exc}.dm_y.stdout.txt") \
        2>> >(tee -a "${log_exc}.dm_y.stderr.txt")

#  Run driver script to generate tracks without 'dep_min': Add '--dry_run' to
#+ view commands without execution
bash "${dir_scr}/${scr_exc}" \
    --verbose \
    --fil_ip "${fil_ip}" \
    --fil_in "${fil_in}" \
    --dir_out "${dir_spk}/dm_n" \
    --typ_out "${typ_out}" \
    --track \
    --scl_fct "${scl_fct}" \
    --err_out "${dir_log}" \
    $(if ${slurm:-false}; then echo "--slurm"; fi) \
         >> >(tee -a "${log_exc}.dm_n.stdout.txt") \
        2>> >(tee -a "${log_exc}.dm_n.stderr.txt")


#  Cleanup: Compress logs and remove empty files ------------------------------
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_log}"
# ls -lhaFG "${dir_log}"  ## Uncomment to check directory for track logs ##
```
</details>
<br />
