
ChIP-seq Protocol Workflow
==========================

**Supporting code and documentation for the *bioRxiv* manuscript "ChIP-seq data processing and relative and quantitative signal normalization for *Saccharomyces cerevisiae*."**

**Author:** *Kris Alavattam*

This notebook provides a guide to the ChIP-seq data processing workflow detailed in the manuscript, including code snippets, explanations, and step-by-step instructions.

Note: If using a high-performance computing cluster (HPCC), request an interactive node to ensure adequate resources for running code in the below chunks. The specific command (e.g., `grabnode` at Fred Hutch Cancer Center) will depend on the job scheduler setup. (This step is unnecessary if running the code on a local machine.)

Note: For detailed instructions on keeping your local version of the [`protocol_chipseq_signal_norm`](https://github.com/kalavattam/protocol_chipseq_signal_norm) repository up-to-date, please see [this GitHub gist](https://gist.github.com/kalavattam/76f123011e8dcd77b445a72d23a64036).

---
<br />

## Table of contents
<details>
<summary><i>Table of contents</i></summary>
<br />
<!-- MarkdownTOC -->

1. [Procedures](#procedures)
    1. [A. Install and configure Miniforge.](#a-install-and-configure-miniforge)
    1. [B. Clone the protocol repository and install project environments.](#b-clone-the-protocol-repository-and-install-project-environments)
    1. [C. Clone the siQ-ChIP repository and install its environment.](#c-clone-the-siq-chip-repository-and-install-its-environment)
1. [Data analysis](#data-analysis)
    1. [A. Prepare and concatenate FASTA and GFF3 files for model and spike-in organisms.](#a-prepare-and-concatenate-fasta-and-gff3-files-for-model-and-spike-in-organisms)
    1. [B. Generate Bowtie 2 indices from the concatenated FASTA file.](#b-generate-bowtie-2-indices-from-the-concatenated-fasta-file)
    1. [C. Obtain and organize ChIP-seq FASTQ files.](#c-obtain-and-organize-chip-seq-fastq-files)
    1. [D. Use Atria to perform adapter and quality trimming of sequenced reads.](#d-use-atria-to-perform-adapter-and-quality-trimming-of-sequenced-reads)
    1. [E. Align sequenced reads with Bowtie 2 and process the read alignments.](#e-align-sequenced-reads-with-bowtie-2-and-process-the-read-alignments)
    1. [F. Construct sample tables recording computed scaling factors for normalization.](#f-construct-sample-tables-recording-computed-scaling-factors-for-normalization)
        1. [1. Construct sample tables recording siQ-ChIP alpha scaling factors.](#1-construct-sample-tables-recording-siq-chip-alpha-scaling-factors)
        1. [2. Construct sample tables recording spike-in scaling factors.](#2-construct-sample-tables-recording-spike-in-scaling-factors)
    1. [G. Compute normalized coverage.](#g-compute-normalized-coverage)
    1. [H. Compute coverage with the sans spike-in quantitative ChIP-seq \(siQ-ChIP\) method.](#h-compute-coverage-with-the-sans-spike-in-quantitative-chip-seq-siq-chip-method)
    1. [I. Compute coverage with the spike-in method.](#i-compute-coverage-with-the-spike-in-method)
    1. [J. Compute log2 ratios of IP to input coverage.](#j-compute-log2-ratios-of-ip-to-input-coverage)

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

`#TODO`
</details>
<br />

<details>
<summary><i>Bash code: Install and configure Miniforge.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
#  Define variable for base directory
dir_bas="${HOME}"  ## WARNING: Change as needed ##

#  Determine the appropriate Miniforge installer to use based on operating
#+ system (OS) and system architecture
case $(uname -s) in
    Darwin) os="MacOSX" ;;
    Linux)  os="Linux"  ;;
    *) echo "Error: Unsupported operating system: '$(uname -s)'." >&2 ;;
esac

ar=$(uname -m)  # e.g., "x86_64" for Intel/AMD, "arm64" for ARM

#  Set Miniforge installer URL and script name
https="https://github.com/conda-forge/miniforge/releases/latest/download"
script="Miniforge3-${os}-${ar}.sh"


#  Download and install Miniforge ---------------------------------------------
#  Move to base directory
cd "${dir_bas}" || echo "Error: Failed to change directory: '${dir_bas}'." >&2

#  Download Miniforge installer
curl -L -O "${https}/${script}"

#  Run Miniforge installer
bash "${script}"
```
</details>
<br />

<a id="b-clone-the-protocol-repository-and-install-project-environments"></a>
### B. Clone the protocol repository and install project environments.
<details>
<summary><i>Text: Clone the protocol repository and install project environments.</i></summary>
<br />

`#TODO`
</details>
<br />

<details>
<summary><i>Bash code: Clone the protocol repository and install project environments.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
repo="protocol_chipseq_signal_norm"
https="https://github.com/kalavattam/${repo}.git"


#  Do the main work -----------------------------------------------------------
#  Make the base directory
mkdir -p "${dir_bas}"

#  Move to base directory
cd "${dir_bas}" || echo "Error: Failed to change directory: '${dir_bas}'." >&2

#  Clone the protocol repository
git clone "${https}"

#  Move to repository directory
cd "${repo}" || echo "Error: Failed to change directory: '${repo}'." >&2

#  Install Conda/Mamba environments with install_envs.sh
bash "scripts/install_envs.sh" --env_nam "env_align" --yes
bash "scripts/install_envs.sh" --env_nam "env_analyze" --yes
```
</details>
<br />

<a id="c-clone-the-siq-chip-repository-and-install-its-environment"></a>
### C. Clone the siQ-ChIP repository and install its environment.
<details>
<summary><i>Text: Clone the siQ-ChIP repository and install its environment.</i></summary>
<br />

`#TODO`
</details>
<br />

<details>
<summary><i>Bash code: Clone the siQ-ChIP repository and install its environment.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/siQ-ChIP"
https="https://github.com/kalavattam/siQ-ChIP.git"
branch="protocol"

dir_scr="${dir_bas}/protocol_chipseq_signal_norm/scripts"


#  Do the main work -----------------------------------------------------------
#  Move to base directory
cd "${dir_bas}" || echo "Error: Failed to change directory: '${dir_bas}'." >&2

#  Clone siQ-ChIP repository
git clone "${https}"

#  Move to siQ-ChIP repository directory
cd "${dir_rep}" || echo "Error: Failed to change directory: '${dir_rep}'." >&2

#  Switch to specified branch
git checkout "${branch}"

#  Install Conda/Mamba environment for siQ-ChIP
bash "${dir_scr}/install_envs.sh" --env_nam "env_siq" --yes
```
</details>
<br />
<br />

<a id="data-analysis"></a>
## Data analysis
<a id="a-prepare-and-concatenate-fasta-and-gff3-files-for-model-and-spike-in-organisms"></a>
### A. Prepare and concatenate FASTA and GFF3 files for model and spike-in organisms.
<details>
<summary><i>Text: Prepare and concatenate FASTA and GFF3 files for model and spike-in organisms.</i></summary>
<br />

`#TODO`
</details>
<br />

<details>
<summary><i>Bash code: Prepare and concatenate FASTA and GFF3 files for model and spike-in organisms.</i></summary>
<br />

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##
```
</details>
<br />

<a id="b-generate-bowtie-2-indices-from-the-concatenated-fasta-file"></a>
### B. Generate Bowtie 2 indices from the concatenated FASTA file.
<details>
<summary><i>Text: Generate Bowtie 2 indices from the concatenated FASTA file.</i></summary>
<br />

To align ChIP-seq reads against both *S. cerevisiae* and *S. pombe* genomes, we first generate Bowtie 2 indices from a concatenated FASTA file. This ensures efficient and accurate alignment for, e.g., spike-in normalization (described below).

**Steps overview:**
1. *Define directories and files:* Set paths for inputs and outputs.
2. *Activate environment:* Load necessary tools and dependencies.
3. *Run Bowtie 2 index creation:* Use the concatenated FASTA file to generate indices, logging output for troubleshooting.
4. *Optional cleanup:* Remove the decompressed FASTA file to save space.
</details>
<br />

<details>
<summary><i>Bash code: Generate Bowtie 2 indices from the concatenated FASTA file.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
#  Define variables for directory paths, etc.
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"
dir_dat="${dir_rep}/data"
dir_gen="${dir_dat}/genomes"
dir_cat="${dir_gen}/concat"
dir_fas="${dir_cat}/fasta/proc"
fil_fas="sc_sp_proc.fasta"
pth_fas="${dir_fas}/${fil_fas}"
dir_idx="${dir_cat}/index/bowtie2"
env_nam="env_align"
day="$(date '+%Y-%m%d')"


#  Do the main work -----------------------------------------------------------
#  Source utility functions
source "${dir_fnc}/check_program_path.sh"
source "${dir_fnc}/handle_env.sh"

#  Activate the required environment
handle_env "${env_nam}"

#  Ensure access to bowtie2-build
check_program_path "bowtie2-build"

#  Create output directory structure for Bowtie 2 index files and logs
mkdir -p ${dir_idx}/{docs,logs}

#  If necessary, decompress the FASTA file
if [[ ! -f "${pth_fas}" && -f "${pth_fas}.gz" ]]; then
    gunzip -c "${pth_fas}.gz" > "${pth_fas}"
fi

#  "Build" the Bowtie 2 index using the decompressed FASTA file
bowtie2-build "${pth_fas}" "${dir_idx}/${fil_fas%.fasta}" \
     > >(tee -a "${dir_idx}/logs/${day}.execute.stdout.txt") \
    2> >(tee -a "${dir_idx}/logs/${day}.execute.stderr.txt")

#  Optional cleanup: Once the index is built, delete the decompressed FASTA
#+ file
if [[ -f "${pth_fas}" ]]; then rm "${pth_fas}"; fi

#  Cleanup: Compress large stdout and stderr files, and remove files with size
#+ 0
bash "${dir_scr}/compress_remove_files.sh" \
    --pattern "*.txt" \
    --dir_fnd "${dir_idx}/logs"

#  Optional: Check the contents of the logs directory
# ls -lhaFG "${dir_idx}/logs"
```
</details>
<br />

<a id="c-obtain-and-organize-chip-seq-fastq-files"></a>
### C. Obtain and organize ChIP-seq FASTQ files.
<details>
<summary><i>Bash code: Obtain and organize ChIP-seq FASTQ files.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
#  Define variables for directory paths, etc.
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"
dir_raw="${dir_rep}/data/raw"
dir_doc="${dir_raw}/docs"
fil_tbl="datasets.tsv"  ## WARNING: Change as needed ##
pth_tsv="${dir_doc}/${fil_tbl}"
dir_log="${dir_raw}/logs"
dir_sym="${dir_rep}/data/symlinked"
env_nam="env_align"
threads=6
time="6:00:00"
day="$(date '+%Y-%m%d')"

#  Source utility functions
source "${dir_fnc}/check_program_path.sh"
source "${dir_fnc}/echo_warning.sh"
source "${dir_fnc}/handle_env.sh"

#  Activate the required environment
handle_env "${env_nam}"

#  Ensure access to necessary dependencies such as GNU Parallel, SLURM sbatch
check_program_path parallel
check_program_path sbatch ||
    echo_warning \
        "SLURM is not available on this system. Do not use the '--slurm'" \
        "flag with the driver script."

#  If needed, create directory structure for raw and symlinked FASTQ files and
#+ logs
mkdir -p ${dir_raw}/{docs,logs}
mkdir -p ${dir_sym}/{docs,logs}

#  Download and symlink FASTQ files 
bash "${dir_scr}/execute_download_fastqs.sh" \
    --threads "${threads}" \
    --infile "${pth_tsv}" \
    --dir_out "${dir_raw}" \
    --dir_sym "${dir_sym}" \
    --err_out "${dir_raw}/logs" \
    --slurm \
    --time "${time}" \
         > >(tee -a "${dir_raw}/logs/${day}.execute.stdout.txt") \
        2> >(tee -a "${dir_raw}/logs/${day}.execute.stderr.txt")

#  Cleanup: Compress large stdout and stderr files, and remove files with size
#+ 0
bash "${dir_scr}/compress_remove_files.sh" \
    --dir_fnd "${dir_raw}/logs"
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
#  Define variables for directory paths, environment, threads, and infiles
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"
dir_dat="${dir_rep}/data"
dir_sym="${dir_dat}/symlinked"
dir_pro="${dir_dat}/processed"
dir_trm="${dir_pro}/trim_atria"
env_nam="env_analyze"
threads=4
infiles="$(  ## WARNING: Change search parameters as needed ##
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_sym}" \
        --pattern "*.fastq.gz" \
        --depth 1 \
        --follow \
        --fastqs
)"
day="$(date '+%Y-%m%d')"

#  Source utility functions
source "${dir_fnc}/check_program_path.sh"
source "${dir_fnc}/echo_warning.sh"
source "${dir_fnc}/handle_env.sh"

#  Activate the required environment
handle_env "${env_nam}"

#  Check availability of Atria and other necessary tools
check_program_path atria
check_program_path pbzip2
check_program_path pigz
check_program_path sbatch ||
    echo_warning \
        "SLURM is not available on this system. Do not use the '--slurm'" \
        "flag with the driver script."

#  Create output directory structure for trimmed FASTQ files and logs
mkdir -p ${dir_trm}/{docs,logs}

#  Run the driver script to trim FASTQ files with Atria
bash "${dir_scr}/execute_trim_fastqs.sh" \
    --verbose \
    --threads ${threads} \
    --infiles "${infiles}" \
    --dir_out "${dir_trm}" \
    --err_out "${dir_trm}/logs" \
    --slurm \
         >> >(tee -a "${dir_trm}/logs/${day}.execute.stdout.txt") \
        2>> >(tee -a "${dir_trm}/logs/${day}.execute.stderr.txt")

#  Cleanup: Move Atria LOG and JSON files to the logs directory
mv ${dir_trm}/*.{log,json} "${dir_trm}/logs"

#  Cleanup: Compress large stdout, stderr, LOG, and JSON files, and remove
#+ files with size 0
bash "${dir_scr}/compress_remove_files.sh" \
    --dir_fnd "${dir_trm}/logs"

bash "${dir_scr}/compress_remove_files.sh" \
    --dir_fnd "${dir_trm}/logs" \
    --pattern "*.log"

bash "${dir_scr}/compress_remove_files.sh" \
    --dir_fnd "${dir_trm}/logs" \
    --pattern "*.json"

#  Optional: Check the contents of the logs directory
# ls -lhaFG "${dir_trm}/logs"
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
#  Define variables for directory paths, environment, driver script arguments,
#+ and so on
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"
dir_dat="${dir_rep}/data"
dir_idx="${dir_dat}/genomes/concat/index"
dir_pro="${dir_dat}/processed"
dir_trm="${dir_pro}/trim_atria"
env_nam="env_align"
threads=8
aligner="bowtie2"
a_type="global"
req_flg=true
flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
mapq=1
str_idx="sc_sp_proc"
pth_idx="${dir_idx}/${aligner}/${str_idx}"
infiles="$(  ## WARNING: Change search parameters as needed ##
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_trm}" \
        --pattern "*.atria.fastq.gz" \
        --depth 1 \
        --fastqs
)"
dir_aln="${dir_pro}/align_${aligner}_${a_type}"
dir_out="${dir_aln}/flag-${flg}_mapq-${mapq}"
nam_job="align_fastqs"
max_job=6
time="1:00:00"
day="$(date '+%Y-%m%d')"

#  Create output directory structure for trimmed FASTQ files and logs
mkdir -p ${dir_out}/{init,sc,sp}/{docs,logs}

#  Source utility functions
source "${dir_fnc}/check_program_path.sh"
source "${dir_fnc}/echo_warning.sh"
source "${dir_fnc}/handle_env.sh"

#  Activate the required environment
handle_env "${env_nam}"

#  Check the availability of Bowtie 2, Samtools, and SLURM sbatch
check_program_path bowtie2
check_program_path samtools
check_program_path sbatch ||
    echo_warning \
        "SLURM is not available on this system. Do not use the '--slurm'" \
        "flag with the driver script."

#  Run the driver script to align and post-process FASTQ files
bash "${dir_scr}/execute_align_fastqs.sh" \
    --verbose \
    --threads "${threads}" \
    --aligner "${aligner}" \
    --a_type "${a_type}" \
    --mapq "${mapq}" \
    --req_flg \
    --index "${pth_idx}" \
    --infiles "${infiles}" \
    --dir_out "${dir_out}/init" \
    --err_out "${dir_out}/init/logs" \
    --slurm \
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
bash "${dir_scr}/execute_filter_bams.sh" \
    --verbose \
    --threads "${threads}" \
    --infiles "${infiles}" \
    --dir_out "${dir_out}/sc" \
    --err_out "${dir_out}/sc/logs" \
    --retain "sc" \
    --slurm \
         >> >(tee -a "${dir_out}/sc/logs/${day}.execute.stdout.txt") \
        2>> >(tee -a "${dir_out}/sc/logs/${day}.execute.stderr.txt")

#  Run the driver script to filter BAM files for S. pombe ("sp") alignments
#+ (i.e., the "spike-in" alignments)
bash "${dir_scr}/execute_filter_bams.sh" \
    --verbose \
    --threads "${threads}" \
    --infiles "${infiles}" \
    --dir_out "${dir_out}/sp" \
    --err_out "${dir_out}/sp/logs" \
    --retain "sp" \
    --slurm \
         >> >(tee -a "${dir_out}/sp/logs/${day}.execute.stdout.txt") \
        2>> >(tee -a "${dir_out}/sp/logs/${day}.execute.stderr.txt")

#  Cleanup: Compress large stdout and stderr files, and remove files with size
#+ 0
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_out}/init/logs"

bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_out}/sc/logs"

bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_out}/sp/logs"

#  Optional: Check the contents of the logs directory
# ls -lhaFG "${dir_trm}/init/logs"
# ls -lhaFG "${dir_trm}/sc/logs"
# ls -lhaFG "${dir_trm}/sp/logs"
```
</details>
<br />

<a id="f-construct-sample-tables-recording-computed-scaling-factors-for-normalization"></a>
### F. Construct sample tables recording computed scaling factors for normalization.
<details>
<summary><i>Text: Construct sample tables recording computed scaling factors for normalization.</i></summary>
<br />

The following Bash code chunk provides a structured approach to generating sample tables that record scaling factors derived from the siQ-ChIP and spike-in methods. These factors are stored in tables for use in downstream processing, where they are applied to ["normalized coverage" data](https://www.nature.com/articles/s41598-023-34430-2#Sec2) to produce scaled coverage tracks. This facilitates ***quantitative*** (siQ-ChIP) and ***semiquantitative*** (spike-in) comparisons across samples.

To generate siQ-ChIP- and spike-in-scaled coverage tracks, we normalize IP (immunoprecipitation) data by dividing it by input (control) data and multiplicatively applying the corresponding scaling factors. *(Note: In siQ-ChIP calculations, the scaling factor is referred to as "alpha".)*

As part of this process, we derive ***minimum input depth factors*** to prevent erroneous (e.g., division by zero) or extreme (e.g., division by anomalously small input values) calculations, promoting numerical stability in resulting scaled coverage tracks. Specifically, when a binned input value falls below the minimum input depth factor, it is replaced by the factor to prevent instability.

To support different analysis resolutions, we compute minimum input depth factors for a range of common genomic bin sizes: 1 (i.e., no binning), 5, 10, 20, 30, 40, and 50 bp. *(Note: The application of minimum input depth factors during scaled coverage track generation is optional. See ... for more details.)*
</details>
<br />

<a id="1-construct-sample-tables-recording-siq-chip-alpha-scaling-factors"></a>
#### 1. Construct sample tables recording siQ-ChIP alpha scaling factors.
<details>
<summary><i>Bash code: Construct sample tables recording siQ-ChIP alpha scaling factors.</i></summary>

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

dir_aln="${dir_pro}/align_reads/${details}"
dir_bam="${dir_aln}/sc"
dir_cvg="${dir_pro}/compute_coverage/${details}"
dir_tbl="${dir_cvg}/tables"
typ_cvg="alpha"

pattern="*.bam"
exclude="*Brn1*"

env_nam="env_analyze"
execute="${dir_scr}/execute_calculate_scaling_factor_${typ_cvg}.sh"

day="$(date '+%Y-%m%d')"
exc_log="${dir_trk}/logs/${day}.$(basename "${execute}" ".sh")"

#  Set hardcoded argument assignments
threads=8
ser_ip="$(  ## WARNING: Change search parameters as needed ##
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_bam}" \
        --pattern "${pattern}" \
        --include "IP*" \
        --exclude "${exclude}"
)"
ser_in="$(sed 's:IP:in:g' < <(echo "${ser_ip}"))"
tbl_met="${dir_dat}/raw/docs/measurements_siqchip.tsv"
eqn="6nd"
fil_out="${dir_tbl}/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_${typ_cvg}.tsv"
err_out="${dir_tbl}/logs"


#  Create required directories if necessary -----------------------------------
mkdir -p "${dir_tbl}/logs"


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
check_program_path parallel
check_program_path python
check_program_path samtools
check_program_path sbatch ||
    echo_warning \
        "SLURM is not available on this system. Do not use the '--slurm'" \
        "flag with the driver script."


#  Construct table ------------------------------------------------------------
if ${debug}; then
    echo "###########################"
    echo "## Call to driver script ##"
    echo "###########################"
    echo ""
    echo "bash \"${execute}\" \\"
    echo "    --verbose \\"
    echo "    --threads ${threads} \\"
    echo "    --ser_ip \"${ser_ip}\" \\"
    echo "    --ser_in \"${ser_in}\" \\"
    echo "    --tbl_met \"${tbl_met}\" \\"
    echo "    --eqn \"${eqn}\" \\"
    echo "    --fil_out \"${fil_out}\" \\"
    echo "    --err_out \"${err_out}\" \\"
    echo "    --slurm"
    echo "         >> >(tee -a ${exc_log}.stdout.txt) \\"
    echo "        2>> >(tee -a ${exc_log}.stderr.txt)"
    echo ""
    echo ""
fi

#  Run the driver script: Change '--verbose' to '--dry_run' to perform a 
#+ test run without executing commands
bash "${execute}" \
    --verbose \
    --threads ${threads} \
    --ser_ip "${ser_ip}" \
    --ser_in "${ser_in}" \
    --tbl_met "${tbl_met}" \
    --eqn "${eqn}" \
    --fil_out "${fil_out}" \
    --err_out "${err_out}" \
    --slurm \
         >> >(tee -a "${exc_log}.stdout.txt") \
        2>> >(tee -a "${exc_log}.stderr.txt")

if [[ -f "${fil_out}" ]]; then
    #  Sort the table of input minimume depths by rows
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

<a id="2-construct-sample-tables-recording-spike-in-scaling-factors"></a>
#### 2. Construct sample tables recording spike-in scaling factors.
<details>
<summary><i>Bash code: Construct sample tables recording spike-in scaling factors.</i></summary>

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

dir_aln="${dir_pro}/align_reads/${details}"
dir_bam="${dir_aln}/sc"
dir_cvg="${dir_pro}/compute_coverage/${details}"
dir_tbl="${dir_cvg}/tables"
typ_cvg="spike"

pattern="*.bam"
exclude="*Brn1*"

env_nam="env_analyze"
execute="${dir_scr}/execute_calculate_scaling_factor_${typ_cvg}.sh"

day="$(date '+%Y-%m%d')"
exc_log="${dir_trk}/logs/${day}.$(basename "${execute}" ".sh")"

#  Set hardcoded argument assignments
threads=8
ser_mip="$(  ## WARNING: Change search parameters as needed ##
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_bam}" \
        --pattern "${pattern}" \
        --include "IP*" \
        --exclude "${exclude}"
)"
ser_sip="$(sed 's:sc:sp:g' < <(echo "${ser_mip}"))"
ser_min="$(sed 's:IP:in:g' < <(echo "${ser_mip}"))"
ser_sin="$(sed 's:IP:in:g' < <(echo "${ser_sip}"))"
fil_out="${dir_tbl}/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_${typ_cvg}.tsv"
err_out="${dir_tbl}/logs"


#  Create required directories if necessary -----------------------------------
mkdir -p "${dir_tbl}/logs"


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
check_program_path parallel
check_program_path python
check_program_path samtools
check_program_path sbatch ||
    echo_warning \
        "SLURM is not available on this system. Do not use the '--slurm'" \
        "flag with the driver script."


#  Construct table ------------------------------------------------------------
if ${debug}; then
    echo "###########################"
    echo "## Call to driver script ##"
    echo "###########################"
    echo ""
    echo "bash \"${execute}\" \\"
    echo "    --verbose \\"
    echo "    --threads ${threads} \\"
    echo "    --ser_mip \"${ser_mip}\" \\"
    echo "    --ser_sip \"${ser_sip}\" \\"
    echo "    --ser_min \"${ser_min}\" \\"
    echo "    --ser_sin \"${ser_sin}\" \\"
    echo "    --fil_out \"${fil_out}\" \\"
    echo "    --err_out \"${err_out}\" \\"
    echo "    --slurm"
    echo "         >> >(tee -a ${exc_log}.stdout.txt) \\"
    echo "        2>> >(tee -a ${exc_log}.stderr.txt)"
    echo ""
    echo ""
fi

#  Run the driver script: Change '--verbose' to '--dry_run' to perform a 
#+ test run without executing commands
bash "${execute}" \
    --verbose \
    --threads ${threads} \
    --ser_mip "${ser_mip}" \
    --ser_sip "${ser_sip}" \
    --ser_min "${ser_min}" \
    --ser_sin "${ser_sin}" \
    --fil_out "${fil_out}" \
    --err_out "${err_out}" \
    --slurm \
         >> >(tee -a "${exc_log}.stdout.txt") \
        2>> >(tee -a "${exc_log}.stderr.txt")

if [[ -f "${fil_out}" ]]; then
    #  Sort the table of input minimume depths by rows
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

<a id="g-compute-normalized-coverage"></a>
### G. Compute normalized coverage.
<details>
<summary><i>Text: Compute normalized raw coverage..</i></summary>
<br />

This following Bash code chunk provides an example of how to compute ChIP-seq normalized coverage per [Dickson et al., *Sci Rep*, 2023](https://www.nature.com/articles/s41598-023-34430-2). The coverage type is determined by setting the variable `typ_cvg` to `"norm"` (other options are `raw` and `frag`; see the [`compute_coverage.py`](https://github.com/kalavattam/protocol_chipseq_signal_norm/blob/main/scripts/compute_coverage.py) documentation for more details). BEDGRAPH and log output files will be saved to separate directories based on the selected coverage type.

<b>Workflow overview</b>  
1. Set paths and environment: Define paths for input BAM files, output directories, and required dependencies.
2. Activate environment and check dependencies: Ensure required software (GNU Parallel, Python, SLURM, etc.) is available.
3. Compute coverage: Calls the driver script (`execute_compute_coverage.sh`) to generate BEDGRAPH coverage tracks.
4. Verify normalization: Optionally perform a check to confirm the output normalized coverage sums to ~1.
5. Cleanup: Compress logs and remove size-0 files.

<b>Output files</b>
1. BEDGRAPH coverage tracks (`*.bdg.gz`) are saved in structured directories based on the chosen normalization type.
2. Log files tracking the execution process are stored separately.
</details>
<br />

<details>
<summary><i>Bash code: Compute normalized coverage.</i></summary>

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

dir_aln="${dir_pro}/align_reads/${details}"
dir_bam="${dir_aln}/sc"
dir_cvg="${dir_pro}/compute_coverage/${details}"
typ_trk="norm"
dir_trk="${dir_cvg}/${typ_trk}"

pattern="*.bam"

env_nam="env_analyze"
execute="${dir_scr}/execute_compute_coverage.sh"

day="$(date '+%Y-%m%d')"
exc_log="${dir_trk}/logs/${day}.$(basename "${execute}" ".sh")"

#  Set hardcoded argument assignments
threads=8
infiles="$(  ## WARNING: Change search parameters as needed ##
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_bam}" \
            --pattern "${pattern}"
    )"
dir_out="${dir_trk}"
typ_out="bdg.gz"
siz_bin=30
typ_cvg="${typ_trk}"
err_out="${dir_trk}/logs"
nam_job="compute_coverage_${typ_cvg}"

if ${debug}; then
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
    echo "flg=${flg}"
    echo "mapq=${mapq}"
    echo "details=${details}"
    echo ""
    echo "dir_aln=${dir_aln}"
    echo "dir_bam=${dir_bam}"
    echo "dir_cvg=${dir_cvg}"
    echo "typ_trk=${typ_trk}"
    echo "dir_trk=${dir_trk}"
    echo ""
    echo "pattern=${pattern}"
    echo ""
    echo "env_nam=${env_nam}"
    echo "execute=${execute}"
    echo ""
    echo "day=${day}"
    echo "exc_log=${exc_log}"
    echo ""
    echo "threads=${threads}"
    echo "infiles=${infiles}"
    echo "dir_out=${dir_out}"
    echo "typ_out=${typ_out}"
    echo "siz_bin=${siz_bin}"
    echo "typ_cvg=${typ_cvg}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo ""
    echo ""
fi


#  Create required directories if necessary -----------------------------------
mkdir -p "${dir_out}/logs"


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
check_program_path parallel
check_program_path python
check_program_path sbatch ||
    echo_warning \
        "SLURM is not available on this system. Do not use the '--slurm'" \
        "flag with the driver script."


#  Compute coverage -----------------------------------------------------------
if ${debug}; then
    echo "###########################"
    echo "## Call to driver script ##"
    echo "###########################"
    echo ""
    echo "bash ${execute} \\"
    echo "    --verbose \\"
    echo "    --threads ${threads} \\"
    echo "    --infiles ${infiles} \\"
    echo "    --dir_out ${dir_trk} \\"
    echo "    --typ_out ${typ_out} \\"
    echo "    --siz_bin ${siz_bin} \\"
    echo "    --typ_cvg ${typ_cvg} \\"
    echo "    --err_out ${err_out} \\"
    echo "    --nam_job ${nam_job} \\"
    echo "    --slurm \\"
    echo "         >> >(tee -a ${exc_log}.stdout.txt) \\"
    echo "        2>> >(tee -a ${exc_log}.stderr.txt)"
    echo ""
    echo ""
fi

#  Run driver script
bash "${execute}" \
    --verbose \
    --threads "${threads}" \
    --infiles "${infiles}" \
    --dir_out "${dir_trk}" \
    --typ_out "${typ_out}" \
    --siz_bin "${siz_bin}" \
    --typ_cvg "${typ_cvg}" \
    --err_out "${err_out}" \
    --nam_job "${nam_job}" \
    --slurm \
         >> >(tee -a "${exc_log}.stdout.txt") \
        2>> >(tee -a "${exc_log}.stderr.txt")

#  Check that each normalized coverage BEDGRAPH file sums to unity
if ${debug}; then
    echo "##############################################"
    echo "## Check unity of normalized coverage files ##"
    echo "##############################################"
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

<a id="h-compute-coverage-with-the-sans-spike-in-quantitative-chip-seq-siq-chip-method"></a>
### H. Compute coverage with the sans spike-in quantitative ChIP-seq (siQ-ChIP) method.
<details>
<summary><i>Text: Compute coverage with the sans spike-in quantitative ChIP-seq (siQ-ChIP) method.</i></summary>
<br />

This section describes the steps to compute ChIP-seq coverage normalized using the siQ-ChIP method. The approach involves... `#TODO`. The procedure makes use of utility scripts and functions, environment handling, and parallel processing where applicable.

**Steps overview:**
1. *Set up directories and paths:* Define variables for key directories, data locations, and output destinations.
2. *Activate environment and check dependencies:* Load the necessary computational environment and ensure essential dependencies are available.
3. *Calculate alpha scaling factors:* Use the driver script to compute siQ-ChIP alpha scaling factors and save the sample-specific values to a TSV file. The script can utilize SLURM for job scheduling if available; otherwise, it will use GNU Parallel for parallel processing.
4. *Sort and update output:* Sort the generated output file, replacing it with the sorted version.
5. *Optional cleanup:* Compress large log files, and remove empty log files.

**Important note:**
- The [`execute_calculate_scaling_factor_alpha.sh`](https://github.com/kalavattam/protocol_chipseq_signal_norm/blob/main/scripts/execute_calculate_scaling_factor_alpha.sh) script in this code chunk requires that *S. cerevisiae* IP BAM files follow a specific naming convention as described in the protocol manuscript. The expected filename format:
    ```txt
    assay_genotype_state_treatment_factor_strain/replicate.
    ```
    + Required name components:
        - *assay:* Must be 'IP' or 'in', and must be followed by an underscore.
        - *factor:* A required component preceded by an underscore.
        - *strain/replicate:* A required component preceded by an underscore; it marks the end of the pattern.
    + Optional name components:
        - *genotype:* If present, must be preceded by an underscore.
        - *state:* An optional component with preferred values (e.g., 'G1', 'G2M', 'log', or 'Q') but can also be flexible; if present, it must be preceded by an underscore.
        - *treatment:* If present, must be preceded by an underscore.
- Failure to adhere to this naming convention may cause the script to fail.
</details>
<br />

<details>
<summary><i>Bash code: Compute coverage with the sans spike-in quantitative ChIP-seq (siQ-ChIP) method.</i></summary>

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

dir_aln="${dir_pro}/align_reads/${details}"
dir_bam="${dir_aln}/sc"
dir_cvg="${dir_pro}/compute_coverage/${details}"
dir_nrm="${dir_cvg}/norm"
dir_tbl="${dir_cvg}/tables"
dir_alf="${dir_cvg}/alpha"

pattern="*.bam"
include="IP*"
exclude="*Brn1*"

env_nam="env_analyze"
execute="execute_compute_coverage_ratio.sh"

day="$(date '+%Y-%m%d')"
exc_log="${dir_alf}/logs/${day}.$(basename "${execute}" ".sh")"

typ_cvg="alpha"
tbl_alf="${dir_tbl}/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_${typ_cvg}.tsv"

#  Set hardcoded argument assignments
# shellcheck disable=SC1091
source "${dir_fnc}/extract_fld_str.sh"

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
typ_out="bdg.gz"
scl_fct="$(extract_fld_str  3 "${tbl_alf}")"
dep_min="$(extract_fld_str 24 "${tbl_alf}")"  ## WARNING: See above note ##
err_out="${dir_alf}/logs"
nam_job="compute_coverage_ratio_${typ_cvg}"

if ${debug}; then
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
    echo ""
    echo "pattern=${pattern}"
    echo "include=${include}"
    echo "exclude=${exclude}"
    echo ""
    echo "env_nam=${env_nam}"
    echo "execute=${execute}"
    echo ""
    echo "day=${day}"
    echo "exc_log=${exc_log}"
    echo ""
    echo "typ_cvg=${typ_cvg}"
    echo "tbl_alf=${tbl_alf}"
    echo ""
    echo "fil_ip=${fil_ip}"
    echo "fil_in=${fil_in}"
    echo "dir_out=${dir_out}"
    echo "typ_out=${typ_out}"
    echo "scl_fct=${scl_fct}"
    echo "dep_min=${dep_min}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo ""
    echo ""
fi


#  Create required directories if necessary -----------------------------------
mkdir -p ${dir_alf}/{dm_y,dm_n,logs}

#  Debug outdirectory paths
if ${debug}; then
    echo "#####################################"
    echo "## Outdirectory paths and contents ##"
    echo "#####################################"
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
    ls -lhaFG "${dir_alf}/logs"
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
check_program_path parallel
check_program_path python
check_program_path samtools
check_program_path sbatch ||
    echo_warning \
        "SLURM is not available on this system. Do not use the '--slurm'" \
        "flag with the driver script."


#  Generate siQ-scaled coverage tracks ----------------------------------------
if ${debug}; then
    echo "############################################"
    echo "## Call to driver script with '--dep_min' ##"
    echo "############################################"
    echo ""
    echo "bash ${dir_scr}/${execute} \\"
    echo "    --verbose \\"
    echo "    --fil_ip ${fil_ip} \\"
    echo "    --fil_in ${fil_in} \\"
    echo "    --dir_out ${dir_alf}/dm_y \\"
    echo "    --typ_out ${typ_out} \\"
    echo "    --track \\"
    echo "    --scl_fct ${scl_fct} \\"
    echo "    --dep_min ${dep_min} \\"
    echo "    --err_out ${dir_alf}/logs \\"
    echo "    --slurm \\"
    echo "         >> >(tee -a ${exc_log}.stdout.txt) \\"
    echo "        2>> >(tee -a ${exc_log}.stderr.txt)"
    echo ""
    echo ""
fi

#  Run driver script to generate tracks with 'dep_min'
bash "${dir_scr}/${execute}" \
    --verbose \
    --fil_ip "${fil_ip}" \
    --fil_in "${fil_in}" \
    --dir_out "${dir_alf}/dm_y" \
    --typ_out "${typ_out}" \
    --track \
    --scl_fct "${scl_fct}" \
    --dep_min "${dep_min}" \
    --err_out "${dir_alf}/logs" \
    --slurm \
         >> >(tee -a "${exc_log}.stdout.txt") \
        2>> >(tee -a "${exc_log}.stderr.txt")

#  Run driver script to generate tracks sans 'dep_min'
bash "${dir_scr}/${execute}" \
    --verbose \
    --fil_ip "${fil_ip}" \
    --fil_in "${fil_in}" \
    --dir_out "${dir_alf}/dm_n" \
    --typ_out "${typ_out}" \
    --track \
    --scl_fct "${scl_fct}" \
    --err_out "${dir_alf}/logs" \
    --slurm \
         >> >(tee -a "${exc_log}.stdout.txt") \
        2>> >(tee -a "${exc_log}.stderr.txt")

#  Cleanup: Compress logs and remove empty files ------------------------------
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_alf}/logs"
# ls -lhaFG "${dir_alf}/logs"  ## Uncomment to check directory for track logs ##
```
</details>
<br />

<a id="i-compute-coverage-with-the-spike-in-method"></a>
### I. Compute coverage with the spike-in method.
<details>
<summary><i>Text: Compute coverage with the spike-in method.</i></summary>
<br />

This section describes the steps to calculate spike-in normalized ChIP-seq coverage. `#TODO`
</details>
<br />

<details>
<summary><i>Bash code: Compute coverage with the spike-in method.</i></summary>

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

dir_aln="${dir_pro}/align_reads/${details}"
dir_bam="${dir_aln}/sc"
dir_cvg="${dir_pro}/compute_coverage/${details}"
dir_nrm="${dir_cvg}/norm"
dir_tbl="${dir_cvg}/tables"
dir_spk="${dir_cvg}/spike"

env_nam="env_analyze"
execute="execute_compute_coverage_ratio.sh"

day="$(date '+%Y-%m%d')"
exc_log="${dir_spk}/logs/${day}.$(basename "${execute}" ".sh")"

typ_cvg="spike"
tbl_spk="${dir_tbl}/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_${typ_cvg}.tsv"

#  Set hardcoded argument assignments
# shellcheck disable=SC1091
source "${dir_fnc}/extract_fld_str.sh"

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
typ_out="bdg.gz"
scl_fct="$(extract_fld_str  5 "${tbl_spk}")"
dep_min="$(extract_fld_str 21 "${tbl_spk}")"  ## WARNING: See above note ##
err_out="${dir_spk}/logs"
nam_job="compute_coverage_ratio_${typ_cvg}"

if ${debug}; then
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
    echo ""
    echo "pattern=${pattern}"
    echo "include=${include}"
    echo "exclude=${exclude}"
    echo ""
    echo "env_nam=${env_nam}"
    echo "execute=${execute}"
    echo ""
    echo "day=${day}"
    echo "exc_log=${exc_log}"
    echo ""
    echo "typ_cvg=${typ_cvg}"
    echo "tbl_spk=${tbl_spk}"
    echo ""
    echo "fil_ip=${fil_ip}"
    echo "fil_in=${fil_in}"
    echo "dir_out=${dir_out}"
    echo "typ_out=${typ_out}"
    echo "scl_fct=${scl_fct}"
    echo "dep_min=${dep_min}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo ""
    echo ""
fi


#  Create required directories if necessary -----------------------------------
mkdir -p ${dir_spk}/{dm_y,dm_n,logs}

#  Debug outdirectory paths
if ${debug}; then
    echo "#####################################"
    echo "## Outdirectory paths and contents ##"
    echo "#####################################"
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
    ls -lhaFG "${dir_spk}/logs/"
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
check_program_path parallel
check_program_path python
check_program_path samtools
check_program_path sbatch ||
    echo_warning \
        "SLURM is not available on this system. Do not use the '--slurm'" \
        "flag with the driver script."


#  Generate siQ-scaled coverage tracks ----------------------------------------
if ${debug}; then
    echo "############################################"
    echo "## Call to driver script with '--dep_min' ##"
    echo "############################################"
    echo ""
    echo "bash ${dir_scr}/${execute} \\"
    echo "    --verbose \\"
    echo "    --fil_ip ${fil_ip} \\"
    echo "    --fil_in ${fil_in} \\"
    echo "    --dir_out ${dir_spk}/dm_y \\"
    echo "    --typ_out ${typ_out} \\"
    echo "    --track \\"
    echo "    --scl_fct ${scl_fct} \\"
    echo "    --dep_min ${dep_min} \\"
    echo "    --err_out ${dir_spk}/logs \\"
    echo "    --slurm \\"
    echo "         >> >(tee -a ${exc_log}.stdout.txt) \\"
    echo "        2>> >(tee -a ${exc_log}.stderr.txt)"
    echo ""
    echo ""
fi

#  Run driver script to generate tracks with 'dep_min'
bash "${dir_scr}/${execute}" \
    --verbose \
    --fil_ip "${fil_ip}" \
    --fil_in "${fil_in}" \
    --dir_out "${dir_spk}/dm_y" \
    --typ_out "${typ_out}" \
    --track \
    --scl_fct "${scl_fct}" \
    --dep_min "${dep_min}" \
    --err_out "${dir_spk}/logs" \
    --slurm \
         >> >(tee -a "${exc_log}.stdout.txt") \
        2>> >(tee -a "${exc_log}.stderr.txt")

#  Run driver script to generate tracks sans 'dep_min'
bash "${dir_scr}/${execute}" \
    --verbose \
    --fil_ip "${fil_ip}" \
    --fil_in "${fil_in}" \
    --dir_out "${dir_spk}/dm_n" \
    --typ_out "${typ_out}" \
    --track \
    --scl_fct "${scl_fct}" \
    --err_out "${dir_spk}/logs" \
    --slurm \
         >> >(tee -a "${exc_log}.stdout.txt") \
        2>> >(tee -a "${exc_log}.stderr.txt")

#  Cleanup: Compress logs and remove empty files ------------------------------
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_spk}/logs"
# ls -lhaFG "${dir_spk}/logs"  ## Uncomment to check directory for track logs ##
```
</details>
<br />

<a id="j-compute-log2-ratios-of-ip-to-input-coverage"></a>
### J. Compute log2 ratios of IP to input coverage.
<details>
<summary><i>Text: Compute log2 ratios of IP to input coverage.</i></summary>
<br />

`#TODO`
</details>
<br />

<details>
<summary><i>Bash code: Compute log2 ratios of IP to input coverage.</i></summary>

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

dir_aln="${dir_pro}/align_reads/${details}"
dir_bam="${dir_aln}/sc"
dir_cvg="${dir_pro}/compute_coverage/${details}"
dir_nrm="${dir_cvg}/norm"
dir_tbl="${dir_cvg}/tables"
dir_lg2="${dir_cvg}/log2"

env_nam="env_analyze"
execute="execute_compute_coverage_ratio.sh"

day="$(date '+%Y-%m%d')"
exc_log="${dir_lg2}/logs/${day}.$(basename "${execute}" ".sh")"

typ_cvg="log2"
tbl_lg2="${dir_tbl}/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_alpha.tsv"

#  Set hardcoded argument assignments
# shellcheck disable=SC1091
source "${dir_fnc}/extract_fld_str.sh"

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
typ_out="bdg.gz"
dep_min="$(extract_fld_str 24 "${tbl_lg2}")"  ## WARNING: See above note ##
err_out="${dir_lg2}/logs"
nam_job="compute_coverage_ratio_${typ_cvg}"

if ${debug}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
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
    echo ""
    echo "env_nam=${env_nam}"
    echo "execute=${execute}"
    echo ""
    echo "day=${day}"
    echo "exc_log=${exc_log}"
    echo ""
    echo "typ_cvg=${typ_cvg}"
    echo "tbl_lg2=${tbl_lg2}"
    echo ""
    echo "fil_ip=${fil_ip}"
    echo "fil_in=${fil_in}"
    echo "dir_out=${dir_out}"
    echo "typ_out=${typ_out}"
    echo "dep_min=${dep_min}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo ""
    echo ""
fi


#  Create required directories if necessary -----------------------------------
# shellcheck disable=SC2086
mkdir -p ${dir_lg2}/{dm_y,dm_n,logs}

#  Debug outdirectory paths
if ${debug}; then
    echo "#####################################"
    echo "## Outdirectory paths and contents ##"
    echo "#####################################"
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
    ls -lhaFG "${dir_lg2}/logs/"
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
check_program_path parallel
check_program_path python
check_program_path samtools
check_program_path sbatch ||
    echo_warning \
        "SLURM is not available on this system. Do not use the '--slurm'" \
        "flag with the driver script."


#  Generate siQ-scaled coverage tracks ----------------------------------------
if ${debug}; then
    echo "############################################"
    echo "## Call to driver script with '--dep_min' ##"
    echo "############################################"
    echo ""
    echo "bash ${dir_scr}/${execute} \\"
    echo "    --verbose \\"
    echo "    --fil_ip ${fil_ip} \\"
    echo "    --fil_in ${fil_in} \\"
    echo "    --dir_out ${dir_lg2}/dm_y \\"
    echo "    --typ_out ${typ_out} \\"
    echo "    --track \\"
    echo "    --dep_min ${dep_min} \\"
    echo "    --log2 \\"
    echo "    --err_out ${dir_lg2}/logs \\"
    echo "    --slurm \\"
    echo "         >> >(tee -a ${exc_log}.stdout.txt) \\"
    echo "        2>> >(tee -a ${exc_log}.stderr.txt)"
    echo ""
    echo ""
fi

#  Generate tracks using 'dep_min'
bash "${dir_scr}/${execute}" \
    --verbose \
    --fil_ip "${fil_ip}" \
    --fil_in "${fil_in}" \
    --dir_out "${dir_lg2}/dm_y" \
    --typ_out "${typ_out}" \
    --track \
    --dep_min "${dep_min}" \
    --log2 \
    --err_out "${dir_lg2}/logs" \
    --slurm \
         >> >(tee -a "${exc_log}.stdout.txt") \
        2>> >(tee -a "${exc_log}.stderr.txt")

#  Generate tracks sans 'dep_min'
bash "${dir_scr}/${execute}" \
    --verbose \
    --fil_ip "${fil_ip}" \
    --fil_in "${fil_in}" \
    --dir_out "${dir_lg2}/dm_n" \
    --typ_out "${typ_out}" \
    --track \
    --log2 \
    --err_out "${dir_lg2}/logs" \
    --slurm \
         >> >(tee -a "${exc_log}.stdout.txt") \
        2>> >(tee -a "${exc_log}.stderr.txt")


#  Cleanup: Compress logs and remove empty files ------------------------------
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_lg2}/logs"

# ls -lhaFG "${dir_lg2}/tracks/logs"  ## Uncomment to check directory for track logs ##
# find . -maxdepth 1 -type f -size 0 -delete -o -type f -not -name "*.gz" -exec gzip {} +
```
</details>
<br />

<details>
<summary><i>Code: Raw</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define functions -----------------------------------------------------------
#  Function to write a compressed bedGraph file
function write_bdg_gz() {
    local fil_in="${1}"     # Input data file (e.g., 'IP.data' or 'IN.data')
    local fil_out="${2}"    # Output compressed bedGraph file: '*.bdg.gz'
    local fct_nrm="${3:-}"  # Norm. factor: Empty or 1 for no norm.

    #  Determine whether to read file with zcat or cat
    if [[ "${fil_in}" == *.gz ]]; then
        reader="zcat"
    else
        reader="cat"
    fi

    #  Use awk to process the file
    # shellcheck disable=SC2002
    ${reader} "${fil_in}" \
        | awk \
            -v OFS="\t" -v nl="${fct_nrm}" '
            {
                if (nl != "") { print $1, $2, $3, $4 / nl }
                else { print $1, $2, $3, $4 }
            }
        '  \
        | gzip \
            > "${fil_out}"
}


#  Function to print error for writing compressed bedGraph coverage files
function print_err_cvg() {
    local type="${1}"
    echo \
        "Error with no exit: Encountered issue writing *.bdg.gz file of" \
        "${type} coverage." >&2
}


#  Helper function to write compressed bedGraph file, handle errors, and
#+ optionally clean up data infile
function write_check_bdg() {
    local fil_src=""
    local fil_out=""
    local typ_cvg=""
    local n_lines=""
    local rmv_src=false

    #  Parse keyword parameters
    while [[ $# -gt 0 ]]; do
        case "${1}" in
             -s|--fil_src) fil_src="${2}"; shift 2 ;;
             -o|--fil_out) fil_out="${2}"; shift 2 ;;
            -tc|--typ_cvg) typ_cvg="${2}"; shift 2 ;;
            -nl|--n_lines) n_lines="${2}"; shift 2 ;;
            -rs|--rmv_src) rmv_src=true;   shift 1 ;;
            *)
                echo "## Unknown parameter passed: ${1} ##" >&2
                return 1
                ;;
        esac
    done

    #  Validate required parameters
    check_supplied_arg "fil_src" "${fil_src}"
    check_exists_file  "fil_src" "${fil_src}"

    check_supplied_arg "fil_out" "${fil_out}"

    check_supplied_arg "typ_cvg" "${typ_cvg}"

    #  If supplied, validate denominator for normalization 
    if [[ -n "${n_lines}" ]]; then check_int_nonneg "n_lines" "${n_lines}"; fi

    #  Write compressed bedGraph file
    # shellcheck disable=SC2086
    if ! write_bdg_gz "${fil_src}" "${fil_out}" ${n_lines}; then
        print_err_cvg "${typ_cvg}"
        return 1
    fi

    #  Optionally remove source file if compressed bedGraph file was
    #+ successfully written
    if ${rmv_src} && [[ -f "${fil_out}" ]]; then
        rm "${fil_src}" || {
            echo \
                "Warning: --rmv_src was specified but could not remove" \
                "source file '${fil_src}'." >&2
        }
    fi

    return 0
}


function check_unity() {
    local fil_in="${1}"
    local rng_gt="${2:-0.98}"
    local rng_lt="${3:-1.02}"
    local reader

    #  Determine whether to read file with zcat or cat
    if [[ "${fil_in}" == *.gz ]]; then
        reader="zcat"
    else
        reader="cat"
    fi

    #  Process the file and check its sum
    ${reader} "${fil_in}" \
        | awk -v rng_gt="${rng_gt}" -v rng_lt="${rng_lt}" '{
            sum += $4
        } END {
            if (sum >= rng_gt && sum <= rng_lt) {
                print "File sums to approximately unity:", sum
            } else {
                print "File does not sum to unity:", sum
            }
        }'
}
```

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
debug=true

#  Define base directory for repository
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##

#  Define paths to protocol repository and its subdirectories
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"
dir_dat="${dir_rep}/data"
dir_raw="${dir_dat}/raw"
dir_pro="${dir_dat}/processed"

#  Define alignment and coverage details
aligner="bowtie2"
a_type="global"
req_flg=true
flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
mapq=1
det_bam="flag-${flg}_mapq-${mapq}"
det_cvg="${aligner}_${a_type}_${det_bam}"

#  Further define directory setup
dir_aln="${dir_pro}/align_${aligner}_${a_type}"
dir_bam="${dir_aln}/${det_bam}/sc"
dir_cvg="${dir_pro}/compute_coverage"
dir_non="${dir_cvg}/${det_cvg}/raw"
dir_nrm="${dir_cvg}/${det_cvg}/norm"
dir_alf="${dir_cvg}/${det_cvg}/alpha"

#  Define environment, resources, and script arguments 
env_nam="env_analyze"
threads=8
tbl_mes="${dir_raw}/docs/measurements_siqchip.tsv"
tbl_col="alpha"
eqn="6nd"
typ_out="bedgraph"
siz_bin=10  # 1

#  Define file search parameters  ## WARNING: Change as needed ##
pattern="*.bam"
include="IP*"
exclude="*Brn1*"
infiles="$(
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_bam}" \
        --pattern "${pattern}" \
        --include "${include}" \
        --exclude "${exclude}"
)"

#  Define scripts and output files
scr_tbl="execute_calculate_scaling_${tbl_col}.sh"
scr_trk="execute_deeptools_coverage.sh"
fil_tbl="${dir_alf}/tables/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_${tbl_col}-${eqn}.tsv"

#  Define information for logging
day="$(date '+%Y-%m%d')"

#  Debug variable assignments
if ${debug}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "\${debug}=${debug}"
    echo ""
    echo "\${dir_bas}=${dir_bas}"
    echo ""
    echo "\${dir_rep}=${dir_rep}"
    echo "\${dir_scr}=${dir_scr}"
    echo "\${dir_fnc}=${dir_fnc}"
    echo "\${dir_dat}=${dir_dat}"
    echo "\${dir_raw}=${dir_raw}"
    echo "\${dir_pro}=${dir_pro}"
    echo ""
    echo "\${aligner}=${aligner}"
    echo "\${a_type}=${a_type}"
    echo "\${req_flg}=${req_flg}"
    echo "\${flg}=${flg}"
    echo "\${mapq}=${mapq}"
    echo "\${det_bam}=${det_bam}"
    echo "\${det_cvg}=${det_cvg}"
    echo ""
    echo "\${dir_aln}=${dir_aln}"
    echo "\${dir_bam}=${dir_bam}"
    echo "\${dir_cvg}=${dir_cvg}"
    echo "\${dir_non}=${dir_non}"
    echo "\${dir_nrm}=${dir_nrm}"
    echo "\${dir_alf}=${dir_alf}"
    echo ""
    echo "\${env_nam}=${env_nam}"
    echo "\${threads}=${threads}"
    echo "\${tbl_mes}=${tbl_mes}"
    echo "\${tbl_col}=${tbl_col}"
    echo "\${eqn}=${eqn}"
    echo "\${typ_out}=${typ_out}"
    echo "\${siz_bin}=${siz_bin}"
    echo ""
    echo "\${pattern}=${pattern}"
    echo "\${include}=${include}"
    echo "\${exclude}=${exclude}"
    echo "\${infiles}=${infiles}"
    echo ""
    echo "\${scr_tbl}=${scr_tbl}"
    echo "\${scr_trk}=${scr_trk}"
    echo "\${fil_tbl}=${fil_tbl}"
    echo ""
    echo "\${day}=${day}"
    echo ""
    echo ""
fi


#  Create required directories if necessary -----------------------------------
mkdir -p {${dir_non},${dir_nrm}}/tracks/{docs,logs}
mkdir -p ${dir_alf}/{tables,tracks}/{docs,logs}

#  Debug outdirectory paths
if ${debug}; then
    echo "#####################################"
    echo "## Outdirectory paths and contents ##"
    echo "#####################################"
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%%%%%"
    echo "%% \${dir_non}/tracks %%"
    echo "%%%%%%%%%%%%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG ${dir_non}/tracks
    echo ""
    ls -lhaFG ${dir_non}/tracks/*
    echo ""
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%%%%%"
    echo "%% \${dir_nrm}/tracks %%"
    echo "%%%%%%%%%%%%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG ${dir_nrm}/tracks
    echo ""
    ls -lhaFG ${dir_nrm}/tracks/*
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%%%%%"
    echo "%% \${dir_alf}/tables %%"
    echo "%%%%%%%%%%%%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG ${dir_alf}/tables
    echo ""
    ls -lhaFG ${dir_alf}/tables/*
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%%%%%"
    echo "%% \${dir_alf}/tracks %%"
    echo "%%%%%%%%%%%%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG ${dir_alf}/tracks
    echo ""
    ls -lhaFG ${dir_alf}/tracks/*
    echo ""
    echo ""
fi


#  Activate the environment and check dependencies ----------------------------
#  Source utility functions
source "${dir_fnc}/check_program_path.sh"
source "${dir_fnc}/echo_warning.sh"
source "${dir_fnc}/handle_env.sh"

#  Activate the required environment
handle_env "${env_nam}"

#  Check the availability of necessary dependencies such as GNU Parallel and
#+ SLURM sbatch
check_program_path awk
check_program_path parallel
check_program_path python
check_program_path samtools
check_program_path sbatch ||
    echo_warning \
        "SLURM is not available on this system. Do not use the '--slurm'" \
        "flag with the driver script."


#  Generate raw, non-normalized signal tracks ---------------------------------
bash "${dir_scr}/execute_deeptools_coverage.sh" \
    --verbose \
    --threads "${threads}" \
    --infiles "${infiles}" \
    --dir_out "${dir_non}/tracks" \
    --typ_out "${typ_out}" \
    --typ_cvg "raw" \
    --siz_bin "${siz_bin}" \
    --err_out "${dir_non}/tracks/logs" \
    --nam_job "raw" \
    --slurm \
         >> >(tee -a "${dir_non}/tracks/logs/${day}_raw.stdout.txt") \
        2>> >(tee -a "${dir_non}/tracks/logs/${day}_raw.stderr.txt")

samtools view \
    -@ ${SLURM_CPUS_ON_NODE} -c -f 64 \
    "${dir_bam}/IP_WT_G1_Hho1_6336.sc.bam"
# 13492920

wc -l "${dir_aln}/${det_bam}/sc_bed/IP_WT_G1_Hho1_6336.sc.bed" \
    | awk '{ print $1 }'
# 13492920

bam="${dir_bam}/IP_WT_G1_Hho1_6336.sc.bam"
# bdg="${bam/.bam/.norm_len}"
bdg="${bam/.bam/.py_norm_len_unity}"
python "${dir_scr}/compute_coverage.py" \
    --verbose \
    --threads ${SLURM_CPUS_ON_NODE} \
    --infile "${bam}" \
    --outfile "${bdg}" \
    --typ_out bedgraph \
    --bin_siz 10 \
    --typ_cvg norm

write_bdg_gz \
    "${bdg}.bdg.gz" \
    "${bdg}_unity.bdg.gz" \
    13492920

check_unity "${bdg}_unity.bdg.gz"

cd ~/repos/protocol_chipseq_signal_norm/data/processed/compute_coverage/bowtie2_global_flag-2_mapq-1/norm/tracks

for bdg in *.bdg.gz; do
    echo "## ${bdg} ##"
    check_unity "${bdg}"
    echo ""
done
```
</details>
<br />

<details>
<summary><i>Bash code: Compute alpha-scaled coverage using the TSV file.</i></summary>

`## #INPROGRESS Draft code ##`
```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##


#  Define variables -----------------------------------------------------------
nam_job="comp_covg_spike"
no_infiles="$(tail -n +2 "${outfile}" | wc -l)"
max_job=6
threads=8
siz_bin=1
outtype="bigwig"
time="0:30:00"
dir_out="${dir_pro}/compute_coverage/${det_cvg}/spike/tracks"
err_out="${dir_out}/logs"
execute="${err_out}/${day}.execute.${nam_job}"

#  Loop over each line (skipping the header) to extract 'sample' and 'scaled' columns
while IFS=$'\t' read -r sample sf scaled main_ip spike_ip main_in spike_in; do
    #  Extract the base name (without directory path) for use as outfile prefix
    outstem="${sample%.bam}"

    echo "sample    ${sample}"
    echo "outstem   ${outstem}"
    echo "sf        ${sf}"
    echo "scaled    ${scaled}"
    echo "main_ip   ${main_ip}"
    echo "spike_ip  ${spike_ip}"
    echo "main_in   ${main_in}"
    echo "spike_in  ${spike_in}"
    echo ""

    ls -lhaFG "${dir_bam}/${sample}"
    echo ""

    #  Run the Python script with parsed values
    cat << EOM
        srun \\
            --job-name=${nam_job} \\
            --nodes=1 \\
            --cpus-per-task=${threads} \\
            --time=${time} \\
            --error=${err_out}/${nam_job}.%A-%a.stderr.txt \\
            --output=${err_out}/${nam_job}.%A-%a.stdout.txt \\
                python "${dir_scr}/compute_coverage.py" \\
                    --infile "${dir_bam}/${sample}" \\
                    --outfile "${dir_cvg}/${det_cvg}/spike/tracks/${outstem}" \\
                    --scl_fct "${scaled}" \\
                    --outtype "${outtype}" \\
                    --threads ${threads} \\
                    --siz_bin 1 \\
                         > >(tee -a "${execute}.${outstem}.stdout.txt") \\
                        2> >(tee -a "${execute}.${outstem}.stderr.txt")
EOM
    echo ""
    echo ""

    srun \
        --job-name=${nam_job} \
        --nodes=1 \
        --cpus-per-task=${threads} \
        --time=${time} \
        --error=${err_out}/${nam_job}.%A-%a.stderr.txt \
        --output=${err_out}/${nam_job}.%A-%a.stdout.txt \
            python "${dir_scr}/compute_coverage.py" \
                --infile "${dir_bam}/${sample}" \
                --outfile "${dir_cvg}/${det_cvg}/spike/tracks/${outstem}" \
                --scl_fct "${scaled}" \
                --outtype "${outtype}" \
                --threads ${threads} \
                --siz_bin 1 \
                     > >(tee -a "${execute}.${outstem}.stdout.txt") \
                    2> >(tee -a "${execute}.${outstem}.stderr.txt")

    sleep 0.3
done < <(tail -n +2 "${outfile}")

# --array=1-${no_infiles}%${max_job} \\


#TODO ...
#  Cleanup: Compress large stdout and stderr files, and remove files with size
#+ 0
bash "${dir_scr}/compress_remove_files.sh" \
    --dir_fnd "${err_out}"

#  Optional: Check the contents of the logs directory
# ls -lhaFG "${err_out}"
```
</details>
<br />
