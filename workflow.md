
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
    1. [F. Compute normalized coverage.](#f-compute-normalized-coverage)
    1. [G. Compute log2 ratios of IP to input coverage.](#g-compute-log2-ratios-of-ip-to-input-coverage)
    1. [G. Compute coverage with the sans spike-in quantitative ChIP-seq \(siQ-ChIP\) method.](#g-compute-coverage-with-the-sans-spike-in-quantitative-chip-seq-siq-chip-method)
    1. [H. Compute coverage using the spike-in method.](#h-compute-coverage-using-the-spike-in-method)

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

<a id="f-compute-normalized-coverage"></a>
### F. Compute normalized coverage.
<details>
<summary><i>Text: Compute normalized raw coverage..</i></summary>
<br />

This following Bash code chunk provides an example of how to compute normalized ChIP-seq coverage per [Dickson et al., *Sci Rep*, 2023](https://www.nature.com/articles/s41598-023-34430-2). The coverage type is determined by setting the variable `typ_cvg` to `"raw"` or `"norm"`. BIGWIG and log output files will be saved to separate directories based on the selected coverage type.
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

#  Define directory paths
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"
dir_dat="${dir_rep}/data"
dir_pro="${dir_dat}/processed"

#  Define alignment and coverage details
aligner="bowtie2"
a_type="global"
req_flg=true
flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
mapq=1
det_bam="flag-${flg}_mapq-${mapq}"
det_cvg="${aligner}_${a_type}_${det_bam}"
typ_cvg="norm"  ## WARNING: "norm" for normalized ##

#  Further define directory setup
dir_aln="${dir_pro}/align_${aligner}_${a_type}"
dir_bam="${dir_aln}/${det_bam}/sc"
dir_cvg="${dir_pro}/compute_coverage"
dir_trk="${dir_cvg}/${det_cvg}/${typ_cvg}/tracks"

#  Define driver script
exc_cvg="${dir_scr}/execute_compute_coverage.sh"

#  Define script arguments, environment, and resources
nam_job="compute_coverage_${typ_cvg}"
typ_out="bdg.gz"
threads=8
siz_bin=10
env_nam="env_analyze"
day="$(date '+%Y-%m%d')"
err_out="${dir_trk}/logs"
exc_pth="${dir_trk}/logs/${day}.$(basename "${exc_cvg}" ".sh")_${typ_cvg}"

#  Define file search parameters
## WARNING: Change search parameters as needed ##
pattern="*.bam"
exclude="*Brn1*"
infiles="$(
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_bam}" \
        --pattern "${pattern}" \
        --exclude "${exclude}"
)"


#  Create required directories if necessary -----------------------------------
mkdir -p echo ${dir_cvg}/${det_cvg}/norm/tracks/{docs,logs}


#  Activate the environment and check dependencies ----------------------------
#  Source utility functions
source "${dir_fnc}/check_program_path.sh"
source "${dir_fnc}/check_unity.sh"
source "${dir_fnc}/echo_warning.sh"
source "${dir_fnc}/handle_env.sh"

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
    echo "bash ${exc_cvg} \\"
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
    echo "         >> >(tee -a ${exc_pth}.stdout.txt) \\"
    echo "        2>> >(tee -a ${exc_pth}.stderr.txt)"
    echo ""
    echo ""
fi

bash "${exc_cvg}" \
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
         >> >(tee -a "${exc_pth}.stdout.txt") \
        2>> >(tee -a "${exc_pth}.stderr.txt")

#  Check that each normalized coverage BEDGRAPH file sums to unity
if ${debug}; then
    echo "##############################################"
    echo "## Check unity of normalized coverage files ##"
    echo "##############################################"
    echo ""
    for bdg in ${dir_trk}/*.${typ_out}; do
        echo "$(basename "${bdg}")"
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

<a id="g-compute-log2-ratios-of-ip-to-input-coverage"></a>
### G. Compute log2 ratios of IP to input coverage.
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

#  Define paths to protocol repository and its subdirectories
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"
dir_dat="${dir_rep}/data"
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
dir_cvg="${dir_pro}/compute_coverage"
dir_nrm="${dir_cvg}/${det_cvg}/norm"
dir_alf="${dir_cvg}/${det_cvg}/alpha"
dir_lg2="${dir_cvg}/${det_cvg}/log2"

#  Define environment, resources, and script arguments 
env_nam="env_analyze"
threads=8
typ_out="bdg.gz"
siz_bin=10

#  Define file search parameters  ## WARNING: Change as needed ##
pattern="${typ_out}"
exclude="*Brn1*"
ser_ip="$(
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_nrm}/tracks" \
        --pattern "*.${typ_out}" \
        --include "IP*" \
        --exclude "${exclude}"
)"
ser_in="$(
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_nrm}/tracks" \
        --pattern "*.${typ_out}" \
        --include "in*" \
        --exclude "${exclude}"
)"

#  Define script and table of minimum input depth values
scr_trk="execute_compute_coverage_ratio.sh"
tbl_min="${dir_alf}/tables/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_depth_min_${siz_bin}.tsv"

#  Define log file prefixes
day="$(date '+%Y-%m%d')"
exc_trk="${dir_lg2}/tracks/logs/${day}.${scr_trk%.sh}"

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
    echo "\${dir_cvg}=${dir_cvg}"
    echo "\${dir_nrm}=${dir_nrm}"
    echo "\${dir_alf}=${dir_alf}"
    echo "\${dir_lg2}=${dir_lg2}"
    echo ""
    echo "\${env_nam}=${env_nam}"
    echo "\${threads}=${threads}"
    echo "\${typ_out}=${typ_out}"
    echo "\${siz_bin}=${siz_bin}"
    echo ""
    echo "\${pattern}=${pattern}"
    echo "\${exclude}=${exclude}"
    echo "\${ser_ip}=${ser_ip}"
    echo "\${ser_in}=${ser_in}"
    echo ""
    echo "\${scr_trk}=${scr_trk}"
    echo "\${tbl_min}=${tbl_min}"
    echo ""
    echo "\${day}=${day}"
    echo "\${exc_trk}=${exc_trk}"
    echo ""
    echo ""
fi


#  Create required directories if necessary -----------------------------------
# shellcheck disable=SC2086
mkdir -p ${dir_lg2}/tracks/{docs,logs}

#  Debug outdirectory paths
if ${debug}; then
    echo "#####################################"
    echo "## Outdirectory paths and contents ##"
    echo "#####################################"
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%%%%"
    echo "%% \${dir_nrm}/tracks %%"
    echo "%%%%%%%%%%%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_nrm}/tracks"
    echo ""
    ls -lhaFG "${dir_nrm}/tracks/logs"
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%%%%"
    echo "%% \${dir_alf}/tables %%"
    echo "%%%%%%%%%%%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_alf}/tables"
    echo ""
    ls -lhaFG "${dir_alf}/tables/logs"
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%%%%"
    echo "%% \${dir_lg2}/tracks %%"
    echo "%%%%%%%%%%%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG "${dir_lg2}/tracks"
    echo ""
    ls -lhaFG "${dir_lg2}/tracks/logs"
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
#  Parse table, assigning serialized string to variables 
dep_min="$(
    awk 'NR > 1 { flt = flt ? flt "," $3 : $3 } END { print flt }' "${tbl_min}"
)"

if ${debug}; then
    echo "###############################################"
    echo "## Call to execute_compute_coverage_ratio.sh ##"
    echo "###############################################"
    echo ""
    echo "bash ${dir_scr}/${scr_trk} \\"
    echo "    --verbose \\"
    echo "    --fil_ip ${ser_ip} \\"
    echo "    --fil_in ${ser_in} \\"
    echo "    --dir_out ${dir_lg2}/tracks \\"
    echo "    --typ_out ${typ_out} \\"
    echo "    --dep_min ${dep_min} \\"
    echo "    --log2 \\"
    echo "    --err_out ${dir_lg2}/tracks/logs \\"
    echo "    --slurm \\"
    echo "         >> >(tee -a ${exc_trk}.stdout.txt) \\"
    echo "        2>> >(tee -a ${exc_trk}.stderr.txt)"
    echo ""
    echo ""
fi

bash "${dir_scr}/${scr_trk}" \
    --verbose \
    --fil_ip "${ser_ip}" \
    --fil_in "${ser_in}" \
    --dir_out "${dir_lg2}/tracks" \
    --typ_out "${typ_out}" \
    --dep_min "${dep_min}" \
    --log2 \
    --err_out "${dir_lg2}/tracks/logs" \
    --slurm \
         >> >(tee -a "${exc_trk}.stdout.txt") \
        2>> >(tee -a "${exc_trk}.stderr.txt")


#  Cleanup: Compress logs and remove empty files ------------------------------
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_lg2}/tracks/logs"

# ls -lhaFG "${dir_lg2}/tracks/logs"  ## Uncomment to check directory for track logs ##


###############
## Debugging ##
###############

vi "${dir_scr}/compute_coverage_ratio.py"

{
    python "${dir_scr}/compute_coverage_ratio.py" \
        -v \
        -fp /home/kalavatt/repos/protocol_chipseq_signal_norm/data/processed/compute_coverage/bowtie2_global_flag-2_mapq-1/norm/tracks/IP_WT_G1_Hho1_6337.sc.bdg.gz \
        -fn /home/kalavatt/repos/protocol_chipseq_signal_norm/data/processed/compute_coverage/bowtie2_global_flag-2_mapq-1/norm/tracks/in_WT_G1_Hho1_6337.sc.bdg.gz \
        -fo /home/kalavatt/repos/protocol_chipseq_signal_norm/data/processed/compute_coverage/bowtie2_global_flag-2_mapq-1/log2/tracks/dm_yes.log2_rat_WT_G1_Hho1_6337.sc.bdg.gz \
        -tr \
        -dm 0.000000822564930190970796 \
        -l2

    python "${dir_scr}/compute_coverage_ratio.py" \
        -v \
        -fp /home/kalavatt/repos/protocol_chipseq_signal_norm/data/processed/compute_coverage/bowtie2_global_flag-2_mapq-1/norm/tracks/IP_WT_G1_Hho1_6337.sc.bdg.gz \
        -fn /home/kalavatt/repos/protocol_chipseq_signal_norm/data/processed/compute_coverage/bowtie2_global_flag-2_mapq-1/norm/tracks/in_WT_G1_Hho1_6337.sc.bdg.gz \
        -fo /home/kalavatt/repos/protocol_chipseq_signal_norm/data/processed/compute_coverage/bowtie2_global_flag-2_mapq-1/log2/tracks/dm_no.log2_rat_WT_G1_Hho1_6337.sc.bdg.gz \
        -tr \
        -l2
}
```
</details>
<br />

<a id="g-compute-coverage-with-the-sans-spike-in-quantitative-chip-seq-siq-chip-method"></a>
### G. Compute coverage with the sans spike-in quantitative ChIP-seq (siQ-ChIP) method.
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
- The [`execute_calculate_scaling_factor_alpha.sh`](https://github.com/kalavattam/protocol_chipseq_signal_norm/blob/main/scripts/execute_calculate_scaling_factor_alpha.sh) script in this code chunk requires that *S. cerevisiae* IP BAM files follow a specific naming convention as outlined in the accompanying manuscript. The expected filename format:
    ```txt
    assay_genotype_state_treatment_factor_strain/replicate.
    ```
    + Required filename components:
        - *assay:* Must be 'IP' or 'in', and must be followed by an underscore.
        - *factor:* A required component preceded by an underscore.
        - *strain/replicate:* A required component preceded by an underscore; it marks the end of the pattern.
    + Optional filename components:
        - *genotype:* If present, must be preceded by an underscore.
        - *state:* An optional component with preferred values (e.g., 'G1', 'G2M', 'log', or 'Q') but can also be flexible; if present, it must be preceded by an underscore.
        - *treatment:* If present, must be preceded by an underscore.
- Failure to adhere to this naming convention may cause the script to fail.
</details>
<br />

<details>
<summary><i>Bash code: Generate a TSV file of sample-specific siQ-ChIP alpha scaling factors.</i></summary>

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
dir_nrm="${dir_cvg}/${det_cvg}/norm"
dir_alf="${dir_cvg}/${det_cvg}/alpha"

#  Define environment, resources, and script arguments 
env_nam="env_analyze"
threads=8
tbl_mes="${dir_raw}/docs/measurements_siqchip.tsv"
typ_cvg="alpha"
eqn="6nd"
typ_out="bdg.gz"
siz_bin=10

#  Define file search parameters  ## WARNING: Change as needed ##
pattern="*.bam"
exclude="*Brn1*"
ser_ip="$(
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_bam}" \
        --pattern "${pattern}" \
        --include "IP*" \
        --exclude "${exclude}"
)"
ser_in="$(
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_bam}" \
        --pattern "${pattern}" \
        --include "in*" \
        --exclude "${exclude}"
)"

#  Define scripts and output files
scr_tbl="execute_calculate_scaling_factor_${typ_cvg}.sh"
scr_trk="execute_compute_coverage_ratio.sh"
tbl_alf="${dir_alf}/tables/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_${typ_cvg}_${eqn}.tsv"
tbl_min="${dir_alf}/tables/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_depth_min_${siz_bin}.tsv"

#  Define log file prefixes
day="$(date '+%Y-%m%d')"
exc_tbl="${dir_alf}/tables/logs/${day}.${scr_tbl%.sh}_${eqn}"
exc_trk="${dir_alf}/tracks/logs/${day}.${scr_trk%.sh}"

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
    echo "\${dir_nrm}=${dir_nrm}"
    echo "\${dir_alf}=${dir_alf}"
    echo ""
    echo "\${env_nam}=${env_nam}"
    echo "\${threads}=${threads}"
    echo "\${eqn}=${eqn}"
    echo "\${tbl_mes}=${tbl_mes}"
    echo "\${typ_cvg}=${typ_cvg}"
    echo "\${typ_out}=${typ_out}"
    echo "\${siz_bin}=${siz_bin}"
    echo ""
    echo "\${pattern}=${pattern}"
    echo "\${exclude}=${exclude}"
    echo "\${ser_ip}=${ser_ip}"
    echo "\${ser_in}=${ser_in}"
    echo ""
    echo "\${scr_tbl}=${scr_tbl}"
    echo "\${scr_trk}=${scr_trk}"
    echo "\${tbl_min}=${tbl_min}"
    echo "\${tbl_alf}=${tbl_alf}"
    echo ""
    echo "\${day}=${day}"
    echo "\${exc_tbl}=${exc_tbl}"
    echo "\${exc_trk}=${exc_trk}"
    echo ""
    echo ""
fi


#  Create required directories if necessary -----------------------------------
mkdir -p ${dir_alf}/{tables,tracks}/{docs,logs}

#  Debug outdirectory paths
if ${debug}; then
    echo "#####################################"
    echo "## Outdirectory paths and contents ##"
    echo "#####################################"
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%%%%"
    echo "%% \${dir_nrm}/tracks %%"
    echo "%%%%%%%%%%%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG ${dir_nrm}/tracks
    echo ""
    ls -lhaFG ${dir_nrm}/tracks/logs
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%%%%"
    echo "%% \${dir_alf}/tables %%"
    echo "%%%%%%%%%%%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG ${dir_alf}/tables
    echo ""
    ls -lhaFG ${dir_alf}/tables/logs
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%%%%"
    echo "%% \${dir_trk}/tracks %%"
    echo "%%%%%%%%%%%%%%%%%%%%%%%"
    echo ""
    ls -lhaFG ${dir_alf}/tracks
    echo ""
    ls -lhaFG ${dir_alf}/tracks/logs
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


#  Calculate and write out input minimum depths -------------------------------
if [[ ! -f "${tbl_min}" ]]; then
    unset arr_ser_in && typeset -a arr_ser_in
    IFS=',' read -r -a arr_ser_in <<< "${ser_in}"

    #TODO: Parallelize this
    {
        echo -e "file\tfrag\tnorm"
        for bam in "${arr_ser_in[@]}"; do
            frag=$(
                python "${dir_scr}/calculate_factor_depth.py" \
                    -i "${bam}" \
                    -m "frag" \
                    -sb ${siz_bin}
            )
            norm=$(
                python "${dir_scr}/calculate_factor_depth.py" \
                    -i "${bam}" \
                    -m "norm" \
                    -sb ${siz_bin}
            )
            echo -e "$(basename ${bam})\t${frag}\t${norm}"
        done
    } \
        | tee -a "${tbl_min}"

    unset arr_ser_in
fi

if [[ -f "${tbl_min}" ]]; then
    #  Sort the table of input minimume depths by rows
    awk 'NR == 1; NR > 1 && NF { print | "sort" }' "${tbl_min}" \
        > "${dir_alf}/tables/tmp.tsv"

    #  Replace the original table with the sorted version
    mv -f "${dir_alf}/tables/tmp.tsv" "${tbl_min}"
fi

# cat "${tbl_min}"  ## Uncomment to check the table contents ##


#  Calculate siQ-ChIP alpha scaling factors -----------------------------------
#  Debug call to alpha computation driver script, etc.
if ${debug}; then
    echo "###################################################"
    echo "## Call to alpha computation driver script, etc. ##"
    echo "###################################################"
    echo ""
    echo "if [[ ! -f ${tbl_alf} ]]; then"
    echo "    bash ${dir_scr}/${scr_tbl} \\"
    echo "        --verbose \\"
    echo "        --threads ${threads} \\"
    echo "        --ser_ip ${ser_ip} \\"
    echo "        --ser_in ${ser_in} \\"
    echo "        --table ${tbl_mes} \\"
    echo "        --eqn ${eqn} \\"
    echo "        --outfile ${tbl_alf} \\"
    echo "        --err_out ${dir_alf}/tables/logs \\"
    echo "        --flg_dep \\"
    echo "        --flg_len \\"
    echo "        --flg_mc \\"
    echo "        --slurm \\"
    echo "             >> >(tee -a ${exc_tbl}.stdout.txt) \\"
    echo "            2>> >(tee -a ${exc_tbl}.stderr.txt)"
    echo "fi"
    echo ""
    echo "if [[ -f ${tbl_alf} ]]; then"
    echo "    awk 'NR == 1; NR > 1 && NF { print | "sort" }' ${tbl_alf} \\"
    echo "        > ${dir_alf}/tables/tmp.tsv"
    echo ""
    echo "    mv -f ${dir_alf}/tables/tmp.tsv ${tbl_alf}"
    echo "fi"
    echo ""
    echo ""
fi

#  Run the driver script to generate a TSV file of sample-specific siQ-ChIP
#+ alpha scaling factors
if [[ ! -f "${tbl_alf}" ]]; then
    bash "${dir_scr}/${scr_tbl}" \
        --verbose \
        --threads "${threads}" \
        --ser_ip "${ser_ip}" \
        --ser_in "${ser_in}" \
        --table "${tbl_mes}" \
        --eqn "${eqn}" \
        --outfile "${tbl_alf}" \
        --err_out "${dir_alf}/tables/logs" \
        --flg_dep \
        --flg_len \
        --flg_mc \
        --slurm \
             >> >(tee -a "${exc_tbl}.stdout.txt") \
            2>> >(tee -a "${exc_tbl}.stderr.txt")
fi

if [[ -f "${tbl_alf}" ]]; then
    #  Sort the table of scaling factors by rows
    awk 'NR == 1; NR > 1 && NF { print | "sort" }' "${tbl_alf}" \
        > "${dir_alf}/tables/tmp.tsv"

    #  Replace the original table with the sorted version
    mv -f "${dir_alf}/tables/tmp.tsv" "${tbl_alf}"
fi
# cat "${tbl_alf}"  ## Uncomment to check the table contents ##

#  Check that 'tbl_min' column 1 matches 'tbl_alf' column 2 (basenames)
if ${debug}; then
    if ! \
        diff -q \
            <(awk 'NR > 1 { print $1 }' "${tbl_min}") \
            <(
                awk 'NR > 1 { print $2 }' "${tbl_alf}" \
                    | xargs -I {} basename {}
            ) > /dev/null
    then
        echo \
            "Mismatch detected between \${tbl_min} column 1 and \${tbl_alf}" \
            "column 2. Differences:"
        diff <(awk 'NR > 1 { print $1 }' "${tbl_min}") \
             <(
                awk 'NR > 1 { print $2 }' "${tbl_alf}" \
                    | xargs -I {} basename {}
            )
    fi
fi


#  Generate siQ-scaled coverage tracks ----------------------------------------
#  Parse tables, assigning serialized strings to variables 
ser_num=$(
    sed \
        -e "s:${dir_bam}:${dir_nrm}/tracks:g" \
        -e "s:.bam:.${typ_out}:g" \
        < <(awk 'NR > 1 { print $1 }' "${tbl_alf}" | paste -sd ',' -)
)
ser_den=$(
    sed \
        -e "s:${dir_bam}:${dir_nrm}/tracks:g" \
        -e "s:.bam:.${typ_out}:g" \
        < <(awk 'NR > 1 { print $2 }' "${tbl_alf}" | paste -sd ',' -)
)
scl_fct="$(awk 'NR > 1 { print $3 }' "${tbl_alf}" | paste -sd ',' -)"
dep_min="$(awk 'NR > 1 { print $3 }' "${tbl_min}" | paste -sd ',' -)"

if ${debug}; then
    echo "###############################################"
    echo "## Call to execute_compute_coverage_ratio.sh ##"
    echo "###############################################"
    echo ""
    echo "bash ${dir_scr}/${scr_trk} \\"
    echo "    --verbose \\"
    echo "    --fil_ip ${ser_num} \\"
    echo "    --fil_in ${ser_den} \\"
    echo "    --dir_out ${dir_alf}/tracks \\"
    echo "    --typ_out ${typ_out} \\"
    echo "    --scl_fct ${scl_fct} \\"
    echo "    --dep_min ${dep_min} \\"
    echo "    --err_out ${dir_alf}/tracks/logs \\"
    echo "    --slurm \\"
    echo "         >> >(tee -a ${exc_trk}.stdout.txt) \\"
    echo "        2>> >(tee -a ${exc_trk}.stderr.txt)"
    echo ""
    echo ""
fi

bash "${dir_scr}/${scr_trk}" \
    --verbose \
    --fil_ip "${ser_num}" \
    --fil_in "${ser_den}" \
    --dir_out "${dir_alf}/tracks" \
    --typ_out "${typ_out}" \
    --scl_fct "${scl_fct}" \
    --dep_min "${dep_min}" \
    --err_out "${dir_alf}/tracks/logs" \
    --slurm \
         >> >(tee -a "${exc_trk}.stdout.txt") \
        2>> >(tee -a "${exc_trk}.stderr.txt")


#  Cleanup: Compress logs and remove empty files ------------------------------
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_alf}/tables/logs"
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${dir_alf}/tracks/logs"

# ls -lhaFG "${dir_alf}/tables/logs"  ## Uncomment to check directory for table logs ##
# ls -lhaFG "${dir_alf}/tracks/logs"  ## Uncomment to check directory for track logs ##
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

<a id="h-compute-coverage-using-the-spike-in-method"></a>
### H. Compute coverage using the spike-in method.
<details>
<summary><i>Text: Compute coverage using the spike-in method.</i></summary>
<br />

This section describes the steps to calculate spike-in normalized ChIP-seq coverage.
</details>
<br />

<details>
<summary><i>Bash code: Generate a TSV file of sample-specific spike-in-derived scaling factors.</i></summary>

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
dir_pro="${dir_dat}/processed"
{
    aligner="bowtie2"
    a_type="global"
    req_flg=true
    flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
    mapq=1
    det_bam="flag-${flg}_mapq-${mapq}"
    det_cvg="${aligner}_${a_type}_${det_bam}"
}
dir_aln="${dir_pro}/align_${aligner}_${a_type}"
dir_bam="${dir_aln}/${det_bam}/sc"
dir_cvg="${dir_pro}/compute_coverage"
dir_out="${dir_cvg}/${det_cvg}/spike/tables"
env_nam="env_analyze"
day="$(date '+%Y-%m%d')"

#  Set hardcoded argument assignments, etc.
threads=8
include="IP*,*Brn1*"  # include="IP*,*Hmo1*"  # include="IP*,*Hho1*"
infiles="$(  ## WARNING: Change search parameters as needed ##
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_bam}" \
        --pattern "*.bam" \
        --include "${include}"
)"
case "${include}" in
    IP*Hho1*) outfile="${dir_out}/IP_WT_G1-G2M-Q_Hho1_6336-6337.tsv" ;;
    IP*Hmo1*) outfile="${dir_out}/IP_WT_G1-G2M-Q_Hmo1_7750-7751.tsv" ;;
    IP*Brn1*) outfile="${dir_out}/IP_WT_log-Q_Brn1_rep1-rep2-rep3.tsv" ;;
    *) echo "Error: No matching pattern found for '${include}'" >&2 ;;
esac
err_out="${dir_out}/logs"
scr_mng="${HOME}/miniforge3/etc/profile.d/conda.sh"

#  Using the date and outfile, set path and prefix for driver script logs
exc_tbl="${dir_out}/logs/${day}.execute.$(basename "${outfile}" .tsv)"

#  Create directory structure for storing output tables and tracks associated
#+ with different normalization methods (alpha, spike, norm, raw)
mkdir -p ${dir_cvg}/${det_cvg}/{alpha,spike}/tables/{docs,logs}
mkdir -p ${dir_cvg}/${det_cvg}/{alpha,norm,raw,spike}/tracks/{docs,logs}

#  Source utility functions
source "${dir_fnc}/check_program_path.sh"
source "${dir_fnc}/echo_warning.sh"
source "${dir_fnc}/handle_env.sh"

#  Activate the required environment
handle_env "${env_nam}"

#  Check the availability of necessary dependencies such as GNU Parallel,
#+ Python, and SLURM sbatch
check_program_path awk
check_program_path parallel
check_program_path python
check_program_path samtools
check_program_path sbatch ||
    echo_warning \
        "SLURM is not available on this system. Do not use the '--slurm'" \
        "flag with the driver script."

#  Run the driver script to calculate per-sample spike-in-derived scaling
#+ factors
bash "${dir_scr}/execute_calculate_scaling_factor_spike.sh" \
    --verbose \
    --threads "${threads}" \
    --infiles "${infiles}" \
    --outfile "${outfile}" \
    --err_out "${err_out}" \
    --flg_mc \
    --slurm \
         > >(tee -a "${exc_tbl}.stdout.txt") \
        2> >(tee -a "${exc_tbl}.stderr.txt")

#  Relativize the scaling factors to the maximum value, and sort the outfile
#+ rows
python "${dir_scr}/relativize_scaling_factors.py" --infile "${outfile}" \
    | awk 'NR == 1; NR > 1 && NF { print | "sort" }' \
        > "${dir_out}/tmp.txt"

#  Replace the original outfile with the newly relativized and sorted version
mv -f "${dir_out}/tmp.txt" "${outfile}"

#  Optional: Check the contents of the outfile
# cat "${outfile}"
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
exc_tbl="${err_out}/${day}.execute.${nam_job}"

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
                         > >(tee -a "${exc_tbl}.${outstem}.stdout.txt") \\
                        2> >(tee -a "${exc_tbl}.${outstem}.stderr.txt")
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
                     > >(tee -a "${exc_tbl}.${outstem}.stdout.txt") \
                    2> >(tee -a "${exc_tbl}.${outstem}.stderr.txt")

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
