
ChIP-seq Protocol Workflow
==========================

**Supporting code and documentation for the *Bio-protocol* manuscript "An Introduction to ChIP-seq Data Analysis Using *Saccharomyces cerevisiae* with a Focus on Relative and Quantitative Signal Normalization."**

**Author:** *Kris Alavattam*

This notebook provides a guide to the ChIP-seq data analysis workflow detailed in the manuscript, including code snippets, explanations, and step-by-step instructions.

Note: If using a high-performance computing cluster (HPCC), request an interactive node to ensure adequate resources for running code in the below chunks. The specific command (e.g., `grabnode` at Fred Hutch Cancer Center) will depend on the job scheduler setup. This step is unnecessary if running the code on a local machine.

Note: For detailed instructions on keeping your local version of the 202X_protocol_ChIP repository up-to-date, please see [this GitHub gist](https://gist.github.com/kalavattam/76f123011e8dcd77b445a72d23a64036).

---
<br />

## F. Generate Bowtie 2 indices from the concatenated FASTA file.
<details>
<summary><i>Text: Generate Bowtie 2 indices from the concatenated FASTA file.</i></summary>
<br />

To align ChIP-seq reads against both *S. cerevisiae* and *S. pombe* genomes, we first generate Bowtie 2 indices from a concatenated FASTA file. This ensures efficient and accurate alignment for, e.g., spike-in normalization (described below).

Steps overview:
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

#  Optional: Request an interactive node
grabnode  # Request 1 core, 20 GB memory, 1 day, no GPU

#  Define variables for directory paths, etc.
## WARNING: Change path if you're not Kris ##
dir_bas="${HOME}/tsukiyamalab/Kris"  ## WARNING: Change if not Kris ##
dir_rep="${dir_bas}/202X_protocol_ChIP"
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

## G. Obtain and organize ChIP-seq FASTQ files.
<details>
<summary><i>Bash code: Obtain and organize ChIP-seq FASTQ files.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node
grabnode  # Request 1 core, 20 GB memory, 1 day, no GPU

#  Define variables for directory paths, etc.
dir_bas="${HOME}/tsukiyamalab/Kris"  ## WARNING: Change if not Kris ##
dir_rep="${dir_bas}/202X_protocol_ChIP"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"
dir_raw="${dir_rep}/data/raw"
dir_doc="${dir_raw}/docs"
fil_tsv="datasets.tsv"  ## WARNING: Change as needed ##
pth_tsv="${dir_doc}/${fil_tsv}"
dir_log="${dir_raw}/logs"
dir_sym="${dir_rep}/data/symlinked"
env_nam="env_align"
threads=6
time="6:00:00"
day="$(date '+%Y-%m%d')"

#  Source utility functions
source "${dir_fnc}/check_program_path.sh"
source "${dir_fnc}/handle_env.sh"

#  Activate the required environment
handle_env "${env_nam}"

#  Ensure access to bowtie2-build
check_program_path "parallel"

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

## H. Use Atria to perform adapter and quality trimming of sequenced reads.
<details>
<summary><i>Bash code: Use Atria to perform adapter and quality trimming of sequenced reads.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node
grabnode  # Request 1 core, 20 GB memory, 1 day, no GPU

#  Define variables for directory paths, environment, threads, and infiles
dir_bas="${HOME}/tsukiyamalab/Kris"  ## WARNING: Change if not Kris ##
dir_rep="${dir_bas}/202X_protocol_ChIP"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"
dir_dat="${dir_rep}/data"
dir_sym="${dir_dat}/symlinked"
dir_pro="${dir_dat}/processed"
dir_trm="${dir_pro}/trim_atria_FASTQ"
env_nam="env_analyze"
threads=4
infiles="$(  ## WARNING: Change the search parameters as needed ##
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
source "${dir_fnc}/handle_env.sh"

#  Activate the required environment
handle_env "${env_nam}"

#  Check availability of Atria and other necessary tools
check_program_path "atria"
check_program_path "pbzip2"
check_program_path "pigz"

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
         > >(tee -a "${dir_trm}/logs/${day}.execute.stdout.txt") \
        2> >(tee -a "${dir_trm}/logs/${day}.execute.stderr.txt")

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

## I. Align sequenced reads with Bowtie 2 and process the read alignments.
<details>
<summary><i>Bash code: Align sequenced reads with Bowtie 2 and process the read alignments.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node
grabnode  # Request 1 core, 20 GB memory, 1 day, no GPU

#  Define variables for directory paths, environment, driver script parameters,
#+ and so on
## WARNING: Change path if you're not Kris ##
dir_bas="${HOME}/tsukiyamalab/Kris"  ## WARNING: Change if not Kris ##
dir_rep="${dir_bas}/202X_protocol_ChIP"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"
dir_dat="${dir_rep}/data"
dir_idx="${dir_dat}/genomes/concat/index"
dir_pro="${dir_dat}/processed"
dir_trm="${dir_pro}/trim_atria_FASTQ"
env_nam="env_align"
threads=8
aligner="bowtie2"
a_type="global"
mapq=1
req_flg=true
flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
str_idx="sc_sp_proc"
pth_idx="${dir_idx}/${aligner}/${str_idx}"
infiles="$(  ## WARNING: Change the search parameters as needed ##
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_trm}" \
        --pattern "*.atria.fastq.gz" \
        --depth 1 \
        --fastqs \
        --include "*Q*,*Brn1*" \
        --exclude "*msn2*,*SMC4off*,*merged*"
);$(  ## TODO: Delete semicolon and subsequent command after testing ##
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_trm}" \
        --pattern "*.atria.fastq.gz" \
        --depth 1 \
        --fastqs \
        --include "*Hho1*"
)"
dir_aln="${dir_pro}/align_${aligner}_${a_type}_BAM"
dir_out="${dir_aln}/flag-${flg}_mapq-${mapq}"
nam_job="align_fastqs"
max_job=6
time="1:00:00"
day="$(date '+%Y-%m%d')"

#  Create output directory structure for trimmed FASTQ files and logs
mkdir -p ${dir_out}/{init,sc,sp}/{docs,logs}

#  Source utility functions
source "${dir_fnc}/check_program_path.sh"
source "${dir_fnc}/handle_env.sh"

#  Activate the required environment
handle_env "${env_nam}"

#  Check availability of Atria and other necessary tools
check_program_path "bowtie2"
check_program_path "samtools"

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
         > >(tee -a "${dir_out}/init/logs/${day}.execute.stdout.txt") \
        2> >(tee -a "${dir_out}/init/logs/${day}.execute.stderr.txt")

#  Adjust variable 'infiles' assignment for filtering BAM files
infiles="$(  ## WARNING: Change the search parameters as needed ##
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
         > >(tee -a "${dir_out}/sc/logs/${day}.execute.stdout.txt") \
        2> >(tee -a "${dir_out}/sc/logs/${day}.execute.stderr.txt")

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
         > >(tee -a "${dir_out}/sp/logs/${day}.execute.stdout.txt") \
        2> >(tee -a "${dir_out}/sp/logs/${day}.execute.stderr.txt")

#  Cleanup: Compress large stdout and stderr files, and remove files with size
#+ 0
bash "${dir_scr}/compress_remove_files.sh" \
    --dir_fnd "${dir_out}/init/logs"

bash "${dir_scr}/compress_remove_files.sh" \
    --dir_fnd "${dir_out}/sc/logs"

bash "${dir_scr}/compress_remove_files.sh" \
    --dir_fnd "${dir_out}/sp/logs"

#  Optional: Check the contents of the logs directory
# ls -lhaFG "${dir_trm}/init/logs"
# ls -lhaFG "${dir_trm}/sc/logs"
# ls -lhaFG "${dir_trm}/sp/logs"
```
</details>
<br />
