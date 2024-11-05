
ChIP-seq Protocol Workflow
==========================

**Supporting code and documentation for the *Bio-protocol* manuscript "An Introduction to ChIP-seq Data Analysis Using *Saccharomyces cerevisiae* with a Focus on Relative and Quantitative Signal Normalization."**

**Author:** *Kris Alavattam*

This notebook provides a guide to the ChIP-seq data analysis workflow detailed in the manuscript, including code snippets, explanations, and step-by-step instructions.

Note: If using a high-performance computing cluster (HPCC), request an interactive node to ensure adequate resources for running code in the below chunks. The specific command (e.g., `grabnode` at Fred Hutch Cancer Center) will depend on the job scheduler setup. This step is unnecessary if running the code on a local machine.

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
dir_rep="${HOME}/tsukiyamalab/Kris/202X_protocol_ChIP"  
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

#  Source utility functions
source "${dir_fnc}/check_program_path.sh"
source "${dir_fnc}/handle_env_activate.sh"

#  Activate the environment for alignment tools, etc.
handle_env_activate "${env_nam}"

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
     > >(tee -a "${dir_idx}/logs/stdout.txt") \
    2> >(tee -a "${dir_idx}/logs/stderr.txt")

#  Optional: Once the index is built, delete the decompressed FASTA file
if [[ -f "${pth_fas}" ]]; then rm "${pth_fas}"; fi
```
</details>
<br />

<details>
<summary><i>Bash code: Write accompanying documentation.</i></summary>

```bash
#!/bin/bash

#  This is a template example: Be sure to update paths, file names, and other
#+ details specific to your work before running it
cat << EOF > "${dir_idx}/docs/2024-1105.txt"
# Bowtie 2 index generation documentation
# Date: 2024-1105
# Author: Kris Alavattam

## Purpose
Generate Bowtie 2 indices from a concatenated FASTA file (sc_sp_proc.fasta) for simultaneous alignment of ChIP-seq reads against S. cerevisiae and S. pombe genomes.

## Environment setup
- Activated environment: "${env_nam}"
- Required tool: bowtie2-build

## Directory structure
- Genome data directory: "${dir_gen}"
- Concatenated file directory: "${dir_cat}"
- Concatenated FASTA directory: "${dir_fas}"
- Index output directory: "${dir_idx}"

## Steps
1. Interactive node request: grabnode (1 core, 20 GB memory, 1 day, no GPU).
2. Environment activation: Activated "${env_nam}" using custom function handle_env_activate.
3. Tool check: Ensured bowtie2-build is in PATH using custom function check_program_path.
4. Directory setup: Created necessary directories for indices, docs, and logs.
5. Decompression: Checked and decompressed the FASTA file if needed.
6. Index building: Ran bowtie2-build on the decompressed FASTA
    A. Output index files to "${dir_idx}".
    B. Output logs to "${dir_idx}/logs".
    C. After index building, contents of "${dir_idx}":
    \`\`\`
    ❯ ls -lhaFG "\${dir_idx}"
    total 49M
    drwxrws--- 4 kalavatt  256 Nov  5 05:43 ./
    drwxrws--- 3 kalavatt   25 Nov  5 05:27 ../
    drwxrws--- 2 kalavatt   31 Nov  5 05:59 docs/
    drwxrws--- 2 kalavatt   56 Nov  5 05:42 logs/
    -rw-rw---- 1 kalavatt  12M Nov  5 05:42 sc_sp_proc.1.bt2
    -rw-rw---- 1 kalavatt 6.0M Nov  5 05:42 sc_sp_proc.2.bt2
    -rw-rw---- 1 kalavatt  269 Nov  5 05:42 sc_sp_proc.3.bt2
    -rw-rw---- 1 kalavatt 6.0M Nov  5 05:42 sc_sp_proc.4.bt2
    -rw-rw---- 1 kalavatt  12M Nov  5 05:43 sc_sp_proc.rev.1.bt2
    -rw-rw---- 1 kalavatt 6.0M Nov  5 05:43 sc_sp_proc.rev.2.bt2
    \`\`\`
7. Cleanup: Optionally removed the decompressed FASTA after index generation.
    A. Upon cleanup, contents of "${dir_fas}":
    \`\`\`
    ❯ ls -lhaFG "\${dir_fas}"
    total 7.5M
    drwxrws--- 2 kalavatt   37 Nov  5 06:03 ./
    drwxrws--- 3 kalavatt   22 Oct 22 10:56 ../
    -rw-rw---- 1 kalavatt 7.4M Oct 22 11:06 sc_sp_proc.fasta.gz
    \`\`\`

## Output
- Bowtie 2 indices stored in "${dir_idx}"
- Logs: stdout and stderr logged to "${dir_idx}/logs"

EOF

# cat "${dir_idx}/docs/2024-1105.txt"
# rm "${dir_idx}/docs/2024-1105.txt"
```
</details>
<br />

## G. Obtain and organize ChIP-seq FASTQ files.
<br />

## H. Use Atria to perform adapter and quality trimming of sequenced reads.
<details>
<summary><i>Bash code: Use Atria to perform adapter and quality trimming of sequenced reads.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node
grabnode  # Request 1 core, 20 GB memory, 1 day, no GPU

#  Define variables for directory paths, environment, and threads
## WARNING: Change path if you're not Kris ##
dir_rep="${HOME}/tsukiyamalab/Kris/202X_protocol_ChIP"  
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"
dir_dat="${dir_rep}/data"
dir_sym="${dir_dat}/symlinked"
dir_pro="${dir_dat}/processed"
dir_trm="${dir_pro}/trim_atria_FASTQ"
env_nam="env_analyze"
threads=4

#  Source utility functions
source "${dir_fnc}/check_program_path.sh"
source "${dir_fnc}/handle_env_activate.sh"

#  Activate the required environment
handle_env_activate "${env_nam}"

#  Check availability of Atria and other necessary tools
check_program_path "atria"
check_program_path "pbzip2"
check_program_path "pigz"

#  Create output directory structure for trimmed FASTQ files and logs
mkdir -p ${dir_trm}/{docs,logs}

#  Find and format FASTQ files into a semicolon- and comma-delimited string
infiles="$(
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_sym}" \
        --pattern "*.fastq.gz" \
        --depth 1 \
        --follow \
        --fastqs
)"

#  Run the driver script to trim FASTQ files with Atria
bash "${dir_scr}/execute_trim_fastqs.sh" \
    --verbose \
    --threads ${threads} \
    --infiles "${infiles}" \
    --dir_out "${dir_trm}" \
    --err_out "${dir_trm}/logs" \
    --slurm \
         > >(tee -a "${dir_trm}/logs/2024-1105.execute.stdout.txt") \
        2> >(tee -a "${dir_trm}/logs/2024-1105.execute.stderr.txt")


#  Move Atria LOG and JSON files to the logs directory
mv ${dir_trm}/*.{log,json} "${dir_trm}/logs"

#  Optional: Compress large stdout, stderr, LOG, and JSON files, and remove
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
