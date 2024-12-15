
ChIP-seq Protocol Workflow
==========================

**Supporting code and documentation for the *bioRxiv* manuscript "ChIP-seq data processing and relative and quantitative signal normalization for *Saccharomyces cerevisiae*."**

**Author:** *Kris Alavattam*

This notebook provides a guide to the ChIP-seq data processing workflow detailed in the manuscript, including code snippets, explanations, and step-by-step instructions.

Note: If using a high-performance computing cluster (HPCC), request an interactive node to ensure adequate resources for running code in the below chunks. The specific command (e.g., `grabnode` at Fred Hutch Cancer Center) will depend on the job scheduler setup. This step is unnecessary if running the code on a local machine.

Note: For detailed instructions on keeping your local version of the [`protocol_chipseq_signal_norm`](https://github.com/kalavattam/protocol_chipseq_signal_norm) repository up-to-date, please see [this GitHub gist](https://gist.github.com/kalavattam/76f123011e8dcd77b445a72d23a64036).

---
<br />

## Procedures
`#TODO`
<br />
<br />

## Data analysis
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

`#TODO`
</details>
<br />

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

#  Optional: Request an interactive node
grabnode  # Request 1 core, 20 GB memory, 1 day, no GPU

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

### C. Obtain and organize ChIP-seq FASTQ files.
<details>
<summary><i>Bash code: Obtain and organize ChIP-seq FASTQ files.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node
grabnode  # Request 1 core, 20 GB memory, 1 day, no GPU

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

### D. Use Atria to perform adapter and quality trimming of sequenced reads.
<details>
<summary><i>Bash code: Use Atria to perform adapter and quality trimming of sequenced reads.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node
grabnode  # Request 1 core, 20 GB memory, 1 day, no GPU

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

### E. Align sequenced reads with Bowtie 2 and process the read alignments.
<details>
<summary><i>Bash code: Align sequenced reads with Bowtie 2 and process the read alignments.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node
grabnode  # Request 1 core, 20 GB memory, 1 day, no GPU

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

### F. Compute normalized (or raw) coverage.
<details>
<summary><i>Text: Compute raw or normalized (fractional) coverage.</i></summary>
<br />

This following Bash code chunk provides an example of how to compute ChIP-seq coverage, either as raw (unadjusted for sequencing depth or other technical biases) or normalized (fractional) per [Dickson et al., *Sci Rep*, 2023](https://www.nature.com/articles/s41598-023-34430-2). The coverage type is determined by setting the variable `typ_cov` to `"raw"` or `"norm"`. BIGWIG and log output files will be saved to separate directories based on the selected coverage type.
</details>
<br />

<details>
<summary><i>Bash code: Compute raw or normalized (fractional) coverage.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node
grabnode  # Request 1 core, 20 GB memory, 1 day, no GPU

#  Define variables for directory paths, environment, driver script arguments,
#+ etc.
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
det_bam="flag-${flg}_mapq-${mapq}"
det_cov="${aligner}_${a_type}_${det_bam}"
typ_cov="norm"  ## WARNING: "raw" for unadjusted, "norm" for normalized ##

dir_aln="${dir_pro}/align_${aligner}_${a_type}"
dir_bam="${dir_aln}/${det_bam}/sc"
dir_cov="${dir_pro}/compute_coverage"
dir_trk="${dir_cov}/${det_cov}/${typ_cov}/tracks"

env_nam="env_analyze"
day="$(date '+%Y-%m%d')"
nam_job="compute_coverage_${typ_cov}"
typ_out="bigwig"
bin_siz=1

#  Set hardcoded argument assignments, etc.
threads=8
infiles="$(  ## WARNING: Change search parameters as needed ##
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_bam}" \
        --pattern "*.bam"
)"
err_out="${dir_trk}/logs"

#  Using the date and outfile, set path and prefix for driver script logs
exc_pth="${dir_trk}/logs/${day}.execute.${nam_job}"

#  Create directory structure for storing output tables and tracks associated
#+ with different normalization methods (alpha, spike, norm, raw)
mkdir -p ${dir_cov}/${det_cov}/{alpha,spike}/tables/{docs,logs}
mkdir -p ${dir_cov}/${det_cov}/{alpha,norm,raw,spike}/tracks/{docs,logs}

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
check_program_path sbatch ||
    echo_warning \
        "SLURM is not available on this system. Do not use the '--slurm'" \
        "flag with the driver script."

#  Run the driver script to calculate per-sample spike-in-derived scaling
#+ factors
bash "${dir_scr}/execute_compute_coverage.sh" \
    --verbose \
    --threads "${threads}" \
    --infiles "${infiles}" \
    --dir_out "${dir_trk}" \
    --typ_out "${typ_out}" \
    --bin_siz "${bin_siz}" \
    $(
        if [[ "${typ_cov}" == "norm" ]]; then
            echo "--norm"
        fi
    ) \
    --err_out "${err_out}" \
    --nam_job "${nam_job}" \
    --slurm \
         >> >(tee -a "${exc_pth}.stdout.txt") \
        2>> >(tee -a "${exc_pth}.stderr.txt")

#  Cleanup: Compress large stdout, stderr, LOG, and JSON files, and remove
#+ files with size 0
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${err_out}"

#  Optional: Check the contents of the logs directory
# ls -lhaFG "${err_out}"
```
</details>
<br />

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
#  Define directory paths
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
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
det_cov="${aligner}_${a_type}_${det_bam}"
typ_cov="alpha"

#  Further define directory setup
dir_aln="${dir_pro}/align_${aligner}_${a_type}"
dir_bam="${dir_aln}/${det_bam}/sc"
dir_cov="${dir_pro}/compute_coverage"
dir_det="${dir_cov}/${det_cov}/${typ_cov}"
dir_tbl="${dir_det}/tables"
dir_trk="${dir_det}/tracks"
eo_tbl="${dir_tbl}/logs"
eo_trk="${dir_trk}/logs"

#  Define environment, resources, and script arguments 
env_nam="env_analyze"
threads=8
mes_tbl="${dir_raw}/docs/measurements_siq_chip.tsv"
bin_siz=1

#  Define file search parameters
## WARNING: Change search parameters as needed ##
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
scr_tbl="execute_calculate_scaling_factor_${typ_cov}.sh"
scr_trk="execute_deeptools_coverage.sh"
fil_tbl="${dir_tbl}/IP_WT_G1-G2M-Q_Hho1-Hmo1_6336-6337_7750-7751.tsv"

#  Define log file prefixes
day="$(date '+%Y-%m%d')"
exc_tbl="${eo_tbl}/${day}.execute.${scr_tbl%.sh}.$(
    basename "${fil_tbl}" .tsv
)"
exc_trk="${eo_trk}/${day}.${scr_trk%.sh}"


#  Create required directories if necessary -----------------------------------
mkdir -p ${dir_tbl}/{docs,logs}
mkdir -p ${dir_trk}/{docs,logs}


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


#  Calculate siQ-ChIP alpha scaling factors -----------------------------------
if [[ ! "${fil_tbl}" ]]; then
    #  Run the driver script to generate a TSV file of sample-specific siQ-ChIP
    #+ alpha scaling factors
    bash "${dir_scr}/execute_calculate_scaling_factor_${typ_cov}.sh" \
        --verbose \
        --threads "${threads}" \
        --infiles "${infiles}" \
        --table "${mes_tbl}" \
        --outfile "${fil_tbl}" \
        --err_out "${eo_tbl}" \
        --flg_dep \
        --flg_len \
        --flg_mc \
        --slurm \
             >> >(tee -a "${exc_tbl}.stdout.txt") \
            2>> >(tee -a "${exc_tbl}.stderr.txt")

    #  Sort the table of scaling factors by rows
    awk 'NR == 1; NR > 1 && NF { print | "sort" }' "${fil_tbl}" \
        > "${dir_tbl}/tmp.tsv"

    #  Replace the original table with the sorted version
    mv -f "${dir_tbl}/tmp.tsv" "${fil_tbl}"

    # cat "${fil_tbl}"  ## Uncomment to check the table contents ##
fi

#  Generate alpha-scaled signal tracks ----------------------------------------
#  Use the TSV file to generate alpha-scaled signal tracks
bash "${dir_scr}/execute_deeptools_coverage.sh" \
    --verbose \
    --threads "${threads}" \
    --table "${fil_tbl}" \
    --tbl_col "${typ_cov}" \
    --dir_out "${dir_trk}" \
    --bin_siz "${bin_siz}" \
    --err_out "${eo_trk}" \
    --slurm \
         >> >(tee -a "${exc_trk}.stdout.txt") \
        2>> >(tee -a "${exc_trk}.stderr.txt")


#  Cleanup: Compress logs and remove empty files ------------------------------
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${eo_tbl}"
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${eo_trk}"

# ls -lhaFG "${eo_tbl}"  ## Uncomment to check directory for table logs ##
# ls -lhaFG "${eo_trk}"  ## Uncomment to check directory for track logs ##
```
</details>
<br />

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

#  Optional: Request an interactive node
grabnode  # Request 1 core, 20 GB memory, 1 day, no GPU

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
    det_cov="${aligner}_${a_type}_${det_bam}"
}
dir_aln="${dir_pro}/align_${aligner}_${a_type}"
dir_bam="${dir_aln}/${det_bam}/sc"
dir_cov="${dir_pro}/compute_coverage"
dir_out="${dir_cov}/${det_cov}/spike/tables"
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
mkdir -p ${dir_cov}/${det_cov}/{alpha,spike}/tables/{docs,logs}
mkdir -p ${dir_cov}/${det_cov}/{alpha,norm,raw,spike}/tracks/{docs,logs}

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
nam_job="comp_covg_spike"
no_infiles="$(tail -n +2 "${outfile}" | wc -l)"
max_job=6
threads=8
bin_siz=1
outtype="bigwig"
time="0:30:00"
dir_out="${dir_pro}/compute_coverage/${det_cov}/spike/tracks"
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
                    --outfile "${dir_cov}/${det_cov}/spike/tracks/${outstem}" \\
                    --scl_fct "${scaled}" \\
                    --outtype "${outtype}" \\
                    --threads ${threads} \\
                    --bin_siz 1 \\
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
                --outfile "${dir_cov}/${det_cov}/spike/tracks/${outstem}" \
                --scl_fct "${scaled}" \
                --outtype "${outtype}" \
                --threads ${threads} \
                --bin_siz 1 \
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
