
Workflow: Validate New siQ-ChIP Implementation with the Original
================================================================

**Instructions to install and run the [original siQ-ChIP implementation](https://github.com/kalavattam/siQ-ChIP/tree/cerevisiae) by Brad Dickson. Outputs from the original implementation (e.g., $\alpha$ values, normalized coverage, siQ-ChIP-scaled coverage) serve as validation standards to confirm the accuracy/reproducibility of the [newly implemented siQ-ChIP computations](https://github.com/kalavattam/protocol_chipseq_signal_norm).**

**Author:** *Kris Alavattam*  

---
<br />

## Table of contents
<details>
<summary><i>Table of contents</i></summary>
<br />
<!-- MarkdownTOC -->

1. [Procedures](#procedures)
    1. [A. Clone the forked siQ-ChIP repository.](#a-clone-the-forked-siq-chip-repository)
    1. [B. Create an environment for running siQ-ChIP.](#b-create-an-environment-for-running-siq-chip)
1. [Data analysis](#data-analysis)
    1. [A. QNAME-sort BAM files filtered for *S. cerevisiae* alignments](#a-qname-sort-bam-files-filtered-for-s-cerevisiae-alignments)
    1. [B. Generate BED files from QNAME-sorted BAM files](#b-generate-bed-files-from-qname-sorted-bam-files)
    1. [C. Run the initial implementation of siQ-ChIP](#c-run-the-initial-implementation-of-siq-chip)

<!-- /MarkdownTOC -->
</details>
<br />

<a id="procedures"></a>
## Procedures
<a id="a-clone-the-forked-siq-chip-repository"></a>
### A. Clone the forked siQ-ChIP repository.
<details>
<summary><i>Text: Clone the forked siQ-ChIP repository.</i></summary>
<br />

`#TODO`
</details>
<br />

<details>
<summary><i>Code: Clone the forked siQ-ChIP repository.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node
grabnode  # Request 1 core, 20 GB memory, 1 day, no GPU

#  Define function to write a warning message to stderr and return code 1
function echo_warning() { echo "Warning: $*" >&2; return 1; }


#  Define variables for directory paths and environment name
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/siQ-ChIP"
git_rep="kalavattam/siQ-ChIP"
branch="cerevisiae"

#  Go to base repository directory
cd "${dir_bas}" || echo_warning "cd'ing failed; check on this"

#  Clone forked siQ-ChIP repository
if [[ ! -d "${dir_rep}" ]]; then gh repo clone "${git_rep}"; fi

#  Checkout the 'cerevisiae' branch
cd "${dir_rep}" || echo_warning "cd'ing failed; check on this"

#  Switch to the specified branch
git checkout -b "${branch}"

echo "Repository is ready and on branch '${branch}'."
```
</details>
<br />

<a id="b-create-an-environment-for-running-siq-chip"></a>
### B. Create an environment for running siQ-ChIP.
<details>
<summary><i>Text: Create an environment for running siQ-ChIP.</i></summary>
<br />

`#TODO`
</details>
<br />

<details>
<summary><i>Bash code: Create an environment for running siQ-ChIP.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node
grabnode  # Request 1 core, 20 GB memory, 1 day, no GPU

#  Define variables for directory paths and environment name
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
env_nam="env_siqchip"

#  Create an environment containing siQ-ChIP dependencies
bash "${dir_scr}/install_envs.sh" \
    --env_nam "${env_nam}" \
    --yes
```
</details>
<br />

<a id="data-analysis"></a>
## Data analysis
<a id="a-qname-sort-bam-files-filtered-for-s-cerevisiae-alignments"></a>
### A. QNAME-sort BAM files filtered for *S. cerevisiae* alignments
<details>
<summary><i>Bash code: QNAME-sort BAM files filtered for </i>S. cerevisiae<i> alignments</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node
grabnode  # Request 1 core, 20 GB memory, 1 day, no GPU

#  Define variables for directory paths, environment, submission script
#+ arguments, metadata, and so on
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
species="sc"

dir_aln="${dir_pro}/align_${aligner}_${a_type}"
dir_bam="${dir_aln}/flag-${flg}_mapq-${mapq}/${species}"

pattern="*.bam"
exclude="*Brn1*"
depth=1
infiles="$(  ## WARNING: Change search parameters as needed ##
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_bam}" \
        --pattern "${pattern}" \
        --exclude "${exclude}" \
        --depth ${depth}
)"

scr_sub="${dir_scr}/submit_qsort_bam_slurm.sh"
env_nam="env_siqchip"
dir_out="${dir_bam}_qnam"
err_out="${dir_out}/logs"

day="$(date '+%Y-%m%d')"
nam_job="qsort_bams"
threads=8
time="2:00:00"
num_job=$(awk -F "," '{ print NF }' <<< "${infiles}")
max_job=6

#  Create output directory structure for trimmed FASTQ files and logs
mkdir -p ${dir_out}/{docs,logs}

#  Debug variable assignments
if ${debug:-true}; then
    {
        echo "####################################"
        echo "## Hardcoded variable assignments ##"
        echo "####################################"
        echo ""
        echo "\${dir_bas}=${dir_bas}"
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
        echo "\${species}=${species}"
        echo ""
        echo "\${dir_aln}=${dir_aln}"
        echo "\${dir_bam}=${dir_bam}"
        echo ""
        echo "\${pattern}=${pattern}"
        echo "\${exclude}=${exclude}"
        echo "\${depth}=${depth}"
        echo "\${infiles}=${infiles}"
        echo ""
        echo "\${scr_sub}=${scr_sub}"
        echo "\${env_nam}=${env_nam}"
        echo "\${dir_out}=${dir_out}"
        echo "\${err_out}=${err_out}"
        echo ""
        echo "\${day}=${day}"
        echo "\${nam_job}=${nam_job}"
        echo "\${threads}=${threads}"
        echo "\${time}=${time}"
        echo "\${num_job}=${num_job}"
        echo "\${max_job}=${max_job}"
        echo ""
        echo ""
    } >> >(tee -a "${err_out}/${day}.execute.stdout.txt")
fi

#  Source utility functions
source "${dir_fnc}/check_program_path.sh"
source "${dir_fnc}/echo_warning.sh"
source "${dir_fnc}/handle_env.sh"

#  Activate the required environment
handle_env "${env_nam}"

#  Check the availability of Samtools and SLURM sbatch
check_program_path samtools
check_program_path sbatch

#  Debug call to sbatch with submission script
if ${debug:-true}; then
    {
        echo "####################"
        echo "## Call to sbatch ##"
        echo "####################"
        echo ""
        echo "sbatch \\"
        echo "    --job-name=${nam_job} \\"
        echo "    --nodes=1 \\"
        echo "    --cpus-per-task=${threads} \\"
        echo "    --time=${time} \\"
        echo "    --error=${err_out}/${nam_job}.%A-%a.stderr.txt \\"
        echo "    --output=${err_out}/${nam_job}.%A-%a.stdout.txt \\"
        echo "    --array=1-${num_job}%${max_job} \\"
        echo "    ${scr_sub} \\"
        echo "        ${env_nam} \\"
        echo "        ${threads} \\"
        echo "        ${infiles} \\"
        echo "        ${dir_out} \\"
        echo "        ${err_out} \\"
        echo "        ${nam_job}"
        echo ""
        echo ""
        echo "#########################################"
        echo "## Contents of SLURM submission script ##"
        echo "#########################################"
        echo ""
        echo "## ${scr_sub} ##"
        echo ""
        cat "${scr_sub}"
        echo ""
    } >> >(tee -a "${err_out}/${day}.execute.stdout.txt")
fi

#  Run SLURM submission script to QNAME-sort BAM files 
# shellcheck disable=SC2046,SC2086
sbatch \
    --job-name=${nam_job} \
    --nodes=1 \
    --cpus-per-task=${threads} \
    --time=${time} \
    --error=${err_out}/${nam_job}.%A-%a.stderr.txt \
    --output=${err_out}/${nam_job}.%A-%a.stdout.txt \
    --array=1-${num_job}%${max_job} \
    ${scr_sub} \
        ${env_nam} \
        ${threads} \
        ${infiles} \
        ${dir_out} \
        ${err_out} \
        ${nam_job}

#  Compress large stdout and stderr files, and remove files with size 0
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${err_out}"
```
</details>
<br />

<a id="b-generate-bed-files-from-qname-sorted-bam-files"></a>
### B. Generate BED files from QNAME-sorted BAM files
<details>
<summary><i>Bash code: Generate BED files from QNAME-sorted BAM files</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node
grabnode  # Request 1 core, 20 GB memory, 1 day, no GPU

#  Define variables for directory paths, environment, submission script
#+ arguments, metadata, and so on
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
spe_typ="sc_qnam"

dir_aln="${dir_pro}/align_${aligner}_${a_type}"
dir_bam="${dir_aln}/flag-${flg}_mapq-${mapq}/${spe_typ}"

pattern="*.bam"
depth=1
infiles="$(  ## WARNING: Change search parameters as needed ##
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_bam}" \
        --pattern "${pattern}" \
        --depth ${depth}
)"

scr_sub="${dir_scr}/submit_convert_bam_bed_slurm.sh"
env_nam="env_siqchip"
dir_out="${dir_bam%_qnam}_bed"
err_out="${dir_out}/logs"

day="$(date '+%Y-%m%d')"
nam_job="convert_bam_bed"
threads=8
time="2:00:00"
num_job=$(awk -F "," '{ print NF }' <<< "${infiles}")
max_job=6

#  Create output directory structure for trimmed FASTQ files and logs
mkdir -p ${dir_out}/{docs,logs}

#  Debug variable assignments
if ${debug:-true}; then
    {
        echo "####################################"
        echo "## Hardcoded variable assignments ##"
        echo "####################################"
        echo ""
        echo "\${dir_bas}=${dir_bas}"
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
        echo "\${spe_typ}=${spe_typ}"
        echo ""
        echo "\${dir_aln}=${dir_aln}"
        echo "\${dir_bam}=${dir_bam}"
        echo ""
        echo "\${pattern}=${pattern}"
        echo "\${depth}=${depth}"
        echo "\${infiles}=${infiles}"
        echo ""
        echo "\${scr_sub}=${scr_sub}"
        echo "\${env_nam}=${env_nam}"
        echo "\${dir_out}=${dir_out}"
        echo "\${err_out}=${err_out}"
        echo ""
        echo "\${day}=${day}"
        echo "\${nam_job}=${nam_job}"
        echo "\${threads}=${threads}"
        echo "\${time}=${time}"
        echo "\${num_job}=${num_job}"
        echo "\${max_job}=${max_job}"
        echo ""
        echo ""
    } >> >(tee -a rm "${err_out}/${day}.execute.stdout.txt")
fi

#  Source utility functions
source "${dir_fnc}/check_program_path.sh"
source "${dir_fnc}/echo_warning.sh"
source "${dir_fnc}/handle_env.sh"

#  Activate the required environment
handle_env "${env_nam}"

#  Check the availability of Samtools and SLURM sbatch
check_program_path samtools
check_program_path sbatch

#  Debug call to sbatch with submission script
if ${debug:-true}; then
    {
        echo "####################"
        echo "## Call to sbatch ##"
        echo "####################"
        echo ""
        echo "sbatch \\"
        echo "    --job-name=${nam_job} \\"
        echo "    --nodes=1 \\"
        echo "    --cpus-per-task=${threads} \\"
        echo "    --time=${time} \\"
        echo "    --error=${err_out}/${nam_job}.%A-%a.stderr.txt \\"
        echo "    --output=${err_out}/${nam_job}.%A-%a.stdout.txt \\"
        echo "    --array=1-${num_job}%${max_job} \\"
        echo "    ${scr_sub} \\"
        echo "        ${env_nam} \\"
        echo "        ${threads} \\"
        echo "        ${infiles} \\"
        echo "        ${dir_out} \\"
        echo "        ${err_out} \\"
        echo "        ${nam_job}"
        echo ""
        echo ""
        echo "#########################################"
        echo "## Contents of SLURM submission script ##"
        echo "#########################################"
        echo ""
        echo "## ${scr_sub} ##"
        echo ""
        cat "${scr_sub}"
        echo ""
    } >> >(tee -a less "${err_out}/${day}.execute.stdout.txt")
fi

#  Run SLURM submission script to generate BED files from QNAME-sorted BAM
#+ files
# shellcheck disable=SC2046,SC2086
sbatch \
    --job-name=${nam_job} \
    --nodes=1 \
    --cpus-per-task=${threads} \
    --time=${time} \
    --error=${err_out}/${nam_job}.%A-%a.stderr.txt \
    --output=${err_out}/${nam_job}.%A-%a.stdout.txt \
    --array=1-${num_job}%${max_job} \
    ${scr_sub} \
        ${env_nam} \
        ${threads} \
        ${infiles} \
        ${dir_out} \
        ${err_out} \
        ${nam_job}

#  Compress large stdout and stderr files, and remove files with size 0
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${err_out}"
```
</details>
<br />

<a id="c-run-the-initial-implementation-of-siq-chip"></a>
### C. Run the initial implementation of siQ-ChIP
<details>
<summary><i>Code: Run the initial implementation of siQ-ChIP</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node
grabnode  # Request 1 core, 20 GB memory, 1 day, no GPU

#  Define variables for directory paths, environment, submission script
#+ arguments, metadata, and so on
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
spe_typ="sc_bed"

dir_aln="${dir_pro}/align_${aligner}_${a_type}"
dir_bed="${dir_aln}/flag-${flg}_mapq-${mapq}/${spe_typ}"

depth=1
IFS=',' read -r -a arr_IP <<< "$(
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_bed}" \
        --pattern "IP_*.bed" \
        --depth ${depth}
)"
IFS=',' read -r -a arr_in <<< "$(
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_bed}" \
        --pattern "in_*.bed" \
        --depth ${depth}
)"

dir_ini="${dir_pro}/compute_coverage"
dir_cov="${dir_ini}/${aligner}_${a_type}_flag-${flg}_mapq-${mapq}"
typ_cov="siqchip"
dir_out="${dir_cov}/${typ_cov}"
dir_doc="${dir_out}/docs"
dir_log="${dir_out}/logs"

dir_raw="${dir_dat}/raw/docs"
fil_raw="measurements_siqchip_initial.tsv"
pth_raw="${dir_raw}/${fil_raw}"

dir_siq="${dir_bas}/siQ-ChIP"
fil_siq="get_siq.sh"
scr_siq="${dir_siq}/${fil_siq}"

#  Create output directory structure for trimmed FASTQ files and logs
mkdir -p ${dir_out}/{docs,logs}

if ${debug:-true}; then
    {
        echo "####################################"
        echo "## Hardcoded variable assignments ##"
        echo "####################################"
        echo ""
        echo "\${dir_bas}=${dir_bas}"
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
        echo "\${spe_typ}=${spe_typ}"
        echo ""
        echo "\${dir_aln}=${dir_aln}"
        echo "\${dir_bed}=${dir_bed}"
        echo ""
        echo "\${depth}=${depth}"
        echo "\${arr_IP[@]}=( ${arr_IP[*]} )"
        echo "\${arr_in[@]}=( ${arr_in[*]} )"
        echo ""
        echo "\${dir_ini}=${dir_ini}"
        echo "\${dir_cov}=${dir_cov}"
        echo "\${typ_cov}=${typ_cov}"
        echo "\${dir_out}=${dir_out}"
        echo "\${dir_doc}=${dir_doc}"
        echo "\${dir_log}=${dir_log}"
        echo ""
        echo "\${dir_raw}=${dir_raw}"
        echo "\${fil_raw}=${fil_raw}"
        echo "\${pth_raw}=${pth_raw}"
        echo ""
        echo "\${dir_siq}=${dir_siq}"
        echo "\${fil_siq}=${fil_siq}"
        echo "\${scr_siq}=${scr_siq}"
        echo ""
        echo ""
    }
fi

if ${need_to_do:-false}; then
    #  Create (temporary) decompressed versions of the BED files
    for file in ${dir_out}/*.bed.gz; do gunzip -k "${file}"; done
fi

if ${need_to_do:-false}; then
    #  Compute the average fragment length per sample
    for file in ${dir_out}/*.bed; do
        name="$(basename "${file}" ".sc.qnam.bed")"
        cat "${file}" \
            | awk -v name="${name}" '
                BEGIN {
                    OFS = "\t"
                } {
                    sum += $4
                    count++
                } END {
                    print name, (count ? sum / count : 0)
                }
        ' \
            > "${dir_doc}/"  #TODO: Pick up with this tomorrow...
    done
fi

if ${need_to_do:-false}; then
    #  Use awk to process a TSV file of per-sample siQ-ChIP measurements to
    #+ create parameter files for each sample
    cat "${pth_raw}" \
        | awk -v dir_doc="${dir_doc}" '
            BEGIN {
                FS = "\t"   # Set input field separator to tab
                OFS = "\n"  # Set output field separator to newline
            } NR == 1 {
                #  Store header values (column names) in an array for outfile
                #+ names
                for (i = 2; i <= NF; i++) {
                    file_names[i] = dir_doc "/params_" $i ".txt";
                }
                next  # Skip processing the header row
            } {
                #  Append rows 2 to 7 to respective files
                for (i = 2; i <= NF; i++) {
                    print $i > file_names[i]
                }
            }
        '
fi

if ${debug:-true}; then
    #  Debug the contents of the TXT parameters
    for file in ${dir_doc}/*.txt; do
        echo "## $(basename "${file}") ##"
        cat "${file}"
        echo ""
    done
fi

if ${debug:-true}; then
    #  Debug lists of the IP and input BED files, and the TXT parameter files
    for file in "${arr_IP[@]}"; do echo "${file}"; done && echo ""
    for file in "${arr_in[@]}"; do echo "${file}"; done && echo ""
    for file in ${dir_doc}/*.txt; do echo "${file}"; done
fi
```
</details>
<br />

