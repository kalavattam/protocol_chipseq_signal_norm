
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
    1. [A. QNAME-sort BAM files filtered for *S. cerevisiae* alignments.](#a-qname-sort-bam-files-filtered-for-s-cerevisiae-alignments)
    1. [B. Generate BED files from QNAME-sorted BAM files.](#b-generate-bed-files-from-qname-sorted-bam-files)
    1. [C. Organize and generate metadata for and run the initial implementation of siQ-ChIP.](#c-organize-and-generate-metadata-for-and-run-the-initial-implementation-of-siq-chip)

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
<summary><i>Bash code: Clone the forked siQ-ChIP repository.</i></summary>

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
### A. QNAME-sort BAM files filtered for *S. cerevisiae* alignments.
<details>
<summary><i>Bash code: QNAME-sort BAM files filtered for </i>S. cerevisiae<i> alignments.</i></summary>

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
if ${debug}; then
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
if ${debug}; then
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
### B. Generate BED files from QNAME-sorted BAM files.
<details>
<summary><i>Bash code: Generate BED files from QNAME-sorted BAM files.</i></summary>

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
if ${debug}; then
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
if ${debug}; then
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

<a id="c-organize-and-generate-metadata-for-and-run-the-initial-implementation-of-siq-chip"></a>
### C. Organize and generate metadata for and run the initial implementation of siQ-ChIP.
<details>
<summary><i>Bash code: Organize and generate metadata for and run the initial implementation of siQ-ChIP.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  # Uncomment to request 1 core,   20 GB memory, 1 day, no GPU
# grabnode  # Uncomment to request 8 cores, 160 GB memory, 1 day, no GPU


#  Define functions -----------------------------------------------------------
#  Populate an array with sorted file names that match a specified pattern
#+ (e.g., a file glob), optionally excluding paths
function populate_array_glob() {
    local -n arr_ref=${1}  # Name of array to populate (passed by reference)
    local search="${2}"    # Directory to search
    local pattern="${3}"   # File name pattern (e.g., "IP_*.bed")
    local path=${4:-true}  # Retain paths (true/false); default is true
    
    mapfile -t arr_ref < <(
        find "${search}" -maxdepth 1 -type f -name "${pattern}" \
            | { 
                if ${path}; then 
                    cat                         # Retain path in file name
                else 
                    awk -F '/' '{ print $NF }'  # Strip path from file name
                fi 
            } \
            | sort
    )
}


#  Compute the average fragment length for a sample
function compute_avg_frag_len() {
    local infile="${1}"
    local strip="${2:-".sc.qnam.bed"}"
    local name

    #  Extract sample name by removing the specified suffix from the file name
    name="$(basename "${infile}" "${strip}")"

    #  Compute the average fragment length (column 4) for the BED file
    awk \
        -v name="${name}" \
        'BEGIN {
            OFS = "\t"
        } {                                        # For each line, sum lengths
            sum += $4; count++                     # ...and increment count
        } END {                                    # Compute, print average
            print name, (count ? sum / count : 0)  # Avoid division by zero
        }' \
        "${infile}"
}
export -f compute_avg_frag_len


#  Split a TSV file into per-sample parameter files; each sample gets its own
#+ TXT file containing rows of measurements
function split_tsv_params() {
    local infile="${1}"  # Path to input TSV file of siQ-ChIP measurements
    local outdir="${2}"  # Directory to save output parameter files

    #  Extract sample-specific measurements and write to separate files
    awk -v outdir="${outdir}" \
        'BEGIN {
            FS = "\t"   # Set field separator to tab
            OFS = "\n"  # Set output field separator to newline
        } NR == 1 {
            #  From the header row (row 1), create outfile paths for each
            #+ sample
            for (i = 2; i <= NF; i++) {
                file_names[i] = outdir "/params_" $i ".txt";
            }
            next  # Skip the header row during data processing
        } {
            #  Write rows of measurements (2-7) to corresponding sample files
            for (i = 2; i <= NF; i++) {
                print $i > file_names[i]
            }
        }' \
        "${infile}"
}


#  Rename and move globbed files
function move_rename_globbed_files() {
    local pattern="${1}"  # File search pattern (e.g., "Norm*.bed")
    local target="${2}"   # Target directory for moved files
    local renamer="${3}"  # 'sed' expression(s) for renaming files

    #  Expand the shell glob before processing
    local files=( ${pattern} )
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "Error: No files matching pattern '${pattern}' found." >&2
        return 1
    fi

    #  Process each file using xargs and pass 'target' and 'renamer' explicitly
    printf '%s\n' "${files[@]}" \
        | sort \
        | env \
            target="${target}" renamer="${renamer}" \
            xargs -I {} bash -c '
                #  Generate new name and validate it
                nam_new=$(echo "{}" | sed -e "${renamer}")
                if [[ -z "${nam_new}" || "${nam_new}" == "{}" ]]; then
                    echo "Error: Failed to assign a valid new name for file: {}" >&2
                    exit 1
                fi

                #  Attempt to move the file
                mv "{}" "${target}/${nam_new}" || {
                    echo "Error: Failed to move {} → ${target}/${nam_new}" >&2
                    exit 1
                }
            '
}


#  Move globbed files without renaming
function move_globbed_files() {
    local pattern="${1}"  # File search pattern (e.g., "siqchip*.bed")
    local target="${2}"   # Target directory for moved files

    #  Expand the shell glob before processing
    local files=( ${pattern} )
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "Error: No files matching pattern '${pattern}' found." >&2
        return 1
    fi

    #  Process each file using xargs and pass 'target' explicitly
    printf '%s\n' "${files[@]}" \
        | sort \
        | env target="${target}" xargs -I {} bash -c '
            #  Attempt to move the file
            mv "{}" "${target}/{}" || {
                echo "Error: Failed to move {} → ${target}/{}" >&2
                exit 1
            }
        '
}


#  Convert a BED file to a headered bedgraph file
function convert_bed_bedgraph() {
    local pair="${1}"     # Key-value pair: "sample.bed;rgb_color"
    local dir_out="${2}"  # Output directory for bedgraph file

    #  Extract BED file name and color information from pair
    local bed=$(awk -F ';' '{ print $1 }' <<< "${pair}")
    local color=$(awk -F ';' '{ print $2 }' <<< "${pair}")
    
    #  Define bedgraph output path
    local bg="$(basename ${bed} ".bed").bedGraph"
    local pth_out="${dir_out}/${bg}"

    #  Add track definition line to bedgraph file
    echo "track type=bedgraph name=\"${bed##*/}\" color=${color}" \
        > "${pth_out}"

    #  Append BED data as bedgraph format, assuming BED columns are 1-4
    awk 'BEGIN { OFS="\t" } { print $1, $2, $3, $4 }' "${bed}" >> "${pth_out}"
}
export -f convert_bed_bedgraph


#TODO: Description
function covert_bedgraph_bigwig() {
    :  #TODO
}


#  Define variables -----------------------------------------------------------
#  Define variables for directory paths, environment, submission script
#+ arguments, metadata, and so on
debug=true
siqchip=false
convert=true
threads=${SLURM_CPUS_ON_NODE:-1}

dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_scr="${dir_rep}/scripts"
dir_fnc="${dir_scr}/functions"
dir_dat="${dir_rep}/data"
dir_pro="${dir_dat}/processed"
dir_gen="${dir_dat}/genomes"

aligner="bowtie2"
a_type="global"
req_flg=true
flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
mapq=1
spe_typ="sc_bed"

dir_aln="${dir_pro}/align_${aligner}_${a_type}"
dir_bed="${dir_aln}/flag-${flg}_mapq-${mapq}/${spe_typ}"

unset arr_IP arr_in
populate_array_glob arr_IP "${dir_bed}" "IP_*.bed" true
populate_array_glob arr_in "${dir_bed}" "in_*.bed" true

dir_ini="${dir_pro}/compute_coverage"
dir_cov="${dir_ini}/${aligner}_${a_type}_flag-${flg}_mapq-${mapq}"
dir_exp="${dir_cov}/siqchip"
dir_doc="${dir_exp}/docs"
dir_log="${dir_exp}/logs"
day="2024-1217"  ## WARNING: Change as needed ##

dir_raw="${dir_dat}/raw/docs"
fil_raw="measurements_siqchip_initial.tsv"
pth_raw="${dir_raw}/${fil_raw}"

dir_siq="${dir_bas}/siQ-ChIP"

dir_fas="${dir_gen}/cerevisiae/fasta/proc"
fil_fas="S288C_R64-5-1_proc.fasta.gz"
fil_chr="S288C_R64-5-1_proc.chrom-info.tsv"

#  Create output directory structure for trimmed FASTQ files and logs
mkdir -p ${dir_exp}/{docs,logs}
mkdir -p ${dir_exp}/{alpha,norm}/{bed,bg,bw}

#  Debug hardcoded variable assignments
if ${debug}; then
    {
        echo "####################################"
        echo "## Hardcoded variable assignments ##"
        echo "####################################"
        echo ""
        echo "\${debug}=${debug}"
        echo "\${siqchip}=${siqchip}"
        echo "\${convert}=${convert}"
        echo "\${threads}=${threads}"
        echo ""
        echo "\${dir_bas}=${dir_bas}"
        echo "\${dir_rep}=${dir_rep}"
        echo "\${dir_scr}=${dir_scr}"
        echo "\${dir_fnc}=${dir_fnc}"
        echo "\${dir_dat}=${dir_dat}"
        echo "\${dir_pro}=${dir_pro}"
        echo "\${dir_gen}=${dir_gen}"
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
        echo "\${arr_IP[@]}=( ${arr_IP[*]} )"
        echo "\${arr_in[@]}=( ${arr_in[*]} )"
        echo ""
        echo "\${dir_ini}=${dir_ini}"
        echo "\${dir_cov}=${dir_cov}"
        echo "\${dir_exp}=${dir_exp}"
        echo "\${dir_doc}=${dir_doc}"
        echo "\${dir_log}=${dir_log}"
        echo "\${day}=${day}"
        echo ""
        echo "\${dir_raw}=${dir_raw}"
        echo "\${fil_raw}=${fil_raw}"
        echo "\${pth_raw}=${pth_raw}"
        echo ""
        echo "\${dir_siq}=${dir_siq}"
        echo ""
        echo "\${dir_fas}=${dir_fas}"
        echo "\${fil_fas}=${fil_fas}"
        echo "\${fil_chr}=${fil_chr}"
        echo ""
        echo ""
    } # \
        # >> #TODO
fi

#  Source utility functions to activate environment and check dependencies
source "${dir_fnc}/check_program_path.sh"
source "${dir_fnc}/handle_env.sh"

#  Activate the required environment
handle_env env_siqchip

#  Ensure access to dependencies
arr_dep=(
    "bc"                # Precision calculator
    "bedClip"           # Program for clipping BED file entries
    "bedtools"          # Program for genome arithmetic
    "bedGraphToBigWig"  # Convert bedgraph files to bigWig format
    "gfortran"          # GNU Fortran compiler
    "gnuplot"           # Program for data visualization
    "parallel"          # GNU Parallel for parallel processing
)
for dep in "${arr_dep[@]}"; do check_program_path "${dep}"; done
unset arr_dep dep

#  Create (temporary) decompressed versions of the BED files
for file in ${dir_bed}/*.bed.gz; do
    if [[ ! -f "${file%.gz}" ]]; then gunzip -k "${file}"; fi
done

#  Compute the average fragment length per sample
if [[ ! -f "${dir_doc}/avg_frag_len.txt" ]]; then
    if [[ ${threads} -gt 1 ]]; then
        #  Run in parallel if threads > 1
        parallel --jobs ${threads} \
            'compute_avg_frag_len {1} > {2}/{#}.tmp' \
                ::: "${dir_bed}"/*.bed \
                ::: "${dir_doc}"

        #  Avoid race conditions by individually writing to and then combining
        #+ temporary files
        cat "${dir_doc}"/*.tmp | sort -k1,1 > "${dir_doc}/avg_frag_len.txt" \
            && rm "${dir_doc}"/*.tmp
    else
        for file in "${dir_bed}"/*.bed; do
            compute_avg_frag_len "${file}" >> "${dir_doc}/avg_frag_len.txt"
        done
    fi
fi

#  Check: Do all files already exist in the output directory? First, define
#+ array of expected parameter files for all samples
arr_param=(
    "params_WT_G1_Hho1_6336.txt"
    "params_WT_G1_Hho1_6337.txt"
    "params_WT_G1_Hmo1_7750.txt"
    "params_WT_G1_Hmo1_7751.txt"
    "params_WT_G2M_Hho1_6336.txt"
    "params_WT_G2M_Hho1_6337.txt"
    "params_WT_G2M_Hmo1_7750.txt"
    "params_WT_G2M_Hmo1_7751.txt"
    "params_WT_Q_Hho1_6336.txt"
    "params_WT_Q_Hho1_6337.txt"
    "params_WT_Q_Hmo1_7750.txt"
    "params_WT_Q_Hmo1_7751.txt"
)

#  Perform the check
all_exist=true
for file in "${arr_param[@]}"; do
    if [[ ! -f "${dir_doc}/${file}" ]]; then
        all_exist=false
        break
    fi
done

#  If any file is missing, split a TSV file of siQ-ChIP measurements into
#+ individual parameter TXT files, one for each sample
if ! ${all_exist}; then split_tsv_params "${pth_raw}" "${dir_doc}"; fi

# #  Less strict check verification of any 'params_*.txt' files exist in the
# #+ documentation directory
# if ! compgen -G "${dir_doc}/params_*.txt" > /dev/null; then
#     split_tsv_params "${pth_raw}" "${dir_doc}"
# fi

#  Unset variables used to check and generate parameter TXT files
unset arr_param file all_exist

#  Debug the contents of the parameter TXT files
if ${debug}; then
    {
        echo "#################################"
        echo "## Parameter TXT file contents ##"
        echo "#################################"
        echo ""

        for file in ${dir_doc}/params*.txt; do
            echo "## $(basename "${file}") ##"
            cat "${file}"
            echo ""
        done
        
        echo ""
    } # \
        # >> #TODO
fi

#  Assign array of parameter TXT files
unset arr_param
populate_array_glob arr_param "${dir_doc}" "params_*.txt" true

#  Copy siQ-ChIP IP BED, input BED, and parameter TXT files into the siQ-ChIP
#+ repository directory
if ${siqchip}; then
    if [[ ${threads} -gt 1 ]]; then
        #  If threads > 1, use GNU Parallel to copy files in parallel
        parallel --jobs ${threads} 'cp {1} {2}' \
            ::: "${arr_IP[@]}" "${arr_in[@]}" "${arr_param[@]}" \
            ::: "${dir_siq}"
    else
        #  Copy files sequentially if threads <= 1
        cp "${arr_IP[*]}" "${arr_in[*]}" "${arr_param[*]}" "${dir_siq}"
    fi
fi

#  Assign arrays of path-free IP BED files, input BED files, parameter TXT
#+ files, and experiment stems
#  Populate arrays using the helper function
unset arr_IP arr_in arr_param
populate_array_glob arr_IP    "${dir_bed}" "IP_*.bed"     false
populate_array_glob arr_in    "${dir_bed}" "in_*.bed"     false
populate_array_glob arr_param "${dir_doc}" "params_*.txt" false

unset arr_stem && typeset -a arr_stem+=( $(
    for file in "${arr_param[@]}"; do
        init="${file##*params_}"
        stem="siqchip_${init%.txt}"
        echo "${stem}"
    done    
) )
unset file init stem

#  Debug lists of the IP BED files, input BED files, and TXT parameter files
if ${debug}; then
    {
        echo "##############"
        echo "## IP files ##"
        echo "##############"
        for file in "${arr_IP[@]}"; do echo "${file}"; done
        echo ""
        echo ""

        echo "#################"
        echo "## Input files ##"
        echo "#################"
        for file in "${arr_in[@]}"; do echo "${file}"; done
        echo ""
        echo ""

        echo "##############################"
        echo "## siQ-ChIP parameter files ##"
        echo "##############################"
        for file in "${arr_param[@]}"; do echo "${file}"; done
        echo ""
        echo ""

        echo "######################"
        echo "## Experiment stems ##"
        echo "######################"
        for file in "${arr_stem[@]}"; do echo "${file}"; done
        echo ""
        echo ""
    } # \
        # >> #TODO
fi

if ${siqchip}; then
    #  Make experiment layout file, 'EXPlayout', for siQ-ChIP
    if [[ ! -f "${dir_siq}/EXPlayout" ]]; then
        echo "#getTracks: IP.bed input.bed params output_name" \
            >> "${dir_siq}/EXPlayout"

        for idx in ${!arr_IP[@]}; do
            echo \
                "${arr_IP[idx]} ${arr_in[idx]} ${arr_param[idx]}" \
                "${arr_stem[idx]}" \
                    >> "${dir_siq}/EXPlayout"
        done

        echo "#getResponse: CNTR.bed EXP.bed output_name" \
            >> "${dir_siq}/EXPlayout"
        echo "#getFracts: data any order, last is output_name" \
            >> "${dir_siq}/EXPlayout"
        echo "#END" >> "${dir_siq}/EXPlayout"
    fi

    #  An 'Annotations.bed' file is required in the siQ-ChIP repository directory
    if [[ ! -f "${dir_siq}/Annotations.bed" ]]; then
        touch "${dir_siq}/Annotations.bed"
    fi

    #  Change to the siQ-ChIP repository directory and execute the siQ-ChIP
    #+ implementation
    cd "${dir_siq}" \
        || echo "Error: Failed to change to directory '${dir_siq}'."
    bash getsiq.sh \
         >> >(tee -a "${dir_log}/${day}.execute_siqchip.stdout.txt") \
        2>> >(tee -a "${dir_log}/${day}.execute_siqchip.stderr.txt")

    #  Clean up/organize siQ-ChIP outfiles
    #  Step 1: Remove copied BED and parameter files
    rm IP*.bed in*.bed params*.txt

    #  Step 2: Rename and move normalized coverage BED files to 'norm/bed'
    #+ directory
    move_rename_globbed_files "Norm*.bed" "${dir_exp}/norm/bed" "
        s/^NormCovIP_siqchip_/norm_IP_/;
        s/^NormCovIN_siqchip_/norm_in_/
    "

    #  Step 3: Move siQ-ChIP-scaled coverage BED files to 'alpha/bed' directory
    move_globbed_files "siqchip*.bed" "${dir_exp}/alpha/bed"

    #  Step 4: Move siQ-ChIP alpha calculation TXT files to 'alpha' directory
    move_globbed_files "siqchip*.alpha" "${dir_exp}/alpha"

    #  Step 5: Move the 'a.out' file to 'log' directory
    mv "a.out" "${dir_log}"

    #  Step 6: Move all other recent files to experiment 'docs' directory
    ## WARNING: Change date and time as needed ##
    find * -type f -newermt "2024-12-17 15:25:00" \
        | sort \
        | env dir_doc="${dir_doc}" \
            xargs -I {} bash -c 'mv "{}" "${dir_doc}/{}"'
fi

if ${convert}; then
    #  Define an array containing key-value pairs for BED files and (IBM)
    #+ colors
    arr_bed=(
        #  red 40; #ff8389; 255,131,137
        "${dir_exp}/norm/bed/norm_in_WT_G1_Hho1_6336.bed;255,131,137"

        #  red 30; #ffb3b8; 255,179,184
        "${dir_exp}/norm/bed/norm_in_WT_G1_Hho1_6337.bed;255,179,184"

        #  magenta 40; #ff7eb6; 255,126,182
        "${dir_exp}/norm/bed/norm_in_WT_G2M_Hho1_6336.bed;255,126,182"

        #  magenta 30; #ffafd2; 255,175,210
        "${dir_exp}/norm/bed/norm_in_WT_G2M_Hho1_6337.bed;255,175,210"

        #  purple 40; #be95ff; 190,149,255
        "${dir_exp}/norm/bed/norm_in_WT_Q_Hho1_6336.bed;190,149,255"

        #  purple 30; #d4bbff; 212,187,255
        "${dir_exp}/norm/bed/norm_in_WT_Q_Hho1_6337.bed;212,187,255"

        #  red 70; #da1e28; 162,25,31
        "${dir_exp}/norm/bed/norm_IP_WT_G1_Hho1_6336.bed;162,25,31" 

        #  red 60; #a2191f; 218,30,40
        "${dir_exp}/norm/bed/norm_IP_WT_G1_Hho1_6337.bed;218,30,40" 

        #  magenta 70; #9f1853; 159,24,83
        "${dir_exp}/norm/bed/norm_IP_WT_G2M_Hho1_6336.bed;159,24,83" 

        #  magenta 60; #d02670; 208,38,112
        "${dir_exp}/norm/bed/norm_IP_WT_G2M_Hho1_6337.bed;208,38,112"

        #  purple 70; #6929c4; 105,41,196
        "${dir_exp}/norm/bed/norm_IP_WT_Q_Hho1_6336.bed;105,41,196"

        #  purple 60; #8a3ffc; 138,63,252
        "${dir_exp}/norm/bed/norm_IP_WT_Q_Hho1_6337.bed;138,63,252"

        #  blue 40; #78a9ff; 120,169,255
        "${dir_exp}/norm/bed/norm_in_WT_G1_Hmo1_7750.bed;120,169,255"

        #  blue 30; #a6c8ff; 166,200,255
        "${dir_exp}/norm/bed/norm_in_WT_G1_Hmo1_7751.bed;166,200,255"

        #  cyan 40; #33b1ff; 51,177,255
        "${dir_exp}/norm/bed/norm_in_WT_G2M_Hmo1_7750.bed;51,177,255"

        #  cyan 30; #82cfff; 130,207,255
        "${dir_exp}/norm/bed/norm_in_WT_G2M_Hmo1_7751.bed;130,207,255"

        #  teal 40; #08bdba; 8,189,186
        "${dir_exp}/norm/bed/norm_in_WT_Q_Hmo1_7750.bed;8,189,186"

        #  teal 30; #3ddbd9; 61,219,217
        "${dir_exp}/norm/bed/norm_in_WT_Q_Hmo1_7751.bed;61,219,217"

        #  blue 70; #0043ce; 0,67,206
        "${dir_exp}/norm/bed/norm_IP_WT_G1_Hmo1_7750.bed;0,67,206"

        #  blue 60; #0f62fe; 15,98,254
        "${dir_exp}/norm/bed/norm_IP_WT_G1_Hmo1_7751.bed;15,98,254"

        #  cyan 70; #00539a; 0,83,154
        "${dir_exp}/norm/bed/norm_IP_WT_G2M_Hmo1_7750.bed;0,83,154"

        #  cyan 60; #0072c3; 0,114,195
        "${dir_exp}/norm/bed/norm_IP_WT_G2M_Hmo1_7751.bed;0,114,195"

        #  teal 70; #005d5d; 0,93,93
        "${dir_exp}/norm/bed/norm_IP_WT_Q_Hmo1_7750.bed;0,93,93" 

        #  teal 60; #007d79; 0,125,121
        "${dir_exp}/norm/bed/norm_IP_WT_Q_Hmo1_7751.bed;0,125,121"
    )

    #  Debug the array key (BED file)-value (color) pairs
    {
        echo "########################################"
        echo "## siQ-ChIP BED-color key-value pairs ##"
        echo "########################################"
        echo ""

        for pair in "${arr_bed[@]}"; do
              key="$(echo "${pair}" | awk -F ';' '{ print $1 }')"
            value="$(echo "${pair}" | awk -F ';' '{ print $2 }')"
            
            echo "  key .......... ${key}"
            echo "value .......... ${value}"
            echo ""
        done
        echo ""
    } # \
        # >> #TODO

    #  Perform BED-to-bedgraph file conversions
    if compgen -G "${dir_exp}/norm/bg/*.bedGraph" > /dev/null; then
        echo \
            "BED-to-bedGraph conversion skipped: Files already exist in" \
            "'${dir_exp}/norm/bg'"
    else
        if [[ ${threads} -gt 1 ]]; then
            #  Run in parallel if threads > 1
            parallel --jobs ${threads} \
                'convert_bed_bedgraph {1} {2}' \
                    ::: "${pair}" \
                    ::: "${dir_exp}/norm/bg"
        else
            #  Convert files sequentially if threads <= 1
            for pair in "${arr_bed[@]}"; do
                convert_bed_bedgraph "${pair}" "${dir_exp}/norm/bg"
            done
        fi
    fi

    #  Create arrays of bedGraph infiles and bigWig outfiles  #TODO #PICKUPHERE
    populate_array_glob arr_bg "${dir_exp}/norm/bg" "*.bedGraph" true
    # for bg in "${arr_bg[@]}"; do head "${bg}"; done

    #  Create chromosome info TSV file if it does not exist
    if [[ ! -f "${dir_fas}/${fil_chr}" ]]; then
        if [[ -f "${dir_fas}/${fil_fas}" ]]; then
            zcat "${dir_fas}/${fil_fas}" \
                | awk '
                    /^>/ {
                        #  Print name and length of previous chromosome
                        if (seq) { print name"\t"length(seq); seq="" }

                        #  Extract chromosome name
                        name = substr($0, 2)
                    } index($0, ">") != 1 {
                        #  Accumulate sequence for length calculation
                        seq = seq $0
                    } END {
                        #  Print name and length of final chromosome
                        if (seq) { print name"\t"length(seq) }
                    }
                ' \
                    > "${dir_fas}/${fil_chr}"
        else
            echo "Error: FASTA file not found: '${dir_fas}/${fil_fas}'" >&2
        fi
    fi

    if [[ ${threads} -gt 1 ]]; then
        #  Run parallelized version of pipeline
        #  Step 1: Exclude headers from bedGraph files
        parallel --jobs ${threads} '
            if [[ -f {.}.bedGraph && ! -f {.}.noheader.bedGraph ]]; then
                if ! tail -n +2 {.}.bedGraph > {.}.noheader.bedGraph; then
                    echo \
                        "Error: Failed to remove headers from file:" \
                        "{.}.bedGraph. Expected outfile:" \
                        "{.}.noheader.bedGraph." >&2
                    exit 1
                fi
            else
                echo \
                    "Error: Missing infile {.}.bedGraph or outfile" \
                    "{.}.noheader.bedGraph already exists." >&2
                exit 1
            fi
        ' \
            ::: "${arr_bg[@]}"

        #  Step 2: Sort the bedGraph files
        parallel --jobs ${threads} '
            if [[ -f {.}.noheader.bedGraph && ! -f {.}.sorted.bedGraph ]]; then
                if ! \
                    LC_ALL=C sort -k1,1 -k2,2n {.}.noheader.bedGraph \
                        > {.}.sorted.bedGraph
                then
                    echo \
                        "Error: Failed to sort file: {.}.noheader.bedGraph." \
                        "Expected outfile: {.}.sorted.bedGraph." >&2
                    exit 1
                fi
            else
                echo \
                    "Error: Missing infile {.}.noheader.bedGraph or outfile" \
                    "{.}.sorted.bedGraph already exists." >&2
                exit 1
            fi
        ' \
            ::: "${arr_bg[@]}"

        #  Step 3: Clip records spanning chromosome ends
        parallel --jobs ${threads} '
            if [[
                     -f {1.}.sorted.bedGraph \
                &&   -f {2} \
                && ! -f {1.}.clipped.bedGraph
            ]]; then
                if ! \
                    bedClip -verbose=2 \
                        {1.}.sorted.bedGraph \
                        {2} \
                        {1.}.clipped.bedGraph
                then
                    echo \
                        "Error: Failed to clip file: {1.}.sorted.bedGraph" \
                        "using chromosome information from {2}. Expected" \
                        "outfile: {1.}.clipped.bedGraph." >&2
                    exit 1
                fi
            else
                echo \
                    "Error: Missing infiles {1.}.sorted.bedGraph or" \
                    "chromosome information file {2}, or outfile" \
                    "{1.}.clipped.bedGraph already exists." >&2
                exit 1
            fi
        ' \
            ::: "${arr_bg[@]}" \
            ::: "${dir_fas}/${fil_chr}"
        #TODO: Need to capture verbose output from bedClip

        #  Convert the above-processed bedgraph files to bigWig files
        parallel --jobs ${threads} '
            if [[
                     -f {1.}.clipped.bedGraph \
                &&   -f {2} \
                && ! -f {3}/{1/.}.bigWig
            ]]; then
                if ! \
                    bedGraphToBigWig \
                        {1.}.clipped.bedGraph \
                        {2} \
                        {3}/{1/.}.bigWig
                then
                    echo \
                        "Error: Failed to convert file {1.}.clipped.bedGraph" \
                        "to bigWig using chromosome information from {2}." \
                        "Expected outfile: {3}/{1/.}.bigWig." >&2
                    exit 1
                fi
            else
                echo \
                    "Error: Missing infiles {1.}.clipped.bedGraph or" \
                    "chromosome information file {2}, or outfile" \
                    "{3}/{1/.}.bigWig already exists." >&2
                exit 1
            fi
        ' \
            ::: "${arr_bg[@]}" \
            ::: "${dir_fas}/${fil_chr}" \
            ::: "${dir_exp}/norm/bw"
    else
        #  Run serial version of the pipeline
        for file in "${arr_bg[@]}"; do
            base="${file%%.bedGraph}"

            #  Step 1: Exclude headers from bedGraph files
            if [[
                     -f "${base}.bedGraph" \
                && ! -f "${base}.noheader.bedGraph"
            ]]; then
                if ! \
                    tail -n +2 "${base}.bedGraph" > "${base}.noheader.bedGraph"
                then
                    echo \
                        "Error: Failed to remove headers from file:" \
                        "${base}.bedGraph. Expected outfile:" \
                        "${base}.noheader.bedGraph." >&2
                    exit 1
                fi
            else
                echo \
                    "Error: Missing infile ${base}.bedGraph or outfile" \
                    "${base}.noheader.bedGraph already exists." >&2
                exit 1
            fi

            #  Step 2: Sort the bedGraph files
            if [[
                     -f "${base}.noheader.bedGraph" \
                && ! -f "${base}.sorted.bedGraph"
            ]]; then
                if ! \
                    LC_ALL=C sort -k1,1 -k2,2n "${base}.noheader.bedGraph" \
                        > "${base}.sorted.bedGraph"
                then
                    echo \
                        "Error: Failed to sort file:" \
                        "${base}.noheader.bedGraph. Expected outfile:" \
                        "${base}.sorted.bedGraph." >&2
                    exit 1
                fi
            else
                echo \
                    "Error: Missing infile ${base}.noheader.bedGraph or" \
                    "outfile ${base}.sorted.bedGraph already exists." >&2
                exit 1
            fi

            #  Step 3: Clip records spanning chromosome ends
            if [[
                     -f "${base}.sorted.bedGraph" \
                &&   -f "${dir_fas}/${fil_chr}" \
                && ! -f "${base}.clipped.bedGraph"
            ]]; then
                if ! \
                    bedClip -verbose=2 \
                        "${base}.sorted.bedGraph" \
                        "${dir_fas}/${fil_chr}" \
                        "${base}.clipped.bedGraph"
                then
                    echo \
                        "Error: Failed to clip file: ${base}.sorted.bedGraph" \
                        "using chromosome information from" \
                        "${dir_fas}/${fil_chr}. Expected outfile:" \
                        "${base}.clipped.bedGraph." >&2
                    exit 1
                fi
            else
                echo \
                    "Error: Missing infile ${base}.sorted.bedGraph or" \
                    "chromosome information file ${dir_fas}/${fil_chr}, or" \
                    "outfile ${base}.clipped.bedGraph already exists." >&2
                exit 1
            fi

            #  Step 4: Convert clipped bedGraph files to bigWig format
            if [[
                     -f "${base}.clipped.bedGraph" \
                &&   -f "${dir_fas}/${fil_chr}" \
                && ! -f "${dir_exp}/norm/bw/${base##*/}.bigWig"
            ]]; then
                if ! \
                    bedGraphToBigWig \
                        "${base}.clipped.bedGraph" \
                        "${dir_fas}/${fil_chr}" \
                        "${dir_exp}/norm/bw/${base##*/}.bigWig"
                then
                    echo \
                        "Error: Failed to convert file" \
                        "${base}.clipped.bedGraph to bigWig using chromosome" \
                        "information from ${dir_fas}/${fil_chr}. Expected" \
                        "outfile: ${dir_exp}/norm/bw/${base##*/}.bigWig." >&2
                    exit 1
                fi
            else
                echo \
                    "Error: Missing infile ${base}.clipped.bedGraph or" \
                    "chromosome information file ${dir_fas}/${fil_chr}, or" \
                    "outfile ${dir_exp}/norm/bw/${base##*/}.bigWig already" \
                    "exists." >&2
                exit 1
            fi
        done
    fi

    #  Clean up intermediate files
    rm *.{noheader,sorted,clipped}.bedGraph
fi
```
</details>
<br />

