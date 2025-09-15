
Workflow: Validate New siQ-ChIP Implementation with the Original
================================================================

**Instructions to install and run the [original siQ-ChIP implementation](https://github.com/kalavattam/siQ-ChIP/tree/cerevisiae) by Brad Dickson. Outputs from this implementation (e.g., $\alpha$ values, normalized coverage, siQ-ChIP-scaled coverage) serve as validation standards to confirm the accuracy and reproducibility of the [Python siQ-ChIP implementation used in this protcol](https://github.com/kalavattam/protocol_chipseq_signal_norm).**

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
    1. [A. Generate BED files from QNAME-sorted BAM files.](#a-generate-bed-files-from-qname-sorted-bam-files)
    1. [B. Organize, generate metadata for, and then run the initial implementation of siQ-ChIP.](#b-organize-generate-metadata-for-and-then-run-the-initial-implementation-of-siq-chip)
1. [Test refactored bespoke code](#test-refactored-bespoke-code)
    1. [F. Compute normalized \(proportional\) or raw coverage.](#f-compute-normalized-proportional-or-raw-coverage)
    1. [G. Compute coverage with the sans spike-in quantitative ChIP-seq \(siQ-ChIP\) method.](#g-compute-coverage-with-the-sans-spike-in-quantitative-chip-seq-siq-chip-method)

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
branch="devel"

#  Go to base repository directory
cd "${dir_bas}" || echo_warning "cd'ing failed; check on this"

#  Clone forked siQ-ChIP repository
if [[ ! -d "${dir_rep}" ]]; then gh repo clone "${git_rep}"; fi

#  Checkout the 'devel' branch
cd "${dir_rep}" || echo_warning "cd'ing failed; check on this"

#  Switch to the specified branch
git checkout "${branch}"

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
<a id="a-generate-bed-files-from-qname-sorted-bam-files"></a>
### A. Generate BED files from QNAME-sorted BAM files.
<details>
<summary><i>Text: Generate BED files from QNAME-sorted BAM files.</i></summary>
<br />

`#TODO`
</details>
<br />

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
spe_typ="sc"

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
scr_cnv="${dir_scr}/compute_coverage.py"
env_nam="env_analyze"
dir_out="${dir_bam}_bed"
err_out="${dir_out}/logs"

day="$(date '+%Y-%m%d')"
nam_job="convert_bam_bed"
threads=1
time="2:00:00"
num_job=$(awk -F "," '{ print NF }' <<< "${infiles}")
max_job=6

#  Create output directory structure for trimmed FASTQ files and logs
mkdir -p ${dir_out}/{docs,logs}

#  Debug variable assignments
if [[ -f "${err_out}/${day}.execute.stdout.txt" ]]; then
    rm "${err_out}/${day}.execute."std???".txt"
fi

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
        echo "\${scr_cnv}=${scr_cnv}"
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
check_program_path python
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
        echo "        ${scr_cnv} \\"
        echo "        ${dir_out} \\"
        echo "        ${err_out} \\"
        echo "        ${nam_job}"
        echo ""
        echo ""
        echo "##############################################"
        echo "## Contents of BAM-to-BED conversion script ##"
        echo "##############################################"
        echo ""
        echo "## ${scr_cnv} ##"
        echo ""
        cat "${scr_cnv}"
        echo ""
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
        ${scr_cnv} \
        ${dir_out} \
        ${err_out} \
        ${nam_job}

#  Compress large stdout and stderr files, and remove files with size 0
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${err_out}"
```
</details>
<br />

<a id="b-organize-generate-metadata-for-and-then-run-the-initial-implementation-of-siq-chip"></a>
### B. Organize, generate metadata for, and then run the initial implementation of siQ-ChIP.
<details>
<summary><i>Bash code: Organize, generate metadata for, and then run the initial implementation of siQ-ChIP.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node --------------------------------------
# grabnode  # Uncomment to request 1 core,   20 GB memory, 1 day, no GPU
# grabnode  # Uncomment to request 8 cores, 160 GB memory, 1 day, no GPU


#  Define functions -----------------------------------------------------------
#  Populate an array with sorted file names that match a specified pattern
#+ (e.g., a Shell file glob), optionally excluding paths
function populate_array_glob() {
    local -n arr_ref=${1}     # Name of array to populate (passed by reference)
    local search="${2}"       # Directory to search
    local pattern="${3}"      # File name pattern (e.g., "IP_*.bed")
    local pth_ret=${4:-true}  # Retain paths (true/false); default is true
    
    mapfile -t arr_ref < <(
        find "${search}" -maxdepth 1 -type f -name "${pattern}" \
            | { 
                if ${pth_ret}; then 
                    cat                         # true: Retain path
                else 
                    awk -F '/' '{ print $NF }'  # false: Strip path
                fi 
            } \
            | sort
    )
}


#  Compute the average fragment length for a sample
function compute_frag_len_avg() {
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
export -f compute_frag_len_avg


#  Split a TSV file into per-sample parameter files; each sample gets its own
#+ TXT file containing rows of measurements
function split_tsv_params() {
    local infile="${1}"  # Path to input TSV file of siQ-ChIP measurements
    local outdir="${2}"  # Directory to save output parameter files

    #  Extract sample-specific measurements and write to separate files
    awk -v outdir="${outdir}" \
        'BEGIN {
            FS = "\t"   # Field separator: tab
            OFS = "\n"  # Output field separator: newline
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
function move_rename_files_glob() {
    local pattern="${1}"  # File search pattern (e.g., "Norm*.bed")
    local target="${2}"   # Target directory for moved files
    local renamer="${3}"  # 'sed' expression(s) for renaming files

    #  Expand the shell glob before processing
    local files=( ${pattern} )
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "Error: No files matching pattern '${pattern}' found." >&2
        return 1
    fi

    #  Process each file using xargs, passing 'target' and 'renamer' explicitly
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
function move_files_glob() {
    local pattern="${1}"  # File search pattern (e.g., "siqchip*.bed")
    local target="${2}"   # Target directory for moved files

    #  Expand the shell glob before processing
    local files=( ${pattern} )
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "Error: No files matching pattern '${pattern}' found." >&2
        return 1
    fi

    #  Process each file using xargs, passing 'target' explicitly
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


#  Function to compile Fortran scripts if necessary
function compile_fortran() {
    local fil_src="${1}"  # Path to Fortran source file
    local fil_bin="${2}"  # Path to output binary file

    if [[ ! -f "${fil_bin}" || "${fil_src}" -nt "${fil_bin}" ]]; then
        if ${dry_run:-false}; then
            echo "Dry run: Would compile '${fil_src}' -> '${fil_bin}'."
        else
            if ! \
                gfortran -O3 -fbounds-check -o "${fil_bin}" "${fil_src}"
            then
                echo "Error: Failed to compile '${fil_src}'." >&2
                return 1
            fi
        fi
    fi
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


#  Check that normalized coverage files sum to approximately unity
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


#  Define variables -----------------------------------------------------------
#  Define variables for directory paths, environment, submission script
#+ arguments, metadata, and so on
debug=true
siqchip=true
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
dir_cvg="${dir_ini}/${aligner}_${a_type}_flag-${flg}_mapq-${mapq}"
dir_exp="${dir_cvg}/siqchip"
dir_doc="${dir_exp}/docs"
dir_log="${dir_exp}/logs"
day="2025-0108"  ## WARNING: Change as needed ##
pth_doc="${dir_doc}/${day}.documentation.txt"

dir_raw="${dir_dat}/raw/docs"
fil_raw="measurements_siqchip_initial.tsv"
pth_raw="${dir_raw}/${fil_raw}"

dir_siq="${dir_bas}/siQ-ChIP"

dir_fas="${dir_gen}/cerevisiae/fasta/proc"
fil_fas="S288C_R64-5-1_proc.fasta.gz"
fil_chr="S288C_R64-5-1_proc.chrom-info.tsv"

siz_bin=5  # 10  # 1  # 30
siz_gen=12157105

submit="gnu"  # "slurm"  # "serial"
env_nam="env_siqchip"

#  Create output directory structure for trimmed FASTQ files and logs
mkdir -p ${dir_exp}/{docs,logs,norm,raw,siq}

#  To avoid repeated writing in testing, remove documentation TXT file
if [[ -f "${pth_doc}" ]]; then rm "${pth_doc}"; fi

if ${debug}; then
    {
        echo ""
        echo "#  siQ-ChIP experiment"
        echo "#  KA"
        echo "#  ${day}"
        echo ""
        echo ""
    } \
        >> >(tee -a "${pth_doc}")
fi  # cat "${pth_doc}"

#  Debug hardcoded variable assignments
if ${debug}; then
    {
        echo "####################################"
        echo "## Hardcoded variable assignments ##"
        echo "####################################"
        echo ""
        echo "\${debug}=${debug}"
        echo "\${siqchip}=${siqchip}"
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
        echo "\${dir_cvg}=${dir_cvg}"
        echo "\${dir_exp}=${dir_exp}"
        echo "\${dir_doc}=${dir_doc}"
        echo "\${dir_log}=${dir_log}"
        echo "\${day}=${day}"
        echo "\${pth_doc}=${pth_doc}"
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
        echo "\${siz_bin}=${siz_bin}"
        echo "\${siz_gen}=${siz_gen}"
        echo ""
        echo "\${submit}=${submit}"
        echo "\${env_nam}=${env_nam}"
        echo ""
        echo ""
    } \
        >> >(tee -a "${pth_doc}")
fi  # cat "${pth_doc}"


#  Source utility functions to activate environment and check dependencies
source "${dir_fnc}/check_program_path.sh"
source "${dir_fnc}/handle_env.sh"

#  Activate the required environment
handle_env "${env_nam}"

#  Ensure access to dependencies
arr_dep=(
    "bc"        # Precision calculator
    "bedtools"  # Program for genome arithmetic
    "gfortran"  # GNU Fortran compiler
    "gnuplot"   # Program for data visualization
    "parallel"  # GNU Parallel for parallel processing
    "sbatch"    #TODO: Description
)
for dep in "${arr_dep[@]}"; do check_program_path "${dep}"; done
unset arr_dep dep

#  If necessary, create (temporary) decompressed versions of the BED files
if compgen -G "${dir_bed}/*.bed.gz" > /dev/null; then
    for file in ${dir_bed}/*.bed.gz; do
        if [[ ! -f "${file%.gz}" ]]; then gunzip -k "${file}"; fi
    done
fi

#  Compute the average fragment length per sample
if [[ ! -f "${dir_doc}/frag_len_avg.txt" ]]; then
    if [[ ${threads} -gt 1 ]]; then
        #  Run in parallel if threads > 1
        parallel --jobs ${threads} \
            'compute_frag_len_avg {1} > {2}/{#}.tmp' \
                ::: "${dir_bed}"/*.bed \
                ::: "${dir_doc}"

        #  Avoid race conditions by individually writing to and then combining
        #+ temporary files
        cat "${dir_doc}"/*.tmp | sort -k1,1 > "${dir_doc}/frag_len_avg.txt" \
            && rm "${dir_doc}"/*.tmp
    else
        for file in "${dir_bed}"/*.bed; do
            compute_frag_len_avg "${file}" >> "${dir_doc}/frag_len_avg.txt"
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
    } \
        >> >(tee -a "${pth_doc}")
fi  # cat "${pth_doc}"

# #  Copy siQ-ChIP IP BED, input BED, and parameter TXT files into the siQ-ChIP
# #+ repository directory
# if ${siqchip}; then
#     if [[ ${threads} -gt 1 ]]; then
#         #  If threads > 1, use GNU Parallel to copy files in parallel
#         parallel --jobs ${threads} 'cp {1} {2}' \
#             ::: "${arr_IP[@]}" "${arr_in[@]}" "${arr_param[@]}" \
#             ::: "${dir_siq}"
#     else
#         #  Copy files sequentially if threads <= 1
#         cp "${arr_IP[*]}" "${arr_in[*]}" "${arr_param[*]}" "${dir_siq}"
#     fi
# fi

#  Assign arrays of path-free IP BED files, input BED files, parameter TXT
#+ files, and experiment stems
#  Populate arrays using the helper function 'populate_array_glob'
unset arr_IP arr_in arr_param
populate_array_glob arr_IP    "${dir_bed}" "IP_*.bed"     true  # false
populate_array_glob arr_in    "${dir_bed}" "in_*.bed"     true  # false
populate_array_glob arr_param "${dir_doc}" "params_*.txt" true  # false

unset arr_stem && typeset -a arr_stem+=( $(
    for file in "${arr_param[@]}"; do
        init="${file##*params_}"
        echo "${init%.txt}"
    done
) )
unset file init

#  Debug lists of the IP BED files, input BED files, and TXT parameter files
if ${debug}; then
    {
        echo "##################"
        echo "## IP BED files ##"
        echo "##################"
        echo ""
        for file in "${arr_IP[@]}"; do echo "${file}"; done
        echo ""
        echo ""

        echo "#####################"
        echo "## Input BED files ##"
        echo "#####################"
        echo ""
        for file in "${arr_in[@]}"; do echo "${file}"; done
        echo ""
        echo ""

        echo "##############################"
        echo "## siQ-ChIP parameter files ##"
        echo "##############################"
        echo ""
        for file in "${arr_param[@]}"; do echo "${file}"; done
        echo ""
        echo ""

        echo "######################"
        echo "## Experiment stems ##"
        echo "######################"
        echo ""
        for file in "${arr_stem[@]}"; do echo "${file}"; done
        echo ""
        echo ""
    } \
        >> >(tee -a "${pth_doc}")
fi  # cat "${pth_doc}"

if ${siqchip}; then
    #  Make experiment layout file for siQ-ChIP
    #TODO: Test files with paths (i.e., not in siQ repo directory)
    if [[ ! -f "${dir_doc}/${day}.layout_exp.txt" ]]; then
        echo "#getTracks: IP.bed input.bed params output_name" \
            >> "${dir_doc}/${day}.layout_exp.txt"

        for idx in ${!arr_IP[@]}; do
            echo \
                "${arr_IP[idx]} ${arr_in[idx]} ${arr_param[idx]}" \
                "${arr_stem[idx]}" \
                    >> "${dir_doc}/${day}.layout_exp.txt"
        done

        echo "#getResponse: CNTR.bed EXP.bed output_name" \
            >> "${dir_doc}/${day}.layout_exp.txt"
        echo "#getFracts: data any order, last is output_name" \
            >> "${dir_doc}/${day}.layout_exp.txt"
        echo "#END" >> "${dir_doc}/${day}.layout_exp.txt"
    fi  # cat "${dir_doc}/${day}.layout_exp.txt"

    #  Change to the siQ-ChIP repository directory and execute the siQ-ChIP
    #+ implementation
    cd "${dir_siq}" \
        || echo "Error: Failed to change to directory '${dir_siq}'."

    #  If necessary, compile required Fortran binaries 
    compile_fortran "${dir_siq}/tracks.f90" "${dir_siq}/tracks"
    compile_fortran "${dir_siq}/get_alpha.f90" "${dir_siq}/get_alpha"
    compile_fortran "${dir_siq}/merge_tracks.f90" "${dir_siq}/merge_tracks"

    #  If doing repeated testing, remove output and log files from experiment
    #+ directory
    if compgen -G "${dir_exp}/*.gz" > /dev/null; then
        rm "${dir_exp}/"*".gz"
    fi

    if compgen -G "${dir_log}/${day}.execute_siqchip.*.txt" > /dev/null; then
        rm "${dir_log}/${day}.execute_siqchip.stdout.txt"
        rm "${dir_log}/${day}.execute_siqchip.stderr.txt"
    fi

    #  Debug call to get_siq.sh
    if ${debug}; then
        {
            echo "######################################"
            echo "## Call to driver script get_siq.sh ##"
            echo "######################################"
            echo ""
            echo "bash ${dir_siq}/get_siq.sh \\"
            echo "    --verbose \\"
            echo "    --exp_lay ${dir_doc}/${day}.layout_exp.txt \\"
            echo "    --siz_bin ${siz_bin} \\"
            echo "    --siz_gen ${siz_gen} \\"
            echo "    --dir_out ${dir_exp} \\"
            echo "    --raw \\"
            echo "    --submit ${submit} \\"

            if [[ "${submit}" == "slurm" ]]; then
                echo "    --env_nam ${env_nam} \\"
            fi

            if [[ "${submit}" == "gnu" ]]; then
                echo "    --max_job ${SLURM_CPUS_ON_NODE} \\"
            fi

            if [[ "${submit}" != "gnu" ]]; then
                echo "         >> >(tee -a ${dir_log}/${day}.execute_siqchip.stdout.txt) \\"
                echo "        2>> >(tee -a ${dir_log}/${day}.execute_siqchip.stderr.txt)"
            fi

            echo ""
            echo ""
        } \
            >> >(tee -a "${pth_doc}")

        {
            echo "################################"
            echo "## Dry-run call to get_siq.sh ##"
            echo "################################"
            echo ""
            bash "${dir_siq}/get_siq.sh" \
                --dry_run \
                --verbose \
                --exp_lay "${dir_doc}/${day}.layout_exp.txt" \
                --siz_bin ${siz_bin} \
                --siz_gen ${siz_gen} \
                --dir_out "${dir_exp}" \
                --raw \
                --submit "${submit}" \
                $(
                    if [[ "${submit}" == "slurm" ]]; then
                        echo "--env_nam ${env_nam}"
                    fi
                ) \
                $(
                    if [[ "${submit}" == "gnu" ]]; then
                        echo "--max_job ${SLURM_CPUS_ON_NODE}"
                    fi
                )
            echo ""
        } \
            >> >(tee -a "${pth_doc}")
    fi  # cat "${pth_doc}"

    #  Run get_siq.sh
    bash "${dir_siq}/get_siq.sh" \
        --verbose \
        --exp_lay "${dir_doc}/${day}.layout_exp.txt" \
        --siz_bin ${siz_bin} \
        --siz_gen ${siz_gen} \
        --dir_out "${dir_exp}" \
        --raw \
        --submit "${submit}" \
        $(
            if [[ "${submit}" == "slurm" ]]; then
                echo "--env_nam ${env_nam}"
            fi
        ) \
        $(
            if [[ "${submit}" == "gnu" ]]; then
                echo "--max_job ${SLURM_CPUS_ON_NODE}"
            fi
        ) \
             >> >(tee -a "${dir_log}/${day}.execute_siqchip.stdout.txt")
            2>> >(tee -a "${dir_log}/${day}.execute_siqchip.stderr.txt")

    #  Clean up/organize siQ-ChIP outfiles: Move siQ-scaled, raw, and
    #+ normalized coverage bedGraph files, and metadata TXT files, to
    #+ respective directories
    cd "${dir_exp}" || echo "Failed to cd; check on this"
    move_files_glob "siq_*.bdg.gz" "${dir_exp}/siq"
    move_files_glob "norm_*.bdg.gz"  "${dir_exp}/norm"
    move_files_glob "raw_*.bdg.gz" "${dir_exp}/raw"
    move_files_glob "meta_*.txt" "${dir_exp}/logs"
fi
```
</details>
<br />

<a id="test-refactored-bespoke-code"></a>
## Test refactored bespoke code
<a id="f-compute-normalized-proportional-or-raw-coverage"></a>
### F. Compute normalized (proportional) or raw coverage.
<details>
<summary><i>Code: Compute normalized (proportional) or raw coverage.</i></summary>

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
dir_pro="${dir_dat}/processed"

#  Define alignment and coverage details
aligner="bowtie2"
a_type="global"
req_flg=true
flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
mapq=1
det_bam="flag-${flg}_mapq-${mapq}"
det_cvg="${aligner}_${a_type}_${det_bam}"
typ_cvg="norm"  ## WARNING: "raw" for unadjusted, "norm" for normalized ##

#  Further define directory setup
dir_aln="${dir_pro}/align_${aligner}_${a_type}"
dir_bam="${dir_aln}/${det_bam}/sc"
dir_cvg="${dir_pro}/compute_coverage"
dir_trk="${dir_cvg}/${det_cvg}/${typ_cvg}/tracks"

#  Define file search parameters
## WARNING: Change search parameters as needed ##
pattern="*.bam"
infiles="$(
    bash "${dir_scr}/find_files.sh" \
        --dir_fnd "${dir_bam}" \
        --pattern "${pattern}"
)"

#  Define script
exc_cvg="${dir_scr}/execute_compute_coverage.sh"

#  Define environment and script arguments, including resources
env_nam="env_analyze"
nam_job="compute_coverage_${typ_cvg}"
typ_out="bedgraph"  # "bigwig"
threads=8
bin_siz=30  # 1

#  Define log files
err_out="${dir_trk}/logs"
day="$(date '+%Y-%m%d')"
exc_pth="${dir_trk}/logs/${day}.execute.${nam_job}"


#  Create required directories if necessary -----------------------------------
mkdir -p ${dir_cvg}/${det_cvg}/norm/tracks/{docs,logs}


#  Activate the environment and check dependencies ----------------------------
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


#  Compute coverage -----------------------------------------------------------
bash "${exc_cvg}" \
    --verbose \
    --threads "${threads}" \
    --infiles "${infiles}" \
    --dir_out "${dir_trk}" \
    --typ_out "${typ_out}" \
    --bin_siz "1" \
    --typ_cvg "${typ_cvg}" \
    --err_out "${err_out}" \
    --nam_job "${nam_job}" \
    --slurm \
         >> >(tee -a "${exc_pth}.stdout.txt") \
        2>> >(tee -a "${exc_pth}.stderr.txt")

    # --bin_siz "${bin_siz}" \

#  Cleanup: Compress logs and remove empty files ------------------------------
bash "${dir_scr}/compress_remove_files.sh" --dir_fnd "${err_out}"

# ls -lhaFG "${err_out}"  ## Uncomment to check directory for logs ##
```
</details>
<br />

<a id="g-compute-coverage-with-the-sans-spike-in-quantitative-chip-seq-siq-chip-method"></a>
### G. Compute coverage with the sans spike-in quantitative ChIP-seq (siQ-ChIP) method.
<details>
<summary><i>Code: Compute coverage with the sans spike-in quantitative ChIP-seq (siQ-ChIP) method.</i></summary>

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
det_cvg="${aligner}_${a_type}_${det_bam}"
typ_cvg="alpha"

#  Further define directory setup
dir_aln="${dir_pro}/align_${aligner}_${a_type}"
dir_bam="${dir_aln}/${det_bam}/sc"
dir_cvg="${dir_pro}/compute_coverage"
dir_det="${dir_cvg}/${det_cvg}/${typ_cvg}"
dir_tbl="${dir_det}/tables"
dir_trk="${dir_det}/tracks"
eo_tbl="${dir_tbl}/logs"
eo_trk="${dir_trk}/logs"

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
scr_tbl="execute_calculate_scaling_factor_${typ_cvg}.sh"
scr_trk="execute_deeptools_coverage.sh"
fil_tbl="${dir_tbl}/IP_WT_G1-G2M-Q_Hho1-Hmo1_6336-6337_7750-7751.tsv"

#  Define log file prefixes
day="$(date '+%Y-%m%d')"
exc_tbl="${eo_tbl}/${day}.execute.${scr_tbl%.sh}.$(
    basename "${fil_tbl}" .tsv
)"
exc_trk="${eo_trk}/${day}.${scr_trk%.sh}"

#  Define environment and script arguments, including resources
env_nam="env_analyze"
threads=8
mes_tbl="${dir_raw}/docs/measurements_siqchip.tsv"
eqn=6  # "6nd"
bin_siz=1

#  Debug hardcoded variable assignments
if ${debug:-true}; then
    {
        echo "####################################"
        echo "## Hardcoded variable assignments ##"
        echo "####################################"
        echo ""
        echo "#  Define directory paths"
        echo "\${dir_bas}=${dir_bas}"
        echo "\${dir_rep}=${dir_rep}"
        echo "\${dir_scr}=${dir_scr}"
        echo "\${dir_fnc}=${dir_fnc}"
        echo "\${dir_dat}=${dir_dat}"
        echo "\${dir_raw}=${dir_raw}"
        echo "\${dir_pro}=${dir_pro}"
        echo ""
        echo "#  Define alignment and coverage details"
        echo "\${aligner}=${aligner}"
        echo "\${a_type}=${a_type}"
        echo "\${req_flg}=${req_flg}"
        echo "\${flg}=${flg}"
        echo "\${mapq}=${mapq}"
        echo "\${det_bam}=${det_bam}"
        echo "\${det_cvg}=${det_cvg}"
        echo "\${typ_cvg}=${typ_cvg}"
        echo ""
        echo "#  Further define directory setup"
        echo "\${dir_aln}=${dir_aln}"
        echo "\${dir_bam}=${dir_bam}"
        echo "\${dir_cvg}=${dir_cvg}"
        echo "\${dir_det}=${dir_det}"
        echo "\${dir_tbl}=${dir_tbl}"
        echo "\${dir_trk}=${dir_trk}"
        echo "\${eo_tbl}=${eo_tbl}"
        echo "\${eo_trk}=${eo_trk}"
        echo ""
        echo "#  Define environment, resources, and script arguments "
        echo "\${env_nam}=${env_nam}"
        echo "\${threads}=${threads}"
        echo "\${mes_tbl}=${mes_tbl}"
        echo "\${eqn}=${eqn}"
        echo "\${bin_siz}=${bin_siz}"
        echo ""
        echo "#  Define file search parameters"
        echo "\${pattern}=${pattern}"
        echo "\${include}=${include}"
        echo "\${exclude}=${exclude}"
        echo "\${infiles}=${infiles}"
        echo ""
        echo "#  Define scripts and output files"
        echo "\${scr_tbl}=${scr_tbl}"
        echo "\${scr_trk}=${scr_trk}"
        echo "\${fil_tbl}=${fil_tbl}"
        echo ""
        echo "#  Define log file prefixes"
        echo "\${day}=${day}"
        echo "\${exc_tbl}=${exc_tbl}"
        echo "\${exc_trk}=${exc_trk}"
        echo ""
        echo ""
    }
fi


#  Create required directories if necessary -----------------------------------
mkdir -p ${dir_tbl}/{docs,logs}
mkdir -p ${dir_trk}/{docs,logs}


#  Activate the environment and check dependencies ----------------------------
#  Source utility functions
source "${dir_fnc}/check_program_path.sh"
source "${dir_fnc}/echo_error.sh"
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
if [[ ! -f "${fil_tbl}" ]]; then
    #  Run the driver script to generate a TSV file of sample-specific siQ-ChIP
    #+ alpha scaling factors
    bash "${dir_scr}/execute_calculate_scaling_factor_${typ_cvg}.sh" \
        --verbose \
        --threads "${threads}" \
        --infiles "${infiles}" \
        --table "${mes_tbl}" \
        --eqn "${eqn}" \
        --outfile "${fil_tbl}" \
        --err_out "${eo_tbl}" \
        --flg_dep \
        --flg_len \
        --flg_mc \
        --slurm \
             >> >(tee -a "${exc_tbl}.stdout.txt") \
            2>> >(tee -a "${exc_tbl}.stderr.txt")

    if [[ $? -ne 0 ]]; then
        echo_error "The driver script for alpha computation failed."
    fi
fi

if [[ -f "${fil_tbl}" ]]; then
    #  Sort the table of scaling factors by rows
    awk 'NR == 1; NR > 1 && NF { print | "sort" }' "${fil_tbl}" \
        > "${dir_tbl}/tmp.tsv"

    #  Replace the original table with the sorted version
    mv -f "${dir_tbl}/tmp.tsv" "${fil_tbl}"

    # cat "${fil_tbl}"  ## Uncomment to check table contents ##
fi


#  Generate alpha-scaled signal tracks ----------------------------------------
#  Use the TSV file to generate alpha-scaled signal tracks
bash "${dir_scr}/execute_deeptools_coverage.sh" \
    --verbose \
    --threads "${threads}" \
    --table "${fil_tbl}" \
    --tbl_col "${typ_cvg}" \
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

<details>
<summary><i>Code: Scraps, 2025-0117</i></summary>

```bash
#!/bin/bash

#################################
## START: Material for testing ##
#################################
siz_gen=12157105
smp_tst="in_WT_G1_Hho1_6336"
dir_siq="${dir_cvg}/bowtie2_global_flag-2_mapq-1/siqchip"

function tally_aln_bam() {
    local bam="${1}"

    samtools view -F 4 "${bam}" \
        | awk '$2 == 99 || $2 == 1123 || $2 == 163 || $2 == 1187' \
        | wc -l
}


tally_aln_bam "${dir_bam}/${smp_tst}.sc.bam"  # 12851824
wc -l "${dir_bam}_bed/${smp_tst}.sc.bed" | awk '{ print $1 }'  # 12851824

python "${dir_scr}/compute_factor_depth.py" \
    -i "${dir_bam}/${smp_tst}.sc.bam" \
    -sb ${siz_bin} \
    -sg ${siz_gen}
    # -i "${dir_bam}_bed/${smp_tst}.sc.bed" \

gzip -dc "${dir_siq}/bin_10_slurm/raw/raw_${smp_tst}.bdg.gz" | head
gzip -dc "${dir_siq}/bin_10_slurm/raw/raw_${smp_tst}.bdg.gz" \
    | awk '
        NR == 1 {
            min = max = $4
        } {
            sum += $4; count++
            if ($4 < min) min = $4
            if ($4 > max) max = $4
        } END {
            if (count > 0) print "Min:", min, "Max:", max, "Avg:", sum / count
        }'
# Min: 5.6497175141242938E-003 Max: 2351.3904606954566 Avg: 10.7354

gzip -dc "${dir_siq}/bin_10_slurm/norm/norm_${smp_tst}.bdg.gz" | head
gzip -dc "${dir_siq}/bin_10_slurm/norm/norm_${smp_tst}.bdg.gz" \
    | awk '
        NR == 1 {
            min = max = $4
        } {
            sum += $4; count++
            if ($4 < min) min = $4
            if ($4 > max) max = $4
        } END {
            if (count > 0) print "Min:", min, "Max:", max, "Avg:", sum / count
        }'
#  Min: 4.39604e-10 Max: 0.000182962 Avg: 8.35317e-07

python "${dir_scr}/compute_factor_depth.py" \
    -i "${dir_bam}/${smp_tst}.sc.bam" \
    -bg "${dir_siq}/bin_10_slurm/norm/norm_${smp_tst}.bdg.gz" \
    -sb ${siz_bin} \
    -sg ${siz_gen}
# 10.58973100570373127027

python "${dir_scr}/compute_factor_depth.py" \
    -i "${dir_bam}/${smp_tst}.sc.bam" \
    -bg "${dir_siq}/bin_10_slurm/raw/raw_${smp_tst}.bdg.gz" \
    -sb ${siz_bin} \
    -sg ${siz_gen}
# 136097359.12745553255081176758

python "${dir_scr}/compute_factor_depth.py" \
    -i "${dir_bam}/${smp_tst}.sc.bam" \
    -m "frag" \
    -sb ${siz_bin} \
    -sg ${siz_gen}
# 10.57145971138664286570

gzip -dc "${dir_siq}/bin_10_slurm/raw/raw_${smp_tst}.bdg.gz" \
    | awk '
        NR == 1 {
            min = max = $4
        } {
            sum += $4; count++
            if ($4 < min) min = $4
            if ($4 > max) max = $4
        } END {
            if (count > 0) print "Sum:", sum
        }'
# Sum: 1.2874e+07

awk 'BEGIN { print 10.57145971138664286570 / 1.2874e+07 }'
# 8.21148e-07  # Not correct

awk 'BEGIN { print 10.57145971138664286570 / 12851824 }'
# 8.22565e-07  # Correct, I think?

#  Min: 5.6497175141242938E-003, 4.39604e-10
awk 'BEGIN { print 5.6497175141242938E-003 / 12851824 }'
# 4.39604e-10

#  Max: 2351.3904606954566, 0.000182962
awk 'BEGIN { print 2351.3904606954566 / 12851824 }'
# 0.000182962

#  Avg: 10.7354, 8.35317e-07
awk 'BEGIN { print 10.7354 / 12851824 }'
# 8.35321e-07

python "${dir_scr}/compute_factor_depth.py" \
    -i "${dir_bam}/${smp_tst}.sc.bam" \
    -m "frag" \
    -sb ${siz_bin} \
    -sg ${siz_gen}
# 10.57145971138664286570

python "${dir_scr}/compute_factor_depth.py" \
    -i "${dir_bam}/${smp_tst}.sc.bam" \
    -m "norm" \
    -sb ${siz_bin} \
    -sg ${siz_gen}
# 0.00000082256493019097

/usr/bin/time -v \
    python "${dir_scr}/compute_coverage_ratio.py" \
        -fp "${dir_nrm}/tracks/IP_WT_G1_Hho1_6336.sc.bedgraph" \
        -fn "${dir_nrm}/tracks/in_WT_G1_Hho1_6336.sc.bedgraph" \
        -fo "${dir_alf}/tracks/siq_WT_G1_Hho1_6336.sc.bedgraph" \
        -sf 0.0023106316806475436 \
        -de 0.00000082256493019097

less "${dir_nrm}/tracks/IP_WT_G1_Hho1_6336.sc.bedgraph"
# ...
# I       139020  139030  0.000000587195438972421814
# I       139030  139040  0.000000618895019319668699
# I       139040  139050  0.000000645281853688476867
# ...
# I       229140  229150  0.000000003042622662004834
# II      0       10      0.000000007616689624089380

less "${dir_nrm}/tracks/in_WT_G1_Hho1_6336.sc.bedgraph"
# I       139020  139030  0.000000829794333179596676
# I       139030  139040  0.000000775919333803290887
# I       139040  139050  0.000000759633062811069910
# ...
# I       229140  229150  0.000000006281134336350769
# I       229150  229160  0.000000004116929465424744
# I       229160  229170  0.000000004116929465424744
# I       229170  229180  0.000000004116929465424744
# I       229180  229190  0.000000004116929465424744
# I       229190  229200  0.000000004116929465424744
# I       229200  229210  0.000000004116929465424744
# I       229210  229220  0.000000004116929465424744
# I       229220  229230  0.000000004116929465424744
# I       229230  229240  0.000000004116929465424744
# I       229240  229250  0.000000004116929465424744
# I       229250  229260  0.000000004116929465424744
# I       229260  229270  0.000000004116929465424744
# I       229270  229280  0.000000002058464732712373
# I       230110  230120  0.000000001556199337930554
# I       230120  230130  0.000000007780996689652769
# I       230130  230140  0.000000007780996689652769
# I       230140  230150  0.000000007780996689652769
# I       230150  230160  0.000000007780996689652769
# I       230160  230170  0.000000007780996689652769
# I       230170  230180  0.000000007780996689652769
# I       230180  230190  0.000000007780996689652769
# I       230190  230200  0.000000007780996689652769
# I       230200  230210  0.000000007780996689652769
# I       230210  230220  0.000000007002897020687492
# II      0       10      0.000000037074198027690852

#  Try again after updating function to handle missing bins
/usr/bin/time -v \
    python "${dir_scr}/compute_coverage_ratio.py" \
        -fp "${dir_nrm}/tracks/IP_WT_G1_Hho1_6336.sc.bedgraph" \
        -fn "${dir_nrm}/tracks/in_WT_G1_Hho1_6336.sc.bedgraph" \
        -fo "${dir_alf}/tracks/siq_WT_G1_Hho1_6336.sc.bedgraph" \
        -sf 0.0023106316806475436 \
        -de 0.00000082256493019097

#  Comparing signal tracks from Brad and Kris's implementations in IGV:
#+ - Brad's: Max signal of...
#+     + 0.003155848  for I:28,298-61,871
#+     + 0.004347577  for I:167,388-191,608
#+     + 0.004283995  for II:191,406-204,401
#+     + 0.00886315   for XII:356,615-541,444
#+     + 0.004360575  for XIV:174,841-421,578
#+ - Kris's: Max signal of...
#+     + 0.0031328022 for I:28,298-61,871
#+     + 0.0043440936 for I:167,388-191,608
#+     + 0.0042771483 for II:191,406-204,401
#+     + 0.008773096  for XII:356,615-541,444
#+     + 0.0043579685 for XIV:174,841-421,578

###############################
## END: Material for testing ##
###############################
```
</details>
<br />
