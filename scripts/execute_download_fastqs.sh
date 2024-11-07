#!/bin/bash

#  execute_download_fastqs.sh
#  KA


#  Run script in interactive/test mode (true) or command-line mode (false)
interactive=false

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive}; then set -euo pipefail; fi

#  Set the path to the "scripts" directory
# shellcheck disable=SC1091
if ${interactive}; then
    ## WARNING: Change path if you're not Kris and `interactive=true` ##
    dir_scr="${HOME}/tsukiyamalab/Kris/202X_protocol_ChIP/scripts"
else
    dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi


#  Source and define functions ================================================
#  Set the path to the "functions" directory
dir_fnc="${dir_scr}/functions"

# shellcheck disable=SC1091
{
    source "${dir_fnc}/check_exists_file_dir.sh"
    source "${dir_fnc}/check_format_time.sh"
    source "${dir_fnc}/check_int_pos.sh"
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/check_supplied_arg.sh"
    source "${dir_fnc}/echo_error.sh"
    source "${dir_fnc}/echo_warning.sh"
    source "${dir_fnc}/exit_0.sh"
    source "${dir_fnc}/exit_1.sh"
    source "${dir_fnc}/handle_env.sh"
    source "${dir_fnc}/handle_env_activate.sh"
    source "${dir_fnc}/handle_env_deactivate.sh"
}


#  Set up paths, values, and parameters for interactive mode
function set_interactive() {
    ## WARNING: Change the values if you're not Kris and `interactive=true` ##
    #  Set hardcoded paths, values, etc.
    dir_rep="${HOME}/tsukiyamalab/Kris/202X_protocol_ChIP"
    dir_raw="${dir_rep}/data/raw"
    dir_doc="${dir_raw}/docs"
    dir_sym="${dir_rep}/data/symlinked"
    dir_log="${dir_raw}/logs"
    pth_tsv="${dir_doc}/test_3.tsv"

    #  Set hardcoded argument assignments
    # shellcheck disable=SC2269
    {
        verbose=true
        threads=4
        infile="${pth_tsv}"
        dir_out="${dir_raw}"
        dir_sym="${dir_sym}"
        nam_job="download_fastqs"
        err_out="${dir_log}"
        slurm=true
        time="3:00:00"
    }
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_align"
scr_sub="${dir_scr}/submit_download_fastqs.sh"

#  Initialize argument variables, assigning default values where applicable
verbose=false
threads=1
infile=""
dir_out=""
dir_sym=""
nam_job="download_fastqs"
err_out=""
slurm=false
time="3:00:00"

#  Define help message
show_help=$(
    cat << EOM
Usage: 
  execute_download_fastqs.sh
    [--verbose] --threads <int> --infile <str> --dir_out <str> --dir_sym <str>
    --nam_job <str> --err_out <str> --slurm --time <str>

Description:
  execute_download_fastqs.sh downloads FASTQ files listed in a TSV file and
  creates symbolic links to them using custom names provided in the TSV file.
  The script supports single- and paired-end sequenced reads, as well as
  downloading from both FTP and HTTPS addresses.

Arguments:
   -h, --help     Display this help message and exit.
   -v, --verbose  Run script in 'verbose' mode (optional).
   -t, --threads  Number of threads to use (default: ${threads}; see 'Notes'
                  below).
   -i, --infile   Input TSV file.
  -do, --dir_out  Output directory for downloaded FASTQ files
  -ds, --dir_sym  Output directory for symlinked FASTQ files.
  -nj, --nam_job  The name of the job (default: ${nam_job}).
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (default: \${dir_out}/err_out).
  -sl, --slurm    Submit jobs to the SLURM scheduler (optional; see 'Notes'
                  below).
  -tm, --time     The length of time (e.g., h:mm:ss) for the SLURM job
                  (required if --slurm is specified, ignored if not; default:
                  ${time}).

Dependencies:
  - Programs
    + Bash or Zsh
    + ln
    + SLURM (if --slurm is specified)
    + wget
  - Functions
    + check_exists_file_dir
    + check_format_time
    + check_int_pos
    + check_program_path
    + check_supplied_arg
    + echo_error
    + echo_warning
    + exit_0
    + exit_1
    + handle_env
    + handle_env_activate
    + handle_env_deactivate

Notes:
  - The script requires a properly formatted TSV (tab-separated value) file
    with a header and columns for run accession numbers, custom file names, and
    URLs (FTP or HTTPS). For paired-end files, URLs in the TSV should be
    separated by semicolons. See TSV files in 202X_protocol_ChIP/data/raw/docs
    for examples.
  - Symbolic links are created in 'dir_sym' with names specified by the 
    'custom_name' column in the TSV file.
  - If 'threads' is a positive integer greater than 1, the job submission
    script will be executed using GNU Parallel with the --jobs option set to
    \${threads}.
  - If the --slurm flag is specified and 'threads' is greater than 1, the job
    submission script will run via GNU Parallel within a SLURM job submission.
  - If --slurm is specified and 'threads' is set to 1, the execution script
    will terminate with an error, as serial job submission to SLURM is not
    permitted (and array job submission code has not been implemented).

Example:
  \`\`\`
  #  Run with SLURM, allowing a maximum of four downloads to run concurrently
  bash execute_download_fastqs.sh
      --threads 4
      --infile \${HOME}/path/to/PRJNA471802.tsv
      --dir_out \${HOME}/path/to/dir_downloaded_files
      --dir_sym \${HOME}/path/to/dir_symlinked_files
      --slurm
  \`\`\`
EOM
)

#  Parse arguments
if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    exit_0
fi

if ${interactive}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -v|--verbose) verbose=true;   shift 1 ;;
             -t|--threads) threads="${2}"; shift 2 ;;
             -i|--infile)  infile="${2}";  shift 2 ;;
            -do|--dir_out) dir_out="${2}"; shift 2 ;;
            -ds|--dir_sym) dir_sym="${2}"; shift 2 ;;
            -nj|--nam_job) nam_job="${2}"; shift 2 ;;
            -eo|--err_out) err_out="${2}"; shift 2 ;;
            -sl|--slurm)   slurm=true;     shift 1 ;;
            -tm|--time)    time="${2}";    shift 2 ;;
            *)
                echo "## Unknown argument passed: ${1} ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                exit_1
                ;;
        esac
    done
fi

#  Check arguments
check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

check_supplied_arg -a "${infile}" -n "infile"
check_exists_file_dir "f" "${infile}" "infile"

check_supplied_arg -a "${dir_out}" -n "dir_out"
check_exists_file_dir "d" "${dir_out}" "dir_out"

check_supplied_arg -a "${dir_sym}" -n "dir_sym"
check_exists_file_dir "d" "${dir_sym}" "dir_sym"

check_supplied_arg -a "${nam_job}" -n "nam_job"

if [[ -n "${err_out}" ]]; then
    check_exists_file_dir "d" "${err_out}" "err_out"
elif [[ -z "${err_out}" ]]; then
    err_out="${dir_out}/err_out"
    check_exists_file_dir "d" "${err_out}" "err_out"
fi

if "${slurm}"; then
    check_supplied_arg -a "${time}" -n "time"
    check_format_time "${time}"
fi

#  Activate environment and check that dependencies are in PATH
handle_env "${env_nam}"

check_program_path cut
check_program_path parallel
check_program_path wget


#  Do the main work ===========================================================
#  Report argument variable assignments if in "verbose mode"
if ${verbose}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "env_nam=${env_nam}"
    echo "scr_sub=${scr_sub}"
    echo ""
    echo ""
    echo "###################################"
    echo "## Argument variable assignments ##"
    echo "###################################"
    echo ""
    echo "verbose=${verbose}"
    echo "threads=${threads}"
    echo "infile=${infile}"
    echo "dir_out=${dir_out}"
    echo "dir_sym=${dir_sym}"
    echo "nam_job=${nam_job}"
    echo "err_out=${err_out}"
    echo "slurm=${slurm}"
    echo "time=${time}"
    echo ""
fi

#  Create arrays to store SRR entries, URLs, and custom names
unset      list_acc list_url_1 list_url_2 list_cus
typeset -a list_acc list_url_1 list_url_2 list_cus

#  Read the TSV file, processing each line to extract SRR accessions (if
#+ available), URLs, and custom names
iter=0
while IFS=$'\t' read -r line; do
    (( iter++ )) || true  # Prevent script exit if `set -e` 
    
    if ${verbose}; then echo "Processing line #${iter}: ${line}"; fi
    
    #  Parse the header and detect available columns
    if [[ "${iter}" -eq 1 ]]; then
        IFS=$'\t' read -r -a headers <<< "${line}"

        #  Determine the index of the required columns dynamically
        for i in "${!headers[@]}"; do
            case "${headers[i]}" in
                "run_accession") run_acc_idx=${i} ;;
                "custom_name") custom_name_idx=${i} ;;
                "fastq_ftp") url_col_idx=${i}; url_col='fastq_ftp'   ;;
                "fastq_https") url_col_idx=${i}; url_col='fastq_https' ;;
            esac
        done
        
        #  Ensure the URL column was found
        if [[ -z "${url_col}" ]]; then
            echo_error "No valid URL column found in header."
            exit_1
        fi
        
        continue
    fi

    #  Read each column based on detected header indices
    run_acc=$(echo "${line}" | cut -f $(( run_acc_idx + 1 )))
    custom_name=$(echo "${line}" | cut -f $(( custom_name_idx + 1 )))
    urls=$(echo "${line}" | cut -f $(( url_col_idx + 1 )))

    #  Handle missing run_acc entries
    if [[ -z "${run_acc}" || "${run_acc}" == "#N/A" ]]; then
        run_acc="SRR_undefined_${iter}"
    fi
    
    #  Parse the FASTQ URLs (paired-end or single-end)
    IFS=';' read -r -a fastq_urls <<< "${urls}"

    #  Add to arrays
    list_acc+=( "${run_acc}" )
    list_cus+=( "${custom_name}" )

    #  For single-end data
    if [[ ${#fastq_urls[@]} -eq 1 ]]; then
        list_url_1+=( "${fastq_urls[0]}" )
        list_url_2+=( "#N/A" )  # No second URL for single-end data
    #  For paired-end data
    elif [[ ${#fastq_urls[@]} -eq 2 ]]; then
        list_url_1+=( "${fastq_urls[0]}" )
        list_url_2+=( "${fastq_urls[1]}" )
    else
        echo_error "Unexpected number of FASTQ URLs for ${run_acc}"
        exit_1
    fi
done < "${infile}"

#  Report array element assignments if in "verbose mode"
if ${verbose}; then
    echo ""
    echo ""
    echo "#######################"
    echo "## list_acc elements ##"
    echo "#######################"
    echo ""
    for el in "${list_acc[@]}"; do echo "${el}"; done
    echo ""
    echo ""

    echo "#######################"
    echo "## list_cus elements ##"
    echo "#######################"
    echo ""
    for el in "${list_cus[@]}"; do echo "${el}"; done
    echo ""
    echo ""

    echo "#########################"
    echo "## list_url_1 elements ##"
    echo "#########################"
    echo ""
    for el in "${list_url_1[@]}"; do echo "${el}"; done
    echo ""
    echo ""

    echo "#########################"
    echo "## list_url_2 elements ##"
    echo "#########################"
    echo ""
    for el in "${list_url_2[@]}"; do echo "${el}"; done
    echo ""
    echo ""

    unset el
fi

#  Create a temporary GNU Parallel configuration file if threads > 1
if [[ "${threads}" -gt 1 ]]; then
    #  Generate a temporary file and check for issues
    # config=$(mktmp)
    config="${err_out}/${nam_job}.${RANDOM}.txt" ||
        {
            echo_error \
                "Failed to create a temporary file for GNU Parallel" \
                "configuration."
            exit_1
        }

    #  Ensure the temporary file is removed upon script exit
    trap 'rm -f "${config}"' EXIT

    #  Populate the GNU Parallel configuration file with parameters for each
    #+ job
    for i in "${!list_acc[@]}"; do
        echo \
            "${list_acc[i]}" \
            "${list_url_1[i]}" \
            "${list_url_2[i]}" \
            "${dir_out}" \
            "${dir_sym}" \
            "${list_cus[i]}" \
            "${err_out}" \
            "${nam_job}" \
                >> "${config}"
    done
fi

#  Execute download and symlink creation based on SLURM flag and thread count
if ${slurm}; then
    if [[ "${threads}" -gt 1 ]]; then
        #  Submit jobs to SLURM using GNU Parallel within a single SLURM job
        sbatch \
            --job-name="${nam_job}" \
            --nodes=1 \
            --cpus-per-task="${threads}" \
            --time="${time}" \
            --error="${config%.txt}.%A.stderr.txt" \
            --output="${config%.txt}.%A.stdout.txt" \
            --wrap="parallel --colsep ' ' --jobs ${threads} \
            'bash \"${scr_sub}\" {1} {2} {3} {4} {5} {6} {7} {8}' :::: ${config}"
    else
        #  Display error if SLURM submission is attempted with threads=1
        echo_error \
            "SLURM submissions require threads > 1; current value:" \
            "threads=${threads}."
        exit_1
    fi
else
    if [[ "${threads}" -gt 1 ]]; then
        #  Use GNU Parallel for multi-threaded execution without SLURM
        parallel --colsep ' ' --jobs "${threads}" \
            "bash \"${scr_sub}\" {1} {2} {3} {4} {5} {6} {7} {8}" \
            :::: "${config}"
    else
        #  Run jobs serially if threads=1
        for i in "${!list_acc[@]}"; do
            bash "${scr_sub}" \
                "${list_acc[i]}" \
                "${list_url_1[i]}" \
                "${list_url_2[i]}" \
                "${dir_out}" \
                "${dir_sym}" \
                "${list_cus[i]}" \
                "${err_out}" \
                "${nam_job}"
        done
    fi
fi

# sra_explorer_prjna471802.sh = Swygert et al., Mol Cell 2019 = ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114566
# sra_explorer_prjna702745.sh = Swygert et al., Elife 2021 = ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167017
# sra_explorer_prjna549445.sh = Dickson et al., J Biol Chem 2020 = ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132906
# sra_explorer_prjna857063.sh = Dickson et al., Sci Rep 2023 = ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207783
