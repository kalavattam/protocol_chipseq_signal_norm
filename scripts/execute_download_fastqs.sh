#!/bin/bash

#  execute_download_fastqs.sh
#  KA


#  Run script in interactive/test mode (true) or command-line mode (false)
interactive=false

#  If not in interactive/test mode, then set script to exit if non-0 exit codes
#+ are encountered
if ! ${interactive}; then set -e; fi

#  Set the path to the "scripts" directory
# shellcheck disable=SC1091
if ${interactive}; then
    ## WARNING: Change path if you're not Kris and `interactive=true` ##
    dir_sc="${HOME}/tsukiyamalab/Kris/202X_protocol_ChIP/scripts"
else
    dir_sc="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi


#  Source and define functions ================================================
#  Set the path to the "functions" directory
dir_fn="${dir_sc}/functions"

# shellcheck disable=SC1091
{
    source "${dir_fn}/check_supplied_arg.sh"
    source "${dir_fn}/check_exists_file_dir.sh"
    source "${dir_fn}/check_format_time.sh"
    source "${dir_fn}/check_int_pos.sh"
    source "${dir_fn}/check_program_in_path.sh"
    source "${dir_fn}/echo_error.sh"
    source "${dir_fn}/echo_warning.sh"
    source "${dir_fn}/exit_0.sh"
    source "${dir_fn}/exit_1.sh"
}


#  Set up paths, values, and parameters for interactive mode
function set_interactive() {
    ## WARNING: Change the values if you're not Kris and `interactive=true` ##
    #  Set hardcoded paths, values, etc.
    dir_bas="${HOME}/tsukiyamalab/Kris"  # ls -lhaFG "${dir_bas}"
    nam_rep="202X_protocol_ChIP"         # echo "${nam_rep}"
    dir_rep="${dir_bas}/${nam_rep}"      # ls -lhaFG "${dir_rep}"
    nam_raw="data/raw"                   # echo "${nam_raw}"
    dir_raw="${dir_rep}/${nam_raw}"      # ls -lhaFG "${dir_raw}"
    nam_doc="docs"                       # echo "${nam_doc}"
    dir_doc="${dir_raw}/${nam_doc}"      # ls -lhaFG "${dir_doc}"
    nam_sym="data/symlinked"             # echo "${nam_sym}"
    dir_sym="${dir_rep}/${nam_sym}"      # ls -lhaFG "${dir_sym}"
    nam_log="logs"                       # echo "${nam_log}"
    dir_log="${dir_raw}/${nam_log}"      # ls -lhaFG "${dir_log}"
    nam_tsv="test_1.tsv"                 # echo "${nam_tsv}"
    pth_tsv="${dir_doc}/${nam_tsv}"      # ls -lhaFG "${pth_tsv}"

    #  Set hardcoded argument assignments
    # shellcheck disable=SC2269
    {
        verbose=true
        dry_run=true
        infile="${pth_tsv}"
        dir_out="${dir_raw}"
        dir_sym="${dir_sym}"
        par="slurm"
        nam_job="download_fastqs"
        err_out="${dir_log}"
        max_job=6
        time="3:00:00"
    }
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
sc_sub="${dir_sc}/submit_download_fastqs.sh"

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
infile=""
dir_out=""
dir_sym=""
par=""
nam_job="download_fastqs"
err_out=""
max_job=6
time="3:00:00"

#  Define help message
show_help=$(
    cat << EOM
Usage: 
  execute_download_fastqs.sh
    [--verbose] [--dry_run] --infile <str> --dir_out <str> --dir_sym <str>
    --par <str> --err_out <str> --nam_job <str> --max_job <int> --time <str>

Description:
  execute_download_fastqs.sh downloads FASTQ files listed in a TSV file and
  creates symbolic links to them using custom names provided in the TSV file.
  The script supports single- and paired-end sequenced reads, as well as
  downloading from both FTP and HTTPS addresses.

Arguments:
   -h, --help     Display this help message and exit.
   -v, --verbose  Run script in 'verbose' mode (optional).
  -dr, --dry_run  Run the command in 'check' mode (optional).
   -i, --infile   Input TSV file.
  -do, --dir_out  Output directory for downloaded FASTQ files
  -ds, --dir_sym  Output directory for symlinked FASTQ files.
   -p, --par      Program for parallelization ('slurm', 'gnu', or 'serial').
  -nj, --nam_job  The name of the job (default: ${nam_job}).
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (required; default: \${dir_out}/err_out).
  -mj, --max_job  The maximum number of jobs to run at one time (default: ${max_job};
                  ignored if par='gnu' or par='serial').
  -tm, --time     The length of time, in 'h:mm:ss' format, for the SLURM job
                  (default: ${time}; ignored if par='gnu' or par='serial').

Notes:
  - The script requires a properly formatted TSV (tab-separated value) file
    with a header and columns for run accession numbers, custom file names, and
    URLs (FTP or HTTPS).
  - For paired-end files, URLs in the TSV should be separated by semicolons.
  - Symbolic links are created in 'dir_sym' with names specified by the 
    'custom_name' column in the TSV file.
  - Use '--dry_run' to validate inputs without downloading files.

Example:
  \`\`\`
  #  Run with SLURM, allowing a maximum of six jobs to run concurrently
  bash execute_download_fastqs.sh
      --infile \${HOME}/path/to/PRJNA471802.tsv
      --dir_out \${HOME}/path/to/dir_downloaded_files
      --dir_sym \${HOME}/path/to/dir_symlinked_files
      --par slurm
      --max_job 6
  \`\`\`
EOM
)

#  Parse arguments
if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    exit_0
fi

if ${interactive}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -v|--verbose) verbose=true;   shift 1 ;;
            -dr|--dry_run) dry_run=true;   shift 1 ;;
             -i|--infile)  infile="${2}";  shift 2 ;;
            -do|--dir_out) dir_out="${2}"; shift 2 ;;
            -ds|--dir_sym) dir_sym="${2}"; shift 2 ;;
             -p|--par)     par="${2,,}";   shift 2 ;;
            -nj|--nam_job) nam_job="${2}"; shift 2 ;;
            -eo|--err_out) err_out="${2}"; shift 2 ;;
            -mj|--max_job) max_job="${2}"; shift 2 ;;
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
check_arg_supplied -a "${infile}" -n "infile"
check_exists_file_dir "f" "${infile}" "infile"

check_arg_supplied -a "${dir_out}" -n "dir_out"
check_exists_file_dir "d" "${dir_out}" "dir_out"

check_arg_supplied -a "${dir_sym}" -n "dir_sym"
check_exists_file_dir "d" "${dir_sym}" "dir_sym"

check_arg_supplied -a "${par}" -n "par"

check_arg_supplied -a "${nam_job}" -n "nam_job"

if [[ -n "${err_out}" ]]; then
    check_exists_file_dir "d" "${err_out}" "err_out"
elif [[ -z "${err_out}" ]]; then
    err_out="${dir_out}/err_out"
    check_exists_file_dir "d" "${err_out}" "err_out"
fi

case "${par}" in
    slurm)
        check_arg_supplied -a "${max_job}" -n "max_job"
        check_int_pos "${max_job}" "max_job"

        check_arg_supplied -a "${time}" -n "time"
        check_format_time "${time}"
        ;;

    serial)
        :
        ;;

    *)
        # Invalid option for --par
        echo_error \
            "Invalid selection for --par: '${par}'. Choose 'slurm', 'gnu'," \
            "or 'serial'."
        exit_1
        ;;
esac


#  Do the main work ===========================================================
#  Report argument variable assignments if in "verbose mode"
if ${verbose}; then
    echo "###################################"
    echo "## Argument variable assignments ##"
    echo "###################################"
    echo ""
    echo "verbose=${verbose}"
    echo "dry_run=${dry_run}"
    echo "infile=${infile}"
    echo "dir_out=${dir_out}"
    echo "dir_sym=${dir_sym}"
    echo "par=${par}"
    echo "nam_job=${nam_job}"
    echo "err_out=${err_out}"
    echo "max_job=${max_job}"
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
    (( iter++ ))
    
    #  Parse the header and detect available columns
    if [[ "${iter}" -eq 1 ]]; then
        IFS=$'\t' read -r -a headers <<< "${line}"

        #  Determine the index of the required columns dynamically
        for i in "${!headers[@]}"; do
            case "${headers[i]}" in
                "run_accession") run_acc_idx=${i} ;;
                "custom_name") custom_name_idx=${i} ;;
                "fastq_ftp") url_column_idx=${i}; url_column='fastq_ftp'   ;;
                "fastq_https") url_column_idx=${i}; url_column='fastq_https' ;;
            esac
        done
        
        #  Ensure the URL column was found
        if [[ -z "${url_column}" ]]; then
            echo_error "No valid URL column found in header."
            exit_1
        fi
        
        continue
    fi

    #  Read each column based on detected header indices
    run_acc=$(echo "${line}" | cut -f $(( run_acc_idx + 1 )))
    custom_name=$(echo "${line}" | cut -f $(( custom_name_idx + 1 )))
    urls=$(echo "${line}" | cut -f $(( url_column_idx + 1 )))

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

#  Limit job concurrency to the number of files, showing a warning if adjusted
case "${par}" in
    slurm)
        if [[ ${#list_acc[@]} -lt ${max_job} ]]; then
            max_job=${#list_acc[@]}
            echo_warning \
                "The maximum number of SLURM jobs to run at a time, ${max_job}," \
                "is greater than the number of FASTQ files, ${#list_acc[@]}." \
                "Adjusting max_job to ${#list_acc[@]}."
        fi
        ;;
esac

#  Download files and create symlinks
if [[ "${par}" == "slurm" ]]; then
    #  Join array elements with a delimiter (e.g., comma) to handle spaces correctly
    str_list_acc=$(IFS=','; echo "${list_acc[*]}")      # echo "${str_list_acc}"
    str_list_url_1=$(IFS=','; echo "${list_url_1[*]}")  # echo "${str_list_url_1}"
    str_list_url_2=$(IFS=','; echo "${list_url_2[*]}")  # echo "${str_list_url_2}"
    str_list_cus=$(IFS=','; echo "${list_cus[*]}")      # echo "${str_list_cus}"

    #  Submit each job via a SLURM job array for parallel downloads
    # shellcheck disable=SC2086
    sbatch \
        --job-name=${nam_job} \
        --nodes=1 \
        --cpus-per-task=1 \
        --time=${time} \
        --error=${err_out}/${nam_job}.%A-%a.stderr.txt \
        --output=${err_out}/${nam_job}.%A-%a.stdout.txt \
        --array=1-${#list_acc[@]}%${max_job} \
        --export=\
"str_list_acc=${str_list_acc},\
str_list_url_1=${str_list_url_1},\
str_list_url_2=${str_list_url_2},\
dir_out=${dir_out},\
dir_sym=${dir_sym},\
str_list_cus=${str_list_cus}" \
            ${sc_sub}
else
    #  Run serially
    for i in "${!list_acc[@]}"; do
        bash "${sc_sub}" \
            "${list_acc[i]}" \
            "${list_url_1[i]}" \
            "${list_url_2[i]}" \
            "${dir_out}" \
            "${dir_sym}" \
            "${list_cus[i]}" \
                 > "${err_out}/${nam_job}.${list_acc[i]}.stdout.txt" \
                2> "${err_out}/${nam_job}.${list_acc[i]}.stderr.txt"
    done
fi

# sra_explorer_prjna471802.sh = Swygert et al., Mol Cell 2019 = ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114566
# sra_explorer_prjna702745.sh = Swygert et al., Elife 2021 = ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167017
# sra_explorer_prjna549445.sh = Dickson et al., J Biol Chem 2020 = ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132906
# sra_explorer_prjna857063.sh = Dickson et al., Sci Rep 2023 = ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207783
