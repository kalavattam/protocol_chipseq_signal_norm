#!/bin/bash

#  execute_trim_fastqs.sh
#  KA


#  Run script in interactive/test mode (true) or command-line mode (false)
interactive=true

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive}; then set -euo pipefail; fi

#  Set the path to the "scripts" directory
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
    source "${dir_fn}/check_exists_file_dir.sh"
    source "${dir_fn}/check_format_time.sh"
    source "${dir_fn}/check_int_pos.sh"
    source "${dir_fn}/check_program_path.sh"
    source "${dir_fn}/check_supplied_arg.sh"
    source "${dir_fn}/echo_error.sh"
    source "${dir_fn}/echo_warning.sh"
    source "${dir_fn}/exit_0.sh"
    source "${dir_fn}/exit_1.sh"
    source "${dir_fn}/handle_env_activate.sh"
    source "${dir_fn}/handle_env_deactivate.sh"
}


function set_interactive() {
    #  Set hardcoded paths, values, etc.
    ## WARNING: Change the values if you're not Kris and `interactive=true` ##
    dir_bas="${HOME}/tsukiyamalab/Kris"   # ls -lhaFG "${dir_bas}"
    nam_rep="202X_protocol_ChIP"          # echo "${nam_rep}"
    dir_rep="${dir_bas}/${nam_rep}"       # ls -lhaFG "${dir_rep}"
    nam_dat="data"                        # echo "${nam_dat}"
    dir_dat="${dir_rep}/${nam_dat}"       # ls -lhaFG "${dir_dat}"
    nam_sym="symlinked"                   # echo "${nam_sym}"
    dir_sym="${dir_dat}/${nam_sym}"       # ls -lhaFG "${dir_sym}"
    nam_pro="processed"                   # echo "${nam_pro}"
    dir_pro="${dir_dat}/${nam_pro}"       # ls -lhaFG "${dir_pro}"
    nam_trm="2024-1104_trim_atria_FASTQ"  # echo "${nam_trm}"
    dir_trm="${dir_pro}/${nam_trm}"       # ls -lhaFG "${dir_trm}"  # mkdir -p ${dir_trm}/{docs,logs}

    #  Set hardcoded argument assignments
    verbose=true                          # echo "${verbose}"
    dry_run=true                          # echo "${dry_run}"
    threads=8                             # echo "${threads}"
    infiles="$(
        bash "${dir_sc}/find_files.sh" \
            --dir_fnd "${dir_sym}" \
            --pattern "*.fastq.gz" \
            --depth 1 \
            --follow \
            --fastqs
    )"                                    # echo "${infiles}"
    dir_out="${dir_trm}"                  # ls -lhaFG "${dir_out}"
    sfx_se=".fastq.gz"                    # echo "${sfx_se}"
    sfx_pe="_R1.fastq.gz"                 # echo "${sfx_pe}"
    err_out="${dir_out}/logs"             # ls -lhaFG "${err_out}"
    nam_job="trim_fastqs"                 # echo "${nam_job}"
    slurm=true                            # echo "${slurm}"
    max_job=6                             # echo "${max_job}"
    time="0:45:00"                        # echo "${time}"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_analyze"                     # echo "${env_nam}"
sc_sub="${dir_sc}/submit_trim_fastqs.sh"  # ls -lhaFG "${sc_sub}"

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=8
infiles=""
dir_out=""
sfx_se=".fastq.gz"
sfx_pe="_R1.fastq.gz"
err_out=""
nam_job="trim_fastqs"
slurm=false
max_job=6
time="0:45:00"

show_help=$(
cat << EOM
execute_trim_fastqs.sh
  [--verbose] [--dry_run] --threads <int> --infiles <str> --dir_out <str>
  --sfx_se <str> --sfx_pe <str> --err_out <str> --nam_job <str> [--slurm]
  [--max_job <int>] [--time <str>]

Description:
  execute_trim_fastqs.sh performs read-adapter and quality trimming with the
  program Atria. It work with either single- or paired-end FASTQ files.

Arguments:
   -h, --help     Display this help message and exit.
   -v, --verbose  Run script in 'verbose' mode (optional).
  -dr, --dry_run  Run the command in check mode (optional).
   -t, --threads  Number of threads to use (required; default: ${threads}).
   -i, --infiles  Semicolon-separated string vector of FASTQ infiles
                  (required). Within the semicolon delimiters, sample-specific
                  corresponding paired-end files must be together separated by
                  commas, e.g.,
                  "\${HOME}/path/samp_1.fastq.gz;\${HOME}/path/samp_2_R1.fastq.gz,\${HOME}/path/samp_2_R2.fastq.gz;\${HOME}/path/samp_3.fastq.gz".
   -o, --dir_out  Directory for Atria-trimmed FASTQ outfile(s) (required).
  -ss, --sfx_se   Suffix to strip from single-end sequenced FASTQ files
                  (required; default: ${sfx_se}).
  -sp, --sfx_pe   Suffix to strip from paired-end sequenced FASTQ files
                  (required; default: ${sfx_pe}).
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (required; default: \${dir_out}/err_out).
  -nj, --nam_job  The name of the job, which is used when writing stderr and
                  stdout (required; default: ${nam_job}).
  -sl, --slurm    Submit jobs to the SLURM scheduler; otherwise, run them in
                  serial (optional).
  -mj, --max_job  The maximum number of jobs to run at one time (required if
                  --slurm is specified, ignored if not; default: ${max_job}).
  -tm, --time     The length of time (e.g., h:mm:ss) for the SLURM job
                  (required if --slurm is specified, ignored if not; default:
                  ${time}).

Dependencies:
  - Programs
    + Atria
    + Bash or Zsh
    + SLURM (if --slurm is specified)
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
    + handle_env_activate.sh
    + handle_env_deactivate.sh

Note:
  Atria is set to not allow read lengths less than 35 bp. It's also set to
  search for and trim known adapters, among other things. For more details, see
  the Atria documentation:
  github.com/cihga39871/Atria/blob/master/docs/2.Atria_trimming_methods_and_usages.md

Example:
  \`\`\`
  bash execute_trim_fastqs.sh
      --verbose
      --dry_run
      --threads 4
      --fq_1 "\${HOME}/path/samp_1_R1.fastq.gz"
      --fq_2 "\${HOME}/path/samp_1_R2.fastq.gz"
      --dir_out "\${HOME}/path/output/dir"
      --slurm
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
             -t|--threads) threads="${2}"; shift 2 ;;
             -i|--infiles) infiles="${2}"; shift 2 ;;
            -do|--dir_out) dir_out="${2}"; shift 2 ;;
            -ss|--sfx_se)  sfx_se="${2}";  shift 2 ;;
            -sp|--sfx_pe)  sfx_pe="${2}";  shift 2 ;;
            -eo|--err_out) err_out="${2}"; shift 2 ;;
            -nj|--nam_job) nam_job="${2}"; shift 2 ;;
            -sl|--slurm)   slurm=true;     shift 1 ;;
            -mj|--max_job) max_job="${2}"; shift 2 ;;
            -tm|--time)    time="${2}";    shift 2 ;;
            *) 
                echo "## Unknown parameter passed: ${1} ##" >&2
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

check_supplied_arg -a "${infiles}" -n "infiles"
check_exists_file_dir "d" "$(dirname "${infiles%%[,;]*}")" "infiles"

check_supplied_arg -a "${dir_out}" -n "dir_out"
check_exists_file_dir "d" "${dir_out}" "dir_out"

check_supplied_arg -a "${sfx_se}" -n "sfx_se"
check_supplied_arg -a "${sfx_pe}" -n "sfx_pe"

if [[ -n "${err_out}" ]]; then
    check_exists_file_dir "d" "${err_out}" "err_out"
elif [[ -z "${err_out}" ]]; then
    err_out="${dir_out}/err_out"
    check_exists_file_dir "d" "${err_out}" "err_out"
fi

check_supplied_arg -a "${nam_job}" -n "nam_job"

if ${slurm}; then
    check_supplied_arg -a "${sc_sub}" -n "sc_sub"
    check_exists_file_dir "f" "${sc_sub}" "sc_sub"

    check_supplied_arg -a "${max_job}" -n "max_job"
    check_int_pos "${max_job}" "max_job"
    
    check_supplied_arg -a "${time}" -n "time"
    check_format_time "${time}"
fi

#  Check dependencies
check_program_path atria
check_program_path check_supplied_arg
check_program_path check_exists_file_dir
check_program_path check_format_time
check_program_path check_int_pos
check_program_path check_program_path
check_program_path echo_error
check_program_path echo_warning
check_program_path exit_0
check_program_path exit_1


#  Parse the --infiles argument -----------------------------------------------
#+ ...into an array, then validate the infile value assignments
IFS=';' read -r -a arr_infiles <<< "${infiles}"  # unset arr_infiles
# for infile in "${arr_infiles[@]}"; do echo "${infile}"; done  # unset infile

#  Check that each infile exists; if not, exit
for infile in "${arr_infiles[@]}"; do
    # infile="${arr_infiles[0]}"  # echo "${infile}"
    
    #  Check that the infile contains a comma
    if [[ "${infile}" == *,* ]]; then
        fq_1="${infile%%,*}"  # echo "${fq_1}"
        fq_2="${infile#*,}"   # echo "${fq_2}"

        #  Ensure there is only one comma
        if [[ "${fq_2}" == *,* ]]; then
            echo_error \
                "It appears there is a problem with the serialized string" \
                "output from running find_files.sh in '--fastqs' mode. More" \
                "than one comma was found in a specific substring:" \
                "'${infile}'."
            exit_1
        fi

        #  Validate sample-specific paired-end FASTQ files exist
        check_exists_file_dir "f" "${fq_1}"
        check_exists_file_dir "f" "${fq_2}"

        #  Validate presence of file suffix in fq_1 assignment
        if [[ "${fq_1}" != *"${sfx_pe}" ]]; then
            echo_error \
                "Suffix '${sfx_pe}' not found in file name '${fq_1}'. Check" \
                "--sfx_pe assignment."
            exit_1
        fi
    else
        fq_1="${infile}"
        unset fq_2
        
        #  Validate sample-specific single-end FASTQ file exists
        check_exists_file_dir "f" "${fq_1}"

        #  Validate presence of file suffix in fq_1 assignment
        if [[ "${fq_1}" != *"${sfx_se}" ]]; then
            echo_error \
                "Suffix '${sfx_se}' not found in file name '${fq_1}'. Check" \
                "--sfx_se assignment."
            exit_1
        fi
    fi
done

#  Check that max_job is not greater than the number of elements in arr_infiles
if ${slurm}; then
    if [[ "${max_job}" -gt "${#arr_infiles[@]}" ]]; then
        echo_warning \
            "The maximum number of SLURM jobs to run at a time, ${max_job}," \
            "is greater than the number of infiles, ${#arr_infiles[@]}." \
            "Adjusting max_job from ${max_job} to ${#arr_infiles[@]}."
        max_job="${#arr_infiles[@]}"
    fi
fi


#  Do the main work ===========================================================
if ${verbose}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "env_nam=${env_nam}"
    echo "sc_sub=${sc_sub}"
    echo ""
    echo ""
    echo "###################################"
    echo "## Argument variable assignments ##"
    echo "###################################"
    echo ""
    echo "verbose=${verbose}"
    echo "dry_run=${dry_run}"
    echo "threads=${threads}"
    echo "infiles=${infiles}"
    echo "dir_out=${dir_out}"
    echo "sfx_se=${sfx_se}"
    echo "sfx_pe=${sfx_pe}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo "slurm=${slurm}"
    echo "sc_sub=${sc_sub}"
    echo "max_job=${max_job}"
    echo "time=${time}"
    echo ""
    echo ""
fi

# shellcheck disable=SC1083,SC2157
if ${slurm}; then
    #  If --slurm was specified, run jobs in parallel via individual job
    #+ submissions to SLURM
    if ${dry_run} || ${verbose}; then
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
        echo "    --array=1-${#arr_infiles[@]}%${max_job} \\"
        echo "    ${sc_sub} \\"
        echo "        ${env_nam} \\"
        echo "        ${threads} \\"
        echo "        ${infiles} \\"
        echo "        ${dir_out} \\"
        echo "        ${sfx_se} \\"
        echo "        ${sfx_pe} \\"
        echo "        ${err_out} \\"
        echo "        ${nam_job}"
        echo ""
        echo ""
        echo "#########################################"
        echo "## Contents of SLURM submission script ##"
        echo "#########################################"
        echo ""
        echo "## ${sc_sub} ##"
        echo ""
        cat "${sc_sub}"
        echo ""
    fi

    if ! ${dry_run}; then
        # shellcheck disable=SC2086
        sbatch \
            --job-name=${nam_job} \
            --nodes=1 \
            --cpus-per-task=${threads} \
            --time=${time} \
            --error=${err_out}/${nam_job}.%A-%a.stderr.txt \
            --output=${err_out}/${nam_job}.%A-%a.stdout.txt \
            --array=1-${#arr_infiles[@]}%${max_job} \
            ${sc_sub} \
                ${env_nam} \
                ${threads} \
                ${infiles} \
                ${dir_out} \
                ${sfx_se} \
                ${sfx_pe} \
                ${err_out} \
                ${nam_job}
    fi
else
    #  If --slurm was not specified, run jobs in serial
    if [[ "${CONDA_DEFAULT_ENV}" != "env_align" ]]; then
        eval "$(conda shell.bash hook)"
        conda activate env_align
    fi
    
    if ${dry_run} || ${verbose}; then
        echo "#############################"
        echo "## Serial call(s) to Atria ##"
        echo "#############################"
        echo ""
    fi
    
    for idx in "${!arr_infiles[@]}"; do
        n_call=$(( idx + 1 ))
        infile="${arr_infiles[idx]}"
        
        if [[ -z "${infile}" ]]; then
            echo_error \
                "Failed to retrieve infile for idx=${idx}:" \
                "\${arr_infiles[idx]}."
            exit_1
        fi

        if [[ "${infile}" == *,* ]]; then
            fq_1="${infile%%,*}"
            fq_2="${infile#*,}"
            samp="$(basename "${fq_1%%"${sfx_pe}"}")"
        else
            fq_1="${infile}"
            unset fq_2
            samp="$(basename "${fq_1%%"${sfx_se}"}")"
        fi

        if ${dry_run} || ${verbose}; then
            echo "## Call no. ${n_call} ##"
            echo "atria \\"
            echo "    -t ${threads} \\"
            echo "    -r ${fq_1} \\"
            echo "    $(if [[ -n ${fq_2} ]]; then echo "-R ${fq_2}"; fi) \\"
            echo "    -o ${dir_out} \\"
            echo "    --length-range 35:500 \\"
            echo "         > ${err_out}/${nam_job}.${samp}.stdout.txt \\"
            echo "        2> ${err_out}/${nam_job}.${samp}.stderr.txt"
            echo ""
        fi

        if ! ${dry_run}; then
            # shellcheck disable=SC2046,SC2086
            atria \
                -t ${threads} \
                -r ${fq_1} \
                $(if [[ -n ${fq_2} ]]; then echo "-R ${fq_2}"; fi) \
                -o ${dir_out} \
                --length-range 35:500 \
                     > ${err_out}/${nam_job}.${samp}.stdout.txt \
                    2> ${err_out}/${nam_job}.${samp}.stderr.txt
        fi
    done
fi
