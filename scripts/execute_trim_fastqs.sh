#!/bin/bash

#  execute_trim_fastqs.sh
#  KA


#  Run script in interactive/test mode (true) or command-line mode (false)
interactive=false

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive}; then set -euo pipefail; fi

#  Set the path to the "scripts" directory
if ${interactive}; then
    ## WARNING: If interactive=true, change path as needed ##
    dir_scr="${HOME}/repos/protocol_chipseq_signal_norm/scripts"
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
    source "${dir_fnc}/debug_array_contents.sh"
    source "${dir_fnc}/echo_error.sh"
    source "${dir_fnc}/echo_warning.sh"
    source "${dir_fnc}/exit_0.sh"
    source "${dir_fnc}/exit_1.sh"
    source "${dir_fnc}/handle_env.sh"
    source "${dir_fnc}/reset_max_job.sh"
    source "${dir_fnc}/set_params_parallel.sh"
}


#  Set up paths, values, and parameters for interactive mode
function set_interactive() {
    #  Set base repository paths
    dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
    dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
    dir_scr="${dir_rep}/scripts"
    dir_fnc="${dir_scr}/functions"

    #  Define data directories
    dir_dat="${dir_rep}/data"
    dir_sym="${dir_dat}/symlinked"
    dir_pro="${dir_dat}/processed"
    dir_trm="${dir_pro}/trim_fastqs"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=true
    threads=8
    infiles="$(  ## WARNING: Change search parameters as needed ##
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_sym}" \
            --pattern "*.fastq.gz" \
            --include "*Q_Hmo1*" \
            --depth 1 \
            --follow \
            --fastqs
    )"
    dir_out="${dir_trm}"
    sfx_se=".fastq.gz"
    sfx_pe="_R1.fastq.gz"
    err_out="${dir_out}/logs"
    nam_job="trim_fastqs"
    slurm=false
    max_job=2
    time="0:45:00"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_protocol"
scr_sub="${dir_scr}/submit_trim_fastqs.sh"
par_job=""

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
max_job=8
time="0:45:00"

show_help=$(cat << EOM
Usage:
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
  -mj, --max_job  Maximum number of jobs to run concurrently (default: ${max_job}).
                    - Required if '--slurm' is specified; controls SLURM array
                      tasks.
                    - If '--slurm' is not specified:
                      + If 'max_job' is greater than 1, jobs run in parallel
                        via GNU Parallel.
                      + If 'max_job' is 1, jobs run sequentially (serial mode).
  -sl, --slurm    Submit jobs to the SLURM scheduler; otherwise, run them in
                  serial (optional).
  -tm, --time     The length of time (e.g., h:mm:ss) for the SLURM job
                  (required if --slurm is specified, ignored if not; default:
                  ${time}).

Dependencies:
  - Programs
    + Atria
    + Bash or Zsh
    + pbzip2
    + pigz
    + SLURM (if --slurm is specified)
  - Functions
    + check_exists_file_dir
    + check_format_time
    + check_int_pos
    + check_program_path
    + check_supplied_arg
    + debug_array_contents
    + echo_error
    + echo_warning
    + exit_0
    + exit_1
    + handle_env
    + handle_env_activate
    + handle_env_deactivate

Note:
  Atria is set to not allow read lengths less than 35 bp. It is also set to
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
    check_supplied_arg -a "${scr_sub}" -n "scr_sub"
    check_exists_file_dir "f" "${scr_sub}" "scr_sub"

    check_supplied_arg -a "${max_job}" -n "max_job"
    check_int_pos "${max_job}" "max_job"
    
    check_supplied_arg -a "${time}" -n "time"
    check_format_time "${time}"
else
    read -r threads max_job par_job < <(
        set_params_parallel "${threads}" "${max_job}" "${par_job}"
    )
    # echo "${threads}"
    # echo "${max_job}"
    # echo "${par_job}"

    check_supplied_arg -a "${par_job}" -n "par_job"
    check_int_pos "${par_job}" "par_job"
fi


#  Activate environment and check that dependencies are in PATH ---------------
handle_env "${env_nam}" > /dev/null

check_program_path atria
if ! ${slurm} && [[ ${par_job} -gt 1 ]]; then check_program_path parallel; fi
check_program_path pbzip2
check_program_path pigz


#  Parse the --infiles argument -----------------------------------------------
#+ ...into an array, then validate the infile value assignments
IFS=';' read -r -a arr_infiles <<< "${infiles}"
# for infile in "${arr_infiles[@]}"; do echo "${infile}"; done  # unset infile

#  Check that each infile exists; if not, exit
for infile in "${arr_infiles[@]}"; do
    # infile="${arr_infiles[1]}"  # echo "${infile}"
    
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
        check_exists_file_dir "f" "${fq_1}" "fq_1"
        check_exists_file_dir "f" "${fq_2}" "fq_2"

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
        check_exists_file_dir "f" "${fq_1}" "fq_1"

        #  Validate presence of file suffix in fq_1 assignment
        if [[ "${fq_1}" != *"${sfx_se}" ]]; then
            echo_error \
                "Suffix '${sfx_se}' not found in file name '${fq_1}'. Check" \
                "--sfx_se assignment."
            exit_1
        fi
    fi
done

if ${slurm}; then
    max_job=$(reset_max_job "${max_job}" "${#arr_infiles[@]}")
fi

if ${verbose}; then
    echo "#######################################"
    echo "## Parsed vector for parallelization ##"
    echo "#######################################"
    echo ""
    debug_array_contents "arr_infiles"
    if ${slurm}; then
        echo "  - Max no. jobs to run at a time (SLURM): ${max_job}"
    elif [[ ${par_job} -gt 1 ]]; then
        echo "  - Max no. jobs to run at a time (GNU Parallel): ${par_job}"
    else
        echo "  - Max no. jobs to run at a time: ${par_job}"
    fi
    echo "  - No. threads per job: ${threads}"
    echo ""
    echo ""
fi


#  Do the main work ===========================================================
if ${verbose}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "env_nam=${env_nam}"
    echo "scr_sub=${scr_sub}"
    echo "par_job=${par_job:-#N/A}"
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
    echo "scr_sub=${scr_sub}"
    echo "max_job=${max_job}"
    echo "time=${time}"
    echo ""
    echo ""
    echo "#################################"
    echo "## Array derived from variable ##"
    echo "#################################"
    echo ""
    echo "arr_infiles=( ${arr_infiles[*]} )"
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
        echo "    ${scr_sub} \\"
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
        echo "## ${scr_sub} ##"
        echo ""
        cat "${scr_sub}"
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
            ${scr_sub} \
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
    #  If --slurm was not specified, run jobs in parallel with GNU Parallel
    if [[ ${par_job} -gt 1 ]]; then
        config="${err_out}/${nam_job}.config_parallel.txt"

        if [[ -f "${config}" ]]; then rm "${config}"; fi
        touch "${config}" || {
            echo_error "Failed to create a GNU Parallel configuration file."
            exit_1
        }

        for idx in "${!arr_infiles[@]}"; do
            echo \
                "${env_nam}" \
                "${threads}" \
                "${arr_infiles[idx]}" \
                "${dir_out}" \
                "${sfx_se}" \
                "${sfx_pe}" \
                "${err_out}" \
                "${nam_job}" \
                    >> "${config}"
        done

        if ${dry_run} || ${verbose}; then
            echo "#######################################"
            echo "## Parallel call(s) to trim fastq(s) ##"
            echo "#######################################"
            echo ""

            parallel --colsep ' ' --jobs "${par_job}" --dryrun \
                "bash \"${scr_sub}\" {1} {2} {3} {4} {5} {6} {7} {8}" \
                :::: "${config}"

            echo ""
            echo ""
        fi

        if ! ${dry_run}; then
            parallel --colsep ' ' --jobs "${par_job}" \
                "bash \"${scr_sub}\" {1} {2} {3} {4} {5} {6} {7} {8}" \
                :::: "${config}"
        fi
    else
        #  If --slurm was not specified and 'par_job=1', then run jobs in
        #+ serial
        if ${dry_run} || ${verbose}; then
            echo "##################################################"
            echo "## Serialized job execution via 'submit' script ##"
            echo "##################################################"
            echo ""
            echo "bash \"${scr_sub}\" \\"
            echo "    \"${env_nam}\" \\"
            echo "    \"${threads}\" \\"
            echo "    \"${infiles}\" \\"
            echo "    \"${dir_out}\" \\"
            echo "    \"${sfx_se}\" \\"
            echo "    \"${sfx_pe}\" \\"
            echo "    \"${err_out}\" \\"
            echo "    \"${nam_job}\""
            echo ""
            echo ""
        fi
        
        if ! ${dry_run}; then
            bash "${scr_sub}" \
                "${env_nam}" \
                "${threads}" \
                "${infiles}" \
                "${dir_out}" \
                "${sfx_se}" \
                "${sfx_pe}" \
                "${err_out}" \
                "${nam_job}"
        fi
    fi
fi
