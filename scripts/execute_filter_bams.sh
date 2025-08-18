#!/bin/bash

#  execute_filter_bams.sh
#  KA


#  Run script in interactive mode (true) or command-line mode (false)
interactive=false

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive:-false}; then set -euo pipefail; fi

#  Set the path to the "scripts" directory
if ${interactive:-false}; then
    ## WARNING: If 'interactive=true', change path as needed ##
    dir_scr="${HOME}/repos/protocol_chipseq_signal_norm/scripts"
else
    dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi


#  Source and define functions ================================================
#  Set the path to the "functions" directory
dir_fnc="${dir_scr}/functions"

# shellcheck disable=SC1090
for fnc in \
    check_exists_file_dir \
    check_format_time \
    check_int_pos \
    check_program_path \
    check_str_delim \
    check_supplied_arg \
    echo_error \
    echo_warning \
    exit_0 \
    exit_1 \
    filter_bam_sc \
    filter_bam_sp \
    handle_env \
    print_parallel_info \
    set_params_parallel
do
    source "${dir_fnc}/${fnc}.sh"
done
unset fnc


#  Set paths, values, arguments, etc. for interactive mode
function set_interactive() {
    #  Set base repository paths
    dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
    dir_rep="${dir_bas}/protocol_chipseq_signal_norm"

    #  Define data directories and related variables
    dir_dat="${dir_rep}/data"
    dir_pro="${dir_dat}/processed"
    dir_aln="${dir_pro}/align_fastqs"
    aligner="bowtie2"
    a_type="global"
    req_flg=true
    flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)" 
    mapq=1
    str_det="${aligner}_${a_type}_flag-${flg}_mapq-${mapq}"
    dir_det="${dir_aln}/${str_det}"
    dir_ini="${dir_det}/init"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=true
    retain="sc"
    threads=8
    infiles="$(  ## WARNING: Change search parameters as needed ##
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_ini}" \
            --pattern "*.bam" \
            --depth 1
    )"
    dir_out="${dir_det}/${retain}"
    mito=false
    tg=false
    mtr=false
    chk_chr=false
    err_out="${dir_out}/logs"
    nam_job="filter_bams"
    max_job=2
    slurm=false
    time="1:00:00"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_protocol"
scr_sub="${dir_scr}/submit_filter_bams.sh"  ## 'scr_fnc' is defined below ##
par_job=""

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=1
infiles=""
dir_out=""
retain="sc"
mito=false
tg=false
mtr=false
chk_chr=false
err_out=""
nam_job="filter_bams"
max_job=6
slurm=false
time="0:30:00"

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  execute_filter_bams.sh
    [--verbose] [--dry_run] --threads <int> --infiles <str> --dir_out <str>
    --retain <str> [--mito] [--tg] [--mtr] [--chk_chr] --err_out <str>
    --nam_job <str> --max_job <int> [--slurm] [--time <str>]

Description:
  The driver script 'execute_filter_bams.sh' filters BAM infiles to retain
  species-specific chromosomes for S. cerevisiae ("main" alignments) or S.
  pombe ("spike-in" alignments).

  Optional features include retaining mitochondrial (S. cerevisiae and S.
  pombe: '--mito') and additional chromosomes (S. pombe: '--tg', '--mtr'), and
  performing checks on chromosomes in filtered BAM outfiles ('--chk_chr').

  The script supports parallel execution via SLURM or GNU Parallel, or can run
  serially.

Arguments:
   -h, --help     Print this help message and exit.
   -v, --verbose  Run script in "verbose" mode (optional).
  -dr, --dry_run  Run the command in "check" mode (optional).
   -t, --threads  Number of threads to use (default: '${threads}').
   -i, --infiles  Comma-separated serialized string of coordinate-sorted BAM
                  input files.
  -do, --dir_out  The directory to store species-filtered and -reheadered BAM
                  output files.
   -r, --retain   Specify species chromosomes to retain: S. cerevisiae, "sc";
                  S. pombe, "sp" (default: '${retain}').
   -m, --mito     Retain mitochondrial chromosome (optional).
  -tg, --tg       Retain SP_II_TG chromosome (optional, sp only).
  -mr, --mtr      Retain SP_MTR chromosome (optional, sp only).
  -cc, --chk_chr  Check chromosomes in filtered BAM outfile (optional)
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (default: '\${dir_out}/err_out').
  -nj, --nam_job  The name of the job, which is used when writing stderr and
                  stdout TXT files (default: '${nam_job}').
  -mj, --max_job  Maximum number of jobs to run concurrently (default: '${max_job}').
                    - If '--slurm' is specified, controls SLURM array tasks.
                    - If '--slurm' is not specified:
                      + If 'max_job' is greater than 1, jobs run in parallel
                        via GNU Parallel.
                      + If 'max_job' is 1, jobs run sequentially (serial mode).
  -sl, --slurm    Submit jobs to the SLURM scheduler; otherwise, run them in
                  serial (optional).
  -tm, --time     The length of time, in 'h:mm:ss' format, for the SLURM job
                  (required if '--slurm' is specified, ignored if not; default:
                  '${time}').

Dependencies:
  - AWK
  - Bash or Zsh
  - GNU Parallel (when '--slurm' is not specified but multiple jobs are)
  - grep
  - mv
  - rm
  - Samtools
  - SLURM (when '--slurm' is specified)

Notes:
  - When the '--slurm' flag is used, jobs are parallelized via SLURM array
    tasks; otherwise... #TODO
  - BAM infiles must be coordinate-sorted.
  - Flag '--mito' applies to either S. cerevisiae or S. pombe data.
  - Flags '--tg' and '--mtr' apply only to S. pombe data.

Examples:
  \`\`\`
  #  Use SLURM to filter BAM files for S. cerevisiae ("sc") chromosomes (i.e.,
  #+ "main" alignments)
  retain="sc"
  bash "\${dir_scr}/execute_filter_bams.sh"
      --verbose
      --threads "\${threads}"
      --infiles "\${infiles}"
      --dir_out "\${dir_out}/\${retain}"
      --err_out "\${dir_out}/\${retain}/logs"
      --retain  "\${retain}"
      --slurm

  #  Use GNU Parallel to filter BAM files for S. pombe ("sp") chromosomes
  #+ (i.e., "spike-in" alignments)
  retain="sp"
  bash "\${dir_scr}/execute_filter_bams.sh"
      --verbose
      --threads "\${threads}"
      --infiles "\${infiles}"
      --dir_out "\${dir_out}/\${retain}"
      --err_out "\${dir_out}/\${retain}/logs"
      --retain  "\${retain}"
  \`\`\`
EOM
)

#  Parse arguments
if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    exit_0
fi

if ${interactive:-false}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -v|--verbose) verbose=true;   shift 1 ;;
            -dr|--dry_run) dry_run=true;   shift 1 ;;
             -t|--threads) threads="${2}"; shift 2 ;;
             -i|--infiles) infiles="${2}"; shift 2 ;;
            -do|--dir_out) dir_out="${2}"; shift 2 ;;
             -r|--retain)
                retain="$(echo "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;
             -m|--mito)    mito=true;      shift 1 ;;
            -tg|--tg)      tg=true;        shift 1 ;;
            -mr|--mtr)     mtr=true;       shift 1 ;;
            -cc|--chk_chr) chk_chr=true;   shift 1 ;;
            -eo|--err_out) err_out="${2}"; shift 2 ;;
            -nj|--nam_job) nam_job="${2}"; shift 2 ;;
            -mj|--max_job) max_job="${2}"; shift 2 ;;
            -sl|--slurm)   slurm=true;     shift 1 ;;
            -tm|--time)    time="${2}";    shift 2 ;;
            *)
                echo "## Unknown parameter passed: '${1}' ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                exit_1
                ;;
        esac
    done
fi

#  Check arguments
check_supplied_arg -a "${env_nam}" -n "env_nam"

check_supplied_arg -a "${scr_sub}" -n "scr_sub"
check_exists_file_dir "f" "${scr_sub}" "scr_sub"

check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

check_supplied_arg -a "${infiles}" -n "infiles"
check_exists_file_dir "d" "$(dirname "${infiles%%,*}")" "infiles"
check_str_delim "infiles" "${infiles}"

check_supplied_arg -a "${dir_out}" -n "dir_out"
check_exists_file_dir "d" "${dir_out}" "dir_out"

case "${retain}" in
    sc|sp) 
        #  Based on '--retain' assignment, define and validate function script
        scr_fnc="${dir_fnc}/filter_bam_${retain}.sh"
        
        check_supplied_arg -a "${scr_fnc}" -n "scr_fnc"
        check_exists_file_dir "f" "${scr_fnc}" "scr_fnc"
        ;;
    *)
        echo_error \
            "Selection associated with '--retain' is not valid: '${retain}'." \
            "Selection must be 'sc' or 'sp'."
        exit_1
        ;;
esac

if [[ -z "${err_out}" ]]; then
    err_out="${dir_out}/err_out"
fi
check_exists_file_dir "d" "${err_out}" "err_out"

check_supplied_arg -a "${nam_job}" -n "nam_job"


#  Parse and validate infiles -------------------------------------------------
IFS=',' read -r -a arr_infile <<< "${infiles}" && unset IFS

#  Check that each infile exists; if not, exit
for infile in "${arr_infile[@]}"; do
    check_exists_file_dir "f" "${infile}" "infile"
done


#  Parse job execution parameters ---------------------------------------------
if ${slurm:-false}; then
    check_supplied_arg -a "${max_job}" -n "max_job"
    check_int_pos "${max_job}" "max_job"

    max_job=$(reset_max_job "${max_job}" "${#arr_infile[@]}")
    
    check_supplied_arg -a "${time}" -n "time"
    check_format_time "${time}"
else
    IFS=';' read -r threads par_job < <(
        set_params_parallel "${threads}" "${max_job}" "${par_job}"
    )
    unset IFS max_job time

    check_supplied_arg -a "${par_job}" -n "par_job"
    check_int_pos "${par_job}" "par_job"
fi

#  Debug parallelization information
if ${verbose}; then
    print_parallel_info \
        "${slurm}" "${max_job:-UNSET}" "${par_job}" "${threads}" "arr_infile"
fi


#  Activate environment and check that dependencies are in PATH ---------------
handle_env "${env_nam}" > /dev/null

check_program_path awk
check_program_path grep
check_program_path mv
check_program_path rm
check_program_path samtools

if ${slurm:-false}; then
    check_program_path sbatch
elif [[ ${par_job} -gt 1 ]]; then
    check_program_path parallel
fi


#  Do the main work ===========================================================
if ${verbose}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "env_nam=${env_nam}"
    echo "scr_sub=${scr_sub}"
    echo "scr_fnc=${scr_fnc}"
    echo "par_job=${par_job:-UNSET}"
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
    echo "retain=${retain}"
    echo "mito=${mito}"
    echo "tg=${tg}"
    echo "mtr=${mtr}"
    echo "chk_chr=${chk_chr}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo "max_job=${max_job:-UNSET}"
    echo "slurm=${slurm}"
    echo "time=${time:-UNSET}"
    echo ""
    echo ""
    echo "#################################"
    echo "## Array derived from variable ##"
    echo "#################################"
    echo ""
    echo "arr_infile=( ${arr_infile[*]} )"
    echo ""
    echo ""
fi

# shellcheck disable=SC1083,SC2157
if ${slurm:-false}; then
    #  If --slurm was specified, run jobs in parallel via SLURM
    if ${dry_run} || ${verbose}; then
        echo "######################"
        echo "## Call to 'sbatch' ##"
        echo "######################"
        echo ""
        echo "sbatch \\"
        echo "    --job-name=${nam_job} \\"
        echo "    --nodes=1 \\"
        echo "    --cpus-per-task=${threads} \\"
        echo "    --time=${time} \\"
        echo "    --error=${err_out}/${nam_job}.%A-%a.stderr.txt \\"
        echo "    --output=${err_out}/${nam_job}.%A-%a.stdout.txt \\"
        echo "    --array=1-${#arr_infile[@]}%${max_job} \\"
        echo "    ${scr_sub} \\"
        echo "        --scr_fnc ${scr_fnc} \\"
        echo "        --threads ${threads} \\"
        echo "        --infiles ${infiles} \\"
        echo "        --dir_out ${dir_out} \\"
        echo "        $(if ${mito}; then echo "--mito"; fi) \\"
        echo "        $(if ${tg}; then echo "--tg"; fi) \\"
        echo "        $(if ${mtr}; then echo "--mtr"; fi) \\"
        echo "        $(if ${chk_chr}; then echo "--chk_chr"; fi) \\"
        echo "        --err_out ${err_out} \\"
        echo "        --nam_job ${nam_job}"
        echo ""
        echo ""
    fi

    if ! ${dry_run}; then
        # shellcheck disable=SC2046,SC2086
        sbatch \
            --job-name=${nam_job} \
            --nodes=1 \
            --cpus-per-task=${threads} \
            --time=${time} \
            --error=${err_out}/${nam_job}.%A-%a.stderr.txt \
            --output=${err_out}/${nam_job}.%A-%a.stdout.txt \
            --array=1-${#arr_infile[@]}%${max_job} \
            ${scr_sub} \
                -en ${env_nam} \
                -sf ${scr_fnc} \
                 -t ${threads} \
                 -i ${infiles} \
                -do ${dir_out} \
                $(if ${mito}; then echo "-m"; fi) \
                $(if ${tg}; then echo "-tg"; fi) \
                $(if ${mtr}; then echo "-mr"; fi) \
                $(if ${chk_chr}; then echo "-cc"; fi) \
                -eo ${err_out} \
                -nj ${nam_job}
    fi
else
    #  If '--slurm' was not specified, run jobs in parallel with GNU Parallel
    if [[ "${par_job}" -gt 1 ]]; then
        config="${err_out}/${nam_job}.config_parallel.txt"

        if [[ -f "${config}" ]]; then rm "${config}"; fi
        touch "${config}" || {
            echo_error "Failed to create a GNU Parallel configuration file."
            exit_1
        }

        #  Populate GNU Parallel configuration file
        for idx in "${!arr_infile[@]}"; do
            echo \
                "${env_nam}" \
                "${scr_fnc}" \
                "${threads}" \
                "${arr_infile[idx]}" \
                "${dir_out}" \
                "${err_out}" \
                "${nam_job}" \
                    >> "${config}"
        done

        #  Construct command to be passed to GNU Parallel
        cmd="bash ${scr_sub}"
        cmd+=" -en {1}"  # Environment name
        cmd+=" -sf {2}"  # Function script
        cmd+="  -t {3}"  # Number of threads
        cmd+="  -i {4}"  # Input BAM file
        cmd+=" -do {5}"  # Output BAM directory
        
        if ${mito};    then cmd+=" -m";  fi  # Flag: Retain mito. alignments
        if ${tg};      then cmd+=" -tg"; fi  # Flag: Retain SP_II_TG alignments
        if ${mtr};     then cmd+=" -mr"; fi  # Flag: Retain SP_MTR alignments
        if ${chk_chr}; then cmd+=" -cc"; fi  # Flag: Check chr. in outfile
        
        cmd+=" -eo {6}"  # Directory for stderr/stdout logs
        cmd+=" -nj {7}"  # Job name
        cmd+="  > {6}/{7}_par.{4/.}.stdout.txt"  # 'scr_sub' stdout log
        cmd+=" 2> {6}/{7}_par.{4/.}.stderr.txt"  # 'scr_sub' stderr log

        if ${dry_run} || ${verbose}; then
            echo "####################################################"
            echo "## Parallelized job execution via 'submit' script ##"
            echo "####################################################"
            echo ""

            parallel --colsep ' ' --jobs "${par_job}" --dryrun \
                "${cmd}" \
                :::: "${config}"

            echo ""
            echo ""
        fi

        if ! ${dry_run}; then
            parallel --colsep ' ' --jobs "${par_job}" \
                "${cmd}" \
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
            echo "    -en ${env_nam} \\"
            echo "    -sf ${scr_fnc} \\"
            echo "     -t ${threads} \\"
            echo "     -i ${infiles} \\"
            echo "    -do ${dir_out} \\"
            if ${mito};    then echo "     -m \\"; fi
            if ${tg};      then echo "    -tg \\"; fi
            if ${mtr};     then echo "    -mr \\"; fi
            if ${chk_chr}; then echo "    -cc \\"; fi
            echo "    -eo ${err_out} \\"
            echo "    -nj ${nam_job} \\"
            echo "         > ${err_out}/${nam_job}_ser.stdout.txt \\"
            echo "        2> ${err_out}/${nam_job}_ser.stderr.txt"
            echo ""
            echo ""
        fi

        # shellcheck disable=SC2046
        bash "${scr_sub}" \
            -en "${env_nam}" \
            -sf "${scr_fnc}" \
             -t "${threads}" \
             -i "${infiles}" \
            -do "${dir_out}" \
            $(if ${mito};    then echo  "-m"; fi) \
            $(if ${tg};      then echo "-tg"; fi) \
            $(if ${mtr};     then echo "-mr"; fi) \
            $(if ${chk_chr}; then echo "-cc"; fi) \
            -eo "${err_out}" \
            -nj "${nam_job}" \
                 > "${err_out}/${nam_job}_ser.stdout.txt" \
                2> "${err_out}/${nam_job}_ser.stderr.txt"
    fi
fi
