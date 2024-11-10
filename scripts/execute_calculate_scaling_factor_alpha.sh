#!/bin/bash

#  execute_calculate_scaling_factor_alpha.sh
#  KA


#  Run script in interactive/test mode (true) or command-line mode (false)
interactive=false

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive}; then set -euo pipefail; fi

#  Set the path to the "scripts" directory
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
}


#  Set up paths, values, and parameters for interactive mode
function set_interactive() {
    #  Set hardcoded paths, values, etc.
    ## WARNING: Change the values if you're not Kris and `interactive=true` ##
    dir_rep="${HOME}/tsukiyamalab/Kris/202X_protocol_ChIP"
    dir_bam="${dir_rep}/03_bam"
    dir_scf="${dir_rep}/06_sf"
    {
        aligner="bowtie2" 
        a_type="global"
        flg=2  # "NA"
        mapq=1
        details="${aligner}/${a_type}/flag-${flg}_mapq-${mapq}"
    }
    dir_ip="${dir_bam}/${details}/sc"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=true
    threads=8
    infiles=$(
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_ip}" \
            --pattern "*.bam" \
            --depth 1 \
            --include "IP*,*Q_*,*Esa1*"
    )
    table="${dir_rep}/data/measurements_siq_chip.tsv"
    {
        dir_out="${dir_scf}/${details}/alpha"
    }
    outfile="${dir_out}/test.txt"  # "/dev/stdout"
    flg_dep=true
    flg_len=true
    flg_in=true
    flg_mc=true
    err_out="${dir_out}/err_out"
    nam_job="calc_sf_alpha"
    slurm=true
    serial=false
    max_job=6
    time="0:30:00"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_analyze"
scr_sub="${dir_scr}/submit_calculate_scaling_factor_alpha.sh"
scr_par="${dir_scr}/parse_metadata_siq_chip.py"
scr_alf="${dir_scr}/calculate_scaling_factor_alpha.py"

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=1
infiles=""
table=""
outfile=""
flg_dep=false
flg_len=false
flg_in=false
flg_mc=false
err_out=""
nam_job="calc_sf_alpha"
serial=false
slurm=false
max_job=6
time="0:30:00"

#  Assign variable for help message
show_help=$(
cat << EOM
execute_calculate_scaling_factor_alpha.sh
  [--verbose] [--dry_run] --threads <int> --infiles <str> --table <str>
  --outfile <str> [--flg_dep] [--flg_len] [--flg_in] [--flg_mc] --err_out <str>
  --nam_job <str> [--serial] [--slurm] [--max_job <int>] [--time <str>]

Description:
  execute_calculate_scaling_factor_alpha.sh performs... #TODO

Arguments:
   -h, --help     Display this help message and exit (0).
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Perform a dry run without executing commands (optional).
   -t, --threads  Number of threads to use (required; default: ${threads}).
   -i, --infiles  Comma-separated serialized list of IP coordinate-sorted BAM
                  infiles, including paths (required).
  -tb, --table    #TODO
   -o, --outfile  #TODO
  -fd, --flg_dep  Use Samtools to calculate the number of alignments.
  -fl, --flg_len  Use Samtools and awk to calculate the mean fragment length.
  -fi, --flg_in   Include sample input alpha values in outfile.
  -fm, --flg_mc   Include additional measurements and calculations in outfile.
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (required; default: \${dir_out}/err_out).
  -nj, --nam_job  The name of the job (required; default: ${nam_job}).
  -se, --serial   Run jobs in serial on the current system.
  -sl, --slurm    Submit jobs to the SLURM scheduler.
  -mj, --max_job  The maximum number of jobs to run at one time (required if
                  --slurm is specified, ignored if not; default: ${max_job}).
  -tm, --time     The length of time, in 'h:mm:ss' format, for the SLURM job
                  (required if --slurm is specified, ignored if not; default:
                  ${time}).

#TODO Notes, examples, etc.
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
            -tb|--table)   table="${2}";   shift 2 ;;
             -o|--outfile) outfile="${2}"; shift 2 ;;
            -fd|--flg_dep) flg_dep=true;   shift 1 ;;
            -fl|--flg_len) flg_len=true;   shift 1 ;;
            -fi|--flg_in)  flg_in=true;    shift 1 ;;
            -fm|--flg_mc)  flg_mc=true;    shift 1 ;;
            -eo|--err_out) err_out="${2}"; shift 2 ;;
            -nj|--nam_job) nam_job="${2}"; shift 2 ;;
            -sl|--slurm)   slurm=true;     shift 1 ;;
            -se|--serial)  serial=true;    shift 1 ;;
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
check_supplied_arg -a "${env_nam}" -n "env_nam"

check_supplied_arg -a "${scr_sub}" -n "scr_sub"
check_exists_file_dir "f" "${scr_sub}" "scr_sub"

check_supplied_arg -a "${scr_par}" -n "scr_par"
check_exists_file_dir "f" "${scr_par}" "scr_par"

check_supplied_arg -a "${scr_alf}" -n "scr_alf"
check_exists_file_dir "f" "${scr_alf}" "scr_alf"

check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

check_supplied_arg -a "${infiles}" -n "infiles"
check_exists_file_dir "d" "$(dirname "${infiles%%[,;]*}")" "infiles"

check_supplied_arg -a "${table}" -n "table"
check_exists_file_dir "f" "${table}" "table"

check_supplied_arg -a "${outfile}" -n "outfile"
check_exists_file_dir "d" "$(dirname "${outfile}")" "outfile"

if [[ -n "${err_out}" ]]; then
    check_exists_file_dir "d" "${err_out}" "err_out"
elif [[ -z "${err_out}" ]]; then
    err_out="${dir_out}/err_out"
    check_exists_file_dir "d" "${err_out}" "err_out"
fi

if [[ -z ${serial} ]] || [[ -z ${slurm} ]]; then
    echo_error "Either --serial or --slurm must be specified."
    exit_1
elif ${serial} && ${slurm}; then
    echo_error "Only one of --serial or --slurm should be provided"
    exit_1
fi

if ${slurm}; then
    check_supplied_arg -a "${nam_job}" -n "nam_job"
    
    check_supplied_arg -a "${max_job}" -n "max_job"
    check_int_pos "${max_job}" "max_job"
    
    check_supplied_arg -a "${time}" -n "time"
    check_format_time "${time}"
fi

#  Activate environment and check that dependencies are in PATH
handle_env "${env_nam}"

check_program_path awk
check_program_path python
check_program_path samtools

#  Parse the --infiles argument into an array, then validate the infile value
#+ assignments
IFS=',' read -r -a arr_infiles <<< "${infiles}"  # unset arr_infiles
# for infile in "${arr_infiles[@]}"; do echo "${infile}"; done  # unset infile

#  Check that each infile exists; if not, exit
for infile in "${arr_infiles[@]}"; do
    check_exists_file_dir "f" "${infile}"
done

if ${slurm}; then
    if [[ "${max_job}" -gt "${#arr_infiles[@]}" ]]; then
        echo_warning \
            "The maximum number of SLURM jobs to run at a time, ${max_job}," \
            "is greater than the number of infiles, ${#arr_infiles[@]}." \
            "Adjusting max_job to ${#arr_infiles[@]}."
        max_job="${#arr_infiles[@]}"
    fi
fi


#  Do the main work ===========================================================
#  Report argument variable assignments if in "verbose mode"
if ${verbose}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "env_nam=${env_nam}"
    echo "scr_sub=${scr_sub}"
    echo "scr_par=${scr_par}"
    echo "scr_alf=${scr_alf}"
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
    echo "table=${table}"
    echo "outfile=${outfile}"
    echo "flg_dep=${flg_dep}"
    echo "flg_len=${flg_len}"
    echo "flg_in=${flg_in}"
    echo "flg_mc=${flg_mc}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo "slurm=${slurm}"
    echo "serial=${serial}"
    echo "max_job=${max_job}"
    echo "time=${time}"
    echo ""
    echo ""
fi

# shellcheck disable=SC1083,SC2157,SC2046,SC2086
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
        echo "    --export=dir_scr=${dir_scr} \\"
        echo "    ${scr_sub} \\"
        echo "        --threads ${threads} \\"
        echo "        --infiles ${infiles} \\"
        echo "        --table ${table} \\"
        echo "        --outfile ${outfile} \\"
        echo "        $(if ${flg_dep}; then echo "--flg_dep"; fi) \\"
        echo "        $(if ${flg_len}; then echo "--flg_len"; fi) \\"
        echo "        $(if ${flg_in}; then echo "--flg_in"; fi) \\"
        echo "        $(if ${flg_mc}; then echo "--flg_mc"; fi) \\"
        echo "        --err_out ${err_out} \\"
        echo "        --nam_job ${nam_job} \\"
        echo "        --env_nam ${env_nam} \\"
        echo "        --scr_par ${scr_par} \\"
        echo "        --scr_alf ${scr_alf}"
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
        sbatch \
            --job-name=${nam_job} \
            --nodes=1 \
            --cpus-per-task=${threads} \
            --time=${time} \
            --error=${err_out}/${nam_job}.%A-%a.stderr.txt \
            --output=${err_out}/${nam_job}.%A-%a.stdout.txt \
            --array=1-${#arr_infiles[@]}%${max_job} \
            --export=dir_scr="${dir_scr}" \
            ${scr_sub} \
                --threads ${threads} \
                --infiles ${infiles} \
                --table ${table} \
                --outfile ${outfile} \
                $(if ${flg_dep}; then echo "--flg_dep"; fi) \
                $(if ${flg_len}; then echo "--flg_len"; fi) \
                $(if ${flg_in}; then echo "--flg_in"; fi) \
                $(if ${flg_mc}; then echo "--flg_mc"; fi) \
                --err_out ${err_out} \
                --nam_job ${nam_job} \
                --env_nam ${env_nam} \
                --scr_par ${scr_par} \
                --scr_alf ${scr_alf}
    fi
fi

#TODO Write code for local, serial job submission
