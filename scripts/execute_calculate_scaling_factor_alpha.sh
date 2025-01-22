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
    source "${dir_fnc}/echo_error.sh"
    source "${dir_fnc}/echo_warning.sh"
    source "${dir_fnc}/exit_0.sh"
    source "${dir_fnc}/exit_1.sh"
    source "${dir_fnc}/handle_env.sh"
    source "${dir_fnc}/reset_max_job.sh"
}


#  Set up paths, values, and parameters for interactive mode
function set_interactive() {
    #  Set hardcoded paths, values, etc.
    ## WARNING: If interactive=true, change values as needed ##
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
    det_cvg="${aligner}_${a_type}_${det_bam}"
    typ_cvg="alpha"

    dir_aln="${dir_pro}/align_${aligner}_${a_type}"
    dir_bam="${dir_aln}/${det_bam}/sc"
    dir_cvg="${dir_pro}/compute_coverage"
    dir_out="${dir_cvg}/${det_cvg}/${typ_cvg}/tables"

    pattern="*.bam"
    exclude="*Brn1*"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=false
    threads=8
    ser_ip="$(  ## WARNING: Change search parameters as needed ##
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_bam}" \
            --pattern "${pattern}" \
            --include "IP*" \
            --exclude "${exclude}"
    )"
    ser_in="$(  ## WARNING: Change search parameters as needed ##
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_bam}" \
            --pattern "${pattern}" \
            --include "in*" \
            --exclude "${exclude}"
    )"
    table="${dir_dat}/raw/docs/measurements_siqchip.tsv"
    eqn="6nd"  # "5"
    outfile="${dir_out}/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_${typ_cvg}_${eqn}.tsv"
    flg_dep=true
    flg_len=true
    flg_mc=true
    err_out="${dir_out}/logs"
    nam_job="calc_sf_alpha_${eqn}"
    slurm=true
    max_job=6
    time="0:30:00"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_analyze"
scr_sub="${dir_scr}/submit_calculate_scaling_factor_alpha.sh"
scr_met="${dir_scr}/parse_metadata_siq_chip.py"
scr_alf="${dir_scr}/calculate_scaling_factor_alpha.py"
denom=4

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=1
ser_ip=""
ser_in=""
table=""
eqn="6nd"
outfile=""
flg_dep=false
flg_len=false
flg_mc=false
err_out=""
nam_job="calc_sf_alpha_${eqn}"
slurm=false
max_job=6
time="0:30:00"

#TODO: Add '--rnd' argument

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  execute_calculate_scaling_factor_alpha.sh
    [--verbose] [--dry_run] --threads <int> --ser_ip <str> --ser_in <str>
    --table <str> --eqn <str> --outfile <str> [--flg_dep] [--flg_len]
    [--flg_mc] --err_out <str> --nam_job <str> [--slurm] [--max_job <int>]
    [--time <str>]

Description:
  execute_calculate_scaling_factor_alpha.sh performs... #TODO

Arguments:
   -h, --help     Display this help message and exit (0).
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Perform a dry run without executing commands (optional).
   -t, --threads  Number of threads to use (required; default: ${threads}).
  -sp, --ser_ip   Comma-separated serialized list of IP coordinate-sorted BAM
                  infiles, including paths (required).
  -sn, --ser_in   Comma-separated serialized list of input coordinate-sorted
                  BAM infiles, including paths (required).
  -tb, --table    #TODO
  -eq, --eqn      #TODO
   -o, --outfile  #TODO
  -fd, --flg_dep  Use Samtools to calculate the number of alignments.
  -fl, --flg_len  Use Samtools and awk to calculate the mean fragment length.
  -fm, --flg_mc   Include additional measurements and calculations in outfile.
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (required; default: \${dir_out}/err_out).
  -nj, --nam_job  The name of the job (required; default: ${nam_job}).
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
            -sp|--ser_ip)  ser_ip="${2}";  shift 2 ;;
            -sn|--ser_in)  ser_in="${2}";  shift 2 ;;
            -tb|--table)   table="${2}";   shift 2 ;;
            -eq|--eqn)     eqn="${2}";     shift 2 ;;
             -o|--outfile) outfile="${2}"; shift 2 ;;
            -fd|--flg_dep) flg_dep=true;   shift 1 ;;
            -fl|--flg_len) flg_len=true;   shift 1 ;;
            -fm|--flg_mc)  flg_mc=true;    shift 1 ;;
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
check_supplied_arg -a "${env_nam}" -n "env_nam"

check_supplied_arg -a "${scr_sub}" -n "scr_sub"
check_exists_file_dir "f" "${scr_sub}" "scr_sub"

check_supplied_arg -a "${scr_met}" -n "scr_met"
check_exists_file_dir "f" "${scr_met}" "scr_met"

check_supplied_arg -a "${scr_alf}" -n "scr_alf"
check_exists_file_dir "f" "${scr_alf}" "scr_alf"

check_supplied_arg -a "${denom}" -n "denom"
check_int_pos "${denom}" "denom"

check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

check_supplied_arg -a "${ser_ip}" -n "ser_ip"
check_exists_file_dir "d" "$(dirname "${ser_ip%%[,;]*}")" "ser_ip"

check_supplied_arg -a "${ser_in}" -n "ser_in"
check_exists_file_dir "d" "$(dirname "${ser_in%%[,;]*}")" "ser_in"

check_supplied_arg -a "${table}" -n "table"
check_exists_file_dir "f" "${table}" "table"

check_supplied_arg -a "${eqn}" -n "eqn"
case "${eqn}" in
    5|5nd|6|6nd) : ;;
    *)
        echo_error \
            "Equation (--eqn) was assigned '${eqn}' but must be '5', '5nd'," \
            "'6', or '6nd'."
        exit_1
    ;;
esac

check_supplied_arg -a "${outfile}" -n "outfile"
check_exists_file_dir "d" "$(dirname "${outfile}")" "outfile"

# check_supplied_arg -a "${rnd}" -n "rnd"
# check_int_pos "${rnd}" "rnd"

if [[ -n "${err_out}" ]]; then
    check_exists_file_dir "d" "${err_out}" "err_out"
elif [[ -z "${err_out}" ]]; then
    err_out="${dir_out}/err_out"
    check_exists_file_dir "d" "${err_out}" "err_out"
fi

check_supplied_arg -a "${nam_job}" -n "nam_job"

if ${slurm}; then
    check_supplied_arg -a "${max_job}" -n "max_job"
    check_int_pos "${max_job}" "max_job"
    
    check_supplied_arg -a "${time}" -n "time"
    check_format_time "${time}"
else
    par_job=$(( threads / denom ))
    threads=${denom}

    check_supplied_arg -a "${par_job}" -n "par_job"
    check_int_pos "${par_job}" "par_job"

    check_supplied_arg -a "${threads}" -n "threads"
    check_int_pos "${threads}" "threads"
fi

#  Activate environment and check that dependencies are in PATH
handle_env "${env_nam}" > /dev/null

check_program_path awk
if ! ${slurm} && [[ ${threads} -gt 1 ]]; then check_program_path parallel; fi
check_program_path python
check_program_path samtools
if ${slurm}; then check_program_path sbatch; fi

#  Parse the '--ser_ip' and '--ser_in' arguments into separate arrays, then
#+ validate the infile value assignments
IFS=',' read -r -a arr_ip <<< "${ser_ip}"
IFS=',' read -r -a arr_in <<< "${ser_in}"

for file in "${arr_ip[@]}" "${arr_in[@]}"; do
    check_exists_file_dir "f" "${file}" "file"
done
unset file

#  Reset 'max_job' if it is greater than the number of infiles
if ${slurm}; then
    max_job=$(reset_max_job "${max_job}" "${#arr_ip[@]}")
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
    echo "scr_met=${scr_met}"
    echo "scr_alf=${scr_alf}"
    echo "denom=${denom}"
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
    echo "ser_ip=${ser_ip}"
    echo "ser_in=${ser_in}"
    echo "table=${table}"
    echo "eqn=${eqn}"
    echo "outfile=${outfile}"
    echo "flg_dep=${flg_dep}"
    echo "flg_len=${flg_len}"
    echo "flg_mc=${flg_mc}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo "slurm=${slurm}"
    echo "max_job=${max_job}"
    echo "time=${time}"
    echo ""
    echo ""
    echo "###################################"
    echo "## Arrays derived from variables ##"
    echo "###################################"
    echo ""
    echo "arr_ip=( ${arr_ip[*]} )"
    echo "arr_in=( ${arr_in[*]} )"
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
        echo "    --array=1-${#arr_ip[@]}%${max_job} \\"
        echo "    ${scr_sub} \\"
        echo "        --threads ${threads} \\"
        echo "        --ser_ip ${ser_ip} \\"
        echo "        --ser_in ${ser_in} \\"
        echo "        --table ${table} \\"
        echo "        --eqn ${eqn} \\"
        echo "        --outfile ${outfile} \\"
        echo "        $(if ${flg_dep}; then echo "--flg_dep"; fi) \\"
        echo "        $(if ${flg_len}; then echo "--flg_len"; fi) \\"
        echo "        $(if ${flg_mc}; then echo "--flg_mc"; fi) \\"
        echo "        --err_out ${err_out} \\"
        echo "        --nam_job ${nam_job} \\"
        echo "        --env_nam ${env_nam} \\"
        echo "        --scr_met ${scr_met} \\"
        echo "        --scr_alf ${scr_alf}"
        echo ""
        echo ""
        # echo "#########################################"
        # echo "## Contents of SLURM submission script ##"
        # echo "#########################################"
        # echo ""
        # echo "## ${scr_sub} ##"
        # echo ""
        # cat "${scr_sub}"
        # echo ""
        # echo ""
    fi

    if ! ${dry_run}; then
        #  To prevent potential race conditions from concurrent writes,
        #+ pre-write the header to the outfile before running SLURM jobs
        #+ (for more details, see related comment for GNU Parallel jobs below)
        if ${flg_mc}; then
            echo -e \
                "samp_ip\tsamp_in\talpha\teqn\tmass_ip\tmass_in\tvol_all\tvol_in\tdep_ip\tdep_in\tlen_ip\tlen_in" \
                    >> "${outfile}"
        else
            echo -e "samp_ip\tsamp_in\talpha\teqn" >> "${outfile}"
        fi

        sbatch \
            --job-name=${nam_job} \
            --nodes=1 \
            --cpus-per-task=${threads} \
            --time=${time} \
            --error=${err_out}/${nam_job}.%A-%a.stderr.txt \
            --output=${err_out}/${nam_job}.%A-%a.stdout.txt \
            --array=1-${#arr_ip[@]}%${max_job} \
            ${scr_sub} \
                --threads ${threads} \
                --ser_ip ${ser_ip} \
                --ser_in ${ser_in} \
                --table ${table} \
                --eqn ${eqn} \
                --outfile ${outfile} \
                $(if ${flg_dep}; then echo "--flg_dep"; fi) \
                $(if ${flg_len}; then echo "--flg_len"; fi) \
                $(if ${flg_mc}; then echo "--flg_mc"; fi) \
                --err_out ${err_out} \
                --nam_job ${nam_job} \
                --env_nam ${env_nam} \
                --scr_met ${scr_met} \
                --scr_alf ${scr_alf}
    fi
else
    :  #TODO
fi
