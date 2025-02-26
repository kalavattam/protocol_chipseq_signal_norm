#!/bin/bash

#  execute_calculate_scaling_factor_spike.sh
#  KA


#  Run script in interactive/test mode (true) or command-line mode (false)
interactive=false

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive}; then set -euo pipefail; fi

#  Set the path to the "scripts" directory
if ${interactive}; then
    ## WARNING: If 'interactive=true', change path as needed ##
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
    details="${aligner}_${a_type}_flag-${flg}_mapq-${mapq}"

    dir_aln="${dir_pro}/align_reads/${details}"
    dir_bam="${dir_aln}/sc"
    dir_cvg="${dir_pro}/compute_coverage/${details}"
    dir_out="${dir_cvg}/tables"

    pattern="*.bam"
    exclude="*Brn1*"
    
    typ_cvg="spike"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=true
    threads=8
    ser_mip="$(  ## WARNING: Change search parameters as needed ##
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_bam}" \
            --pattern "${pattern}" \
            --include "IP*" \
            --exclude "${exclude}"
    )"
    ser_sip="$(sed 's:sc:sp:g' < <(echo "${ser_mip}"))"
    ser_min="$(sed 's:IP:in:g' < <(echo "${ser_mip}"))"
    ser_sin="$(sed 's:IP:in:g' < <(echo "${ser_sip}"))"
    fil_out="${dir_out}/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_${typ_cvg}.tsv"
    rnd=24
    err_out="${dir_out}/logs"
    nam_job="calc_sf_spike"
    slurm=true
    max_job=6
    time="0:30:00"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
# shellcheck disable=SC2269
{
    env_nam="env_protocol"
    dir_scr="${dir_scr}"
    scr_sub="${dir_scr}/submit_calculate_scaling_factor_spike.sh"
    denom=4
    par_job=""
}

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=1
ser_mip=""
ser_sip=""
ser_min=""
ser_sin=""
fil_out=""
rnd=24
err_out=""
nam_job="calc_sf_spike"
slurm=false
max_job=6
time="0:30:00"

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  execute_calculate_scaling_factor_spike.sh
    [--verbose] [--dry_run] --threads <int> --ser_mip <str> --ser_sip <str>
    --ser_min <str> --ser_sin <str> --fil_out <str> --rnd <int> --err_out <str>
    --nam_job <str> [--slurm] [--max_job <int>] [--time <str>]

Description:
  The driver script 'execute_calculate_scaling_factor_spike.sh' coordinates the
  calculation of spike-in-derived scaling factors for ChIP-seq data. This
  involves counting alignments in "main" (S. cerevisiae) and spike-in
  (S. pombe) BAM files, calculating scaling coefficients using a custom Python
  script, and outputting results to a specified table file.

  [#TODO: Note on calculating 'dep_min' for modes "frag" and "norm" at common
  bin sizes.]

  The script facilitates serial processing or parallelized processing with
  SLURM or GNU Parallel.

Arguments:
   -h, --help     Display this help message and exit (0).
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Perform a dry run without executing commands (optional).
   -t, --threads  Number of threads to use (default: ${threads}).
  -mp, --ser_mip  Comma-separated list of "main" model organism IP BAM files
                  that are coordinate-sorted.
  -sp, --ser_sip  Comma-separated serialized string of "spike-in" organism IP
                  BAM files.
  -mn, --ser_min  Comma-separated serialized string of "main" organism input
                  BAM files.
  -sn, --ser_sin  Comma-separated serialized string of "spike-in" organism
                  input BAM files.
   -o, --fil_out  Tab-delimited text outfile in which calculated
                  spike-in-derived scaling factors and, optionally, additional
                  metrics are saved. If the file already exists, new data will
                  be appended.
  -fd, --flg_dep  #TODO
   -r, --rnd      Number of decimal places for rounding the spike-in scaling
                  factor (default: ${rnd}).
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (default: \${dir_out}/err_out).
  -nj, --nam_job  The name of the job (default: '${nam_job}').
  -sl, --slurm    Submit jobs to the SLURM scheduler.
  -mj, --max_job  The maximum number of jobs to run at one time (required if
                  --slurm is specified, ignored if not; default: ${max_job}).
  -tm, --time     The length of time, in 'h:mm:ss' format, for the SLURM job
                  (required if '--slurm' is specified, ignored if not; default:
                  '${time}').

Dependencies:
  - Programs
    + awk
    + Bash or Zsh
    + Conda
    + GNU Parallel (if '--slurm' is not specified)
    + Python
    + Samtools
    + SLURM (if '--slurm' is specified)
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
  - If '--dry_run' is enabled, the script will output commands that would be
    executed without actually running them.
  - For non-SLURM job submissions, GNU Parallel is used. This requires
    specifying a set number of parallel jobs ('par_job') and the number of
    threads per job ('threads'). The number of parallel jobs is determined by
    dividing the user-specified 'threads' value by a denominator ('denom'),
    which is set (hardcoded) to 4. The value of 'threads' is then reset to
    'denom'. For example, if the user-specified 'threads' value is 8, 'par_job'
    would be set to 2 (8 divided by 4) and 'threads' would be reset to 4. Thus,
    2 jobs would be run by GNU Parallel (because 'par_job=2'), each with 4
    threads (because 'threads=4').
  - #TODO

Examples:
  \`\`\`
  #TODO
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
            -mp|--ser_mip) ser_mip="${2}"; shift 2 ;;
            -sp|--ser_sip) ser_sip="${2}"; shift 2 ;;
            -mn|--ser_min) ser_min="${2}"; shift 2 ;;
            -sn|--ser_sin) ser_sin="${2}"; shift 2 ;;
            -fo|--fil_out) fil_out="${2}"; shift 2 ;;
             -r|--rnd)     rnd="${2}";     shift 2 ;;
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

check_supplied_arg -a "${dir_scr}" -n "dir_scr"
check_exists_file_dir "d" "${dir_scr}" "dir_scr"

check_supplied_arg -a "${scr_sub}" -n "scr_sub"
check_exists_file_dir "f" "${scr_sub}" "scr_sub"

check_supplied_arg -a "${denom}" -n "denom"
check_int_pos "${denom}" "denom"

check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

check_supplied_arg -a "${ser_mip}" -n "ser_mip"
check_exists_file_dir "d" "$(dirname "${ser_mip%%[,;]*}")" "ser_mip"

check_supplied_arg -a "${ser_sip}" -n "ser_sip"
check_exists_file_dir "d" "$(dirname "${ser_sip%%[,;]*}")" "ser_sip"

check_supplied_arg -a "${ser_min}" -n "ser_min"
check_exists_file_dir "d" "$(dirname "${ser_min%%[,;]*}")" "ser_min"

check_supplied_arg -a "${ser_sin}" -n "ser_sin"
check_exists_file_dir "d" "$(dirname "${ser_sin%%[,;]*}")" "ser_sin"

check_supplied_arg -a "${fil_out}" -n "fil_out"
check_exists_file_dir "d" "$(dirname "${fil_out}")" "fil_out"

check_supplied_arg -a "${rnd}" -n "rnd"
check_int_pos "${rnd}" "rnd"

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

#  Parse the '--ser_mip', etc. arguments into arrays, then validate the file
#+ value assignments
IFS=',' read -r -a arr_mip <<< "${ser_mip}"  # unset arr_mip
IFS=',' read -r -a arr_sip <<< "${ser_sip}"  # unset arr_sip
IFS=',' read -r -a arr_min <<< "${ser_min}"  # unset arr_min
IFS=',' read -r -a arr_sin <<< "${ser_sin}"  # unset arr_sin

#  Check that each file exists; if not, exit
for file in \
    "${arr_mip[@]}" "${arr_sip[@]}" \
    "${arr_min[@]}" "${arr_sin[@]}"
do
    check_exists_file_dir "f" "${file}" "file"
done
unset file

#  Reset 'max_job' if it is greater than the array lengths
if ${slurm}; then
    max_job=$(reset_max_job "${max_job}" "${#arr_mip[@]}")
    #FIXME: Error when 'max_job' is unset or empty
fi


#  Do the main work ===========================================================
#  Report argument variable assignments if in "verbose mode"
if ${verbose}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "env_nam=${env_nam}"
    echo "dir_scr=${dir_scr}"
    echo "scr_sub=${scr_sub}"
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
    echo "ser_mip=${ser_mip}"
    echo "ser_sip=${ser_sip}"
    echo "ser_min=${ser_min}"
    echo "ser_sin=${ser_sin}"
    echo "fil_out=${fil_out}"
    echo "rnd=${rnd}"
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
    echo "arr_mip=( ${arr_mip[*]} )"
    echo ""
    echo "arr_sip=( ${arr_sip[*]} )"
    echo ""
    echo "arr_min=( ${arr_min[*]} )"
    echo ""
    echo "arr_sin=( ${arr_sin[*]} )"
    echo ""
    echo ""
fi

#TODO: Not needed, right?
ser_mip=$(echo "${arr_mip[*]}" | tr ' ' ',')
ser_sip=$(echo "${arr_sip[*]}" | tr ' ' ',')
ser_min=$(echo "${arr_min[*]}" | tr ' ' ',')
ser_sin=$(echo "${arr_sin[*]}" | tr ' ' ',')

if ${verbose}; then
    echo "###############################################################"
    echo "## Variable assignments re- or newly constructed from arrays ##"
    echo "###############################################################"
    echo ""
    echo "ser_mip=\"${ser_mip}\""
    echo ""
    echo "ser_sip=\"${ser_sip}\""
    echo ""
    echo "ser_min=\"${ser_min}\""
    echo ""
    echo "ser_sin=\"${ser_sin}\""
    echo ""
    echo ""
fi

#  To prevent potential race conditions from concurrent writes,
#+ pre-write the header to the outfile before running SLURM jobs
prt_1="main_ip\tspike_ip\tmain_in\tspike_in\tsf\t"
prt_2="num_mp\tnum_sp\tnum_mn\tnum_sn\t"
prt_3="dm_fr_1\tdm_fr_5\tdm_fr_10\tdm_fr_20\tdm_fr_30\t"
prt_4="dm_fr_40\tdm_fr_50\t"
prt_5="dm_nm_1\tdm_nm_5\tdm_nm_10\tdm_nm_20\tdm_nm_30\t"
prt_6="dm_nm_40\tdm_nm_50"

if ${dry_run} || ${verbose}; then
    echo "##################"
    echo "## Table header ##"
    echo "##################"
    echo ""
    echo -e "${prt_1}${prt_2}${prt_3}${prt_4}${prt_5}${prt_6}"
    echo ""
    echo ""
fi

if ! ${dry_run}; then
    echo -e \
        "${prt_1}${prt_2}${prt_3}${prt_4}${prt_5}${prt_6}" \
            >> "${fil_out}"
fi


# shellcheck disable=SC1083,SC2157,SC2046,SC2086
if ${slurm}; then
    #  If --slurm was specified, run jobs in parallel via individual job
    #+ submissions to SLURM
    if ${dry_run} || ${verbose}; then
        echo "#####################"
        echo "## SLURM execution ##"
        echo "#####################"
        echo ""
        echo "sbatch \\"
        echo "    --job-name=${nam_job} \\"
        echo "    --nodes=1 \\"
        echo "    --cpus-per-task=${threads} \\"
        echo "    --time=${time} \\"
        echo "    --error=${err_out}/${nam_job}.%A-%a.stderr.txt \\"
        echo "    --output=${err_out}/${nam_job}.%A-%a.stdout.txt \\"
        echo "    --array=1-${#arr_mip[@]}%${max_job} \\"
        echo "    ${scr_sub} \\"
        echo "        --threads ${threads} \\"
        echo "        --ser_mip ${ser_mip} \\"
        echo "        --ser_sip ${ser_sip} \\"
        echo "        --ser_min ${ser_min} \\"
        echo "        --ser_sin ${ser_sin} \\"
        echo "        --fil_out ${fil_out} \\"
        echo "        --rnd ${rnd} \\"
        echo "        --err_out ${err_out} \\"
        echo "        --nam_job ${nam_job} \\"
        echo "        --env_nam ${env_nam} \\"
        echo "        --dir_scr ${dir_scr}"
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
        sbatch \
            --job-name=${nam_job} \
            --nodes=1 \
            --cpus-per-task=${threads} \
            --time=${time} \
            --error=${err_out}/${nam_job}.%A-%a.stderr.txt \
            --output=${err_out}/${nam_job}.%A-%a.stdout.txt \
            --array=1-${#arr_mip[@]}%${max_job} \
            ${scr_sub} \
                --threads ${threads} \
                --ser_mip ${ser_mip} \
                --ser_sip ${ser_sip} \
                --ser_min ${ser_min} \
                --ser_sin ${ser_sin} \
                --fil_out ${fil_out} \
                --rnd ${rnd} \
                --err_out ${err_out} \
                --nam_job ${nam_job} \
                --env_nam ${env_nam} \
                --dir_scr ${dir_scr}
    fi
else
    #  GNU Parallel execution
    if [[ ${threads} -gt 1 ]]; then
        config="${err_out}/${nam_job}.config_parallel.txt"

        if [[ -f "${config}" ]]; then rm "${config}"; fi
        touch "${config}" || {
            echo_error "Failed to create a GNU Parallel configuration file."
            exit_1
        }

        for idx in "${!arr_mip[@]}"; do
            echo \
                "${threads}" \
                "${arr_mip[idx]}" \
                "${arr_sip[idx]}" \
                "${arr_min[idx]}" \
                "${arr_sin[idx]}" \
                "${fil_out}" \
                "${rnd}" \
                "${err_out}" \
                "${nam_job}" \
                "${env_nam}" \
                "${dir_scr}" \
                    >> "${config}"
        done
        # cat "${config}"

        cmd="bash \"${scr_sub}\" -t {1} -mp {2} -sp {3} -mn {4} -sn {5} -fo {6} -r {7} -eo {8} -nj {9} -en {10} -ds {11}"

        if ${dry_run} || ${verbose}; then
            echo "############################"
            echo "## GNU Parallel execution ##"
            echo "############################"
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
        #  Serial execution
        pth_std="${err_out}/${nam_job}"
        if ${dry_run} || ${verbose}; then
            echo "######################"
            echo "## Serial execution ##"
            echo "######################"
            echo ""
            echo "bash ${scr_sub} \\"
            echo "     -t ${threads} \\"
            echo "    -mp ${ser_mip} \\"
            echo "    -sp ${ser_sip} \\"
            echo "    -mn ${ser_min} \\"
            echo "    -sn ${ser_sin} \\"
            echo "    -fo ${fil_out} \\"
            echo "     -r ${rnd} \\"
            echo "    -eo ${err_out} \\"
            echo "    -nj ${nam_job} \\"
            echo "    -en ${env_nam} \\"
            echo "    -ds ${dir_scr} \\"
            echo "         > ${pth_std}.stdout.txt \\"
            echo "        2> ${pth_std}.stderr.txt"
            echo ""
            echo ""
        fi

        if ! ${dry_run}; then
            bash "${scr_sub}" \
                 -t "${threads}" \
                -mp "${ser_mip}" \
                -sp "${ser_sip}" \
                -mn "${ser_min}" \
                -sn "${ser_sin}" \
                -fo "${fil_out}" \
                 -r "${rnd}" \
                -eo "${err_out}" \
                -nj "${nam_job}" \
                -en "${env_nam}" \
                -ds "${dir_scr}" \
                     > "${pth_std}.stdout.txt" \
                    2> "${pth_std}.stderr.txt"
        fi
    fi
fi
