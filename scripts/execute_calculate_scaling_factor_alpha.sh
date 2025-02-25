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
    
    typ_cvg="alpha"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=true
    threads=8
    ser_ip="$(  ## WARNING: Change search parameters as needed ##
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_bam}" \
            --pattern "${pattern}" \
            --include "IP*" \
            --exclude "${exclude}"
    )"
    ser_in="$(sed 's:IP:in:g' < <(echo "${ser_ip}"))"
    tbl_met="${dir_dat}/raw/docs/measurements_siqchip.tsv"
    eqn="6nd"  # "5"
    fil_out="${dir_out}/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_${typ_cvg}_${eqn}.tsv"
    rnd=24
    err_out="${dir_out}/logs"
    nam_job="calc_sf_alpha_${eqn}"
    slurm=false
    max_job=6
    time="0:30:00"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
# shellcheck disable=SC2269
{
    env_nam="env_protocol"
    dir_scr="${dir_scr}"
    scr_sub="${dir_scr}/submit_calculate_scaling_factor_alpha.sh"
    denom=4
    par_job=""
}

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=1
ser_ip=""
ser_in=""
tbl_met="$(dirname "${dir_scr}")/data/raw/docs/measurements_siqchip.tsv"
eqn="6nd"
fil_out=""
rnd=24
err_out=""
nam_job="calc_sf_alpha_${eqn}"
slurm=false
max_job=6
time="0:30:00"

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  execute_calculate_scaling_factor_alpha.sh
    [--verbose] [--dry_run] --threads <int> --ser_ip <str> --ser_in <str>
    --tbl_met <str> --eqn <str> --fil_out <str> --rnd <int> --err_out <str>
    --nam_job <str> [--slurm] [--max_job <int>] [--time <str>]

Description:
  The driver script 'execute_calculate_scaling_factor_alpha.sh' performs...
  #TODO

Arguments:
   -h, --help     Display this help message and exit (0).
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Perform a dry run without executing commands (optional).
   -t, --threads  Number of threads to use (default: ${threads}).
  -sp, --ser_ip   Comma-separated serialized list of IP coordinate-sorted BAM
                  infiles, including paths.
  -sn, --ser_in   Comma-separated serialized list of input coordinate-sorted
                  BAM infiles, including paths.
  -tb, --tbl_met  #TODO Write description (default: '${tbl_met}')
  -eq, --eqn      #TODO Write description (default: '${eqn}').
  -fo, --fil_out  #TODO Write description
  -fl, --flg_len  #TODO Implement
  -fd, --flg_dep  #TODO Implement
   -r, --rnd      #TODO Write description (default: ${rnd}).
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (required; default: \${dir_out}/err_out).
  -nj, --nam_job  The name of the job (required; default: '${nam_job}').
  -sl, --slurm    Submit jobs to the SLURM scheduler.
  -mj, --max_job  The maximum number of jobs to run at one time (required if
                  '--slurm' is specified, ignored if not; default: ${max_job}).
  -tm, --time     The length of time, in 'h:mm:ss' format, for the SLURM job
                  (required if '--slurm' is specified, ignored if not; default:
                  '${time}').

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
            -tb|--tbl_met) tbl_met="${2}"; shift 2 ;;
            -eq|--eqn)     eqn="${2}";     shift 2 ;;
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

check_supplied_arg -a "${ser_ip}" -n "ser_ip"
check_exists_file_dir "d" "$(dirname "${ser_ip%%[,;]*}")" "ser_ip"

check_supplied_arg -a "${ser_in}" -n "ser_in"
check_exists_file_dir "d" "$(dirname "${ser_in%%[,;]*}")" "ser_in"

check_supplied_arg -a "${tbl_met}" -n "tbl_met"
check_exists_file_dir "f" "${tbl_met}" "tbl_met"

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

#  Parse the '--ser_ip' and '--ser_in' arguments into separate arrays, then
#+ validate the infile value assignments
IFS=',' read -r -a arr_ip <<< "${ser_ip}"
IFS=',' read -r -a arr_in <<< "${ser_in}"

#  Check that each file exists; if not, exit
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
    echo "ser_ip=${ser_ip}"
    echo "ser_in=${ser_in}"
    echo "tbl_met=${tbl_met}"
    echo "eqn=${eqn}"
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
    echo "arr_ip=( ${arr_ip[*]} )"
    echo ""
    echo "arr_in=( ${arr_in[*]} )"
    echo ""
    echo ""
fi

#TODO: Not needed, right?
ser_ip=$(echo "${arr_ip[*]}" | tr ' ' ',')
ser_in=$(echo "${arr_in[*]}" | tr ' ' ',')

if ${verbose}; then
    echo "###############################################################"
    echo "## Variable assignments re- or newly constructed from arrays ##"
    echo "###############################################################"
    echo ""
    echo "ser_ip=\"${ser_ip}\""
    echo ""
    echo "ser_in=\"${ser_in}\""
    echo ""
    echo ""
fi

#  To prevent potential race conditions from concurrent writes,
#+ pre-write the header to the outfile before running SLURM jobs
prt_1="fil_ip\tfil_in\talpha\teqn\t"
prt_2="mass_ip\tmass_in\tvol_all\tvol_in\t"
prt_3="dep_ip\tdep_in\tlen_ip\tlen_in\t"
prt_4="dm_fr_1\tdm_fr_5\tdm_fr_10\tdm_fr_20\tdm_fr_30\t"
prt_5="dm_fr_40\tdm_fr_50\t"
prt_6="dm_nm_1\tdm_nm_5\tdm_nm_10\tdm_nm_20\tdm_nm_30\t"
prt_7="dm_nm_40\tdm_nm_50"

if ${dry_run} || ${verbose}; then
    echo "##################"
    echo "## Table header ##"
    echo "##################"
    echo ""
    echo -e "${prt_1}${prt_2}${prt_3}${prt_4}${prt_5}${prt_6}${prt_7}"
    echo ""
    echo ""
fi

if ! ${dry_run}; then
    echo -e \
        "${prt_1}${prt_2}${prt_3}${prt_4}${prt_5}${prt_6}${prt_7}" \
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
        echo "    --array=1-${#arr_ip[@]}%${max_job} \\"
        echo "    ${scr_sub} \\"
        echo "        --threads ${threads} \\"
        echo "        --ser_ip ${ser_ip} \\"
        echo "        --ser_in ${ser_in} \\"
        echo "        --tbl_met ${tbl_met} \\"
        echo "        --eqn ${eqn} \\"
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
    fi

    if ! ${dry_run}; then        
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
                --tbl_met ${tbl_met} \
                --eqn ${eqn} \
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

        for idx in "${!arr_ip[@]}"; do
            echo \
                "${threads}" \
                "${arr_ip[idx]}" \
                "${arr_in[idx]}" \
                "${tbl_met}" \
                "${eqn}" \
                "${fil_out}" \
                "${rnd}" \
                "${err_out}" \
                "${nam_job}" \
                "${env_nam}" \
                "${dir_scr}" \
                    >> "${config}"
        done

        cmd="bash \"${scr_sub}\" -t {1} -sp {2} -sn {3} -tm {4} -eq {5} -fo {6} -r {7} -eo {8} -nj {9} -en {10} -ds {11}"

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
        if ${dry_run} || ${verbose}; then
            echo "######################"
            echo "## Serial execution ##"
            echo "######################"
            echo ""
            echo "bash ${scr_sub} \\"
            echo "     -t ${threads} \\"
            echo "    -sp ${ser_ip} \\"
            echo "    -sn ${ser_in} \\"
            echo "    -tm ${tbl_met} \\"
            echo "    -eq ${eqn} \\"
            echo "    -fo ${fil_out} \\"
            echo "     -r ${rnd} \\"
            echo "    -eo ${err_out} \\"
            echo "    -nj ${nam_job} \\"
            echo "    -en ${env_nam} \\"
            echo "    -ds ${dir_scr} \\"
            echo "         > ${err_out}/${nam_job}.stdout.txt \\"
            echo "        2> ${err_out}/${nam_job}.stderr.txt"
            echo ""
            echo ""
        fi

        if ! ${dry_run}; then
            bash "${scr_sub}" \
                 -t "${threads}" \
                -sp "${ser_ip}" \
                -sn "${ser_in}" \
                -tm "${tbl_met}" \
                -eq "${eqn}" \
                -fo "${fil_out}" \
                 -r "${rnd}" \
                -eo "${err_out}" \
                -nj "${nam_job}" \
                -en "${env_nam}" \
                -ds "${dir_scr}" \
                     > "${err_out}/${nam_job}.stdout.txt" \
                    2> "${err_out}/${nam_job}.stderr.txt"
        fi
    fi
fi
