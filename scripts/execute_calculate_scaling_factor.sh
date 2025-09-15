#!/bin/bash

#  execute_calculate_scaling_factor.sh
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
    check_arrays_lengths \
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
    handle_env \
    print_parallel_info \
    reset_max_job \
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
    dir_scr="${dir_rep}/scripts"

    #  Set alignment parameters
    aligner="bowtie2"
    a_type="global"
    req_flg=true
    flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
    mapq=1
    str_det="${aligner}_${a_type}_flag-${flg}_mapq-${mapq}"

    #  Define data directories
    dir_dat="${dir_rep}/data"
    dir_pro="${dir_dat}/processed"
    dir_aln="${dir_pro}/align_reads"  # align_fastqs
    dir_det="${dir_aln}/${str_det}"
    dir_bam="${dir_det}/sc"

    #  Set output directories
    dir_sig="${dir_pro}/compute_signal/${str_det}"
    dir_out="${dir_sig}/tables"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=true
    threads=6
    mode="alpha"
    ser_mip="$(  ## WARNING: Change search parameters as needed ##
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_bam}" \
            --pattern "*.bam" \
            --include "IP*" \
            --exclude "*Brn1*"
    )"
    ser_min="${ser_mip//IP_/in_}"

    fil_out="${dir_out}/test_${mode}"  # ChIP_WT_G1-G2M-Q_Hho1-Hmo1_${mode}
    nam_job="calc_sf_${mode}"

    if [[ "${mode}" == "spike" ]]; then
        ser_sip="${ser_mip//sc/sp}"
        ser_sin="${ser_sip//IP_/in_}"
        fil_out="${fil_out}.tsv"
    elif [[ "${mode}" == "alpha" ]]; then
        tbl_met="${dir_dat}/raw/docs/measurements_siqchip.tsv"
        eqn="6nd"
        fil_out="${fil_out}_${eqn}.tsv"
        nam_job="calc_sf_${eqn}"
    fi

    rnd=24
    err_out="${dir_out}/logs"
    max_job=2
    slurm=false
    time="0:30:00"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
# shellcheck disable=SC2269
{
    env_nam="env_protocol"
    dir_scr="${dir_scr}"
    scr_hdr="${dir_scr}/write_header.sh"
    scr_sub="${dir_scr}/submit_calculate_scaling_factor.sh"
    par_job=""
}

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=1
mode="alpha"
ser_mip=""
ser_min=""
ser_sip=""
ser_sin=""
tbl_met="$(dirname "${dir_scr}")/data/raw/docs/measurements_siqchip.tsv"
eqn="6nd"
fil_out=""
rnd=24
err_out=""
nam_job="calc_sf_${mode}_${eqn}"
max_job=6
slurm=false
time="0:30:00"

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  execute_calculate_scaling_factor.sh
    [--verbose] [--dry_run] --threads <int> --mode <str> --ser_mip <str>
    --ser_min <str> [--ser_sip <str>] [--ser_sin <str>] --tbl_met <str>
    --eqn <str> --fil_out <str> --rnd <int> --err_out <str> --nam_job <str>
    --max_job <int> [--slurm] [--time <str>]

Description:
  Coordinates the calculation of siQ-ChIP alpha or spike-in scaling factors for
  ChIP-seq data.

  For alpha scaling factors, this involves... #TODO

  For spike-in scaling factors, this involves counting alignments in "main" or
  "experimental" (e.g., S. cerevisiae) and "spike-in" (e.g., S. pombe) BAM
  files, calculating scaling coefficients using custom Python scripts, and
  outputting results to a specified tab-delimited table file.

  Facilitates parallelized job submission with either SLURM or GNU Parallel, or
  serialized job submission.

Arguments:
   -h, --help     Display this help message and exit
   -v, --verbose  Run script in 'verbose mode' (optional)
  -dr, --dry_run  Perform a dry run without executing commands (optional)
   -t, --threads  Number of threads to use (default: '${threads}')
   -m, --mode     Type of scaling factor to compute: 'alpha' or 'spike'
                  (default: '${mode}')
  -mp, --ser_mip  Comma-separated string of "main" organism IP BAM files
  -mn, --ser_min  Comma-separated string of "main" organism input BAM files
  -sp, --ser_sip  Comma-separated string of "spike-in" organism IP BAM files
                  (required if '--mode spike', ignored if not)
  -sn, --ser_sin  Comma-separated string of "spike-in" organism input BAM files
                  (required if '--mode spike', ignored if not)
  -tb, --tbl_met  Tab-delimited input file of siQ-ChIP metadata metrics
                  (required if '--mode alpha', ignored if not)
  -eq, --eqn      Alpha equation to compute: '5', '5nd', '6', '6nd' (required
                  if '--mode alpha', ignored if not; default: '${eqn}')
  -fo, --fil_out  Tab-delimited text output file in which scaling factors and
                  additional metrics and values are written.
  -fl, --flg_len  #TODO: Implement this
  -fd, --flg_dep  #TODO: Implement this
   -r, --rnd      Number of decimal places for rounding scaling factors and
                  minimum input depth factors (default: '${rnd}')
  -eo, --err_out  Directory for stderr and stdout output files (default:
                  '\${dir_out}/err_out')
  -nj, --nam_job  Name of job (default: '${nam_job}')
  -mj, --max_job  Maximum number of jobs to run concurrently (default: '${max_job}')
                    - If '--slurm' is specified, controls SLURM array tasks
                    - If '--slurm' is not specified:
                      + If 'max_job' is greater than 1, jobs run in parallel
                        via GNU Parallel
                      + If 'max_job' is 1, jobs run sequentially (serial mode)
  -sl, --slurm    Submit jobs to SLURM scheduler
  -tm, --time     Length of time, in 'h:mm:ss' format, for SLURM job (required
                  if '--slurm' is specified, ignored if not; default: '${time}')

Dependencies:
  - awk
  - Bash or Zsh
  - Conda
  - GNU Parallel (when '--slurm' is not specified and 'threads' > 1)
  - Python
  - Samtools
  - SLURM (when '--slurm' is specified)

Notes:
  - If '--dry_run' is enabled, the script will output commands that would be
    executed without actually running them.
  - For non-SLURM job submissions, GNU Parallel is used. This requires
    specifying a set number of parallel jobs ('par_job') and the number of
    threads per job ('threads'). The number of parallel jobs is determined by
    dividing the user-specified 'threads' value by a denominator,
    which is determined... #TODO: Explain this. The value of 'threads' is then
    reset to... #TODO: Explain this. For example, if the user-specified...
    #TODO: Explain this.
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
    echo ""
    if [[ -z "${1:-}" ]]; then exit_1; else exit_0; fi
fi

if ${interactive:-false}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -v|--verbose) verbose=true;   shift 1 ;;
            -dr|--dry_run) dry_run=true;   shift 1 ;;
             -t|--threads) threads="${2}"; shift 2 ;;
             -m|--mode)    mode="${2}";    shift 2 ;;
            -mp|--ser_mip) ser_mip="${2}"; shift 2 ;;
            -mn|--ser_min) ser_min="${2}"; shift 2 ;;
            -sp|--ser_sip) ser_sip="${2}"; shift 2 ;;
            -sn|--ser_sin) ser_sin="${2}"; shift 2 ;;
            -tb|--tbl_met) tbl_met="${2}"; shift 2 ;;
            -eq|--eqn)     eqn="${2}";     shift 2 ;;
            -fo|--fil_out) fil_out="${2}"; shift 2 ;;
             -r|--rnd)     rnd="${2}";     shift 2 ;;
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

check_supplied_arg -a "${dir_scr}" -n "dir_scr"
check_exists_file_dir "d" "${dir_scr}" "dir_scr"

check_supplied_arg -a "${scr_hdr}" -n "scr_hdr"
check_exists_file_dir "f" "${scr_hdr}" "scr_hdr"

check_supplied_arg -a "${scr_sub}" -n "scr_sub"
check_exists_file_dir "f" "${scr_sub}" "scr_sub"

check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

case "${mode}" in
    alpha|spike) : ;;
    *)
        echo_error \
            "Scaling factor mode ('--mode') was assigned '${mode}' but must" \
            "be 'alpha' or 'spike'."
        exit_1
        ;;
esac

check_supplied_arg -a "${ser_mip}" -n "ser_mip"
check_exists_file_dir "d" "$(dirname "${ser_mip%%[,;]*}")" "ser_mip"
check_str_delim "ser_mip" "${ser_mip}"

check_supplied_arg -a "${ser_min}" -n "ser_min"
check_exists_file_dir "d" "$(dirname "${ser_min%%[,;]*}")" "ser_min"
check_str_delim "ser_min" "${ser_min}"

if [[ "${mode}" == "spike" ]]; then
    check_supplied_arg -a "${ser_sip}" -n "ser_sip"
    check_exists_file_dir "d" "$(dirname "${ser_sip%%[,;]*}")" "ser_sip"
    check_str_delim "ser_sip" "${ser_sip}"

    check_supplied_arg -a "${ser_sin}" -n "ser_sin"
    check_exists_file_dir "d" "$(dirname "${ser_sin%%[,;]*}")" "ser_sin"
    check_str_delim "ser_sin" "${ser_sin}"

    unset tbl_met eqn
elif [[ "${mode}" == "alpha" ]]; then
    check_supplied_arg -a "${tbl_met}" -n "tbl_met"
    check_exists_file_dir "f" "${tbl_met}" "tbl_met"

    check_supplied_arg -a "${eqn}" -n "eqn"
    case "${eqn}" in
        5|5nd|6|6nd) : ;;
        *)
            echo_error \
                "Equation ('--eqn') was assigned '${eqn}' but must be '5'," \
                "'5nd', '6', or '6nd'."
            exit_1
        ;;
    esac

    unset ser_sip ser_sin
fi

check_supplied_arg -a "${fil_out}" -n "fil_out"
check_exists_file_dir "d" "$(dirname "${fil_out}")" "fil_out"

check_supplied_arg -a "${rnd}" -n "rnd"
check_int_pos "${rnd}" "rnd"

if [[ -z "${err_out}" ]]; then
    err_out="${dir_out}/err_out"
fi
check_exists_file_dir "d" "${err_out}" "err_out"

check_supplied_arg -a "${nam_job}" -n "nam_job"


#  Parse and validate input vector elements -----------------------------------
IFS=',' read -r -a arr_mip <<< "${ser_mip}"
IFS=',' read -r -a arr_min <<< "${ser_min}"

for file in "${arr_mip[@]}" "${arr_min[@]}"; do
    check_exists_file_dir "f" "${file}" "file"
done
unset file

check_arrays_lengths "arr_mip" "arr_min"

#  Handle spike-in files if specified
if [[ "${mode}" == "spike" ]]; then
    IFS=',' read -r -a arr_sip <<< "${ser_sip}"
    IFS=',' read -r -a arr_sin <<< "${ser_sin}"

    for file in "${arr_sip[@]}" "${arr_sin[@]}"; do
        check_exists_file_dir "f" "${file}" "file"
    done
    unset file

    check_arrays_lengths "arr_sip" "arr_sin"
    check_arrays_lengths "arr_sip" "arr_mip"
fi


#  Parse job execution parameters ---------------------------------------------
if ${slurm:-false}; then
    check_supplied_arg -a "${max_job}" -n "max_job"
    check_int_pos "${max_job}" "max_job"

    max_job=$(reset_max_job "${max_job}" "${#arr_mip[@]}")
    
    check_supplied_arg -a "${time}" -n "time"
    check_format_time "${time}"
else
    IFS=';' read -r threads par_job < <(
        set_params_parallel "${threads}" "${max_job}" "${par_job}"
    )
    unset max_job time

    check_supplied_arg -a "${par_job}" -n "par_job"
    check_int_pos "${par_job}" "par_job"
fi

#  Debug parallelization information
if [[ "${mode}" == "spike" ]]; then
    unset extra && typeset -a extra
    extra=( "arr_sip" "arr_sin" )
fi

print_parallel_info \
    "${slurm}" "${max_job:-UNSET}" "${par_job}" "${threads}" \
    "arr_mip" "arr_min" "${extra[@]}"

if [[ "${mode}" == "spike" ]]; then unset extra; fi


#  Activate environment and check that dependencies are in PATH ---------------
handle_env "${env_nam}" > /dev/null

check_program_path awk
check_program_path python
check_program_path samtools

if ${slurm:-false}; then
    check_program_path sbatch
elif [[ ${threads} -gt 1 ]]; then
    check_program_path parallel
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
    echo "scr_hdr=${scr_hdr}"
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
    echo "ser_mip=${ser_mip}"
    echo "ser_min=${ser_min}"
    echo "ser_sip=${ser_sip:-UNSET}"
    echo "ser_sin=${ser_sin:-UNSET}"
    echo "tbl_met=${tbl_met:-UNSET}"
    echo "eqn=${eqn:-UNSET}"
    echo "fil_out=${fil_out}"
    echo "rnd=${rnd}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo "max_job=${max_job:-UNSET}"
    echo "slurm=${slurm}"
    echo "time=${time:-UNSET}"
    echo ""
    echo ""
    echo "###################################"
    echo "## Arrays derived from variables ##"
    echo "###################################"
    echo ""
    echo "arr_mip=( ${arr_mip[*]} )"
    echo ""
    echo "arr_min=( ${arr_min[*]} )"

    if [[ "${mode}" == "spike" ]]; then
        echo ""
        echo "arr_sip=( ${arr_sip[*]} )"
        echo ""
        echo "arr_sin=( ${arr_sin[*]} )"
    fi

    echo ""
    echo ""
fi

# shellcheck disable=SC1083,SC2157,SC2046,SC2086
if ${slurm:-false}; then
    #  If '--slurm' was specified, run jobs in parallel via individual job
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
        echo "        --dir_scr ${dir_scr} \\"
        echo "        --env_nam ${env_nam} \\"
        echo "        --threads ${threads} \\"
        echo "        --mode ${mode} \\"
        echo "        --ser_mip ${ser_mip} \\"
        echo "        --ser_min ${ser_min} \\"

        if [[ "${mode}" == "spike" ]]; then
            echo "        --ser_sip ${ser_sip} \\"
            echo "        --ser_sin ${ser_sin} \\"
        else
            echo "        --tbl_met ${tbl_met} \\"
            echo "        --eqn ${eqn} \\"
        fi

        echo "        --fil_out ${fil_out} \\"
        echo "        --rnd ${rnd} \\"
        echo "        --err_out ${err_out} \\"
        echo "        --nam_job ${nam_job} \\"
        echo ""
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
            --array=1-${#arr_mip[@]}%${max_job} \
            ${scr_sub} \
                --dir_scr ${dir_scr} \
                --env_nam ${env_nam} \
                --threads ${threads} \
                --mode ${mode} \
                --ser_mip ${ser_mip} \
                --ser_min ${ser_min} \
                $(
                    if [[ "${mode}" == "spike" ]]; then
                        echo "--ser_sip ${ser_sip} --ser_sin ${ser_sin}"
                    else
                        echo "--tbl_met ${tbl_met} --eqn ${eqn}"
                    fi
                ) \
                --fil_out ${fil_out} \
                --rnd ${rnd} \
                --err_out ${err_out} \
                --nam_job ${nam_job}
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
                "${dir_scr}" \
                "${env_nam}" \
                "${threads}" \
                "${mode}" \
                "${arr_mip[idx]}" \
                "${arr_min[idx]}" \
                $(
                    if [[ "${mode}" == "spike" ]]; then
                        echo "${arr_sip[idx]}"
                        echo "${arr_sin[idx]}"
                    else
                        echo "${tbl_met}"
                        echo "${eqn}"
                    fi
                ) \
                "${fil_out}" \
                "${rnd}" \
                "${err_out}" \
                "${nam_job}" \
                    >> "${config}"
        done

        cmd="bash ${scr_sub}"
        cmd+=" -ds {1}"
        cmd+=" -en {2}"
        cmd+="  -t {3}"
        cmd+="  -m {4}"
        cmd+=" -mp {5}"
        cmd+=" -mn {6}"

        if [[ "${mode}" == "spike" ]]; then
            cmd+=" -sp {7}"
            cmd+=" -sn {8}"
        else
            cmd+=" -tm {7}"
            cmd+=" -eq {8}"
        fi

        cmd+=" -fo {9}"
        cmd+="  -r {10}"
        cmd+=" -eo {11}"
        cmd+=" -nj {12}"
        cmd+="  > {11}/{12}_par.{5/.}.stdout.txt"  # 'scr_sub' stdout log
        cmd+=" 2> {11}/{12}_par.{5/.}.stderr.txt"  # 'scr_sub' stderr log

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
            echo "    -sp ${ser_mip} \\"
            echo "    -sn ${ser_min} \\"
            echo "    -tm ${tbl_met} \\"
            echo "    -eq ${eqn} \\"
            echo "    -fo ${fil_out} \\"
            echo "     -r ${rnd} \\"
            echo "    -eo ${err_out} \\"
            echo "    -nj ${nam_job} \\"
            echo "    -en ${env_nam} \\"
            echo "    -ds ${dir_scr} \\"
            echo "         > ${err_out}/${nam_job}_ser.stdout.txt \\"
            echo "        2> ${err_out}/${nam_job}_ser.stderr.txt"
            echo ""
            echo ""
        fi

        if ! ${dry_run}; then
            bash "${scr_sub}" \
                 -t "${threads}" \
                -sp "${ser_mip}" \
                -sn "${ser_min}" \
                -tm "${tbl_met}" \
                -eq "${eqn}" \
                -fo "${fil_out}" \
                 -r "${rnd}" \
                -eo "${err_out}" \
                -nj "${nam_job}" \
                -en "${env_nam}" \
                -ds "${dir_scr}" \
                     > "${err_out}/${nam_job}_ser.stdout.txt" \
                    2> "${err_out}/${nam_job}_ser.stderr.txt"
        fi
    fi
fi
