#!/bin/bash

#  execute_compute_coverage_ratio.sh
#  KA


#  Run script in interactive/test mode (true) or command-line mode (false)
interactive=true

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
    source "${dir_fnc}/check_array_files.sh"
    source "${dir_fnc}/check_arrays_lengths.sh"
    source "${dir_fnc}/check_exists_file_dir.sh"
    source "${dir_fnc}/check_flt_pos.sh"
    source "${dir_fnc}/check_format_time.sh"
    source "${dir_fnc}/check_int_pos.sh"
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/check_str_delim.sh"
    source "${dir_fnc}/check_supplied_arg.sh"
    source "${dir_fnc}/debug_array_contents.sh"
    source "${dir_fnc}/echo_error.sh"
    source "${dir_fnc}/echo_warning.sh"
    source "${dir_fnc}/exit_0.sh"
    source "${dir_fnc}/exit_1.sh"
    source "${dir_fnc}/handle_env.sh"
    source "${dir_fnc}/populate_array_empty.sh"
    source "${dir_fnc}/reset_max_job.sh"
}


#  Set up paths, values, and parameters for interactive mode
function set_interactive() {
    #  Set hardcoded paths, values, etc.
    ## WARNING: If interactive=true, change values as needed ##
    dir_bas="${HOME}/repos"
    dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
    dir_scr="${dir_rep}/scripts"
    dir_dat="${dir_rep}/data"
    dir_pro="${dir_dat}/processed"

    aligner="bowtie2"
    a_type="global"
    req_flg=true
    flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
    mapq=1
    det_bam="flag-${flg}_mapq-${mapq}"
    det_cov="${aligner}_${a_type}_${det_bam}"

    dir_cov="${dir_pro}/compute_coverage"
    dir_nrm="${dir_cov}/${det_cov}/norm/tracks"
    dir_alf="${dir_cov}/${det_cov}/alpha/tracks"

    pattern="*.bdg.gz"
    exclude="*Brn1*"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=true
    fil_ip="$(
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_nrm}" \
            --pattern "${pattern}" \
            --include "IP*" \
            --exclude "${exclude}"
    )"
    fil_in="$(
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_nrm}" \
            --pattern "${pattern}" \
            --include "in*" \
            --exclude "${exclude}"
    )"
    dir_out="${dir_alf}"
    typ_out="bdg.gz"
    scl_fct=""
    dep_min=""  #TODO
    log2=true
    rnd=24
    err_out="${dir_alf}/logs"
    nam_job="compute_coverage_ratio"
    slurm=true
    max_job=6
    time="0:30:00"
}

#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_analyze"
scr_sub="${dir_scr}/submit_compute_coverage_ratio.sh"
scr_cvg="${dir_scr}/compute_coverage_ratio.py"
denom=4

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
fil_ip=""
fil_in=""
dir_out=""
typ_out="bdg.gz"
scl_fct="none"
dep_min=""
log2=false
rnd=24
err_out=""
nam_job="compute_coverage_ratio"
slurm=false
max_job=6
time="0:30:00"

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  execute_compute_coverage_ratio.sh
    [--verbose] [--dry_run] --fil_ip <str> --fil_in <str> --dir_out <str>
    --typ_out <str> --scl_fct <flt> --dep_min <flt> [--log2] --rnd <int>
    --err_out <str> --nam_job <str> --max_job <int> [--slurm] [--time <str>]

Description:
  The driver script 'execute_compute_coverage_ratio.sh' automates the
  computation of ratio coverage tracks for ChIP-seq data using BEDGRAPH signal
  (coverage score) tracks. It can be used to generate either siQ- or spike-in-
  scaled coverage values, or log2(IP/input) enrichment tracks. When used with
  '--scl_fct', the script applies a user-specified scaling factor, allowing for
  normalization based on the siQ-ChIP method or external spike-in controls.
  Alternatively, enabling flag '--log2' computes the log2-transformed ratio of
  IP signal to input signal, which is useful for evaluating relative enrichment
  across genomic regions. The script supports parallel processing via SLURM or
  GNU Parallel, or it run submit jobs in serial.

Arguments:
   -h, --help     Print this help message and exit.
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Run script in 'dry-run mode' (optional).
  -ip, --fil_ip   Comma-separated string of ChIP-seq IP BEDGRAPH files.
  -in, --fil_in   Comma-separated string of ChIP-seq input BEDGRAPH files.
  -do, --dir_out  Directory to write coverage signal ratio track outfiles.
  -to, --typ_out  Format of coverage signal ratio track outfiles: 'bedgraph',
                  'bedgraph.gz', 'bdg', 'bdg.gz', 'bg', or 'bg.gz' (default:
                  '${typ_out}').
                    - 'bedgraph', 'bdg', 'bg': BEDGRAPH format (coverage).
                    - 'bedgraph.gz', 'bdg.gz', 'bg.gz': gzip-compressed
                       BEDGRAPH format.
  -sf, --scl_fct  Comma-separated string of scaling factors.
  -dm, --dep_min  Comma-separated string of minimum input depth values used to
                  avoid extreme/erroneous divisions.
  -l2, --log2     Compute log2(sig_ip/sig_in) with a pseudocount of 1 to
                  prevent log(0) errors (optional). The scaling factor
                  ('--scl_fct') is applied after the log2 transformation:
                  'scl_fct * log2(sig_ip/sig_in + 1)'.
   -r, --rnd      Number of decimal places for rounding binned signal (coverage
                  score) ratio values (default: ${rnd}).
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (default: "\${dir_out}/err_out").
  -nj, --nam_job  Prefix for the names of jobs (default: '${nam_job}').
  -mj, --max_job  Maximum number of jobs to run concurrently (default: ${max_job}).
                    - Required if '--slurm' is specified; controls SLURM array
                      tasks.
                    - If '--slurm' is not specified:
                      + If 'max_job' is greater than 1, jobs run in parallel
                        via GNU Parallel.
                      + If 'max_job' is 1, jobs run sequentially (serial mode).
                  - Default: ${max_job}.
  -sl, --slurm    Submit jobs to the SLURM scheduler (optional; otherwise, if
                  'max_job' > 1, run them via GNU Parallel, and if 'max_job' is
                  1, run them in serial.)
  -tm, --time     The length of time, in 'h:mm:ss' format, for the SLURM job
                  (required if '--slurm' is specified, ignored if not; default:
                  '${time}').

Dependencies:
  - Programs
    + awk
    + Bash or Zsh
    + GNU Parallel (if not using SLURM and '--max_job' > 1)
    + Python
    + SLURM (if using '--slurm')
  - Functions
    + check_array_files
    + check_arrays_lengths
    + check_exists_file_dir
    + check_flt_pos
    + check_format_time
    + check_int_pos
    + check_program_path
    + check_str_delim
    + check_supplied_arg
    + debug_array_contents
    + echo_error
    + echo_warning
    + exit_0
    + exit_1
    + handle_env
    + populate_array_empty
    + reset_max_job

Notes:
  - When the '--slurm' flag is used, jobs are parallelized via SLURM array
    tasks; otherwise, jobs are parallelized with GNU Parallel.
  - Outfile names are derived from BEDGRAPH IP infiles and the value(s)
    associated with '--typ_out'.
  - #TODO

Example:
  bash execute_compute_coverage_ratio.sh
      --fil_ip "/path/to/fil_ip_1.bdg.gz,/path/to/fil_ip_2.bdg.gz"
      --fil_in "/path/to/fil_in_1.bdg.gz,/path/to/fil_in_2.bdg.gz"
      --dir_out "/path/to/write/output/files"
      --typ_out "bdg.gz"
      --scl_fct "1.0,1.0"
      --dep_min "0.0035,0.0041"
      --err_out "/path/to/write/output/files/logs"
      --nam_job "compute_ratio"
      --slurm
EOM
)

#  Parse arguments
if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    exit 0
fi

if ${interactive}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -v|--verbose) verbose=true;   shift 1 ;;
            -dr|--dry_run) dry_run=true;   shift 1 ;;
            -fp|--fil_ip)  fil_ip="${2}";  shift 2 ;;
            -fn|--fil_in)  fil_in="${2}";  shift 2 ;;
            -do|--dir_out) dir_out="${2}"; shift 2 ;;
            -to|--typ_out) typ_out="${2}"; shift 2 ;;
            -sf|--scl_fct) scl_fct="${2}"; shift 2 ;;
            -dm|--dep_min) dep_min="${2}"; shift 2 ;;
             -l|--log2)    log2=true;      shift 1 ;;
             -r|--rnd)     rnd="${2}";     shift 2 ;;
            -eo|--err_out) err_out="${2}"; shift 2 ;;
            -nj|--nam_job) nam_job="${2}"; shift 2 ;;
            -sl|--slurm)   slurm=true;     shift 1 ;;
            -mj|--max_job) max_job="${2}"; shift 2 ;;
            -tm|--time)    time="${2}";    shift 2 ;;
            *)
                echo "## Unknown parameter: ${1} ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                exit 1
                ;;
        esac
    done
fi

#  Check arguments ------------------------------------------------------------
check_supplied_arg -a "${env_nam}" -n "env_nam"

check_supplied_arg -a "${scr_sub}" -n "scr_sub"
check_exists_file_dir "f" "${scr_sub}" "scr_sub"

check_supplied_arg -a "${scr_cvg}" -n "scr_cvg"
check_exists_file_dir "f" "${scr_cvg}" "scr_cvg"

check_supplied_arg -a "${denom}" -n "denom"
check_int_pos "${denom}" "denom"

check_supplied_arg -a "${fil_ip}" -n "fil_ip"
check_exists_file_dir "d" "$(dirname "${fil_ip%%[,;]*}")" "fil_ip"
check_str_delim "fil_ip" "${fil_ip}"

check_supplied_arg -a "${fil_in}" -n "fil_in"
check_exists_file_dir "d" "$(dirname "${fil_in%%[,;]*}")" "fil_in"
check_str_delim "fil_in" "${fil_in}"

check_supplied_arg -a "${dir_out}" -n "dir_out"
check_exists_file_dir "d" "${dir_out}" "dir_out"

check_supplied_arg -a "${typ_out}" -n "typ_out"
case "${typ_out}" in
    bedgraph|bedgraph.gz|bdg|bdg.gz|bg|bg.gz) : ;;  # Valid options
    *)
        echo_error \
            "Selection associated with '--typ_out' is not valid:" \
            "'${typ_out}'. Expected 'bedgraph', 'bedgraph.gz', 'bdg'," \
            "'bdg.gz', 'bg', or 'bg.gz'."
        exit_1
        ;;
esac

if [[ -n "${scl_fct}" ]]; then check_str_delim "scl_fct" "${scl_fct}"; fi

check_supplied_arg -a "${dep_min}" -n "dep_min"
check_str_delim "dep_min" "${dep_min}"

check_supplied_arg -a "${rnd}" -n "rnd"
check_int_pos "${rnd}" "rnd"

if [[ -n "${err_out}" ]]; then
    check_exists_file_dir "d" "${err_out}" "err_out"
elif [[ -z "${err_out}" ]]; then
    err_out="${dir_out}/err_out"
    check_exists_file_dir "d" "${err_out}" "err_out"
fi

check_supplied_arg -a "${nam_job}" -n "nam_job"

check_supplied_arg -a "${max_job}" -n "max_job"
check_int_pos "${max_job}" "max_job"

if ${slurm}; then    
    check_supplied_arg -a "${time}" -n "time"
    check_format_time "${time}"
else
    par_job=${max_job}

    check_supplied_arg -a "${par_job}" -n "par_job"
    check_int_pos "${par_job}" "par_job"
fi


#  Activate environment and check that dependencies are in PATH ---------------
handle_env "${env_nam}" > /dev/null

if ! ${slurm} && [[ ${par_job} -gt 1 ]]; then check_program_path parallel; fi
check_program_path python
if ${slurm}; then check_program_path sbatch; fi


#  Parse and validate vector elements -----------------------------------------
IFS=',' read -r -a arr_fil_ip <<< "${fil_ip}"
check_array_files "fil_ip" "${arr_fil_ip[@]}"

IFS=',' read -r -a arr_fil_in <<< "${fil_in}"
check_array_files "fil_in" "${arr_fil_in[@]}"

unset arr_fil_out && typeset -a arr_fil_out
for i in "${arr_fil_ip[@]}"; do
    base=$(basename "${i}")

    base="${base%.bedgraph.gz}"
    base="${base%.bdg.gz}"
    base="${base%.bg.gz}"
    base="${base%.bedgraph}"
    base="${base%.bdg}"
    base="${base%.bg}"

    arr_fil_out+=( "${dir_out}/${base}.${typ_out}" )
done
unset base

if [[ -z "${scl_fct}" ]]; then
    unset arr_scl_fct && typeset -a arr_scl_fct
    populate_array_empty arr_scl_fct "${#arr_fil_ip[@]}"
else
    IFS=',' read -r -a arr_scl_fct <<< "${scl_fct}"
fi

for s in "${arr_scl_fct[@]}"; do
    if [[ "${s}" != "#N/A" ]]; then check_flt_pos "${s}" "scl_fct"; fi
done

#TODO: Need processing for dep_min, which is an array...



#  Submit job via SLURM or GNU Parallel
if ${slurm}; then
    sbatch --job-name="${nam_job}" --array=1-${max_job}%${max_job} "${scr_sub}" \
        "${env_nam}" "${scr_cvg}" "${fil_ip}" "${fil_in}" "${fil_out}" \
        "${scl_fct}" "${dep_min}" "${log2}" "${rnd}" "${err_out}" "${nam_job}"
else
    bash "${scr_sub}" "${env_nam}" "${scr_cvg}" "${fil_ip}" "${fil_in}" \
        "${fil_out}" "${scl_fct}" "${dep_min}" "${log2}" "${rnd}" "${err_out}" "${nam_job}"
fi
















