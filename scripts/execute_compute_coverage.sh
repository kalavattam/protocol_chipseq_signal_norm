#!/bin/bash

#  execute_compute_coverage.sh
#  KA


#  Run script in interactive/test mode (true) or command-line mode (false)
interactive=true

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
    source "${dir_fnc}/check_flt_pos.sh"
    source "${dir_fnc}/check_format_time.sh"
    source "${dir_fnc}/check_int_pos.sh"
    source "${dir_fnc}/check_mut_excl_args.sh"
    source "${dir_fnc}/check_mut_excl_flags.sh"
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/check_str_delim.sh"
    source "${dir_fnc}/check_supplied_arg.sh"
    source "${dir_fnc}/echo_error.sh"
    source "${dir_fnc}/echo_warning.sh"
    source "${dir_fnc}/exit_0.sh"
    source "${dir_fnc}/exit_1.sh"
    source "${dir_fnc}/handle_env.sh"
    source "${dir_fnc}/reset_max_job.sh"
}

#  Helper function to validate table is not empty
function check_table() {
    local table="${1}"

    if [[ $(wc -l < "${table}") -le 1 ]]; then
        echo \
            "Error: Table file '${table}' is empty or contains only a" \
            "header." >&2
        return 1
    fi
}


#  Helper function to vaildate existence of column in table
function check_table_column() {
    local table="${1}"
    local column="${2}"

    if \
        ! awk -F '\t' -v col="${column}" '
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    if ($i == col) exit 0
                } exit 1
            }
        ' "${table}"
    then
        echo \
            "Error: Column '${column}' not found in table header:" \
            "$(awk 'NR==1' "${table}" | tr '\t' ' ')"
        return 1
    fi
}


#  Helper function to warn that user-supplied scaling factors override those in
#+ table
function check_scaling_factor_table() {
    local scl_fct="${1}"
    local table="${2}"

    if [[ -n "${scl_fct}" && -n "${table}" ]]; then
        echo \
            "Warning: --scl_fct will override scaling factors from" \
            "--table." >&2
        return 1
    fi
}


#  Helper function to warn that normalization flags override scaling factors in
#+ table
function check_normalization_table() {
    local name="${1}"
    local norm="${2}"
    local table="${3}"

    if ${norm} && [[ -n "${table}" ]]; then
        echo \
            "Warning: --${name} will override scaling factors from" \
            "--table." >&2
        return 1
    fi
}


#  Helper function to populate array with '#N/A' if it is empty
function populate_array_empty() {
    local -n arr="${1}"
    local target="${2}"

    if [[ -z "${arr[*]}" ]]; then
        for ((i = 0; i < target; i++)); do arr+=( "#N/A" ); done
    fi
}


#  Helper function to validate that two arrays have matching lengths
function validate_length_arrays() {
    local arr_nam_1="${1}"
    local -n  arr_1="${2}"
    local arr_nam_2="${3}"
    local -n  arr_2="${4}"

    if [[ "${#arr_1[@]}" -ne "${#arr_2[@]}" ]]; then
        echo \
            "Error: '${arr_nam_1}' must match the number of '${arr_nam_2}'." \
            "Got ${#arr_1[@]} for '${arr_nam_1}' but ${#arr_2[@]} for" \
            "'${arr_nam_2}'." >&2
        return 1
    fi
}


#  Helper function to validate existence of files in arrays
function validate_files_array() {
    local desc="${1}"
    shift

    #  Check that files are supplied
    if [[ "$#" -eq 0 ]]; then
        echo "Error: No files supplied to validate for ${desc}." >&2
        return 1
    fi
    
    for file in "$@"; do
        check_exists_file_dir "f" "${file}" "${desc}"
    done
}


#  Helper function to debug array contents
function debug_contents_array() {
    for arr_nam in "$@"; do
        #  Access the array indirectly using eval
        eval "arr=( \"\${${arr_nam}[@]}\" )"
        if [[ -n "${arr[*]}" ]]; then
            echo "  - ${arr_nam}: ( ${arr[*]} )"
        fi
    done
}


#  Set up paths, values, and parameters for interactive mode
function set_interactive() {
    #  Set hardcoded paths, values, etc.
    ## WARNING: Change values if you're not Kris and `interactive=true` ##
    dir_bas="${HOME}/tsukiyamalab/Kris"
    dir_rep="${dir_bas}/202X_protocol_ChIP"
    dir_scr="${dir_rep}/scripts"
    dir_dat="${dir_rep}/data"
    dir_pro="${dir_dat}/processed"
    {
        aligner="bowtie2"
        a_type="global"
        req_flg=true
        flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
        mapq=1
        det_bam="flag-${flg}_mapq-${mapq}"
        det_cov="${aligner}_${a_type}_${det_bam}"
    }
    dir_aln="${dir_pro}/align_${aligner}_${a_type}"
    dir_bam="${dir_aln}/${det_bam}/sc"
    dir_cov="${dir_pro}/compute_coverage"
    # dir_tbl="${dir_cov}/${det_cov}/alpha/tables"  # "${dir_cov}/${det_cov}/spike/tables"
    dir_trk="${dir_cov}/${det_cov}/alpha/tracks"  # "${dir_cov}/${det_cov}/spike/tracks"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=true
    threads=8
    # infiles=""
    infiles="$(  ## WARNING: Change the search parameters as needed ##
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_bam}" \
            --pattern "*.bam" \
            --include "IP*,*Hho1*"  # "IP*,*Hmo1*"  # "IP*,*Brn1*"
    )"
    # table="${dir_tbl}/IP_WT_log-Q_Brn1_rep1-rep2-rep3.tsv"  # "${dir_tbl}/IP_WT_G1-G2M-Q_Hho1_6336-6337.tsv"  # "${dir_tbl}/IP_WT_G1-G2M-Q_Hmo1_7750-7751.tsv"  # ""  # "IP_WT_G1-G2M-Q_Hho1-Hmo1_6336-6337_7750-7751.tsv"
    tbl_col="alpha"  # "scaled"  # ""
    dir_out="${dir_trk}"
    typ_out="bigwig"
    bin_siz=1
    scl_fct="0.002054,0.003138,0.003127,0.003522,0.056611,0.02906"  # ""
    norm=false
    raw=false
    usr_frg=""
    err_out="${dir_trk}/logs"
    nam_job="compute_coverage"
    slurm=true
    max_job=4  # 6
    time="0:30:00"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_analyze"
scr_slm="${dir_scr}/submit_compute_coverage_slurm.sh"
scr_par="${dir_scr}/submit_compute_coverage_parallel.sh"
scr_cvg="${dir_scr}/compute_coverage.py"
denom=4

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=1
infiles=""
table=""
tbl_col="alpha"
dir_out=""
typ_out="bigwig"
bin_siz=10
scl_fct=""
norm=false
raw=false
usr_frg=""
err_out=""
nam_job="compute_coverage"
slurm=false
max_job=6
time="0:30:00"

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  execute_compute_coverage.sh
    [--verbose] [--dry_run] --threads <int> --infiles <str> [--table <str>]
    --dir_out <str> --typ_out <str> --bin_siz <int> [--scl_fct <flt>] [--norm]
    [--raw] [--usr_frg <flt>] --err_out <str> --nam_job <str> [--slurm]
    [--max_job <int>] [--time <str>]

Description:
  execute_compute_coverage.sh automates the calculation of coverage signal
  tracks from BAM files. It supports various normalization methods, including
  user-supplied scaling factors (e.g., unadjusted or relativized spike-in
  scaling, or siQ-ChIP scaling), normalized (fractional) coverage (per Dickson
  et al., Sci Rep 2023), and raw (unadjusted) coverage. The script integrates
  with GNU Parallel or SLURM for parallel processing.

Arguments:
   -h, --help     Print this help message and exit.
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Run script in 'dry-run mode' (optional).
   -t, --threads  Number of threads to use (default: ${threads}).
   -i, --infiles  Comma-separated string vector of coordinate-sorted BAM
                  infiles. Cannot be used with --table.
  -tb, --table    Path to TSV infile listing BAM files with calculated spike-in
                  (relativized or not) or alpha scaling factors. Columns must
                  include 'sample' and, depending on the scaling factor type,
                  'sf', 'scaled', or 'alpha'. Cannot be used with --infiles.
  -tc, --tbl_col  Column in --table containing scaling factors: 'sf', 'scaled',
                  or 'alpha' (required if --table is invoked, ignored if not;
                  default: ${tbl_col}).
  -do, --dir_out  The directory to write coverage signal track outfiles.
  -to, --typ_out  Format of coverage signal track outfiles: 'bedgraph',
                  'bigwig', or 'both' (default: ${typ_out}).
  -bs, --bin_siz  Bin size for coverage calculation in base pairs (default:
                  ${bin_siz}).
  -sf, --scl_fct  Comma-separated string vector of scaling factors to apply to
                  coverage (optional). Must match the number of infiles via
                  --infiles or --table. If used with --table, these values
                  override those in the specified column (via --tbl_col) with a
                  warning. Cannot be used with --norm or --raw. 
  -no, --norm     Normalize coverage by fragment length and total alignments,
                  generating 'normalized coverage' per Dickson et al., Sci Rep
                  2023 (optional). If used with --table, these values override
                  those in the specified column (via --tbl_col) with a warning.
                  Cannot be used with --scl_fct or --raw.
   -r, --raw      Calculate raw, unadjusted coverage using BAM infiles as-is,
                  without applying any normalization or scaling factors
                  (optional). This is equivalent to specifying --infiles
                  without any additional coverage modification options. If used
                  with --table, these values override those in the specified
                  column (via --tbl_col) with a warning. Cannot be used with
                  --scl_fct or --norm.
  -uf, --usr_frg  Comma-separated string vector of fragment lengths to use
                  instead of read lengths (single-end alignments) or template
                  lengths (paired-end alignments; optional). Must match the
                  number of infiles via --infiles or --table.
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (default: \${dir_out}/err_out).
  -nj, --nam_job  Prefix for the names of jobs (default: ${nam_job}).
  -sl, --slurm    Submit jobs to the SLURM scheduler; otherwise, run them via
                  GNU Parallel (optional).
  -mj, --max_job  The maximum number of jobs to run at one time (required if
                  --slurm is specified, ignored if not; default: ${max_job}).
  -tm, --time     The length of time, in 'h:mm:ss' format, for the SLURM job
                  (required if '--slurm' is specified, ignored if not; default:
                  ${time}).

Dependencies:
  - Programs
    + awk
    + Bash or Zsh
    + GNU Parallel (if not using SLURM and threads > 1)
    + Python
    + SLURM (if using --slurm)
  - Functions
    + check_exists_file_dir
    + check_flt_pos
    + check_format_time
    + check_int_pos
    + check_mut_excl_args
    + check_mut_excl_flags
    + check_normalization_table
    + check_program_path
    + check_scaling_factor_table
    + check_str_delim
    + check_supplied_arg
    + check_table
    + check_table_column
    + debug_contents_array
    + echo_error
    + echo_warning
    + exit_0
    + exit_1
    + handle_env
    + populate_array_empty
    + reset_max_job
    + validate_files_array
    + validate_length_arrays

Notes:
  - When the '--slurm' flag is used, jobs are parallelized via SLURM array
    tasks; otherwise, jobs are parallelized with GNU Parallel.
  - BAM infiles must be coordinate-sorted.
  - Outfile names are derived BAM infiles and the value(s) associated with
    --typ_out.
  - #TODO

Examples:
  \`\`\`
  #  Example 1: Run with GNU Parallel and --infiles
  bash "\${HOME}/path/to/scripts/execute_compute_coverage.sh" \
      --threads 8
      --infiles "\${HOME}/path/to/sample1.bam,\${HOME}/path/to/sample2.bam"
      --dir_out "\${HOME}/path/to/outfiles"
      --typ_out "bigwig"
      --bin_siz 50
      --norm
      --err_out "\${HOME}/path/to/outfiles/logs"
      --nam_job "coverage_calculation"

  #  Example 2: Run with SLURM and --table
  bash "\${HOME}/path/to/scripts/execute_compute_coverage.sh" \
      --threads 8
      --table "\${HOME}/path/to/table.tsv"
      --tbl_col "scaled"
      --dir_out "\${HOME}/path/to/outfiles"
      --typ_out "both"
      --bin_siz 1
      --slurm
      --max_job 10
      --time "0:20:00"
      --err_out "\${HOME}/path/to/outfiles/logs"
      --nam_job "slurm_coverage"
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
             -v|--verbose) verbose=true;     shift 1 ;;
            -dr|--dry_run) dry_run=true;     shift 1 ;;
             -t|--threads) threads="${2}";   shift 2 ;;
             -i|--infiles) infiles="${2}";   shift 2 ;;
            -tb|--table)   table="${2}";     shift 2 ;;
            -tc|--tbl_col) tbl_col="${2,,}"; shift 2 ;;
            -do|--dir_out) dir_out="${2}";   shift 2 ;;
            -to|--typ_out) typ_out="${2,,}"; shift 2 ;;
            -bs|--bin_siz) bin_siz="${2}";   shift 2 ;;
            -sf|--scl_fct) scl_fct="${2}";   shift 2 ;;
            -no|--norm)    norm=true;        shift 1 ;;
             -r|--raw)     raw=true;         shift 1 ;;
            -uf|--usr_frg) usr_frg="${2}";   shift 2 ;;
            -eo|--err_out) err_out="${2}";   shift 2 ;;
            -nj|--nam_job) nam_job="${2}";   shift 2 ;;
            -sl|--slurm)   slurm=true;       shift 1 ;;
            -mj|--max_job) max_job="${2}";   shift 2 ;;
            -tm|--time)    time="${2}";      shift 2 ;;
            *)
                echo "## Unknown parameter passed: ${1} ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                exit_1
                ;;
        esac
    done
fi


#  Check arguments ------------------------------------------------------------
#  Check that environment name is supplied
check_supplied_arg -a "${env_nam}" -n "env_nam"

#  Check that submission scripts are supplied and exists as files
check_supplied_arg -a "${scr_slm}" -n "scr_slm"
check_exists_file_dir "f" "${scr_slm}" "scr_slm"

check_supplied_arg -a "${scr_par}" -n "scr_par"
check_exists_file_dir "f" "${scr_par}" "scr_par"

#  Check that coverage script is supplied and exists as a file
check_supplied_arg -a "${scr_cvg}" -n "scr_cvg"
check_exists_file_dir "f" "${scr_cvg}" "scr_cvg"

#  Check that denominator is supplied
check_supplied_arg -a "${denom}" -n "denom"

#  Check that threads value is supplied and is a positive integer
check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

#  Check that --infiles and --table are mutually exclusive; one must be
#+ specified
check_mut_excl_args "infiles" "${infiles}" "table" "${table}"

#  Validate --infiles or --table
if [[ -n "${infiles}" ]]; then
    #  Validate existence of directory containing infiles 
    check_exists_file_dir "d" "$(dirname "${infiles%%[,;]*}")" "infiles"

    #  Ensure --infiles is not empty or improperly formatted
    check_str_delim "infiles" "${infiles}"
elif [[ -n "${table}" ]]; then
    #  Validate existence of table file
    check_exists_file_dir "f" "${table}" "table"

    #  Validate table isn't empty by checking its number of lines
    check_table "${table}"

    #  When --table is specified, validate --tbl_col
    check_supplied_arg -a "${tbl_col}" -n "tbl_col"
    case "${tbl_col}" in
        alpha|scaled|sf) : ;;  # Valid options
        *)
            echo_error \
                "Selection associated with --tbl_col is not valid:" \
                "${tbl_col}. Selection must be 'alpha', 'scaled', or 'sf'."
            exit_1
            ;;
    esac

    #  If present, trim leading/trailing whitespace from tbl_col
    tbl_col=$(echo "${tbl_col}" | xargs)

    #  Validate tbl_col exists in the table header
    check_table_column "${table}" "${tbl_col}"
fi

#  Check that output directory is supplied and exists
check_supplied_arg -a "${dir_out}" -n "dir_out"
check_exists_file_dir "d" "${dir_out}" "dir_out"

#  Check that output type is supplied and valid
check_supplied_arg -a "${typ_out}" -n "typ_out"
case "${typ_out}" in
    bedgraph|bigwig|both) : ;;  # Valid options
    *)
        echo_error \
            "Selection associated with --typ_out is not valid: ${typ_out}." \
            "Selection must be 'bedgraph', 'bigwig', or 'both'."
        exit_1
        ;;
esac

#  Check that bin size is supplied and is a positive integer
check_supplied_arg -a "${bin_siz}" -n "bin_siz"
check_int_pos "${bin_siz}" "bin_siz"

#  Ensure --scl_fct and --table compatibility
check_scaling_factor_table "${scl_fct}" "${table}"

#  If applicable, validate the mutual exclusivity of --scl_fct with --norm or
#+ --raw, and check that --scl_fct is properly formatted
if [[ -n "${scl_fct}" ]]; then
    if ${norm}; then
        check_mut_excl_args "scl_fct" "${scl_fct}" "norm" "${norm}"
    elif ${raw:-false}; then
        check_mut_excl_args "scl_fct" "${scl_fct}" "raw" "${raw:-false}"
    fi

    check_str_delim "scl_fct" "${scl_fct}"
fi

#  Ensure mutual exclusivity: --norm and --raw
check_mut_excl_flags "${norm}" "norm" "${raw:-false}" "raw"

#  Warn about --norm or --raw overriding --table
if [[ -n "${table}" ]]; then
    check_normalization_table "norm" ${norm} "${table}"
    check_normalization_table "raw" ${raw:-false} "${table}"
fi

#  Default to --raw if no scaling option is provided
if [[ -z "${scl_fct}" && -z "${table}" ]] && ! ${norm} && ! ${raw:-false}; then
    echo_warning \
        "No scaling factor or normalization option specified. Defaulting to" \
        "non-normalized/unadjusted coverage (--raw)."
    # raw=true
fi

#  If applicable, ensure --usr_frg is not empty or improperly formatted
if [[ -n "${usr_frg}" ]]; then check_str_delim "usr_frg" "${usr_frg}"; fi

#  Check that stderr/stdout directory is supplied and exists, or assign default
if [[ -n "${err_out}" ]]; then
    check_exists_file_dir "d" "${err_out}" "err_out"
elif [[ -z "${err_out}" ]]; then
    err_out="${dir_out}/err_out"
    check_exists_file_dir "d" "${err_out}" "err_out"
fi

#  Check that job name is supplied
check_supplied_arg -a "${nam_job}" -n "nam_job"

#  If submitting jobs to SLURM, perform checks for maximum jobs allowed at a
#  time (max_job) and the maximum permitted runtime (time)
if ${slurm}; then
    #  Ensure max_job is provided and is a positive integer
    check_supplied_arg -a "${max_job}" -n "max_job"
    check_int_pos "${max_job}" "max_job"
    
    #  Ensure time is provided and has a valid time format
    check_supplied_arg -a "${time}" -n "time"
    check_format_time "${time}"
else
    #  For non-SLURM jobs, calculate the number of parallel jobs to run
    #  based on the number of threads, setting the thread count accordingly
    par_job=$(( threads / denom ))
    threads=${denom}

    #  Ensure the calculated number of parallel jobs (par_job) is provided and
    #  is a positive integer
    check_supplied_arg -a "${par_job}" -n "par_job"
    check_int_pos "${par_job}" "par_job"

    #  Ensure the number of threads is provided and is a positive integer
    check_supplied_arg -a "${threads}" -n "threads"
    check_int_pos "${threads}" "threads"
fi


#  Debug summary output of resolved argument states ---------------------------
#+ ...i.e., for --table, --tbl_col, --scl_fct, --norm, and --raw
if ${verbose}; then
    #  Determine the normalization method
    if ${norm}; then
        mth_nrm="Normalized coverage (Dickson et al., Sci Rep 2023)"
        if [[ -n "${table}" && -n "${tbl_col}" ]]; then
            src_scl="--table (--tbl_col: '${tbl_col}') overridden by --norm"
        else
            src_scl="--norm"
        fi
    elif ${raw:-false}; then
        mth_nrm="Raw (unadjusted) coverage"
        if [[ -n "${table}" && -n "${tbl_col}" ]]; then
            src_scl="--table (--tbl_col: '${tbl_col}') overridden by --raw"
        else
            src_scl="--raw"
        fi
    else
        #  Default to --raw if no scaling or normalization method is specified
        if [[ -n "${scl_fct}" ]]; then
            src_scl="--scl_fct"
            mth_nrm="User-supplied scaling via --scl_fct"
        elif [[ -n "${table}" && -n "${tbl_col}" ]]; then
            src_scl="--table via --tbl_col '${tbl_col}'"
            mth_nrm="User-supplied scaling via --table and --tbl_col"
        else
            src_scl="None"
            mth_nrm="None; default to raw (unadjusted) coverage"
        fi
    fi

    echo "##############################"
    echo "## Resolved argument states ##"
    echo "##############################"
    echo ""
    echo "- Scaling factor source: ${src_scl}"
    echo "- Normalization method: ${mth_nrm}"
    echo ""
    echo ""
fi


#  Activate environment and check that dependencies are in PATH ---------------
handle_env "${env_nam}" > /dev/null

check_program_path awk
if ! ${slurm}; then check_program_path parallel; fi
check_program_path python
if ${slurm}; then check_program_path sbatch; fi


#  Parse and validate table and/or vector elements ----------------------------
#  If applicable, parse values associated with --table
if [[ -n "${table}" ]]; then
    #  Read the table's header to determine available columns
    header=$(awk 'NR == 1' "${table}")
    IFS=$'\t' read -r -a arr_header <<< "${header}"

    #  Determine indices of relevant columns
    idx_smp=-1
    idx_sf=-1
    idx_scl=-1
    idx_alf=-1
    for i in "${!arr_header[@]}"; do
        case "${arr_header[i]}" in
            sample) idx_smp=${i} ;;
            sf)     idx_sf=${i}  ;;
            scaled) idx_scl=${i} ;;
            alpha)  idx_alf=${i} ;;
        esac
    done

    #  Initialize arrays to store values
    unset arr_infiles && typeset -a arr_infiles
    unset arr_sf      && typeset -a arr_sf
    unset arr_scaled  && typeset -a arr_scaled
    unset arr_alpha   && typeset -a arr_alpha

    #  Parse the table rows, skipping the header
    while IFS=$'\t' read -r row; do
        IFS=$'\t' read -r -a fields <<< "${row}"

        #  Extract the 'sample' column (always present)
        arr_infiles+=( "${fields[idx_smp]}" )

        #  Extract scaling factor columns if present
        [[ ${idx_sf}  -ne -1 ]] && arr_sf+=( "${fields[idx_sf]}" )
        [[ ${idx_scl} -ne -1 ]] && arr_scaled+=( "${fields[idx_scl]}" )
        [[ ${idx_alf} -ne -1 ]] && arr_alpha+=( "${fields[idx_alf]}" )
    done < <(awk 'NR > 1' "${table}")

    #  Dynamically assign arr_scl_fct based on tbl_col
    unset arr_scl_fct && typeset -a arr_scl_fct
    if [[ -z "${scl_fct}" ]] && ! ${norm} && ! ${raw:-false}; then
        eval "arr_scl_fct=( \"\${arr_${tbl_col}[@]}\" )"
    fi
fi

#  If applicable, parse --infiles
if [[ -n "${infiles}" ]]; then
    IFS=',' read -r -a arr_infiles <<< "${infiles}"
fi

#  Having parsed --table or --infiles, validate each infile exists
validate_files_array "infiles/table" "${arr_infiles[@]}"

#  Based on arr_infiles and dir_out, construct outfile paths, using them to
#+ populate arr_outfiles
unset arr_outfiles && typeset -a arr_outfiles
for i in "${arr_infiles[@]}"; do
    if [[ "${tbl_col}" == "scaled" || "${tbl_col}" == "sf" ]]; then
        arr_outfiles+=( "${dir_out}/$(basename "${i}" .bam).${tbl_col}" )
    else
        arr_outfiles+=( "${dir_out}/$(basename "${i}" .bam)" )
    fi
done

#  If applicable, parse --scl_fct
if [[ -n "${scl_fct}" ]]; then
    IFS=',' read -r -a arr_scl_fct <<< "${scl_fct}"
fi

#  If empty, populate --scl_fct array with "#N/A"
populate_array_empty arr_scl_fct "${#arr_infiles[@]}"

#  Having parsed --table or --scl_fct, validate each scaling factor as a
#+ positive float (if applicable)
if [[ "${arr_scl_fct[0]}" != "#N/A" ]]; then
    for s in "${arr_scl_fct[@]}"; do check_flt_pos "${s}" "scl_fct"; done
fi

#  If applicable, parse and validate --usr_frg; if empty, populate with "#N/A"
if [[ -n "${usr_frg}" ]]; then
    IFS=',' read -r -a arr_usr_frg <<< "${usr_frg}"
fi

#  If empty, populate --usr_frg array with "#N/A"
populate_array_empty arr_usr_frg "${#arr_infiles[@]}"

#  Having parsed --usr_frg, validate each scaling factor as a positive float
#+ (if applicable)
if [[ "${arr_usr_frg[0]}" != "#N/A" ]]; then
    for u in "${arr_usr_frg[@]}"; do check_flt_pos "${u}" "usr_frg"; done
fi

#  Ensure matching element counts between arr_infiles and the other arrays
validate_length_arrays "outfiles" arr_outfiles "infiles" arr_infiles
validate_length_arrays "scl_fct"  arr_scl_fct  "infiles" arr_infiles
validate_length_arrays "usr_frg"  arr_usr_frg  "infiles" arr_infiles

#  Reset max_job if it is greater than the number of infiles
if ${slurm}; then
    max_job=$(reset_max_job "${max_job}" "${#arr_infiles[@]}")
fi

#  Debug output for infiles and other values
if ${verbose}; then
    echo "## Parsed vectors for parallelization ##"
    debug_contents_array "arr_infiles" "arr_scl_fct" "arr_usr_frg"
    if ${slurm}; then
        echo "  - Max no. jobs to run at a time (SLURM): ${max_job}"
    else
        echo "  - Max no. jobs to run at a time (GNU Parallel): ${par_job}"
    fi
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
    echo "scr_slm=${scr_slm}"
    echo "scr_par=${scr_par}"
    echo "scr_cvg=${scr_cvg}"
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
    echo "infiles=${infiles:-#N/A}"
    echo "table=${table:-#N/A}"
    echo "tbl_col=${tbl_col:-#N/A}"
    echo "dir_out=${dir_out}"
    echo "typ_out=${typ_out}"
    echo "bin_siz=${bin_siz}"
    echo "scl_fct=${scl_fct:-#N/A}"
    echo "norm=${norm}"
    echo "raw=${raw:-false}"
    echo "usr_frg=${usr_frg:-#N/A}"
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
    echo "arr_infiles=( ${arr_infiles[*]} )"
    echo "arr_outfiles=( ${arr_outfiles[*]} )"
    echo "arr_scl_fct=( ${arr_scl_fct[*]:-#N/A} )"
    echo "arr_usr_frg=( ${arr_usr_frg[*]:-#N/A} )"
    echo ""
    echo ""
fi


# shellcheck disable=SC2016
if ${slurm}; then
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
        echo "        ${scr_slm} \\"
        echo "            ${env_nam} \\"
        echo "            ${scr_cvg} \\"
        echo "            ${threads} \\"
        echo "            $(echo "${arr_infiles[*]}" | tr ' ' ',') \\"
        echo "            $(echo "${arr_outfiles[*]}" | tr ' ' ',') \\"
        echo "            ${typ_out} \\"
        echo "            ${bin_siz} \\"
        echo "            ${norm} \\"
        echo "            ${raw:-false} \\"
        echo "            $(echo "${arr_scl_fct[*]}" | tr ' ' ',') \\"
        echo "            $(echo "${arr_usr_frg[*]}" | tr ' ' ',') \\"
        echo "            ${err_out} \\"
        echo "            ${nam_job}"
        echo ""
        echo ""
        echo "#########################################"
        echo "## Contents of SLURM submission script ##"
        echo "#########################################"
        echo ""
        echo "## ${scr_slm} ##"
        echo ""
        cat "${scr_slm}"
        echo ""
        echo ""
        echo "###########################################"
        echo "## Contents of Python computation script ##"
        echo "###########################################"
        echo ""
        echo "## ${scr_cvg} ##"
        echo ""
        cat "${scr_cvg}"
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
            --array=1-${#arr_infiles[@]}%${max_job} \
                ${scr_slm} \
                    ${env_nam} \
                    ${scr_cvg} \
                    ${threads} \
                    $(echo "${arr_infiles[*]}" | tr ' ' ',') \
                    $(echo "${arr_outfiles[*]}" | tr ' ' ',') \
                    ${typ_out} \
                    ${bin_siz} \
                    ${norm} \
                    ${raw:-false} \
                    $(echo "${arr_scl_fct[*]}" | tr ' ' ',') \
                    $(echo "${arr_usr_frg[*]}" | tr ' ' ',') \
                    ${err_out} \
                    ${nam_job}
    fi
else
    #  Create a temporary GNU Parallel configuration file if threads > 1
    if [[ "${threads}" -gt 1 ]]; then
        #  Generate a configuration file
        config="${err_out}/${nam_job}.${RANDOM}.txt"
        touch "${config}" ||
            {
                echo_error \
                    "Failed to create a GNU Parallel configuration file."
                exit_1
            }

        #  Populate the GNU Parallel configuration file with parameters for
        #+ each job
        for i in "${!arr_infiles[@]}"; do
            echo \
                "${env_nam}" \
                "${scr_cvg}" \
                "${threads}" \
                "${arr_infiles[i]}" \
                "${arr_outfiles[i]}" \
                "${typ_out}" \
                "${bin_siz}" \
                ${norm} \
                ${raw:-false} \
                "${arr_scl_fct[i]}" \
                "${arr_usr_frg[i]}" \
                "${err_out}" \
                "${nam_job}" \
                    >> "${config}"
        done

        #  Perform a dry-run to debug and/or verify commands
        if ${dry_run} || ${verbose}; then
            # shellcheck disable=SC2090
            parallel --colsep ' ' --jobs "${par_job}" --dryrun \
                "bash \"${scr_par}\" {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13}" \
                :::: "${config}"
        fi

        #  Use GNU Parallel for multi-threaded execution without SLURM
        if ! ${dry_run}; then
            # shellcheck disable=SC2090
            parallel --colsep ' ' --jobs "${par_job}" \
                "bash \"${scr_par}\" {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13}" \
                :::: "${config}"
        fi
    else
        if ${dry_run} || ${verbose}; then
            echo "###########################################"
            echo "## Serial call(s) for computing coverage ##"
            echo "###########################################"
            echo ""

            iter=0
            for i in "${!arr_infiles[@]}"; do
                (( iter++ )) || true
                pth_std="${err_out}/${nam_job}.$(basename "${arr_infiles[i]}").${iter}"

                echo "## Call no. ${iter} ##"
                echo "bash ${scr_par} \\"
                echo "    ${env_nam} \\"
                echo "    ${scr_cvg} \\"
                echo "    ${threads} \\"
                echo "    ${arr_infiles[i]} \\"
                echo "    ${arr_outfiles[i]} \\"
                echo "    ${typ_out} \\"
                echo "    ${bin_siz} \\"
                echo "    ${norm} \\"
                echo "    ${raw:-false} \\"
                echo "    ${arr_scl_fct[i]} \\"
                echo "    ${arr_usr_frg[i]} \\"
                echo "    ${err_out} \\"
                echo "    ${nam_job} \\"
                echo "         > ${pth_std}.stdout.txt \\"
                echo "        2> ${pth_std}.stderr.txt"
                echo ""
            done
            echo ""
        fi

        if ! ${dry_run}; then
            for i in "${!arr_infiles[@]}"; do
                pth_std="${err_out}/${nam_job}.$(basename "${arr_infiles[i]}").${iter}"
                
                bash "${scr_par}" \
                    "${env_nam}" \
                    "${scr_cvg}" \
                    "${threads}" \
                    "${arr_infiles[i]}" \
                    "${arr_outfiles[i]}" \
                    "${typ_out}" \
                    "${bin_siz}" \
                    ${norm} \
                    ${raw:-false} \
                    "${arr_scl_fct[i]}" \
                    "${arr_usr_frg[i]}" \
                    "${err_out}" \
                    "${nam_job}" \
                         > "${pth_std}.stdout.txt" \
                        2> "${pth_std}.stderr.txt"
            done
        fi
    fi
fi
