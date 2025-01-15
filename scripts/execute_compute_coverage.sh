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
    source "${dir_fnc}/check_mut_excl_args.sh"
    source "${dir_fnc}/check_mut_excl_flags.sh"
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/check_str_delim.sh"
    source "${dir_fnc}/check_supplied_arg.sh"
    source "${dir_fnc}/check_table.sh"
    source "${dir_fnc}/check_table_column.sh"
    source "${dir_fnc}/check_table_scaling_factor.sh"
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

    dir_aln="${dir_pro}/align_${aligner}_${a_type}"
    dir_bam="${dir_aln}/${det_bam}/sc"
    dir_cov="${dir_pro}/compute_coverage"
    # dir_tbl="${dir_cov}/${det_cov}/alpha/tables"
    typ_trk="norm"
    dir_trk="${dir_cov}/${det_cov}/${typ_trk}/tracks"

    pattern="*.bam"
    include=""  # "IP*"
    exclude="*Brn1*"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=true
    threads=8
    # infiles=""
    infiles="$(  ## WARNING: Change search parameters as needed ##
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_bam}" \
            --pattern "${pattern}" \
            --include "${include}" \
            --exclude "${exclude}"
    )"
    # table=""
    # table="IP_WT_G1-G2M-Q_Hho1-Hmo1_6336-6337_7750-7751.tsv"
    tbl_col="none"  # "alpha"  # "scaled"
    dir_out="${dir_trk}"
    typ_out="bdg"
    siz_bin=10  # 30  # 1
    scl_fct="none"  # "0.002054,0.003138,0.003127,0.003522,0.056611,0.02906"
    typ_cvg="${typ_trk}"
    usr_frg="none"
    err_out="${dir_trk}/logs"
    nam_job="compute_coverage"
    slurm=true
    max_job=6  # 4
    time="0:30:00"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_analyze"
scr_sub="${dir_scr}/submit_compute_coverage.sh"
scr_cvg="${dir_scr}/compute_coverage.py"
denom=4

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=1
infiles=""
table=""
tbl_col="none"  # "alpha"
dir_out=""
typ_out="bigwig"
siz_bin=10
scl_fct="none"
typ_cvg="norm"
usr_frg="none"
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
    [--tbl_col <str>] --dir_out <str> --typ_out <str> --siz_bin <int>
    [--scl_fct <flt>] --typ_cvg <str> [--usr_frg <flt>] --err_out <str>
    --nam_job <str> [--slurm] [--max_job <int>] [--time <str>]

Description:
  The driver script 'execute_compute_coverage.sh' automates the calculation of
  coverage signal tracks from BAM files. It supports various normalization
  methods, including user-supplied scaling factors, normalized (probability
  distribution) coverage (per Dickson et al., Sci Rep 2023), and raw
  (unadjusted) coverage. The script integrates with GNU Parallel or SLURM for
  parallel processing.

Arguments:
   -h, --help     Print this help message and exit.
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Run script in 'dry-run mode' (optional).
   -t, --threads  Number of threads to use (default: ${threads}).
   -i, --infiles  Comma-separated string vector of coordinate-sorted BAM
                  infiles. Cannot be used with '--table'.
  -tb, --table    Path to TSV infile listing BAM files with or without scaling
                  factors. Columns must include 'sample' and, depending on the
                  scaling factor type, 'sf', 'scaled', or 'alpha'. Cannot be
                  used with '--infiles'.
  -tc, --tbl_col  Specify the column in '--table' that contains scaling
                  factors ('sf', 'scaled', or 'alpha'), or use 'none' to ignore
                  scaling factors in the table (default: '${tbl_col}'). This
                  argument is required if '--table' is specified. If '--scl_fct'
                  is set to any value other than 'none', it takes precedence
                  over '--tbl_col'. Ignored if '--table' is not specified.
  -do, --dir_out  Directory to write coverage signal track outfiles.
  -to, --typ_out  Format of coverage signal track outfiles: 'bedgraph' (or
                  'bdg' or 'bg'), 'bigwig' (or 'bw'), or 'both'. Alternatively,
                  one can set '--typ_out bed' to output alignments processed by
                  function parse_bam() BED-like format (note: these are not
                  coverage values; default: '${typ_out}').
  -sb, --siz_bin  Bin size for coverage calculation in base pairs (default:
                  ${siz_bin}).
  -sf, --scl_fct  Comma-separated string vector of scaling factors to apply to
                  coverage (optional; default: '${scl_fct}'). Must match the
                  number of infiles via '--infiles' or '--table'. If used with
                  '--table', these values override those in the specified
                  column (via '--tbl_col') with a warning. To explicitly avoid
                  applying a scaling factor with a table, set '--scl_fct none'.
  -tv, --typ_cvg  Specify coverage calculation type (default: 'norm'). Options:
                    - 'raw', 'unadj', 'unadjusted': Compute unadjusted
                      coverage.
                    - 'len', 'len_frag': Normalize coverage by fragment length.
                    - 'norm', 'normalized': Per Dickson et al., Sci Rep 2023
                      (PMID: 37160995), normalize coverage by fragment length
                      and total fragments/unity.
  -uf, --usr_frg  Comma-separated string vector of fragment lengths to use
                  instead of read lengths (single-end alignments) or template
                  lengths (paired-end alignments; optional; default:
                  '${usr_frg}'). Must match the number of infiles via
                  '--infiles' or '--table'. To explicitly avoid applying a
                  fragment length, set '--usr_frg none'.
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (default: "\${dir_out}/err_out").
  -nj, --nam_job  Prefix for the names of jobs (default: '${nam_job}').
  -sl, --slurm    Submit jobs to the SLURM scheduler (optional; otherwise, if
                  'threads' > 1, run them via GNU Parallel, and if 'threads' is
                  1, run them in serial.)
  -mj, --max_job  The maximum number of jobs to run at one time (required if
                  '--slurm' is specified, ignored if not; default: ${max_job}).
  -tm, --time     The length of time, in 'h:mm:ss' format, for the SLURM job
                  (required if '--slurm' is specified, ignored if not; default:
                  '${time}').

Dependencies:
  - Programs
    + awk
    + Bash or Zsh
    + GNU Parallel (if not using SLURM and threads > 1)
    + Python
    + SLURM (if using --slurm)
  - Functions
    + check_array_files
    + check_arrays_lengths
    + check_exists_file_dir
    + check_flt_pos
    + check_format_time
    + check_int_pos
    + check_mut_excl_args
    + check_mut_excl_flags
    + check_program_path
    + check_str_delim
    + check_supplied_arg
    + check_table
    + check_table_column
    + check_table_scaling_factor
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
  - BAM infiles must be coordinate-sorted.
  - Outfile names are derived from BAM infiles and the value(s) associated with
    '--typ_out'.
  - #TODO

Examples:
  \`\`\`
  #  Example 1: Run with GNU Parallel and '--infiles'
  bash "\${HOME}/path/to/scripts/execute_compute_coverage.sh"
      --threads 8
      --infiles "\${HOME}/path/to/sample1.bam,\${HOME}/path/to/sample2.bam"
      --dir_out "\${HOME}/path/to/outfiles"
      --typ_out "bigwig"
      --siz_bin 50
      --typ_cvg "norm"
      --err_out "\${HOME}/path/to/outfiles/logs"
      --nam_job "coverage_calculation"

  #  Example 2: Run with SLURM and '--table' and '--tbl_col'
  bash "\${HOME}/path/to/scripts/execute_compute_coverage.sh"
      --threads 8
      --table "\${HOME}/path/to/table.tsv"
      --tbl_col "scaled"
      --dir_out "\${HOME}/path/to/outfiles"
      --typ_out "both"
      --siz_bin 5
      --slurm
      --max_job 10
      --time "0:20:00"
      --err_out "\${HOME}/path/to/outfiles/logs"
      --nam_job "slurm_coverage"

  #  Example 3: Run with SLURM and '--table' and '--typ_cvg'
  bash "\${HOME}/path/to/scripts/execute_compute_coverage.sh"
      --threads 8
      --table "\${HOME}/path/to/table.tsv"
      --typ_cvg "norm"
      --dir_out "\${HOME}/path/to/outfiles"
      --typ_out "bedgraph"
      --siz_bin 10
      --slurm
      --err_out "\${HOME}/path/to/outfiles/logs"
      --nam_job "compute_coverage"
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
            -sb|--siz_bin) siz_bin="${2}";   shift 2 ;;
            -sf|--scl_fct) scl_fct="${2}";   shift 2 ;;
            -tv|--typ_cvg) typ_cvg="${2}";   shift 2 ;;
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

#  Check that submission script is supplied and exists as a file
check_supplied_arg -a "${scr_sub}" -n "scr_sub"
check_exists_file_dir "f" "${scr_sub}" "scr_sub"

#  Check that coverage script is supplied and exists as a file
check_supplied_arg -a "${scr_cvg}" -n "scr_cvg"
check_exists_file_dir "f" "${scr_cvg}" "scr_cvg"

#  Check that denominator is supplied
check_supplied_arg -a "${denom}" -n "denom"

#  Check that threads value is supplied and is a positive integer
check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

#  Check that '--infiles' and '--table' are mutually exclusive
check_mut_excl_args "infiles" "${infiles}" "table" "${table}"

#  Validate '--infiles' or '--table'
if [[ -n "${infiles}" ]]; then
    #  Validate existence of directory containing infiles 
    check_exists_file_dir "d" "$(dirname "${infiles%%[,;]*}")" "infiles"

    #  Ensure '--infiles' is not empty or improperly formatted
    check_str_delim "infiles" "${infiles}"
elif [[ -n "${table}" ]]; then
    #  Validate existence of table file
    check_exists_file_dir "f" "${table}" "table"

    #  Validate table isn't empty by checking its number of lines
    check_table "${table}"

    #  Ensure '--tbl_col' is specified and valid
    check_supplied_arg -a "${tbl_col}" -n "tbl_col"
    case "${tbl_col}" in
        alpha|scaled|sf|none) : ;;  # Valid options
        *)
            echo_error \
                "Selection associated with '--tbl_col' is not valid:" \
                "${tbl_col}. Selection must be 'alpha', 'scaled', 'sf', or" \
                "'none'."
            exit_1
            ;;
    esac

    #  Skip further checks for '--tbl_col none'
    if [[ "${tbl_col}" != "none" ]]; then
        #  If present, trim leading/trailing whitespace from 'tbl_col'
        tbl_col=$(echo "${tbl_col}" | xargs)

        #  Validate 'tbl_col' exists in the table header
        check_table_column "${table}" "${tbl_col}"
    fi

    #  Warn if '--scl_fct' overrides '--tbl_col'
    if [[ -n "${scl_fct}" && "${scl_fct}" != "none" ]]; then
        if [[ "${tbl_col}" == "none" ]]; then
            echo_warning \
                "Scaling factors provided via '--scl_fct' override the" \
                "absence or non-specification a column in the table:" \
                "'--table ${table} --tbl_col none'."
        else
            echo_warning \
                "Scaling factors provided via '--scl_fct' override those in the" \
                "table: '--table ${table} --tbl_col ${tbl_col}'."
        fi
    fi
fi

#  Check that output directory is supplied and exists
check_supplied_arg -a "${dir_out}" -n "dir_out"
check_exists_file_dir "d" "${dir_out}" "dir_out"

#  Check that output type is supplied and valid
check_supplied_arg -a "${typ_out}" -n "typ_out"
case "${typ_out}" in
    bedgraph|bdg|bg|bigwig|bw|both|bed) : ;;  # Valid options
    *)
        echo_error \
            "Selection associated with '--typ_out' is not valid:" \
            "'${typ_out}'. Selection must be 'bedgraph', 'bdg', 'bg'," \
            "'bigwig', 'bw', 'both', or 'bed'."
        exit_1
        ;;
esac

#  Check that bin size is supplied and is a positive integer
check_supplied_arg -a "${siz_bin}" -n "siz_bin"
check_int_pos "${siz_bin}" "siz_bin"

#  If supplied, check that '--scl_fct' is properly formatted
if [[ -n "${scl_fct}" ]]; then
    if [[ "${scl_fct}" != "none" ]]; then
        check_str_delim "scl_fct" "${scl_fct}"
    fi
fi

#  If applicable, warn about '--scl_fct' overriding scaling factors in the table
if [[ -n "${table}" && -n "${scl_fct}" && "${scl_fct}" != "none" ]]; then
    if [[ "${tbl_col}" != "none" ]]; then
        echo_warning \
            "Scaling factors provided via '--scl_fct' override those in the" \
            "table: '--table ${table} --tbl_col ${tbl_col}'."
    else
        echo_warning \
            "Scaling factors provided via '--scl_fct' are being used even" \
            "though '--tbl_col none' is specified for '--table ${table}'."
    fi
fi

#  Default to '--typ_cvg raw' if no normalization is provided
if [[ -z "${typ_cvg}" ]]; then
    typ_cvg="raw"
    echo_warning \
        "No coverage normalization provided. Defaulting to 'raw' (i.e.,
        'non-normalized' or 'unadjusted') coverage ('--typ_cvg raw')."
else
    case "${typ_cvg}" in
        raw|unadj|unadjusted) : ;;  # Valid options for unadjusted coverage
        len|len_frag) : ;;          # Valid options for fragment-length normalization
        norm|normalized) : ;;       # Valid options for "normalized" coverage
        *)
            echo_error \
                "Invalid value for '--typ_cvg': '${typ_cvg}'. Expected one" \
                "of 'raw', 'unadj', 'unadjusted', 'len', 'len_frag', 'norm'," \
                "or 'normalized'."
            exit_1
            ;;
    esac
fi

#  If applicable, ensure '--usr_frg' is not empty or improperly formatted
if [[ -n "${usr_frg}" ]]; then
    if [[ "${usr_frg}" != "none" ]]; then
        check_str_delim "usr_frg" "${usr_frg}"
    fi
fi

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
if ${verbose}; then
    case "${typ_cvg}" in
        raw|unadj|unadjusted)
            mth_nrm="No normalization; returning raw (unadjusted) coverage:"
            mth_nrm+=" '--typ_cvg ${typ_cvg}'."
            ;;
        len|len_frag)
            mth_nrm="Performing fragment-length normalization: '--typ_cvg"
            mth_nrm+=" ${typ_cvg}'."
            ;;
        norm|normalized)
            mth_nrm="Generating normalized coverage (Dickson et al., Sci Rep"
            mth_nrm+=" 2023) via fragment-length and unity normalization:"
            mth_nrm+=" '--typ_cvg ${typ_cvg}'."
            ;;
        *)
            #  Should not be possible to see this
            mth_nrm="Unknown normalization method: '--typ_cvg ${typ_cvg}'."
            ;;
    esac

    if [[ "${scl_fct}" != "none" ]]; then
        src_scl="Custom multiplicative scaling factor(s):"
        src_scl+=" '--scl_fct ${scl_fct}'."
    elif [[
        -n "${table}" && -n "${tbl_col}" && "${tbl_col}" != "none"
    ]]; then
        src_scl="Custom multiplicative scaling factor(s) via table column:"
        src_scl+=" '--tbl_col ${tbl_col}'."
    else
        src_scl="No multiplicative scaling factor(s)."
    fi

    #  Print debug summary
    echo "##############################"
    echo "## Resolved Argument States ##"
    echo "##############################"
    echo ""
    echo "- Normalization method: ${mth_nrm}"
    echo "- Scaling factor source: ${src_scl}"
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
#  If applicable, parse values associated with '--table'  #TODO: Modularize
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
    unset arr_infile && typeset -a arr_infile
    unset arr_sf      && typeset -a arr_sf
    unset arr_scaled  && typeset -a arr_scaled
    unset arr_alpha   && typeset -a arr_alpha

    #  Parse the table rows, skipping the header
    while IFS=$'\t' read -r row; do
        IFS=$'\t' read -r -a fields <<< "${row}"

        #  Extract the 'sample' column (always present)
        arr_infile+=( "${fields[idx_smp]}" )

        #  Extract scaling factor columns if present
        [[ ${idx_sf}  -ne -1 ]] && arr_sf+=( "${fields[idx_sf]}" )
        [[ ${idx_scl} -ne -1 ]] && arr_scaled+=( "${fields[idx_scl]}" )
        [[ ${idx_alf} -ne -1 ]] && arr_alpha+=( "${fields[idx_alf]}" )
    done < <(awk 'NR > 1' "${table}")

    #  Dynamically assign arr_scl_fct based on tbl_col
    unset arr_scl_fct && typeset -a arr_scl_fct
    if [[ -z "${scl_fct}" && -z "${typ_cvg}" ]]; then
        eval "arr_scl_fct=( \"\${arr_${tbl_col}[@]}\" )"
    fi
fi

#  If applicable, parse '--infiles'
if [[ -n "${infiles}" ]]; then
    IFS=',' read -r -a arr_infile <<< "${infiles}"
fi

#  Having parsed '--table' or '--infiles', validate each infile exists
check_array_files "infiles/table" "${arr_infile[@]}"

#  Based on 'arr_infile' and 'dir_out', construct outfile paths, using them to
#+ populate 'arr_outfile'
unset arr_outfile && typeset -a arr_outfile
for i in "${arr_infile[@]}"; do
    if [[ "${tbl_col}" == "scaled" || "${tbl_col}" == "sf" ]]; then
        arr_outfile+=( "${dir_out}/$(basename "${i}" .bam).${tbl_col}" )
    else
        arr_outfile+=( "${dir_out}/$(basename "${i}" .bam)" )
    fi
done

#  If applicable, parse '--scl_fct'
if [[ -n "${scl_fct}" && "${scl_fct}" != "none" ]]; then
    IFS=',' read -r -a arr_scl_fct <<< "${scl_fct}"
fi

#  If empty or "none", populate '--scl_fct' array with "#N/A"
populate_array_empty arr_scl_fct "${#arr_infile[@]}"

#  Having parsed '--table' or '--scl_fct', validate each scaling factor as a
#+ positive float (if applicable)
for s in "${arr_scl_fct[@]}"; do
    if [[ "${s}" != "#N/A" ]]; then check_flt_pos "${s}" "scl_fct"; fi
done

#  If applicable, parse and validate '--usr_frg'
if [[ -n "${usr_frg}" && "${usr_frg}" != "none" ]]; then
    IFS=',' read -r -a arr_usr_frg <<< "${usr_frg}"
fi

#  If empty or "none", populate '--usr_frg' array with "#N/A"
populate_array_empty arr_usr_frg "${#arr_infile[@]}"

#  Having parsed '--usr_frg', validate each scaling factor as a positive float
#+ (if applicable)
for u in "${arr_usr_frg[@]}"; do
    if [[ "${u}" != "#N/A" ]]; then check_flt_pos "${u}" "usr_frg"; fi
done


#  Ensure matching element counts between 'arr_infile' and the other arrays
check_arrays_lengths "outfiles" arr_outfile "infiles" arr_infile
check_arrays_lengths "scl_fct"  arr_scl_fct  "infiles" arr_infile
check_arrays_lengths "usr_frg"  arr_usr_frg  "infiles" arr_infile

#  Reset 'max_job' if it is greater than the number of infiles
if ${slurm}; then
    max_job=$(reset_max_job "${max_job}" "${#arr_infile[@]}")
fi

#  Debug output for infiles and other values
if ${verbose}; then
    echo "########################################"
    echo "## Parsed vectors for parallelization ##"
    echo "########################################"
    echo ""
    debug_array_contents \
        "arr_infile" "arr_outfile" "arr_scl_fct" "arr_usr_frg"
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
    echo "scr_sub=${scr_sub}"
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
    echo "siz_bin=${siz_bin}"
    echo "scl_fct=${scl_fct:-#N/A}"
    echo "typ_cvg=${typ_cvg}"
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
    echo "arr_infile=( ${arr_infile[*]} )"
    echo "arr_outfile=( ${arr_outfile[*]} )"
    echo "arr_scl_fct=( ${arr_scl_fct[*]} )"
    echo "arr_usr_frg=( ${arr_usr_frg[*]} )"
    echo ""
    echo ""
fi

#  Having validated individual array elements, re- or newly assign variables
#+ with comma-separated string values from array elements
 infiles=$(echo "${arr_infile[*]}"  | tr ' ' ',')
outfiles=$(echo "${arr_outfile[*]}" | tr ' ' ',')
 scl_fct=$(echo "${arr_scl_fct[*]}" | tr ' ' ',')
 usr_frg=$(echo "${arr_usr_frg[*]}" | tr ' ' ',')

if ${verbose}; then
    echo "##################################################"
    echo "## Variable assignments constructed from arrays ##"
    echo "##################################################"
    echo ""
    echo "infiles=\"${infiles}\""
    echo "outfiles=\"${outfiles}\""
    echo "scl_fct=\"${scl_fct}\""
    echo "usr_frg=\"${usr_frg}\""
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
        echo "    --array=1-${#arr_infile[@]}%${max_job} \\"
        echo "        ${scr_sub} \\"
        echo "            ${env_nam} \\"
        echo "            ${scr_cvg} \\"
        echo "            ${threads} \\"
        echo "            ${infiles} \\"
        echo "            ${outfiles} \\"
        echo "            ${typ_out} \\"
        echo "            ${siz_bin} \\"
        echo "            ${typ_cvg} \\"
        echo "            ${scl_fct} \\"
        echo "            ${usr_frg} \\"
        echo "            ${err_out} \\"
        echo "            ${nam_job}"
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
        # echo "###########################################"
        # echo "## Contents of Python computation script ##"
        # echo "###########################################"
        # echo ""
        # echo "## ${scr_cvg} ##"
        # echo ""
        # cat "${scr_cvg}"
        # echo ""
        # echo ""
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
                    ${env_nam} \
                    ${scr_cvg} \
                    ${threads} \
                    ${infiles} \
                    ${outfiles} \
                    ${typ_out} \
                    ${siz_bin} \
                    ${typ_cvg} \
                    ${scl_fct} \
                    ${usr_frg} \
                    ${err_out} \
                    ${nam_job}
    fi
else
    #  GNU Parallel execution
    if [[ "${threads}" -gt 1 ]]; then
        #  Generate a configuration file
        config="${err_out}/${nam_job}.${RANDOM}.txt"
        touch "${config}" || {
            echo_error "Failed to create a GNU Parallel configuration file."
            exit_1
        }

        #  Populate the GNU Parallel configuration file with parameters for
        #+ each job
        for i in "${!arr_infile[@]}"; do
            echo \
                "${env_nam}" \
                "${scr_cvg}" \
                "${threads}" \
                "${arr_infile[i]}" \
                "${arr_outfile[i]}" \
                "${typ_out}" \
                "${siz_bin}" \
                "${typ_cvg}" \
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
                "bash \"${scr_sub}\" {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12}" \
                :::: "${config}"
        fi

        #  Use GNU Parallel for multi-threaded execution without SLURM
        if ! ${dry_run}; then
            # shellcheck disable=SC2090
            parallel --colsep ' ' --jobs "${par_job}" \
                "bash \"${scr_sub}\" {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12}" \
                :::: "${config}"
        fi
    else
        pth_std="${err_out}/${nam_job}"
        if ${dry_run} || ${verbose}; then
            echo "###########################################"
            echo "## Serial call(s) for computing coverage ##"
            echo "###########################################"
            echo ""
            echo "bash ${scr_sub} \\"
            echo "    ${env_nam} \\"
            echo "    ${scr_cvg} \\"
            echo "    ${threads} \\"
            echo "    ${infiles} \\"
            echo "    ${outfiles} \\"
            echo "    ${typ_out} \\"
            echo "    ${siz_bin} \\"
            echo "    ${typ_cvg} \\"
            echo "    ${scl_fct} \\"
            echo "    ${usr_frg} \\"
            echo "    ${err_out} \\"
            echo "    ${nam_job} \\"
            echo "         > ${pth_std}.stdout.txt \\"
            echo "        2> ${pth_std}.stderr.txt"
            echo ""
            echo ""
        fi

        if ! ${dry_run}; then
            bash "${scr_sub}" \
                "${env_nam}" \
                "${scr_cvg}" \
                "${threads}" \
                "${infiles}" \
                "${outfiles}" \
                "${typ_out}" \
                "${siz_bin}" \
                "${typ_cvg}" \
                "${scl_fct}" \
                "${usr_frg}" \
                "${err_out}" \
                "${nam_job}" \
                     > "${pth_std}.stdout.txt" \
                    2> "${pth_std}.stderr.txt"

            # for i in "${!arr_infile[@]}"; do
            #     bash "${scr_sub}" \
            #         "${env_nam}" \
            #         "${scr_cvg}" \
            #         "${threads}" \
            #         "${arr_infile[i]}" \
            #         "${arr_outfile[i]}" \
            #         "${typ_out}" \
            #         "${siz_bin}" \
            #         "${typ_cvg}" \
            #         "${arr_scl_fct[i]}" \
            #         "${arr_usr_frg[i]}" \
            #         "${err_out}" \
            #         "${nam_job}" \
            #              > "${err_out}/${nam_job}.$(basename "${arr_outfile[i]}").stdout.txt" \
            #             2> "${err_out}/${nam_job}.$(basename "${arr_outfile[i]}").stderr.txt"
            # done
        fi
    fi
fi
