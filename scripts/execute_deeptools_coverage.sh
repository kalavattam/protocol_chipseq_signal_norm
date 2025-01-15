#!/bin/bash

#  execute_deeptools_coverage.sh
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
    source "${dir_fnc}/check_array_files.sh"
    source "${dir_fnc}/check_arrays_lengths.sh"
    source "${dir_fnc}/check_exists_file_dir.sh"
    source "${dir_fnc}/check_flt_pos.sh"
    source "${dir_fnc}/check_format_time.sh"
    source "${dir_fnc}/check_int_pos.sh"
    source "${dir_fnc}/check_mut_excl_args.sh"
    source "${dir_fnc}/check_mut_excl_flags.sh"
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/check_region.sh"
    source "${dir_fnc}/check_region_bam.sh"  #TODO
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
    # dir_bam="${dir_aln}/${det_bam}/sc"
    dir_cov="${dir_pro}/compute_coverage"
    dir_tbl="${dir_cov}/${det_cov}/alpha/tables"  # "spike"
    dir_trk="${dir_cov}/${det_cov}/alpha/tracks"  # "spike"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=false
    threads=8
    # infiles="$(  ## WARNING: Change search parameters as needed ##
    #     bash "${dir_scr}/find_files.sh" \
    #         --dir_fnd "${dir_bam}" \
    #         --pattern "*.bam" \
    #         --exclude "in*,*Brn1*"
    # )"
    table="${dir_tbl}/IP_WT_G1-G2M-Q_Hho1-Hmo1_6336-6337_7750-7751.tsv"  # ""  # "${dir_tbl}/IP_WT_log-Q_Brn1_rep1-rep2-rep3.tsv"  # "${dir_tbl}/IP_WT_G1-G2M-Q_Hho1_6336-6337.tsv"  # "${dir_tbl}/IP_WT_G1-G2M-Q_Hmo1_7750-7751.tsv"
    tbl_col="alpha"  # "scaled"  # ""
    dir_out="${dir_trk}"
    typ_out="bedgraph"  # "bigwig"
    siz_bin=10  # 1
    siz_gen=12157105
    region=""
    scl_fct=""  # "0.002054,0.003138,0.003127,0.003522,0.056611,0.02906"
    typ_cvg=""  # "raw"  # "none"  # "rpkm"  # "fpkm"  # "cpm"  # "bpm"  # "rpgc"
    exact=false
    usr_frg=""
    err_out="${dir_trk}/logs"
    nam_job="run_deeptools_coverage"
    slurm=true
    max_job=8  # 4  # 6
    time="0:30:00"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_analyze"
scr_slm="${dir_scr}/submit_deeptools_coverage_slurm.sh"
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
siz_bin=10
siz_gen=12157105
region=""
scl_fct=""
typ_cvg=""
exact=false
usr_frg=""
err_out=""
nam_job="run_deeptools_coverage"
slurm=false
max_job=6
time="0:30:00"

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  execute_deeptools_coverage.sh
    [--verbose] [--dry_run] --threads <int> --infiles <str> [--table <str>]
    [--tbl_col <str>] --dir_out <str> --typ_out <str> --siz_bin <int>
    [--siz_gen <int>] [--region <str>] --scl_fct <flt> --typ_cvg <str>
    [--exact] [--usr_frg <flt>] --err_out <str> --nam_job <str> [--slurm]
    [--max_job <int>] [--time <str>]

Description:
  This driver script, 'execute_deeptools_coverage.sh', automates the creation
  of coverage signal tracks (e.g., BIGWIG or BEDGRAPH files) from coordinate-
  sorted BAM input files. It supports scaling normalization using spike-in,
  siQ-ChIP (alpha), or other methods based on aligned read depth (for more
  details, refer to the documentation for deepTools utility 'bamCoverage'). The
  script can execute jobs in serial or parallel using SLURM or GNU Parallel.

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
                  or 'alpha' (required if '--table' is invoked, ignored if not;
                  default: "${tbl_col}").
  -do, --dir_out  The directory to write coverage signal track outfiles.
  -to, --typ_out  Format of coverage signal track outfiles: 'bedgraph',
                  'bigwig', or 'both' (default: "${typ_out}").
  -sb, --siz_bin  Bin size for coverage calculation in base pairs (default:
                  ${siz_bin}).
  -sg, --siz_gen  Effective genome size of model organism (required if
                  '--typ_cvg rpgc', ignored if not; default: ${siz_gen}).
   -r, --region   Region of the genome to plot in 'chr' or 'chr:start-stop'
                  format (optional).
  -sf, --scl_fct  Comma-separated string vector of scaling factors to apply to
                  coverage (optional). Must match the number of infiles via
                  '--infiles' or '--table'. If used with '--table', these
                  values override those in the specified column (via
                  '--tbl_col') with a warning.
  -tv, --typ_cvg  Specify a type of coverage normalizatio to in coverage
                  computation; available options are: 'raw', 'none', 'rpkm',
                  'fpkm', 'cpm', 'bpm', or 'rpgc'. If used with '--table', the
                  normalization scaling factors override values in the column
                  specified by '--tbl_col', and a warning is issued.
  -ex, --exact    Compute scaling factors using all alignments instead of a
                  sampled subset (optional). Only applicable if --typ_cvg is
                  specified, and ignored if --typ_cvg is not invoked.
  -uf, --usr_frg  Comma-separated string vector of fragment lengths to use
                  instead of read lengths (single-end alignments) or template
                  lengths (paired-end alignments; optional). Must match the
                  number of infiles via --infiles or --table.
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (default: \${dir_out}/err_out).
  -nj, --nam_job  Prefix for the names of jobs (default: "${nam_job}").
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
    + deepTools bamCoverage
    + GNU Parallel (if not using SLURM and threads > 1)
    + Python
    + Samtools
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
    + check_region
    + check_region_bam
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
  - Running the script with '--exact' significantly increases the time required
    for coverage computation.
  - #TODO

Examples:
  \`\`\`
  #  Run with a table file and column
  execute_compute_coverage.sh
    --threads 8
    --table sample_scaling_factors.tsv
    --tbl_col alpha
    --dir_out /path/to/outdir
    --typ_out bigwig
    --siz_bin 10
    --slurm
    --max_job 10
    --time 1:00:00

  #  Run with specific BAM infiles and scaling factors
  execute_compute_coverage.sh
    -t 4
    -i file1.bam,file2.bam
    -sf 1.0,0.8
    -do /path/to/output
    -to bedgraph
    -sb 5
    -tv raw
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
            -sg|--siz_gen) siz_gen="${2}";   shift 2 ;;
             -r|--region)  region="${2}";    shift 2 ;;
            -sf|--scl_fct) scl_fct="${2}";   shift 2 ;;
            -tv|--typ_cvg) typ_cvg="${2,,}"; shift 2 ;;
            -ex|--exact)   exact=true;       shift 1 ;;
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

#  Check that denominator is supplied
check_supplied_arg -a "${denom}" -n "denom"

#  Check that a positive integer value is supplied to 'threads'
check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

#  Check that '--infiles' and '--table' are mutually exclusive; one must be
#+ specified
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

    #  When '--table' is specified, validate '--tbl_col'
    check_supplied_arg -a "${tbl_col}" -n "tbl_col"
    case "${tbl_col}" in
        alpha|scaled|sf) : ;;
        *)
            echo_error \
                "Selection associated with --tbl_col is not valid:" \
                "'${tbl_col}'. Selection must be 'alpha', 'scaled', or 'sf'."
            exit_1
            ;;
    esac

    #  If present, trim leading/trailing whitespace from 'tbl_col'
    tbl_col=$(echo "${tbl_col}" | xargs)

    #  Validate 'tbl_col' exists in the table header
    check_table_column "${table}" "${tbl_col}"
fi

#  Check that output directory is supplied and exists
check_supplied_arg -a "${dir_out}" -n "dir_out"
check_exists_file_dir "d" "${dir_out}" "dir_out"

#  Check that output type is supplied and valid
check_supplied_arg -a "${typ_out}" -n "typ_out"
case "${typ_out}" in
    bedgraph|bdg|bg|bigwig|bw|both) : ;;
    *)
        echo_error \
            "Selection associated with '--typ_out' is not valid:" \
            "'${typ_out}'. Selection must be 'bedgraph', 'bdg', or 'bg';" \
            "'bigwig' or 'bw'; or 'both'."
        exit_1
        ;;
esac

#  Check that bin size is supplied and is a positive integer
check_supplied_arg -a "${siz_bin}" -n "siz_bin"
check_int_pos "${siz_bin}" "siz_bin"

#  Check that effective genome size is supplied and a positive integer
check_supplied_arg -a "${siz_gen}" -n "siz_gen"
check_int_pos "${siz_gen}" "siz_gen"

#  If supplied, check that region is properly formatted
if [[ -n "${region}" ]]; then check_region "${region}"; fi

#  Validate '--scl_fct'
if [[ -n "${scl_fct}" ]]; then
    #  Ensure compatibility of '--table' and '--scl_fct'
    check_table_scaling_factor "str" "${table}" "${scl_fct}" "scl_fct"
    
    #  Check '--scl_fct' is properly formatted
    check_str_delim "scl_fct" "${scl_fct}"
fi

#  Validate '--typ_cvg'
if [[ -n "${typ_cvg}" ]]; then
    #  Parse and validate '--typ_cvg' assignment
    case "${typ_cvg}" in
        raw|none)  typ_cvg="None" ;;
        rpkm|fpkm) typ_cvg="RPKM" ;;
        cpm)       typ_cvg="CPM"  ;;
        bpm)       typ_cvg="BPM"  ;;
        rpgc)      typ_cvg="RPGC" ;;
        *)
            echo_error \
                "Selection associated with --typ_cvg is not valid: '${typ_cvg}'." \
                "Assignment must be 'raw', 'none', 'rpkm', 'fpkm', 'cpm'," \
                "'bpm', or 'rpgc'."
            exit_1
            ;;
    esac

    #  If table is provided, ensure compatibility of '--table' and '--typ_cvg'
    if [[ -n "${table}" ]]; then
        check_table_scaling_factor "str" "${table}" "${typ_cvg}" "typ_cvg"
    fi
else
    typ_cvg="None"
fi

#  If applicable, ensure '--usr_frg' is not empty or improperly formatted
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
#  time ('max_job') and the maximum permitted runtime ('time')
if ${slurm}; then
    #  Ensure 'max_job' is provided and is a positive integer
    check_supplied_arg -a "${max_job}" -n "max_job"
    check_int_pos "${max_job}" "max_job"
    
    #  Ensure 'time' is provided and has a valid time format
    check_supplied_arg -a "${time}" -n "time"
    check_format_time "${time}"
else
    #  For non-SLURM jobs, calculate the number of parallel jobs to run
    #+ based on the number of threads, setting the thread count accordingly
    par_job=$(( threads / denom ))
    threads=${denom}

    #  Ensure the calculated number of parallel jobs (par_job) is provided and
    #+ is a positive integer
    check_supplied_arg -a "${par_job}" -n "par_job"
    check_int_pos "${par_job}" "par_job"

    #  Ensure the number of threads is provided and is a positive integer
    check_supplied_arg -a "${threads}" -n "threads"
    check_int_pos "${threads}" "threads"
fi


#  Debug summary output of resolved argument states ---------------------------
#+ ...i.e., for '--typ_cvg', '--table', '--tbl_col', and '--scl_fct'
if ${verbose}; then
    if [[ -n "${typ_cvg}" ]]; then
        if [[ -n "${scl_fct}" ]]; then
            src_scl="--scl_fct '${scl_fct}' and --typ_cvg"
            src_scl+=" '${typ_cvg/None/raw}'"
            mth_nrm="User-supplied scaling factors × '${typ_cvg/None/raw}'"

            #  --scl_fct takes precedence when both --scl_fct and --tbl_col are
            #+ provided
            if [[ -n "${table}" && -n "${tbl_col}" ]]; then
                src_scl+=" (overriding '--tbl_col ${tbl_col}')"
                mth_nrm+=" (overriding table column '${tbl_col}')"
            fi
        elif [[ -n "${table}" && -n "${tbl_col}" ]]; then
            src_scl="--table (--tbl_col '${tbl_col}') and --typ_cvg"
            src_scl+=" '${typ_cvg/None/raw}'"
            mth_nrm="'${tbl_col}' × '${typ_cvg/None/raw}'"
        else
            src_scl="--typ_cvg '${typ_cvg/None/raw}'"
            mth_nrm="Depth-based normalization using '${typ_cvg/None/raw}'"
        fi
    elif [[ -n "${scl_fct}" ]]; then
        src_scl="--scl_fct '${scl_fct}'"
        mth_nrm="User-supplied scaling factors via --scl_fct"

        if [[ -n "${table}" && -n "${tbl_col}" ]]; then
            src_scl+=" (overriding '--tbl_col ${tbl_col}')"
            mth_nrm+=" (overriding table column '${tbl_col}')"
        fi
    elif [[ -n "${table}" && -n "${tbl_col}" ]]; then
        src_scl="--table (--tbl_col '${tbl_col}')"
        mth_nrm="Sample-specific scaling factors via table column '${tbl_col}'"
    else
        echo_error \
            "Neither --typ_cvg nor --scl_fct nor --table (with --tbl_col)" \
            "was specified. This violates input requirements."
        exit 1
    fi

    echo "##############################"
    echo "## Resolved Argument States ##"
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
check_program_path bamCoverage
if ! ${slurm}; then check_program_path parallel; fi
check_program_path python
check_program_path samtools
if ${slurm}; then check_program_path sbatch; fi


#  Parse and validate table and/or vector elements ----------------------------
#  If applicable, parse values associated with --table
if [[ -n "${table}" ]]; then  #TODO: Modularize this
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
fi

#  If applicable, parse --infiles
if [[ -n "${infiles}" ]]; then
    IFS=',' read -r -a arr_infiles <<< "${infiles}"
fi

#  Having parsed --table or --infiles, validate each infile exists
check_array_files "infiles/table" "${arr_infiles[@]}"

#  Based on arr_infiles and dir_out, construct outfile paths, using them to
#+ populate arr_outfiles
unset arr_outfiles && typeset -a arr_outfiles
for i in "${arr_infiles[@]}"; do
    arr_outfiles+=( "${dir_out}/$(basename "${i}" .bam)" )
done

#  Dynamically assign arr_scl_fct
unset arr_scl_fct && typeset -a arr_scl_fct
if [[ -n "${table}" ]]; then
    if [[ -n "${scl_fct}" ]]; then
        #  Assign user-supplied scaling factors, overriding table scaling
        #+ factors if present
        IFS=',' read -r -a arr_scl_fct <<< "${scl_fct}"
    else
        #  Use scaling factors from table column
        if [[ -n "${tbl_col}" ]]; then
            eval "arr_scl_fct=( \"\${arr_${tbl_col}[@]}\" )"
        fi
    fi
elif [[ -n "${infiles}" ]]; then
    if [[ -n "${scl_fct}" ]]; then
        #  Assign user-supplied scaling factors
        IFS=',' read -r -a arr_scl_fct <<< "${scl_fct}"
    fi
fi

#  Populate 'arr_scl_fct' with placeholders if it is empty
populate_array_empty arr_scl_fct "${#arr_infiles[@]}"

#  Having parsed --table or --scl_fct, validate each scaling factor as a
#+ positive float (if applicable)
if [[ "${arr_scl_fct[0]}" != "#N/A" ]]; then
    for s in "${arr_scl_fct[@]}"; do check_flt_pos "${s}" "scl_fct"; done
fi

#  If applicable, parse and validate --usr_frg
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
check_arrays_lengths "outfiles" arr_outfiles "infiles" arr_infiles
check_arrays_lengths "scl_fct"  arr_scl_fct  "infiles" arr_infiles
check_arrays_lengths "usr_frg"  arr_usr_frg  "infiles" arr_infiles

#  Reset max_job if it is greater than the number of infiles
if ${slurm}; then
    max_job=$(reset_max_job "${max_job}" "${#arr_infiles[@]}")
fi

#  Debug output for infiles and other values
if ${verbose}; then
    echo "########################################"
    echo "## Parsed vectors for parallelization ##"
    echo "########################################"
    echo ""
    debug_array_contents \
        "arr_infiles" "arr_outfiles" "arr_scl_fct" "arr_usr_frg"
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
    echo "region=${region:-#N/A}"
    echo "scl_fct=${scl_fct:-#N/A}"
    echo "typ_cvg=${typ_cvg}"
    echo "exact=${exact}"
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
    echo "arr_scl_fct=( ${arr_scl_fct[*]} )"
    echo "arr_usr_frg=( ${arr_usr_frg[*]} )"
    echo ""
    echo ""
fi

#TODO: Used the wrong equation to compute alpha, but test anyway
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
        echo "             -v ${verbose} \\"
        echo "             -t ${threads} \\"
        echo "            -si $(echo "${arr_infiles[*]}" | tr ' ' ',') \\"
        echo "            -so $(echo "${arr_outfiles[*]}" | tr ' ' ',') \\"
        echo "            -to ${typ_out} \\"
        echo "            -sb ${siz_bin} \\"
        echo "             -r ${region:-#N/A} \\"
        echo "            -sf $(echo "${arr_scl_fct[*]}" | tr ' ' ',') \\"
        echo "            -tv ${typ_cvg} \\"
        echo "            -ex ${exact} \\"
        echo "            -su $(echo "${arr_usr_frg[*]}" | tr ' ' ',') \\"
        echo "            -eo ${err_out} \\"
        echo "            -nj ${nam_job} \\"
        echo "            -en ${env_nam}"
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
                     -v ${verbose} \
                     -t ${threads} \
                    -si $(echo "${arr_infiles[*]}" | tr ' ' ',') \
                    -so $(echo "${arr_outfiles[*]}" | tr ' ' ',') \
                    -to ${typ_out} \
                    -sb ${siz_bin} \
                     -r ${region:-#N/A} \
                    -sf $(echo "${arr_scl_fct[*]}" | tr ' ' ',') \
                    -tv ${typ_cvg} \
                    -ex ${exact} \
                    -su $(echo "${arr_usr_frg[*]}" | tr ' ' ',') \
                    -eo ${err_out} \
                    -nj ${nam_job} \
                    -en ${env_nam}
    fi
fi

#TODO: Implement GNU Parallel and serial job submissions
#TODO: Add siz_gen where appropriate
