#!/bin/bash

#  execute_deeptools_compare.sh
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
    source "${dir_fnc}/check_array_files.sh"
    source "${dir_fnc}/check_arrays_lengths.sh"
    source "${dir_fnc}/check_exists_file_dir.sh"
    source "${dir_fnc}/check_flt_pos.sh"
    source "${dir_fnc}/check_format_time.sh"
    source "${dir_fnc}/check_int_pos.sh"
    source "${dir_fnc}/check_match.sh"
    source "${dir_fnc}/check_mut_excl_args.sh"
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/check_region.sh"
    source "${dir_fnc}/check_region_bam.sh"      #TODO
    source "${dir_fnc}/check_seq_type.sh"
    source "${dir_fnc}/check_str_delim.sh"
    source "${dir_fnc}/check_supplied_arg.sh"
    source "${dir_fnc}/debug_array_contents.sh"
    source "${dir_fnc}/echo_error.sh"
    source "${dir_fnc}/echo_warning.sh"          #TODO
    source "${dir_fnc}/exit_0.sh"
    source "${dir_fnc}/exit_1.sh"
    source "${dir_fnc}/handle_env.sh"
    source "${dir_fnc}/populate_array_empty.sh"
    source "${dir_fnc}/reset_max_job.sh"
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
    aligner="bowtie2"
    a_type="global"
    req_flg=true
    flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
    mapq=1
    det_bam="flag-${flg}_mapq-${mapq}"
    det_cov="${aligner}_${a_type}_${det_bam}"
    typ_fil="bigwig"
    case "${typ_fil}" in
        bedgraph|bdg|bg|bigwig|bw)
            dir_cov="${dir_pro}/compute_coverage/${det_cov}"
            ;;
        bam)
            dir_aln="${dir_pro}/align_${aligner}_${a_type}"
            dir_bam="${dir_aln}/${det_bam}/sc"
            ;;
    esac
    cov_nrm="norm"
    dir_bdg="${dir_cov}/${cov_nrm}/tracks"
    dir_bwg="${dir_cov}/${cov_nrm}/tracks"
    dir_trk="${dir_cov}/log2/${cov_nrm}/tracks"
    include="*Hho1*"  # "*Hmo1*"  # "*Brn1*"
    case "${typ_fil}" in
        bedgraph|bdg|bg)
            dir_fnd="${dir_bdg}"
            pattern="*.bg"
            ;;
        bigwig|bw)
            dir_fnd="${dir_bwg}"
            pattern="*.bw"
            ;;
        bam)
            dir_fnd="${dir_bam}"
            pattern="*.bam"
            ;;
    esac

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=true
    threads=8
    fil_num="$(
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_fnd}" \
            --pattern "${pattern}" \
            --include "IP*,${include}"
    )"
    fil_den="$(
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_fnd}" \
            --pattern "${pattern}" \
            --include "in*,${include}"
    )"
    fil_stm=$(
        echo "${fil_num}" \
            | sed -e "s:${dir_fnd}:${dir_trk}:g" \
                  -e 's:/IP:/IP-in:g' \
                  -e 's:\.\(bigwig\|bw\|bedgraph\|bdg\|bg\|bam\)::g'
    )
    typ_out="bigwig"
    bin_siz=1  # 10
    region=""
    oper="log2"
    scl_fct=""
    norm=""
    scl_pre=""
    exact=true
    usr_frg=""
    err_out="${dir_trk}/logs"
    nam_job="run_deeptools_compare"
    slurm=true
    max_job=4  # 6
    time="0:30:00"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_analyze"
scr_slm="${dir_scr}/submit_deeptools_compare_slurm.sh"
scr_par="${dir_scr}/submit_deeptools_compare_parallel.sh"
denom=4

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=1
fil_num=""
fil_den=""
fil_stm=""
typ_out="bigwig"
bin_siz=10
region=""
oper="log2"
scl_fct=""
norm=""
scl_pre=""
exact=false
usr_frg=""
err_out=""
nam_job="run_deeptools_compare"
slurm=false
max_job=6
time="0:30:00"

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  execute_compute_compare.sh
    [--verbose] [--dry_run] --threads <int> --fil_num <str> --fil_den <str>
    --fil_stm <str> --typ_out <str> --bin_siz <int> [--region <str>]
    --scl_fct <flt> --norm <str> --scl_pre <str> [--exact] [--usr_frg <flt>]
    --err_out <str> --nam_job <str> [--slurm] [--max_job <int>] [--time <str>]

Description:
  execute_deeptools_compare.sh automates...

Arguments:
   -h, --help     Print this help message and exit.
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Run script in 'dry-run mode' (optional).
   -t, --threads  Number of threads to use (default: ${threads}).
  -fn, --fil_num  Comma-separated string of BAM, BIGWIG, or BEDGRAPH infiles to
                  be used as the first file/numerator in the comparisons.
  -fd, --fil_den  Comma-separated string of BAM, BIGWIG, or BEDGRAPH infiles to
                  be used as the second file/denominator in the comparisons.
  -fs, --fil_stm  Comma-separated string of outstem sans extensions.
  -to, --typ_out  Format of comparison signal track outfiles: 'bedgraph' or
                  'bigwig' (default: '${typ_out}').
  -bs, --bin_siz  Bin size for comparison calculation in base pairs (default:
                  ${bin_siz}).
   -r, --region   Specify a genomic region to limit the operation (optional).
                  The format can be either: 'chr' (entire chromosome) or
                  'chr:start-stop' (specific range within a chromosome).
  -op, --oper     Operation to perform for comparisons of first/numerator and
                  second/denominator BAM infiles: 'log2', 'ratio', 'subtract',
                  'add', 'mean', 'reciprocal_ratio', 'first', 'second'
                  (default: '${oper}').
  -sf, --scl_fct  Comma-separated string of scaling factors for numerator and
                  denominator files. Each scaling factor pair should be in the
                  format 'num:den', where 'num' is the scaling factor for the
                  first file in a pair/numerator and 'den' is the scaling
                  factor for the second file in a pair/denominator. Multiple
                  pairs should be separated by commas. For example,
                  '0.7:1,0.3:0.6' scales the first numerator file by 0.7 and
                  the first denominator file by 1.0, and the second numerator
                  file by 0.3 and the second denominator file by 0.6. This
                  argument cannot be used with '--norm <str>' or
                  '--scl_pre <str>'.
  -no, --norm     Specify a normalization method to use when computing
                  comparisons; available options are: 'raw', 'none', 'rpkm',
                  'fpkm', 'cpm', 'bpm', or 'rpgc'. This argument cannot be used
                  with '--scl_fct <str>' or '--scl_pre <str>'.
  -sp, --scl_pre  Predefined scaling method for compensating sequencing depth
                  differences in comparison operations: 'depth', 'ses', or
                  'none'. This argument cannot be used with '--norm <str>' or
                  '--scl_fct <str>'.
  -ex, --exact    Compute scaling factors using all alignments instead of a
                  sampled subset. Only applicable if '--norm <str>' is
                  specified, and ignored if '--norm "none"' or '--norm' is not
                  invoked. Significantly increases the time required for
                  comparison computation.
  -uf, --usr_frg  Comma-separated string vector of fragment lengths to use
                  instead of read lengths (single-end alignments) or template
                  lengths (paired-end alignments; optional). Must match the
                  number of infiles and outstems via '--fil_num <str>',
                  '--fil_den <str>', and '--fil_stm <str>'.
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (default: '\$(dirname "\${fil_stm%%[,;]*}")/err_out').
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
    + deepTools bamCoverage
    + GNU Parallel (if not using SLURM and threads > 1)
    + Python
    + Samtools
    + SLURM (if using --slurm)
  - Functions
    + #TODO

Notes:
  - When the '--slurm' flag is used, jobs are parallelized via SLURM array
    tasks; otherwise, jobs are parallelized with GNU Parallel.
  - BAM infiles must be coordinate-sorted.
  - Outfile names are derived BAM infiles and the value(s) associated with
    --typ_out.
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
             -v|--verbose) verbose=true;     shift 1 ;;
            -dr|--dry_run) dry_run=true;     shift 1 ;;
             -t|--threads) threads="${2}";   shift 2 ;;
            -fn|--fil_num) fil_num="${2}";   shift 2 ;;
            -fd|--fil_den) fil_den="${2}";   shift 2 ;;
            -fs|--fil_stm) fil_stm="${2}";   shift 2 ;;
            -to|--typ_out) typ_out="${2,,}"; shift 2 ;;
            -bs|--bin_siz) bin_siz="${2}";   shift 2 ;;
             -r|--region)  region="${2}";    shift 2 ;;
            -op|--oper)    oper="${2}";      shift 2 ;;
            -sf|--scl_fct) scl_fct="${2}";   shift 2 ;;
            -no|--norm)    norm="${2,,}";    shift 2 ;;
            -sp|--scl_pre) scl_pre="${2,,}"; shift 2 ;;
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
check_supplied_arg -a "${env_nam}" -n "env_nam"

check_supplied_arg -a "${scr_slm}" -n "scr_slm"
check_exists_file_dir "f" "${scr_slm}" "scr_slm"

# check_supplied_arg -a "${scr_par}" -n "scr_par"   #TODO: Implement "scr_par"
# check_exists_file_dir "f" "${scr_par}" "scr_par"

check_supplied_arg -a "${denom}" -n "denom"

check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

check_supplied_arg -a "${fil_num}" -n "fil_num"
check_exists_file_dir "d" "$(dirname "${fil_num%%[,;]*}")" "fil_num"

check_supplied_arg -a "${fil_den}" -n "fil_den"
check_exists_file_dir "d" "$(dirname "${fil_den%%[,;]*}")" "fil_den"

check_supplied_arg -a "${fil_stm}" -n "fil_stm"
check_exists_file_dir "d" "$(dirname "${fil_stm%%[,;]*}")" "fil_stm"

check_supplied_arg -a "${typ_out}" -n "typ_out"
case "${typ_out}" in
    bedgraph|bdg|bg) typ_out="bedgraph" ;;
    bigwig|bw)       typ_out="bigwig"   ;;
    *)
        echo_error \
            "Selection associated with --typ_out is not valid: ${typ_out}." \
            "Selection must be 'bedgraph' or 'bigwig'."
        exit_1
        ;;
esac

check_supplied_arg -a "${bin_siz}" -n "bin_siz"
check_int_pos "${bin_siz}" "bin_siz"

if [[ -n "${region}" ]]; then check_region "${region}"; fi

case "${oper}" in
    log*)         oper="log2"             ;;
    1|1st|first)  oper="first"            ;;
    2|2nd|second) oper="second"           ;;
    mean|mn)      oper="mean"             ;;
    rat*)         oper="ratio"            ;;
    sub*)         oper="subtract"         ;;
    add*)         oper="add"              ;;
    recip*)       oper="reciprocal_ratio" ;;
    *)
        echo_error \
            "Invalid operation '${oper}'. Supported operations are 'log2'," \
            "'ratio', 'subtract', 'add', 'mean', 'reciprocal_ratio'," \
            "'first', or 'second'." >&2
        exit_1
        ;;
esac

#  Check that only one of --scl_fct, --norm, or --scl_pre is assigned
#MAYBE: Allow combinations, i.e., 'asmgt > 1'?
asmgt=0
if [[ -n "${scl_fct}" ]]; then (( asmgt++ )) || true; fi
if [[ -n "${norm}" ]];    then (( asmgt++ )) || true; fi
if [[ -n "${scl_pre}" ]]; then (( asmgt++ )) || true; fi

if (( asmgt > 1 )); then
    echo \
        "Error: Only one of --scl_fct, --norm, or --scl_pre can be" \
        "specified." >&2
    exit_1
fi

if [[ -n "${scl_fct}" ]]; then
    check_str_delim "scl_fct" "${scl_fct}"
elif [[ -n "${norm}" ]]; then
    case "${norm}" in
         raw|none) norm="None" ;;  #MAYBE: Disallow this selection
        rpkm|fpkm) norm="RPKM" ;;
              cpm) norm="CPM"  ;;
              bpm) norm="BPM"  ;;
             rpgc) norm="RPGC" ;;
        *)
            echo_error \
                "Selection associated with '--norm' is not valid: '${norm}'." \
                "Assignment must be 'raw', 'none', 'rpkm', 'fpkm', 'cpm'," \
                "'bpm', or 'rpgc'."
            exit_1
            ;;
    esac
elif [[ -n "${scl_pre}" ]]; then
    case "${scl_pre}" in
        depth) scl_pre="readCount" ;;
          ses) scl_pre="SES"       ;;
         none) scl_pre="None"      ;;  #MAYBE: Disallow this selection
        *)
            echo_error \
                "Selection associated with '--scl_pre' is not valid:" \
                "'${scl_pre}'. Assignment must be 'depth', 'ses', or 'none'."
            exit_1
            ;;
    esac
fi

if [[ -n "${usr_frg}" ]]; then check_str_delim "usr_frg" "${usr_frg}"; fi

if [[ -n "${err_out}" ]]; then
    check_exists_file_dir "d" "${err_out}" "err_out"
elif [[ -z "${err_out}" ]]; then
    err_out="$(dirname "${fil_stm%%[,;]*}")/err_out"
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


#  Debug summary output of resolved argument states ---------------------------
#+ ...i.e., for --oper, --scl_fct, --norm, and --scl_pre
#MAYBE: Allow combinations, i.e., 'asmgt > 1'?
if ${verbose}; then
    if [[ -n "${norm}" ]]; then
        src_scl="--norm"
        if [[ "${norm}" != "none" ]]; then
            mth_nrm="Coverage normalization via --norm '${norm}'"
        else
            mth_nrm="Raw (unadjusted) coverage: --norm '${norm}'"
        fi
    elif [[ -n "${scl_fct}" ]]; then
        src_scl="--scl_fct"
        mth_nrm="User-supplied scaling via --scl_fct '${scl_fct}'"
    elif [[ -n "${scl_pre}" ]]; then
        src_scl="--scl_pre"
        if [[ "${scl_pre}" != "none" ]]; then
            mth_nrm="deepTools-calculated scaling via --scl_pre '${scl_pre}'"
        else
            mth_nrm="Unscaled data: --scl_pre '${scl_pre}'"
        fi
    fi

    echo "##############################"
    echo "## Resolved argument states ##"
    echo "##############################"
    echo ""
    echo "- Comparison operation: ${oper}"
    echo "- Scaling factor source: ${src_scl:-#N/A}"
    echo "- Normalization method: ${mth_nrm:-#N/A}"
    echo ""
    echo ""
fi


#  Activate environment and check that dependencies are in PATH ---------------
handle_env "${env_nam}" > /dev/null

check_program_path awk
check_program_path bamCompare
if ! ${slurm}; then check_program_path parallel; fi
check_program_path python
check_program_path samtools
if ${slurm}; then check_program_path sbatch; fi


#  Parse and validate vector elements -----------------------------------------
IFS=',' read -r -a arr_fil_num <<< "${fil_num}"
IFS=',' read -r -a arr_fil_den <<< "${fil_den}"
IFS=',' read -r -a arr_fil_stm <<< "${fil_stm}"

#  Validate file existence for file "numerator" and "denominator" arrays
check_array_files "fil_num" "${arr_fil_num[@]}"
check_array_files "fil_den" "${arr_fil_den[@]}"

#  Validate directory existence for outfile stem array
for fil in "${arr_fil_stm[@]}"; do
    dir="$(dirname "${fil}")"
    if [[ ! -d "${dir}" ]]; then
        echo \
            "Error: --fil_stm directory '${dir}' does not exist. Please" \
            "create it before proceeding." >&2
        exit_1
    fi
done

#  Validate argument compatibility with file extensions, populating an array of
#+ file extensions
unset arr_typ_fil && typeset -a arr_typ_fil
for i in "${!arr_fil_num[@]}"; do
    check_match \
        "File extensions for numerator and denominator at index ${i}" \
        "${arr_fil_num[i]##*.}" \
        "${arr_fil_den[i]##*.}" \
        "arr_typ_fil[i]"
done

#  Validate compatibility of --norm and --scl_pre with file types
for ext in "${arr_typ_fil[@]}"; do
    if [[ "${ext}" != "bam" ]]; then
        if [[ -n "${norm}" || -n "${scl_pre}" ]]; then
            echo \
                "Error: The --norm and --scl_pre arguments are only" \
                "applicable for BAM files. Detected file type '${ext}' is" \
                "not supported for these arguments. Please remove --norm" \
                "and --scl_pre if using non-BAM files (e.g., for deepTools" \
                "bigwigCompare)." >&2
            exit_1
        fi
    fi
done

#  Determine and validate sequencing types for numerator-denominator file pairs
unset arr_typ_seq && typeset -a arr_typ_seq
for i in "${!arr_fil_num[@]}"; do
    if [[ "${arr_typ_fil[i]}" == "bam" ]]; then
        typ_num=$(check_seq_type "${arr_fil_num[i]}")
        typ_den=$(check_seq_type "${arr_fil_den[i]}")

        check_match \
            "Sequencing type for file pair at index ${i}" \
            "${typ_num}" \
            "${typ_den}" \
            "arr_typ_seq[i]"
    else
        arr_typ_seq[i]="#N/A"
    fi
done

#  Handle scaling factors (--scl_fct)
if [[ -n "${scl_fct}" ]]; then
    IFS=',' read -r -a arr_scl_fct <<< "${scl_fct}"

    for s in "${arr_scl_fct[@]}"; do
        IFS=':' read -r num den <<< "${s}"
        check_flt_pos "${num}" "scl_fct_num"
        check_flt_pos "${den}" "scl_fct_den"
    done
else
    populate_array_empty arr_scl_fct "${#arr_fil_num[@]}"
fi

#  Handle user-defined fragment sizes (--usr_frg)
if [[ -n "${usr_frg}" ]]; then
    IFS=',' read -r -a arr_usr_frg <<< "${usr_frg}"
    for u in "${arr_usr_frg[@]}"; do check_flt_pos "${u}" "usr_frg"; done
else
    populate_array_empty arr_usr_frg "${#arr_fil_num[@]}"
fi

#  Ensure array lengths match across all vectors
check_arrays_lengths "arr_fil_num" arr_fil_num  "arr_fil_den" arr_fil_den
check_arrays_lengths "arr_fil_num" arr_fil_num  "arr_fil_stm" arr_fil_stm
check_arrays_lengths "arr_fil_num" arr_fil_num  "arr_typ_fil" arr_typ_fil
check_arrays_lengths "arr_fil_num" arr_fil_num  "arr_typ_seq" arr_typ_seq
check_arrays_lengths "arr_fil_num" arr_fil_num  "arr_scl_fct" arr_scl_fct
check_arrays_lengths "arr_fil_num" arr_fil_num  "arr_usr_frg" arr_usr_frg

#  Adjust max jobs for SLURM parallelization
if ${slurm}; then
    max_job=$(reset_max_job "${max_job}" "${#arr_fil_num[@]}")
fi

#  Debug output for parsed values
if ${verbose}; then
    echo "########################################"
    echo "## Parsed vectors for parallelization ##"
    echo "########################################"
    echo ""
    debug_array_contents \
        "arr_fil_num" "arr_fil_den" "arr_fil_stm" "arr_typ_fil" "arr_typ_seq" \
        "arr_scl_fct" "arr_usr_frg"
        
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
    echo "fil_num=${fil_num}"
    echo "fil_den=${fil_den}"
    echo "fil_stm=${fil_stm}"
    echo "typ_out=${typ_out}"
    echo "bin_siz=${bin_siz}"
    echo "region=${region:-#N/A}"
    echo "scl_fct=${scl_fct:-#N/A}"
    echo "norm=${norm:-#N/A}"
    echo "scl_pre=${scl_pre:-#N/A}"
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
    echo "arr_fil_num=( ${arr_fil_num[*]} )"
    echo "arr_fil_den=( ${arr_fil_den[*]} )"
    echo "arr_fil_stm=( ${arr_fil_stm[*]} )"
    echo "arr_typ_fil=( ${arr_typ_fil[*]} )"
    echo "arr_typ_seq=( ${arr_typ_seq[*]} )"
    echo "arr_scl_fct=( ${arr_scl_fct[*]} )"
    echo "arr_usr_frg=( ${arr_usr_frg[*]} )"
    echo ""
    echo ""
fi

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
        echo "    --array=1-${#arr_fil_num[@]}%${max_job} \\"
        echo "        ${scr_slm} \\"
        echo "             -v ${verbose} \\"
        echo "             -t ${threads} \\"
        echo "            -sn $(echo "${arr_fil_num[*]}" | tr ' ' ',') \\"
        echo "            -sd $(echo "${arr_fil_den[*]}" | tr ' ' ',') \\"
        echo "            -st $(echo "${arr_typ_fil[*]}" | tr ' ' ',') \\"
        echo "            -ss $(echo "${arr_fil_stm[*]}" | tr ' ' ',') \\"
        echo "            -to ${typ_out} \\"
        echo "            -op ${oper} \\"
        echo "            -bs ${bin_siz} \\"
        echo "             -r ${region:-#N/A} \\"
        echo "            -sf $(echo "${arr_scl_fct[*]}" | tr ' ' ',') \\"
        echo "            -no ${norm:-#N/A} \\"
        echo "            -sp ${scl_pre:-#N/A} \\"
        echo "            -ex ${exact} \\"
        echo "            -su $(echo "${arr_usr_frg[*]}" | tr ' ' ',') \\"
        echo "            -eo ${err_out} \\"
        echo "            -nj ${nam_job} \\"
        echo "            -en ${env_nam}"
        echo ""
        echo ""
        # echo "#########################################"
        # echo "## Contents of SLURM submission script ##"
        # echo "#########################################"
        # echo ""
        # echo "## ${scr_slm} ##"
        # echo ""
        # cat "${scr_slm}"
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
            --array=1-${#arr_fil_num[@]}%${max_job} \
                ${scr_slm} \
                     -v ${verbose} \
                     -t ${threads} \
                    -sn $(echo "${arr_fil_num[*]}" | tr ' ' ',') \
                    -sd $(echo "${arr_fil_den[*]}" | tr ' ' ',') \
                    -st $(echo "${arr_typ_fil[*]}" | tr ' ' ',') \
                    -ss $(echo "${arr_fil_stm[*]}" | tr ' ' ',') \
                    -to ${typ_out} \
                    -op ${oper} \
                    -bs ${bin_siz} \
                     -r ${region:-#N/A} \
                    -sf $(echo "${arr_scl_fct[*]}" | tr ' ' ',') \
                    -no ${norm:-#N/A} \
                    -sp ${scl_pre:-#N/A} \
                    -ex ${exact} \
                    -su $(echo "${arr_usr_frg[*]}" | tr ' ' ',') \
                    -eo ${err_out} \
                    -nj ${nam_job} \
                    -en ${env_nam}
    fi
fi
