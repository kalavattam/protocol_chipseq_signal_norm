#!/bin/bash

#  execute_deeptools_compare.sh
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
    source "${dir_fnc}/check_match.sh"
    source "${dir_fnc}/check_mut_excl_args.sh"
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/check_region.sh"
    source "${dir_fnc}/check_region_bam.sh"  #TODO
    source "${dir_fnc}/check_scl_fct.sh"
    source "${dir_fnc}/check_seq_type.sh"
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
    ## WARNING: If interactive=true, change values as needed ##
    #  Define base directories
    dir_bas="${HOME}/repos"
    dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
    dir_scr="${dir_rep}/scripts"
    dir_dat="${dir_rep}/data"
    dir_pro="${dir_dat}/processed"

    #  Define alignment settings
    aligner="bowtie2"
    a_type="global"
    req_flg=true
    flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
    mapq=1
    det_bam="flag-${flg}_mapq-${mapq}"

    #  Define coverage and normalization settings
    det_cov="${aligner}_${a_type}_${det_bam}"
    dir_cov="${dir_pro}/compute_coverage/${det_cov}"
    # cov_nrm="depth"  # "norm"

    #  Define variables for table of samples, corresponding alpha values
    dir_tbl="${dir_cov}/alpha/tables"
    fil_tbl="ChIP_WT_G1-G2M-Q_Hho1-Hmo1_alpha_6nd.tsv"
    pth_tbl="${dir_tbl}/${fil_tbl}"

    # #  Define file type and location
    # typ_fil="bam"  # "bigwig"
    # case "${typ_fil}" in
    #     bedgraph|bdg|bg|bigwig|bw)
    #         dir_bdg="${dir_cov}/${cov_nrm}/tracks"
    #         dir_bwg="${dir_cov}/${cov_nrm}/tracks"
    #         ;;
    #     bam)
    #         dir_aln="${dir_pro}/align_${aligner}_${a_type}"
    #         dir_bam="${dir_aln}/${det_bam}/sc"
    #         ;;
    # esac

    # #  Define file search settings
    # case "${typ_fil}" in
    #     bedgraph|bdg|bg)
    #         dir_fnd="${dir_bdg}"
    #         pattern="*.bg"
    #         ;;
    #     bigwig|bw)
    #         dir_fnd="${dir_bwg}"
    #         pattern="*.bw"
    #         ;;
    #     bam)
    #         dir_fnd="${dir_bam}"
    #         pattern="*.bam"
    #         ;;
    # esac
    # exclude="*Brn1*"

    #  Define track outfile track location
    dir_trk="${dir_cov}/alpha/tracks"
    # dir_trk="${dir_cov}/log2/${cov_nrm}/tracks"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=false
    threads=8
    ser_num="$(awk 'NR > 1 { print $1 }' "${pth_tbl}" | paste -sd ',' -)"
    # ser_num="$(
    #     bash "${dir_scr}/find_files.sh" \
    #         --dir_fnd "${dir_fnd}" \
    #         --pattern "${pattern}" \
    #         --include "IP*" \
    #         --exclude "${exclude}"
    # )"
    ser_den="$(awk 'NR > 1 { print $2 }' "${pth_tbl}" | paste -sd ',' -)"
    # ser_den="$(
    #     bash "${dir_scr}/find_files.sh" \
    #         --dir_fnd "${dir_fnd}" \
    #         --pattern "${pattern}" \
    #         --include "in*" \
    #         --exclude "${exclude}"
    # )"
    ser_stm=$(
        echo "${ser_num}" \
            | sed -e "s|$(dirname "${ser_num%%,*}")|${dir_trk}|g" \
                  -e 's:/IP:/siq:g' \
                  -e 's:\.\(bigwig\|bw\|bedgraph\|bdg\|bg\|bam\)::g'
    )
    # ser_stm=$(
    #     echo "${ser_num}" \
    #         | sed -e "s|${dir_fnd}|${dir_trk}|g" \
    #               -e 's:/IP:/IP-in:g' \
    #               -e 's:\.\(bigwig\|bw\|bedgraph\|bdg\|bg\|bam\)::g'
    # )
    typ_out="bedgraph"  # "bigwig"
    siz_bin=10  # 1
    region=""
    oper="ratio"
    scl_fct=$(awk 'NR > 1 { print $3 }' "${pth_tbl}" | paste -sd ',' -)  # ""
    typ_cvg="none"
    scl_pre="none"  # "depth"
    exact=false
    usr_frg=""
    err_out="${dir_trk}/logs"
    nam_job="run_deeptools_compare"
    slurm=true
    max_job=6  # 4
    time="0:30:00"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_analyze"
scr_sub="${dir_scr}/submit_deeptools_compare_slurm.sh"  #TODO: Lose "_slurm"
denom=4

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=1
ser_num=""
ser_den=""
ser_stm=""
typ_out="bedgraph"  # "bigwig"
siz_bin=10
siz_gen=12157105
region=""
oper="ratio"  # "log2"
scl_fct="none"
typ_cvg="none"
scl_pre="none"
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
    [--verbose] [--dry_run] --threads <int> --ser_num <str> --ser_den <str>
    --ser_stm <str> --typ_out <str> --siz_bin <int> [--region <str>]
    --scl_fct <flt> --typ_cvg <str> --scl_pre <str> [--exact] [--usr_frg <flt>]
    --err_out <str> --nam_job <str> [--slurm] [--max_job <int>] [--time <str>]

Description:
  This driver script, 'execute_deeptools_compare.sh', automates the comparison
  of genomic signal tracks, such as ChIP-seq IP/input signal, using deepTools
  utilities 'bamCompare' and 'bigwigCompare'. The script supports BAM,
  BEDGRAPH, or BIGWIG input files, generating normalized and/or scaled
  comparison tracks in either BEDGRAPH or BIGWIG formats.

Arguments:
   -h, --help     Print this help message and exit.
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Run script in 'dry-run mode' (optional).
   -t, --threads  Number of threads to use (default: ${threads}).
  -sn, --ser_num  Comma-separated serialized string of BAM, BIGWIG, or BEDGRAPH
                  infiles to be used as the first file (e.g., numerator) in the
                  comparisons.
  -sd, --ser_den  Comma-separated serialized string of BAM, BIGWIG, or BEDGRAPH
                  infiles to be used as the second file (e.g., denominator) in
                  the comparisons.
  -ss, --ser_stm  Comma-separated serialized string of outstems without
                  extensions.
  -to, --typ_out  Format of comparison signal track outfiles: 'bedgraph' or
                  'bigwig' (default: "${typ_out}").
  -sb, --siz_bin  Bin size for comparison calculation(s) in base pairs
                  (default: ${siz_bin}).
  -sg, --siz_gen  Effective genome size of model organism (required if
                  '--typ_cvg rpgc', ignored if not; default: ${siz_gen}).
   -r, --region   Specify a genomic region to limit the operation (optional).
                  The format can be either: 'chr' (entire chromosome) or
                  'chr:start-stop' (specific range within a chromosome).
  -op, --oper     Operation to perform for comparison(s) of first/numerator and
                  second/denominator BAM infiles: 'log2', 'ratio', 'subtract',
                  'add', 'mean', 'reciprocal_ratio', 'first', 'second'
                  (default: "${oper}").
  -sf, --scl_fct  Comma-separated serialized string of precomputed scaling
                  factor pairs for file comparisons. Each pair should be in the
                  format  'first:second' (e.g., 'num:den'), where 'first'
                  ('num') is the scaling factor for the first file and 'second'
                  ('den') is for the second. If only 'num' is provided (e.g.,
                  '0.5'), it defaults to 'num:1' ('0.5:1'). Multiple pairs
                  should be comma-separated. Use 'none' to skip precomputed
                  scaling factors (default: "${scl_fct}").
  -tv, --typ_cvg  Specify a normalization method to use when computing
                  comparison(s); available options are 'none', 'raw', 'rpkm',
                  'fpkm', 'cpm', 'bpm', or 'rpgc' (default: "${typ_cvg}").
  -sp, --scl_pre  Predefined scaling method for compensating sequencing depth
                  differences in comparison operations: 'depth', 'ses', or
                  'none' (default: "${scl_pre}").
  -ex, --exact    Compute scaling factors using all alignments instead of a
                  sampled subset. Only applicable if '--typ_cvg <str>' is
                  specified and ignored if '--typ_cvg none', '--typ_cvg #N/A',
                  or '--typ_cvg' is not invoked. Significantly increases the
                  time required for comparison computation.
  -uf, --usr_frg  Comma-separated string vector of fragment lengths to use
                  instead of read lengths (for alignments of single-end reads)
                  or template lengths (for alignments of paired-end reads;
                  optional). Must match the number of infiles and outstems via
                  '--ser_num <str>', '--ser_den <str>', and '--ser_stm <str>'.
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (default: "\$(dirname "\${ser_stm%%[,;]*}")/err_out").
  -nj, --nam_job  Prefix for the names of jobs (default: "${nam_job}").
  -sl, --slurm    Submit jobs to the SLURM scheduler; otherwise, run them via
                  GNU Parallel (optional).
  -mj, --max_job  The maximum number of jobs to run at one time (required if
                  '--slurm' is specified, ignored if not; default: ${max_job}).
  -tm, --time     The length of time, in 'h:mm:ss' format, for the SLURM job
                  (required if '--slurm' is specified, ignored if not; default:
                  ${time}).

Dependencies:
  - Programs
    + awk
    + Bash or Zsh
    + deepTools 'bamCompare' and 'bigwigCompare'
    + GNU Parallel (if not using SLURM and 'threads' > 1)
    + Python
    + Samtools
    + SLURM (if using '--slurm')
  - Functions
    + #TODO

Notes:
  - When the '--slurm' flag is used, jobs are parallelized via SLURM array
    tasks; otherwise, jobs are parallelized with GNU Parallel.
  - BAM infiles must be coordinate-sorted.
  - On '--scl_fct': For example, assigning '0.7:1,0.3:0.6' scales the first
    file in the first pair (i.e., the numerator file for pair #1) by 0.7 and
    the second file in the first pair (i.e., the denominator file for pair #1)
    by 1.0, and the first file in the second pair (i.e., the numerator file for
    pair #2) by 0.3 and the second file in the second pair (i.e., the
    denominator file for pair #2) by 0.6.
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
            -sn|--ser_num) ser_num="${2}";   shift 2 ;;
            -sd|--ser_den) ser_den="${2}";   shift 2 ;;
            -ss|--ser_stm) ser_stm="${2}";   shift 2 ;;
            -to|--typ_out) typ_out="${2,,}"; shift 2 ;;
            -sb|--siz_bin) siz_bin="${2}";   shift 2 ;;
            -sg|--siz_gen) siz_gen="${2}";   shift 2 ;;
             -r|--region)  region="${2}";    shift 2 ;;
            -op|--oper)    oper="${2,,}";    shift 2 ;;
            -sf|--scl_fct) scl_fct="${2,,}"; shift 2 ;;
            -no|--typ_cvg) typ_cvg="${2,,}"; shift 2 ;;
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

check_supplied_arg -a "${scr_sub}" -n "scr_sub"
check_exists_file_dir "f" "${scr_sub}" "scr_sub"

check_supplied_arg -a "${denom}" -n "denom"

check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

check_supplied_arg -a "${ser_num}" -n "ser_num"
check_exists_file_dir "d" "$(dirname "${ser_num%%[,;]*}")" "ser_num"

check_supplied_arg -a "${ser_den}" -n "ser_den"
check_exists_file_dir "d" "$(dirname "${ser_den%%[,;]*}")" "ser_den"

check_supplied_arg -a "${ser_stm}" -n "ser_stm"
check_exists_file_dir "d" "$(dirname "${ser_stm%%[,;]*}")" "ser_stm"

check_supplied_arg -a "${typ_out}" -n "typ_out"
case "${typ_out}" in
    bedgraph|bdg|bg) typ_out="bedgraph" ;;
    bigwig|bw)       typ_out="bigwig"   ;;
    *)
        echo_error \
            "Selection associated with '--typ_out' is not valid:" \
            "'${typ_out}'. Selection must be 'bedgraph' or 'bigwig'."
        exit_1
        ;;
esac

check_supplied_arg -a "${siz_bin}" -n "siz_bin"
check_int_pos "${siz_bin}" "siz_bin"

check_supplied_arg -a "${siz_gen}" -n "siz_gen"
check_int_pos "${siz_gen}" "siz_gen"

if [[ -n "${region}" ]]; then check_region "${region}"; fi

case "${oper}" in
    log*)        oper="log2"             ;;
    1*|first)    oper="first"            ;;
    2*|second)   oper="second"           ;;
    av*|mean|mn) oper="mean"             ;;
    rat*)        oper="ratio"            ;;
    sub*)        oper="subtract"         ;;
    add*)        oper="add"              ;;
    recip*)      oper="reciprocal_ratio" ;;
    *)
        echo_error \
            "Invalid operation '${oper}'. Supported operations are 'log2'," \
            "'ratio', 'subtract', 'add', 'mean', 'reciprocal_ratio'," \
            "'first', or 'second'." >&2
        exit_1
        ;;
esac

#  Process '--scl_fct'
if [[ -z "${scl_fct}" || "${scl_fct}" == "none" ]]; then
    scl_fct="None"
else
    #  Validate and format 'scl_fct'
    if ! scl_fct_pro=$(check_scl_fct "${scl_fct}" 2> /dev/null); then
        echo "Error: Invalid value for '--scl_fct': '${scl_fct}'." >&2
        exit 1
    fi
    scl_fct="${scl_fct_pro}"
    unset scl_fct_pro
fi

#  Process '--typ_cvg'
if [[ -z "${typ_cvg}" ]]; then typ_cvg="none"; fi
case "${typ_cvg}" in
    none|raw)  typ_cvg="None" ;;
    rpkm|fpkm) typ_cvg="RPKM" ;;
    cpm)       typ_cvg="CPM"  ;;
    bpm)       typ_cvg="BPM"  ;;
    rpgc)      typ_cvg="RPGC" ;;
    *)
        echo_error \
            "Invalid '--typ_cvg' selection: '${typ_cvg}'. Valid options are" \
            "'none', 'raw', 'rpkm', 'fpkm', 'cpm', 'bpm', or 'rpgc'."
        exit_1
        ;;
esac

#  Process '--scl_pre'
if [[ -z "${scl_pre}" ]]; then scl_pre="none"; fi
case "${scl_pre}" in
    depth) scl_pre="readCount" ;;
    ses)   scl_pre="SES"       ;;
    none)
        # if [[ "${typ_cvg}" == "None" ]]; then
        #     echo_warning \
        #         "Specified both '--scl_fct none' and '--typ_cvg none'. The" \
        #         "'bamCompare' or 'bigwigCompare' operation ('--oper" \
        #         "${oper}')" "will proceed without scaling, allowing" \
        #         "sequencing depth to potentially affect results."
        # fi
        scl_pre="None"
        ;;
    *)
        echo_error \
            "Invalid '--scl_pre' selection: '${scl_pre}'. Valid options are" \
            "'depth', 'ses', or 'none'."
        exit_1
        ;;
esac

#  Issue a warning for the dual selection of non-"none" '--scl_fct' and
#+ '--typ_cvg'
if [[ "${scl_fct}" != "None" && "${typ_cvg}" != "None" ]]; then
    echo_warning \
        "Using '--scl_fct ${scl_fct}' with '--typ_cvg ${typ_cvg}' (i.e., when" \
        "neither '--typ_cvg' is nor '--scl_fct' are 'none') is not" \
        "recommended. Options associated with '--typ_cvg' are designed to" \
        "adjust for differences in sequenced read alignment counts (i.e.," \
        "sequencing depth after alignment and processing), while values" \
        "supplied to '--scl_fct <str>' typically account for compositional" \
        "differences between libraries, such as those derived from exogenous" \
        "spike-in DNA or scaling factors computed by tools such as siQ-ChIP," \
        "DESeq2, or edgeR. Combining these options can skew y-axis values," \
        "leading to misleading or unreliable semiquantitative results."
fi

if [[ -n "${usr_frg}" ]]; then check_str_delim "usr_frg" "${usr_frg}"; fi

if [[ -n "${err_out}" ]]; then
    check_exists_file_dir "d" "${err_out}" "err_out"
elif [[ -z "${err_out}" ]]; then
    err_out="$(dirname "${ser_stm%%[,;]*}")/err_out"
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
if ${verbose}; then
    echo "##############################"
    echo "## Resolved argument states ##"
    echo "##############################"
    echo ""
    echo "- Comparison operation: '--oper ${oper}'"
    echo "- Scaling factor source(s):"
    echo "    + '--scl_fct ${scl_fct}'"
    echo "    + '--typ_cvg ${typ_cvg}'"
    echo "    + '--scl_pre ${scl_pre}'"
    echo ""
    echo ""
fi


#  Activate environment and check that dependencies are in PATH ---------------
handle_env "${env_nam}" > /dev/null

check_program_path awk
check_program_path bamCompare
check_program_path bigwigCompare
if ! ${slurm}; then check_program_path parallel; fi
check_program_path python
check_program_path samtools
if ${slurm}; then check_program_path sbatch; fi


#  Parse and validate vector elements -----------------------------------------
IFS=',' read -r -a arr_fil_num <<< "${ser_num}"
IFS=',' read -r -a arr_fil_den <<< "${ser_den}"
IFS=',' read -r -a arr_fil_stm <<< "${ser_stm}"

#  Validate file existence for file "numerator" and "denominator" arrays
check_array_files "ser_num" "${arr_fil_num[@]}"
check_array_files "ser_den" "${arr_fil_den[@]}"

#  Validate directory existence for outfile stem array
for fil in "${arr_fil_stm[@]}"; do
    # echo "${fil}"
    dir="$(dirname "${fil}")"
    if [[ ! -d "${dir}" ]]; then
        echo \
            "Error: --ser_stm directory '${dir}' does not exist. Please" \
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

#  Validate compatibility of --typ_cvg and --scl_pre with file types
for ext in "${arr_typ_fil[@]}"; do
    if [[ "${ext}" != "bam" ]]; then
        if [[ "${typ_cvg}" != "None" || "${scl_pre}" != "None" ]]; then
            echo_error \
                "The '--typ_cvg' and '--scl_pre' arguments are only" \
                "applicable for BAM files; you specified '--typ_cvg" \
                "${typ_cvg}' and '--scl_pre ${scl_pre}'. Detected file type" \
                "'${ext}' is not supported for these arguments. Please" \
                "specify '--typ_cvg none' and '--scl_pre none' if using" \
                "BEDGRAPH or BIGWIG files."
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

#  Handle scaling factors ('--scl_fct')
if [[ -n "${scl_fct}" && "${scl_fct}" != "None" ]]; then
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
    echo "ser_num=${ser_num}"
    echo "ser_den=${ser_den}"
    echo "ser_stm=${ser_stm}"
    echo "typ_out=${typ_out}"
    echo "siz_bin=${siz_bin}"
    echo "region=${region:-#N/A}"
    echo "scl_fct=${scl_fct}"
    echo "typ_cvg=${typ_cvg}"
    echo "scl_pre=${scl_pre}"
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
        echo "        ${scr_sub} \\"
        echo "             -v ${verbose} \\"
        echo "             -t ${threads} \\"
        echo "            -sn $(echo "${arr_fil_num[*]}" | tr ' ' ',') \\"
        echo "            -sd $(echo "${arr_fil_den[*]}" | tr ' ' ',') \\"
        echo "            -st $(echo "${arr_typ_fil[*]}" | tr ' ' ',') \\"
        echo "            -ss $(echo "${arr_fil_stm[*]}" | tr ' ' ',') \\"
        echo "            -to ${typ_out} \\"
        echo "            -op ${oper} \\"
        echo "            -sb ${siz_bin} \\"
        echo "            -sg ${siz_gen} \\"
        echo "             -r ${region:-#N/A} \\"
        echo "            -sf $(echo "${arr_scl_fct[*]}" | tr ' ' ',') \\"
        echo "            -tv ${typ_cvg} \\"
        echo "            -sp ${scl_pre} \\"
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
        # echo "## ${scr_sub} ##"
        # echo ""
        # cat "${scr_sub}"
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
                ${scr_sub} \
                    $(if ${verbose}; then echo "-v"; fi) \
                     -t ${threads} \
                    -sn $(echo "${arr_fil_num[*]}" | tr ' ' ',') \
                    -sd $(echo "${arr_fil_den[*]}" | tr ' ' ',') \
                    -st $(echo "${arr_typ_fil[*]}" | tr ' ' ',') \
                    -ss $(echo "${arr_fil_stm[*]}" | tr ' ' ',') \
                    -to ${typ_out} \
                    -op ${oper} \
                    -sb ${siz_bin} \
                    $(
                        if [[ -n "${siz_gen}" ]]; then
                            echo "-sg ${siz_gen}"
                        fi
                    ) \
                     -r ${region:-#N/A} \
                    -sf $(echo "${arr_scl_fct[*]}" | tr ' ' ',') \
                    -tv ${typ_cvg} \
                    -sp ${scl_pre} \
                    $(if ${exact}; then echo "-ex"; fi) \
                    -su $(echo "${arr_usr_frg[*]}" | tr ' ' ',') \
                    -eo ${err_out} \
                    -nj ${nam_job} \
                    -en ${env_nam}
    fi
fi

#TODO: Implement GNU Parallel and serial job submissions
#TODO: Correct siz_gen echo statement
