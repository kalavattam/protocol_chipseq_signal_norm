#!/bin/bash

#  execute_compute_signal.sh
#  KA


#  Run script in interactive mode (true) or command-line mode (false)
interactive=false

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive}; then set -euo pipefail; fi

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
    check_array_files \
    check_arrays_lengths \
    check_exists_file_dir \
    check_flt_pos \
    check_format_time \
    check_int_pos \
    check_program_path \
    check_str_delim \
    check_supplied_arg \
    debug_array_contents \
    echo_error \
    echo_warning \
    exit_0 \
    exit_1 \
    handle_env \
    populate_array_empty \
    print_parallel_info \
    reset_max_job \
    set_params_parallel \
    summarize_sig_nrm
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
    dir_aln="${dir_pro}/align_fastqs"
    dir_det="${dir_aln}/${str_det}"
    dir_bam="${dir_det}/sc"

    #  Set output directories
    dir_sig="${dir_pro}/compute_signal/${str_det}"
    dir_trk="${dir_sig}/norm"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=true
    threads=6
    infiles="$(  ## WARNING: Change search parameters as needed ##
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_bam}" \
            --pattern "*.bam"
    )"
    dir_out="${dir_trk}"
    typ_out="bdg.gz"
    siz_bin=30
    scl_fct=""
    typ_sig="$(basename "${dir_trk}")"
    usr_frg=""
    rnd=24
    err_out="${dir_trk}/logs"
    nam_job="compute_signal"
    max_job=2
    slurm=false
    time="0:30:00"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_protocol"
scr_sub="${dir_scr}/submit_compute_signal.sh"
scr_sig="${dir_scr}/compute_signal.py"
par_job=""

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=4
infiles=""
dir_out=""
typ_out="bdg.gz"
siz_bin=10
scl_fct=""
typ_sig="norm"
usr_frg=""
rnd=24
err_out=""
nam_job="compute_signal"
max_job=6
slurm=false
time="0:30:00"

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  execute_compute_signal.sh
    [--verbose] [--dry_run] --threads <int> --infiles <str> --dir_out <str>
    --typ_out <str> --siz_bin <int> [--scl_fct <flt>] --typ_sig <str>
    [--usr_frg <flt>] --rnd <int> --err_out <str> --nam_job <str> --max_job <int>
    [--slurm] [--time <str>]

Description:
  The driver script 'execute_compute_signal.sh' automates the calculation of
  signal tracks from BAM files. It supports various normalization methods,
  including user-supplied scaling factors, normalized coverage (Dickson et al.,
  Sci Rep 2023), and raw (unadjusted) coverage.

  The script supports parallel execution via SLURM or GNU Parallel, or can run
  serially.

Arguments:
   -h, --help     Print this help message and exit.
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Run script in 'dry-run mode' (optional).
   -t, --threads  Number of threads to use (default: ${threads}).
   -i, --infiles  Comma-separated string of coordinate-sorted BAM infiles.
  -do, --dir_out  Directory to write coverage signal track outfiles.
  -to, --typ_out  Format of coverage signal track outfiles: 'bedgraph',
                  'bedgraph.gz', 'bdg', 'bdg.gz', 'bg', 'bg.gz', 'bed', or
                  'bed.gz' (default: '${typ_out}').
                    - 'bedgraph', 'bdg', 'bg': BEDGRAPH format (coverage).
                    - 'bedgraph.gz', 'bdg.gz', 'bg.gz': gzip-compressed
                       BEDGRAPH format.
                    - 'bed', 'bed.gz': BED-like format (fragment coordinates,
                      not coverage signal). Using this option will cause
                      arguments such as '--siz_bin' and '--typ_sig' to be
                      ignored.
  -sb, --siz_bin  Bin size for coverage calculation in base pairs (default:
                  ${siz_bin}).
  -sf, --scl_fct  Comma-separated string of scaling factors to apply to
                  coverage (optional). Must match the number of infiles via
                  '--infiles'.
  -tv, --typ_sig  Specify coverage calculation type (default: '${typ_sig}').
                  Options:
                    - 'unadj', 'unadjusted': Compute unadjusted signal.
                    - 'len', 'len_frag': Adjust signal by fragment length.
                    - 'norm', 'normalized': Per Dickson et al., Sci Rep 2023
                      (PMID: 37160995), adjust signal by fragment length
                      and total fragments/unity, generating "normalized
                      coverage".
  -uf, --usr_frg  Comma-separated string of fragment lengths to use instead of
                  read lengths (single-end alignments) or template lengths
                  (paired-end alignments; optional). Must match the number of
                  infiles via '--infiles'.
   -r, --rnd      Number of decimal places for rounding signal (coverage score)
                  values (default: '${rnd}').
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (default: "\${dir_out}/err_out").
  -nj, --nam_job  Prefix for the names of jobs (default: '${nam_job}').
  -mj, --max_job  Maximum number of jobs to run concurrently (default: '${max_job}').
                    - If '--slurm' is specified, controls SLURM array tasks.
                    - If '--slurm' is not specified:
                      + If 'max_job' is greater than 1, jobs run in parallel
                        via GNU Parallel.
                      + If 'max_job' is 1, jobs run sequentially (serial mode).
  -sl, --slurm    Submit jobs to the SLURM scheduler (optional; otherwise, if
                  'threads' > 1, run them via GNU Parallel, and if 'threads' is
                  1, run them in serial.)
  -tm, --time     The length of time, in 'h:mm:ss' format, for the SLURM job
                  (required if '--slurm' is specified, ignored if not; default:
                  '${time}').

Dependencies:
  - Bash or Zsh
  - GNU Parallel (if not using SLURM and threads > 1)
  - Python
  - SLURM (if using --slurm)

Notes:
  - When the '--slurm' flag is used, jobs are parallelized via SLURM array
    tasks; otherwise... #TODO
  - BAM infiles must be coordinate-sorted.
  - Outfile names are derived from BAM infiles and the value(s) associated with
    '--typ_out'.

Examples:
  \`\`\`
  #  Example 1: Run with GNU Parallel
  bash "\${HOME}/path/to/scripts/execute_compute_signal.sh"
      --threads 8
      --infiles "\${HOME}/path/to/sample1.bam,\${HOME}/path/to/sample2.bam"
      --dir_out "\${HOME}/path/to/outfiles"
      --typ_out "bdg.gz"
      --siz_bin 50
      --typ_sig "norm"
      --err_out "\${HOME}/path/to/outfiles/logs"
      --nam_job "coverage_calculation"

  #  Example 2: #TODO
  \`\`\`
EOM
)

#  Parse arguments
if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    exit_0
fi

if ${interactive:-false}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -v|--verbose) verbose=true;     shift 1 ;;
            -dr|--dry_run) dry_run=true;     shift 1 ;;
             -t|--threads) threads="${2}";   shift 2 ;;
             -i|--infiles) infiles="${2}";   shift 2 ;;
            -do|--dir_out) dir_out="${2}";   shift 2 ;;
            -to|--typ_out)
                typ_out="$(echo "${typ_out}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;
            -sb|--siz_bin) siz_bin="${2}";   shift 2 ;;
            -sf|--scl_fct) scl_fct="${2}";   shift 2 ;;
            -tv|--typ_sig) typ_sig="${2}";   shift 2 ;;
            -uf|--usr_frg) usr_frg="${2}";   shift 2 ;;
             -r|--rnd)     rnd="${2}";       shift 2 ;;
            -eo|--err_out) err_out="${2}";   shift 2 ;;
            -nj|--nam_job) nam_job="${2}";   shift 2 ;;
            -mj|--max_job) max_job="${2}";   shift 2 ;;
            -sl|--slurm)   slurm=true;       shift 1 ;;
            -tm|--time)    time="${2}";      shift 2 ;;
            *)
                echo "## Unknown parameter passed: '${1}' ##" >&2
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

check_supplied_arg -a "${scr_sig}" -n "scr_sig"
check_exists_file_dir "f" "${scr_sig}" "scr_sig"

check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

check_supplied_arg -a "${infiles}" -n "infiles"
check_exists_file_dir "d" "$(dirname "${infiles%%[,;]*}")" "infiles"
check_str_delim "infiles" "${infiles}"

check_supplied_arg -a "${dir_out}" -n "dir_out"
check_exists_file_dir "d" "${dir_out}" "dir_out"

check_supplied_arg -a "${typ_out}" -n "typ_out"
case "${typ_out}" in
    bedgraph|bedgraph.gz|bdg|bdg.gz|bg|bg.gz|bed|bed.gz) : ;;
    *)
        echo_error \
            "Selection associated with '--typ_out' is not valid:" \
            "'${typ_out}'. Expected 'bedgraph', 'bedgraph.gz', 'bdg'," \
            "'bdg.gz', 'bg', 'bg.gz', 'bed', or 'bed.gz'."
        exit_1
        ;;
esac

check_supplied_arg -a "${siz_bin}" -n "siz_bin"
check_int_pos "${siz_bin}" "siz_bin"

if [[ -n "${scl_fct}" ]]; then check_str_delim "scl_fct" "${scl_fct}"; fi

if [[ -z "${typ_sig}" ]]; then
    typ_sig="unadj"
    echo_warning \
        "No coverage normalization provided. Defaulting to '${typ_sig}' (i.e.,
        'non-normalized') coverage ('--typ_sig ${typ_sig}')."
else
    case "${typ_sig}" in
        raw|unadj|unadjusted) : ;;  # Valid options for unadjusted coverage
        len|len_frag)         : ;;  # Valid options for frag.-length norm.
        norm|normalized)      : ;;  # Valid options for normalized coverage
        *)
            echo_error \
                "Invalid value for '--typ_sig': '${typ_sig}'. Expected" \
                "'raw', 'unadj', 'unadjusted', 'len', 'len_frag', 'norm', or" \
                "'normalized'."
            exit_1
            ;;
    esac
fi

if [[ -n "${usr_frg}" ]]; then check_str_delim "usr_frg" "${usr_frg}"; fi

check_supplied_arg -a "${rnd}" -n "rnd"
check_int_pos "${rnd}" "rnd"

if [[ -n "${err_out}" ]]; then
    check_exists_file_dir "d" "${err_out}" "err_out"
elif [[ -z "${err_out}" ]]; then
    err_out="${dir_out}/err_out"
    check_exists_file_dir "d" "${err_out}" "err_out"
fi

check_supplied_arg -a "${nam_job}" -n "nam_job"


#  Parse and validate input vector elements -----------------------------------
IFS=',' read -r -a arr_infile <<< "${infiles}"
check_array_files "infiles" "${arr_infile[@]}"

unset arr_outfile && typeset -a arr_outfile
for i in "${arr_infile[@]}"; do
    arr_outfile+=( "${dir_out}/$(basename "${i}" .bam).${typ_out}" )
done

if [[ -z "${scl_fct}" ]]; then
    unset arr_scl_fct && typeset -a arr_scl_fct
    populate_array_empty arr_scl_fct "${#arr_infile[@]}"
else
    IFS=',' read -r -a arr_scl_fct <<< "${scl_fct}"
fi

for s in "${arr_scl_fct[@]}"; do
    if [[ "${s}" != "#N/A" ]]; then check_flt_pos "${s}" "scl_fct"; fi
done

if [[ -z "${usr_frg}" ]]; then
    unset arr_usr_frg && typeset -a arr_usr_frg
    populate_array_empty arr_usr_frg "${#arr_infile[@]}"
else
    IFS=',' read -r -a arr_usr_frg <<< "${usr_frg}"
fi

for u in "${arr_usr_frg[@]}"; do
    if [[ "${u}" != "#N/A" ]]; then check_flt_pos "${u}" "usr_frg"; fi
done

check_arrays_lengths "arr_outfile" "arr_infile"
check_arrays_lengths "arr_scl_fct" "arr_infile"
check_arrays_lengths "arr_usr_frg" "arr_infile"


#  Parse job execution parameters ---------------------------------------------
if ${slurm:-false}; then
    check_supplied_arg -a "${max_job}" -n "max_job"
    check_int_pos "${max_job}" "max_job"

    max_job=$(reset_max_job "${max_job}" "${#arr_infile[@]}")
    
    check_supplied_arg -a "${time}" -n "time"
    check_format_time "${time}"
else
    IFS=';' read -r threads par_job < <(
        set_params_parallel "${threads}" "${max_job}" "${par_job}"
    )
    unset IFS max_job time

    check_supplied_arg -a "${par_job}" -n "par_job"
    check_int_pos "${par_job}" "par_job"
fi

#  Debug parallelization information and summary output of resolved states
print_parallel_info \
    "${slurm}" "${max_job:-#N/A}" "${par_job}" "${threads}" \
    "arr_infile" "arr_outfile" "arr_scl_fct" "arr_usr_frg"

summarize_sig_nrm "${typ_sig}" "${scl_fct}"


#  Activate environment and check that dependencies are in PATH ---------------
handle_env "${env_nam}" > /dev/null

check_program_path python

if ${slurm:-false}; then
    check_program_path sbatch
elif [[ "${threads}" -gt 1 ]]; then
    check_program_path parallel
fi


#  Do the main work ===========================================================
IFS=" "
if ${verbose}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "env_nam=${env_nam}"
    echo "scr_sub=${scr_sub}"
    echo "scr_sig=${scr_sig}"
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
    echo "infiles=${infiles}"
    echo "dir_out=${dir_out}"
    echo "typ_out=${typ_out}"
    echo "siz_bin=${siz_bin}"
    echo "scl_fct=${scl_fct:-#N/A}"
    echo "typ_sig=${typ_sig}"
    echo "usr_frg=${usr_frg:-#N/A}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo "max_job=${max_job:-#N/A}"
    echo "slurm=${slurm}"
    echo "time=${time:-#N/A}"
    echo ""
    echo ""
    echo "###################################"
    echo "## Arrays derived from variables ##"
    echo "###################################"
    echo ""
    echo "arr_infile=( ${arr_infile[*]} )"
    echo ""
    echo "arr_outfile=( ${arr_outfile[*]} )"
    echo ""
    echo "arr_scl_fct=( ${arr_scl_fct[*]} )"
    echo ""
    echo "arr_usr_frg=( ${arr_usr_frg[*]} )"
    echo ""
    echo ""
fi
unset IFS

infiles=$(echo "${arr_infile[*]}"  | tr ' ' ',')
outfiles=$(echo "${arr_outfile[*]}" | tr ' ' ',')
scl_fct=$(echo "${arr_scl_fct[*]}" | tr ' ' ',')
usr_frg=$(echo "${arr_usr_frg[*]}" | tr ' ' ',')

if ${verbose}; then
    echo "###############################################################"
    echo "## Variable assignments re- or newly constructed from arrays ##"
    echo "###############################################################"
    echo ""
    echo "infiles=\"${infiles}\""
    echo ""
    echo "outfiles=\"${outfiles}\""
    echo ""
    echo "scl_fct=\"${scl_fct}\""
    echo ""
    echo "usr_frg=\"${usr_frg}\""
    echo ""
    echo ""
fi

# shellcheck disable=SC2016,SC2090
if ${slurm:-false}; then
    #  SLURM execution
    if ${dry_run} || ${verbose}; then
        echo "######################"
        echo "## Call to 'sbatch' ##"
        echo "######################"
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
        echo "            ${scr_sig} \\"
        echo "            ${threads} \\"
        echo "            ${infiles} \\"
        echo "            ${outfiles} \\"
        echo "            ${siz_bin} \\"
        echo "            ${typ_sig} \\"
        echo "            ${scl_fct} \\"
        echo "            ${usr_frg} \\"
        echo "            ${rnd} \\"
        echo "            ${err_out} \\"
        echo "            ${nam_job}"
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
            --array=1-${#arr_infile[@]}%${max_job} \
                ${scr_sub} \
                    ${env_nam} \
                    ${scr_sig} \
                    ${threads} \
                    ${infiles} \
                    ${outfiles} \
                    ${siz_bin} \
                    ${typ_sig} \
                    ${scl_fct} \
                    ${usr_frg} \
                    ${rnd} \
                    ${err_out} \
                    ${nam_job}
    fi
else
    #  GNU Parallel execution
    if [[ ${threads} -gt 1 ]]; then
        #  Create and populate GNU Parallel configuration file
        config="${err_out}/${nam_job}.config_parallel.txt"

        if [[ -f "${config}" ]]; then rm "${config}"; fi
        
        if ! \
            touch "${config}"
        then
            echo_error "Failed to create GNU Parallel configuration file."
            exit_1
        fi

        if ! \
            for idx in "${!arr_infile[@]}"; do
                echo \
                    "${env_nam}" \
                    "${scr_sig}" \
                    "${threads}" \
                    "${arr_infile[idx]}" \
                    "${arr_outfile[idx]}" \
                    "${siz_bin}" \
                    "${typ_sig}" \
                    "${arr_scl_fct[idx]}" \
                    "${arr_usr_frg[idx]}" \
                    "${rnd}" \
                    "${err_out}" \
                    "${nam_job}" \
                    "$(basename "${arr_infile[idx]%.*}")" \
                        >> "${config}"
            done
        then
            echo_error "Failed to populate GNU Parallel configuration file."
        fi

        #  Construct GNU Parallel command with argument placeholders and log
        #+ redirection
        cmd="bash ${scr_sub}"
        cmd+=" {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12}"
        cmd+="  > {11}/{12}_par.{13}.stdout.txt"  # 'scr_sub' stdout log
        cmd+=" 2> {11}/{12}_par.{13}.stderr.txt"  # 'scr_sub' stderr log

        #  Execute jobs with GNU Parallel
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
        pth_std="${err_out}/${nam_job}_ser"
        if ${dry_run} || ${verbose}; then
            echo "######################"
            echo "## Serial execution ##"
            echo "######################"
            echo ""
            echo "bash ${scr_sub} \\"
            echo "    ${env_nam} \\"
            echo "    ${scr_sig} \\"
            echo "    ${threads} \\"
            echo "    ${infiles} \\"
            echo "    ${outfiles} \\"
            echo "    ${siz_bin} \\"
            echo "    ${typ_sig} \\"
            echo "    ${scl_fct} \\"
            echo "    ${usr_frg} \\"
            echo "    ${rnd} \\"
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
                "${scr_sig}" \
                "${threads}" \
                "${infiles}" \
                "${outfiles}" \
                "${siz_bin}" \
                "${typ_sig}" \
                "${scl_fct}" \
                "${usr_frg}" \
                "${rnd}" \
                "${err_out}" \
                "${nam_job}" \
                     > "${pth_std}.stdout.txt" \
                    2> "${pth_std}.stderr.txt"
        fi
    fi
fi
