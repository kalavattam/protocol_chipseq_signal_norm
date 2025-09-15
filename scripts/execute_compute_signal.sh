#!/bin/bash

#  execute_compute_signal.sh
#  KA


#  Run script in interactive mode (true) or command-line mode (false)
interactive=false  ## WARNING: User should change needed ##

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if [[ "${interactive}" == "false" ]]; then set -euo pipefail; fi

#  Set the path to the "scripts" directory
if [[ "${interactive}" == "true" ]]; then
    ## WARNING: If 'interactive=true', user should change needed ##
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
    extract_field_str \
    format_print_cmd \
    handle_env \
    make_lower \
    populate_array_empty \
    print_parallel_info \
    reset_max_job \
    set_params_parallel \
    summarize_sig_norm
do
    source "${dir_fnc}/${fnc}.sh"
done
unset fnc


#  Determine prefix for ratio track filenames based on provided arguments
function generate_prefix() {
    local method="${1}"
    local scl_fct="${2}"

    if [[ "${method}" == "log2" ]]; then
        if [[ -n "${scl_fct}" ]]; then
            echo "scl_log2_rat"
        else
            echo "log2_rat"
        fi
    else
        if [[ -n "${scl_fct}" ]]; then
            echo "scl_rat"
        else
            echo "rat"
        fi
    fi
}


#  Construct the call to 'scr_sub' with appropriate flags based on context
function build_cmd() {
    local idx="${1:-UNSET}"  # Optional index for GNU Parallel or serial modes
    local infile fil_ip fil_in outfile scl_fct usr_frg dep_min
    local cmd wrap base

    #  Assign default local values from global variables or arrays, but only
    #+ for the active mode
    if [[ "${mode}" == "signal" || "${mode}" == "coord" ]]; then
        infile="${ser_infile-}"
        outfile="${ser_outfile-}"
        scl_fct="${ser_scl_fct-}"
        usr_frg="${ser_usr_frg-}"
    elif [[ "${mode}" == "ratio" ]]; then
        fil_ip="${ser_fil_ip-}"
        fil_in="${ser_fil_in-}"
        outfile="${ser_outfile-}"
        scl_fct="${ser_scl_fct-}"
        dep_min="${ser_dep_min-}"
    fi

    if [[ "${idx}" != "UNSET" ]]; then
        infile="${arr_infile[idx]}"
        outfile="${arr_outfile[idx]}"
        scl_fct="${arr_scl_fct[idx]}"
        usr_frg="${arr_usr_frg[idx]}"
        dep_min="${arr_dep_min[idx]}"

        if [[ "${mode}" == "ratio" ]]; then
            fil_ip="${arr_fil_ip[idx]}"
            fil_in="${arr_fil_in[idx]}"
        fi
    fi

    #  Build command
    cmd="${scr_sub}"
    cmd+=" --env_nam ${env_nam}"
    cmd+=" --dir_scr ${dir_scr}"
    cmd+=" --threads ${threads}"
    cmd+=" --mode ${mode}"

    if [[ "${mode}" != "coord" ]]; then
        cmd+=" --method ${method}"
    fi

    if [[ "${mode}" == "signal" || "${mode}" == "coord" ]]; then
        cmd+=" --ser_infile ${infile}"
    elif [[ "${mode}" == "ratio" ]]; then
        cmd+=" --ser_fil_ip ${fil_ip}"
        cmd+=" --ser_fil_in ${fil_in}"
    fi

    cmd+=" --ser_outfile ${outfile}"

    if [[ "${mode}" == "signal" ]]; then
        cmd+=" --siz_bin ${siz_bin}"
        cmd+=" --ser_scl_fct ${scl_fct}"
        cmd+=" --ser_usr_frg ${usr_frg}"
    elif [[ "${mode}" == "ratio" ]]; then
        cmd+=" --ser_scl_fct ${scl_fct}"
        cmd+=" --ser_dep_min ${dep_min}"
        if [[ "${track}" == "true" ]]; then
            cmd+=" --track"
        fi
    fi

    cmd+=" --rnd ${rnd}"
    cmd+=" --err_out ${err_out}"
    cmd+=" --nam_job ${nam_job}"

    if [[ "${slurm:-false}" == "true" ]]; then
        wrap="sbatch --job-name=${nam_job}"
        wrap+=" --nodes=1"
        wrap+=" --cpus-per-task=${threads}"
        wrap+=" --time=${time}"
        wrap+=" --output=${err_out}/${nam_job}.%A-%a.stdout.txt"
        wrap+=" --error=${err_out}/${nam_job}.%A-%a.stderr.txt"
        wrap+=" --array=1-${#arr_outfile[@]}%${max_job}"
    else
        base=$(basename "${outfile}")
        base="${base%.gz}"
        base="${base%.bedgraph}"
        base="${base%.bdg}"
        base="${base%.bg}"

        cmd+=" >> ${err_out}/${nam_job}.${base}.stdout.txt"
        cmd+=" 2>> ${err_out}/${nam_job}.${base}.stderr.txt"
        wrap="bash"
    fi

    #  Return command
    echo "${wrap} ${cmd}"
}


#  Set paths, values, arguments, etc. for interactive mode
function set_interactive() {
    #  Set base repository paths
    dir_bas="${HOME}/repos"  ## WARNING: User should change as needed ##
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

    #  Set hardcoded argument assignments (and related assignments)
    verbose=true
    dry_run=true
    threads=6
    sf=""  # "spike"  # "alpha"
    mode="ratio"  # "signal"
    method="log2"  # "unadj"  # "norm"

    if [[ "${sf}" =~ ^(alpha|spike)$ ]]; then
        tbl="${dir_pro}/compute_signal/${str_det}/tables/${sf}_test.tsv"
    elif [[ "${mode}" == "ratio" && "${method}" == "log2" && -z "${sf}" ]]; then
        use_dm=true
        tbl="${dir_pro}/compute_signal/${str_det}/tables/alpha_test.tsv"
    fi

    if [[ "${mode}" != "ratio" ]]; then
        infiles="$(  ## WARNING: User should change as needed ##
            bash "${dir_scr}/find_files.sh" \
                --dir_fnd "${dir_pro}/align_fastqs/${str_det}/sc" \
                --pattern "*.bam" \
                --exclude "*Brn1*"
        )"
    else
        fil_ip="$(  ## WARNING: User should change as needed ##
            bash "${dir_scr}/find_files.sh" \
                --dir_fnd "${dir_pro}/compute_signal/${str_det}/norm" \
                --pattern "*.bdg.gz" \
                --include "IP*" \
                --exclude "*Brn1*"
        )"
        fil_in="$(sed 's:\/IP:\/in:g' < <(echo "${fil_ip}"))"
    fi

    dir_out="${dir_pro}/compute_${mode}/${str_det}/${method}"
    typ_out="bdg.gz"

    if [[ "${sf}" =~ ^(alpha|spike)$ ]]; then
        prefix="${sf}"
    else
        prefix=""
    fi

    if [[ "${mode}" != "ratio" ]]; then
        track=false
    else
        track=true
        scl_fct="$(
            if [[ "${sf}" == "alpha" ]]; then
                extract_field_str "${tbl}" 3 false
            elif [[ "${sf}" == "spike" ]]; then
                extract_field_str "${tbl}" 5 false
            fi
        )"
        dep_min="$(
            if [[ "${sf}" == "alpha" || "${use_dm}" == "true" ]]; then
                extract_field_str "${tbl}" 24 false
            elif [[ "${sf}" == "spike" ]]; then
                extract_field_str "${tbl}" 21 false
            fi
        )"
    fi

    siz_bin=10
    rnd=24
    err_out="${dir_out}/logs"
    nam_job="compute_${mode}_${method}"
    max_job=2
    slurm=false
    time="0:30:00"

    if true; then
        echo "#######################################################"
        echo "## Debug variable assignments for 'interactive mode' ##"
        echo "#######################################################"
        echo ""
        echo "verbose=${verbose}"
        echo "dry_run=${dry_run}"
        echo "threads=${threads}"
        echo "mode=${mode}"
        echo "method=${method}"
        echo "infiles=${infiles:-UNSET}"
        echo "fil_ip=${fil_ip:-UNSET}"
        echo "fil_in=${fil_in:-UNSET}"
        echo "dir_out=${dir_out}"
        echo "typ_out=${typ_out}"
        echo "prefix=${prefix:-UNSET}"
        echo "track=${track:-UNSET}"
        echo "scl_fct=${scl_fct:-UNSET}"
        echo "dep_min=${dep_min:-UNSET}"
        echo "siz_bin=${siz_bin:-UNSET}"
        echo "rnd=${rnd:-UNSET}"
        echo "err_out=${err_out}"
        echo "nam_job=${nam_job}"
        echo "max_job=${max_job}"
        echo "slurm=${slurm}"
        echo "time=${time}"
    fi
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_protocol"
scr_sub="${dir_scr}/submit_compute_signal.sh"
par_job=""

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=4
mode="signal"
method="norm"
infiles=""
fil_ip=""
fil_in=""
dir_out=""
typ_out="bdg.gz"
prefix=""
siz_bin=10
scl_fct=""
usr_frg=""
dep_min=""
rnd=24
err_out=""
nam_job="compute_${mode}_${method}"
max_job=6
slurm=false
time="0:30:00"

#  Assign variables for help messages
show_help=$(cat << EOM
Usage:
  execute_compute_signal.sh
    [--help] [--details] [--verbose] [--dry_run] --threads <int> --mode <str> --method <str> --infiles <str> --fil_ip <str> --fil_in <str> --dir_out <str>
    --typ_out <str> [--prefix <str>] --siz_bin <int> [--scl_fct <flt>] [--usr_frg <flt>] [--dep_min <str>] --rnd <int> --err_out <str> --nam_job <str>
    --max_job <int> [--slurm] [--time <str>]

Description:
  This driver script automates the computation of signal tracks, ratio tracks, or fragment coordinate files from BAM or bedGraph input files. It supports
  multiple normalization strategies and can run in serial or parallel via GNU Parallel or SLURM.

  For more details on what this script can do, including copious notes and example usage, run the following:
  \`\`\`
  bash ${0} --details
  \`\`\`

Arguments:
   -h, --help     Print this short help message and exit.
   -d, --details  Print full documentation with notes and examples, then exit.
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Run script in 'dry-run mode' (optional).
   -t, --threads  Number of threads to use (default: "${threads}").
   -m, --mode     Type of computation to perform: 'signal', 'ratio', or 'coord' (default: "${mode}").
  -me, --method   Signal or ratio computation subtype (default: "${method}"). Used only with '--mode signal' or '--mode ratio'.
   -i, --infiles  Comma-separated list of coordinate-sorted BAM files (used only with '--mode signal' or '--mode coord').
  -ip, --fil_ip   Comma-separated list of coordinate-sorted bedGraph files for ChIP IP signal (used only with '--mode ratio').
  -in, --fil_in   Comma-separated list of coordinate-sorted bedGraph files for input signal (used only with '--mode ratio').
  -do, --dir_out  Output directory for generated files.
  -to, --typ_out  Format of output signal, ratio, or coordinate files (default: "${typ_out}").
  -pr, --prefix   Custom prefix to prepend to output filenames (used only with '--mode signal' or '--mode ratio').
  -tr, --track    If '--mode ratio', also output a companion bedGraph file with all '-inf' and 'nan' rows removed (optional).
  -sb, --siz_bin  Bin size in base pairs for signal computation (default: "${siz_bin}"; used only with '--mode signal').
  -sf, --scl_fct  Comma-separated list of scaling factors for signal or ratio values (optional; used only with '--mode signal' or '--mode ratio').
  -uf, --usr_frg  Comma-separated list of fragment lengths to use instead of read/template lengths (optional; used with '--mode signal').
  -dm, --dep_min  Comma-separated list of minimum input depth values to avoid extreme division (optional; used only with '--mode ratio').
   -r, --rnd      Decimal precision for signal scores (default: "${rnd}"; sed only with '--mode signal' or '--mode ratio').
  -eo, --err_out  Directory for stderr and stdout TXT output files (default: "\${dir_out}/err_out").
  -nj, --nam_job  Prefix for job names (default: "${nam_job}").
  -mj, --max_job  Maximum concurrent jobs to run (default: "${max_job}").
  -sl, --slurm    Submit jobs to the SLURM scheduler (optional).
  -tm, --time     SLURM job time in 'h:mm:ss' format (required if '--slurm'; default: "${time}").
EOM
)

show_details=$(cat << EOM
Usage:
  execute_compute_signal.sh
    [--help] [--details] [--verbose] [--dry_run] --threads <int> --mode <str> --method <str> --infiles <str> --fil_ip <str> --fil_in <str> --dir_out <str>
    --typ_out <str> [--prefix <str>] --siz_bin <int> [--scl_fct <flt>] [--usr_frg <flt>] [--dep_min <str>] --rnd <int> --err_out <str> --nam_job <str>
    --max_job <int> [--slurm] [--time <str>]

Description:
  The driver script 'execute_compute_signal.sh' automates the computations of bedGraph signal or ratio tracks, or BED-like fragment coordinate files, from BAM
  (for signal tracks or fragment coordinate files) or bedGraph (for ratio tracks) inputs.

  It supports multiple signal normalization strategies, including
    - unadjusted (raw) signal
    - fragment-length adjusted signal (Dickson et al., JBC 2020; Dickson et al., Sci Rep 2023)
    - normalized coverage (Dickson et al., Sci Rep 2023)
    - input-normalized signal ratios (siQ-ChIP; Dickson et al., JBC 2020; Dickson et al., Sci Rep 2023)
    - log2-transformed input-normalized signal ratios (e.g., as described in "Data Analysis G" of Alavattam et al., Bio-protocol 2025)
    - spike-in-normalized signal (i.e., unmodified ChIP-Rx; Orlando et al., Cell Rep 2014)
    - input- and spike-in-normalized signal ratios (modified ChIP-Rx as described in "Data Analysis I" of Alavattam et al., Bio-protocol 2025)

  The script supports parallel job execution via SLURM or GNU Parallel, or can execute jobs in serial.

Arguments:
   -h, --help     Print a short help message and exit.

   -d, --details  Print this full documentation with notes and examples, then exit.

   -v, --verbose  Run script in 'verbose mode' (optional).

  -dr, --dry_run  Run script in 'dry-run mode' (optional).

   -t, --threads  Number of threads to use (default: "${threads}").

   -m, --mode     Type of computation to perform: 'signal', 'ratio', or 'coord' (default: "${mode}"). Available options:
                    - 's', 'sig', 'signal':
                      + Compute signal tracks directly from BAM input files.
                      + Supports unadjusted signal, fragment-length-adjusted signal, or normalized coverage.
                      + See '--method' for calculation styles.
                      + Use this option with appropriately computed scaling factors ('--scl_fct') to compute unmodified ChIP-Rx (Orlando et al., Cell Rep 2014)
                        spike-in-normalized signal.
                      + If 's' or 'sig' are supplied, variable 'mode' is set to "signal".
                    - 'r', 'rat', 'ratio':
                      + Compute IP/input signal ratios from bedGraph files.
                      + Supports both untransformed and log2-transformed ratios.
                      + Untransformed ratios are used in the modified ChIP-Rx spike-in and siQ-ChIP normalizations described in Alavattam et al., Bio-protocol
                        2025.
                      + Use this option with '--method log2' to compute log2(IP/input) ratios.
                      + If 'r' or 'rat' are supplied, variable 'mode' is set to "ratio".
                    - 'c', 'coord', 'coordinates':
                      + Instead of computing signal or ratio tracks, output fragment coordinates in BED-like format from BAM input files.
                      + Use this option to prepare input files for the original siQ-ChIP implementation (Dickson et al., JBC 2020; Dickson et al.,
                        Sci Rep 2023). For more details, see github.com/BradleyDickson/siQ-ChIP or github.com/kalavattam/siQ-ChIP.
                      + This mode disables '--siz_bin' and '--method', and sets '--typ_out' to 'bed.gz' by default (or 'bed' if '--typ_out bed' is specified).
                      + If 'c' or 'coordinates' are supplied, variable 'mode' is set to "coord".

  -me, --method   Signal or ratio computation subtype (default: "${method}"), used only with '--mode signal' or '--mode ratio'.
                    - If '--mode signal', then the available options are
                      + 'u', 'unadj', 'unadjusted':
                        - Compute unadjusted signal.
                        - If 'u' or 'unadjusted', variable 'method' is set to "unadj".
                      + 'l', 'len', 'len_frag', 'frag':
                        - Adjust signal by fragment length.
                        - If 'l', 'len', or 'len_frag', variable 'method' is set to "frag".
                        - For example, one might use this option to compute siQ-ChIP-scaled signal using the initial siQ-ChIP equation described in Dickson et
                          al., JBC 2020.
                      + 'n', 'norm', 'normalized':
                        - Compute normalized coverage per Dickson et al., Sci Rep 2023, adjusting by fragment length and total fragments (unity scaling).
                        - If 'n' or 'normalized', variable 'method' is set to "norm".
                    - If '--mode ratio', then the available options are
                      + 'u', 'unadj', 'unadjusted':
                        - Compute unadjusted (non-log2) IP/input ratio.
                        - If 'u' or 'unadjusted', variable 'method' is set to "unadj".
                      + '2', 'l2', 'log2':
                        - Compute log2(IP/input) ratio.
                        - If '2' or 'l2', variable 'method' is set to "log2".

   -i, --infiles  Comma-separated list of coordinate-sorted BAM infiles.
                    - Required when '--mode signal' or '--mode coord'.
                    - Ignored for '--mode ratio'.

  -ip, --fil_ip   Comma-separated list of coordinate-sorted bedGraph files representing ChIP IP signal tracks.
                    - Required when '--mode ratio'; ignored otherwise.
                    - The list order must match that of '--fil_in' files.

  -in, --fil_in   Comma-separated list of coordinate-sorted bedGraph files representing input signal tracks.
                    - Required when '--mode ratio'; ignored otherwise.
                    - The list order must match that of '--fil_ip' files.

  -do, --dir_out  Output directory for generated files:
                    - Signal tracks if '--mode signal'.
                    - Ratio tracks if '--mode ratio'.
                    - BED-like files of fragment coordinates if '--mode coord'.

  -to, --typ_out  Format of signal track output files (default: "${typ_out}"). Available options:
                    - 'bedgraph', 'bdg', 'bg':
                      + Signal/ratio in bedGraph format.
                      + Cannot be used with '--mode coord'.
                    - 'bedgraph.gz', 'bdg.gz', 'bg.gz':
                      + Signal/ratio in gzip-compressed bedGraph format.
                      + Cannot be used with '--mode coord'.
                    - 'bed', 'bed.gz':
                      + BED-like format for fragment coordinates instead of signal.
                      + Can only be used with '--mode coord'.

  -pr, --prefix   Custom prefix to prepend to output filenames.
                    - When '--mode signal':
                        + If not specified, no prefix is added.
                        + If specified, the prefix is prepended to the base filename, and any leading 'IP_' or 'in_' string in the base filename is stripped
                          before applying it.
                    - When '--mode ratio':
                        + If not specified, a default prefix is automatically constructed based on '--method' and '--scl_fct':
                            - 'rat' (default)
                            - 'log2_rat' (if '--method log2')
                            - 'scl_rat' (if '--scl_fct' is supplied)
                            - 'scl_log2_rat' (if '--method log2' and '--scl_fct' is supplied)
                        + If specified, the custom prefix replaces the default.
                        + Whether specified or not, any leading 'IP_' string in the base name is stripped before the prefix.
                    - When '--mode coordinates', this argument is ignored.

  -tr, --track    If '--mode ratio', also output a companion bedGraph file with all '-inf' and 'nan' rows removed (optional).
                    - The new file will include '.track' before the extension.
                    - This cleaned version is ideal for visualization in genome browsers such as IGV, avoiding issues caused by '-inf' or 'nan' values.

  -sb, --siz_bin  Bin size in base pairs for signal computation (default: "${siz_bin}").
                    - Used only with '--mode signal'; ignored otherwise.

  -sf, --scl_fct  Comma-separated list of scaling factors to apply to signal or ratio values.
                    - Used with either '--mode signal' or '--mode ratio'; ignored otherwise.
                    - List size must match the number of input files via '--infiles' or '--fil_ip'/'--fil_in'.

  -uf, --usr_frg  Comma-separated list of fragment lengths to use instead of read lengths (single-end alignments) or template lengths (paired-end alignments;
                  optional).
                    - Used with either '--mode signal' or '--mode coord'; ignored otherwise.
                    - List size must match the number of input files via '--infiles'.

  -dm, --dep_min  Comma-separated list of minimum input depth values used to avoid extreme division operations (optional).
                    - Used only with '--mode ratio'; ignored otherwise.
                    - List size must match the number of input files via '--fil_ip'/'--fil_in'.

   -r, --rnd      Decimal precision for signal scores (default: "${rnd}").

  -eo, --err_out  Directory for stderr and stdout TXT output files (default: "\${dir_out}/err_out").

  -nj, --nam_job  Prefix for job names (default: "${nam_job}").

  -mj, --max_job  Maximum concurrent jobs to run (default: "${max_job}").
                    - With '--slurm': Number of SLURM array tasks.
                    - Without '--slurm':
                      + If 'max_job' is greater than 1, jobs run in parallel via GNU Parallel.
                      + If 'max_job' is 1, jobs run in serial, i.e., sequentially.

  -sl, --slurm    Submit jobs to the SLURM scheduler (optional).
                    - If '--slurm' is not specified and 'threads' is greater than 1, then run jobs in parallel with GNU Parallel.
                    - If '--slurm' is not specified and 'threads' is 1, run jobs in serial.

  -tm, --time     The length of time, in 'h:mm:ss' format, for the SLURM job (required if '--slurm' is specified, ignored if not; default: "${time}").

Dependencies:
  - Bash or Zsh
  - GNU Parallel (if not using SLURM and 'threads' is greater than 1)
  - Python
  - SLURM (if using '--slurm')

Notes:
  - BAM and bedGraph input files must be coordinate-sorted.
  - If applicable, use consistent file ordering between IP and input files.
  - '--typ_out' must be compatible with selected '--mode'; otherwise, an error is thrown.  ## #TODO: Is this true? Check this... ##
  - Output filenames are derived from BAM or bedGraph input files and the value associated with '--typ_out'.
  - BED-like files of fragment coordinates are, e.g., used as input to the original siQ-ChIP implementation (Dickson et al., JBC 2020; Dickson et al., Sci Rep
    2023).
  - Job execution mode (serial, GNU Parallel, SLURM) is chosen automatically based on '--threads' and '--slurm'.
    ## #TODO: More explanation. ##

Examples:
  \`\`\`
  #  Example 1: Compute normalized coverage using GNU Parallel
  bash "\${HOME}/scripts/execute_compute_signal.sh" \\
      --threads 8 \\
      --mode "signal" \\
      --method "norm" \\
      --infiles "\${HOME}/project/samples/sample_1.bam,\${HOME}/project/samples/sample_2.bam" \\
      --dir_out "\${HOME}/project/tracks" \\
      --typ_out "bdg.gz" \\
      --siz_bin 50 \\
      --err_out "\${HOME}/project/logs" \\
      --nam_job "norm_sig"

  #  Example 2: Compute log2 IP/input ratios from bedGraph files in serial
  bash "\${HOME}/scripts/execute_compute_signal.sh" \\
      --threads 1 \\
      --mode "ratio" \\
      --method "log2" \\
      --fil_ip "\${HOME}/project/norm/IP_1.bdg,\${HOME}/project/norm/IP_2.bdg" \\
      --fil_in "\${HOME}/project/norm/in_1.bdg,\${HOME}/project/norm/in_2.bdg" \\
      --dir_out "\${HOME}/project/ratios" \\
      --typ_out "bg"
  \`\`\`
EOM
)

#  Parse arguments
if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    exit_0
elif [[ "${1}" == "-d" || "${1}" == "--details" ]]; then
    echo "${show_details}"
    exit_0
fi

if [[ "${interactive}" == "true" ]]; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -v|--verbose) verbose=true;                   shift 1 ;;
            -dr|--dry_run) dry_run=true;                   shift 1 ;;
             -t|--threads) threads="${2}";                 shift 2 ;;
             -m|--mode)    mode="$(make_lower "${2}")";    shift 2 ;;
            -me|--method)  method="$(make_lower "${2}")";  shift 2 ;;
             -i|--infiles) infiles="${2}";                 shift 2 ;;
            -ip|--fil_ip)  fil_ip="${2}" ;                 shift 2 ;;
            -in|--fil_in)  fil_in="${2}" ;                 shift 2 ;;
            -do|--dir_out) dir_out="${2}";                 shift 2 ;;
            -to|--typ_out) typ_out="$(make_lower "${2}")"; shift 2 ;;
            -pr|--prefix)  prefix="${2}";                  shift 2 ;;
            -tr|--track)   track=true;                     shift 1 ;;
            -sb|--siz_bin) siz_bin="${2}";                 shift 2 ;;
            -sf|--scl_fct) scl_fct="${2}";                 shift 2 ;;
            -uf|--usr_frg) usr_frg="${2}";                 shift 2 ;;
            -dm|--dep_min) dep_min="${2}";                 shift 2 ;;
             -r|--rnd)     rnd="${2}";                     shift 2 ;;
            -eo|--err_out) err_out="${2}";                 shift 2 ;;
            -nj|--nam_job) nam_job="${2}";                 shift 2 ;;
            -mj|--max_job) max_job="${2}";                 shift 2 ;;
            -sl|--slurm)   slurm=true;                     shift 1 ;;
            -tm|--time)    time="${2}";                    shift 2 ;;
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

check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

case "${mode}" in
    s|sig|signal)
        mode="signal"

        case "${method}" in
            u|unadj|unadjusted)  method="unadj" ;;
            l|len|len_frag|frag) method="frag"   ;;
            n|norm|normalized)   method="norm"  ;;
            *)
                echo_error \
                    "Invalid value for '--method': '${method}'. Expected" \
                    "'u', 'unadj', 'unadjusted', 'l', 'len', 'len_frag'," \
                    "'frag', 'n', 'norm', or 'normalized'."
                exit_1
                ;;
        esac

        check_supplied_arg -a "${infiles}" -n "infiles"
        check_exists_file_dir "d" "$(dirname "${infiles%%[,;]*}")" "infiles"
        check_str_delim "infiles" "${infiles}"

        check_supplied_arg -a "${typ_out}" -n "typ_out"
        case "${typ_out}" in
            bedgraph|bedgraph.gz|bdg|bdg.gz|bg|bg.gz) : ;;
            *)
                echo_error \
                    "Invalid value for '--typ_out': '${typ_out}'. Expected" \
                    "'bedgraph', 'bedgraph.gz', 'bdg', 'bdg.gz', 'bg', or" \
                    "'bg.gz'."
                exit_1
                ;;
        esac
        ;;

    r|rat|ratio)
        mode="ratio"

        case "${method}" in
            u|unadj|unadjusted) method="unadj" ;;
            2|l2|log2)          method="log2"  ;;
            *)
                echo_error \
                    "Invalid value for '--method': '${method}'. Expected" \
                    "'u', 'unadj', 'unadjusted', '2', 'l2', or 'log2'."
                exit_1
                ;;

        esac

        check_supplied_arg -a "${fil_ip}" -n "fil_ip"
        check_exists_file_dir "d" "$(dirname "${fil_ip%%[,;]*}")" "fil_ip"
        check_str_delim "fil_ip" "${fil_ip}"

        check_supplied_arg -a "${fil_in}" -n "fil_in"
        check_exists_file_dir "d" "$(dirname "${fil_in%%[,;]*}")" "fil_in"
        check_str_delim "fil_in" "${fil_in}"
        ;;

    c|coord|coordinates)
        mode="coord"

        # #MAYBE: 'siz_bin' not applicable in coord mode; silently unset
        # unset siz_bin
        
        if [[ -n "${method}" ]]; then
            echo_warning \
                "Argument '--method' is not applicable with '--mode" \
                "coordinates'. Ignoring/unsetting '--method ${method}'."
        fi
        unset method

        check_supplied_arg -a "${infiles}" -n "infiles"
        check_exists_file_dir "d" "$(dirname "${infiles%%[,;]*}")" "infiles"
        check_str_delim "infiles" "${infiles}"
        ;;

    *)
        echo_error \
            "Invalid value for '--mode': '${mode}'. Expected 's', 'sig'," \
            "'signal', 'r', 'rat', 'ratio', 'c', 'coord', or 'coordinates'."
        exit_1
        ;;
esac

check_supplied_arg -a "${dir_out}" -n "dir_out"
check_exists_file_dir "d" "${dir_out}" "dir_out"

if [[ "${mode}" =~ ^(s|sig|signal|r|rat|ratio)$ ]]; then
    check_supplied_arg -a "${typ_out}" -n "typ_out"
    case "${typ_out}" in
        bedgraph|bedgraph.gz|bdg|bdg.gz|bg|bg.gz) : ;;
        *)
            echo_error \
                "Invalid value for '--typ_out': '${typ_out}'. Expected" \
                "'bedgraph', 'bedgraph.gz', 'bdg', 'bdg.gz', 'bg', or" \
                "'bg.gz'."
            exit_1
            ;;
    esac

    if [[ -n "${prefix}" && ! "${prefix}" =~ ^[a-zA-Z0-9._-]+$ ]]; then
        echo_warning \
            "User-supplied '--prefix' contains unusual characters:" \
            "'${prefix}'. Proceeding, but this may result in malformed" \
            "filenames or other issues."
    fi

    if [[ "${mode}" == "signal" ]]; then
        check_supplied_arg -a "${siz_bin}" -n "siz_bin"
        check_int_pos "${siz_bin}" "siz_bin"
    else
        #FIXME: May want to change this, b/c of default 'siz_bin=10'
        if [[ -n "${siz_bin}" ]]; then
            echo_warning \
                "Argument '--siz_bin' is not applicable with '--mode ${mode}'." \
                "Ignoring/unsetting '--siz_bin ${siz_bin}'."
            unset siz_bin
        fi
    fi
else
    check_supplied_arg -a "${typ_out}" -n "typ_out"
    case "${typ_out}" in
        bedgraph|bedgraph.gz|bdg|bdg.gz|bg|bg.gz)
            echo_warning \
                "Unsupported value for '--typ_out' with '--mode ${mode}':" \
                "'${typ_out}'. Resetting '--typ_out' to 'bed.gz'."
            typ_out="bed.gz"
            ;;

        bed|bed.gz) : ;;

        *)
            echo_error \
                "Invalid value for '--typ_out': '${typ_out}'. Expected" \
                "'bed' or 'bed.gz'."
            exit_1
            ;;
    esac
fi

if [[ -n "${scl_fct}" ]]; then check_str_delim "scl_fct" "${scl_fct}"; fi
if [[ -n "${usr_frg}" ]]; then check_str_delim "usr_frg" "${usr_frg}"; fi
if [[ -n "${dep_min}" ]]; then check_str_delim "dep_min" "${dep_min}"; fi

check_supplied_arg -a "${rnd}" -n "rnd"
check_int_pos "${rnd}" "rnd"

if [[ -z "${err_out}" ]]; then err_out="${dir_out}/err_out"; fi
check_exists_file_dir "d" "${err_out}" "err_out"

check_supplied_arg -a "${nam_job}" -n "nam_job"

check_supplied_arg -a "${max_job}" -n "max_job"
check_int_pos "${max_job}" "max_job"


#  Parse and validate input vector elements -----------------------------------
if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
    if [[ -n "${infiles}" ]]; then
        IFS=',' read -r -a arr_infile <<< "${infiles}"
        check_array_files "infiles" "${arr_infile[@]}"
    else
        echo_error "Variable 'infiles' unexpectedly empty."
        exit_1
    fi

    unset arr_outfile && typeset -a arr_outfile
    for i in "${arr_infile[@]}"; do
        base="$(basename "${i}" .bam)"

        if [[ -n "${prefix}" ]]; then
            if [[ "${base}" =~ ^IP_ ]]; then
                echo_warning \
                    "Filename '${base}' starts with 'IP_'. Stripping this" \
                    "prefix before applying custom '--prefix ${prefix}'."
                base="${base#IP_}"
            elif [[ "${base}" =~ ^in_ ]]; then
                echo_warning \
                    "Filename '${base}' starts with 'in_'. Stripping this" \
                    "prefix before applying custom '--prefix ${prefix}'."
                base="${base#in_}"
            fi
            arr_outfile+=( "${dir_out}/${prefix}.${base}.${typ_out}" )
        else
            arr_outfile+=( "${dir_out}/${base}.${typ_out}" )
        fi
    done

    check_arrays_lengths "arr_outfile" "arr_infile"

    if [[ "${mode}" == "signal" ]]; then
        if [[ -z "${scl_fct}" ]]; then
            unset arr_scl_fct && typeset -a arr_scl_fct
            populate_array_empty arr_scl_fct "${#arr_infile[@]}"
        else
            IFS=',' read -r -a arr_scl_fct <<< "${scl_fct}"
        fi

        for s in "${arr_scl_fct[@]}"; do
            if [[ "${s}" != "NA" ]]; then check_flt_pos "${s}" "scl_fct"; fi
        done
        unset s

        if [[ -z "${usr_frg}" ]]; then
            unset arr_usr_frg && typeset -a arr_usr_frg
            populate_array_empty arr_usr_frg "${#arr_infile[@]}"
        else
            IFS=',' read -r -a arr_usr_frg <<< "${usr_frg}"
        fi

        for u in "${arr_usr_frg[@]}"; do
            if [[ "${u}" != "NA" ]]; then check_flt_pos "${u}" "usr_frg"; fi
        done
        unset u

        check_arrays_lengths "arr_scl_fct" "arr_infile"
        check_arrays_lengths "arr_usr_frg" "arr_infile"
    fi
else
    if [[ -n "${fil_ip}" ]]; then
        IFS=',' read -r -a arr_fil_ip <<< "${fil_ip}"
        check_array_files "fil_ip" "${arr_fil_ip[@]}"
    else
        echo_error "Variable 'fil_ip' unexpectedly empty."
        exit_1
    fi

    if [[ -n "${fil_in}" ]]; then
        IFS=',' read -r -a arr_fil_in <<< "${fil_in}"
        check_array_files "fil_in" "${arr_fil_in[@]}"
    else
        echo_error "Variable 'fil_in' unexpectedly empty."
        exit_1
    fi

    unset arr_outfile && typeset -a arr_outfile
    for i in "${arr_fil_ip[@]}"; do
        base=$(basename "${i}")

        #  Strip 'IP_' prefix (if present)
        base="${base#IP_}"

        #  If not user-assigned, determine filename 'prefix' based on provided
        #+ arguments
        if [[ -z "${prefix}" ]]; then
            prefix="$(generate_prefix "${method}" "${scl_fct}")"
        fi

        #  Remove file extensions
        for ext in bedgraph.gz bedgraph bdg.gz bdg bg.gz bg; do
            base="${base%."${ext}"}"
        done
        unset ext

        arr_outfile+=( "${dir_out}/${prefix}_${base}.${typ_out}" )
    done
    unset i base prefix

    if [[ -z "${scl_fct}" ]]; then
        unset arr_scl_fct && typeset -a arr_scl_fct
        populate_array_empty arr_scl_fct "${#arr_fil_ip[@]}"
    else
        IFS=',' read -r -a arr_scl_fct <<< "${scl_fct}"
    fi

    for s in "${arr_scl_fct[@]}"; do
        if [[ "${s}" != "NA" ]]; then check_flt_pos "${s}" "scl_fct"; fi
    done
    unset s

    if [[ -z "${dep_min}" ]]; then
        unset arr_dep_min && typeset -a arr_dep_min
        populate_array_empty arr_dep_min "${#arr_fil_ip[@]}"
    else
        IFS=',' read -r -a arr_dep_min <<< "${dep_min}"
    fi

    for d in "${arr_dep_min[@]}"; do
        if [[ "${d}" != "NA" ]]; then check_flt_pos "${d}" "dep_min"; fi
    done
    unset d

    check_arrays_lengths "arr_fil_in"  "arr_fil_ip"
    check_arrays_lengths "arr_outfile" "arr_fil_ip"
    check_arrays_lengths "arr_scl_fct" "arr_fil_ip"
    check_arrays_lengths "arr_dep_min" "arr_fil_ip"
fi


#  Parse job execution parameters ---------------------------------------------
if [[ "${slurm}" == "true" ]]; then
    if [[ "${mode}" == "ratio" ]]; then
        max_job=$(reset_max_job "${max_job}" "${#arr_outfile[@]}")
    else
        max_job=$(reset_max_job "${max_job}" "${#arr_infile[@]}")
    fi
    
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

#  Debug parallelization information and summary output of resolved states
print_parallel_info \
    "${slurm}" "${max_job:-UNSET}" "${par_job}" "${threads}" \
    "arr_infile" "arr_fil_ip" "arr_fil_in" "arr_outfile" \
    "arr_scl_fct" "arr_usr_frg" "arr_dep_min"
#TODO: Change "parameter" to "argument" in message

if [[ "${mode}" == "signal" ]]; then
    summarize_sig_norm "${method}" "${scl_fct}"
fi


#  Activate environment and check that dependencies are in PATH ---------------
handle_env "${env_nam}" > /dev/null

check_program_path python

if [[ "${slurm}" == "true" ]]; then
    check_program_path sbatch
elif [[ "${threads}" -gt 1 ]]; then
    check_program_path parallel
fi


#  Do the main work ===========================================================
if [[ "${verbose}" == "true" ]]; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "env_nam=${env_nam}"
    echo "scr_sub=${scr_sub}"
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
    echo "mode=${mode}"
    echo "method=${method}"
    echo "infiles=${infiles:-UNSET}"
    echo "fil_ip=${fil_ip:-UNSET}"
    echo "fil_in=${fil_in:-UNSET}"
    echo "dir_out=${dir_out}"
    echo "typ_out=${typ_out}"
    echo "prefix=${prefix:-UNSET}"
    echo "siz_bin=${siz_bin:-UNSET}"
    echo "scl_fct=${scl_fct:-UNSET}"
    echo "usr_frg=${usr_frg:-UNSET}"
    echo "dep_min=${dep_min:-UNSET}"
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

    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        echo "arr_infile=( ${arr_infile[*]} )"
        echo ""
    elif [[ "${mode}" == "ratio" ]]; then
        echo "arr_fil_ip=( ${arr_fil_ip[*]} )"
        echo ""
        echo "arr_fil_in=( ${arr_fil_in[*]} )"
        echo ""
    fi

    echo "arr_outfile=( ${arr_outfile[*]} )"
    echo ""

    if [[ "${mode}" != "coord" ]]; then
        echo "arr_scl_fct=( ${arr_scl_fct[*]} )"
        echo ""
    fi

    if [[ "${mode}" == "signal" ]]; then
        echo "arr_usr_frg=( ${arr_usr_frg[*]} )"
        echo ""
    elif [[ "${mode}" == "ratio" ]]; then
        echo "arr_dep_min=( ${arr_dep_min[*]} )"
        echo ""
    fi

    echo ""
fi

if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
    ser_infile=$(echo "${arr_infile[*]}"  | tr ' ' ',')
    ser_outfile=$(echo "${arr_outfile[*]}" | tr ' ' ',')
    ser_scl_fct=$(echo "${arr_scl_fct[*]}" | tr ' ' ',')
    ser_usr_frg=$(echo "${arr_usr_frg[*]}" | tr ' ' ',')
elif [[ "${mode}" == "ratio" ]]; then
    ser_fil_ip=$(echo "${arr_fil_ip[*]}"  | tr ' ' ',')
    ser_fil_in=$(echo "${arr_fil_in[*]}"  | tr ' ' ',')
    ser_outfile=$(echo "${arr_outfile[*]}" | tr ' ' ',')
    ser_scl_fct=$(echo "${arr_scl_fct[*]}" | tr ' ' ',')
    ser_dep_min=$(echo "${arr_dep_min[*]}" | tr ' ' ',')
fi

if [[ "${verbose}" == "true" ]]; then
    echo "##################################################"
    echo "## Variable assignments constructed from arrays ##"
    echo "##################################################"
    echo ""

    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        echo "ser_infile=\"${ser_infile}\""
        echo ""
    elif [[ "${mode}" == "ratio" ]]; then
        echo "ser_fil_ip=\"${ser_fil_ip}\""
        echo ""
        echo "ser_fil_in=\"${ser_fil_in}\""
        echo ""
    fi

    echo "ser_outfile=\"${ser_outfile}\""
    echo ""

    if [[ "${mode}" != "coord" ]]; then
        echo "ser_scl_fct=\"${ser_scl_fct}\""
        echo ""
    fi

    if [[ "${mode}" == "signal" ]]; then
        echo "ser_usr_frg=\"${ser_usr_frg}\""
        echo ""
    elif [[ "${mode}" == "ratio" ]]; then
        echo "ser_dep_min=\"${ser_dep_min}\""
        echo ""
    fi

    echo ""
fi

# shellcheck disable=SC2016,SC2090
if [[ "${slurm}" == "true" ]]; then
    #  SLURM execution
    if [[ "${dry_run}" == "true" || "${verbose}" == "true" ]]; then
        echo "######################"
        echo "## Call to 'sbatch' ##"
        echo "######################"
        echo ""
        build_cmd "UNSET" | format_print_cmd "${slurm}" "${scr_sub}"
        echo ""
        echo ""
    fi

    if [[ "${dry_run}" == "false" ]]; then eval "$(build_cmd)"; fi
else
    #  GNU Parallel execution
    if [[ ${threads} -gt 1 ]]; then
        #  Create and populate GNU Parallel configuration file
        config="${err_out}/${nam_job}.config_parallel.txt"

        #  Clean old config file
        if [[ -f "${config}" ]]; then rm "${config}"; fi

        #  Populate config file with full job commands
        for idx in "${!arr_outfile[@]}"; do
            build_cmd "${idx}" >> "${config}" || {
                echo_error "Failed to write command for index no. '${idx}'."
                exit_1
            }
        done  # cat "${config}"

        #  Execute jobs with GNU Parallel
        if [[ "${dry_run}" == "true" || "${verbose}" == "true" ]]; then
            echo "############################"
            echo "## GNU Parallel execution ##"
            echo "############################"
            echo ""
            parallel --jobs "${par_job}" --dryrun < "${config}"
            echo ""
            echo ""
        fi

        if [[ "${dry_run}" == "false" ]]; then
            parallel --jobs "${par_job}" < "${config}"
        fi
    else
        #  Serial execution
        pth_std="${err_out}/${nam_job}_ser"
        if [[ "${dry_run}" == "true" || "${verbose}" == "true" ]]; then
            echo "######################"
            echo "## Serial execution ##"
            echo "######################"
            echo ""
            echo "bash ${scr_sub} \\"
            echo "    ${env_nam} \\"
            echo "    ${dir_scr} \\"
            echo "    ${threads} \\"
            echo "    ${ser_infile} \\"
            echo "    ${ser_outfile} \\"
            echo "    ${siz_bin} \\"
            echo "    ${method} \\"
            echo "    ${ser_scl_fct} \\"
            echo "    ${ser_usr_frg} \\"
            echo "    ${rnd} \\"
            echo "    ${err_out} \\"
            echo "    ${nam_job} \\"
            echo "         > ${pth_std}.stdout.txt \\"
            echo "        2> ${pth_std}.stderr.txt"
            echo ""
            echo ""
        fi

        if [[ "${dry_run}" == "false" ]]; then
            bash "${scr_sub}" \
                "${env_nam}" \
                "${dir_scr}" \
                "${threads}" \
                "${ser_infile}" \
                "${ser_outfile}" \
                "${siz_bin}" \
                "${method}" \
                "${ser_scl_fct}" \
                "${ser_usr_frg}" \
                "${rnd}" \
                "${err_out}" \
                "${nam_job}" \
                     > "${pth_std}.stdout.txt" \
                    2> "${pth_std}.stderr.txt"
        fi
    fi
fi
