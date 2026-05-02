#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: submit_compute_signal.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


#  Require Bash >= 4.4 before doing any work
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be run under Bash >= 4.4." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

#  Run in safe mode, exiting on errors, unset variables, and pipe failures
set -euo pipefail


#  Set hardcoded arguments
## WARNING: Do not change unless testing/stepping through ##
#  If true, print verbose/debug Bash-level logging
debug=true

#  If true, dry-run script
dry_run=false

#  If true, parse arguments and exit before validation or job submission
p_only=false

#  If true, parse and check arguments, then exit before execution
pc_only=false


#  Define script-specific functions
function process_io() {
    local mode=""
    local infile=""
    local fil_A=""
    local fil_B=""
    local outfile=""
    local scl_fct=""
    local opt_var=""
    local show_help      # Help text
    local samp dsc ext   # Sample name, output descriptor, and loop variable
    local -a exts        # Recognized output/input extensions

    show_help=$(cat << EOM
Description:
  Check and parse input/output file arguments for downstream processing, returning a sample name and output descriptor.

Usage:
  process_io
    [--help] --mode <str> (--infile <str> | --fil_A <str> --fil_B <str>) --outfile <str> [--scl_fct <mlt>] [--opt_var <mlt>]

Arguments:
  -h, --hlp, --help  <flg>
    Print this help message and return 0.

  -md, --mode  <str>
    Computation mode: 'signal', 'ratio', or 'coord'.

  -i, --infile  <str>
    Scalar BAM input file ('mode=signal' or 'mode=coord').  #TODO: SAM and CRAM, or at least just CRAM, support

  -fA, --fil_A  <str>
    Scalar numerator bedGraph[.gz] file ('mode=ratio'; e.g., IP).

  -fB, --fil_B  <str>
    Scalar denominator bedGraph[.gz] file ('mode=ratio'; e.g., input).

  -o, --outfile  <str>
    Output file (bedGraph[.gz] for 'mode={signal,ratio}', BED[.gz] for 'mode=coord').

  -sf, --scl_fct  <mlt>
    Optional scaling factor (<flt>) or sentinel (NA) ('mode=signal' or 'mode=ratio').

  -ov, --opt_var  <mlt>
    Optional variable: fragment length (<int>, 'mode=signal') or minimum input depth (<flt>, 'mode=ratio').
EOM
    )

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -md|--mode)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                mode="${2}"
                shift 2
                ;;

            -i|--infile)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                infile="${2}"
                shift 2
                ;;

            -fA|--fil[_-]A)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                fil_A="${2}"
                shift 2
                ;;

            -fB|--fil[_-]B)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                fil_B="${2}"
                shift 2
                ;;

            -o|--outfile)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                outfile="${2}"
                shift 2
                ;;

            -sf|--scl[_-]fct)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                scl_fct="${2}"
                shift 2
                ;;

            -ov|--opt[_-]var)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                opt_var="${2}"
                shift 2
                ;;

            *)
                echo_err_func "${FUNCNAME[0]}" \
                    "unknown option/parameter passed: '${1}'."
                echo >&2
                echo "${show_help}" >&2
                return 1
                ;;
        esac
    done

    if [[ -z "${mode}" ]]; then
        echo_err_func "${FUNCNAME[0]}" "'--mode' is required."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var "outfile" "${outfile}" || return 1

    if [[ "${mode}" == "ratio" ]]; then
        validate_var "fil_A" "${fil_A}" || return 1
        validate_var "fil_B" "${fil_B}" || return 1
    else
        validate_var "infile" "${infile}" || return 1
    fi

    if [[ "${mode}" != "coord" ]]; then
        validate_var "scl_fct" "${scl_fct}" || return 1
        validate_var "opt_var" "${opt_var}" || return 1
    fi

    exts=( bedGraph bedGraph.gz bedgraph bedgraph.gz bdg bdg.gz bg bg.gz )
    if [[ "${mode}" == "ratio" ]]; then
        samp="${outfile##*/}"
        for ext in "${exts[@]}"; do
            samp="${samp%."${ext}"}"
        done
    else
        samp="$(basename "${infile}" ".bam")"
    fi

    exts+=( bed bed.gz )
    dsc="$(basename "${outfile}")"
    for ext in "${exts[@]}"; do
        dsc="${dsc%."${ext}"}"
    done
    unset ext

    #  Return sample name and output descriptor (comma-delimited)
    echo "${samp},${dsc}"
}


function set_args_opt() {
    local mode="${1}"
    local scl_fct="${2}"
    local opt_var="${3}"
    local rnd="${4}"
    local track="${5:-false}"
    local pseudo="${6:-NA}"
    local eps="${7:-NA}"
    local skip_00="${8:-NA}"
    local drp_nan="${9:-false}"
    local -a optional             # Optional CLI arguments
    local show_help               # Help text

    show_help=$(cat << EOM
Description:
  Build a comma-delimited list of optional CLI arguments for 'compute_signal.py' ('mode=signal') or 'compute_signal_ratio.py' ('mode=ratio').

Usage:
  set_args_opt
    [-h|--hlp|--help] mode scl_fct opt_var rnd [track] [pseudo] [eps] [skip_00] [drp_nan]

Positional arguments:
  1  mode     <str>  Mode: 'signal' or 'ratio'.
  2  scl_fct  <mlt>  Scaling factor/coefficient (<flt>) or sentinel (NA).
  3  opt_var  <mlt>  'usr_frg' ('compute_signal', <int>) or 'dep_min' ('compute_signal_ratio', <flt>) or sentinel (NA).
  4  rnd      <int>  Maximum number of decimal places retained for finite emitted values.
  5  track    <bol>  Mode 'ratio': return track sans '-inf', 'nan' (default: false).
  6  pseudo   <mlt>  Mode 'ratio': per-sample pseudocount spec 'A[:B]' (<str>) or sentinel (NA; default: NA).
  7  eps      <mlt>  Mode 'ratio': zero-tolerance epsilon (<flt>) or sentinel (NA; default: NA).
  8  skip_00  <mlt>  Mode 'ratio': zero-zero skip mode ('pre_scale' or 'post_scale') or sentinel (NA; default: NA).
  9  drp_nan  <bol>  Mode 'ratio': drop non-finite values (default: false).

Notes:
  - Returns optional arguments as a single comma-delimited line on stdout.
  - For 'mode=signal', 'opt_var' maps to '--usr_frg'.
  - For 'mode=ratio',  'opt_var' maps to '--dep_min'.
EOM
    )

    if [[ "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ $# -lt 4 || $# -gt 9 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "expected 4-9 arguments, but got $#."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    optional=()

    if [[ "${scl_fct}" != "NA" ]]; then
        optional+=( --scl_fct "${scl_fct}" )
    fi

    if [[ "${opt_var}" != "NA" ]]; then
        if [[ "${mode}" == "signal" ]]; then
            optional+=( --usr_frg "${opt_var}" )
        elif [[ "${mode}" == "ratio" ]]; then
            optional+=( --dep_min "${opt_var}" )
        fi
    fi

    if [[ "${rnd}" != "NA" ]]; then
        optional+=( --rnd "${rnd}" )
    fi

    if [[ "${mode}" == "ratio" ]]; then
        if [[ "${track}" == "true" ]]; then optional+=( --track ); fi

        if [[ "${pseudo}" != "NA" ]]; then
            optional+=( --pseudo "${pseudo}" )
        fi

        if [[ "${eps}" != "NA" ]]; then optional+=( --eps "${eps}" ); fi

        if [[ "${skip_00}" != "NA" ]]; then
            optional+=( --skip_00 "${skip_00}" )
        fi

        if [[ "${drp_nan}" == "true" ]]; then optional+=( --drp_nan ); fi
    fi

    #  Return values as comma-separated list
    ( IFS=','; echo "${optional[*]}" )
}


function run_dry_or_wet() {
    local arr_nam="${1:-}"
    local log_out="${2:-}"
    local log_err="${3:-}"
    local dir_out dir_err decl  # Derived dirs and declaration metadata
    local -a cmd_cpy            # Local copy of the command array
    local show_help             # Help text


    show_help=$(cat << EOM
Description:
  Print or execute a command stored in an array variable, with stdout/stderr redirected to log files.

Usage:
  run_dry_or_wet
    [-h|--hlp|--help] arr_nam log_out log_err

Positional arguments:
  1  arr_nam  <str>  Name of the command array.
  2  log_out  <str>  Path to file for stdout redirection.
  3  log_err  <str>  Path to file for stderr redirection.

Notes:
  - In debug or dry-run mode, prints the fully quoted command plus redirections.
  - In non-dry-run mode, executes the command and returns its exit status.
EOM
)

    if [[ "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ $# -ne 3 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "expected 3 arguments, but got $#."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! [[ "${arr_nam}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "invalid array name '${arr_nam}'."
        return 1
    fi

    if ! decl="$(declare -p "${arr_nam}" 2> /dev/null)"; then
        echo_err_func "${FUNCNAME[0]}" \
            "command array '${arr_nam}' is unset."
        return 1
    elif [[ "${decl}" != declare\ -a* ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'${arr_nam}' is not an indexed array."
        return 1
    fi

    local -n cmd_ref="${arr_nam}"

    if (( ${#cmd_ref[@]} == 0 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "received an empty command array '${arr_nam}'."
        return 1
    fi

    cmd_cpy=( "${cmd_ref[@]}" )

    #  Refuse to run if log dirs are neither existent nor writable
    dir_out="$(dirname "${log_out}")"
    dir_err="$(dirname "${log_err}")"

    if [[ ! -d "${dir_out}" || ! -d "${dir_err}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "log directory missing: '${dir_out}' or '${dir_err}'."
        return 1
    fi

    if [[ ! -w "${dir_out}" || ! -w "${dir_err}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "log directory not writable: '${dir_out}' or '${dir_err}'."
        return 1
    fi

    #  In "debug" or "dry-run mode", show the exact command and redirections
    if [[ "${debug}" == "true" || "${dry_run}" == "true" ]]; then
        printf '%q ' "${cmd_cpy[@]}" >&2
        echo ">> ${log_out} 2>> ${log_err}" >&2
        echo >&2
        echo >&2
    fi

    #  Execute the command with logging when not in dry-run mode
    if [[ "${dry_run}" == "false" ]]; then
        "${cmd_cpy[@]}" >> "${log_out}" 2>> "${log_err}"
        return $?
    fi

    return 0
}


function run_comp_sig() {
    local debug="${1}"
    local threads="${2}"
    local infile="${3}"
    local outfile="${4}"
    local siz_bin="${5}"
    local method="${6}"
    local scl_fct="${7}"
    local usr_frg="${8}"
    local rnd="${9}"
    local err_out="${10}"
    local nam_job="${11}"
    local dsc="${12}"
    local log_out log_err  # Explicit local variable declarations
    local -a optional cmd  # Optional arguments and command array
    local show_help        # Help text

    show_help=$(cat << EOM
Description:
  Build and run the per-sample call to 'compute_signal.py'.

Usage:
  run_comp_sig
    [-h|--hlp|--help] debug threads infile outfile siz_bin method scl_fct usr_frg rnd err_out nam_job dsc

Positional arguments:
   1  debug    <bol>  Print debug messages or not.
   2  threads  <int>  Number of threads to use.
   3  infile   <str>  Input BAM file.
   4  outfile  <str>  Output filename.
   5  siz_bin  <int>  Bin size in base pairs.
   6  method   <mlt>  Type of signal computation (<str>) or empty sentinel ("").
   7  scl_fct  <mlt>  Scaling factor/coefficient (<flt>) or sentinel (NA).
   8  usr_frg  <mlt>  Fragment length (<int>) or sentinel (NA).
   9  rnd      <int>  Maximum number of decimal places retained for finite emitted values.
  10  err_out  <str>  Directory for stdout and stderr.
  11  nam_job  <str>  Job name.
  12  dsc      <str>  Descriptor for log file naming.

Notes:
  - Optional CLI arguments are derived via 'set_args_opt'.
  - Logging is handled through 'run_dry_or_wet'.
EOM
    )

    if [[ "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ $# -ne 12 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'run_comp_sig()' expects 12 arguments, but got $#."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Generate optional arguments array dynamically
    IFS="," read -r -a optional <<< "$(
        set_args_opt "signal" "${scl_fct}" "${usr_frg}" "${rnd}"
    )"

    #  Debug array of optional arguments
    if [[ "${debug}" == "true" ]]; then
        echo "set_args_opt(): optional=( ${optional[*]} )" >&2
        echo >&2
    fi

    #  Define paths for log output files
    log_out="${err_out}/${nam_job}.${dsc}.stdout.txt"
    log_err="${err_out}/${nam_job}.${dsc}.stderr.txt"

    #  Build call to 'compute_signal.py' via 'run_py'
    cmd=(
        run_py "${scr_sig}"
            --verbose
            --threads "${threads}"
            --infile "${infile}"
            --outfile "${outfile}"
            --siz_bin "${siz_bin}"
    )

    if [[ -n "${method}" ]]; then
        cmd+=( --method "${method}" )
    fi

    if [[ "${#optional[@]}" -gt 0 && -n "${optional[0]}" ]]; then
        cmd+=( "${optional[@]}" )
    fi

    #  Debug or execute call to 'compute_signal.py'
    run_dry_or_wet cmd "${log_out}" "${log_err}" || return 1
}


function run_comp_rat() {
    local debug="${1}"
    local fil_A="${2}"
    local fil_B="${3}"
    local outfile="${4}"
    local method="${5}"
    local scl_fct="${6}"
    local dep_min="${7}"
    local rnd="${8}"
    local track="${9}"
    local pseudo="${10}"
    local eps="${11}"
    local skip_00="${12}"
    local drp_nan="${13}"
    local skp_pfx="${14}"
    local err_out="${15}"
    local nam_job="${16}"
    local dsc="${17}"
    local log_out log_err  # Local variable declarations
    local -a optional cmd  # Local array declarations
    local show_help        # Help text

    show_help=$(cat << EOM
Description:
  Build and run the per-sample call to 'compute_signal_ratio.py'.

Usage:
  run_comp_rat
    [-h|--hlp|--help] debug fil_A fil_B outfile method scl_fct dep_min rnd track pseudo eps skip_00 drp_nan skp_pfx err_out nam_job dsc

Positional arguments:
   1  debug    <bol>  Print debug messages or not.
   2  fil_A    <str>  Numerator bedGraph file.
   3  fil_B    <str>  Denominator bedGraph file.
   4  outfile  <str>  Output ratio bedGraph file.
   5  method   <str>  Ratio method: 'unadj', 'log2', 'unadj_r', 'log2_r'.
   6  scl_fct  <mlt>  Optional scaling factor/coefficient (<flt>) or sentinel (NA).
   7  dep_min  <mlt>  Optional minimum input depth (<flt>) or sentinel (NA).
   8  rnd      <int>  Maximum number of decimal places retained for finite emitted values.
   9  track    <bol>  Emit companion '.track' file.
  10  pseudo   <mlt>  Per-sample pseudocount spec 'A[:B]' (<str>) or sentinel (NA).
  11  eps      <mlt>  Shared epsilon (<flt>) or sentinel (NA).
  12  skip_00  <mlt>  Shared zero-zero skip mode ('pre_scale' or 'post_scale') or sentinel (NA).
  13  drp_nan  <bol>  Drop non-finite values from main output.
  14  skp_pfx  <mlt>  Shared skip-prefix string (<str>) or sentinel (NA).
  15  err_out  <str>  Directory for log file output.
  16  nam_job  <str>  Job name used in log filenames.
  17  dsc      <str>  Descriptor string for logs.

Notes:
  - Optional CLI arguments are derived via 'set_args_opt'.
  - Shared header-skip prefixes are appended via '--skp_pfx' when set.
  - Logging is handled through 'run_dry_or_wet'.
EOM
    )

    if [[ "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ $# -ne 17 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'run_comp_rat()' expects 17 arguments, but got $#."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Generate optional arguments array dynamically
    IFS="," read -r -a optional <<< "$(
        set_args_opt \
            "ratio" \
            "${scl_fct}" \
            "${dep_min}" \
            "${rnd}" \
            "${track}" \
            "${pseudo}" \
            "${eps}" \
            "${skip_00}" \
            "${drp_nan}"
    )"

    #  Debug array of optional arguments
    if [[ "${debug}" == "true" ]]; then
        echo "set_args_opt(): optional=( ${optional[*]} )" >&2
        echo >&2
    fi

    #  Define paths for log output files
    log_out="${err_out}/${nam_job}.${dsc}.stdout.txt"
    log_err="${err_out}/${nam_job}.${dsc}.stderr.txt"

    #  Build call to 'compute_signal_ratio.py' via 'run_py'
    cmd=(
        run_py "${scr_rat}"
            --verbose
            --fil_A "${fil_A}"
            --fil_B "${fil_B}"
            --fil_out "${outfile}"
            --method "${method}"
    )

    #  If any, then append optional flags
    if [[ "${#optional[@]}" -gt 0 && -n "${optional[0]}" ]]; then
        cmd+=( "${optional[@]}" )
    fi

    if [[ "${skp_pfx}" != "NA" ]]; then
        cmd+=( --skp_pfx "${skp_pfx}" )
    fi

    #  Debug or execute call to 'compute_signal_ratio.py'
    run_dry_or_wet cmd "${log_out}" "${log_err}" || return 1
}


function get_arr_elem() {
    local arr_nam="${1:-}"  # Name of indexed array
    local idx="${2:-}"      # Array index
    local decl              # Output from 'declare -p'
    local show_help         # Help message

    show_help=$(cat << EOM
Usage:
  get_arr_elem [-h|--hlp|--help] arr_nam idx

Description:
  Print one element from a named indexed array.

Positional arguments:
  1  arr_nam  <str>  Name of indexed array.
  2  idx      <int>  Zero-based array index.
EOM
    )

    if [[ "${arr_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${arr_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'arr_nam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${idx}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'idx', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif ! [[ "${arr_nam}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "invalid array name '${arr_nam}'."
        return 1
    elif ! [[ "${idx}" =~ ^[0-9]+$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "array index must be a non-negative integer: '${idx}'."
        return 1
    fi

    if ! decl="$(declare -p "${arr_nam}" 2> /dev/null)"; then
        echo_err_func "${FUNCNAME[0]}" \
            "array '${arr_nam}' is unset."
        return 1
    elif [[ "${decl}" != declare\ -a* ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'${arr_nam}' is not an indexed array."
        return 1
    fi

    local -n arr_ref="${arr_nam}"

    if (( idx >= ${#arr_ref[@]} )); then
        echo_err_func "${FUNCNAME[0]}" \
            "array index '${idx}' is out of range for '${arr_nam}'" \
            "with length '${#arr_ref[@]}'."
        return 1
    fi

    printf "%s\n" "${arr_ref[idx]}"
}


function task_pro() {
    local mode="${1}"       # Mode: 'signal', 'ratio', or 'coord'
    local idx="${2}"        # Array index (integer ≥ 0)
    local arr_fil_A="${3}"  # Array name for scalar infile or file A
    local arr_fil_B="${4}"  # Array name for scalar file B; empty if unused
    local arr_out="${5}"    # Array name for scalar outfile
    local arr_scl="${6}"    # Array name for scalar scl_fct; empty if unused
    local arr_opt="${7}"    # Array name for scalar opt_var; empty if unused
    local fil_A fil_B outfile scl_fct opt_var samp dsc
    local err_ini out_ini err_dsc out_dsc
    local show_help

    show_help=$(cat << EOM
Description:
  Prepare per-task inputs and logging metadata for the 'run_task_*' helpers ("pro" is "prologue").

Usage:
  task_pro
    [-h|--hlp|--help] mode idx arr_fil_A arr_fil_B arr_out arr_scl arr_opt

Positional arguments:
  1  mode       <str>  'signal', 'ratio', or 'coord'.
  2  idx        <int>  Zero-based task index.
  3  arr_fil_A  <str>  Array name for scalar infile or file A.
  4  arr_fil_B  <str>  Array name for scalar file B; "" if unused.
  5  arr_out    <str>  Array name for scalar outfile.
  6  arr_scl    <str>  Array name for scalar scl_fct; "" if unused.
  7  arr_opt    <str>  Array name for scalar opt_var; "" if unused.

Behavior:
  - Pulls per-task values from arrays by name.
  - Emits debug snapshots when 'debug=true'.
  - Calls 'process_io' to derive sample ('samp') and descriptor ('dsc').
  - Under Slurm, also derives initial log paths via 'set_logs_slurm'.

Returns one comma-delimited line on stdout:
  fil_A,fil_B?,outfile,samp,dsc,err_ini?,out_ini?,err_dsc?,out_dsc?
EOM
    )

    if [[ "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ $# -ne 7 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'task_pro()' expects 7 arguments, but got $#."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    fil_A="$(get_arr_elem "${arr_fil_A}" "${idx}")" || return 1

    if [[ -n "${arr_fil_B}" ]]; then
        fil_B="$(get_arr_elem "${arr_fil_B}" "${idx}")" || return 1
    fi

    outfile="$(get_arr_elem "${arr_out}" "${idx}")" || return 1

    if [[ -n "${arr_scl}" ]]; then
        scl_fct="$(get_arr_elem "${arr_scl}" "${idx}")" || return 1
    fi

    if [[ -n "${arr_opt}" ]]; then
        opt_var="$(get_arr_elem "${arr_opt}" "${idx}")" || return 1
    fi

    #  Debug inputs
    if [[ "${debug}" == "true" ]]; then
        if [[ "${mode}" == "signal" ]]; then
            debug_var \
                "infile=${fil_A}" \
                "outfile=${outfile}" \
                "scl_fct=${scl_fct}" \
                "usr_frg=${opt_var}"
        elif [[ "${mode}" == "ratio" ]]; then
            debug_var \
                "fil_A=${fil_A}" \
                "fil_B=${fil_B}" \
                "outfile=${outfile}" \
                "scl_fct=${scl_fct}" \
                "dep_min=${opt_var}"
        else
            debug_var \
                "infile=${fil_A}" \
                "outfile=${outfile}" \
                "usr_frg=${opt_var}"
        fi
    fi

    #  Derive sample ('samp') and descriptor ('dsc') via 'process_io'
    if [[ "${mode}" == "ratio" ]]; then
        IFS=',' read -r samp dsc < <(
            process_io \
                -md "${mode}" \
                -fA "${fil_A}" \
                -fB "${fil_B}" \
                 -o "${outfile}" \
                -sf "${scl_fct:-NA}" \
                -ov "${opt_var:-NA}"
        ) || return 1
    elif [[ "${mode}" == "signal" ]]; then
        IFS=',' read -r samp dsc < <(
            process_io \
                -md "${mode}" \
                 -i "${fil_A}" \
                 -o "${outfile}" \
                -sf "${scl_fct:-NA}" \
                -ov "${opt_var:-NA}"
        ) || return 1
    else
        IFS=',' read -r samp dsc < <(
            process_io \
                -md "${mode}" \
                 -i "${fil_A}" \
                 -o "${outfile}"
        ) || return 1
    fi

    #  Debug 'samp' and 'dsc' output by process_io
    if [[ "${debug}" == "true" ]]; then
        debug_var "samp=${samp}" "dsc=${dsc}"
    fi

    #  If using Slurm, request initial and descriptor log paths
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        IFS=',' read -r err_ini out_ini err_dsc out_dsc < <(
            set_logs_slurm \
                "${id_job}" "${id_tsk}" "${samp}" "${err_out}" "${nam_job}"
        ) || return 1

        if [[ "${debug}" == "true" ]]; then
            debug_var \
                "err_ini=${err_ini}" \
                "out_ini=${out_ini}" \
                "err_dsc=${err_dsc}" \
                "out_dsc=${out_dsc}"
        fi
    fi

    #  Return all values for one of the task callers to parse
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "${fil_A}" "${fil_B:-}" "${outfile}" "${samp}" "${dsc}" \
        "${err_ini:-}" "${out_ini:-}" "${err_dsc:-}" "${out_dsc:-}"
}


function task_epi() {
    local show_help

    show_help=$(cat << EOM
Description:
  Remove initial Slurm-generated log files after a task completes ("epi" is "epilogue").

Usage:
  task_epi
    [-h|--hlp|--help] err_ini out_ini

Positional arguments:
  1  err_ini  <str>  Initial stderr log path.
  2  out_ini  <str>  Initial stdout log path.

Notes:
  - No action is taken outside Slurm array tasks.
  - No action is taken in dry-run mode.
EOM
)
    if [[ "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ $# -ne 2 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'task_epi()' expects 2 arguments, but got $#."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" && "${dry_run}" == "false" ]]; then
        local err_ini="${1}"
        local out_ini="${2}"

        if [[ -n "${err_ini}" && -f "${err_ini}" ]]; then
            rm -f "${err_ini}" \
                || echo_warn \
                    "failed to remove initial stderr log: '${err_ini}'."
        fi

        if [[ -n "${out_ini}" && -f "${out_ini}" ]]; then
            rm -f "${out_ini}" \
                || echo_warn \
                    "failed to remove initial stdout log: '${out_ini}'."
        fi
    fi
}


function run_task_sig() {
    local idx="${1}"
    local infile unused outfile samp dsc err_ini out_ini err_dsc out_dsc
    local rc=0
    local show_help

    show_help=$(cat << EOM
Description:
  Run one signal-computation task, either under Slurm-array execution or in local GNU Parallel / serial iteration.

Usage:
  run_task_sig
    [-h|--hlp|--help] idx

Positional argument:
  1  idx  <int>  Zero-based task index into 'arr_infile', 'arr_outfile', 'arr_scl_fct', and 'arr_usr_frg'.

Notes:
  - Uses 'task_pro' to derive task-specific inputs and logging metadata.
  - Preserves and returns the exit status from 'run_comp_sig'.
  - Uses 'task_epi' to remove initial Slurm logs after execution.
EOM
    )

    if [[ "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ $# -ne 1 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'run_task_sig()' expects 1 argument, but got $#."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    # shellcheck disable=SC2034
    IFS=',' read -r \
        infile unused outfile samp dsc \
        err_ini out_ini err_dsc out_dsc \
    < <(
        task_pro "signal" "${idx}" \
            "arr_infile" "" "arr_outfile" "arr_scl_fct" "arr_usr_frg"
    ) || return 1

    run_comp_sig \
        "${debug}" \
        "${threads}" \
        "${infile}" \
        "${outfile}" \
        "${siz_bin}" \
        "${method}" \
        "$(get_arr_elem arr_scl_fct "${idx}")" \
        "$(get_arr_elem arr_usr_frg "${idx}")" \
        "${rnd}" \
        "${err_out}" \
        "${nam_job}" \
        "${dsc}"
    rc=$?

    task_epi "${err_ini}" "${out_ini}"
    return "${rc}"
}


function run_task_rat() {
    local idx="${1}"
    local fil_A fil_B outfile samp dsc err_ini out_ini err_dsc out_dsc
    local rc=0
    local show_help

    show_help=$(cat << EOM
Description:
  Run one ratio-computation task, either under Slurm-array execution or in local GNU Parallel / serial iteration.

Usage:
  run_task_rat
    [-h|--hlp|--help] idx

Positional argument:
  1  idx  <int>  Zero-based task index into 'arr_fil_A', 'arr_fil_B', 'arr_outfile', 'arr_scl_fct', 'arr_dep_min', and 'arr_pseudo'.

Notes:
  - Uses 'task_pro' to derive task-specific inputs and logging metadata.
  - Preserves and returns the exit status from 'run_comp_rat'.
  - Uses 'task_epi' to remove initial Slurm logs after execution.
EOM
    )

    if [[ "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ $# -ne 1 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'run_task_rat()' expects 1 argument, but got $#."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    IFS=',' read -r \
        fil_A fil_B outfile samp dsc \
        err_ini out_ini err_dsc out_dsc \
    < <(
        task_pro "ratio" "${idx}" \
            "arr_fil_A" "arr_fil_B" "arr_outfile" "arr_scl_fct" "arr_dep_min"
    ) || return 1

    run_comp_rat \
        "${debug}" \
        "${fil_A}" \
        "${fil_B}" \
        "${outfile}" \
        "${method}" \
        "$(get_arr_elem arr_scl_fct "${idx}")" \
        "$(get_arr_elem arr_dep_min "${idx}")" \
        "${rnd}" \
        "${track}" \
        "$(get_arr_elem arr_pseudo "${idx}")" \
        "${eps}" \
        "${skip_00}" \
        "${drp_nan}" \
        "${skp_pfx}" \
        "${err_out}" \
        "${nam_job}" \
        "${dsc}"
    rc=$?

    task_epi "${err_ini}" "${out_ini}"
    return "${rc}"
}


function run_task_coord() {
    local idx="${1}"
    local infile unused outfile samp dsc err_ini out_ini err_dsc out_dsc
    local rc=0
    local show_help

    show_help=$(cat << EOM
Description:
  Run one fragment-coordinate extraction task, either under Slurm-array execution or in local GNU Parallel / serial iteration.

Usage:
  run_task_coord
    [-h|--hlp|--help] idx

Positional argument:
  1  idx  <int>  Zero-based task index into 'arr_infile', 'arr_outfile', and 'arr_usr_frg'.

Notes:
  - Uses 'task_pro' to derive task-specific inputs and logging metadata.
  - Reuses 'run_comp_sig' with stub signal-mode arguments to preserve the original coordinate-extraction pathway.
  - Preserves and returns the exit status from 'run_comp_sig'.
  - Uses 'task_epi' to remove initial Slurm logs after execution.
EOM
    )

    if [[ "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ $# -ne 1 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'run_task_coord()' expects 1 argument, but got $#."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    # shellcheck disable=SC2034
    IFS=',' read -r \
        infile unused outfile samp dsc \
        err_ini out_ini err_dsc out_dsc \
    < <(
        task_pro "coord" "${idx}" \
            "arr_infile" "" "arr_outfile" "" "arr_usr_frg"
    ) || return 1

    #  (Use stub parameters per original 'coord' behavior)
    run_comp_sig \
        "${debug}" \
        1 \
        "${infile}" \
        "${outfile}" \
        1 \
        "" \
        "NA" \
        "$(get_arr_elem arr_usr_frg "${idx}")" \
        1 \
        "${err_out}" \
        "${nam_job}" \
        "${dsc}"
    rc=$?

    task_epi "${err_ini}" "${out_ini}"
    return "${rc}"
}


#  Resolve '--dir_scr' before sourced parser helpers are available
function resolve_dir_scr() {
    local script="${1:-}"
    shift

    local -a args=( "$@" )
    local i=0

    if [[ -z "${script}" ]]; then
        script="unknown_script"
    fi

    for (( i = 0; i < ${#args[@]}; i++ )); do
        case "${args[i]}" in
            -ds|--dir[_-]scr)
                if (( i + 1 >= ${#args[@]} )) \
                    || [[ -z "${args[i + 1]:-}" || "${args[i + 1]}" == -* ]]
                then
                    echo "error(${script}):" \
                        "option '${args[i]}' requires a value." >&2
                    echo >&2
                    show_help_main
                    return 1
                fi

                printf "%s\n" "${args[i + 1]}"
                return 0
                ;;
        esac
    done

    echo "error(${script}):" \
        "required option '--dir_scr' was not supplied." >&2
    echo >&2
    show_help_main
    return 1
}


#  Source 'source_helpers.sh' and requested helper scripts from 'dir_scr'
function source_submit_helpers() {
    local script="${1:-}"
    local dir_scr_arg="${2:-}"
    local fnc_src

    if (( $# < 2 )); then
        echo "error(${script:-unknown_script}):" \
            "expected at least two arguments: 'script' and 'dir_scr_arg'." >&2
        return 1
    fi

    shift 2

    if [[ -z "${script}" ]]; then
        script="unknown_script"
    fi

    if [[ -z "${dir_scr_arg}" ]]; then
        echo "error(${script}):" \
            "positional argument 2, 'dir_scr_arg', is missing." >&2
        return 1
    elif [[ ! -d "${dir_scr_arg}" ]]; then
        echo "error(${script}):" \
            "script directory not found: '${dir_scr_arg}'." >&2
        return 1
    elif (( $# < 1 )); then
        echo "error(${script}):" \
            "at least one helper script name must be supplied." >&2
        return 1
    fi

    dir_fnc="${dir_scr_arg}/functions"
    fnc_src="${dir_fnc}/source_helpers.sh"

    if [[ ! -f "${fnc_src}" ]]; then
        echo "error(${script}):" \
            "script not found: '${fnc_src}'." >&2
        return 1
    fi

    # shellcheck disable=SC1090
    source "${fnc_src}" || {
        echo "error(${script}):" \
            "failed to source '${fnc_src}'." >&2
        return 1
    }

    source_helpers "${dir_fnc}" "$@" || {
        echo "error(${script}):" \
            "failed to source required helper scripts." >&2
        return 1
    }
}


#  Parse keyword arguments after helper scripts have been sourced
function parse_args() {
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -en|--env|--env[_-]nam)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                env_nam="${2}"
                shift 2
                ;;

            -ds|--dir[_-]scr)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                dir_scr="${2}"
                shift 2
                ;;

            -t|--thr|--threads)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                threads="${2}"
                shift 2
                ;;

            -md|--mode)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                mode="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -me|--method)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                method="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -i|-fi|-ci|--infile|--infiles|--fil[_-]in|--csv[_-]infile|--csv[_-]infiles)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                csv_infile="${2}"
                shift 2
                ;;

            -fA|-f1|-cA|-c1|--fil[_-]A|--fil[_-]1|--csv[_-]A|--csv[_-]1|--csv[_-]fil[_-]A|--csv[_-]fil[_-]IP)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                csv_fil_A="${2}"
                shift 2
                ;;

            -fB|-f2|-cB|-c2|--fil[_-]B|--fil[_-]2|--csv[_-]B|--csv[_-]2|--csv[_-]fil[_-]B|--csv[_-]fil[_-]input)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                csv_fil_B="${2}"
                shift 2
                ;;

            -o|-fo|-co|--outfile|--outfiles|--fil[_-]out|--csv[_-]outfile|--csv[_-]outfiles)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                csv_outfile="${2}"
                shift 2
                ;;

            -tr|--track)
                track=true
                shift 1
                ;;

            -sb|--siz[_-]bin)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                siz_bin="${2}"
                shift 2
                ;;

            -sf|--scale|--scl[_-]fct|--csv[_-]scl[_-]fct)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                csv_scl_fct="${2}"
                shift 2
                ;;

            -uf|--usr[_-]frg|--csv[_-]usr[_-]frg)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                csv_usr_frg="${2}"
                shift 2
                ;;

            -dm|--dep[_-]min|--depth[_-]min|--csv[_-]dep[_-]min)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                csv_dep_min="${2}"
                shift 2
                ;;

            -ps|--pseudo|--pseudocount|--csv[_-]pseudo)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                csv_pseudo="${2}"
                shift 2
                ;;

            -e|--eps)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                eps="${2}"
                shift 2
                ;;

            -s0|--skip[_-]00)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                skip_00="${2}"
                shift 2
                ;;

            -dn|--drp[_-]nan)
                drp_nan=true
                shift 1
                ;;

            -sk|--skp[_-]pfx|--skip[_-]pfx|--skip[_-]prefix)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                skp_pfx="${2}"
                shift 2
                ;;

            -dp|--dp|--rnd|--round|--decimals|--digits)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                rnd="${2}"
                shift 2
                ;;

            -eo|--err[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                err_out="${2}"
                shift 2
                ;;

            -nj|--nam[_-]job)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                nam_job="${2}"
                shift 2
                ;;

            *)
                echo_err "unknown parameter passed: '${1}'."
                echo >&2
                show_help_main
                return 1
                ;;
        esac
    done
}


#  Canonicalize mode, method, and derived job name
function canonicalize_args() {
    case "${mode}" in
        signal)
            if [[ -z "${method}" ]]; then
                method="norm"
            fi

            case "${method}" in
                u|unadj|unadjusted|s|smp|simple|r|raw)
                    method="unadj"
                    ;;
                f|frg|frag|frg[_-]len|frag[_-]len|l|len|len[_-]frg|len[_-]frag)
                    method="frag"
                    ;;
                n|nrm|norm|normalized)
                    method="norm"
                    ;;
                *)
                    echo_err \
                        "invalid value for '--method': '${method}'." \
                        "Expected a signal method alias for 'unadj', 'frag'," \
                        "or 'norm'."
                    return 1
                    ;;
            esac
            ;;

        ratio)
            if [[ -z "${method}" ]]; then
                method="unadj"
            fi

            case "${method}" in
                u|unadj|unadjusted|s|smp|simple|r|raw)
                    method="unadj"
                    ;;
                2|l2|lg2|log2)
                    method="log2"
                    ;;
                ur|unadj[_-]r|unadjusted[_-]r|sr|smp[_-]r|simple[_-]r|rr|raw[_-]r)
                    method="unadj_r"
                    ;;
                2r|l2r|l2[_-]r|lg2[_-]r|log2[_-]r)
                    method="log2_r"
                    ;;
                *)
                    echo_err \
                        "invalid value for '--method': '${method}'." \
                        "Expected a ratio method alias for 'unadj', 'log2'," \
                        "'unadj_r', or 'log2_r'."
                    return 1
                    ;;
            esac
            ;;

        coord)
            method=""
            ;;

        *)
            echo_err \
                "invalid value for '--mode': '${mode}'." \
                "Expected 'signal', 'ratio', or 'coord'."
            return 1
            ;;
    esac

    if [[ -z "${nam_job}" ]]; then
        if [[ "${mode}" == "coord" ]]; then
            nam_job="compute_${mode}"
        else
            nam_job="compute_${mode}_${method}"
        fi
    fi
}


#  Validate required arguments and simple numeric values
function validate_args() {
    validate_var     "env_nam" "${env_nam}"         || return 1
    validate_var_dir "dir_scr" "${dir_scr}" 0 false || return 1
    validate_var     "threads" "${threads}"         || return 1

    check_int_pos "${threads}" "threads" || return 1
    check_int_pos "${rnd}"     "rnd"     || return 1

    if [[ "${mode}" == "ratio" ]]; then
        case "${skip_00}" in
            NA|pre_scale|post_scale) : ;;
            *)
                echo_err \
                    "'--skip_00' must be 'NA', 'pre_scale', or" \
                    "'post_scale': '${skip_00}'."
                return 1
                ;;
        esac
    fi

    if [[ "${mode}" == "signal" ]]; then
        validate_var "csv_infile"  "${csv_infile}"  || return 1
        validate_var "csv_outfile" "${csv_outfile}" || return 1
        validate_var "siz_bin"     "${siz_bin}"     || return 1
        check_int_pos "${siz_bin}" "siz_bin"        || return 1
    elif [[ "${mode}" == "coord" ]]; then
        validate_var "csv_infile"  "${csv_infile}"  || return 1
        validate_var "csv_outfile" "${csv_outfile}" || return 1
    elif [[ "${mode}" == "ratio" ]]; then
        validate_var "csv_fil_A"   "${csv_fil_A}"   || return 1
        validate_var "csv_fil_B"   "${csv_fil_B}"   || return 1
        validate_var "csv_outfile" "${csv_outfile}" || return 1
    fi

    validate_var "nam_job" "${nam_job}" || return 1

    if [[ -z "${err_out}" ]]; then
        echo_err "'--err_out' is required."
        return 1
    fi

    validate_var_dir "err_out" "${err_out}" || return 1
}


#  Resolve script paths after 'dir_scr' has been validated and helpers sourced
function resolve_script_paths() {
    validate_var_dir "dir_scr" "${dir_scr}" 0 false || return 1

    dir_rep="$(cd "${dir_scr}/.." && pwd)"
    scr_sig="${dir_scr}/compute_signal.py"
    scr_rat="${dir_scr}/compute_signal_ratio.py"

    validate_var_file "scr_sig" "${scr_sig}" || return 1
    validate_var_file "scr_rat" "${scr_rat}" || return 1
}


#  Initialize argument variables, assigning default values where applicable
env_nam="env_protocol"
dir_scr=""
threads=4
mode="signal"
method=""
csv_infile=""
csv_fil_A=""
csv_fil_B=""
csv_outfile=""
track=false
siz_bin=10
csv_scl_fct=""
csv_usr_frg=""
csv_dep_min=""
csv_pseudo=""
eps="NA"
skip_00="NA"
drp_nan=false
skp_pfx="NA"
rnd=24
err_out=""
nam_job=""

#  Define the help message
function show_help_main() {
    cat << EOM >&2
Usage:
  submit_compute_signal.sh
    [--help] [--env_nam <str>] --dir_scr <str> [--threads <int>] [--mode {signal,ratio,coord}]
    [--method {unadj,frag,norm,log2,unadj_r,log2_r,...}]
    (--csv_infile <bam1,bam2,...> | --csv_fil_A <IP1.bdg[.gz],IP2.bdg[.gz],...> --csv_fil_B <in1.bdg[.gz],in2.bdg[.gz],...>)
    --csv_outfile <out1.bdg[.gz],out2.bdg[.gz],...> [--track] [--drp_nan]
    [--siz_bin <int>] [--csv_scl_fct <mlt1,mlt2,...>] [--csv_usr_frg <mlt1,mlt2,...>] [--csv_dep_min <mlt1,mlt2,...>] [--csv_pseudo <mlt1,mlt2,...>]
    [--eps <mlt>] [--skip_00 <mlt>] [--skp_pfx <mlt>] [--rnd <int>] --err_out <str> [--nam_job <str>]

Description:
  Submit per-sample signal, ratio, or fragment-coordinate jobs from comma-separated file lists to 'compute_signal.py' or 'compute_signal_ratio.py' under Slurm, GNU Parallel, or serial execution.

Keyword arguments:
  -h, --hlp, --help  <flg>
      Print this help message and exit.

  -en, --env_nam  <str>
      Mamba environment to activate (default: ${env_nam}).

  -ds, --dir_scr  <str>
      Directory containing scripts and functions.

  -t, --thr, --threads  <int>
      Number of threads to use per job (default: ${threads}).

  -md, --mode  <str>
      Type of computation: 'signal', 'ratio', or 'coord' (default: ${mode});

  -me, --method <str>
      Computation subtype (ignored if '--mode coord'; default if '--mode signal': norm; default if '--mode ratio': unadj):
        - For '--mode signal':
          + 'u', 'unadj', 'unadjusted', 's', 'smp', 'simple', 'r', 'raw'
          + 'f', 'frg', 'frag', 'frg_len', 'frag_len', 'l', 'len', 'len_frg', 'len_frag'
          + 'n', 'nrm', 'norm', 'normalized'
        - For '--mode ratio':
          + 'u', 'unadj', 'unadjusted', 's', 'smp', 'simple', 'r', 'raw'
          + '2', 'l2', 'lg2', 'log2'
          + 'ur', 'unadj_r', 'unadjusted_r', 'sr', 'smp_r', 'simple_r', 'rr', 'raw_r'
          + '2r', 'l2r', 'l2_r', 'lg2_r', 'log2_r'

  -i, -fi, -ci, --infile, --infiles, --fil_in, --csv_infile, --csv_infiles  <mlt>  <element: str>
      Comma-separated list of BAM files ('--mode signal', '--mode coord').

  -fA, -f1, -cA, -c1, --fil_A, --fil_1, --csv_A, --csv_1, --csv_fil_A, --csv_fil_IP  <mlt>  <element: str>
      Comma-separated list of numerator (i.e., "file A"; e.g., IP) bedGraph files ('--mode ratio').

  -fB, -f2, -cB, -c2, --fil_B, --fil_2, --csv_B, --csv_2, --csv_fil_B, --csv_fil_input  <mlt>  <element: str>
      Comma-separated list of denominator (i.e., "file B"; e.g., input) bedGraph files ('--mode ratio').

  -o, -fo, -co, --outfile, --outfiles, --fil_out, --csv_outfile, --csv_outfiles  <mlt>  <element: str>
      Comma-separated list of output files (e.g., full bedGraph[.gz] or BED[.gz] paths).

  -tr, --track  <flg>
      Output a companion bedGraph without '-inf' or 'nan' rows ('--mode ratio').

  -sb, --siz_bin  <int>
      Bin size in base pairs for signal computation ('--mode signal'; default: ${siz_bin}).

  -sf, --csv_scl_fct  <mlt>  <element: str>
      Optional comma-separated list of scaling factors or sentinels ('--mode signal', '--mode ratio').

      Compatibility aliases include '--scale' and '--scl_fct'.

  -uf, --csv_usr_frg  <mlt>  <element: str>
      Optional comma-separated list of fragment lengths or sentinels ('--mode signal', '--mode coord').

      Compatibility alias: '--usr_frg'.

  -dm, --csv_dep_min  <mlt>  <element: str>
      Optional comma-separated list of minimum input depth values or sentinels ('--mode ratio').

      Compatibility aliases include '--dep_min' and '--depth_min'.

  -ps, --csv_pseudo  <mlt>  <element: str>
      Optional comma-separated list of per-sample pseudocount specs 'A[:B]' or sentinels ('--mode ratio').

      Compatibility aliases include '--pseudo' and '--pseudocount'.

   -e, --eps  <mlt>
      Shared epsilon for zero checks in ratio mode (<flt>) or sentinel (NA).

  -s0, --skip_00  <str>
      Shared zero-zero skip mode in ratio mode ('pre_scale' or 'post_scale') or sentinel (NA).

  -dn, --drp_nan  <str>
      Drop 'nan', 'inf', and '-inf' rows from the main ratio output.

  -sk, --skp_pfx, --skip_pfx, --skip_prefix  <mlt>
      Shared comma-separated bedGraph header prefixes to skip in ratio mode (<str>) or sentinel (NA).

  -dp, --dp, --rnd, --round, --decimals, --digits  <int>
      Maximum number of decimal places retained for finite output signal or ratio values (default: ${rnd}).

  -eo, --err_out  <str>
      Directory for stderr and stdout logs.

  -nj, --nam_job  <str>
      Prefix for job names (default depends on resolved '--mode' and '--method'; e.g., 'compute_signal_norm', 'compute_ratio_unadj', or 'compute_coord').

Notes:
  - BAM and bedGraph input files must be coordinate-sorted.
  - Input and output paths supplied to this wrapper interface must not contain spaces, commas, or semicolons.
  - This wrapper does not support '-' for stdin/stdout. Use the underlying Python scripts directly for streaming input/output workflows.
  - If and where applicable, use consistent file ordering in- and outfiles, and between IP (file A) and input (file B) files.
  - To run in "debug mode", set hardcoded variable 'debug=true'                   [debug=${debug:-UNSET}]
  - To run in "dry-run mode", set hardcoded variable 'dry_run=true'               [dry_run=${dry_run:-UNSET}]
  - To run in "parse-only mode", set hardcoded variable 'p_only=true'             [p_only=${p_only:-UNSET}]
  - To run in "parse-and-check-only mode", set hardcoded variable 'pc_only=true'  [pc_only=${pc_only:-UNSET}]
EOM
}
#TODO: the following needs to be incorporated into the help docs
#+   Sourced function scripts:
#+     - check_args.sh
#+     - check_inputs.sh
#+     - check_numbers.sh
#+     - format_outputs.sh
#+     - handle_env.sh
#+     - manage_slurm.sh
#+     - populate_array_empty.sh
#+     - run_python.sh


#  Main script execution
function main() {
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        show_help_main
        echo >&2
        exit 0
    fi

    dir_scr="$(resolve_dir_scr "${0##*/}" "$@")" || exit 1

    source_submit_helpers "${0##*/}" "${dir_scr}" \
        check_args \
        check_inputs \
        check_numbers \
        format_outputs \
        handle_env \
        manage_slurm \
        populate_array_empty \
        run_python \
        || exit 1

    parse_args "$@"       || exit 1

    if [[ "${p_only}" == "true" ]]; then
        if [[ "${debug}" == "true" ]]; then debug_var "p_only=true"; fi
        exit 0
    fi

    canonicalize_args     || exit 1
    validate_args         || exit 1
    resolve_script_paths  || exit 1

    #  Debug argument variable assignments
    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "env_nam=${env_nam}" \
            "dir_scr=${dir_scr}" \
            "threads=${threads}" \
            "mode=${mode}" \
            "method=${method:-UNSET}"

        if [[ "${mode}" != "ratio" ]]; then
            debug_var \
                "csv_infile=${csv_infile}"
        elif [[ "${mode}" == "ratio" ]]; then
            debug_var \
                "csv_fil_A=${csv_fil_A}" \
                "csv_fil_B=${csv_fil_B}"
        fi

        debug_var "csv_outfile=${csv_outfile}"

        if [[ "${mode}" == "ratio" ]]; then
            debug_var \
                "track=${track}" \
                "csv_dep_min=${csv_dep_min}" \
                "csv_pseudo=${csv_pseudo}" \
                "eps=${eps}" \
                "skip_00=${skip_00}" \
                "drp_nan=${drp_nan}" \
                "skp_pfx=${skp_pfx}"
        fi

        if [[ "${mode}" == "signal" ]]; then
            debug_var \
                "siz_bin=${siz_bin}" \
                "csv_usr_frg=${csv_usr_frg}"
        fi

        if [[ "${mode}" != "coord" ]]; then
            debug_var \
                "csv_scl_fct=${csv_scl_fct}" \
                "rnd=${rnd}"
        fi

        debug_var \
            "err_out=${err_out}" \
            "nam_job=${nam_job}"
    fi

    #  Perform mode-dependent array reconstruction from serialized strings
    #+
    #+ Reconstruct required arrays first, then auto-populate omitted optional
    #+ per-sample arrays with 'NA' sentinels so this script is standalone-ready
    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        IFS=',' read -r -a arr_infile <<< "${csv_infile}"
        IFS=',' read -r -a arr_outfile <<< "${csv_outfile}"

        check_arr_nonempty "arr_infile"  "csv_infile"  || exit 1
        check_arr_nonempty "arr_outfile" "csv_outfile" || exit 1

        if [[ -n "${csv_usr_frg}" ]]; then
            IFS=',' read -r -a arr_usr_frg <<< "${csv_usr_frg}"
        else
            unset arr_usr_frg && declare -a arr_usr_frg
            populate_array_empty arr_usr_frg "${#arr_infile[@]}"
        fi

        if [[ "${mode}" == "signal" ]]; then
            if [[ -n "${csv_scl_fct}" ]]; then
                IFS=',' read -r -a arr_scl_fct <<< "${csv_scl_fct}"
            else
                unset arr_scl_fct && declare -a arr_scl_fct
                populate_array_empty arr_scl_fct "${#arr_infile[@]}"
            fi
        fi
    elif [[ "${mode}" == "ratio" ]]; then
        IFS=',' read -r -a arr_fil_A   <<< "${csv_fil_A}"
        IFS=',' read -r -a arr_fil_B   <<< "${csv_fil_B}"
        IFS=',' read -r -a arr_outfile <<< "${csv_outfile}"

        check_arr_nonempty "arr_fil_A"   "csv_fil_A"   || exit 1
        check_arr_nonempty "arr_fil_B"   "csv_fil_B"   || exit 1
        check_arr_nonempty "arr_outfile" "csv_outfile" || exit 1

        if [[ -n "${csv_scl_fct}" ]]; then
            IFS=',' read -r -a arr_scl_fct <<< "${csv_scl_fct}"
        else
            unset arr_scl_fct && declare -a arr_scl_fct
            populate_array_empty arr_scl_fct "${#arr_fil_A[@]}"
        fi

        if [[ -n "${csv_dep_min}" ]]; then
            IFS=',' read -r -a arr_dep_min <<< "${csv_dep_min}"
        else
            unset arr_dep_min && declare -a arr_dep_min
            populate_array_empty arr_dep_min "${#arr_fil_A[@]}"
        fi

        if [[ -n "${csv_pseudo}" ]]; then
            IFS=',' read -r -a arr_pseudo <<< "${csv_pseudo}"
        else
            unset arr_pseudo && declare -a arr_pseudo
            populate_array_empty arr_pseudo "${#arr_fil_A[@]}"
        fi
    fi

    #  Check final array lengths after any auto-population of optional vectors
    if [[ "${mode}" == "signal" ]]; then
        check_arr_lengths "arr_infile" "arr_outfile" || exit 1
        check_arr_lengths "arr_infile" "arr_scl_fct" || exit 1
        check_arr_lengths "arr_infile" "arr_usr_frg" || exit 1
    elif [[ "${mode}" == "ratio" ]]; then
        check_arr_lengths "arr_fil_A" "arr_fil_B"    || exit 1
        check_arr_lengths "arr_fil_A" "arr_outfile"  || exit 1
        check_arr_lengths "arr_fil_A" "arr_scl_fct"  || exit 1
        check_arr_lengths "arr_fil_A" "arr_dep_min"  || exit 1
        check_arr_lengths "arr_fil_A" "arr_pseudo"   || exit 1
    elif [[ "${mode}" == "coord" ]]; then
        check_arr_lengths "arr_infile" "arr_outfile" || exit 1
        check_arr_lengths "arr_infile" "arr_usr_frg" || exit 1
    fi

    #  Reject stdin input in this wrapper; use the Python scripts directly
    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        for infile in "${arr_infile[@]}"; do
            if [[ "${infile}" == "-" ]]; then
                echo_err \
                    "'-' is not allowed in '--csv_infile' for" \
                    "'submit_compute_signal.sh'. Use the underlying Python" \
                    "script directly for stdin input."
                exit 1
            fi
        done
        unset infile
    elif [[ "${mode}" == "ratio" ]]; then
        for infile in "${arr_fil_A[@]}"; do
            if [[ "${infile}" == "-" ]]; then
                echo_err \
                    "'-' is not allowed in '--csv_fil_A' for" \
                    "'submit_compute_signal.sh'. Use the underlying Python" \
                    "script directly for stdin input."
                exit 1
            fi
        done

        for infile in "${arr_fil_B[@]}"; do
            if [[ "${infile}" == "-" ]]; then
                echo_err \
                    "'-' is not allowed in '--csv_fil_B' for" \
                    "'submit_compute_signal.sh'. Use the underlying Python" \
                    "script directly for stdin input."
                exit 1
            fi
        done
        unset infile
    fi

    #  Reject stdout output in this wrapper; use the Python scripts directly
    for outfile in "${arr_outfile[@]}"; do
        if [[ "${outfile}" == "-" ]]; then
            echo_err \
                "'-' is not allowed in '--csv_outfile' for" \
                "'submit_compute_signal.sh'. Use the underlying Python" \
                "script directly for stdout output."
            exit 1
        fi
    done
    unset outfile

    #  Check input files
    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        for idx in "${!arr_infile[@]}"; do
            validate_var_file "arr_infile" "${arr_infile[${idx}]}" "${idx}" \
                || exit 1
        done
    elif [[ "${mode}" == "ratio" ]]; then
        for idx in "${!arr_fil_A[@]}"; do
            validate_var_file "arr_fil_A" "${arr_fil_A[${idx}]}" "${idx}" \
                || exit 1
        done

        for idx in "${!arr_fil_B[@]}"; do
            validate_var_file "arr_fil_B" "${arr_fil_B[${idx}]}" "${idx}" \
                || exit 1
        done
    fi

    #  Debug number of array elements and array element values
    if [[ "${debug}" == "true" ]]; then
        if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
            debug_var \
                "\${#arr_infile[@]}=${#arr_infile[@]}" \
                "\${#arr_usr_frg[@]}=${#arr_usr_frg[@]}"
            echo "arr_infile=( ${arr_infile[*]} )"   >&2 && echo >&2
            echo "arr_usr_frg=( ${arr_usr_frg[*]} )" >&2 && echo >&2
        fi

        if [[ "${mode}" == "ratio" ]]; then
            debug_var \
                "\${#arr_fil_A[@]}=${#arr_fil_A[@]}" \
                "\${#arr_fil_B[@]}=${#arr_fil_B[@]}" \
                "\${#arr_dep_min[@]}=${#arr_dep_min[@]}" \
                "\${#arr_pseudo[@]}=${#arr_pseudo[@]}"
            echo "arr_fil_A=( ${arr_fil_A[*]} )"     >&2 && echo >&2
            echo "arr_fil_B=( ${arr_fil_B[*]} )"     >&2 && echo >&2
            echo "arr_dep_min=( ${arr_dep_min[*]} )" >&2 && echo >&2
            echo "arr_pseudo=( ${arr_pseudo[*]} )"   >&2 && echo >&2
        fi

        if [[ "${mode}" != "coord" ]]; then
            debug_var "\${#arr_scl_fct[@]}=${#arr_scl_fct[@]}"
            echo "arr_scl_fct=( ${arr_scl_fct[*]} )" >&2 && echo >&2
        fi

        debug_var "\${#arr_outfile[@]}=${#arr_outfile[@]}"
        echo "arr_outfile=( ${arr_outfile[*]} )"     >&2 && echo >&2
    fi

    #  Exit after parsing and validation/checking, but before execution
    if [[ "${pc_only}" == "true" ]]; then
        if [[ "${debug}" == "true" ]]; then debug_var "pc_only=true"; fi
        exit 0
    fi


    #  Do the main work -----------------------------------------------------------
    #  Activate environment
    handle_env "${env_nam}" || exit 1

    #  For 'scripts.*' imports, ensure project root is on 'PYTHONPATH' (exported
    #+ after environment activation to prevent clobbering)
    #+
    #+ Note: 'run_py' handles 'PYTHONPATH' when 'PY_INVOKE=module', so this is
    #+ redundant but also harmless
    if [[ -z "${PYTHONPATH:-}" ]]; then
        export PYTHONPATH="${dir_rep}"
    else
        export PYTHONPATH="${dir_rep}:${PYTHONPATH}"
    fi

    #  Determine and run mode: Slurm or GNU Parallel/serial
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        #  Job submission type: Slurm
        id_job=${SLURM_ARRAY_JOB_ID}
        id_tsk=${SLURM_ARRAY_TASK_ID}

        if ! [[ "${id_tsk}" =~ ^[1-9][0-9]*$ ]]; then
            echo_err "Slurm task ID is invalid: '${id_tsk}'."
            exit 1
        elif [[
            "${mode}" == "ratio" && "${id_tsk}" -gt "${#arr_fil_A[@]}"
        ]]; then
            echo_err \
                "Slurm task ID '${id_tsk}' exceeds number of ratio entries:" \
                "'${#arr_fil_A[@]}'."
            exit 1
        elif [[
            "${mode}" != "ratio" && "${id_tsk}" -gt "${#arr_infile[@]}"
        ]]; then
            echo_err \
                "Slurm task ID '${id_tsk}' exceeds number of input entries:" \
                "'${#arr_infile[@]}'."
            exit 1
        else
            idx=$(( id_tsk - 1 ))
        fi

        if [[ "${mode}" == "signal" ]]; then
            run_task_sig "${idx}"   || exit 1
        elif [[ "${mode}" == "ratio" ]]; then
            run_task_rat "${idx}"   || exit 1
        else
            run_task_coord "${idx}" || exit 1
        fi
    else
        #  Job submission type: GNU Parallel or serial
        if [[ "${mode}" == "signal" ]]; then
            for idx in "${!arr_infile[@]}"; do
                run_task_sig "${idx}"   || exit 1
            done
        elif [[ "${mode}" == "ratio" ]]; then
            for idx in "${!arr_fil_A[@]}"; do
                run_task_rat "${idx}"   || exit 1
            done
        else
            for idx in "${!arr_infile[@]}"; do
                run_task_coord "${idx}" || exit 1
            done
        fi
    fi
}


main "$@"
