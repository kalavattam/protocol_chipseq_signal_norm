#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: submit_compute_signal.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
#   GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


# Require Bash >= 4.4 before doing any work.
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

# Run in safe mode, exiting on errors, unset variables, and pipe failures.
set -euo pipefail

# 'sbatch' copies this script, so '--dir_scr' must outrank 'BASH_SOURCE'.
function resolve_dir_scr() {
    local -a args=( "$@" )
    local i

    for (( i = 0; i < ${#args[@]}; i++ )); do
        case "${args[i]}" in
            -ds|--dir[_-]scr)
                if \
                       (( i + 1 < ${#args[@]} )) \
                    && [[ -n "${args[i + 1]}" ]] \
                    && [[ "${args[i + 1]}" != -* ]]
                then
                    printf '%s\n' "${args[i + 1]}"
                    return 0
                fi
                ;;
        esac
    done

    cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
}

dir_scr="$(resolve_dir_scr "$@")" || {
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to resolve the script directory." >&2
    exit 1
}

# shellcheck source=lib/bash/help/help_submit_compute_signal.sh
source "${dir_scr}/../lib/bash/help/help_submit_compute_signal.sh" || {
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to source main help." >&2
    exit 1
}


# Define script-specific functions.
function process_io() {
    local mode=""
    local fil_in=""
    local fil_A=""
    local fil_B=""
    local fil_out=""
    local scl_fct=""
    local opt_var=""
    local show_help     # Help text
    local samp dsc ext  # Sample name, output descriptor, and loop variable
    local -a exts       # Recognized output/input extensions

    show_help=$(cat << EOM
Usage
-----
  process_io
    [--help] --mode <mode> (--fil_in <file> | --fil_A <file> --fil_B <file>) --fil_out <file> [--scl_fct <num>] [--opt_var <num>]

  Check and parse input/output file arguments for downstream processing, returning a sample name and output descriptor.

Parameters
----------
  -h, --help : flag
    Print this help message and return 0.

  -md, --mode : {'signal', 'ratio', 'coord'}
    Workflow mode: 'signal', 'ratio', or 'coord'.

  -fi, --fil_in : file
    Input file path. Scalar BAM or CRAM input file ('mode=signal' or 'mode=coord').

  -fA, --fil_A : file
    First bedGraph input file, file A ('mode=ratio'; e.g., IP).

  -fB, --fil_B : file
    Second bedGraph input file, file B ('mode=ratio'; e.g., input).

  -fo, --fil_out : file
    Output file path. Output file (bedGraph[.gz] for 'mode={signal,ratio}', BED[.gz] for 'mode=coord').

  -sf, --scl_fct : number
    Scaling factor (<flt>) or sentinel (NA) ('mode=signal' or 'mode=ratio').

  -ov, --opt_var : number
    Optional variable: fragment length (<int>, 'mode=signal') or minimum input depth (<flt>, 'mode=ratio').

Returns
-------
  Prints a comma-delimited sample name and output descriptor to stdout.

Notes
-----
  Runtime requirements:
    - basename
    - bash >= 4.4

Examples
--------
  1. Derive names for one signal-mode BAM and bedGraph output.
    '''bash
    process_io \\
        --mode signal \\
        --fil_in sample.bam \\
        --fil_out sample.norm.bedGraph \\
        --scl_fct NA \\
        --opt_var NA
    '''

  2. Derive names for a ratio-mode bedGraph pair with scaling metadata.
    '''bash
    process_io \\
        --mode ratio \\
        --fil_A IP.bedGraph \\
        --fil_B input.bedGraph \\
        --fil_out IP_over_input.bedGraph \\
        --scl_fct 2:1 \\
        --opt_var 0.5
    '''
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

            -fi|--fil_in)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                fil_in="${2}"
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

            -fo|--fil_out)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                fil_out="${2}"
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

    validate_var "fil_out" "${fil_out}" || return 1

    if [[ "${mode}" == "ratio" ]]; then
        validate_var "fil_A" "${fil_A}" || return 1
        validate_var "fil_B" "${fil_B}" || return 1
    else
        validate_var "fil_in" "${fil_in}" || return 1
    fi

    if [[ "${mode}" != "coord" ]]; then
        validate_var "scl_fct" "${scl_fct}" || return 1
        validate_var "opt_var" "${opt_var}" || return 1
    fi

    exts=( bedGraph bedGraph.gz bedgraph bedgraph.gz bdg bdg.gz bg bg.gz )
    if [[ "${mode}" == "ratio" ]]; then
        samp="${fil_out##*/}"
        for ext in "${exts[@]}"; do
            samp="${samp%."${ext}"}"
        done
    else
        samp="$(basename "${fil_in}")"
        samp="${samp%.bam}"
        samp="${samp%.cram}"
    fi

    exts+=( bed bed.gz )
    dsc="$(basename "${fil_out}")"
    for ext in "${exts[@]}"; do
        dsc="${dsc%."${ext}"}"
    done
    unset ext

    # Return sample name and output descriptor (comma-delimited).
    echo "${samp},${dsc}"
}


function set_args_opt() {
    local mode="${1:-}"
    local scl_fct="${2:-}"
    local opt_var="${3:-}"
    local dp="${4:-}"
    local track="${5:-false}"
    local pseudo="${6:-NA}"
    local eps="${7:-NA}"
    local skip_00="${8:-NA}"
    local drp_nan="${9:-false}"
    local -a optional  # Optional CLI arguments
    local show_help    # Help text

    show_help=$(cat << EOM
Usage
-----
  set_args_opt
    [--help] mode scl_fct opt_var dp [track] [pseudo] [eps] [skip_00] [drp_nan]

  Build a comma-delimited list of optional CLI arguments for 'compute_signal.py' ('mode=signal') or 'compute_signal_ratio.py' ('mode=ratio').

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  mode : {'signal', 'ratio'}
    Workflow mode. Mode: 'signal' or 'ratio'.

  2  scl_fct : number
    Scaling factor. Use sentinel 'NA' to omit scaling.

  3  opt_var : number
    Optional mode-specific value: fixed fragment length for 'compute_signal.py' or minimum input depth for 'compute_signal_ratio.py'. Use sentinel 'NA' to omit.

  4  dp : int
    Maximum number of decimal places retained for finite emitted values.

  5  track : bool
    Write a companion track file. Mode 'ratio': return track sans '-inf', 'nan' (default: false).

  6  pseudo : structured string
    Per-file pseudocount spec 'A[:B]'. Mode 'ratio' only; use sentinel 'NA' to omit (default: 'NA').

  7  eps : float
    Zero tolerance epsilon or sentinel 'NA'. Mode 'ratio' only (default: 'NA').

  8  skip_00 : {'pre_scale', 'post_scale'}
    Skip rows where both compared values are zero. Mode 'ratio': zero-zero skip mode ('pre_scale' or 'post_scale') or sentinel (NA; default: NA).

  9  drp_nan : bool
    Drop non-finite values from main output. Mode 'ratio' only (default: false).

Returns
-------
  - Prints optional arguments as a single comma-delimited line to stdout.
  - Returns 0 when optional arguments are built successfully; 1 otherwise.

Notes
-----
  Runtime requirements:
    bash >= 4.4

  - For 'mode=signal', 'opt_var' maps to '--usr_frg'.
  - For 'mode=ratio',  'opt_var' maps to '--dep_min'.

Examples
--------
  1. Build signal options with scaling and a fixed fragment length.
    '''bash
    set_args_opt signal 2 150 3
    '''

  2. Build ratio options with pseudocount, zero handling, and track output.
    '''bash
    set_args_opt ratio 2:1 0.5 4 true 1:2 0.001 post_scale true
    '''
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

    if [[ "${dp}" != "NA" ]]; then
        optional+=( --dp "${dp}" )
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

    # Return values as comma-separated list.
    ( IFS=','; echo "${optional[*]}" )
}


function run_dry_or_wet() {
    local arr_nam="${1:-}"
    local log_out="${2:-}"
    local log_err="${3:-}"
    local dir_out dir_err decl  # Derived directories and declaration metadata
    local -a cmd_cpy            # Local copy of the command array
    local show_help             # Help text


    show_help=$(cat << EOM
Usage
-----
  run_dry_or_wet
    [--help] arr_nam log_out log_err

  Print or execute a command stored in an array variable, with stdout/stderr redirected to log files.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  arr_nam : str
    Name of the command array.

  2  log_out : file
    Stdout log file.

  3  log_err : file
    Stderr log file.

Returns
-------
  0 for successful dry runs and the command exit status for wet runs; 1 if command-array or log-path validation fails.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - Command-array executable (when 'dry_run' is false)
    - dirname

  - In debug or dry-run mode, prints the fully quoted command plus redirections.
  - In non-dry-run mode, executes the command and returns its exit status.

Examples
--------
  1. Print a command in dry-run mode without executing it.
    '''bash
    tmp="\$(mktemp -d)"; trap 'rm -r -- "\${tmp}"' EXIT
    declare -a cmd=(printf '%s\n' dry-run)
    debug=false; dry_run=true
    run_dry_or_wet cmd "\${tmp}/stdout.txt" "\${tmp}/stderr.txt"
    '''

  2. Execute a command and append its output to writable log files.
    '''bash
    tmp="\$(mktemp -d)"; trap 'rm -r -- "\${tmp}"' EXIT
    declare -a cmd=(printf '%s\n' executed)
    debug=false; dry_run=false
    if run_dry_or_wet cmd "\${tmp}/stdout.txt" "\${tmp}/stderr.txt"; then
        printf '%s\n' 'command completed'
    fi
    '''
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

    # Refuse to run if log dirs are neither existent nor writable.
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

    # In "debug" or "dry-run mode", show the exact command and redirections.
    if [[ "${debug}" == "true" || "${dry_run}" == "true" ]]; then
        printf '%q ' "${cmd_cpy[@]}" >&2
        echo ">> ${log_out} 2>> ${log_err}" >&2
        echo >&2
        echo >&2
    fi

    # Execute the command with logging when not in dry-run mode.
    if [[ "${dry_run}" == "false" ]]; then
        "${cmd_cpy[@]}" >> "${log_out}" 2>> "${log_err}"
        return $?
    fi

    return 0
}


function run_comp_sig() {
    local debug="${1:-}"
    local threads="${2:-}"
    local fil_in="${3:-}"
    local fil_out="${4:-}"
    local siz_bin="${5:-}"
    local method="${6:-}"
    local scl_fct="${7:-}"
    local usr_frg="${8:-}"
    local dp="${9:-}"
    local ref_fa="${10:-}"
    local chr_sizes="${11:-}"
    local chunk_size="${12:-}"
    local engine="${13:-}"
    local dir_eo="${14:-}"
    local nam_job="${15:-}"
    local dsc="${16:-}"
    local log_out log_err  # Explicit local variable declarations
    local -a optional cmd  # Optional arguments and command array
    local show_help        # Help text

    show_help=$(cat << EOM
Usage
-----
  run_comp_sig
    [--help] debug threads fil_in fil_out siz_bin method scl_fct usr_frg dp ref_fa chr_sizes chunk_size engine dir_eo nam_job dsc

  Build and run the per-sample call to 'compute_signal.py'.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  01  debug : bool
    Print debug messages or not.

  02  threads : int
    Number of threads to use.

  03  fil_in : file
    Input file path. Input BAM or CRAM file.

  04  fil_out : file
    Output file path.

  05  siz_bin : int
    Bin size in base pairs.

  06  method : {'unadj', 'frag', 'norm'}
    Workflow method. Type of signal computation or empty sentinel ("").

  07  scl_fct : number
    Scaling factor. Use sentinel 'NA' to omit scaling.

  08  usr_frg : int
    Fixed fragment length or sentinel 'NA'.

  09  dp : int
    Maximum number of decimal places retained for finite emitted values.

  10  ref_fa : file
    Reference FASTA file for CRAM input or empty string.

  11  chr_sizes : file
    Chromosome sizes file or empty string.

  12  chunk_size : int
    Number of records to process per chunk. Used by compute-signal engines.

  13  engine : {'chrom', 'window'}
    Processing engine.

  14  dir_eo : dir
    Directory for stderr and stdout log files.

  15  nam_job : str
    Job name.

  16  dsc : str
    Descriptor for log file naming.

Returns
-------
  Returns 0 when the per-sample signal command succeeds; 1 otherwise.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - bash >= 4.4
    - dirname
    - Input BAM or CRAM index (when the selected signal engine uses indexed access)
    - python >= 3.11
    - Reference FASTA and required index (when processing CRAM)

  - Optional CLI arguments are derived via 'set_args_opt'.
  - Logging is handled through 'run_dry_or_wet'.

Examples
--------
  1. Display the signal-dispatch positional contract.
    '''bash
    run_comp_sig --help
    '''

  2. Build and print one signal command safely in dry-run mode.
    '''bash
    tmp="\$(mktemp -d)"; trap 'rm -r -- "\${tmp}"' EXIT
    scr_sig=compute_signal; debug=false; dry_run=true
    run_comp_sig \\
        false \\
        1 \\
        sample.bam \\
        sample.norm.bedGraph \\
        10 \\
        norm \\
        NA \\
        NA \\
        3 \\
        '' \\
        '' \\
        100000 \\
        chrom \\
        "\${tmp}" \\
        compute_signal \\
        sample_norm
    '''
EOM
    )

    if [[ "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ $# -ne 16 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'run_comp_sig()' expects 16 arguments, but got $#."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    # Generate optional arguments array dynamically.
    IFS="," read -r -a optional <<< "$(
        set_args_opt "signal" "${scl_fct}" "${usr_frg}" "${dp}"
    )"

    # Debug array of optional arguments.
    if [[ "${debug}" == "true" ]]; then
        echo "set_args_opt(): optional=( ${optional[*]} )" >&2
        echo >&2
    fi

    # Define paths for log output files.
    log_out="${dir_eo}/${nam_job}.${dsc}.stdout.txt"
    log_err="${dir_eo}/${nam_job}.${dsc}.stderr.txt"

    # Build call to 'compute_signal.py' via 'run_py'.
    cmd=(
        run_py "${scr_sig}"
            --verbose
            --threads "${threads}"
            --fil_in "${fil_in}"
            --fil_out "${fil_out}"
            --siz_bin "${siz_bin}"
    )

    if [[ -n "${method}" ]]; then
        cmd+=( --method "${method}" )
    fi

    if [[ -n "${ref_fa}" ]]; then
        cmd+=( --ref_fa "${ref_fa}" )
    fi

    if [[ -n "${chr_sizes}" ]]; then
        cmd+=( --chr_sizes "${chr_sizes}" )
    fi

    if [[ -n "${chunk_size}" ]]; then
        cmd+=( --chunk_size "${chunk_size}" )
    fi

    if [[ -n "${engine}" ]]; then
        cmd+=( --engine "${engine}" )
    fi

    if [[ "${#optional[@]}" -gt 0 && -n "${optional[0]}" ]]; then
        cmd+=( "${optional[@]}" )
    fi

    # Debug or execute call to 'compute_signal.py'.
    run_dry_or_wet cmd "${log_out}" "${log_err}" || return 1
}


function run_comp_rat() {
    local debug="${1:-}"
    local fil_A="${2:-}"
    local fil_B="${3:-}"
    local fil_out="${4:-}"
    local method="${5:-}"
    local scl_fct="${6:-}"
    local dep_min="${7:-}"
    local dp="${8:-}"
    local track="${9:-}"
    local pseudo="${10:-}"
    local eps="${11:-}"
    local skip_00="${12:-}"
    local drp_nan="${13:-}"
    local skp_pfx="${14:-}"
    local chr_sizes="${15:-}"
    local strict_bins="${16:-false}"
    local dir_eo="${17:-}"
    local nam_job="${18:-}"
    local dsc="${19:-}"
    local log_out log_err  # Local variable declarations
    local -a optional cmd  # Local array declarations
    local show_help        # Help text

    show_help=$(cat << EOM
Usage
-----
  run_comp_rat
    [--help] debug fil_A fil_B fil_out method scl_fct dep_min dp track pseudo eps skip_00 drp_nan skp_pfx chr_sizes strict_bins dir_eo nam_job dsc

  Build and run the per-sample call to 'compute_signal_ratio.py'.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  01  debug : bool
    Print debug messages or not.

  02  fil_A : file
    First bedGraph input file, file A.

  03  fil_B : file
    Second bedGraph input file, file B.

  04  fil_out : file
    Output file path. Output ratio bedGraph file.

  05  method : {'unadj', 'log2', 'unadj_r', 'log2_r'}
    Workflow method. Ratio method: 'unadj', 'log2', 'unadj_r', 'log2_r'.

  06  scl_fct : number
    Scaling factor. Use sentinel 'NA' to omit scaling.

  07  dep_min : float
    Optional minimum input depth or sentinel 'NA'.

  08  dp : int
    Maximum number of decimal places retained for finite emitted values.

  09  track : bool
    Write a companion track file. The name carries a '.track' suffix.

  10  pseudo : structured string
    Per-file pseudocount spec 'A[:B]' or sentinel 'NA'.

  11  eps : float
    Zero tolerance epsilon or sentinel 'NA'.

  12  skip_00 : {'pre_scale', 'post_scale'}
    Skip rows where both compared values are zero. Shared zero-zero skip mode ('pre_scale' or 'post_scale') or sentinel (NA).

  13  drp_nan : bool
    Drop non-finite values from main output.

  14  skp_pfx : str
    Comma-separated list of header prefixes to skip. Use sentinel 'NA' to omit.

  15  chr_sizes : file
    Chromosome sizes file or empty string.

  16  strict_bins : bool
    Require strict bin compatibility across identical ordered input bin grids.

  17  dir_eo : dir
    Directory for stderr and stdout log files.

  18  nam_job : str
    Job name used in log filenames.

  19  dsc : str
    Descriptor string for logs.

Returns
-------
  Returns 0 when the per-sample ratio command succeeds; 1 otherwise.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - bash >= 4.4
    - dirname
    - python >= 3.11

  - Optional CLI arguments are derived via 'set_args_opt'.
  - Shared header-skip prefixes are appended via '--skp_pfx' when set.
  - Logging is handled through 'run_dry_or_wet'.

Examples
--------
  1. Display the ratio-dispatch positional contract.
    '''bash
    run_comp_rat --help
    '''

  2. Build and print one ratio command safely in dry-run mode.
    '''bash
    tmp="\$(mktemp -d)"; trap 'rm -r -- "\${tmp}"' EXIT
    scr_rat=compute_signal_ratio; debug=false; dry_run=true
    run_comp_rat \\
        false \\
        IP.bedGraph \\
        input.bedGraph \\
        ratio.bedGraph \\
        unadj \\
        NA \\
        NA \\
        3 \\
        false \\
        NA \\
        0 \\
        NA \\
        false \\
        NA \\
        '' \\
        false \\
        "\${tmp}" \\
        compute_ratio \\
        IP_over_input
    '''
EOM
    )

    if [[ "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ $# -ne 19 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'run_comp_rat()' expects 17 arguments, but got $#."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    # Generate optional arguments array dynamically.
    IFS="," read -r -a optional <<< "$(
        set_args_opt \
            "ratio" \
            "${scl_fct}" \
            "${dep_min}" \
            "${dp}" \
            "${track}" \
            "${pseudo}" \
            "${eps}" \
            "${skip_00}" \
            "${drp_nan}"
    )"

    # Debug array of optional arguments.
    if [[ "${debug}" == "true" ]]; then
        echo "set_args_opt(): optional=( ${optional[*]} )" >&2
        echo >&2
    fi

    # Define paths for log output files.
    log_out="${dir_eo}/${nam_job}.${dsc}.stdout.txt"
    log_err="${dir_eo}/${nam_job}.${dsc}.stderr.txt"

    # Build call to 'compute_signal_ratio.py' via 'run_py'.
    cmd=(
        run_py "${scr_rat}"
            --verbose
            --fil_A "${fil_A}"
            --fil_B "${fil_B}"
            --fil_out "${fil_out}"
            --method "${method}"
    )

    # If any, then append optional flags.
    if [[ "${#optional[@]}" -gt 0 && -n "${optional[0]}" ]]; then
        cmd+=( "${optional[@]}" )
    fi

    if [[ "${skp_pfx}" != "NA" ]]; then
        cmd+=( --skp_pfx "${skp_pfx}" )
    fi

    if [[ -n "${chr_sizes}" ]]; then
        cmd+=( --chr_sizes "${chr_sizes}" )
    fi

    if [[ "${strict_bins}" == "true" ]]; then
        cmd+=( --strict_bins )
    fi

    # Debug or execute call to 'compute_signal_ratio.py'.
    run_dry_or_wet cmd "${log_out}" "${log_err}" || return 1
}


function get_arr_elem() {
    local arr_nam="${1:-}"  # Name of indexed array.
    local idx="${2:-}"      # Array index.
    local decl              # Output from 'declare -p'.
    local show_help         # Help message.

    show_help=$(cat << EOM
Usage
-----
  get_arr_elem
    [--help] arr_nam idx

  Print one element from a named indexed array.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  arr_nam : str
    Name of indexed array.

  2  idx : int
    Zero-based array index.

Returns
-------
  - Prints the selected array element to stdout.
  - Returns 0 when the element is found; 1 otherwise.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Print the second element of a prepared input array.
    '''bash
    declare -a arr_fil_in=(sample_A.bam sample_B.bam)
    get_arr_elem arr_fil_in 1
    '''

  2. Confirm that an out-of-range element request is rejected.
    '''bash
    declare -a arr_fil_out=(sample_A.bedGraph)
    if ! get_arr_elem arr_fil_out 2; then
        printf '%s\n' 'out-of-range index rejected as expected'
    fi
    '''
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
    local mode="${1:-}"
    local idx="${2:-}"
    local nam_arr_fil_A="${3:-}"
    local nam_arr_fil_B="${4:-}"
    local nam_arr_out="${5:-}"
    local nam_arr_scl="${6:-}"
    local nam_arr_opt="${7:-}"
    local fil_A fil_B fil_out scl_fct opt_var samp dsc
    local err_ini out_ini err_dsc out_dsc
    local show_help

    show_help=$(cat << EOM
Usage
-----
  task_pro
    [--help] mode idx arr_fil_A arr_fil_B arr_out arr_scl arr_opt

  Prepare per-task inputs and logging metadata for the 'run_task_*' helpers ("pro" is "prologue").

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  mode : {'signal', 'ratio', 'coord'}
    Workflow mode. 'signal', 'ratio', or 'coord'.

  2  idx : int
    Zero-based task index.

  3  nam_arr_fil_A : str
    Array name for scalar fil_in or file A.

  4  nam_arr_fil_B : str
    Array name for scalar file B; "" if unused.

  5  nam_arr_out : str
    Array name for scalar fil_out.

  6  nam_arr_scl : str
    Array name for scalar scl_fct; "" if unused.

  7  nam_arr_opt : str
    Array name for scalar opt_var; "" if unused.

Returns
-------
  - Prints one comma-delimited task metadata line to stdout:
    fil_A,fil_B?,fil_out,samp,dsc,err_ini?,out_ini?,err_dsc?,out_dsc?
  - Returns 0 when task metadata are prepared successfully; 1 otherwise.

Notes
-----
  Runtime requirements:
    - basename
    - bash >= 4.4
    - ln (when run as a Slurm array task)

  - Pulls per-task values from arrays by name.
  - Emits debug snapshots when 'debug=true'.
  - Calls 'process_io' to derive sample ('samp') and descriptor ('dsc').
  - Under Slurm, also derives initial log paths via 'set_logs_slurm'.

Examples
--------
  1. Prepare local metadata for one signal task.
    '''bash
    debug=false; unset SLURM_ARRAY_TASK_ID
    arr_fil_in=(sample.bam); arr_fil_out=(sample.norm.bedGraph)
    arr_scl_fct=(NA); arr_usr_frg=(NA)
    task_pro signal 0 arr_fil_in '' arr_fil_out arr_scl_fct arr_usr_frg
    '''

  2. Capture local metadata for one ratio task.
    '''bash
    debug=false; unset SLURM_ARRAY_TASK_ID
    arr_fil_A=(IP.bedGraph); arr_fil_B=(input.bedGraph)
    arr_fil_out=(ratio.bedGraph); arr_scl_fct=(2:1); arr_dep_min=(0.5)
    if \\
        metadata="\$(
            task_pro \\
                ratio 0 arr_fil_A arr_fil_B arr_fil_out arr_scl_fct arr_dep_min
        )"
    then
        printf '%s\n' "\${metadata}"
    fi
    '''
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

    fil_A="$(get_arr_elem "${nam_arr_fil_A}" "${idx}")" || return 1

    if [[ -n "${nam_arr_fil_B}" ]]; then
        fil_B="$(get_arr_elem "${nam_arr_fil_B}" "${idx}")" || return 1
    fi

    fil_out="$(get_arr_elem "${nam_arr_out}" "${idx}")" || return 1

    if [[ -n "${nam_arr_scl}" ]]; then
        scl_fct="$(get_arr_elem "${nam_arr_scl}" "${idx}")" || return 1
    fi

    if [[ -n "${nam_arr_opt}" ]]; then
        opt_var="$(get_arr_elem "${nam_arr_opt}" "${idx}")" || return 1
    fi

    # Debug inputs.
    if [[ "${debug}" == "true" ]]; then
        if [[ "${mode}" == "signal" ]]; then
            debug_var \
                "fil_in=${fil_A}" \
                "fil_out=${fil_out}" \
                "scl_fct=${scl_fct}" \
                "usr_frg=${opt_var}"
        elif [[ "${mode}" == "ratio" ]]; then
            debug_var \
                "fil_A=${fil_A}" \
                "fil_B=${fil_B}" \
                "fil_out=${fil_out}" \
                "scl_fct=${scl_fct}" \
                "dep_min=${opt_var}"
        else
            debug_var \
                "fil_in=${fil_A}" \
                "fil_out=${fil_out}" \
                "usr_frg=${opt_var}"
        fi
    fi

    # Derive sample ('samp') and descriptor ('dsc') via 'process_io'.
    if [[ "${mode}" == "ratio" ]]; then
        IFS=',' read -r samp dsc < <(
            process_io \
                -md "${mode}" \
                -fA "${fil_A}" \
                -fB "${fil_B}" \
                -fo "${fil_out}" \
                -sf "${scl_fct:-NA}" \
                -ov "${opt_var:-NA}"
        ) || return 1
    elif [[ "${mode}" == "signal" ]]; then
        IFS=',' read -r samp dsc < <(
            process_io \
                -md "${mode}" \
                -fi "${fil_A}" \
                -fo "${fil_out}" \
                -sf "${scl_fct:-NA}" \
                -ov "${opt_var:-NA}"
        ) || return 1
    else
        IFS=',' read -r samp dsc < <(
            process_io \
                -md "${mode}" \
                -fi "${fil_A}" \
                -fo "${fil_out}"
        ) || return 1
    fi

    # Debug 'samp' and 'dsc' output by process_io.
    if [[ "${debug}" == "true" ]]; then
        debug_var "samp=${samp}" "dsc=${dsc}"
    fi

    # If using Slurm, request initial and descriptor log paths.
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        IFS=',' read -r err_ini out_ini err_dsc out_dsc < <(
            set_logs_slurm \
                "${id_job}" "${id_tsk}" "${samp}" "${dir_eo}" "${nam_job}"
        ) || return 1

        if [[ "${debug}" == "true" ]]; then
            debug_var \
                "err_ini=${err_ini}" \
                "out_ini=${out_ini}" \
                "err_dsc=${err_dsc}" \
                "out_dsc=${out_dsc}"
        fi
    fi

    # Return all values for one of the task callers to parse.
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "${fil_A}" "${fil_B:-}" "${fil_out}" "${samp}" "${dsc}" \
        "${err_ini:-}" "${out_ini:-}" "${err_dsc:-}" "${out_dsc:-}"
}


function task_epi() {
    local show_help

    show_help=$(cat << EOM
Usage
-----
  task_epi
    [--help] err_ini out_ini

  Remove initial Slurm-generated log files after a task completes ("epi" is "epilogue").

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  err_ini : str
    Initial stderr log path.

  2  out_ini : str
    Initial stdout log path.

Returns
-------
  Returns 0 when initial log cleanup succeeds or no cleanup is needed; 1 when argument validation fails.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - rm (when run as a wet Slurm array task)

  - No action is taken outside Slurm array tasks.
  - No action is taken in dry-run mode.

Examples
--------
  1. Perform the no-op cleanup path outside a Slurm array task.
    '''bash
    unset SLURM_ARRAY_TASK_ID; dry_run=false
    task_epi '' ''
    '''

  2. Remove existing initial logs after a wet Slurm-array task.
    '''bash
    tmp="\$(mktemp -d)"; trap 'rm -r -- "\${tmp}"' EXIT
    touch "\${tmp}/initial.stderr.txt" "\${tmp}/initial.stdout.txt"
    SLURM_ARRAY_TASK_ID=0; dry_run=false
    task_epi "\${tmp}/initial.stderr.txt" "\${tmp}/initial.stdout.txt"
    '''
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
    local idx="${1:-}"
    local fil_in unused fil_out samp dsc err_ini out_ini err_dsc out_dsc
    local rc=0
    local show_help

    show_help=$(cat << EOM
Usage
-----
  run_task_sig
    [--help] idx

  Run one signal-computation task, either under Slurm-array execution or in GNU Parallel or in a serial iteration.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  idx : int
    Zero-based task index into 'arr_fil_in', 'arr_fil_out', 'arr_scl_fct', and 'arr_usr_frg'.

Returns
-------
  Returns the exit status from 'run_comp_sig', or 1 if task setup fails.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - basename
    - bash >= 4.4
    - dirname
    - Input BAM or CRAM index (when the selected signal engine uses indexed access)
    - ln (when run as a Slurm array task)
    - python >= 3.11
    - Reference FASTA and required index (when processing CRAM)
    - rm (when run as a wet Slurm array task)

  - Uses 'task_pro' to derive task-specific inputs and logging metadata.
  - Preserves and returns the exit status from 'run_comp_sig'.
  - Uses 'task_epi' to remove initial Slurm logs after execution.

Examples
--------
  1. Display the signal-task contract before wrapper setup.
    '''bash
    run_task_sig --help
    '''

  2. Run the first signal task after 'main' has prepared globals and arrays.
    '''bash
    run_task_sig 0
    '''
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
        fil_in unused fil_out samp dsc \
        err_ini out_ini err_dsc out_dsc \
    < <(
        task_pro "signal" "${idx}" \
            "arr_fil_in" "" "arr_fil_out" "arr_scl_fct" "arr_usr_frg"
    ) || return 1

    run_comp_sig \
        "${debug}" \
        "${threads}" \
        "${fil_in}" \
        "${fil_out}" \
        "${siz_bin}" \
        "${method}" \
        "$(get_arr_elem arr_scl_fct "${idx}")" \
        "$(get_arr_elem arr_usr_frg "${idx}")" \
        "${dp}" \
        "${ref_fa}" \
        "${chr_sizes}" \
        "${chunk_size}" \
        "${engine}" \
        "${dir_eo}" \
        "${nam_job}" \
        "${dsc}"
    rc=$?

    task_epi "${err_ini}" "${out_ini}"
    return "${rc}"
}


function run_task_rat() {
    local idx="${1:-}"
    local fil_A fil_B fil_out samp dsc err_ini out_ini err_dsc out_dsc
    local rc=0
    local show_help

    show_help=$(cat << EOM
Usage
-----
  run_task_rat
    [--help] idx

  Run one ratio-computation task, either under Slurm-array execution or in GNU Parallel or in a serial iteration.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  idx : int
    Zero-based task index into 'arr_fil_A', 'arr_fil_B', 'arr_fil_out', 'arr_scl_fct', 'arr_dep_min', and 'arr_pseudo'.

Returns
-------
  Returns the exit status from 'run_comp_rat', or 1 if task setup fails.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - basename
    - bash >= 4.4
    - dirname
    - ln (when run as a Slurm array task)
    - python >= 3.11
    - rm (when run as a wet Slurm array task)

  - Uses 'task_pro' to derive task-specific inputs and logging metadata.
  - Preserves and returns the exit status from 'run_comp_rat'.
  - Uses 'task_epi' to remove initial Slurm logs after execution.

Examples
--------
  1. Display the ratio-task contract before wrapper setup.
    '''bash
    run_task_rat --help
    '''

  2. Run the first ratio task after 'main' has prepared globals and arrays.
    '''bash
    run_task_rat 0
    '''
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
        fil_A fil_B fil_out samp dsc \
        err_ini out_ini err_dsc out_dsc \
    < <(
        task_pro "ratio" "${idx}" \
            "arr_fil_A" "arr_fil_B" "arr_fil_out" "arr_scl_fct" "arr_dep_min"
    ) || return 1

    run_comp_rat \
        "${debug}" \
        "${fil_A}" \
        "${fil_B}" \
        "${fil_out}" \
        "${method}" \
        "$(get_arr_elem arr_scl_fct "${idx}")" \
        "$(get_arr_elem arr_dep_min "${idx}")" \
        "${dp}" \
        "${track}" \
        "$(get_arr_elem arr_pseudo "${idx}")" \
        "${eps}" \
        "${skip_00}" \
        "${drp_nan}" \
        "${skp_pfx}" \
        "${chr_sizes}" \
        "${strict_bins}" \
        "${dir_eo}" \
        "${nam_job}" \
        "${dsc}"
    rc=$?

    task_epi "${err_ini}" "${out_ini}"
    return "${rc}"
}


function run_task_coord() {
    local idx="${1:-}"
    local fil_in unused fil_out samp dsc err_ini out_ini err_dsc out_dsc
    local rc=0
    local show_help

    show_help=$(cat << EOM
Usage
-----
  run_task_coord
    [--help] idx

  Run one fragment-coordinate extraction task, either under Slurm-array execution or in GNU Parallel or in a serial iteration.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  idx : int
    Zero-based task index into 'arr_fil_in', 'arr_fil_out', and 'arr_usr_frg'.

Returns
-------
  Returns the exit status from the coordinate extraction command, or 1 if task setup fails.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - basename
    - bash >= 4.4
    - dirname
    - ln (when run as a Slurm array task)
    - python >= 3.11
    - Reference FASTA and required index (when processing CRAM)
    - rm (when run as a wet Slurm array task)

  - Uses 'task_pro' to derive task-specific inputs and logging metadata.
  - Reuses 'run_comp_sig' with stub signal-mode arguments to preserve the original coordinate-extraction pathway.
  - Preserves and returns the exit status from 'run_comp_sig'.
  - Uses 'task_epi' to remove initial Slurm logs after execution.

Examples
--------
  1. Display the coordinate-task contract before wrapper setup.
    '''bash
    run_task_coord --help
    '''

  2. Run the first coordinate task after 'main' has prepared globals and arrays.
    '''bash
    run_task_coord 0
    '''
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
        fil_in unused fil_out samp dsc \
        err_ini out_ini err_dsc out_dsc \
    < <(
        task_pro "coord" "${idx}" \
            "arr_fil_in" "" "arr_fil_out" "" "arr_usr_frg"
    ) || return 1

    # (Use stub parameters per original 'coord' behavior).
    run_comp_sig \
        "${debug}" \
        1 \
        "${fil_in}" \
        "${fil_out}" \
        1 \
        "" \
        "NA" \
        "$(get_arr_elem arr_usr_frg "${idx}")" \
        1 \
        "${ref_fa}" \
        "${chr_sizes}" \
        "" \
        "" \
        "${dir_eo}" \
        "${nam_job}" \
        "${dsc}"
    rc=$?

    task_epi "${err_ini}" "${out_ini}"
    return "${rc}"
}



# Source 'source_helpers.sh' and requested helper scripts from 'dir_scr'.
function source_helpers_submit() {
    local scr="${1:-}"
    local dir_scr_arg="${2:-}"
    local fnc_src

    if (( $# < 2 )); then
        echo "error(${scr:-unknown_script}):" \
            "expected at least two arguments: 'scr' and 'dir_scr_arg'." >&2
        return 1
    fi

    shift 2

    if [[ -z "${scr}" ]]; then
        scr="unknown_script"
    fi

    if [[ -z "${dir_scr_arg}" ]]; then
        echo "error(${scr}):" \
            "positional argument 2, 'dir_scr_arg', is missing." >&2
        return 1
    elif [[ ! -d "${dir_scr_arg}" ]]; then
        echo "error(${scr}):" \
            "script directory not found: '${dir_scr_arg}'." >&2
        return 1
    elif (( $# < 1 )); then
        echo "error(${scr}):" \
            "at least one helper script name must be supplied." >&2
        return 1
    fi

    dir_fnc="${dir_scr_arg}/../lib/bash"
    fnc_src="${dir_fnc}/core/source_helpers.sh"

    if [[ ! -f "${fnc_src}" ]]; then
        echo "error(${scr}):" \
            "script not found: '${fnc_src}'." >&2
        return 1
    fi

    # shellcheck disable=SC1090
    source "${fnc_src}" || {
        echo "error(${scr}):" \
            "failed to source '${fnc_src}'." >&2
        return 1
    }

    source_helpers "${dir_fnc}" "$@" || {
        echo "error(${scr}):" \
            "failed to source required helper scripts." >&2
        return 1
    }
}


#  Initialize hardcoded arguments
function init_args_hardcoded() {
    # WARNING: Do not change unless testing/stepping through. If true, print
    # verbose/debug Bash-level logging.
    debug=true

    # If true, dry-run script.
    dry_run=false

    # If true, parse arguments and exit before validation or job submission.
    p_only=false

    # If true, parse and check arguments, then exit before execution.
    pc_only=false
}


# Initialize argument variables, assigning default values where applicable.
function init_arg_defs() {
    env_nam="env_protocol"
    threads=4
    mode="signal"
    method=""
    csv_fil_in=""
    csv_fil_A=""
    csv_fil_B=""
    csv_fil_out=""
    ref_fa=""
    chr_sizes=""
    track=false
    siz_bin=10
    chunk_size=100000
    engine="chrom"
    csv_scl_fct=""
    csv_usr_frg=""
    csv_dep_min=""
    csv_pseudo=""
    eps="NA"
    skip_00="NA"
    strict_bins=false
    drp_nan=false
    skp_pfx="NA"
    dp=24
    dir_eo=""
    nam_job=""
}


# Initialize hardcoded arguments and user-facing argument defaults.
function init_defs() {
    init_args_hardcoded
    init_arg_defs
}


# Parse keyword arguments after helper scripts have been sourced.
function parse_args() {
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -en|--env|--env[_-]nam)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                env_nam="${2}"
                shift 2
                ;;

            -ds|--dir[_-]scr)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                dir_scr="${2}"
                shift 2
                ;;

            -t|--thr|--threads)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                threads="${2}"
                shift 2
                ;;

            -md|--mode)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                mode="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -me|--method)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                method="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -ci|--csv[_-]fil_in)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                csv_fil_in="${2}"
                shift 2
                ;;

            -rf|--ref[_-]fa)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                ref_fa="${2}"
                shift 2
                ;;

            -cs|--chr[_-]sizes|--chrom[_-]sizes)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                chr_sizes="${2}"
                shift 2
                ;;

            -cA|--csv[_-]fil[_-]A)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                csv_fil_A="${2}"
                shift 2
                ;;

            -cB|--csv[_-]fil[_-]B)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                csv_fil_B="${2}"
                shift 2
                ;;

            -co|--csv[_-]fil_out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                csv_fil_out="${2}"
                shift 2
                ;;

            -tr|--track)
                track=true
                shift 1
                ;;

            -sb|--siz[_-]bin)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                siz_bin="${2}"
                shift 2
                ;;

            -ck|--chunk[_-]size|--chnk[_-]size)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                chunk_size="${2}"
                shift 2
                ;;

            -eg|--engine)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                engine="${2,,}"
                shift 2
                ;;

            -csf|--csv[_-]scl[_-]fct)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                csv_scl_fct="${2}"
                shift 2
                ;;

            -cuf|--csv[_-]usr[_-]frg)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                csv_usr_frg="${2}"
                shift 2
                ;;

            -cdm|--csv[_-]dep[_-]min)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                csv_dep_min="${2}"
                shift 2
                ;;

            -cps|--csv[_-]pseudo)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                csv_pseudo="${2}"
                shift 2
                ;;

            -e|--eps)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                eps="${2}"
                shift 2
                ;;

            -s0|--skip[_-]00)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                skip_00="${2}"
                shift 2
                ;;

            --strict[_-]bins)
                strict_bins=true
                shift 1
                ;;

            -dn|--drp[_-]nan)
                drp_nan=true
                shift 1
                ;;

            -sp|--skp[_-]pfx)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                skp_pfx="${2}"
                shift 2
                ;;

            -dp|--dp)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                dp="${2}"
                shift 2
                ;;

            -deo|--dir[_-]eo)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                dir_eo="${2}"
                shift 2
                ;;

            -nj|--nam[_-]job)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_compute_signal
                    return 1
                }
                nam_job="${2}"
                shift 2
                ;;

            *)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                help_submit_compute_signal
                return 1
                ;;
        esac
    done
}


# Canonicalize mode, method, and derived job name.
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
                        "invalid value for '--method': '${method}'. Expected" \
                        "a ratio method alias for 'unadj', 'log2'," \
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
                "Expected Workflow mode. 'signal', 'ratio', or 'coord'."
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


# Validate required arguments and simple numeric values.
function validate_args() {
    validate_var     "env_nam" "${env_nam}"         || return 1
    validate_var_dir "dir_scr" "${dir_scr}" 0 false || return 1
    validate_var     "threads" "${threads}"         || return 1

    check_int_pos "${threads}" "threads" || return 1
    check_int_pos "${dp}"     "dp"     || return 1

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
        validate_var "csv_fil_in"  "${csv_fil_in}"  || return 1
        validate_var "csv_fil_out" "${csv_fil_out}" || return 1
        validate_var "siz_bin"     "${siz_bin}"     || return 1
        check_int_pos "${siz_bin}" "siz_bin"        || return 1
        check_int_pos "${chunk_size}" "chunk_size"  || return 1
        case "${engine}" in
            chrom|window) : ;;
            *)
                echo_err \
                    "'--engine' must be 'chrom' or 'window': '${engine}'."
                return 1
                ;;
        esac
    elif [[ "${mode}" == "coord" ]]; then
        validate_var "csv_fil_in"  "${csv_fil_in}"  || return 1
        validate_var "csv_fil_out" "${csv_fil_out}" || return 1
    elif [[ "${mode}" == "ratio" ]]; then
        validate_var "csv_fil_A"   "${csv_fil_A}"   || return 1
        validate_var "csv_fil_B"   "${csv_fil_B}"   || return 1
        validate_var "csv_fil_out" "${csv_fil_out}" || return 1
    fi

    validate_var "nam_job" "${nam_job}" || return 1

    if [[ -n "${chr_sizes}" ]]; then
        validate_var_file "chr_sizes" "${chr_sizes}" 0 true || return 1
    fi

    if [[ -z "${dir_eo}" ]]; then
        echo_err "'--dir_eo' is required."
        return 1
    fi

    validate_var_dir "dir_eo" "${dir_eo}" || return 1
}


# Resolve installed Python command names after helper sourcing.
function resolve_paths_scrs() {
    validate_var_dir "dir_scr" "${dir_scr}" 0 false || return 1

    # 'run_python.sh' reads the caller-owned repository root dynamically.
    # shellcheck disable=SC2034
    dir_rep="$(cd "${dir_scr}/.." && pwd)"
    scr_sig=compute_signal
    scr_rat=compute_signal_ratio

    if [[ "${PY_INVOKE}" == "console" ]]; then
        check_pgrm_path "${scr_sig}" || return 1
        check_pgrm_path "${scr_rat}" || return 1
    fi
}


# Print parsed scalar arguments when debugging is enabled.
function print_state_debug() {
    if [[ "${debug}" != "true" ]]; then
        return 0
    fi

    debug_var \
        "env_nam=${env_nam}" \
        "dir_scr=${dir_scr}" \
        "threads=${threads}" \
        "mode=${mode}" \
        "method=${method:-UNSET}"

    if [[ "${mode}" != "ratio" ]]; then
        debug_var \
            "csv_fil_in=${csv_fil_in}" \
            "ref_fa=${ref_fa:-UNSET}"
    elif [[ "${mode}" == "ratio" ]]; then
        debug_var \
            "csv_fil_A=${csv_fil_A}" \
            "csv_fil_B=${csv_fil_B}"
    fi

    debug_var "csv_fil_out=${csv_fil_out}"
    debug_var "chr_sizes=${chr_sizes:-UNSET}"

    if [[ "${mode}" == "ratio" ]]; then
        debug_var \
            "track=${track}" \
            "csv_dep_min=${csv_dep_min}" \
            "csv_pseudo=${csv_pseudo}" \
            "eps=${eps}" \
            "skip_00=${skip_00}" \
            "strict_bins=${strict_bins}" \
            "drp_nan=${drp_nan}" \
            "skp_pfx=${skp_pfx}"
    fi

    if [[ "${mode}" == "signal" ]]; then
        debug_var \
            "siz_bin=${siz_bin}" \
            "chunk_size=${chunk_size}" \
            "engine=${engine}" \
            "csv_usr_frg=${csv_usr_frg}"
    fi

    if [[ "${mode}" != "coord" ]]; then
        debug_var \
            "csv_scl_fct=${csv_scl_fct}" \
            "dp=${dp}"
    fi

    debug_var \
        "dir_eo=${dir_eo}" \
        "nam_job=${nam_job}"
}


# Reconstruct mode-dependent arrays from serialized argument strings.
function prepare_vecs() {
    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        IFS=',' read -r -a arr_fil_in <<< "${csv_fil_in}"
        IFS=',' read -r -a arr_fil_out <<< "${csv_fil_out}"

        check_arr_nonempty "arr_fil_in"  "csv_fil_in"  || return 1
        check_arr_nonempty "arr_fil_out" "csv_fil_out" || return 1

        if [[ -n "${csv_usr_frg}" ]]; then
            IFS=',' read -r -a arr_usr_frg <<< "${csv_usr_frg}"
        else
            unset arr_usr_frg && declare -ga arr_usr_frg
            populate_array_empty arr_usr_frg "${#arr_fil_in[@]}"
        fi

        if [[ "${mode}" == "signal" ]]; then
            if [[ -n "${csv_scl_fct}" ]]; then
                IFS=',' read -r -a arr_scl_fct <<< "${csv_scl_fct}"
            else
                unset arr_scl_fct && declare -ga arr_scl_fct
                populate_array_empty arr_scl_fct "${#arr_fil_in[@]}"
            fi
        fi
    elif [[ "${mode}" == "ratio" ]]; then
        IFS=',' read -r -a arr_fil_A   <<< "${csv_fil_A}"
        IFS=',' read -r -a arr_fil_B   <<< "${csv_fil_B}"
        IFS=',' read -r -a arr_fil_out <<< "${csv_fil_out}"

        check_arr_nonempty "arr_fil_A"   "csv_fil_A"   || return 1
        check_arr_nonempty "arr_fil_B"   "csv_fil_B"   || return 1
        check_arr_nonempty "arr_fil_out" "csv_fil_out" || return 1

        if [[ -n "${csv_scl_fct}" ]]; then
            IFS=',' read -r -a arr_scl_fct <<< "${csv_scl_fct}"
        else
            unset arr_scl_fct && declare -ga arr_scl_fct
            populate_array_empty arr_scl_fct "${#arr_fil_A[@]}"
        fi

        if [[ -n "${csv_dep_min}" ]]; then
            IFS=',' read -r -a arr_dep_min <<< "${csv_dep_min}"
        else
            unset arr_dep_min && declare -ga arr_dep_min
            populate_array_empty arr_dep_min "${#arr_fil_A[@]}"
        fi

        if [[ -n "${csv_pseudo}" ]]; then
            IFS=',' read -r -a arr_pseudo <<< "${csv_pseudo}"
        else
            unset arr_pseudo && declare -ga arr_pseudo
            populate_array_empty arr_pseudo "${#arr_fil_A[@]}"
        fi
    fi
}


# Validate vector lengths, sentinels, CRAM reference needs, and input files.
function validate_vecs() {
    local idx fil_in fil_out need_ref

    if [[ "${mode}" == "signal" ]]; then
        check_arr_lengths "arr_fil_in" "arr_fil_out" || return 1
        check_arr_lengths "arr_fil_in" "arr_scl_fct" || return 1
        check_arr_lengths "arr_fil_in" "arr_usr_frg" || return 1
    elif [[ "${mode}" == "ratio" ]]; then
        check_arr_lengths "arr_fil_A"  "arr_fil_B"   || return 1
        check_arr_lengths "arr_fil_A"  "arr_fil_out" || return 1
        check_arr_lengths "arr_fil_A"  "arr_scl_fct" || return 1
        check_arr_lengths "arr_fil_A"  "arr_dep_min" || return 1
        check_arr_lengths "arr_fil_A"  "arr_pseudo"  || return 1
    elif [[ "${mode}" == "coord" ]]; then
        check_arr_lengths "arr_fil_in" "arr_fil_out" || return 1
        check_arr_lengths "arr_fil_in" "arr_usr_frg" || return 1
    fi

    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        for fil_in in "${arr_fil_in[@]}"; do
            if [[ "${fil_in}" == "-" ]]; then
                echo_err \
                    "'-' is not allowed in '--csv_fil_in' for" \
                    "'submit_compute_signal.sh'. Provide file paths."
                return 1
            fi
        done
    elif [[ "${mode}" == "ratio" ]]; then
        for fil_in in "${arr_fil_A[@]}"; do
            if [[ "${fil_in}" == "-" ]]; then
                echo_err \
                    "'-' is not allowed in '--csv_fil_A' for" \
                    "'submit_compute_signal.sh'. Provide file paths."
                return 1
            fi
        done

        for fil_in in "${arr_fil_B[@]}"; do
            if [[ "${fil_in}" == "-" ]]; then
                echo_err \
                    "'-' is not allowed in '--csv_fil_B' for" \
                    "'submit_compute_signal.sh'. Provide file paths."
                return 1
            fi
        done
    fi

    for fil_out in "${arr_fil_out[@]}"; do
        if [[ "${fil_out}" == "-" ]]; then
            echo_err \
                "'-' is not allowed in '--csv_fil_out' for" \
                "'submit_compute_signal.sh'. Provide output file paths."
            return 1
        fi
    done

    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        need_ref=false

        for fil_in in "${arr_fil_in[@]}"; do
            if [[ "${fil_in,,}" == *.cram ]]; then
                need_ref=true
                break
            fi
        done

        if [[ "${need_ref}" == "true" && -z "${ref_fa}" ]]; then
            echo_err \
                "'--ref_fa' is required when '--csv_fil_in' contains CRAM" \
                "input."
            return 1
        fi

        if [[ -n "${ref_fa}" ]]; then
            validate_var_file "ref_fa" "${ref_fa}" 0 true || return 1
        fi
    fi

    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        for idx in "${!arr_fil_in[@]}"; do
            validate_var_file "arr_fil_in" "${arr_fil_in[${idx}]}" "${idx}" \
                || return 1
        done
    elif [[ "${mode}" == "ratio" ]]; then
        for idx in "${!arr_fil_A[@]}"; do
            validate_var_file "arr_fil_A" "${arr_fil_A[${idx}]}" "${idx}" \
                || return 1
        done

        for idx in "${!arr_fil_B[@]}"; do
            validate_var_file "arr_fil_B" "${arr_fil_B[${idx}]}" "${idx}" \
                || return 1
        done
    fi
}


# Print reconstructed array values when debugging is enabled.
function print_vecs_debug() {
    if [[ "${debug}" != "true" ]]; then
        return 0
    fi

    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        debug_var \
            "\${#arr_fil_in[@]}=${#arr_fil_in[@]}" \
            "\${#arr_usr_frg[@]}=${#arr_usr_frg[@]}"
        echo "arr_fil_in=( ${arr_fil_in[*]} )"   >&2 && echo >&2
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

    debug_var "\${#arr_fil_out[@]}=${#arr_fil_out[@]}"
    echo "arr_fil_out=( ${arr_fil_out[*]} )"     >&2 && echo >&2
}


# Activate the environment containing the managed editable installation.
function setup_env() {
    handle_env "${env_nam}" || return 1
}


# Dispatch Slurm-array, GNU Parallel, or serial work.
function run_jobs() {
    local id_job id_tsk idx

    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        id_job=${SLURM_ARRAY_JOB_ID}
        id_tsk=${SLURM_ARRAY_TASK_ID}

        if ! [[ "${id_tsk}" =~ ^[1-9][0-9]*$ ]]; then
            echo_err "Slurm task ID is invalid: '${id_tsk}'."
            return 1
        elif [[
            "${mode}" == "ratio" && "${id_tsk}" -gt "${#arr_fil_A[@]}"
        ]]; then
            echo_err \
                "Slurm task ID '${id_tsk}' exceeds number of ratio entries:" \
                "'${#arr_fil_A[@]}'."
            return 1
        elif [[
            "${mode}" != "ratio" && "${id_tsk}" -gt "${#arr_fil_in[@]}"
        ]]; then
            echo_err \
                "Slurm task ID '${id_tsk}' exceeds number of input entries:" \
                "'${#arr_fil_in[@]}'."
            return 1
        else
            idx=$(( id_tsk - 1 ))
        fi

        if [[ "${mode}" == "signal" ]]; then
            run_task_sig "${idx}"   || return 1
        elif [[ "${mode}" == "ratio" ]]; then
            run_task_rat "${idx}"   || return 1
        else
            run_task_coord "${idx}" || return 1
        fi
    else
        if [[ "${mode}" == "signal" ]]; then
            for idx in "${!arr_fil_in[@]}"; do
                run_task_sig "${idx}"   || return 1
            done
        elif [[ "${mode}" == "ratio" ]]; then
            for idx in "${!arr_fil_A[@]}"; do
                run_task_rat "${idx}"   || return 1
            done
        else
            for idx in "${!arr_fil_in[@]}"; do
                run_task_coord "${idx}" || return 1
            done
        fi
    fi
}


# Main script execution.
function main() {
    init_defs

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        help_submit_compute_signal
        echo >&2
        return 0
    fi

    source_helpers_submit "${0##*/}" "${dir_scr}" \
        check_args \
        check_env \
        check_inputs \
        check_numbers \
        format_outputs \
        handle_env \
        manage_slurm \
        populate_array_empty \
        run_python \
        help/help_submit_compute_signal \
        || return 1

    parse_args "$@" || return 1

    if [[ "${p_only}" == "true" ]]; then
        if [[ "${debug}" == "true" ]]; then debug_var "p_only=true"; fi
        return 0
    fi

    canonicalize_args  || return 1
    validate_args      || return 1
    resolve_paths_scrs || return 1
    print_state_debug  || return 1
    prepare_vecs       || return 1
    validate_vecs      || return 1
    print_vecs_debug   || return 1

    if [[ "${pc_only}" == "true" ]]; then
        if [[ "${debug}" == "true" ]]; then debug_var "pc_only=true"; fi
        return 0
    fi

    setup_env || return 1
    run_jobs  || return 1
}


main "$@"
