#!/bin/bash

#  submit_compute_signal.sh
#  KA


#  Set hardcoded arguments ----------------------------------------------------
#  Run script in interactive mode (true) or command-line mode (false)
interactive=false

#  If true, run script in debug mode
debug=true

#  If true, dry-run script
dry_run=false

#  If true, run script in unit-test mode
unit=false


#  Define script-specific functions -------------------------------------------
#  Function to debug, validate, and parse input and output files, returning
#+ a sample name and output file descriptor
function process_io() {
    local mode=""        # Computation mode: 'signal', 'ratio', or 'coord'
    local infile_1=""    # 1° input: BAM (signal/coord) or IP bedGraph (ratio)
    local infile_2=""    # 2° input: input bedGraph (ratio only)
    local outfile=""     # Output file path (bedGraph/bedGraph.gz)
    local scl_fct=""     # Optional scaling factor, e.g., 'alpha' or 'spike'
    local opt_var=""     # Optional: 'frg_len' (signal), 'dep_min' (ratio)
    local help samp dsc  # Help text; sample name; output descriptor

    help=$(cat << EOM
Usage:
  process_io [--help] --mode <str> --infile_1 <str> [--infile_2 <str>] --outfile <str> [--scl_fct <flt>] [--opt_var <mult>]

Description:
  Validate and parse input/output file arguments for downstream processing, returning a sample name and output descriptor.

Arguments:
   -h, --help      <flg>  Print this help message and return 0.
   -m, --mode      <str>  Computation mode: 'signal', 'ratio', or 'coord'.
  -i1, --infile_1  <str>  Input file 1: BAM ('mode=signal' or 'mode=coord') or IP bedGraph ('mode=ratio').
  -i2, --infile_2  <str>  Input file 2: input bedGraph ('mode=ratio'; ignored otherwise).
   -o, --outfile   <str>  Output file.
  -sf, --scl_fct   <flt>  Optional scaling factor ('mode=signal' or 'mode=ratio').
  -ov, --opt_var   <mlt>  Optional variable: fragment length (<int>, 'mode=signal') or minimum input depth (<flt>, 'mode=ratio').
EOM
    )

    if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${help}"
        return 0
    fi

    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -m|--mode)     mode="${2}";      shift 2  ;;
            -i1|--infile_1) infile_1="${2}";  shift 2  ;;
            -i2|--infile_2) infile_2="${2}" ; shift 2  ;;
             -o|--outfile)  outfile="${2}";   shift 2  ;;
            -sf|--scl_fct)  scl_fct="${2}";   shift 2  ;;
            -ov|--opt_var)  opt_var="${2}";   shift 2  ;;
            *)
                echo "## Unknown parameter passed: '${1}' ##" >&2
                echo "" >&2
                echo "${help}" >&2
                return 1
                ;;
        esac
    done

    if [[ -z "${mode}" ]]; then
        echo "Error: '--mode' is required." >&2
        echo "" >&2
        echo "${help}" >&2
        return 1
    fi

    validate_var "infile_1" "${infile_1}" || return 1
    validate_var "outfile"  "${outfile}"  || return 1
    
    if [[ "${mode}" != "coord" ]]; then
        validate_var "scl_fct"  "${scl_fct}"  || return 1
        validate_var "opt_var"  "${opt_var}"  || return 1
    fi

    if [[ "${mode}" == "ratio" ]]; then
        validate_var "infile_2" "${infile_2}" || return 1
        samp="${outfile##*/}"
        for ext in bedgraph.gz bedgraph bdg.gz bdg bg.gz bg; do
            samp="${samp%."${ext}"}"
        done
    else
        samp="$(basename "${infile_1}" ".bam")"
    fi

    dsc="$(basename "${outfile}")"
    for ext in bedgraph.gz bedgraph bdg.gz bdg bg.gz bg; do
        dsc="${dsc%."${ext}"}"
    done
    unset ext

    echo "${samp},${dsc}"
}


#  Function to set optional arguments for the call to 'scr_sig'
function set_args_opt() {
    local mode="${1}"          # Mode: 'signal' or 'ratio'
    local scl_fct="${2}"       # Scaling factor/coefficient
    local opt_var="${3}"       # 'usr_frg' ('signal') or 'dep_min' ('ratio')
    local rnd="${4}"           # Number of decimal places for rounding
    local track="${5:-false}"  # Mode 'ratio': return track sans '-inf', 'nan'
    local log2="${6:-false}"   # Mode 'ratio': log2-transform ratio values
    
    unset optional && typeset -a optional

    if [[ $# -lt 4 ]]; then
        echo \
            "Error: set_args_opt() expects at least 4 arguments, but got" \
            "$#." >&2
        return 1
    elif [[ $# -gt 6 ]]; then
        echo \
            "Error: set_args_opt() expects at most 6 arguments, but got" \
            "$#." >&2
        return 1
    fi

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

    if [[ "${mode}" == "ratio" && "${track}" == "true" ]]; then
        optional+=( --track )
    fi

    if [[ "${mode}" == "ratio" && "${log2}" == "true" ]]; then
        optional+=( --log2 )
    fi

    #  Return values (comma-delimited)
    IFS=','; echo "${optional[*]}"
}


######################
## Work in progress ##
######################

#  Helper function: When debugging or dry-running, print an array-built
#+ command; otherwise execute it
function run_dry_or_wet() {
    local arr_nam="${1}"   # Name of the command array
    local log_out="${2}"   # Path to file for stdout redirection
    local log_err="${3}"   # Path to file for stderr redirection
    local dir_out dir_err  # Derived dirs for std(out|err) files
    local -a commnd        # Local copy of the command array
    local arr_len=0        # Preset local array length

    #  Check that array exists
    if ! eval 'declare -p '"${arr_nam}"' >/dev/null 2>&1'; then
        echo "Error: command array '${arr_nam}' is unset." >&2
        return 1
    fi

    #  Get array length and refuse to run if empty
    eval 'arr_len=${#'"${arr_nam}"'[@]}'
    if (( arr_len == 0 )); then
        echo "## WARNING: run_dry_or_wet received an empty command array ##" >&2
        return 1
    fi

    #  Make a safe local copy of the array, preserving element boundaries
    eval 'commnd=( "${'"${arr_nam}"'[@]}" )'

    #  Refuse to run if log dirs are neither existent nor writable
    dir_out="$(dirname "${log_out}")"
    dir_err="$(dirname "${log_err}")"

    if [[ ! -d "${dir_out}" || ! -d "${dir_err}" ]]; then
        echo "Error: Log directory missing: '${dir_out}' or '${dir_err}'." >&2
        return 1
    fi

    if [[ ! -w "${dir_out}" || ! -w "${dir_err}" ]]; then
        echo "Error: Log directory not writable: '${dir_out}' or '${dir_err}'." >&2
        return 1
    fi

    #  In debug or dry-run mode, show the exact command and redirections
    if [[ "${debug}" == "true" || "${dry_run}" == "true" ]]; then
        printf '%q ' "${commnd[@]}" >&2
        echo ">> ${log_out} 2>> ${log_err}" >&2
        echo "" >&2
        echo "" >&2
    fi

    #  Execute the command with logging when not in dry-run mode
    if [[ "${dry_run}" == "false" ]]; then
        "${commnd[@]}" >> "${log_out}" 2>> "${log_err}"
        return $?
    fi

    return 0
}

######################
## Work in progress ##
######################


#  Helper function to run 'compute_signal.py'
function run_comp_signal() {
    local debug="${1}"     # Print 'debug' statements: 'true' or 'false'
    local threads="${2}"   # Number of threads to use
    local infile="${3}"    # Input BAM file
    local outfile="${4}"   # Output filename
    local siz_bin="${5}"   # Bin size in base pairs
    local method="${6}"    # Type of signal computation
    local scl_fct="${7}"   # Scaling factor/coefficient
    local usr_frg="${8}"   # Fragment length
    local rnd="${9}"       # Number of decimal places for rounding
    local err_out="${10}"  # Directory for stdout and stderr
    local nam_job="${11}"  # Job name
    local dsc="${12}"      # Descriptor for log file naming
    local log_out log_err  # Explicit local variable declarations
    local -a optional cmd  # Explicit local array declarations

    #  Generate optional arguments array dynamically
    IFS="," read -r -a optional <<< "$(
        set_args_opt "signal" "${scl_fct}" "${usr_frg}" "${rnd}"
    )"

    #  Debug array of optional arguments
    if [[ "${debug}" == "true" ]]; then
        echo "set_args_opt(): optional=( ${optional[*]} )" >&2
        echo "" >&2
    fi

    #  Define paths for log output files
    log_out="${err_out}/${nam_job}.${dsc}.stdout.txt"
    log_err="${err_out}/${nam_job}.${dsc}.stderr.txt"

    #  Build call to 'compute_signal.py'
    cmd=(
        python "${scr_sig}"
            --verbose
            --threads "${threads}"
            --infile "${infile}"
            --outfile "${outfile}"
            --siz_bin "${siz_bin}"
    )
    if [[ -n "${method}" ]]; then cmd+=( --method "${method}" ); fi
    cmd+=( "${optional[@]}" )

    #  Debug or execute call to 'compute_signal.py'
    run_dry_or_wet cmd "${log_out}" "${log_err}" || return 1  # Fail fast
}


#  Helper function to run 'compute_signal_ratio.py'
function run_comp_ratio() {
    local debug="${1}"     # Print debug messages if 'true'; suppress if 'false'
    local fil_ip="${2}"    # Input IP bedGraph file (numerator in ratio)
    local fil_in="${3}"    # Input input bedGraph file (denominator in ratio)
    local outfile="${4}"   # Output ratio bedGraph file
    local scl_fct="${5}"   # Optional scaling factor applied to the IP signal before ratio computation
    local dep_min="${6}"   # Optional minimum input depth to avoid division by near-zero values
    local rnd="${7}"       # Number of decimal places for rounding ratio values
    local track="${8}"     # If 'true', write an additional bedGraph without '-inf' and 'nan' rows
    local log2="${9}"      # If 'true', compute log2(IP/input) instead of untransformed ratio
    local err_out="${10}"  # Directory for log file output (stderr, stdout)
    local nam_job="${11}"  # Job name used in log filenames
    local dsc="${12}"      # Descriptor string for logs
    local log_out log_err  # Explicit local variable declarations
    local -a optional cmd  # Explicit local array declarations

    #  Generate optional arguments array dynamically
    IFS="," read -r -a optional <<< "$(
        set_args_opt "ratio" \
            "${scl_fct}" "${dep_min}" "${rnd}" "${track}" "${log2}"
    )"

    #  Debug array of optional arguments
    if [[ "${debug}" == "true" ]]; then
        echo "set_args_opt(): optional=( ${optional[*]} )" >&2
        echo "" >&2
    fi

    #  Define paths for log output files
    log_out="${err_out}/${nam_job}.${dsc}.stdout.txt"
    log_err="${err_out}/${nam_job}.${dsc}.stderr.txt"

    #  Build call to 'compute_signal_ratio.py'
    cmd=(
        python "${scr_rat}"
            --verbose
            --fil_ip "${fil_ip}"
            --fil_in "${fil_in}"
            --fil_out "${outfile}"
            "${optional[@]}"
    )

    #  Debug or execute call to 'compute_signal_ratio.py'
    run_dry_or_wet cmd "${log_out}" "${log_err}" || return 1  # Fail fast
}


######################
## Work in progress ##
######################

#  Function to print an element from an indexed array by its name and index to
#+ stdout
function get_arr_elem() {
    # ${1}: array name
    # ${2}: index
    eval 'declare -p '"${1}"' >/dev/null 2>&1' || {
        echo "Error: array '${1}' is unset." >&2
        return 1
    }
    eval 'printf "%s\n" "${'"${1}"'['"${2}"']}"'
}


#  Prepare per-task inputs and logging info used by all run_task_* helper
#+ functions
#+   - Pulls array elements by name (infile(s), outfile, optional vars)
#+   - Emits debug snapshots
#+   - Calls process_io to derive sample (samp) and descriptor (dsc)
#+   - If under SLURM, obtains init/descriptor log paths via set_logs_slurm
#+ 
#+ Returns a comma-separated line:
#+   infile_1,infile_2?,outfile,samp,dsc,err_ini?,out_ini?,err_dsc?,out_dsc?
function task_prelude() {
    local mode="${1}"     # Mode: 'signal', 'ratio', or 'coord'
    local idx="${2}"      # Array index (integer ≥ 0)
    local arr_in1="${3}"  # Array name: 'infile_1'
    local arr_in2="${4}"  # Array name: 'infile_2'
    local arr_out="${5}"  # Array name: 'outfile'
    local arr_scl="${6}"  # Array name: 'scl_fct'
    local arr_opt="${7}"  # Array name: 'opt_var' ("" if unused)
    local infile_1 infile_2 outfile scl_fct opt_var samp dsc 
    local err_ini out_ini err_dsc out_dsc

    infile_1="$(get_arr_elem "${arr_in1}"   "${idx}")"

    if [[ -n "${arr_in2}"  ]]; then
       infile_2="$(get_arr_elem "${arr_in2}"  "${idx}")"
    fi

    outfile="$(get_arr_elem "${arr_out}" "${idx}")"

    if [[ -n "${arr_scl}" ]] ; then
       scl_fct="$(get_arr_elem "${arr_scl}"  "${idx}")"
    fi

    if [[ -n "${arr_opt}" ]] ; then
       opt_var="$(get_arr_elem "${arr_opt}" "${idx}")"
    fi

    #  Debug inputs
    if [[ "${debug}" == "true" ]]; then
        if [[ "${mode}" == "signal" ]]; then
            debug_var \
                "infile=${infile_1}" \
                "outfile=${outfile}" \
                "scl_fct=${scl_fct}" \
                "usr_frg=${opt_var}"
        elif [[ "${mode}" == "ratio" ]]; then
            debug_var \
                "fil_ip=${infile_1}" \
                "fil_in=${infile_2}" \
                "outfile=${outfile}" \
                "scl_fct=${scl_fct}" \
                "dep_min=${opt_var}"
        else
            debug_var \
                "infile=${infile_1}" \
                "outfile=${outfile}"
        fi
    fi

    #  Derive sample ('samp') and descriptor ('dsc') via process_io
    if [[ "${mode}" == "ratio" ]]; then
        IFS=',' read -r samp dsc < <(
            process_io \
                 -m "${mode}" \
                -i1 "${infile_1}" \
                -i2 "${infile_2}" \
                 -o "${outfile}" \
                -sf "${scl_fct:-NA}" \
                -ov "${opt_var:-NA}"
        ) || return 1
    elif [[ "${mode}" == "signal" ]]; then
        IFS=',' read -r samp dsc < <(
            process_io \
                 -m "${mode}" \
                -i1 "${infile_1}" \
                 -o "${outfile}" \
                -sf "${scl_fct:-NA}" \
                -ov "${opt_var:-NA}"
        ) || return 1
    else
        IFS=',' read -r samp dsc < <(
            process_io \
                 -m "${mode}" \
                -i1 "${infile_1}" \
                 -o "${outfile}"
        ) || return 1
    fi

    #  Debug 'samp' and 'dsc' output by process_io
    if [[ "${debug}" == "true" ]]; then
        debug_var "samp=${samp}" "dsc=${dsc}"
    fi

    #  If using SLURM, request initial and descriptor log paths
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
    echo "${infile_1},${infile_2:-},${outfile},${samp},${dsc},${err_ini:-},${out_ini:-},${err_dsc:-},${out_dsc:-}"
}


#  Remove initial SLURM logs only if not in dry-run mode; used by all
#+ run_task_* helper functions
function task_epilogue() {
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        local err_ini="${1}"
        local out_ini="${2}"

        if [[ -n "${err_ini}" && -f "${err_ini}" ]]; then
            rm -f "${err_ini}"
        fi

        if [[ -n "${out_ini}" && -f "${out_ini}" ]]; then
            rm -f "${out_ini}"
        fi
    fi
}


#  Helper function to run SLURM or non-SLURM signal computation jobs
function run_task_signal() {
    local idx="${1}"
    local infile unused outfile samp dsc err_ini out_ini err_dsc out_dsc

    # shellcheck disable=SC2034
    IFS=',' read -r \
        infile unused outfile samp dsc \
        err_ini out_ini err_dsc out_dsc \
    < <(
        task_prelude "signal" "${idx}" \
            "arr_infile" "" "arr_outfile" "arr_scl_fct" "arr_usr_frg"
    ) || return 1

    run_comp_signal \
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

    task_epilogue "${err_ini}" "${out_ini}"
}


#  Helper function to run SLURM or non-SLURM ratio computation jobs
function run_task_ratio() {
    local idx="${1}"
    local fil_ip fil_in outfile samp dsc err_ini out_ini err_dsc out_dsc

    IFS=',' read -r \
        fil_ip fil_in outfile samp dsc \
        err_ini out_ini err_dsc out_dsc \
    < <(
        task_prelude "ratio" "${idx}" \
            "arr_fil_ip" "arr_fil_in" "arr_outfile" "arr_scl_fct" "arr_dep_min"
    ) || return 1

    run_comp_ratio \
        "${debug}" \
        "${fil_ip}" \
        "${fil_in}" \
        "${outfile}" \
        "$(get_arr_elem arr_scl_fct "${idx}")" \
        "$(get_arr_elem arr_dep_min "${idx}")" \
        "${rnd}" \
        "${track}" \
        "${log2}" \
        "${err_out}" \
        "${nam_job}" \
        "${dsc}"

    task_epilogue "${err_ini}" "${out_ini}"
}


#  Helper function to run SLURM or non-SLURM fragment coordinate extraction
#+ jobs
function run_task_coord() {
    local idx="${1}"
    local infile unused outfile samp dsc err_ini out_ini err_dsc out_dsc

    # shellcheck disable=SC2034
    IFS=',' read -r \
        infile unused outfile samp dsc \
        err_ini out_ini err_dsc out_dsc \
    < <(
        task_prelude "coord" "${idx}" \
            "arr_infile" "" "arr_outfile" "" ""
    ) || return 1

    #  Use stub parameters per your original coord behavior
    run_comp_signal \
        "${debug}" \
        1 \
        "${infile}" \
        "${outfile}" \
        1 \
        "" \
        "NA" \
        "NA" \
        1 \
        "${err_out}" \
        "${nam_job}" \
        "${dsc}"

    task_epilogue "${err_ini}" "${out_ini}"
}

######################
## Work in progress ##
######################


#  Extract values from a specified field from a delimited text table (e.g.,
#+ CSV or TSV), optionally skipping a header line; returns values as a comma-
#+ separated string (e.g., "val1,val2,val3")
function extract_field_str() {
    tbl="${1}"         # Input table file path
    fld="${2}"         # 1-based field index to extract
    hdr="${3:-false}"  # Skip header line if 'true'; default: 'false'

    awk \
        -v fld="${fld}" \
        -v hdr="${hdr}" \
        '(hdr == "true" && NR == 1) {
            next
        } {
            val = val ? val "," $fld : $fld
        } END {
            print val
        }' \
            "${tbl}"
}


#  Set paths, values, arguments, etc. for interactive mode
function set_interactive() {
    #  Set base repository paths
    dir_bas="${HOME}/repos"
    dir_rep="${dir_bas}/protocol_chipseq_signal_norm"

    #  Set alignment parameters
    aligner="bowtie2"
    a_type="global"
    req_flg=true
    flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
    mapq=1
    str_det="${aligner}_${a_type}_flag-${flg}_mapq-${mapq}"

    #  Define data directories and data table
    dir_dat="${dir_rep}/data"
    dir_pro="${dir_dat}/processed"

    #  Set hardcoded argument assignments (and related assignments), including
    #+ input and output directories
    env_nam="env_protocol"
    dir_scr="${dir_rep}/scripts"
    threads=4  # 1  # 8
    sf=""  # "spike"  # "alpha"
    mode="ratio"  # "coord"  # "signal"
    method="log2"  # "norm"  # "unadj"

    if [[ "${sf}" =~ ^(alpha|spike)$ ]]; then
        tbl="${dir_pro}/compute_signal/${str_det}/tables/${sf}_test.tsv"
    elif [[ "${mode}" == "ratio" && "${method}" == "log2" && -z "${sf}" ]]; then
        use_dm=true
        tbl="${dir_pro}/compute_signal/${str_det}/tables/alpha_test.tsv"
    fi

    if [[ "${mode}" == "coord" ]]; then
        unset method
        dir_out="${dir_pro}/compute_${mode}/${str_det}"
    else
        if [[ "${sf}" =~ ^(alpha|spike)$ ]]; then
            dir_out="${dir_pro}/compute_${mode}/${str_det}/${method}/${sf}"
        else
            dir_out="${dir_pro}/compute_${mode}/${str_det}/${method}"
        fi
    fi

    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        dir_in="${dir_pro}/align_fastqs/${str_det}/sc"
        ser_infile="$(
            if [[ "${sf}" != "spike" ]]; then
                bash "${dir_scr}/find_files.sh" \
                    --dir_fnd "${dir_in}" \
                    --pattern "*.bam" \
                    --exclude "*Brn1*"
            else
                bash "${dir_scr}/find_files.sh" \
                    --dir_fnd "${dir_in}" \
                    --pattern "*.bam" \
                    --include "IP*" \
                    --exclude "*Brn1*"
            fi
        )"
        ser_outfile="$(
            if [[ "${mode}" == "coord" ]]; then
                sed \
                    -e "s|.bam|.bed.gz|g" \
                    -e "s|${dir_in}|${dir_out}|g" \
                        < <(echo "${ser_infile}")
            else
                if [[ "${sf}" != "spike" ]]; then
                    sed \
                        -e "s|.bam|.bdg.gz|g" \
                        -e "s|${dir_in}|${dir_out}|g" \
                            < <(echo "${ser_infile}")
                else
                    sed \
                        -e "s|/IP|/${sf}_IP|g" \
                        -e "s|.bam|.bdg.gz|g" \
                        -e "s|${dir_in}|${dir_out}|g" \
                            < <(echo "${ser_infile}")
                fi
            fi
        )"
    elif [[ "${mode}" == "ratio" ]]; then
        dir_in="${dir_pro}/compute_signal/${str_det}/norm"
        ser_fil_ip="$(
            bash "${dir_scr}/find_files.sh" \
                --dir_fnd "${dir_in}" \
                --pattern "*.bdg.gz" \
                --include "IP*" \
                --exclude "*Brn1*"
        )"
        ser_fil_in="$(sed 's|/IP|/in|g' < <(echo "${ser_fil_ip}"))"
        ser_outfile="$(
            if [[ "${sf}" =~ ^(alpha|spike)$ ]]; then
                if [[ "${mode}" == "signal" ]]; then
                    sed \
                        -e "s|/IP|/${sf}_IP|g" \
                        -e "s|${dir_in}|${dir_out}|g" \
                            < <(echo "${ser_fil_ip}")
                else
                    sed \
                        -e "s|/IP|/${sf}|g" \
                        -e "s|${dir_in}|${dir_out}|g" \
                            < <(echo "${ser_fil_ip}")
                fi
            else
                sed \
                    -e "s|/IP|/${mode}_${method}|g" \
                    -e "s|${dir_in}|${dir_out}|g" \
                        < <(echo "${ser_fil_ip}")
            fi
        )"
    fi

    if [[ "${mode}" != "ratio" ]]; then
        track=false
    else
        track=true
        ser_scl_fct="$(
            if [[ "${sf}" == "alpha" ]]; then
                extract_field_str "${tbl}" 3 false
            elif [[ "${sf}" == "spike" ]]; then
                extract_field_str "${tbl}" 5 false
            else
                echo "NA,NA"
            fi
        )"
        ser_dep_min="$(
            if [[ "${sf}" == "alpha" || "${use_dm}" == "true" ]]; then
                extract_field_str "${tbl}" 24 false
            elif [[ "${sf}" == "spike" ]]; then
                extract_field_str "${tbl}" 21 false
            else
                echo "NA,NA"
            fi
        )"
    fi

    if [[ "${mode}" == "signal" ]]; then
        siz_bin=30
        ser_scl_fct="$(
            if [[ "${sf}" == "spike" ]]; then
                extract_field_str "${tbl}" 5 false
            else
                echo "NA,NA,NA,NA"
            fi
        )"
        ser_usr_frg="$(
            if [[ "${sf}" == "spike" ]]; then
                echo "NA,NA"
            else
                echo "NA,NA,NA,NA"
            fi
        )"
    fi

    if [[ "${mode}" != "coord" ]]; then
        rnd=24
    fi

    err_out="${dir_out}/logs"
    if [[ -n "${method}" ]]; then
        nam_job="compute_${mode}_${method}"
    else
        nam_job="compute_${mode}"
    fi

    #  If in "debug mode", print variable assignments 
    if [[ "${debug}" == "true" ]]; then
        echo "#######################################################"
        echo "## Debug variable assignments for 'interactive mode' ##"
        echo "#######################################################"
        echo ""
        echo "#  Hardcoded arguments"
        echo "interactive=${interactive}"
        echo "debug=${debug}"
        echo "dry_run=${dry_run}"
        echo "unit=${unit}"
        echo ""
        echo "#  Set base repository paths"
        echo "dir_bas=${dir_bas}"
        echo "dir_rep=${dir_rep}"
        echo ""
        echo "sf=${sf:-UNSET}"
        echo "tbl=${tbl:-UNSET}"
        echo ""
        echo "dir_in=${dir_in}"
        echo "dir_out=${dir_out}"
        echo ""
        echo "#  Set alignment parameters"
        echo "aligner=${aligner}"
        echo "a_type=${a_type}"
        echo "req_flg=${req_flg}"
        echo "flg=${flg}"
        echo "mapq=${mapq}"
        echo "str_det=${str_det}"
        echo ""
        echo "#  Define data directories"
        echo "dir_dat=${dir_dat}"
        echo "dir_pro=${dir_pro}"
        echo ""
        echo "#  Set hardcoded argument assignments, etc."
        echo "env_nam=${env_nam}"
        echo "dir_scr=${dir_scr}"
        echo "threads=${threads}"
        echo "mode=${mode}"
        echo "method=${method:-UNSET}"
        echo "dir_in=${dir_in}"
        echo "dir_out=${dir_out}"
        echo "ser_infile=${ser_infile:-UNSET}"
        echo "ser_fil_ip=${ser_fil_ip:-UNSET}"
        echo "ser_fil_in=${ser_fil_in:-UNSET}"
        echo "ser_outfile=${ser_outfile}"
        echo "track=${track}"
        echo "siz_bin=${siz_bin:-UNSET}"
        echo "ser_scl_fct=${ser_scl_fct:-UNSET}"
        echo "ser_usr_frg=${ser_usr_frg:-UNSET}"
        echo "ser_dep_min=${ser_dep_min:-UNSET}"
        echo "rnd=${rnd:-UNSET}"
        echo "err_out=${err_out}"
        echo "nam_job=${nam_job}"
        echo ""
        echo ""
    fi
}


#  Initialize hardcoded argument variable (if '--method log2', then later
#+ overwritten with 'true')
log2=false

#  Initialize argument variables, assigning default values where applicable
env_nam="env_protocol"
dir_scr=""
threads=4
mode="signal"
method="norm"
ser_infile=""
ser_fil_ip=""
ser_fil_in=""
ser_outfile=""
track=false
siz_bin=10
ser_scl_fct=""
ser_usr_frg=""
ser_dep_min=""
rnd=24
err_out=""
nam_job="compute_${mode}_${method}"

#  Define the help message
show_help=$(cat << EOM
Usage:
  submit_compute_signal.sh
    [--help] --env_nam <str> --dir_scr <str> --threads <int> --mode <str> --method <str> --ser_infile <str> --ser_fil_ip <str> --ser_fil_in <str>
    --ser_outfile <str> [--track] --siz_bin <int> [--ser_scl_fct <str>] [--ser_usr_frg <str>] [--ser_dep_min <str>] --rnd <int> --err_out <str> --nam_job <str>

Arguments:
   -h, --help         <flg>  Print this help message and exit
  -en, --env_nam      <str>  Mamba environment to activate (default: "${env_nam}")
  -ds, --dir_scr      <str>  Directory containing scripts and functions
   -t, --threads      <int>  Number of threads to use per job (default: "${threads}")
   -m, --mode         <str>  Type of computation: 'signal', 'ratio', or 'coord' (default: "${mode}")
  -me, --method       <str>  Computation subtype (ignored if '--mode coord'; default: "${method}"):
                               - 'unadj', 'frag', 'norm' ('--mode signal')
                               - 'unadj', 'log2' ('--mode ratio')
   -i, --ser_infile   <str>  Comma-separated list of BAM files ('--mode signal', '--mode coord') <element: str>
  -ip, --ser_fil_ip   <str>  Comma-separated list of IP bedGraph files ('--mode ratio') <element: str>
  -in, --ser_fil_in   <str>  Comma-separated list of input bedGraph files ('--mode ratio') <element: str>
   -o, --ser_outfile  <str>  Comma-separated list of output file stems <element: str>
  -tr, --track        <flg>  Output a companion bedGraph without '-inf' or 'nan' rows ('--mode ratio')
  -sb, --siz_bin      <int>  Bin size in base pairs for signal computation ('--mode signal'; default: "${siz_bin}")
  -sf, --ser_scl_fct  <str>  Comma-separated list of scaling factors ('--mode signal', '--mode ratio') <element: flt>
  -uf, --ser_usr_frg  <str>  Comma-separated list of fragment lengths to override auto-detection ('--mode signal') <element: int>
  -dm, --ser_dep_min  <str>  Comma-separated list of minimum input depth values ('--mode ratio') <element: flt>
   -r, --rnd          <int>  Decimal precision for output signal or ratio values (default: "${rnd}")
  -eo, --err_out      <str>  Directory for stderr and stdout logs
  -nj, --nam_job      <str>  Job name prefix (default: "${nam_job}")

Notes:
  - To run in "interactive mode", set hardcoded variable 'interactive=true'  [interactive=${interactive:-UNSET}]
  - To run in "debug mode", set hardcoded variable 'debug=true'              [debug=${debug:-UNSET}]
  - To run in "dry-run mode", set hardcoded variable 'dry_run=true'          [dry_run=${dry_run:-UNSET}]
  - To run in "unit-test mode", set hardcoded variable 'unit=true'           [unit=${unit:-UNSET}]
EOM
)

#  Parse arguments, assigning them to variables
#+ 
#+ Most of the argument inputs are not checked, as this is performed by
#+ 'execute_*.sh' and, to a certain extent, corresponding Python scripts
if [[ "${interactive}" == "false" ]]; then
    if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        exit 0
    fi
fi

if [[ "${interactive}" == "true" ]]; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -en|--env_nam)     env_nam="${2}";     shift 2 ;;
            -ds|--dir_scr)     dir_scr="${2}";     shift 2 ;;
             -t|--threads)     threads="${2}";     shift 2 ;;
             -m|--mode)        mode="${2}";        shift 2 ;;
            -me|--method)      method="${2}";      shift 2 ;;
             -i|--ser_infile)  ser_infile="${2}";  shift 2 ;;
            -ip|--ser_fil_ip)  ser_fil_ip="${2}" ; shift 2 ;;
            -in|--ser_fil_in)  ser_fil_in="${2}" ; shift 2 ;;
             -o|--ser_outfile) ser_outfile="${2}"; shift 2 ;;
            -tr|--track)       track=true;         shift 1 ;;
            -sb|--siz_bin)     siz_bin="${2}";     shift 2 ;;
            -sf|--ser_scl_fct) ser_scl_fct="${2}"; shift 2 ;;
            -uf|--ser_usr_frg) ser_usr_frg="${2}"; shift 2 ;;
            -dm|--ser_dep_min) ser_dep_min="${2}"; shift 2 ;;
             -r|--rnd)         rnd="${2}";         shift 2 ;;
            -eo|--err_out)     err_out="${2}";     shift 2 ;;
            -nj|--nam_job)     nam_job="${2}";     shift 2 ;;
            *)
                echo "## Unknown parameter passed: '${1}' ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                exit 1
                ;;
        esac
    done
fi

if [[ "${unit}" == "true" ]]; then exit 0; fi

#  Validate 'mode' and 'method'
case "${mode}" in
    signal)
        case "${method}" in
            unadj|frag|norm) : ;;
            *)
                echo \
                    "Error: Invalid value for '--method': '${method}'." \
                    "Expected 'unadj', 'frag', or 'norm' for '--mode ${mode}'."
                if [[ "${interactive}" == "false" ]]; then exit 1; fi
                ;;
        esac
        ;;
    ratio)
        case "${method}" in
            unadj) : ;;
            log2)  log2=true ;;
            *)
                echo \
                    "Error: Invalid value for '--method': '${method}'." \
                    "Expected 'unadj' or 'log2' for '--mode ${mode}'."
                if [[ "${interactive}" == "false" ]]; then exit 1; fi
                ;;
        esac
        ;;
    coord) unset method ;;
    *)
        echo \
            "Error: Invalid value for '--mode': '${mode}'. Expected" \
            "'signal', 'ratio', or 'coord'."
        if [[ "${interactive}" == "false" ]]; then exit 1; fi
        ;;
esac

#  Assign and validate variables for scripts and functions
if [[ ! -d "${dir_scr}" ]]; then
    echo "Error: Script directory not found: '${dir_scr}'." >&2
    if [[ "${interactive}" == "false" ]]; then exit 1; fi
else
    scr_sig="${dir_scr}/compute_signal.py"
    scr_rat="${dir_scr}/compute_signal_ratio.py"
    scr_sub="${dir_scr}/functions/submit.sh"

    for fil in "${scr_sig}" "${scr_rat}" "${scr_sub}"; do
        if [[ ! -f "${fil}" ]]; then
            echo "Error: Script not found: '${fil}'." >&2
            if [[ "${interactive}" == "false" ]]; then exit 1; fi
        fi
    done
    unset fil

    # shellcheck disable=SC1090
    source "${scr_sub}" || {
        echo "Error: Failed to source '${scr_sub}'." >&2
        if [[ "${interactive}" == "false" ]]; then exit 1; fi
    }
fi

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
            "ser_infile=${ser_infile}"
    elif [[ "${mode}" == "ratio" ]]; then
        debug_var \
            "ser_fil_ip=${ser_fil_ip}" \
            "ser_fil_in=${ser_fil_in}"
    fi

    debug_var \
        "ser_outfile=${ser_outfile}"

    if [[ "${mode}" == "ratio" ]]; then
        debug_var \
            "track=${track}" \
            "log2=${log2}" \
            "ser_dep_min=${ser_dep_min}"
    fi

    if [[ "${mode}" == "signal" ]]; then
        debug_var \
            "siz_bin=${siz_bin}" \
            "ser_usr_frg=${ser_usr_frg}"
    fi

    if [[ "${mode}" != "coord" ]]; then
        debug_var \
            "ser_scl_fct=${ser_scl_fct}" \
            "rnd=${rnd}"
    fi

    debug_var \
        "err_out=${err_out}" \
        "nam_job=${nam_job}"
fi

#  Activate environment
activate_env "${env_nam}" || exit 1

#  Perform mode-dependent array reconstruction from serialized strings
IFS=','
if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
    read -r -a arr_infile  <<< "${ser_infile}"
    
    if [[ "${mode}" == "signal" ]]; then
        read -r -a arr_usr_frg <<< "${ser_usr_frg}"
    fi
fi

if [[ "${mode}" == "ratio" ]]; then
    read -r -a arr_fil_ip  <<< "${ser_fil_ip}"
    read -r -a arr_fil_in  <<< "${ser_fil_in}"
    read -r -a arr_dep_min <<< "${ser_dep_min}"
fi

if [[ "${mode}" != "coord" ]]; then
    read -r -a arr_scl_fct <<< "${ser_scl_fct}"
fi

read -r -a arr_outfile <<< "${ser_outfile}"
IFS=$' \t\n'

#  Debug output to check number of array elements and array element values
if [[ "${debug}" == "true" ]]; then
    if [[ "${mode}" =~ ^(signal|coord)$ ]]; then
        echo "\${#arr_infile[@]}=${#arr_infile[@]}" >&2       && echo "" >&2
        echo "arr_infile=( ${arr_infile[*]} )" >&2            && echo "" >&2
        
        if [[ "${mode}" == "signal" ]]; then
            echo "\${#arr_usr_frg[@]}=${#arr_usr_frg[@]}" >&2 && echo "" >&2
            echo "arr_usr_frg=( ${arr_usr_frg[*]} )" >&2      && echo "" >&2
        fi
    fi

    if [[ "${mode}" == "ratio" ]]; then
        echo "\${#arr_fil_ip[@]}=${#arr_fil_ip[@]}" >&2       && echo "" >&2
        echo "arr_fil_ip=( ${arr_fil_ip[*]} )" >&2            && echo "" >&2
        echo "\${#arr_fil_in[@]}=${#arr_fil_in[@]}" >&2       && echo "" >&2
        echo "arr_fil_in=( ${arr_fil_in[*]} )" >&2            && echo "" >&2
        echo "\${#arr_dep_min[@]}=${#arr_dep_min[@]}" >&2     && echo "" >&2
        echo "arr_dep_min=( ${arr_dep_min[*]} )" >&2          && echo "" >&2
    fi

    if [[ "${mode}" != "coord" ]]; then
        echo "\${#arr_scl_fct[@]}=${#arr_scl_fct[@]}" >&2     && echo "" >&2
        echo "arr_scl_fct=( ${arr_scl_fct[*]} )" >&2          && echo "" >&2
    fi

    echo "\${#arr_outfile[@]}=${#arr_outfile[@]}" >&2         && echo "" >&2
    echo "arr_outfile=( ${arr_outfile[*]} )" >&2              && echo "" >&2
fi

#  Determine and run mode: SLURM or GNU Parallel/serial
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    #  Job submission type: SLURM
    id_job=${SLURM_ARRAY_JOB_ID}
    id_tsk=${SLURM_ARRAY_TASK_ID}

    if [[ "${id_tsk}" -lt 1 ]]; then
        echo "Error: SLURM task ID is invalid: '${id_tsk}'." >&2
        exit 1
    else
        idx=$(( id_tsk - 1 ))
    fi

    if [[ "${mode}" == "signal" ]]; then
        run_task_signal "${idx}" || exit 1
    elif [[ "${mode}" == "ratio" ]]; then
        run_task_ratio "${idx}" || exit 1
    else
        run_task_coord "${idx}" || exit 1
    fi
else
    #  Job submission type: GNU Parallel or serial
    if [[ "${mode}" == "signal" ]]; then
        for idx in "${!arr_infile[@]}"; do
            run_task_signal "${idx}" || exit 1
        done
    elif [[ "${mode}" == "ratio" ]]; then
        for idx in "${!arr_fil_ip[@]}"; do
            run_task_ratio "${idx}" || exit 1
        done
    else
        for idx in "${!arr_infile[@]}"; do
            run_task_coord "${idx}" || exit 1
        done
    fi
fi
