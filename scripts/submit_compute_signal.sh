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
    local mode=""
    local infile_1=""
    local infile_2=""
    local outfile=""
    local scl_fct=""
    local opt_var=""
    local help samp dsc

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

    echo "${samp};${dsc}"
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
    local log_out log_err cmd

    #  Generate optional arguments array dynamically
    IFS="," read -r -a optional <<< "$(
        set_args_opt "signal" "${scl_fct}" "${usr_frg}" "${rnd}"
    )"

    #  Debug array of optional arguments
    if [[ "${debug}" == "true" ]]; then
        echo "set_args_opt(): optional=( ${optional[*]} )"
        echo ""
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

    #  Debug call to 'compute_signal.py'
    if [[ "${debug}" == "true" || "${dry_run}" == "true" ]]; then
        printf '%q ' "${cmd[@]}"
        echo ">> ${log_out} 2>> ${log_err}"
        echo ""
        echo ""
    fi

    #  Execute 'compute_signal.py', writing stdout and stderr to log files
    if [[ "${dry_run}" == "false" ]]; then
        "${cmd[@]}" >> "${log_out}" 2>> "${log_err}"
    fi
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
    
    #  Generate optional arguments array dynamically
    IFS="," read -r -a optional <<< "$(
        set_args_opt "ratio" \
            "${scl_fct}" "${dep_min}" "${rnd}" "${track}" "${log2}"
    )"

    #  Debug array of optional arguments
    if [[ "${debug}" == "true" ]]; then
        echo "set_args_opt(): optional=( ${optional[*]} )"
        echo ""
    fi

    #  Define paths for log output files
    log_out="${err_out}/${nam_job}.${dsc}.stdout.txt"
    log_err="${err_out}/${nam_job}.${dsc}.stderr.txt"

    #  Debug call to 'compute_signal_ratio.py'
    if [[ "${debug}" == "true" || ${dry_run} == "true" ]]; then
        echo "python ${scr_rat} \\"
        echo "    --verbose \\"
        echo "    --fil_ip ${fil_ip} \\"
        echo "    --fil_in ${fil_in} \\"
        echo "    --fil_out ${outfile} \\"
        echo "    ${optional[*]} \\"
        echo "        >> ${log_out} \\"
        echo "       2>> ${log_err}"
        echo ""
    fi

    #  Execute 'compute_signal_ratio.py', writing stdout and stderr to log
    #+ files
    if [[ ${dry_run} == "false" ]]; then
        python "${scr_rat}" \
            --verbose \
            --fil_ip "${fil_ip}" \
            --fil_in "${fil_in}" \
            --fil_out "${outfile}" \
            "${optional[@]}" \
                 >> "${log_out}" \
                2>> "${log_err}"
    fi
}


#  Helper function to run SLURM or non-SLURM signal computation jobs
function run_task_signal() {
    local idx="${1}"

    infile="${arr_infile[idx]}"
    outfile="${arr_outfile[idx]}"
    scl_fct="${arr_scl_fct[idx]}"
    usr_frg="${arr_usr_frg[idx]}"

    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "infile=${infile}" \
            "outfile=${outfile}" \
            "scl_fct=${scl_fct}" \
            "usr_frg=${usr_frg}"
    fi

    IFS=';' read -r samp dsc < <(
        process_io \
             -m "${mode}"    -i1 "${infile}"  \
             -o "${outfile}" -sf "${scl_fct}" \
            -ov "${usr_frg}"
    ) || return 1

    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "samp=${samp}" \
            "dsc=${dsc}"
    fi

    #  Set up log files if using SLURM
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        IFS=';' read -r err_ini out_ini err_dsc out_dsc < <(
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

    #  Run 'compute_signal.py'
    run_comp_signal \
        "${debug}" \
        "${threads}" \
        "${infile}" \
        "${outfile}" \
        "${siz_bin}" \
        "${method}" \
        "${scl_fct}" \
        "${usr_frg}" \
        "${rnd}" \
        "${err_out}" \
        "${nam_job}" \
        "${dsc}"

    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" && "${dry_run}" == "false" ]]; then
        rm "${err_ini}" "${out_ini}"
    fi
}


#  Helper function to run SLURM or non-SLURM ratio computation jobs
function run_task_ratio() {
    local idx="${1}"

    fil_ip="${arr_fil_ip[idx]}"
    fil_in="${arr_fil_in[idx]}"
    outfile="${arr_outfile[idx]}"
    scl_fct="${arr_scl_fct[idx]}"
    dep_min="${arr_dep_min[idx]}"

    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "fil_ip=${fil_ip}" \
            "fil_in=${fil_in}" \
            "outfile=${outfile}" \
            "scl_fct=${scl_fct}" \
            "dep_min=${dep_min}"
    fi

    IFS=';' read -r samp dsc < <(
        process_io \
             -m "${mode}"    -i1 "${fil_ip}"  \
            -i2 "${fil_in}"   -o "${outfile}" \
            -sf "${scl_fct}" -ov "${dep_min}"
    ) || return 1

    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "samp=${samp}" \
            "dsc=${dsc}"
    fi

    #  Set up log files if using SLURM
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        IFS=';' read -r err_ini out_ini err_dsc out_dsc < <(
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

    #  Run 'compute_signal_ratio.py'
    run_comp_ratio \
        "${debug}" \
        "${fil_ip}" \
        "${fil_in}" \
        "${outfile}" \
        "${scl_fct}" \
        "${dep_min}" \
        "${rnd}" \
        "${track}" \
        "${log2}" \
        "${err_out}" \
        "${nam_job}" \
        "${dsc}"

    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" && "${dry_run}" == "false" ]]; then
        rm "${err_ini}" "${out_ini}"
    fi
}


#  Helper function to run SLURM or non-SLURM fragment coordinate extraction
#+ jobs
function run_task_coord() {
    local idx="${1}"

    infile="${arr_infile[idx]}"
    outfile="${arr_outfile[idx]}"

    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "infile=${infile}" \
            "outfile=${outfile}"
    fi

    IFS=';' read -r samp dsc < <(
        process_io \
            -m "${mode}" -i1 "${infile}" -o "${outfile}"
    ) || return 1

    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "samp=${samp}" \
            "dsc=${dsc}"
    fi

    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        IFS=';' read -r err_ini out_ini err_dsc out_dsc < <(
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

    #  Call 'run_comp_signal' with stub arguments
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

    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" && "${dry_run}" == "false" ]]; then
        rm "${err_ini}" "${out_ini}"
    fi
}


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
  -si, --ser_infile   <str>  Comma-separated list of BAM files ('--mode signal', '--mode coord') <element: str>
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
        echo "\${#arr_infile[@]}=${#arr_infile[@]}"       && echo ""
        echo "arr_infile=( ${arr_infile[*]} )"            && echo ""
        
        if [[ "${mode}" == "signal" ]]; then
            echo "\${#arr_usr_frg[@]}=${#arr_usr_frg[@]}" && echo ""
            echo "arr_usr_frg=( ${arr_usr_frg[*]} )"      && echo ""
        fi
    fi

    if [[ "${mode}" == "ratio" ]]; then
        echo "\${#arr_fil_ip[@]}=${#arr_fil_ip[@]}"       && echo ""
        echo "arr_fil_ip=( ${arr_fil_ip[*]} )"            && echo ""
        echo "\${#arr_fil_in[@]}=${#arr_fil_in[@]}"       && echo ""
        echo "arr_fil_in=( ${arr_fil_in[*]} )"            && echo ""
        echo "\${#arr_dep_min[@]}=${#arr_dep_min[@]}"     && echo ""
        echo "arr_dep_min=( ${arr_dep_min[*]} )"          && echo ""
    fi

    if [[ "${mode}" != "coord" ]]; then
        echo "\${#arr_scl_fct[@]}=${#arr_scl_fct[@]}"     && echo ""
        echo "arr_scl_fct=( ${arr_scl_fct[*]} )"          && echo ""
    fi

    echo "\${#arr_outfile[@]}=${#arr_outfile[@]}"         && echo ""
    echo "arr_outfile=( ${arr_outfile[*]} )"              && echo ""
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
