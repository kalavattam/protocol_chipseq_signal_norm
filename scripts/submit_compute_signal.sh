#!/bin/bash

#  submit_compute_signal.sh
#  KA


#  If true, run script in debug mode
debug=true


#  Define functions
#  Function to print variable assignment(s)
function debug_var() { printf "%s\n\n" "$@"; }


#  Function to activate user-supplied Conda/Mamba environment
function activate_env() {
    local env_nam="${1}"  # Name of environment to activate
    if [[ "${CONDA_DEFAULT_ENV}" != "${env_nam}" ]]; then
        # shellcheck disable=SC1091
        source "$(conda info --base)/etc/profile.d/conda.sh"
        if ! \
            conda activate "${env_nam}"
        then
            echo "Error: Failed to activate environment: '${env_nam}'." >&2
            return 1
        fi
    fi
}


#  Function to validate that a variable is not empty or unset
function validate_var() {
    local var_nam="${1}"   # Variable name (for error messages)
    local var_val="${2-}"  # Value to check for emptiness/unset state
    local idx="${3:-0}"    # Optional index (for arrays); defaults to '0'

    if [[ -z "${var_val}" ]]; then
        if [[ "${idx}" -ne 0 ]]; then
            echo \
                "Error: '${var_nam}' is empty or unset for array index" \
                "'${idx}'." >&2
        else
            echo "Error: '${var_nam}' is empty or unset." >&2
        fi
        return 1
    fi
}


#  Function to debug, validate, and parse 'infile', 'outfile', 'scl_fct', and
#+ 'usr_frg', assigning 'samp' (sample name) from 'infile' and 'dsc' (output
#+ sample description) from 'outfile'
function process_input() {
    local infile="${1}"   # Input BAM file
    local outfile="${2}"  # Output file stem
    local scl_fct="${3}"  # Scaling factor/coefficient
    local usr_frg="${4}"  # Fragment length
    local samp  # Sample name derived from input BAM file

    #  Validate input arguments
    validate_var "infile"  "${infile}"  || return 1
    validate_var "outfile" "${outfile}" || return 1
    validate_var "scl_fct" "${scl_fct}" || return 1
    validate_var "usr_frg" "${usr_frg}" || return 1

    #  Extract the sample name by stripping the ".bam" extension
    samp="$(basename "${infile}" ".bam")"
    dsc=$(basename "${outfile}")

    #  Return value
    echo "${samp};${dsc}"
}


#  Function to set up SLURM log file paths and create symlinked log files
function set_logs_slurm() {
    local id_job="${1}"   # SLURM job ID
    local id_tsk="${2}"   # SLURM task ID within the job array
    local samp="${3}"     # Sample name associated with log files
    local err_out="${4}"  # Directory for stderr and stdout log files
    local nam_job="${5}"  # Name of the job (used in log file naming)
    local err_ini  # SLURM initial stderr log 
    local out_ini  # SLURM initial stdout log 
    local err_dsc  # Symlinked stderr log with dscriptive name
    local out_dsc  # Symlinked stdout log with dscriptive name

    #  Set log paths
    err_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stderr.txt"
    out_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stdout.txt"
    err_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stderr.txt"
    out_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stdout.txt"

    #  Create symlinked log files
    ln -f "${err_ini}" "${err_dsc}"
    ln -f "${out_ini}" "${out_dsc}"

    #  Return values
    echo "${err_ini};${out_ini};${err_dsc};${out_dsc}"
}


#  Function to set optional arguments for the call to 'scr_sig'
function set_args_opt() {
    unset optional
    typeset -a optional
    # typeset -g -a optional  # No go for BSD Bash v3.2

    if [[ -n "${typ_cvg}" ]]; then
        optional+=( --typ_cvg "${typ_cvg}" )
    fi

    if [[ -n "${scl_fct}" && "${scl_fct}" != "#N/A" ]]; then
        optional+=( --scl_fct "${scl_fct}" )
    fi

    if [[ "${usr_frg}" != "#N/A" ]]; then
        optional+=( --usr_frg "${usr_frg}" )
    fi

    if [[ -n "${rnd}" && "${rnd}" != "#N/A" ]]; then
        optional+=( --rnd "${rnd}" )
    fi
}


#  Function to run 'compute_signal.py'
function run_comp_sig() {
    local threads="${1}"  # Number of threads to use
    local infile="${2}"   # Input BAM file
    local outfile="${3}"  # Output file stem
    local siz_bin="${4}"  # Bin size in base pairs
    local str_opt="${5}"  # Serialized optional arguments (comma-separated)
    local err_out="${6}"  # Directory for stdout and stderr
    local nam_job="${7}"  # Job name
    local dsc="${8}"      # Descriptor for log file naming

    #  Reconstruct optional arguments array
    unset optional
    IFS="," read -r -a optional <<< "${str_opt}"
    unset IFS

    #  Execute the 'compute_signal.py', writing stdout and stderr to log files
    python "${scr_sig}" \
        --verbose \
        --threads "${threads}" \
        --infile "${infile}" \
        --outfile "${outfile}" \
        --siz_bin "${siz_bin}" \
        "${optional[@]}" \
             > "${err_out}/${nam_job}.${dsc}.stdout.txt" \
            2> "${err_out}/${nam_job}.${dsc}.stderr.txt"
}


#  Define the help message
show_help=$(cat << EOM
\${1}=env_nam       # str: Name of Conda/Mamba environment to activate
\${2}=scr_sig       # str: Path to signal script
\${3}=threads       # int: Number of threads to use
\${4}=str_infile    # str: Comma-separated string of infiles (str)
\${5}=str_outfile   # str: Comma-separated string of outfile stems (str)
\${6}=siz_bin       # int: Bin size in base pairs
\${7}=typ_cvg       # str: Type of signal to compute
\${8}=str_scl_fct   # str: Comma-separated string of scaling factors (flt)
\${9}=str_usr_frg   # str: Comma-separated string of fragment lengths (int)
\${10}=rnd          # int: No. decimal places for rounding signal values
\${11}=err_out      # str: Directory for stdout and stderr files
\${12}=nam_job      # str: Name of job
EOM
)

#  Display help message if a help option or no arguments are given
if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    cat << EOM
$(basename "${0}") requires 12 positional arguments:
${show_help}
EOM
    exit 0
fi

#  Check for exactly 12 positional arguments
if [[ $# -ne 12 ]]; then
    msg="but $# were supplied."
    [[ $# -eq 1 ]] && msg="but only $# was supplied."
    cat << EOM
Error: $(basename "${0}") requires 12 positional arguments, ${msg}

The necessary positional arguments:
${show_help}
EOM
    exit 1
fi

#  Parse positional arguments, assigning them to variables
#+ 
#+ Most of the argument inputs are not checked, as this is performed by
#+ 'execute_*.sh' and, to a certain extent, the corresponding Python script
env_nam="${1}"
scr_sig="${2}"
threads="${3}"
str_infile="${4}"
str_outfile="${5}"
siz_bin="${6}"
typ_cvg="${7}"
str_scl_fct="${8}"
str_usr_frg="${9}"
rnd="${10}"
err_out="${11}"
nam_job="${12}"

#  Debug argument variable assignments
if ${debug:-false}; then
    debug_var \
        "env_nam=${env_nam}" \
        "scr_sig=${scr_sig}" \
        "threads=${threads}" \
        "str_infile=${str_infile}" \
        "str_outfile=${str_outfile}" \
        "siz_bin=${siz_bin}" \
        "typ_cvg=${typ_cvg}" \
        "str_scl_fct=${str_scl_fct}" \
        "str_usr_frg=${str_usr_frg}" \
        "rnd=${rnd}" \
        "err_out=${err_out}" \
        "nam_job=${nam_job}"
fi

#  Activate environment
activate_env "${env_nam}" || exit 1

#  Reconstruct arrays from serialized strings
IFS=","
read -r -a arr_infile  <<< "${str_infile}"
read -r -a arr_outfile <<< "${str_outfile}"
read -r -a arr_scl_fct <<< "${str_scl_fct}"
read -r -a arr_usr_frg <<< "${str_usr_frg}"
unset IFS

#  Debug output to check number of array elements and array element values
if ${debug:-false}; then
    echo "\${#arr_infile[@]}=${#arr_infile[@]}"   && echo ""
    echo "arr_infile=( ${arr_infile[*]} )"        && echo ""
    echo "\${#arr_outfile[@]}=${#arr_outfile[@]}" && echo ""
    echo "arr_outfile=( ${arr_outfile[*]} )"      && echo ""
    echo "\${#arr_scl_fct[@]}=${#arr_scl_fct[@]}" && echo ""
    echo "arr_scl_fct=( ${arr_scl_fct[*]} )"      && echo ""
    echo "\${#arr_usr_frg[@]}=${#arr_usr_frg[@]}" && echo ""
    echo "arr_usr_frg=( ${arr_usr_frg[*]} )"      && echo ""
fi

#  Determine and run mode: SLURM or GNU Parallel/serial
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    #  Mode: SLURM
    id_job=${SLURM_ARRAY_JOB_ID}
    id_tsk=${SLURM_ARRAY_TASK_ID}

    if [[ "${id_tsk}" -lt 1 ]]; then
        echo "Error: SLURM task ID is invalid: ${id_tsk}" >&2
        exit 1
    else
        idx=$(( id_tsk - 1 ))
    fi

    #  Retrieve input and output files, and related values, by indexing into
    #+ reconstructed arrays
    infile="${arr_infile[idx]}"
    outfile="${arr_outfile[idx]}"
    scl_fct="${arr_scl_fct[idx]}"
    usr_frg="${arr_usr_frg[idx]}"

    if ${debug:-false}; then
        debug_var \
            "infile=${infile}" \
            "outfile=${outfile}" \
            "scl_fct=${scl_fct}" \
            "usr_frg=${usr_frg}"
    fi

    #  Run function to debug and validate 'infile', 'outfile', 'scl_fct', and
    #+ 'usr_frg', using it to assign 'samp'
    IFS=';' read -r samp dsc < <(
        process_input "${infile}" "${outfile}" "${scl_fct}" "${usr_frg}"
    ) || exit 1
    unset IFS

    if ${debug:-false}; then echo "samp=${samp}" && echo ""; fi

    #  Run function to set SLURM and symlinked log files
    IFS=';' read -r err_ini out_ini err_dsc out_dsc < <(
        set_logs_slurm \
        "${id_job}" "${id_tsk}" "${samp}" "${err_out}" "${nam_job}"
    ) || exit 1
    unset IFS

    if ${debug:-false}; then
        debug_var \
            "err_ini=${err_ini}" \
            "out_ini=${out_ini}" \
            "err_dsc=${err_dsc}" \
            "out_dsc=${out_dsc}"
    fi

    #  Set optional arguments if applicable
    set_args_opt  # Defines and assigns array 'optional'
    echo "set_args_opt(): optional=( ${optional[*]} )"
    echo ""

    #  Run Python script to compute signal
    if ${debug:-false}; then
        echo "python ${scr_sig} \\"
        echo "    --verbose \\"
        echo "    --threads ${threads} \\"
        echo "    --infile ${infile} \\"
        echo "    --outfile ${outfile} \\"
        echo "    --siz_bin ${siz_bin} \\"
        echo "    ${optional[*]}"
        echo ""
    fi

    #TODO: Modularize this
    python "${scr_sig}" \
        --verbose \
        --threads "${threads}" \
        --infile "${infile}" \
        --outfile "${outfile}" \
        --siz_bin "${siz_bin}" \
        "${optional[@]}"

    rm "${err_ini}" "${out_ini}"
else
    #  Mode: GNU Parallel/serial
    for idx in "${!arr_infile[@]}"; do
        infile="${arr_infile[idx]}"
        outfile="${arr_outfile[idx]}"
        scl_fct="${arr_scl_fct[idx]}"
        usr_frg="${arr_usr_frg[idx]}"

        #  Debug and validate variable assignments
        if ${debug:-false}; then
            debug_var \
                "infile=${infile}" \
                "outfile=${outfile}" \
                "scl_fct=${scl_fct}" \
                "usr_frg=${usr_frg}"
        fi

        #  Run function to debug and validate 'infile', 'outfile', 'scl_fct', and
        #+ 'usr_frg', using it to assign 'samp'
        IFS=';' read -r samp dsc < <(
            process_input "${infile}" "${outfile}" "${scl_fct}" "${usr_frg}"
        ) || exit 1
        unset IFS

        if ${debug:-false}; then echo "samp=${samp}" && echo ""; fi

        #  Set optional arguments if applicable
        set_args_opt  # Defines and assigns array 'optional'
        echo "set_args_opt(): optional=( ${optional[*]} )"
        echo ""

        #  Run Python script to compute signal
        if ${debug:-false}; then
            echo "python ${scr_sig} \\"
            echo "    --verbose \\"
            echo "    --threads ${threads} \\"
            echo "    --infile ${infile} \\"
            echo "    --outfile ${outfile} \\"
            echo "    --siz_bin ${siz_bin} \\"
            echo "    ${optional[*]} \\"
            echo "         >> ${err_out}/${nam_job}.${dsc}.stdout.txt \\"
            echo "        2>> ${err_out}/${nam_job}.${dsc}.stderr.txt"
            echo ""
        fi

        #TODO: Modularize this
        python "${scr_sig}" \
            --verbose \
            --threads "${threads}" \
            --infile "${infile}" \
            --outfile "${outfile}" \
            --siz_bin "${siz_bin}" \
            "${optional[@]}" \
                 >> "${err_out}/${nam_job}.${dsc}.stdout.txt" \
                2>> "${err_out}/${nam_job}.${dsc}.stderr.txt"
    done
fi
