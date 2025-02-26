#!/bin/bash

#  submit_filter_bams.sh
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


#  Function to safely source a function script required for SLURM submission
function source_function() {
    local fil_scr="${1}"  # Function script to source

    if [[ ! -f "${fil_scr}" ]]; then
        echo "Error: Function script not found: '${fil_scr}'." >&2
        return 1
    fi

    # shellcheck disable=SC1090
    if ! \
        source "${fil_scr}"
    then
        echo "Error: Failed to source function: '${fil_scr}'." >&2
        return 1
    fi
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
    local err_dsc  # Symlinked stderr log with descriptive name
    local out_dsc  # Symlinked stdout log with descriptive name

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


#  Function to parse 'infile' and 'scr_fnc' into 'samp' (sample name),
#+ 'nam_fnc' (name of function to run), and 'outfile'
function process_infile_function() {
    local infile="${1}"   # Input BAM file
    local scr_fnc="${2}"  # Function script
    local dir_out="${3}"  # Directory for output BAM files
    local samp     # Sample name derived from infile
    local nam_fnc  # Function name derived from 'scr_fnc'
    local outfile  # Output BAM file

    #  Validate input arguments
    validate_var "infile"   "${infile}"  || return 1
    validate_var "scr_fnc"  "${scr_fnc}" || return 1
    validate_var "dir_out"  "${dir_out}" || return 1

    #  Extract sample and function names from input values
    samp="$(basename "${infile}" ".bam")"
    nam_fnc="$(basename "${scr_fnc}" ".sh")"

    #  Assign outfile name based on function type: 'sc' or 'sp'
    if [[ ${nam_fnc} =~ "sc" ]]; then
        outfile="${dir_out}/${samp}.sc.bam"
    elif [[ ${nam_fnc} =~ "sp" ]]; then
        outfile="${dir_out}/${samp}.sp.bam"
    else
        echo "Error: Sample name could not be processed." >&2
        return 1
    fi

    #  Return values
    echo "${samp};${nam_fnc};${outfile}"
}


#  Function to execute filtering using the specified function
function run_filtering() {
    local nam_fnc="${1}"   # Name of function to run
    local threads="${2}"   # Number of threads
    local infile="${3}"    # Input BAM file
    local outfile="${4}"   # Output BAM file
    local mito="${5}"      # Retain mito. chr. (true/false)
    local tg="${6}"        # Retain SP_II_TG chr. (true/false)
    local mtr="${7}"       # Retain SP_MTR chr. (true/false)
    local chk_chr="${8}"   # Check chr. in output (true/false)
    local err_out="${9}"   # Directory for stderr and stdout logs
    local nam_job="${10}"  # Job name for log file naming
    local samp="${11}"     # Sample name for log file naming
    local log_out  # 'nam_fnc' stdout log file
    local log_err  # 'nam_fnc' stderr log file

    #  Define paths for log output files
    log_out="${err_out}/${nam_job}.${samp}.stdout.txt"
    log_err="${err_out}/${nam_job}.${samp}.stderr.txt"

    #  Run the filtering function and capture logs
    # shellcheck disable=SC2046
    "${nam_fnc}" \
        --threads "${threads}" \
        --infile "${infile}" \
        --outfile "${outfile}" \
        $(if ${mito}; then echo "--mito"; fi) \
        $(if ${tg}; then echo "--tg"; fi) \
        $(if ${mtr}; then echo "--mtr"; fi) \
        $(if ${chk_chr}; then echo "--chk_chr"; fi) \
             > "${log_out}" \
            2> "${log_err}"
}


#  Parse keyword arguments, assigning them to variables
#+ 
#+ Most of the argument inputs are not checked, as this is performed by
#+ 'execute_*.sh' and the function submitted to SLURM, 'scr_fnc'
env_nam="env_protocol"
scr_fnc=""
threads=4
str_infile=""
dir_out=""
mito=false
tg=false
mtr=false
chk_chr=false
err_out=""
nam_job=""

#  Define the help message
show_help=$(cat << EOM
$(basename "${0}") takes the following keyword arguments:
  -en, --env_nam     Conda/Mamba environment to activate.
  -sf, --scr_fnc     Path and name of function to source and submit to SLURM.
   -t, --threads     Number of threads to use.
   -i, --str_infile  Comma-separated serialized string of BAM infiles.
  -do, --dir_out     Directory to write BAM outfiles.
   -m, --mito        (flag) Retain mitochondrial chromosome in BAM outfiles.
  -tg, --tg          (flag) Retain SP_II_TG chromosome in BAM outfiles.
  -mr, --mtr         (flag) Retain SP_MTR chromosome in BAM outfiles.
  -cc, --chk_chr     (flag) Check chromosomes in BAM outfiles.
  -eo, --err_out     Directory to store stderr and stdout outfiles.
  -nj, --nam_job     Name of job.

All arguments are required with the following notes and exceptions:
  - '--env_nam' defaults to 'env_nam=${env_nam}' if not specified.
  - '--threads' default to 'threads=${threads}' if not specified.
  - '--mito', '--tg', '--mtr', and --chk_chr are optional.
EOM
)

#  Display help message if a help option or no arguments are given
if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    exit 0
fi

#  Parse keyword arguments
while [[ "$#" -gt 0 ]]; do
    case "${1}" in
        -en|--env_nam)    env_nam="${2}";    shift 2 ;;
        -sf|--scr_fnc)    scr_fnc="${2}";    shift 2 ;;
         -t|--threads)    threads="${2}";    shift 2 ;;
         -i|--str_infile) str_infile="${2}"; shift 2 ;;
        -do|--dir_out)    dir_out="${2}";    shift 2 ;;
         -m|--mito)       mito=true;         shift 1 ;;
        -tg|--tg)         tg=true;           shift 1 ;;
        -mr|--mtr)        mtr=true;          shift 1 ;;
        -cc|--chk_chr)    chk_chr=true;      shift 1 ;;
        -eo|--err_out)    err_out="${2}";    shift 2 ;;
        -nj|--nam_job)    nam_job="${2}";    shift 2 ;;
        *)
            echo "## Unknown parameter passed: ${1} ##" >&2
            echo "" >&2
            echo "${show_help}" >&2
            exit 1
            ;;
    esac
done

#  Debug argument variable assignments
if ${debug:-false}; then
    debug_var \
        "env_nam=${env_nam}" \
        "scr_fnc=${scr_fnc}" \
        "threads=${threads}" \
        "str_infile=${str_infile}" \
        "dir_out=${dir_out}" \
        "mito=${mito}" \
        "tg=${tg}" \
        "mtr=${mtr}" \
        "chk_chr=${chk_chr}" \
        "err_out=${err_out}" \
        "nam_job=${nam_job}"
fi

#  Activate environment
activate_env "${env_nam}" || exit 1

#  Reconstruct array from serialized string
validate_var "str_infile" "${str_infile}" || exit 1
IFS=';' read -r -a arr_infile <<< "${str_infile}"
unset IFS

#  Debug output to check number of array elements and array element values
if ${debug:-false}; then
    echo "\${#arr_infile[@]}=${#arr_infile[@]}" && echo ""
    echo "arr_infile=( ${arr_infile[*]} )"      && echo ""
fi

#  Source function to submit to SLURM
if [[ -z "${scr_fnc}" ]]; then
    echo "Error: '--scr_fnc' is required but not set." >&2
    exit 1
else
    source_function "${scr_fnc}" || exit 1
fi

#  Determine and run mode: SLURM or GNU Parallel/serial
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    #  Mode: SLURM
    id_job="${SLURM_ARRAY_JOB_ID}"
    id_tsk="${SLURM_ARRAY_TASK_ID}"

    if [[ "${id_tsk}" -lt 1 ]]; then
        echo "Error: SLURM task ID is invalid: ${id_tsk}" >&2
        exit 1
    else
        idx=$(( id_tsk - 1 ))
    fi

    #  Retrieve the input file by indexing into the reconstructed input file
    #+ array
    infile="${arr_infile[idx]}"

    if ${debug:-false}; then debug_var "infile=${infile}"; fi

    #  Run function to debug and validate 'infile', using it with 'scr_fnc' to
    #+ assign values to variables 'samp', 'nam_fnc', and 'outfile'
    IFS=';' read -r samp nam_fnc outfile < <(
        process_infile_function "${infile}" "${scr_fnc}" "${dir_out}"
    ) || exit 1
    unset IFS

    if ${debug:-false}; then
        debug_var \
            "samp=${samp}" \
            "nam_fnc=${nam_fnc}" \
            "outfile=${outfile}"
    fi

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

    #  Perform filtering
    if ! \
        run_filtering \
            "${nam_fnc}" "${threads}" "${infile}" "${outfile}" "${mito}" \
            "${tg}" "${mtr}" "${chk_chr}" "${err_out}" "${nam_job}" "${samp}"
    then
        echo "Error: Failed to filter BAM file: '${infile}'." >&2
        exit 1
    fi

    #  Remove the initial SLURM stderr and stdout TXT outfiles
    rm "${err_ini}" "${out_ini}"
else
    #  Mode: GNU Parallel/serial
    for idx in "${!arr_infile[@]}"; do
        #  Retrieve the input file by indexing into the reconstructed input
        #+ file array
        infile="${arr_infile[idx]}"

        if ${debug:-false}; then debug_var "infile=${infile}"; fi

        #  Run function to debug and validate 'infile', using it with 'scr_fnc'
        #+ to assign values to variables 'samp', 'nam_fnc', and 'outfile'
        IFS=';' read -r samp nam_fnc outfile < <(
            process_infile_function "${infile}" "${scr_fnc}" "${dir_out}"
        ) || exit 1
        unset IFS

        if ${debug:-false}; then
            debug_var \
                "samp=${samp}" \
                "nam_fnc=${nam_fnc}" \
                "outfile=${outfile}"
        fi

        #  Perform filtering
        if ! \
            run_filtering \
                "${nam_fnc}" "${threads}" "${infile}" "${outfile}" "${mito}" \
                "${tg}" "${mtr}" "${chk_chr}" "${err_out}" "${nam_job}" \
                "${samp}"
        then
            echo "Error: Failed to filter BAM file: '${infile}'." >&2
            exit 1
        fi
    done
fi
