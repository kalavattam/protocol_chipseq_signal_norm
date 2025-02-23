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
    local env_nam="${1}"
    if [[ "${CONDA_DEFAULT_ENV}" != "${env_nam}" ]]; then
        eval "$(conda shell.bash hook)"
        conda activate "${env_nam}" || 
            {
                echo "Error: Failed to activate environment '${env_nam}'." >&2
                return 1
            }
    fi
}


#  Function to validate that a variable is not empty or unset
function validate_var() {
    local var_nam="${1}"
    local var_val="${2}"
    local idx="${3:-0}"  # Default to '0' if no index is given

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


#  Function to source a function script for SLURM submission
function source_function() {
    local fil_scr="${1}"

    if [[ ! -f "${fil_scr}" ]]; then
        echo "Error: Function script not found: '${fil_scr}'." >&2
        return 1
    fi

    # shellcheck disable=SC1090
    if ! source "${fil_scr}"; then
        echo "Error: Failed to source function: '${fil_scr}'." >&2
        return 1
    fi
}


#  Function to set up SLURM log file paths and create symlinked log files
function setup_slurm_logs() {
    local id_job="${1}"
    local id_tsk="${2}"
    local samp="${3}"
    local err_out="${4}"
    local nam_job="${5}"
    local err_ini out_ini err_dsc out_dsc

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
    local infile="${1}"
    local scr_fnc="${2}"
    local dir_out="${3}"
    local samp nam_fnc outfile

    validate_var "infile"   "${infile}"  || return 1
    validate_var "scr_fnc"  "${scr_fnc}" || return 1

    #  Derive sample and function names from argument variable assignments
    samp="$(basename "${infile}" ".bam")"
    nam_fnc="$(basename "${scr_fnc}" ".sh")"

    #  Perform pattern matching to assign an outfile name
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
    local nam_fnc="${1}"
    local threads="${2}"
    local infile="${3}"
    local outfile="${4}"
    local mito="${5}"
    local tg="${6}"
    local mtr="${7}"
    local chk_chr="${8}"

    # shellcheck disable=SC2046
    "${nam_fnc}" \
        --threads "${threads}" \
        --infile "${infile}" \
        --outfile "${outfile}" \
        $(if ${mito}; then echo "--mito"; fi) \
        $(if ${tg}; then echo "--tg"; fi) \
        $(if ${mtr}; then echo "--mtr"; fi) \
        $(if ${chk_chr}; then echo "--chk_chr"; fi)
}


#  Parse keyword arguments, assigning them to variables
#+ 
#+ Most of the argument inputs are not checked, as this is performed by
#+ execute_*.sh and the function submitted to SLURM, "${scr_fnc}"
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

if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    exit 0
fi

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
if ${debug}; then
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
if ${debug}; then
    echo "\${#arr_infile[@]}=${#arr_infile[@]}" && echo ""
    echo "arr_infile=( ${arr_infile[*]} )"      && echo ""
fi

#  Source function to submit to SLURM
source_function "${scr_fnc}" || exit 1

#  Determine and run mode: SLURM or GNU Parallel/serial
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    #  Mode: SLURM
    id_job="${SLURM_ARRAY_JOB_ID}"
    id_tsk="${SLURM_ARRAY_TASK_ID}"
    idx=$(( id_tsk - 1 ))

    #  Retrieve the input file by indexing into the reconstructed input file
    #+ array
    infile="${arr_infile[idx]}"

    if ${debug}; then debug_var "infile=${infile}"; fi

    #  Run function to debug and validate 'infile', using it with 'scr_fnc' to
    #+ assign values to variables 'samp', 'nam_fnc', and 'outfile'
    IFS=';' read -r samp nam_fnc outfile < <(
        process_infile_function "${infile}" "${scr_fnc}" "${dir_out}"
    ) || exit 1
    unset IFS

    if ${debug}; then
        debug_var \
            "samp=${samp}" \
            "nam_fnc=${nam_fnc}" \
            "outfile=${outfile}"
    fi

    #  Run function to set SLURM and symlinked log files
    IFS=';' read -r err_ini out_ini err_dsc out_dsc < <(
        setup_slurm_logs \
        "${id_job}" "${id_tsk}" "${samp}" "${err_out}" "${nam_job}"
    ) || exit 1
    unset IFS

    if ${debug}; then
        debug_var \
            "err_ini=${err_ini}" \
            "out_ini=${out_ini}" \
            "err_dsc=${err_dsc}" \
            "out_dsc=${out_dsc}"
    fi

    if ! \
        run_filtering \
            "${nam_fnc}" "${threads}" "${infile}" "${outfile}" "${mito}" \
            "${tg}" "${mtr}" "${chk_chr}"
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

        if ${debug}; then debug_var "infile=${infile}"; fi

        #  Run function to debug and validate 'infile', using it with 'scr_fnc'
        #+ to assign values to variables 'samp', 'nam_fnc', and 'outfile'
        IFS=';' read -r samp nam_fnc outfile < <(
            process_infile_function "${infile}" "${scr_fnc}" "${dir_out}"
        ) || exit 1
        unset IFS

        if ${debug}; then
            debug_var \
                "samp=${samp}" \
                "nam_fnc=${nam_fnc}" \
                "outfile=${outfile}"
        fi

        if ! \
            run_filtering \
                "${nam_fnc}" "${threads}" "${infile}" "${outfile}" "${mito}" \
                "${tg}" "${mtr}" "${chk_chr}"
        then
            echo "Error: Failed to filter BAM file: '${infile}'." >&2
            exit 1
        fi
    done
fi
