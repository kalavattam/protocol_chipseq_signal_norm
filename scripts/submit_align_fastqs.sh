#!/bin/bash

#  submit_align_fastqs.sh
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


#  Function to debug, validate, and parse 'infile' to output variable
#+ assignments: 'fq_1', 'fq_2', 'samp' (sample name), and 'outfile'
function process_infile() {
    local infile="${1}"
    local sfx_se="${2}"
    local sfx_pe="${3}"
    local dir_out="${4}"
    local fq_1 fq_2 samp outfile

    validate_var "infile"  "${infile}" || return 1

    if [[ "${infile}" == *,* ]]; then
        fq_1="${infile%%,*}"
        fq_2="${infile#*,}"
        samp="$(basename "${fq_1%%"${sfx_pe}"}")"
    else
        fq_1="${infile}"
        unset fq_2
        samp="$(basename "${fq_1%%"${sfx_se}"}")"
    fi

    validate_var "fq_1" "${fq_1}" || return 1
    if [[ -n "${fq_2}" ]]; then validate_var "fq_2" "${fq_2}" || return 1; fi

    #  Assign the outfile
    outfile="${dir_out}/${samp}.bam"

    #  Return values
    echo "${fq_1};${fq_2};${samp};${outfile}"
}


#  Function to execute alignment using the specified function 'align_fastqs'
function run_alignment() {
    local scr_fnc="${1}"
    local threads="${2}"
    local aligner="${3}"
    local a_type="${4}"
    local mapq="${5}"
    local req_flg="${6}"
    local index="${7}"
    local fq_1="${8}"
    local fq_2="${9}"
    local outfile="${10}"
    local qname="${11}"
    local nam_fnc

    nam_fnc="$(basename "${scr_fnc}" .sh)"
    
    # shellcheck disable=SC2046
    "${nam_fnc}" \
        --threads "${threads}" \
        --aligner "${aligner}" \
        $(if [[ -n "${a_type}" ]]; then echo "--a_type ${a_type}"; fi) \
        $(if [[ "${mapq}" -gt 0 ]]; then echo "--mapq ${mapq}"; fi) \
        $(if ${req_flg:-false}; then echo "--req_flg"; fi) \
        --index "${index}" \
        --fq_1 "${fq_1}" \
        $(if [[ -n "${fq_2}" ]]; then echo "--fq_2 ${fq_2}"; fi) \
        --outfile "${outfile}" \
        $(if ${qname:-false}; then echo "--qname"; fi)
}


#  Parse keyword arguments, assigning them to variables
#+ 
#+ Most of the argument inputs are not checked, as this is performed by
#+ execute_*.sh and the function submitted to SLURM, "${scr_fnc}"
env_nam="env_protocol"
scr_fnc=""
threads=4
aligner="bowtie2"
a_type="end-to-end"
mapq=""
req_flg=false
index=""
str_infile=""
dir_out=""
qname=false
sfx_se=""
sfx_pe=""
err_out=""
nam_job="align_fastqs"

#  Define the help message
show_help=$(cat << EOM
$(basename "${0}") takes the following keyword arguments:
  -en, --env_nam     Conda/Mamba environment to activate.
  -sf, --scr_fnc     Path and name of function to source and submit to SLURM.
   -t, --threads     Number of threads to use.
   -a, --aligner     Alignment program to use.
  -at, --a_type      Alignment type.
  -mq, --mapq        MAPQ threshold for filtering BAM outfiles.
  -rf, --req_flg     (flag) Require flag bit 2 for filtering BAM outfiles.
  -ix, --index       Path to directory containing aligner index.
   -i, --str_infile  Semicolon-delimited serialized string of FASTQ infiles.
  -do, --dir_out     Directory to write BAM outfiles.
  -qn, --qname       (flag) Retain queryname-sorted intermediate BAM files.
  -ss, --sfx_se      Suffix to strip from SE FASTQ files.
  -sp, --sfx_pe      Suffix to strip from PE FASTQ files.
  -eo, --err_out     Directory to store stderr and stdout outfiles.
  -nj, --nam_job     Name of job.

All arguments are required with the following notes and exceptions:
  - '--env_nam' defaults to 'env_nam=${env_nam}' if not specified.
  - '--threads' defaults to 'threads=${threads}' if not specified.
  - '--aligner' defaults to 'aligner=${aligner}' if not specified.
  - '--a_type' defaults to 'a_type=${a_type}' if not specified.
  - '--a_type' is only required if 'aligner=bowtie2'; if 'aligner=bwa', then
    the assignment is ignored.
  - '--req_flg' and '--qname' are optional.
  - '--nam_job' defaults to 'nam_job=${nam_job}' if not specified.
EOM
)

#  Display help message if no arguments or help option is given
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
         -a|--aligner)
            aligner="$(echo "${2}" | tr '[:upper:]' '[:lower:]')"
            shift 2
            ;;
        -at|--a_type)
            a_type="$(echo "${2}" | tr '[:upper:]' '[:lower:]')"
            shift 2
            ;;
        -mq|--mapq)       mapq="${2}";       shift 2 ;;
        -rf|--req_flg)    req_flg=true;      shift 1 ;;
        -ix|--index)      index="${2}";      shift 2 ;;
         -i|--str_infile) str_infile="${2}"; shift 2 ;;
        -do|--dir_out)    dir_out="${2}";    shift 2 ;;
        -qn|--qname)      qname=true;        shift 1 ;;
        -ss|--sfx_se)     sfx_se="${2}";     shift 2 ;;
        -sp|--sfx_pe)     sfx_pe="${2}";     shift 2 ;;
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
        "aligner=${aligner}" \
        "a_type=${a_type}" \
        "mapq=${mapq}" \
        "req_flg=${req_flg}" \
        "index=${index}" \
        "str_infile=${str_infile}" \
        "dir_out=${dir_out}" \
        "qname=${qname}" \
        "sfx_se=${sfx_se}" \
        "sfx_pe=${sfx_pe}" \
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

    #  Run function to validate 'infile', using it to assign values to
    #+ variables 'fq_1', 'fq_2', 'samp', and 'outfile'
    IFS=';' read -r fq_1 fq_2 samp outfile < <(
        process_infile "${infile}" "${sfx_se}" "${sfx_pe}" "${dir_out}"
    ) || exit 1
    unset IFS

    if ${debug}; then
        debug_var \
            "fq_1=${fq_1}" \
            "fq_2=${fq_2}" \
            "samp=${samp}" \
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

    #  Run 'align_fastqs'
    if ! \
        run_alignment \
            "${scr_fnc}" "${threads}" "${aligner}" "${a_type}" "${mapq}" \
            "${req_flg}" "${index}" "${fq_1}" "${fq_2}" "${outfile}" "${qname}"
    then
        echo "Error: Failed to perform alignment, etc." >&2
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

        #  Run function to debug and validate 'infile', using it to assign
        #+ values to variables 'fq_1', 'fq_2', 'samp', and 'outfile'
        IFS=';' read -r fq_1 fq_2 samp outfile < <(
            process_infile "${infile}" "${sfx_se}" "${sfx_pe}" "${dir_out}"
        ) || exit 1
        unset IFS

        if ${debug}; then
            debug_var \
                "fq_1=${fq_1}" \
                "fq_2=${fq_2}" \
                "samp=${samp}" \
                "outfile=${outfile}"
        fi

        #  Run 'align_fastqs'
        if ! \
            run_alignment \
                "${scr_fnc}" "${threads}" "${aligner}" "${a_type}" "${mapq}" \
                "${req_flg}" "${index}" "${fq_1}" "${fq_2}" "${outfile}" "${qname}"
        then
            echo "Error: Failed to perform alignment, etc." >&2
            exit 1
        fi
    done
fi
