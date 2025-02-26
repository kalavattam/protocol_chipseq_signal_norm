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
    local var_nam="${1}"  # Variable name (for error messages)
    local var_val="${2}"  # Value to check for emptiness/unset state
    local idx="${3:-0}"   # Optional index (for arrays); defaults to '0'

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


#  Function to debug, validate, and parse 'infile' into 'fq_1', 'fq_2', 'samp'
#+ (sample name), and 'outfile'
function process_infile() {
    local infile="${1}"   # Input FASTQ file(s)
    local sfx_se="${2}"   # Suffix for SE FASTQ files
    local sfx_pe="${3}"   # Suffix for PE FASTQ files
    local dir_out="${4}"  # Directory for output files
    local fq_1     # FASTQ file #1 (SE or PE read 1)
    local fq_2     # FASTQ file #2 (PE read 2; '#N/A' if SE)
    local samp     # Sample name derived from FASTQ filename
    local outfile  # Output file

    #  Validate input arguments
    validate_var "infile"   "${infile}"  || return 1
    validate_var "sfx_se"   "${sfx_se}"  || return 1
    validate_var "sfx_pe"   "${sfx_pe}"  || return 1
    validate_var "dir_out"  "${dir_out}" || return 1

    #  Parse input FASTQ file(s) to assign 'fq_1', 'fq_2', and 'samp'
    if [[ "${infile}" == *,* ]]; then
        fq_1="${infile%%,*}"
        fq_2="${infile#*,}"
        samp="$(basename "${fq_1%%"${sfx_pe}"}")"
    else
        fq_1="${infile}"
        fq_2="#N/A"
        samp="$(basename "${fq_1%%"${sfx_se}"}")"
    fi

    #  Validate parsed FASTQ file paths
    validate_var "fq_1" "${fq_1}" || return 1
    
    if [[ "${fq_2}" != "#N/A" ]]; then
        validate_var "fq_2" "${fq_2}" || return 1
    fi

    #  Assign outfile
    outfile="${dir_out}/${samp}.bam"

    #  Return values
    echo "${fq_1};${fq_2};${samp};${outfile}"
}


#  Function to execute alignment using the specified function 'align_fastqs'
function run_alignment() {
    local scr_fnc="${1}"   # Script containing the function to execute
    local threads="${2}"   # Number of threads
    local aligner="${3}"   # Aligner program: 'bowtie2' or 'bwa'
    local a_type="${4}"    # Alignment type: 'local', 'global', '#N/A'
    local mapq="${5}"      # MAPQ threshold
    local req_flg="${6}"   # Flag: Require flag bit 2 (properly paired reads)
    local index="${7}"     # Index path
    local fq_1="${8}"      # First FASTQ file
    local fq_2="${9}"      # Second FASTQ file (optional)
    local outfile="${10}"  # Output BAM file
    local qname="${11}"    # Flag: Retain queryname-sorted BAM
    local err_out="${12}"  # Directory for logging stderr/stdout
    local nam_job="${13}"  # Job name for log file naming
    local samp="${14}"     # Sample name for log file naming
    local nam_fnc  # Name of function to execute (derived from 'scr_fnc')
    local log_out  # 'nam_fnc' stdout log file
    local log_err  # 'nam_fnc' stderr log file

    #  Extract function name from 'scr_fnc', and define paths for log output
    #+ files
    nam_fnc="$(basename "${scr_fnc}" .sh)"
    log_out="${err_out}/${nam_job}.${samp}.stdout.txt"
    log_err="${err_out}/${nam_job}.${samp}.stderr.txt"

    #  Run alignment function with specified parameters
    # shellcheck disable=SC2046
    if ! \
        "${nam_fnc}" \
            --threads "${threads}" \
            --aligner "${aligner}" \
            --a_type "${a_type}" \
            $(if [[ "${mapq}" -gt 0 ]]; then echo "--mapq ${mapq}"; fi) \
            $(if ${req_flg:-false}; then echo "--req_flg"; fi) \
            --index "${index}" \
            --fq_1 "${fq_1}" \
            $(
                if [[ "${fq_2}" != "#N/A" ]]; then
                    echo "--fq_2 ${fq_2}"
                fi
            ) \
            --outfile "${outfile}" \
            $(if ${qname:-false}; then echo "--qname"; fi) \
                 > "${log_out}" \
                2> "${log_err}"
    then
        echo \
            "Error: Alignment failed for sample '${samp}'. See log:" \
            "'${log_err}'" >&2
        return 1
    fi
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
  - '--req_flg' and '--qname' are optional.
  - '--nam_job' defaults to 'nam_job=${nam_job}' if not specified.
EOM
)

#  Display help message if no arguments, or help option, is given
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

    if ${debug}; then debug_var "infile=${infile}"; fi

    #  Run function to debug and validate 'infile', using it to assign values
    #+ to variables 'fq_1', 'fq_2', 'samp', and 'outfile'
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
        set_logs_slurm \
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
            "${req_flg}" "${index}" "${fq_1}" "${fq_2}" "${outfile}" \
            "${qname}" "${err_out}" "${nam_job}" "${samp}"
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
                "${req_flg}" "${index}" "${fq_1}" "${fq_2}" "${outfile}" \
                "${qname}" "${err_out}" "${nam_job}" "${samp}"
        then
            echo "Error: Failed to perform alignment, etc." >&2
            exit 1
        fi
    done
fi
