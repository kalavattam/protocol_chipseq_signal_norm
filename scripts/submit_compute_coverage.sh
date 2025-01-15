#!/bin/bash

#  If true, run script in debug mode
debug=true

#  Define functions
function debug_var() {
    for var in "$@"; do
        echo "${var}"
        echo ""
    done
}


function validate_var() {
    local var_nam=${1}
    local var_val=${2}
    local     idx=${3}  # Pass index explicitly
    if [[ -z "${var_val}" ]]; then
        echo \
            "Error: '${var_nam}' is empty or unset for array index" \
            "'${idx}'." >&2
        return 1
    fi
}


function set_logs() {
    local  job_id=${1}
    local task_id=${2}
    local nam_smp=${3}
    local dir_log=${4}
    err_ini="${dir_log}/${nam_job}.${job_id}-${task_id}.stderr.txt"
    out_ini="${dir_log}/${nam_job}.${job_id}-${task_id}.stdout.txt"
    err_dsc="${dir_log}/${nam_job}.${nam_smp}.${job_id}-${task_id}.stderr.txt"
    out_dsc="${dir_log}/${nam_job}.${nam_smp}.${job_id}-${task_id}.stdout.txt"
}


function set_args_opt() {
    unset optional && typeset -a optional
    if [[ -n "${typ_cvg}" ]]; then
        optional+=( --typ_cvg "${typ_cvg}" )
    fi

    if [[ -n "${scl_fct}" && "${scl_fct}" != "#N/A" ]]; then
        optional+=( --scl_fct "${scl_fct}" )
    fi

    if [[ "${usr_frg}" != "#N/A" ]]; then
        optional+=( --usr_frg "${usr_frg}" )
    fi
}


#  Display help message if no arguments or help option is given
show_help=$(cat << EOM
\${1}=env_nam       # str: Name of Conda/Mamba environment to activate
\${2}=scr_cvg       # str: Path to coverage script
\${3}=threads       # int: Number of threads to use
\${4}=str_infile    # str: Comma-separated string of infiles
\${5}=str_outfile   # str: Comma-separated string of outfile stems
\${6}=typ_out       # str: Outfile type: 'bedgraph', 'bigwig', or 'both'
\${7}=siz_bin       # int: Bin size in base pairs
\${8}=typ_cvg       # str: Type of coverage to compute
\${9}=str_scl_fct   # flt: Comma-separated string of scaling factors
\${10}=str_usr_frg  # int: Comma-separated string of fragment lengths
\${11}=err_out      # str: Directory for stdout and stderr files
\${12}=nam_job      # str: Name of job
EOM
)

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
#+ execute_*.sh and, to a certain extent, the corresponding Python script
env_nam="${1}"
scr_cvg="${2}"
threads="${3}"
str_infile="${4}"
str_outfile="${5}"
typ_out="${6}"
siz_bin="${7}"
typ_cvg="${8}"
str_scl_fct="${9}"
str_usr_frg="${10}"
err_out="${11}"
nam_job="${12}"

#  Activate environment
if [[ "${CONDA_DEFAULT_ENV}" != "${env_nam}" ]]; then
    eval "$(conda shell.bash hook)"
    conda activate "${env_nam}" ||
        {
            echo "Error: Failed to activate environment '${env_nam}'." >&2
            exit 1
        }
fi

#  Reconstruct arrays from serialized strings
IFS=',' read -r -a arr_infile <<< "${str_infile}"
IFS=',' read -r -a arr_outfile <<< "${str_outfile}"
IFS=',' read -r -a arr_scl_fct <<< "${str_scl_fct}"
IFS=',' read -r -a arr_usr_frg <<< "${str_usr_frg}"

#  Debug output to check number of array elements and array element values
if ${debug}; then
    echo "\${#arr_infile[@]}=${#arr_infile[@]}"
    echo ""
    echo "arr_infile=( ${arr_infile[*]} )"
    echo ""
    echo "\${#arr_outfile[@]}=${#arr_outfile[@]}"
    echo ""
    echo "arr_outfile=( ${arr_outfile[*]} )"
    echo ""
    echo "\${#arr_scl_fct[@]}=${#arr_scl_fct[@]}"
    echo ""
    echo "arr_scl_fct=( ${arr_scl_fct[*]} )"
    echo ""
    echo "\${#arr_usr_frg[@]}=${#arr_usr_frg[@]}"
    echo ""
    echo "arr_usr_frg=( ${arr_usr_frg[*]} )"
    echo ""
fi

#  Determine and run mode: SLURM or GNU Parallel/serial
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    #  Mode: SLURM
    id_job=${SLURM_ARRAY_JOB_ID}
    id_tsk=${SLURM_ARRAY_TASK_ID}
    idx=$(( id_tsk - 1 ))

    infile="${arr_infile[idx]}"
    outfile="${arr_outfile[idx]}"
    scl_fct="${arr_scl_fct[idx]}"
    usr_frg="${arr_usr_frg[idx]}"

    #  Debug and validate variable assignments
    if ${debug}; then
        debug_var \
            "infile=${infile}" "outfile=${outfile}" \
            "scl_fct=${scl_fct}" "usr_frg=${usr_frg}"
    fi
    validate_var "infile"  "${infile}"  "${idx}" || exit 1
    validate_var "outfile" "${outfile}" "${idx}" || exit 1
    validate_var "scl_fct" "${scl_fct}" "${idx}" || exit 1
    validate_var "usr_frg" "${usr_frg}" "${idx}" || exit 1

    samp="${outfile##*/}"
    samp="${samp%.bam}"

    #  Debug sample name
    if ${debug}; then echo "samp=${samp}" && echo ""; fi

    #  Run subroutine to set SLURM and symlinked/better-named log files
    set_logs "${id_job}" "${id_tsk}" "${samp}" "${err_out}"
    ln -f "${err_ini}" "${err_dsc}"
    ln -f "${out_ini}" "${out_dsc}"

    #  Set optional arguments if applicable
    set_args_opt  # Subroutine defines/assigns array 'optional'

    python "${scr_cvg}" \
        --verbose \
        --threads "${threads}" \
        --infile "${infile}" \
        --outfile "${outfile}" \
        --typ_out "${typ_out}" \
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
        if ${debug}; then
            debug_var \
                "infile=${infile}" "outfile=${outfile}" \
                "scl_fct=${scl_fct}" "usr_frg=${usr_frg}"
        fi
        validate_var "infile"  "${infile}"  "${idx}" || exit 1
        validate_var "outfile" "${outfile}" "${idx}" || exit 1
        validate_var "scl_fct" "${scl_fct}" "${idx}" || exit 1
        validate_var "usr_frg" "${usr_frg}" "${idx}" || exit 1

        #  Set optional arguments if applicable
        set_args_opt  # Subroutine defines/assigns array 'optional'

        python "${scr_cvg}" \
            --verbose \
            --threads "${threads}" \
            --infile "${infile}" \
            --outfile "${outfile}" \
            --typ_out "${typ_out}" \
            --siz_bin "${siz_bin}" \
            "${optional[@]}" \
                 >> "${err_out}/${nam_job}.$(basename "${outfile}").stdout.txt" \
                2>> "${err_out}/${nam_job}.$(basename "${outfile}").stderr.txt"
    done
fi
