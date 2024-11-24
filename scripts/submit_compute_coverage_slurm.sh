#!/bin/bash

#  If true, run script in debug mode
debug=true

#  Display help message if no arguments or help option is given
show_help=$(cat << EOM
\${1}=env_nam       # str: Name of Conda/Mamba environment to activate
\${2}=scr_cvg       # str: Path to coverage script
\${3}=threads       # int: Number of threads to use
\${4}=str_infile    # str: Comma-separated string of infiles
\${5}=str_outfile   # str: Comma-separated string of outfile stems
\${6}=typ_out       # str: Outfile type: 'bedgraph', 'bigwig', or 'both'
\${7}=bin_siz       # int: Bin size in base pairs
\${8}=norm          # bol: Calculate 'normalized coverage' (PMID: 37160995)
\${9}=raw           # bol: Calculate 'raw coverage'
\${10}=str_scl_fct  # flt: Comma-separated string of scaling factors
\${11}=str_usr_frg  # int: Comma-separated string of fragment lengths
\${12}=err_out      # str: Directory for stdout and stderr files
\${13}=nam_job      # str: Name of job
EOM
)

if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    cat << EOM
$(basename "${0}") requires 13 positional arguments:
${show_help}
EOM
    exit 0
fi

#  Check for exactly 13 positional arguments
if [[ $# -ne 13 ]]; then
    msg="but $# were supplied."
    [[ $# -eq 1 ]] && msg="but only $# was supplied."
    cat << EOM
Error: $(basename "${0}") requires 13 positional arguments, ${msg}

The necessary positional arguments:
${show_help}
EOM
    exit 1
fi

#  Parse positional arguments, assigning them to variables; most of the
#+ argument inputs are not checked, as this is performed by execute_*.sh and to
#+ a certain extent by the Python script submitted to SLURM
env_nam="${1}"
scr_cvg="${2}"
threads="${3}"
str_infile="${4}"
str_outfile="${5}"
typ_out="${6}"
bin_siz="${7}"
norm="${8}"
raw="${9}"
str_scl_fct="${10}"
str_usr_frg="${11}"
err_out="${12}"
nam_job="${13}"

#  Activate environment
if [[ "${CONDA_DEFAULT_ENV}" != "${env_nam}" ]]; then
    eval "$(conda shell.bash hook)"
    conda activate "${env_nam}" ||
        {
            echo "Error: Failed to activate environment '${env_nam}'." >&2
            exit 1
        }
fi

#  Check that SLURM environment variables are set
if [[ -z "${SLURM_ARRAY_JOB_ID}" ]]; then
    echo "Error: SLURM_ARRAY_JOB_ID is not set." >&2
    exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID}" ]]; then
    echo "Error: SLURM_ARRAY_TASK_ID is not set." >&2
    exit 1
fi

#  Give important SLURM environmental variables shorter names
id_job=${SLURM_ARRAY_JOB_ID}
id_tsk=${SLURM_ARRAY_TASK_ID}

#  Reconstruct arrays from serialized strings
IFS=',' read -r -a arr_infiles  <<< "${str_infile}"
IFS=',' read -r -a arr_outfiles <<< "${str_outfile}"
IFS=',' read -r -a arr_scl_fct  <<< "${str_scl_fct}"
IFS=',' read -r -a arr_usr_frg  <<< "${str_usr_frg}"

#  Debug output to check the task ID/array index, number of array elements, and
#+ array element values
if ${debug}; then
    echo "SLURM_ARRAY_TASK_ID=${id_tsk}"
    echo ""
    echo "\${#arr_infiles[@]}=${#arr_infiles[@]}"
    echo ""
    echo "arr_infiles=( ${arr_infiles[*]} )"
    echo ""
    echo "\${#arr_outfiles[@]}=${#arr_outfiles[@]}"
    echo ""
    echo "arr_outfiles=( ${arr_outfiles[*]} )"
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

#  Determine variable assignments from reconstructed arrays based on
#+ SLURM_ARRAY_TASK_ID
infile="${arr_infiles[$(( id_tsk - 1 ))]}"
outfile="${arr_outfiles[$(( id_tsk - 1 ))]}"
scl_fct="${arr_scl_fct[$(( id_tsk - 1 ))]}"
usr_frg="${arr_usr_frg[$(( id_tsk - 1 ))]}"

#  Debug variable assignments from reconstructed arrays
if ${debug}; then
    echo "infile=${infile}"
    echo ""
    echo "outfile=${outfile}"
    echo ""
    echo "scl_fct=${scl_fct}"
    echo ""
    echo "usr_frg=${usr_frg}"
    echo ""
fi

#  Exit if any variable assignment is empty
if [[ -z "${infile}" ]]; then
    echo \
        "Error: Failed to retrieve infile for id_tsk=${id_tsk}:" \
        "\${arr_infiles[$(( id_tsk - 1 ))]}." >&2
    exit 1
fi

if [[ -z "${outfile}" ]]; then
    echo \
        "Error: Failed to retrieve outfile for id_tsk=${id_tsk}:" \
        "\${arr_outfiles[$(( id_tsk - 1 ))]}." >&2
    exit 1
fi

if [[ -z "${scl_fct}" ]]; then
    echo \
        "Error: Failed to retrieve scl_fct for id_tsk=${id_tsk}:" \
        "\${arr_scl_fct[$(( id_tsk - 1 ))]}." >&2
    exit 1
fi

if [[ -z "${usr_frg}" ]]; then
    echo \
        "Error: Failed to retrieve usr_frg for id_tsk=${id_tsk}:" \
        "\${arr_usr_frg[$(( id_tsk - 1 ))]}." >&2
    exit 1
fi

#  Derive sample name from infile assignment
samp="${infile##*/}"
samp="${samp%.bam}"

#  Debug sample name
if ${debug}; then
    echo "samp=${samp}"
    echo ""
fi

#  Assign stdout and stderr outfiles to variables
err_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stderr.txt"
out_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stdout.txt"
err_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stderr.txt"
out_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stdout.txt"

#  Give the initial stderr and stdout TXT outfiles more descriptive names
ln -f "${err_ini}" "${err_dsc}"
ln -f "${out_ini}" "${out_dsc}"

#  Execute Python script to compute coverage
# shellcheck disable=SC2046
python "${scr_cvg}" \
    --verbose \
    --threads "${threads}" \
    --infile "${infile}" \
    --outfile "${outfile}" \
    --typ_out "${typ_out}" \
    --bin_siz "${bin_siz}" \
    $(
        if ! ${norm} && ! ${raw} && [[ "${scl_fct}" != "#N/A" ]]; then
            echo "--scl_fct ${scl_fct}"
        fi
    ) \
    $(
        if "${norm}" && ! "${raw}" && [[ "${scl_fct}" == "#N/A" ]]; then
            echo "--norm"
        fi
    ) \
    $(
        if [[ "${usr_frg}" != "#N/A" ]]; then
            echo "--usr_frg ${usr_frg}"
        fi
    )

#  Remove the initial SLURM stderr and stdout TXT outfiles
rm "${err_ini}" "${out_ini}"
