#!/bin/bash

#  Display help message if no arguments or help option is given
show_help=$(cat << EOM
\${1}=env_nam   # str: Name of Conda/Mamba environment to activate
\${2}=scr_cvg   # str: Path to coverage script
\${3}=threads   # int: Number of threads to use
\${4}=infile    # str: Infile, including path
\${5}=outfile   # str: Outfile stem, including path
\${6}=typ_out   # str: Outfile type: 'bedgraph', 'bigwig', or 'both'
\${7}=bin_siz   # int: Bin size in base pairs
\${8}=norm      # bol: Calculate 'normalized coverage' (PMID: 37160995)
\${9}=raw       # bol: Calculate 'raw coverage'
\${10}=scl_fct  # flt: User-supplied scaling factor for coverage calculation
\${11}=usr_frg  # int: User-supplied fragment length in base pairs
\${12}=err_out  # str: Directory, including path, for stdout and stderr files
\${13}=nam_job  # str: Name of job
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
infile="${4}"
outfile="${5}"
typ_out="${6}"
bin_siz="${7}"
norm="${8}"
raw="${9}"
scl_fct="${10}"
usr_frg="${11}"
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
    ) \
         >> "${err_out}/${nam_job}.$(basename "${outfile}").stdout.txt" \
        2>> "${err_out}/${nam_job}.$(basename "${outfile}").stderr.txt"
