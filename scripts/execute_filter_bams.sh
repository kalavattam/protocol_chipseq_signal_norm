#!/bin/bash

#  execute_filter_bams.sh
#  KA


#  Run script in interactive/test mode (true) or command-line mode (false)
interactive=false

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive}; then set -euo pipefail; fi

#  Set the path to the "scripts" directory
# shellcheck disable=SC1091
if ${interactive}; then
    ## WARNING: Change path if you're not Kris and `interactive=true` ##
    dir_scr="${HOME}/tsukiyamalab/Kris/202X_protocol_ChIP/scripts"
else
    dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi


#  Source and define functions ================================================
#  Set the path to the "functions" directory
dir_fnc="${dir_scr}/functions"

# shellcheck disable=SC1091
{
    source "${dir_fnc}/check_supplied_arg.sh"
    source "${dir_fnc}/check_exists_file_dir.sh"
    source "${dir_fnc}/check_format_time.sh"
    source "${dir_fnc}/check_int_pos.sh"
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/echo_error.sh"
    source "${dir_fnc}/echo_warning.sh"
    source "${dir_fnc}/exit_0.sh"
    source "${dir_fnc}/exit_1.sh"
    source "${dir_fnc}/filter_bam_sc.sh"
    source "${dir_fnc}/filter_bam_sp.sh"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=1
infiles=""
dir_out=""
retain="sc"
mito=false
tg=false
mtr=false
chk_chr=false
err_out=""
slurm=false
max_job=6
time="0:30:00"

#  Initialize hardcoded argument variables
env_nam="env_align"
nam_job="filter_bams_${retain}"
scr_fnc="${dir_fnc}/${nam_job}.sh"
scr_sub="${dir_scr}/submit_filter_bams.sh"

#  Assign variable for help message
show_help=$(cat << EOM
execute_filter_bams.sh
  [--verbose] [--dry_run] --threads <int> --infiles <str> --dir_out <str>
  --retain <str> [--mito] [--tg] [--mtr] [--chk_chr] --err_out <str>
  --nam_job <str> [--slurm] [--max_job <int>] [--time <str>]

Description:
  execute_filter_bams.sh performs... #TODO

Arguments:
   -h, --help     Print this help message and exit.
   -v, --verbose  Run script in 'verbose' mode (optional).
  -dr, --dry_run  Run the command in 'check' mode (optional).
   -t, --threads  Number of threads to use (required; default: ${threads}).
   -i, --infiles  Comma-separated string vector of coordinate-sorted BAM
                  infiles (required).
  -do, --dir_out  The directory to store species-filtered and -reheadered BAM
                  outfiles (required).
   -r, --retain   Specify species chromosomes to retain: S. cerevisiae, "sc";
                  S. pombe, "sp" (required; default: ${retain}).
   -m, --mito     Retain mitochondrial chromosome (optional).
  -tg, --tg       Retain SP_II_TG chromosome (optional, sp only).
  -mr, --mtr      Retain SP_MTR chromosome (optional, sp only).
  -cc, --chk_chr  Check chromosomes in filtered BAM outfile (optional)
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (required; default: \${dir_out}/err_out).
  -nj, --nam_job  The name of the job, which is used when writing stderr and
                  stdout (required; default: ${nam_job}).
  -sl, --slurm    Submit jobs to the SLURM scheduler; otherwise, run them in
                  serial (optional).
  -mj, --max_job  The maximum number of jobs to run at one time (required if
                  --slurm is specified, ignored if not; default: ${max_job}).
  -tm, --time     The length of time, in 'h:mm:ss' format, for the SLURM job
                  (required if --slurm is specified, ignored if not; default:
                  ${time}).

Dependencies:
  - Programs
    + awk
    + Bash or Zsh
    + grep
    + mv
    + rm
    + Samtools
    + SLURM (if --slurm is specified)
  - Functions
    + check_supplied_arg
    + check_exists_file_dir
    + check_format_time
    + check_int_pos
    + check_program_path
    + echo_error
    + echo_warning
    + exit_0
    + exit_1
    + filter_bam_sc
    + filter_bam_sp

Note:
  #TODO

Example:
  \`\`\`
  #TODO
  \`\`\`
EOM
)

#  Parse arguments
if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    exit_0
fi

if ${interactive}; then
    #  Hardcoded paths, etc.
    dir_bas="${HOME}/tsukiyamalab/Kris"   # ls -lhaFG "${dir_bas}"
    nam_rep="202X_protocol_ChIP"          # echo "${nam_rep}"
    dir_rep="${dir_bas}/${nam_rep}"       # ls -lhaFG "${dir_rep}"
    nam_bam="03_bam"                      # echo "${nam_bam}"
    dir_bgn="${dir_rep}/${nam_bam}"       # ls -lhaFG "${dir_bgn}"
    {
        aligner="bowtie2"                                          # echo "${aligner}"
        a_type="global"                                            # echo "${a_type}"
        req_flg=true                                               # echo "${req_flg}"
        flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"  # echo "${flg}"
        mapq=1                                                     # echo "${mapq}"
        details="${aligner}/${a_type}/flag-${flg}_mapq-${mapq}"    # echo "${details}"
    }
    dir_bam="${dir_bgn}/${details}/init"  # ls -lhaFG "${dir_bam}"

    #  Hardcoded argument assignments
    verbose=true                               # echo "${verbose}"
    dry_run=true                               # echo "${dry_run}"
    retain="sp"  # "sc"                        # echo "${retain}"
    threads=4                                  # echo "${threads}"
    # infiles="$(
    #     bash "${dir_scr}/find_files.sh" \
    #         --dir_fnd "${dir_bam}" \
    #         --pattern "*.bam" \
    #         --depth 1 \
    #         --include "*G1-1st*" \
    #         --exclude "bak.*"
    # )"                                       # echo "${infiles}"
    infiles="${dir_bam}/IP_log_Brn1_rep1.bam"  # echo "${infiles}"
    dir_out="${dir_bgn}/${details}/${retain}"  # ls -lhaFG "${dir_out}"  # mkdir -p "${dir_out}/err_out"
    mito=false                                 # echo "${mito}"
    tg=false                                   # echo "${tg}"
    mtr=false                                  # echo "${mtr}"
    chk_chr=false                              # echo "${chk_chr}"
    err_out="${dir_out}/err_out"               # ls -lhaFG "${err_out}"
    slurm=true                                 # echo "${slurm}"
    max_job=12                                 # echo "${max_job}"
    time="0:30:00"                             # echo "${time}"
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -v|--verbose) verbose=true;    shift 1 ;;
            -dr|--dry_run) dry_run=true;    shift 1 ;;
             -t|--threads) threads="${2}";  shift 2 ;;
             -i|--infiles) infiles="${2}";  shift 2 ;;
            -do|--dir_out) dir_out="${2}";  shift 2 ;;
             -r|--retain)  retain="${2,,}"; shift 2 ;;
             -m|--mito)    mito=true;       shift 1 ;;
            -tg|--tg)      tg=true;         shift 1 ;;
            -mr|--mtr)     mtr=true;        shift 1 ;;
            -cc|--chk_chr) chk_chr=true;    shift 1 ;;
            -eo|--err_out) err_out="${2}";  shift 2 ;;
            -sl|--slurm)   slurm=true;      shift 1 ;;
            -mj|--max_job) max_job="${2}";  shift 2 ;;
            -tm|--time)    time="${2}";     shift 2 ;;
            *)
                echo "## Unknown parameter passed: ${1} ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                exit_1
                ;;
        esac
    done
fi

#  Check arguments
check_supplied_arg -a "${scr_fnc}" -n "scr_fnc"
check_exists_file_dir "f" "${scr_fnc}" "scr_fnc"

check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

check_supplied_arg -a "${infiles}" -n "infiles"
check_exists_file_dir "d" "$(dirname "${infiles%%,*}")" "infiles"

check_supplied_arg -a "${dir_out}" -n "dir_out"
check_exists_file_dir "d" "${dir_out}" "dir_out"

case "${retain}" in
    sc|sp) : ;;
    *)
        echo_error \
            "Selection associated with --retain is not valid: ${retain}." \
            "Selection must be 'sc' or 'sp'."
        exit_1
        ;;
esac

if [[ -n "${err_out}" ]]; then
    check_exists_file_dir "d" "${err_out}" "err_out"
elif [[ -z "${err_out}" ]]; then
    err_out="${dir_out}/err_out"
    check_exists_file_dir "d" "${err_out}" "err_out"
fi

check_supplied_arg -a "${nam_job}" -n "nam_job"

if ${slurm}; then
    check_supplied_arg -a "${scr_sub}" -n "scr_sub"
    check_exists_file_dir "f" "${scr_sub}" "scr_sub"

    check_supplied_arg -a "${max_job}" -n "max_job"
    check_int_pos "${max_job}" "max_job"
    
    check_supplied_arg -a "${time}" -n "time"
    check_format_time "${time}"
fi

#  Check that dependencies are in PATH
check_program_path awk
check_program_path grep
check_program_path mv
check_program_path rm
check_program_path samtools


#  Parse the --infiles argument -----------------------------------------------
#+ ...into an array, then validate the infile value assignments
IFS=',' read -r -a arr_infiles <<< "${infiles}"  # unset arr_infiles
# for infile in "${arr_infiles[@]}"; do echo "${infile}"; done  # unset infile

#  Check that each infile exists; if not, exit
for infile in "${arr_infiles[@]}"; do
    check_exists_file_dir "f" "${infile}"
done

if ${slurm}; then
    if [[ "${max_job}" -gt "${#arr_infiles[@]}" ]]; then
        echo_warning \
            "The maximum number of SLURM jobs to run at a time, ${max_job}," \
            "is greater than the number of infiles, ${#arr_infiles[@]}." \
            "Adjusting max_job to ${#arr_infiles[@]}."
        max_job="${#arr_infiles[@]}"
    fi
fi


#  Do the main work ===========================================================
#  Report argument variable assignments if in "verbose mode"
if ${verbose}; then
cat << EOM
#TODO
env_nam=${env_nam}
nam_job=${nam_job}
scr_fnc=${scr_fnc}
scr_sub=${scr_sub}

###################################
## Argument variable assignments ##
###################################

verbose=${verbose}
dry_run=${dry_run}
threads=${threads}
infiles=${infiles}
dir_out=${dir_out}
retain=${retain}
mito=${mito}
tg=${tg}
mtr=${mtr}
chk_chr=${chk_chr}
err_out=${err_out}
slurm=${slurm}
max_job=${max_job}
time=${time}


EOM
fi

# shellcheck disable=SC1083,SC2157
if ${slurm}; then
    #  If --slurm was specified, run jobs in parallel via individual job
    #+ submissions to SLURM
    if ${dry_run} || ${verbose}; then
        echo "####################"
        echo "## Call to sbatch ##"
        echo "####################"
        echo ""
        echo "sbatch \\"
        echo "    --job-name=${nam_job} \\"
        echo "    --nodes=1 \\"
        echo "    --cpus-per-task=${threads} \\"
        echo "    --time=${time} \\"
        echo "    --error=${err_out}/${nam_job}.%A-%a.stderr.txt \\"
        echo "    --output=${err_out}/${nam_job}.%A-%a.stdout.txt \\"
        echo "    --array=1-${#arr_infiles[@]}%${max_job} \\"
        echo "    ${scr_sub} \\"
        echo "        --scr_fnc ${scr_fnc} \\"
        echo "        --threads ${threads} \\"
        echo "        --infiles ${infiles} \\"
        echo "        --dir_out ${dir_out} \\"
        echo "        $(if ${mito}; then echo "--mito"; fi) \\"
        echo "        $(if ${tg}; then echo "--tg"; fi) \\"
        echo "        $(if ${mtr}; then echo "--mtr"; fi) \\"
        echo "        $(if ${chk_chr}; then echo "--chk_chr"; fi) \\"
        echo "        --err_out ${err_out} \\"
        echo "        --nam_job ${nam_job}"
        echo ""
        echo ""
        echo "#########################################"
        echo "## Contents of SLURM submission script ##"
        echo "#########################################"
        echo ""
        echo "## ${scr_sub} ##"
        echo ""
        cat "${scr_sub}"
        echo ""
    fi

    if ! ${dry_run}; then
        # shellcheck disable=SC2046,SC2086
        sbatch \
            --job-name=${nam_job} \
            --nodes=1 \
            --cpus-per-task=${threads} \
            --time=${time} \
            --error=${err_out}/${nam_job}.%A-%a.stderr.txt \
            --output=${err_out}/${nam_job}.%A-%a.stdout.txt \
            --array=1-${#arr_infiles[@]}%${max_job} \
            ${scr_sub} \
                --scr_fnc ${scr_fnc} \
                --threads ${threads} \
                --infiles ${infiles} \
                --dir_out ${dir_out} \
                $(if ${mito}; then echo "--mito"; fi) \
                $(if ${tg}; then echo "--tg"; fi) \
                $(if ${mtr}; then echo "--mtr"; fi) \
                $(if ${chk_chr}; then echo "--chk_chr"; fi) \
                --err_out ${err_out} \
                --nam_job ${nam_job}
    fi
else
    #  If --slurm was not specified, run jobs in serial
    if [[ "${CONDA_DEFAULT_ENV}" != "env_align" ]]; then
        eval "$(conda shell.bash hook)"
        conda activate env_align
    fi
    
    #  Derive function name from function script
    nam_fnc="$(basename "${scr_fnc}" ".sh")"

    if ${dry_run} || ${verbose}; then
        echo "###################################################"
        echo "## Serial call(s) for species-specific filtering ##"
        echo "###################################################"
        echo ""
    fi

    for idx in "${!arr_infiles[@]}"; do
        #  Assign infile via array index
        n_call=$(( idx + 1 ))
        infile="${arr_infiles[idx]}"
        
        if [[ -z "${infile}" ]]; then
            echo_error \
                "Failed to retrieve infile for idx=${idx}:" \
                "\${arr_infiles[idx]}."
            exit_1
        fi

        #  Derive sample name from infile assignment
        samp="${infile##*/}"
        samp="${samp%.bam}"

        #  Perform pattern matching to assign the outfile name
        # shellcheck disable=SC2086
        if [[ ${nam_fnc} =~ "sc" ]]; then
            outfile="${dir_out}/${samp}.sc.bam"
        elif [[ ${nam_fnc} =~ "sp" ]]; then
            outfile="${dir_out}/${samp}.sp.bam"
        else
            echo "Error: Sample name could not be processed." >&2
            exit 1
        fi

        if ${dry_run} || ${verbose}; then
            echo "## Call no. ${n_call} ##"
            echo "${nam_fnc} \\"
            echo "    --threads ${threads} \\"
            echo "    --infile ${infile} \\"
            echo "    --outfile ${outfile} \\"
            echo "    $(if ${mito}; then echo "--mito"; fi) \\"
            echo "    $(if ${tg}; then echo "--tg"; fi) \\"
            echo "    $(if ${mtr}; then echo "--mtr"; fi) \\"
            echo "    $(if ${chk_chr}; then echo "--chk_chr"; fi) \\"
            echo "         > ${err_out}/${nam_job}.${samp}.stdout.txt \\"
            echo "        2> ${err_out}/${nam_job}.${samp}.stderr.txt"
            echo ""
        fi

        if ! ${dry_run}; then
            # shellcheck disable=SC2046,SC2086
            ${nam_fnc} \
                --threads ${threads} \
                --infile ${infile} \
                --outfile ${outfile} \
                $(if ${mito}; then echo "--mito"; fi) \
                $(if ${tg}; then echo "--tg"; fi) \
                $(if ${mtr}; then echo "--mtr"; fi) \
                $(if ${chk_chr}; then echo "--chk_chr"; fi) \
                     > ${err_out}/${nam_job}.${samp}.stdout.txt \
                    2> ${err_out}/${nam_job}.${samp}.stderr.txt
        fi
    done
fi
