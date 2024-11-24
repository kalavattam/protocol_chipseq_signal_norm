#!/bin/bash

#  execute_align_fastqs.sh
#  KA


#  Run script in interactive/test mode (true) or command-line mode (false)
interactive=false

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive}; then set -euo pipefail; fi

#  Set the path to the "scripts" directory
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
    source "${dir_fnc}/align_fastqs.sh"
    source "${dir_fnc}/check_exists_file_dir.sh"
    source "${dir_fnc}/check_format_time.sh"
    source "${dir_fnc}/check_int_nonneg.sh"
    source "${dir_fnc}/check_int_pos.sh"
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/check_supplied_arg.sh"
    source "${dir_fnc}/echo_error.sh"
    source "${dir_fnc}/echo_warning.sh"
    source "${dir_fnc}/exit_0.sh"
    source "${dir_fnc}/exit_1.sh"
    source "${dir_fnc}/handle_env.sh"
}


#  Set up paths, values, and parameters for interactive mode
function set_interactive() {
    #  Set hardcoded paths, values, etc.
    ## WARNING: Change values if you're not Kris and `interactive=true` ##
    dir_bas="${HOME}/tsukiyamalab/Kris"
    dir_rep="${dir_bas}/202X_protocol_ChIP"
    dir_dat="${dir_rep}/data"
    dir_idx="${dir_dat}/genomes/concat/index"
    str_idx="sc_sp_proc"
    dir_pro="${dir_dat}/processed"
    dir_trm="${dir_pro}/trim_atria"
    dir_aln="${dir_pro}/align"

    #  Set hardcoded argument assignments, etc.
    verbose=true
    dry_run=true
    threads=8
    aligner="bowtie2"
    a_type="global"
    mapq=1
    req_flg=true
    {
        flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
    }
    index="${dir_idx}/${str_idx}"
    if [[ ${aligner} == "bwa" ]]; then index="${index}.fa"; fi
    infiles="$(  ## WARNING: Change the search parameters as needed ##
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_trm}" \
            --pattern "*.atria.fastq.gz" \
            --depth 1 \
            --fastqs
    )"
    {
        dir_aln="${dir_aln}_${aligner}_${a_type}"
    }
    dir_out="${dir_aln}/flag-${flg}_mapq-${mapq}/init"
    qname=false  # true
    sfx_se=".atria.fastq.gz"
    sfx_pe="_R1.atria.fastq.gz"
    err_out="${dir_out}/logs"
    nam_job="align_fastqs"
    slurm=true  # false
    max_job=6
    time="1:00:00"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_align"
scr_sub="${dir_scr}/submit_align_fastqs.sh"
scr_fnc="${dir_fnc}/align_fastqs.sh"

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=1
aligner="bowtie2"
a_type="global"
mapq=1
req_flg=false
index=""
infiles=""
dir_out=""
qname=false
sfx_se=".atria.fastq.gz"
sfx_pe="_R1.atria.fastq.gz"
err_out=""
nam_job="align_fastqs"
slurm=false
max_job=6
time="1:00:00"

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  execute_align_fastqs.sh
    [--verbose] [--dry_run] --threads <int> --aligner <str> --a_type <str>
    --mapq <int> [--req_flg] --index <str> --infiles <str> --dir_out <str>
    [--qname] --sfx_se <str> --sfx_pe <str> --err_out <str> --nam_job <str>
    [--slurm] [--max_job <int>] [--time <str>]

Description:
  execute_align_fastqs.sh aligns single- or paired-end FASTQ files using Bowtie
  2 or BWA, followed by processing steps with Samtools, including filtering,
  sorting, mate fixing (for paired-end reads), duplicate marking, and indexing.
  The script supports both parallel execution via SLURM or serial execution.
  The '--qname' flag optionally retains an intermediate queryname-sorted BAM
  file.

Arguments:
   -h, --help     Display this help message and exit (0).
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Perform a dry run without executing commands (optional).
   -t, --threads  Number of threads to use (default: ${threads}).
   -a, --aligner  Alignment program to use ('bowtie2' or 'bwa'; default:
                  ${aligner}).
  -at, --a_type   Alignment type (required if using Bowtie 2, ignored if using
                  BWA: 'local', 'global', or 'end-to-end'; default:
                  ${a_type}).
  -mq, --mapq     MAPQ threshold for filtering BAM outfiles (default: ${mapq}).
                  To perform no filtering, specify 0.
  -rf, --req_flg  Require flag bit 2, signifying that paired-end alignments are
                  properly paired, for filtering BAM outfiles (optional;
                  ignored if working with single-end sequenced reads).
  -ix, --index    Path to the directory containing the aligner index. If using
                  Bowtie 2, the path should end with the index stem:
                  "\${HOME}/path/stem"; if using BWA, the path should end with
                  the stem and a FASTA file extension (.fa):
                  "\${HOME}/path/stem.fa".
   -i, --infiles  Semicolon-separated string vector of FASTQ infiles. Within
                  the semicolon delimiters, sample-specific paired-end files
                  must be together separated by commas, e.g.,
                  "\${HOME}/path/samp_1.fastq.gz;\${HOME}/path/samp_2_R1.fastq.gz,\${HOME}/path/samp_2_R2.fastq.gz;\${HOME}/path/samp_3.fastq.gz".
  -do, --dir_out  Directory to write BAM outfiles.
  -qn, --qname    Retain queryname-sorted intermediate BAM files (optional).
  -ss, --sfx_se   Suffix to strip from single-end sequenced FASTQ files
                  (default: ${sfx_se}).
  -sp, --sfx_pe   Suffix to strip from paired-end sequenced FASTQ files
                  (default: ${sfx_pe}).
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (default: \${dir_out}/err_out).
  -nj, --nam_job  The name of the job, which is used when writing stderr and
                  stdout TXT files (default: ${nam_job}).
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
    + mv
    + rm
    + Samtools
    + SLURM (if --slurm is specified)
  - Functions
    + align_fastqs
    + check_exists_file_dir
    + check_format_time
    + check_int_nonneg
    + check_int_pos
    + check_program_path
    + check_supplied_arg
    + echo_error
    + echo_warning
    + exit_0
    + exit_1
    + handle_env
    + handle_env_activate
    + handle_env_deactivate

Notes:
  - When the '--slurm' flag is used, jobs are parallelized via SLURM array
    tasks; otherwise, jobs run serially.
  - Warning: Running alignment and post-alignment processing serially with many
    files is not recommended, as it can be extremely time-consuming.
  - If using Bowtie 2, ensure the path to index files ends with the index stem,
    such as "\${HOME}/path/stem" or "\${HOME}/path/sc_sp_proc". For BWA, the
    path should include the stem and a FASTA file extension, e.g.,
    "\${HOME}/path/stem.fa" or "\${HOME}/path/sc_sp_proc.fa".
  - Queryname-sorted BAM files will share the same path and stem as specified
    for the '--outfile' argument, but with the extension '.qnam.bam' instead of
    '.bam'.

Example:
  \`\`\`
  bash "\${dir_scr}/execute_align_fastqs.sh"
      --verbose
      --threads "\${threads}"
      --aligner "\${aligner}"
      --a_type "\${a_type}"
      --mapq "\${mapq}"
      --req_flg
      --index "\${pth_idx}"
      --infiles "\${infiles}"
      --dir_out "\${dir_out}/init"
      --err_out "\${dir_out}/init/logs"
      --slurm
  \`\`\`
EOM
)

#  Parse arguments
if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    exit_0
fi

if ${interactive}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -v|--verbose) verbose=true;     shift 1 ;;
            -dr|--dry_run) dry_run=true;     shift 1 ;;
             -t|--threads) threads="${2}";   shift 2 ;;
             -a|--aligner) aligner="${2,,}"; shift 2 ;;
            -at|--a_type)  a_type="${2,,}";  shift 2 ;;
            -mq|--mapq)    mapq="${2}";      shift 2 ;;
            -rf|--req_flg) req_flg=true;     shift 1 ;;
            -ix|--index)   index="${2}";     shift 2 ;;
             -i|--infiles) infiles="${2}";   shift 2 ;;
            -do|--dir_out) dir_out="${2}";   shift 2 ;;
            -qn|--qname)   qname=true;       shift 1 ;;
            -ss|--sfx_se)  sfx_se="${2}";    shift 2 ;;
            -sp|--sfx_pe)  sfx_pe="${2}";    shift 2 ;;
            -eo|--err_out) err_out="${2}";   shift 2 ;;
            -nj|--nam_job) nam_job="${2}";   shift 2 ;;
            -sl|--slurm)   slurm=true;       shift 1 ;;
            -mj|--max_job) max_job="${2}";   shift 2 ;;
            -tm|--time)    time="${2}";      shift 2 ;;
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
check_supplied_arg -a "${env_nam}" -n "env_nam"

check_supplied_arg -a "${scr_sub}" -n "scr_sub"
check_exists_file_dir "f" "${scr_sub}" "scr_sub"

check_supplied_arg -a "${scr_fnc}" -n "scr_fnc"
check_exists_file_dir "f" "${scr_fnc}" "scr_fnc"

check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

case "${aligner}" in
    bowtie2)
        case "${a_type}" in
            local|global|end-to-end) : ;;
            *)
                echo_error \
                    "Selection associated with --a_type is not valid:" \
                    "${a_type}. Selection must be 'local', 'global', or" \
                    "'end-to-end'."
                exit_1
                ;;
        esac
        ;;
    bwa) unset a_type ;;
    *)
        echo_error \
            "Selection associated with --aligner is not valid: ${aligner}." \
            "Selection must be 'bowtie2' or 'bwa'."
        exit_1
        ;;
esac

check_supplied_arg -a "${mapq}" -n "mapq"
check_int_nonneg "${mapq}" "mapq"

check_supplied_arg -a "${index}" -n "index"
check_exists_file_dir "d" "$(dirname "${index}")" "index"

check_supplied_arg -a "${infiles}" -n "infiles"
check_exists_file_dir "d" "$(dirname "${infiles%%[,;]*}")" "infiles"

check_supplied_arg -a "${dir_out}" -n "dir_out"
check_exists_file_dir "d" "${dir_out}" "dir_out"

check_supplied_arg -a "${sfx_se}" -n "sfx_se"
check_supplied_arg -a "${sfx_pe}" -n "sfx_pe"

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

    check_supplied_arg -a "${scr_fnc}" -n "scr_fnc"
    check_exists_file_dir "f" "${scr_fnc}" "scr_fnc"

    check_supplied_arg -a "${max_job}" -n "max_job"
    check_int_pos "${max_job}" "max_job"
    
    check_supplied_arg -a "${time}" -n "time"
    check_format_time "${time}"
fi

#  Activate environment and check that dependencies are in PATH
handle_env "${env_nam}" > /dev/null

case "${aligner}" in
    bowtie2) check_program_path bowtie2 ;;
    bwa)     check_program_path bwa     ;;
esac

check_program_path samtools


#  Parse the --infiles argument -----------------------------------------------
#+ ...into an array, then validate the infile value assignments
IFS=';' read -r -a arr_infiles <<< "${infiles}"  # unset arr_infiles
# for infile in "${arr_infiles[@]}"; do echo "${infile}"; done  # unset infile

#  Check that each infile exists; if not, exit
for infile in "${arr_infiles[@]}"; do
    # infile="${arr_infiles[0]}"
    
    #  Check that the infile contains a comma
    if [[ "${infile}" == *,* ]]; then
        fq_1="${infile%%,*}"
        fq_2="${infile#*,}"

        #  Ensure there is only one comma
        if [[ "${fq_2}" == *,* ]]; then
            echo_error \
                "It appears an error is present in the string output from" \
                "running find_files.sh in '--fastqs' mode. More than one" \
                "comma was found in a specific substring: '${infile}'." \
                "Exiting now."
            exit_1
        fi

        #  Validate sample-specific paired-end FASTQ files exist
        check_exists_file_dir "f" "${fq_1}" "fq_1"
        check_exists_file_dir "f" "${fq_2}" "fq_2"

        #  Validate presence of file suffix in fq_1 assignment
        if [[ "${fq_1}" != *"${sfx_pe}" ]]; then
            echo_error \
                "Suffix '${sfx_pe}' not found in file name '${fq_1}'. Check" \
                "--sfx_pe assignment."
            exit_1
        fi
    else
        fq_1="${infile}"
        unset fq_2
        
        #  Validate sample-specific single-end FASTQ file exists
        check_exists_file_dir "f" "${fq_1}" "fq_1"

        #  Validate presence of file suffix in fq_1 assignment
        if [[ "${fq_1}" != *"${sfx_se}" ]]; then
            echo_error \
                "Suffix '${sfx_se}' not found in file name '${fq_1}'. Check" \
                "--sfx_se assignment."
            exit_1
        fi
    fi
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
if ${verbose}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "env_nam=${env_nam}"
    echo "scr_sub=${scr_sub}"
    echo "scr_fnc=${scr_fnc}"
    echo ""
    echo ""
    echo "###################################"
    echo "## Argument variable assignments ##"
    echo "###################################"
    echo ""
    echo "verbose=${verbose}"
    echo "dry_run=${dry_run}"
    echo "threads=${threads}"
    echo "aligner=${aligner}"
    echo "a_type=${a_type}"
    echo "mapq=${mapq}"
    echo "req_flg=${req_flg}"
    echo "index=${index}"
    echo "infiles=${infiles}"
    echo "dir_out=${dir_out}"
    echo "qname=${qname}"
    echo "sfx_se=${sfx_se}"
    echo "sfx_pe=${sfx_pe}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo "slurm=${slurm}"
    echo "max_job=${max_job}"
    echo "time=${time}"
    echo ""
    echo ""
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
        echo "        --aligner ${aligner} \\"
        echo "        $(if [[ -n ${a_type} ]]; then echo "--a_type ${a_type}"; fi) \\"
        echo "        --mapq ${mapq} \\"
        echo "        $(if ${req_flg}; then echo "--req_flg"; fi) \\"
        echo "        --index ${index} \\"
        echo "        --infiles ${infiles} \\"
        echo "        --dir_out ${dir_out} \\"
        echo "        $(if ${qname}; then echo "--qname"; fi) \\"
        echo "        --sfx_se ${sfx_se} \\"
        echo "        --sfx_pe ${sfx_pe} \\"
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
                --aligner ${aligner} \
                $(if [[ -n ${a_type} ]]; then echo "--a_type ${a_type}"; fi) \
                --mapq ${mapq} \
                $(if ${req_flg}; then echo "--req_flg"; fi) \
                --index ${index} \
                --infiles ${infiles} \
                --dir_out ${dir_out} \
                $(if ${qname}; then echo "--qname"; fi) \
                --sfx_se ${sfx_se} \
                --sfx_pe ${sfx_pe} \
                --err_out ${err_out} \
                --nam_job ${nam_job}
    fi
else
    #  If --slurm was not specified, run jobs in serial
    if ${dry_run} || ${verbose}; then
        echo "#################################################"
        echo "## Serial call(s) for alignment and processing ##"
        echo "#################################################"
        echo ""
    fi

    for idx in "${!arr_infiles[@]}"; do
        n_call=$(( idx + 1 ))
        infile="${arr_infiles[idx]}"
        
        if [[ -z "${infile}" ]]; then
            echo_error \
                "Failed to retrieve infile for idx=${idx}:" \
                "\${arr_infiles[idx]}."
            exit_1
        fi

        if [[ "${infile}" == *,* ]]; then
            fq_1="${infile%%,*}"
            fq_2="${infile#*,}"
            samp="$(basename "${fq_1%%"${sfx_pe}"}")"
        else
            fq_1="${infile}"
            unset fq_2
            samp="$(basename "${fq_1%%"${sfx_se}"}")"
        fi

        outfile="${dir_out}/${samp}.bam"

        if ${dry_run} || ${verbose}; then
            echo "## Call no. ${n_call} ##"
            echo "align_fastqs \\"
            echo "    --threads ${threads} \\"
            echo "    --aligner ${aligner} \\"
            echo "    $(if [[ -n ${a_type} ]]; then echo "--a_type ${a_type}"; fi) \\"
            echo "    $(if [[ ${mapq} -gt 0 ]]; then echo "--mapq ${mapq}"; fi) \\"
            echo "    $(if ${req_flg}; then echo "--req_flg"; fi) \\"
            echo "    --index ${index} \\"
            echo "    --fq_1 ${fq_1} \\"
            echo "    $(if [[ -n ${fq_2} ]]; then echo "--fq_2 ${fq_2}"; fi) \\"
            echo "    --outfile ${outfile} \\"
            echo "    $(if ${qname}; then echo "--qname"; fi) \\"
            echo "         > ${err_out}/${nam_job}.${samp}.stdout.txt \\"
            echo "        2> ${err_out}/${nam_job}.${samp}.stderr.txt"
            echo ""
        fi

        if ! ${dry_run}; then
            # shellcheck disable=SC2046,SC2086
            align_fastqs \
                --threads ${threads} \
                --aligner ${aligner} \
                $(if [[ -n ${a_type} ]]; then echo "--a_type ${a_type}"; fi) \
                $(if [[ ${mapq} -gt 0 ]]; then echo "--mapq ${mapq}"; fi) \
                $(if ${req_flg}; then echo "--req_flg"; fi) \
                --index ${index} \
                --fq_1 ${fq_1} \
                $(if [[ -n ${fq_2} ]]; then echo "--fq_2 ${fq_2}"; fi) \
                --outfile ${outfile} \
                $(if ${qname}; then echo "--qname"; fi) \
                     > ${err_out}/${nam_job}.${samp}.stdout.txt \
                    2> ${err_out}/${nam_job}.${samp}.stderr.txt
        fi
    done
fi
