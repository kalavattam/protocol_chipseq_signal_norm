#!/bin/bash

#  execute_align_fastqs.sh
#  KA


#  Run script in interactive mode (true) or command-line mode (false)
interactive=false

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive}; then set -euo pipefail; fi

#  Set the path to the "scripts" directory
if ${interactive:-false}; then
    ## WARNING: If 'interactive=true', change path as needed ##
    dir_scr="${HOME}/repos/protocol_chipseq_signal_norm/scripts"
else
    dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi


#  Source and define functions ================================================
#  Set the path to the "functions" directory
dir_fnc="${dir_scr}/functions"

# shellcheck disable=SC1090
for fnc in \
    align_fastqs \
    check_exists_file_dir \
    check_format_time \
    check_int_nonneg \
    check_int_pos \
    check_program_path \
    check_string_fastqs \
    check_supplied_arg \
    echo_error \
    echo_warning \
    exit_0 \
    exit_1 \
    handle_env \
    print_parallel_info \
    reset_max_job \
    set_params_parallel
do
    source "${dir_fnc}/${fnc}.sh"
done
unset fnc


#  Set paths, values, arguments, etc. for interactive mode
function set_interactive() {
    #  Set base repository paths
    dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
    dir_rep="${dir_bas}/protocol_chipseq_signal_norm"

    #  Define data directories
    dir_dat="${dir_rep}/data"
    dir_idx="${dir_dat}/genomes/concat/index"
    dir_pro="${dir_dat}/processed"
    dir_trm="${dir_pro}/trim_fastqs"
    dir_aln="${dir_pro}/align_fastqs"
    
    #  Set hardcoded argument assignments and related variables
    verbose=true
    dry_run=true
    threads=8
    aligner=bowtie2
    a_type="global"
    mapq=1
    req_flg=true
    flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
    str_idx="sc_sp_proc"
    if [[ ${aligner} != "bwa" ]]; then
        index="${dir_idx}/${aligner}/${str_idx}"
    else
        index="${dir_idx}/${aligner}/${str_idx}.fa"
    fi
    infiles="$(  ## WARNING: Change search parameters as needed ##
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_trm}" \
            --pattern "*.atria.fastq.gz" \
            --depth 1 \
            --fastqs
    )"
    str_det="${aligner}_${a_type}_flag-${flg}_mapq-${mapq}"
    dir_out="${dir_aln}/${str_det}/init"
    qname=false
    sfx_se=".atria.fastq.gz"
    sfx_pe="_R1.atria.fastq.gz"
    err_out="${dir_out}/logs"
    nam_job="align_fastqs"
    slurm=false
    max_job=2
    time="1:00:00"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_protocol"
scr_sub="${dir_scr}/submit_align_fastqs.sh"
scr_fnc="${dir_fnc}/align_fastqs.sh"
par_job=""

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
max_job=6
slurm=false
time="1:00:00"

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  execute_align_fastqs.sh
    [--verbose] [--dry_run] --threads <int> --aligner <str> [--a_type <str>]
    --mapq <int> [--req_flg] --index <str> --infiles <str> --dir_out <str>
    [--qname] --sfx_se <str> --sfx_pe <str> --err_out <str> --nam_job <str>
    --max_job <int> [--slurm] [--time <str>]

Description:
  The driver script 'execute_align_fastqs.sh' aligns single- or paired-end
  FASTQ files using Bowtie 2 or BWA, followed by processing steps with
  Samtools, including filtering, sorting, mate fixing (for paired-end reads),
  duplicate marking, and indexing.

  The script supports both parallel execution via SLURM ('--slurm') or GNU
  Parallel, or serial execution. 

Arguments:
   -h, --help     Display this help message and exit (0).
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Perform a dry run without executing commands (optional).
   -t, --threads  Number of threads to use (default: ${threads}).
   -a, --aligner  Alignment program to use ('bowtie2' or 'bwa'; default:
                  '${aligner}').
  -at, --a_type   Alignment type (required if using Bowtie 2, ignored if using
                  BWA: 'local', 'global', or 'end-to-end'; default:
                  '${a_type}').
  -mq, --mapq     MAPQ threshold for filtering BAM outfiles (default: '${mapq}').
                  To perform no filtering, specify 0.
  -rf, --req_flg  Require flag bit 2, signifying that paired-end alignments are
                  properly paired, for filtering BAM outfiles (optional;
                  ignored if working with single-end sequenced reads).
  -ix, --index    Path to the directory containing the aligner index. If using
                  Bowtie 2, the path should end with the index stem:
                  "\${HOME}/path/stem"; if using BWA, the stem should end with
                  a FASTA file extension, e.g., "\${HOME}/path/stem.fa".
   -i, --infiles  Semicolon-separated string vector of FASTQ infiles. Within
                  the semicolon delimiters, sample-specific paired-end files
                  must be together separated by commas, e.g.,
                  "\${HOME}/path/samp_1.fastq.gz;\${HOME}/path/samp_2_R1.fastq.gz,\${HOME}/path/samp_2_R2.fastq.gz;\${HOME}/path/samp_3.fastq.gz".
  -do, --dir_out  Directory to write BAM outfiles.
  -qn, --qname    Retain queryname-sorted intermediate BAM files (optional).
  -ss, --sfx_se   Suffix to strip from single-end sequenced FASTQ files
                  (default: '${sfx_se}').
  -sp, --sfx_pe   Suffix to strip from paired-end sequenced FASTQ files
                  (default: '${sfx_pe}').
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (default: \${dir_out}/err_out).
  -nj, --nam_job  The name of the job, which is used when writing stderr and
                  stdout TXT files (default: '${nam_job}').
  -mj, --max_job  Maximum number of jobs to run concurrently (default: '${max_job}').
                    - If '--slurm' is specified, controls SLURM array tasks.
                    - If '--slurm' is not specified:
                      + If 'max_job' is greater than 1, jobs run in parallel
                        via GNU Parallel.
                      + If 'max_job' is 1, jobs run sequentially (serial mode).
  -sl, --slurm    Submit jobs to the SLURM scheduler; otherwise, run them in
                  serial (optional).
  -tm, --time     The length of time, in 'h:mm:ss' format, for the SLURM job
                  (required if '--slurm' is specified, ignored if not; default:
                  '${time}').

Dependencies:
  - AWK
  - Bash or Zsh
  - GNU Parallel (when '--slurm' is not specified but multiple jobs are)
  - mv
  - rm
  - Samtools
  - SLURM (when '--slurm' is specified)

Notes:
  - When the '--slurm' flag is used, jobs are parallelized via SLURM array
    tasks; otherwise... #TODO
  - Warning: Running alignment and post-alignment processing serially (i.e., in
    a sequence one after the other) with many files is not recommended, as it
    can be extremely time-consuming.
  - If using Bowtie 2, ensure the path to index files ends with the index stem,
    e.g., "\${HOME}/path/stem" or "\${HOME}/path/sc_sp_proc". For BWA, the path
    should include the stem and a FASTA file extension, e.g.,
    "\${HOME}/path/stem.fa" or "\${HOME}/path/sc_sp_proc.fa".
  - Calling the script with the optional '--qname' flag retains an intermediate
    queryname-sorted BAM file used for mate fixing.
  - Retained queryname-sorted BAM files will share the same path and stem as
    specified for the '--outfile' argument, but with the extension ".qnam.bam"
    in place of ".bam".

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

if ${interactive:-false}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -v|--verbose) verbose=true;   shift 1 ;;
            -dr|--dry_run) dry_run=true;   shift 1 ;;
             -t|--threads) threads="${2}"; shift 2 ;;
             -a|--aligner)
                aligner="$(echo "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;
            -at|--a_type)
                a_type="$(echo "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;
            -mq|--mapq)    mapq="${2}";    shift 2 ;;
            -rf|--req_flg) req_flg=true;   shift 1 ;;
            -ix|--index)   index="${2}";   shift 2 ;;
             -i|--infiles) infiles="${2}"; shift 2 ;;
            -do|--dir_out) dir_out="${2}"; shift 2 ;;
            -qn|--qname)   qname=true;     shift 1 ;;
            -ss|--sfx_se)  sfx_se="${2}";  shift 2 ;;
            -sp|--sfx_pe)  sfx_pe="${2}";  shift 2 ;;
            -eo|--err_out) err_out="${2}"; shift 2 ;;
            -nj|--nam_job) nam_job="${2}"; shift 2 ;;
            -mj|--max_job) max_job="${2}"; shift 2 ;;
            -sl|--slurm)   slurm=true;     shift 1 ;;
            -tm|--time)    time="${2}";    shift 2 ;;
            *)
                echo "## Unknown parameter passed: '${1}' ##" >&2
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
                    "Selection associated with '--a_type' is not valid:" \
                    "'${a_type}'. Selection must be 'local', 'global', or" \
                    "'end-to-end'."
                exit_1
                ;;
        esac
        ;;
    bwa) a_type="#N/A" ;;
    *)
        echo_error \
            "Selection associated with '--aligner' is not valid: "\
            "'${aligner}'. Selection must be 'bowtie2' or 'bwa'."
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


#  Parse and validate infiles -------------------------------------------------
IFS=';' read -r -a arr_infile <<< "${infiles}" && unset IFS

check_string_fastqs "${infiles}"derive "${sfx_se}" "${sfx_pe}" || exit_1


#  Parse job execution parameters ---------------------------------------------
if ${slurm:-false}; then
    check_supplied_arg -a "${max_job}" -n "max_job"
    check_int_pos "${max_job}" "max_job"
    
    max_job=$(reset_max_job "${max_job}" "${#arr_infile[@]}")
    
    check_supplied_arg -a "${time}" -n "time"
    check_format_time "${time}"
else
    IFS=';' read -r threads par_job < <(
        set_params_parallel "${threads}" "${max_job}" "${par_job}"
    )
    unset IFS max_job time

    check_supplied_arg -a "${par_job}" -n "par_job"
    check_int_pos "${par_job}" "par_job"
fi

#  Debug parallelization information
if ${verbose}; then
    print_parallel_info \
        "${slurm}" "${max_job:-#N/A}" "${par_job}" "${threads}" "arr_infile"
fi


#  Activate environment and check that dependencies are in PATH ---------------
handle_env "${env_nam}" > /dev/null

case "${aligner}" in
    bowtie2) check_program_path bowtie2 ;;
    bwa)     check_program_path bwa     ;;
esac

check_program_path samtools

if ${slurm:-false}; then
    check_program_path sbatch
elif [[ ${par_job} -gt 1 ]]; then
    check_program_path parallel
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
    echo "par_job=${par_job:-#N/A}"
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
    echo "max_job=${max_job:-#N/A}"
    echo "slurm=${slurm}"
    echo "time=${time:-#N/A}"
    echo ""
    echo ""
    echo "#################################"
    echo "## Array derived from variable ##"
    echo "#################################"
    echo ""
    echo "arr_infile=( ${arr_infile[*]} )"
    echo ""
    echo ""
fi

# shellcheck disable=SC1083,SC2157
if ${slurm:-false}; then
    #  If --slurm was specified, run jobs in parallel via SLURM
    if ${dry_run} || ${verbose}; then
        echo "######################"
        echo "## Call to 'sbatch' ##"
        echo "######################"
        echo ""
        echo "sbatch \\"
        echo "    --job-name=${nam_job} \\"
        echo "    --nodes=1 \\"
        echo "    --cpus-per-task=${threads} \\"
        echo "    --time=${time} \\"
        echo "    --error=${err_out}/${nam_job}.%A-%a.stderr.txt \\"
        echo "    --output=${err_out}/${nam_job}.%A-%a.stdout.txt \\"
        echo "    --array=1-${#arr_infile[@]}%${max_job} \\"
        echo "    ${scr_sub} \\"
        echo "        -en ${env_nam} \\"
        echo "        -sf ${scr_fnc} \\"
        echo "         -t ${threads} \\"
        echo "         -a ${aligner} \\"
        echo "        -at ${a_type} \\"
        echo "        -mq ${mapq} \\"
        echo "        $(if ${req_flg}; then echo "-rf"; fi) \\"
        echo "        -ix ${index} \\"
        echo "         -i ${infiles} \\"
        echo "        -do ${dir_out} \\"
        echo "        $(if ${qname}; then echo "-qn"; fi) \\"
        echo "        -ss ${sfx_se} \\"
        echo "        -sp ${sfx_pe} \\"
        echo "        -eo ${err_out} \\"
        echo "        -nj ${nam_job}"
        echo ""
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
            --array=1-${#arr_infile[@]}%${max_job} \
            ${scr_sub} \
                -en ${env_nam} \
                -sf ${scr_fnc} \
                 -t ${threads} \
                 -a ${aligner} \
                -at ${a_type} \
                -mq ${mapq} \
                $(if ${req_flg}; then echo "-rf"; fi) \
                -ix ${index} \
                 -i ${infiles} \
                -do ${dir_out} \
                $(if ${qname}; then echo "-qn"; fi) \
                -ss ${sfx_se} \
                -sp ${sfx_pe} \
                -eo ${err_out} \
                -nj ${nam_job}
    fi
else
    #  If '--slurm' was not specified, run jobs in parallel with GNU Parallel
    if [[ "${par_job}" -gt 1 ]]; then
        config="${err_out}/${nam_job}.config_parallel.txt"

        if [[ -f "${config}" ]]; then rm "${config}"; fi
        touch "${config}" || {
            echo_error "Failed to create a GNU Parallel configuration file."
            exit_1
        }

        #  Populate GNU Parallel configuration file
        for idx in "${!arr_infile[@]}"; do
            echo \
                "${env_nam}" \
                "${scr_fnc}" \
                "${threads}" \
                "${aligner}" \
                "${a_type}" \
                "${mapq}" \
                "${index}" \
                "${arr_infile[idx]}" \
                "${dir_out}" \
                "${sfx_se}" \
                "${sfx_pe}" \
                "${err_out}" \
                "${nam_job}" \
                    >> "${config}"
        done

        #  Construct command to be passed to GNU Parallel
        cmd="bash ${scr_sub}"
        cmd+=" -en {1}"   # Environment name
        cmd+=" -sf {2}"   # Function script
        cmd+="  -t {3}"   # Number of threads
        cmd+="  -a {4}"   # Aligner
        cmd+=" -at {5}"   # Alignment type
        cmd+=" -mq {6}"   # MAPQ threshold (for filtering BAM files)
        
        if ${req_flg}; then 
            cmd+=" -rf"   # Flag: Require flag bit 2
        fi  

        cmd+=" -ix {7}"   # Index path
        cmd+="  -i {8}"   # Input FASTQ file(s)
        cmd+=" -do {9}"   # Output BAM directory
        
        if ${qname}; then 
            cmd+=" -qn"   # Flag: Retain queryname-sorted BAM files
        fi  

        cmd+=" -ss {10}"  # SE FASTQ file suffix
        cmd+=" -sp {11}"  # PE FASTQ file suffix
        cmd+=" -eo {12}"  # Directory for stderr/stdout logs
        cmd+=" -nj {13}"  # Job name
        cmd+="  > {12}/{13}_par.{8/.}.stdout.txt"  # 'scr_sub' stdout log
        cmd+=" 2> {12}/{13}_par.{8/.}.stderr.txt"  # 'scr_sub' stderr log
        #TODO: Not enough is stripped from {8}; may want to pass 'samp' instead

        if ${dry_run} || ${verbose}; then
            echo "####################################################"
            echo "## Parallelized job execution via 'submit' script ##"
            echo "####################################################"
            echo ""

            parallel --colsep ' ' --jobs "${par_job}" --dryrun \
                "${cmd}" \
                :::: "${config}"

            echo ""
            echo ""
        fi

        if ! ${dry_run}; then
            parallel --colsep ' ' --jobs "${par_job}" \
                "${cmd}" \
                :::: "${config}"
        fi
    else
        #  If --slurm was not specified and 'par_job=1', then run jobs in
        #+ serial
        if ${dry_run} || ${verbose}; then
            echo "##################################################"
            echo "## Serialized job execution via 'submit' script ##"
            echo "##################################################"
            echo ""
            echo "bash ${scr_sub} \\"
            echo "    -en ${env_nam} \\"
            echo "    -sf ${scr_fnc} \\"
            echo "     -t ${threads} \\"
            echo "     -a ${aligner} \\"
            echo "    -at ${a_type} \\"
            echo "    -mq ${mapq} \\"
            if ${req_flg}; then echo "    -rf \\"; fi
            echo "    -ix ${index} \\"
            echo "     -i ${infiles} \\"
            echo "    -do ${dir_out} \\"
            if ${qname}; then echo "    -qn \\"; fi
            echo "    -ss ${sfx_se} \\"
            echo "    -sp ${sfx_pe} \\"
            echo "    -eo ${err_out} \\"
            echo "    -nj ${nam_job} \\"
            echo "         > ${err_out}/${nam_job}_ser.stdout.txt \\"
            echo "        2> ${err_out}/${nam_job}_ser.stderr.txt"
            echo ""
            echo ""
        fi

        # shellcheck disable=SC2046
        bash "${scr_sub}" \
            -en "${env_nam}" \
            -sf "${scr_fnc}" \
             -t "${threads}" \
             -a "${aligner}" \
            -at "${a_type}" \
            -mq "${mapq}" \
            $(if ${req_flg}; then echo "-rf"; fi) \
            -ix "${index}" \
             -i "${infiles}" \
            -do "${dir_out}" \
            $(if ${qname}; then echo "-qn"; fi) \
            -ss "${sfx_se}" \
            -sp "${sfx_pe}" \
            -eo "${err_out}" \
            -nj "${nam_job}" \
                 > "${err_out}/${nam_job}_ser.stdout.txt" \
                2> "${err_out}/${nam_job}_ser.stderr.txt"
    fi
fi
