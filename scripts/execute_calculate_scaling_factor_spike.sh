#!/bin/bash

#  execute_calculate_scaling_factor_spike.sh
#  KA


#  Run script in interactive/test mode (true) or command-line mode (false)
interactive=false

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive}; then set -euo pipefail; fi

#  Set the path to the "scripts" directory
if ${interactive}; then
    ## WARNING: If interactive=true, change path as needed ##
    dir_scr="${HOME}/repos/protocol_chipseq_signal_norm/scripts"
else
    dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi


#  Source and define functions ================================================
#  Set the path to the "functions" directory
dir_fnc="${dir_scr}/functions"

# shellcheck disable=SC1091
{
    source "${dir_fnc}/check_exists_file_dir.sh"
    source "${dir_fnc}/check_format_time.sh"
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
    ## WARNING: If interactive=true, change values as needed ##
    dir_bas="${HOME}/repos"
    dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
    dir_scr="${dir_rep}/scripts"
    dir_dat="${dir_rep}/data"
    dir_pro="${dir_dat}/processed"
    {
        aligner="bowtie2"
        a_type="global"
        req_flg=true
        flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
        mapq=1
    }
    dir_aln="${dir_pro}/align_${aligner}_${a_type}"
    dir_bam="${dir_aln}/flag-${flg}_mapq-${mapq}/sc"
    dir_out="${dir_pro}/calculate_scaling_factor_spike"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=true
    threads=8
    infiles="$(  ## WARNING: Change search parameters as needed ##
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_bam}" \
            --pattern "*.bam" \
            --include "IP*,*Hho1*"
    )"
    outfile="${dir_out}/IP-in_WT_G1-G2M-Q_Hho1_6336-6337.mc.txt"
    flg_in=true
    flg_mc=true
    err_out="${dir_out}/logs"
    nam_job="calc_sf_spike"
    slurm=false
    max_job=6
    time="0:30:00"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_analyze"
scr_sub="${dir_scr}/submit_calculate_scaling_factor_spike.sh"
scr_spk="${dir_scr}/calculate_scaling_factor_spike.py"
denom=4

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=1
infiles=""
outfile=""
flg_in=false
flg_mc=false
err_out=""
nam_job="calc_sf_spike"
slurm=false
max_job=6
time="0:30:00"

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  execute_calculate_scaling_factor_spike.sh
    [--verbose] [--dry_run] --threads <int> --infiles <str> --outfile <str>
    [--flg_in] [--flg_mc] --err_out <str> --nam_job <str> [--slurm]
    [--max_job <int>] [--time <str>]

Description:
  execute_calculate_scaling_factor_spike.sh orchestrates the calculation of
  spike-in-derived scaling factors for ChIP-seq data. The main workflow
  involves counting alignments in main and spike-in S. cerevisiae BAM files,
  calculating scaling coefficients using a custom Python script, and outputting
  results to a specified table file. To handle groups of sample BAM data
  efficiently, the script uses parallelized processing via SLURM scheduling or
  GNU Parallel.

Arguments:
   -h, --help     Display this help message and exit (0).
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Perform a dry run without executing commands (optional).
   -t, --threads  Number of threads to use (default: ${threads}).
   -i, --infiles  Comma-separated list of S. cerevisiae IP BAM files that are 
                  coordinate-sorted. (See Notes for file name and path
                  structure requirements.)
   -o, --outfile  Tab-delimited text outfile in which calculated
                  spike-in-derived scaling factors and, optionally, additional
                  metrics are saved. If the file already exists, new data will
                  be appended.
  -fi, --flg_in   Include sample input coefficients in outfile.
  -fm, --flg_mc   Include additional measurements and calculations in outfile.
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (default: \${dir_out}/err_out).
  -nj, --nam_job  The name of the job (default: ${nam_job}).
  -sl, --slurm    Submit jobs to the SLURM scheduler.
  -mj, --max_job  The maximum number of jobs to run at one time (required if
                  --slurm is specified, ignored if not; default: ${max_job}).
  -tm, --time     The length of time, in 'h:mm:ss' format, for the SLURM job
                  (required if --slurm is specified, ignored if not; default:
                  ${time}).

Dependencies:
  - Programs
    + awk
    + Bash or Zsh
    + Conda
    + GNU Parallel (if '--slurm' is not specified)
    + Python
    + Samtools
    + SLURM (if '--slurm' is specified)
  - Functions
    + check_exists_file_dir
    + check_format_time
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
  - If '--dry_run' is enabled, the script will output commands that would be
    executed without actually running them.
  - For non-SLURM job submissions, GNU Parallel is used. This requires
    specifying a set number of parallel jobs ('par_job') and the number of
    threads per job ('threads'). The number of parallel jobs is determined by
    dividing the user-specified 'threads' value by a denominator ('denom'),
    which is set (hardcoded) to 4. The value of 'threads' is then reset to
    'denom'. For example, if the user-specified 'threads' value is 8, 'par_job'
    would be set to 2 (8 divided by 4) and 'threads' would be reset to 4. Thus,
    2 jobs would be run by GNU Parallel (because 'par_job=2'), each with 4
    threads (because 'threads=4').
  - This and the accompanying submission script assume a specific directory and
    naming structure for the infiles:
      + The primary infiles are coordinate-sorted S. cerevisiae ("sc") IP BAM
        files.
      + Based on the paths of these files, the script derives paths to
        additional files: S. cerevisiae input files and S. pombe ("sp") IP and
        input files:
        \`\`\`
        align_bowtie2_global/flag-2_mapq-1/
        ├── init
        │   ├── docs
        │   └── logs
        ├── sc
        │   ├── docs
        │   └── logs
        └── sp
            ├── docs
            └── logs
        \`\`\`
      + The paths are generated by making systematic substitutions to the file
        path and names, saved in variables 'mp', 'mn', 'sp', and 'sn':
          - sp: Derived from 'mp' (main S. cerevisiae IP file) by replacing
            subdirectory /sc/ with /sp/ and string .sc. with .sp.
          - mn: Derived from 'mp' by replacing string IP_ with in_.
          - sn: Derived from 'sp' by replacing string IP_ with in_.
      + The following logic is used for file name and path manipulation:
        \`\`\`
        mp={arr_infiles}
        sp=\$(echo "\${mp}" | sed "s:/sc/:/sp/:g; s:.sc.:.sp.:g")
        mn=\$(echo "\${mp}" | sed "s:/IP_:/in_:/g")
        sn=\$(echo "\${sp}" | sed "s:/IP_:/in_:/g")
        \`\`\`
      + If any of these files are missing, the script will print an error
        message and exit.

Examples:
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
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -v|--verbose) verbose=true;   shift 1 ;;
            -dr|--dry_run) dry_run=true;   shift 1 ;;
             -t|--threads) threads="${2}"; shift 2 ;;
             -i|--infiles) infiles="${2}"; shift 2 ;;
             -o|--outfile) outfile="${2}"; shift 2 ;;
            -fi|--flg_in)  flg_in=true;    shift 1 ;;
            -fm|--flg_mc)  flg_mc=true;    shift 1 ;;
            -eo|--err_out) err_out="${2}"; shift 2 ;;
            -nj|--nam_job) nam_job="${2}"; shift 2 ;;
            -sl|--slurm)   slurm=true;     shift 1 ;;
            -mj|--max_job) max_job="${2}"; shift 2 ;;
            -tm|--time)    time="${2}";    shift 2 ;;
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

check_supplied_arg -a "${scr_spk}" -n "scr_spk"
check_exists_file_dir "f" "${scr_spk}" "scr_spk"

check_supplied_arg -a "${denom}" -n "denom"

check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

check_supplied_arg -a "${infiles}" -n "infiles"
check_exists_file_dir "d" "$(dirname "${infiles%%[,;]*}")" "infiles"

check_supplied_arg -a "${outfile}" -n "outfile"
check_exists_file_dir "d" "$(dirname "${outfile}")" "outfile"

if [[ -n "${err_out}" ]]; then
    check_exists_file_dir "d" "${err_out}" "err_out"
elif [[ -z "${err_out}" ]]; then
    err_out="${dir_out}/err_out"
    check_exists_file_dir "d" "${err_out}" "err_out"
fi

check_supplied_arg -a "${nam_job}" -n "nam_job"

if ${slurm}; then
    check_supplied_arg -a "${max_job}" -n "max_job"
    check_int_pos "${max_job}" "max_job"
    
    check_supplied_arg -a "${time}" -n "time"
    check_format_time "${time}"
else
    #  For non-SLURM submissions, calculate the number of parallel jobs
    #+ ('par_job') by dividing the user-defined 'threads' value by the
    #+ hardcoded denominator ('denom'); then, reset 'threads' to the value of
    #+ 'denom' for consistency in processing
    par_job=$(( threads / denom ))
    threads=${denom}

    check_supplied_arg -a "${par_job}" -n "par_job"
    check_int_pos "${par_job}" "par_job"

    check_supplied_arg -a "${threads}" -n "threads"
    check_int_pos "${threads}" "threads"
fi

#  Activate environment and check that dependencies are in PATH
handle_env "${env_nam}" > /dev/null

check_program_path awk
if ! ${slurm}; then check_program_path parallel; fi
check_program_path python
check_program_path samtools
if ${slurm}; then check_program_path sbatch; fi

#  Parse the --infiles argument into an array, then validate the infile value
#+ assignments
IFS=',' read -r -a arr_infiles <<< "${infiles}"  # unset arr_infiles
# for infile in "${arr_infiles[@]}"; do echo "${infile}"; done  # unset infile

#  Check that each infile exists; if not, exit
for infile in "${arr_infiles[@]}"; do
    check_exists_file_dir "f" "${infile}" "infile"
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
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "env_nam=${env_nam}"
    echo "scr_sub=${scr_sub}"
    echo "scr_spk=${scr_spk}"
    echo "denom=${denom}"
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
    echo "infiles=${infiles}"
    echo "outfile=${outfile}"
    echo "flg_in=${flg_in}"
    echo "flg_mc=${flg_mc}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo "slurm=${slurm}"
    echo "max_job=${max_job}"
    echo "time=${time}"
    echo ""
    echo ""
fi

# shellcheck disable=SC1083,SC2157,SC2046,SC2086
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
        echo "    --export=dir_scr=${dir_scr} \\"
        echo "    ${scr_sub} \\"
        echo "        --threads ${threads} \\"
        echo "        --infiles ${infiles} \\"
        echo "        --outfile ${outfile} \\"
        echo "        $(if ${flg_in}; then echo "--flg_in"; fi) \\"
        echo "        $(if ${flg_mc}; then echo "--flg_mc"; fi) \\"
        echo "        --err_out ${err_out} \\"
        echo "        --nam_job ${nam_job} \\"
        echo "        --env_nam ${env_nam} \\"
        echo "        --scr_spk ${scr_spk}"
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
        #  To prevent potential race conditions from concurrent writes,
        #+ pre-write the header to the outfile before running SLURM jobs
        #+ (for more details, see related comment for GNU Parallel jobs below)
        if ! ${flg_mc}; then
            echo -e "sample\tsf" >> "${outfile}"
        else
            echo -e \
                "sample\tsf\tmain_ip\tspike_ip\tmain_in\tspike_in" \
                    >> "${outfile}"
        fi

        sbatch \
            --job-name=${nam_job} \
            --nodes=1 \
            --cpus-per-task=${threads} \
            --time=${time} \
            --error=${err_out}/${nam_job}.%A-%a.stderr.txt \
            --output=${err_out}/${nam_job}.%A-%a.stdout.txt \
            --array=1-${#arr_infiles[@]}%${max_job} \
            --export=dir_scr=${dir_scr} \
            ${scr_sub} \
                --threads ${threads} \
                --infiles ${infiles} \
                --outfile ${outfile} \
                $(if ${flg_in}; then echo "--flg_in"; fi) \
                $(if ${flg_mc}; then echo "--flg_mc"; fi) \
                --err_out ${err_out} \
                --nam_job ${nam_job} \
                --env_nam ${env_nam} \
                --scr_spk ${scr_spk}
    fi
else
    #  Define the base GNU Parallel command, specifying the number of jobs to
    #+ run, using named headers for input variables, and setting a comma as the
    #+ column separator for input data
    cmd_par="parallel --jobs ${par_job} --header : --colsep ','"
    
    #  Define the main processing logic to be run by GNU Parallel; this block
    #+ handles input file manipulation, alignment counting, scaling factor
    #+ calculation, and output writing
    # shellcheck disable=SC2016,SC2089
    logic='
        #  Set variables for main and spike-in input files based on the input
        #+ filename pattern
        mp={arr_infiles}
        sp=$(echo "${mp}" | sed "s:\/sc\/:\/sp\/:g; s:\.sc\.:\.sp\.:g")
        mn=$(echo "${mp}" | sed "s:\/IP_:\/in_:g")
        sn=$(echo "${sp}" | sed "s:\/IP_:\/in_:g")

        #  Check for required input files; if missing, print an error message
        #+ and exit
        if [[
            ! -f "${mp}" || ! -f "${sp}" || ! -f "${mn}" || ! -f "${sn}"
        ]]; then
            echo \
                "Error: Missing input file(s) for mp=${mp}, sp=${sp}," \
                "mn=${mn}, sn=${sn}." >&2
            exit 1
        fi

        #  Count the number of alignments in each file using Samtools
        num_mp=$(samtools view -c -@ {threads} "${mp}")
        num_sp=$(samtools view -c -@ {threads} "${sp}")
        num_mn=$(samtools view -c -@ {threads} "${mn}")
        num_sn=$(samtools view -c -@ {threads} "${sn}")

        #  Extract the sample name from the main input file for naming output
        #+ entries
        samp=$(basename "${mp}" | sed "s:IP_::; s:.bam::; s:\\.:_:g")

        #  Calculate the scaling factor (coefficient) using a specified Python
        #+ script
        sf=$(
            python {scr_spk} \
                --main_ip "${num_mp}" \
                --spike_ip "${num_sp}" \
                --main_in "${num_mn}" \
                --spike_in "${num_sn}"
        )

        # #  Write a header to the output file if it is not already present
        # #+ using file locking to prevent concurrent writes from multiple
        # #+ processes
        # {
        #     flock -n 200 || exit 1
        #     if ! \
        #         grep -q "$(printf "^sample\tsf")" {outfile} 2> /dev/null
        #     then
        #         if ! {flg_mc}; then
        #             echo -e "sample\tsf" >> {outfile}
        #         else
        #             echo -e \
        #                 "sample\tsf\tmain_ip\tspike_ip\tmain_in\tspike_in" \
        #                     >> {outfile}
        #         fi
        #     fi
        # } 200> "{outfile}.lock"

        #  Write the calculated scaling factor and optional metrics to the
        #+ output file
        if ! {flg_mc}; then
            echo -e "${mp}\t${sf}" >> {outfile}
            if {flg_in}; then
                echo -e "${mn}\t1" >> {outfile}
            fi
        else
            echo -e \
                "${mp}\t${sf}\t${num_mp}\t${num_sp}\t${num_mn}\t${num_sn}" \
                    >> {outfile}
            if {flg_in}; then
                echo -e \
                    "${mn}\t1\t${num_mp}\t${num_sp}\t${num_mn}\t${num_sn}" \
                        >> {outfile}
            fi
        fi
    '
    #NOTE: 2024-1119, replaced ${mp##*/} with ${mp} and ${mn##*/} with ${mn}

    #  Define a helper function to run GNU Parallel with the specified logic
    #+ and variables
    function run_parallel() {
        local cmd="${1}"

        ${cmd} "${logic}" \
            ::: arr_infiles "${arr_infiles[@]}" \
            ::: threads "${threads}" \
            ::: scr_spk "${scr_spk}" \
            ::: flg_mc "${flg_mc}" \
            ::: flg_in "${flg_in}" \
            ::: outfile "${outfile}"
    }

    #  If `dry_run` or `verbose` is enabled, run GNU Parallel in dry-run mode
    #+ to display commands without executing them, capturing output and error
    #+ logs
    if ${dry_run} || ${verbose}; then
        # #  Based on the presence of `flg_mc`, output the appropriate header
        # if ! ${flg_mc}; then
        #     echo -e "sample\tsf"
        # else
        #     echo -e "sample\tsf\tmain_ip\tspike_ip\tmain_in\tspike_in"
        # fi

        #  Execute a dry-run of the parallel processing logic, capturing stdout
        #+ and stderr logs for each job
        run_parallel "${cmd_par} --dryrun" \
             > >(
                tee -a "${err_out}/${nam_job}.$(
                    basename "${outfile}" ".txt"
                ).stdout.txt"
            ) \
            2> >(
                tee -a "${err_out}/${nam_job}.$(
                    basename "${outfile}" ".txt"
                ).stderr.txt"
            )
    fi

    #  Run GNU Parallel for actual processing if `dry_run` is not enabled,
    #+ capturing output and error logs
    if ! ${dry_run}; then
        #  To prevent potential race conditions from concurrent writes,
        #+ pre-write the header to the outfile before running GNU Parallel jobs
        #+ (this avoids relying on `flock`, which seems to exhibit buggy
        #+ behavior in GNU Parallel job submissions); if `flg_mc=true`, include
        #+ additional column names in the header
        if ! ${flg_mc}; then
            echo -e "sample\tsf" >> "${outfile}"
        else
            echo -e \
                "sample\tsf\tmain_ip\tspike_ip\tmain_in\tspike_in" \
                    >> "${outfile}"
        fi

        #  Execute the parallel processing logic, capturing stdout and stderr
        #+ logs for each job
        run_parallel "${cmd_par}" \
             >> >(
                tee -a "${err_out}/${nam_job}.$(
                    basename "${outfile}" ".txt"
                ).stdout.txt"
            ) \
            2>> >(
                tee -a "${err_out}/${nam_job}.$(
                    basename "${outfile}" ".txt"
                ).stderr.txt"
            )
    fi
fi
