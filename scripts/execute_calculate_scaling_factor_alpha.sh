#!/bin/bash

#  execute_calculate_scaling_factor_alpha.sh
#  KA


#  Run script in interactive/test mode (true) or command-line mode (false)
interactive=true

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
    ## WARNING: Change values if you're not Kris and `interactive=true` ##
    dir_bas="${HOME}/tsukiyamalab/Kris"
    dir_rep="${dir_bas}/202X_protocol_ChIP"
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
    dir_aln="${dir_pro}/align_${aligner}_${a_type}_BAM"
    dir_bam="${dir_aln}/flag-${flg}_mapq-${mapq}/sc"
    dir_out="${dir_pro}/calculate_scaling_factor_alpha"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=true
    threads=8
    infiles=$(  ## WARNING: Change the search parameters as needed ##
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_bam}" \
            --pattern "*.bam" \
            --include "IP*,*Hho1*"
    )
    table="${dir_dat}/raw/docs/measurements_siq_chip.tsv"
    outfile="${dir_out}/IP-in_WT_G1-G2M-Q_Hho1_6336-6337.mc.txt"
    flg_dep=true
    flg_len=true
    flg_in=true
    flg_mc=true
    err_out="${dir_out}/logs"
    nam_job="calc_sf_alpha"
    slurm=true
    scr_mng="${HOME}/miniforge3/etc/profile.d/conda.sh"
    max_job=6
    time="0:30:00"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_analyze"
scr_sub="${dir_scr}/submit_calculate_scaling_factor_alpha.sh"
scr_par="${dir_scr}/parse_metadata_siq_chip.py"
scr_alf="${dir_scr}/calculate_scaling_factor_alpha.py"
denom=4

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
threads=1
infiles=""
table=""
outfile=""
flg_dep=false
flg_len=false
flg_in=false
flg_mc=false
err_out=""
nam_job="calc_sf_alpha"
slurm=false
scr_mng="${HOME}/miniforge3/etc/profile.d/conda.sh"
max_job=6
time="0:30:00"

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  execute_calculate_scaling_factor_alpha.sh
    [--verbose] [--dry_run] --threads <int> --infiles <str> --table <str>
    --outfile <str> [--flg_dep] [--flg_len] [--flg_in] [--flg_mc]
    --err_out <str> --nam_job <str> [--slurm] [--scr_mng <str>]
    [--max_job <int>] [--time <str>]

Description:
  execute_calculate_scaling_factor_alpha.sh performs... #TODO

Arguments:
   -h, --help     Display this help message and exit (0).
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Perform a dry run without executing commands (optional).
   -t, --threads  Number of threads to use (required; default: ${threads}).
   -i, --infiles  Comma-separated serialized list of IP coordinate-sorted BAM
                  infiles, including paths (required).
  -tb, --table    #TODO
   -o, --outfile  #TODO
  -fd, --flg_dep  Use Samtools to calculate the number of alignments.
  -fl, --flg_len  Use Samtools and awk to calculate the mean fragment length.
  -fi, --flg_in   Include sample input alpha values in outfile.
  -fm, --flg_mc   Include additional measurements and calculations in outfile.
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (required; default: \${dir_out}/err_out).
  -nj, --nam_job  The name of the job (required; default: ${nam_job}).
  -sl, --slurm    Submit jobs to the SLURM scheduler.
  -sm, --scr_mng  Conda package manager shell script, conda.sh (required if
                  --slurm is specified, ignored if not; default:
                  ${scr_mng}).
  -mj, --max_job  The maximum number of jobs to run at one time (required if
                  --slurm is specified, ignored if not; default: ${max_job}).
  -tm, --time     The length of time, in 'h:mm:ss' format, for the SLURM job
                  (required if --slurm is specified, ignored if not; default:
                  ${time}).

#TODO Notes, examples, etc.
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
            -tb|--table)   table="${2}";   shift 2 ;;
             -o|--outfile) outfile="${2}"; shift 2 ;;
            -fd|--flg_dep) flg_dep=true;   shift 1 ;;
            -fl|--flg_len) flg_len=true;   shift 1 ;;
            -fi|--flg_in)  flg_in=true;    shift 1 ;;
            -fm|--flg_mc)  flg_mc=true;    shift 1 ;;
            -eo|--err_out) err_out="${2}"; shift 2 ;;
            -nj|--nam_job) nam_job="${2}"; shift 2 ;;
            -sl|--slurm)   slurm=true;     shift 1 ;;
            -sm|--scr_mng) scr_mng="${2}"; shift 2 ;;
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

check_supplied_arg -a "${scr_par}" -n "scr_par"
check_exists_file_dir "f" "${scr_par}" "scr_par"

check_supplied_arg -a "${scr_alf}" -n "scr_alf"
check_exists_file_dir "f" "${scr_alf}" "scr_alf"

check_supplied_arg -a "${denom}" -n "denom"

check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

check_supplied_arg -a "${infiles}" -n "infiles"
check_exists_file_dir "d" "$(dirname "${infiles%%[,;]*}")" "infiles"

check_supplied_arg -a "${table}" -n "table"
check_exists_file_dir "f" "${table}" "table"

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
    check_supplied_arg -a "${scr_mng}" -n "scr_mng"
    check_exists_file_dir "f" "${scr_mng}" "scr_mng"

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
    echo "scr_par=${scr_par}"
    echo "scr_alf=${scr_alf}"
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
    echo "table=${table}"
    echo "outfile=${outfile}"
    echo "flg_dep=${flg_dep}"
    echo "flg_len=${flg_len}"
    echo "flg_in=${flg_in}"
    echo "flg_mc=${flg_mc}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo "slurm=${slurm}"
    echo "scr_mng=${scr_mng}"
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
        echo "        --table ${table} \\"
        echo "        --outfile ${outfile} \\"
        echo "        $(if ${flg_dep}; then echo "--flg_dep"; fi) \\"
        echo "        $(if ${flg_len}; then echo "--flg_len"; fi) \\"
        echo "        $(if ${flg_in}; then echo "--flg_in"; fi) \\"
        echo "        $(if ${flg_mc}; then echo "--flg_mc"; fi) \\"
        echo "        --err_out ${err_out} \\"
        echo "        --nam_job ${nam_job} \\"
        echo "        --scr_mng ${scr_mng} \\"
        echo "        --fnc_env ${dir_fnc}/handle_env.sh \\"
        echo "        --env_nam ${env_nam} \\"
        echo "        --scr_par ${scr_par} \\"
        echo "        --scr_alf ${scr_alf}"
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
            echo -e "sample\talpha" >> "${outfile}"
        else
            echo -e \
                "sample\talpha\tmass_ip\tmass_in\tvolume_ip\tvolume_in\tdepth_ip\tdepth_in\tlength_ip\tlength_in" \
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
            --export=dir_scr="${dir_scr}" \
            ${scr_sub} \
                --threads ${threads} \
                --infiles ${infiles} \
                --table ${table} \
                --outfile ${outfile} \
                $(if ${flg_dep}; then echo "--flg_dep"; fi) \
                $(if ${flg_len}; then echo "--flg_len"; fi) \
                $(if ${flg_in}; then echo "--flg_in"; fi) \
                $(if ${flg_mc}; then echo "--flg_mc"; fi) \
                --err_out ${err_out} \
                --nam_job ${nam_job} \
                --scr_mng ${scr_mng} \
                --fnc_env ${dir_fnc}/handle_env.sh \
                --env_nam ${env_nam} \
                --scr_par ${scr_par} \
                --scr_alf ${scr_alf}
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
        #  Set variables for IP and input files based on the input filename
        file_ip={arr_infiles}
        file_in=$(echo "${file_ip}" | sed "s:\/IP_:\/in_:g")

        #  Check for required input files; if missing, print an error message
        #+ and exit
        if [[ ! -f "${file_ip}" || ! -f "${file_in}" ]]; then
            echo \
                "Error: Missing input file(s) for file_ip=${file_ip}," \
                "file_in=${file_in}." >&2
            exit 1
        fi

        #  Run the Python script to parse siQ-ChIP metadata and assign
        #+ variables for measurements
        source <(
            python {scr_par} \
                --text {table} \
                --bam "${file_ip}" \
                --shell
        )

        #  If --flg_dep, calculate sequencing depth (number of alignments) for
        #+ IP and input samples (note the weird single-quoting here)
        if {flg_dep}; then
            depth_ip=$(samtools view -@ {threads} -c "${file_ip}")
            depth_in=$(samtools view -@ {threads} -c "${file_in}")
        fi

        #  If --flg_len, calculate mean fragment length for IP and input samples
        if {flg_len}; then
            length_ip="$(
                samtools view -@ {threads} "${file_ip}" \
                    | awk '\''{
                        if ($9 > 0) { sum += $9; count++ }
                    } END {
                        if (count > 0) { print sum / count }
                    }'\'' 
            )"
            length_in="$(
                samtools view -@ {threads} "${file_in}" \
                    | awk '\''{
                        if ($9 > 0) { sum += $9; count++ }
                    } END {
                        if (count > 0) { print sum / count }
                    }'\'' 
            )"
        fi

        # #  If --flg_len, calculate mean fragment length for IP and input
        # #+ samples (here, the embedded awk scripts are passed as strings using
        # #+ double quotes)
        # if {flg_len}; then
        #     length_ip="$(
        #         samtools view -@ {threads} "${file_ip}" \
        #             | awk "{ \
        #                 if (\$9 > 0) { sum += \$9; count++ } \
        #             } END { \
        #                 if (count > 0) { print sum / count } \
        #             }"
        #     )"
        #     length_in="$(
        #         samtools view -@ {threads} "${file_in}" \
        #             | awk "{ \
        #                 if (\$9 > 0) { sum += \$9; count++ } \
        #             } END { \
        #                 if (count > 0) { print sum / count } \
        #             }"
        #     )"
        # fi

        #  Extract the sample name from the IP file for naming output entries
        samp=$(basename "${file_ip}" | sed "s:IP_::; s:.bam::; s:\\.:_:g")

        #  Calculate the siQ-ChIP alpha value using a specified Python script
        alpha=$(
            python {scr_alf} \
                --mass_ip ${mass_ip} \
                --mass_in ${mass_in} \
                --volume_ip ${volume_ip} \
                --volume_in ${volume_in} \
                --depth_ip ${depth_ip} \
                --depth_in ${depth_in} \
                --length_ip ${length_ip} \
                --length_in ${length_in}
        )

        #  Write the calculated alpha value and optional metrics to the output file
        if ! {flg_mc}; then
            echo -e "${file_ip##*/}\t${alpha}" >> {outfile}
            if {flg_in}; then
                echo -e "${file_in##*/}\t#N/A" >> {outfile}
            fi
        else
            echo -e "${file_ip##*/}\t${alpha}\t${mass_ip}\t${mass_in}\t${volume_ip}\t${volume_in}\t${depth_ip}\t${depth_in}\t${length_ip}\t${length_in}" \
                >> {outfile}
            if {flg_in}; then
                echo -e "${file_in##*/}\t#N/A\t${mass_ip}\t${mass_in}\t${volume_ip}\t${volume_in}\t${depth_ip}\t${depth_in}\t${length_ip}\t${length_in}" \
                    >> {outfile}
            fi
        fi
    '

    #  Define a helper function to run GNU Parallel with the specified logic
    #+ and variables
    function run_parallel() {
        local cmd="${1}"

        ${cmd} "${logic}" \
            ::: arr_infiles "${arr_infiles[@]}" \
            ::: threads "${threads}" \
            ::: scr_par "${scr_par}" \
            ::: scr_alf "${scr_alf}" \
            ::: flg_dep "${flg_dep}" \
            ::: flg_len "${flg_len}" \
            ::: flg_mc "${flg_mc}" \
            ::: flg_in "${flg_in}" \
            ::: outfile "${outfile}"
    }

    #  If `dry_run` or `verbose` is enabled, run GNU Parallel in dry-run mode
    #+ to display commands without executing them, capturing output and error
    #+ logs
    if ${dry_run} || ${verbose}; then
        #  Based on the presence of `flg_mc`, output the appropriate header
        if ! ${flg_mc}; then
            echo -e "sample\talpha"
        else
            echo -e "sample\talpha\tmass_ip\tmass_in\tvolume_ip\tvolume_in\tdepth_ip\tdepth_in\tlength_ip\tlength_in"
        fi

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
            echo -e "sample\talpha" >> "${outfile}"
        else
            echo -e \
                "sample\talpha\tmass_ip\tmass_in\tvolume_ip\tvolume_in\tdepth_ip\tdepth_in\tlength_ip\tlength_in" \
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
