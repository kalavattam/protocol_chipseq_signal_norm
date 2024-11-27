#!/bin/bash

#  execute_deeptools_compare.sh
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
    source "${dir_fnc}/check_array_files.sh"
    source "${dir_fnc}/check_arrays_lengths.sh"
    source "${dir_fnc}/check_exists_file_dir.sh"
    source "${dir_fnc}/check_flt_pos.sh"
    source "${dir_fnc}/check_format_time.sh"
    source "${dir_fnc}/check_int_pos.sh"
    source "${dir_fnc}/check_mut_excl_args.sh"
    source "${dir_fnc}/check_mut_excl_flags.sh"
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/check_region.sh"
    source "${dir_fnc}/check_region_bam.sh"
    source "${dir_fnc}/check_str_delim.sh"
    source "${dir_fnc}/check_supplied_arg.sh"
    source "${dir_fnc}/check_table.sh"
    source "${dir_fnc}/check_table_column.sh"
    source "${dir_fnc}/check_table_scaling_factor.sh"
    source "${dir_fnc}/debug_array_contents.sh"
    source "${dir_fnc}/echo_error.sh"
    source "${dir_fnc}/echo_warning.sh"
    source "${dir_fnc}/exit_0.sh"
    source "${dir_fnc}/exit_1.sh"
    source "${dir_fnc}/handle_env.sh"
    source "${dir_fnc}/populate_array_empty.sh"
    source "${dir_fnc}/reset_max_job.sh"
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
        det_bam="flag-${flg}_mapq-${mapq}"
        det_cov="${aligner}_${a_type}_${det_bam}"
    }
    # dir_aln="${dir_pro}/align_${aligner}_${a_type}"
    # dir_bam="${dir_aln}/${det_bam}/sc"
    dir_cov="${dir_pro}/compute_coverage/${det_cov}"
    cov_nrm="norm"
    dir_bwg="${dir_cov}/${cov_nrm}/tracks"
    # dir_tbl="${dir_cov}/${cov_nrm}/tables"
    dir_trk="${dir_cov}/log2/${cov_nrm}/tracks"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=true
    threads=8
    fil_num="$(  ## WARNING: Change the search parameters as needed ##
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_bwg}" \
            --pattern "*.bw" \
            --include "IP*,*Hho1*"  # "IP*,*Hmo1*"  # "IP*,*Brn1*"
    )"
    fil_den="$(  ## WARNING: Change the search parameters as needed ##
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_bwg}" \
            --pattern "*.bw" \
            --include "in*,*Hho1*"  # "in*,*Hmo1*"  # "in*,*Brn1*"
    )"
    dir_out="${dir_trk}"
    typ_out="bigwig"
    bin_siz=1  # 10
    region=""
    scl_fct=""  # "0.002054,0.003138,0.003127,0.003522,0.056611,0.02906"  
    norm=""  # "raw"  # "none"  # "rpkm"  # "fpkm"  # "cpm"  # "bpm"  # "rpgc"
    exact=true
    usr_frg=""
    err_out="${dir_trk}/logs"
    nam_job="run_deeptools_compare"
    slurm=true
    max_job=4  # 6
    time="0:30:00"
}


source activate env_analyze

#  For BIGWIG
typ_fil="bw"
verbose=true
threads=8
fil_num="/home/kalavatt/tsukiyamalab/Kris/202X_protocol_ChIP/data/processed/compute_coverage/bowtie2_global_flag-2_mapq-1/norm/tracks/IP_WT_Q_Hho1_6337.sc"
fil_den="/home/kalavatt/tsukiyamalab/Kris/202X_protocol_ChIP/data/processed/compute_coverage/bowtie2_global_flag-2_mapq-1/norm/tracks/in_WT_Q_Hho1_6337.sc"
dir_out="/home/kalavatt/tsukiyamalab/Kris/202X_protocol_ChIP/data/processed/compute_coverage/bowtie2_global_flag-2_mapq-1/log2/norm/tracks"
outfile="${dir_out}/log2_IP_in_Q_Hho1_6337.sc"
typ_out="bigwig"
operation="log2"
bin_siz=1
region="#N/A"
scl_fct="#N/A"

#  For BAM (#TODO: Change 'raw' to 'depth' because alignment count-based scaling is still happening)
typ_fil="bam"
verbose=true
threads=8
fil_num="/home/kalavatt/tsukiyamalab/Kris/202X_protocol_ChIP/data/processed/align_bowtie2_global/flag-2_mapq-1/sc/IP_WT_Q_Hho1_6337.sc"
fil_den="/home/kalavatt/tsukiyamalab/Kris/202X_protocol_ChIP/data/processed/align_bowtie2_global/flag-2_mapq-1/sc/in_WT_Q_Hho1_6337.sc"
dir_out="/home/kalavatt/tsukiyamalab/Kris/202X_protocol_ChIP/data/processed/compute_coverage/bowtie2_global_flag-2_mapq-1/log2/depth/tracks"
outfile="${dir_out}/log2_IP_in_Q_Hho1_6337.sc"
typ_out="bigwig"
operation="log2"
bin_siz=1
region="#N/A"
scl_fct="#N/A"
norm="#N/A"
method="readCount"
exact=true
typ_seq="paired"
usr_frg="#N/A"

# shellcheck disable=SC2046
case "${typ_fil}" in
    bedgraph|bg|bigwig|bw)
        bigwigCompare \
            $(if ${verbose}; then echo "--verbose"; fi) \
            --numberOfProcessors "${threads}" \
            --bigwig1 "${fil_num}.${typ_fil}" \
            --bigwig2 "${fil_den}.${typ_fil}" \
            --outFileName "${outfile}.${typ_out}" \
            --outFileFormat "${typ_out}" \
            --operation "${operation}" \
            $(
                if [[ "${operation}" =~ ^(log2|ratio)$ ]]; then
                    echo "--pseudocount 1 1"
                fi
            ) \
            --binSize "${bin_siz}" \
            $(
                if [[ "${region}" != "#N/A" ]]; then
                    echo "--region ${region}"
                fi
            ) \
            $(
                if [[ "${scl_fct}" != "#N/A" ]]; then
                    echo "--scaleFactor ${scl_fct}"
                fi
            ) \
            --skipZeroOverZero \
            --skipNonCoveredRegions
        ;;
    bam)
        bamCompare \
            $(if ${verbose}; then echo "--verbose"; fi) \
            --numberOfProcessors "${threads}" \
            --bamfile1 "${fil_num}.${typ_fil}" \
            --bamfile2 "${fil_den}.${typ_fil}" \
            --outFileName "${outfile}.${typ_out}" \
            --outFileFormat "${typ_out}" \
            --operation "${operation}" \
            $(
                if [[ "${operation}" =~ ^(log2|ratio)$ ]]; then
                    echo "--pseudocount 1 1"
                fi
            ) \
            --binSize "${bin_siz}" \
            $(
                if [[ "${region}" != "#N/A" ]]; then
                    echo "--region ${region}"
                fi
            ) \
            --skipZeroOverZero \
            --skipNonCoveredRegions \
            $(
                if [[ "${scl_fct}" != "#N/A" ]]; then
                    echo "--scaleFactor ${scl_fct}"
                    echo "--scaleFactorsMethod None"
                elif [[ "${norm}" != "#N/A" ]]; then
                    echo "--normalizeUsing ${norm}"
                    if [[ "${norm}" == "RPGC" ]]; then
                        echo "--effectiveGenomeSize ${gen_siz:-11624332}"
                    fi
                else
                    echo "--scaleFactorsMethod ${method:-readCount}"
                fi
            ) \
            $(if ${exact}; then echo "--exactScaling"; fi) \
            $(
                if [[ "${typ_seq}" == "paired" ]]; then
                    echo "--samFlagInclude 64"
                fi
            ) \
            $(
                if [[
                       "${typ_seq}" == "paired" && "${usr_frg}" == "#N/A"
                ]]; then
                    echo "--extendReads"
                elif [[ "${usr_frg}" != "#N/A" ]]; then
                    echo "--extendReads ${usr_frg}"
                fi
            )
        ;;
esac
