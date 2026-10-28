#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_execute_align_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
#   GPT-5.6);
# - Anthropic Claude Code (Opus 5, Opus 5.5).
#
# Distributed under the MIT license.


# TODO: add support for aligner 'minibwa' to 'align_fastqs' chain.
# shellcheck disable=SC2154
function help_execute_align_fastqs() {
    cat >&2 << EOM
Usage
-----
  execute_align_fastqs.sh
    [--help] [--verbose] [--dry_run]
    [--env_nam <str>] [--threads <int>]
    [--aligner <spec>] [--bt2_mode <spec>] [--bt2_X <int>] [--bwa_alg <spec>] [--mapq <int>] [--req_flg] [--unal_mode <mode>]
    --index <path> --csv_fil_in <csv> [--ref_fa <file>]
    --dir_out <dir> [--out_ext <format>] [--qname] [--sfx_se <str>] [--sfx_pe <str>]
    [--dir_eo <dir>] [--nam_job <str>]
    [--max_job <int>] [--slurm] [--time <time>]

  Align single- or paired-end short-read FASTQ files using Bowtie 2, BWA, or BWA-MEM2, followed by post-alignment processing with Samtools, including MAPQ-based filtering, sorting, mate fixing (for paired-end alignments), duplicate marking, and indexing.

  Jobs may be run through Slurm ('--slurm'), GNU Parallel, or serial execution, depending on user arguments and the resolved number of jobs.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -v, --verbose : flag
    Run script in verbose mode.

  -dr, --dry_run : flag
    Run script in dry-run mode.

    Print the Slurm, GNU Parallel, or serial command that would run, then call 'submit_align_fastqs.sh' locally with '--dry_run' to print the commands for each input entry, without submitting or running them or writing logs.

  -en, --env, --env_nam : str
    Conda environment to activate (default: '${env_nam}').

  -t, --thr, --threads : int
    Number of threads to use (default: ${threads}).

  -a, --aln, --aligner : {'bowtie2', 'bwa', 'bwa-mem2'}
    Alignment program to use: 'bowtie2', 'bwa', or 'bwa-mem2' (default: '${aligner}').

  -2m, --bt2_mode : {'local', 'global', 'end-to-end'}
    Bowtie 2 alignment type when '--aligner bowtie2': 'local', 'global', or 'end-to-end' (default: '${bt2_mode}'); refused otherwise.

  -2X, --bt2_X : int
    Maximum fragment length for a Bowtie 2 paired-end alignment, passed as 'bowtie2 -X' (default: ${bt2_X}); refused unless '--aligner bowtie2'; ignored with a warning for single-end data.

  -ba, --bwa_alg : {'mem', 'aln'}
    BWA algorithm when '--aligner bwa': 'mem' or 'aln' (default: '${bwa_alg}').

    With '--aligner bwa-mem2', which has only 'mem', 'mem' is ignored with a warning and 'aln' is refused; with '--aligner bowtie2', it is refused.

  -mq, --mapq : int
    MAPQ threshold for filtering alignment output files (default: ${mapq}).

    To disable MAPQ-based filtering, specify 0.

    With paired-end data, a pair is kept only when both mates pass. Secondary (SAM flag 256) and supplementary (SAM flag 2048) alignments are kept or dropped with their read, whatever their own MAPQ.

    With Bowtie 2, a MAPQ of 1 or more does not mean a unique alignment, as a read alignment with multiple equally good placements can get MAPQ 1, and its placement is chosen at random.

  -rq, --req_flg : flag
    Require SAM flag bit 2 for properly paired alignments.

    The flag is checked as the aligner set it and again after 'samtools fixmate', which clears it on a pair whose mates map to different references or align in the same orientation (both forward or both reverse).

    Ignored with a warning when every entry is single-end.

  -um, --unal_mode : {'remove', 'separate'}
    What to do with the records that the alignment file leaves out: 'remove' or 'separate' (default: '${unal_mode}'). These are the unaligned reads and, in paired-end data, their mates (aligned or not), plus the templates that '--req_flg' and '--mapq' reject.

    With 'remove', they are dropped. With 'separate', nothing is dropped: each goes to a diagnostic file beside the alignment file, and the alignment file holds the same records as with 'remove' (see Notes).

    Refused as 'separate' with '--bt2_inherited', whose call writes no unaligned reads.

  -ix, --index : path
    Path to the aligner index/reference.

    If using Bowtie 2, the path should end with the index stem, e.g., "\${HOME}/path/stem".

    If using BWA or BWA-MEM2, the path should be the indexed reference FASTA path, e.g., "\${HOME}/path/stem.fa".

  -ci, --csv_fil_in : list of structured string
    Semicolon-delimited serialized string of FASTQ input entries.

    For single-end data, each entry is one FASTQ file. For paired-end data, each entry contains a comma-delimited FASTQ pair.

    E.g., "\${HOME}/path/samp_1.fastq.gz;\${HOME}/path/samp_2_R1.fastq.gz,\${HOME}/path/samp_2_R2.fastq.gz;\${HOME}/path/samp_3.fastq.gz".

  -rf, --ref_fa : file
    Reference FASTA file. Required when '--out_ext cram'; otherwise, ignored with a warning.

  -do, --dir_out : dir
    Output directory for alignment files.

  -ox, --out_ext : {'bam', 'cram'}
    Final output extension for alignment files: 'bam' or 'cram' (default: '${out_ext}').

  -qn, --qname : flag
    Also write a queryname-sorted copy of the final alignment file.

  -sxs, --sfx_se, --suffix_se : str
    Suffix to strip from single-end FASTQ filenames (default: '${sfx_se}').

    Ignored with a warning when every entry is paired-end.

  -sxp, --sfx_pe, --suffix_pe : str
    Suffix to strip from paired-end FASTQ read-1 filenames (default: '${sfx_pe}').

    Ignored with a warning when every entry is single-end.

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files in TXT format (default: \${dir_out}/logs).

  -nj, --nam_job : str
    Job name used when writing stderr and stdout TXT files (default: '${nam_job}').

  -mj, --max_job : int
    Maximum number of jobs to run concurrently (default: ${max_job}).
      - If '--slurm' is specified, controls Slurm array-task concurrency.
      - If '--slurm' is not specified:
        + If 'max_job' is greater than 1, jobs run in parallel via GNU Parallel.
        + If 'max_job' is 1, jobs run in serial.

  -sl, --slurm : flag
    Submit jobs to the Slurm scheduler; otherwise, run them through GNU Parallel or in serial, as resolved locally.

  -tm, --time : time
    Slurm job time limit, in 'h:mm:ss' format (default: '${time}'); required if '--slurm' is specified, ignored with a warning if not.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - bash >= 4.4
    - bowtie2 (when '--aligner bowtie2' is specified)
    - bwa (when '--aligner bwa' is specified)
    - bwa-mem2 (when '--aligner bwa-mem2' is specified)
    - conda (when the requested environment is not active)
    - md5sum
    - parallel (when '--slurm' is not specified and multiple jobs are run)
    - samtools
    - sbatch (when '--slurm' is specified)

  - When the '--slurm' flag is used, jobs are parallelized via Slurm array tasks. Otherwise, if multiple jobs are to be run, they are parallelized locally via GNU Parallel; if only one job is to be run, execution proceeds in serial.
  - If using Bowtie 2, ensure the path to index files ends with the index stem, e.g., "\${HOME}/path/stem" or "\${HOME}/path/sc_sp_proc". If using BWA or BWA-MEM2, the path should include the reference FASTA filename, e.g., "\${HOME}/path/stem.fa" or "\${HOME}/path/sc_sp_proc.fa".
  - Calling the script with '--qname' also writes a queryname-sorted copy of each final alignment file, with the same records and flags.
  - Retained queryname-sorted output files will share the same path and stem as the final alignment output file, but with '.qnam' inserted before the final extension (for example, '.qnam.bam' or '.qnam.cram').
  - When '--out_ext cram' is used, '--ref_fa' must also be supplied. Although intermediate work files are processed in BAM format, the final output, and any retained queryname-sorted output, are written as CRAM, which requires a reference FASTA.
  - With '--unal_mode separate', every record the aligner writes lands in exactly one file ('docs/design/align_fastqs_separate.md'):
    + '<stem>.unal.bam' (or '.cram'), the side file, holds the unaligned reads and their mates:
      - With single-end data, these are the reads flagged unaligned (SAM flag 4).
      - With Bowtie 2 and paired-end data, they are the pairs that did not align as concordant pairs. Bowtie 2 writes both mates of such a pair unaligned even when one mate could align alone.
      - With BWA or BWA-MEM2, they are the records flagged unaligned or with an unaligned mate.
      - When BWA-backtrack flags one mate unaligned but leaves its partner's mate flag clear, the partner goes beside its mate.
    + The templates that '--req_flg' and '--mapq' reject go to their own files:
      - '<stem>.improp.bam' holds the pairs that '--req_flg' rejects. Single-end data has no improper file, and without '--req_flg' the improper file is empty.
      - '<stem>.mapq.bam' holds the templates that '--mapq' rejects.
      - A template that fails both is improper.
      - '--mapq' and '--req_flg' do not apply to the side file, as both judge aligned pairs and would remove every pair it holds.
    + Each template is kept whole: all of a read pair's records go to the same file.
    + For each diagnostic file:
      - It is coordinate-sorted and indexed.
      - Its duplicates are marked against all the files together.
      - It carries the pipeline record and follows the alignment file into CRAM and into a '--qname' copy.
    + The diagnostic files are not for signal or scaling.
    + '<stem>.diag.tsv.gz' counts the templates and records each file holds, by category.
  - Duplicates are marked (SAM flag 1024) after MAPQ filtering, and kept.
    + Like any other read, '--mapq' filters them by the MAPQ the aligner gave them.
    + When every read name carries Illumina flow-cell coordinates ('...:tile:x:y'), a duplicate is tagged as optical ('dt:Z:SQ') if it shares its original's tile and lies within 2500 pixels of it in x and y. Other duplicates are tagged as library duplicates ('dt:Z:LB').
    + Otherwise, as with most SRA-dumped names, no 'dt' tags are written and a warning says why.
    + Duplicate statistics are written to '<stem>.markdup.txt.gz' beside the alignment file; Samtools estimates library size for paired-end data only. The file's last line, 'OPTICAL DISTANCE:', gives 2500 when optical duplicates were told apart, or 'not used' and the reason.
    + With '--unal_mode separate', duplicates are marked twice:
      - Markdup 1 runs on the alignment file alone, as with 'remove', so its flags and '<stem>.markdup.txt.gz' are the same in both modes.
      - Markdup 2 runs on all files together, so a read in a diagnostic file can be marked as a duplicate of a read in the alignment file. Its statistics are written to '<stem>.all.markdup.txt.gz'.
  - Bowtie 2 runs with '--reorder', so its output order does not depend on the thread count. It also runs with '--seed 0', its default, which the record gives as 'bt2_seed=0'.
    + Which of several equally good alignments a read gets depends on the seed and on the read's sequence, qualities, and whole FASTQ header line, so renaming reads can move them.
    + For paired-end data, the two mates' seeds are combined, which cancels the seed and any header text identical in both mates; only header differences between the mates move pairs.
  - BWA-MEM and BWA-MEM2 run with '-K 100000000', which fixes the batch of read bases over which paired-end insert sizes are estimated. Their output therefore does not depend on the thread count; the record gives the setting as 'bwa_K=100000000'.
  - Each alignment file's header, the '.qnam' copy's included, carries a pipeline record as '@CO' lines:
    + 'pipeline:' lines give the commit, the options in effect, the input file names, and any skipped steps;
    + 'pipeline-run:' lines give the thread count, the start time, the host, the scripts and functions that ran, and the full paths.
  - The record names each input and gives the MD5 checksum of its bytes, as stored: 'index_md5' for the index files, 'fq_1_md5' and 'fq_2_md5' for the FASTQs, and 'ref_fa_md5' for the reference with CRAM output.
  - It also gives the Conda package version and build of the aligner and Samtools ('pkg_bowtie2=2.5.5-h9e91881_0', for example), since a program's own '@PG' version can be wrong: BWA-MEM2 2.3 reports 2.2.1 (bwa-mem2 issue #276).

Examples
--------
  1. Run local Bowtie 2 alignment and write BAM files.
    '''bash
    bash "\${dir_scr}/execute_align_fastqs.sh" \\
        --threads "\${threads}" \\
        --aligner bowtie2 \\
        --bt2_mode global \\
        --index "\${idx_bt2}" \\
        --csv_fil_in "\${csv_fil_in}" \\
        --dir_out "\${dir_out}/bam" \\
        --out_ext bam \\
        --dir_eo "\${dir_out}/logs"
    '''

  2. Run local BWA-MEM2 alignment and write CRAM files.
    '''bash
    bash "\${dir_scr}/execute_align_fastqs.sh" \\
        --threads "\${threads}" \\
        --aligner bwa-mem2 \\
        --index "\${idx_bwa_mem2}" \\
        --csv_fil_in "\${csv_fil_in}" \\
        --ref_fa "\${ref_fa}" \\
        --dir_out "\${dir_out}/cram" \\
        --out_ext cram \\
        --dir_eo "\${dir_out}/logs"
    '''

  3. Submit a Slurm array job.
    '''bash
    bash "\${dir_scr}/execute_align_fastqs.sh" \\
        --threads "\${threads}" \\
        --aligner bowtie2 \\
        --index "\${idx_bt2}" \\
        --csv_fil_in "\${csv_fil_in}" \\
        --dir_out "\${dir_out}/bam" \\
        --dir_eo "\${dir_out}/logs" \\
        --max_job 12 \\
        --slurm \\
        --time 1:00:00
    '''
EOM
}
