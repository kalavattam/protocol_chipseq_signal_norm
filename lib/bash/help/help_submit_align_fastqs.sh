#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_submit_align_fastqs.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5, Opus 5.5).
#
# Distributed under the MIT license.


# TODO: do we need '--dir_scr' in the examples? Only if using 'sbatch'.
# TODO: include Conda, aligner, and Samtools versions.
# shellcheck disable=SC2154
function help_submit_align_fastqs() {
    cat >&2 << EOM
Usage
-----
  submit_align_fastqs.sh
    [--help] [--dry_run]
    [--env_nam <str>] [--dir_scr <dir>] [--threads <int>]
    [--aligner <spec>] [--bt2_mode <spec>] [--bt2_X <int>] [--bwa_alg <spec>] [--mapq <int>] [--req_flg] [--unal_mode <mode>]
    --index <path> --csv_fil_in <csv> [--ref_fa <file>]
    --dir_out <dir> [--out_ext <format>]
    [--qname] [--sfx_se <str>] [--sfx_pe <str>]
    --dir_eo <dir> [--nam_job <str>]

  Submit or execute one or more FASTQ-alignment jobs by calling the downstream function 'align_fastqs'.

  This wrapper
    - parses a semicolon-delimited list of FASTQ input entries,
    - derives sample names and alignment fil_out paths,
    - activates the requested Conda environment, and then
    - runs alignment either under Slurm array execution or by serial/local iteration, depending on how the script is invoked.

  For each input entry, this script writes log files, except in a dry run, to:

    \${dir_eo}/\${nam_job}.\${samp}.stdout.txt
    \${dir_eo}/\${nam_job}.\${samp}.stderr.txt

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -dr, --dry_run : flag
    Run script in dry-run mode.

    For each input entry, print the 'align_fastqs' call with its log redirections, then the commands it would run, to stdout, without running them or writing files or logs.

    A dry run always iterates over the entries locally, even inside a Slurm array task.

  -en, --env, --env_nam : str
    Conda environment to activate (default: '${env_nam}').

  -ds, --dir_scr : dir
    Directory containing scripts and functions.

    Passed by the 'execute_*.sh' wrappers, and needed when this script is run from a copy, as 'sbatch <script>' does, rather than from its real path.

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
    Directory to write alignment output files.

  -ox, --out_ext : {'bam', 'cram'}
    Final output extension for alignment files: 'bam' or 'cram' (default: '${out_ext}').

  -qn, --qnam, --qname : flag
    Also write a queryname-sorted copy of the final alignment file.

  -sxs, --sfx_se, --suffix_se : str
    Suffix to strip from single-end FASTQ filenames. Required when an entry is single-end; otherwise, ignored with a warning.

  -sxp, --sfx_pe, --suffix_pe : str
    Suffix to strip from paired-end FASTQ read-1 filenames. Required when an entry is paired-end; otherwise, ignored with a warning.

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files.

  -nj, --nam_job : str
    Job name (default: '${nam_job}').

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
    - samtools
    - Slurm allocation (when run as a Slurm array task)

  - '--ref_fa' is required when '--out_ext cram'.
  - '--req_flg' and '--qname' are optional flags.
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
    + '<stem>.diag.tsv' counts the templates and records each file holds, by category.
  - Duplicates are marked (SAM flag 1024) after MAPQ filtering, and kept.
    + Like any other read, '--mapq' filters them by the MAPQ the aligner gave them.
    + When every read name carries Illumina flow-cell coordinates ('...:tile:x:y') in one format, duplicates on the same tile as their original and within 2500 pixels of it in x and y are tagged as optical ('dt:Z:SQ'); the rest are tagged as library duplicates ('dt:Z:LB').
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
  1. Align one single-end library with Bowtie 2 and write BAM output.
    '''bash
    bash "\${dir_scr}/submit_align_fastqs.sh" \\
        --env_nam "env_protocol" \\
        --dir_scr "\${dir_scr}" \\
        --threads 4 \\
        --aligner "bowtie2" \\
        --bt2_mode "global" \\
        --mapq 1 \\
        --index "\${dir_ref}/bowtie2/genome" \\
        --csv_fil_in "\${dir_fastq}/sample.fastq.gz" \\
        --dir_out "\${dir_out}" \\
        --out_ext "bam" \\
        --sfx_se ".fastq.gz" \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "align_bowtie2"
    '''

  2. Align paired-end reads with BWA-MEM2 and write CRAM output.
    '''bash
    bash "\${dir_scr}/submit_align_fastqs.sh" \\
        --env_nam "env_protocol" \\
        --dir_scr "\${dir_scr}" \\
        --threads 8 \\
        --aligner "bwa-mem2" \\
        --mapq 1 \\
        --index "\${dir_ref}/bwa-mem2/genome.fa" \\
        --csv_fil_in "\${dir_fastq}/sample_R1.fastq.gz,\${dir_fastq}/sample_R2.fastq.gz" \\
        --ref_fa "\${dir_ref}/genome.fa" \\
        --dir_out "\${dir_out}" \\
        --out_ext "cram" \\
        --req_flg \\
        --sfx_pe "_R1.fastq.gz" \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "align_bwa_mem2"
    '''
EOM
}
