#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_execute_compute_signal.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
#   GPT-5.6);
# - Anthropic Claude Code (Opus 5, Fable 5).
#
# Distributed under the MIT license.


# TODO: more examples, studies, and PMIDs in detailed help: e.g., Bressan et
# al., 2024 (PMID: 39211331).
# TODO: check the ordering and semantic paragraphs, particularly around
# '... [--csv_pseudo <csv>] [--typ_sig <str>] ...'.
usage=$(cat << EOM
Usage
-----
  execute_compute_signal.sh
    [--help] [--details] [--all_help]
    [--verbose] [--dry_run]
    [--env_nam <str>] [--threads <int>]
    [--mode <mode>] [--method <method>]
    (--csv_fil_in <csv> [--ref_fa <file>] | --csv_fil_A <csv> --csv_fil_B <csv> [--chr_siz <file>])
    --dir_out <dir> [--typ_out <format>] [--prefix <str>]
    [--siz_bin <int>] [--engine <engine>] [--siz_win <int>] [--csv_usr_frg <csv>] [--csv_scl_fct <csv>]
    [--csv_dep_min <csv>] [--csv_pseudo <csv>] [--typ_sig <str>] [--eps <num>] [--skip_00 <choice>] [--strict_bins] [--drp_nan] [--skp_pfx <csv>]
    [--report_only] [--no_report]
    [--track] [--dp <int>]
    [--dir_eo <dir>] [--nam_job <str>] [--max_job <int>] [--slurm] [--time <time>]

EOM
)

# shellcheck disable=SC2154
function help_execute_compute_signal() {
cat >&2 << EOM
${usage}

  Coordinate and automate the computation of signal tracks, ratio tracks, or fragment-coordinate files from BAM/CRAM or bedGraph input files. Supports multiple normalization strategies and runs computations in serial or parallel via GNU Parallel or Slurm.

  For more details on what this script can do, including notes and usage examples, run either of the following:
  '''bash
  bash path/to/execute_compute_signal.sh --details
  bash path/to/execute_compute_signal.sh --all_help
  '''

Parameters
----------
  -h, --help : flag
    Print this short help message and exit.

  -d, --details : flag
    Print full documentation with notes and examples, then exit.

  -ah, --all_help : flag
    Print both the short and full documentation, then exit.

  -v, --verbose : flag
    Run script in verbose mode.

  -dr, --dry, --dry_run : flag
    Run script in dry-run mode.

  -en, --env, --env_nam : str
    Conda environment to activate (default: '${env_nam}').

  -t, --thr, --threads : int
    Number of threads to use (default: ${threads}).

  -md, --mode : {'signal', 'ratio', 'coord'}
    Workflow mode: 'signal', 'ratio', or 'coord' (default: '${mode}').

    See '--details' for accepted spellings.

  -me, --method : {'unadj', 'frag', 'norm', 'count', 'cpm', 'linear', 'log2', 'linear_r', 'log2_r'}
    Workflow method. Signal or ratio computation subtype.

    Used with '--mode signal' or '--mode ratio' (default if '--mode signal': 'norm'; default if '--mode ratio': 'linear').

    For '--mode signal', the main choices are
      - 'unadj' (base-pair overlap),
      - 'frag' (base-pair overlap, divided by fragment length),
      - 'norm' (base-pair overlap, divided by fragment length and total fragments, so coverage integrates to 1),
      - 'count' (one whole count per bin a fragment touches), and
      - 'cpm' (those whole counts rescaled to sum to one million).

    For '--mode ratio', the main choices are
      - 'linear' ('file A / file B'),
      - 'log2' ['log2(file A / file B)'],
      - 'linear_r' ('file B / file A'), and
      - 'log2_r' ['log2(file B / file A)' = '-log2(file A / file B)'].

    See '--details' for accepted spellings and references.

  -ci, --csv_fil_in : list of file
    Comma-separated list of input file paths for coordinate-sorted BAM/CRAM files.

    Used with '--mode signal' or '--mode coord'.

  -rf, --ref_fa : file
    Reference FASTA file for CRAM input files (required if any '--csv_fil_in' element ends in '.cram').

    Used with '--mode signal' or '--mode coord'.

  -cs, --chr_siz : file
    Chromosome sizes file in UCSC-style TSV format.

    Used with '--mode ratio'.

    In '--mode ratio', '--chr_siz <file>' validates bedGraph interval bounds. Signal and coordinate modes, on the other hand, take chromosome sizes from the BAM/CRAM header.

  -cA, --csv_fil_A : list of file
    Comma-separated list of file A paths for coordinate-sorted numerator bedGraph files (e.g., IP signal).

    Used with '--mode ratio'.

  -cB, --csv_fil_B : list of file
    Comma-separated list of file B paths for coordinate-sorted denominator bedGraph files (e.g., input signal).

    Used with '--mode ratio'.

  -do, --dir_out : dir
    Output directory for generated files.

  -to, --typ_out : {'bedGraph', 'bedGraph.gz', 'bedgraph', 'bedgraph.gz', 'bdg', 'bdg.gz', 'bg', 'bg.gz', 'bed', 'bed.gz'}
    Output file format (default: '${typ_out}').

    For '--mode signal' or '--mode ratio', the typical choice is a bedGraph-style track (e.g., 'bedGraph.gz').

    For '--mode coord', the typical choice is a BED-like coordinate file (e.g., 'bed.gz').

    If an incompatible combination is requested, the value is coerced to a sensible default with a warning; see '--details' for the full list and mode-specific behavior.

  -px, --pfx, --prefix : str
    Custom prefix to prepend to output filenames.

  -sb, --siz_bin : int
    Bin size in base pairs for signal computation (default: 10).

    Used with '--mode signal'

  -eg, --engine : {'chrom', 'window'}
    Processing engine for signal computation (default: '${engine}').

    Both engines dispatch indexed fetch tasks and produce the same signal; they differ only in how fetch work is divided among threads.
      - 'chrom': one fetch task per chromosome. Task size tracks chromosome size, so the longest chromosomes dominate wall time.
      - 'window': each chromosome is split into fixed-size coordinate windows, with one fetch task per window. Task sizes are uniform, giving finer load balance across threads, but at the cost of more fetch calls. Window size is set by '--siz_win'.

    Tentative recommendation for S. cerevisiae datasets: keep 'chrom' as the general choice and the current best choice for CRAM input; try 'window' for large BAM inputs.

  -sw, --siz_win : int
    Window size in base pairs for the 'window' engine's indexed fetch tasks (default: ${siz_win}). Ignored by the 'chrom' engine.

    Used with '--mode signal'.

  -csf, --csv_scl_fct : list of structured string
    Comma-separated list of scaling factors or sentinels.

    Used with '--mode signal' or '--mode ratio'.

    For '--mode signal', each element must be 'NA' or a positive scalar float.

    For '--mode ratio', each element may be 'NA', a positive scalar float, or a positive 'A:B' spec, where A scales '--csv_fil_A' and B scales '--csv_fil_B' before ratio calculation.

    A factor multiplies the finished track verbatim, whatever '--method' produced it.

  -cuf, --csv_usr_frg : list of int
    Comma-separated list of fixed fragment-length values or sentinels.

    Used with '--mode signal' or '--mode coord'.

  -cdm, --csv_dep_min : list of number
    Comma-separated list of minimum-depth values or sentinels; here 'min' abbreviates minimum.

    Used with '--mode ratio'.

  -cps, --csv_pseudo : list of structured string
    Comma-separated list of pseudocount values as per-sample specs 'A[:B]', the element 'edger' to derive one from the counts beside each track, or sentinels.

    Used with '--mode ratio'.

  -ts, --typ_sig : str
    Signal type the input tracks carry, required by '--csv_pseudo edger'.

    Used with '--mode ratio'.

  -e, --eps : number
    Zero tolerance epsilon or sentinel used for ratio-mode zero checks.

    Used with '--mode ratio'.

  -s0, --skp_00, --skip_00 : {'pre_scale', 'post_scale'}
    Skip rows where both compared values are zero. Shared zero-zero skip mode or sentinel for ratio computation: 'pre_scale' or 'post_scale'.

    Used with '--mode ratio'.

  -stn, --strict_bins : flag
    Require strict bin compatibility. If '--mode ratio', require both input bedGraph files to have the same ordered '(chrom, start, end)' grid across all data rows.

  -dn, --drp_nan, --drop_nan : flag
    Drop non-finite ratio rows ('inf', '-inf', and 'nan') from the main ratio output.

    Used with '--mode ratio'.

  -sp, --skp_pfx : list of str
    Comma-separated list of header prefixes to skip. Shared comma-separated list of bedGraph header prefixes or sentinel to skip while parsing ratio-mode input files.

    Used with '--mode ratio'.

  -ro, --report_only : flag
    Write only the reports, with no signal track. By default, the fragment count 'N' / 'n_frg' and the fragment-bin overlap count 'L' / 'n_ovlp' are written beside each output track as '<track>.n_frg.txt' and '<track>.n_ovlp.txt', both of which are needed for pseudocount regularization per (or adapted from) edgeR.

    Used with '--mode signal'.

  -nr, --no_report : flag
    Suppress both counts. Cannot be combined with '--report_only'.

    Used with '--mode signal'.

  -tr, --trk, --track : flag
    Write a companion track file. If '--mode ratio', also write a companion bedGraph with all non-finite rows ('inf', '-inf', and 'nan') removed.

  -dp, --dp : int
    Maximum number of decimal places retained for finite emitted values (default: ${dp}).

    After rounding, non-informative trailing zeros are stripped.

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files in TXT format (default: '\${dir_out}/logs').

  -nj, --nam_job : str
    Job name. Prefix for job names (default depends on resolved '--mode' and '--method'; e.g., 'compute_signal_norm', 'compute_ratio_linear', or 'compute_coord').

  -mj, --max_job : int
    Maximum number of jobs to run concurrently (default: ${max_job}).

  -sl, --slurm : flag
    Submit jobs to the Slurm scheduler.

  -tm, --time : time
    Slurm job time limit in 'h:mm:ss' format (required if '--slurm'; default: '${time}').

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - basename
    - bash >= 4.4
    - dirname
    - parallel (when '--slurm' is not specified and multiple jobs are run)
    - python >= 3.11
    - rm (when '--slurm' is not specified and multiple jobs are run)
    - sbatch (when '--slurm' is specified)
    - tr

  - BAM/CRAM and bedGraph input files must be coordinate-sorted.
  - CRAM inputs require '--ref_fa'.
  - Input and output paths supplied to this wrapper interface must not contain spaces, commas, or semicolons.
  - See 'execute_compute_signal.sh --details' for more notes.

Examples
--------
  1. Compute normalized coverage from two BAM files locally.
    '''bash
    bash "\${dir_scr}/execute_compute_signal.sh" \\
        --mode signal \\
        --method norm \\
        --csv_fil_in IP1.bam,IP2.bam \\
        --dir_out tracks \\
        --typ_out bedGraph.gz \\
        --siz_bin 50
    '''

  2. Compute log2 IP/input ratios from bedGraph files.
    '''bash
    bash "\${dir_scr}/execute_compute_signal.sh" \\
        --mode ratio \\
        --method log2 \\
        --csv_fil_A IP1.bedGraph,IP2.bedGraph \\
        --csv_fil_B in1.bedGraph,in2.bedGraph \\
        --csv_pseudo 0.5,0.4 \\
        --dir_out ratios \\
        --typ_out bedGraph.gz \\
        --track
    '''
EOM
}


function detail_execute_compute_signal() {
    local mode="${1:-}"

    # Only print the top-level 'Usage' block if '--no-usage' is not invoked.
    if [[ ! "${mode}" =~ ^--no[_-]usage$ ]]; then
cat >&2 << EOM
${usage}

EOM
    fi

cat >&2 << EOM
  Driver script automating the computations of bedGraph signal or ratio tracks, or BED-like fragment coordinate files, from BAM/CRAM (for signal tracks or fragment coordinate files) or bedGraph (for ratio tracks) input files.

  Supports multiple signal normalization strategies, including the following:
    - unadjusted (raw) signal (i.e., base-pair overlap with no fragment-length or total-fragment adjustment)
    - fragment-length adjusted signal (Dickson et al., JBC 2020 [PMID: 32994221]; Dickson et al., Sci Rep 2023 [PMID: 37160995])
    - normalized coverage (Dickson et al., Sci Rep 2023 [PMID: 37160995])
    - siQ-ChIP IP efficiency (input-normalized ratio-based analyses; Dickson et al., JBC 2020 [PMID: 32994221]; Dickson et al., Sci Rep 2023 [PMID: 37160995])
    - log2-transformed input-normalized signal ratios (e.g., as described in "Data Analysis G" of Alavattam et al., Bio-protocol 2025 [PMID: 40364978])
    - spike-in-normalized signal (i.e., unmodified ChIP-Rx; Orlando et al., Cell Rep 2014 [PMID: 25437568])
    - input- and spike-in-normalized signal ratios (i.e., ChIP-Rx alpha ratio as described in "Data Analysis I" of Alavattam et al., Bio-protocol 2025 [PMID: 40364978])

  Supports serial job execution and parallel job execution via Slurm or GNU Parallel.

Parameters
----------
  -h, --help : flag
    Print a short help message and exit.

  -d, --details : flag
    Print this full documentation with notes and examples, then exit.

  -ah, --all_help : flag
    Print both the short and full documentation, then exit.

  -v, --verbose : flag
    Run script in verbose mode.

  -dr, --dry, --dry_run : flag
    Run script in dry-run mode.

  -en, --env, --env_nam : str
    Conda environment to activate (default: '${env_nam}').

  -t, --thr, --threads : int
    Number of threads to use (default: ${threads}).

  -md, --mode : {'signal', 'ratio', 'coord'}
    Workflow mode: 'signal', 'ratio', or 'coord' (default: '${mode}'). Available options:
      - 's', 'sig', 'signal':
        + Compute signal tracks directly from BAM/CRAM input files.
        + Supports unadjusted signal, fragment-length-adjusted signal, or normalized coverage.
        + See '--method' for calculation styles.
        + Use this option with appropriately computed scaling factors ('--csv_scl_fct') to compute unmodified ChIP-Rx (Orlando et al., Cell Rep 2014 [PMID: 25437568]) spike-in-normalized signal.
        + If 's' or 'sig' are supplied, variable 'mode' is set to "signal".

      - 'r', 'rat', 'ratio':
        + Compute IP/input signal ratios from bedGraph files.
        + Supports linear ratios, log2 ratios, and reciprocal variants of both.
        + Untransformed ratios are used in the modified ChIP-Rx spike-in and siQ-ChIP normalizations described in Alavattam et al., Bio-protocol 2025.
        + Use this option with '--method log2' to compute log2(IP/input) ratios.
        + If 'r' or 'rat' are supplied, variable 'mode' is set to "ratio".

      - 'c', 'coord', 'coordinates':
        + Instead of computing signal or ratio tracks, output fragment coordinates in BED-like format from BAM/CRAM input files.
        + Use this option to prepare input files for the original siQ-ChIP implementation ((Dickson et al., JBC 2020 [PMID: 32994221]; Dickson et al., Sci Rep 2023 [PMID: 37160995])).
        + For more details, see github.com/BradleyDickson/siQ-ChIP or github.com/kalavattam/siQ-ChIP.
        + This mode disables '--siz_bin' and '--method', and sets '--typ_out' to 'bed.gz' by default (or 'bed' if '--typ_out bed' is specified).
        + If 'c' or 'coordinates' are supplied, variable 'mode' is set to "coord".

  -me, --method : {'unadj', 'frag', 'norm', 'count', 'cpm', 'linear', 'log2', 'linear_r', 'log2_r'}
    Workflow method. Signal or ratio computation subtype used with '--mode signal' or '--mode ratio' (default if '--mode signal': norm; default if '--mode ratio': linear).
      - If '--mode signal', then the available options are
        + 'unadj':
          - Compute unadjusted signal (base-pair overlap with no fragment-length or total-fragment adjustment).

        + 'frag':
          - Adjust the above signal by fragment length.
          - For example, use this option to compute siQ-ChIP-scaled signal with the initial equation described in Dickson et al., JBC 2020 (PMID: 32994221) or Equation 5 in Dickson et al., Sci Rep 2023 (PMID: 37160995).

        + 'nc', 'norm':
          - Compute normalized coverage per Dickson et al., Sci Rep 2023, adjusting base-pair overlap signal by both fragment length and the total number of fragments so that the genome-wide coverage sums to 1.
            + That is, the coverage integrates to unity and can be interpreted as a probability distribution over the genome.
          - Internally, both of these values are standardized to 'method=norm'.

        + 'count':
          - Deposit one whole count in every bin a fragment touches, whether the fragment covers the whole bin or a single base.
          - The track's value column sums to the overlap count 'L' / 'n_ovlp' reported beside the track on the same run.

        + 'cpm':
          - Rescale the whole counts so the track sums to one million, following the edgeR construction.

      - If '--mode ratio', then the available options are
        + 'linear':
          - Compute the linear (i.e., non-log2) fil_A/fil_B ratio (e.g., IP/input): 'ratio = fil_A / fil_B'.

        + 'log2', 'l2':
          - Compute log2(fil_A/fil_B) ratio [e.g., log2(IP/input)]: 'ratio = log2(fil_A / fil_B)'.
          - Internally, both of these values are standardized to 'method=log2'.

        + 'linear_r':
          - Compute the reciprocal of the linear ratio: 'ratio = fil_B / fil_A = 1 / (fil_A / fil_B)'.

        + 'log2_r', 'l2_r':
          - Compute the reciprocal of the log2(fil_A/fil_B) ratio: 'ratio = log2(fil_B / fil_A) = -log2(fil_A / fil_B)'.
          - Internally, both of these values are standardized to 'method=log2_r'.

  -ci, --csv_fil_in : list of file
    Comma-separated list of input file paths for coordinate-sorted BAM/CRAM files.

    Required when '--mode signal' or '--mode coord'. Ignored for '--mode ratio'.

  -rf, --ref_fa : file
    Reference FASTA file for CRAM input files.

    Required when any '--csv_fil_in' element ends in '.cram'.

    Used with '--mode signal' or '--mode coord'.

  -cs, --chr_siz : file
    Chromosome sizes file in UCSC-style TSV format, with chromosome names in the first column and positive integer sizes in the second.

    Used with '--mode ratio'.

    In '--mode ratio', it validates that every bedGraph interval falls within its chromosome's bounds. This validation is independent of '--strict_bins', which compares the two inputs to each other rather than to a reference.

    On the other hand, '--mode signal' and '--mode coord' do not accept '--chr_siz <file>': they take chromosome sizes from the BAM/CRAM header, which is authoritative for the alignments being read.

  -cA, --csv_fil_A : list of file
    Comma-separated list of file A paths for coordinate-sorted numerator bedGraph files (e.g., ChIP IP signal tracks).

    Used with '--mode ratio'.

    The list order must match that of '--csv_fil_B' files.

  -cB, --csv_fil_B : list of file
    Comma-separated list of file B paths for coordinate-sorted denominator bedGraph files (e.g., input signal tracks).

    Used with '--mode ratio'.

    The list order must match that of '--csv_fil_A' files.

  -do, --dir_out : dir
    Output directory for generated files:
      - Signal tracks if '--mode signal'.
      - Ratio tracks if '--mode ratio'.
      - BED-like files of fragment coordinates if '--mode coord'.

  -to, --typ_out : {'bedGraph', 'bedGraph.gz', 'bedgraph', 'bedgraph.gz', 'bdg', 'bdg.gz', 'bg', 'bg.gz', 'bed', 'bed.gz'}
    Output file format for signal track output files (default: '${typ_out}'). Available options:
      - 'bedGraph', 'bedgraph', 'bdg', 'bg':
        + Signal/ratio in bedGraph format.
        + Intended for '--mode signal' or '--mode ratio'; with '--mode coord' these values are accepted but coerced to 'bed.gz' (see Notes).

      - 'bedGraph.gz', 'bedgraph.gz', 'bdg.gz', 'bg.gz':
        + Signal/ratio in gzip-compressed bedGraph format.
        + Intended for '--mode signal' or '--mode ratio'; with '--mode coord' these values are accepted but coerced to 'bed.gz' (see Notes).

      - 'bed', 'bed.gz':
        + BED-like format for fragment coordinates instead of signal.
        + Intended for '--mode coord'; with '--mode signal' or '--mode ratio' these values are accepted but coerced to 'bedGraph.gz' (see Notes).

  -px, --pfx, --prefix : str
    Custom prefix to prepend to output filenames.

    When '--mode signal' or '--mode coord':
      - If not specified, no prefix is added.
      - If specified, the prefix is prepended to the base filename, and any leading 'IP_' or 'in_' in that base filename is stripped before applying it.
      - Note: If you use the same '--prefix' for both IP and input BAMs/CRAMs in the same output directory, stripping 'IP_'/ 'in_' may cause filename collisions.
        + In that case, prefer distinct prefixes (e.g., 'IP', 'in').

    When '--mode ratio':
      - If not specified, a default prefix is automatically constructed based on '--method' and '--csv_scl_fct'; for example:
        + 'rat' (default)
        + 'log2_rat' (if '--method log2')
        + 'recip_rat' (if '--method linear_r')
        + 'log2_recip_rat' (if '--method log2_r')
        + 'scl_rat' (if '--method linear' and '--csv_scl_fct' is supplied)
        + 'scl_log2_rat' (if '--method log2' and '--csv_scl_fct' is supplied)
        + 'scl_recip_rat' (if '--method linear_r' and '--csv_scl_fct' is supplied)
        + 'scl_log2_recip_rat' (if '--method log2_r' and '--csv_scl_fct' is supplied)
      - If specified, the custom prefix replaces the default.
      - Whether specified or not, any leading 'IP_' string in the base name is stripped before the prefix.

  -sb, --siz_bin : int
    Bin size in base pairs for signal computation (default: 10).

    Used with '--mode signal'.

  -eg, --engine : {'chrom', 'window'}
    Processing engine for signal computation (default: '${engine}').

    Both engines dispatch indexed fetch tasks and produce the same signal; they differ only in how fetch work is divided among threads. Available options:
      - 'chrom':
        + One fetch task per chromosome.
        + Task size tracks chromosome size, so the longest chromosomes dominate wall time.
        + The general choice, and the current best choice for CRAM input.

      - 'window':
        + Each chromosome is split into fixed-size coordinate windows, with one fetch task per window.
        + Task sizes are uniform, giving finer load balance across threads, at the cost of more fetch calls.
        + Worth trying for large BAM inputs, particularly with '--threads' above one.
        + Window size is set by '--siz_win'.

    Used with '--mode signal'; ignored otherwise.

  -sw, --siz_win : int
    Window size in base pairs for the 'window' engine's indexed fetch tasks (default: ${siz_win}).

    Smaller windows give finer load balance across threads at the cost of more fetch calls; larger windows do the reverse. Ignored by the 'chrom' engine, which uses whole chromosomes as its unit of work.

    Used with '--mode signal'; ignored otherwise.

  -csf, --csv_scl_fct : list of structured string
    Comma-separated list of scaling factors or sentinels to apply to signal or ratio values.

    Used with '--mode signal' or '--mode ratio'; ignored otherwise.

    List size must match the number of input files via '--csv_fil_in' or '--csv_fil_A'/'--csv_fil_B'.

    For '--mode signal', each non-sentinel element must be a positive scalar float.

    For '--mode ratio', each non-sentinel element may be either:
      - 'A'    Scale file A by A and file B by 1.0.
      - 'A:B'  Scale file A by A and file B by B.

    A factor multiplies the finished track verbatim, whatever '--method' produced it.

    A siQ-ChIP alpha is derived against a ratio of signal types: equations 5 and 6 for 'frag', 5nd and 6nd for 'norm'.

    A spike-in alpha is likewise derived against a ratio of signal types: either 'count' or 'unadj', since the per-bin deposition difference largely cancels in the division. Earlier spike-in work used 'count', so prefer it where continuity with those results matters.

  -cuf, --csv_usr_frg : list of int
    Comma-separated list of fixed fragment-length values or sentinels to use instead of read lengths (single-end alignments) or template lengths (paired-end alignments).

    Used with either '--mode signal' or '--mode coord'; ignored otherwise.

    List size must match the number of input files via '--csv_fil_in'.

  -cdm, --csv_dep_min : list of number
    Comma-separated list of minimum-depth values or sentinels used to avoid extreme division operations; here 'min' abbreviates minimum.

    Used with '--mode ratio'; ignored otherwise.

    List size must match the number of input files via '--csv_fil_A'/'--csv_fil_B'.

    Although allowed, using '--csv_dep_min' together with '--csv_pseudo' is usually harder to interpret, since both stabilize low-depth ratio behavior in different ways.

  -cps, --csv_pseudo : list of structured string
    Comma-separated list of pseudocount values or sentinels used during ratio computation.

    Each element is one of:
      - 'A' or 'A:B', a literal pseudocount for the pair.
      - 'edger', which is a request rather than a value: the counts a signal run wrote beside each track are read back and passed to 'compute_pseudo --method edger', whose result is used for that pair. Needs the following:
        + '--typ_sig <spec>' and
        + both signal bedGraph tracks must have their '<track>.n_frg.txt' and '<track>.n_ovlp.txt' beside them.
      - 'NA', no pseudocount for that pair.

    Elements may be mixed, so one pair can be derived while another takes a literal.

    Used with '--mode ratio'; ignored otherwise.

  -ts, --typ_sig : str
    Signal type the input tracks carry, required by '--csv_pseudo edger'.

    Used with '--mode ratio'; ignored otherwise.

    The '--mode ratio --typ_sig <spec>' value is the '--mode signal --method <spec>' that produced the track ('--method' chooses what to do; '--typ_sig' declares what an existing input is). Rejected with '--mode signal', where '--method' already says what is being written.

    List size must match the number of input files via '--csv_fil_A'/'--csv_fil_B'.

    Each non-sentinel element may be either:
      - 'A': add pseudocount A symmetrically.
      - 'A:B': add pseudocount A to file A and B to file B.

    Although allowed, using '--csv_pseudo' together with '--csv_dep_min' is usually harder to interpret, since both stabilize low-depth ratio behavior in different ways.

  -e, --eps : number
    Zero tolerance epsilon or sentinel used for ratio-mode zero checks.

    Used with '--mode ratio'; ignored otherwise.

    Non-sentinel values must be non-negative floats.

  -s0, --skp_00, --skip_00 : {'pre_scale', 'post_scale'}
    Skip rows where both compared values are zero. Shared zero-zero skip mode or sentinel for ratio computation.

    Used with '--mode ratio'; ignored otherwise.

    Non-sentinel values must be one of 'pre_scale' or 'post_scale'.

  -stn, --strict_bins : flag
    Require strict bin compatibility between the two ratio inputs.

    With this flag, both input bedGraph files must share the same ordered '(chrom, start, end)' grid across all data rows. Without it, only the first few paired rows are checked for equal bin width, which catches gross mismatches but not divergence later in the file.

    Used with '--mode ratio'; ignored otherwise.

  -dn, --drp_nan, --drop_nan : flag
    Drop non-finite values from the main ratio output.

    Used only with '--mode ratio'; ignored otherwise.

    If set, rows yielding 'inf', '-inf', or 'nan' are omitted from the main ratio output.

  -sp, --skp_pfx : list of str
    Comma-separated list of header prefixes to skip. Shared comma-separated list of bedGraph header prefixes or sentinel to skip.

    Used with '--mode ratio'; ignored otherwise.

    Passed through to 'submit_compute_signal.sh' and then to 'compute_signal_ratio.py'.

  -ro, --report_only : flag
    Write the two per-sample counts without writing a signal track.

    By default, the following two per-sample counts are written beside each output track, whether or not this flag is given:
      - 'N' / 'n_frg', the fragment count, written as '<track>.n_frg.txt':
        + Counted from the same iterator that builds the track, under the same '--csv_usr_frg' and alignment-filter settings.
        + Thus, it's the number a '--method norm' run divides by, not an estimate of it.

      - 'L' / 'n_ovlp', the fragment-bin overlap count, written as '<track>.n_ovlp.txt':
        + The total number of bins the counted fragments span / overlap, counting a fragment once per bin it touches.
        + Not the summed base pairs an unadjusted track reports: the bin count is what puts 'k = L / N' in bins, which is the unit needed to compute per-bin pseudocounts.
        + In addition to alignment-filter settings, it depends on '--siz_bin' and, if supplied, '--csv_usr_frg', since fragment extension changes how far each fragment reaches.

    Both names follow the track's own derived name, so under '--prefix', '.sample.bedGraph' yields 'run_1.sample.n_frg.txt' and 'run_1.sample.n_ovlp.txt'. The pair is what 'compute_pseudo --method edger' consumes.

    Counting happens before the output branch, so the reported values are identical to those a track-writing run would produce on the same input. Use this to obtain prior counts without outputting a bedGraph track.

    Cannot be combined with '--no_report', as doing so would write neither a track nor a report. Such a run is rejected.

    Used with '--mode signal'; ignored otherwise.

  -nr, --no_report : flag
    Suppress the two per-sample counts described under '--report_only' above, which are otherwise written beside each output track.

    They are the inputs 'compute_pseudo --method edger' consumes, so suppressing them leaves a later ratio run without the numbers it needs for pseudocount regularization per (or adapted from) edgeR.

    Used with '--mode signal'; ignored otherwise.

  -tr, --trk, --track : flag
    Write a companion track file. If '--mode ratio', also write a companion bedGraph with all non-finite rows ('inf', '-inf', and 'nan') removed.

    The new file will include '.track' before the extension.

    This cleaned version is ideal for visualization in genome browsers such as IGV, avoiding issues caused by 'inf', '-inf' or 'nan' values.

  -dp, --dp : int
    Maximum number of decimal places retained for finite emitted values (default: ${dp}).

    After rounding, non-informative trailing zeros are stripped.

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files in TXT format (default: '\${dir_out}/logs').

  -nj, --nam_job : str
    Job name. Prefix for job names (default depends on resolved '--mode' and '--method'; e.g., 'compute_signal_norm', 'compute_ratio_log2', or 'compute_coord').

  -mj, --max_job : int
    Maximum number of jobs to run concurrently (default: ${max_job}).

    With '--slurm': maximum number of Slurm array tasks allowed to run concurrently.

    Without '--slurm':
      - If the resolved number of parallel jobs is greater than 1, jobs run in parallel via GNU Parallel.
      - If the resolved number of parallel jobs is 1, jobs run serially.

  -sl, --slurm : flag
    Submit jobs to the Slurm scheduler.

    If '--slurm' is not specified, this script uses '--threads', '--max_job', and the detected CPU core count (via 'set_params_parallel') to decide whether to run jobs in parallel with GNU Parallel or in serial.
      - In this non-Slurm path, '--threads' is treated as a total CPU/thread budget for the local machine, and 'set_params_parallel' converts that into threads per job plus number of parallel jobs.
      - If the resolved number of parallel jobs is greater than 1, jobs are run with GNU Parallel.
      - If the resolved number of parallel jobs is 1, jobs are run serially.

  -tm, --time : time
    Slurm job time limit. The length of time, in 'h:mm:ss' format, for the Slurm job (required if '--slurm' is specified, ignored if not; default: '${time}').

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - basename
    - bash >= 4.4
    - conda (when the requested environment is not active)
    - dirname
    - parallel (when '--slurm' is not specified and multiple jobs are run)
    - python >= 3.11
    - rm (when '--slurm' is not specified and multiple jobs are run)
    - sbatch (when '--slurm' is specified)
    - tr

  - BAM/CRAM and bedGraph input files must be coordinate-sorted.
  - CRAM inputs require '--ref_fa'.
  - Input and output paths supplied to this wrapper interface must not contain spaces, commas, or semicolons.
    + Commas are used internally as list delimiters.
    + Semicolons are also considered unsafe in this wrapper workflow.
    + Spaces in paths are not supported by the current list-serialization and reconstruction logic.
  - Use consistent file ordering between IP and input files.
  - '--typ_out' must be compatible with the selected '--mode'.
    + With '--mode signal' or '--mode ratio', bedGraph-style values ('bedGraph', 'bedgraph', 'bdg', and 'bg', and their '.gz' variants) are allowed; 'bed'/'bed.gz' are accepted but are automatically converted to 'bedGraph.gz' with a warning.
    + With '--mode coord', 'bed'/'bed.gz' are allowed; bedGraph-style values are accepted but are automatically converted to 'bed.gz' with a warning.
  - Output file path. Output filenames are derived from BAM/CRAM or bedGraph input files and the value associated with '--typ_out'.
  - For bedGraph-style output, '--dp' sets the maximum number of decimal places retained for finite emitted values; after rounding, non-informative trailing zeros and any trailing decimal point are stripped.
  - BED-like files of fragment coordinates are, e.g., used as input to the original siQ-ChIP implementation (Dickson et al., JBC 2020 [PMID: 32994221]; Dickson et al., Sci Rep 2023 [PMID: 37160995]).
  - Job execution mode (serial, GNU Parallel, or Slurm array) is chosen automatically from '--slurm', '--threads', and '--max_job':
    + If '--slurm' is specified:
      - Jobs are submitted as a Slurm array.
      - '--max_job' (after adjustment) sets the maximum number of array tasks running concurrently.
      - '--threads' controls '--cpus-per-task' for each array element.
    + If '--slurm' is not specified:
      - Helper function 'set_params_parallel' uses '--threads', '--max_job', and the detected CPU core count to determine a safe combination of threads per job and the number of parallel jobs.
      - If the resulting number of parallel jobs is greater than 1, commands are written to a configuration file and executed with GNU Parallel.
      - If the resulting number of parallel jobs is 1, all jobs are run serially in a single Bash process (neither GNU Parallel nor Slurm).

Examples
--------
  1. Compute normalized coverage using GNU Parallel.
    '''bash
    bash "\${HOME}/bin/execute_compute_signal.sh" \\
        --threads 8 \\
        --mode "signal" \\
        --method "norm" \\
        --csv_fil_in "\${HOME}/project/samples/sample_1.bam,\${HOME}/project/samples/sample_2.bam" \\
        --dir_out "\${HOME}/project/tracks" \\
        --typ_out "bedGraph.gz" \\
        --siz_bin 50 \\
        --dir_eo "\${HOME}/project/logs" \\
        --nam_job "norm_sig"
    '''

  2. Compute log2 IP/input ratios from bedGraph files in serial with user-supplied pseudocounts.
    '''bash
    bash "\${HOME}/bin/execute_compute_signal.sh" \\
        --threads 1 \\
        --mode "ratio" \\
        --method "log2" \\
        --csv_fil_A "\${HOME}/project/norm/IP_1.bedGraph,\${HOME}/project/norm/IP_2.bedGraph" \\
        --csv_fil_B "\${HOME}/project/norm/in_1.bedGraph,\${HOME}/project/norm/in_2.bedGraph" \\
        --csv_pseudo "0.5,0.4" \\
        --dir_out "\${HOME}/project/ratios" \\
        --typ_out "bedGraph.gz"
    '''

  3. Compute normalized coverage without the per-sample counts of 'N' / 'n_frg' and 'L' / 'n_ovlp'.
    '''bash
    bash "\${HOME}/bin/execute_compute_signal.sh" \\
        --threads 8 \\
        --mode "signal" \\
        --method "norm" \\
        --csv_fil_in "\${HOME}/project/samples/sample_1.bam,\${HOME}/project/samples/sample_2.bam" \\
        --dir_out "\${HOME}/project/tracks" \\
        --typ_out "bedGraph.gz" \\
        --siz_bin 50 \\
        --no_report \\
        --dir_eo "\${HOME}/project/logs" \\
        --nam_job "norm_sig"
    '''

    Writes 'sample_1.bedGraph.gz' and 'sample_2.bedGraph.gz' and nothing else. Example 1, which omits '--no_report', writes 'sample_1.n_frg.txt' and 'sample_1.n_ovlp.txt' beside each track. Suppress them only where no pseudocount is wanted later: 'compute_pseudo --method edger' consumes these two counts.

  4. Write only the fragment and overlap counts, with no signal track.
    '''bash
    bash "\${HOME}/bin/execute_compute_signal.sh" \\
        --threads 8 \\
        --mode "signal" \\
        --csv_fil_in "\${HOME}/project/samples/sample_1.bam,\${HOME}/project/samples/sample_2.bam" \\
        --dir_out "\${HOME}/project/counts" \\
        --siz_bin 50 \\
        --report_only \\
        --dir_eo "\${HOME}/project/logs" \\
        --nam_job "counts"
    '''

    Each report is named from the track that would have been written, so this writes 'sample_1.n_frg.txt' and 'sample_1.n_ovlp.txt' to '--dir_out', and the same pair for 'sample_2'. 'L' / 'n_ovlp' counts fragment-bin overlaps, so it depends on '--siz_bin'.

  5. Compute counts-per-million signal from whole-count deposition.
    '''bash
    bash "\${HOME}/bin/execute_compute_signal.sh" \\
        --threads 8 \\
        --mode "signal" \\
        --method "cpm" \\
        --csv_fil_in "\${HOME}/project/samples/sample_1.bam,\${HOME}/project/samples/sample_2.bam" \\
        --dir_out "\${HOME}/project/tracks" \\
        --typ_out "bedGraph.gz" \\
        --siz_bin 20 \\
        --prefix "cpm" \\
        --dir_eo "\${HOME}/project/logs" \\
        --nam_job "cpm_sig"
    '''

    Each fragment deposits one whole count per touched bin, and the track is rescaled to sum to one million.

    The divisor 'L' / 'n_ovlp' is written beside it as 'cpm.sample_1.n_ovlp.txt', following the track's prefixed name, so a CPM run can share an output directory with the normalized run of example 1.

  6. Compute log2 ratios and write a browser-ready companion track.
    '''bash
    bash "\${HOME}/bin/execute_compute_signal.sh" \\
        --threads 4 \\
        --mode "ratio" \\
        --method "log2" \\
        --csv_fil_A "\${HOME}/project/norm/IP_1.bedGraph" \\
        --csv_fil_B "\${HOME}/project/norm/in_1.bedGraph" \\
        --dir_out "\${HOME}/project/ratios" \\
        --typ_out "bedGraph.gz" \\
        --track
    '''

    Writes the ratio track and, beside it, a companion carrying '.track' before the extension with 'inf', '-inf', and 'nan' rows removed. Load that copy in a genome browser, where non-finite values otherwise cause rendering trouble.

  7. Derive an edgeR-styled pseudocount from the counts the signal run wrote.
    '''bash
    bash "\${HOME}/bin/execute_compute_signal.sh" \\
        --threads 4 \\
        --mode "ratio" \\
        --method "log2" \\
        --csv_fil_A "\${HOME}/project/norm/IP_1.bedGraph" \\
        --csv_fil_B "\${HOME}/project/norm/in_1.bedGraph" \\
        --dir_out "\${HOME}/project/ratios" \\
        --typ_out "bedGraph.gz" \\
        --typ_sig "norm" \\
        --csv_pseudo "edger"
    '''

    Reads 'IP_1.n_frg.txt', 'IP_1.n_ovlp.txt' and the matching pair beside 'in_1.bedGraph', passes all four to 'compute_pseudo --method edger', and applies the result to the IP-input pair.

    The signal runs that wrote those tracks wrote the counts beside them, so nothing extra is needed unless they were given '--no_report'. '--typ_sig norm' says the tracks carry normalized coverage, which is what '--method norm' produced.
EOM
}
