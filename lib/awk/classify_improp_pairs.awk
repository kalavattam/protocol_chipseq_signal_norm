# -*- coding: utf-8 -*-
#
# Script: classify_improp_pairs.awk
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


# Usage
# -----
#   gawk -f classify_improp_pairs.awk \
#       -v aligner=<aligner> -v max_frag=<int> -v min_frag=<int> -v mapq=<int>
#
#   Count the improper pairs of an '--unal_mode separate' run by category, for
#   the 'improp' rows of '<stem>.diag.tsv.gz' (Step 5 of 'align_fastqs()').
#   Read the records of the improper file from stdin as SAM text without a
#   header, in any order.
#
# Parameters
# ----------
#   aligner : str
#     Aligner that wrote the records: 'bowtie2', 'bwa', or 'bwa-mem2'.
#
#   max_frag : int
#     Bowtie 2's '-X', the longest fragment of a proper pair.
#
#   min_frag : int
#     Bowtie 2's '-I', the shortest fragment of a proper pair.
#
#   mapq : int
#     MAPQ threshold ('--mapq') below which a primary record makes its
#     template low-MAPQ.
#
# Returns
# -------
#   Exits 0 on success. Prints an error and exits 1 on a template without a
#   primary first-mate record or on a primary record without 'MC'.
#
# Output
# ------
#   One tab-separated row per category to stdout, in the order of the rows in
#   'docs/design/align_fastqs_separate.md', with six columns: 'improp', the
#   category, its templates and records, and its low-MAPQ templates and
#   records. 'too_far' and 'too_close' are NA for the BWA aligners, since only
#   Bowtie 2's fragment bounds are known here, and 'failed_est' is NA for
#   every aligner.
#
# Notes
# -----
#   Every template in the input is a pair with both mates aligned, and
#   'samtools fixmate -m' has given its primary records the mate CIGAR ('MC').
#   The design page defines each category and how a pair is judged.
#
#   Runtime requirements:
#     - gawk
#
# Examples
# --------
#   1. Count the improper pairs of a Bowtie 2 run at the default '-X' and
#      '--mapq'.
#     '''bash
#     samtools view sample.improp.bam \
#         | gawk -f lib/awk/classify_improp_pairs.awk \
#             -v aligner=bowtie2 -v max_frag=1000 -v min_frag=0 -v mapq=1
#     '''


BEGIN {
    FS = "\t"
    OFS = "\t"

    # The categories that need fragment bounds, which only Bowtie 2 gives.
    bt2 = (aligner == "bowtie2")
    max_frag += 0
    min_frag += 0
    mapq += 0

    # Bowtie 2 judges a fragment shorter than 1 base too close whatever '-I'
    # says ('pe.cpp', 'PairedEndPolicy::peClassifyPair').
    if (min_frag < 1) min_frag = 1

    n_cat = split( \
        "all diff_ref orient_ff orient_rr orient_rf dovetailed too_far " \
        "too_close failed_est other unjudged duplicate_all", \
        arr_cat, " " \
    )
    na["failed_est"] = 1
    if (!bt2) {
        na["too_far"] = 1
        na["too_close"] = 1
    }
}


# Return the reference bases a CIGAR string covers: its 'M', 'D', 'N', '=', and
# 'X' operations.
function get_ref_len(cig,    len, n, op) {
    len = 0
    while (match(cig, /^[0-9]+[MIDNSHP=X]/)) {
        n = substr(cig, 1, RLENGTH - 1) + 0
        op = substr(cig, RLENGTH, 1)
        if (op ~ /[MDN=X]/) len += n
        cig = substr(cig, RLENGTH + 1)
    }

    return len
}


# Return the value of a record's SAM tag, or "" if it has none.
function get_tag(tag,    fld) {
    for (fld = 12; fld <= NF; fld++) {
        if (substr($fld, 1, 5) == tag ":Z:") return substr($fld, 6)
    }

    return ""
}


# Return the first reason that applies to a first mate's primary record.
function get_reason(    mc, s1, e1, s2, e2, rev1, rev2, fs, fe, rs, re, \
    lo, hi, span, mx) {
    if ($7 != "=" && $7 != $3) return "diff_ref"

    rev1 = int($2 / 16) % 2
    rev2 = int($2 / 32) % 2
    if (!rev1 && !rev2) return "orient_ff"
    if (rev1 && rev2) return "orient_rr"

    mc = get_tag("MC")
    if (mc == "") {
        printf "error(classify_improp_pairs.awk): '%s' has no 'MC' tag.\n", \
            $1 > "/dev/stderr"
        bad = 1
        exit 1
    }

    s1 = $4 + 0
    e1 = s1 + get_ref_len($6) - 1
    s2 = $8 + 0
    e2 = s2 + get_ref_len(mc) - 1

    # The forward mate ('fs' to 'fe') and the reverse mate ('rs' to 're').
    if (rev1) {
        fs = s2; fe = e2; rs = s1; re = e1
    } else {
        fs = s1; fe = e1; rs = s2; re = e2
    }

    # Facing outward: the reverse mate ends before the forward mate starts.
    if (re < fs) return "orient_rf"

    # One mate extends past the other: the reverse mate starts before the
    # forward mate, or the forward mate ends after it. Since the reverse mate
    # does not end before the forward mate starts, the two overlap.
    if (rs < fs || fe > re) return "dovetailed"

    if (bt2) {
        lo = (fs < rs) ? fs : rs
        hi = (fe > re) ? fe : re
        span = hi - lo + 1

        # Bowtie 2 raises '-X' to a mate's length when the mate is longer.
        mx = max_frag
        if (e1 - s1 + 1 > mx) mx = e1 - s1 + 1
        if (e2 - s2 + 1 > mx) mx = e2 - s2 + 1

        if (span > mx) return "too_far"
        if (span < min_frag) return "too_close"
    }

    return "other"
}


{
    nam = $1
    flg = $2 + 0
    n_rec[nam]++

    # 'markdup' skips unaligned (0x4), secondary (0x100), and supplementary
    # (0x800) records.
    if (int(flg / 4) % 2 || int(flg / 256) % 2 || int(flg / 2048) % 2) {
        n_unj[nam]++
    } else {
        # A primary record below the threshold makes the template low-MAPQ.
        if ($5 + 0 < mapq) low[nam] = 1
        if (int(flg / 64) % 2) rsn[nam] = get_reason()
    }

    if (int(flg / 1024) % 2) n_dup[nam]++
}


END {
    if (bad) exit 1

    for (nam in n_rec) {
        if (!(nam in rsn)) {
            printf "error(classify_improp_pairs.awk): '%s' has no primary " \
                "first-mate record.\n", nam > "/dev/stderr"
            exit 1
        }

        add_count("all", nam, n_rec[nam])
        add_count(rsn[nam], nam, n_rec[nam])
        if (n_unj[nam]) add_count("unjudged", nam, n_unj[nam])
        if (n_dup[nam]) add_count("duplicate_all", nam, n_dup[nam])
    }

    for (i = 1; i <= n_cat; i++) {
        cat = arr_cat[i]
        if (cat in na) {
            print "improp", cat, "NA", "NA", "NA", "NA"
        } else {
            print "improp", cat, tpl[cat] + 0, rec[cat] + 0, \
                tpl_low[cat] + 0, rec_low[cat] + 0
        }
    }
}


# Add one template and its records in a category to the counts.
function add_count(category, qname, n) {
    tpl[category]++
    rec[category] += n
    if (qname in low) {
        tpl_low[category]++
        rec_low[category] += n
    }
}
