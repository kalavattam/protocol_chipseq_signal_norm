#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_help_style.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="help style"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


#  Scan a caller-supplied file list for a fixed pattern
function scan_fixed_files() {
    local lbl="${1:-}"
    local ptn="${2:-}"
    local sev="${3:-fail}"

    shift 3

    local fil_scn=( "$@" )
    local found=0
    local file line

    for file in "${fil_scn[@]}"; do
        while IFS= read -r line; do
            found=1
            if [[ "${sev}" == "warn" ]]; then
                record_warn "${lbl}: $(print_relpath "${file}"):${line}"
            else
                record_fail "${lbl}: $(print_relpath "${file}"):${line}"
            fi
        done < <(grep -n -- "${ptn}" "${file}" || true)
    done

    if (( found == 0 )); then
        record_pass "no ${lbl} findings"
    fi
}


#  Scan the default style-check file set for a fixed pattern
function scan_fixed() {
    local lbl="${1:-}"
    local ptn="${2:-}"
    local sev="${3:-fail}"

    scan_fixed_files "${lbl}" "${ptn}" "${sev}" "${files[@]}"
}


#  Fail when retired reference-FASTA aliases are used in shell/docs/test text
function scan_retired_ref_fa_aliases() {
    local found=0
    local file line
    local -a fil_scn=()

    while IFS= read -r file; do
        fil_scn+=( "${file}" )
    done < <(
        {
            find \
                "${ROOT_REPO}/bin" \
                -type f \
                \( -name '*.sh' -o -name '*.py' \) \
                -print

            find \
                "${ROOT_REPO}/install/scripts" \
                -type f \
                -name '*.sh' \
                -print

            find \
                "${ROOT_REPO}/tests" \
                -type f \
                -name '*.sh' \
                -print

            find \
                "${ROOT_REPO}/docs/standards" \
                -type f \
                -name '*.md' \
                -print
        } | sort
    )

    for file in "${fil_scn[@]}"; do
        case "$(basename "${file}")" in
            test_help_style.sh|test_python_style.sh)
                continue
                ;;
        esac
        while IFS= read -r line; do
            found=1
            record_fail "retired ref_fa alias: $(print_relpath "${file}"):${line}"
        done < <(
            perl -ne '
                next if /retired ref_fa alias/;
                next if /--ref\(?!/;
                if (
                    /--ref(?![A-Za-z0-9_-]|\[)/ ||
                    /-r,\s+--ref_fa\b/ ||
                    /-r\|--ref\[_-\]fa\)/ ||
                    /"-r",\s*"--ref_fa"/
                ) {
                    print "$.:$_";
                }
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no retired ref_fa alias findings"
    fi
}


#  Warn on hard tabs in user-facing help text
function scan_tabs() {
    local found=0
    local file line

    for file in "${files[@]}"; do
        while IFS= read -r line; do
            found=1
            record_warn "hard tab: $(print_relpath "${file}"):${line}"
        done < <(grep -n $'\t' "${file}" || true)
    done

    if (( found == 0 )); then
        record_pass "no hard tab findings"
    fi
}


#  Fail when stale precision aliases appear as the canonical Usage spelling
function scan_usage_legacy_precision() {
    local found=0
    local file line
    local ptn="--r""nd"

    for file in "${files[@]}"; do
        while IFS= read -r line; do
            found=1
            record_fail \
                "stale precision alias in Usage block:" \
                "$(print_relpath "${file}"):${line}"
        done < <(
            awk -v ptn="${ptn}" '
                /^Usage$/ { in_usage = 1 }
                in_usage && /^$/ { in_usage = 0 }
                in_usage && index($0, ptn) { print FNR ":" $0 }
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no stale precision alias Usage findings"
    fi
}


#  Fail when a fixed pattern appears outside parser case-pattern lines
function scan_fixed_nonparser() {
    local lbl="${1:-}"
    local ptn="${2:-}"

    local found=0
    local file line

    for file in "${files[@]}"; do
        while IFS= read -r line; do
            found=1
            record_fail "${lbl}: $(print_relpath "${file}"):${line}"
        done < <(
            awk -v ptn="${ptn}" '
                index($0, ptn) && $0 !~ /\|/ { print FNR ":" $0 }
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no ${lbl} findings"
    fi
}


#  Fail when legacy colon-terminated top-level help headings are used
function scan_legacy_colon_headings() {
    local found=0
    local file line

    for file in "${files[@]}"; do
        while IFS= read -r line; do
            found=1
            record_fail \
                "legacy colon help heading: $(print_relpath "${file}"):${line}"
        done < <(
            awk '
                BEGIN {
                    heading = "Usage|Description|Argument|Arguments"
                    heading = heading "|Keyword argument|Keyword arguments"
                    heading = heading "|Positional argument|Positional arguments"
                    heading = heading "|Parameters"
                    heading = heading "|Dependencies|Returns|Output|Notes"
                    heading = heading "|References|See Also"
                    heading = heading "|Examples|Example|Testing"
                    heading = heading "|Expected global|Expected global variable|Expected global variables"
                    heading = heading "|Expected globals"
                    heading = heading "|Generated global|Generated global variable|Generated global variables"
                    heading = heading "|Generated globals|#TODO"
                }
                $0 ~ "^(" heading "):$" {
                    print FNR ":" $0
                }
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no legacy colon help heading findings"
    fi
}


#  Fail when major underlined help headings have malformed spacing
function scan_major_section_format() {
    local found=0
    local file line

    for file in "${files[@]}"; do
        while IFS= read -r line; do
            found=1
            record_fail \
                "major section format: $(print_relpath "${file}"):${line}"
        done < <(
            awk '
                BEGIN {
                    heading = "Usage|Parameters|Returns|Notes"
                    heading = heading "|References|See Also|Examples"
                    heading = heading "|Expected globals"
                    heading = heading "|Generated globals|#TODO"
                }
                function is_heading(line) {
                    return line ~ "^(" heading ")$"
                }
                function is_underline(line) {
                    return line ~ /^-{3,}$/
                }
                function missing_preceding_blank(line, prev_line, ok) {
                    ok = (line != "Usage")
                    ok = ok && (FNR > 2)
                    ok = ok && (prev_line != "")
                    ok = ok && (prev_line !~ /<< ?'\''?EOM'\''?/)
                    return ok
                }
                is_heading($0) {
                    if (getline underline <= 0 || ! is_underline(underline)) {
                        print FNR ":" $0 " [missing underline]"
                    } else if (missing_preceding_blank($0, prev)) {
                        print (FNR - 1) ":" $0 " [missing preceding blank line]"
                    }
                }
                { prev = $0 }
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no major section format findings"
    fi
}


#  Fail when retired user-facing headings appear in shell help
function scan_retired_user_headings() {
    local found=0
    local file line
    local allow_output=0

    for file in "${files[@]}"; do
        allow_output=0
        if [[ "${file}" == "${ROOT_REPO}/lib/bash/"* ]]; then
            allow_output=1
        fi
        while IFS= read -r line; do
            found=1
            record_fail \
                "retired user-facing help heading:" \
                "$(print_relpath "${file}"):${line}"
        done < <(
            awk -v allow_output="${allow_output}" '
                /^(Description|Dependencies|Testing|#TODO)$/ {
                    print FNR ":" $0
                }
                $0 == "Output" && ! allow_output {
                    print FNR ":" $0
                }
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no retired user-facing help heading findings"
    fi
}


#  Fail when shell Parameters rows still use angle-bracket type syntax
function scan_legacy_parameter_rows() {
    local found=0
    local file line

    for file in "${files[@]}"; do
        while IFS= read -r line; do
            found=1
            record_fail \
                "legacy shell parameter row:" \
                "$(print_relpath "${file}"):${line}"
        done < <(
            perl -ne '
                if (/^Parameters\s*$/) {
                    $in_params = 1;
                    next;
                }
                if (/^(Usage|Returns|Notes|References|See Also|Examples|Expected globals|Generated globals)\s*$/) {
                    $in_params = 0;
                }
                next unless $in_params;
                if (/^\s+\d+\+?\s+\S+\s+<[^>]+>/ ||
                    /^\s*-[A-Za-z0-9][^#\n]*<[^>]+>/) {
                    print "$.:$_";
                }
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no legacy shell parameter row findings"
    fi
}


#  Fail when positional shell Parameters rows put the ordinal inside the type pair
function scan_malformed_positional_parameter_rows() {
    local found=0
    local file line

    for file in "${files[@]}"; do
        while IFS= read -r line; do
            found=1
            record_fail \
                "malformed positional parameter row:" \
                "$(print_relpath "${file}"):${line}"
        done < <(
            perl -ne '
                if (/^Parameters\s*$/) {
                    $in_params = 1;
                    next;
                }
                if (/^(Usage|Returns|Notes|References|See Also|Examples|Expected globals|Generated globals)\s*$/) {
                    $in_params = 0;
                }
                next unless $in_params;
                if (/^\s*\d+\+?\s*:\s*\S+\s*:/) {
                    print "$.:$_";
                }
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no malformed positional parameter row findings"
    fi
}


#  Fail on noncanonical <...> placeholders in argument documentation
function scan_argument_placeholders() {
    local found=0
    local file line
    local -a fil_scn=( "${files[@]}" "${files_py[@]}" )

    for file in "${fil_scn[@]}"; do
        while IFS= read -r line; do
            found=1
            record_fail \
                "noncanonical argument placeholder:" \
                "$(print_relpath "${file}"):${line}"
        done < <(
            perl -ne '
                my $heading = qr/
                    ^(?:Script|Description|Parameters|Arguments|
                    Keyword\ arguments|Positional\ arguments|
                    Required|Optional|Dependencies|Returns|Output|
                    Notes|Examples|Testing|Expected\ globals|
                    Generated\ globals|\#TODO):?\s*$
                /x;

                if (/^Usage\s*$/) {
                    $in_usage = 1;
                    next;
                }
                if (/$heading/) {
                    $in_usage = 0 unless /^Usage\s*$/;
                }

                my $is_arg_row = (
                    /^\s+\d+\+?\s+\S+\s+<[^>]+>/ ||
                    /^\s*-[A-Za-z0-9][^#\n]*<[^>]+>/
                );
                my $is_usage_row = (
                    $in_usage &&
                    /<[^>]+>/ &&
                    m{^\s*(?:[A-Za-z0-9_./-]+|\[?\-\S+)}
                );
                next unless ($is_arg_row || $is_usage_row);

                while (/<([^>]+)>/g) {
                    my $tok = $1;
                    next if $tok =~ /^(flag|str|bool|int|flt|num|path|file|dir|csv|mode|method|format|engine|layout|equation|aligner|algorithm|choice|spec|time|size)$/;
                    print "$.:<$tok>: $_";
                    last;
                }
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no noncanonical argument placeholder findings"
    fi
}


#  Fail when Parameters rows use compact mini-grammar instead of readable types
function scan_parameter_mini_grammar_types() {
    local found=0
    local file line

    for file in "${files[@]}"; do
        while IFS= read -r line; do
            found=1
            record_fail \
                "compact mini-grammar parameter type:" \
                "$(print_relpath "${file}"):${line}"
        done < <(
            perl -ne '
                if (/^Parameters\s*$/) {
                    $in_params = 1;
                    next;
                }
                if (/^(Usage|Returns|Notes|References|See Also|Examples|Expected globals|Generated globals)\s*$/) {
                    $in_params = 0;
                }
                next unless $in_params;
                next unless /^\s+(?:\d+\+?\s+\S+|-{1,2}\S.*?)\s*:\s*(.+)$/;
                my $typ = $1;
                if ($typ =~ /\b(?:enum|csv):/ || $typ =~ /^spec(?:\s|$)/) {
                    print "$.:$_";
                }
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no compact mini-grammar parameter type findings"
    fi
}


#  Fail when choice-set parameter types are not quoted and comma-spaced
function scan_choice_set_format() {
    local found=0
    local file line

    for file in "${files[@]}"; do
        while IFS= read -r line; do
            found=1
            record_fail \
                "choice-set parameter type format:" \
                "$(print_relpath "${file}"):${line}"
        done < <(
            perl -ne '
                if (/^Parameters\s*$/) {
                    $in_params = 1;
                    next;
                }
                if (/^(Usage|Returns|Notes|References|See Also|Examples|Expected globals|Generated globals)\s*$/) {
                    $in_params = 0;
                }
                next unless $in_params;
                next unless /^\s+(?:\d+\+?\s+\S+|-{1,2}\S.*?)\s*:\s*(\{.*\})\s*$/;
                my $typ = $1;
                if ($typ !~ /^\{\x27[^\x27]+\x27(, \x27[^\x27]+\x27)*\}$/) {
                    print "$.:$_";
                }
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no choice-set parameter type format findings"
    fi
}


#  Defer Parameter alias cardinality to exact parser/public-set comparison
function scan_parameter_alias_counts() {
    record_pass \
        "Parameter rows use the complete public alias set without a shell" \
        "alias-count cap"
}


#  Fail on vague user-facing optionality markers
function scan_optional_marker() {
    local found=0
    local file line

    for file in "${files[@]}"; do
        while IFS= read -r line; do
            found=1
            record_fail \
                "vague optionality marker:" \
                "$(print_relpath "${file}"):${line}"
        done < <(grep -n -i -- '(optional)' "${file}" || true)
    done

    if (( found == 0 )); then
        record_pass "no vague optionality marker findings"
    fi
}


#  Warn when prose suggests an option should be bracketed in Usage
function scan_usage_bracket_advisory() {
    local found=0
    local file line

    for file in "${files[@]}"; do
        while IFS= read -r line; do
            found=1
            record_warn \
                "possible Usage bracket mismatch:" \
                "$(print_relpath "${file}"):${line}"
        done < <(
            perl -ne '
                if (/^Usage\s*$/) {
                    $in_usage = 1;
                    @usage = ();
                    next;
                }
                if ($in_usage && /^Parameters\s*$/) {
                    $in_usage = 0;
                    $usage_txt = join(" ", @usage);
                    next;
                }
                if ($in_usage) {
                    push @usage, $_;
                    next;
                }

                if (/^Parameters\s*$/) {
                    $in_params = 1;
                    next;
                }
                if (/^(Usage|Returns|Notes|References|See Also|Examples|Expected globals|Generated globals)\s*$/) {
                    $in_params = 0;
                    $opt = "";
                    $row = "";
                }
                next unless $in_params;

                if (/^\s+(?:\d+\+?\s+\S+|-{1,2}\S.*?)\s*:\s*/) {
                    check_row();
                    $row = $_;
                    $opt = "";
                    while (/--[A-Za-z0-9_]+/g) {
                        $opt = $&;
                        last;
                    }
                    next;
                }
                $row .= $_ if $row;

                END { check_row(); }

                sub check_row {
                    return unless $opt && $usage_txt;
                    return unless $row =~ /(default:|Required (?:when|if)|required (?:when|if)|Used only with|used only with|ignored otherwise)/;
                    my $quoted = quotemeta($opt);
                    my $bad = 0;
                    while ($usage_txt =~ /$quoted\b/g) {
                        my $pos = $-[0];
                        my $prev = substr($usage_txt, $pos - 2, 2);
                        my $prev1 = substr($usage_txt, $pos - 1, 1);
                        next if $prev1 eq "[";
                        next if $prev1 eq "(";
                        next if $prev =~ /\|\s/;
                        $bad = 1;
                    }
                    if ($bad) {
                        print "$.:$opt may need brackets in Usage: $row";
                    }
                }
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no possible Usage bracket mismatch findings"
    fi
}


#  Warn when canonical placeholders are probably too vague for path-like args
function scan_under_specific_placeholders() {
    local found=0
    local file line
    local -a fil_scn=( "${files[@]}" "${files_py[@]}" )

    for file in "${fil_scn[@]}"; do
        while IFS= read -r line; do
            found=1
            record_warn \
                "possibly under-specific argument placeholder:" \
                "$(print_relpath "${file}"):${line}"
        done < <(
            perl -ne '
                my $heading = qr/
                    ^(?:Script|Description|Parameters|Arguments|
                    Keyword\ arguments|Positional\ arguments|
                    Required|Optional|Dependencies|Returns|Output|
                    Notes|Examples|Testing|Expected\ globals|
                    Generated\ globals|\#TODO):?\s*$
                /x;

                if (/^Usage\s*$/) {
                    $in_usage = 1;
                    next;
                }
                if (/$heading/) {
                    $in_usage = 0 unless /^Usage\s*$/;
                }

                my $is_arg_row = (
                    /^\s+\d+\+?\s+\S+\s+<[^>]+>/ ||
                    /^\s*-[A-Za-z0-9][^#\n]*<[^>]+>/
                );
                my $is_usage_row = (
                    $in_usage &&
                    /<str>/ &&
                    m{^\s*(?:[A-Za-z0-9_./-]+|\[?\-\S+)}
                );
                next unless ($is_arg_row || $is_usage_row);
                next unless /<str>/;
                next if /(?:sfx|suffix|pfx|prefix|nam_job|job|smp|sample|method|mode|env_nam|aligner|qname|eqn|typ_out|out_ext|pattern|arr_nam|arr_name|nam_arr|desc|fill|type|tbl_col|\barr\d?\b)/;
                next unless /(?:fil|file|fil_in|fil_out|tbl|table|cfg|pth|path|dir|ref_fa|chr_sizes|log_out|log_err|script|scr_|bam|fastq|bdg|fasta|yaml|index)/i;
                print "$.:$_";
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no possibly under-specific argument placeholder findings"
    fi
}


#  Fail when string-like defaults are shown without required single quotes
function scan_unquoted_string_defaults() {
    local found=0
    local file line
    local -a fil_scn=( "${files[@]}" )

    for file in "${fil_scn[@]}"; do
        while IFS= read -r line; do
            found=1
            record_fail \
                "unquoted string-like default:" \
                "$(print_relpath "${file}"):${line}"
        done < <(
            perl -ne '
                my %string_vars = map { $_ => 1 } qw(
                    aligner aln_typ align_typ bt2_mode bwa_alg cfg_met
                    chr_sizes dir_eo dir_out dir_scr engine env_nam eqn
                    fil_in fil_out if_exists index method mode nam_job
                    out_ext prefix ref_fa retain sfx_pe sfx_se suffix_pe
                    suffix_se tbl_met time typ_out
                );

                if (/\(default:\s*\$\{([A-Za-z_][A-Za-z0-9_]*)\}\)/) {
                    print "$.:$_" if $string_vars{$1};
                    next;
                }

                if (/\(default:\s*([A-Za-z_][A-Za-z0-9_.-]*)\)/) {
                    my $val = $1;
                    next if $val =~ /^(?:true|false)$/;
                    print "$.:$_";
                }
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no unquoted string-like default findings"
    fi
}


#  Fail when prefix uses the retired -pr short alias instead of -px
function scan_retired_prefix_short_alias() {
    local found=0
    local file line

    for file in "${files[@]}"; do
        while IFS= read -r line; do
            found=1
            record_fail \
                "retired prefix short alias:" \
                "$(print_relpath "${file}"):${line}"
        done < <(
            perl -ne '
                next unless /--(?:pfx|prefix)\b/;
                if (/(^|[\s,|])\-pr(?:[,|\s]|$)/) {
                    print "$.:$_";
                }
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no retired prefix short alias findings"
    fi
}


#  Fail when a generic --suffix option omits the canonical -sx/--sfx aliases
function scan_generic_suffix_aliases() {
    local found=0
    local file line

    for file in "${files[@]}"; do
        while IFS= read -r line; do
            found=1
            record_fail \
                "generic suffix alias set:" \
                "$(print_relpath "${file}"):${line}"
        done < <(
            perl -ne '
                next unless /--suffix\b/;
                next if /--suffix\[_-\]/;
                next if /-sx\b/ && /--sfx\b/;
                print "$.:$_";
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no generic suffix alias set findings"
    fi
}


#  Fail when user-facing help contains sourced helper/function inventories
function scan_user_facing_helper_inventories() {
    local found=0
    local file line

    for file in "${files[@]}"; do
        while IFS= read -r line; do
            found=1
            record_fail \
                "user-facing helper/function inventory:" \
                "$(print_relpath "${file}"):${line}"
        done < <(
            perl -ne '
                next if /^\s*#/;

                if (/^\s*(?:-\s*)?(?:Functions|Function scripts|Sourced functions|Shell functions|Helpers|Helper functions|Sourced helpers)\s*:?\s*$/) {
                    print "$.:$_";
                    next;
                }

                if (/^\s*[+-]\s+[A-Za-z_][A-Za-z0-9_.\/-]*\s+##\s*[A-Za-z_][A-Za-z0-9_.\/-]*\s*##\s*$/) {
                    print "$.:$_";
                    next;
                }

                if (/^\s*[+-]\s+help\/[A-Za-z0-9_.\/-]+\s*$/) {
                    print "$.:$_";
                    next;
                }
            ' "${file}"
        )
    done

    if (( found == 0 )); then
        record_pass "no user-facing helper/function inventory findings"
    fi
}


#  Run the authoritative changed-content structural help checker
function scan_structural_help_style() {
    local output=""

    if output="$(
        PYTHONDONTWRITEBYTECODE=1 \
            python3 -m dev.audit.help_style \
                --root "${ROOT_REPO}" 2>&1
    )"; then
        record_pass "${output}"
        return 0
    fi

    while IFS= read -r line; do
        record_fail "${line}"
    done <<< "${output}"
}


#  Enforce changed prose and repository-wide physical Examples commands
function scan_help_heredoc_source_reflow() {
    local output=""

    if output="$(
        PYTHONDONTWRITEBYTECODE=1 \
            python3 -m dev.audit.help_heredoc_reflow \
                --root "${ROOT_REPO}" 2>&1
    )"; then
        record_pass "${output}"
        return 0
    fi

    while IFS= read -r line; do
        record_fail "${line}"
    done <<< "${output}"
}


#  Run the authoritative full-help Examples checker
function scan_help_examples() {
    local output=""

    if output="$(
        PYTHONDONTWRITEBYTECODE=1 \
            python3 -m dev.audit.help_examples \
                --root "${ROOT_REPO}" 2>&1
    )"; then
        record_pass "${output}"
        return 0
    fi

    while IFS= read -r line; do
        record_fail "${line}"
    done <<< "${output}"
}


#  Enforce repository-wide Runtime requirements for every help owner
function scan_runtime_requirements() {
    local output=""

    if output="$(
        PYTHONDONTWRITEBYTECODE=1 \
            python3 -m dev.audit.help_runtime_requirements \
                --root "${ROOT_REPO}" \
                --strict 2>&1
    )"; then
        record_pass "${output}"
        return 0
    fi

    while IFS= read -r line; do
        record_fail "${line}"
    done <<< "${output}"
}


#  Run the authoritative Usage/Parameters alias-set checker
function scan_parameter_alias_sets() {
    local output=""

    if output="$(
        PYTHONDONTWRITEBYTECODE=1 \
            python3 -m dev.audit.help_aliases \
                --root "${ROOT_REPO}" "${files[@]}" 2>&1
    )"; then
        record_pass "${output}"
        return 0
    fi

    while IFS= read -r line; do
        record_fail "${line}"
    done <<< "${output}"
}


#  Check current-diff ShellCheck source mappings
function scan_shellcheck_source_policy() {
    local output=""

    if output="$(
        PYTHONDONTWRITEBYTECODE=1 \
            python3 -m dev.audit.source_policy \
                --root "${ROOT_REPO}" 2>&1
    )"; then
        record_pass "${output}"
        return 0
    fi

    while IFS= read -r line; do
        record_fail "${line}"
    done <<< "${output}"
}


#  Check interface-owner selection for recognized unknown-option branches
function scan_unknown_option_helpers() {
    local output=""

    if output="$(
        PYTHONDONTWRITEBYTECODE=1 \
            python3 -m dev.audit.unknown_option_helpers \
                --root "${ROOT_REPO}" 2>&1
    )"; then
        record_pass "${output}"
        return 0
    fi

    while IFS= read -r line; do
        record_fail "${line}"
    done <<< "${output}"
}


#  Report bounded submit-bootstrap spacing as an advisory
function warn_submit_help_source_spacing() {
    local output=""

    output="$(
        PYTHONDONTWRITEBYTECODE=1 \
            python3 -m dev.audit.source_policy \
                --root "${ROOT_REPO}" --help-source-spacing
    )"
    record_warn "${output}"
}


print_section "${TEST_NAME}"

#  Scan top-level script help text and extracted help files for style drift
files=()
while IFS= read -r file; do
    files+=( "${file}" )
done < <(
    {
        find "${ROOT_REPO}/bin" -maxdepth 1 -type f -name '*.sh' -print
        find "${ROOT_REPO}/lib/bash" -type f -name '*.sh' -print
        find "${ROOT_REPO}/install/scripts" -maxdepth 1 -type f -name '*.sh' -print
        find "${ROOT_REPO}/tests" -maxdepth 1 -type f -name '*.sh' -print
    } | sort
)

#  Scan user-facing Python CLIs for parser and entrypoint consistency
files_py=()
while IFS= read -r file; do
    files_py+=( "${file}" )
done < <(
    {
        find "${ROOT_REPO}/src/protocol_chipseq_signal_norm/cli" \
            -maxdepth 1 -type f -name '*.py' ! -name '__init__.py' -print
    } | sort
)

#  Run lightweight help-style scans
ptn_flg="<fl""g>"
ptn_mlt="<ml""t>"
ptn_dry_run="--dry""-run"
ptn_old_log_var="err""_out"
ptn_old_log_opt="--err""_out"
ptn_old_dir_eo_row="-e""o, --dir_eo"
ptn_old_dir_eo_usage="-e""o|--dir_eo"
ptn_old_dir_eo_parser="-e""o|--dir[_-]eo"
ptn_old_dir_eo_compat="-d""eo|-e""o"
ptn_old_dir_eo_example=" -e""o "
ptn_old_req_flg_row="-r""f, --req_flg"
ptn_old_req_flg_parser="-r""f|--req[_-]flg"

scan_tabs

scan_fixed "literal ${ptn_flg}" "${ptn_flg}"

scan_fixed "literal ${ptn_mlt}" "${ptn_mlt}"

scan_usage_legacy_precision

scan_fixed \
    "stale Options heading" \
    'Options:'

scan_fixed_nonparser \
    "stale dry-run help spelling" \
    "${ptn_dry_run}"

scan_fixed \
    "compatibility alias prose" \
    'Compatibility ali''ases'

scan_fixed \
    "user-facing sourced-helper inventory" \
    'Sourced function scr''ipts:'

scan_fixed \
    "retired stderr/stdout log variable" \
    "${ptn_old_log_var}"

scan_fixed \
    "retired stderr/stdout log option" \
    "${ptn_old_log_opt}"

scan_fixed \
    "retired dir_eo short help row" \
    "${ptn_old_dir_eo_row}"

scan_fixed \
    "retired dir_eo short Usage alias" \
    "${ptn_old_dir_eo_usage}"

scan_fixed \
    "retired dir_eo short parser alias" \
    "${ptn_old_dir_eo_parser}"

scan_fixed \
    "retired dir_eo short compatibility alias" \
    "${ptn_old_dir_eo_compat}"

scan_fixed \
    "retired dir_eo short example" \
    "${ptn_old_dir_eo_example}"

scan_retired_ref_fa_aliases

scan_fixed \
    "retired req_flg short help row" \
    "${ptn_old_req_flg_row}"

scan_fixed \
    "retired req_flg short parser alias" \
    "${ptn_old_req_flg_parser}"

scan_legacy_colon_headings

scan_major_section_format

scan_retired_user_headings

scan_legacy_parameter_rows

scan_malformed_positional_parameter_rows

scan_argument_placeholders

scan_parameter_mini_grammar_types

scan_choice_set_format

scan_parameter_alias_counts

scan_optional_marker

scan_usage_bracket_advisory

scan_under_specific_placeholders

scan_unquoted_string_defaults

scan_retired_prefix_short_alias

scan_generic_suffix_aliases

scan_user_facing_helper_inventories

scan_structural_help_style

scan_help_examples

scan_runtime_requirements

scan_parameter_alias_sets

scan_shellcheck_source_policy

scan_unknown_option_helpers

warn_submit_help_source_spacing

scan_help_heredoc_source_reflow

scan_fixed \
    "stale helper name echo_error" \
    'echo_error'

scan_fixed \
    "stale helper name echo_warning" \
    'echo_warning'

scan_fixed \
    "stale helper dependency submit.sh" \
    'submit\.sh'

finish
