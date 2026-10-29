#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: make.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


# Require Bash >= 4.4 before doing any work.
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(make.sh):" \
        "this script must be run under Bash >= 4.4." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

# Run in safe mode, exiting on errors, unset variables, and pipe failures.
set -euo pipefail


# Resolve paths relative to 'tests/fixtures'.
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_fix="${dir_scr}"

# Source shared fixture-generation helpers.
# shellcheck source=tests/support/fixture_helpers.sh
source "${dir_scr}/../../support/fixture_helpers.sh"


# Write one fixture file, one argument per line.
function write_fixture() {
    local fil="${1}"
    shift 1

    printf '%s\n' "$@" > "${fil}"
}


# Remove stale outputs so regeneration is idempotent, then create the
# directories.
rm_files \
    "${dir_fix}" \
    "${dir_fix}/bdg/bin_10.bedGraph" \
    "${dir_fix}/bdg/bin_10_two_rows.bedGraph" \
    "${dir_fix}/bdg/bin_20.bedGraph" \
    "${dir_fix}/bdg/bin_50.bedGraph" \
    "${dir_fix}/bdg/clamped_terminal_bin.bedGraph" \
    "${dir_fix}/bdg/grid_A.bedGraph" \
    "${dir_fix}/bdg/grid_late_mismatch_B.bedGraph" \
    "${dir_fix}/bdg/no_interval.bedGraph" \
    "${dir_fix}/bdg/one_row.bedGraph" \
    "${dir_fix}/bdg/out_of_bounds.bedGraph" \
    "${dir_fix}/bdg/partial_terminal_bin.bedGraph" \
    "${dir_fix}/bdg/runs.bedGraph" \
    "${dir_fix}/bdg/shared_factor.bedGraph" \
    "${dir_fix}/bdg/sizes_duplicate.tsv" \
    "${dir_fix}/bdg/sizes_mixed_delimiters.tsv" \
    "${dir_fix}/bdg/two_chromosomes.bedGraph" \
    "${dir_fix}/input_floor/five_bins.bedGraph" \
    "${dir_fix}/input_floor/headers_and_blanks.bed" \
    "${dir_fix}/input_floor/one_bin_2.bedGraph" \
    "${dir_fix}/input_floor/one_bin_nonround.bedGraph" \
    "${dir_fix}/input_floor/one_interval.bed" \
    "${dir_fix}/input_floor/two_intervals.bed" \
    "${dir_fix}/merge_bins/epsilon.bedGraph" \
    "${dir_fix}/merge_bins/exact_text.bedGraph" \
    "${dir_fix}/merge_bins/malformed.bedGraph" \
    "${dir_fix}/merge_bins/midstream_header.bedGraph" \
    "${dir_fix}/merge_bins/nonfinite.bedGraph" \
    "${dir_fix}/merge_bins/one_row.bedGraph" \
    "${dir_fix}/merge_bins/rounding.bedGraph" \
    "${dir_fix}/merge_bins/rounding_epsilon.bedGraph" \
    "${dir_fix}/merge_bins/zero_width.bedGraph" \
    "${dir_fix}/pseudo/all_malformed_rows.bedGraph" \
    "${dir_fix}/pseudo/grid_A.bedGraph" \
    "${dir_fix}/pseudo/grid_B.bedGraph" \
    "${dir_fix}/pseudo/one_bin_2.bedGraph" \
    "${dir_fix}/pseudo/one_bin_3.bedGraph" \
    "${dir_fix}/pseudo/one_malformed_row.bedGraph" \
    "${dir_fix}/ratio/bare_pseudocount_A.bedGraph" \
    "${dir_fix}/ratio/bare_pseudocount_B.bedGraph" \
    "${dir_fix}/ratio/dash_io_B.bedGraph" \
    "${dir_fix}/ratio/end_mismatch_A.bedGraph" \
    "${dir_fix}/ratio/end_mismatch_B.bedGraph" \
    "${dir_fix}/ratio/out_of_bounds_A.bedGraph" \
    "${dir_fix}/ratio/out_of_bounds_B.bedGraph" \
    "${dir_fix}/ratio/paired_pseudocount_A.bedGraph" \
    "${dir_fix}/ratio/paired_pseudocount_B.bedGraph" \
    "${dir_fix}/ratio/small_A.bedGraph" \
    "${dir_fix}/ratio/small_B.bedGraph" \
    "${dir_fix}/signal/sizes.tsv" \
    "${dir_fix}/sum/header_and_nan.bedGraph" \
    "${dir_fix}/sum/two_widths.bedGraph"

mkdirs \
    "${dir_fix}/bdg" \
    "${dir_fix}/input_floor" \
    "${dir_fix}/merge_bins" \
    "${dir_fix}/pseudo" \
    "${dir_fix}/ratio" \
    "${dir_fix}/signal" \
    "${dir_fix}/sum"


# Inputs for 'utils_bdg' (tests/unit/utilities/test_bdg.py): bin-size and grid
# checks, chromosome sizes, bounds, and bin counting.
write_fixture \
    "${dir_fix}/bdg/bin_10.bedGraph" \
    'chrI 0 10 1'

write_fixture \
    "${dir_fix}/bdg/bin_10_two_rows.bedGraph" \
    $'chrI\t0\t10\t1' \
    $'chrI\t10\t20\t2'

write_fixture \
    "${dir_fix}/bdg/bin_20.bedGraph" \
    'chrI 0 20 1'

write_fixture \
    "${dir_fix}/bdg/bin_50.bedGraph" \
    $'chrI\t0\t50\t1' \
    $'chrI\t50\t150\t2'

write_fixture \
    "${dir_fix}/bdg/clamped_terminal_bin.bedGraph" \
    $'chrI\t0\t10\t2' \
    $'chrI\t10\t70\t3' \
    $'chrI\t70\t75\t4'

write_fixture \
    "${dir_fix}/bdg/grid_A.bedGraph" \
    'chrI 0 10 1' \
    'chrI 10 20 1' \
    'chrI 20 30 1'

write_fixture \
    "${dir_fix}/bdg/grid_late_mismatch_B.bedGraph" \
    'chrI 0 10 1' \
    'chrI 10 20 1' \
    'chrI 20 31 1'

write_fixture \
    "${dir_fix}/bdg/no_interval.bedGraph" \
    'track name=empty'

write_fixture \
    "${dir_fix}/bdg/one_row.bedGraph" \
    $'chrI\t0\t10\t1'

write_fixture \
    "${dir_fix}/bdg/out_of_bounds.bedGraph" \
    $'chrI\t70\t81\t1'

write_fixture \
    "${dir_fix}/bdg/partial_terminal_bin.bedGraph" \
    $'chrI\t0\t10\t1' \
    $'chrI\t10\t14\t5'

write_fixture \
    "${dir_fix}/bdg/runs.bedGraph" \
    $'chrI\t0\t10\t2' \
    $'chrI\t10\t40\t3'

write_fixture \
    "${dir_fix}/bdg/shared_factor.bedGraph" \
    $'chrI\t0\t20\t1' \
    $'chrI\t20\t40\t2' \
    $'chrI\t40\t60\t3'

write_fixture \
    "${dir_fix}/bdg/sizes_duplicate.tsv" \
    $'chrI\t80' \
    $'chrI\t90'

write_fixture \
    "${dir_fix}/bdg/sizes_mixed_delimiters.tsv" \
    '# comment' \
    $'chrI\t80' \
    'chrII 120'

write_fixture \
    "${dir_fix}/bdg/two_chromosomes.bedGraph" \
    $'chrI\t0\t10\t1' \
    $'chrI\t10\t15\t1' \
    $'chrII\t0\t10\t1' \
    $'chrII\t10\t15\t1'


# Inputs for 'compute_input_floor', in 'test_input_floor.py' and
# 'test_input_floor_opt_applicability.py' (tests/unit/compute_signal/).
write_fixture \
    "${dir_fix}/input_floor/five_bins.bedGraph" \
    'chrI 0 10 0' \
    'chrI 10 20 -1' \
    'chrI 20 30 1' \
    'chrI 30 40 2' \
    'chrI 40 50 10'

write_fixture \
    "${dir_fix}/input_floor/headers_and_blanks.bed" \
    '# h' \
    'track name=x' \
    $'chrI\t0\t10' \
    '' \
    '   ' \
    $'chrI\t10\t20'

write_fixture \
    "${dir_fix}/input_floor/one_bin_2.bedGraph" \
    $'chrI\t0\t10\t2'

write_fixture \
    "${dir_fix}/input_floor/one_bin_nonround.bedGraph" \
    $'chrI\t0\t10\t1.23'

write_fixture \
    "${dir_fix}/input_floor/one_interval.bed" \
    $'chrI\t0\t10'

write_fixture \
    "${dir_fix}/input_floor/two_intervals.bed" \
    $'chrI\t0\t10' \
    $'chrI\t5\t15'


# Inputs for 'merge_bins' (tests/unit/bedgraph/test_merge_bins.py),
# space-delimited as the script accepts.
write_fixture \
    "${dir_fix}/merge_bins/epsilon.bedGraph" \
    'chrI 0 10 1.0000' \
    'chrI 10 20 1.0004' \
    'chrI 20 30 1.0020'

write_fixture \
    "${dir_fix}/merge_bins/exact_text.bedGraph" \
    'chrI 0 10 1' \
    'chrI 10 20 1' \
    'chrI 20 30 1.0' \
    'chrI 30 40 1.0'

write_fixture \
    "${dir_fix}/merge_bins/malformed.bedGraph" \
    'chrI 0 10'

write_fixture \
    "${dir_fix}/merge_bins/midstream_header.bedGraph" \
    'chrI 0 10 2' \
    'chrI 10 20 2' \
    'track type=bedGraph' \
    'chrI 20 30 2' \
    'chrI 30 40 2'

write_fixture \
    "${dir_fix}/merge_bins/nonfinite.bedGraph" \
    'chrI 0 10 NaN' \
    'chrI 10 20 nan' \
    'chrI 20 30 inf' \
    'chrI 30 40 inf'

write_fixture \
    "${dir_fix}/merge_bins/one_row.bedGraph" \
    'chrI 0 10 1'

write_fixture \
    "${dir_fix}/merge_bins/rounding.bedGraph" \
    'chrI 0 10 1.001' \
    'chrI 10 20 1.002' \
    'chrI 20 30 1.014'

write_fixture \
    "${dir_fix}/merge_bins/rounding_epsilon.bedGraph" \
    'chrI 0 10 1.001' \
    'chrI 10 20 1.009' \
    'chrI 20 30 1.021'

write_fixture \
    "${dir_fix}/merge_bins/zero_width.bedGraph" \
    'chrI 10 10 1'


# Inputs for 'compute_pseudo' (tests/unit/compute_signal/test_pseudo.py).
write_fixture \
    "${dir_fix}/pseudo/all_malformed_rows.bedGraph" \
    $'chrI\t0\t10\tnot-a-number'

write_fixture \
    "${dir_fix}/pseudo/grid_A.bedGraph" \
    $'chrI\t0\t10\t2' \
    $'chrI\t10\t20\t4'

write_fixture \
    "${dir_fix}/pseudo/grid_B.bedGraph" \
    $'chrI\t0\t10\t3' \
    $'chrI\t10\t20\t5'

write_fixture \
    "${dir_fix}/pseudo/one_bin_2.bedGraph" \
    $'chrI\t0\t10\t2'

write_fixture \
    "${dir_fix}/pseudo/one_bin_3.bedGraph" \
    $'chrI\t0\t10\t3'

write_fixture \
    "${dir_fix}/pseudo/one_malformed_row.bedGraph" \
    $'chrI\t0\t10\tnot-a-number' \
    $'chrI\t10\t20\t2'


# Inputs for 'comp_sig_rat' (tests/unit/compute_signal/test_ratio.py), as A and
# B pairs, space-delimited.
write_fixture \
    "${dir_fix}/ratio/bare_pseudocount_A.bedGraph" \
    'chrI 0 10 4' \
    'chrI 10 20 5'

write_fixture \
    "${dir_fix}/ratio/bare_pseudocount_B.bedGraph" \
    'chrI 0 10 2' \
    'chrI 10 20 0'

write_fixture \
    "${dir_fix}/ratio/dash_io_B.bedGraph" \
    'chrI 0 10 2'

write_fixture \
    "${dir_fix}/ratio/end_mismatch_A.bedGraph" \
    'chrI 0 10 4' \
    'chrI 10 20 2'

write_fixture \
    "${dir_fix}/ratio/end_mismatch_B.bedGraph" \
    'chrI 0 10 2' \
    'chrI 10 21 1'

write_fixture \
    "${dir_fix}/ratio/out_of_bounds_A.bedGraph" \
    'chrI 70 81 4'

write_fixture \
    "${dir_fix}/ratio/out_of_bounds_B.bedGraph" \
    'chrI 70 80 2'

write_fixture \
    "${dir_fix}/ratio/paired_pseudocount_A.bedGraph" \
    'chrI 0 10 5'

write_fixture \
    "${dir_fix}/ratio/paired_pseudocount_B.bedGraph" \
    'chrI 0 10 0'

write_fixture \
    "${dir_fix}/ratio/small_A.bedGraph" \
    'chrI 0 10 4' \
    'chrI 10 20 0'

write_fixture \
    "${dir_fix}/ratio/small_B.bedGraph" \
    'chrI 0 10 2' \
    'chrI 10 20 0'


# Inputs for 'compute_signal' (tests/unit/compute_signal/test_signal.py): a
# chromosome-size table.
write_fixture \
    "${dir_fix}/signal/sizes.tsv" \
    $'chrI\t80'


# Inputs for 'sum_bdg' (tests/unit/bedgraph/test_sum.py), space-delimited.
write_fixture \
    "${dir_fix}/sum/header_and_nan.bedGraph" \
    'track name=x' \
    'chrI 0 10 1' \
    'chrI 10 20 2' \
    'chrI 20 30 nan'

write_fixture \
    "${dir_fix}/sum/two_widths.bedGraph" \
    'chrI 0 10 1' \
    'chrI 10 30 2'


succeed "generated bedGraph unit-test fixtures under ${dir_fix}"
