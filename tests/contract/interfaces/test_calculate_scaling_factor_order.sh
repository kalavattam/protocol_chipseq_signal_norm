#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_calculate_scaling_factor_order.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="calculate_scaling_factor option order"

# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


function assert_order_in_region() {
    local file="${1:?}"
    local label="${2:?}"
    local begin="${3:?}"
    local end="${4:?}"
    local output=""
    shift 4

    if output="$(
        "${TEST_PYTHON}" - "${file}" "${begin}" "${end}" "$@" << 'PY'
import sys
from pathlib import Path

path = Path(sys.argv[1])
begin = sys.argv[2]
end = sys.argv[3]
tokens = sys.argv[4:]
text = path.read_text(encoding="utf-8")
start = text.find(begin)
if start < 0:
    raise SystemExit(f"missing region start {begin!r}")
stop = text.find(end, start + len(begin))
if stop < 0:
    raise SystemExit(f"missing region end {end!r}")
region = text[start:stop]
cursor = -1
for token in tokens:
    position = region.find(token, cursor + 1)
    if position < 0:
        raise SystemExit(f"missing ordered token {token!r}")
    cursor = position
PY
    )"; then
        record_pass "${label}"
    else
        record_fail "${label}: ${output}"
    fi
}


function assert_fixture_order() {
    local output=""

    if output="$(
        "${TEST_PYTHON}" - "$@" << 'PY'
import sys
from pathlib import Path

failures = []
checked = 0
for name in sys.argv[1:]:
    path = Path(name)
    lines = path.read_text(encoding="utf-8").splitlines()
    for index, line in enumerate(lines):
        if "--csv_mip" not in line:
            continue
        start = max(0, index - 12)
        prior = "\n".join(lines[start:index])
        uses_base_array = "arr_cmd" in prior
        if (
            "--aln_typ" not in prior
            and "--aln_typ" not in line
            and not uses_base_array
        ):
            failures.append(f"{path}:{index + 1}: --csv_mip lacks prior --aln_typ")
        checked += 1
if failures:
    raise SystemExit("\n".join(failures))
if checked == 0:
    raise SystemExit("no calculate-scaling-factor fixture commands checked")
print(checked)
PY
    )"; then
        record_pass \
            "${output} fixture command(s) place --aln_typ before --csv_mip"
    else
        while IFS= read -r line; do
            record_fail "fixture option order: ${line}"
        done <<< "${output}"
    fi
}


execute="${ROOT_REPO}/bin/execute_calculate_scaling_factor.sh"
submit="${ROOT_REPO}/bin/submit_calculate_scaling_factor.sh"
help_execute="$(
    printf '%s' \
        "${ROOT_REPO}/lib/bash/help/" \
        "help_execute_calculate_scaling_factor.sh"
)"
help_submit="$(
    printf '%s' \
        "${ROOT_REPO}/lib/bash/help/" \
        "help_submit_calculate_scaling_factor.sh"
)"
log_execute="${TEST_DIR_LOG}/calculate_scaling_factor_execute.help.log"
log_submit="${TEST_DIR_LOG}/calculate_scaling_factor_submit.help.log"
order=(
    --aln_typ
    --ref_fa
    --csv_mip
    --csv_min
    --csv_sip
    --csv_sin
    --fil_out
)


print_section "${TEST_NAME}"

run_capture \
    "execute calculate-scaling-factor help" \
    "${log_execute}" \
    "${TEST_BASH}" "${execute}" --help || true
run_capture \
    "submit calculate-scaling-factor help" \
    "${log_submit}" \
    "${TEST_BASH}" "${submit}" --help || true

for spec in \
    "execute source help:${help_execute}" \
    "submit source help:${help_submit}" \
    "execute rendered help:${log_execute}" \
    "submit rendered help:${log_submit}"
do
    label="${spec%%:*}"
    file="${spec#*:}"
    assert_order_in_region \
        "${file}" "${label} Usage follows canonical option order" \
        $'Usage\n-----' $'\nParameters\n----------' \
        "${order[@]}"
    assert_order_in_region \
        "${file}" "${label} Parameters follows canonical option order" \
        $'Parameters\n----------' $'\nNotes\n-----' \
        "${order[@]}"
done

assert_order_in_region \
    "${execute}" "execute build_cmd forwards canonical option order" \
    "function build_cmd() {" "function build_cmd_cmb() {" \
    "${order[@]}"
assert_order_in_region \
    "${execute}" "execute defaults follow canonical option order" \
    "function init_arg_defs() {" "function init_defs() {" \
    aln_typ= ref_fa= csv_mip= csv_min= csv_sip= csv_sin= fil_out=
assert_order_in_region \
    "${execute}" "execute parser follows canonical option order" \
    "function parse_args() {" "function canonicalize_args() {" \
    '-at|' '-rf|' '-cmip|' '-cmin|' '-csip|' '-csin|' '-fo|'
assert_order_in_region \
    "${execute}" "execute debug report follows canonical option order" \
    "function print_state_debug() {" "function run_jobs() {" \
    aln_typ= ref_fa= csv_mip= csv_min= csv_sip= csv_sin= fil_out=

assert_order_in_region \
    "${submit}" "submit defaults follow canonical option order" \
    "function init_arg_defs() {" "function init_defs() {" \
    aln_typ= ref_fa= csv_mip= csv_min= csv_sip= csv_sin= fil_out=
assert_order_in_region \
    "${submit}" "submit parser follows canonical option order" \
    "function parse_args() {" "function canonicalize_args() {" \
    '-at|' '-rf|' '-cmip|' '-cmin|' '-csip|' '-csin|' '-fo|'
assert_order_in_region \
    "${submit}" "submit debug report follows canonical option order" \
    "function print_state_debug() {" "function prepare_vecs() {" \
    aln_typ= ref_fa= csv_mip= csv_min= csv_sip= csv_sin= fil_out=

assert_fixture_order \
    "${ROOT_REPO}/tests/support/test_helpers.sh" \
    "${ROOT_REPO}/tests/integration/local/calculate_scaling_factor/test_execute_calculate_scaling_factor_siq.sh" \
    "${ROOT_REPO}/tests/integration/local/calculate_scaling_factor/test_execute_calculate_scaling_factor_spike.sh" \
    "${ROOT_REPO}/tests/integration/local/calculate_scaling_factor/test_submit_calculate_scaling_factor_siq.sh" \
    "${ROOT_REPO}/tests/integration/local/calculate_scaling_factor/test_submit_calculate_scaling_factor_spike.sh" \
    "${ROOT_REPO}/tests/integration/parallel/calculate_scaling_factor/test_execute_calculate_scaling_factor_parallel.sh"

finish "$@"
