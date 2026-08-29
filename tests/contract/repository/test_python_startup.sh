#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_python_startup.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="python startup"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


log_env_rslv="${TEST_DIR_LOG}/python_startup/resolve_env_python.log"
log_pol_py="${TEST_DIR_LOG}/python_startup/python_version_policy.log"
log_def_spk="${TEST_DIR_LOG}/python_startup/calculate_scaling_factor_spike_default.log"
log_ret_in="${TEST_DIR_LOG}/python_startup/compute_signal_retired_input.log"
log_ret_out="${TEST_DIR_LOG}/python_startup/compute_signal_retired_output.log"
log_ret_i="${TEST_DIR_LOG}/python_startup/compute_signal_retired_short_i.log"
log_tup_hlp="${TEST_DIR_LOG}/python_startup/help_literal_tuples.log"
ret_alias_in="--in""file"
ret_alias_out="--out""file"
ret_alias_short="-""i"

hlp_scr=(
    "src/protocol_chipseq_signal_norm/cli/add_coeffs_namespaced.py"
    "src/protocol_chipseq_signal_norm/cli/calculate_scaling_factor_siqchip.py"
    "src/protocol_chipseq_signal_norm/cli/calculate_scaling_factor_spike.py"
    "src/protocol_chipseq_signal_norm/cli/compute_input_floor.py"
    "src/protocol_chipseq_signal_norm/cli/compute_pseudo.py"
    "src/protocol_chipseq_signal_norm/cli/compute_signal.py"
    "src/protocol_chipseq_signal_norm/cli/compute_signal_ratio.py"
    "src/protocol_chipseq_signal_norm/cli/merge_bins_bdg.py"
    "src/protocol_chipseq_signal_norm/cli/parse_metadata_siqchip.py"
    "src/protocol_chipseq_signal_norm/cli/relativize_scaling_factors.py"
    "src/protocol_chipseq_signal_norm/cli/sum_bdg.py"
)


print_section "${TEST_NAME}"

# Resolve the project environment locally for dependency-backed Python checks.
require_env_project env_nam || {
    finish
    exit $?
}

if [[ -n "${CONDA_DEFAULT_ENV:-}" && "${CONDA_DEFAULT_ENV}" != "base" ]]; then
    if ! \
        py="$(find_python)"
    then
        record_fail "active project environment has no python/python3 on PATH"
        finish
        exit $?
    fi

    py_cmd=( "${py}" )
else
    # shellcheck disable=SC2154
    if \
        run_capture \
            "resolve env python" \
            "${log_env_rslv}" \
            conda run -n "${env_nam}" python -c \
                'import sys; print(sys.executable)'
    then
        IFS= read -r py < "${log_env_rslv}"
    else
        record_fail \
            "failed to resolve python from '${env_nam}'; see" \
            "$(print_relpath "${log_env_rslv}")"
        finish
        exit $?
    fi

    if [[ -z "${py}" || ! -x "${py}" ]]; then
        record_fail \
            "resolved python is not executable; see" \
            "$(print_relpath "${log_env_rslv}")"
        finish
        exit $?
    fi

    py_cmd=( "${py}" )
fi

if ! \
    "${py_cmd[@]}" - << PY
import sys
raise SystemExit(0 if sys.version_info >= (3, 11) else 1)
PY
then
    record_fail \
        "$("${py_cmd[@]}" --version 2>&1) is older than Python 3.11"
    finish
    exit $?
fi

if \
    PYTHONDONTWRITEBYTECODE=1 \
    run_capture \
        "python version policy" \
        "${log_pol_py}" \
        "${py_cmd[@]}" -m dev.audit.python_version_policy \
            --root "${ROOT_REPO}" \
            --strict
then
    record_pass "repository Python >= 3.11 policy"
else
    record_fail \
        "repository Python >= 3.11 policy; see" \
        "$(print_relpath "${log_pol_py}")"
fi

# Run lightweight Python syntax checks without writing pycache in-place.
while IFS= read -r file; do
    rel="$(print_relpath "${file}")"
    log="${TEST_DIR_LOG}/python_compile/${rel//\//__}.log"

    if \
        PYTHONDONTWRITEBYTECODE=1 \
        PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
        PYTHONPATH="${ROOT_REPO}/src" \
        run_capture \
            "python compile ${rel}" \
            "${log}" \
            "${py_cmd[@]}" -m py_compile "${file}"
    then
        record_pass "python syntax ${rel}"
    else
        record_fail "python syntax ${rel}; see $(print_relpath "${log}")"
    fi
done < <(
    find "${ROOT_REPO}/src/protocol_chipseq_signal_norm" \
        -type f -name '*.py' -print \
        | sort
)

# Reconcile the declared roster against what the package ships: a hardcoded
# list cannot notice a twelfth CLI, and a missing entry degraded to a skip.
cli_found=()

while IFS= read -r file; do
    cli_found+=( "$(print_relpath "${file}")" )
done < <(
    find "${ROOT_REPO}/src/protocol_chipseq_signal_norm/cli" \
        -type f -name '*.py' ! -name '__init__.py' -print \
        | LC_ALL=C sort
)

if [[ "${cli_found[*]}" == "${hlp_scr[*]}" ]]; then
    record_pass \
        "declared help roster matches the shipped CLI modules" \
        "(${#cli_found[@]})"
else
    # Name the differing entries rather than dumping both rosters.
    cnt_dif=0

    while IFS= read -r rel; do
        cnt_dif=$(( cnt_dif + 1 ))

        record_fail \
            "shipped CLI module absent from the declared help roster: ${rel}"
    done < <(
        comm -13 \
            <(printf '%s\n' "${hlp_scr[@]}"   | LC_ALL=C sort) \
            <(printf '%s\n' "${cli_found[@]}" | LC_ALL=C sort)
    )

    while IFS= read -r rel; do
        cnt_dif=$(( cnt_dif + 1 ))

        record_fail \
            "declared help roster entry is not a shipped CLI module:" \
            "${rel}"
    done < <(
        comm -23 \
            <(printf '%s\n' "${hlp_scr[@]}"   | LC_ALL=C sort) \
            <(printf '%s\n' "${cli_found[@]}" | LC_ALL=C sort)
    )

    # Equal sets in a different order reach this branch with no difference.
    if (( cnt_dif == 0 )); then
        record_fail \
            "declared help roster holds the shipped CLI modules but is not" \
            "in sorted order"
    fi
fi

# 'CapArgumentParser' sets 'add_help=False', so a module that omits
# 'add_help_cap()' carries no help flag at all. Render both forms.
cnt_hlp=0

for rel in "${cli_found[@]}"; do
    file="${ROOT_REPO}/${rel}"

    if [[ ! -f "${file}" ]]; then
        record_fail "python help ${rel}: file not present"
        continue
    fi

    for flg in "--help" "-h"; do
        log="${TEST_DIR_LOG}/python_help/${rel//\//__}${flg//-/_}.log"

        if \
            PYTHONDONTWRITEBYTECODE=1 \
            PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
            PYTHONPATH="${ROOT_REPO}/src" \
            run_capture \
                "python help ${flg} ${rel}" \
                "${log}" \
                "${py_cmd[@]}" "${file}" "${flg}"
        then
            assert_file_nonempty "${log}" "python ${flg} ${rel} help"

            cnt_hlp=$(( cnt_hlp + 1 ))
        else
            record_fail \
                "python ${flg} ${rel} failed; see $(print_relpath "${log}")"
        fi
    done
done

# Assert the inspected count, not the finding count alone.
if (( cnt_hlp == 2 * ${#cli_found[@]} )); then
    record_pass \
        "help rendered in both forms for ${#cli_found[@]} CLI modules"
else
    record_fail \
        "help render coverage is ${cnt_hlp} of $(( 2 * ${#cli_found[@]} ))" \
        "expected renders"
fi

# Stray commas in a parenthesized concatenation make 'help=' a tuple, which
# 'argparse' raises on only at render time. 'PY.CLI.HELP.LAYOUT' scopes itself
# to adjacent string literals, so its checker cannot see this shape.
cnt_py="$(
    find "${ROOT_REPO}/src" -type f -name '*.py' -print | wc -l | tr -d ' '
)"

if \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
    run_capture \
        "help literal tuple sweep" \
        "${log_tup_hlp}" \
        "${py_cmd[@]}" - "${ROOT_REPO}/src" << 'PY'
import ast
import pathlib
import sys

# 'metavar' takes a tuple when 'nargs' is set; the rest are always strings.
KEY_STR = ("help", "description", "epilog")

root = pathlib.Path(sys.argv[1])
paths = sorted(root.rglob("*.py"))
found = []

for path in paths:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue

        keys = {kwd.arg for kwd in node.keywords}

        for kwd in node.keywords:
            if not isinstance(kwd.value, ast.Tuple):
                continue

            if kwd.arg in KEY_STR or (
                kwd.arg == "metavar" and "nargs" not in keys
            ):
                found.append(
                    f"{path}:{kwd.value.lineno}: {kwd.arg}= is a tuple"
                )

print(f"inspected {len(paths)}")

for line in found:
    print(line)

raise SystemExit(1 if found else 0)
PY
then
    assert_pattern_found \
        "${log_tup_hlp}" \
        "^inspected ${cnt_py}$" \
        "help literal tuple sweep inspected ${cnt_py} Python files"
else
    record_fail \
        "tuple-valued help keyword in maintained Python, or the sweep failed" \
        "to run; see $(print_relpath "${log_tup_hlp}")"
fi

# Confirm retired compute_signal.py input/output aliases are not accepted.
if \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
    PYTHONPATH="${ROOT_REPO}/src" \
    run_capture \
        "python compute_signal rejects retired input alias" \
        "${log_ret_in}" \
        "${py_cmd[@]}" "${ROOT_REPO}/src/protocol_chipseq_signal_norm/cli/compute_signal.py" \
            "${ret_alias_in}" x \
            --fil_in canonical.bam \
            --fil_out y
then
    record_fail "compute_signal.py unexpectedly accepts retired input alias"
else
    assert_pattern_found \
        "${log_ret_in}" \
        "unrecognized arguments: ${ret_alias_in}" \
        "compute_signal.py rejects retired input alias"
fi

if \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
    PYTHONPATH="${ROOT_REPO}/src" \
    run_capture \
        "python compute_signal rejects retired output alias" \
        "${log_ret_out}" \
        "${py_cmd[@]}" "${ROOT_REPO}/src/protocol_chipseq_signal_norm/cli/compute_signal.py" \
            --fil_in x \
            --fil_out canonical.bdg \
            "${ret_alias_out}" y
then
    record_fail "compute_signal.py unexpectedly accepts retired output alias"
else
    assert_pattern_found \
        "${log_ret_out}" \
        "unrecognized arguments: ${ret_alias_out}" \
        "compute_signal.py rejects retired output alias"
fi

if \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
    PYTHONPATH="${ROOT_REPO}/src" \
    run_capture \
        "python compute_signal rejects retired short input alias" \
        "${log_ret_i}" \
        "${py_cmd[@]}" "${ROOT_REPO}/src/protocol_chipseq_signal_norm/cli/compute_signal.py" \
            "${ret_alias_short}" x \
            --fil_in canonical.bam \
            --fil_out y
then
    record_fail \
        "compute_signal.py unexpectedly accepts retired short input alias"
else
    assert_pattern_found \
        "${log_ret_i}" \
        "unrecognized arguments: ${ret_alias_short}" \
        "compute_signal.py rejects retired short input alias"
fi

# Confirm the direct spike-in calculator defaults to the ChIP-Rx ratio.
if \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
    PYTHONPATH="${ROOT_REPO}/src" \
    run_capture \
        "python calculate-scaling-factor spike default" \
        "${log_def_spk}" \
        "${py_cmd[@]}" -m protocol_chipseq_signal_norm.cli.calculate_scaling_factor_spike \
            --main_ip 100 \
            --spike_ip 10 \
            --main_in 100 \
            --spike_in 20
then
    assert_pattern_found \
        "${log_def_spk}" \
        '^2$' \
        "calculate_scaling_factor_spike.py defaults to chiprx_alpha_ratio"
else
    record_fail \
        "calculate_scaling_factor_spike.py default calculation failed; see" \
        "$(print_relpath "${log_def_spk}")"
fi

finish
