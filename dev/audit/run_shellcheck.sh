#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: run_shellcheck.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


#  Run repository ShellCheck with explicit managed-tool provenance
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." > /dev/null 2>&1 && pwd)"
posix_bootstrap="${root}/install/scripts/install_envs_entrypoint.sh"
cd "${root}"


function die_infrastructure() {
    echo "error(run_shellcheck.sh):" "$@" >&2
    exit 2
}


function absolute_file() {
    local supplied="${1}"
    local parent=""

    if [[ "${supplied}" != /* ]]; then
        supplied="${root}/${supplied}"
    fi
    parent="$(cd "$(dirname "${supplied}")" > /dev/null 2>&1 && pwd -P)" ||
        return 1
    printf '%s/%s\n' "${parent}" "$(basename "${supplied}")"
}


provenance_only=false
list_only=false
output_dir="${TEST_ARTIFACT_ROOT:-${root}/artifacts/tests}/shellcheck"
declare -a arr_requested=()
while (( "$#" > 0 )); do
    case "${1}" in
        --provenance-only)
            provenance_only=true
            shift
            ;;
        --list)
            list_only=true
            shift
            ;;
        --output-dir)
            (( "$#" >= 2 )) || die_infrastructure "--output-dir needs a path."
            output_dir="${2}"
            shift 2
            ;;
        --)
            shift
            arr_requested+=( "$@" )
            break
            ;;
        -*)
            die_infrastructure "unknown option: ${1}"
            ;;
        *)
            arr_requested+=( "${1}" )
            shift
            ;;
    esac
done

# Require the managed environment before deriving any tool path.
[[ -n "${CONDA_PREFIX:-}" ]] || \
    die_infrastructure "CONDA_PREFIX must identify env_protocol."

# Resolve managed tools directly, never through PATH.
shellcheck_path="${CONDA_PREFIX}/bin/shellcheck"
python_path="${CONDA_PREFIX}/bin/python"

[[ -x "${shellcheck_path}" ]] || \
    die_infrastructure \
        "managed ShellCheck is not executable: ${shellcheck_path}"
[[ -x "${python_path}" ]] || \
    die_infrastructure "managed Python is not executable: ${python_path}"

# Verify the pinned ShellCheck version from its own reported provenance.
version=""

while IFS=':' read -r key value; do
    if [[ "${key}" == "version" ]]; then
        version="${value# }"
        break
    fi
done < <("${shellcheck_path}" --version)

[[ "${version}" == "0.10.0" ]] || \
    die_infrastructure \
        "expected ShellCheck 0.10.0; found ${version:-unknown}."

printf 'ShellCheck environment prefix: %s\n' "${CONDA_PREFIX}"
printf 'ShellCheck executable: %s\n' "${shellcheck_path}"
printf 'ShellCheck version: %s\n' "${version}"

if [[ "${provenance_only}" == "true" ]]; then
    (( ${#arr_requested[@]} == 0 )) || \
        die_infrastructure "--provenance-only accepts no paths."

    exit 0
fi

# Split discovered or supplied paths into Bash and POSIX dialect sets.
declare -a arr_bash_paths=()
declare -a arr_sh_paths=()

if (( ${#arr_requested[@]} > 0 )); then
    for supplied in "${arr_requested[@]}"; do
        resolved="$(absolute_file "${supplied}")" || \
            die_infrastructure "cannot resolve path: ${supplied}"
        [[ -f "${resolved}" ]] || \
            die_infrastructure "ShellCheck path is not a file: ${resolved}"

        if [[ "${resolved}" == "${posix_bootstrap}" ]]; then
            arr_sh_paths+=( "${resolved}" )
        else
            arr_bash_paths+=( "${resolved}" )
        fi
    done
else
    while IFS= read -r -d '' relative; do
        [[ -f "${root}/${relative}" ]] || continue

        case "${relative}" in
            bin/*.sh|\
            lib/bash/*.sh|\
            lib/bash/**/*.sh|\
            install/scripts/*.sh|\
            tests/*.sh|\
            tests/**/*.sh)
                resolved="${root}/${relative}"
                if [[ "${resolved}" == "${posix_bootstrap}" ]]; then
                    arr_sh_paths+=( "${resolved}" )
                else
                    arr_bash_paths+=( "${resolved}" )
                fi
                ;;
        esac
    done < <(
        git -C "${root}" ls-files -z --cached --others --exclude-standard -- \
            '*.sh'
    )
fi

if [[ "${list_only}" == "true" ]]; then
    if (( ${#arr_bash_paths[@]} > 0 )); then
        for resolved in "${arr_bash_paths[@]}"; do
            printf 'bash\t%s\n' "${resolved}"
        done
    fi

    if (( ${#arr_sh_paths[@]} > 0 )); then
        for resolved in "${arr_sh_paths[@]}"; do
            printf 'sh\t%s\n' "${resolved}"
        done
    fi

    exit 0
fi

if [[ "${output_dir}" != /* ]]; then
    output_dir="${root}/${output_dir}"
fi

mkdir -p "${output_dir}"
bash_raw="${output_dir}/bash_findings.json"
sh_raw="${output_dir}/sh_findings.json"
summary="${output_dir}/summary.json"


function run_language() {
    local shell_mode="${1}"
    local raw_output="${2}"
    shift 2
    local -a arr_files=( "$@" )
    local -a arr_batch=()
    local -a arr_parts=()
    local batch_size=5
    local batch_status=0
    local end=0
    local index=0
    local language_status=0
    local merge_status=0
    local part=""
    local parts_dir=""

    # Emit an empty but well-formed inventory when a dialect has no files.
    if (( ${#arr_files[@]} == 0 )); then
        printf '{"comments":[]}\n' > "${raw_output}"

        return 0
    fi

    parts_dir="$(mktemp -d "${output_dir}/.${shell_mode}.parts.XXXXXX")" || \
        return 2

    # Scan in bounded batches, keeping the worst observed status.
    for (( index = 0; index < ${#arr_files[@]}; index += batch_size )); do
        end=$(( index + batch_size ))
        (( end > ${#arr_files[@]} )) && end="${#arr_files[@]}"
        arr_batch=( "${arr_files[@]:index:end-index}" )
        part="${parts_dir}/part_$(( index / batch_size )).json"
        arr_parts+=( "${part}" )

        "${shellcheck_path}" --format=json1 --shell="${shell_mode}" -x \
            -P "${root}" "${arr_batch[@]}" > "${part}"

        batch_status="$?"

        if (( batch_status > 1 )); then
            language_status="${batch_status}"
        elif (( batch_status == 1 && language_status == 0 )); then
            language_status=1
        fi
    done

    PYTHONDONTWRITEBYTECODE=1 "${python_path}" - \
        "${raw_output}" "${arr_parts[@]}" << 'PY'
import json
import sys
from pathlib import Path


output = Path(sys.argv[1])
comments: list[object] = []
for supplied in sys.argv[2:]:
    value = json.loads(Path(supplied).read_text(encoding="utf-8"))
    if isinstance(value, dict) and isinstance(value.get("comments"), list):
        comments.extend(value["comments"])
    elif isinstance(value, list):
        comments.extend(value)
    else:
        raise ValueError(f"unexpected ShellCheck JSON in {supplied}")
output.write_text(
    json.dumps({"comments": comments}, sort_keys=True) + "\n",
    encoding="utf-8",
)
PY
    merge_status="$?"

    rm -rf -- "${parts_dir}"

    (( merge_status == 0 )) || return 2

    return "${language_status}"
}


# Scan each dialect separately so empty-language evidence is preserved.
set +e

if (( ${#arr_bash_paths[@]} > 0 )); then
    run_language bash "${bash_raw}" "${arr_bash_paths[@]}"
else
    run_language bash "${bash_raw}"
fi

bash_status="$?"

if (( ${#arr_sh_paths[@]} > 0 )); then
    run_language sh "${sh_raw}" "${arr_sh_paths[@]}"
else
    run_language sh "${sh_raw}"
fi

sh_status="$?"

set -e

if (( bash_status > 1 )); then
    status="${bash_status}"
elif (( sh_status > 1 )); then
    status="${sh_status}"
elif (( bash_status == 1 || sh_status == 1 )); then
    status=1
else
    status=0
fi

set +e
report_values="$(
    PYTHONDONTWRITEBYTECODE=1 "${python_path}" - \
        "${summary}" "${shellcheck_path}" "${version}" \
        "${bash_raw}" "${#arr_bash_paths[@]}" "${bash_status}" \
        "${sh_raw}" "${#arr_sh_paths[@]}" "${sh_status}" \
        "${status}" << 'PY'
import json
import sys
from pathlib import Path


def load_findings(path: str) -> tuple[list[object], bool]:
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return [], False
    if isinstance(value, dict) and isinstance(value.get("comments"), list):
        return value["comments"], True
    if isinstance(value, list):
        return value, True
    return [], False


(
    summary_path,
    executable,
    version,
    bash_path,
    bash_files,
    bash_status,
    sh_path,
    sh_files,
    sh_status,
    requested_status,
) = sys.argv[1:]
bash_findings, bash_valid = load_findings(bash_path)
sh_findings, sh_valid = load_findings(sh_path)
status = int(requested_status)
if not bash_valid or not sh_valid:
    status = 2
disposition = (
    "clean"
    if status == 0
    else "completed_with_findings"
    if status == 1
    else "infrastructure_failure"
)
summary = {
    "schema_version": 1,
    "executable": executable,
    "version": version,
    "source_path": str(Path.cwd()),
    "flags": ["-x", "-P", str(Path.cwd())],
    "file_count": int(bash_files) + int(sh_files),
    "finding_count": len(bash_findings) + len(sh_findings),
    "status": status,
    "disposition": disposition,
    "languages": {
        "bash": {
            "shell": "bash",
            "file_count": int(bash_files),
            "finding_count": len(bash_findings),
            "status": int(bash_status),
            "raw_findings": Path(bash_path).name,
            "raw_valid": bash_valid,
        },
        "sh": {
            "shell": "sh",
            "file_count": int(sh_files),
            "finding_count": len(sh_findings),
            "status": int(sh_status),
            "raw_findings": Path(sh_path).name,
            "raw_valid": sh_valid,
        },
    },
}
Path(summary_path).write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n",
    encoding="utf-8",
)
print(f"{summary['finding_count']}\t{status}")
PY
)"
report_status="$?"
set -e
(( report_status == 0 )) || \
    die_infrastructure "failed to produce structured ShellCheck evidence."
IFS=$'\t' read -r finding_count status <<< "${report_values}"

printf 'ShellCheck file count: %d\n' \
    "$(( ${#arr_bash_paths[@]} + ${#arr_sh_paths[@]} ))"
printf 'ShellCheck finding count: %d\n' "${finding_count}"
printf 'ShellCheck evidence: %s\n' "${summary}"
case "${status}" in
    0) echo "ShellCheck status: clean" ;;
    1) echo "ShellCheck status: completed with findings" >&2 ;;
    *) echo "ShellCheck status: infrastructure failure (${status})" >&2 ;;
esac
exit "${status}"
