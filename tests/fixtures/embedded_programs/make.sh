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

# Declare every generated path up front. Each verdict directory is a small
# repository root, since the checker judges a file by its path from the root.
dir_acc="${dir_fix}/accepted"
dir_rej="${dir_fix}/rejected"
dir_bnd="${dir_fix}/boundary"
dir_nap="${dir_fix}/nonapplicable"

fil_acc_sh="${dir_acc}/bin/short_programs.sh"
fil_acc_tst="${dir_acc}/tests/contract/test_marked.sh"
fil_acc_py="${dir_acc}/tests/unit/test_short.py"
fil_acc_rcp="${dir_acc}/tests/fixtures/demo/make.sh"
fil_rej_sh="${dir_rej}/bin/long_programs.sh"
fil_rej_tst="${dir_rej}/tests/contract/test_writes.sh"
fil_rej_py="${dir_rej}/tests/unit/test_literal.py"
fil_rej_rcp="${dir_rej}/tests/fixtures/demo/make.sh"
fil_bnd_sh="${dir_bnd}/bin/limits.sh"
fil_bnd_tst="${dir_bnd}/tests/contract/test_strings.sh"
fil_bnd_py="${dir_bnd}/tests/unit/test_three.py"
fil_bnd_rcp="${dir_bnd}/tests/fixtures/demo/make.sh"
fil_nap_rcp="${dir_nap}/tests/fixtures/demo/make.sh"
fil_nap_gen="${dir_nap}/tests/fixtures/demo/generated/test_generated.sh"
fil_nap_hlp="${dir_nap}/tests/integration/helper.py"


# Remove stale outputs so regeneration is idempotent.
rm_files \
    "${dir_fix}" \
    "${fil_acc_sh}" \
    "${fil_acc_tst}" \
    "${fil_acc_py}" \
    "${fil_acc_rcp}" \
    "${fil_rej_sh}" \
    "${fil_rej_tst}" \
    "${fil_rej_py}" \
    "${fil_rej_rcp}" \
    "${fil_bnd_sh}" \
    "${fil_bnd_tst}" \
    "${fil_bnd_py}" \
    "${fil_bnd_rcp}" \
    "${fil_nap_rcp}" \
    "${fil_nap_gen}" \
    "${fil_nap_hlp}"

mkdirs \
    "${dir_acc}/bin" \
    "${dir_acc}/tests/contract" \
    "${dir_acc}/tests/unit" \
    "${dir_acc}/tests/fixtures/demo" \
    "${dir_rej}/bin" \
    "${dir_rej}/tests/contract" \
    "${dir_rej}/tests/unit" \
    "${dir_rej}/tests/fixtures/demo" \
    "${dir_bnd}/bin" \
    "${dir_bnd}/tests/contract" \
    "${dir_bnd}/tests/unit" \
    "${dir_bnd}/tests/fixtures/demo" \
    "${dir_nap}/tests/fixtures/demo/generated" \
    "${dir_nap}/tests/integration"


# Write programs within their limits. Help text in a heredoc may show longer
# programs, which are not run and so are not counted.
cat << 'EOM' > "${fil_acc_sh}"
#!/usr/bin/env bash
python3 -c 'import sys; print(sys.argv)'
bash -c 'echo one'
awk '{ print $1 }' input.txt
sed 's/a/b/' input.txt
sed -E 's/a+/b/' input.txt
function show_help() {
    cat << HELP
awk '
    { a = 1 }
    { b = 2 }
    { c = 3 }
    { d = 4 }
'
python3 - << 'PY'
HELP
}
EOM

# Write a test whose writes are marked copies, two of them continued commands,
# and whose other redirections capture output rather than write inputs.
cat << 'EOM' > "${fil_acc_tst}"
#!/usr/bin/env bash
cp "${src}" "${tmp}/input.txt"  # fixture-perturbation: the code rewrites it
cp \
    "${src}" \
    "${tmp}/other.txt"  # fixture-perturbation: the code writes beside it
printf '%s\n' "${out}" > "${dir_log}/run.log"
printf '%s\n' "${out}" >&2
mkdir -p "${tmp}"
cp -R \
    "${src}/a" \
    "${src}/b" \
    "${src}/c" \
    "${tmp}/"  # fixture-perturbation: the build writes beside it
EOM

# Write a Python test with a long docstring, a three-line literal, and a copy
# marked on its statement's last line.
cat << 'EOM' > "${fil_acc_py}"
"""
A docstring
of four
or more lines.
"""

import shutil


def test_short() -> None:
    lines = """a
b
c"""
    assert lines


def test_copy(tmp_path) -> None:
    shutil.copyfile(
        "input.txt", tmp_path / "input.txt"
    )  # fixture-perturbation: the code rewrites it
EOM

# Write one program over each limit: Python through a heredoc and as a two-line
# '-c' program, two-line 'bash -c' and 'sed', four-line 'awk', a GNU-only 'sed'
# option, and 'perl'.
cat << 'EOM' > "${fil_rej_sh}"
#!/usr/bin/env bash
python3 - << 'PY'
print("heredoc")
PY
python3 -c 'import sys
print(sys.argv)'
bash -c 'echo one
echo two'
awk '
    { a = 1 }
    { b = 2 }
    { c = 3 }
    { print }
' input.txt
sed '
s/a/b/
s/c/d/
' input.txt
sed -i 's/a/b/' input.txt
perl -e 'print 1'
EOM

# Write a test with unmarked writes and a marker that gives no reason.
cat << 'EOM' > "${fil_rej_tst}"
#!/usr/bin/env bash
printf 'x\n' > "${tmp}/input.txt"
cp "${src}" "${tmp}/copy.txt"
cp "${src}" "${tmp}/kept.txt"  # fixture-perturbation:
: > "${tmp}/empty.txt"
cat > "${tmp}/doc.txt" << 'DOC'
body
DOC
EOM

# Write a Python test that writes a four-line literal as its input, and an
# f-string with one short and one long piece between its fields, which counts
# once as a whole.
cat << 'EOM' > "${fil_rej_py}"
from pathlib import Path


def test_reads(tmp_path: Path) -> None:
    text = """a
b
c
d
"""
    (tmp_path / "x.txt").write_text(text)


def test_formats() -> None:
    name = "x"
    text = f"""{name}
a
{name}
b
c
d
"""
    assert text
EOM

# Write programs at their limits: three 'awk' lines, and a one-line 'bash -c'
# inside a command substitution.
cat << 'EOM' > "${fil_bnd_sh}"
#!/usr/bin/env bash
awk '
    { a = 1 }
    { b = 2 }
    { c = 3 }
' input.txt
out="$(
    "${TEST_BASH}" -c 'echo one'
)"
EOM

# Write a test whose writer text sits inside a multi-line string and inside an
# 'awk' program, where the test writes nothing.
cat << 'EOM' > "${fil_bnd_tst}"
#!/usr/bin/env bash
exp="first line
printf 'x' > not_a_file
last line"
got="$(awk '{ printf "%s", (NR > 1 ? " " : "") }' input.txt)"
EOM

# Write a Python test with a literal and an f-string of exactly three lines.
cat << 'EOM' > "${fil_bnd_py}"
def test_three() -> None:
    text = """a
b
c"""
    assert text


def test_three_formatted() -> None:
    name = "x"
    text = f"""{name}
b
{name}"""
    assert text
EOM

# Write recipes that install stubs: one that copies a stub and makes that copy
# executable, one that makes files executable without copying them from a stub,
# once at the start of a line and once indented in a loop, and one whose
# 'chmod' is continued across two lines.
cat << 'EOM' > "${fil_acc_rcp}"
#!/usr/bin/env bash
cp "${dir_stb}/stub_tool.sh" "${dir_tool}/tool"
chmod +x "${dir_tool}/tool"
EOM

cat << 'EOM' > "${fil_rej_rcp}"
#!/usr/bin/env bash
cat << 'STUB' > "${dir_tool}/stub.sh"
python3 - << 'PY'
print("stub")
PY
STUB
chmod +x "${dir_tool}/stub.sh"
for name in a b; do
    chmod +x "${dir_tool}/${name}"
done
EOM

cat << 'EOM' > "${fil_bnd_rcp}"
#!/usr/bin/env bash
cp "${dir_stb}/stub_ok.sh" "${dir_tool}/ok"
cp "${dir_stb}/stub_fail.sh" "${dir_tool}/fail"
chmod +x "${dir_tool}/ok" \
    "${dir_tool}/fail"
EOM

# Write a recipe whose 'chmod' sits inside a script it writes as data, a
# generated fixture that holds a Python heredoc, and a helper that is not a
# test.
cat << 'EOM' > "${fil_nap_rcp}"
#!/usr/bin/env bash
cat << 'DATA' > "${dir_script}/install.sh"
chmod +x tool
DATA
EOM

cat << 'EOM' > "${fil_nap_gen}"
#!/usr/bin/env bash
python3 - << 'PY'
print("generated")
PY
EOM

cat << 'EOM' > "${fil_nap_hlp}"
from pathlib import Path


def write_bundle(path: Path) -> None:
    path.write_text("bundle\n", encoding="utf-8")
EOM


succeed "generated embedded-program fixtures under ${dir_fix}"
