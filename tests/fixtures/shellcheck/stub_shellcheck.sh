#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: stub_shellcheck.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5, Opus 5.5).
#
# Distributed under the MIT license.


set -u

if [[ "${1:-}" == "--version" ]]; then
    cat << 'VERSION'
ShellCheck - shell script analysis tool
version: 0.10.0
license: GNU General Public License, version 3
website: https://www.shellcheck.net
VERSION
    exit 0
fi

if [[ -n "${FAKE_SHELLCHECK_LOG:-}" ]]; then
    printf '%s\n' "$*" >> "${FAKE_SHELLCHECK_LOG}"
fi

comment='{"file":"fixture.sh","line":1,"column":1,'
comment+='"level":"warning","code":9999,'
comment+='"message":"fixture finding","fix":null}'

status="${FAKE_SHELLCHECK_STATUS:-0}"
case "${status}" in
    0)
        printf '{"comments":[]}\n'
        ;;
    1)
        printf '{"comments":[%s]}\n' "${comment}"
        ;;
    *)
        echo "fixture infrastructure failure" >&2
        ;;
esac
exit "${status}"
