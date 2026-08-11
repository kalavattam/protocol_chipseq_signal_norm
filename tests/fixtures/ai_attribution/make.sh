#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: make.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


# Require Bash >= 4.4 before doing any work.
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
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


# Resolve paths relative to 'tests/fixtures'
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_fix="${dir_scr}"

# Source shared fixture-generation helpers.
# shellcheck source=tests/support/fixture_helpers.sh
source "${dir_scr}/../../support/fixture_helpers.sh"

# Declare every generated path up front. The directory names the verdict the
# checker must return for the header inside it: the two credited forms are
# accepted, and a coherent no-AI profile is outside what the rule reaches
# rather than a header it accepts.
dir_acc="${dir_fix}/accepted"
dir_nap="${dir_fix}/non_applicable"
fil_single="${dir_acc}/single_vendor.sh"
fil_multi="${dir_acc}/multi_vendor.sh"
fil_none="${dir_nap}/no_attribution.sh"

# Remove stale outputs so regeneration is idempotent.
rm_files "${dir_fix}" "${fil_single}" "${fil_multi}" "${fil_none}"

mkdirs "${dir_acc}" "${dir_nap}"


# Author every fixture literally. They share one bounded header and differ
# only in the attribution block, which is the subject under test. The
# delimiter is quoted so the '$' in no body reaches the shell.

# One vendor credited in the bounded prose form.
cat << 'EOM' > "${fil_single}"
#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: single_vendor.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.
EOM

# Two vendors credited in the lead-in and semicolon-list form.
cat << 'EOM' > "${fil_multi}"
#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: multi_vendor.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.
EOM

# A coherent no-AI profile, which declares no attribution at all.
cat << 'EOM' > "${fil_none}"
#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: no_attribution.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Distributed under the MIT license.
EOM


succeed "generated ai-attribution fixtures under ${dir_fix}"
