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


# Resolve paths relative to 'tests/fixtures'
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_fix="${dir_scr}"

# Source shared fixture-generation helpers.
# shellcheck source=tests/support/fixture_helpers.sh
source "${dir_scr}/../../support/fixture_helpers.sh"

# Declare every generated path up front. This fixture set is one record of
# every verdict at once rather than a set of inputs each carrying one, so it
# sits at the fixture root with no verdict directories beneath it.
fil_cases="${dir_fix}/cases.json"


# Remove stale outputs so regeneration is idempotent.
rm_file "${dir_fix}" "${fil_cases}"


# Author the record literally, in the canonical form owned by
# 'JSON.SOURCE.FORM'. The delimiter is quoted, so no '$' in the body reaches
# the shell.
cat << 'EOM' > "${fil_cases}"
{
  "accepted": {
    "schema_version": 1,
    "capture_status": "reconstructed_after_fact",
    "sources": [
      {
        "kind": "immutable_accepted_artifact",
        "locator": "S3-STD-008",
        "sha256": "7c8b5424ae9a0e05f2ed37a4636e859f2b3369ade60c24489be55eb3986b447e"
      }
    ],
    "record_id": "option_order",
    "rule_ids": ["S3-STD-008"],
    "old_owner": "HELP.SOURCE_STYLE",
    "old_obligation": "semantic option-order relationships",
    "new_owner": "HELP.OPTION.ORDER",
    "source_anchor": "docs/standards/help.md#source-help-structure",
    "destination_anchor": "docs/standards/help.md#semantic-option-order",
    "source_fingerprint": "70618498e2bd37b13a820068d640732bbbea92d92a966379b33ea0f8df5d79c4",
    "destination_fingerprint": "b1b587467856e82ede919cc7d60c62d1deb43a50704e000cbc35d5ccb2a5942c",
    "old_to_new_complete": true,
    "old_to_new_dispositions": [
      {
        "obligation": "shared semantic precedence",
        "disposition": "moved",
        "evidence": "S3-STD-008"
      }
    ],
    "new_to_old_sources": ["HELP.SOURCE_STYLE", "S3-STD-008"],
    "new_to_old_provenance": [
      {"new_fact": "default category precedence", "source": "S3-STD-008"}
    ],
    "remaining_old_scope": "Shell source and heredoc realization",
    "audience_before": "public help",
    "audience_after": "public help",
    "availability_reduced": false,
    "behavior_changed": false,
    "consequential_delta": false,
    "deltas": {
      "behavior": "none",
      "interface": "none",
      "scientific": "none",
      "exception": "ownership_boundary_only",
      "information": "owner_movement_with_bidirectional_record"
    },
    "evidence_roles": [
      {"source": "S3-STD-008", "role": "immutable_accepted_evidence"},
      {"source": "LLM comparison", "role": "supporting_comparison"}
    ],
    "llm_role": "supporting_comparison",
    "human_authorization": {
      "required": false,
      "status": "not_required",
      "reference": "approved Phase A/B"
    }
  },
  "rejected": {"authorized_by": "LLM review", "approved": true},
  "boundary": {
    "consequential_delta": true,
    "human_authorization": {
      "required": true,
      "status": "pending",
      "reference": ""
    }
  },
  "non_applicable": {"matching_hash": "evidence_only"}
}
EOM


succeed "generated semantic-movement fixtures under ${dir_fix}"
