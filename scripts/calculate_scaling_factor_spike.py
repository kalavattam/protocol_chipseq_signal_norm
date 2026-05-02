#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.

"""
############
## Script ##
############

calculate_scaling_factor_spike.py


#################
## Description ##
#################

Calculates spike-in-based scaling factor(s) for ChIP-seq datasets from counts
of fragments(†) mapped to the main (experimental) and spike-in genomes for a
matched IP (immunoprecipitate) and input pair.

The script’s primary output is the Bio-protocol spike-in scaling coefficient
(‡). It can also output the ChIP-Rx coefficients α^{IP} and α^{in} (Orlando et
al.), the ratio α^{IP} / α^{in}, and the Rx-input coefficient α_{Rx} (Bressan
et al., Ma et al., Niu et al., Fursova et al.), each of which is an algebraic
function of the same counts.

(†) Here, “fragments” can mean either (i) read alignments or (ii) fragments
inferred from alignments (preferred). Either is acceptable as long as the same
counting rules (e.g., filters, secondary/supplementary handling, duplicate
handling, etc.) are applied consistently across IP and input, and across main
and spike-in references.

(‡) 'fractional', 'bioprotocol' and 'bio_protocol' are used here as convenient,
informal labels for the coefficient used in the Bio-protocol manuscript.


###########
## Usage ##
###########

python -m scripts.calculate_scaling_factor_spike \\
    [--help] [--verbose] [--coef <str>] [--fmt <str>] \\
    --main_ip <int> --spike_ip <int> \\
    --main_in <int> --spike_in <int> \\
    [--rnd <int>]


###############
## Arguments ##
###############

 -h, --help
    Return help message and exit.

 -v, --verbose
    Increase output verbosity.

 -c, --coef, --coefficient
    Which coefficient to output (default: fractional).

-ft, --fmt, --format
    Output format (default: plain).

-mp, --main_ip
    Number of “main” alignments in IP sample.

-sp, --spike_ip
    Number of spike-in alignments in IP sample.

-mn, --main_in
    Number of “main” alignments in input sample.

-sn, --spike_in
    Number of spike-in alignments in input sample.

-dp, --dp, --rnd, --round, --decimals, --digits
    Decimal precision for rounding coefficient(s) (default: 24).


%%%%%%%%%%%%%%%%%%%%%%
%% Argument details %%
%%%%%%%%%%%%%%%%%%%%%%

- -c, --coef, --coefficient {fractional,chiprx_alpha_ip,chiprx_alpha_in,
  chiprx_alpha_ratio,rxinput_alpha,all}
    Select which coefficient to output (default: fractional).

    Accepted values and aliases (case-insensitive):
        + fractional | bioprotocol | bio_protocol
            (N_s^{in} / T^{in}) / (N_s^{IP} / T^{IP})
        + chiprx_alpha_ip | alpha_chiprx_ip | chiprx_ip
            10^6 / N_s^{IP}
        + chiprx_alpha_in | alpha_chiprx_in | chiprx_in
            10^6 / N_s^{in}
        + chiprx_alpha_ratio | alpha_chiprx_ratio | chiprx_ratio
            N_s^{in} / N_s^{IP}
        + rxinput_alpha | alpha_rxinput | rxi_alpha | alpha_rxi | rxinput | rxi
            (10^6 * N_s^{in}) / (N_s^{IP} * T^{in})
        + all
            Compute and emit all of the above.

    Notes:
        + The canonical coefficient names emitted by the program are
            fractional, chiprx_alpha_ip, chiprx_alpha_in, chiprx_alpha_ratio,
            rxinput_alpha.
        + Aliases are accepted only for '--coef' input.
        + Additional legacy aliases are accepted for backward compatibility but
          are intentionally not exposed in user-facing documentation.

- -ft, --fmt, --format {plain,tsv,json}
    Select output format (default: plain).

    Contracts:
    + 'fmt=plain'
        stdout: a single floating-point number with no label and a trailing
        newline, rounded using --rnd.

        Allowed only when '--coef' is not 'all'.

    + 'fmt=tsv'
        stdout: a header line followed by one or more data lines.

        Columns (always in this order):
            coef    value

        If '--coef' is not 'all', emit one row. If '--coef' is 'all', emit one
        row per coefficient, in the following order:
            fractional, chiprx_alpha_ip, chiprx_alpha_in, chiprx_alpha_ratio,
            rxinput_alpha.

        Values are rounded using '--rnd'.

    + 'fmt=json'
        stdout: a JSON object mapping coefficient names to values.

        If '--coef' is not 'all', the object contains a single key. If '--coef'
        is 'all', the object contains keys:
            "fractional", "chiprx_alpha_ip", "chiprx_alpha_in",
            "chiprx_alpha_ratio", "rxinput_alpha"

        Values are rounded using '--rnd'.

- Error handling for 'format'/'coef' mismatches
    If '--format plain' is requested with '--coef all', exit non-zero with an
    error:
        Format 'plain' requires a single coefficient; use '--format tsv|json'
        or a specific '--coef' value.

- Rounding
    '--rnd' applies to all numeric outputs in all formats and controls the
    number of decimal places used when printing.


############
## Output ##
############

Output depends on '--coef' and '--format':
    - If '--coef' is a single coefficient (not 'all'):
        + 'fmt=plain':
            stdout: a single floating-point number with no label and a trailing
            newline (rounded to '--rnd' decimals).
        + 'fmt=tsv':
            stdout: a header line followed by one row:
                coef    value
        + 'fmt=json':
            stdout: a JSON object with one key/value pair.
    - If '--coef all':
        + 'fmt=tsv':
            stdout: a header line followed by one row per coefficient, in the
            following order:
                fractional, chiprx_alpha_ip, chiprx_alpha_in,
                chiprx_alpha_ratio, rxinput_alpha
        + 'fmt=json':
            stdout: a JSON object with the following keys:
                "fractional", "chiprx_alpha_ip", "chiprx_alpha_in",
                "chiprx_alpha_ratio", "rxinput_alpha"
    - Format/coef mismatch:
        + '--format plain' is not allowed with '--coef all' and exits non-zero.
    - stderr:
        + With '--verbose', a banner of parsed arguments and computed values is
          printed to stderr.
        + The script may emit warning messages to stderr if spike-in fractions
          are extremely small or large (see “General notes”).

Non-zero exit on error:
    - invalid counts
    - zero totals
    - division by zero
    - invalid '--coef' alias/value
    - invalid '--format' value or '--format plain' with '--coef all'


##############
## Examples ##
##############

1a. Compute the “fractional” coefficient from the Bio-protocol manuscript
    (defaults: '--coef fractional', '--format plain') as a single plain number;
    use short-form arguments
'''bash
python -m scripts.calculate_scaling_factor_spike \
    -mp 100000 -sp 5000 -mn 90000 -sn 4500
'''

1b. Same as the above, but use long-form arguments
'''bash
python -m scripts.calculate_scaling_factor_spike \
    --main_ip 100000 --spike_ip 5000 \
    --main_in 90000  --spike_in 4500
'''

2. Emit a machine-readable JSON object (for a single coefficient)
'''bash
python -m scripts.calculate_scaling_factor_spike \
    --coef rxinput_alpha --format json \
    --main_ip 100000 --spike_ip 5000 \
    --main_in 90000  --spike_in 4500
'''

3. Emit the ratio of ChIP-Rx coefficients
'''bash
python -m scripts.calculate_scaling_factor_spike \
    --coef chiprx_alpha_ratio --format tsv \
    --main_ip 100000 --spike_ip 5000 \
    --main_in 90000  --spike_in 4500
'''

4. Emit all coefficients in TSV (easy to inspect or paste into a sheet)
'''bash
python -m scripts.calculate_scaling_factor_spike \
    --coef all --format tsv \
    --main_ip 100000 --spike_ip 5000 \
    --main_in 90000  --spike_in 4500 \
    --rnd 12
'''

5. Debug: show alias canonicalization (including case insensitivity) in the
   verbose banner
'''bash
python -m scripts.calculate_scaling_factor_spike \
    -v --coef aLpHa_RxInPuT \
    --main_ip 100000 --spike_ip 5000 \
    --main_in 90000  --spike_in 4500 \
    --rnd 12
'''


##########################################
## Scaling factor details, calculations ##
##########################################

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fractional scaling coefficient from Alavattam et al., Bio-protocol 2025 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

For a single ChIP-seq sample, let
    - N_m be the number of fragments from the main genome, and
    - N_s be the number of fragments from the spike-in genome.

Define the total number of fragments in that sample as

    T = N_m + N_s

and define the spike-in fraction (the fraction of fragments from the spike-in
genome) as

    f_s = N_s / T.

For a matched IP/input pair, use superscripts ^{IP} and ^{in}:
    - N_m^{IP}, N_s^{IP}: main and spike-in fragments in the IP sample;
    - N_m^{in}, N_s^{in}: main and spike-in fragments in the input sample.

so that
    - T^{IP} = N_m^{IP} + N_s^{IP},
    - T^{in} = N_m^{in} + N_s^{in},
    - f_s^{IP} = N_s^{IP} / T^{IP}, and
    - f_s^{in} = N_s^{in} / T^{in}.

In this notation, the Bio-protocol spike-in coefficient is the ratio of spike-
in fractions:

    S = f_s^{in} / f_s^{IP}

which expands to

    S = (N_s^{in} / T^{in}) / (N_s^{IP} / T^{IP})
      = (N_s^{in} × T^{IP}) / (N_s^{IP} × T^{in}).

Interpretation:
    - When applied as a multiplicative scale factor to an IP-derived signal
      track, S multiplies every IP bin by the same constant.
    - Because S is a ratio of spike-in fractions, it up-weights IP tracks whose
      spike-in fraction is smaller than input’s, and down-weights IP tracks
      whose spike-in fraction is larger than input’s.

For more details, see the following:
    - https://pubmed.ncbi.nlm.nih.gov/40364978/
    - https://biostars.org/p/9572653/#9572655


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ChIP-Rx coefficient α (Orlando et al., Cell Rep 2014) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Orlando et al. define N_d as the number of spike-in fragments in “millions”
(subscript _d reflects D. melanogaster for the spike-in reference in their
implementation). In the notation above:

    N_d = N_s / 10^6.

They then define the ChIP-Rx normalization factor α as the reciprocal of N_d:

    α = 1 / N_d.

Substituting N_d = N_s / 10^6 gives the commonly used form

    α = 1 / (N_s / 10^6) = 10^6 / N_s.

Applied separately to a matched IP/input pair,
    - α^{IP} = 10^6 / N_s^{IP} and
    - α^{in} = 10^6 / N_s^{in}.

Interpretation:
    α is a “per-million spike-in fragments” scaling: if N_s is small, α is
    large; if N_s is large, α is small.

Relationship of S to α:
    Because α = 10^6 / N_s, the ratio of ChIP-Rx coefficients between IP and
    input is

        α^{IP} / α^{in}
            = (10^6 / N_s^{IP}) / (10^6 / N_s^{in})
            = N_s^{in} / N_s^{IP}.

    Using S = (N_s^{in} × T^{IP}) / (N_s^{IP} × T^{in}), it follows that

         α^{IP} / α^{in} = S × (T^{in} / T^{IP}).

    Interpretation:
        - α^{IP} / α^{in} captures only the relative spike-in depths
          (N_s^{in} / N_s^{IP}).
        - S equals (decomposes to) the relative spike-in depth multiplied by
          the paired-sample depth factor (T^{IP} / T^{in}), which anchors the
          result to the matched input baseline.

For more details, see the following:
    - https://pubmed.ncbi.nlm.nih.gov/25437568/


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ratio of ChIP-Rx coefficients α^{IP} / α^{in} %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

From the ChIP-Rx definitions

    α^{IP} = 10^6 / N_s^{IP}
    α^{in} = 10^6 / N_s^{in},

their ratio is

    α^{IP} / α^{in}
      = (10^6 / N_s^{IP}) / (10^6 / N_s^{in})
      = N_s^{in} / N_s^{IP}.

This script refers to that quantity as 'chiprx_alpha_ratio'.

Relationship to the Bio-protocol coefficient:
    Recall

        fractional
            = (N_s^{in} / T^{in}) / (N_s^{IP} / T^{IP})
            = (N_s^{in} / N_s^{IP}) × (T^{IP} / T^{in}).

    Therefore,

        chiprx_alpha_ratio = fractional × (T^{in} / T^{IP})

    and equivalently,

        fractional = chiprx_alpha_ratio × (T^{IP} / T^{in}).

Interpretation:
    - 'chiprx_alpha_ratio' captures only the relative spike-in depths
      (N_s^{in} / N_s^{IP}).
    - Unlike the Bio-protocol coefficient, it does not include the paired-
      sample depth factor (T^{IP} / T^{in}).
    - This distinction matters when the coefficient is used together with a raw
      per-bin IP/input ratio, because the raw bin ratio already carries the
      paired-sample depth term.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rx-input coefficient α_{Rx} (Bressan et al., NAR Genom Bioinform 2024) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Based on the work of Ma et al., Niu et al., and Fursova et al., Bressan et al.
summarize an “Rx-input” style coefficient that uses:
    - the spike-in count in IP (N_s^{IP}) to set a ChIP-Rx-like scaling, and
    - the spike-in fraction in input (f_s^{in} = N_s^{in} / T^{in}) to anchor
      that scaling to the matched input baseline.

In the notation above, the resulting coefficient can be written as

    α_{Rx} = (10^6 / N_s^{IP}) × (N_s^{in} / T^{in})
           = (10^6 × N_s^{in}) / (N_s^{IP} × T^{in})
           = α^{IP} × f_s^{in}.

Interpretation:
    - The first factor (10^6 / N_s^{IP}) is ChIP-Rx-style “per-million spike-in
      fragments” scaling based on the IP spike-in depth;
    - The second factor (N_s^{in} / T^{in}) uses the matched input’s spike-in
      fraction as an anchor, down-weighting the IP-based scale when the input
      spike-in fraction is small and up-weighting it when that fraction is
      large.

Relationship of α_{Rx} to S:
    Start from the definitions:

        S = (N_s^{in} / T^{in}) / (N_s^{IP} / T^{IP})
          = (N_s^{in} × T^{IP}) / (N_s^{IP} × T^{in}),

    and

        α_{Rx} = (10^6 / N_s^{IP}) × (N_s^{in} / T^{in})
               = (10^6 × N_s^{in}) / (N_s^{IP} × T^{in}).

    Compare the right-hand sides:
        - S has numerator (N_s^{in} × T^{IP}), and
        - α_{Rx} has numerator (10^6 × N_s^{in})

    and both share denominator (N_s^{IP} × T^{in}).

    Therefore,

        α_{Rx} = (10^6 / T^{IP}) × S
               = 10^6 × (S / T^{IP}),

    and after rearranging,

        S = (α_{Rx} / 10^6) × T^{IP}.

    Interpretation:
        - S is dimensionless and scales an IP-derived track directly.
        - α_{Rx} is “S expressed in per-million units,” after additionally
          dividing by the IP library depth T^{IP}.  In other words, α_{Rx}
          encodes the same spike-in/input anchoring as S, but with an explicit
          per-million scaling.

Relationship of α_{Rx} to α:
    Recall:
        - α^{IP} = 10^6 / N_s^{IP}
        - α^{in} = 10^6 / N_s^{in}
        - f_s^{in} = N_s^{in} / T^{in}

    Then,

        α_{Rx} = (10^6 / N_s^{IP}) × (N_s^{in} / T^{in})
               = α^{IP} × f_s^{in}.

    Interpretation:
        α_{Rx} starts with a ChIP-Rx-style scale set by the IP spike-in depth
        (α^{IP}), then down-weights that scale by the spike-in fraction
        observed in the matched input (f_s^{in}).

For more details, see the following:
    - #TODO (Kris will do so)


###################
## General notes ##
###################

- Counts should reflect the same filtering applied to downstream signal (e.g.,
  processed main alignments, consistent FLAG/secondary handling) and should be
  computed the same way for IP and input.
- Prefer read alignment-inferred fragments over raw alignment counts when
  possible; more importantly, keep whatever choice consistent across IP and
  input.
- 'fractional' and 'chiprx_alpha_ratio' are dimensionless. By contrast,
  'chiprx_alpha_ip', 'chiprx_alpha_in', and 'rxinput_alpha' retain the
  conventional per-million scaling inherited from ChIP-Rx/Rx-input
  formulations; this is related to Orlando et al.’s RRPM (“reference reads
  per million”) framing.
- If IP and input have identical spike-in fractions, then 'fractional' = 1. If
  IP’s spike-in fraction is larger than input’s, then 'fractional' < 1 (and
  vice versa).
- The script may emit non-fatal warnings to stderr when spike-in fractions are
  extremely small or large; computation still succeeds when denominators are
  valid.
- Some coefficients require non-zero spike-in counts. In particular,
  N_s^{IP} > 0 is required for 'fractional', 'chiprx_alpha_ip',
  'chiprx_alpha_ratio', and 'rxinput_alpha'; N_s^{in} > 0 is required for
  'chiprx_alpha_in' and 'chiprx_alpha_ratio'. The program exits non-zero when a
  requested coefficient would require division by zero.
- Some extreme count compositions (e.g., 'main == 0' with 'spike > 0', or vice
  versa) are mathematically valid and are allowed, but the script emits a
  warning because they often indicate a counting problem or a main/spike-in
  genome partitioning problem.


#######################
## Performance notes ##
#######################

- O(1) arithmetic.
- Since there’s no file I/O, input sizes do not affect runtime.
    + Runtime: argument parsing, arithmetic, and printing.
- Floating-point is double precision.
- Formatting controls final string rounding.
- Given the same inputs, output is deterministic.


##########
## TODO ##
##########

- Future experiment: refine the current spike-in fraction warnings by
  empirically deriving “expected range” bounds, and consider warning on
  spike/main ratios instead of or in addition to spike/total fractions (or let
  the user decide depending on analysis approach, counting definition, etc.).
    + For example,
        '''python
        #  Check for very high or low ratios
        if not (0.1 <= ratio_ip <= 10):
            print(
                f'Warning: The ratio of spike-in to "main" alignments for IP '
                f'data ({ratio_ip:.2f}) may be outside the expected range.'
            )

        if not (0.1 <= ratio_in <= 10):
            print(
                'Warning: The ratio of spike-in to "main" alignments for '
                f'input data ({ratio_in:.2f}) may be outside the expected '
                'range.'
            )
        '''
    + Rough, brief plan:
        - collect several IP/input pairs with known spike-in rates spanning low
          to moderate (e.g., ~0.05%–5%),
        - compute 'ratio_ip = spike_ip / (main_ip + spike_ip)' and
          'ratio_in = spike_in / (main_in + spike_in)',
        - summarize distributions across runs/replicates, and
        - choose non-blocking warning thresholds (e.g., lower/upper quantiles).
        - Decide whether to warn on fractions (spike/total; 0–1) or on ratios
          (spike/main; unbounded), then
        - implement optional '--warn-bounds ip_low ip_high' in the CLI to allow
          per-lab overrides.
"""

from __future__ import annotations

import argparse
import json
import re
import signal
import sys

from contextlib import redirect_stdout

from scripts.functions.utils_check import check_cmp
from scripts.functions.utils_cli import (
    add_help_cap,
    CapArgumentParser
)

try:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except (AttributeError, ValueError):
    pass

assert sys.version_info >= (3, 10), "Python >= 3.10 required."


#  Map user-facing '--coef' aliases to canonical coefficient names
COEF_ALIAS_CANON = {
    "fractional": "fractional",
    "bioprotocol": "fractional",
    "bio_protocol": "fractional",
    "alavattam": "fractional",                   # Not exposed in docs
    "tsukiyama": "fractional",                   # Not exposed in docs
    "s": "fractional",                           # For b/c; not exposed in docs

    "chiprx_alpha_ip": "chiprx_alpha_ip",
    "alpha_chiprx_ip": "chiprx_alpha_ip",
    "alpha_ip": "chiprx_alpha_ip",               # For b/c; not exposed in docs
    "chiprx_ip": "chiprx_alpha_ip",
    "orlando_ip": "chiprx_alpha_ip",             # Not exposed in docs

    "chiprx_alpha_in": "chiprx_alpha_in",
    "alpha_chiprx_in": "chiprx_alpha_in",
    "alpha_in": "chiprx_alpha_in",               # For b/c; not exposed in docs
    "chiprx_in": "chiprx_alpha_in",
    "orlando_in": "chiprx_alpha_in",             # Not exposed in docs

    "chiprx_alpha_ratio": "chiprx_alpha_ratio",
    "alpha_chiprx_ratio": "chiprx_alpha_ratio",
    "alpha_ratio": "chiprx_alpha_ratio",         # For consistency; not exposed
    "chiprx_ratio": "chiprx_alpha_ratio",
    "orlando_ratio": "chiprx_alpha_ratio",       # Not exposed in docs
    "r": "chiprx_alpha_ratio",                   # Not exposed in docs

    "rxinput_alpha": "rxinput_alpha",
    "alpha_rxinput": "rxinput_alpha",
    "rxi_alpha": "rxinput_alpha",
    "alpha_rxi": "rxinput_alpha",
    "rxinput": "rxinput_alpha",
    "rxi": "rxinput_alpha",
    "bressan": "rxinput_alpha",                  # Not exposed in docs
    "ma": "rxinput_alpha",                       # Not exposed in docs
    "niu": "rxinput_alpha",                      # Not exposed in docs
    "fursova": "rxinput_alpha",                  # Not exposed in docs

    "all": "all",
}

#  Canonical output order used when '--coef' all
COEF_ORDER = (
    "fractional", "chiprx_alpha_ip", "chiprx_alpha_in", "chiprx_alpha_ratio",
    "rxinput_alpha"
)


def normalize_coef(raw: str) -> str:
    """
    Normalize a user-supplied coefficient name or alias to its canonical name.
    """
    key = raw.strip().lower()

    #  Treat hyphens like underscores (and collapse runs)
    key = re.sub(r"[-_]+", "_", key)

    if key in COEF_ALIAS_CANON:
        return COEF_ALIAS_CANON[key]

    msg = (
        f"Invalid --coef '{raw}'.\n\n"
        "Accepted values and aliases (case-insensitive):\n"
        "    - fractional | bioprotocol | bio_protocol\n"
        "    - chiprx_alpha_ip | alpha_chiprx_ip | chiprx_ip\n"
        "    - chiprx_alpha_in | alpha_chiprx_in | chiprx_in\n"
        "    - chiprx_alpha_ratio | alpha_chiprx_ratio | chiprx_ratio\n"
        "    - rxinput_alpha | alpha_rxinput | rxi_alpha | alpha_rxi | "
        "rxinput | rxi\n"
        "    - all\n"
    )
    raise ValueError(msg)


def round_value(x: float, rnd: int) -> float:
    """
    Round a value to the requested number of decimal places, avoiding -0.0.
    """
    y = round(x, rnd)
    if y == 0.0:
        return 0.0
    return y


def emit_output(
    vals: dict[str, float],
    coef: str,
    fmt: str,
    rnd: int
) -> None:
    """
    Emit one or more computed coefficients in plain, TSV, or JSON format.
    """
    if fmt == "plain" and coef == "all":
        raise ValueError(
            "Format 'plain' requires a single coefficient; use "
            "'--format tsv|json' or a specific '--coef' value."
        )

    if coef != "all" and coef not in vals:
        raise ValueError(
            f"Requested coefficient '{coef}' was not computed. "
            f"Computed keys: {', '.join(sorted(vals))}."
        )

    if coef != "all":
        key = coef
        v = round_value(vals[key], rnd)
        if fmt == "plain":
            print(f"{v:.{rnd}f}")
            return
        if fmt == "tsv":
            print("coef\tvalue")
            print(f"{key}\t{v:.{rnd}f}")
            return
        if fmt == "json":
            print(json.dumps({key: v}, sort_keys=False))
            return

    # coef == all
    if fmt == "tsv":
        print("coef\tvalue")
        for k in COEF_ORDER:
            v = round_value(vals[k], rnd)
            print(f"{k}\t{v:.{rnd}f}")
        return

    if fmt == "json":
        obj = {k: round_value(vals[k], rnd) for k in COEF_ORDER}
        print(json.dumps(obj, sort_keys=False))
        return

    raise ValueError(f"Unknown format '{fmt}'.")


def validate_counts(args: argparse.Namespace) -> None:
    """
    Ensure counts are valid and totals are non-zero.

    Args:
        args : argparse.Namespace
            main_ip, spike_ip, main_in, spike_in

    Raises:
        ValueError
            If any count < 0 or if rnd < 0.
        ZeroDivisionError
            If a per-sample total is zero.
    """
    for name in ("main_ip", "spike_ip", "main_in", "spike_in"):
        check_cmp(getattr(args, name), "ge", 0, name, allow_none=False)

    check_cmp(args.rnd, "ge", 0, "rnd", allow_none=False)

    #  Preserve raw user input before canonicalizing
    args.coef_raw = args.coef

    #  Canonicalize coef (aliases to canonical)
    args.coef = normalize_coef(args.coef)

    #  (Can’t be triggered, but keep it anyway)
    if args.fmt not in ("plain", "tsv", "json"):
        raise ValueError(f"Invalid --format '{args.fmt}'.")

    if args.fmt == "plain" and args.coef == "all":
        raise ValueError(
            "Format 'plain' requires a single coefficient; use "
            "'--format tsv|json' or a specific '--coef' value."
        )

    if (args.main_ip + args.spike_ip) == 0:
        raise ZeroDivisionError(
            "IP totals are zero ('main_ip' + 'spike_ip' = 0)."
        )
    if (args.main_in + args.spike_in) == 0:
        raise ZeroDivisionError(
            "Input totals are zero ('main_in' + 'spike_in' = 0)."
        )


def calculate_scaling_factors(
    main_ip: int, spike_ip: int,
    main_in: int, spike_in: int,
    required: tuple[str, ...] = COEF_ORDER
) -> dict[str, float]:
    """
    Compute spike-count algebra coefficients from main/spike counts in
    IP/input.

    Args:
        main_ip, spike_ip, main_in, spike_in
            Non-negative integer counts.
        required
            Tuple of canonical coefficient names to compute. Canonical names:
            fractional, chiprx_alpha_ip, chiprx_alpha_in, chiprx_alpha_ratio,
            rxinput_alpha.

    Returns:
        dict[str, float]
            Keys are exactly those listed in 'required' (in any order).
    """
    #  Validate that all inputs are integers and are non-negative
    for name, count in (
        ("main_ip", main_ip), ("spike_ip", spike_ip),
        ("main_in", main_in), ("spike_in", spike_in)
    ):
        if not isinstance(count, int):
            raise TypeError(
                f"Expected type 'int' for '{name}', but got "
                f"'{type(count).__name__}'."
            )
        if count < 0:
            raise ValueError(
                f"Count for '{name}' must be >= 0, but got '{count}'."
            )

    t_ip = main_ip + spike_ip
    t_in = main_in + spike_in

    if t_ip == 0:
        raise ZeroDivisionError(
            "IP totals are zero ('main_ip' + 'spike_ip' = 0)."
        )
    if t_in == 0:
        raise ZeroDivisionError(
            "Input totals are zero ('main_in' + 'spike_in' = 0)."
        )

    #  Validate 'required' against the allowed canonical coefficient names
    req = set(required)
    unknown = req - set(COEF_ORDER)
    if unknown:
        raise ValueError(
            "Invalid coefficient name(s) in 'required': "
            f"{', '.join(sorted(unknown))}. Allowed: "
            f"{', '.join(COEF_ORDER)}."
        )

    #  Denominator checks only for requested coefficients
    if req & {
        "fractional", "chiprx_alpha_ip", "chiprx_alpha_ratio",
        "rxinput_alpha"
    }:
        if spike_ip == 0:
            raise ZeroDivisionError(
                "'spike_ip' is 0; cannot compute requested coefficient(s) "
                "that require division by N_s^IP."
            )

    if req & {"chiprx_alpha_in", "chiprx_alpha_ratio"}:
        if spike_in == 0:
            raise ZeroDivisionError(
                "'spike_in' is 0; cannot compute requested coefficient(s) "
                "that require division by N_s^in."
            )

    vals: dict[str, float] = {}

    #  Bio-protocol coefficient (Alavattam et al.)
    if "fractional" in req:
        #  fractional = (spike_in / t_in) / (spike_ip / t_ip)
        vals["fractional"] = (spike_in / t_in) / (spike_ip / t_ip)

    #  ChIP-Rx alpha for IP (Orlando et al.)
    if "chiprx_alpha_ip" in req:
        vals["chiprx_alpha_ip"] = 1e6 / spike_ip

    #  ChIP-Rx alpha for input (Orlando et al.)
    if "chiprx_alpha_in" in req:
        vals["chiprx_alpha_in"] = 1e6 / spike_in

    #  Ratio of ChIP-Rx alpha coefficients
    if "chiprx_alpha_ratio" in req:
        #  chiprx_alpha_ratio = (10^6 / spike_ip) / (10^6 / spike_in)
        #                     = spike_in / spike_ip
        vals["chiprx_alpha_ratio"] = spike_in / spike_ip

    #  Rx-input alpha (Bressan et al.)
    if "rxinput_alpha" in req:
        #  rxinput_alpha = (1e6 × spike_in) / (spike_ip × t_in)
        vals["rxinput_alpha"] = (1e6 * spike_in) / (spike_ip * t_in)

    return vals


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command line arguments.
    """
    parser = CapArgumentParser(description=(
        "Calculate one or more spike-in-derived uniform scaling coefficients "
        "for a ChIP-seq IP/input pair."
    ))
    add_help_cap(parser)
    parser.add_argument(
        "-v", "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Increase output verbosity.\n\n"
    )
    parser.add_argument(
        "-c", "--coef", "--coefficient",
        dest="coef",
        type=str,
        required=False,
        default="fractional",
        help=(
            "Which coefficient to output (case-insensitive; default: "
            "%(default)s).\n\n"
            "Aliases:\n"
            "    - fractional | bioprotocol | bio_protocol\n"
            "        (N_s^{in} / T^{in}) / (N_s^{IP} / T^{IP})\n"
            "    - chiprx_alpha_ip | alpha_chiprx_ip | chiprx_ip\n"
            "        10^6 / N_s^{IP}\n"
            "    - chiprx_alpha_in | alpha_chiprx_in | chiprx_in\n"
            "        10^6 / N_s^{in}\n"
            "    - chiprx_alpha_ratio | alpha_chiprx_ratio | chiprx_ratio\n"
            "        N_s^{in} / N_s^{IP}\n"
            "    - rxinput_alpha | alpha_rxinput | rxi_alpha |\n"
            "      alpha_rxi | rxinput | rxi\n"
            "        (10^6 * N_s^{in}) / (N_s^{IP} * T^{in})\n"
            "    - all\n"
            "\n"
            "Associated publications:\n"
            "    - fractional\tAlavattam et al.\n"
            "    - chiprx_alpha_ip\tOrlando et al.\n"
            "    - chiprx_alpha_in\tOrlando et al.\n"
            "    - rxinput_alpha\tNiu et al., Ma et al., Fursova et al., "
            "Bressan et al.\n\n"
        )
    )
    parser.add_argument(
        "-ft", "--fmt", "--format",
        dest="fmt",
        type=str,
        choices=("plain", "tsv", "json"),
        default="plain",
        help=(
            "Output format (default: %(default)s). 'plain' requires a single "
            "coefficient (i.e., not 'all').\n\n"
        )
    )
    parser.add_argument(
        "-mp", "-mip", "--mip", "--main_ip",
        dest="main_ip",
        type=int,
        required=True,
        help=(
            "Number of “main” alignments (reads or read-inferred fragments) "
            "for the ChIP-seq IP data.\n\n"
        )
    )
    parser.add_argument(
        "-sp", "-sip", "--sip", "--spike_ip",
        dest="spike_ip",
        type=int,
        required=True,
        help="Number of spike-in alignments for the ChIP-seq IP data.\n\n"
    )
    parser.add_argument(
        "-mn", "-min", "--min", "--main_in",
        dest="main_in",
        type=int,
        required=True,
        help=(
            "Number of “main” alignments for the corresponding ChIP-seq input "
            "data.\n\n"
        )
    )
    parser.add_argument(
        "-sn", "-sin", "--sin", "--spike_in",
        dest="spike_in",
        type=int,
        required=True,
        help=(
            "Number of spike-in alignments for the corresponding ChIP-seq "
            "input data.\n\n"
        )
    )
    parser.add_argument(
        "-dp", "--dp", "--rnd", "--round", "--decimals", "--digits",
        dest="rnd",
        type=int,
        default=24,
        required=False,
        help=(
            "Number of decimal places for rounding emitted coefficient "
            "values (default: %(default)s).\n\n"
        )
    )

    #  If no arguments are provided, use 'argv' to display help and exit
    argv_parse = sys.argv[1:] if argv is None else argv
    if not argv_parse:
        parser.print_help(sys.stderr)
        sys.exit(0)

    return parser.parse_args(argv_parse)


def main(argv: list[str] | None = None) -> int:
    """
    Execute the primary control flow for the script.

    Parse arguments, validate counts, compute one or more spike-in scaling
    coefficients, and print the result to stdout.

    Args:
        argv: list[str] | None
            Optional command-line argument list. If None, arguments are read
            from 'sys.argv[1:]'.

    Returns:
        int
            Exit status code. On success, prints the requested coefficient
            output and returns 0.

    Side effects:
        Emits warnings (e.g., unusually small or large IP or input ratios) to
        stderr. Prints human-readable error messages to stderr on failure.

    Exits:
        0 on success or when showing help with no arguments, 1 on validation or
        computation errors (e.g., negative counts, zero per-sample totals, or a
        requested division by zero).
    """
    #  Parse CLI arguments
    args = parse_args(argv)

    #  Compute requested coefficient(s), handling exceptions as necessary
    try:
        #  Check inputs
        validate_counts(args)

        #  Return “verbose banner” to stderr
        if args.verbose:
            with redirect_stdout(sys.stderr):
                print("####################################################")
                print("## Arguments for 'calculate_scaling_factor_spike' ##")
                print("####################################################")
                print("")
                print("--verbose")
                print(f"--coef     {args.coef_raw}  (canon: {args.coef})")
                print(f"--format   {args.fmt}")
                print(f"--main_ip  {args.main_ip}")
                print(f"--spike_ip {args.spike_ip}")
                print(f"--main_in  {args.main_in}")
                print(f"--spike_in {args.spike_in}")
                print(f"--rnd      {args.rnd}")
                print("")
                print("")

        #  Compute fractions for optional warnings and verbose echo
        ratio_ip = args.spike_ip / (args.main_ip + args.spike_ip)
        ratio_in = args.spike_in / (args.main_in + args.spike_in)

        #  Additional “implausible composition” warnings (non-fatal); these are
        #  mathematically valid but usually indicate a counting or main/spike-
        #  in partitioning problem
        for label, main, spike, frac in (
            ("IP", args.main_ip, args.spike_ip, ratio_ip),
            ("input", args.main_in, args.spike_in, ratio_in),
        ):
            if main == 0 and spike > 0:
                print(
                    f"Warning: '{label}' has 'main == 0' and 'spike > 0' "
                    f"(spike-in fraction = {frac:.3f}). This is unusual and "
                    "may indicate a counting issue, main/spike-in "
                    "partitioning problem, etc.",
                    file=sys.stderr,
                )
            if main > 0 and spike == 0:
                print(
                    f"Warning: '{label}' has 'spike == 0' and 'main > 0' "
                    "(no spike-in signal). This is unusual for spike-in "
                    "experiments and may indicate a counting issue, "
                    "main/spike-in partitioning problem, etc.",
                    file=sys.stderr,
                )

        #  Warn if spike-in fraction is vanishingly small or oddly large
        for label, frac in (("IP", ratio_ip), ("input", ratio_in)):
            if frac < 1e-6:
                print(
                    f"Warning: '{label}' spike-in fraction may be unusually "
                    f"low ({frac:.3e}).",
                    file=sys.stderr
                )
            elif frac > 0.5:
                print(
                    f"Warning: '{label}' spike-in fraction may be unusually "
                    f"high ({frac:.3f}).",
                    file=sys.stderr
                )

        required = COEF_ORDER if args.coef == "all" else (args.coef,)
        vals = calculate_scaling_factors(
            args.main_ip, args.spike_ip, args.main_in, args.spike_in,
            required=required
        )

        if args.verbose:
            with redirect_stdout(sys.stderr):
                print("###############################################")
                print("## Computed spike fractions and coefficients ##")
                print("###############################################")
                print("")

                #  Keep one shared “value column” start for everything in this
                #  block; ratio labels are “IP spike-in fraction” and ”input
                #  spike-in fraction”
                base_lbls = ("IP spike-in fraction", "input spike-in fraction")
                coef_lbls = COEF_ORDER if args.coef == "all" else (args.coef,)

                all_lbls = base_lbls + tuple(coef_lbls)
                width_lbl = max(len(s) for s in all_lbls)

                def print_dotted(lbl: str, value: str) -> None:
                    dots = "." * (width_lbl - len(lbl) + 8)
                    print(f"{lbl} {dots} {value}")

                print_dotted(
                    "IP spike-in fraction", f"{ratio_ip:.{args.rnd}f}"
                )
                print_dotted(
                    "input spike-in fraction", f"{ratio_in:.{args.rnd}f}"
                )

                for k in coef_lbls:
                    print_dotted(k, f"{vals[k]:.{args.rnd}f}")

                print("")
                print("")

        emit_output(vals, args.coef, args.fmt, args.rnd)
        return 0

    except (ValueError, TypeError, ZeroDivisionError) as e:
        raise SystemExit(str(e))


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:
        try:
            sys.stdout.close()
        except Exception:
            pass
        try:
            sys.stderr.close()
        except Exception:
            pass
        sys.exit(0)
