
# Input-floor command reference and design
`compute_input_floor` computes `dep_min`, a floor for positive input denominators that helps avoid extreme or erroneous IP/input ratios. The current repository consumer is `compute_signal_ratio`, which applies the result as `denominator := max(denominator, dep_min)`.

This document is maintained repository documentation for command-reference, scientific, design, compatibility, and provenance detail. It is not included in the installed Python distribution. Installed command users can rely on `compute_input_floor --help`, which independently documents the complete runtime contract.

Modes `frag` and `norm` compute the fragment-normalization and normalized-coverage floors used in the Dickson/siQ-ChIP and *Bio-protocol* workflows. Those floors were designed for large, sparse genomes such as human and mouse and are likely suitable for most of them&mdash;but not for small, dense genomes such as *S. cerevisiae* and *S. pombe*. Mode `dist` (the default mode) derives the floor from the input itself and suits both regimes.

<br />

## Inputs and modes
The callable accepts `fil_in`, `siz_bin`, `siz_gen`, `mode`, `method`, `qntl_nz`, `coef`, `floor`, `eps`, `mode_nz`, `paired_flags`, `single_flags`, `skp_pfx`, `fmt_in`, and `ref_fa`. The CLI exposes the corresponding options without changing the callable calculation.

| Mode   | Input and calculation                                                                                                     |
| :---   | :---                                                                                                                      |
| `dist` | Requires bedGraph, bdg, or bg input, optionally with `.gz`, and uses column four. It does not use `siz_bin` or `siz_gen`. |
| `frag` | Requires BAM, CRAM, or BED/BED.GZ alignment records and uses their counted-record total with `siz_bin` and `siz_gen`.     |
| `norm` | Uses only `siz_bin` and `siz_gen`; it does not use `fil_in` or the input-format options.                                  |

For `dist` and `frag`, `fil_in='-'` reads standard input and requires `fmt_in`. The public explicit choices are `bam`, `cram`, `bed`, `bedGraph`, `bdg`, and `bg`. Hints accept any letter case and resolve once to the canonical internal values `bam`, `cram`, `bed`, and `bedGraph`. Thus existing aliases such as a lowercase `bedgraph` hint remain accepted even though help displays only `bedGraph`. A named path determines its format case-insensitively from its suffix, so `fmt_in` does nothing for it. CRAM decoding requires `ref_fa`; other formats do not use it.

`siz_bin` is the target signal-bin width in base pairs and `siz_gen` is the effective genome size in base pairs. Both dimensions must be positive and `siz_bin` must be smaller than `siz_gen` in `frag` and `norm`. Dimension validation does not apply to `dist`, because neither dimension enters its statistic and the command does not infer or validate bedGraph interval widths. The callable takes both dimensions positionally and does not use them in `dist`.

For BAM and CRAM counting in `frag`, `flags_pe` and `flags_se` select paired-end and single-end main alignments. `None` uses the defaults. The canonical public CLI spellings are `--flags_pe` and `--flags_se`; the former `--flags-pe` and `--flags-se` spellings remain accepted as hidden compatibility aliases. Help, diagnostics, verbose reporting, documentation, and new examples use only the canonical underscore spellings. The CLI accepts decimal and hexadecimal FLAGs and defaults to `99,1123,163,1187` for paired-end data and `0,16,1024,1040` for single-end data. BED input does not use these FLAG selections, so they are refused with it.

Under [`HELP.PARAMETER.APPLICABILITY`](../standards/help.md), the CLI never ignores a given option silently. It refuses an option that would change the result where it does not apply: `--method`, `--qntl_nz`, `--coef`, `--eps`, `--mode_nz`, and `--floor` outside `dist`; `--siz_bin` and `--siz_gen` in `dist`; `--flags_pe` and `--flags_se` outside `frag` or with BED input; `--fil_in` in `norm`; `--qntl_nz` with a method other than `qntl_nz`; and `--coef` with `qntl_nz`. It ignores a supporting input with a warning, without checking it or passing it on: `--fmt_in` with a named path or in `norm`, `--ref_fa` without CRAM input, and `--skp_pfx` in `norm` or with BAM or CRAM input. The callable raises `ValueError` for `coef` outside `dist` or with `qntl_nz` and for FLAG allowlists outside `frag` or with BED input, and warns about `fmt_in` with a named path and `ref_fa` without CRAM input; parameters with non-`None` defaults cannot show whether they were given, so it does not check them.

`skp_pfx` contains prefixes skipped as header or metadata rows in BED and bedGraph-like inputs. The CLI default is `#,track,browser`; an empty string disables prefix skipping.

<br />

## Distribution calculation
Distribution mode reads bedGraph column-four values and applies this order:
1. Skip configured header and metadata rows and non-data rows.
2. Retain finite values.
3. Apply the selected `eps` and `mode_nz` rule.
4. Retain only positive values, `v_i > 0`.
5. Apply the selected `method`.
6. Apply the lower bound as `dep_min := max(dep_min, floor)`.

The epsilon rules are exact. `closed` drops values with `abs(v_i) <= eps`, `open` drops values with `abs(v_i) < eps`, and `off` disables epsilon filtering. Positive-only filtering follows in every case because ratio denominators occupy a positive domain; retaining positive bins keeps `dep_min` in the same conceptual space as the denominators it floors.

The four methods operate on the filtered positive values:

| Method       | Calculation                                                                                                                                                                   |
| :---         | :---                                                                                                                                                                          |
| `qntl_nz`    | `dep_min = Q_q({v_i : v_i > 0})`, where `q = qntl_nz / 100`. Sort the `N` values and select `i = floor(q * (N - 1))`, clamped to `[0, N - 1]`, so `dep_min = sorted_vals[i]`. |
| `frc_mdn_nz` | `dep_min = coef * median({v_i : v_i > 0})`.                                                                                                                                   |
| `frc_avg_nz` | `dep_min = coef * mean({v_i : v_i > 0})`.                                                                                                                                     |
| `min_nz`     | `dep_min = coef * min({v_i : v_i > 0})`.                                                                                                                                      |

`qntl_nz` is a percentage in `[0, 100]` and is used only by `qntl_nz`. `coef` is nonnegative and is used only by `frc_mdn_nz`, `frc_avg_nz`, and `min_nz`. When `coef` is omitted, its defaults match `compute_pseudo.py`: `0.01` for `frc_*` and `1.0` for `min_nz`. The `qntl_nz` method does not use `coef`. `floor` and `eps` are nonnegative and are used only by `dist`.

<br />

## Fragment and normalized calculations
Fragment mode computes:
```text
dep_min = ((n * b) / g) / [1 - (b / g)]
```
Here, `n` is the number of counted BAM, CRAM, or BED/BED.GZ fragment (or alignment) records, `b = siz_bin`, and `g = siz_gen`. This matches the fragment-normalized signal derivation used in the siQ-ChIP code paths.

Normalized mode computes:
```text
dep_min = (b / g) / [1 - (b / g)]
```
This is the per-record fragment expression after dividing by `n`:
```text
([(n * b) / g] / [1 - (b / g)]) / n = (b / g) / [1 - (b / g)]
```
This matches the normalized-coverage track that `compute_signal --method norm` writes. Each fragment's overlap with a bin is divided by the fragment's length and by the total fragment count, so each fragment contributes `1 / n` in all and the bin values sum to about `1`: the track is a distribution of mass over bins (a probability mass function over bins; its running sum along the genome, not the track itself, is the cumulative distribution function). With `g / b` bins, the mean bin value is then about `b / g`, whatever the depth, which is why `norm` needs no input file. The sum falls slightly short of `1` through clipping at chromosome ends and floating-point accumulation, and any later scaling or decimal rounding changes it further.

This normalized-coverage derivation and its siQ-ChIP provenance come from Brad Dickson's workflow, designed for large, sparse genomes such as human and mouse. The `frag` and `norm` modes keep those calculations; they are likely suitable for most large, sparse genomes, but not for small, dense ones.

The formulas compute the mean per-bin input depth correctly for any genome, depth, or bin width. What does not transfer is using that mean as the floor, which assumes a sparse regime. A bin's relative sampling noise is roughly `1 / sqrt(λ)`, where `λ ≈ N (ℓ + b) / g` is the expected number of fragments overlapping it and `ℓ` is the fragment length. In a large genome at typical sequencing depth, `λ` is small, a bin well below the mean is within sampling noise of it, and clamping it to the mean mostly removes noise. In a small, dense genome sequenced to comparable depth, `λ` is one to two orders of magnitude larger, so a bin well below the mean is more likely a real depletion, and clamping overwrites it. `dist` instead takes the floor from the measured input distribution, so it adapts to either regime.

Limited initial testing supports this, and the figures below are not yet general. On a few *S. cerevisiae* libraries at `b = 10`, `λ` was about 240. For comparison, Dickson's code is set up for human data ([`README.md`](https://github.com/BradleyDickson/siQ-ChIP/blob/24e601beb373479e1e78ff907850a3527c4d9e3e/README.md#L274)) with `b = 30` ([`runCrunch.sh`, line 1](https://github.com/BradleyDickson/siQ-ChIP/blob/24e601beb373479e1e78ff907850a3527c4d9e3e/runCrunch.sh#L1)) and `g = 3.2e9` ([line 21](https://github.com/BradleyDickson/siQ-ChIP/blob/24e601beb373479e1e78ff907850a3527c4d9e3e/runCrunch.sh#L21), where the floor is commented "average layer on input"); a comment in [`mergetracks.f90`, line 62](https://github.com/BradleyDickson/siQ-ChIP/blob/24e601beb373479e1e78ff907850a3527c4d9e3e/mergetracks.f90#L62) gives a floor of 7 as the estimate for a depth of 100 million, which implies about 225 bp of coverage per fragment, and these give `λ` of about 8. In the *S. cerevisiae* libraries, the mean-layer floor clamped a large share of bins, roughly 40-93% depending on the library, in normalized-coverage ratios. These come from a small number of datasets and need confirmation across many datasets and organisms, including *S. pombe*, before they are treated as general.

<br />

## Examples
Compute a first-percentile distribution floor from the positive values in a bedGraph track:
```bash
compute_input_floor \
    --mode dist \
    --fil_in signal.bdg \
    --method qntl_nz \
    --qntl_nz 1
```

Compute the normalized floor from explicit dimensions:
```bash
compute_input_floor \
    --mode norm \
    --siz_bin 30 \
    --siz_gen 12157105
```

<br />

## Defaults, output, and failures
The CLI defaults are `mode='dist'` and `dp=24`. The `dist` options default to `method='qntl_nz'`, `qntl_nz=1.0`, `floor=0.0`, `eps=0.0`, and `mode_nz='closed'`, and the `frag` and `norm` dimensions default to `siz_bin=10` and `siz_gen=12157105`. A default takes effect only in the modes its option applies to, and only an option given explicitly is refused or warned about. The `siz_gen` default is appropriate for *Saccharomyces cerevisiae* when retaining multi-mapping alignments.

`compute_input_floor()` returns one unrounded floating-point `dep_min`. Distribution mode returns the selected statistic after its lower bound; fragment and normalized modes return their bin-to-genome-size formulas. The CLI uses the shared `utilities.utils_format.format_value()` behavior: it emits finite values with at most `dp` decimal places, removes non-informative trailing zeros and a trailing decimal point, and normalizes negative zero to `0`.

The callable raises `InputFloorValidationError` for an invalid mode, invalid `frag`/`norm` dimensions, an invalid input format, an empty filtered positive distribution, or a non-finite computed floor. It raises `AlignmentReadError` when an alignment input cannot be read. CLI validation also requires a finite `qntl_nz` in `[0, 100]`, nonnegative `coef`, `floor`, and `eps` values when applicable, and a nonnegative `dp`. Anticipated validation and computation failures retain their actionable stderr diagnostics and status `1`; parser-controlled usage errors retain argparse behavior.
