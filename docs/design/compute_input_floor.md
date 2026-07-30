
# Input-floor command reference and design
`compute_input_floor` computes `dep_min`, a floor for positive input denominators that helps avoid extreme or erroneous IP/input ratios. The current repository consumer is `compute_signal_ratio`, which applies the result as `denominator := max(denominator, dep_min)`.

This document is maintained repository documentation for command-reference, scientific, design, compatibility, and provenance detail. It is not included in the installed Python distribution. Installed command users can rely on `compute_input_floor --help`, which independently documents the complete runtime contract.

Use `dist` for new analyses. Use `frag` or `norm` when reproducing the fragment-normalization or normalized-coverage floor calculations used in the Dickson/siQ-ChIP and *Bio-protocol* workflows.

<br />

## Inputs and modes
The callable accepts `fil_in`, `siz_bin`, `siz_gen`, `mode`, `method`, `qntl_nz`, `coef`, `floor`, `eps`, `mode_nz`, `paired_flags`, `single_flags`, `skp_pfx`, `infmt`, and `ref_fa`. The CLI exposes the corresponding options without changing the callable calculation.

| Mode   | Input and calculation                                                                                                 |
| :---   | :---                                                                                                                  |
| `dist` | Requires bedGraph, bdg, or bg input, optionally with `.gz`, and uses column four. It ignores `siz_bin` and `siz_gen`. |
| `frag` | Requires BAM, CRAM, or BED/BED.GZ alignment records and uses their counted-record total with `siz_bin` and `siz_gen`. |
| `norm` | Uses only `siz_bin` and `siz_gen`; it ignores `fil_in` and all input-format options.                                  |

For `dist` and `frag`, `fil_in='-'` reads standard input and requires `infmt`. The public explicit choices are `bam`, `cram`, `bed`, `bedGraph`, `bdg`, and `bg`. Hints accept any letter case and resolve once to the canonical internal values `bam`, `cram`, `bed`, and `bedgraph`. Thus existing aliases such as a lowercase `bedgraph` hint remain accepted even though help displays only `bedGraph`. A named path determines its format case-insensitively from its suffix, so `infmt` is ignored. CRAM decoding requires `ref_fa`; other formats ignore it.

`siz_bin` is the target signal-bin width in base pairs and `siz_gen` is the effective genome size in base pairs. Both dimensions must be positive and `siz_bin` must be smaller than `siz_gen` in `frag` and `norm`. Dimension validation does not apply to `dist`, because neither dimension enters its statistic and the command does not infer or validate bedGraph interval widths.

For BAM and CRAM counting in `frag`, `flags_pe` and `flags_se` select paired-end and single-end main alignments. `None` uses the defaults. The canonical public CLI spellings are `--flags_pe` and `--flags_se`; the former `--flags-pe` and `--flags-se` spellings remain accepted as hidden compatibility aliases. Help, diagnostics, verbose reporting, documentation, and new examples use only the canonical underscore spellings. The CLI accepts decimal and hexadecimal FLAGs and defaults to `99,1123,163,1187` for paired-end data and `0,16,1024,1040` for single-end data. BED input ignores these FLAG selections.

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

`qntl_nz` is a percentage in `[0, 100]` and is used only by `qntl_nz`. `coef` is nonnegative and is used only by the other three methods. When `coef` is omitted, its defaults match `compute_pseudo.py`: `0.01` for `frc_*` and `1.0` for `min_nz`. The `qntl_nz` method ignores `coef`. `floor` and `eps` are nonnegative and are used only by `dist`.

<br />

## Fragment and normalized calculations
Fragment mode computes:
```text
dep_min = ((n * b) / g) / [1 - (b / g)]
```
Here, `n` is the number of counted BAM, CRAM, or BED/BED.GZ alignment records, `b = siz_bin`, and `g = siz_gen`. This matches the fragment-normalized signal derivation used in the siQ-ChIP code paths.

Normalized mode computes:
```text
dep_min = (b / g) / [1 - (b / g)]
```
This is the per-record fragment expression after dividing by `n`:
```text
([(n * b) / g] / [1 - (b / g)]) / n = (b / g) / [1 - (b / g)]
```
Under the corresponding downstream normalized-coverage signal convention, each accepted fragment contributes total mass one across the genomic bins it overlaps, and the aggregate is divided by the accepted fragment count. Therefore, before optional downstream scaling and rendered rounding, the nonnegative bin masses sum approximately to `1` and form a probability mass function. Their cumulative sum, not the signal itself, is a cumulative distribution function. The sum is approximate because chromosome-boundary clipping can shorten terminal bins and floating-point accumulation can introduce small numerical differences. Optional downstream scaling changes the mass sum, and decimal rendering can further change the displayed sum without changing the pre-rendering values.

Under this convention, expected per-bin depth depends only on `b` and `g`; no input file is needed. This normalized-coverage derivation and its siQ-ChIP provenance come from Brad Dickson's legacy workflow. The `frag` and `norm` modes preserve those calculations for reproducibility but are not recommended for new analyses.

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

Reproduce the legacy normalized floor from explicit dimensions:
```bash
compute_input_floor \
    --mode norm \
    --siz_bin 30 \
    --siz_gen 12157105
```

<br />

## Defaults, output, and failures
The CLI defaults are `mode='dist'`, `method='qntl_nz'`, `qntl_nz=1.0`, `floor=0.0`, `eps=0.0`, `mode_nz='closed'`, `siz_bin=10`, `siz_gen=12157105`, and `dp=24`. The `siz_gen` default is appropriate for *Saccharomyces cerevisiae* when retaining multi-mapping alignments.

`compute_input_floor()` returns one unrounded floating-point `dep_min`. Distribution mode returns the selected statistic after its lower bound; fragment and normalized modes return their bin-to-genome-size formulas. The CLI uses the shared `utilities.utils_format.format_value()` behavior: it emits finite values with at most `dp` decimal places, removes non-informative trailing zeros and a trailing decimal point, and normalizes negative zero to `0`.

The callable raises `InputFloorValidationError` for an invalid mode, invalid `frag`/`norm` dimensions, an invalid input format, an empty filtered positive distribution, or a non-finite computed floor. It raises `AlignmentReadError` when an alignment input cannot be read. CLI validation also requires a finite `qntl_nz` in `[0, 100]`, nonnegative `coef`, `floor`, and `eps` values when applicable, and a nonnegative `dp`. Anticipated validation and computation failures retain their actionable stderr diagnostics and status `1`; parser-controlled usage errors retain argparse behavior.
