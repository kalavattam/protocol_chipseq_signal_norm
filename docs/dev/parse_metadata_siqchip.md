
# siQ-ChIP Metadata Parser Design
## Overview
`parse_metadata_siqchip.py` connects one alignment filename to exactly one row in a siQ-ChIP metadata TSV. It is intentionally deterministic: it parses filename tokens according to YAML, matches only configured metadata fields, and fails on zero or multiple matching rows.

The parser retrieves metadata values only. It does not calculate scaling factors, count alignments, estimate fragment lengths, choose equations, or infer missing measurements.

<br />

## Data model
The parser takes three inputs:
- an alignment filename, such as a BAM or CRAM path;
- a metadata TSV containing sample identifiers and calculator input values; and
- a YAML file describing filename tokens, row-matching fields, and exported calculator inputs.

<br />

## Data flow
The parser has a narrow role in the scaling-factor workflow:
```text
alignment filename + metadata TSV + parser YAML
        ↓
metadata parser
        ↓
one matched metadata row + configured metadata-derived calculator inputs
        ↓
downstream siQ-ChIP wrapper and calculator
        ↓
scaling coefficient
```

The parser only exports values present in the matched metadata row. Some exported values may be precomputed quantities, such as alignment depths or fragment lengths, that allow the workflow to skip downstream calculations. The parser does not compute those values itself; it only retrieves them from the matched metadata row.

<br />

### Important terms
| Term                     | Definition                                                                                                                |
| :----------------------- | :------------------------------------------------------------------------------------------------------------------------ |
| Parsed fields            | Filename-derived values assigned positionally using `filename.fields`.                                                    |
| Match fields             | Parsed fields listed under `matching.fields`; these become the metadata-row match criteria.                               |
| Metadata columns         | Literal column headers in the metadata TSV or CSV table.                                                                  |
| Field-to-column mappings | Explicit mappings from parsed field names to metadata column names when the names differ.                                 |
| Calculator inputs        | Configured metadata columns exported under fixed wrapper-facing keys for downstream scaling-factor calculation.           |
| Precomputed values       | Optional metadata values, such as `dep_*` or `len_*`, that may satisfy workflow inputs before alignment-derived fallback. |

<br />

## YAML scheme
The YAML specifies how the filename is parsed, which parsed fields are used for metadata matching, and which metadata columns are exported for downstream calculators. Parsed fields not listed under `matching.fields` remain filename annotations and do not constrain metadata-row matching.

<br />

### `filename`
The `filename` section strips one alignment extension, strips one configured suffix, splits the basename, and assigns tokens to declared fields.
```yaml
filename:
  delimiter: "_"
  strip_extensions:
    - ".bam"
    - ".cram"
    - ".sam"
  strip_suffixes:
    - ".sc"
    - ".sp"
  fields:
    - assay
    - genotype
    - state
    - factor
    - id
```

<br />

### `matching`
The `matching` section declares which parsed fields identify a metadata row. Each field must exist in `filename.fields`.
```yaml
matching:
  fields:
    - genotype
    - state
    - factor
    - id
```

<br />

### `field_to_column`
The `field_to_column` section maps parsed field names to metadata TSV columns when they differ.
```yaml
field_to_column:
  id: strain  # 'id' in filename maps to 'strain' in metadata TSV
```

<br />

### `calculator_inputs`
The `calculator_inputs` section maps fixed parser-to-wrapper export keys to user-editable metadata TSV column names. The left-side keys are workflow interface names and should not be renamed (unless the shell wrappers are updated too).
```yaml
calculator_inputs:
  siqchip:
    required:
      vol_in: volume_in
      vol_all: volume_all
      mass_in: mass_in
      mass_ip: mass_ip
    optional:
      len_in: length_in
      len_ip: length_ip
      dep_in: dep_in
      dep_ip: dep_ip
      lib_vol_in: lib_vol_in
      lib_vol_ip: lib_vol_ip
```

For custom metadata headers, keep the fixed workflow export keys on the left and change only the mapped TSV column names on the right:
```yaml
calculator_inputs:
  siqchip:
    optional:
      len_ip: IP_fragment_length_bp
      len_in: input_fragment_length_bp
```

<br />

### Takeaways
Users should customize `filename.fields`, `matching.fields`, and `field_to_column` to match their filename scheme and metadata table. In `calculator_inputs`, users should customize only the right-side TSV column names.

<br />

## Processing Rules
For one alignment filename and metadata TSV, the parser:
1. strips one configured alignment-file extension,
2. strips one configured sample suffix,
3. splits the remaining basename using `filename.delimiter`,
4. assigns tokens to `filename.fields`,
5. builds match criteria from `matching.fields`,
6. maps match fields through `field_to_column`,
7. requires exactly one matching TSV row, and
8. exports configured calculator inputs from that row.

The filename must split into exactly the declared number of fields; extra or missing tokens are errors.

Parsed fields not listed under `matching.fields` are retained as filename annotations and are not used to filter the metadata table.

Required calculator inputs must be present and non-missing in the matched row. Optional inputs export as `NA` when absent. Paired optional values must be both present or both absent:
- `dep_ip` and `dep_in`;
- `len_ip` and `len_in`; and
- `lib_vol_ip` and `lib_vol_in`.

<br />

## Examples
Default Tsukiyama-style yeast (*S. cerevisiae*) names use `id` as the filename field and map it to the `strain` metadata column:
```yaml
filename:
  delimiter: "_"
  strip_extensions: [".bam", ".cram", ".sam"]
  strip_suffixes: [".sc", ".sp"]
  fields: [assay, genotype, state, factor, id]

matching:
  fields: [genotype, state, factor, id]

field_to_column:
  id: strain
```

For `IP_WT_Q_Hmo1_7750.sc.cram`, the parser can parse `assay=IP`, `genotype=WT`, `state=Q`, `factor=Hmo1`, and `id=7750`. Because `id` is listed in `matching.fields` and `field_to_column` maps `id` to the `strain` metadata column, the metadata-row match uses `genotype=WT`, `state=Q`, `factor=Hmo1`, and `strain=7750`.

PRJNA857063-style HeLa (*H. sapiens*) names include an `id` token mapped to `replicate`, but do not use it for matching because this study has no biological replication:
```yaml
filename:
  delimiter: "_"
  strip_extensions: [".bam", ".cram", ".sam"]
  strip_suffixes: [".hs", ".dm"]
  fields: [assay, cell, treatment, factor, id]

matching:
  fields: [cell, treatment, factor]

field_to_column:
  id: replicate
```

For `IP_HeLa_DMSO_H3K18ac_rep1.hs.bam`, the parser can parse `assay=IP`, `cell=HeLa`, `treatment=DMSO`, `factor=H3K18ac`, and `id=rep1`. With the PRJNA857063 config, row matching uses `cell`, `treatment`, and `factor`; if `id` were added to `matching.fields`, it would match through `field_to_column` as `replicate=rep1`.

<br />

## Source precedence
The parser does not compute final depth or fragment-length values. It exports metadata-table candidates, and the workflow resolves final values later.

For siQ-ChIP depth values:
```text
dep_ip / dep_in:
  1. command-line depth override, if provided
  2. metadata TSV value, if present
  3. count from alignment file
```

For siQ-ChIP fragment-length values:
```text
len_ip / len_in:
  1. command-line length override, if provided
  2. metadata TSV value, if present
  3. estimate from alignment file
  4. single-end default length, if applicable
```

In brief:
- `mass_*` and `vol_*` must come from metadata.
- `dep_*` may come from metadata or be counted from alignment files.
- `len_*` may come from metadata or be estimated from alignment files.
- `lib_vol_*` is optional paired metadata.

<br />

## Fail-fast rules
The parser and calling workflow intentionally fail rather than guess:
- After extension and suffix stripping, the filename must split into exactly the number of fields declared in `filename.fields`.
- `matching.fields` must be a subset of `filename.fields`.
- Each matching field must correspond to a literal metadata column, either by the same name or through `field_to_column`.
- Required calculator inputs must be present and non-missing in the matched metadata row.
- Matching must return exactly one metadata row.
- Paired optional fields must be supplied together or omitted together.
- Values with source precedence must resolve to exactly one final value.

Unsupported fallbacks are not attempted:
- no fuzzy aliases,
- no wildcard row matching,
- no regex mode,
- no optional-token backtracking, and
- no study-specific fallback behavior.

<br />

## Non-Goals
The parser does not perform fuzzy matching, alias resolution, regex matching, optional-token backtracking, wildcard row matching, or study-specific fallback behavior. Those behaviors were intentionally removed to keep metadata matching predictable, auditable, and easier to maintain.
