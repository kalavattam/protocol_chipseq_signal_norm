
# Python source-policy test fixtures
These tracked, text-only fixtures exercise deterministic portions of the bounded Python source-policy checker. They are intentionally named `*.py.fixture`, so maintained-Python discovery does not treat rejected syntax or source forms as migration targets.

Run the fixture recipe from the repository root in `env_protocol`:
```bash
conda run --no-capture-output -n env_protocol \
    bash tests/fixtures/python_source_policy/make.sh
```

The recipe is nonmutating and idempotent. It validates the tracked inputs in place rather than generating ignored data: canonical, boundary, and literal-exception fixtures must pass, while the negative fixture must fail under every expected deterministic owner.

<br />

## Provenance and roles
The fixtures are repository-authored micro-cases derived from the reconciled Python source-policy workshop and the four-file pilot requirements. They contain no downloaded, private, biological, scheduler, or workflow data.

- `positive.py.fixture` covers canonical docstrings, annotations, comments, strings, definition topology, multiline parameter and construction delimiters, a multiline NumPy type, compact nested tuples, and a greedily wrapped `parse_args()` help literal.
- `negative.py.fixture` contains bounded near misses for every deterministic checker owner, including docstring alignment, a misaligned multiline NumPy type closer, multiline parameters, and an overlong single static help literal.
- `boundary.py.fixture` supplies exact source-form boundary material for parameter markers and generated 79- and 80-column unit cases.
- `exceptions.py.fixture` covers governed source-header and directive exclusions, a necessary raw docstring, embedded-double-quote selection, and the multiline triple-single literal-content exception.
- `make.sh` validates cohort status and expected rule-ID coverage without rewriting any fixture.

The focused unit tests add generated in-memory cases for exact 79-column help and comment boundaries, rendered help equivalence, `self` receiver scope, generator-expression exclusions, X/Y/Z thresholds, candidate-schema stability, and checker or producer idempotence.

<br />

## Deterministic expectations
The positive, boundary, and exception fixtures must produce zero deterministic findings. The negative fixture must produce at least one finding for each of:
- `HELP.PROSE.SENTENCES`
- `PY.CLI.HELP.LAYOUT`
- `PY.COMMENT.FORM`
- `PY.DOCSTRING.LAYOUT`
- `PY.DOCSTRING.NUMPY`
- `PY.NAMING.IDENTIFIERS`
- `PY.SOURCE.LAYOUT`
- `PY.STRING.QUOTES`
- `PY.TYPE.ANNOTATIONS`
- `SOURCE.DELIMITED.MULTILINE`

Semantic readability, docstring usefulness, naming quality, abbreviation clarity, paragraph boundaries, and X/Y/Z candidate dispositions are deliberately absent from the recipe’s pass/fail decision.
