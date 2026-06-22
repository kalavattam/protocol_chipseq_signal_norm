
# Python Style Guide
This file owns Python structure rules for this repository. `AGENTS.md` is the top-level routing guide.

<br />

## Script Structure
Python command-line scripts should use an explicit `main()` function and a final entrypoint guard:
```python
def main() -> int:
    ...


if __name__ == "__main__":
    raise SystemExit(main())
```

Keep imports, constants, small helper functions, parser construction, and `main()` in a predictable order. Use `main()` for orchestration, not for large blocks of domain logic.

Utility modules under `scripts/functions/` do not need a `main()` function unless they are directly executable. Keep reusable formatting, parsing, and calculation helpers importable without side effects.

<br />

## Formatting and Output
Prefer clear `snake_case` names and preserve existing short project idioms where they are already established.

Format emitted finite numeric values without useless trailing decimal zeros when that is the existing workflow convention. Keep output formatting helpers shared rather than reimplementing them in multiple scripts.

Keep source code, code comments, and diagnostics near 80 characters where practical. Markdown prose and shell help heredocs have their own rules in `docs/style/help.md` and `docs/style/shell.md`.

<br />

## Docstrings
Use concise docstrings that state what a function does, not a step-by-step restatement of the implementation. For nontrivial helpers, mention important assumptions about input shape, emitted values, or failure behavior.
