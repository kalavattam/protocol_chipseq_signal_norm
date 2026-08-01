#!/usr/bin/env python3
"""Run numeric applicability with a temporary, deterministic inventory."""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

from dev.audit import numeric_emission_inventory, numeric_output_applicability


def main(argv: list[str] | None = None) -> int:
    """Build only an ephemeral inventory, then run the registered checker."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path)
    parser.add_argument("--schema", type=Path)
    parser.add_argument("--inventory", type=Path)
    args = parser.parse_args(argv)
    forwarded: list[str] = ["--root", "."]
    if args.config is not None:
        forwarded.extend(["--config", str(args.config)])
    if args.schema is not None:
        forwarded.extend(["--schema", str(args.schema)])
    if args.inventory is not None:
        forwarded.extend(["--inventory", str(args.inventory)])
        return numeric_output_applicability.main(forwarded)
    with tempfile.TemporaryDirectory(
        prefix="numeric_output_applicability."
    ) as tmp:
        inventory = Path(tmp) / "inventory.json"
        result = numeric_emission_inventory.main(
            ["--root", ".", "--json-out", str(inventory)]
        )
        if result:
            return result
        forwarded.extend(["--inventory", str(inventory)])
        return numeric_output_applicability.main(forwarded)


if __name__ == "__main__":
    raise SystemExit(main())
